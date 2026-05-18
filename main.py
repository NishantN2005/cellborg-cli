# dropin_cellborg_ui.py
from nicegui import ui, app
import os
import shutil
import asyncio
from pathlib import Path
import json
from datetime import datetime

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from qc_functions import (
    read_10x_mtx,
    find_species,
    calculate_qc_metrics,
    voilin_plot,
    gate_adata,
    scrublet,
    MT_TO_SPECIES,
    SPECIES_TO_MT,
)

from pa_functions import (
    read_adata,
    init_project,
    do_clustering,
    gene_expression,
    annotations,
    run_differential_expression,
)

from heatmap_functions import compute_cluster_gene_heatmap, prepare_heatmap_matrix

PROJECTS_DIR = Path("projects")
SELECTED_PROJECT: dict | None = None


# -----------------------------
# filesystem / metadata helpers
# -----------------------------
def ensure_projects_dir() -> None:
    PROJECTS_DIR.mkdir(parents=True, exist_ok=True)


def get_projects() -> list[Path]:
    ensure_projects_dir()
    return [Path(e.path) for e in os.scandir(PROJECTS_DIR) if e.is_dir()]


def unique_dest_folder(name: str) -> Path:
    base = PROJECTS_DIR / name
    if not base.exists():
        return base
    i = 1
    while True:
        cand = PROJECTS_DIR / f"{name}-{i}"
        if not cand.exists():
            return cand
        i += 1


def load_project_metadata(project_dir: Path) -> dict | None:
    meta_path = project_dir / "cellborg-cli" / "metadata.json"
    if not meta_path.exists():
        return None
    try:
        data = json.loads(meta_path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else None
    except Exception:
        return None


def existing_titles() -> set[str]:
    titles = set()
    for p in get_projects():
        md = load_project_metadata(p)
        if md and md.get("project_title"):
            titles.add(str(md["project_title"]).strip())
    return titles


# -----------------------------
# Annotation helpers
# -----------------------------
def _get_cluster_series_from_adata(adata) -> pd.Series:
    for key in ("leiden", "cluster", "clusters", "louvain"):
        if key in adata.obs.columns:
            return adata.obs[key].astype(str)
    for key in adata.obs.columns:
        lk = key.lower()
        if "cluster" in lk or "leiden" in lk or "louvain" in lk:
            return adata.obs[key].astype(str)
    raise KeyError("Could not find cluster labels in adata.obs (expected 'leiden' or 'cluster').")


def _cluster_sort_key(x: str):
    try:
        return (0, int(x))
    except Exception:
        return (1, str(x))


def _load_existing_annotations(project_path: str | Path) -> dict:
    p = Path(project_path) / "cellborg-cli" / "cluster_annotations.json"
    if not p.exists():
        return {}
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
        return d if isinstance(d, dict) else {}
    except Exception:
        return {}


# ------------------------------
# QC functions
# ------------------------------
def load_violin_arrays(project_path: str | Path):
    hc_path = Path(project_path) / "cellborg-cli" / "highcharts_data.json"
    if not hc_path.exists():
        raise FileNotFoundError(f"{hc_path} not found")

    hc = json.loads(hc_path.read_text(encoding="utf-8"))

    a1, a2, a3 = [], [], []

    if isinstance(hc, dict):
        vals = list(hc.values())

        # Layout 2: {"CELLID": {"n_genes":..,"total_counts":..,"pct_counts_mt":..}, ...}
        if vals and isinstance(vals[0], dict):
            for v in vals:
                if not isinstance(v, dict):
                    continue
                ng = v.get("n_genes")
                tc = v.get("total_counts")
                mt = v.get("pct_counts_mt")
                if ng is not None:
                    try:
                        a1.append(float(ng))
                    except Exception:
                        pass
                if tc is not None:
                    try:
                        a2.append(float(tc))
                    except Exception:
                        pass
                if mt is not None:
                    try:
                        a3.append(float(mt))
                    except Exception:
                        pass

        # Layout 1: {"n_genes_by_counts":[...], "total_counts":[...], "pct_counts_mt":[...]}
        else:
            k1 = "n_genes_by_counts" if "n_genes_by_counts" in hc else "n_genes"
            k2 = "total_counts"
            k3 = "pct_counts_mt"
            for x in hc.get(k1, []) or []:
                try:
                    a1.append(float(x))
                except Exception:
                    pass
            for x in hc.get(k2, []) or []:
                try:
                    a2.append(float(x))
                except Exception:
                    pass
            for x in hc.get(k3, []) or []:
                try:
                    a3.append(float(x))
                except Exception:
                    pass

    return a1, a2, a3


def violin_figure(values, title: str, x_label: str = ""):
    v = np.array(values, dtype=float)
    v = v[np.isfinite(v)]
    fig = go.Figure()
    fig.add_trace(
        go.Violin(
            y=v,
            box_visible=True,
            meanline_visible=True,
            points="outliers",
            name=title,
        )
    )
    fig.update_layout(
        title=title,
        height=320,
        margin=dict(l=10, r=10, t=40, b=10),
    )
    fig.update_xaxes(title=x_label)
    return fig


# ------------------------------
# QC SAVE (FIXED)
# ------------------------------
def _as_float(x, default=0.0) -> float:
    try:
        if x is None:
            return float(default)
        return float(x)
    except Exception:
        return float(default)


def thresholds_elements_to_values(qc_thresholds: dict) -> dict:
    """
    Converts QC_THRESHOLDS (ui.number elements) into plain floats
    so background threads + gate_adata don't see UI objects.
    """
    out: dict = {}
    for metric, bundle in (qc_thresholds or {}).items():
        min_el = bundle.get("min_input")
        max_el = bundle.get("max_input")
        out[metric] = {
            "min": _as_float(getattr(min_el, "value", None)),
            "max": _as_float(getattr(max_el, "value", None)),
        }
    return out


async def save_qc_async(*, adata, qc_thresholds_elements: dict, project_path: str, ui_slot, save_button=None):
    """
    Runs gating + scrublet + write in a background thread (so websocket stays alive),
    then notifies + navigates on the UI thread inside ui_slot.
    """
    try:
        thresholds_values = thresholds_elements_to_values(qc_thresholds_elements)

        def work():
            a = gate_adata(adata, thresholds_values)
            scrublet(a, 0.065)
            out_dir = Path(project_path) / "cellborg-cli"
            out_dir.mkdir(parents=True, exist_ok=True)
            out_file = out_dir / "adata_qc.h5ad"
            a.write(out_file)
            return str(out_file)

        out_file = await asyncio.to_thread(work)

        with ui_slot:
            ui.notify(f"Saved QC dataset → {out_file}", color="positive")
            ui.navigate.to("/")

    except Exception as e:
        with ui_slot:
            ui.notify(f"QC save failed: {e}", color="negative")
        if save_button is not None:
            try:
                save_button.enable()
            except Exception:
                pass


# ------------------------------
# PA: load UMAP JSON + plot (cluster colors)
# ------------------------------
def load_umap_dict(project_path: str | Path, resolution: float) -> dict:
    umap_path = (
        Path(project_path)
        / "cellborg-cli"
        / "figures"
        / "umap"
        / f"umap_clustering_res_{int(resolution * 100)}.json"
    )
    if not umap_path.exists():
        raise FileNotFoundError(f"UMAP json not found: {umap_path}")
    return json.loads(umap_path.read_text(encoding="utf-8"))


def umap_scatter_figure_from_dict(umap_dict: dict, title="UMAP (Leiden clusters)", annotations: dict | None = None):
    """
    umap_dict format:
    {
        "CELLID": {"UMAP1": float, "UMAP2": float, "cluster": "0"},
        ...
    }
    One trace per cluster => Plotly assigns different colors automatically.
    
    annotations: optional dict mapping cluster IDs to cell type names
    """
    df = pd.DataFrame.from_dict(umap_dict, orient="index")

    df["UMAP1"] = pd.to_numeric(df.get("UMAP1"), errors="coerce")
    df["UMAP2"] = pd.to_numeric(df.get("UMAP2"), errors="coerce")
    df["cluster"] = df.get("cluster").astype(str)

    df = df.dropna(subset=["UMAP1", "UMAP2"])

    clusters = sorted(df["cluster"].unique(), key=_cluster_sort_key)

    # Build cluster display names from annotations
    cluster_display_names = {}
    for c in clusters:
        if annotations and str(c) in annotations and annotations[str(c)]:
            cluster_display_names[c] = f"{c}: {annotations[str(c)]}"
        else:
            cluster_display_names[c] = c

    fig = go.Figure()
    for c in clusters:
        sub = df[df["cluster"] == c]
        fig.add_trace(
            go.Scattergl(
                x=sub["UMAP1"],
                y=sub["UMAP2"],
                mode="markers",
                name=cluster_display_names[c],
                marker=dict(size=3),
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        title=title,
        height=700,
        margin=dict(l=20, r=20, t=50, b=20),
        legend=dict(itemsizing="constant"),
    )
    fig.update_xaxes(title="UMAP1")
    fig.update_yaxes(title="UMAP2")
    return fig


# -----------------------------
# native folder picker
# -----------------------------
async def choose_folder_via_native_file_dialog() -> Path | None:
    if not getattr(app, "native", None) or app.native.main_window is None:
        ui.notify("Native window not available. Run with ui.run(native=True).", color="negative")
        return None

    try:
        import webview  # pip install pywebview
    except Exception as e:
        ui.notify(f"pywebview not installed/available: {e}", color="negative")
        return None

    dialog_type = webview.FileDialog.OPEN if hasattr(webview, "FileDialog") else webview.OPEN_DIALOG

    files = await app.native.main_window.create_file_dialog(
        dialog_type=dialog_type,
        allow_multiple=False,
        file_types=("All files (*.*)",),
    )
    if not files:
        return None

    file_path = Path(files[0])
    if not file_path.exists():
        ui.notify(f"File not found: {file_path}", color="negative")
        return None

    folder = file_path.parent
    if not folder.is_dir():
        ui.notify(f"Not a folder: {folder}", color="negative")
        return None

    return folder


async def copy_folder_to_projects(src_folder: Path) -> Path:
    ensure_projects_dir()
    dest = unique_dest_folder(src_folder.name)

    def do_copy():
        shutil.copytree(str(src_folder), str(dest))

    await asyncio.to_thread(do_copy)
    return dest


# -----------------------------
# UI callbacks
# -----------------------------
def select_project(md: dict) -> None:
    global SELECTED_PROJECT
    SELECTED_PROJECT = md
    project_details.refresh()


@ui.refreshable
def projects_list():
    for project in get_projects():
        md = load_project_metadata(project)
        if not md:
            ui.button(f"{project.name} (no metadata)", on_click=lambda p=project: ui.notify(str(p)))
            continue

        title = str(md.get("project_title") or project.name).strip() or project.name
        ui.button(title, on_click=lambda m=md: select_project(m))


@ui.refreshable
def project_details():
    ui.label("Project Details").style("font-size: 28px; font-weight: 700; color: #4ecda4; margin: 6px 0;")

    if not SELECTED_PROJECT:
        ui.label("Select a project to see details here.").style("color: #555;")
        return

    ui.label(f"Title: {SELECTED_PROJECT.get('project_title', '')}").style("color:#ddd; font-size: 18px;")
    ui.label(f"Description: {SELECTED_PROJECT.get('project_description', '')}").style("color:#bbb;")

    status = str(SELECTED_PROJECT.get("status", "unknown")).upper()
    with ui.row().classes("w-full").style("gap:12px; margin-top:12px;"):
        if status == "QC":
            ui.button("Run QC", on_click=lambda: ui.navigate.to("/qc")).style("flex:1;")
        elif status in ("PROC_ANNO", "PA"):
            ui.button("Run PA", on_click=lambda: ui.navigate.to("/pca")).style("flex:1;")
        ui.button("Run Analysis", on_click=lambda: ui.navigate.to("/analysis")).style("flex:1;")


async def add_project() -> None:
    ui.notify("Select any file inside the folder you want to add.", color="info")
    folder = await choose_folder_via_native_file_dialog()
    if folder is None:
        return

    try:
        dest = await copy_folder_to_projects(folder)
        species = find_species(dest / "features.tsv.gz")
    except Exception as e:
        ui.notify(f"Failed to copy: {e}", color="negative")
        return

    dialog.clear()

    with dialog, ui.card().classes("w-96"):
        ui.label("Project Title")
        title_input = ui.input(
            validation=lambda v: "Title already exists" if v and v.strip() in existing_titles() else None
        ).props("autofocus")

        ui.label("Project description")
        desc_input = ui.textarea().style("height: 100px; width: 100%;")

        def save():
            md = {
                "original_folder_name": folder.name,
                "project_path": str(dest),
                "species": species,
                "project_title": (title_input.value or dest.name).strip(),
                "project_description": desc_input.value or "",
                "added_at": datetime.now().isoformat(),
                "status": "QC",
            }
            (dest / "cellborg-cli").mkdir(parents=True, exist_ok=True)
            (dest / "cellborg-cli" / "metadata.json").write_text(json.dumps(md, indent=4), encoding="utf-8")

            ui.notify("Saved project metadata", color="positive")
            dialog.close()
            projects_list.refresh()

        with ui.row().classes("justify-end w-full"):
            ui.button("Cancel", on_click=dialog.close)
            ui.button("Save", on_click=save)

    dialog.open()


# -----------------------------
# Pages
# -----------------------------
@ui.page("/")
def dashboard_page():
    ensure_projects_dir()

    with ui.dialog() as d:
        global dialog
        dialog = d

    ui.dark_mode().enable()

    with ui.row().style(
        """
        padding:12px;
        width:100%;
        height:95vh;
        box-sizing:border-box;
        gap:12px;
        overflow:hidden;
        """
    ):
        with ui.column().style(
            """
            width:520px;
            height:100%;
            padding:20px;
            border:1px solid #e5e7eb;
            border-radius:8px;
            box-shadow:0 1px 3px rgba(0,0,0,0.06);
            overflow-y:auto;
            """
        ):
            ui.label("Projects").style("font-size: 36px; font-weight: 800; color: #4ecda4;")
            ui.button("Add New Project", on_click=add_project)
            projects_list()

        with ui.column().style(
            """
            flex:1;
            height:100%;
            padding:20px;
            border:1px solid #e5e7eb;
            border-radius:8px;
            box-shadow:0 1px 3px rgba(0,0,0,0.06);
            overflow-y:auto;
            """
        ):
            project_details()


@ui.page("/qc")
def qc_page():
    SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path") if SELECTED_PROJECT else None
    QC_THRESHOLDS = {}

    if not SELECTED_PROJECT_PATH:
        ui.label("No project selected. Go back and select one.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return

    ui.dark_mode().enable()

    ui.label("QC Runner").style("font-size: 32px; font-weight: 800; color: #4ecda4; margin: 12px 0;")
    ui.label(f"Project: {SELECTED_PROJECT.get('project_title', '')}").style("color:#ddd;")
    ui.label(SELECTED_PROJECT.get("project_description", "")).style("color:#bbb;")

    ui.separator().style("opacity:0.25; margin: 12px 0;")

    try:
        adata = read_10x_mtx(SELECTED_PROJECT_PATH)
    except Exception as e:
        ui.notify(f"Failed to read project data: {e}", color="negative")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return

    mt = "MT-"
    try:
        species = SELECTED_PROJECT.get("species")
        if species:
            mt = SPECIES_TO_MT.get(species, mt)
    except Exception:
        pass

    try:
        calculate_qc_metrics(mt, adata)
        voilin_plot(adata, SELECTED_PROJECT_PATH)
    except Exception as e:
        ui.notify(f"QC metrics failed: {e}", color="negative")
        return

    try:
        a1, a2, a3 = load_violin_arrays(SELECTED_PROJECT_PATH)
    except Exception as e:
        ui.notify(f"Could not load violin plot data: {e}", color="negative")
        return

    def finite_min_max(arr):
        v = np.array(arr, dtype=float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            return 0.0, 1.0
        lo, hi = float(v.min()), float(v.max())
        if hi <= lo:
            hi = lo + 1.0
        return lo, hi

    def make_violin_with_lines(values, title, y_min, y_max):
        v = np.array(values, dtype=float)
        v = v[np.isfinite(v)]

        fig = go.Figure()
        fig.add_trace(
            go.Violin(
                y=v,
                box_visible=True,
                meanline_visible=True,
                points="outliers",
                name=title,
            )
        )

        fig.update_layout(
            title=dict(text=title, x=0.5),
            height=520,
            margin=dict(l=25, r=25, t=55, b=20),
            xaxis=dict(showticklabels=False),
            dragmode="pan",
            newshape=dict(line=dict(width=2)),
            shapes=[
                dict(
                    type="line",
                    xref="paper",
                    x0=0,
                    x1=1,
                    yref="y",
                    y0=y_min,
                    y1=y_min,
                    line=dict(width=3),
                    name="min_line",
                ),
                dict(
                    type="line",
                    xref="paper",
                    x0=0,
                    x1=1,
                    yref="y",
                    y0=y_max,
                    y1=y_max,
                    line=dict(width=3),
                    name="max_line",
                ),
            ],
        )
        return fig

    def set_line(fig, idx: int, y: float):
        fig.layout.shapes[idx].y0 = y
        fig.layout.shapes[idx].y1 = y

    def qc_tile(metric: str, values: list[float]):
        data_min, data_max = finite_min_max(values)
        y_min = data_min
        y_max = data_max

        with ui.card().style("width: 560px; padding: 16px;"):
            fig = make_violin_with_lines(values, metric, y_min, y_max)
            plot = ui.plotly(fig).style("width:100%;")

            with ui.row().classes("w-full items-center").style("gap: 12px; margin-top: 10px;"):
                min_in = ui.number("Min line", value=y_min, format="%.3f").style("flex:1;")
                max_in = ui.number("Max line", value=y_max, format="%.3f").style("flex:1;")

            QC_THRESHOLDS[metric] = {
                "min_input": min_in,
                "max_input": max_in,
                "plot": plot,
            }

            def apply_from_inputs():
                nonlocal y_min, y_max
                y_min = float(min_in.value)
                y_max = float(max_in.value)

                if y_min > y_max:
                    y_min = y_max
                    min_in.value = y_min

                set_line(plot.figure, 0, y_min)
                set_line(plot.figure, 1, y_max)

                plot.figure.update_yaxes(range=[y_min, y_max])
                plot.update()

            min_in.on("change", lambda e: apply_from_inputs())
            max_in.on("change", lambda e: apply_from_inputs())

            def on_relayout(e):
                payload = getattr(e, "args", None)
                if not isinstance(payload, dict):
                    return

                new_min = None
                new_max = None

                for k, v in payload.items():
                    if k in ("shapes[0].y0", "shapes[0].y1"):
                        try:
                            new_min = float(v)
                        except Exception:
                            pass
                    if k in ("shapes[1].y0", "shapes[1].y1"):
                        try:
                            new_max = float(v)
                        except Exception:
                            pass

                changed = False
                if new_min is not None:
                    min_in.value = new_min
                    changed = True
                if new_max is not None:
                    max_in.value = new_max
                    changed = True

                if changed:
                    apply_from_inputs()

            plot.on("plotly_relayout", on_relayout)
            plot.on("plotly_relayouting", on_relayout)

    with ui.row().classes("w-full justify-center").style("gap:24px; padding-top:16px; flex-wrap:wrap;"):
        qc_tile("n_genes", a1)
        qc_tile("total_counts", a2)
        qc_tile("pct_counts_mt", a3)

    ui.separator().style("opacity:0.25; margin: 12px 0;")

    actions = ui.row().classes("w-full").style("gap:12px;")
    with actions:
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat").style("flex:1")
        save_btn = ui.button("Save").props("flat").style("flex:1")

        def on_save_click():
            save_btn.disable()
            ui.notify("Saving QC…", color="info")

            # update metadata status (kept as you had it)
            md_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "metadata.json"
            if md_path.exists():
                try:
                    md = json.loads(md_path.read_text(encoding="utf-8"))
                    md["status"] = "PROC_ANNO"
                    md_path.write_text(json.dumps(md, indent=4), encoding="utf-8")
                    SELECTED_PROJECT.update(md)
                except Exception as e:
                    print(f"Failed to update metadata.json: {e}")

            asyncio.create_task(
                save_qc_async(
                    adata=adata,
                    qc_thresholds_elements=QC_THRESHOLDS,
                    project_path=SELECTED_PROJECT_PATH,
                    ui_slot=actions,
                    save_button=save_btn,
                )
            )

        save_btn.on("click", lambda e: on_save_click())


@ui.page("/pca")
def pca_page():
    PCA_DONE = False
    SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path") if SELECTED_PROJECT else None

    if not SELECTED_PROJECT_PATH:
        ui.label("No project selected. Go back and select one.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return

    ui.dark_mode().enable()

    ui.label("Processing and Annotation").style("font-size: 32px; font-weight: 800; color: #4ecda4; margin: 12px 0;")
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    adata_qc_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_qc.h5ad"
    if not adata_qc_path.exists():
        ui.label("QC dataset not found. Please run QC first.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return

    adata = read_adata(adata_qc_path)
    init_project(SELECTED_PROJECT_PATH, adata)

    with ui.row().classes("w-full items-center").style("gap:16px;"):
        ui.label("Leiden resolution").style("color:#ddd; min-width:160px;")
        res_slider = ui.slider(min=0.1, max=2.0, value=0.8, step=0.05).props("label-always").style("flex:1;")
        res_label = ui.label("0.80").style("color:#ddd; width:70px; text-align:right;")
        run_btn = ui.button("Run clustering").style("min-width:170px;")

    def sync_res_label():
        res_label.text = f"{float(res_slider.value):.2f}"

    res_slider.on("change", lambda e: sync_res_label())
    sync_res_label()

    ui.separator().style("opacity:0.25; margin: 12px 0;")

    with ui.card().classes("w-full").style("padding:16px;"):
        status_lbl = ui.label("Adjust resolution and click “Run clustering” to render UMAP.").style("color:#bbb;")
        umap_plot = ui.plotly(go.Figure()).style("width:100%;")

    next_btn = ui.button("Next →", on_click=lambda: ui.navigate.to("/anno")).props("unelevated")
    next_btn.disable()

    async def run_and_render():
        nonlocal PCA_DONE
        PCA_DONE = False
        run_btn.disable()
        next_btn.disable()
        res = float(res_slider.value)
        status_lbl.text = f"Clustering at resolution {res:.2f}…"
        try:
            await asyncio.to_thread(do_clustering, adata, res)

            umap_dict = load_umap_dict(SELECTED_PROJECT_PATH, res)
            # Load annotations for display
            try:
                annotations = _load_existing_annotations(SELECTED_PROJECT_PATH)
            except Exception:
                annotations = None
            fig = umap_scatter_figure_from_dict(umap_dict, title=f"UMAP (Leiden res={res:.2f})", annotations=annotations)

            umap_plot.figure = fig
            umap_plot.update()
            status_lbl.text = f"Rendered UMAP for resolution {res:.2f}"
            PCA_DONE = True
            next_btn.enable()
        except Exception as e:
            status_lbl.text = f"Failed: {e}"
            ui.notify(f"PA failed: {e}", color="negative")
        finally:
            run_btn.enable()

    run_btn.on("click", lambda e: asyncio.create_task(run_and_render()))

    with ui.row().classes("w-full justify-end mt-6"):
        next_btn


@ui.page("/anno")
def annotation_page():
    ui.dark_mode().enable()

    SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path") if SELECTED_PROJECT else None
    if not SELECTED_PROJECT_PATH:
        ui.label("No project selected. Go back and select one.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return

    pv_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "project_values.json"
    if not pv_path.exists():
        ui.label("project_values.json not found. Run PA first.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return

    try:
        pv = json.loads(pv_path.read_text(encoding="utf-8"))
        gene_list = pv.get("gene_list", [])
        if not isinstance(gene_list, list):
            gene_list = []
        gene_list = [str(g) for g in gene_list]
    except Exception as e:
        ui.label(f"Failed to read project_values.json: {e}").style("color:#fbbf24;")
        return

    genes_upper = [(g, g.upper()) for g in gene_list]

    adata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_clustered.h5ad"
    if not adata_path.exists():
        ui.label("Clustered dataset not found. Please run PA first.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return
    adata = read_adata(adata_path)

    ui.label("Annotations").style("font-size: 32px; font-weight: 800; color: #4ecda4; margin: 12px 0;")
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    with ui.row().classes("w-full").style("gap:18px; padding:12px; height: calc(100vh - 120px);"):
        # LEFT COLUMN
        with ui.column().style("width:520px; gap:10px; position:relative;"):
            ui.label("Gene search").style("color:#e5e7eb; font-weight:700;")

            chips_row = ui.row().classes("items-center w-full").style(
                "gap:10px; padding:10px; border:1px solid #374151; border-radius:8px; background:#0b1220;"
            )

            selected_genes: list[str] = []

            def render_chips():
                chips_row.clear()
                with chips_row:
                    if not selected_genes:
                        ui.label("Selected genes will appear here.").style("color:#9ca3af;")
                        return
                    for g in selected_genes:
                        with ui.row().classes("items-center").style(
                            "gap:8px; padding:7px 10px; border:1px solid #4b5563; border-radius:8px; background:#111827;"
                        ):
                            ui.label(g).style("color:#e5e7eb; font-size:16px;")
                            ui.button(
                                icon="close",
                                on_click=lambda gg=g: (selected_genes.remove(gg), render_chips()),
                            ).props("flat dense").style("color:#ef4444;")

            search_wrap = ui.column().classes("w-full").style("position:relative;")
            with search_wrap:
                search_in = (
                    ui.input(placeholder="Search gene (e.g., Dnpep, Rnpepl1, Sept2)")
                    .props("outlined clearable")
                    .style("width:100%;")
                )

                dropdown = ui.card().classes("w-full").style(
                    """
                    position:absolute;
                    top:52px;
                    left:0;
                    right:0;
                    z-index:9999;
                    padding:6px;
                    border:1px solid #374151;
                    border-radius:10px;
                    background:#0b1220;
                    box-shadow: 0 8px 20px rgba(0,0,0,0.35);
                    max-height: 260px;
                    overflow-y: auto;
                    display:none;
                    """
                )

            def hide_dropdown():
                dropdown.style("display:none;")

            def show_dropdown():
                dropdown.style("display:block;")

            def add_gene(g: str):
                if g not in selected_genes:
                    selected_genes.append(g)
                    render_chips()
                search_in.value = ""
                hide_dropdown()
                ui.notify(f"Selected {g}", color="positive")

            def render_dropdown(matches: list[str], query: str):
                dropdown.clear()
                q = (query or "").strip()
                if not q:
                    hide_dropdown()
                    return

                if not matches:
                    with dropdown:
                        ui.label(f'No matches for "{q}".').style("color:#9ca3af; padding:6px;")
                    show_dropdown()
                    return

                with dropdown:
                    for g in matches[:20]:
                        ui.button(
                            g,
                            on_click=lambda gg=g: add_gene(gg),
                        ).props("flat dense").style(
                            "width:100%; justify-content:flex-start; color:#e5e7eb; text-transform:none;"
                        )
                show_dropdown()

            def do_search(query: str):
                q = (query or "").strip()
                if not q:
                    render_dropdown([], "")
                    return
                q_up = q.upper()
                matches = [g for (g, gu) in genes_upper if q_up in gu]
                render_dropdown(matches, q)

            search_in.on("input", lambda e: do_search(search_in.value))
            search_in.on("change", lambda e: do_search(search_in.value))

            ui.separator().style("opacity:0.25; margin: 10px 0;")
            ui.label("Plot type").style("color:#e5e7eb; font-weight:700;")

            plot_type = ui.radio(
                options=["Feature Plot", "Violin Plot", "Dot Plot"],
                value="Feature Plot",
            ).props("inline").style("color:#e5e7eb;")

            load_btn = ui.button("Load Plot").style(
                "width:100%; height:46px; background:#41c99a; color:#0b1220; font-weight:800; border-radius:8px;"
            )

            render_chips()

            # ---- Annotate clusters section ----
            ui.separator().style("opacity:0.25; margin: 10px 0;")
            ui.label("Annotate Clusters").style("color:#e5e7eb; font-weight:800;")

            try:
                cluster_series = _get_cluster_series_from_adata(adata)
                cluster_ids = sorted(cluster_series.unique().tolist(), key=_cluster_sort_key)
            except Exception as e:
                cluster_ids = []
                ui.label(f"Could not load clusters: {e}").style("color:#fbbf24;")

            existing_ann = _load_existing_annotations(SELECTED_PROJECT_PATH)

            annot_card = ui.card().classes("w-full").style(
                "padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
            )

            cluster_inputs: dict[str, any] = {}

            with annot_card:
                if not cluster_ids:
                    ui.label("No clusters found. Run PA first.").style("color:#9ca3af;")
                else:
                    with ui.column().style("max-height: 280px; overflow-y: auto; padding-right:6px; gap:10px;"):
                        for cid in cluster_ids:
                            with ui.row().classes("w-full items-center").style("gap:12px;"):
                                ui.label(f"Cluster {cid}:").style("color:#e5e7eb; width:120px;")
                                inp = ui.input(
                                    value=str(existing_ann.get(str(cid), "")),
                                    placeholder="e.g., T cells, Fibroblasts, ..."
                                ).props("outlined dense").style("flex:1;")
                                cluster_inputs[str(cid)] = inp

            with ui.row().classes("w-full").style("gap:12px; margin-top:10px;"):
                save_ann_btn = ui.button("Save").style(
                    "flex:1; height:44px; border-radius:8px; background:#41c99a; color:#0b1220; font-weight:900;"
                )
                next_page_btn = ui.button("Next Page").style(
                    "flex:1; height:44px; border-radius:8px; background:#111827; color:#e5e7eb; font-weight:900;"
                )

            def _save_cluster_annotations():
                if not cluster_inputs:
                    ui.notify("No cluster inputs to save.", color="negative")
                    return
                try:
                    annotations_dict: dict[str, str] = {}
                    for cid, inp in cluster_inputs.items():
                        annotations_dict[str(cid)] = (inp.value or "").strip()

                    # persist mapping for UI reload
                    ann_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "cluster_annotations.json"
                    ann_path.parent.mkdir(parents=True, exist_ok=True)
                    ann_path.write_text(json.dumps(annotations_dict, indent=2), encoding="utf-8")

                    # call your backend writer (umap_annotations.json etc.)
                    annotations(adata, annotations_dict)

                    # save updated adata
                    adata_out = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_annotated.h5ad"
                    adata.write(adata_out)

                    ui.notify("Saved annotations successfully", color="positive")
                except Exception as e:
                    ui.notify(f"Annotation save failed: {e}", color="negative")

            save_ann_btn.on("click", lambda e: _save_cluster_annotations())

            def finish_annotations():
                metadata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "metadata.json"
                if metadata_path.exists():
                    try:
                        md = json.loads(metadata_path.read_text(encoding="utf-8"))
                        md["status"] = "ANNO_DONE"
                        metadata_path.write_text(json.dumps(md, indent=4), encoding="utf-8")
                        SELECTED_PROJECT.update(md)
                    except Exception as e:
                        print(f"Failed to update metadata.json: {e}")
                ui.navigate.to("/")
            next_page_btn.on("click", lambda e: finish_annotations())
            hide_dropdown()

        # RIGHT COLUMN
        with ui.column().style("flex:1; gap:12px;"):
            ui.label("Single-cell RNAseq Gene Expression").style(
                "font-size: 20px; font-weight: 800; color:#e5e7eb; text-align:center; margin-top:6px;"
            )

            plot_card = ui.card().classes("w-full").style(
                "flex:1; padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
            )
            with plot_card:
                plot_status = ui.label("Choose a plot type and click Load Plot.").style("color:#9ca3af;")
                plot_el = ui.plotly(go.Figure()).style("width:100%; height:100%;")

            with ui.row().classes("w-full").style("gap:12px;"):
                cluster_btn = ui.button("Cluster Plot").style(
                    "flex:1; height:48px; border-radius:8px; background:#41c99a; color:#0b1220; font-weight:800;"
                )

                plot_name_btn = ui.button(plot_type.value).style(
                    "flex:1; height:48px; border-radius:8px; background:#41c99a; color:#0b1220; font-weight:800;"
                )

            async def load_cluster_plot():
                with plot_card:
                    cluster_btn.disable()
                    try:
                        res_raw = pv.get("clust_resolution", None)
                        if res_raw is None:
                            plot_status.text = 'Missing "clust_resolution" in project_values.json'
                            ui.notify('Missing "clust_resolution" in project_values.json', color="negative")
                            return

                        try:
                            res = float(res_raw)
                        except Exception:
                            plot_status.text = f'Invalid clust_resolution: {res_raw!r}'
                            ui.notify(f'Invalid "clust_resolution": {res_raw!r}', color="negative")
                            return

                        plot_status.text = f"Loading cluster UMAP (res={res:.2f})…"
                        umap_dict = load_umap_dict(SELECTED_PROJECT_PATH, res)
                        # Load annotations for display
                        try:
                            annotations = _load_existing_annotations(SELECTED_PROJECT_PATH)
                        except Exception:
                            annotations = None
                        fig = umap_scatter_figure_from_dict(umap_dict, title=f"UMAP (Leiden res={res:.2f})", annotations=annotations)

                        plot_el.figure = fig
                        plot_el.update()

                        plot_status.text = f"Loaded cluster UMAP (res={res:.2f})"
                        ui.notify("Loaded cluster plot", color="positive")
                    except Exception as e:
                        plot_status.text = f"Failed to load cluster plot: {e}"
                        ui.notify(f"Cluster plot failed: {e}", color="negative")
                    finally:
                        cluster_btn.enable()

            # IMPORTANT: pass async function directly (do NOT wrap in lambda)
            cluster_btn.on("click", load_cluster_plot)

            def load_gene_expression_df(project_path: str | Path) -> pd.DataFrame:
                p = Path(project_path) / "cellborg-cli" / "figures" / "gene_expression_per_cell_with_clusters.json"
                if not p.exists():
                    raise FileNotFoundError(f"Gene expression json not found: {p}")
                data = json.loads(p.read_text(encoding="utf-8"))
                df = pd.DataFrame.from_dict(data, orient="index")
                df["UMAP1"] = pd.to_numeric(df.get("UMAP1"), errors="coerce")
                df["UMAP2"] = pd.to_numeric(df.get("UMAP2"), errors="coerce")
                df["cluster"] = df.get("cluster").astype(str)
                df = df.dropna(subset=["UMAP1", "UMAP2", "cluster"])
                return df

            def feature_plot_multi_gene(df: pd.DataFrame, genes: list[str]) -> go.Figure:
                genes = [g for g in genes if g in df.columns]
                if not genes:
                    raise KeyError("None of the selected genes exist in the gene-expression JSON.")

                for g in genes:
                    df[g] = pd.to_numeric(df[g], errors="coerce").fillna(0.0)

                clusters = sorted(df["cluster"].unique(), key=_cluster_sort_key)
                fig = go.Figure()

                traces_per_gene = len(clusters) + 1  # clusters + colorbar trace

                for gi, gene in enumerate(genes):
                    for c in clusters:
                        sub = df[df["cluster"] == c]
                        fig.add_trace(go.Scattergl(
                            x=sub["UMAP1"],
                            y=sub["UMAP2"],
                            mode="markers",
                            name=f"Cluster {c}",
                            marker=dict(size=3, color=sub[gene], colorscale="Viridis", showscale=False),
                            hovertemplate=f"cluster={c}<br>{gene}=%{{marker.color:.3f}}<extra></extra>",
                            visible=(gi == 0),
                        ))

                    gmin = float(df[gene].min())
                    gmax = float(df[gene].max())
                    fig.add_trace(go.Scattergl(
                        x=[None], y=[None],
                        mode="markers",
                        marker=dict(
                            color=[gmin, gmax],
                            colorscale="Viridis",
                            showscale=True,
                            colorbar=dict(title=f"{gene} expr"),
                        ),
                        hoverinfo="skip",
                        showlegend=False,
                        visible=(gi == 0),
                    ))

                buttons = []
                for gi, gene in enumerate(genes):
                    vis = [False] * (traces_per_gene * len(genes))
                    start = gi * traces_per_gene
                    for k in range(traces_per_gene):
                        vis[start + k] = True
                    buttons.append(dict(
                        label=gene,
                        method="update",
                        args=[{"visible": vis}, {"title": f"Feature Plot: {gene} (per cluster)"}],
                    ))

                fig.update_layout(
                    title=f"Feature Plot: {genes[0]} (per cluster)",
                    height=650,
                    margin=dict(l=20, r=20, t=70, b=20),
                    legend=dict(itemsizing="constant"),
                    updatemenus=[dict(
                        type="dropdown",
                        x=0.01, y=1.15,
                        xanchor="left", yanchor="top",
                        buttons=buttons,
                    )],
                )
                fig.update_xaxes(title="UMAP1")
                fig.update_yaxes(title="UMAP2")
                return fig

            def violin_plot_multi_gene(df: pd.DataFrame, genes: list[str]) -> go.Figure:
                genes = [g for g in genes if g in df.columns]
                if not genes:
                    raise KeyError("None of the selected genes exist in the gene-expression JSON.")

                for g in genes:
                    df[g] = pd.to_numeric(df[g], errors="coerce").fillna(0.0)

                clusters = sorted(df["cluster"].unique(), key=_cluster_sort_key)
                fig = go.Figure()

                traces_per_gene = len(clusters)

                for gi, gene in enumerate(genes):
                    for c in clusters:
                        sub = df[df["cluster"] == c]
                        fig.add_trace(go.Violin(
                            y=sub[gene],
                            name=f"{c}",
                            box_visible=True,
                            meanline_visible=True,
                            points=False,
                            visible=(gi == 0),
                        ))

                buttons = []
                for gi, gene in enumerate(genes):
                    vis = [False] * (traces_per_gene * len(genes))
                    start = gi * traces_per_gene
                    for k in range(traces_per_gene):
                        vis[start + k] = True
                    buttons.append(dict(
                        label=gene,
                        method="update",
                        args=[{"visible": vis}, {"title": f"Violin: {gene} by cluster"}],
                    ))

                fig.update_layout(
                    title=f"Violin: {genes[0]} by cluster",
                    height=650,
                    margin=dict(l=20, r=20, t=70, b=20),
                    xaxis_title="Cluster",
                    yaxis_title="Expression",
                    violinmode="group",
                    updatemenus=[dict(
                        type="dropdown",
                        x=0.01, y=1.15,
                        xanchor="left", yanchor="top",
                        buttons=buttons,
                    )],
                )
                return fig

            def dot_plot_per_cluster(df: pd.DataFrame, genes: list[str]) -> go.Figure:
                genes = [g for g in genes if g in df.columns]
                if not genes:
                    raise KeyError("None of the selected genes exist in the gene-expression JSON.")

                clusters = sorted(df["cluster"].unique(), key=_cluster_sort_key)

                rows = []
                for c in clusters:
                    sub = df[df["cluster"] == c]
                    for g in genes:
                        v = pd.to_numeric(sub[g], errors="coerce").fillna(0.0)
                        mean_expr = float(v.mean())
                        pct_expr = float((v > 0).mean())
                        rows.append((c, g, mean_expr, pct_expr))

                stats = pd.DataFrame(rows, columns=["cluster", "gene", "mean_expr", "pct_expr"])
                x = stats["gene"].tolist()
                y = stats["cluster"].tolist()
                color = stats["mean_expr"].to_numpy()
                size = (stats["pct_expr"].to_numpy() * 30.0) + 3.0

                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=x,
                    y=y,
                    mode="markers",
                    marker=dict(
                        size=size,
                        color=color,
                        colorscale="Viridis",
                        showscale=True,
                        colorbar=dict(title="Mean expr"),
                        sizemode="diameter",
                    ),
                    hovertemplate="cluster=%{y}<br>gene=%{x}<br>mean=%{marker.color:.3f}<extra></extra>",
                    showlegend=False,
                ))

                fig.update_layout(
                    title="Dot Plot (per cluster)",
                    height=650,
                    margin=dict(l=20, r=20, t=50, b=20),
                    xaxis_title="Gene",
                    yaxis_title="Cluster",
                )
                return fig

            def load_plot():
                gene_expression(SELECTED_PROJECT_PATH,adata, selected_genes)

                pt = plot_type.value
                genes_label = ", ".join(selected_genes[:5]) if selected_genes else "(no genes selected)"
                plot_status.text = f"Loaded: {pt} | Genes: {genes_label}"
                ui.notify(f"Loading {pt}", color="info")

                try:
                    df = load_gene_expression_df(SELECTED_PROJECT_PATH)

                    if pt == "Feature Plot":
                        if not selected_genes:
                            raise ValueError("Select at least one gene for Feature Plot.")
                        fig = feature_plot_multi_gene(df, selected_genes)

                    elif pt == "Violin Plot":
                        if not selected_genes:
                            raise ValueError("Select at least one gene for Violin Plot.")
                        fig = violin_plot_multi_gene(df, selected_genes)

                    else:
                        if not selected_genes:
                            raise ValueError("Select at least one gene for Dot Plot.")
                        fig = dot_plot_per_cluster(df, selected_genes[:12])

                    plot_el.figure = fig
                    plot_el.update()

                except Exception as e:
                    plot_status.text = f"Failed to load plot: {e}"
                    ui.notify(f"Plot failed: {e}", color="negative")
                    return

                plot_name_btn.text = pt

            def on_plot_type_change():
                plot_name_btn.text = plot_type.value

            plot_type.on("change", lambda e: on_plot_type_change())
            on_plot_type_change()

            load_btn.on("click", lambda e: load_plot())

@ui.page("/analysis")
def analysis_page():
    ui.dark_mode().enable()

    ui.label("Analysis").style(
        "font-size: 32px; font-weight: 800; color: #4ecda4; margin: 12px 0;"
    )
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    # --- helper to build a "pop-up" tile card ---
    def analysis_tile(title: str, subtitle: str, icon_name: str, accent: str, on_click=None):
        # nice hover pop + glow
        base = (
            "w-full h-full p-6 rounded-2xl "
            "bg-gray-900/60 border border-white/10 "
            "transition-all duration-200 ease-out "
            "hover:-translate-y-1 hover:scale-[1.01] hover:border-white/20 "
            "hover:shadow-2xl cursor-pointer select-none"
        )

        card = ui.card().classes(base)
        if on_click:
            card.on("click", on_click)

        with card:
            # icon (upper-left)
            ui.icon(icon_name).style(
                f"font-size: 28px; color: {accent}; margin-bottom: 12px;"
            )

            # title + subtitle
            ui.label(title).style(
                "font-size: 22px; font-weight: 800; color: #e5e7eb; line-height: 1.1;"
            )
            ui.label(subtitle).style(
                "margin-top: 8px; color: #94a3b8; font-size: 14px;"
            )

            # little accent bar at bottom
            ui.separator().style(f"opacity:0.25; margin: 16px 0 0 0;")
            ui.element("div").style(
                f"height: 6px; border-radius: 999px; background: {accent}; opacity: 0.9; margin-top: 14px;"
            )

        return card

    # --- layout: 2x2 big tiles (responsive) ---
    with ui.element("div").classes("w-full"):
        with ui.element("div").classes("grid grid-cols-1 md:grid-cols-2 gap-6"):
            analysis_tile(
                "Heatmap",
                "Cluster × gene (or module) expression overview.",
                "grid_on",
                "#4ecda4",
                on_click=lambda: ui.navigate.to("/heatmap"),
            )

            analysis_tile(
                "Pseudotime",
                "Trajectory inference and lineage progression plots.",
                "timeline",
                "#60a5fa",
                on_click=lambda: ui.navigate.to("/pseudotime"),
            )

            analysis_tile(
                "Feature Expression",
                "Marker / gene expression across embeddings and groups.",
                "scatter_plot",
                "#fbbf24",
                on_click=lambda: ui.navigate.to("/feature"),
            )

            analysis_tile(
                "Receptor–Ligand",
                "Cell–cell communication and interaction scoring.",
                "link",
                "#fb7185",
                on_click=lambda: ui.navigate.to("/receptor_ligand"),
            )

            analysis_tile(
                "Differential Expression",
                "Marker genes per cluster via Wilcoxon or t-test.",
                "bar_chart",
                "#a78bfa",
                on_click=lambda: ui.navigate.to("/de"),
            )

@ui.page("/heatmap")
def heatmap_page():
    ui.dark_mode().enable()

    ui.label("Heatmap").style(
        "font-size: 32px; font-weight: 800; color: #4ecda4; margin: 12px 0;"
    )
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    # ----------------------------
    # Resolve project + load gene list
    # ----------------------------
    SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path") if SELECTED_PROJECT else None
    if not SELECTED_PROJECT_PATH:
        ui.label("No project selected. Go back and select one.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    pv_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "project_values.json"
    if not pv_path.exists():
        ui.label("project_values.json not found. Run PA first.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    try:
        pv = json.loads(pv_path.read_text(encoding="utf-8"))
        gene_list = pv.get("gene_list", [])
        if not isinstance(gene_list, list):
            gene_list = []
        gene_list = [str(g) for g in gene_list]
    except Exception as e:
        ui.label(f"Failed to read project_values.json: {e}").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    genes_upper = [(g, g.upper()) for g in gene_list]

    # Load clustered adata (needed for heatmap compute)
    adata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_annotated.h5ad"
    if not adata_path.exists():
        ui.label("adata_annotated.h5ad not found. Run PA first.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    adata = read_adata(adata_path)

    # ----------------------------
    # Gene search UI (same as /anno)
    # ----------------------------
    ui.label("Gene search").style("color:#e5e7eb; font-weight:700; margin-top: 6px;")

    chips_row = ui.row().classes("items-center w-full").style(
        "gap:10px; padding:10px; border:1px solid #374151; border-radius:8px; background:#0b1220;"
    )
    selected_genes: list[str] = []

    def render_chips():
        chips_row.clear()
        with chips_row:
            if not selected_genes:
                ui.label("Selected genes will appear here.").style("color:#9ca3af;")
                return
            for g in selected_genes:
                with ui.row().classes("items-center").style(
                    "gap:8px; padding:7px 10px; border:1px solid #4b5563; border-radius:8px; background:#111827;"
                ):
                    ui.label(g).style("color:#e5e7eb; font-size:16px;")
                    ui.button(
                        icon="close",
                        on_click=lambda gg=g: (selected_genes.remove(gg), render_chips()),
                    ).props("flat dense").style("color:#ef4444;")

    search_wrap = ui.column().classes("w-full").style("position:relative; max-width: 820px;")
    with search_wrap:
        search_in = (
            ui.input(placeholder="Search gene (e.g., Dnpep, Rnpepl1, Sept2)")
            .props("outlined clearable")
            .style("width:100%;")
        )

        dropdown = ui.card().classes("w-full").style(
            """
            position:absolute;
            top:52px;
            left:0;
            right:0;
            z-index:9999;
            padding:6px;
            border:1px solid #374151;
            border-radius:10px;
            background:#0b1220;
            box-shadow: 0 8px 20px rgba(0,0,0,0.35);
            max-height: 260px;
            overflow-y: auto;
            display:none;
            """
        )

    def hide_dropdown():
        dropdown.style("display:none;")

    def show_dropdown():
        dropdown.style("display:block;")

    def add_gene(g: str):
        if g not in selected_genes:
            selected_genes.append(g)
            render_chips()
        search_in.value = ""
        hide_dropdown()
        ui.notify(f"Selected {g}", color="positive")

    def render_dropdown(matches: list[str], query: str):
        dropdown.clear()
        q = (query or "").strip()
        if not q:
            hide_dropdown()
            return
        if not matches:
            with dropdown:
                ui.label(f'No matches for "{q}".').style("color:#9ca3af; padding:6px;")
            show_dropdown()
            return

        with dropdown:
            for g in matches[:20]:
                ui.button(g, on_click=lambda gg=g: add_gene(gg)).props("flat dense").style(
                    "width:100%; justify-content:flex-start; color:#e5e7eb; text-transform:none;"
                )
        show_dropdown()

    def do_search(query: str):
        q = (query or "").strip()
        if not q:
            render_dropdown([], "")
            return
        q_up = q.upper()
        matches = [g for (g, gu) in genes_upper if q_up in gu]
        render_dropdown(matches, q)

    search_in.on("input", lambda e: do_search(search_in.value))
    search_in.on("change", lambda e: do_search(search_in.value))

    render_chips()
    hide_dropdown()

    # ----------------------------
    # Controls (metric/scale/threshold) + Render button
    # ----------------------------
    ui.separator().style("opacity:0.25; margin: 14px 0;")

    with ui.row().classes("w-full items-center").style("gap:14px; flex-wrap:wrap;"):
        metric_sel = ui.select(
            options={
                "mean": "Mean expression",
                "pct": "% expressing",
            },
            value="mean",
            label="Metric",
        ).props("outlined dense").style("min-width:220px;")

        scale_sel = ui.select(
            options={
                "raw": "Raw",
                "zscore": "Z-score (per gene)",
            },
            value="raw",
            label="Scale",
        ).props("outlined dense").style("min-width:240px;")



        expr_thr = ui.number("Expr threshold (for % expressing)", value=0.0, format="%.3f").props("outlined dense").style(
            "min-width:260px;"
        )

        render_btn = ui.button("Render heatmap").props("unelevated").style(
            "height:40px; background:#41c99a; color:#0b1220; font-weight:900; border-radius:10px;"
        )

    # ----------------------------
    # Plot area
    # ----------------------------
    plot_card = ui.card().classes("w-full").style(
        "margin-top: 12px; padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
    )
    with plot_card:
        status = ui.label("Select genes and click Render heatmap.").style("color:#9ca3af;")
        heatmap_el = ui.plotly(go.Figure()).style("width:100%; height: 720px;")

    async def render_heatmap():
        render_btn.disable()
        with plot_card:
            try:
                if not selected_genes:
                    ui.notify("Select at least one gene.", color="negative")
                    return

                status.text = "Computing heatmap…"
                metric = str(metric_sel.value)
                scale = str(scale_sel.value)
                thr = float(expr_thr.value or 0.0)

                # Heavy compute in background thread
                def _work():
                    res = compute_cluster_gene_heatmap(
                        adata,
                        selected_genes,
                        cluster_key="leiden",     # change if needed
                        use_layer=None,           # or "log1p" if you store it
                        use_raw=False,            # set True if you want adata.raw
                        expr_threshold=thr,
                    )
                    mat = prepare_heatmap_matrix(res, metric=metric, scale=scale)
                    return res, mat

                res, mat = await asyncio.to_thread(_work)

                if mat.shape[1] == 0:
                    status.text = "No selected genes were found in this dataset."
                    ui.notify("None of the selected genes exist in adata.var_names.", color="negative")
                    return

                # Load cluster annotations and map cluster IDs to names
                cluster_labels = mat.index.tolist()
                try:
                    existing_ann = _load_existing_annotations(SELECTED_PROJECT_PATH)
                    # Map cluster ID -> annotation (or keep ID if no annotation)
                    cluster_labels = [
                        f"{cid}: {existing_ann.get(str(cid), '')}" if existing_ann.get(str(cid)) else cid
                        for cid in cluster_labels
                    ]
                except Exception:
                    pass  # If annotation loading fails, use raw cluster IDs

                # Plotly heatmap
                fig = go.Figure()
                fig.add_trace(
                    go.Heatmap(
                        z=mat.to_numpy(dtype=float),
                        x=mat.columns.tolist(),
                        y=cluster_labels,
                        hovertemplate="Cluster=%{y}<br>Gene=%{x}<br>Value=%{z:.4f}<extra></extra>",
                        colorbar=dict(title=("Z" if scale == "zscore" else metric)),
                    )
                )
                fig.update_layout(
                    title=f"Heatmap ({'Mean expr' if metric=='mean' else '% expressing'}; {scale})",
                    height=700,
                    margin=dict(l=60, r=20, t=50, b=80),
                    xaxis=dict(tickangle=45),
                )

                heatmap_el.figure = fig
                heatmap_el.update()

                missing = [g for g in selected_genes if g not in res.genes]
                if missing:
                    ui.notify(f"Missing genes (not found): {', '.join(missing[:8])}" + ("…" if len(missing) > 8 else ""), color="warning")

                status.text = f"Rendered heatmap for {len(res.genes)} genes across {len(res.clusters)} clusters."
                ui.notify("Heatmap rendered", color="positive")

            except Exception as e:
                status.text = f"Failed: {e}"
                ui.notify(f"Heatmap failed: {e}", color="negative")
            finally:
                render_btn.enable()

    render_btn.on("click", lambda e: asyncio.create_task(render_heatmap()))

    ui.separator().style("opacity:0.25; margin: 14px 0;")
    ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")

@ui.page("/pseudotime")
def pseudotime_page():
    ui.dark_mode().enable()

    ui.label("Pseudotime").style(
        "font-size: 32px; font-weight: 800; color: #60a5fa; margin: 12px 0;"
    )
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path") if SELECTED_PROJECT else None
    if not SELECTED_PROJECT_PATH:
        ui.label("No project selected. Go back and select one.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    adata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_annotated.h5ad"
    if not adata_path.exists():
        ui.label("adata_annotated.h5ad not found. Run PA/Annotations first.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    # Load adata
    try:
        adata = read_adata(adata_path)
    except Exception as e:
        ui.notify(f"Failed to load adata: {e}", color="negative")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    # Cluster key
    cluster_key = "leiden"
    if cluster_key not in getattr(adata, "obs", {}):
        # fallback attempt
        try:
            _ = _get_cluster_series_from_adata(adata)
            # _get_cluster_series_from_adata returns a series but we still want a key name;
            # if leiden missing, just pick first matching key
            for k in adata.obs.columns:
                lk = k.lower()
                if "leiden" in lk or "cluster" in lk or "louvain" in lk:
                    cluster_key = k
                    break
        except Exception:
            cluster_key = None

    # Root cluster options
    try:
        if cluster_key:
            cluster_ids = sorted(adata.obs[cluster_key].astype(str).unique().tolist(), key=_cluster_sort_key)
        else:
            cluster_ids = []
    except Exception:
        cluster_ids = []

    # Layout
    with ui.row().classes("w-full").style("gap:18px; padding:12px; height: calc(100vh - 140px);"):
        # LEFT: controls
        with ui.column().style("width:520px; gap:12px;"):
            ui.label("Controls").style("color:#e5e7eb; font-weight:800; font-size:18px;")

            ctrl_card = ui.card().classes("w-full").style(
                "padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
            )
            with ctrl_card:
                ui.label("Root selection").style("color:#9ca3af; font-weight:700;")

                root_mode = ui.radio(
                    options=["Root cluster", "Root cell id"],
                    value="Root cluster",
                ).props("inline").style("color:#e5e7eb;")

                root_cluster_sel = ui.select(
                    options=cluster_ids,
                    value=(cluster_ids[0] if cluster_ids else None),
                    label=f"Root cluster ({cluster_key or 'missing'})",
                ).props("outlined dense").style("width:100%;")

                root_cell_in = ui.input(
                    label="Root cell id (index in adata.obs_names)",
                    placeholder="Paste a cell id exactly (e.g., AAACCTG...)",
                ).props("outlined clearable dense").style("width:100%;")

                ui.separator().style("opacity:0.25; margin: 8px 0;")

                neighbors_k = ui.number("Neighbors (k)", value=15, format="%.0f").props("outlined dense").style("width:100%;")
                n_dcs = ui.number("Diffusion components", value=10, format="%.0f").props("outlined dense").style("width:100%;")

                ui.separator().style("opacity:0.25; margin: 8px 0;")

                color_mode = ui.select(
                    options={
                        "pseudotime": "Pseudotime",
                        "cluster": "Cluster labels",
                    },
                    value="pseudotime",
                    label="Color by",
                ).props("outlined dense").style("width:100%;")

                run_btn = ui.button("Compute pseudotime").props("unelevated").style(
                    "width:100%; height:44px; border-radius:10px; background:#60a5fa; color:#0b1220; font-weight:900;"
                )

                ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat").style("width:100%;")

            # enable/disable controls based on root mode
            def sync_root_mode():
                if root_mode.value == "Root cluster":
                    root_cluster_sel.enable()
                    root_cell_in.disable()
                else:
                    root_cluster_sel.disable()
                    root_cell_in.enable()

            root_mode.on("change", lambda e: sync_root_mode())
            sync_root_mode()

        # RIGHT: plot
        with ui.column().style("flex:1; gap:12px;"):
            ui.label("Trajectory view").style(
                "font-size: 20px; font-weight: 800; color:#e5e7eb; text-align:center; margin-top:6px;"
            )

            plot_card = ui.card().classes("w-full").style(
                "flex:1; padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
            )
            with plot_card:
                status = ui.label("Pick a root and click Compute pseudotime.").style("color:#9ca3af;")
                plot_el = ui.plotly(go.Figure()).style("width:100%; height: 100%;")

    def _plot_from_adata():
        # needs UMAP coords
        if "X_umap" in getattr(adata, "obsm", {}):
            X = np.array(adata.obsm["X_umap"])
        elif "umap" in getattr(adata, "obsm", {}):
            X = np.array(adata.obsm["umap"])
        else:
            raise KeyError("UMAP coordinates not found in adata.obsm (expected 'X_umap'). Run UMAP first in PA pipeline.")

        df = pd.DataFrame({
            "UMAP1": X[:, 0],
            "UMAP2": X[:, 1],
        }, index=adata.obs_names.astype(str))

        # color
        if color_mode.value == "cluster":
            if not cluster_key:
                raise KeyError("No cluster key found in adata.obs to color by cluster.")
            df["color"] = adata.obs[cluster_key].astype(str).values
            # one trace per cluster for categorical colors
            clusters = sorted(df["color"].unique().tolist(), key=_cluster_sort_key)
            
            # Load cluster annotations for display names
            try:
                existing_ann = _load_existing_annotations(SELECTED_PROJECT_PATH)
                cluster_display_names = {
                    cid: f"{cid}: {existing_ann.get(str(cid), '')}" if existing_ann.get(str(cid)) else cid
                    for cid in clusters
                }
            except Exception:
                cluster_display_names = {cid: cid for cid in clusters}
            
            fig = go.Figure()
            for c in clusters:
                sub = df[df["color"] == c]
                fig.add_trace(go.Scattergl(
                    x=sub["UMAP1"], y=sub["UMAP2"],
                    mode="markers",
                    name=cluster_display_names.get(c, c),
                    marker=dict(size=3),
                    hoverinfo="skip",
                ))
            title = f"Pseudotime view (colored by {cluster_key})"
        else:
            if "dpt_pseudotime" not in adata.obs.columns:
                raise KeyError("dpt_pseudotime not found. Compute pseudotime first.")
            pt = pd.to_numeric(adata.obs["dpt_pseudotime"], errors="coerce").fillna(0.0).to_numpy()
            df["pt"] = pt
            fig = go.Figure()
            fig.add_trace(go.Scattergl(
                x=df["UMAP1"], y=df["UMAP2"],
                mode="markers",
                marker=dict(
                    size=3,
                    color=df["pt"],
                    colorscale="Viridis",
                    showscale=True,
                    colorbar=dict(title="DPT"),
                ),
                hovertemplate="pt=%{marker.color:.4f}<extra></extra>",
                showlegend=False,
            ))
            title = "Pseudotime (Scanpy DPT)"

        fig.update_layout(
            title=title,
            height=760,
            margin=dict(l=20, r=20, t=55, b=20),
        )
        fig.update_xaxes(title="UMAP1")
        fig.update_yaxes(title="UMAP2")
        return fig

    async def _compute_pseudotime():
        run_btn.disable()
        try:
            # Import scanpy lazily
            try:
                import scanpy as sc
            except Exception as e:
                status.text = f"Scanpy not available: {e}"
                ui.notify("Scanpy is required for DPT pseudotime. Install: pip install scanpy", color="negative")
                return

            if cluster_key is None and root_mode.value == "Root cluster":
                ui.notify("No clusters found in adata.obs to use as a root cluster.", color="negative")
                return

            status.text = "Computing neighbors / diffusion map / DPT…"

            k = int(float(neighbors_k.value or 15))
            ndc = int(float(n_dcs.value or 10))

            def _work():
                # neighbors on X_pca if exists else default
                sc.pp.neighbors(adata, n_neighbors=max(2, k), use_rep=("X_pca" if "X_pca" in adata.obsm else None))
                sc.tl.diffmap(adata, n_comps=max(2, ndc))

                # set root
                if root_mode.value == "Root cell id":
                    rid = (root_cell_in.value or "").strip()
                    if not rid:
                        raise ValueError("Root cell id is empty.")
                    if rid not in adata.obs_names:
                        raise KeyError(f"Root cell id not found in adata.obs_names: {rid}")
                    adata.uns["iroot"] = int(np.where(adata.obs_names == rid)[0][0])
                else:
                    rc = str(root_cluster_sel.value)
                    # choose first cell in that cluster as root
                    mask = (adata.obs[cluster_key].astype(str) == rc).to_numpy()
                    idxs = np.where(mask)[0]
                    if idxs.size == 0:
                        raise ValueError(f"No cells found in root cluster {rc}.")
                    adata.uns["iroot"] = int(idxs[0])

                sc.tl.dpt(adata)  # creates adata.obs['dpt_pseudotime']
                return True

            await asyncio.to_thread(_work)

            # plot
            fig = _plot_from_adata()
            plot_el.figure = fig
            plot_el.update()
            status.text = "Done. Pseudotime computed and rendered."
            ui.notify("Pseudotime rendered", color="positive")

        except Exception as e:
            status.text = f"Failed: {e}"
            ui.notify(f"Pseudotime failed: {e}", color="negative")
        finally:
            run_btn.enable()

    run_btn.on("click", _compute_pseudotime)

@ui.page("/feature")
def feature_expression_page():
    ui.dark_mode().enable()

    ui.label("Feature Expression").style(
        "font-size: 32px; font-weight: 800; color: #fbbf24; margin: 12px 0;"
    )
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    # ----------------------------
    # Resolve project + load gene list
    # ----------------------------
    SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path") if SELECTED_PROJECT else None
    if not SELECTED_PROJECT_PATH:
        ui.label("No project selected. Go back and select one.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    pv_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "project_values.json"
    if not pv_path.exists():
        ui.label("project_values.json not found. Run PA first.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    try:
        pv = json.loads(pv_path.read_text(encoding="utf-8"))
        gene_list = pv.get("gene_list", [])
        if not isinstance(gene_list, list):
            gene_list = []
        gene_list = [str(g) for g in gene_list]
    except Exception as e:
        ui.label(f"Failed to read project_values.json: {e}").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    genes_upper = [(g, g.upper()) for g in gene_list]

    # Optional: load adata (mostly for cluster key detection / sanity)
    adata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_annotated.h5ad"
    if not adata_path.exists():
        ui.label("adata_annotated.h5ad not found. Run PA/Annotations first.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return
    adata = read_adata(adata_path)

    # ----------------------------
    # Load gene expression DF from JSON written by your backend
    # ----------------------------
    def load_gene_expression_df(project_path: str | Path) -> pd.DataFrame:
        print('------------------ project path ------------------', project_path)
        p = Path(project_path) / "cellborg-cli" / "figures" / "gene_expression_per_cell_with_clusters.json"
        if not p.exists():
            raise FileNotFoundError(f"Gene expression json not found: {p}")
        data = json.loads(p.read_text(encoding="utf-8"))
        df = pd.DataFrame.from_dict(data, orient="index")
        df["UMAP1"] = pd.to_numeric(df.get("UMAP1"), errors="coerce")
        df["UMAP2"] = pd.to_numeric(df.get("UMAP2"), errors="coerce")
        df["cluster"] = df.get("cluster").astype(str)
        df = df.dropna(subset=["UMAP1", "UMAP2", "cluster"])
        return df

    # ----------------------------
    # Gene search UI (same style as /anno & /heatmap)
    # ----------------------------
    ui.label("Gene search").style("color:#e5e7eb; font-weight:700; margin-top: 6px;")

    chips_row = ui.row().classes("items-center w-full").style(
        "gap:10px; padding:10px; border:1px solid #374151; border-radius:8px; background:#0b1220;"
    )
    selected_genes: list[str] = []

    def render_chips():
        chips_row.clear()
        with chips_row:
            if not selected_genes:
                ui.label("Selected genes will appear here.").style("color:#9ca3af;")
                return
            for g in selected_genes:
                with ui.row().classes("items-center").style(
                    "gap:8px; padding:7px 10px; border:1px solid #4b5563; border-radius:8px; background:#111827;"
                ):
                    ui.label(g).style("color:#e5e7eb; font-size:16px;")
                    ui.button(
                        icon="close",
                        on_click=lambda gg=g: (selected_genes.remove(gg), render_chips()),
                    ).props("flat dense").style("color:#ef4444;")

    search_wrap = ui.column().classes("w-full").style("position:relative; max-width: 820px;")
    with search_wrap:
        search_in = (
            ui.input(placeholder="Search gene (e.g., Dnpep, Rnpepl1, Sept2)")
            .props("outlined clearable")
            .style("width:100%;")
        )

        dropdown = ui.card().classes("w-full").style(
            """
            position:absolute;
            top:52px;
            left:0;
            right:0;
            z-index:9999;
            padding:6px;
            border:1px solid #374151;
            border-radius:10px;
            background:#0b1220;
            box-shadow: 0 8px 20px rgba(0,0,0,0.35);
            max-height: 260px;
            overflow-y: auto;
            display:none;
            """
        )

    def hide_dropdown():
        dropdown.style("display:none;")

    def show_dropdown():
        dropdown.style("display:block;")

    def add_gene(g: str):
        if g not in selected_genes:
            selected_genes.append(g)
            render_chips()
        search_in.value = ""
        hide_dropdown()
        ui.notify(f"Selected {g}", color="positive")

    def render_dropdown(matches: list[str], query: str):
        dropdown.clear()
        q = (query or "").strip()
        if not q:
            hide_dropdown()
            return
        if not matches:
            with dropdown:
                ui.label(f'No matches for "{q}".').style("color:#9ca3af; padding:6px;")
            show_dropdown()
            return

        with dropdown:
            for g in matches[:20]:
                ui.button(g, on_click=lambda gg=g: add_gene(gg)).props("flat dense").style(
                    "width:100%; justify-content:flex-start; color:#e5e7eb; text-transform:none;"
                )
        show_dropdown()

    def do_search(query: str):
        q = (query or "").strip()
        if not q:
            render_dropdown([], "")
            return
        q_up = q.upper()
        matches = [g for (g, gu) in genes_upper if q_up in gu]
        render_dropdown(matches, q)

    search_in.on("input", lambda e: do_search(search_in.value))
    search_in.on("change", lambda e: do_search(search_in.value))

    render_chips()
    hide_dropdown()

    # ----------------------------
    # Controls
    # ----------------------------
    ui.separator().style("opacity:0.25; margin: 14px 0;")

    with ui.row().classes("w-full items-center").style("gap:14px; flex-wrap:wrap;"):
        # behavior toggles
        split_by_cluster = ui.switch("Split by cluster", value=False).style("color:#e5e7eb;")
        log1p_toggle = ui.switch("log1p color", value=False).style("color:#e5e7eb;")

        pt_size = ui.number("Point size", value=3, format="%.0f").props("outlined dense").style("width:180px;")
        max_genes = ui.number("Max genes (dropdown)", value=25, format="%.0f").props("outlined dense").style("width:220px;")

        render_btn = ui.button("Render feature plot").props("unelevated").style(
            "height:40px; background:#fbbf24; color:#0b1220; font-weight:900; border-radius:10px;"
        )

    # ----------------------------
    # Plot area
    # ----------------------------
    plot_card = ui.card().classes("w-full").style(
        "margin-top: 12px; padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
    )
    with plot_card:
        status = ui.label("Select genes and click Render feature plot.").style("color:#9ca3af;")
        plot_el = ui.plotly(go.Figure()).style("width:100%; height: 760px;")

    # ----------------------------
    # Plot builder
    # ----------------------------
    def build_feature_plot(df: pd.DataFrame, genes: list[str]) -> go.Figure:
        genes = [g for g in genes if g in df.columns]
        if not genes:
            raise KeyError("None of the selected genes exist in the gene-expression JSON.")

        # sanitize genes to numeric
        for g in genes:
            df[g] = pd.to_numeric(df[g], errors="coerce").fillna(0.0)

        clusters = sorted(df["cluster"].unique().tolist(), key=_cluster_sort_key)
        
        # Load cluster annotations for display names
        try:
            existing_ann = _load_existing_annotations(SELECTED_PROJECT_PATH)
            cluster_display_names = {
                cid: f"{cid}: {existing_ann.get(str(cid), '')}" if existing_ann.get(str(cid)) else cid
                for cid in clusters
            }
        except Exception:
            cluster_display_names = {cid: cid for cid in clusters}
        
        size = int(float(pt_size.value or 3))
        do_log = bool(log1p_toggle.value)
        split = bool(split_by_cluster.value)

        # one "gene group" = (traces + optional colorbar trace)
        # If split_by_cluster: traces per gene = len(clusters) + 1(colorbar)
        # Else: traces per gene = 1 + 1(colorbar)
        if split:
            traces_per_gene = len(clusters) + 1
        else:
            traces_per_gene = 1 + 1

        fig = go.Figure()

        # helper to compute color array
        def color_arr(series: pd.Series) -> np.ndarray:
            v = series.to_numpy(dtype=float)
            if do_log:
                v = np.log1p(np.clip(v, 0.0, None))
            return v

        for gi, gene in enumerate(genes):
            if split:
                for c in clusters:
                    sub = df[df["cluster"] == c]
                    fig.add_trace(go.Scattergl(
                        x=sub["UMAP1"],
                        y=sub["UMAP2"],
                        mode="markers",
                        name=cluster_display_names.get(c, c),
                        marker=dict(
                            size=size,
                            color=color_arr(sub[gene]),
                            colorscale="Viridis",
                            showscale=False,
                        ),
                        hovertemplate=f"cluster={c}<br>{gene}=%{{marker.color:.4f}}<extra></extra>",
                        visible=(gi == 0),
                    ))
            else:
                fig.add_trace(go.Scattergl(
                    x=df["UMAP1"],
                    y=df["UMAP2"],
                    mode="markers",
                    name=gene,
                    marker=dict(
                        size=size,
                        color=color_arr(df[gene]),
                        colorscale="Viridis",
                        showscale=False,
                    ),
                    hovertemplate=f"{gene}=%{{marker.color:.4f}}<extra></extra>",
                    showlegend=False,
                    visible=(gi == 0),
                ))

            # add a shared colorbar for this gene via dummy trace
            vals = color_arr(df[gene])
            gmin = float(np.nanmin(vals)) if vals.size else 0.0
            gmax = float(np.nanmax(vals)) if vals.size else 1.0
            fig.add_trace(go.Scattergl(
                x=[None], y=[None],
                mode="markers",
                marker=dict(
                    color=[gmin, gmax],
                    colorscale="Viridis",
                    showscale=True,
                    colorbar=dict(title=f"{gene} ({'log1p' if do_log else 'raw'})"),
                ),
                hoverinfo="skip",
                showlegend=False,
                visible=(gi == 0),
            ))

        # dropdown to switch gene
        buttons = []
        for gi, gene in enumerate(genes):
            vis = [False] * (traces_per_gene * len(genes))
            start = gi * traces_per_gene
            for k in range(traces_per_gene):
                vis[start + k] = True
            buttons.append(dict(
                label=gene,
                method="update",
                args=[
                    {"visible": vis},
                    {"title": f"Feature Expression: {gene} ({'split clusters' if split else 'all cells'})"},
                ],
            ))

        fig.update_layout(
            title=f"Feature Expression: {genes[0]} ({'split clusters' if split else 'all cells'})",
            height=760,
            margin=dict(l=20, r=20, t=70, b=20),
            legend=dict(itemsizing="constant"),
            updatemenus=[dict(
                type="dropdown",
                x=0.01, y=1.15,
                xanchor="left", yanchor="top",
                buttons=buttons,
            )] if len(genes) > 1 else [],
        )
        fig.update_xaxes(title="UMAP1")
        fig.update_yaxes(title="UMAP2")
        return fig

    # ----------------------------
    # Render handler (slot-safe: no create_task)
    # ----------------------------
    async def render_feature(e=None):
        with plot_card:
            render_btn.disable()
            try:
                if not selected_genes:
                    ui.notify("Select at least one gene.", color="negative")
                    return

                # cap genes to keep plot responsive
                cap = int(float(max_genes.value or 25))
                genes = selected_genes[: max(1, cap)]

                status.text = "Computing gene expression + loading JSON…"

                # Ensure the JSON exists/updated (your backend)
                # This can be heavy, so run in a thread.
                def _work():
                    gene_expression(SELECTED_PROJECT_PATH, adata, genes)
                    print('inside _work', SELECTED_PROJECT_PATH)
                    return load_gene_expression_df(SELECTED_PROJECT_PATH)

                df = await asyncio.to_thread(_work)

                status.text = "Rendering…"
                fig = build_feature_plot(df, genes)

                plot_el.figure = fig
                plot_el.update()

                status.text = f"Rendered feature plot for {len([g for g in genes if g in df.columns])} gene(s)."
                ui.notify("Feature plot rendered", color="positive")

            except Exception as ex:
                status.text = f"Failed: {ex}"
                ui.notify(f"Feature plot failed: {ex}", color="negative")
            finally:
                render_btn.enable()

    render_btn.on("click", render_feature)

    ui.separator().style("opacity:0.25; margin: 14px 0;")
    ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")

@ui.page("/de")
def de_page():
    ui.dark_mode().enable()

    ui.label("Differential Expression").style(
        "font-size: 32px; font-weight: 800; color: #a78bfa; margin: 12px 0;"
    )
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path") if SELECTED_PROJECT else None
    if not SELECTED_PROJECT_PATH:
        ui.label("No project selected. Go back and select one.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    adata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_annotated.h5ad"
    if not adata_path.exists():
        adata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_clustered.h5ad"
    if not adata_path.exists():
        ui.label("Clustered dataset not found. Run PA/Annotations first.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    try:
        adata = read_adata(adata_path)
    except Exception as e:
        ui.notify(f"Failed to load adata: {e}", color="negative")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    # Mutable dict — mutated in-place so closures always see the latest data
    de_results: dict = {}
    de_json_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "de_results.json"
    if de_json_path.exists():
        try:
            de_results.update(json.loads(de_json_path.read_text(encoding="utf-8")))
        except Exception:
            pass

    with ui.row().classes("w-full").style("gap:18px; padding:12px; height: calc(100vh - 140px);"):
        # LEFT: controls
        with ui.column().style("width:520px; gap:12px;"):
            ui.label("Controls").style("color:#e5e7eb; font-weight:800; font-size:18px;")

            ctrl = ui.card().classes("w-full").style(
                "padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
            )
            with ctrl:
                method_sel = ui.select(
                    options={"wilcoxon": "Wilcoxon rank-sum", "t-test": "t-test"},
                    value="wilcoxon",
                    label="Test method",
                ).props("outlined dense").style("width:100%;")

                n_genes_in = ui.number("Top N genes per cluster", value=25, format="%.0f").props("outlined dense").style("width:100%;")

                run_btn = ui.button("Run DE").props("unelevated").style(
                    "width:100%; height:44px; border-radius:10px; background:#a78bfa; color:#0b1220; font-weight:900; margin-top:8px;"
                )

                ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat").style("width:100%;")

        # RIGHT: results
        with ui.column().style("flex:1; gap:12px;"):
            ui.label("Marker genes per cluster").style(
                "font-size: 20px; font-weight: 800; color:#e5e7eb; text-align:center; margin-top:6px;"
            )

            out_card = ui.card().classes("w-full").style(
                "flex:1; padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
            )
            with out_card:
                init_msg = "Previous results loaded — click Run DE to recompute." if de_results else "Click Run DE to compute marker genes."
                status = ui.label(init_msg).style("color:#9ca3af;")

                controls_row = ui.row().classes("w-full items-center").style("gap:12px; margin-bottom:10px;")
                plot_el = ui.plotly(go.Figure()).style("width:100%; height: 440px;")

                table_el = ui.table(
                    columns=[
                        {"name": "rank", "label": "#", "field": "rank", "align": "left", "sortable": True},
                        {"name": "gene", "label": "Gene", "field": "gene", "sortable": True},
                        {"name": "score", "label": "Score", "field": "score", "sortable": True},
                        {"name": "logfc", "label": "Log2FC", "field": "logfc", "sortable": True},
                        {"name": "pval_adj", "label": "Adj. p-val", "field": "pval_adj", "sortable": True},
                    ],
                    rows=[],
                    row_key="rank",
                ).classes("w-full").style("margin-top:10px;")

                export_btn = ui.button("Export CSV").props("flat").style("margin-top:8px; color:#a78bfa;")

    def _build_bar_chart(cluster: str) -> go.Figure:
        rows = de_results.get(cluster, [])
        if not rows:
            return go.Figure()
        genes = [r["gene"] for r in rows]
        scores = [r["score"] for r in rows]
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=scores[::-1],
            y=genes[::-1],
            orientation="h",
            marker=dict(color="#a78bfa"),
        ))
        fig.update_layout(
            title=f"Top marker genes — Cluster {cluster}",
            height=440,
            margin=dict(l=10, r=20, t=50, b=20),
            xaxis_title="Score",
        )
        return fig

    def _update_table(cluster: str):
        rows = de_results.get(cluster, [])
        table_el.rows = [
            {
                "rank": i + 1,
                "gene": r["gene"],
                "score": f"{r['score']:.4f}",
                "logfc": f"{r['logfc']:.4f}",
                "pval_adj": f"{r['pval_adj']:.2e}",
            }
            for i, r in enumerate(rows)
        ]
        table_el.update()

    def _render_for_cluster(cluster: str):
        plot_el.figure = _build_bar_chart(cluster)
        plot_el.update()
        _update_table(cluster)

    def _setup_results_ui():
        clusters = sorted(de_results.keys(), key=_cluster_sort_key)
        if not clusters:
            return

        controls_row.clear()
        with controls_row:
            cluster_sel = ui.select(
                options=clusters,
                value=clusters[0],
                label="View cluster",
            ).props("outlined dense").style("min-width:200px;")
            cluster_sel.on("change", lambda e: _render_for_cluster(str(cluster_sel.value)))

        _render_for_cluster(clusters[0])

    if de_results:
        _setup_results_ui()

    def _export_csv():
        if not de_results:
            ui.notify("Run DE first.", color="negative")
            return
        rows = []
        for cluster, genes in de_results.items():
            for r in genes:
                rows.append({"cluster": cluster, **r})
        csv_bytes = pd.DataFrame(rows).to_csv(index=False).encode("utf-8")
        ui.download(csv_bytes, filename="de_results.csv")

    export_btn.on("click", lambda e: _export_csv())

    async def _run_de():
        run_btn.disable()
        status.text = "Running DE analysis…"
        try:
            method = str(method_sel.value)
            n = max(5, min(int(float(n_genes_in.value or 25)), 200))

            def _work():
                return run_differential_expression(adata, SELECTED_PROJECT_PATH, method=method, n_genes=n)

            new_results = await asyncio.to_thread(_work)
            de_results.clear()
            de_results.update(new_results)

            _setup_results_ui()
            status.text = f"Done. {len(de_results)} clusters, top {n} genes each."
            ui.notify("Differential expression complete", color="positive")
        except Exception as e:
            status.text = f"Failed: {e}"
            ui.notify(f"DE failed: {e}", color="negative")
        finally:
            run_btn.enable()

    run_btn.on("click", lambda e: asyncio.create_task(_run_de()))


@ui.page("/receptor_ligand")
def receptor_ligand_page():
    ui.dark_mode().enable()

    ui.label("Receptor–Ligand").style(
        "font-size: 32px; font-weight: 800; color: #fb7185; margin: 12px 0;"
    )
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    # ----------------------------
    # Resolve project
    # ----------------------------
    SELECTED_PROJECT_PATH = None
    if SELECTED_PROJECT and isinstance(SELECTED_PROJECT, dict):
        SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path")
    if not SELECTED_PROJECT_PATH:
        ui.notify("No project selected. Returning to dashboard.", color="negative")
        ui.navigate.to("/")
        return
    SELECTED_PROJECT_PATH = str(SELECTED_PROJECT_PATH)

    adata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_annotated.h5ad"
    if not adata_path.exists():
        ui.label("adata_annotated.h5ad not found. Run PA/Annotations first.").style("color:#fbbf24;")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    # Load adata
    try:
        adata = read_adata(adata_path)
    except Exception as e:
        ui.notify(f"Failed to load adata: {e}", color="negative")
        ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat")
        return

    # Detect cluster key
    try:
        cluster_series = _get_cluster_series_from_adata(adata)
        # find the actual key name if possible
        cluster_key = None
        for k in ("leiden", "cluster", "clusters", "louvain"):
            if k in adata.obs.columns:
                cluster_key = k
                break
        if cluster_key is None:
            # fallback: first matching column
            for k in adata.obs.columns:
                lk = k.lower()
                if "leiden" in lk or "cluster" in lk or "louvain" in lk:
                    cluster_key = k
                    break
        clusters = sorted(cluster_series.unique().tolist(), key=_cluster_sort_key)
    except Exception as e:
        cluster_key = None
        clusters = []
        ui.notify(f"Could not detect clusters: {e}", color="negative")

    # ----------------------------
    # Find LR database files (optional)
    # ----------------------------
    # Expected formats:
    # - CSV/TSV with columns: ligand, receptor (case-insensitive)
    # Put one here: <project>/cellborg-cli/lr_db.csv  (or .tsv)
    lr_candidates = []
    for fn in ("lr_db.csv", "lr_db.tsv", "ligand_receptor.csv", "ligand_receptor.tsv"):
        p = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / fn
        if p.exists():
            lr_candidates.append(p)

    # ----------------------------
    # UI layout
    # ----------------------------
    with ui.row().classes("w-full").style("gap:18px; padding:12px; height: calc(100vh - 140px);"):
        # LEFT: controls
        with ui.column().style("width:520px; gap:12px;"):
            ui.label("Controls").style("color:#e5e7eb; font-weight:800; font-size:18px;")

            ctrl = ui.card().classes("w-full").style(
                "padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
            )
            with ctrl:
                ui.label("Groups").style("color:#9ca3af; font-weight:700;")
                if not clusters:
                    ui.label("No clusters found in adata.obs.").style("color:#fbbf24;")
                sender_sel = ui.select(
                    options=clusters,
                    value=(clusters[0] if clusters else None),
                    label=f"Sender cluster ({cluster_key or 'missing'})",
                ).props("outlined dense").style("width:100%;")

                receiver_sel = ui.select(
                    options=clusters,
                    value=(clusters[1] if len(clusters) > 1 else (clusters[0] if clusters else None)),
                    label=f"Receiver cluster ({cluster_key or 'missing'})",
                ).props("outlined dense").style("width:100%;")

                ui.separator().style("opacity:0.25; margin: 8px 0;")

                ui.label("LR database").style("color:#9ca3af; font-weight:700;")
                if lr_candidates:
                    lr_map = {str(p): p.name for p in lr_candidates}
                    lr_sel = ui.select(
                        options=lr_map,
                        value=str(lr_candidates[0]),
                        label="Ligand–Receptor file",
                    ).props("outlined dense").style("width:100%;")
                else:
                    lr_sel = None
                    ui.label(
                        "No LR DB file found.\nAdd one at: cellborg-cli/lr_db.csv (columns: ligand,receptor)"
                    ).style("white-space:pre-line; color:#fbbf24; font-size:14px;")

                ui.separator().style("opacity:0.25; margin: 8px 0;")

                top_n = ui.number("Top N interactions", value=50, format="%.0f").props("outlined dense").style("width:100%;")
                min_mean = ui.number("Min mean expr (ligand & receptor)", value=0.0, format="%.3f").props("outlined dense").style("width:100%;")
                use_log1p = ui.switch("log1p transform", value=False).style("color:#e5e7eb;")

                run_btn = ui.button("Compute interactions").props("unelevated").style(
                    "width:100%; height:44px; border-radius:10px; background:#fb7185; color:#0b1220; font-weight:900;"
                )

                ui.button("Back to Analysis", on_click=lambda: ui.navigate.to("/analysis")).props("flat").style("width:100%;")

        # RIGHT: outputs
        with ui.column().style("flex:1; gap:12px;"):
            ui.label("Predicted interactions").style(
                "font-size: 20px; font-weight: 800; color:#e5e7eb; text-align:center; margin-top:6px;"
            )

            out_card = ui.card().classes("w-full").style(
                "flex:1; padding:12px; border:1px solid #374151; border-radius:10px; background:#0b1220;"
            )
            with out_card:
                status = ui.label("Pick sender/receiver and click Compute interactions.").style("color:#9ca3af;")
                plot_el = ui.plotly(go.Figure()).style("width:100%; height: 520px;")
                table_el = ui.table(
                    columns=[
                        {"name": "rank", "label": "#", "field": "rank", "align": "left"},
                        {"name": "ligand", "label": "Ligand", "field": "ligand"},
                        {"name": "receptor", "label": "Receptor", "field": "receptor"},
                        {"name": "lig_mean", "label": "Lig mean (sender)", "field": "lig_mean"},
                        {"name": "rec_mean", "label": "Rec mean (receiver)", "field": "rec_mean"},
                        {"name": "score", "label": "Score", "field": "score"},
                    ],
                    rows=[],
                    row_key="rank",
                ).classes("w-full").style("margin-top:10px;")

    # ----------------------------
    # Helpers
    # ----------------------------
    def _read_lr_db(path: Path) -> pd.DataFrame:
        # Accept csv/tsv. Must include ligand/receptor cols.
        if path.suffix.lower() == ".tsv":
            df = pd.read_csv(path, sep="\t")
        else:
            df = pd.read_csv(path)

        cols = {c.lower(): c for c in df.columns}
        if "ligand" not in cols or "receptor" not in cols:
            raise ValueError("LR DB must have columns named ligand and receptor (case-insensitive).")

        out = df[[cols["ligand"], cols["receptor"]]].copy()
        out.columns = ["ligand", "receptor"]
        out["ligand"] = out["ligand"].astype(str)
        out["receptor"] = out["receptor"].astype(str)
        out = out.dropna().drop_duplicates()
        return out

    def _mean_expr_for_genes(adata_sub, genes: list[str]) -> dict[str, float]:
        """
        Return mean expression per gene for adata_sub.
        Works with sparse/dense X.
        """
        genes = [g for g in genes if g in adata_sub.var_names]
        if not genes:
            return {}

        # slice in var dimension once
        sub = adata_sub[:, genes]
        X = sub.X

        # mean per column
        try:
            import scipy.sparse as sp  # type: ignore
            if sp.issparse(X):
                m = np.asarray(X.mean(axis=0)).ravel()
            else:
                m = np.asarray(X).mean(axis=0)
        except Exception:
            m = np.asarray(X).mean(axis=0)

        if bool(use_log1p.value):
            m = np.log1p(np.clip(m, 0.0, None))

        return {g: float(v) for g, v in zip(genes, m)}

    def _build_scatter(df: pd.DataFrame, title: str) -> go.Figure:
        # Scatter: lig_mean vs rec_mean colored by score
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df["lig_mean"],
            y=df["rec_mean"],
            mode="markers",
            marker=dict(
                size=10,
                color=df["score"],
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title="Score"),
            ),
            text=df["ligand"] + " → " + df["receptor"],
            hovertemplate="%{text}<br>lig=%{x:.4f}<br>rec=%{y:.4f}<br>score=%{marker.color:.4f}<extra></extra>",
            showlegend=False,
        ))
        fig.update_layout(
            title=title,
            height=520,
            margin=dict(l=40, r=20, t=60, b=40),
            xaxis_title="Ligand mean (sender)",
            yaxis_title="Receptor mean (receiver)",
        )
        return fig

    # ----------------------------
    # Compute handler (slot-safe)
    # ----------------------------
    async def run_rl(e=None):
        with out_card:
            run_btn.disable()
            try:
                if not cluster_key or not clusters:
                    ui.notify("Clusters not available. Run Leiden/PA first.", color="negative")
                    return

                sender = str(sender_sel.value)
                receiver = str(receiver_sel.value)
                if not sender or not receiver:
                    ui.notify("Select sender and receiver.", color="negative")
                    return

                if lr_sel is None:
                    ui.notify("Add an LR DB file first (cellborg-cli/lr_db.csv).", color="negative")
                    return

                lr_path = Path(str(lr_sel.value))
                if not lr_path.exists():
                    ui.notify(f"LR DB not found: {lr_path}", color="negative")
                    return

                N = int(float(top_n.value or 50))
                N = max(5, min(N, 500))
                thr = float(min_mean.value or 0.0)

                status.text = "Loading LR DB + computing means…"

                def _work():
                    lr = _read_lr_db(lr_path)

                    # restrict to genes that exist in var_names
                    ligands = sorted(set(lr["ligand"].tolist()))
                    receptors = sorted(set(lr["receptor"].tolist()))

                    # subset sender/receiver cells
                    s_mask = (adata.obs[cluster_key].astype(str) == sender).to_numpy()
                    r_mask = (adata.obs[cluster_key].astype(str) == receiver).to_numpy()

                    if s_mask.sum() == 0 or r_mask.sum() == 0:
                        raise ValueError("Sender or receiver cluster has 0 cells.")

                    ad_s = adata[s_mask, :]
                    ad_r = adata[r_mask, :]

                    lig_means = _mean_expr_for_genes(ad_s, ligands)
                    rec_means = _mean_expr_for_genes(ad_r, receptors)

                    rows = []
                    for _, row in lr.iterrows():
                        lig = row["ligand"]
                        rec = row["receptor"]
                        lm = lig_means.get(lig)
                        rm = rec_means.get(rec)
                        if lm is None or rm is None:
                            continue
                        if lm < thr or rm < thr:
                            continue
                        score = lm * rm
                        rows.append((lig, rec, lm, rm, score))

                    out = pd.DataFrame(rows, columns=["ligand", "receptor", "lig_mean", "rec_mean", "score"])
                    if out.empty:
                        return out

                    out = out.sort_values("score", ascending=False).head(N).reset_index(drop=True)
                    out.insert(0, "rank", np.arange(1, len(out) + 1))
                    return out

                df = await asyncio.to_thread(_work)

                if df is None or df.empty:
                    status.text = "No interactions passed filters (genes missing or below threshold)."
                    table_el.rows = []
                    table_el.update()
                    plot_el.figure = go.Figure()
                    plot_el.update()
                    ui.notify("No interactions found.", color="warning")
                    return

                # render
                status.text = f"Top {len(df)} interactions: sender {sender} → receiver {receiver}"
                table_el.rows = [
                    {
                        "rank": int(r.rank),
                        "ligand": r.ligand,
                        "receptor": r.receptor,
                        "lig_mean": f"{float(r.lig_mean):.4f}",
                        "rec_mean": f"{float(r.rec_mean):.4f}",
                        "score": f"{float(r.score):.4f}",
                    }
                    for r in df.itertuples(index=False)
                ]
                table_el.update()

                fig = _build_scatter(df, title=f"LR interactions (sender {sender} → receiver {receiver})")
                plot_el.figure = fig
                plot_el.update()

                ui.notify("Computed receptor–ligand interactions", color="positive")

            except Exception as ex:
                status.text = f"Failed: {ex}"
                ui.notify(f"Receptor–ligand failed: {ex}", color="negative")
            finally:
                run_btn.enable()

    run_btn.on("click", run_rl)


# -----------------------------
# run
# -----------------------------
if __name__ in {"__main__", "__mp_main__"}:
    ui.run(host="127.0.0.1", port=8080, native=True)
