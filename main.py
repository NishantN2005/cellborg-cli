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
    gene_expression
)

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


def umap_scatter_figure_from_dict(umap_dict: dict, title="UMAP (Leiden clusters)"):
    """
    umap_dict format:
    {
        "CELLID": {"UMAP1": float, "UMAP2": float, "cluster": "0"},
        ...
    }
    One trace per cluster => Plotly assigns different colors automatically.
    """
    df = pd.DataFrame.from_dict(umap_dict, orient="index")

    df["UMAP1"] = pd.to_numeric(df.get("UMAP1"), errors="coerce")
    df["UMAP2"] = pd.to_numeric(df.get("UMAP2"), errors="coerce")
    df["cluster"] = df.get("cluster").astype(str)

    df = df.dropna(subset=["UMAP1", "UMAP2"])

    def _cluster_key(x: str):
        try:
            return (0, int(x))
        except Exception:
            return (1, x)

    clusters = sorted(df["cluster"].unique(), key=_cluster_key)

    fig = go.Figure()
    for c in clusters:
        sub = df[df["cluster"] == c]
        fig.add_trace(
            go.Scattergl(
                x=sub["UMAP1"],
                y=sub["UMAP2"],
                mode="markers",
                name=f"Cluster {c}",
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
        ui.button("Run Analysis").style("flex:1;")


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
            # NOTE: fixed key "project_path" (your old code used "projzect_path")
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

    dark = ui.dark_mode()
    dark.enable()

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

    dark = ui.dark_mode()
    dark.enable()

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

            # update metadata status
            md_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "metadata.json"
            if md_path.exists():
                try:
                    md = json.loads(md_path.read_text(encoding="utf-8"))
                    md["status"] = "PROC_ANNO"
                    md_path.write_text(json.dumps(md, indent=4), encoding="utf-8")
                    # update selected project in memory too
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
    PCA_DONE = False  # global flag to track if PCA has been done (enables "Next" button)
    SELECTED_PROJECT_PATH = SELECTED_PROJECT.get("project_path") if SELECTED_PROJECT else None

    if not SELECTED_PROJECT_PATH:
        ui.label("No project selected. Go back and select one.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return

    dark = ui.dark_mode()
    dark.enable()

    ui.label("Processing and Annotation").style("font-size: 32px; font-weight: 800; color: #4ecda4; margin: 12px 0;")
    #ui.label(f"Project: {SELECTED_PROJECT.get('project_title', '')}").style("color:#ddd;")
    #ui.label(SELECTED_PROJECT.get("project_description", "")).style("color:#bbb;")
    ui.separator().style("opacity:0.25; margin: 12px 0;")

    adata_qc_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_qc.h5ad"
    if not adata_qc_path.exists():
        ui.label("QC dataset not found. Please run QC first.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return

    # Load + init once per page load
    adata = read_adata(adata_qc_path)
    init_project(SELECTED_PROJECT_PATH, adata)

    # Controls: resolution slider + run button
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

    async def run_and_render():
        nonlocal PCA_DONE
        PCA_DONE = False
        run_btn.disable()
        res = float(res_slider.value)
        status_lbl.text = f"Clustering at resolution {res:.2f}…"
        try:
            await asyncio.to_thread(do_clustering, adata, res)

            # load UMAP json + plot (each cluster different color via per-trace split)
            umap_dict = load_umap_dict(SELECTED_PROJECT_PATH, res)
            fig = umap_scatter_figure_from_dict(umap_dict, title=f"UMAP (Leiden res={res:.2f})")

            umap_plot.figure = fig
            umap_plot.update()
            status_lbl.text = f"Rendered UMAP for resolution {res:.2f}"
            PCA_DONE = True
        except Exception as e:
            status_lbl.text = f"Failed: {e}"
            ui.notify(f"PA failed: {e}", color="negative")
        finally:             
            next_btn.enable()
            run_btn.enable() 

    run_btn.on("click", lambda e: asyncio.create_task(run_and_render()))

    next_btn = ui.button(
        'Next →',
        on_click=lambda: ui.navigate.to("/anno"),
    ).props('unelevated')
    next_btn.disable()  # default locked

    # enable if PCA already completed
    if PCA_DONE:
        print("PCA already done, enabling Next button")
        next_btn.enable()

    with ui.row().classes('w-full justify-end mt-6'):
        next_btn


@ui.page("/anno")
def annotation_page():
    ui.dark_mode().enable()

    # ---- locate project + gene list ----
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

    #------ premptively load adata -------
    adata_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "adata_clustered.h5ad"
    if not adata_path.exists():
        ui.label("Clustered dataset not found. Please run PA first.").style("color:#fbbf24;")
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat")
        return 
    adata = read_adata(adata_path)
    
    ui.label("Annotations").style(
        "font-size: 32px; font-weight: 800; color: #4ecda4; margin: 12px 0;"
    )
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

            # ---- buttons under chart ----
            with ui.row().classes("w-full").style("gap:12px;"):
                cluster_btn = ui.button("Cluster Plot").style(
                    "flex:1; height:48px; border-radius:8px; background:#41c99a; color:#0b1220; font-weight:800;"
                )

                plot_name_btn = ui.button(
                    plot_type.value,
                    on_click=lambda: (
                        print(f"{plot_type.value} clicked"),
                        ui.notify(f"{plot_type.value} clicked", color="info"),
                    ),
                ).style(
                    "flex:1; height:48px; border-radius:8px; background:#41c99a; color:#0b1220; font-weight:800;"
                )

            async def load_cluster_plot():
                # IMPORTANT: always enter a slot when updating UI from async handlers
                with plot_card:
                    cluster_btn.disable()
                    try:
                        print('------------ Loading cluster plot ---', pv.get("clust_resolution"))
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
                        fig = umap_scatter_figure_from_dict(umap_dict, title=f"UMAP (Leiden res={res:.2f})")

                        plot_el.figure = fig
                        plot_el.update()

                        plot_status.text = f"Loaded cluster UMAP (res={res:.2f})"
                        ui.notify("Loaded cluster plot", color="positive")

                    except Exception as e:
                        plot_status.text = f"Failed to load cluster plot: {e}"
                        ui.notify(f"Cluster plot failed: {e}", color="negative")
                    finally:
                        cluster_btn.enable()

            # DON’T use asyncio.create_task here; keep NiceGUI context intact
            cluster_btn.on("click", lambda e: load_cluster_plot())

            def _cluster_sort_key(x: str):
                try:
                    return (0, int(x))
                except Exception:
                    return (1, str(x))

            def load_gene_expression_df(project_path: str | Path) -> pd.DataFrame:
                p = Path(project_path) / "cellborg-cli" / "figures" / "gene_expression_per_cell_with_clusters.json"
                if not p.exists():
                    raise FileNotFoundError(f"Gene expression json not found: {p}")
                data = json.loads(p.read_text(encoding="utf-8"))
                df = pd.DataFrame.from_dict(data, orient="index")
                # normalize types
                df["UMAP1"] = pd.to_numeric(df.get("UMAP1"), errors="coerce")
                df["UMAP2"] = pd.to_numeric(df.get("UMAP2"), errors="coerce")
                df["cluster"] = df.get("cluster").astype(str)
                df = df.dropna(subset=["UMAP1", "UMAP2", "cluster"])
                return df

            def feature_plot_multi_gene(df: pd.DataFrame, genes: list[str]) -> go.Figure:
                genes = [g for g in genes if g in df.columns]
                if not genes:
                    raise KeyError("None of the selected genes exist in the gene-expression JSON.")

                # numeric + fill
                for g in genes:
                    df[g] = pd.to_numeric(df[g], errors="coerce").fillna(0.0)

                clusters = sorted(df["cluster"].unique(), key=_cluster_sort_key)

                fig = go.Figure()

                # trace bookkeeping: for each gene, we add:
                # - one Scattergl per cluster
                # - one hidden "colorbar" trace
                traces_per_gene = len(clusters) + 1

                for gi, gene in enumerate(genes):
                    # cluster traces
                    for c in clusters:
                        sub = df[df["cluster"] == c]
                        fig.add_trace(go.Scattergl(
                            x=sub["UMAP1"],
                            y=sub["UMAP2"],
                            mode="markers",
                            name=f"Cluster {c}",
                            marker=dict(
                                size=3,
                                color=sub[gene],
                                colorscale="Viridis",
                                showscale=False,   # only show on the gene colorbar trace
                            ),
                            hovertemplate=f"cluster={c}<br>{gene}=%{{marker.color:.3f}}<extra></extra>",
                            visible=(gi == 0),
                        ))

                    # shared colorbar trace for this gene
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

                # dropdown: toggle visibility gene-by-gene
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
                            {"title": f"Feature Plot: {gene} (per cluster)"}
                        ],
                    ))

                fig.update_layout(
                    title=f"Feature Plot: {genes[0]} (per cluster)",
                    height=650,
                    margin=dict(l=20, r=20, t=70, b=20),
                    legend=dict(itemsizing="constant"),
                    updatemenus=[dict(
                        type="dropdown",
                        x=0.01,
                        y=1.15,
                        xanchor="left",
                        yanchor="top",
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
                        args=[
                            {"visible": vis},
                            {"title": f"Violin: {gene} by cluster"}
                        ],
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
                        x=0.01,
                        y=1.15,
                        xanchor="left",
                        yanchor="top",
                        buttons=buttons,
                    )],
                )
                return fig

            def dot_plot_per_cluster(df: pd.DataFrame, genes: list[str]) -> go.Figure:
                genes = [g for g in genes if g in df.columns]
                if not genes:
                    raise KeyError("None of the selected genes exist in the gene-expression JSON.")

                clusters = sorted(df["cluster"].unique(), key=_cluster_sort_key)

                # stats: mean expr + pct expressing (>0)
                rows = []
                for c in clusters:
                    sub = df[df["cluster"] == c]
                    for g in genes:
                        v = pd.to_numeric(sub[g], errors="coerce").fillna(0.0)
                        mean_expr = float(v.mean())
                        pct_expr = float((v > 0).mean())  # 0..1
                        rows.append((c, g, mean_expr, pct_expr))

                stats = pd.DataFrame(rows, columns=["cluster", "gene", "mean_expr", "pct_expr"])

                # convert to plotly arrays
                x = stats["gene"].tolist()
                y = stats["cluster"].tolist()
                color = stats["mean_expr"].to_numpy()
                size = (stats["pct_expr"].to_numpy() * 30.0) + 3.0  # scale sizes

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

            def load_plot_shim():
                # 1) write JSON for current selected genes (your function)
                gene_expression(adata, selected_genes)  # make sure this writes to the figures json

                pt = plot_type.value
                genes_label = ", ".join(selected_genes[:5]) if selected_genes else "(no genes selected)"
                plot_status.text = f"Loaded: {pt} | Genes: {genes_label}"
                ui.notify(f"Loading {pt}", color="info")

                try:
                    # 2) load JSON -> df
                    df = load_gene_expression_df(SELECTED_PROJECT_PATH)

                    # 3) plot per cluster
                    if pt == "Feature Plot":
                        if not selected_genes:
                            raise ValueError("Select at least one gene for Feature Plot.")
                        gene = selected_genes[0]
                        fig = feature_plot_multi_gene(df, selected_genes)

                    elif pt == "Violin Plot":
                        if not selected_genes:
                            raise ValueError("Select at least one gene for Violin Plot.")
                        gene = selected_genes[0]
                        fig = violin_plot_multi_gene(df, selected_genes)

                    else:  # Dot Plot
                        if not selected_genes:
                            raise ValueError("Select at least one gene for Dot Plot.")
                        # keep it readable
                        genes = selected_genes[:12]
                        fig = dot_plot_per_cluster(df, genes)

                    plot_el.figure = fig
                    plot_el.update()

                except Exception as e:
                    plot_status.text = f"Failed to load plot: {e}"
                    ui.notify(f"Plot failed: {e}", color="negative")
                    return

                plot_name_btn.text = pt


            # keep the label button synced even before loading (when user changes radio)
            def on_plot_type_change():
                plot_name_btn.text = plot_type.value

            plot_type.on("change", lambda e: on_plot_type_change())
            on_plot_type_change()

            load_btn.on("click", lambda e: load_plot_shim())

# -----------------------------
# run
# -----------------------------
if __name__ in {"__main__", "__mp_main__"}:
    ui.run(host="127.0.0.1", port=8080, native=True)
