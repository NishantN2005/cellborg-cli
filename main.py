from nicegui import ui, app
import os
import shutil
import asyncio
from pathlib import Path
import json
from datetime import datetime
import numpy as np
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
    Converts your QC_THRESHOLDS (which currently stores ui.number elements)
    into plain floats so background threads + gate_adata don't see UI objects.
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

    if SELECTED_PROJECT:
        status = SELECTED_PROJECT.get("status", "unknown").upper()
        with ui.row().classes("w-full").style("gap:12px; margin-top:12px;"):
            if status == "QC":
                ui.button("Run QC", on_click=lambda: ui.navigate.to("/qc")).style("flex:1;")
            elif status == "PROC_ANNO":
                ui.button("Run PA", on_click=lambda: ui.notify("Not implemented yet")).style("flex:1;")
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
        title_input = ui.input(validation=lambda v: "Title already exists" if v and v.strip() in existing_titles() else None).props(
            "autofocus"
        )

        ui.label("Project description")
        desc_input = ui.textarea().style("height: 100px; width: 100%;")

        def save():
            md = {
                "original_folder_name": folder.name,
                "projzect_path": str(dest),
                "species": species,
                "project_title": (title_input.value or dest.name).strip(),
                "project_description": desc_input.value or "",
                "added_at": datetime.now().isoformat(),
                "status": "qc"
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

            # (left exactly as you had it)
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

    def get_line(fig, idx: int) -> float:
        return float(fig.layout.shapes[idx].y0)

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
                print(f"[{metric}] min={y_min} max={y_max}")

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

    # IMPORTANT: define 'actions' container (you were using it but never created it)
    actions = ui.row().classes("w-full").style("gap:12px;")
    with actions:
        ui.button("Back to Projects", on_click=lambda: ui.navigate.to("/")).props("flat").style("flex:1")

        save_btn = ui.button("Save").props("flat").style("flex:1")

        def on_save_click():
            save_btn.disable()
            with actions:
                ui.notify("Saving QC…", color="info")

            #change status in metadata.json file to qc complete
            md_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "metadata.json"
            if md_path.exists():
                try:
                    md = json.loads(md_path.read_text(encoding="utf-8"))
                    md["status"] = "proc_anno"
                    md_path.write_text(json.dumps(md, indent=4), encoding="utf-8")
                except Exception as e:
                    print(f"Failed to update metadata.json: {e}")

            asyncio.create_task(
                save_qc_async(
                    adata=adata,
                    qc_thresholds_elements=QC_THRESHOLDS,  # <- UI elements dict (we convert to floats inside)
                    project_path=SELECTED_PROJECT["project_path"],
                    ui_slot=actions,  # <- explicit UI slot
                    save_button=save_btn,
                )
            )

        save_btn.on("click", lambda e: on_save_click())


# -----------------------------
# run
# -----------------------------
if __name__ in {"__main__", "__mp_main__"}:
    ui.run(host="127.0.0.1", port=8080, native=True)
