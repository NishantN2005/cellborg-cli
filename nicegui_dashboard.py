from nicegui import ui, app
import os
import shutil
import asyncio
from pathlib import Path
import json
from datetime import datetime

PROJECTS_DIR = Path('projects')
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
        cand = PROJECTS_DIR / f'{name}-{i}'
        if not cand.exists():
            return cand
        i += 1


def load_project_metadata(project_dir: Path) -> dict | None:
    meta_path = project_dir / 'cellborg-cli' / 'metadata.json'
    if not meta_path.exists():
        return None
    try:
        data = json.loads(meta_path.read_text(encoding='utf-8'))
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
# native folder picker
# -----------------------------
async def choose_folder_via_native_file_dialog() -> Path | None:
    if not getattr(app, 'native', None) or app.native.main_window is None:
        ui.notify('Native window not available. Run with ui.run(native=True).', color='negative')
        return None

    try:
        import webview  # pip install pywebview
    except Exception as e:
        ui.notify(f'pywebview not installed/available: {e}', color='negative')
        return None

    dialog_type = webview.FileDialog.OPEN if hasattr(webview, 'FileDialog') else webview.OPEN_DIALOG

    files = await app.native.main_window.create_file_dialog(
        dialog_type=dialog_type,
        allow_multiple=False,
        file_types=('All files (*.*)',),
    )
    if not files:
        return None

    file_path = Path(files[0])
    if not file_path.exists():
        ui.notify(f'File not found: {file_path}', color='negative')
        return None

    folder = file_path.parent
    if not folder.is_dir():
        ui.notify(f'Not a folder: {folder}', color='negative')
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
            ui.button(f'{project.name} (no metadata)', on_click=lambda p=project: ui.notify(str(p)))
            continue

        title = str(md.get('project_title') or project.name).strip() or project.name
        ui.button(title, on_click=lambda m=md: select_project(m))


@ui.refreshable
def project_details():
    ui.label('Project Details').style('font-size: 28px; font-weight: 700; color: #4ecda4; margin: 6px 0;')

    if not SELECTED_PROJECT:
        ui.label('Select a project to see details here.').style('color: #555;')
        return

    ui.label(f"Title: {SELECTED_PROJECT.get('project_title', '')}").style('color:#ddd; font-size: 18px;')
    ui.label(f"Description: {SELECTED_PROJECT.get('project_description', '')}").style('color:#bbb;')

    with ui.row().classes('w-full').style('gap:12px; margin-top:12px;'):
        ui.button('Run QC', on_click=lambda: ui.navigate.to('/qc')).style('flex:1;')
        ui.button('Run Analysis').style('flex:1;')


async def add_project() -> None:
    ui.notify('Select any file inside the folder you want to add.', color='info')
    folder = await choose_folder_via_native_file_dialog()
    if folder is None:
        return

    try:
        dest = await copy_folder_to_projects(folder)
    except Exception as e:
        ui.notify(f'Failed to copy: {e}', color='negative')
        return

    dialog.clear()

    with dialog, ui.card().classes('w-96'):
        ui.label('Project Title')
        title_input = ui.input(
            validation=lambda v: 'Title already exists' if v and v.strip() in existing_titles() else None
        ).props('autofocus')

        ui.label('Project description')
        desc_input = ui.textarea().style('height: 100px; width: 100%;')

        def save():
            md = {
                "original_folder_name": folder.name,
                "project_path": str(dest),
                "project_title": (title_input.value or dest.name).strip(),
                "project_description": desc_input.value or "",
                "added_at": datetime.now().isoformat(),
            }
            (dest / "cellborg-cli").mkdir(parents=True, exist_ok=True)
            (dest / "cellborg-cli" / "metadata.json").write_text(json.dumps(md, indent=4), encoding='utf-8')

            ui.notify('Saved project metadata', color='positive')
            dialog.close()
            projects_list.refresh()

        with ui.row().classes('justify-end w-full'):
            ui.button('Cancel', on_click=dialog.close)
            ui.button('Save', on_click=save)

    dialog.open()


# -----------------------------
# Pages
# -----------------------------
@ui.page('/')
def dashboard_page():
    ensure_projects_dir()

    with ui.dialog() as d:
        # keep a global ref so add_project can use it
        global dialog
        dialog = d

    dark = ui.dark_mode()
    dark.enable()

    with ui.row().style(
        '''
        padding:12px;
        width:100%;
        height:95vh;
        box-sizing:border-box;
        gap:12px;
        overflow:hidden;
        '''
    ):
        with ui.column().style(
            '''
            width:520px;
            height:100%;
            padding:20px;
            border:1px solid #e5e7eb;
            border-radius:8px;
            box-shadow:0 1px 3px rgba(0,0,0,0.06);
            overflow-y:auto;
            '''
        ):
            ui.label('Projects').style('font-size: 36px; font-weight: 800; color: #4ecda4;')
            ui.button('Add New Project', on_click=add_project)
            projects_list()

        with ui.column().style(
            '''
            flex:1;
            height:100%;
            padding:20px;
            border:1px solid #e5e7eb;
            border-radius:8px;
            box-shadow:0 1px 3px rgba(0,0,0,0.06);
            overflow-y:auto;
            '''
        ):
            project_details()


@ui.page('/qc')
def qc_page():
    dark = ui.dark_mode()
    dark.enable()

    ui.label('QC Runner').style('font-size: 32px; font-weight: 800; color: #4ecda4; margin: 12px 0;')

    if SELECTED_PROJECT:
        ui.label(f"Project: {SELECTED_PROJECT.get('project_title', '')}").style('color:#ddd;')
        ui.label(SELECTED_PROJECT.get('project_description', '')).style('color:#bbb;')
    else:
        ui.label('No project selected. Go back and select one.').style('color:#fbbf24;')

    ui.button('Back to Projects', on_click=lambda: ui.navigate.to('/')).props('flat')


# -----------------------------
# run
# -----------------------------
if __name__ in {"__main__", "__mp_main__"}:
    ui.run(host='127.0.0.1', port=8080, native=True)
