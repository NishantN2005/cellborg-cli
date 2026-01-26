#!/usr/bin/env python3
import sys
import os
import shutil
from pathlib import Path

from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QLabel, QFileDialog,
    QMessageBox, QScrollArea, QFrame, QHBoxLayout, QDialog, QLineEdit, QTextEdit,
    QFormLayout, QDialogButtonBox, QSpinBox, QDoubleSpinBox
)
from PySide6.QtCore import Qt
import json
import webbrowser
import tempfile
from PySide6.QtCore import QUrl

# Prefer to embed Plotly in a Qt WebEngine view if available, otherwise open in browser
try:
    from PySide6.QtWebEngineWidgets import QWebEngineView
    WEBENGINE_AVAILABLE = True
except Exception:
    QWebEngineView = None
    WEBENGINE_AVAILABLE = False

import plotly.graph_objs as go
import plotly.io as pio
from plotly.subplots import make_subplots

from qc_functions import read_10x_mtx, calculate_qc_metrics, find_species, voilin_plot, SPECIES_TO_MT

BASE_DIR = Path(__file__).resolve().parent
PROJECTS_DIR = BASE_DIR / "projects"
PROJECTS_DIR.mkdir(exist_ok=True)

# Module-level selected project attributes (updated when a project is selected)
SELECTED_PROJECT_NAME = None
SELECTED_PROJECT_PATH = None
SELECTED_PROJECT_METADATA = None


class ProjectEntry(QFrame):
    def __init__(self, folder_name: str, title: str | None = None, on_click=None, parent=None):
        super().__init__(parent)
        self.folder = folder_name
        self.on_click = on_click
        display = title or folder_name
        self.setFrameShape(QFrame.StyledPanel)
        self.setObjectName("projectEntry")
        layout = QHBoxLayout(self)
        layout.setContentsMargins(10, 6, 10, 6)
        label = QLabel(display)
        label.setStyleSheet("font-weight:600")
        layout.addWidget(label)
        layout.addStretch()
        # make it appear clickable
        self.setCursor(Qt.PointingHandCursor)

    def mousePressEvent(self, event):
        if callable(self.on_click):
            try:
                self.on_click(self.folder)
            except Exception:
                pass
        super().mousePressEvent(event)


class Dashboard(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Cellborg Projects Dashboard")
        self.resize(600, 700)

        root = QVBoxLayout(self)
        # Top controls
        controls = QHBoxLayout()
        self.new_btn = QPushButton("+ New Project")
        self.new_btn.clicked.connect(self.add_project)
        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.clicked.connect(self.refresh)
        controls.addWidget(self.new_btn)
        controls.addWidget(self.refresh_btn)
        controls.addStretch()
        root.addLayout(controls)

        # Main area: left=project list, right=detail view
        main_h = QHBoxLayout()

        # Left: Scroll area for stacked project widgets
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.list_container = QWidget()
        self.list_layout = QVBoxLayout(self.list_container)
        self.list_layout.setAlignment(Qt.AlignTop)
        self.scroll.setWidget(self.list_container)
        main_h.addWidget(self.scroll, 2)

        # Right: project detail view (initially empty)
        self.detail = QFrame()
        self.detail.setFrameShape(QFrame.StyledPanel)
        detail_layout = QVBoxLayout(self.detail)
        detail_layout.setContentsMargins(12, 12, 12, 12)
        self.title_label = QLabel("")
        self.title_label.setStyleSheet("font-size:16px; font-weight:700")
        self.desc_label = QLabel("")
        self.desc_label.setWordWrap(True)
        detail_layout.addWidget(self.title_label)
        detail_layout.addWidget(self.desc_label)
        detail_layout.addStretch()

        # Action buttons
        self.run_qc_btn = QPushButton("Run QC")
        self.run_analysis_btn = QPushButton("Run Analysis")
        # not wired yet (do nothing)
        self.run_qc_btn.clicked.connect(self.run_qc)
        self.run_analysis_btn.clicked.connect(lambda: None)
        btns = QHBoxLayout()
        btns.addWidget(self.run_qc_btn)
        btns.addWidget(self.run_analysis_btn)
        detail_layout.addLayout(btns)

        main_h.addWidget(self.detail, 3)

        root.addLayout(main_h)

        self.info = QLabel("")
        root.addWidget(self.info)

        self.refresh()

    def run_qc(self):
        # Guard: ensure a project is selected
        if not SELECTED_PROJECT_PATH:
            QMessageBox.warning(self, "No project selected", "Please select a project on the left first.")
            return

        try:
            adata = read_10x_mtx(SELECTED_PROJECT_PATH)
            print("Read adata with shape:", adata.shape)
        except Exception as e:
            QMessageBox.critical(self, "Read failed", f"Failed to read project data: {e}")
            return

        # determine mt prefix from metadata species when available
        mt = "MT-"
        try:
            species = SELECTED_PROJECT_METADATA.get("species") if SELECTED_PROJECT_METADATA else None
            if species:
                mt = SPECIES_TO_MT.get(species, mt)
        except Exception:
            mt = "MT-"

        try:
            calculate_qc_metrics(mt, adata)
        except Exception as e:
            QMessageBox.critical(self, "QC failed", f"QC calculation failed: {e}")
        try:
            voilin_plot(adata, SELECTED_PROJECT_PATH)
        except Exception as e:
            QMessageBox.critical(self, "Plot failed", f"Failed to create violin plot: {e}")
            return

        # After voilin_plot runs, attempt to read generated highcharts_data.json
        try:
            hc_path = Path(SELECTED_PROJECT_PATH) / "cellborg-cli" / "highcharts_data.json"
            if not hc_path.exists():
                raise FileNotFoundError(f"{hc_path} not found")
            with hc_path.open("r", encoding="utf-8") as fh:
                hc = json.load(fh)

            # extract arrays for the three violin metrics
            # Support two JSON layouts:
            # 1) {"n_genes_by_counts": [...], "total_counts": [...], ...}
            # 2) {"CELLID1": {"n_genes": ..., "total_counts": ..., "pct_counts_mt": ...}, ...}
            a1 = []
            a2 = []
            a3 = []
            print("hc: ", hc)
            if isinstance(hc, dict):
                # assume  values are dicts per cell id
                vals = list(hc.values())
                if vals and isinstance(vals[0], dict):
                    for v in vals:
                        try:
                            if v is None:
                                continue
                            a1.append(float(v.get("n_genes"))) if v.get("n_genes") is not None else None
                        except Exception:
                            pass
                        try:
                            a2.append(float(v.get("total_counts"))) if v.get("total_counts") is not None else None
                        except Exception:
                            pass
                        try:
                            a3.append(float(v.get("pct_counts_mt"))) if v.get("pct_counts_mt") is not None else None
                        except Exception:
                            pass

            # build row tuples for plotting and filtering
                rows = list(zip(a1, a2, a3))
                dlg = ViolinDialog(rows, project_path=SELECTED_PROJECT_PATH, parent=self)
                dlg.exec()
        except Exception as e:
            QMessageBox.warning(self, "Plot data", f"Could not load violin plot data: {e}")
            return
        
    def refresh(self):
        # Clear existing entries
        for i in reversed(range(self.list_layout.count())):
            w = self.list_layout.itemAt(i).widget()
            if w:
                w.setParent(None)

        dirs = [p for p in PROJECTS_DIR.iterdir() if p.is_dir()]
        dirs.sort(key=lambda p: p.name)
        if not dirs:
            empty = QLabel('No projects found. Use "New Project" to add one.')
            empty.setStyleSheet('color: #6b7280; padding: 12px')
            self.list_layout.addWidget(empty)
        else:
            for p in dirs:
                # try to read metadata.json inside the project folder
                title = None
                meta_file = p / "metadata.json"
                if meta_file.exists():
                    try:
                        with meta_file.open("r", encoding="utf-8") as fh:
                            meta = json.load(fh)
                            title = meta.get("title")
                    except Exception:
                        title = None

                entry = ProjectEntry(p.name, title, on_click=self.show_project)
                self.list_layout.addWidget(entry)

        self.info.setText(f"Projects: {len(dirs)}")

    def add_project(self):
        # Let user pick a folder to move into projects/
        folder = QFileDialog.getExistingDirectory(self, "Select folder to add as project")
        if not folder:
            return
        src = Path(folder)
        dest = PROJECTS_DIR / src.name

        #check if folder is a valid project (contains expected files)
        required_files = ['barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz']
        if not all((src / f).exists() for f in required_files):
            QMessageBox.critical(self, "Invalid Project", f"The selected folder does not contain the required files: {', '.join(required_files)}")
            return
        
        #if dest exists, ask user if they want to overwrite, rename, or cancel
        if dest.exists():
            resp = QMessageBox.question(
                self, "Conflict", f"A project named '{src.name}' already exists. Overwrite?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel,
                QMessageBox.No
            )
            if resp == QMessageBox.Cancel:
                return
            if resp == QMessageBox.No:
                # append suffix
                i = 1
                while True:
                    candidate = PROJECTS_DIR / f"{src.name}-{i}"
                    if not candidate.exists():
                        dest = candidate
                        break
                    i += 1
            else:
                # Yes: remove existing dest
                try:
                    if dest.is_dir():
                        shutil.rmtree(dest)
                    else:
                        dest.unlink()
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to remove existing project: {e}")
                    return

        try:
            # prompt user for metadata (title, description)
            dlg = MetadataDialog(default_title=src.name)
            if dlg.exec() != QDialog.Accepted:
                return
            meta = dlg.get_data()

            shutil.move(str(src), str(dest))
            # create cellborg-cli subdir inside project for metadata/cache
            try:
                cfg_dir = Path(dest) / "cellborg-cli"
                cfg_dir.mkdir(parents=True, exist_ok=True)

                # detect species using qc function (read features file inside moved folder)
                try:
                    features_path = Path(dest) / "features.tsv.gz"
                    species = find_species(str(features_path))
                except Exception as e:
                    print(f"Could not determine species: {e}")
                    species = None
                if species:
                    meta["species"] = species

                # write metadata.json inside moved project's cellborg-cli folder
                with (cfg_dir / "metadata.json").open("w", encoding="utf-8") as fh:
                    json.dump(meta, fh, ensure_ascii=False, indent=2)
            except Exception as e:
                QMessageBox.warning(self, "Warning", f"Saved project but failed to write metadata: {e}")
        except Exception as e:
            QMessageBox.critical(self, "Move failed", f"Could not move folder: {e}")
            return
        
        QMessageBox.information(self, "Added", f"Project '{dest.name}' added to projects/")
        self.refresh()

    def show_project(self, folder_name: str):
        p = PROJECTS_DIR / folder_name
        if not p.exists() or not p.is_dir():
            self.title_label.setText("")
            self.desc_label.setText("")
            return

        # load metadata (from project/cellborg-cli/) and update global selected attributes
        global SELECTED_PROJECT_NAME, SELECTED_PROJECT_PATH, SELECTED_PROJECT_METADATA
        meta_file = p / "cellborg-cli" / "metadata.json"
        title = folder_name
        desc = ""
        meta = {}
        if meta_file.exists():
            try:
                with meta_file.open("r", encoding="utf-8") as fh:
                    meta = json.load(fh) or {}
                    title = meta.get("title") or title
                    desc = meta.get("description") or ""
            except Exception:
                meta = {}

        SELECTED_PROJECT_NAME = title
        SELECTED_PROJECT_PATH = str(p)
        SELECTED_PROJECT_METADATA = meta

        # also update the UI
        self.title_label.setText(title)
        self.desc_label.setText(desc)


class MetadataDialog(QDialog):
    def __init__(self, default_title: str = "", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Project metadata")
        self.setModal(True)
        layout = QFormLayout(self)
        self.title_edit = QLineEdit(default_title)
        self.desc_edit = QTextEdit()
        self.desc_edit.setFixedHeight(120)
        layout.addRow("Title:", self.title_edit)
        layout.addRow("Description:", self.desc_edit)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def get_data(self):
        return {"title": self.title_edit.text().strip() or None, "description": self.desc_edit.toPlainText().strip()}


class ViolinDialog(QDialog):
    """Modal dialog that shows three violin plots using Plotly and provides filters.

    The dialog accepts `rows` as an iterable of tuples: (n_genes, total_counts, pct_counts_mt).
    When Qt WebEngine is available the plot is embedded and updates in-place when the
    filter controls change. Otherwise the plot is opened in the browser and controls
    are disabled.
    """
    def __init__(self, rows, project_path: str | None = None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Violin + Scatter")
        self.resize(1100, 750)
        self.rows = [tuple(r) for r in (rows or [])]
        layout = QVBoxLayout(self)

        # Controls: min nFeature, max nFeature, pct mt less than
        ctrl_layout = QHBoxLayout()
        self.min_feat = QSpinBox()
        self.min_feat.setRange(0, 10_000_000)
        self.min_feat.setValue(0)
        self.max_feat = QSpinBox()
        self.max_feat.setRange(0, 10_000_000)
        # default max = max n_genes in data or 10000
        max_default = max((int(r[0]) for r in self.rows if r and r[0] is not None), default=10000)
        self.max_feat.setValue(max_default)
        self.pct_mt = QDoubleSpinBox()
        self.pct_mt.setRange(0.0, 100.0)
        self.pct_mt.setDecimals(3)
        self.pct_mt.setValue(100.0)

        ctrl_layout.addWidget(QLabel("nFeature_RNA greater than:"))
        ctrl_layout.addWidget(self.min_feat)
        ctrl_layout.addWidget(QLabel("nFeature_RNA less than:"))
        ctrl_layout.addWidget(self.max_feat)
        ctrl_layout.addWidget(QLabel("% mitochondrial genes less than:"))
        ctrl_layout.addWidget(self.pct_mt)
        ctrl_layout.addStretch()
        layout.addLayout(ctrl_layout)

        # connect to update
        self.min_feat.valueChanged.connect(self._on_filter_change)
        self.max_feat.valueChanged.connect(self._on_filter_change)
        self.pct_mt.valueChanged.connect(self._on_filter_change)

        # content area for plot or info
        self._view = None
        self._info_label = None
        self._tmpfile = None

        if WEBENGINE_AVAILABLE:
            self._view = QWebEngineView()
            layout.addWidget(self._view, 1)
        else:
            self._info_label = QLabel("Plots will open in your default browser (Qt WebEngine not available). Filters are disabled.")
            layout.addWidget(self._info_label)
            # disable controls (cannot re-render in browser fallback)
            self.min_feat.setEnabled(False)
            self.max_feat.setEnabled(False)
            self.pct_mt.setEnabled(False)

        # Save button + Close
        btn_row = QHBoxLayout()
        self.save_btn = QPushButton("Save Filters")
        self.save_btn.clicked.connect(self._on_save)
        btn_row.addWidget(self.save_btn)
        btn_row.addStretch()
        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(self.reject)
        btns.accepted.connect(self.accept)
        btn_row.addWidget(btns)
        layout.addLayout(btn_row)

        # store project path for saving metadata
        self._project_path = project_path

        # initial render
        self._render_plot(self.rows)

    def _on_filter_change(self, *_):
        # filter rows according to controls and re-render
        minv = self.min_feat.value()
        maxv = self.max_feat.value()
        pct_max = float(self.pct_mt.value())
        filtered = []
        for r in self.rows:
            try:
                n_genes = float(r[0])
                total = float(r[1])
                pct = float(r[2])
            except Exception:
                continue
            if n_genes >= minv and n_genes <= maxv and pct <= pct_max:
                filtered.append((n_genes, total, pct))

        self._render_plot(filtered)

    def _render_plot(self, rows):
        # build data arrays from rows
        a_n = [r[0] for r in rows]
        a_tot = [r[1] for r in rows]
        a_pct = [r[2] for r in rows]

        titles = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
        specs = [[{}, {}, {}], [{"colspan": 3}, None, None]]
        fig = make_subplots(rows=2, cols=3, specs=specs, subplot_titles=(titles + ["Scatter: total_counts vs n_genes_by_counts (color=pct_counts_mt)"]))

        # violins
        for i, arr in enumerate((a_n, a_tot, a_pct), start=1):
            y = []
            try:
                y = [float(x) for x in arr]
            except Exception:
                y = []
            if y:
                fig.add_trace(go.Violin(y=y, name=titles[i-1], box_visible=True, meanline_visible=True), row=1, col=i)
            else:
                fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers', marker_opacity=0), row=1, col=i)

        # scatter
        x_vals = a_tot
        y_vals = a_n
        c_vals = a_pct
        if x_vals and y_vals:
            scatter = go.Scatter(
                x=x_vals,
                y=y_vals,
                mode='markers',
                marker=dict(color=c_vals if c_vals else None, colorscale='Viridis', showscale=True, colorbar=dict(title='pct_counts_mt')),
                text=[f"pct_counts_mt: {v}" for v in c_vals] if c_vals else None,
            )
            fig.add_trace(scatter, row=2, col=1)
        else:
            fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers', marker_opacity=0), row=2, col=1)

        fig.update_layout(height=780, width=980, showlegend=False)

        html = pio.to_html(fig, full_html=True, include_plotlyjs=True)

        # remove previous tempfile
        try:
            if self._tmpfile:
                os.unlink(self._tmpfile)
        except Exception:
            pass

        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.html')
        tmp.write(html.encode('utf-8'))
        tmp.flush()
        tmp.close()
        self._tmpfile = tmp.name

        if WEBENGINE_AVAILABLE and self._view is not None:
            # load file in web view
            self._view.load(QUrl.fromLocalFile(tmp.name))
        else:
            # open in browser (fallback)
            webbrowser.open(f'file://{tmp.name}')

    def closeEvent(self, event):
        try:
            if hasattr(self, '_tmpfile') and self._tmpfile:
                os.unlink(self._tmpfile)
        except Exception:
            pass
        super().closeEvent(event)

    def _on_save(self):
        """Write current filter values to the project's metadata.json under cellborg-cli/."""
        if not self._project_path:
            QMessageBox.warning(self, "Save failed", "No project path provided; cannot save filters.")
            return

        proj = Path(self._project_path)
        cfg_dir = proj / "cellborg-cli"
        try:
            cfg_dir.mkdir(parents=True, exist_ok=True)
            meta_file = cfg_dir / "metadata.json"
            meta = {}
            if meta_file.exists():
                try:
                    with meta_file.open("r", encoding="utf-8") as fh:
                        meta = json.load(fh) or {}
                except Exception:
                    meta = {}

            # write filters under a `filters` key
            filters = {
                "nFeature_RNA_min": int(self.min_feat.value()),
                "nFeature_RNA_max": int(self.max_feat.value()),
                "pct_counts_mt_max": float(self.pct_mt.value()),
            }
            meta["filters"] = filters

            with meta_file.open("w", encoding="utf-8") as fh:
                json.dump(meta, fh, ensure_ascii=False, indent=2)

            QMessageBox.information(self, "Saved", f"Filters saved to {meta_file}")
        except Exception as e:
            QMessageBox.critical(self, "Save failed", f"Failed to save filters: {e}")


def main():
    app = QApplication(sys.argv)
    win = Dashboard()
    win.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()


def get_selected_project_name():
    """Return the currently selected project's title (or None)."""
    return SELECTED_PROJECT_NAME


def get_selected_project_path():
    """Return the filesystem path to the currently selected project (or None)."""
    return SELECTED_PROJECT_PATH


def get_selected_project_metadata():
    """Return the metadata dict for the currently selected project (or None)."""
    return SELECTED_PROJECT_METADATA
