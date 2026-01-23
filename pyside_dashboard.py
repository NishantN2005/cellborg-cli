#!/usr/bin/env python3
import sys
import os
import shutil
from pathlib import Path

from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QLabel, QFileDialog,
    QMessageBox, QScrollArea, QFrame, QHBoxLayout, QDialog, QLineEdit, QTextEdit,
    QFormLayout, QDialogButtonBox
)
from PySide6.QtCore import Qt
import json

BASE_DIR = Path(__file__).resolve().parent
PROJECTS_DIR = BASE_DIR / "projects"
PROJECTS_DIR.mkdir(exist_ok=True)


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
        self.run_qc_btn.clicked.connect(lambda: None)
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
            # write metadata.json inside moved folder
            try:
                with (Path(dest) / "metadata.json").open("w", encoding="utf-8") as fh:
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

        meta_file = p / "metadata.json"
        title = folder_name
        desc = ""
        if meta_file.exists():
            try:
                with meta_file.open("r", encoding="utf-8") as fh:
                    meta = json.load(fh)
                    title = meta.get("title") or title
                    desc = meta.get("description") or ""
            except Exception:
                pass

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


def main():
    app = QApplication(sys.argv)
    win = Dashboard()
    win.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
