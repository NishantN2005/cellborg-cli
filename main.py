import sys
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QPushButton, QLabel
)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("scRNA-seq Tool")
        self.resize(900, 600)

        central = QWidget()
        layout = QVBoxLayout(central)

        self.status = QLabel("Ready")
        run_btn = QPushButton("Run QC")

        run_btn.clicked.connect(self.run_qc)

        layout.addWidget(self.status)
        layout.addWidget(run_btn)
        self.setCentralWidget(central)

    def run_qc(self):
        self.status.setText("Running QC...")

app = QApplication(sys.argv)
window = MainWindow()
window.show()
sys.exit(app.exec())
