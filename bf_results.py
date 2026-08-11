"""Results browser shown after a solver run.

Eight actions regenerate (via ParaView pvpython, only if missing) and then
display a result: the folds-pattern PNG in a QLabel, the rest as .avi videos in
a Qt-native QMediaPlayer/QVideoWidget. Save/Close absorb the old save.py and the
close-confirmation flow from make_run.py.
"""

from __future__ import annotations

from os.path import exists

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QPixmap
from PySide6.QtMultimedia import QMediaPlayer
from PySide6.QtMultimediaWidgets import QVideoWidget
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton,
    QStackedWidget, QVBoxLayout, QWidget,
)

import bf_runner

# label, output file under Folder_Output, pvpython script stem, is_image
ACTIONS = [
    ("Folds pattern", "Folder_Output/folds_pattern.png", "final_folding_pattren_{c}.py", True),
    ("Cell density distribution", "Folder_Output/Cell_desnity.avi", "cell_density_video_{c}.py", False),
    ("Stiffness distribution", "Folder_Output/Stiffness.avi", "Stiffness_video_{c}.py", False),
    ("Velocity distribution", "Folder_Output/Velocity.avi", "velocity_video_{c}.py", False),
    ("Tangential growth factor", "Folder_Output/growth_factor_t.avi", "growth_factors_video_{c}.py", False),
    ("Radial growth factor", "Folder_Output/growth_factor_r.avi", "growth_factors_video_{c}.py", False),
    ("RGCs proliferation", "Folder_Output/source_vz.avi", "source_terms_video_{c}.py", False),
    ("ORGCs proliferation", "Folder_Output/source_osvz.avi", "source_terms_video_{c}.py", False),
]


class ResultsWindow(QWidget):
    def __init__(self, case: str):
        super().__init__()
        self.case = case
        self._worker: bf_runner.ParaviewWorker | None = None
        self.setWindowTitle("Results Output")
        self.resize(760, 520)

        root = QHBoxLayout(self)
        root.setContentsMargins(16, 16, 16, 16)
        root.setSpacing(14)

        # --- left: action buttons ---
        buttons = QVBoxLayout()
        buttons.setSpacing(8)
        self._buttons: list[QPushButton] = []
        for label, output, script, is_image in ACTIONS:
            btn = QPushButton(label, objectName="actbtn")
            btn.setMinimumWidth(230)
            btn.clicked.connect(
                lambda _=False, o=output, s=script, im=is_image, t=label: self._show(o, s, im, t))
            buttons.addWidget(btn)
            self._buttons.append(btn)

        buttons.addStretch(1)
        save_close = QHBoxLayout()
        self.save_btn = QPushButton("Save", objectName="primary")
        self.save_btn.clicked.connect(self._save)
        self.close_btn = QPushButton("Close", objectName="dangerbtn")
        self.close_btn.clicked.connect(self._close)
        save_close.addWidget(self.save_btn)
        save_close.addWidget(self.close_btn)
        buttons.addLayout(save_close)
        root.addLayout(buttons, 0)

        # --- right: media / image / status pane ---
        self.stack = QStackedWidget()
        self.stack.setObjectName("mediaPane")
        self.status = QLabel("Select a result to view.")
        self.status.setAlignment(Qt.AlignCenter)
        self.status.setWordWrap(True)
        self.image = QLabel()
        self.image.setAlignment(Qt.AlignCenter)
        self.video = QVideoWidget()
        self.stack.addWidget(self.status)   # 0
        self.stack.addWidget(self.image)    # 1
        self.stack.addWidget(self.video)    # 2
        root.addWidget(self.stack, 1)

        self.player = QMediaPlayer(self)
        self.player.setVideoOutput(self.video)
        self.player.mediaStatusChanged.connect(self._loop_video)

    # ---- viewing ----
    def _set_buttons_enabled(self, enabled: bool) -> None:
        for b in self._buttons:
            b.setEnabled(enabled)

    def _show(self, output: str, script_tmpl: str, is_image: bool, title: str) -> None:
        script = script_tmpl.format(c=self.case)
        self.status.setText(f"Preparing “{title}” …")
        self.stack.setCurrentWidget(self.status)
        self._set_buttons_enabled(False)

        self._worker = bf_runner.ParaviewWorker(output, script, self)
        self._worker.ready.connect(lambda p, im=is_image: self._display(p, im))
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(lambda: self._set_buttons_enabled(True))
        self._worker.start()

    def _display(self, path: str, is_image: bool) -> None:
        if is_image:
            self.player.stop()
            pix = QPixmap(path)
            if pix.isNull():
                self.status.setText("Could not load image.")
                self.stack.setCurrentWidget(self.status)
                return
            self.image.setPixmap(pix.scaled(
                self.stack.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation))
            self.stack.setCurrentWidget(self.image)
        else:
            self.stack.setCurrentWidget(self.video)
            self.player.setSource(QUrl.fromLocalFile(path))
            self.player.play()

    def _loop_video(self, status: QMediaPlayer.MediaStatus) -> None:
        if status == QMediaPlayer.EndOfMedia and self.player.source().isValid():
            self.player.setPosition(0)
            self.player.play()

    def _on_failed(self, message: str) -> None:
        self.player.stop()
        self.status.setText(message)
        self.stack.setCurrentWidget(self.status)

    # ---- save / close ----
    def _save(self) -> None:
        target = QFileDialog.getExistingDirectory(self, "Choose a directory to save results")
        if not target:
            return
        self.player.stop()
        bf_runner.save_outputs_to(target)
        self.close()

    def _close(self) -> None:
        reply = QMessageBox.question(
            self, "", "Are you sure you want to continue without saving results?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.player.stop()
            bf_runner.discard_outputs()
            self.close()
        else:
            self._save()

    def closeEvent(self, event):
        # The results window is the last thing on screen after a run — closing it
        # ends the application.
        self.player.stop()
        QApplication.instance().quit()
        super().closeEvent(event)
