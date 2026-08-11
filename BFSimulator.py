#!/usr/bin/env python3
"""Brain-Folding Simulator — Adaptive Bento Scientific Dashboard (PySide6).

A single lightweight PySide6 app: a slim header, a responsive bento grid of the
six parameter categories (reflowing 3 → 2 → 1 columns), a contextual Parameter
Guide (right side on wide windows, docked below on narrow ones), and a sticky
action bar. The whole form is built by iterating the data-driven registry in
``bf_fields.py``; styling lives in ``bf_style.py`` and the reusable widgets in
``bf_widgets.py``. All validation, the ``Parameters.prm`` write, and the
build/run/results flow are preserved from the previous implementation.
"""

from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

from PySide6.QtCore import QObject, QSize, Qt, Signal
from PySide6.QtGui import QFontDatabase, QMovie, QPixmap
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QFileDialog, QFrame, QHBoxLayout, QLabel,
    QMainWindow, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QScrollArea,
    QSplitter, QVBoxLayout, QWidget,
)

import bf_fields as F
import bf_style
import bf_widgets as W
from bf_fields import FIELDS, FIELDS_BY_KEY, SECTION_ORDER
from bf_results import ResultsWindow
from bf_runner import BuildRunWorker, discard_outputs, save_outputs_to

ASSET_DIR = Path(__file__).resolve().parent / "Images"


@lru_cache(maxsize=None)
def _pixmap(name: str) -> QPixmap:
    return QPixmap(str(ASSET_DIR / name))


def _info_pixmap(pair) -> QPixmap | None:
    if not pair:
        return None
    light, dark = pair
    dark_mode = bf_style.CURRENT is bf_style.DARK
    pm = _pixmap(dark if dark_mode else light)
    return pm if not pm.isNull() else None


# --- progress dialog --------------------------------------------------------
class ProgressDialog(QDialog):
    stop_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Running simulation")
        self.setModal(True)
        self.resize(1120, 460)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 16, 16, 16)
        lay.setSpacing(12)

        # header: small spinner + title
        head = QHBoxLayout()
        head.setSpacing(10)
        spinner = QLabel()
        self._movie = QMovie(str(ASSET_DIR / "Layer-80.gif"))
        self._movie.setScaledSize(QSize(28, 28))
        spinner.setMovie(self._movie)
        spinner.setFixedSize(28, 28)
        self._movie.start()
        head.addWidget(spinner)
        title = QLabel("Simulation in progress")
        title.setStyleSheet("font-weight:600; font-size:15px;")
        head.addWidget(title)
        head.addStretch(1)
        lay.addLayout(head)

        # progress bar driven by simulation time / total time
        self.bar = QProgressBar(objectName="progressBar")
        self.bar.setRange(0, 0)          # busy until the first timestep arrives
        self.bar.setTextVisible(False)
        self.bar.setFixedHeight(10)
        self.prog_label = QLabel("Preparing…", objectName="appSub")
        lay.addWidget(self.bar)
        lay.addWidget(self.prog_label)

        # live terminal output
        self.console = QPlainTextEdit(objectName="progressLog")
        self.console.setReadOnly(True)
        self.console.setMaximumBlockCount(5000)
        mono = QFontDatabase.systemFont(QFontDatabase.FixedFont)
        mono.setPointSize(11)
        self.console.setFont(mono)
        lay.addWidget(self.console, 1)

        stop = QPushButton("Stop", objectName="dangerbtn")
        stop.clicked.connect(self.stop_requested.emit)
        lay.addWidget(stop, 0, Qt.AlignRight)

    def append_line(self, text: str) -> None:
        self.console.appendPlainText(text)
        sb = self.console.verticalScrollBar()
        sb.setValue(sb.maximum())

    def set_progress(self, current: float, total: float) -> None:
        if total and total > 0:
            pct = max(0, min(100, int(current / total * 100)))
            if self.bar.maximum() == 0:      # leave "busy" mode
                self.bar.setRange(0, 100)
            self.bar.setValue(pct)
            self.prog_label.setText(
                f"Simulation time: {current:g} / {total:g} s   ·   {pct}%")
        else:
            self.prog_label.setText(f"Simulation time: {current:g} s")

    def closeEvent(self, event):
        self._movie.stop()
        super().closeEvent(event)


# --- main window ------------------------------------------------------------
class ParameterWindow(QMainWindow):
    GUIDE_W = 340

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Brain Model Parameters")
        self.setMinimumSize(720, 560)

        self.controls: dict[str, QWidget] = {}
        self.row_widgets: dict[str, tuple[QLabel, QLabel]] = {}
        self.combo_keys: dict[QComboBox, str] = {}
        self.errors: set[str] = set()
        self._active_key: str | None = None
        self._worker: BuildRunWorker | None = None
        self._progress: ProgressDialog | None = None
        self._results: ResultsWindow | None = None

        self.guide = W.ParameterGuide()
        self._build_ui()
        self._apply_enable_logic()
        self._guide_show("vz_raduis")

    # ---- construction ----
    def _build_ui(self) -> None:
        root = QWidget(objectName="root")
        self.setCentralWidget(root)
        v = QVBoxLayout(root)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        v.addWidget(self._header())

        # workspace (scrollable bento grid)
        self.grid = W.ResponsiveCardGrid()
        self._build_cards()
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.scroll.setWidget(self.grid)

        self.splitter = QSplitter(Qt.Horizontal)
        self.splitter.addWidget(self.scroll)
        self.splitter.addWidget(self.guide)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 0)
        self.splitter.setChildrenCollapsible(False)
        self.splitter.setSizes([900, self.GUIDE_W])
        v.addWidget(self.splitter, 1)

        v.addWidget(self._footer())

    def _header(self) -> QWidget:
        bar = QFrame(objectName="header")
        h = QHBoxLayout(bar)
        h.setContentsMargins(16, 9, 16, 9)
        h.setSpacing(12)

        logo = QLabel()
        bpix = _info_pixmap(F.BACKGROUND)
        if bpix and not bpix.isNull():
            logo.setPixmap(bpix.scaledToHeight(34, Qt.SmoothTransformation))
        h.addWidget(logo)

        titles = QVBoxLayout()
        titles.setSpacing(0)
        titles.addWidget(QLabel("Brain Model Parameters", objectName="appTitle"))
        titles.addWidget(QLabel("BRAINIACS · Cortical-folding simulation", objectName="appSub"))
        h.addLayout(titles)
        h.addStretch(1)

        theme = QPushButton("Theme", objectName="hbtn")
        theme.clicked.connect(self._toggle_theme)
        help_btn = QPushButton("Help", objectName="hbtn")
        help_btn.clicked.connect(self._show_help)
        h.addWidget(theme)
        h.addWidget(help_btn)
        return bar

    def _build_cards(self) -> None:
        watcher = self  # eventFilter host for combos
        for section in SECTION_ORDER:
            card = W.ParameterCard(F.SECTION_SHORT[section], F.SECTION_ICON[section],
                                   F.SECTION_ACCENT[section])
            n = 0
            for field in [f for f in FIELDS if f.section == section]:
                control = self._make_control(field)
                lbl, unit = card.add_row(field.label, control)
                self.row_widgets[field.key] = (lbl, unit)
                self.controls[field.key] = control
                n += 1
            card.set_count(n)
            self.grid.add_card(card)

    def _make_control(self, field) -> QWidget:
        accent = F.SECTION_ACCENT[field.section]
        if field.kind == "entry":
            w = W.NumericField(accent)
            w.focused.connect(lambda k=field.key: self._guide_show(k))
            w.edited.connect(lambda _t, k=field.key: self._on_edit(k))
            return w
        if field.kind == "radio":
            w = W.SegmentedControl(field.options, accent)
            w.changed.connect(lambda _v, k=field.key: self._on_seg_changed(k))
            w.focused.connect(lambda k=field.key: self._guide_show(k))
            return w
        # combo
        w = QComboBox(objectName="pcombo")
        w.setProperty("accent", accent)
        w.addItem("", "")
        for disp, val in field.options:
            w.addItem(disp, val)
        w.installEventFilter(self)
        self.combo_keys[w] = field.key
        w.currentIndexChanged.connect(lambda _i, k=field.key: self._on_combo_changed(k))
        return w

    def _footer(self) -> QWidget:
        bar = QFrame(objectName="footer")
        h = QHBoxLayout(bar)
        h.setContentsMargins(16, 10, 16, 10)
        h.setSpacing(8)
        for text, info in [("About", F.ABOUT_PROGRAM), ("Author", F.ABOUT_AUTHOR),
                           ("Copyright", F.COPYRIGHT)]:
            b = QPushButton(text, objectName="linkbtn")
            b.clicked.connect(lambda _=False, i=info, t=text: self._show_static(t, i))
            h.addWidget(b)
        h.addStretch(1)
        for text, slot in [("Restore Defaults", self._defaults_clicked),
                           ("Save Parameters", self._save_clicked),
                           ("Load Parameters", self._load_clicked)]:
            b = QPushButton(text, objectName="actbtn")
            b.clicked.connect(slot)
            h.addWidget(b)
        run = QPushButton("Run Simulation", objectName="primary")
        run.clicked.connect(self._run_clicked)
        h.addWidget(run)
        return bar

    # ---- values ----
    def value(self, key: str) -> str:
        w = self.controls.get(key)
        if isinstance(w, W.NumericField):
            return w.text()
        if isinstance(w, W.SegmentedControl):
            return w.value()
        if isinstance(w, QComboBox):
            return w.currentData() or ""
        return ""

    def all_values(self) -> dict:
        d = {f.key: self.value(f.key) for f in FIELDS}
        d["case"] = self.value("case")
        return d

    def set_value(self, key: str, v: str) -> None:
        w = self.controls.get(key)
        if isinstance(w, W.NumericField):
            w.setText(str(v))
        elif isinstance(w, W.SegmentedControl):
            w.setValue(str(v))
        elif isinstance(w, QComboBox):
            idx = w.findData(str(v))
            if idx >= 0:
                w.setCurrentIndex(idx)

    # ---- guide ----
    def _guide_show(self, key: str) -> None:
        f = FIELDS_BY_KEY.get(key)
        if f is None:
            return
        self._active_key = key
        val = self.value(key)
        info = f.info.get(val) if isinstance(f.info, dict) else f.info
        expl = info.text if info else ""
        img = _info_pixmap(info.img) if (info and info.img) else None
        invalid = key in self.errors
        warn = ""
        if invalid:
            warn = f"Current value “{val}” is outside the accepted range ({f.hint})."
        self.guide.show_field(
            category=F.SECTION_SHORT[f.section], accent=F.SECTION_ACCENT[f.section],
            name=f.label, symbol=f.symbol, unit=f.unit, hint=f.hint,
            explanation=expl, requirement=(f.error_msg or "No special constraints."),
            value=val, image=img, warn=warn)

    def _show_static(self, title: str, info) -> None:
        self._active_key = None
        self.guide.show_static(
            title=title, explanation=info.text, image=_info_pixmap(info.img),
            link_url=info.url or "", link_text=info.url_text or "")

    def _show_help(self) -> None:
        QMessageBox.information(
            self, "Help",
            "Enter model parameters in the cards. Hover the ⓘ next to any parameter — "
            "or focus its field — to see its guide, unit, range, and validation on the "
            "right. Use Restore Defaults for a 2D or 3D preset, then Run Simulation.")

    # ---- change handlers ----
    def _on_edit(self, key: str) -> None:
        self._revalidate()
        if self._active_key == key or self.controls[key].hasFocus():
            self._guide_show(key)

    def _on_seg_changed(self, key: str) -> None:
        self._apply_enable_logic()
        self._revalidate()
        self._guide_show(key)

    def _on_combo_changed(self, key: str) -> None:
        # A combo change is always a deliberate choice — refresh the guide so the
        # per-option figure/text (e.g. OSVZ distribution) updates immediately.
        self._revalidate()
        self._guide_show(key)

    def eventFilter(self, obj, event):
        if event.type() == event.Type.FocusIn and obj in self.combo_keys:
            self._guide_show(self.combo_keys[obj])
        return False

    # ---- enable/disable ----
    def _apply_enable_logic(self) -> None:
        constant = self.value("stiffness_case") == "Constant"
        self._set_field_enabled("max_density", not constant)
        if constant:
            self.controls["max_density"].clear()

        direct = self.value("solver_type") == "Direct"
        self._set_field_enabled("linear_it", not direct)
        if direct:
            self.controls["linear_it"].clear()

    def _set_field_enabled(self, key: str, on: bool) -> None:
        self.controls[key].setEnabled(on)
        lbl, unit = self.row_widgets[key]
        lbl.setEnabled(on)
        if unit is not None:
            unit.setEnabled(on)

    # ---- validation ----
    def _revalidate(self, *_args) -> None:
        ctx = self.all_values()
        for f in FIELDS:
            w = self.controls.get(f.key)
            if not f.validate or not isinstance(w, W.NumericField):
                continue
            ok = f.validate(ctx.get(f.key, ""), ctx)
            if ok:
                self.errors.discard(f.key)
            else:
                self.errors.add(f.key)
            w.set_error(not ok, f.error_msg)

    # ---- defaults ----
    def _defaults_clicked(self) -> None:
        if QMessageBox.question(
                self, "Restore defaults", "Set all values to the recommended defaults?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No) != QMessageBox.Yes:
            return
        box = QMessageBox(self)
        box.setWindowTitle("Choose case")
        box.setText("Which geometry preset?")
        b2 = box.addButton("2D", QMessageBox.AcceptRole)
        b3 = box.addButton("3D", QMessageBox.AcceptRole)
        box.exec()
        self._apply_defaults("3" if box.clickedButton() is b3 else "2")

    def _apply_defaults(self, case: str) -> None:
        for f in FIELDS:
            self.set_value(f.key, f.default_3d if case == "3" else f.default_2d)
        self._apply_enable_logic()
        self._revalidate()
        if self._active_key:
            self._guide_show(self._active_key)

    # ---- save / load parameters ----
    def _save_clicked(self) -> None:
        if self.errors:
            QMessageBox.critical(self, "Invalid values", "Fix the highlighted fields first.")
            return
        target, _ = QFileDialog.getSaveFileName(
            self, "Save parameters", "Parameters.prm", "Parameter files (*.prm)")
        if not target:
            return
        self._write_prm(Path(target))
        QMessageBox.information(self, "Saved", f"Parameters written to\n{target}")

    def _load_clicked(self) -> None:
        src, _ = QFileDialog.getOpenFileName(
            self, "Load parameters", "", "Parameter files (*.prm)")
        if not src:
            return
        self._read_prm(Path(src))
        self._apply_enable_logic()
        self._revalidate()
        QMessageBox.information(self, "Loaded", "Parameters loaded from file.")

    def _read_prm(self, path: Path) -> None:
        """Populate fields from a .prm file (best-effort, by matching set-keys)."""
        try:
            lines = path.read_text().splitlines()
        except OSError as exc:
            QMessageBox.critical(self, "Load failed", str(exc))
            return
        prm_fields = [f for f in FIELDS if f.prm]
        for line in lines:
            for f in prm_fields:
                if f.prm in line and "=" in line:
                    raw = line.split("=", 1)[1].strip()
                    self.set_value(f.key, raw)
                    break

    # ---- run ----
    def _run_clicked(self) -> None:
        missing = []
        for f in FIELDS:
            if not f.required:
                continue
            if f.key == "max_density" and self.value("stiffness_case") == "Constant":
                continue
            if self.value(f.key) == "":
                missing.append(f.key)
        if missing:
            QMessageBox.critical(self, "Missing values", "One or more fields are empty.")
            self._guide_show(missing[0])
            return
        if self.errors:
            QMessageBox.critical(self, "Invalid values", "One or more entered values are not correct.")
            self._guide_show(sorted(self.errors)[0])
            return

        self._write_prm(Path("Parameters.prm"))
        case = self.value("case") or "2"
        self._progress = ProgressDialog(self)
        self._progress.stop_requested.connect(self._on_stop)
        self._worker = BuildRunWorker(case, self.value("total_time"), self)
        self._worker.output.connect(self._progress.append_line)
        self._worker.progress.connect(self._progress.set_progress)
        self._worker.finished_ok.connect(lambda c=case: self._on_run_ok(c))
        self._worker.failed.connect(self._on_run_failed)
        self.hide()
        self._worker.start()
        self._progress.show()

    def _write_prm(self, path: Path) -> None:
        vals = {f.key: self.value(f.key) for f in FIELDS if f.prm}
        if self.value("stiffness_case") == "Constant":
            vals["max_density"] = "700"
        if self.value("solver_type") == "Direct":
            vals["linear_it"] = "100"
        prm_fields = [f for f in FIELDS if f.prm]
        template = Path("Parameters.prm")
        base = template if template.exists() else path
        lines = base.read_text().splitlines(keepends=True)
        out = []
        for line in lines:
            for f in prm_fields:
                if f.prm in line and "=" in line:
                    idx = line.find("=")
                    line = line[:idx + 1] + " " + vals.get(f.key, "") + " \n"
                    break
            out.append(line)
        path.write_text("".join(out))

    def _on_stop(self) -> None:
        if self._worker:
            self._worker.stop()
        if self._progress:
            self._progress.close()
        self.show()

    def _on_run_failed(self, message: str) -> None:
        if self._progress:
            self._progress.close()
        QMessageBox.critical(self, "Simulation error", message)
        self.show()

    def _on_run_ok(self, case: str) -> None:
        if self._progress:
            self._progress.close()
        see = QMessageBox.question(
            self, "Done", "Simulation finished. Do you want to see the results?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        if see == QMessageBox.Yes:
            self._results = ResultsWindow(case)
            self._results.show()
            return
        save = QMessageBox.question(
            self, "Save results", "Do you want to save the results?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        if save == QMessageBox.Yes:
            target = QFileDialog.getExistingDirectory(self, "Choose a directory to save results")
            if target:
                save_outputs_to(target)
        else:
            discard_outputs()
        QApplication.instance().quit()

    # ---- theme + responsive ----
    def _toggle_theme(self) -> None:
        app = QApplication.instance()
        bf_style.CURRENT = bf_style.LIGHT if bf_style.CURRENT is bf_style.DARK else bf_style.DARK
        app.setPalette(bf_style._palette(bf_style.CURRENT))
        app.setStyleSheet(bf_style._qss(bf_style.CURRENT))
        if self._active_key:
            self._guide_show(self._active_key)

    def closeEvent(self, event):
        # Closing the main window ends the app (unless a results window is open,
        # which manages its own lifetime and quits when the user closes it).
        if self._results is None or not self._results.isVisible():
            QApplication.instance().quit()
        super().closeEvent(event)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        below = self.width() < 1080
        want = Qt.Vertical if below else Qt.Horizontal
        if self.splitter.orientation() != want:
            self.splitter.setOrientation(want)
            if below:
                self.splitter.setSizes([self.height(), 300])
            else:
                self.splitter.setSizes([self.width() - self.GUIDE_W, self.GUIDE_W])


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("Brain Model Parameters")
    # Don't tear the app down when a transient window (progress dialog) closes
    # while the main window is hidden during a run — otherwise the results
    # window would close the instant it opens. We quit explicitly instead.
    app.setQuitOnLastWindowClosed(False)
    bf_style.apply(app)
    win = ParameterWindow()
    win.resize(1480, 750)
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
