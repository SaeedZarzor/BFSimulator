"""Reusable Bento-dashboard components for the Brain-Folding Simulator.

  * ``SegmentedControl`` — replaces radio groups (2D/3D, Varying/Constant, Direct/CG)
  * ``NumericField``     — a spin-box-like line edit with steppers + inline error,
                           preserving scientific notation (1.0e-8, 4.7e-4)
  * ``InfoDot``          — the small "i" that drives the Parameter Guide
  * ``ParameterCard``    — a category card with icon, title, and an aligned row grid
  * ``ParameterGuide``   — the contextual right/bottom guide panel
  * ``ResponsiveCardGrid`` — reflows cards 3 → 2 → 1 columns with width

All styling lives in bf_style.py; these widgets only set object names and the
``accent`` / ``hot`` / ``err`` properties the stylesheet keys off.
"""

from __future__ import annotations

from PySide6.QtCore import QByteArray, QSize, Qt, Signal
from PySide6.QtGui import QPixmap
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (
    QButtonGroup, QComboBox, QFrame, QGridLayout, QHBoxLayout, QLabel,
    QLineEdit, QPushButton, QSizePolicy, QToolButton, QVBoxLayout, QWidget,
)

import bf_fields as F
import bf_style


# --- icon helper ------------------------------------------------------------
def svg_icon(inner_paths: str, color: str, size: int = 18) -> QPixmap:
    """Render inline SVG stroke paths to a crisp (2x) QPixmap tinted `color`."""
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" '
           f'fill="none" stroke="{color}" stroke-width="1.8" '
           f'stroke-linecap="round" stroke-linejoin="round">{inner_paths}</svg>')
    renderer = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    dpr = 2
    pm = QPixmap(size * dpr, size * dpr)
    pm.fill(Qt.transparent)
    from PySide6.QtGui import QPainter
    p = QPainter(pm)
    renderer.render(p)
    p.end()
    pm.setDevicePixelRatio(dpr)
    return pm


def _repolish(w: QWidget) -> None:
    w.style().unpolish(w)
    w.style().polish(w)


# --- info dot ---------------------------------------------------------------
class InfoDot(QToolButton):
    requested = Signal(str)   # field key, on hover or click

    def __init__(self, key: str):
        super().__init__()
        self.setObjectName("infoDot")
        self.key = key
        self.setText("i")
        self.setFixedSize(16, 16)
        self.setCursor(Qt.WhatsThisCursor)
        self.setFocusPolicy(Qt.NoFocus)
        self.clicked.connect(lambda: self.requested.emit(self.key))

    def enterEvent(self, event):
        self.requested.emit(self.key)
        super().enterEvent(event)


# --- numeric field ----------------------------------------------------------
class NumericField(QWidget):
    focused = Signal()
    edited = Signal(str)

    def __init__(self, accent: str = "geo"):
        super().__init__()
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(2)

        self._box = QFrame()
        self._box.setObjectName("numfield")
        self._box.setProperty("accent", accent)
        self._box.setFixedHeight(30)
        row = QHBoxLayout(self._box)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)

        self._edit = QLineEdit()
        self._edit.installEventFilter(self)
        self._edit.textChanged.connect(self.edited.emit)
        row.addWidget(self._edit, 1)

        steps = QVBoxLayout()
        steps.setContentsMargins(0, 0, 0, 0)
        steps.setSpacing(0)
        self._up = QToolButton(objectName="stepUp", text="▲")
        self._dn = QToolButton(objectName="stepDown", text="▼")
        for b, d in ((self._up, 1), (self._dn, -1)):
            b.setFocusPolicy(Qt.NoFocus)
            b.setFixedWidth(20)
            b.clicked.connect(lambda _=False, dd=d: self._step(dd))
        steps.addWidget(self._up)
        steps.addWidget(self._dn)
        row.addLayout(steps)

        lay.addWidget(self._box)
        self._err = QLabel(objectName="fieldError")
        self._err.setWordWrap(True)
        self._err.setVisible(False)
        lay.addWidget(self._err)

    # value api
    def text(self) -> str:
        return self._edit.text()

    def setText(self, v) -> None:
        self._edit.setText(str(v))

    def clear(self) -> None:
        self._edit.clear()

    def setEnabled(self, on: bool) -> None:
        self._edit.setEnabled(on)
        self._up.setEnabled(on)
        self._dn.setEnabled(on)

    def isEnabled(self) -> bool:
        return self._edit.isEnabled()

    def lineEdit(self) -> QLineEdit:
        return self._edit

    def set_accent(self, key: str) -> None:
        self._box.setProperty("accent", key)
        _repolish(self._box)

    def set_error(self, is_error: bool, message: str = "") -> None:
        if self._box.property("err") != is_error:
            self._box.setProperty("err", is_error)
            _repolish(self._box)
        self._err.setText(message)
        self._err.setVisible(bool(is_error and message))

    def eventFilter(self, obj, event):
        if obj is self._edit:
            if event.type() == event.Type.FocusIn:
                self._box.setProperty("hot", True)
                _repolish(self._box)
                self.focused.emit()
            elif event.type() == event.Type.FocusOut:
                self._box.setProperty("hot", False)
                _repolish(self._box)
        return False

    def _step(self, direction: int) -> None:
        t = self._edit.text().strip()
        if not t or "e" in t.lower():
            return
        try:
            v = float(t)
        except ValueError:
            return
        dec = len(t.split(".")[1]) if "." in t else 0
        dec = min(dec, 6)
        step = 10 ** (-dec) if dec else 1
        v += direction * step
        self._edit.setText((f"%.{dec}f" % v) if dec else str(int(round(v))))


# --- segmented control ------------------------------------------------------
class SegmentedControl(QFrame):
    changed = Signal(str)
    focused = Signal()

    def __init__(self, options: list[tuple[str, str]], accent: str = "geo"):
        super().__init__()
        self.setObjectName("seg")
        self.setFixedHeight(30)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.setSpacing(2)
        self._group = QButtonGroup(self)
        self._group.setExclusive(True)
        self._buttons: dict[str, QPushButton] = {}
        for disp, val in options:
            b = QPushButton(disp)
            b.setObjectName("segbtn")
            b.setProperty("accent", accent)
            b.setCheckable(True)
            b.setCursor(Qt.PointingHandCursor)
            b.installEventFilter(self)
            self._group.addButton(b)
            self._buttons[val] = b
            lay.addWidget(b)
            b.clicked.connect(lambda _=False, v=val: self._on_click(v))

    def _on_click(self, val: str) -> None:
        self.changed.emit(val)

    def value(self) -> str:
        for val, b in self._buttons.items():
            if b.isChecked():
                return val
        return ""

    def setValue(self, val: str) -> None:
        b = self._buttons.get(val)
        if b:
            b.setChecked(True)

    def set_accent(self, key: str) -> None:
        for b in self._buttons.values():
            b.setProperty("accent", key)
            _repolish(b)

    def eventFilter(self, obj, event):
        if event.type() == event.Type.FocusIn:
            self.focused.emit()
        return False


# --- category card ----------------------------------------------------------
class ParameterCard(QFrame):
    def __init__(self, title: str, icon_paths: str, accent: str):
        super().__init__()
        self.setObjectName("pcard")
        self.setProperty("accent", accent)
        # Minimum (not Maximum) vertical policy: the card's sizeHint is a floor,
        # so masonry columns never compress rows — they scroll instead.
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(14, 12, 14, 14)
        outer.setSpacing(11)

        head = QHBoxLayout()
        head.setSpacing(9)
        icon = QLabel()
        icon.setPixmap(svg_icon(icon_paths, bf_style.accent_hex(accent), 20))
        icon.setFixedSize(24, 24)
        icon.setAlignment(Qt.AlignCenter)
        head.addWidget(icon)
        title_lbl = QLabel(title, objectName="cardTitle")
        title_lbl.setProperty("accent", accent)   # colored per category via QSS
        head.addWidget(title_lbl)
        head.addStretch(1)
        self._count = QLabel(objectName="cardCount")
        head.addWidget(self._count)
        outer.addLayout(head)

        self._grid = QGridLayout()
        self._grid.setHorizontalSpacing(10)
        self._grid.setVerticalSpacing(9)
        self._grid.setColumnStretch(0, 1)
        self._grid.setColumnMinimumWidth(1, 130)
        outer.addLayout(self._grid)
        self._row = 0

    def add_row(self, label: str, control: QWidget, unit: str = ""):
        # No info dot and no unit column here — the Parameter Guide reports both,
        # updating on field focus.
        lbl = QLabel(label, objectName="pLabel")
        lbl.setWordWrap(True)
        lbl.setTextInteractionFlags(Qt.NoTextInteraction)
        lbl.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self._grid.addWidget(lbl, self._row, 0, Qt.AlignVCenter)
        self._grid.addWidget(control, self._row, 1, Qt.AlignVCenter)
        self._row += 1
        return lbl, None

    def set_count(self, n: int) -> None:
        self._count.setText(str(n))

    def weight(self) -> int:
        """Rough height proxy (row count) used to balance the masonry columns."""
        return max(self._row, 1)


# --- figure that preserves aspect ratio -------------------------------------
class FigureLabel(QLabel):
    def __init__(self):
        super().__init__()
        self.setAlignment(Qt.AlignCenter)
        self.setMinimumHeight(150)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        self._src: QPixmap | None = None

    def set_source(self, pm: QPixmap | None) -> None:
        self._src = pm if (pm and not pm.isNull()) else None
        self._rescale()

    def resizeEvent(self, event):
        self._rescale()
        super().resizeEvent(event)

    def _rescale(self) -> None:
        if self._src is None:
            self.clear()
            return
        target = self.size() * self.devicePixelRatio() if False else self.size()
        scaled = self._src.scaled(
            max(1, self.width() - 4), max(1, self.height() - 4),
            Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.setPixmap(scaled)


# --- parameter guide --------------------------------------------------------
class ParameterGuide(QFrame):
    def __init__(self):
        super().__init__()
        self.setObjectName("guidePanel")
        self.setMinimumWidth(280)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 16, 16, 16)
        lay.setSpacing(12)

        # panel header: open-book icon + "Parameter Guide"
        header = QHBoxLayout()
        header.setSpacing(8)
        book = QLabel()
        book.setPixmap(svg_icon(F.GUIDE_ICON, bf_style.token("muted"), 18))
        book.setFixedSize(20, 20)
        header.addWidget(book)
        header.addWidget(QLabel("Parameter Guide", objectName="gPanelTitle"))
        header.addStretch(1)
        lay.addLayout(header)

        self.eyebrow = QLabel("", objectName="gEyebrow")
        self.name = QLabel("", objectName="gName")
        self.name.setWordWrap(True)
        lay.addWidget(self.eyebrow)
        lay.addWidget(self.name)

        figwrap = QFrame(objectName="gFig")
        figwrap.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        fl = QVBoxLayout(figwrap)
        fl.setContentsMargins(12, 12, 12, 12)
        self.figure = FigureLabel()
        fl.addWidget(self.figure)
        self._figwrap = figwrap
        lay.addWidget(figwrap, 4)   # figure box absorbs most of the vertical slack

        self.expl = QLabel("", objectName="gExpl")
        self.expl.setWordWrap(True)
        lay.addWidget(self.expl)

        self._meta = QWidget()          # holds the Symbol/Unit/Recommended/Current cells
        meta = QGridLayout(self._meta)
        meta.setContentsMargins(0, 0, 0, 0)
        meta.setHorizontalSpacing(8)
        meta.setVerticalSpacing(8)
        self._cells = {}
        for i, (k, key) in enumerate([("Symbol", "sym"), ("Unit", "unit"),
                                      ("Recommended", "range"), ("Current value", "val")]):
            cell = QFrame(objectName="gCell")
            cl = QVBoxLayout(cell)
            cl.setContentsMargins(9, 7, 9, 7)
            cl.setSpacing(1)
            cl.addWidget(QLabel(k, objectName="gCellK"))
            v = QLabel("—", objectName="gCellV")
            v.setTextFormat(Qt.RichText)
            cl.addWidget(v)
            self._cells[key] = v
            meta.addWidget(cell, i // 2, i % 2)
        lay.addWidget(self._meta)

        self.valid = QLabel("", objectName="gValid")
        self.valid.setWordWrap(True)
        lay.addWidget(self.valid)

        self.warn = QFrame(objectName="gWarn")
        wl = QHBoxLayout(self.warn)
        wl.setContentsMargins(11, 9, 11, 9)
        self.warn_text = QLabel("", objectName="gWarnText")
        self.warn_text.setWordWrap(True)
        wl.addWidget(self.warn_text)
        self.warn.setVisible(False)
        lay.addWidget(self.warn)
        lay.addStretch(1)

    def show_field(self, *, category: str, accent: str, name: str, symbol: str,
                   unit: str, hint: str, explanation: str, requirement: str,
                   value: str, image: QPixmap | None, warn: str = "") -> None:
        self.eyebrow.setText(category.upper())
        self.eyebrow.setStyleSheet(f"color:{bf_style.accent_hex(accent)};")
        self.name.setText(name)
        self.figure.set_source(image)
        self._figwrap.setVisible(image is not None)
        self.expl.setText(explanation or "")
        self.expl.setVisible(bool(explanation))
        self._meta.setVisible(True)
        self._cells["sym"].setText(latex_html(symbol))
        self._cells["unit"].setText(unit or "—")
        self._cells["range"].setText(hint or "—")
        self._cells["val"].setText(value if value else "—")
        self.valid.setOpenExternalLinks(False)
        self.valid.setText("✓  " + (requirement or "No special constraints."))
        self.valid.setVisible(True)
        if warn:
            self.warn_text.setText(warn)
            self.warn.setVisible(True)
        else:
            self.warn.setVisible(False)

    def show_static(self, *, title: str, explanation: str, image: QPixmap | None,
                    link_url: str = "", link_text: str = "") -> None:
        """About / Author / Copyright view — no Symbol/Unit/Recommended/Current cells."""
        self.eyebrow.setText("INFORMATION")
        self.eyebrow.setStyleSheet(f"color:{bf_style.token('faint')};")
        self.name.setText(title)
        self.figure.set_source(image)
        self._figwrap.setVisible(image is not None)
        self.expl.setText(explanation or "")
        self.expl.setVisible(bool(explanation))
        self._meta.setVisible(False)
        self.warn.setVisible(False)
        if link_url:
            self.valid.setOpenExternalLinks(True)
            self.valid.setText(f'<a href="{link_url}">{link_text or link_url}</a>')
            self.valid.setVisible(True)
        else:
            self.valid.setVisible(False)


_GREEK = {
    r"\alpha": "α", r"\beta": "β", r"\gamma": "γ", r"\delta": "δ", r"\Delta": "Δ",
    r"\epsilon": "ε", r"\theta": "θ", r"\kappa": "κ", r"\lambda": "λ", r"\mu": "μ",
    r"\nu": "ν", r"\pi": "π", r"\rho": "ρ", r"\sigma": "σ", r"\tau": "τ",
    r"\phi": "φ", r"\omega": "ω",
}


def latex_html(symbol: str) -> str:
    """Render a small LaTeX subset (\\greek, x_{..}, x^{..}) as math-styled rich text.

    Not a full TeX engine — just enough for parameter symbols like ``r_{vz}``,
    ``\\Delta t``, ``c_{max}`` or ``k_s``. Output is italic serif with real
    sub/superscripts, giving a LaTeX-like look without any external dependency.
    """
    if not symbol or symbol in ("—",):
        return "—"
    s = symbol
    for cmd, ch in _GREEK.items():
        s = s.replace(cmd, ch)

    out, i = [], 0
    while i < len(s):
        c = s[i]
        if c in "_^":
            tag = "sub" if c == "_" else "sup"
            i += 1
            if i < len(s) and s[i] == "{":
                j = s.find("}", i)
                content = s[i + 1:j] if j != -1 else s[i + 1:]
                i = j + 1 if j != -1 else len(s)
            elif i < len(s):
                content, i = s[i], i + 1
            else:
                content = ""
            out.append(f"<{tag}>{content}</{tag}>")
        else:
            out.append(c)
            i += 1
    body = "".join(out)
    return (f'<span style="font-family:Georgia,\'Times New Roman\',serif;'
            f'font-style:italic;font-size:14px">{body}</span>')


# --- responsive card grid ---------------------------------------------------
class ResponsiveCardGrid(QWidget):
    """Masonry board: packs cards into 3 → 2 → 1 balanced columns with width.

    Each card is placed in the currently-shortest column (greedy by row count),
    so tall and short cards interleave without the row-height gaps a plain grid
    leaves under shorter cards.
    """

    def __init__(self):
        super().__init__()
        self.setObjectName("workspace")
        self._outer = QHBoxLayout(self)
        self._outer.setContentsMargins(16, 16, 16, 12)
        self._outer.setSpacing(14)
        self._outer.setAlignment(Qt.AlignTop)
        self._cards: list[ParameterCard] = []
        self._cols = 0

    def add_card(self, card: ParameterCard) -> None:
        self._cards.append(card)
        self._relayout(force=True)

    def _cols_for(self, w: int) -> int:
        if w >= 960:
            return 3
        if w >= 640:
            return 2
        return 1

    def resizeEvent(self, event):
        self._relayout()
        super().resizeEvent(event)

    def _relayout(self, force: bool = False) -> None:
        cols = self._cols_for(self.width())
        if cols == self._cols and not force:
            return
        self._cols = cols

        # detach cards, then drop the old column containers
        for card in self._cards:
            card.setParent(None)
        while self._outer.count():
            item = self._outer.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()

        columns: list[QVBoxLayout] = []
        heights = [0] * cols
        for _ in range(cols):
            holder = QWidget()
            col = QVBoxLayout(holder)
            col.setContentsMargins(0, 0, 0, 0)
            col.setSpacing(14)
            self._outer.addWidget(holder, 1)
            columns.append(col)

        for card in self._cards:
            i = heights.index(min(heights))
            columns[i].addWidget(card)
            heights[i] += card.weight()
        for col in columns:
            col.addStretch(1)
