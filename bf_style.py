"""Adaptive Bento Scientific Dashboard theme for the Brain-Folding Simulator.

A single centralized stylesheet (Fusion base + tuned palette + QSS) drives the
whole UI: cool-gray ground, white cards, navy text, deep-teal primary action,
and six restrained category accents used only on card top-borders, icons, focus
rings, and the segmented selection. Components opt in via object names and the
dynamic properties ``accent`` / ``hot`` / ``err`` — no per-widget inline styles.

``apply(app)`` installs everything and stores the active tokens; components read
colors back through ``token()`` and ``accent_hex()``.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont, QPalette
from PySide6.QtWidgets import QApplication

ACCENT_KEYS = ["geo", "adv", "mech", "disc", "solv", "grow"]

LIGHT = {
    "bg": "#f4f6f9", "panel": "#eef1f6", "card": "#ffffff", "card_alt": "#fbfcfe",
    "text": "#172033", "muted": "#64748b", "faint": "#94a3b8",
    "border": "#e4e8ef", "border_strong": "#d3dae6",
    "primary": "#0f766e", "primary_ink": "#ffffff",
    "danger": "#dc2626", "danger_tint": "#fdeeee",
    "seg_ink": "#ffffff",
    "accents": {"geo": "#e11d48", "adv": "#0891b2", "mech": "#7c3aed",
                "disc": "#d97706", "solv": "#64748b", "grow": "#16a34a"},
}
DARK = {
    "bg": "#0e131b", "panel": "#131a24", "card": "#171f2b", "card_alt": "#141c27",
    "text": "#e6ecf4", "muted": "#93a1b5", "faint": "#6b7a90",
    "border": "#26303e", "border_strong": "#38465a",
    "primary": "#2dd4bf", "primary_ink": "#07231f",
    "danger": "#f87171", "danger_tint": "#3a1f22",
    "seg_ink": "#0a1622",
    "accents": {"geo": "#fb7185", "adv": "#22d3ee", "mech": "#a78bfa",
                "disc": "#fbbf24", "solv": "#94a3b8", "grow": "#4ade80"},
}

CURRENT = LIGHT  # replaced on apply()


def token(name: str) -> str:
    return CURRENT.get(name, "#000000")


def accent_hex(key: str) -> str:
    return CURRENT["accents"].get(key, CURRENT["primary"])


def is_dark(app: QApplication) -> bool:
    try:
        return app.styleHints().colorScheme() == Qt.ColorScheme.Dark
    except AttributeError:
        return app.palette().color(QPalette.Window).lightness() < 128


def _accent_rules(c: dict) -> str:
    """Per-category rules for card top-border, focus rings, and segmented selection."""
    out = []
    ink = c["seg_ink"]
    for k, hexv in c["accents"].items():
        out.append(f'QFrame#pcard[accent="{k}"] {{ border-top-color: {hexv}; }}')
        out.append(f'QLabel#cardTitle[accent="{k}"] {{ color: {hexv}; }}')
        out.append(f'QFrame#numfield[accent="{k}"][hot="true"] {{ border-color: {hexv}; }}')
        out.append(f'QComboBox#pcombo[accent="{k}"]:focus {{ border-color: {hexv}; }}')
        out.append(f'QPushButton#segbtn[accent="{k}"]:checked {{ background: {hexv}; color: {ink}; }}')
    return "\n".join(out)


def _qss(c: dict) -> str:
    return f"""
    QWidget {{ color: {c['text']}; font-size: 13px; }}
    QMainWindow, #root {{ background: {c['bg']}; }}
    QToolTip {{ background: {c['card']}; color: {c['text']}; border: 1px solid {c['border']}; padding: 4px 6px; }}

    /* scroll area / workspace */
    QScrollArea {{ border: 0; background: transparent; }}
    #workspace {{ background: {c['bg']}; }}
    QScrollBar:vertical {{ background: transparent; width: 11px; margin: 2px; }}
    QScrollBar::handle:vertical {{ background: {c['border_strong']}; border-radius: 5px; min-height: 32px; }}
    QScrollBar::handle:vertical:hover {{ background: {c['muted']}; }}
    QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; }}
    QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}

    /* header */
    #header {{ background: {c['card']}; border-bottom: 1px solid {c['border']}; }}
    #appTitle {{ font-size: 15px; font-weight: 600; }}
    #appSub {{ font-size: 11px; color: {c['muted']}; }}
    QPushButton#hbtn {{
        background: {c['card_alt']}; color: {c['text']}; border: 1px solid {c['border_strong']};
        border-radius: 8px; padding: 6px 11px; font-size: 12px;
    }}
    QPushButton#hbtn:hover {{ border-color: {c['primary']}; color: {c['primary']}; }}

    /* cards */
    QFrame#pcard {{
        background: {c['card']}; border: 1px solid {c['border']};
        border-top-width: 3px; border-radius: 12px;
    }}
    QLabel#cardTitle {{ font-size: 15px; font-weight: 700; letter-spacing: 0.2px; }}
    QLabel#cardCount {{ color: {c['faint']}; font-size: 11px; }}
    QLabel#pLabel {{ font-size: 12px; color: {c['text']}; }}
    QLabel#pUnit {{ font-size: 11px; color: {c['muted']}; }}
    QLabel#pLabel:disabled, QLabel#pUnit:disabled {{ color: {c['faint']}; }}

    /* info dot */
    QToolButton#infoDot {{
        border: 1px solid {c['border_strong']}; border-radius: 8px;
        background: {c['card_alt']}; color: {c['muted']};
        font-style: italic; font-weight: 700; font-size: 9px; padding: 0;
    }}
    QToolButton#infoDot:hover {{ border-color: {c['primary']}; color: {c['primary']}; }}

    /* numeric field (line edit + steppers) */
    QFrame#numfield {{
        background: {c['card_alt']}; border: 1px solid {c['border_strong']}; border-radius: 7px;
    }}
    QFrame#numfield[err="true"] {{ border-color: {c['danger']}; background: {c['danger_tint']}; }}
    QFrame#numfield QLineEdit {{
        border: 0; background: transparent; color: {c['text']}; padding: 0 8px; font-size: 12px;
    }}
    QFrame#numfield QLineEdit:disabled {{ color: {c['faint']}; }}
    QToolButton#stepUp, QToolButton#stepDown {{
        background: {c['panel']}; color: {c['muted']}; border: 0; border-left: 1px solid {c['border']};
        font-size: 8px; padding: 0; width: 20px;
    }}
    QToolButton#stepDown {{ border-top: 1px solid {c['border']}; }}
    QToolButton#stepUp:hover, QToolButton#stepDown:hover {{ color: {c['primary']}; }}
    QLabel#fieldError {{ color: {c['danger']}; font-size: 10.5px; }}

    /* segmented control */
    QFrame#seg {{ background: {c['card_alt']}; border: 1px solid {c['border_strong']}; border-radius: 7px; }}
    QPushButton#segbtn {{
        border: 0; background: transparent; color: {c['muted']};
        border-radius: 5px; font-size: 12px; padding: 4px 8px;
    }}
    QPushButton#segbtn:hover:!checked {{ color: {c['text']}; }}

    /* combo box */
    QComboBox#pcombo {{
        background: {c['card_alt']}; border: 1px solid {c['border_strong']}; border-radius: 7px;
        padding: 3px 8px; color: {c['text']}; font-size: 12px; min-height: 22px;
    }}
    QComboBox#pcombo::drop-down {{ border: 0; width: 20px; }}
    QComboBox QAbstractItemView {{
        background: {c['card']}; border: 1px solid {c['border']};
        selection-background-color: {c['primary']}; selection-color: #ffffff; outline: none;
    }}

    /* guide panel */
    QFrame#guidePanel {{ background: {c['card']}; border-left: 1px solid {c['border']}; }}
    QLabel#gPanelTitle {{ font-size: 13px; font-weight: 700; color: {c['muted']}; }}
    QLabel#gEyebrow {{ font-size: 10px; font-weight: 600; color: {c['faint']}; }}
    QLabel#gName {{ font-size: 15px; font-weight: 600; }}
    QLabel#gExpl {{ color: {c['muted']}; font-size: 12px; }}
    QFrame#gFig {{ background: {c['card_alt']}; border: 1px solid {c['border']}; border-radius: 10px; }}
    QFrame#gCell {{ background: {c['card_alt']}; border: 1px solid {c['border']}; border-radius: 8px; }}
    QLabel#gCellK {{ font-size: 10px; color: {c['faint']}; }}
    QLabel#gCellV {{ font-size: 12.5px; font-weight: 600; }}
    QLabel#gValid {{ color: {c['muted']}; font-size: 11.5px; }}
    QFrame#gWarn {{ background: {c['danger_tint']}; border: 1px solid {c['danger']}; border-radius: 9px; }}
    QLabel#gWarnText {{ color: {c['danger']}; font-size: 11.5px; }}

    /* footer / action bar */
    #footer {{ background: {c['card']}; border-top: 1px solid {c['border']}; }}
    QPushButton#actbtn {{
        background: {c['card_alt']}; color: {c['text']}; border: 1px solid {c['border_strong']};
        border-radius: 8px; padding: 7px 13px; font-size: 12px;
    }}
    QPushButton#actbtn:hover {{ border-color: {c['primary']}; color: {c['primary']}; }}
    QPushButton#linkbtn {{ background: transparent; border: 0; color: {c['muted']}; padding: 6px 8px; font-size: 12px; }}
    QPushButton#linkbtn:hover {{ color: {c['text']}; }}
    QPushButton#primary {{
        background: {c['primary']}; color: {c['primary_ink']}; border: 0;
        border-radius: 8px; padding: 7px 18px; font-size: 12.5px; font-weight: 600;
    }}
    QPushButton#primary:hover {{ background: {c['primary']}; }}
    QPushButton#dangerbtn {{
        background: {c['danger']}; color: #ffffff; border: 0;
        border-radius: 8px; padding: 7px 16px; font-size: 12.5px; font-weight: 600;
    }}

    QDialog {{ background: {c['bg']}; }}
    QPlainTextEdit#progressLog {{
        background: {c['card_alt']}; color: {c['text']};
        border: 1px solid {c['border']}; border-radius: 8px; padding: 8px;
    }}
    QProgressBar#progressBar {{
        background: {c['card_alt']}; border: 1px solid {c['border']};
        border-radius: 5px;
    }}
    QProgressBar#progressBar::chunk {{
        background: {c['primary']}; border-radius: 5px;
    }}
    QStackedWidget#mediaPane {{
        background: {c['card_alt']}; border: 1px solid {c['border']}; border-radius: 12px;
    }}
    {_accent_rules(c)}
    """


def _palette(c: dict) -> QPalette:
    p = QPalette()
    p.setColor(QPalette.Window, QColor(c["bg"]))
    p.setColor(QPalette.Base, QColor(c["card_alt"]))
    p.setColor(QPalette.AlternateBase, QColor(c["panel"]))
    p.setColor(QPalette.Text, QColor(c["text"]))
    p.setColor(QPalette.WindowText, QColor(c["text"]))
    p.setColor(QPalette.Button, QColor(c["card"]))
    p.setColor(QPalette.ButtonText, QColor(c["text"]))
    p.setColor(QPalette.Highlight, QColor(c["primary"]))
    p.setColor(QPalette.HighlightedText, QColor("#ffffff"))
    p.setColor(QPalette.ToolTipBase, QColor(c["card"]))
    p.setColor(QPalette.ToolTipText, QColor(c["text"]))
    p.setColor(QPalette.PlaceholderText, QColor(c["faint"]))
    return p


def apply(app: QApplication) -> dict:
    """Install the Bento theme; return the active token dict."""
    global CURRENT
    app.setStyle("Fusion")
    font = QFont(app.font())
    font.setPointSizeF(max(font.pointSizeF(), 12.5))
    app.setFont(font)

    CURRENT = DARK if is_dark(app) else LIGHT
    app.setPalette(_palette(CURRENT))
    app.setStyleSheet(_qss(CURRENT))
    return CURRENT
