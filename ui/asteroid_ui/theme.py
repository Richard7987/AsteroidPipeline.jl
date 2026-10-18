"""Tema Gruvbox oscuro — la misma paleta que nezzontli-ctl.

Texto claro (#ebdbb2) sobre fondo #282828: contraste ~11:1, muy por encima del
4.5:1 que piden las pautas de accesibilidad (WCAG AA).
"""

from textual.theme import Theme

GRUVBOX = Theme(
    name="gruvbox-asteroides",
    primary="#83a598",
    secondary="#d3869b",
    accent="#fabd2f",
    warning="#fabd2f",
    error="#fb4934",
    success="#b8bb26",
    foreground="#ebdbb2",
    background="#282828",
    surface="#3c3836",
    panel="#504945",
    dark=True,
)
