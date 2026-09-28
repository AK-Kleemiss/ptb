LANTHANIDES = [
    ("La", 57), ("Ce", 58), ("Pr", 59), ("Nd", 60), ("Pm", 61),
    ("Sm", 62), ("Eu", 63), ("Gd", 64), ("Tb", 65), ("Dy", 66),
    ("Ho", 67), ("Er", 68), ("Tm", 69), ("Yb", 70), ("Lu", 71),
]

SYMBOL_TO_Z = {symbol: z for symbol, z in LANTHANIDES}
Z_TO_SYMBOL = {z: symbol for symbol, z in LANTHANIDES}

# NIST ground-state atomic configurations for La-Lu, reduced to the two
# quantities that matter for ionization: how many 4f electrons the neutral
# atom has, and whether it carries the "5d anomaly" -- La, Ce, Gd, and Lu
# each have an extra unpaired 5d1 electron (half-filled/filled 4f shells at
# Gd/Lu, and an empty/singly-occupied 4f shell at La/Ce, both stabilize a 5d
# electron over doubling up in 4f).
_NEUTRAL_F_COUNT = {
    "La": 0, "Ce": 1, "Pr": 3, "Nd": 4, "Pm": 5, "Sm": 6, "Eu": 7,
    "Gd": 7, "Tb": 9, "Dy": 10, "Ho": 11, "Er": 12, "Tm": 13, "Yb": 14, "Lu": 14,
}
_HAS_5D1 = {"La", "Ce", "Gd", "Lu"}

# Highest oxidation state we're willing to generate structures for, even
# where the electron count would allow more (touching the Xe core beyond
# this is not a meaningful lanthanide oxidation state by any convention).
_MAX_OXIDATION_STATE = 5


def _electron_config_for_oxidation_state(symbol: str, ox: int) -> tuple[int, int, int]:
    """Returns (n_f, d_present, s_unpaired) for Ln^(ox+): electrons are
    removed in the order 6s (2 electrons), then 5d1 if present, then 4f --
    this is the standard picture for lanthanide ionization (4f orbitals are
    core-like and the last to ionize), and reproduces every experimentally
    known oxidation state's electron count (Ce4+ f0, Eu2+ f7, Yb2+ f14,
    Tb4+ f7, ...) as well as the neutral-atom ground-state configurations.
    s_unpaired is 1 only when ox is odd enough to remove exactly one of the
    two 6s electrons (i.e. ox=1) -- every other oxidation state leaves the
    6s shell either full (both paired) or empty (none left), contributing
    no net spin; only the ox=1 case leaves a single unpaired 6s electron,
    which a from a purely n_f/d_present accounting used to silently drop,
    producing electron-count/multiplicity parity violations ORCA rejects
    outright (e.g. Ce+1: 57 explicit electrons, odd, demands an even
    multiplicity, but the old formula gave 3)."""
    remaining = ox
    s_removed = min(remaining, 2)  # 6s
    remaining -= s_removed
    s_unpaired = 1 if s_removed == 1 else 0
    d_removed = 0
    if symbol in _HAS_5D1:
        d_removed = min(remaining, 1)  # 5d1
        remaining -= d_removed
    n0 = _NEUTRAL_F_COUNT[symbol]
    f_removed = min(remaining, n0)  # 4f, only once 6s/5d are exhausted
    n_f = n0 - f_removed
    d_present = (1 if symbol in _HAS_5D1 else 0) - d_removed
    return n_f, d_present, s_unpaired


def _max_oxidation_state(symbol: str) -> int:
    return min(_MAX_OXIDATION_STATE, 2 + (1 if symbol in _HAS_5D1 else 0) + _NEUTRAL_F_COUNT[symbol])


def _multiplicity(n_f: int, d_present: int, s_unpaired: int) -> int:
    # Maximum-multiplicity (high-spin, Hund's-rule) ground state: the first
    # 7 of the 14 4f orbitals fill singly before any pairing starts, the
    # 5d1 electron (if still present at this oxidation state) is an
    # independent unpaired electron on top of that, and a lone remaining
    # 6s electron (only at ox=1) is a third, independent unpaired electron.
    f_unpaired = n_f if n_f <= 7 else 14 - n_f
    return f_unpaired + d_present + s_unpaired + 1


# Every oxidation state from 0 up to each element's maximum (electron-count
# permitting, capped at +5) -- deliberately includes states with no real
# chemical precedent (e.g. Ln0, Ln+1) alongside the well-known ones, since
# PTB parameterization benefits from broad coverage of the computationally
# accessible space, not just experimentally isolable compounds.
COMMON_OXIDATION_STATES: dict[str, list[int]] = {
    symbol: list(range(0, _max_oxidation_state(symbol) + 1))
    for symbol, _ in LANTHANIDES
}

# 2S+1 for every (element, oxidation_state) pair in COMMON_OXIDATION_STATES,
# derived from the electron-removal model above. Validated against known
# spectroscopic ground-state terms at ox=0 (e.g. Gd atom 9D2 -> mult=9,
# La atom 2D3/2 -> mult=2) and every previously-hardcoded ionic state
# (Gd3+ f7 -> mult=8, Ce4+ f0 -> mult=1, Eu2+ f7 -> mult=8, Yb2+ f14 -> mult=1).
DEFAULT_SPIN_MULTIPLICITY: dict[tuple[str, int], int] = {
    (symbol, ox): _multiplicity(*_electron_config_for_oxidation_state(symbol, ox))
    for symbol, _ in LANTHANIDES
    for ox in COMMON_OXIDATION_STATES[symbol]
}
