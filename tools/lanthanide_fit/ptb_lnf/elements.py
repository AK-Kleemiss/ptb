LANTHANIDES = [
    ("La", 57), ("Ce", 58), ("Pr", 59), ("Nd", 60), ("Pm", 61),
    ("Sm", 62), ("Eu", 63), ("Gd", 64), ("Tb", 65), ("Dy", 66),
    ("Ho", 67), ("Er", 68), ("Tm", 69), ("Yb", 70), ("Lu", 71),
]

SYMBOL_TO_Z = {symbol: z for symbol, z in LANTHANIDES}
Z_TO_SYMBOL = {z: symbol for symbol, z in LANTHANIDES}

COMMON_OXIDATION_STATES = {
    "La": [3],
    "Ce": [3, 4],
    "Pr": [3, 4],
    "Nd": [3],
    "Pm": [3],
    "Sm": [2, 3],
    "Eu": [2, 3],
    "Gd": [3],
    "Tb": [3, 4],
    "Dy": [3],
    "Ho": [3],
    "Er": [3],
    "Tm": [2, 3],
    "Yb": [2, 3],
    "Lu": [3],
}

DEFAULT_SPIN_MULTIPLICITY = {
    ("La", 3): 1, ("Ce", 3): 2, ("Ce", 4): 1, ("Pr", 3): 3, ("Pr", 4): 2,
    ("Nd", 3): 4, ("Pm", 3): 5, ("Sm", 2): 7, ("Sm", 3): 6,
    ("Eu", 2): 8, ("Eu", 3): 7, ("Gd", 3): 8, ("Tb", 3): 7,
    ("Tb", 4): 6, ("Dy", 3): 6, ("Ho", 3): 5, ("Er", 3): 4,
    ("Tm", 2): 4, ("Tm", 3): 3, ("Yb", 2): 1, ("Yb", 3): 2,
    ("Lu", 3): 1,
}
