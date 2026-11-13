"""Critical constants of common fluids (SI units).

Values are the usual tabulated ones (for example from the NIST Chemistry WebBook
and Poling, Prausnitz and O'Connell), rounded.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Species", "SPECIES", "get"]

R = 8.314462618  # J / (mol K)


@dataclass(frozen=True)
class Species:
    name: str
    Tc: float      # K
    Pc: float      # Pa
    omega: float   # acentric factor
    M: float       # kg / mol


SPECIES = {
    "N2": Species("nitrogen", 126.20, 3.3958e6, 0.0372, 28.0134e-3),
    "O2": Species("oxygen", 154.58, 5.043e6, 0.0222, 31.9988e-3),
    "Ar": Species("argon", 150.69, 4.863e6, -0.0022, 39.948e-3),
    "CO2": Species("carbon dioxide", 304.13, 7.3773e6, 0.2239, 44.0095e-3),
    "CH4": Species("methane", 190.56, 4.5992e6, 0.0114, 16.0425e-3),
    "C2H6": Species("ethane", 305.32, 4.872e6, 0.0995, 30.069e-3),
    "C3H8": Species("propane", 369.83, 4.248e6, 0.1521, 44.0956e-3),
    "H2O": Species("water", 647.10, 22.064e6, 0.3443, 18.0153e-3),
    "NH3": Species("ammonia", 405.40, 11.333e6, 0.2560, 17.0305e-3),
}


def get(formula: str) -> Species:
    try:
        return SPECIES[formula]
    except KeyError:
        raise KeyError(f"unknown species {formula!r}; available: {sorted(SPECIES)}") from None
