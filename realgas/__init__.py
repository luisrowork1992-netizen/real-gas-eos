"""realgas: cubic equations of state (van der Waals, RK, SRK, Peng-Robinson) for pure fluids."""

from .cubic import CubicEOS, PengRobinson, RedlichKwong, SoaveRedlichKwong, VanDerWaals
from .species import R, SPECIES, Species, get

EOS = {"vdW": VanDerWaals, "RK": RedlichKwong, "SRK": SoaveRedlichKwong, "PR": PengRobinson}

__version__ = "0.1.0"
