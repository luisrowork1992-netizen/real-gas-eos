"""Generic two-parameter cubic equations of state.

All of them share one form::

    P = R T / (V - b) - a(T) / ((V + eps b) (V + sig b))

| EOS            | eps        | sig        | a(T)                          |
|----------------|------------|------------|-------------------------------|
| van der Waals  | 0          | 0          | constant                      |
| Redlich-Kwong  | 0          | 1          | ~ T^-1/2                      |
| Soave (SRK)    | 0          | 1          | alpha(T, omega)               |
| Peng-Robinson  | 1 - sqrt 2 | 1 + sqrt 2 | alpha(T, omega)               |

Units are SI: T in K, P in Pa, V in m^3/mol.
"""

from __future__ import annotations

import numpy as np

from .species import R, Species

__all__ = ["CubicEOS", "VanDerWaals", "RedlichKwong", "SoaveRedlichKwong", "PengRobinson"]


class CubicEOS:
    eps = 0.0
    sig = 0.0
    omega_a = 27 / 64
    omega_b = 1 / 8
    name = "cubic"

    def __init__(self, species: Species):
        self.sp = species
        self.ac = self.omega_a * R**2 * species.Tc**2 / species.Pc
        self.b = self.omega_b * R * species.Tc / species.Pc

    # temperature dependence -----------------------------------------
    def alpha(self, T):
        return 1.0

    def dalpha_dT(self, T):
        return 0.0

    def a(self, T):
        return self.ac * self.alpha(T)

    def da_dT(self, T):
        return self.ac * self.dalpha_dT(T)

    # P-V-T -------------------------------------------------------------
    def pressure(self, T, V):
        V = np.asarray(V, dtype=float)
        return R * T / (V - self.b) - self.a(T) / ((V + self.eps * self.b) * (V + self.sig * self.b))

    def _AB(self, T, P):
        return self.a(T) * P / (R * T) ** 2, self.b * P / (R * T)

    def z_roots(self, T, P):
        """All real compressibility roots with ``Z > B`` (sorted ascending)."""
        A, B = self._AB(T, P)
        e, s = self.eps, self.sig
        coeffs = [
            1.0,
            -(1 + B - (e + s) * B),
            A + e * s * B**2 - (e + s) * B * (B + 1),
            -(A * B + e * s * B**2 * (B + 1)),
        ]
        roots = np.roots(coeffs)
        # near the critical point the three roots coalesce and np.roots returns tiny
        # spurious imaginary parts (the triple root is ill-conditioned), so allow some slack
        real = np.sort(roots[np.abs(roots.imag) < 1e-4 * np.maximum(1, np.abs(roots.real))].real)
        real = real[real > B]
        if real.size == 0:
            raise ArithmeticError("no physical root")
        return real

    def Z(self, T, P, phase: str = "stable"):
        """Compressibility factor. ``phase`` is "vapor", "liquid" or "stable" (lowest Gibbs energy)."""
        roots = self.z_roots(T, P)
        if phase == "vapor":
            return float(roots[-1])
        if phase == "liquid":
            return float(roots[0])
        if phase != "stable":
            raise ValueError("phase must be 'vapor', 'liquid' or 'stable'")
        if roots.size == 1:
            return float(roots[0])
        lnphi = [self._lnphi_from_Z(T, P, z) for z in (roots[0], roots[-1])]
        return float(roots[0] if lnphi[0] < lnphi[1] else roots[-1])

    def volume(self, T, P, phase: str = "stable"):
        return self.Z(T, P, phase) * R * T / P

    def density(self, T, P, phase: str = "stable"):
        """Mass density in kg/m^3."""
        return self.sp.M / self.volume(T, P, phase)

    # fugacity and departure functions ------------------------------------
    def _log_term(self, Z, B):
        if self.sig == self.eps:  # van der Waals
            return B / Z
        return np.log((Z + self.sig * B) / (Z + self.eps * B)) / (self.sig - self.eps)

    def _lnphi_from_Z(self, T, P, Z):
        A, B = self._AB(T, P)
        return Z - 1 - np.log(Z - B) - (A / B) * self._log_term(Z, B)

    def ln_fugacity_coefficient(self, T, P, phase: str = "stable"):
        return float(self._lnphi_from_Z(T, P, self.Z(T, P, phase)))

    def fugacity(self, T, P, phase: str = "stable"):
        return P * np.exp(self.ln_fugacity_coefficient(T, P, phase))

    def enthalpy_departure(self, T, P, phase: str = "stable"):
        """Residual molar enthalpy ``H - H_ideal`` (J/mol)."""
        Z = self.Z(T, P, phase)
        A, B = self._AB(T, P)
        I = self._log_term(Z, B)
        return float(R * T * (Z - 1) + (T * self.da_dT(T) - self.a(T)) / self.b * I)

    def entropy_departure(self, T, P, phase: str = "stable"):
        """Residual molar entropy ``S - S_ideal`` at the same T and P (J/(mol K))."""
        Z = self.Z(T, P, phase)
        A, B = self._AB(T, P)
        I = self._log_term(Z, B)
        return float(R * np.log(Z - B) + self.da_dT(T) / self.b * I)

    def second_virial(self, T):
        """Second virial coefficient implied by the EOS, ``b - a / (R T)`` (m^3/mol)."""
        return self.b - self.a(T) / (R * T)

    # phase equilibrium ---------------------------------------------------
    def saturation_pressure(self, T, P0=None, tol=1e-12, max_iter=500):
        """Vapour pressure from equal liquid and vapour fugacities (successive substitution)."""
        Tr = T / self.sp.Tc
        if Tr >= 1:
            raise ValueError("no saturation above the critical temperature")
        P = P0 or self.sp.Pc * np.exp(5.373 * (1 + self.sp.omega) * (1 - 1 / Tr))  # Wilson estimate
        for _ in range(max_iter):
            roots = self._roots_or_none(T, P)
            if roots is None or roots.size < 2:
                # only one root: step towards the two-phase region
                V = self.z_roots(T, P)[0] * R * T / P
                P *= 1.1 if V > 2.5 * self.b else 0.9  # vapour-like root -> raise P
                continue
            ln_l = self._lnphi_from_Z(T, P, roots[0])
            ln_v = self._lnphi_from_Z(T, P, roots[-1])
            P_new = P * np.exp(ln_l - ln_v)
            if abs(P_new / P - 1) < tol:
                return float(P_new)
            P = P_new
        raise RuntimeError("saturation pressure did not converge")

    def _roots_or_none(self, T, P):
        try:
            return self.z_roots(T, P)
        except ArithmeticError:
            return None

    def __repr__(self):
        return f"{self.name}({self.sp.name})"


class VanDerWaals(CubicEOS):
    name = "vdW"


class RedlichKwong(CubicEOS):
    sig = 1.0
    omega_a = 1 / (9 * (2 ** (1 / 3) - 1))
    omega_b = (2 ** (1 / 3) - 1) / 3
    name = "RK"

    def alpha(self, T):
        return (T / self.sp.Tc) ** -0.5

    def dalpha_dT(self, T):
        return -0.5 * (T / self.sp.Tc) ** -1.5 / self.sp.Tc


class _SoaveAlpha(CubicEOS):
    def m(self):
        raise NotImplementedError

    def alpha(self, T):
        return (1 + self.m() * (1 - np.sqrt(T / self.sp.Tc))) ** 2

    def dalpha_dT(self, T):
        m = self.m()
        return -m * (1 + m * (1 - np.sqrt(T / self.sp.Tc))) / np.sqrt(T * self.sp.Tc)


class SoaveRedlichKwong(_SoaveAlpha):
    sig = 1.0
    omega_a = 1 / (9 * (2 ** (1 / 3) - 1))
    omega_b = (2 ** (1 / 3) - 1) / 3
    name = "SRK"

    def m(self):
        w = self.sp.omega
        return 0.480 + 1.574 * w - 0.176 * w**2


class PengRobinson(_SoaveAlpha):
    eps = 1 - np.sqrt(2)
    sig = 1 + np.sqrt(2)
    omega_a = 0.45723552892138218
    omega_b = 0.07779607390388845
    name = "PR"

    def m(self):
        w = self.sp.omega
        return 0.37464 + 1.54226 * w - 0.26992 * w**2
