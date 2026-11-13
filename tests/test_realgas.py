import numpy as np
import pytest
from scipy.integrate import quad

import realgas as rg

ALL = list(rg.EOS.values())
CO2 = rg.get("CO2")
ZC = {"vdW": 3 / 8, "RK": 1 / 3, "SRK": 1 / 3, "PR": 0.3074013086987}


@pytest.mark.parametrize("E", ALL)
def test_critical_point_is_an_inflection(E):
    eos = E(CO2)
    Vc = ZC[E.name] * rg.R * CO2.Tc / CO2.Pc
    h = 1e-4 * Vc
    P = lambda V: eos.pressure(CO2.Tc, V)
    assert P(Vc) == pytest.approx(CO2.Pc, rel=1e-9)
    dP = (P(Vc + h) - P(Vc - h)) / (2 * h)
    d2P = (P(Vc + h) - 2 * P(Vc) + P(Vc - h)) / h**2
    # compare with the natural scale Pc / Vc (and Pc / Vc^2)
    assert abs(dP) * Vc / CO2.Pc < 1e-6
    assert abs(d2P) * Vc**2 / CO2.Pc < 1e-4


@pytest.mark.parametrize("E", ALL)
def test_critical_compressibility(E):
    eos = E(CO2)
    assert eos.Z(CO2.Tc, CO2.Pc) == pytest.approx(ZC[E.name], abs=2e-3)


@pytest.mark.parametrize("E", ALL)
def test_ideal_gas_limit(E):
    eos = E(rg.get("N2"))
    assert eos.Z(300.0, 1.0) == pytest.approx(1.0, abs=1e-6)
    assert eos.ln_fugacity_coefficient(300.0, 1.0) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("E", ALL)
def test_low_pressure_slope_is_second_virial(E):
    eos = E(rg.get("CH4"))
    T, P = 250.0, 100.0
    slope = (eos.Z(T, P) - 1) * rg.R * T / P
    assert slope == pytest.approx(eos.second_virial(T), rel=1e-3)


@pytest.mark.parametrize("E", ALL)
def test_volume_round_trip(E):
    eos = E(CO2)
    for T, P in [(320.0, 2e6), (350.0, 1.5e7), (280.0, 5e5)]:
        assert eos.pressure(T, eos.volume(T, P)) == pytest.approx(P, rel=1e-9)


@pytest.mark.parametrize("E", [rg.SoaveRedlichKwong, rg.PengRobinson])
@pytest.mark.parametrize("name", ["CH4", "C3H8", "CO2", "N2", "Ar"])
def test_acentric_factor_definition(E, name):
    # omega = -1 - log10(Psat / Pc) at Tr = 0.7; SRK and PR are fitted to honour it
    sp = rg.get(name)
    psat = E(sp).saturation_pressure(0.7 * sp.Tc)
    assert psat == pytest.approx(sp.Pc * 10 ** (-1 - sp.omega), rel=1e-2)


def test_co2_vapour_pressure_at_0C():
    # experimental value 3.485 MPa (NIST)
    assert rg.PengRobinson(CO2).saturation_pressure(273.15) == pytest.approx(3.485e6, rel=1e-2)


@pytest.mark.parametrize("E", ALL)
def test_saturation_equal_fugacities(E):
    eos = E(rg.get("C3H8"))
    T = 300.0
    P = eos.saturation_pressure(T)
    lnl = eos.ln_fugacity_coefficient(T, P, "liquid")
    lnv = eos.ln_fugacity_coefficient(T, P, "vapor")
    assert lnl == pytest.approx(lnv, abs=1e-10)
    assert eos.density(T, P, "liquid") > 5 * eos.density(T, P, "vapor")


def test_vdw_maxwell_equal_area():
    eos = rg.VanDerWaals(CO2)
    T = 0.85 * CO2.Tc
    P = eos.saturation_pressure(T)
    VL, VV = eos.volume(T, P, "liquid"), eos.volume(T, P, "vapor")
    area, _ = quad(lambda V: eos.pressure(T, V), VL, VV, limit=200)
    assert area == pytest.approx(P * (VV - VL), rel=1e-8)


def test_stable_phase_switches_at_saturation():
    eos = rg.PengRobinson(rg.get("C3H8"))
    T = 300.0
    P = eos.saturation_pressure(T)
    assert eos.Z(T, 0.98 * P) == pytest.approx(eos.Z(T, 0.98 * P, "vapor"))
    assert eos.Z(T, 1.02 * P) == pytest.approx(eos.Z(T, 1.02 * P, "liquid"))


@pytest.mark.parametrize("E", ALL)
@pytest.mark.parametrize("phase, T, P", [("vapor", 330.0, 3e6), ("liquid", 250.0, 5e6)])
def test_enthalpy_departure_is_consistent_with_fugacity(E, phase, T, P):
    # (d ln phi / dT)_P = -H_res / (R T^2)
    eos = E(CO2)
    h = 1e-3
    dlnphi = (eos.ln_fugacity_coefficient(T + h, P, phase)
              - eos.ln_fugacity_coefficient(T - h, P, phase)) / (2 * h)
    H = eos.enthalpy_departure(T, P, phase)
    assert dlnphi == pytest.approx(-H / (rg.R * T**2), rel=1e-5)


@pytest.mark.parametrize("E", ALL)
def test_gibbs_identity(E):
    # G_res = H_res - T S_res = R T ln phi
    eos = E(CO2)
    T, P = 310.0, 4e6
    G = eos.enthalpy_departure(T, P) - T * eos.entropy_departure(T, P)
    assert G == pytest.approx(rg.R * T * eos.ln_fugacity_coefficient(T, P), rel=1e-9)


def test_errors():
    with pytest.raises(ValueError):
        rg.PengRobinson(CO2).saturation_pressure(1.1 * CO2.Tc)
    with pytest.raises(ValueError):
        rg.PengRobinson(CO2).Z(300.0, 1e5, phase="gas")
    with pytest.raises(KeyError):
        rg.get("Xe2")
