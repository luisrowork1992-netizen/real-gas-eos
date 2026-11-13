# real-gas-eos

![tests](https://github.com/luisrowork1992-netizen/real-gas-eos/actions/workflows/tests.yml/badge.svg)

`realgas` implements the classic **cubic equations of state** for pure fluids: van der Waals, Redlich–Kwong, Soave–Redlich–Kwong and Peng–Robinson. It computes compressibility factors, densities, fugacities, vapour pressures and residual enthalpy and entropy.

## Features

- One generic cubic engine; each equation of state is a few lines (ε, σ, Ω_a, Ω_b, α(T))
- Critical constants Ω_a and Ω_b in closed form, so every model has an exact critical inflection point
- Picks the liquid, vapour or thermodynamically stable root (the one with the lowest Gibbs energy)
- Fugacity coefficient, residual enthalpy and residual entropy, with analytic dα/dT
- Saturation pressure from equal fugacities
  - starts from the Wilson estimate
  - uses successive substitution
- Second virial coefficient implied by each model
- Built-in critical data (T_c, P_c, ω, M) for N₂, O₂, Ar, CO₂, CH₄, C₂H₆, C₃H₈, H₂O and NH₃

## Install

```bash
git clone https://github.com/luisrowork1992-netizen/real-gas-eos.git
cd real-gas-eos
pip install -e .
```

## Quick start

```python
import realgas as rg

co2 = rg.PengRobinson(rg.get("CO2"))
print(co2.Z(310.0, 5e6))                    # compressibility at 310 K, 5 MPa
print(co2.density(250.0, 5e6))              # liquid CO2 density, kg/m^3
print(co2.saturation_pressure(273.15))      # 3.477e6 Pa (experiment: 3.485e6)
print(co2.enthalpy_departure(310.0, 5e6))   # H - H_ideal, J/mol
```

## Notebook

[`notebooks/real_gases.ipynb`](notebooks/real_gases.ipynb) covers five topics:

- CO₂ isotherms with the saturation lines
- the methane compressibility chart, comparing the four equations of state
- vapour-pressure curves compared with the acentric-factor definition
- the CO₂ coexistence curve
- residual enthalpy of nitrogen

## Tests

```bash
pip install pytest
pytest
```

The tests check the following:

- every model has an exact critical inflection point (P = P_c, ∂P/∂V = ∂²P/∂V² = 0), with Z_c = 3/8, 1/3 or 0.3074
- the ideal-gas limit
- the low-pressure slope equals the model's second virial coefficient
- the P(V) round trip
- the acentric-factor rule at T_r = 0.7 holds within 1 % for SRK and PR (5 fluids)
- the CO₂ vapour pressure at 0 °C is within 1 % of experiment
- liquid and vapour fugacities are equal at saturation
- the van der Waals equal-area (Maxwell) construction
- the stable root switches at P_sat
- residual enthalpy matches the temperature derivative of ln φ
- G_res = H_res − T·S_res

## Scope and limitations

- Pure fluids only; there are no mixing rules.
- Cubic equations of state are weak for strongly polar fluids. PR predicts the vapour pressure of water at 100 °C about 5 % low.
- Liquid densities from SRK and PR are typically 5–15 % off without a volume translation.

## References

- J. M. Smith, H. C. Van Ness and M. M. Abbott, *Introduction to Chemical Engineering Thermodynamics*, 7th ed., McGraw-Hill.
- D.-Y. Peng and D. B. Robinson, *Ind. Eng. Chem. Fundam.* 15, 59 (1976).
- G. Soave, *Chem. Eng. Sci.* 27, 1197 (1972).
- B. E. Poling, J. M. Prausnitz and J. P. O'Connell, *The Properties of Gases and Liquids*, 5th ed.

## License

MIT
