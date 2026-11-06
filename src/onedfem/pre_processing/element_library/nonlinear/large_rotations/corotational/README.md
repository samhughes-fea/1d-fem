# Corotational beams

**Canonical family:** `corotational_section.py` — Battini–Pacoste local map `p(U)`
with cubic/quartic local energy, wrapping any linear core
(`CorotationalEulerBernoulliBeamElement3D`, Timoshenko, Benscoter, GBT, Reddy, Levinson).

**Simple 12-DOF path:** `corotational_3D.py` (`CorotationalBeamElement3D`) — `T.T @ K_local(L) @ T`
with optional finite-difference tangent. Diagnostic / kernel switch (`timoshenko` | `euler_bernoulli`);
not merged with the section family.

Frame, rotation, Battini `p(U)`, local energy and simple-path `K_local` live in `utilities/`.
The two implementations are **not merged**.

**Geometry vs constitutive curvature**

- Discrete curved axes: **kinked nodes** (straight CR chords).
- `precurvature.txt` `k_x0`: pretwist rate (`SectionAnalysis.with_pretwist`); not constitutive `φ_x` in `E_0`.
- `k_y0`, `k_z0`: constitutive precurvature `E_0 = [0, κ_y0, κ_z0, 0, 0, 0]` on the chord (`ε − E_0` in the local core). CR does not build the linear 8-subelement arc.
