# Beam shear correction (\(\kappa\)) and thin-walled sections

This note complements [`JOB_INPUT_BEAM_WARPING.md`](JOB_INPUT_BEAM_WARPING.md) (warping / \(\Gamma\)) and applies to **shear-deformable** beams (`LinearTimoshenkoBeamElement3D`, Total Lagrangian Timoshenko, etc.).

## Default \(\kappa = 5/6\)

When `section.txt` uses the **6-column** tier (area and inertias only), Timoshenko elements use the classical **rectangular** shear correction \(\kappa = 5/6\) for the \(\kappa G A\) shear stiffness diagonal in **`D`**.

That value is **not** universal:

- **Thin-walled open sections** (I-beams, channels, angles) require \(\kappa\) that depends on the profile geometry and loading; \(5/6\) can be **very inaccurate** for shear stiffness.
- **Solid rectangles** are closer to \(5/6\); **solid circles** use a different convention in some texts.

## Authoring guidance

1. Prefer **tier 8+** `section.txt` rows with an explicit **`[kappa]`** (and **`[alpha]`** when used by your pipeline) when using Timoshenko theory on **non-rectangular** or **thin-walled** sections.
2. For **open thin-walled** analysis, treat beam theory shear stiffness as **engineering approximation**; plate/shell or specialized thin-walled beam theories may be needed for production accuracy.
3. Optional runtime hint: set environment variable **`FEM_BEAM_SECTION_HINTS=1`** so elements log a **warning** when the section slice has **no explicit \(\kappa\)** (tier-6-style row) and the implementation falls back to \(5/6\).

## Warping interaction

Vlasov warping stiffness (\(\Gamma\), \(\chi\) DOFs) addresses **constraint-torsion / warping** interaction; it does **not** replace a correct \(\kappa G A\) for transverse shear in Timoshenko theory.
