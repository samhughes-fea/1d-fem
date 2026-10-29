# I-beam large twist with restrained warping

A steel I-beam cantilever (300 × 150 mm, flanges 6 mm, web 4 mm, 2 m long). The root is
clamped with its warping restrained, and a tip torque is applied in eight increments through
the nonlinear static job runner with `CorotationalBenscoterBeamElement3D`.

```
python examples/thin_walled/large_twist/run.py
```

| torque (kN m) | corotational | nonlinear Vlasov | linear Benscoter |
|---|---|---|---|
| 1 | 0.1409 rad | 0.1405 rad | 0.1435 rad |
| 2 | 0.2688 | 0.2682 | 0.2871 |
| 4 | 0.4752 | 0.4742 | 0.5742 |
| 6 | 0.6327 | 0.6313 | 0.8612 |

Twisting stretches the outer fibres by `r²φ'²/2`. The axial stress of that stretch stiffens
the beam (the Wagner effect), so at 6 kN m the twist is 27 % below the linear value. The
reference is the nonlinear Vlasov equation `GJφ' − EΓφ''' + ½EI_nφ'³ = T`, solved by
collocation. The corotational elements stay within 0.25 % of it; the difference is the
warping shear flexibility that Benscoter has and Vlasov leaves out.

Results are in `summary.json` and `twist.png`. Method:
[docs/element_library/corotational_section_beams.md](../../../docs/element_library/corotational_section_beams.md).
