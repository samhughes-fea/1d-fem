# Scope

`1d-fem` (the `onedfem` package) answers *what does the beam do?*: beam finite
elements in 3D for slender composite structures, driven by text job inputs.

## In scope

- Beam, bar and truss elements, from Euler-Bernoulli to geometrically exact,
  with taper, pretwist, warping and GBT distortion.
- Assembly, supports and loads; linear and nonlinear static, eigen, buckling,
  harmonic and transient analyses; rotating structures and follower loads.
- Ply checks along the beam, load-case envelopes and layup sensitivities of beam
  responses.
- Beam models of a blade from a geometry-kernel loft, through the section provider.
- Validation against closed forms and external reference solvers.

## Out of scope

| Topic | Where it belongs |
| --- | --- |
| Blade and section geometry | [implicit-geometry-kernel](https://github.com/samhughes-fea/implicit-geometry-kernel) |
| Ply and laminate constitutive behaviour | [laminate-constitutive-model](https://github.com/samhughes-fea/laminate-constitutive-model) |
| One cross-section's stiffness, mass and stress recovery | [anisotropic-section-model](https://github.com/samhughes-fea/anisotropic-section-model) |
| Shell and solid models of the blade | [shell-fem](https://github.com/samhughes-fea/shell-fem) |
