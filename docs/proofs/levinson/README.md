# Levinson beam — proofs

Levinson is a shear-deformable formulation (like Timoshenko). Proofs specific to Levinson can be added here; for shear force from the constitutive law, the same qualitative conclusion as Timoshenko applies (non-zero \(V_y\), \(V_z\) from \(\sigma = D\,\varepsilon\) when shear strain is present).

Levinson uses **GA** in the D-matrix (no shear correction factor κ). The need for κ is removed by the higher-order approximation of shear deformation in the kinematics: the B-matrix includes terms with α (e.g. α(∂²θ/∂x²) in the shear strain), which better matches the actual shear stress distribution over the cross-section. Unlike Timoshenko (constant shear strain → κGA), the constitutive shear stiffness is therefore GA (κ = 1 implicitly).
