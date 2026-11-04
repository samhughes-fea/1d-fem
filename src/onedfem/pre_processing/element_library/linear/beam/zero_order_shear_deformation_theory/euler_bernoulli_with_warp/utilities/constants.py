# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/constants.py
"""Outer sizes for Vlasov warping EB: ``N_DOF`` = len(U_e), ``N_STRAIN`` = rows of B and D.

``N_DOF = 14`` (7 per node), ``N_STRAIN = 7``; ``B`` is (n_gp, N_STRAIN, N_DOF), ``D`` is (N_STRAIN, N_STRAIN).
"""

N_STANDARD_DOF = 12
N_WARPING_DOF = 2
N_DOF = N_STANDARD_DOF + N_WARPING_DOF
N_STRAIN = 7
