# Linear K_local at current chord length for the simple 12-DOF corotational path.

from __future__ import annotations

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.utilities.B_matrix import (
    StrainDisplacementOperator as TimoshenkoStrain,
)
from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.utilities.D_matrix import (
    MaterialStiffnessOperator as TimoshenkoMaterial,
)
from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.utilities.k0_timoshenko import (
    assemble_timoshenko_K0,
    timoshenko_quadrature_orders_from_element_array,
)
from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.utilities.B_matrix import (
    StrainDisplacementOperator as EulerBernoulliStrain,
)
from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.utilities.D_matrix import (
    MaterialStiffnessOperator as EulerBernoulliMaterial,
)
from onedfem.pre_processing.element_library.shape_function_registry import get_shape_function_operator


def timoshenko_K_local(elem, L: float) -> np.ndarray:
    orders = timoshenko_quadrature_orders_from_element_array(elem.element_array)
    sf = get_shape_function_operator("LinearTimoshenkoBeamElement3D", L)
    strain = TimoshenkoStrain(element_length=L)
    mat = TimoshenkoMaterial(
        youngs_modulus=elem.E,
        shear_modulus=elem.G,
        cross_section_area=elem.A,
        moment_inertia_y=elem.I_y,
        moment_inertia_z=elem.I_z,
        torsion_constant=elem.J_t,
        shear_correction_factor=elem.kappa,
        y_sc=elem.y_sc,
        z_sc=elem.z_sc,
    )
    return assemble_timoshenko_K0(sf, strain, mat.assembly_form(), L / 2.0, orders)


def euler_bernoulli_K_local(elem, L: float) -> np.ndarray:
    axial = max(int(elem.element_array[3]), 1)
    bend_y = max(int(elem.element_array[4]), 2)
    bend_z = max(int(elem.element_array[5]), 2)
    torsion = max(int(elem.element_array[8]), 1)
    n_gauss = max(axial, bend_y, bend_z, torsion, 2)
    sf = get_shape_function_operator("LinearEulerBernoulliBeamElement3D", L)
    strain = EulerBernoulliStrain(element_length=L)
    mat = EulerBernoulliMaterial(
        youngs_modulus=elem.E,
        shear_modulus=elem.G,
        cross_section_area=elem.A,
        moment_inertia_y=elem.I_y,
        moment_inertia_z=elem.I_z,
        torsion_constant=elem.J_t,
    )
    D = mat.assembly_form()
    xi, w = np.polynomial.legendre.leggauss(n_gauss)
    Ke = np.zeros((12, 12), dtype=np.float64)
    detJ = L / 2.0
    for xi_g, w_g in zip(xi, w):
        _, dN_dξ, d2N_dξ2 = sf.natural_coordinate_form(np.array([xi_g]))
        B = strain.physical_coordinate_form(dN_dξ, d2N_dξ2)[0]
        Ke += B.T @ D @ B * w_g * detJ
    return Ke


def K_local(elem, L: float) -> np.ndarray:
    if elem.kernel == "timoshenko":
        return timoshenko_K_local(elem, L)
    return euler_bernoulli_K_local(elem, L)
