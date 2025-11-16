# Investigation Diagnostic Report: Timoshenko and Levinson Element Formulation

## Executive Summary

Investigation of why Timoshenko and Levinson beam elements produce displacements that are 300-500x too small compared to analytical solutions.

**Status**: Integration order configuration implemented, but root cause not yet identified.

## Completed Investigations

### Phase 0: Integration Order Configuration ✅

1. **Updated element.txt files**:
   - Timoshenko (job_0001): Changed shear_y_order and shear_z_order from 0 to 2
   - Levinson (job_0002): Changed shear_y_order and shear_z_order from 0 to 3

2. **Modified element classes**:
   - `timoshenko_3D.py`: Now uses integration orders from element_array, defaults to 2 for shear if 0
   - `levinson_3D.py`: Now uses integration orders from element_array, defaults to 3 for shear if 0

3. **Verification**: Integration orders are now being read and used correctly.

### Phase 1: Component Verification ✅

1. **D-Matrix Verification**:
   - ✅ Timoshenko: D[3,3] = kappa*G*A = 8.842500e+07 N (correct)
   - ✅ Levinson: D[3,3] = G*A (verified in previous work)
   - ✅ Shear stiffness terms are non-zero and correct

2. **B-Matrix Verification**:
   - ✅ B-matrix has non-zero shear terms (B[3, :] norm ~ 4-10)
   - ✅ B-matrix correctly computes γ_xy = du_y/dx - θ_z
   - ⚠️ B[3, 11] sign appears positive at some Gauss points, but this is correct given shape function values

3. **Stiffness Matrix Integration**:
   - ✅ Element stiffness matrix is computed correctly
   - ✅ Integration K^e = ∫ B^T D B |J| dξ is performed correctly
   - ✅ Computed values match logged values:
     - Ke[1,1] = 5.305500e+08
     - Ke[1,5] = 1.768500e+06
     - Ke[5,5] = 1.817556e+04

## Current Results

After implementing integration order fixes:
- **Euler-Bernoulli**: ✅ Perfect match (-3.041259 mm)
- **Timoshenko**: ❌ Still 300x too small (-0.009885 mm vs -3.052568 mm)
- **Levinson**: ❌ Still 500x too small (-0.005844 mm vs -3.050683 mm)

## Remaining Issues

### Potential Root Causes

1. **Selective Integration Not Implemented**:
   - Current implementation uses a single quadrature_order for all terms
   - Timoshenko theory may require different integration orders for bending vs shear terms
   - Reduced integration (order 2) for shear is common to avoid shear locking
   - Full integration (order 3) for bending may be needed

2. **Stiffness Matrix Scaling**:
   - Computed stiffness values appear correct
   - But if stiffness is 300x too high, displacements would be 300x too small
   - Need to verify if stiffness matrix is being scaled incorrectly somewhere

3. **Load Application**:
   - Point loads are applied through element force vectors
   - Need to verify that loads are correctly applied to global force vector
   - Check if loads are being scaled or filtered incorrectly

4. **Global Assembly**:
   - Element stiffness matrices are assembled correctly
   - But need to verify DOF mapping is correct
   - Check if there's any scaling or transformation applied during assembly

## Recommendations

1. **Implement Selective Integration**:
   - Use different quadrature orders for different terms
   - Bending: order 3 (full integration)
   - Shear: order 2 (reduced integration to avoid shear locking)
   - This would require significant code changes to element stiffness computation

2. **Verify Load Application**:
   - Check element force vectors for point loads
   - Verify global force vector contains correct load values
   - Compare with Euler-Bernoulli (which works correctly)

3. **Compare Global Stiffness Matrices**:
   - Extract and compare global K for all three element types
   - Check if Timoshenko/Levinson global K is significantly different from Euler-Bernoulli
   - Verify no unexpected scaling or transformation

4. **Single Element Test**:
   - Create a single-element test case
   - Compare computed vs analytical for isolated element
   - This will isolate whether issue is in element formulation or global assembly

## Files Modified

1. `jobs/job_0001/element.txt` - Updated shear integration orders
2. `jobs/job_0002/element.txt` - Updated shear integration orders
3. `pre_processing/element_library/timoshenko/timoshenko_3D.py` - Added integration order logic
4. `pre_processing/element_library/levinson/levinson_3D.py` - Added integration order logic

## Verification Scripts Created

1. `verify_timoshenko_b_matrix.py` - B-matrix verification
2. `verify_timoshenko_stiffness_integration.py` - Stiffness integration verification
3. `compare_element_stiffness.py` - Stiffness comparison (needs fixing)

## Next Steps

1. Implement selective integration for Timoshenko and Levinson elements
2. Verify load application and global force vector
3. Compare global stiffness matrices between element types
4. Create single-element test cases for isolated verification

