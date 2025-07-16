# Validation Data for Testbench

This directory contains validation data files (.dat) used for testing the implementation of various operators and propagation steps in the deep tissue imaging simulation. Below is a description of the input-output relationships for each set of validation data.

## Individual Operators Validation Data

### ADI X Operator
- **Input**: `individual_ops_in.dat` - Initial complex field matrix
- **Output**: `adi_x_out.dat` - Result after applying the ADI X operator

### ADI Y Operator
- **Input**: `individual_ops_in.dat` - Initial complex field matrix
- **Output**: `adi_y_out.dat` - Result after applying the ADI Y operator

### Combined Nonlinear Operators
- **Input**: `individual_ops_in.dat` - Initial complex field matrix
- **Output**: `half_nonlinear_ops_combined_out.dat` - Result after applying the following operators in sequence:
  1. half_2photon_absorption
  2. half_nonlinear
  3. half_linear_absorption

## B-Vector and Thomas Solver Validation Data

### B-Vector Computation
- **Inputs**:
  - `b_vector_dp_in.dat` - dp parameter
  - `b_vector_dp1_in.dat` - dp1 parameter
  - `b_vector_dp2_in.dat` - dp2 parameter
  - `b_vector_do_in.dat` - do parameter
  - `b_vector_x0_in.dat` - x0 vector
- **Output**: `b_vector_out.dat` - Computed b-vector

### Thomas Solver
- **Inputs**:
  - `b_vector_dp_in.dat` - dp parameter
  - `b_vector_dp1_in.dat` - dp1 parameter
  - `b_vector_dp2_in.dat` - dp2 parameter
  - `b_vector_do_in.dat` - do parameter
  - `b_vector_out.dat` - b-vector
- **Output**: `thomas_solver_out.dat` - Solution from the Thomas solver

## Full Step Propagation Validation Data

### Single Step Propagation
- **Input**: `initial_field.dat` - Initial complex field
- **Output**: `full_step_within_tissue_out.dat` - Field after a single full step propagation

### Multi-Step Propagation
- **Input**: `initial_field.dat` - Initial complex field
- **Intermediate Outputs**:
  - `full_step_within_tissue_step_1.dat` - Field after 1 step
  - `full_step_within_tissue_step_40.dat` - Field after 40 steps
  - `full_step_within_tissue_step_80.dat` - Field after 80 steps
  - `full_step_within_tissue_step_120.dat` - Field after 120 steps
- **Final Output**: `full_step_within_tissue_120_steps.dat` - Field after 120 steps (same as the last intermediate output)

## File Format

All .dat files contain complex matrices stored as text with real and imaginary parts. Each line in the file represents a single complex number with the format:
```
<real_part> <imaginary_part>
```

For 2D matrices, the data is stored in row-major order.