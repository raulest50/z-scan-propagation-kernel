/**
 * @file hls_step_operators.h
 * @brief Header file for HLS implementation of step operators for wave propagation
 * 
 * This file contains the declarations of functions for implementing the step operators
 * in HLS for the Kria KV 260 FPGA. These functions are based on the Python implementations
 * in step_operators.py from the deep_tissue_imaging package.
 */

#ifndef HLS_STEP_OPERATORS_H
#define HLS_STEP_OPERATORS_H

#include <complex>
#include "ap_fixed.h"
#include "hls_stream.h"
#include "hls_math.h"

// Define maximum dimensions for arrays
#define MAX_SIZE 1024
#define MAX_NX 512
#define MAX_NY 512

// Define fixed-point data types for optimized resource usage
typedef ap_fixed<32, 16> fixed_t;  // For real values
typedef std::complex<fixed_t> complex_t;  // For complex values

/**
 * Solves a tridiagonal system with special structure using Thomas algorithm.
 * 
 * Latency estimated: O(N) cycles
 * Resources estimated: ~500 LUTs, ~800 FFs, 4-6 DSPs
 * 
 * @param dp Value for all elements in the main diagonal except first and last
 * @param dp1 Value for the first element in the main diagonal [0,0]
 * @param dp2 Value for the last element in the main diagonal [-1,-1]
 * @param do_val Value for all elements in the off-diagonals
 * @param b Right-hand side vector
 * @param x Solution vector
 * @param n Size of the system
 */
void custom_thomas_solver(
    complex_t dp,        // Diagonal principal (valor común)
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t b[MAX_SIZE],      // Vector del lado derecho
    complex_t x[MAX_SIZE],      // Vector solución
    int n                // Tamaño del sistema
);

/**
 * Multiplies a tridiagonal matrix with special structure by a vector.
 * 
 * Latency estimated: O(N) cycles
 * Resources estimated: ~300 LUTs, ~500 FFs, 2-4 DSPs
 * 
 * @param dp Value for all elements in the main diagonal except first and last
 * @param dp1 Value for the first element in the main diagonal [0,0]
 * @param dp2 Value for the last element in the main diagonal [-1,-1]
 * @param do_val Value for all elements in the off-diagonals
 * @param x0 Input vector to be multiplied
 * @param b Result of the matrix-vector multiplication
 * @param n Size of the vector
 */
void compute_b_vector(
    complex_t dp,        // Diagonal principal (valor común)
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t x0[MAX_SIZE],     // Vector de entrada
    complex_t b[MAX_SIZE],      // Vector resultado
    int n                // Tamaño del vector
);

/**
 * Implements the ADI method in X direction.
 * 
 * Latency estimated: O(Ny*Nx) cycles
 * Resources estimated: ~1200 LUTs, ~2000 FFs, 8-12 DSPs
 * 
 * @param phi Input field
 * @param phi_inter Intermediate field
 * @param Ny Y dimension
 * @param Nx X dimension
 * @param eps Epsilon for stability
 * @param k Wave number
 * @param dz Step in z
 * @param dx Step in x
 */
void adi_x(
    complex_t phi[MAX_NY][MAX_NX],         // Campo de entrada
    complex_t phi_inter[MAX_NY][MAX_NX],   // Campo intermedio
    int Ny,                                // Dimensión Y
    int Nx,                                // Dimensión X
    float eps,                             // Epsilon para estabilidad
    float k,                               // Número de onda
    float dz,                              // Paso en z
    float dx                               // Paso en x
);

/**
 * Implements the ADI method in Y direction.
 * 
 * Latency estimated: O(Nx*Ny) cycles
 * Resources estimated: ~1200 LUTs, ~2000 FFs, 8-12 DSPs
 * 
 * @param phi Input field
 * @param phi_inter Intermediate field
 * @param Nx X dimension
 * @param Ny Y dimension
 * @param eps Epsilon for stability
 * @param k Wave number
 * @param dz Step in z
 * @param dy Step in y
 */
void adi_y(
    complex_t phi[MAX_NY][MAX_NX],         // Campo de entrada
    complex_t phi_inter[MAX_NY][MAX_NX],   // Campo intermedio
    int Nx,                                // Dimensión X
    int Ny,                                // Dimensión Y
    float eps,                             // Epsilon para estabilidad
    float k,                               // Número de onda
    float dz,                              // Paso en z
    float dy                               // Paso en y
);

/**
 * Applies Kerr nonlinear effect.
 * 
 * Latency estimated: O(N) cycles
 * Resources estimated: ~400 LUTs, ~600 FFs, 4-6 DSPs
 * 
 * @param phi Input field
 * @param phi_out Output field
 * @param k_sample Wave number in the sample
 * @param n2_sample Nonlinear index
 * @param dz Step in z
 * @param size Size of the field
 */
void half_nonlinear(
    complex_t phi[MAX_SIZE],      // Campo de entrada
    complex_t phi_out[MAX_SIZE],  // Campo de salida
    float k_sample,               // Número de onda en la muestra
    float n2_sample,              // Índice no lineal
    float dz,                     // Paso en z
    int size                      // Tamaño del campo
);

/**
 * Applies linear absorption.
 * 
 * Latency estimated: O(N) cycles
 * Resources estimated: ~200 LUTs, ~300 FFs, 2-3 DSPs
 * 
 * @param phi Input field
 * @param phi_out Output field
 * @param alpha Absorption coefficient
 * @param dz Step in z
 * @param size Size of the field
 */
void half_linear_absorption(
    complex_t phi[MAX_SIZE],      // Campo de entrada
    complex_t phi_out[MAX_SIZE],  // Campo de salida
    float alpha,                  // Coeficiente de absorción
    float dz,                     // Paso en z
    int size                      // Tamaño del campo
);

/**
 * Applies two-photon absorption.
 * 
 * Latency estimated: O(N) cycles
 * Resources estimated: ~300 LUTs, ~500 FFs, 3-5 DSPs
 * 
 * @param phi Input field
 * @param phi_out Output field
 * @param beta Two-photon absorption coefficient
 * @param dz Step in z
 * @param size Size of the field
 */
void half_2photon_absorption(
    complex_t phi[MAX_SIZE],      // Campo de entrada
    complex_t phi_out[MAX_SIZE],  // Campo de salida
    float beta,                   // Coeficiente de absorción de dos fotones
    float dz,                     // Paso en z
    int size                      // Tamaño del campo
);

/**
 * Complete propagation step combining all operators.
 * 
 * Latency estimated: O(Nx*Ny) cycles
 * Resources estimated: ~3000 LUTs, ~5000 FFs, 20-30 DSPs
 * 
 * @param phi_in Input field
 * @param phi_out Output field
 * @param Nx X dimension
 * @param Ny Y dimension
 * @param eps Epsilon for stability
 * @param k Wave number
 * @param k_sample Wave number in the sample
 * @param n2_sample Nonlinear index
 * @param alpha Linear absorption coefficient
 * @param beta Two-photon absorption coefficient
 * @param dz Step in z
 * @param dx Step in x
 * @param dy Step in y
 */
void propagation_step(
    complex_t phi_in[MAX_NY][MAX_NX],      // Campo de entrada
    complex_t phi_out[MAX_NY][MAX_NX],     // Campo de salida
    int Nx,                                // Dimensión X
    int Ny,                                // Dimensión Y
    float eps,                             // Epsilon para estabilidad
    float k,                               // Número de onda
    float k_sample,                        // Número de onda en la muestra
    float n2_sample,                       // Índice no lineal
    float alpha,                           // Coeficiente de absorción lineal
    float beta,                            // Coeficiente de absorción de dos fotones
    float dz,                              // Paso en z
    float dx,                              // Paso en x
    float dy                               // Paso en y
);

#endif // HLS_STEP_OPERATORS_H