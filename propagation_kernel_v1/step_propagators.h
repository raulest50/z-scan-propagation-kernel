// step_propagators.h
#ifndef STEP_PROPAGATORS_H
#define STEP_PROPAGATORS_H

#include <hls_x_complex.h>
#include "hls_math.h"

#define Nz 120     // Maximum number of z-steps (NDZ)
#define DIM 256
#define MATRIX_SIZE DIM*DIM
#define MAX_NY DIM
#define MAX_NX DIM

typedef hls::x_complex<float> complex_t;

// Variables globales
extern float eps;          // Epsilon para estabilidad
extern float k;            // Número de onda
extern float k0;           // Número de onda
extern float dz;           // Paso en z
extern float dy;           // Paso en y
extern float dx;           // Paso en x
extern float n0;           // Índice de refracción
extern float n2;           // Índice de refracción no lineal
extern float alpha;        // Coeficiente de absorción lineal
extern float beta;         // Coeficiente de absorción de dos fotones
extern complex_t ung;      // Factor para método ADI
extern float phase_const;  // Constante para cálculo de fase (Kerr)
extern float attenuation;  // Factor de atenuación
extern float tpa_const;    // Constante para absorción de dos fotones

/**
 * Solves a tridiagonal system with special structure using Thomas algorithm.
 */
void custom_thomas_solver(
    complex_t dp,        // Diagonal principal (valor común)
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t b[DIM],    // Vector del lado derecho
    complex_t x[DIM]     // Vector solución
);

/**
 * Multiplies a tridiagonal matrix with special structure by a vector.
 */
void compute_b_vector(
    complex_t dp,        // Diagonal principal (valor común)
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t x0[DIM],   // Vector de entrada
    complex_t b[DIM]     // Vector resultado
);

/**
 * Implements the ADI method in X direction.
 */
void adi_x(
    complex_t phi[DIM][DIM],       // Campo de entrada
    complex_t phi_inter[DIM][DIM]  // Campo intermedio
);

/**
 * Implements the ADI method in Y direction.
 */
void adi_y(
    complex_t phi[DIM][DIM],       // Campo de entrada
    complex_t phi_inter[DIM][DIM]  // Campo intermedio
);

/**
 * Applies Kerr nonlinear effect.
 */
void half_nonlinear(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]   // Campo de salida
);

/**
 * Applies linear absorption.
 */
void half_linear_absorption(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]   // Campo de salida
);

/**
 * Applies two-photon absorption.
 */
void half_2photon_absorption(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]   // Campo de salida
);

/**
 * Complete propagation step combining all operators.
 */
void propagation_step(
    complex_t phi_in[DIM][DIM],      // Campo de entrada
    complex_t phi_out[DIM][DIM]      // Campo de salida
);

#endif // STEP_PROPAGATORS_H
