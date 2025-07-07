#ifndef STEP_PROPAGATORS_H
#define STEP_PROPAGATORS_H

#include <hls_x_complex.h>
#include <hls_math.h>

#define DIM 256

typedef hls::x_complex<float> complex_t;

// Parámetros globales (definidos en step_propagators.cpp)
extern const float eps;          // Epsilon para estabilidad
extern const float k;            // Número de onda
extern const float dz;           // Paso en z
extern const float dy;           // Paso en y
extern const float dx;           // Paso en x
extern const float n0;           // Índice de refracción lineal
extern const float n2;           // Índice de refracción no lineal
extern const float alpha;        // Coeficiente de absorción lineal
extern const float beta;         // Coeficiente de absorción de dos fotones
extern complex_t ung;      // Factor para método ADI
extern const float phase_const;  // Constante para cálculo de fase (Kerr)
extern const float attenuation;  // Factor de atenuación lineal
extern const float tpa_const;    // Constante para absorción de dos fotones

/**
 * Resuelve un sistema tridiagonal con estructura especial (Thomas).
 */
void custom_thomas_solver(
    complex_t dp,        // Valor común de la diagonal principal
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t b[DIM],    // Vector del lado derecho
    complex_t x[DIM]     // Vector solución
);

/**
 * Multiplica una matriz tridiagonal con estructura especial por un vector.
 */
void compute_b_vector(
    complex_t dp,        // Valor común de la diagonal principal
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t x0[DIM],   // Vector de entrada
    complex_t b[DIM]     // Vector resultado
);

/**
 * Aplica el método ADI en la dirección X.
 */
void adi_x(
    complex_t phi[DIM][DIM],       // Campo de entrada
    complex_t phi_inter[DIM][DIM]  // Campo intermedio
);

/**
 * Aplica el método ADI en la dirección Y.
 */
void adi_y(
    complex_t phi[DIM][DIM],       // Campo de entrada
    complex_t phi_inter[DIM][DIM]  // Campo intermedio
);

/**
 * Aplica el efecto no lineal Kerr (half-step).
 */
void half_nonlinear(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]   // Campo de salida
);

/**
 * Aplica absorción lineal (half-step).
 */
void half_linear_absorption(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]   // Campo de salida
);

/**
 * Aplica absorción de dos fotones (half-step).
 */
void half_2photon_absorption(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]   // Campo de salida
);

/**
 * Paso completo de propagación combinando todos los operadores:
 * ADI X → half-TPA → half-Kerr → half-linear → ADI Y → half-TPA → half-Kerr → half-linear
 */
void propagation_step(
    complex_t phi_in[DIM][DIM],      // Campo de entrada
    complex_t phi_out[DIM][DIM]      // Campo de salida
);

#endif // STEP_PROPAGATORS_H
