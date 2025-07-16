// step_propagators.h
#ifndef STEP_PROPAGATORS_H
#define STEP_PROPAGATORS_H

#if __has_include(<hls_x_complex.h>)
#  include <hls_x_complex.h>
#  include <hls_math.h>
#  include <hls_stream.h>
#else
#  include "hls_stub.h"
#endif
#include <ap_fixed.h>

#ifndef DIM
#define DIM 256
#endif

using data_t   = ap_fixed<32,16>;
using complex_t = hls::x_complex<data_t>;

// Parámetros globales (definidos en step_propagators.cpp)
extern const data_t eps;          // Epsilon para estabilidad
extern const data_t k;            // Número de onda
extern const data_t dz;           // Paso en z
extern const data_t dy;           // Paso en y
extern const data_t dx;           // Paso en x
extern const data_t n0;           // Índice de refracción lineal
extern const data_t n2;           // Índice de refracción no lineal
extern const data_t alpha;        // Coeficiente de absorción lineal
extern const data_t beta;         // Coeficiente de absorción de dos fotones
extern complex_t ung;      // Factor para método ADI
extern const data_t phase_const;  // Constante para cálculo de fase (Kerr)
extern const data_t attenuation;  // Factor de atenuación lineal
extern const data_t tpa_const;    // Constante para absorción de dos fotones

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
 * Streaming version applying TPA, Kerr and linear attenuation.
 */
void half_nonlin_ops(
    hls::stream<complex_t>& in,
    hls::stream<complex_t>& out
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
