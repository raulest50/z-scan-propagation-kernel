
#include "step_propagators.h"
#include <hls_math.h>


float eps = 1e-12; // Epsilon para estabilidad
float k = 10681416.0;  // Número de onda
float k0 = 7853982.0;  // Número de onda
float dz = 1e-6;  // Paso en z
float dy = 1.7578e-7;  // Paso en y
float dx = 1.7578e-7;  // Paso en x

float n0 = 1.36;
float n2 = 3e-20;
float alpha = 0.3;
float beta = 1e-11;

// Calcular el factor ung
complex_t ung = complex_t ( 0, dz / (4 * k * dy * dy) );

// Constante para el cálculo de fase,  Kerr
float phase_const = k * n2 * dz / 2.0f;

// Calcular el factor de atenuación
float attenuation = hls::exp(-alpha * dz / 4.0f);

// constante para la exp de TPA
float tpa_const = -beta * dz / 4.0f;


/**
 * Solves a tridiagonal system with special structure using Thomas algorithm.
 */
void custom_thomas_solver(
    complex_t dp,        // Diagonal principal (valor común)
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t b[DIM],      // Vector del lado derecho
    complex_t x[DIM]      // Vector solución
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=dp,dp1,dp2,do_val bundle=CTRL
    #pragma HLS INTERFACE m_axi port=b offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=x offset=slave bundle=GMEM1

    // Buffers locales para mejorar el rendimiento
    complex_t c_prime[DIM-1];
    complex_t d_prime[DIM];

    #pragma HLS ARRAY_PARTITION variable=c_prime cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=d_prime cyclic factor=4

    // Eliminación hacia adelante
    // Primera fila
    c_prime[0] = do_val / dp1;
    d_prime[0] = b[0] / dp1;

    // Filas intermedias
    forward_elimination:
    for (int i = 1; i < DIM-1; i++) {
        #pragma HLS PIPELINE II=1
        complex_t denominator = dp - do_val * c_prime[i-1];
        c_prime[i] = do_val / denominator;
        d_prime[i] = (b[i] - do_val * d_prime[i-1]) / denominator;
    }

    // Última fila
    d_prime[DIM-1] = (b[DIM-1] - do_val * d_prime[DIM-2]) / (dp2 - do_val * c_prime[DIM-2]);

    // Sustitución hacia atrás
    x[DIM-1] = d_prime[DIM-1];

    back_substitution:
    for (int i = DIM-2; i >= 0; i--) {
        #pragma HLS PIPELINE II=1
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}


/**
 * Multiplies a tridiagonal matrix with special structure by a vector.
 */
void compute_b_vector(
    complex_t dp,        // Diagonal principal (valor común)
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t x0[DIM],     // Vector de entrada
    complex_t b[DIM]      // Vector resultado
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=dp,dp1,dp2,do_val bundle=CTRL
    #pragma HLS INTERFACE m_axi port=x0 offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=b offset=slave bundle=GMEM1

    // Primera fila
    b[0] = dp1 * x0[0] + do_val * x0[1];

    // Filas intermedias
    middle_rows:
    for (int i = 1; i < DIM-1; i++) {
        #pragma HLS PIPELINE II=1
        b[i] = do_val * x0[i-1] + dp * x0[i] + do_val * x0[i+1];
    }

    // Última fila
    b[DIM-1] = do_val * x0[DIM-2] + dp2 * x0[DIM-1];
}


/**
 * Implements the ADI method in X direction.
 */
void adi_x(
    complex_t phi[DIM][DIM],         // Campo de entrada
    complex_t phi_inter[DIM][DIM]   // Campo intermedio
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_inter offset=slave bundle=GMEM1

    // Buffer local para una columna
    complex_t phi_col[DIM];
    complex_t b_vec[DIM];
    complex_t result[DIM];

    #pragma HLS ARRAY_PARTITION variable=phi_col cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=b_vec cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=result cyclic factor=4


    process_rows:
    for (int j = 0; j < DIM; j++) {
        #pragma HLS DATAFLOW

        // Cargar columna en buffer local
        load_column:
        for (int i = 0; i < DIM; i++) {
            #pragma HLS PIPELINE II=1
            phi_col[i] = phi[j][i];
        }

        // Calcular ratios para condiciones de contorno
        complex_t ratio_x0, ratio_xn;
        if ( hls::hypot(phi_col[1].real(), phi_col[1].imag()) < eps) {
            ratio_x0.real(1.0);
            ratio_x0.imag(0);
        } else {
            ratio_x0 = phi_col[0] / phi_col[1];
        }

        if (hls::hypot(phi_col[DIM-2].real(), phi_col[DIM-2].imag()) < eps) {
            ratio_xn.real(1.0);
            ratio_xn.imag(0);
        } else {
            ratio_xn = phi_col[DIM-1] / phi_col[DIM-2];
        }

        // Calcular coeficientes para matriz B
        // (-ung - ung) == -2 * ung
        complex_t dp1_B = (-ung - ung) + complex_t(1.0, 0) + ung * ratio_x0;
        complex_t dp2_B = (-ung - ung) + complex_t(1.0, 0) + ung * ratio_xn;
        complex_t dp_B = (-ung - ung) + complex_t(1.0, 0);
        complex_t do_B = ung;

        // Calcular vector b
        compute_b_vector(dp_B, dp1_B, dp2_B, do_B, phi_col, b_vec);

        // Calcular coeficientes para matriz A
        // 2 * ung == ung + ung
        complex_t dp1_A = ung + ung + complex_t(1.0, 0) - ung * ratio_x0;
        complex_t dp2_A = ung + ung + complex_t(1.0, 0) - ung * ratio_xn;
        complex_t dp_A = ung + ung + complex_t(1.0, 0);
        complex_t do_A = -ung;

        // Resolver sistema tridiagonal
        custom_thomas_solver(dp_A, dp1_A, dp2_A, do_A, b_vec, result);

        // Guardar resultados
        store_results:
        for (int i = 0; i < DIM; i++) {
            #pragma HLS PIPELINE II=1
            phi_inter[j][i] = result[i];
        }
    }
}


/**
 * Implements the ADI method in Y direction.
 */
void adi_y(
    complex_t phi[DIM][DIM],         // Campo de entrada
    complex_t phi_inter[DIM][DIM]   // Campo intermedio
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_inter offset=slave bundle=GMEM1

    // Buffer local para una fila
    complex_t phi_row[DIM];
    complex_t b_vec[DIM];
    complex_t result[DIM];

    #pragma HLS ARRAY_PARTITION variable=phi_row cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=b_vec cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=result cyclic factor=4



    process_columns:
    for (int i = 0; i < DIM; i++) {
        #pragma HLS DATAFLOW

        // Cargar fila en buffer local
        load_row:
        for (int j = 0; j < DIM; j++) {
            #pragma HLS PIPELINE II=1
            phi_row[j] = phi[j][i];
        }

        // Calcular ratios para condiciones de contorno
        complex_t ratio_y0, ratio_yn;
        if (hls::hypot(phi_row[1].real(), phi_row[1].imag()) < eps) {
            ratio_y0.real(1.0);
            ratio_y0.imag(0);
        } else {
            ratio_y0 = phi_row[0] / phi_row[1];
        }

        if (hls::hypot(phi_row[DIM-2].real(), phi_row[DIM-2].imag()) < eps) {
            ratio_yn.real(1.0);
            ratio_yn.imag(0);
        } else {
            ratio_yn = phi_row[DIM-1] / phi_row[DIM-2];
        }

        // Calcular coeficientes para matriz B
        // -2 * ung == -ung -ung
        complex_t dp1_B = (-ung - ung) + complex_t(1.0, 0) + ung * ratio_y0;
        complex_t dp2_B = (-ung - ung) + complex_t(1.0, 0) + ung * ratio_yn;
        complex_t dp_B = (-ung - ung) + complex_t(1.0, 0);
        complex_t do_B = ung;

        // Calcular vector b
        compute_b_vector(dp_B, dp1_B, dp2_B, do_B, phi_row, b_vec);

        // Calcular coeficientes para matriz A
        // 2 * ung == ung + ung
        complex_t dp1_A = ung + ung + complex_t(1.0, 0) - ung * ratio_y0;
        complex_t dp2_A = ung + ung + complex_t(1.0, 0) - ung * ratio_yn;
        complex_t dp_A = ung + ung + complex_t(1.0, 0);
        complex_t do_A = -ung;

        // Resolver sistema tridiagonal
        custom_thomas_solver(dp_A, dp1_A, dp2_A, do_A, b_vec, result);

        // Guardar resultados
        store_results:
        for (int j = 0; j < DIM; j++) {
            #pragma HLS PIPELINE II=1
            phi_inter[j][i] = result[j];
        }
    }
}


/**
 * Applies Kerr nonlinear effect.
 */
void half_nonlinear(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]  // Campo de salida
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_out offset=slave bundle=GMEM1

    process_nonlinear:
    for (int i = 0; i < DIM; i++) {
        #pragma HLS PIPELINE II=1

        // Calcular el módulo al cuadrado
        float abs_squared = phi[i].real() * phi[i].real() + phi[i].imag() * phi[i].imag();

        // Calcular la fase
        float phase_angle = phase_const * abs_squared;
        float sin_val, cos_val;
        hls::sincos(phase_angle, &sin_val, &cos_val);

        complex_t phase(cos_val, sin_val);

        // Aplicar la fase
        phi_out[i] = phi[i] * phase;
    }
}

/**
 * Applies linear absorption.
 */
void half_linear_absorption(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]  // Campo de salida
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_out offset=slave bundle=GMEM1

    process_absorption:
    for (int i = 0; i < DIM; i++) {
        #pragma HLS PIPELINE II=1
        phi_out[i] = phi[i] * attenuation;
    }
}

/**
 * Applies two-photon absorption.
 */
void half_2photon_absorption(
    complex_t phi[DIM],      // Campo de entrada
    complex_t phi_out[DIM]  // Campo de salida
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_out offset=slave bundle=GMEM1

    process_2photon:
    for (int i = 0; i < DIM; i++) {
        #pragma HLS PIPELINE II=1

        // Calcular el módulo al cuadrado
        float abs_squared = phi[i].real() * phi[i].real() + phi[i].imag() * phi[i].imag();

        // Calcular el factor de atenuación
        float attenuation = hls::exp( tpa_const * abs_squared);

        // Aplicar la atenuación
        phi_out[i] = phi[i] * attenuation;
    }
}


/**
 * Complete propagation step combining all operators.
 */
void propagation_step(
    complex_t phi_in[DIM][DIM],      // Campo de entrada
    complex_t phi_out[DIM][DIM]      // Campo de salida
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi_in offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_out offset=slave bundle=GMEM1

    #pragma HLS DATAFLOW

    // Buffers intermedios
    complex_t phi_after_absorption[DIM][DIM];
    complex_t phi_after_nonlinear[DIM][DIM];
    complex_t phi_after_adi_x[DIM][DIM];

    #pragma HLS ARRAY_PARTITION variable=phi_after_absorption cyclic factor=2 dim=2
    #pragma HLS ARRAY_PARTITION variable=phi_after_nonlinear cyclic factor=2 dim=2
    #pragma HLS ARRAY_PARTITION variable=phi_after_adi_x cyclic factor=2 dim=2

    // Aplicar absorción lineal
    apply_linear_absorption:
    for (int j = 0; j < DIM; j++) {
        for (int i = 0; i < DIM; i++) {
            #pragma HLS PIPELINE II=1
            phi_after_absorption[j][i] = phi_in[j][i] * attenuation;
        }
    }

    // Aplicar absorción de dos fotones
    apply_2photon_absorption:
    for (int j = 0; j < DIM; j++) {
        for (int i = 0; i < DIM; i++) {
            #pragma HLS PIPELINE II=1
            float abs_squared = phi_after_absorption[j][i].real() * phi_after_absorption[j][i].real() + 
                               phi_after_absorption[j][i].imag() * phi_after_absorption[j][i].imag();
            float tpa_attenuation = hls::exp(tpa_const * abs_squared);
            phi_after_absorption[j][i] = phi_after_absorption[j][i] * tpa_attenuation;
        }
    }

    // Aplicar efecto no lineal
    apply_nonlinear:
    for (int j = 0; j < DIM; j++) {
        for (int i = 0; i < DIM; i++) {
            #pragma HLS PIPELINE II=1
            float abs_squared = phi_after_absorption[j][i].real() * phi_after_absorption[j][i].real() + 
                               phi_after_absorption[j][i].imag() * phi_after_absorption[j][i].imag();
            float phase_angle = phase_const * abs_squared;
            float sin_val, cos_val;
            hls::sincos(phase_angle, &sin_val, &cos_val);
            complex_t phase(cos_val, sin_val);
            phi_after_nonlinear[j][i] = phi_after_absorption[j][i] * phase;
        }
    }

    // Aplicar ADI en dirección X
    adi_x(phi_after_nonlinear, phi_after_adi_x);

    // Aplicar ADI en dirección Y
    adi_y(phi_after_adi_x, phi_out);
}
