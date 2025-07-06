/**
 * @file hls_step_operators.cpp
 * @brief Implementation of HLS functions for wave propagation step operators
 * 
 * This file contains the HLS implementations of the step operators for wave propagation
 * based on the Python implementations in step_operators.py from the deep_tissue_imaging package.
 */

#include "hls_step_operators.h"

/**
 * Solves a tridiagonal system with special structure using Thomas algorithm.
 */
void custom_thomas_solver(
    complex_t dp,        // Diagonal principal (valor común)
    complex_t dp1,       // Primer elemento de la diagonal principal
    complex_t dp2,       // Último elemento de la diagonal principal
    complex_t do_val,    // Valor de las diagonales secundarias
    complex_t b[MAX_SIZE],      // Vector del lado derecho
    complex_t x[MAX_SIZE],      // Vector solución
    int n                // Tamaño del sistema
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=dp,dp1,dp2,do_val,n bundle=CTRL
    #pragma HLS INTERFACE m_axi port=b offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=x offset=slave bundle=GMEM1
    
    // Buffers locales para mejorar el rendimiento
    complex_t c_prime[MAX_SIZE-1];
    complex_t d_prime[MAX_SIZE];
    
    #pragma HLS ARRAY_PARTITION variable=c_prime cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=d_prime cyclic factor=4
    
    // Eliminación hacia adelante
    // Primera fila
    c_prime[0] = do_val / dp1;
    d_prime[0] = b[0] / dp1;
    
    // Filas intermedias
    forward_elimination:
    for (int i = 1; i < n-1; i++) {
        #pragma HLS PIPELINE II=1
        complex_t denominator = dp - do_val * c_prime[i-1];
        c_prime[i] = do_val / denominator;
        d_prime[i] = (b[i] - do_val * d_prime[i-1]) / denominator;
    }
    
    // Última fila
    d_prime[n-1] = (b[n-1] - do_val * d_prime[n-2]) / (dp2 - do_val * c_prime[n-2]);
    
    // Sustitución hacia atrás
    x[n-1] = d_prime[n-1];
    
    back_substitution:
    for (int i = n-2; i >= 0; i--) {
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
    complex_t x0[MAX_SIZE],     // Vector de entrada
    complex_t b[MAX_SIZE],      // Vector resultado
    int n                // Tamaño del vector
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=dp,dp1,dp2,do_val,n bundle=CTRL
    #pragma HLS INTERFACE m_axi port=x0 offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=b offset=slave bundle=GMEM1
    
    // Primera fila
    b[0] = dp1 * x0[0] + do_val * x0[1];
    
    // Filas intermedias
    middle_rows:
    for (int i = 1; i < n-1; i++) {
        #pragma HLS PIPELINE II=1
        b[i] = do_val * x0[i-1] + dp * x0[i] + do_val * x0[i+1];
    }
    
    // Última fila
    b[n-1] = do_val * x0[n-2] + dp2 * x0[n-1];
}

/**
 * Implements the ADI method in X direction.
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
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=Ny,Nx,eps,k,dz,dx bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_inter offset=slave bundle=GMEM1
    
    // Buffer local para una columna
    complex_t phi_col[MAX_NX];
    complex_t b_vec[MAX_NX];
    complex_t result[MAX_NX];
    
    #pragma HLS ARRAY_PARTITION variable=phi_col cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=b_vec cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=result cyclic factor=4
    
    // Calcular el factor ung
    complex_t ung;
    ung.real(0);
    ung.imag(dz / (4 * k * dx * dx));
    
    process_rows:
    for (int j = 0; j < Ny; j++) {
        #pragma HLS DATAFLOW
        
        // Cargar columna en buffer local
        load_column:
        for (int i = 0; i < Nx; i++) {
            #pragma HLS PIPELINE II=1
            phi_col[i] = phi[j][i];
        }
        
        // Calcular ratios para condiciones de contorno
        complex_t ratio_x0, ratio_xn;
        if (hls::abs(phi_col[1]) < eps) {
            ratio_x0.real(1.0);
            ratio_x0.imag(0);
        } else {
            ratio_x0 = phi_col[0] / phi_col[1];
        }
        
        if (hls::abs(phi_col[Nx-2]) < eps) {
            ratio_xn.real(1.0);
            ratio_xn.imag(0);
        } else {
            ratio_xn = phi_col[Nx-1] / phi_col[Nx-2];
        }
        
        // Calcular coeficientes para matriz B
        complex_t dp1_B = -2 * ung + complex_t(1.0, 0) + ung * ratio_x0;
        complex_t dp2_B = -2 * ung + complex_t(1.0, 0) + ung * ratio_xn;
        complex_t dp_B = -2 * ung + complex_t(1.0, 0);
        complex_t do_B = ung;
        
        // Calcular vector b
        compute_b_vector(dp_B, dp1_B, dp2_B, do_B, phi_col, b_vec, Nx);
        
        // Calcular coeficientes para matriz A
        complex_t dp1_A = 2 * ung + complex_t(1.0, 0) - ung * ratio_x0;
        complex_t dp2_A = 2 * ung + complex_t(1.0, 0) - ung * ratio_xn;
        complex_t dp_A = 2 * ung + complex_t(1.0, 0);
        complex_t do_A = -ung;
        
        // Resolver sistema tridiagonal
        custom_thomas_solver(dp_A, dp1_A, dp2_A, do_A, b_vec, result, Nx);
        
        // Guardar resultados
        store_results:
        for (int i = 0; i < Nx; i++) {
            #pragma HLS PIPELINE II=1
            phi_inter[j][i] = result[i];
        }
    }
}

/**
 * Implements the ADI method in Y direction.
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
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=Nx,Ny,eps,k,dz,dy bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_inter offset=slave bundle=GMEM1
    
    // Buffer local para una fila
    complex_t phi_row[MAX_NY];
    complex_t b_vec[MAX_NY];
    complex_t result[MAX_NY];
    
    #pragma HLS ARRAY_PARTITION variable=phi_row cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=b_vec cyclic factor=4
    #pragma HLS ARRAY_PARTITION variable=result cyclic factor=4
    
    // Calcular el factor ung
    complex_t ung;
    ung.real(0);
    ung.imag(dz / (4 * k * dy * dy));
    
    process_columns:
    for (int i = 0; i < Nx; i++) {
        #pragma HLS DATAFLOW
        
        // Cargar fila en buffer local
        load_row:
        for (int j = 0; j < Ny; j++) {
            #pragma HLS PIPELINE II=1
            phi_row[j] = phi[j][i];
        }
        
        // Calcular ratios para condiciones de contorno
        complex_t ratio_y0, ratio_yn;
        if (hls::abs(phi_row[1]) < eps) {
            ratio_y0.real(1.0);
            ratio_y0.imag(0);
        } else {
            ratio_y0 = phi_row[0] / phi_row[1];
        }
        
        if (hls::abs(phi_row[Ny-2]) < eps) {
            ratio_yn.real(1.0);
            ratio_yn.imag(0);
        } else {
            ratio_yn = phi_row[Ny-1] / phi_row[Ny-2];
        }
        
        // Calcular coeficientes para matriz B
        complex_t dp1_B = -2 * ung + complex_t(1.0, 0) + ung * ratio_y0;
        complex_t dp2_B = -2 * ung + complex_t(1.0, 0) + ung * ratio_yn;
        complex_t dp_B = -2 * ung + complex_t(1.0, 0);
        complex_t do_B = ung;
        
        // Calcular vector b
        compute_b_vector(dp_B, dp1_B, dp2_B, do_B, phi_row, b_vec, Ny);
        
        // Calcular coeficientes para matriz A
        complex_t dp1_A = 2 * ung + complex_t(1.0, 0) - ung * ratio_y0;
        complex_t dp2_A = 2 * ung + complex_t(1.0, 0) - ung * ratio_yn;
        complex_t dp_A = 2 * ung + complex_t(1.0, 0);
        complex_t do_A = -ung;
        
        // Resolver sistema tridiagonal
        custom_thomas_solver(dp_A, dp1_A, dp2_A, do_A, b_vec, result, Ny);
        
        // Guardar resultados
        store_results:
        for (int j = 0; j < Ny; j++) {
            #pragma HLS PIPELINE II=1
            phi_inter[j][i] = result[j];
        }
    }
}

/**
 * Applies Kerr nonlinear effect.
 */
void half_nonlinear(
    complex_t phi[MAX_SIZE],      // Campo de entrada
    complex_t phi_out[MAX_SIZE],  // Campo de salida
    float k_sample,               // Número de onda en la muestra
    float n2_sample,              // Índice no lineal
    float dz,                     // Paso en z
    int size                      // Tamaño del campo
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=k_sample,n2_sample,dz,size bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_out offset=slave bundle=GMEM1
    
    // Constante para el cálculo de fase
    float phase_const = k_sample * n2_sample * dz / 2.0f;
    
    process_nonlinear:
    for (int i = 0; i < size; i++) {
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
    complex_t phi[MAX_SIZE],      // Campo de entrada
    complex_t phi_out[MAX_SIZE],  // Campo de salida
    float alpha,                  // Coeficiente de absorción
    float dz,                     // Paso en z
    int size                      // Tamaño del campo
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=alpha,dz,size bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_out offset=slave bundle=GMEM1
    
    // Calcular el factor de atenuación
    float attenuation = hls::exp(-alpha * dz / 4.0f);
    
    process_absorption:
    for (int i = 0; i < size; i++) {
        #pragma HLS PIPELINE II=1
        phi_out[i] = phi[i] * attenuation;
    }
}

/**
 * Applies two-photon absorption.
 */
void half_2photon_absorption(
    complex_t phi[MAX_SIZE],      // Campo de entrada
    complex_t phi_out[MAX_SIZE],  // Campo de salida
    float beta,                   // Coeficiente de absorción de dos fotones
    float dz,                     // Paso en z
    int size                      // Tamaño del campo
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=beta,dz,size bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_out offset=slave bundle=GMEM1
    
    process_2photon:
    for (int i = 0; i < size; i++) {
        #pragma HLS PIPELINE II=1
        
        // Calcular el módulo al cuadrado
        float abs_squared = phi[i].real() * phi[i].real() + phi[i].imag() * phi[i].imag();
        
        // Calcular el factor de atenuación
        float attenuation = hls::exp(-beta * dz / 4.0f * abs_squared);
        
        // Aplicar la atenuación
        phi_out[i] = phi[i] * attenuation;
    }
}

/**
 * Complete propagation step combining all operators.
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
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=Nx,Ny,eps,k,k_sample,n2_sample,alpha,beta,dz,dx,dy bundle=CTRL
    #pragma HLS INTERFACE m_axi port=phi_in offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi port=phi_out offset=slave bundle=GMEM1
    
    #pragma HLS DATAFLOW
    
    // Buffers intermedios
    complex_t phi_after_absorption[MAX_NY][MAX_NX];
    complex_t phi_after_nonlinear[MAX_NY][MAX_NX];
    complex_t phi_after_adi_x[MAX_NY][MAX_NX];
    
    #pragma HLS ARRAY_PARTITION variable=phi_after_absorption cyclic factor=2 dim=2
    #pragma HLS ARRAY_PARTITION variable=phi_after_nonlinear cyclic factor=2 dim=2
    #pragma HLS ARRAY_PARTITION variable=phi_after_adi_x cyclic factor=2 dim=2
    
    // Aplicar absorción lineal
    apply_linear_absorption:
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            #pragma HLS PIPELINE II=1
            phi_after_absorption[j][i] = phi_in[j][i] * hls::exp(-alpha * dz / 4.0f);
        }
    }
    
    // Aplicar absorción de dos fotones
    apply_2photon_absorption:
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            #pragma HLS PIPELINE II=1
            float abs_squared = phi_after_absorption[j][i].real() * phi_after_absorption[j][i].real() + 
                               phi_after_absorption[j][i].imag() * phi_after_absorption[j][i].imag();
            float attenuation = hls::exp(-beta * dz / 4.0f * abs_squared);
            phi_after_absorption[j][i] = phi_after_absorption[j][i] * attenuation;
        }
    }
    
    // Aplicar efecto no lineal
    apply_nonlinear:
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            #pragma HLS PIPELINE II=1
            float abs_squared = phi_after_absorption[j][i].real() * phi_after_absorption[j][i].real() + 
                               phi_after_absorption[j][i].imag() * phi_after_absorption[j][i].imag();
            float phase_angle = k_sample * n2_sample * dz / 2.0f * abs_squared;
            float sin_val, cos_val;
            hls::sincos(phase_angle, &sin_val, &cos_val);
            complex_t phase(cos_val, sin_val);
            phi_after_nonlinear[j][i] = phi_after_absorption[j][i] * phase;
        }
    }
    
    // Aplicar ADI en dirección X
    adi_x(phi_after_nonlinear, phi_after_adi_x, Ny, Nx, eps, k, dz, dx);
    
    // Aplicar ADI en dirección Y
    adi_y(phi_after_adi_x, phi_out, Nx, Ny, eps, k, dz, dy);
}