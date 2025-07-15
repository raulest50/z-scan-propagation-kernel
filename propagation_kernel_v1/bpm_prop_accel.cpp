// bpm_prop_accel.cpp

#include "step_propagators.h"
#include <ap_int.h>

// Size of one DIM x DIM complex matrix used by the AXI interfaces
#ifndef MATRIX_SIZE
#define MATRIX_SIZE (DIM*DIM)
#endif

// Number of propagation steps along the Z axis
#ifndef Nz
#define Nz 120
#endif


extern "C" {

// The accelerator function for the propagation loop.
// This function takes as input the initial beam profile (PHI_m0),
// It returns the final propagated beam profile (PHI_m) after 120 propagation steps.
void bpm_prop_accel(
    complex_t PHI_m0[DIM][DIM],
    complex_t PHI_m[DIM][DIM]
    )
{
    // Interface pragmas for AXI access:
    #pragma HLS INTERFACE m_axi port=PHI_m0 depth=MATRIX_SIZE offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=PHI_m0 bundle=control

    #pragma HLS INTERFACE m_axi port=PHI_m depth=MATRIX_SIZE offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=PHI_m bundle=control

    #pragma HLS INTERFACE s_axilite port=return bundle=control

    // Create local arrays for the propagation process
    static complex_t current[DIM][DIM];
    static complex_t next[DIM][DIM];

    // Bind to URAM for better performance
    #pragma HLS bind_storage variable=current type=ram_1p impl=uram
    #pragma HLS bind_storage variable=next type=ram_1p impl=uram

    // Copy initial beam profile to current array
    for (int j = 0; j < DIM; j++) {
        for (int i = 0; i < DIM; i++) {
            #pragma HLS PIPELINE II=1
            current[i][j] = PHI_m0[i][j];
        }
    }

    // Perform 120 propagation steps
    for (int k = 0; k < Nz; k++) {
        // Apply one complete propagation step
        propagation_step(current, next);

        // Copy result back to current for next iteration
        for (int j = 0; j < DIM; j++) {
            for (int i = 0; i < DIM; i++) {
                #pragma HLS PIPELINE II=1
                current[i][j] = next[i][j];
            }
        }
    }

    // Copy final result to output array
    for (int j = 0; j < DIM; j++) {
        for (int i = 0; i < DIM; i++) {
            #pragma HLS PIPELINE II=1
            PHI_m[i][j] = current[i][j];
        }
    }


}



}
