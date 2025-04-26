#include <hls_x_complex.h>
#include <ap_int.h>
#include "bpm_split_step.h"

#define MAX_ND_Z 600     // Maximum number of z-steps (NDZ)
#define MATRIX_SIZE DIM*DIM
#define MAX_STOPS_SIZE 128


extern "C" {

// The accelerator function for the propagation loop.
// This function takes as input the initial beam profile (PHI_m0),
// the propagation constant (k), grid sizes (NDX, NDY, NDZ),
// step sizes (DX, DY, DZ), the number of z-points to save,
// and the arrays nalongZ and n2alongZ (each of length MAX_ND_Z).
// It returns the final propagated beam profile (PHI_m),
// the saved intermediate profiles (PHI_m_alongZ), and the saved z positions (z_to_save).
void bpm_prop_accel(
    hls::x_complex<float> PHI_m0[DIM][DIM],
    float k,
    int NDX, int NDY, int NDZ, float DX, float DY, float DZ,
    float n_medium, float n_sample, float n2_sample,
    int stops_locations[MAX_STOPS_SIZE], int Nstops,
    hls::x_complex<float> PHI_m[DIM][DIM] // 
    )
{
    // Interface pragmas for AXI access:
    #pragma HLS INTERFACE m_axi port=PHI_m0 depth=MATRIX_SIZE offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=PHI_m0 bundle=control

    #pragma HLS INTERFACE s_axilite port=k bundle=control
    
    #pragma HLS INTERFACE s_axilite port=NDX bundle=control
    #pragma HLS INTERFACE s_axilite port=NDY bundle=control
    #pragma HLS INTERFACE s_axilite port=NDZ bundle=control
    #pragma HLS INTERFACE s_axilite port=DX bundle=control
    #pragma HLS INTERFACE s_axilite port=DY bundle=control
    #pragma HLS INTERFACE s_axilite port=DZ bundle=control

    #pragma HLS INTERFACE s_axilite port=n_medium bundle=control
    #pragma HLS INTERFACE s_axilite port=n_sample bundle=control
    #pragma HLS INTERFACE s_axilite port=n2_sample bundle=control

    #pragma HLS INTERFACE m_axi port=stops_locations depth=MAX_STOPS_SIZE offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=stops_locations bundle=control

    #pragma HLS INTERFACE s_axilite port=Nstops bundle=control

    #pragma HLS INTERFACE m_axi port=PHI_m depth=MATRIX_SIZE offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=PHI_m bundle=control

    #pragma HLS INTERFACE s_axilite port=return bundle=control
    


    


    float n0 = 1, n2=1.4;

    bpm_1st_half(
        PHI_m0, PHI_m0, 
        k, 
        n0, n2, 
        NDX, NDY, DX, DY, DZ, 
        PHI_m
        );
    


}



}
