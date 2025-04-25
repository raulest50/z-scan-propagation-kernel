#include <hls_x_complex.h>
#include <ap_int.h>
#include "bpm_split_step.h"

#define MAX_ND_Z 600     // Maximum number of z-steps (NDZ)


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
    float nalongZ[MAX_ND_Z], float n2alongZ[MAX_ND_Z],
    hls::x_complex<float> PHI_m[DIM][DIM] // 
    )
{
    // Interface pragmas for AXI access:
    #pragma HLS INTERFACE m_axi port=PHI_m0 depth=2304 offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=PHI_m0 bundle=control

    #pragma HLS INTERFACE s_axilite port=k bundle=control
    #pragma HLS INTERFACE s_axilite port=NDX bundle=control
    #pragma HLS INTERFACE s_axilite port=NDY bundle=control
    #pragma HLS INTERFACE s_axilite port=NDZ bundle=control
    #pragma HLS INTERFACE s_axilite port=DX bundle=control
    #pragma HLS INTERFACE s_axilite port=DY bundle=control
    #pragma HLS INTERFACE s_axilite port=DZ bundle=control

    #pragma HLS INTERFACE m_axi port=nalongZ depth=600 offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=nalongZ bundle=control
    #pragma HLS INTERFACE m_axi port=n2alongZ depth=600 offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=n2alongZ bundle=control

    #pragma HLS INTERFACE m_axi port=PHI_m depth=2304 offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=PHI_m bundle=control

    #pragma HLS INTERFACE s_axilite port=return bundle=control

    
    


}



}
