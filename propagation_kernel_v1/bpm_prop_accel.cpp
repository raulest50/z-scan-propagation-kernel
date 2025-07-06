#include "step_propagators.h"
#include <ap_int.h>


extern "C" {

// The accelerator function for the propagation loop.
// This function takes as input the initial beam profile (PHI_m0),
// It returns the final propagated beam profile (PHI_m),
// the saved intermediate profiles (PHI_m_alongZ), and the saved z positions (z_to_save).
void bpm_prop_accel(
    hls::x_complex<float> PHI_m0[DIM][DIM],
    hls::x_complex<float> PHI_m[DIM][DIM] // 
    )
{
    // Interface pragmas for AXI access:
    #pragma HLS INTERFACE m_axi port=PHI_m0 depth=MATRIX_SIZE offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=PHI_m0 bundle=control

    #pragma HLS INTERFACE m_axi port=PHI_m depth=MATRIX_SIZE offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=PHI_m bundle=control

    #pragma HLS INTERFACE s_axilite port=return bundle=control

    for( int k = 0; k<Nz; k++ ){
        
        // aplicar Adix
        // aplicar halfSteps

        // aplicar Adiy
        // aplicar HalfSteps
        
    }


}



}
