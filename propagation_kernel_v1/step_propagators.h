
#include <hls_x_complex.h>
#include "hls_math.h"

#define Nz 120     // Maximum number of z-steps (NDZ)
#define DIM 256
#define MATRIX_SIZE DIM*DIM

typedef hls::x_complex<float> complex_t;