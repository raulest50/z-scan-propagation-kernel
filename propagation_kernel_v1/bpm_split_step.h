// bpm_split_step.h
#ifndef BPM_SPLIT_STEP_H
#define BPM_SPLIT_STEP_H

#include <hls_x_complex.h>

#define DIM 46

// First half–step: diffraction + Kerr nonlinearity (with TBC)
void bpm_1st_half(
    hls::x_complex<float> PHI_m[DIM][DIM],
    hls::x_complex<float> PHI_m_auxNL[DIM][DIM],
    float                 k,
    float                 n0,
    float                 n2[DIM][DIM],
    int                   NDX,
    int                   NDY,
    float                 DX,
    float                 DY,
    float                 DZ,
    hls::x_complex<float> PHI_half[DIM][DIM]
);

// Second half–step: Kerr nonlinearity + diffraction (mirror of first)
void bpm_2nd_half(
    hls::x_complex<float> PHI_half[DIM][DIM],
    hls::x_complex<float> PHI_m_auxNL[DIM][DIM],
    float                 k,
    float                 n0,
    float                 n2[DIM][DIM],
    int                   NDX,
    int                   NDY,
    float                 DX,
    float                 DY,
    float                 DZ,
    hls::x_complex<float> PHI_out[DIM][DIM]
);

#endif // BPM_SPLIT_STEP_H
