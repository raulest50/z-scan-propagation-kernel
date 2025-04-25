// bpm_split_step.cpp
#include "bpm_split_step.h"
#include <cmath>

//---------------------------------------------------------------------
// Helper: compute reflection coefficient for transparent BC
//---------------------------------------------------------------------
static inline hls::x_complex<float> computeGamma(
    hls::x_complex<float> u_out,
    hls::x_complex<float> u_in,
    bool                  positive_imag
) {
    if (u_out != hls::x_complex<float>(0,0) && u_in != hls::x_complex<float>(0,0)) {
        auto g = hls::x_complex<float>(1,0) / (u_out / u_in);
        float re = g.real(), im = g.imag();
        return positive_imag
             ? hls::x_complex<float>(re,  std::fabs(im))
             : hls::x_complex<float>(re, -std::fabs(im));
    }
    return hls::x_complex<float>(0,0);
}

//---------------------------------------------------------------------
// Thomas algorithm for tridiagonal systems of size N
//---------------------------------------------------------------------
static void thomasSolver(
    int N,
    hls::x_complex<float> a[],  // lower diag (a[0] unused)
    hls::x_complex<float> b[],  // main diag
    hls::x_complex<float> c[],  // upper diag (c[N-1] unused)
    hls::x_complex<float> d[],  // RHS
    hls::x_complex<float> x[]   // solution
) {
    #pragma HLS INLINE
    hls::x_complex<float> cp[DIM], dp[DIM];

    // forward sweep
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < N; ++i) {
        #pragma HLS PIPELINE II=1
        auto m = b[i] - a[i-1]*cp[i-1];
        cp[i] = (i < N-1 ? c[i] : hls::x_complex<float>(0,0)) / m;
        dp[i] = (d[i] - a[i-1]*dp[i-1]) / m;
    }
    // back substitution
    x[N-1] = dp[N-1];
    for (int i = N-2; i >= 0; --i) {
        #pragma HLS PIPELINE II=1
        x[i] = dp[i] - cp[i]*x[i+1];
    }
}

//---------------------------------------------------------------------
// First half-step: diffraction + Kerr nonlinearity (with TBC)
//---------------------------------------------------------------------
void bpm_1st_half(
    hls::x_complex<float> PHI_m[DIM][DIM],
    hls::x_complex<float> PHI_m_auxNL[DIM][DIM],
    float                 k,
    float                 n0,
    float                 n2,
    int                   NDX,
    int                   NDY,
    float                 DX,
    float                 DY,
    float                 DZ,
    hls::x_complex<float> PHI_half[DIM][DIM]
) {
    #pragma HLS INLINE off
    #pragma HLS PIPELINE

    // define complex forms of needed reals
    hls::x_complex<float> czero(0,0), cone(1,0);
    hls::x_complex<float> j1(0,1), minus_j1(0,-1);
    hls::x_complex<float> invDZ(1.0f/DZ, 0);
    hls::x_complex<float> invDX2(1.0f/(DX*DX), 0);
    hls::x_complex<float> invDY2(1.0f/(DY*DY), 0);

    // split-step coefficients
    // B = -j / (2 k n0)
    hls::x_complex<float> B = minus_j1 * hls::x_complex<float>(1.0f/(2.0f * k * n0), 0);
    hls::x_complex<float> C = B;
    // alpha = -B/(2 DX^2)
    hls::x_complex<float> alpha = B * hls::x_complex<float>(-0.5f,0) * invDX2;
    hls::x_complex<float> gamma = alpha;

    // partition these small arrays
    hls::x_complex<float> a[DIM], b[DIM], c[DIM], d[DIM], x_sol[DIM];
    #pragma HLS ARRAY_PARTITION variable=a     complete
    #pragma HLS ARRAY_PARTITION variable=b     complete
    #pragma HLS ARRAY_PARTITION variable=c     complete
    #pragma HLS ARRAY_PARTITION variable=d     complete
    #pragma HLS ARRAY_PARTITION variable=x_sol complete

    hls::x_complex<float> Dnl_arr[DIM];
    #pragma HLS ARRAY_PARTITION variable=Dnl_arr complete

    // for each row
    for (int l = 0; l <= NDY; ++l) {
        #pragma HLS LOOP_TRIPCOUNT min=1 max=DIM

        // build diagonals and Dnl
        for (int i = 0; i <= NDX; ++i) {
            #pragma HLS UNROLL
            float mag2 = PHI_m_auxNL[i][l].real()*PHI_m_auxNL[i][l].real()
                       + PHI_m_auxNL[i][l].imag()*PHI_m_auxNL[i][l].imag();
            // Dnl = -j * k * n2 * |phi|^2
            Dnl_arr[i] = minus_j1 * hls::x_complex<float>(k * n2 * mag2, 0);

            a[i] = alpha;
            c[i] = gamma;
            // b = invDZ + B*invDX2 - Dnl/4
            b[i] = invDZ + (B * invDX2) - (Dnl_arr[i] * hls::x_complex<float>(0.25f,0));
        }

        // transparent BC on X ends
        b[0]   += alpha * computeGamma(PHI_m[1][l],     PHI_m[0][l],   true);
        b[NDX] += alpha * computeGamma(PHI_m[NDX-1][l], PHI_m[NDX][l], false);

        // build RHS d
        for (int i = 0; i <= NDX; ++i) {
            #pragma HLS UNROLL
            auto up = (l>0)
                    ? PHI_m[i][l-1]
                    : (computeGamma(PHI_m[i][1], PHI_m[i][0], false) * PHI_m[i][0]);
            auto down = (l<NDY)
                    ? PHI_m[i][l+1]
                    : (computeGamma(PHI_m[i][NDY-1], PHI_m[i][NDY], true) * PHI_m[i][NDY]);

            // c/(2 DY2)
            auto c_fac = C * hls::x_complex<float>(0.5f,0) * invDY2;
            // mid = invDZ - C*invDY2 + Dnl/4
            auto mid = invDZ - (C * invDY2) + (Dnl_arr[i] * hls::x_complex<float>(0.25f,0));

            d[i] = c_fac*up
                 + mid*PHI_m[i][l]
                 + c_fac*down;
        }

        // solve and store
        thomasSolver(NDX+1, a, b, c, d, x_sol);
        for (int i = 0; i <= NDX; ++i) {
            PHI_half[i][l] = x_sol[i];
        }
    }
}

//---------------------------------------------------------------------
// Second half-step: Kerr nonlinearity + diffraction (mirror of first)
//---------------------------------------------------------------------
void bpm_2nd_half(
    hls::x_complex<float> PHI_half[DIM][DIM],
    hls::x_complex<float> PHI_m_auxNL[DIM][DIM],
    float                 k,
    float                 n0,
    float                 n2,
    int                   NDX,
    int                   NDY,
    float                 DX,
    float                 DY,
    float                 DZ,
    hls::x_complex<float> PHI_out[DIM][DIM]
) {
    // TODO: implement the reverse order:
    //   1) apply Dnl_arr = -j*k*n2*|phi_half|^2
    //   2) solve tridiagonal in Y-direction (swap loops)
}
