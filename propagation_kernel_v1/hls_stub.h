// hls_stub.h
#ifndef HLS_STUB_H
#define HLS_STUB_H
#include <complex>
#include <cmath>
namespace hls {
template<typename T>
using x_complex = std::complex<T>;

inline float exp(float x) { return std::exp(x); }
inline double exp(double x) { return std::exp(x); }

inline float hypot(float x, float y) { return std::hypot(x, y); }
inline double hypot(double x, double y) { return std::hypot(x, y); }

inline void sincos(float x, float* s, float* c) { *s = std::sin(x); *c = std::cos(x); }
inline void sincos(double x, double* s, double* c) { *s = std::sin(x); *c = std::cos(x); }
}
#endif
