// hls_stub.h
#ifndef HLS_STUB_H
#define HLS_STUB_H
#include <complex>
#include <cmath>
#include <queue>
#include <ap_fixed.h>
namespace hls {
template<typename T>
using x_complex = std::complex<T>;

// ---------------------------------------------------------------------
// Minimal hls::stream stub used for host compilation
// ---------------------------------------------------------------------
template<typename T>
class stream {
    std::queue<T> q;
public:
    inline bool empty() const { return q.empty(); }
    inline void write(const T& v) { q.push(v); }
    inline T read() { T v = q.front(); q.pop(); return v; }
};

inline float exp(float x) { return std::exp(x); }
inline double exp(double x) { return std::exp(x); }
inline ap_fixed<32,16> exp(ap_fixed<32,16> x) {
    return ap_fixed<32,16>(std::exp(static_cast<float>(x)));
}

inline float hypot(float x, float y) { return std::hypot(x, y); }
inline double hypot(double x, double y) { return std::hypot(x, y); }
inline ap_fixed<32,16> hypot(ap_fixed<32,16> x, ap_fixed<32,16> y) {
    return ap_fixed<32,16>(std::hypot(static_cast<float>(x), static_cast<float>(y)));
}

inline void sincos(float x, float* s, float* c) { *s = std::sin(x); *c = std::cos(x); }
inline void sincos(double x, double* s, double* c) { *s = std::sin(x); *c = std::cos(x); }
inline void sincos(ap_fixed<32,16> x, ap_fixed<32,16>* s, ap_fixed<32,16>* c) {
    float sx = std::sin(static_cast<float>(x));
    float cx = std::cos(static_cast<float>(x));
    *s = ap_fixed<32,16>(sx);
    *c = ap_fixed<32,16>(cx);
}
}
#endif
