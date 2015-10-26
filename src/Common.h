#pragma once

#include <Eigen/Core>

#include <tinyformat.h>

#include <cmath>
#include <exception>

namespace pbs {

/* Forward declarations */
template <typename Scalar, int Dimension>  struct TVector;
template <typename Vector>                 struct TBox;

typedef TVector<float, 1>       Vector1f;
typedef TVector<float, 2>       Vector2f;
typedef TVector<float, 3>       Vector3f;
typedef TVector<float, 4>       Vector4f;
typedef TVector<double, 1>      Vector1d;
typedef TVector<double, 2>      Vector2d;
typedef TVector<double, 3>      Vector3d;
typedef TVector<double, 4>      Vector4d;
typedef TVector<int, 1>         Vector1i;
typedef TVector<int, 2>         Vector2i;
typedef TVector<int, 3>         Vector3i;
typedef TVector<int, 4>         Vector4i;
typedef TBox<Vector1f>          Box1f;
typedef TBox<Vector2f>          Box2f;
typedef TBox<Vector3f>          Box3f;
typedef TBox<Vector4f>          Box4f;
typedef TBox<Vector1d>          Box1d;
typedef TBox<Vector2d>          Box2d;
typedef TBox<Vector3d>          Box3d;
typedef TBox<Vector4d>          Box4d;
typedef TBox<Vector1i>          Box1i;
typedef TBox<Vector2i>          Box2i;
typedef TBox<Vector3i>          Box3i;
typedef TBox<Vector4i>          Box4i;

template<typename T>
static constexpr T sqr(T x) { return x*x; }

template<typename T>
static constexpr T cube(T x) { return x*x*x; }

template<typename T>
static inline T clamp(T x, T lo, T hi)  {
    return std::max(lo, std::min(hi, x));
}

template<typename... Args>
static inline void DBG(const char *fmt, const Args&... args) {
    std::cout << tfm::format(fmt, args...) << std::endl;
}

template<typename... Args>
static inline void handleAssert(bool cond, const char *fmt, const Args&... args) {
    if (!cond) {
        std::string msg = tfm::format(fmt, args...);
        std::cout << msg << std::endl;
        throw std::runtime_error(msg);
    }
}

#define ASSERT(_cond_, _fmt_, _args_...) \
handleAssert(_cond_, "%s(%d) ASSERT: " _fmt_, __FILE__, __LINE__, ##_args_)


} // namespace pbs
