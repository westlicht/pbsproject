#pragma once

#include "core/Common.h"
#include "core/Vector.h"

namespace pbs {

// SPH Kernels
// Kernels are split into constant and variable part.
// Note: r is supposed to be within the filter support (e.g. |r| <= h)
// Arguments are as follows:
// r  = displacement vector
// r2 = |r|^2 (squared norm of r)
// rn = |r|   (norm of r)
struct Kernel {
    float h;
    float h2;
    float halfh;

    void init(float h_) {
        h = h_;
        h2 = sqr(h);
        halfh = 0.5f * h;
        poly6Constant = 315.f / (64.f * M_PI * std::pow(h, 9.f));
        poly6GradConstant = -945.f / (32.f * M_PI * std::pow(h, 9.f));
        poly6LaplaceConstant = -945.f / (32.f * M_PI * std::pow(h, 9.f));
        spikyConstant = 15.f / (M_PI * std::pow(h, 6.f));
        spikyGradConstant = -45.f / (M_PI * std::pow(h, 6.f));
        spikyLaplaceConstant = -90.f / (M_PI * std::pow(h, 6.f));
        viscosityLaplaceConstant = 45.f / (M_PI * std::pow(h, 6.f));

        surfaceTensionConstant = 32.f / (M_PI * std::pow(h, 9.f));
        surfaceTensionOffset = -std::pow(h, 6.f) / 64.f;
    }

    float poly6Constant;
    inline float poly6(float r2) const {
        return cube(h2 - r2);
    }

    float poly6GradConstant;
    inline Vector3f poly6Grad(const Vector3f &r, float r2) const {
        return sqr(h2 - r2) * r;
    }

    float poly6LaplaceConstant;
    inline float poly6Laplace(float r2) {
        return (h2 - r2) * (3.f * h2 - 7.f * r2);
    }

    float spikyConstant;
    inline float spiky(float rn) const {
        return cube(h - rn);
    }

    float spikyGradConstant;
    inline Vector3f spikyGrad(const Vector3f &r, float rn) const {
        return sqr(h - rn) * r * (1.f / rn);
    }

    float spikyLaplaceConstant;
    inline float spikyLaplace(float rn) const {
        return (h - rn) * (h - 2.f * rn) / rn;
    }

    float viscosityLaplaceConstant;
    inline float viscosityLaplace(float rn) const {
        return (h - rn);
    }

    float surfaceTensionConstant;
    float surfaceTensionOffset;
    inline float surfaceTension(float rn) const {
        if (rn < halfh) {
            return 2.f * cube(h - rn) * cube(rn) + surfaceTensionOffset;
        } else {
            return cube(h - rn) * cube(rn);
        }
    }
};

} // namespace pbs
