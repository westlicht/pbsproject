#pragma once

#include "core/Common.h"

namespace pbs {

class Mesh;

// Creates an isosurface from a list of particles.
class ParticleMesher {
public:
    struct Parameters {
        int supportParticles;
        int particlesPerUnitVolume;
        float restDensity;
        float restSpacing;
        float particleMass;
        float h;
        float isoLevel = 0.5f;
    };

    static Mesh createMeshIsotropic(const MatrixXf &positions, const Box3f &bounds, const Vector3i &cells, const Parameters &params);
    static Mesh createMeshAnisotropic(MatrixXf &positions, const Box3f &bounds, const Vector3i &cells, const Parameters &params);
};

} // namespace pbs
