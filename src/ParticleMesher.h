#pragma once

#include "Common.h"

namespace pbs {

class Mesh;

// Creates an isosurface from a list of particles.
class ParticleMesher {
public:
    static Mesh createMeshIsotropic(const MatrixXf &positions, const Box3f &bounds, const Vector3i &cells, float smoothRadius, float normalization, float isoLevel = 0.5f);

};

} // namespace pbs
