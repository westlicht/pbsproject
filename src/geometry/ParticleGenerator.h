#pragma once

#include "core/Common.h"

#include <vector>

namespace pbs {

class Mesh;

// Helpers to generate particles from meshes.
class ParticleGenerator {
public:
    // Generates a uniform particle distribution on the surface of a mesh.
    static std::vector<Vector3f> generateSurfaceParticles(const Mesh &mesh, float density, int cells = 100);
};

} // namespace pbs
