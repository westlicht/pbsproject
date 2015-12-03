#pragma once

#include "core/Common.h"

#include <vector>

namespace pbs {

class Mesh;

// Helpers to generate particles from meshes.
class ParticleGenerator {
public:
    struct Result {
        std::vector<Vector3f> positions;
        std::vector<Vector3f> normals;        
    };

    // Generates a uniform particle distribution on the surface of a box.
    static Result generateSurfaceParticles(const Box3f &box, float particleRadius);

    // Generates a uniform particle distribution on the surface of a mesh.
    static Result generateSurfaceParticles(const Mesh &mesh, float particleRadius, int cells = 100);
};

} // namespace pbs
