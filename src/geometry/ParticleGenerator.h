#pragma once

#include "core/Common.h"

#include <vector>

namespace pbs {

class Mesh;

// Helpers to generate particle for boundaries and volumes.
class ParticleGenerator {
public:
    struct Boundary {
        std::vector<Vector3f> positions;
        std::vector<Vector3f> normals;        
    };

    static Boundary generateBoundaryBox(const Box3f &box, float particleRadius, bool flipNormals = false);
    static Boundary generateBoundarySphere(const Vector3f &position, float radius, float particleRadius);
    static Boundary generateBoundaryMesh(const Mesh &mesh, float particleRadius, int cells = 100);


    struct Volume {
        std::vector<Vector3f> positions;
    };

    static Volume generateVolumeBox(const Box3f &box, float particleRadius);
    static Volume generateVolumeSphere(const Vector3f &center, float radius, float particleRadius);
    static Volume generateVolumeMesh(const Mesh &mesh, float particleRadius);
};

} // namespace pbs
