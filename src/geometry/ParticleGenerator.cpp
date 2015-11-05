#include "ParticleGenerator.h"
#include "SDF.h"
#include "Mesh.h"

#include "core/Common.h"
#include "core/Vector.h"
#include "core/Box.h"

#include <pcg32.h>

#include <Eigen/Geometry>

namespace pbs {

std::vector<Vector3f> ParticleGenerator::generateSurfaceParticles(const Mesh &mesh, float density, int cells) {

    DBG("Generating surface particles ...");

    // Compute bounds of mesh and expand by 10%
    Box3f bounds = mesh.computeBounds();
    bounds = bounds.expanded(bounds.extents() * 0.1f);

    // Compute cell and grid size for signed distance field
    float cellSize = bounds.extents()[bounds.majorAxis()] / cells;
    Vector3i size(
        int(std::ceil(bounds.extents().x() / cellSize)),
        int(std::ceil(bounds.extents().y() / cellSize)),
        int(std::ceil(bounds.extents().z() / cellSize))
    );

    VoxelGrid<float> sdf(size);
    sdf.setOrigin(bounds.min);
    sdf.setCellSize(cellSize);

    // Build signed distance field
    DBG("Building signed distance field ...");
    SDF::build(mesh, sdf);

    // Generate initial point distribution
    DBG("Generating initial point distribution ...");
    pcg32 rng;
    std::vector<Vector3f> positions;
    float totalArea = 0.f;
    for (int i = 0; i < mesh.triangles().cols(); ++i) {
        const Vector3f &p0 = mesh.vertices().col(mesh.triangles()(0, i));
        const Vector3f &p1 = mesh.vertices().col(mesh.triangles()(1, i));
        const Vector3f &p2 = mesh.vertices().col(mesh.triangles()(2, i));
        Vector3f e0 = p1 - p0;
        Vector3f e1 = p2 - p0;
        float area = 0.5f * std::abs(e0.cross(e1).norm());
        totalArea += area;

        float n = density * area;
        int ni = int(std::floor(n));

        auto samplePoint = [&] () {
            float s = rng.nextFloat();
            float t = rng.nextFloat();
            float ss = std::sqrt(s);
            positions.emplace_back(p0 + e0 * t * ss + e1 * (1.f - ss));
        };

        for (int j = 0; j < ni; ++j) {
            samplePoint();
        }
        if (rng.nextFloat() < n) {
            samplePoint();
        }
    }
    DBG("Generated %d points", positions.size());

    // Relax point distribution
    DBG("Relaxing point distribution ...");

    // Choose radius to support roughly 10 neighbour particles
    float radius = std::sqrt(totalArea / positions.size() * 10.f / M_PI);
    float radius2 = sqr(radius);

    for (int iteration = 0; iteration < 10; ++iteration) {
        int count = 0;
        std::vector<Vector3f> velocities(positions.size(), Vector3f());
        // Relax positions
        for (size_t i = 0; i < positions.size(); ++i) {
            for (size_t j = i + 1; j < positions.size(); ++j) {
                Vector3f r = positions[j] - positions[i];
                float r2 = r.squaredNorm();
                if (r2 < radius2) {
                    r *= (1.f / std::sqrt(r2));
                    float weight = 0.01f * cube(1.f - r2 / radius2); // TODO use a proper kernel
                    velocities[i] -= weight * r;
                    velocities[j] += weight * r;
                    ++count;
                }
            }
            positions[i] += velocities[i];
        }
        // Reproject to surface
        for (size_t i = 0; i < positions.size(); ++i) {
            Vector3f p = sdf.toVoxelSpace(positions[i]);
            Vector3f n = sdf.gradient(p).normalized();
            positions[i] -= sdf.trilinear(p) * n;
        }
        DBG("avg neighbours = %d", 2 * count / positions.size());
    }

    return positions;
}

} // namespace pbs
