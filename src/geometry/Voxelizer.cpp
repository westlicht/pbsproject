#include "Voxelizer.h"
#include "Mesh.h"

#include "core/Box.h"
#include "core/Vector.h"

#include "Eigen/Geometry"

namespace pbs {

void Voxelizer::voxelize(const Mesh &mesh, float cellSize, Result &result) {
    Box3f &bounds = result.bounds;
    VoxelGrid<bool> &grid = result.grid;
    result.cellSize = cellSize;

    bounds = mesh.computeBounds();

    Vector3i cells(
        int(std::ceil(bounds.extents().x() / cellSize)),
        int(std::ceil(bounds.extents().y() / cellSize)),
        int(std::ceil(bounds.extents().z() / cellSize))
    );

    grid.resize(cells);
    grid.fill(false);

    auto index = [&bounds, &cellSize] (const Vector3f &p) {
        return Vector3i(
            int(std::floor((p.x() - bounds.min.x()) / cellSize)),
            int(std::floor((p.y() - bounds.min.y()) / cellSize)),
            int(std::floor((p.z() - bounds.min.z()) / cellSize))
        );
    };

    auto rasterizeTriangle = [index, &bounds, &cellSize, &grid] (const Vector3f &p0, const Vector3f &p1, const Vector3f &p2) {
        Box3i iBounds;
        iBounds.expandBy(index(p0));
        iBounds.expandBy(index(p1));
        iBounds.expandBy(index(p2));
        Vector3f e0 = p1 - p0;
        Vector3f e1 = p2 - p0;
        Vector3f n(e0.cross(e1).normalized());
        for (int z = iBounds.min.z(); z <= iBounds.max.z(); ++z) {
            for (int y = iBounds.min.y(); y <= iBounds.max.y(); ++y) {
                for (int x = iBounds.min.x(); x <= iBounds.max.x(); ++x) {
                    Vector3f c = bounds.min + Vector3f(x + 0.5f, y + 0.5f, z + 0.5f) * cellSize;
                    float d = (c - p0).dot(n);
                    // check if voxel intersects triangle plane
                    if (std::abs(d) < 0.70710678118655f * cellSize * 2) {
                        // check if projected voxel center lies within triangle
                        Vector3f p = c - d * n;
                        float invTwoArea = 1.f / e0.cross(e1).norm();
                        float s = invTwoArea * (p0.y()*p2.x() - p0.x() * p2.y() + (p2.y() - p0.y())*p.x() + (p0.x() - p2.x()) * p.y());
                        float t = invTwoArea * (p0.x()*p1.y() - p0.y() * p1.x() + (p0.y() - p1.y())*p.x() + (p1.x() - p0.x()) * p.y());
                        // account for counter-clockwise triangles
                        if (s < 0.f && t < 0.f) {
                            s = -s;
                            t = -t;
                        }
                        if (s >= 0.f && t >= 0.f && (1.f - s - t) >= 0.f) {
                            grid.setValue(x, y, z, true);
                        }
                    }
                }
            }
        }
    };

    for (int i = 0; i < mesh.triangles().cols(); ++i) {
        const Vector3u &triangle = mesh.triangles().col(i);
        const Vector3f &p0 = mesh.vertices().col(triangle[0]);
        const Vector3f &p1 = mesh.vertices().col(triangle[1]);
        const Vector3f &p2 = mesh.vertices().col(triangle[2]);
        rasterizeTriangle(p0, p1, p2);
    }
}

void Voxelizer::voxelize(const Mesh &mesh, float cellSize, std::vector<Vector3f> &positions) {
    Result result;
    voxelize(mesh, cellSize, result);
    auto size = result.grid.size();
    for (int z = 0; z < size.z(); ++z) {
        for (int y = 0; y < size.y(); ++y) {
            for (int x = 0; x < size.x(); ++x) {
                if (result.grid.value(x, y, z)) {
                    positions.emplace_back(result.bounds.min + Vector3f(x + 0.5f, y + 0.5f, z + 0.5f) * cellSize);
                }
            }
        }
    }

}

} // namespace pbs
