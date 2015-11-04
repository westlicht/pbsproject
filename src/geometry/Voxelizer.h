#pragma once

#include "VoxelGrid.h"

#include "core/Common.h"
#include "core/Box.h"

namespace pbs {

class Mesh;

class Voxelizer {
public:
    struct Result {
        float cellSize;
        Box3f bounds;
        VoxelGrid<bool> grid;
    };

    static void voxelize(const Mesh &mesh, float cellSize, Result &result);
    static void voxelize(const Mesh &mesh, float cellSize, std::vector<Vector3f> &positions);

};

} // namespace pbs
