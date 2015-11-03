#pragma once

#include "core/Common.h"

namespace pbs {

// Simple 3D voxel grid.
template<typename T>
class VoxelGrid {
public:
    VoxelGrid() : _size(0) {}

    VoxelGrid(const Vector3i &size) {
        resize(size);
    }

    void resize(const Vector3i &size) {
        _size = size;
        _x = _size.x();
        _xy = _size.x() * _size.y();
        _voxels.resize(_size.prod());
    }

    void resize(int x, int y, int z) {
        resize(Vector3i(x, y, z));
    }

    const Vector3i &size() const { return _size; }

    const T &operator()(const Vector3i &index) const { return _voxels[linearize(index)]; }
          T &operator()(const Vector3i &index)       { return _voxels[linearize(index)]; }

    const T &operator()(int x, int y, int z) const { return _voxels[linearize(Vector3i(x, y, z))]; }
          T &operator()(int x, int y, int z)       { return _voxels[linearize(Vector3i(x, y, z))]; }

    const T *data() const { return _voxels.data(); }

private:
    inline size_t linearize(const Vector3i &index) const {
        return index.z() * _xy + index.y() * _x + index.x();
    }

    Vector3i _size;
    int _x;
    int _xy;
    std::vector<T> _voxels;
};

typedef VoxelGrid<float> VoxelGridf;

} // namespace pbs
