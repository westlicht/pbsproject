#pragma once

#include "core/Common.h"
#include "core/Vector.h"

#include <vector>

namespace pbs {

// Simple 3D voxel grid.
template<typename T>
class VoxelGrid {
public:
    VoxelGrid() {}

    VoxelGrid(const Vector3i &size) {
        resize(size);
    }

    VoxelGrid(const Vector3i &size, const T &value) {
        resize(size);
        fill(value);
    }

    // resize the voxel grid
    void resize(const Vector3i &size) {
        _size = size;
        _x = _size.x();
        _xy = _size.x() * _size.y();
        _voxels.resize(_size.prod());
    }

    // resize the voxel grid
    void resize(int x, int y, int z) {
        resize(Vector3i(x, y, z));
    }

    // fills the grid with the given value
    void fill(const T &value) {
        std::fill(_voxels.begin(), _voxels.end(), value);
    }

    // specifies the size of the voxel grid (number of voxels)
    const Vector3i &size() const { return _size; }

    // specifies the origin in world space
    const Vector3f &origin() const { return _origin; }
    void setOrigin(const Vector3f &origin) { _origin = origin; }

    // specifies the cell size in world space
    float cellSize() const { return _cellSize; }
    void setCellSize(float cellSize) { _cellSize = cellSize; }

    // transforms a point in world space to voxel space
    inline Vector3f toVoxelSpace(const Vector3f &vsP) const {
        return (vsP - _origin) * (1.f / _cellSize);
    }

    // transforms a point in voxel space to world space
    inline Vector3f toWorldSpace(const Vector3f &wsP) const {
        return _origin + wsP * _cellSize;
    }

    // voxel data accessor
    const T &operator()(const Vector3i &index) const { return _voxels[linearize(index)]; }
          T &operator()(const Vector3i &index)       { return _voxels[linearize(index)]; }

    const T &operator()(int x, int y, int z) const { return _voxels[linearize(Vector3i(x, y, z))]; }
          T &operator()(int x, int y, int z)       { return _voxels[linearize(Vector3i(x, y, z))]; }

    T value(const Vector3i &index) const { return _voxels[linearize(index)]; }
    T value(int x, int y, int z) const { return _voxels[linearize(Vector3i(x, y, z))]; }
    void setValue(const Vector3i &index, const T &value) { _voxels[linearize(index)] = value; }
    void setValue(int x, int y, int z, const T &value) { _voxels[linearize(Vector3i(x, y, z))] = value; }

    // trilinear filtering
    T trilinear(const Vector3f &vsP) const {
        #if 0
        Vector3f uvw(
            clamp((p.x() - _origin.x()) / _cellSize - 0.5f, 0.f, float(_size.x() - 1)),
            clamp((p.y() - _origin.y()) / _cellSize - 0.5f, 0.f, float(_size.y() - 1)),
            clamp((p.z() - _origin.z()) / _cellSize - 0.5f, 0.f, float(_size.z() - 1))
        );
        #endif

        Vector3f uvw(vsP - Vector3f(0.5f));

        int i0 = std::max(0, int(std::floor(uvw.x())));
        int j0 = std::max(0, int(std::floor(uvw.y())));
        int k0 = std::max(0, int(std::floor(uvw.z())));
        int i1 = std::min(int(_size.x() - 1), i0 + 1);
        int j1 = std::min(int(_size.y() - 1), j0 + 1);
        int k1 = std::min(int(_size.z() - 1), k0 + 1);
        uvw -= Vector3f(float(i0), float(j0), float(k0));

        T temp1, temp2;

        temp1 = (*this)(i0,j0,k0) + T(((*this)(i0,j0,k1) - (*this)(i0,j0,k0)) * uvw.z());
        temp2 = (*this)(i0,j1,k0) + T(((*this)(i0,j1,k1) - (*this)(i0,j1,k0)) * uvw.z());
        T result1 = temp1 + T((temp2-temp1) * uvw.y());

        temp1 = (*this)(i1,j0,k0) + T(((*this)(i1,j0,k1) - (*this)(i1,j0,k0)) * uvw.z());
        temp2 = (*this)(i1,j1,k0) + T(((*this)(i1,j1,k1) - (*this)(i1,j1,k0)) * uvw.z());
        T result2 = temp1 + T((temp2 - temp1) * uvw.y());

        return result1 + T(uvw.x() * (result2 - result1));
    }

    // returns the gradient at the given position using central differences
    TVector<T, 3> gradient(const Vector3f &vsP, float eps = 1e-5f) const {
        return TVector<T, 3>(
            trilinear(vsP + Vector3f(eps, 0.f, 0.f)) - trilinear(vsP - Vector3f(eps, 0.f, 0.f)),
            trilinear(vsP + Vector3f(0.f, eps, 0.f)) - trilinear(vsP - Vector3f(0.f, eps, 0.f)),
            trilinear(vsP + Vector3f(0.f, 0.f, eps)) - trilinear(vsP - Vector3f(0.f, 0.f, eps))
        ) * (0.5f / eps);
    }

    // raw data
    const T *data() const { return _voxels.data(); }

private:
    inline size_t linearize(const Vector3i &index) const {
        return index.z() * _xy + index.y() * _x + index.x();
    }

    Vector3i _size;
    int _x;
    int _xy;
    Vector3f _origin;
    float _cellSize = 1.f;
    std::vector<T> _voxels;

};

typedef VoxelGrid<float> VoxelGridf;
typedef VoxelGrid<bool> VoxelGridb;

} // namespace pbs
