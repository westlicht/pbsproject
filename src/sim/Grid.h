#pragma once

#include "core/Common.h"
#include "core/Vector.h"
#include "core/Box.h"
#include "core/Morton.h"

#include <vector>

namespace pbs {

class Grid {
public:
    void init(const Box3f &bounds, float cellSize) {
        _bounds = bounds;
        _cellSize = cellSize;
        _invCellSize = 1.f / cellSize;

        _size = Vector3i(
            nextPowerOfTwo(int(std::floor(_bounds.extents().x() / _cellSize)) + 1),
            nextPowerOfTwo(int(std::floor(_bounds.extents().y() / _cellSize)) + 1),
            nextPowerOfTwo(int(std::floor(_bounds.extents().z() / _cellSize)) + 1)
        );

        _cellOffset.resize(_size.prod() + 1);

        DBG("Initialized grid: bounds = %s, cellSize = %f, size = %s", _bounds, _cellSize, _size);
    }

    inline Vector3i index(const Vector3f &pos) {
        return Vector3i(
            int(std::floor((pos.x() - _bounds.min.x()) * _invCellSize)),
            int(std::floor((pos.y() - _bounds.min.y()) * _invCellSize)),
            int(std::floor((pos.z() - _bounds.min.z()) * _invCellSize))
        );
    }

    inline size_t indexLinear(const Vector3f &pos) {
        Vector3i i = index(pos);
        return i.z() * (_size.x() * _size.y()) + i.y() * _size.x() + i.x();
    }

    inline uint32_t indexMorton(const Vector3i &index) {
        return Morton3D::morton10bit(index.x(), index.y(), index.z());
    }

    inline uint32_t indexMorton(const Vector3f &pos) {
        return indexMorton(index(pos));
    }

    template<typename SwapFunc>
    void update(const std::vector<Vector3f> &positions, SwapFunc swap) {
        std::vector<uint32_t> cellCount(_size.prod(), 0);
        std::vector<uint32_t> cellIndex(_size.prod(), 0);

        size_t count = positions.size();

        std::vector<uint32_t> indices(count);

        // Update particle index and count number of particles per cell
        for (size_t i = 0; i < count; ++i) {
            uint32_t index = indexLinear(positions[i]);
            //uint32_t index = indexMorton(particles[i].p);
            indices[i] = index;
            ++cellCount[index];
        }

        // Initialize cell indices & offsets
        size_t index = 0;
        for (size_t i = 0; i < cellIndex.size(); ++i) {
            cellIndex[i] = index;
            _cellOffset[i] = index;
            index += cellCount[i];
        }
        _cellOffset.back() = index;

        // Sort particles by index
        for (size_t i = 0; i < count; ++i) {
            while (i >= cellIndex[indices[i]] || i < _cellOffset[indices[i]]) {
                size_t j = cellIndex[indices[i]]++;
                std::swap(indices[i], indices[j]);
                swap(i, j);
            }
        }
    }

    template<typename Func>
    void lookup(const Vector3f &pos, float radius, Func func) {
        Vector3i min = index(pos - Vector3f(radius)).cwiseMax(Vector3i(0));
        Vector3i max = index(pos + Vector3f(radius)).cwiseMin(_size - Vector3i(1));
        for (int z = min.z(); z <= max.z(); ++z) {
            for (int y = min.y(); y <= max.y(); ++y) {
                for (int x = min.x(); x <= max.x(); ++x) {
                    size_t i = z * (_size.x() * _size.y()) + y * _size.x() + x;
                    for (size_t j = _cellOffset[i]; j < _cellOffset[i + 1]; ++j) {
                        func(j);
                    }
                }
            }
        }
    }

private:
    Box3f _bounds;
    float _cellSize;
    float _invCellSize;

    Vector3i _size;
    std::vector<size_t> _cellOffset;
};

} // namespace pbs
