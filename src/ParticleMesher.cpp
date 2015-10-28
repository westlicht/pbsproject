#include "ParticleMesher.h"
#include "Mesh.h"
#include "VoxelGrid.h"
#include "MarchingCubes.h"
#include "Box.h"
#include "Vector.h"
#include "Timer.h"

#include <tbb/tbb.h>

#include <atomic>

namespace pbs {

#define USE_TBB 1

// Uniform grid for indexing particles.
// To allow for particles with different footprints (anisotropic kernels)
// particles may be added to multiple cells. Lookups are fast as only one
// cell needs to be visited.
class Grid {
public:
    template<typename Func>
    Grid(const MatrixXf &positions, const Box3f &bounds, const Vector3i &cells, Func func) :
        _bounds(bounds),
        _cells(cells),
        _cellSize(_bounds.extents().cwiseProduct(Vector3f(_cells.x(), _cells.y(), _cells.z()).cwiseInverse())),
        _invCellSize(_cellSize.cwiseInverse())
    {
        // Compute particle bounds
        std::vector<Box3i> particleBounds(positions.cols());
#if USE_TBB
        tbb::parallel_for(0, int(positions.cols()), 1, [this, func, &positions, &particleBounds] (int i) {
            Box3f box = func(positions.col(i));
            particleBounds[i] = Box3i(index(box.min).cwiseMax(0), index(box.max).cwiseMin(_cells - Vector3i(1)));
        });
#else
        for (int i = 0; i < positions.cols(); ++i) {
            Box3f box = func(positions.col(i));
            particleBounds[i] = Box3i(index(box.min).cwiseMax(0), index(box.max).cwiseMin(_cells - Vector3i(1)));
        }
#endif


#if USE_TBB
        std::atomic<size_t> *cellCount = new std::atomic<size_t>[_cells.prod()];
        std::atomic<size_t> *cellIndex = new std::atomic<size_t>[_cells.prod()];
#else
        std::vector<size_t> cellCount(_cells.prod(), 0);
        std::vector<size_t> cellIndex(_cells.prod(), 0);
#endif
        _cellOffset.resize(_cells.prod() + 1);

        // Count number of particles per cell
#if USE_TBB
        tbb::parallel_for(0, int(particleBounds.size()), 1, [this, &particleBounds, &cellCount] (int i) {
            iterate(particleBounds[i], [&] (const Vector3i &index) {
                ++cellCount[linearize(index)];
            });
        });
#else
        for (const auto &box : particleBounds) {
            iterate(box, [&] (const Vector3i &index) {
                ++cellCount[linearize(index)];
            });
        }
#endif

        // Initialize cell indices & offsets
        size_t index = 0;
        for (size_t i = 0; i < size_t(_cells.prod()); ++i) {
            cellIndex[i] = index;
            _cellOffset[i] = index;
            index += cellCount[i];
        }
        _cellOffset.back() = index;

        // Put particles into cells
        _indices.resize(_cellOffset.back());
#if USE_TBB
        tbb::parallel_for(0, int(particleBounds.size()), 1, [this, &particleBounds, &cellIndex] (int i) {
            iterate(particleBounds[i], [&] (const Vector3i &index) {
                _indices[cellIndex[linearize(index)]++] = i;
            });
        });
#else
        for (size_t i = 0; i < particleBounds.size(); ++i) {
            iterate(particleBounds[i], [&] (const Vector3i &index) {
                _indices[cellIndex[linearize(index)]++] = i;
            });
        }
#endif

#if USE_TBB
        delete [] cellCount;
        delete [] cellIndex;
#endif
    }

    template<typename Func>
    void lookup(const Vector3f &pos, Func func) {
        size_t i = linearize(index(pos).cwiseMax(0).cwiseMin(_cells - Vector3i(1)));
        for (size_t j = _cellOffset[i]; j < _cellOffset[i + 1]; ++j) {
            func(_indices[j]);
        }
    }

private:
    inline size_t linearize(const Vector3i &index) {
        return index.z() * (_cells.x() * _cells.y()) + index.y() * _cells.x() + index.x();
    }

    inline Vector3i index(const Vector3f &pos) {
        return Vector3i(
            int(std::floor((pos.x() - _bounds.min.x()) * _invCellSize.x())),
            int(std::floor((pos.y() - _bounds.min.y()) * _invCellSize.y())),
            int(std::floor((pos.z() - _bounds.min.z()) * _invCellSize.z()))
        );
    }

    inline size_t indexLinear(const Vector3f &pos) {
        return linearize(index(pos));
    }

    template<typename Func>
    inline void iterate(const Box3i &range, Func func) {
        for (int z = range.min.z(); z <= range.max.z(); ++z) {
            for (int y = range.min.y(); y <= range.max.y(); ++y) {
                for (int x = range.min.x(); x <= range.max.x(); ++x) {
                    func(Vector3i(x, y, z));
                }
            }
        }
    }


    Box3f _bounds;
    Vector3i _cells;
    Vector3f _cellSize;
    Vector3f _invCellSize;

    std::vector<size_t> _cellOffset;
    std::vector<size_t> _indices;
};

Mesh ParticleMesher::createMeshIsotropic(const MatrixXf &positions, const Box3f &bounds, const Vector3i &cells, float smoothRadius, float normalization, float isoLevel) {

    DBG("Generating mesh from particles (isotropic kernel)");

    DBG("Building acceleration grid ...");
    Timer timer;
    Grid grid(positions, bounds, cells, [smoothRadius] (const Vector3f &p) { return Box3f(p - Vector3f(smoothRadius), p + Vector3f(smoothRadius)); });
    DBG("Took %s", timer.elapsedString());

    DBG("Building voxel grid ...");
    timer.reset();
    VoxelGridf voxelGrid(cells + Vector3i(1));
    Vector3f min = bounds.min;
    Vector3f extents = bounds.extents();

    //smoothRadius *= 0.5;
    float poly6Constant = 365.f / (64.f * M_PI * std::pow(smoothRadius, 9.f));
    normalization *= poly6Constant;
    float h2 = sqr(smoothRadius);
#if USE_TBB
    tbb::parallel_for(0, (cells + Vector3i(1)).prod(), 1, [cells, min, extents, poly6Constant, h2, normalization, &grid, &positions, &voxelGrid] (int i) {
        int x = i % (cells.x() + 1);
        int y = (i / (cells.x() + 1)) % (cells.y() + 1);
        int z = (i / ((cells.x() + 1) * (cells.y() + 1))) % (cells.z() + 1);
#else
    for (int z = 0; z <= cells.z(); ++z) {
        for (int y = 0; y <= cells.y(); ++y) {
            for (int x = 0; x <= cells.x(); ++x) {
#endif
                Vector3f p = Vector3f(float(x) / cells.x(), float(y) / cells.y(), float(z) / cells.z()).cwiseProduct(extents) + min;
                float c = 0.f;
                grid.lookup(p, [&] (size_t i) {
                    float r2 = (p - positions.col(i)).squaredNorm();
                    if (r2 < h2 && r2 != 0.f) {
                        c += cube(h2 - r2);
                    }
                });
                c *= normalization;
                voxelGrid(x, y, z) = c;
                //DBG("c = %f", c);
#if USE_TBB
    });
#else
            }
        }
    }
#endif

    DBG("Took %s", timer.elapsedString());

    DBG("Building surface ...");
    timer.reset();
    MarchingCubes<float> mc;
    mc.GenerateSurface(voxelGrid.data(), isoLevel, cells.x(), cells.y(), cells.z(), 1.f / cells.x(), 1.f / cells.y(), 1.f / cells.z());
    DBG("Took %s", timer.elapsedString());

    Mesh mesh = mc.mesh();
    return std::move(mesh);
}

} // namespace pbs
