#include "ParticleMesher.h"
#include "Mesh.h"
#include "VoxelGrid.h"
#include "MarchingCubes.h"

#include "core/Common.h"
#include "core/Box.h"
#include "core/Vector.h"
#include "core/Timer.h"

#include <tbb/tbb.h>

#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <atomic>

namespace pbs {

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
        for (size_t i = 0; i < size_t(_cells.prod()); ++i) {
            cellCount[i] = 0;
            cellIndex[i] = 0;
        }
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

Mesh ParticleMesher::createMeshIsotropic(const MatrixXf &positions, const Box3f &bounds, const Vector3i &cells, const Parameters &params) {
    float kernelRadius = params.kernelRadius;
    float kernelRadius2 = sqr(kernelRadius);

    DBG("Generating mesh from particles (isotropic kernel)");
    Timer timer;
    Timer timerTotal;

    Vector3i gridCells(cells / 4);

    DBG("Building acceleration grid (resolution=%s) ...", gridCells);
    timer.reset();
    Grid grid(positions, bounds, gridCells, [kernelRadius] (const Vector3f &p) { return Box3f(p - Vector3f(kernelRadius), p + Vector3f(kernelRadius)); });
    DBG("Took %s", timer.elapsedString());

    DBG("Building voxel grid (resolution=%s) ...", cells);
    timer.reset();
    VoxelGridf voxelGrid(cells + Vector3i(1));
    Vector3f min = bounds.min;
    Vector3f extents = bounds.extents();

    float poly6Constant = 365.f / (64.f * M_PI * std::pow(kernelRadius, 9.f));
    float normalization = params.particleMass / params.restDensity * poly6Constant;
#if USE_TBB
    tbb::parallel_for(0, (cells + Vector3i(1)).prod(), 1, [cells, min, extents, poly6Constant, kernelRadius2, normalization, &grid, &positions, &voxelGrid] (int i) {
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
                    if (r2 < kernelRadius2 && r2 != 0.f) {
                        c += cube(kernelRadius2 - r2);
                    }
                });
                c *= normalization;
                voxelGrid(x, y, z) = c;
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
    Mesh mesh = mc.generateIsoSurface(voxelGrid.data(), params.isoLevel, bounds, cells);
    DBG("Took %s", timer.elapsedString());

    DBG("Building mesh took %s", timerTotal.elapsedString());

    return std::move(mesh);
}

Mesh ParticleMesher::createMeshAnisotropic(MatrixXf &positions, const Box3f &bounds, const Vector3i &cells, const Parameters &params) {
    typedef Eigen::Matrix<float, 3, 3> Matrix3f;
    typedef Eigen::DiagonalMatrix<float, 3> DiagonalMatrix3f;

    DBG("Generating mesh from particles (anisotropic kernel)");
    Timer timer;

#if 1
    // Smooth particle positions
    {
        float kernelRadius = params.kernelRadius;
        float kernelRadius2 = sqr(kernelRadius);
        float lambda = 0.9f; // smoothing strength

        DBG("Building acceleration grid (smoothing) ...");
        timer.reset();
        Grid grid(positions, bounds, cells, [kernelRadius] (const Vector3f &p) { return Box3f(p - Vector3f(kernelRadius), p + Vector3f(kernelRadius)); });
        DBG("Took %s", timer.elapsedString());

        DBG("Smoothing particle positions ...");
        timer.reset();
#if USE_TBB
        tbb::parallel_for(0ul, size_t(positions.cols()), 1ul, [&] (size_t i) {
#else
        for (size_t i = 0; i < size_t(positions.cols()); ++i) {
#endif
            Vector3f x = positions.col(i);
            Vector3f xh;
            float totalWeight = 0.f;
            grid.lookup(positions.col(i), [&] (size_t j) {
                Vector3f r = positions.col(j) - x;
                float r2 = r.squaredNorm();
                if (r2 < kernelRadius2) {
                    float weight = cube(kernelRadius2 - r2);
                    xh += weight * positions.col(j);
                    totalWeight += weight;
                }
            });
            xh *= (1.f / totalWeight);
            positions.col(i) = lerp(lambda, x, xh);
#if USE_TBB
        });
#else
        }
#endif
        DBG("Took %s", timer.elapsedString());
    }

    std::vector<Matrix3f> kernels(positions.cols());
    std::vector<float> determinants(positions.cols());

    // Compute anisotropic kernels
    {
        auto kernel = [] (const Vector3f &r, float kernelRadius) {
            float t = r.norm() / kernelRadius;
            return 1.f - cube(t);
        };

        float kernelRadius = params.kernelRadius;
        float kernelRadius2 = sqr(kernelRadius);
        int Ns = params.kernelSupportParticles;
        int Nsurface = 0;

        DBG("Building acceleration grid (kernel estimation) ...");
        timer.reset();
        Grid grid(positions, bounds, cells, [kernelRadius] (const Vector3f &p) { return Box3f(p - Vector3f(kernelRadius), p + Vector3f(kernelRadius)); });
        DBG("Took %s", timer.elapsedString());

        DBG("Estimating kernels ...");
        timer.reset();
#if USE_TBB
        tbb::parallel_for(0ul, size_t(positions.cols()), 1ul, [&] (size_t i) {
#else
        for (size_t i = 0; i < size_t(positions.cols()); ++i) {
#endif
            // Compute mass center and number of neighbors for early out
            Vector3f x = positions.col(i);
            Vector3f xm;
            int N = 0;
            grid.lookup(positions.col(i), [&] (size_t j) {
                Vector3f r = positions.col(j) - x;
                float r2 = r.squaredNorm();
                if (r2 < kernelRadius2) {
                    xm += x;
                    ++N;
                }
            });
            xm *= (1.f / N);

            // Classify particle as surface or interior
            if (float(N - Ns) / Ns > 0.1f || (x - xm).squaredNorm() / sqr(2.f * params.particleRadius) > 0.1f) {
                // Surface particle -> compute kernel
                ++Nsurface;
                // Compute weighted mean
                Vector3f xw;
                float totalWeight;
                grid.lookup(positions.col(i), [&] (size_t j) {
                    Vector3f r = positions.col(j) - x;
                    if (r.squaredNorm() < kernelRadius2) {
                        float weight = kernel(r, kernelRadius);
                        xw += weight * x;
                        totalWeight += weight;
                    }
                });
                xw *= (1.f / totalWeight);

                // Compute weighted covariance matrix
                Matrix3f C;
                grid.lookup(positions.col(i), [&] (size_t j) {
                    Vector3f r = positions.col(j) - x;
                    if (r.squaredNorm() < kernelRadius2) {
                        float weight = kernel(r, kernelRadius);
                        Vector3f d(positions.col(j) - xw);
                        C += weight * d * d.transpose();
                    }
                });
                C.array() *= (1.f / totalWeight);

                // Compute SVD
                float kr = 4.f;
                float ks = 1400.f;
                //float kn = 0.5f;

                Eigen::JacobiSVD<Matrix3f> svd(C, Eigen::ComputeFullU);
                const auto &R = svd.matrixU();

                Vector3f sigma = svd.singularValues();
                sigma.y() = std::max(sigma.y(), sigma.x() / kr);
                sigma.z() = std::max(sigma.z(), sigma.x() / kr);
                sigma *= ks;

                //DBG("sigma = %s", sigma);
                DiagonalMatrix3f invSigma(sigma.cwiseInverse());
                //DBG("invSigma = %s", Matrix3f(invSigma));
                Matrix3f G(R * invSigma * R.transpose().eval());
                G *= (1.f / kernelRadius);

                //DBG("G = %s", G);

                G = DiagonalMatrix3f(Vector3f(1.f / kernelRadius, 1.f / kernelRadius, 0.01 / kernelRadius));

                kernels[i] = G;
                determinants[i] = G.determinant();

            } else {
                // Interior particle -> use isotropic kernel
                kernels[i] = Matrix3f::Identity() * (1.f / kernelRadius);
                determinants[i] = cube(1.f / kernelRadius);
                //DBG("Gs = %s", kernels[i]);
            }
#if USE_TBB
        });
#else
        }
#endif
        DBG("Took %s", timer.elapsedString());
        DBG("%d of %d (%.1f%%) particles are classified surface particles", Nsurface, positions.cols(), float(Nsurface) / positions.cols() * 100.f);
    }
#endif

    float kernelRadius = params.kernelRadius;
    float kernelRadius2 = sqr(kernelRadius);

    DBG("Building acceleration grid ...");
    timer.reset();
    Grid grid(positions, bounds, cells, [kernelRadius] (const Vector3f &p) { return Box3f(p - Vector3f(kernelRadius), p + Vector3f(kernelRadius)); });
    DBG("Took %s", timer.elapsedString());


    DBG("Building voxel grid ...");
    timer.reset();
    VoxelGridf voxelGrid(cells + Vector3i(1));
    Vector3f min = bounds.min;
    Vector3f extents = bounds.extents();

    //float poly6Constant = 365.f / (64.f * M_PI * std::pow(kernelRadius, 9.f));
    float poly6Constant = 365.f / (64.f * M_PI /* * std::pow(kernelRadius, 3.f)*/);
    float normalization = params.particleMass / params.restDensity * poly6Constant;

    //DiagonalMatrix3f G(Vector3f(1.f / kernelRadius));
    //float det = G.diagonal().prod();


#if USE_TBB
    tbb::parallel_for(0, (cells + Vector3i(1)).prod(), 1, [cells, min, extents, poly6Constant, kernelRadius2, normalization, &grid, &positions, &kernels, &determinants, &voxelGrid] (int i) {
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
                    Vector3f r = p - positions.col(i);
                    float r2 = r.squaredNorm();
                    if (r2 < kernelRadius2 && r2 != 0.f) {
                        const auto &G = kernels[i];
                        const auto &det = determinants[i];
                        c += cube(1.f - sqr((G * r).norm())) * det;
                    }
                });
                c *= normalization;
                voxelGrid(x, y, z) = c;
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
    Mesh mesh = mc.generateIsoSurface(voxelGrid.data(), params.isoLevel, bounds, cells);
    DBG("Took %s", timer.elapsedString());

    return std::move(mesh);
}

} // namespace pbs
