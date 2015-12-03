#pragma once

// [1] Weakly compressible SPH for free surface flows
// [2] Predictive-Corrective Incompressible SPH
// [3] Versatile Surface Tension and Adhesion for SPH Fluids

#include "Scene.h"
#include "Grid.h"
#include "Kernel.h"

#include "core/Common.h"
#include "core/Vector.h"
#include "core/Box.h"
#include "core/AlignedAllocator.h"
#include "core/Timer.h"
#include "core/Profiler.h"

#include "geometry/Mesh.h"
#include "geometry/ObjReader.h"
#include "geometry/Voxelizer.h"
#include "geometry/VoxelGrid.h"
#include "geometry/ParticleGenerator.h"

#include <tbb/tbb.h>

#include <vector>

namespace pbs {

class SPH {
public:
    // Simulation parameters
    struct Parameters {
        float particleRadius;
        float particleDiameter;
        float kernelRadius;
        int kernelSupportParticles;
        float particleMass;
        float restDensity;
    };

    SPH(const Scene &scene);

    void reset();
    void update(float dt);

    const Parameters &parameters() const { return _parameters; }
    const Box3f &bounds() const { return _bounds; }
    float maxTimestep() const { return _maxTimestep; }

    const std::vector<Vector3f> &fluidPositions() const { return _fluidPositions; }
    const std::vector<Vector3f> &boundaryPositions() const { return _boundaryPositions; }
    const std::vector<Vector3f> &boundaryNormals() const { return _boundaryNormals; }

private:
    // iterate i=0..count-1 calling func(i)
    template<typename Func>
    inline void iterate(size_t count, Func func) {
#if USE_TBB
        tbb::parallel_for(0ul, count, 1ul, [func] (size_t i) { func(i); });
#else
        for (size_t i = 0; i < count; ++i) { func(i); }
#endif
    }

    // iterate over all neighbours around p, calling func(j, r, r2)
    template<typename Func>
    inline void iterateNeighbours(const Grid &grid, const std::vector<Vector3f> &positions, const Vector3f &p, Func func) {
        grid.lookup(p, _kernelRadius, [&] (size_t j) {
            Vector3f r = p - positions[j];
            float r2 = r.squaredNorm();
            if (r2 < _kernelRadius2) {
                func(j, r, r2);
            }
            return true;
        });
    }

    // returns true if there are neighbours around p
    inline bool hasNeighbours(const Grid &grid, const std::vector<Vector3f> &positions, const Vector3f &p) {
        bool result = false;
        grid.lookup(p, _kernelRadius, [&] (size_t j) {
            if ((p - positions[j]).squaredNorm() < _kernelRadius2) {
                result = true;
                return false;
            } else {
                return true;
            }
        });
        return result;
    }

    void activateBoundaryParticles();
    void computeDensitiesAndPressures();
    void computeDensity();
    void computeNormals();
    void computeForces();
    void computeCollisions(std::function<void(size_t i, const Vector3f &n, float d)> handler);

    void buildScene(const Scene &scene);
    void addFluidParticles(const ParticleGenerator::Volume &volume);
    void addBoundaryParticles(const ParticleGenerator::Boundary &boundary);

    float _particleRadius = 0.01f;
    float _particleRadius2;
    float _particleDiameter;
    float _kernelRadius;
    float _kernelRadius2;
    int _kernelSupportParticles;
    float _restDensity = 1000.f;            ///< Rest density in kg/m^3
    float _particleMass;                    ///< Particle mass
    float _particleMass2;                   ///< Squared particle mass
    float _maxTimestep;                     ///< Maximum allowed timestep

    Parameters _parameters;

    Vector3f _gravity = Vector3f(0.f, -9.81f, 0.f);

    struct {
        const float gamma = 7.f;
        float cs = 10.f;
        float B;
        float viscosity = 0.005f;
        float dt;
    } wcsph;

    Kernel _kernel;

    Box3f _bounds;

    // Fluid particle buffers
    std::vector<Vector3f> _fluidPositions;
    std::vector<Vector3f> _fluidVelocities;
    std::vector<Vector3f> _fluidNormals;
    std::vector<Vector3f> _fluidForces;
    std::vector<float> _fluidDensities;
    std::vector<float> _fluidPressures;
    Grid _fluidGrid;

    // Boundary particle buffers
    std::vector<Vector3f> _boundaryPositions;
    std::vector<Vector3f> _boundaryNormals;
    std::vector<float> _boundaryDensities;
    std::vector<float> _boundaryPressures;
    std::vector<int> _boundaryActive;
    Grid _boundaryGrid;

    float _t = 0.f;
};

} // namespace pbs
