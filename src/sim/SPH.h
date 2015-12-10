#pragma once

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
    void updateStep();

    const Parameters &parameters() const { return _parameters; }
    const Box3f &bounds() const { return _bounds; }
    float timeStep() const { return _timeStep; }
    float time() const { return _time; }

    const std::vector<Vector3f> &fluidPositions() const { return _fluidPositions; }
          std::vector<Vector3f> &fluidPositions()       { return _fluidPositions; }
    const std::vector<Vector3f> &boundaryPositions() const { return _boundaryPositions; }
    const std::vector<Vector3f> &boundaryNormals() const { return _boundaryNormals; }
    const std::vector<Mesh> &boundaryMeshes() const { return _boundaryMeshes; }

private:

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

    // iterate over all neighbours around p, calling func(j, r, r2)
    template<typename Func>
    inline void iterateNeighbours2(const Grid &grid, const std::vector<Vector3f> &positionsNew, const Vector3f &p, const Vector3f &pNew, Func func) {
        grid.lookup(p, _kernelRadius, [&] (size_t j) {
            Vector3f r = pNew - positionsNew[j];
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

    // Shared update methods
    void activateBoundaryParticles();
    void updateBoundaryGrid();
    void updateBoundaryMasses();
    void updateDensities();
    void updateNormals();
    void computeCollisions(std::function<void(size_t i, const Vector3f &n, float d)> handler);
    void enforceBounds();

    // WCSPH update methods
    void wcsphUpdateDensitiesAndPressures();
    void wcsphUpdateForces();

    void wcsphInit();
    void wcsphUpdate();

    // PCISPH update methods
    void pcisphUpdateGrid();
    void pcisphUpdateDensityVariationScaling();
    void pcisphInitializeForces();
    void pcisphPredictVelocitiesAndPositions();
    void pcisphUpdatePressures();
    void pcisphUpdatePressureForces();
    void pcisphUpdateVelocitiesAndPositions();

    void pcisphInit();
    void pcisphUpdate(int maxIterations = 100);

    void buildScene(const Scene &scene);
    void addFluidParticles(const ParticleGenerator::Volume &volume);
    void addBoundaryParticles(const ParticleGenerator::Boundary &boundary);

    enum Method {
        WCSPH,
        PCISPH,
    };

    static std::string methodToString(Method method);
    static Method stringToMethod(const std::string &str);

    Method _method;
    float _particleRadius = 0.01f;
    float _particleRadius2;
    float _particleDiameter;
    int _kernelRadiusFactor = 4;
    float _kernelRadius;
    float _kernelRadius2;
    int _kernelSupportParticles;
    float _restDensity = 1000.f;            ///< Rest density in kg/m^3
    float _surfaceTension = 1.f;            ///< Surface tension amount
    float _viscosity = 0.f;                 ///< Viscosity
    float _timeStep = 0.001f;
    float _compressionThreshold = 0.02f;

    float _particleMass;                    ///< Particle mass
    float _particleMass2;                   ///< Squared particle mass
    float _invParticleMass;                 ///< Inverse particle mass

    float _maxDensityVariationThreshold;
    float _avgDensityVariationThreshold;
    float _densityVariationScaling;
    float _maxDensityVariation;
    float _prevMaxDensityVariation = 1000.f;
    float _avgDensityVariation;
    float _maxVelocity;
    float _maxForce;

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
    std::vector<Vector3f> _fluidPositionsNew;
    std::vector<Vector3f> _fluidVelocitiesNew;
    std::vector<Vector3f> _fluidPositionsPreShock;
    std::vector<Vector3f> _fluidVelocitiesPreShock;
    std::vector<Vector3f> _fluidNormals;
    std::vector<Vector3f> _fluidForces;
    std::vector<Vector3f> _fluidPressureForces;
    std::vector<float> _fluidDensities;
    std::vector<float> _fluidPressures;
    Grid _fluidGrid;

    // Boundary particle buffers
    std::vector<Vector3f> _boundaryPositions;
    std::vector<Vector3f> _boundaryNormals;
    std::vector<float> _boundaryDensities;
    std::vector<float> _boundaryPressures;
    std::vector<float> _boundaryMasses;
    std::vector<int> _boundaryActive;
    Grid _boundaryGrid;

    std::vector<Mesh> _boundaryMeshes;

    float _time = 0.f;
    float _timePreShock;
};

} // namespace pbs
