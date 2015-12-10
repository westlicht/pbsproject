#include "SPH.h"

#include "core/DebugMonitor.h"

#include <tbb/enumerable_thread_specific.h>

// [1] Weakly compressible SPH for free surface flows
// [2] Predictive-Corrective Incompressible SPH
// [3] Versatile Surface Tension and Adhesion for SPH Fluids
// [4] Versatile Rigid-Fluid Coupling for Incompressible SPH

namespace pbs {

#define HANDLE_BOUNDARIES 1


template<typename T>
static void dumpVector(const std::vector<T> &v) {
    for (const auto &i : v) {
        DBG("%s", i);
    }
}

SPH::SPH(const Scene &scene) {
    // Load scene settings
    _method = stringToMethod(scene.settings.getString("method", "pcisph"));
    _particleRadius = scene.settings.getFloat("particleRadius", _particleRadius);
    _restDensity = scene.settings.getFloat("restDensity", _restDensity);
    _gravity = scene.settings.getVector3("gravity", _gravity);
    _surfaceTension = scene.settings.getFloat("surfaceTension", _surfaceTension);
    _viscosity = scene.settings.getFloat("viscosity", _viscosity);
    _timeStep = scene.settings.getFloat("timeStep", _timeStep);
    _compressionThreshold = scene.settings.getFloat("compressionThreshold", _compressionThreshold);

    // Compute derived constants
    _particleRadius2 = sqr(_particleRadius2);
    _particleDiameter = 2.f * _particleRadius;

    _kernelRadius = _kernelRadiusFactor * _particleRadius;
    _kernelRadius2 = sqr(_kernelRadius);
    _kernelSupportParticles = int(std::ceil((4.f / 3.f * M_PI * cube(_kernelRadius)) / cube(_particleDiameter)));

    //_particleMass = _restDensity / cube(1.f / _particleDiameter);
    _particleMass = _restDensity * cube(_particleDiameter);
    _particleMass /= 1.17f;
    _particleMass2 = sqr(_particleMass);
    _invParticleMass = 1.f / _particleMass;


    _avgDensityVariationThreshold = _compressionThreshold * _restDensity;
    _maxDensityVariationThreshold = _avgDensityVariationThreshold * 10.f;

    wcsph.B = _restDensity * sqr(wcsph.cs) / wcsph.gamma;
    wcsph.dt = std::min(0.25f * _kernelRadius / (_particleMass * 9.81f), 0.4f * _kernelRadius / (wcsph.cs * (1.f + 0.6f * wcsph.viscosity)));

    // Store parameters
    _parameters.particleRadius = _particleRadius;
    _parameters.particleDiameter = _particleDiameter;
    _parameters.kernelRadius = _kernelRadius;
    _parameters.kernelSupportParticles = _kernelSupportParticles;
    _parameters.particleMass = _particleMass;
    _parameters.restDensity = _restDensity;

    buildScene(scene);

    // Compute bounds
    _bounds.reset();
    for (const auto &p : _boundaryPositions) {
        _bounds.expandBy(p);
    }

    _fluidVelocities.resize(_fluidPositions.size());
    _fluidPositionsNew.resize(_fluidPositions.size());
    _fluidVelocitiesNew.resize(_fluidPositions.size());
    _fluidPositionsPreShock.resize(_fluidPositions.size());
    _fluidPositionsPreShock.resize(_fluidPositions.size());
    _fluidNormals.resize(_fluidPositions.size());
    _fluidForces.resize(_fluidPositions.size());
    _fluidPressureForces.resize(_fluidPositions.size());
    _fluidDensities.resize(_fluidPositions.size());
    _fluidPressures.resize(_fluidPositions.size());

    _boundaryDensities.resize(_boundaryPositions.size());
    _boundaryPressures.resize(_boundaryPositions.size());
    _boundaryMasses.resize(_boundaryPositions.size());
    _boundaryActive.resize(_boundaryPositions.size());

    _kernel.init(_kernelRadius);
    _fluidGrid.init(_bounds, _kernelRadius);
    _boundaryGrid.init(_bounds, _kernelRadius);

    // Preprocessing
    updateBoundaryGrid();
    updateBoundaryMasses();

    DBG("method = %s", methodToString(_method));
    DBG("particleRadius = %f", _particleRadius);
    DBG("kernelRadius = %f", _kernelRadius);
    DBG("kernelSupportParticles = %d", _kernelSupportParticles);
    DBG("restDensity = %f", _restDensity);
    DBG("particleMass = %f", _particleMass);
    DBG("gravity = %s", _gravity);
    DBG("surfaceTension = %f", _surfaceTension);
    DBG("viscosity = %f", _viscosity);
    DBG("timeStep = %f", _timeStep);

    DBG("wcsph.gamma = %f", wcsph.gamma);
    DBG("wcsph.cs = %f", wcsph.cs);
    DBG("wcsph.B = %f", wcsph.B);
    DBG("wcsph.viscosity = %f", wcsph.viscosity);
    DBG("wcsph.dt = %f", wcsph.dt);

    DBG("# particles = %d", _fluidPositions.size());
    DBG("# boundary particles = %d", _boundaryPositions.size());

    DBG("Initializing simulation ...");
    switch (_method) {
    case WCSPH: wcsphInit(); break;
    case PCISPH: pcisphInit(); break;
    }
}

void SPH::reset() {

}

void SPH::update(float dt) {
    float targetTime = _time + dt;
    while (_time < targetTime) {
        updateStep();
    }
}

void SPH::updateStep() {
    switch (_method) {
    case WCSPH: wcsphUpdate(); break;
    case PCISPH: pcisphUpdate(); break;
    }
}

// Activate all boundary particles that are nearby fluid particles
void SPH::activateBoundaryParticles() {
    parallelFor(_boundaryPositions.size(), [this] (size_t i) {
        _boundaryActive[i] = hasNeighbours(_fluidGrid, _fluidPositions, _boundaryPositions[i]);
    });
}

void SPH::updateBoundaryGrid() {
    _boundaryGrid.update(_boundaryPositions, [this] (size_t i, size_t j) {
        std::swap(_boundaryPositions[i], _boundaryPositions[j]);
        std::swap(_boundaryNormals[i], _boundaryNormals[j]);
    });
}

// Compute the approximate mass of boundary particles based on [4] equation 4 and 5
void SPH::updateBoundaryMasses() {
    parallelFor(_boundaryPositions.size(), [this] (size_t i) {
        float weight = 0.f;
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _boundaryPositions[i], [this, &weight] (size_t j, const Vector3f &r, float r2) {
            weight += _kernel.poly6(r2);
        });
        _boundaryMasses[i] = _restDensity / (_kernel.poly6Constant * weight);
        _boundaryMasses[i] /= 1.17f;
    });
}

// Computes densities of fluid and boundary particles based on [4] equation 6
void SPH::updateDensities() {
#if HANDLE_BOUNDARIES
    parallelFor(_boundaryPositions.size(), [this] (size_t i) {
        if (!_boundaryActive[i]) {
            return;
        }
        float fluidDensity = 0.f;
        iterateNeighbours(_fluidGrid, _fluidPositions, _boundaryPositions[i], [this, &fluidDensity] (size_t j, const Vector3f &r, float r2) {
            fluidDensity += _kernel.poly6(r2);
        });
        float boundaryDensity = 0.f;
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _boundaryPositions[i], [this, &boundaryDensity] (size_t j, const Vector3f &r, float r2) {
            boundaryDensity += _kernel.poly6(r2) * _boundaryMasses[j];
        });
        float density = _kernel.poly6Constant * _particleMass * fluidDensity;
        density += _kernel.poly6Constant * boundaryDensity;

        _boundaryDensities[i] = density;
    });
#endif

    parallelFor(_fluidPositions.size(), [this] (size_t i) {
        float fluidDensity = 0.f;
        iterateNeighbours(_fluidGrid, _fluidPositions, _fluidPositions[i], [&] (size_t j, const Vector3f &r, float r2) {
            fluidDensity += _kernel.poly6(r2);
        });
        float density = _kernel.poly6Constant * _particleMass * fluidDensity;
#if HANDLE_BOUNDARIES
        float boundaryDensity = 0.f;
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _fluidPositions[i], [&] (size_t j, const Vector3f &r, float r2) {
            boundaryDensity += _kernel.poly6(r2) * _boundaryMasses[j];;
        });
        density += _kernel.poly6Constant * boundaryDensity;
#endif

        _fluidDensities[i] = density;
    });
}

// Compute normals based on [3]
void SPH::updateNormals() {
    parallelFor(_fluidPositions.size(), [this] (size_t i) {
        Vector3f normal;
        iterateNeighbours(_fluidGrid, _fluidPositions, _fluidPositions[i], [&] (size_t j, const Vector3f &r, float r2) {
            normal += _kernel.poly6Grad(r, r2) / _fluidDensities[j];
        });
        normal *= _kernelRadius * _particleMass * _kernel.poly6GradConstant;
        _fluidNormals[i] = normal;
    });
}

void SPH::computeCollisions(std::function<void(size_t i, const Vector3f &n, float d)> handler) {
    for (size_t i = 0; i < _fluidPositions.size(); ++i) {
        const auto &p = _fluidPositions[i];
        if (p.x() < _bounds.min.x()) {
            handler(i, Vector3f(1.f, 0.f, 0.f), _bounds.min.x() - p.x());
        }
        if (p.x() > _bounds.max.x()) {
            handler(i, Vector3f(-1.f, 0.f, 0.f), p.x() - _bounds.max.x());
        }
        if (p.y() < _bounds.min.y()) {
            handler(i, Vector3f(0.f, 1.f, 0.f), _bounds.min.y() - p.y());
        }
        if (p.y() > _bounds.max.y()) {
            handler(i, Vector3f(0.f, -1.f, 0.f), p.y() - _bounds.max.y());
        }
        if (p.z() < _bounds.min.z()) {
            handler(i, Vector3f(0.f, 0.f, 1.f), _bounds.min.z() - p.z());
        }
        if (p.z() > _bounds.max.z()) {
            handler(i, Vector3f(0.f, 0.f, -1.f), p.z() - _bounds.max.z());
        }
    }
}

void SPH::enforceBounds() {
    computeCollisions([&] (size_t i, const Vector3f &n, float d) {
        float c = 0.5f;
        _fluidPositions[i] += n * d;
        _fluidVelocities[i] -= (1 + c) * _fluidVelocities[i].dot(n) * n;
    });
}

void SPH::wcsphUpdateDensitiesAndPressures() {
    parallelFor(_boundaryPositions.size(), [this] (size_t i) {
        if (!_boundaryActive[i]) {
            return;
        }
        float fluidDensity = 0.f;
        iterateNeighbours(_fluidGrid, _fluidPositions, _boundaryPositions[i], [this, &fluidDensity] (size_t j, const Vector3f &r, float r2) {
            fluidDensity += _kernel.poly6(r2);
        });
        float boundaryDensity = 0.f;
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _boundaryPositions[i], [this, &boundaryDensity] (size_t j, const Vector3f &r, float r2) {
            boundaryDensity += _kernel.poly6(r2) * _boundaryMasses[j];
        });
        float density = _kernel.poly6Constant * _particleMass * fluidDensity;
        density += _kernel.poly6Constant * boundaryDensity;

        // Tait pressure (WCSPH)
        float t = density / _restDensity;
        float pressure = wcsph.B * ((t*t)*(t*t)*(t*t)*t - 1.f);

        _boundaryDensities[i] = density;
        _boundaryPressures[i] = pressure;
    });

    parallelFor(_fluidPositions.size(), [this] (size_t i) {
        float fluidDensity = 0.f;
        iterateNeighbours(_fluidGrid, _fluidPositions, _fluidPositions[i], [this, &fluidDensity] (size_t j, const Vector3f &r, float r2) {
            fluidDensity += _kernel.poly6(r2);
        });
        float boundaryDensity = 0.f;
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _fluidPositions[i], [this, &boundaryDensity] (size_t j, const Vector3f &r, float r2) {
            boundaryDensity += _kernel.poly6(r2) * _boundaryMasses[j];;
        });
        float density = _kernel.poly6Constant * _particleMass * fluidDensity;
        density += _kernel.poly6Constant * boundaryDensity;

        // Tait pressure (WCSPH)
        float t = density / _restDensity;
        float pressure = wcsph.B * ((t*t)*(t*t)*(t*t)*t - 1.f);

        _fluidDensities[i] = density;
        _fluidPressures[i] = pressure;
    });
}

void SPH::wcsphUpdateForces() {
    parallelFor(_fluidPositions.size(), [this] (size_t i) {
        Vector3f force(0.f);
        Vector3f forceViscosity;
        Vector3f forceCohesion;
        Vector3f forceCurvature;

        _fluidGrid.lookup(_fluidPositions[i], _kernelRadius, [this, i, &force, &forceCohesion, &forceCurvature, &forceViscosity] (size_t j) {
            const Vector3f &v_i = _fluidVelocities[i];
            const Vector3f &v_j = _fluidVelocities[j];
            const Vector3f &n_i = _fluidNormals[i];
            const Vector3f &n_j = _fluidNormals[j];
            const float &density_i = _fluidDensities[i];
            const float &density_j = _fluidDensities[j];
            const float &pressure_i = _fluidPressures[i];
            const float &pressure_j = _fluidPressures[j];

            if (i != j) {
                Vector3f r = _fluidPositions[i] - _fluidPositions[j];
                float r2 = r.squaredNorm();
                if (r2 < _kernelRadius2 && r2 > 0.00001f) {
                    float rn = std::sqrt(r2);
                    //force -= 0.5f * (pressure_i + pressure_j) * _m / density_j * Kernel::spikyGrad(r);
                    //force -= _particleMass2 * (pressure_i + pressure_j) / (2.f * density_i * density_j) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);

                    // Viscosity force
                    //force += _particleMass2 * _settings.viscosity * (v_j - v_i) / (density_i * density_j) * _kernel.viscosityLaplaceConstant * _kernel.viscosityLaplace(rn);


                    // Pressure force (WCSPH)
                    //if (pressure_i > 0.f || pressure_j > 0.f)
                    force -= _particleMass2 * (pressure_i / sqr(density_i) + pressure_j / sqr(density_j)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
                    //force -= _particleMass2 * (pressure_i / sqr(density_i)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);

                    #if 0
                    // Viscosity force (WCSPH)
                    Vector3f v = (v_i - v_j);
                    if (v.dot(r) < 0.f) {
                        float vu = 2.f * wcsph.viscosity * _kernelRadius * wcsph.cs / (density_i + density_j);
                        force += vu * _particleMass2 * (v.dot(r) / (r2 + 0.001f * sqr(_kernelRadius))) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
                    }
                    #endif

                    // Surface tension force (WCSPH)
                    #if 0
                    float K = 0.1f;
                    Vector3f a = -K * _kernel.poly6Constant * _kernel.poly6(r2) * r / rn;
                    force += _particleMass * a;
                    #endif

                    // Viscosity
                    if (density_j > 0.0001f) {
                        forceViscosity -= (v_i - v_j) * (_kernel.viscosityLaplace(rn) / density_j);
                    }

                    // Surface tension (according to [3])
                    float correctionFactor = 2.f * _restDensity / (density_i + density_j);
                    forceCohesion += correctionFactor * (r / rn) * _kernel.surfaceTension(rn);
                    forceCurvature += correctionFactor * (n_i - n_j);
                } else if (r2 == 0.f) {
                    // Avoid collapsing particles
                    _fluidPositions[j] += Vector3f(1e-5f);
                }
            }
            return true;
        });

#if HANDLE_BOUNDARIES
        _boundaryGrid.lookup(_fluidPositions[i], _kernelRadius, [this, i, &force, &forceCohesion, &forceCurvature, &forceViscosity] (size_t j) {
            const float &density_i = _fluidDensities[i];
            const float &density_j = _boundaryDensities[j];
            const float &pressure_i = _fluidPressures[i];
            const float &pressure_j = _boundaryPressures[j];

            Vector3f r = _fluidPositions[i] - _boundaryPositions[j];
            float r2 = r.squaredNorm();
            if (r2 < _kernelRadius2 && r2 > 0.00001f) {
                float rn = std::sqrt(r2);
                // Pressure force (WCSPH)
                //force -= /*2.f **/ _particleMass2 * (pressure_i / sqr(density_i) + pressure_j / sqr(density_j)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
                force -= _particleMass * _boundaryMasses[j] * (pressure_i / sqr(density_i) + pressure_j / sqr(density_j)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
                //force -= _particleMass * _boundaryMasses[j] * (pressure_i / sqr(density_i)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
            }
            return true;
        });
#endif

        //const float viscosity = 0.0005f;
        //const float viscosity = 0.001f;

        forceViscosity *= _viscosity * _particleMass * _kernel.viscosityLaplaceConstant;

        forceCohesion *= -_surfaceTension * _particleMass2 * _kernel.surfaceTensionConstant;
        forceCurvature *= -_surfaceTension * _particleMass;

        force += forceCohesion + forceCurvature + forceViscosity;
        force += _particleMass * _gravity;

        _fluidForces[i] = force;
    });
}

void SPH::wcsphInit() {

}

void SPH::wcsphUpdate() {
    DebugMonitor::clear();

    Profiler::profile("Grid Update", [&] () {
        _fluidGrid.update(_fluidPositions, [&] (size_t i, size_t j) {
            std::swap(_fluidPositions[i], _fluidPositions[j]);
            std::swap(_fluidVelocities[i], _fluidVelocities[j]);
        });
    });

    Profiler::profile("Activate Boundary", [&] () {
        activateBoundaryParticles();
    });

    Profiler::profile("Density Update", [&] () {
        wcsphUpdateDensitiesAndPressures();
    });

    Profiler::profile("Normal Update", [&] () {
        updateNormals();
    });

    Profiler::profile("Force Update", [&] () {
        wcsphUpdateForces();
    });

    Profiler::profile("Integrate", [&] () {
        parallelFor(_fluidPositions.size(), [&] (size_t i) {
            Vector3f a = _invParticleMass * _fluidForces[i];
            _fluidVelocities[i] += a * _timeStep;
            _fluidPositions[i] += _fluidVelocities[i] * _timeStep;
        });
    });

    Profiler::profile("Collision Update", [&] () {
        enforceBounds();
    });

    _time += _timeStep;
}

void SPH::pcisphUpdateGrid() {
    _fluidGrid.update(_fluidPositions, [&] (size_t i, size_t j) {
        std::swap(_fluidPositions[i], _fluidPositions[j]);
        std::swap(_fluidVelocities[i], _fluidVelocities[j]);
    });
}

void SPH::pcisphUpdateDensityVariationScaling() {
    // Compute density error scaling factor
    Vector3f gradSum;
    float gradDotSum = 0.f;
    for (float x = -_kernelRadius - _particleRadius; x <= _kernelRadius + _particleRadius; x += 2.f * _particleRadius) {
        for (float y = -_kernelRadius - _particleRadius; y <= _kernelRadius + _particleRadius; y += 2.f * _particleRadius) {
            for (float z = -_kernelRadius - _particleRadius; z <= _kernelRadius + _particleRadius; z += 2.f * _particleRadius) {
                Vector3f r = Vector3f(x, y, z);
                float r2 = r.squaredNorm();
                if (r2 < _kernelRadius2) {
                    Vector3f grad = _kernel.poly6GradConstant * _kernel.poly6Grad(r, r2);
                    gradSum += grad;
                    gradDotSum += grad.dot(grad);
                }
            }
        }
    }

    float beta = 2.f * sqr((_particleMass * _timeStep) / _restDensity);
    _densityVariationScaling = -1.f / (beta * (-gradSum.dot(gradSum) - gradDotSum));
#if 0
    DBG("densityVariationScaling = %f", _densityVariationScaling);
#endif
}

// Initialize forces
// - compute all forces that are constant during PCISPH iterations (e.g. viscosity, surface tension, external forces)
// - reset pressures and pressure forces
void SPH::pcisphInitializeForces() {
    parallelFor(_fluidPositions.size(), [&] (size_t i) {
        Vector3f forceViscosity;
        Vector3f forceCohesion;
        Vector3f forceCurvature;

        iterateNeighbours(_fluidGrid, _fluidPositions, _fluidPositions[i], [&] (size_t j, const Vector3f &r, float r2) {
            const Vector3f &v_i = _fluidVelocities[i];
            const Vector3f &v_j = _fluidVelocities[j];
            const Vector3f &n_i = _fluidNormals[i];
            const Vector3f &n_j = _fluidNormals[j];
            const float &density_i = _fluidDensities[i];
            const float &density_j = _fluidDensities[j];

            if (r2 < 1e-7f) {
                return;
            }

            float rn = std::sqrt(r2);

            // Viscosity
            if (density_j > 0.0001f) {
                forceViscosity -= (v_i - v_j) * (_kernel.viscosityLaplace(rn) / density_j);
            }

            // Surface tension (according to [3])
            float correctionFactor = 2.f * _restDensity / (density_i + density_j);
            forceCohesion += correctionFactor * (r / rn) * _kernel.surfaceTension(rn);
            forceCurvature += correctionFactor * (n_i - n_j);
        });

        forceViscosity *= _viscosity * _particleMass * _kernel.viscosityLaplaceConstant;

        forceCohesion *= -_surfaceTension * _particleMass2 * _kernel.surfaceTensionConstant;
        forceCurvature *= -_surfaceTension * _particleMass;

        Vector3f force;
        force += forceCohesion + forceCurvature + forceViscosity;
        force += _particleMass * _gravity;

        _fluidForces[i] = force;
        _fluidPressures[i] = 0.f;
        _fluidPressureForces[i] = Vector3f(0.f);
    });
}

void SPH::pcisphPredictVelocitiesAndPositions() {
    parallelFor(_fluidPositions.size(), [&] (size_t i) {
        Vector3f a = _invParticleMass * (_fluidForces[i] + _fluidPressureForces[i]);
        _fluidVelocitiesNew[i] = _fluidVelocities[i] + a * _timeStep;
        _fluidPositionsNew[i] = _fluidPositions[i] + _fluidVelocitiesNew[i] * _timeStep;
    });
}

void SPH::pcisphUpdatePressures() {
    tbb::enumerable_thread_specific<float> maxDensityVariation(-std::numeric_limits<float>::infinity());
    tbb::enumerable_thread_specific<float> accDensityVariation(0.f);

    parallelFor(_fluidPositions.size(), [&] (size_t i) {
        float fluidDensity = 0.f;
        iterateNeighbours2(_fluidGrid, _fluidPositionsNew, _fluidPositions[i], _fluidPositionsNew[i], [&] (size_t j, const Vector3f &r, float r2) {
            fluidDensity += _kernel.poly6(r2);
        });
        float density = _kernel.poly6Constant * _particleMass * fluidDensity;
#if HANDLE_BOUNDARIES
        float boundaryDensity = 0.f;
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _fluidPositionsNew[i], [&] (size_t j, const Vector3f &r, float r2) {
            boundaryDensity += _kernel.poly6(r2) * _boundaryMasses[j];;
        });
        density += _kernel.poly6Constant * boundaryDensity;
#endif

        float densityVariation = std::max(0.f, density - _restDensity);
        maxDensityVariation.local() = std::max(maxDensityVariation.local(), densityVariation);
        accDensityVariation.local() += densityVariation;

        _fluidPressures[i] += _densityVariationScaling * densityVariation;
    });

    _maxDensityVariation = std::accumulate(maxDensityVariation.begin(), maxDensityVariation.end(), 0.f, [] (float a, float b) { return std::max(a, b); });
    _avgDensityVariation = std::accumulate(accDensityVariation.begin(), accDensityVariation.end(), 0.f) / _fluidPositions.size();

#if 0
    DBG("maxDensityVariation = %f", _maxDensityVariation);
    DBG("avgDensityVariation = %f", _avgDensityVariation);
#endif
}

void SPH::pcisphUpdatePressureForces() {
    parallelFor(_fluidPositions.size(), [&] (size_t i) {
        Vector3f pressureForce;

        iterateNeighbours(_fluidGrid, _fluidPositions, _fluidPositions[i], [&] (size_t j, const Vector3f &r, float r2) {
            if (r2 < 1e-5f) {
                return;
            }

            float rn = std::sqrt(r2);

#if 1
            const float &density_i = _fluidDensities[i];
            const float &density_j = _fluidDensities[j];
            const float &pressure_i = _fluidPressures[i];
            const float &pressure_j = _fluidPressures[j];

            pressureForce -= _particleMass2 * 0.5f * (pressure_i / sqr(density_i) + pressure_j / sqr(density_j)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
#else
            const size_t k = j;//std::min(i, j);
            const float &density_k = _fluidDensities[k];
            const float &pressure_k = _fluidPressures[k];

            pressureForce -= _particleMass2 * (pressure_k / sqr(density_k)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
#endif
        });

#if HANDLE_BOUNDARIES
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _fluidPositions[i], [&] (size_t j, const Vector3f &r, float r2) {
            if (r2 < 1e-5f) {
                return;
            }

            float rn = std::sqrt(r2);

            const float &density_i = _fluidDensities[i];
            const float &density_j = _boundaryDensities[j];
            const float &pressure_i = _fluidPressures[i];
            //const float &pressure_j = _boundaryPressures[j];
            const float &pressure_j = _fluidPressures[i];

            //pressureForce -= _particleMass2 * (pressure_i / sqr(density_i) + pressure_j / sqr(density_j)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
            pressureForce -= _particleMass * _boundaryMasses[j] * 0.5f * (pressure_i / sqr(density_i) + pressure_j / sqr(density_j)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
        });
#endif

        _fluidPressureForces[i] = pressureForce;
    });
}

void SPH::pcisphUpdateVelocitiesAndPositions() {
    tbb::enumerable_thread_specific<float> maxVelocity(0.f);
    tbb::enumerable_thread_specific<float> maxForce(0.f);

    parallelFor(_fluidPositions.size(), [&] (size_t i) {
        Vector3f force = _fluidForces[i] + _fluidPressureForces[i];
        maxForce.local() = std::max(maxForce.local(), force.squaredNorm());
        _fluidVelocitiesNew[i] = _fluidVelocities[i] + _invParticleMass * force * _timeStep;
        _fluidPositionsNew[i] = _fluidPositions[i] + _fluidVelocitiesNew[i] * _timeStep;
        maxVelocity.local() = std::max(maxVelocity.local(), _fluidVelocities[i].squaredNorm());
    });

    _maxVelocity = std::sqrt(std::accumulate(maxVelocity.begin(), maxVelocity.end(), 0.f, [] (float a, float b) { return std::max(a, b); }));
    _maxForce = std::sqrt(std::accumulate(maxForce.begin(), maxForce.end(), 0.f, [] (float a, float b) { return std::max(a, b); }));

    std::swap(_fluidPositionsNew, _fluidPositions);
    std::swap(_fluidVelocitiesNew, _fluidVelocities);
}

void SPH::pcisphInit() {
    // Compute densities
    pcisphUpdateGrid();
    updateDensities();
    float mind = std::numeric_limits<float>::infinity();
    float maxd = -std::numeric_limits<float>::infinity();
    for (auto d : _fluidDensities) {
        mind = std::min(mind, d);
        maxd = std::max(maxd, d);
    }
    DBG("min/max densities = %f/%f", mind, maxd);

    _fluidPositionsPreShock = _fluidPositions;
    _fluidVelocitiesPreShock = _fluidVelocities;

    // Relax initial particle distribution and reset velocities
    pcisphUpdate(10000);
    for (auto &v : _fluidVelocities) {
        v = Vector3f(0.f);
    }

    _time = 0.f;
    _timePreShock = 0.f;
}

void SPH::pcisphUpdate(int maxIterations) {
    DebugMonitor::clear();

    Profiler::profile("Updage Grid", [&] () {
        pcisphUpdateGrid();
    });

    Profiler::profile("Activate Boundary", [&] () {
        activateBoundaryParticles();
    });

    Profiler::profile("Update Densities", [&] () {
        updateDensities();
    });

    Profiler::profile("Update Normals", [&] () {
        updateNormals();
    });

    Profiler::profile("Initialize Forces", [&] () {
        pcisphInitializeForces();
    });

    int k = 0;
    while (k < maxIterations) {
        Profiler::profile("Predict velocities/positions", [&] () {
            pcisphPredictVelocitiesAndPositions();
        });
        Profiler::profile("Update pressures", [&] () {
            pcisphUpdateDensityVariationScaling();
            pcisphUpdatePressures();
        });
        Profiler::profile("Update pressure forces", [&] () {
            pcisphUpdatePressureForces();
        });
        ++k;
        if (k >= 3 && _maxDensityVariation < _maxDensityVariationThreshold) {
            break;
        }
    }
    if (k > 3) {
        DBG("Computed %d pressure iterations!", k);
    }
    DebugMonitor::addItem("pressureIterations", "%d", k);

    Profiler::profile("Update velocities/positions", [&] () {
        pcisphUpdateVelocitiesAndPositions();
    });

    Profiler::profile("Collision Update", [&] () {
        enforceBounds();
    });

    DebugMonitor::addItem("fluidParticles", "%d", _fluidPositions.size());
    DebugMonitor::addItem("boundaryParticles", "%d", _boundaryPositions.size());

    DebugMonitor::addItem("maxDensityVariation", "%.1f", _maxDensityVariation);
    DebugMonitor::addItem("avgDensityVariation", "%.1f", _avgDensityVariation);
    DebugMonitor::addItem("maxVelocity", "%.3f", _maxVelocity);
    DebugMonitor::addItem("maxForce", "%.3f", _maxForce);

    // Prevent divison by zero
    _maxVelocity = std::max(1e-8f, _maxVelocity);
    _maxForce = std::max(1e-8f, _maxForce);

    // Adjust timestep
    if ((0.19f * std::sqrt(_kernelRadius / _maxForce) > _timeStep) &&
        (_maxDensityVariation < 4.5f * _avgDensityVariationThreshold) &&
        (_avgDensityVariation < 0.9f * _avgDensityVariationThreshold) &&
        (0.39f * _kernelRadius / _maxVelocity > _timeStep)) {
        _timeStep *= 1.002f;
    }
    if ((0.2f * std::sqrt(_kernelRadius / _maxForce) < _timeStep) ||
        (_maxDensityVariation > 5.5f * _avgDensityVariationThreshold) ||
        (_avgDensityVariation >= _avgDensityVariationThreshold) ||
        (0.4f * _kernelRadius / _maxVelocity <= _timeStep)) {
        _timeStep *= 0.998f;
    }

    // Detect shock
    //if ((_maxDensityVariation - _prevMaxDensityVariation > 0.5f * (_avgDensityVariationThreshold + _maxDensityVariationThreshold)) ||
    if ((_maxDensityVariation - _prevMaxDensityVariation > _maxDensityVariationThreshold) ||
        (_maxDensityVariation > _maxDensityVariationThreshold) ||
        (0.45f * _kernelRadius / _maxVelocity < _timeStep)) {
        //if (_maxDensityVariation - _prevMaxDensityVariation > 0.5f * (_avgDensityVariationThreshold + _maxDensityVariationThreshold)) {
        if (_maxDensityVariation - _prevMaxDensityVariation > _maxDensityVariationThreshold) {
            DBG("shock due to [1]");
        }
        if (_maxDensityVariation > _maxDensityVariationThreshold) {
            DBG("shock due to [2]");
        }
        if (0.45f * _kernelRadius / _maxVelocity < _timeStep) {
            DBG("shock due to [3]");
        }
        _timeStep = std::min(0.2f * std::sqrt(_kernelRadius / _maxForce), 0.25f * _kernelRadius / _maxVelocity);

        // Go back two timesteps
        _time = _timePreShock;
        _fluidPositions = _fluidPositionsPreShock;
        _fluidVelocities = _fluidVelocitiesPreShock;

        DebugMonitor::addItem("shock", "yes");
    } else {
        _prevMaxDensityVariation = _maxDensityVariation;
        DebugMonitor::addItem("shock", "no");
    }

    // store time, positions and velocities of two timesteps back
    // "new" buffer holds positions/velocities before current integration step!
    _timePreShock = _time;
    std::swap(_fluidPositionsNew, _fluidPositionsPreShock);
    std::swap(_fluidVelocitiesNew, _fluidVelocitiesPreShock);

    _time += _timeStep;

    DebugMonitor::addItem("timeStep", "%.5f", _timeStep);
    DebugMonitor::addItem("time", "%.5f", _time);
}


void SPH::buildScene(const Scene &scene) {
    for (const auto &sceneBox : scene.boxes) {
        switch (sceneBox.type) {
        case Scene::Fluid:
            addFluidParticles(ParticleGenerator::generateVolumeBox(sceneBox.bounds, _particleRadius));
            break;
        case Scene::Boundary:
            addBoundaryParticles(ParticleGenerator::generateBoundaryBox(sceneBox.bounds, _particleRadius));
            _boundaryMeshes.emplace_back(Mesh::createBox(sceneBox.bounds));
            break;
        }
    }
    for (const auto &sceneSphere : scene.spheres) {
        switch (sceneSphere.type) {
        case Scene::Fluid:
            addFluidParticles(ParticleGenerator::generateVolumeSphere(sceneSphere.position, sceneSphere.radius, _particleRadius));
            break;
        case Scene::Boundary:
            addBoundaryParticles(ParticleGenerator::generateBoundarySphere(sceneSphere.position, sceneSphere.radius, _particleRadius));
            _boundaryMeshes.emplace_back(Mesh::createSphere(sceneSphere.position, sceneSphere.radius));
            break;
        }
    }
    for (const auto &sceneMesh : scene.meshes) {
        Mesh mesh = ObjReader::load(sceneMesh.filename);
        switch (sceneMesh.type) {
        case Scene::Fluid:
            addFluidParticles(ParticleGenerator::generateVolumeMesh(mesh, _particleRadius));
            break;
        case Scene::Boundary:
            addBoundaryParticles(ParticleGenerator::generateBoundaryMesh(mesh, _particleRadius));
            _boundaryMeshes.emplace_back(mesh);
            break;
        }
    }

    addBoundaryParticles(ParticleGenerator::generateBoundaryBox(scene.world.bounds, _particleRadius, true));
}

void SPH::addFluidParticles(const ParticleGenerator::Volume &volume) {
    _fluidPositions.insert(_fluidPositions.end(), volume.positions.begin(), volume.positions.end());
}

void SPH::addBoundaryParticles(const ParticleGenerator::Boundary &boundary) {
    _boundaryPositions.insert(_boundaryPositions.end(), boundary.positions.begin(), boundary.positions.end());
    _boundaryNormals.insert(_boundaryNormals.end(), boundary.normals.begin(), boundary.normals.end());
}

std::string SPH::methodToString(Method method) {
    switch (method) {
    case WCSPH: return "wcsph";
    case PCISPH: return "pcisph";
    }
    return "unknown";
}

SPH::Method SPH::stringToMethod(const std::string &str) {
    if (str == "wcsph") {
        return WCSPH;
    } else if (str == "pcisph") {
        return PCISPH;
    } else {
        return PCISPH;
    }
}


} // namespace pbs
