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

    SPH(const Scene &scene) {
        _particleRadius = scene.settings.getFloat("particleRadius", _particleRadius);
        _particleRadius2 = sqr(_particleRadius2);
        _particleDiameter = 2.f * _particleRadius;

        _kernelRadius = 4.f * _particleRadius;
        _kernelRadius2 = sqr(_kernelRadius);
        _kernelSupportParticles = int(std::ceil((4.f / 3.f * M_PI * cube(_kernelRadius)) / cube(_particleDiameter)));

        //_supportParticles = scene.settings.getInteger("supportParticles", _supportParticles);
        //_particlesPerUnitVolume = scene.settings.getInteger("particlesPerUnitVolume", _particlesPerUnitVolume);
        _restDensity = scene.settings.getFloat("restDensity", _restDensity);

        //_restSpacing = 1.f / std::pow(_particlesPerUnitVolume, 1.f / 3.f);
        //_particleMass = _restDensity / _particlesPerUnitVolume;
        _particleMass = _restDensity / cube(1.f / _particleDiameter);
        _particleMass2 = sqr(_particleMass);
        //_h = std::pow((3.f * _supportParticles) / (4.f * M_PI * _particlesPerUnitVolume), 1.f / 3.f);
        //_h = _restSpacing * 2.f;
        //_kernelRadius2 = sqr(_h);

        _gravity = scene.settings.getVector3("gravity", _gravity);

        wcsph.B = _restDensity * sqr(wcsph.cs) / wcsph.gamma;
        wcsph.dt = std::min(0.25f * _kernelRadius / (_particleMass * 9.81f), 0.4f * _kernelRadius / (wcsph.cs * (1.f + 0.6f * wcsph.viscosity)));

        _maxTimestep = 1e-3f;

        _bounds = scene.world.bounds;
        _kernel.init(_kernelRadius);
        _grid.init(_bounds, _kernelRadius);

        DBG("particleRadius = %f", _particleRadius);
        DBG("kernelRadius = %f", _kernelRadius);
        DBG("kernelSupportParticles = %d", _kernelSupportParticles);
        DBG("restDensity = %f", _restDensity);
        DBG("particleMass = %f", _particleMass);
        DBG("gravity = %s", _gravity);

        DBG("wcsph.gamma = %f", wcsph.gamma);
        DBG("wcsph.cs = %f", wcsph.cs);
        DBG("wcsph.B = %f", wcsph.B);
        DBG("wcsph.viscosity = %f", wcsph.viscosity);
        DBG("wcsph.dt = %f", wcsph.dt);

        for (const auto &box : scene.boxes) {
            voxelizeBox(box.bounds);
        }
        for (const auto &sphere : scene.spheres) {
            voxelizeSphere(sphere.position, sphere.radius);
        }
        for (const auto &mesh : scene.meshes) {
            voxelizeMesh(mesh);
        }

        //voxelizeBox(Box3f(Vector3f(0.25f), Vector3f(0.75f)));
        //voxelizeBox(Box3f(Vector3f(0.1f, 0.5f), Vector3f(0.9f, 0.9f)));
        //voxelizeBox(Box3f(Vector3f(0.3f, 0.5f), Vector3f(0.7f, 0.9f)));
        //voxelizeBox(Box3f(Vector3f(0.4f), Vector3f(0.6f)));


        DBG("simulating %d particles ...", _positions.size());

        _velocities.resize(_positions.size());
        _normals.resize(_positions.size());
        _forces.resize(_positions.size());
        _densities.resize(_positions.size());
        _pressures.resize(_positions.size());
    }

    // iterate i=0..count-1 calling func(i)
    template<typename Func>
    void iterate(size_t count, Func func) {
#if USE_TBB
        tbb::parallel_for(0ul, count, 1ul, [func] (size_t i) { func(i); });
#else
        for (size_t i = 0; i < count; ++i) { func(i); }
#endif
    }

    // iterate over all neighbours around p, calling func(j, r, r2)
    template<typename Func>
    void iterateNeighbours(const Grid &grid, const std::vector<Vector3f> &positions, const Vector3f &p, Func func) {
        grid.lookup(p, _kernelRadius, [&] (size_t j) {
            Vector3f r = p - positions[j];
            float r2 = r.squaredNorm();
            if (r2 < _kernelRadius2) {
                func(j, r, r2);
            }
        });
    }

    void computeDensity() {
        iterate(_positions.size(), [this] (size_t i) {
            float density = 0.f;
            iterateNeighbours(_grid, _positions, _positions[i], [this, &density] (size_t j, const Vector3f &r, float r2) {
                density += _kernel.poly6(r2);
            });
            density *= _particleMass * _kernel.poly6Constant;

            // Tait pressure (WCSPH)
            float t = density / _restDensity;
            float pressure = wcsph.B * ((t*t)*(t*t)*(t*t)*t - 1.f);

            _densities[i] = density;
            _pressures[i] = pressure;
        });
    }

    // Compute normals based on [3]
    void computeNormals() {
        iterate(_positions.size(), [this] (size_t i) {
            Vector3f normal;
            iterateNeighbours(_grid, _positions, _positions[i], [this, &normal] (size_t j, const Vector3f &r, float r2) {
                normal += _kernel.poly6Grad(r, r2) / _densities[j];
            });
            normal *= _kernelRadius * _particleMass * _kernel.poly6GradConstant;
            _normals[i] = normal;
        });
    }

    void computeForces() {
        iterate(_positions.size(), [this] (size_t i) {
            Vector3f force(0.f);
            Vector3f forceViscosity;
            Vector3f forceCohesion;
            Vector3f forceCurvature;

            _grid.lookup(_positions[i], _kernelRadius, [this, i, &force, &forceCohesion, &forceCurvature, &forceViscosity] (size_t j) {
                const Vector3f &v_i = _velocities[i];
                const Vector3f &v_j = _velocities[j];
                const Vector3f &n_i = _normals[i];
                const Vector3f &n_j = _normals[j];
                const float &density_i = _densities[i];
                const float &density_j = _densities[j];
                const float &pressure_i = _pressures[i];
                const float &pressure_j = _pressures[j];

                if (i != j) {
                    Vector3f r = _positions[i] - _positions[j];
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
                        _positions[j] += Vector3f(1e-5f);
                    }
                }
            });

            const float viscosity = 0.0005f;
            forceViscosity *= viscosity * _particleMass * _kernel.viscosityLaplaceConstant;

            const float surfaceTension = 2.f;
            forceCohesion *= -surfaceTension * _particleMass2 * _kernel.surfaceTensionConstant;
            forceCurvature *= -surfaceTension * _particleMass;

            force += forceCohesion + forceCurvature + forceViscosity;
            force += _particleMass * _gravity;

            _forces[i] = force;
        });
    }

    void computeCollisions(std::function<void(size_t i, const Vector3f &n, float d)> handler) {
        for (size_t i = 0; i < _positions.size(); ++i) {
            const auto &p = _positions[i];
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

    void update(float dt) {
        _t += dt;

        Vector3f gd;
        float t = std::fmod(_t * 0.5f, 4.f);
        if (t < 1.f) {
            _gravity = Vector3f(0.f, -9.81f, 0.f);
        } else if (t < 2.f) {
            _gravity = Vector3f(9.81f, 0.f, 0.f);
        } else if (t < 3.f) {
            _gravity = Vector3f(0.f, 9.81f, 0.f);
        } else {
            _gravity = Vector3f(-9.81f, 0.f, 0.f);
        }

        //_gravity = Vector3f(0.f);
        _gravity = Vector3f(0.f, -9.81f, 0.f);

        {
            ProfileScope profile("Grid Update");
            _grid.update(_positions, [this] (size_t i, size_t j) {
                std::swap(_positions[i], _positions[j]);
                std::swap(_velocities[i], _velocities[j]);
            });
        }

        {
            ProfileScope profile("Density Update");
            computeDensity();
        }

        {
            ProfileScope profile("Normal Update");
            computeNormals();
        }

        {
            ProfileScope profile("Force Update");
            computeForces();
        }

        {
            ProfileScope profile("Integrate");
            float invM = 1.f / _particleMass;
            iterate(_positions.size(), [this, invM, dt] (size_t i) {
                Vector3f a = _forces[i] * invM;
                _velocities[i] += a * dt;
                _positions[i] += _velocities[i] * dt;
            });
        }

        {
            ProfileScope profile("Collision Update");

            // Collision handling
            computeCollisions([this] (size_t i, const Vector3f &n, float d) {
                float c = 0.5f;
                _positions[i] += n * d;
                _velocities[i] -= (1 + c) * _velocities[i].dot(n) * n;
            });
        }

        //Profiler::dump();
    }


    void voxelizeBox(const Box3f &box) {
        Vector3i min(
            int(std::ceil(box.min.x() / _particleDiameter)),
            int(std::ceil(box.min.y() / _particleDiameter)),
            int(std::ceil(box.min.z() / _particleDiameter))
        );
        Vector3i max(
            int(std::floor(box.max.x() / _particleDiameter)),
            int(std::floor(box.max.y() / _particleDiameter)),
            int(std::floor(box.max.z() / _particleDiameter))
        );
        for (int z = min.z(); z <= max.z(); ++z) {
            for (int y = min.y(); y <= max.y(); ++y) {
                for (int x = min.x(); x <= max.x(); ++x) {
                    Vector3f p(x * _particleDiameter, y * _particleDiameter, z * _particleDiameter);
                    _positions.emplace_back(p);
                }
            }
        }
    }

    void voxelizeSphere(const Vector3f &pos, float radius) {
        Vector3i min(
            int(std::ceil((pos.x() - radius) / _particleDiameter)),
            int(std::ceil((pos.y() - radius) / _particleDiameter)),
            int(std::ceil((pos.z() - radius) / _particleDiameter))
        );
        Vector3i max(
            int(std::floor((pos.x() + radius) / _particleDiameter)),
            int(std::floor((pos.y() + radius) / _particleDiameter)),
            int(std::floor((pos.z() + radius) / _particleDiameter))
        );
        float r2 = sqr(radius);
        for (int z = min.z(); z <= max.z(); ++z) {
            for (int y = min.y(); y <= max.y(); ++y) {
                for (int x = min.x(); x <= max.x(); ++x) {
                    Vector3f p(x * _particleDiameter, y * _particleDiameter, z * _particleDiameter);
                    if ((p - pos).squaredNorm() <= r2) {
                        _positions.emplace_back(p);
                    }
                }
            }
        }
    }

    void voxelizeMesh(const Scene::Mesh &sceneMesh) {
        Mesh mesh = ObjReader::load(sceneMesh.filename);
        if (sceneMesh.type == Scene::Liquid) {
            Voxelizer::voxelize(mesh, _particleDiameter, _positions);
        } else {
            auto positions = ParticleGenerator::generateSurfaceParticles(mesh, 1000.f);
            _blockerPositions.insert(_blockerPositions.end(), positions.begin(), positions.end());
        }
    }

    const Box3f &bounds() const { return _bounds; }

    // Returns a set of simulation parameters
    Parameters parameters() const {
        Parameters params;
        params.particleRadius = _particleRadius;
        params.particleDiameter = _particleDiameter;
        params.kernelRadius = _kernelRadius;
        params.kernelSupportParticles = _kernelSupportParticles;
        params.particleMass = _particleMass;
        params.restDensity = _restDensity;
        return params;
    }

    float maxTimestep() const { return _maxTimestep; }

    // Returns particle positions in matrix form
    MatrixXf positions() const {
        MatrixXf positions;
        positions.resize(3, _positions.size());
        for (size_t i = 0; i < _positions.size(); ++i) {
            positions.col(i) = _positions[i];
        }
        return std::move(positions);
    }

    // Returns particle positions in matrix form
    MatrixXf blockerPositions() const {
        MatrixXf positions;
        positions.resize(3, _blockerPositions.size());
        for (size_t i = 0; i < _blockerPositions.size(); ++i) {
            positions.col(i) = _blockerPositions[i];
        }
        return std::move(positions);
    }

private:
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
    Grid _grid;

    // Fluid particle buffers
    std::vector<Vector3f> _positions;
    std::vector<Vector3f> _velocities;
    std::vector<Vector3f> _normals;
    std::vector<Vector3f> _forces;
    std::vector<float> _densities;
    std::vector<float> _pressures;

    // Blocker particle buffers
    std::vector<Vector3f> _blockerPositions;

    float _t = 0.f;
};

} // namespace pbs
