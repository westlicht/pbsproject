#pragma once

#include "Common.h"
#include "Vector.h"
#include "Box.h"

#include <vector>

namespace pbs {
namespace sph3d {

struct Particle {
    Vector3f p;
    Vector3f v;
    Vector3f force;
    float density;

    Particle(const Vector3f &p) : p(p), v(0.f) {}
};

class Grid {
public:
    Grid(const Box3f &bounds, float cellSize) :
        _bounds(bounds),
        _cellSize(cellSize),
        _invCellSize(1.f / cellSize)
    {
        _size = Vector3i(
            int(std::floor(_bounds.extents().x() / _cellSize)) + 1,
            int(std::floor(_bounds.extents().y() / _cellSize)) + 1,
            int(std::floor(_bounds.extents().z() / _cellSize)) + 1
        );

        _cellOffset.resize(_size.prod() + 1);

        DBG("Grid(): bounds = %s, cellSize = %f, size = %s", _bounds, _cellSize, _size);
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

    void update(const std::vector<Particle> &particles) {
        std::vector<size_t> cellCount(_size.prod(), 0);
        std::vector<size_t> cellIndex(_size.prod(), 0);
        // Count number of particles per cell
        for (size_t i = 0; i < particles.size(); ++i) {
            size_t index = indexLinear(particles[i].p);
            ASSERT(index < size_t(_size.prod()), "particle out of bounds (pos=%s, bounds=%s, index=%d)", particles[i].p, _bounds, index);
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
        // Put particles into cells
        _indices.resize(particles.size());
        for (size_t i = 0; i < particles.size(); ++i) {
            size_t index = indexLinear(particles[i].p);
            _indices[cellIndex[index]++] = i;
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
                        func(_indices[j]);
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
    std::vector<size_t> _indices;
};

class SPH {
public:
    // Simulation constants
    struct Constants {
        // Number of particles expected to be within smoothing kernel support
        static constexpr int supportParticles = 30;
        // Number of particles per unit volume
        static constexpr int particlesPerUnitVolume = 10000;
        // Rest density in kg/m^2
        static constexpr float restDensity = 100.f;
    };

    // Simulation settings
    struct Settings {
        // Stiffness constant
        float stiffness = 3.f;
        // Viscosity
        float viscosity = 1.f;
        // Gravity force
        Vector3f gravity = Vector3f(0.f, -9.81f, 0.f);
    };

    // Kernels
    struct Kernel {
        float h;
        float h2;

        Kernel(float h) : h(h), h2(sqr(h)) {
            poly6Constant = 4.f / (M_PI * std::pow(h, 8.f));
            spikyConstant = 10.f / (M_PI * std::pow(h, 5.f));
            spikyGradConstant = -30.f / (M_PI * std::pow(h, 5.f));
            viscosityLaplaceConstant = 360.f / (29.f * M_PI * std::pow(h, 5.f));
        }

        // Kernels are split into constant and variable part. Arguments are as follows:
        // r: displacement vector
        // r2 = |r|^2 (squared norm of r)
        // rn = |r|   (norm of r)

        float poly6Constant;
        inline float poly6(float r2) const {
            return cube(h2 - r2);
        }

        float spikyConstant;
        inline float spiky(float rn) const {
            return cube(h - rn);
        }

        float spikyGradConstant;
        inline Vector3f spikyGrad(const Vector3f &r, float rn) {
            return sqr(h - rn) * r * (1.f / rn);
        }

        float viscosityLaplaceConstant;
        inline float viscosityLaplace(float rn) {
            return (h - rn);
        }

    };

    SPH() :
        _restSpacing(1.f / std::sqrt(Constants::particlesPerUnitVolume)),
        _particleMass(Constants::restDensity / Constants::particlesPerUnitVolume),
        _particleMass2(sqr(_particleMass)),
        _h(std::sqrt(Constants::supportParticles * (1.f / Constants::particlesPerUnitVolume) / M_PI)),
        _h2(sqr(_h)),
        _kernel(_h),
        _bounds(Box3f(Vector3f(0.f), Vector3f(1.f))),
        _grid(_bounds, _h)
    {
        DBG("Constants::supportParticles = %d", Constants::supportParticles);
        DBG("Constants::particlesPerUnitVolume = %d", Constants::particlesPerUnitVolume);
        DBG("Constants::restDensity = %f", Constants::restDensity);
        DBG("settings.stiffness = %f", _settings.stiffness);
        DBG("settings.viscosity = %f", _settings.viscosity);
        DBG("restSpacing = %f", _restSpacing);
        DBG("particleMass = %f", _particleMass);
        DBG("h = %f", _h);

        voxelizeBox(Box3f(Vector3f(0.25f), Vector3f(0.75f)));
        //voxelizeBox(Box3f(Vector3f(0.1f, 0.5f), Vector3f(0.9f, 0.9f)));
        //voxelizeBox(Box3f(Vector3f(0.3f, 0.5f), Vector3f(0.7f, 0.9f)));
        //voxelizeBox(Box3f(Vector3f(0.4f), Vector3f(0.6f)));

        DBG("simulating %d particles ...", _particles.size());
    }

    const Settings &settings() const { return _settings; }
          Settings &settings()       { return _settings; }

    template<typename Func>
    void iterate(Func func) {
        for (size_t i = 0; i < _particles.size(); ++i) {
            for (size_t j = i + 1; j < _particles.size(); ++j) {
                Func(_particles[i], _particles[j]);
            }
        }
    }

    void computeDensity() {
        for (size_t i = 0; i < _particles.size(); ++i) {
            float density = 0.f;
            _grid.lookup(_particles[i].p, _h, [this, i, &density] (size_t j) {
                Vector3f r = _particles[i].p - _particles[j].p;
                float r2 = r.squaredNorm();
                if (r2 < _h2) {
                    density += _particleMass * _kernel.poly6Constant * _kernel.poly6(r2);
                }
            });
            _particles[i].density = density;
        }
    }

    void computeForces() {
        for (size_t i = 0; i < _particles.size(); ++i) {
            Vector3f force(0.f);
            _grid.lookup(_particles[i].p, _h, [this, i, &force] (size_t j) {
                const float &density_i = _particles[i].density;
                const float &density_j = _particles[j].density;
                const float &p_i = _settings.stiffness * (density_i - Constants::restDensity);
                const float &p_j = _settings.stiffness * (density_j - Constants::restDensity);
                const Vector3f &v_i = _particles[i].v;
                const Vector3f &v_j = _particles[j].v;
                if (i != j) {
                    Vector3f r = _particles[i].p - _particles[j].p;
                    float r2 = r.squaredNorm();
                    if (r2 < _h2 && r2 > 0) {
                        float rn = std::sqrt(r2);
                        // Pressure force
                        //force -= density_i * (p_i / sqr(density_i) + p_j / sqr(density_j)) * Constants::m * Kernel::spikyGrad(r);
                        //force -= 0.5f * (p_i + p_j) * Constants::m / density_j * Kernel::spikyGrad(r);
                        force -= _particleMass2 * (p_i + p_j) / (2.f * density_i * density_j) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);

                        // Viscosity force
                        force += _particleMass2 * _settings.viscosity * (v_j - v_i) / (density_i * density_j) * _kernel.viscosityLaplaceConstant * _kernel.viscosityLaplace(rn);

                        // Surface tension force
                        // according to "Weakly compressible SPH for free surface flows"
                        //float K = 0.001f;
                        //Vector3f a = -K * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
                        //force += _particleMass * a;

                    } else if (r2 == 0.f) {
                        // Avoid collapsing particles
                        _particles[j].p += Vector3f(1e-5f);
                    }
                }
            });

            force += _particleMass * _settings.gravity;

            _particles[i].force = force;
        }
    }

    void computeCollisions(std::function<void(Particle &particle, const Vector3f &n, float d)> handler) {
        for (auto &particle : _particles) {
            if (particle.p.x() < _bounds.min.x()) {
                handler(particle, Vector3f(1.f, 0.f, 0.f), _bounds.min.x() - particle.p.x());
            }
            if (particle.p.x() > _bounds.max.x()) {
                handler(particle, Vector3f(-1.f, 0.f, 0.f), particle.p.x() - _bounds.max.x());
            }
            if (particle.p.y() < _bounds.min.y()) {
                handler(particle, Vector3f(0.f, 1.f, 0.f), _bounds.min.y() - particle.p.y());
            }
            if (particle.p.y() > _bounds.max.y()) {
                handler(particle, Vector3f(0.f, -1.f, 0.f), particle.p.y() - _bounds.max.y());
            }
            if (particle.p.z() < _bounds.min.z()) {
                handler(particle, Vector3f(0.f, 0.f, 1.f), _bounds.min.z() - particle.p.z());
            }
            if (particle.p.z() > _bounds.max.z()) {
                handler(particle, Vector3f(0.f, 0.f, -1.f), particle.p.z() - _bounds.max.z());
            }
        }
    }

    void update(float dt) {
        _t += dt;

        Vector3f gd;
        float t = std::fmod(_t * 0.5f, 4.f);
        if (t < 1.f) {
            _settings.gravity = Vector3f(0.f, -9.81f, 0.f);
        } else if (t < 2.f) {
            _settings.gravity = Vector3f(9.81f, 0.f, 0.f);
        } else if (t < 3.f) {
            _settings.gravity = Vector3f(0.f, 9.81f, 0.f);
        } else {
            _settings.gravity = Vector3f(-9.81f, 0.f, 0.f);
        }

        //_settings.gravity = Vector3f(0.f);

        _grid.update(_particles);

        computeDensity();
        computeForces();


        float invM = 1.f / _particleMass;
        for (auto &particle : _particles) {
            Vector3f a = particle.force * invM;
            particle.v += a * dt;
            particle.p += particle.v * dt;
        }

        // Collision handling
        computeCollisions([] (Particle &particle, const Vector3f &n, float d) {
            float c = 0.5f;
            particle.p += n * d;
            particle.v = particle.v - (1 + c) * particle.v.dot(n) * n;

        });
    }


    void voxelizeBox(const Box3f &box) {
        Vector3i min(
            int(std::floor(box.min.x() / _restSpacing)),
            int(std::floor(box.min.y() / _restSpacing)),
            int(std::floor(box.min.z() / _restSpacing))
        );
        Vector3i max(
            int(std::floor(box.max.x() / _restSpacing)),
            int(std::floor(box.max.y() / _restSpacing)),
            int(std::floor(box.max.z() / _restSpacing))
        );
        for (int z = min.z(); z <= max.z(); ++z) {
            for (int y = min.y(); y <= max.y(); ++y) {
                for (int x = min.x(); x <= max.x(); ++x) {
                    Vector3f p(x * _restSpacing, y * _restSpacing, z * _restSpacing) ;
                    _particles.emplace_back(Particle(p));
                }
            }
        }
    }

    const Box3f &bounds() const { return _bounds; }
    const std::vector<Particle> &particles() const { return _particles; }

private:
    float _restSpacing;
    float _particleMass;
    float _particleMass2;
    float _h;
    float _h2;

    Settings _settings;

    Kernel _kernel;

    Box3f _bounds;
    std::vector<Particle> _particles;
    Grid _grid;

    float _t = 0.f;

};

} // namespace sph3d
} // namespace pbs
