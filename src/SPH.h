#pragma once

#include "Common.h"
#include "Vector.h"
#include "Box.h"

#include <vector>

namespace pbs {

class SPH {
    static constexpr int Size = 50;
    static constexpr float CellSize = 2.f;

public:

    struct Particle {
        Vector2f p;
        Vector2f v;
        Vector2f force;
        float density;

        Particle(const Vector2f &p) : p(p), v(0.f) {}
    };

    struct Constants {
        // water:
        // rest density (3D) = 998.29 kg/m^3
        // rest density (2D) = 99.89 kg/m^2
        // viscosity = 3.5 Pa s
        // surface tension = 0.0728 N/m
        // gas stiffness = 3 J

        static constexpr int supportedParticles = 20; // number of particles supported within radius h

        static constexpr int ppv = 50*50; // particles per unit volume (1 m^2)
        static constexpr float particleSpacing = 1.f / 50.f; // 1.f / std::sqrt(ppv);

        static constexpr float restDensity = 99.89; // rest density (kg/m^2)
        static constexpr float m = restDensity / ppv; // particle mass (kg)

        static constexpr float h = 0.0504f; // std::sqrt(supportedParticles * (1.f / ppv) / M_PI);

        //static constexpr float h = 2.f; // smoothing radius
        //static constexpr float m = 1.f; // particle mass
        static constexpr float k = 3.f; // stiffness constant
        //static constexpr float restDensity = 1.f; // density in rest state

        static constexpr float viscosity = 3.5f * 0.5;

        static constexpr float g = 9.81f; // gravity
    };

    // Optimized kernel evaluation
    struct Kernel {
        static constexpr float h = Constants::h;
        static constexpr float h2 = sqr(h);
        static constexpr float h4 = sqr(h2);
        static constexpr float h5 = h4*h;
        static constexpr float h6 = h4*h2;
        static constexpr float h8 = sqr(h4);
        static constexpr float h9 = h8*h;
        //static constexpr float normPoly6 = 315.f / (64.f * M_PI * h9);
        static constexpr float normPoly6 = 4.f / (M_PI * h8);
        static constexpr float normPoly6Grad = -945.f / (32.f * M_PI * h9);
        static constexpr float normPoly6Laplace = -945.f / (32.f * M_PI * h9);
//        static constexpr float normSpiky = 15.f / (M_PI * h6);
        static constexpr float normSpiky = 10.f / (M_PI * h5);
//        static constexpr float normSpikyGrad = -45.f / (M_PI * h6);
        static constexpr float normSpikyGrad = -30.f / (M_PI * h5);
        static constexpr float normSpikyLaplace = -90.f / (M_PI * h6);
        static constexpr float normViscosityLaplace = 360.f / (29.f * M_PI * h5);

        static inline float poly6(const Vector2f &r) {
            float r2 = r.squaredNorm();
            return normPoly6 * cube(h2 - r2);
        }

        static inline Vector2f poly6Grad(const Vector2f &r) {
            float r2 = r.squaredNorm();
            return normPoly6Grad * sqr(h2 - r2) * r;
        }

        static inline float poly6Laplace(const Vector2f &r) {
            float r2 = r.squaredNorm();
            return normPoly6Laplace * (h2 - r2) * (3.f * h2 - 7.f * r2);
        }

        static inline float spiky(const Vector2f &r) {
            float rn = r.norm();
            return normSpiky * cube(h - rn);
        }

        static inline Vector2f spikyGrad(const Vector2f &r) {
            float rn = r.norm();
            return normSpikyGrad * sqr(h - rn) * r * (1.f / rn);
        }

        static inline float spikyLaplace(const Vector2f &r) {
            float rn = r.norm();
            return normSpikyLaplace * (h - rn) * (h - 2.f * rn) / rn;
        }

        static inline float viscosityLaplace(const Vector2f &r) {
            float rn = r.norm();
            return normViscosityLaplace * (h - rn);
        }
    };

    class Grid {
    public:
        Grid(const Box2f &bounds, float cellSize) :
            _bounds(bounds),
            _cellSize(cellSize),
            _invCellSize(1.f / cellSize)
        {
            _size = Vector2i(
                int(std::floor(_bounds.extents().x() / _cellSize)) + 1,
                int(std::floor(_bounds.extents().y() / _cellSize)) + 1
            );

            _cellOffset.resize(_size.prod() + 1);

            DBG("Grid(): bounds = %s, cellSize = %f, size = %s", _bounds, _cellSize, _size);
        }

        inline Vector2i index(const Vector2f &pos) {
            return Vector2i(
                int(std::floor((pos.x() - _bounds.min.x()) * _invCellSize)),
                int(std::floor((pos.y() - _bounds.min.y()) * _invCellSize))
            );
        }

        inline size_t indexLinear(const Vector2f &pos) {
            Vector2i i = index(pos);
            return i.y() * _size.x() + i.x();
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
        void lookup(const Vector2f &pos, float radius, Func func) {
            Vector2i min = index(pos - Vector2f(radius)).cwiseMax(Vector2i(0));
            Vector2i max = index(pos + Vector2f(radius)).cwiseMin(_size - Vector2i(1));
            for (int y = min.y(); y <= max.y(); ++y) {
                for (int x = min.x(); x <= max.x(); ++x) {
                    size_t i = y * _size.x() + x;
                    for (size_t j = _cellOffset[i]; j < _cellOffset[i + 1]; ++j) {
                        func(_indices[j]);
                    }
                }
            }
        }

    private:
        Box2f _bounds;
        float _cellSize;
        float _invCellSize;

        Vector2i _size;
        std::vector<size_t> _cellOffset;
        std::vector<size_t> _indices;
    };

    SPH() :
        _bounds(Box2f(Vector2f(0.f), Vector2f(1.f))),
        //_bounds(Box2f(Vector2f(0.f, 0.f), Vector2f(Size, Size))),
        _grid(_bounds, Constants::h)
    {
        #if 0
        _particles.reserve(sqr(Size));
        for (int x = 0; x < Size; ++x) {
            for (int y = Size / 2; y < Size; ++y) {
                _particles.emplace_back(Particle({Vector2f(x + 0.5f, y + 0.5f)}));
            }
        }
        #endif
        //voxelizeBox(Box2f(Vector2f(0.25f), Vector2f(0.75f)));
        voxelizeBox(Box2f(Vector2f(0.1f, 0.5f), Vector2f(0.9f, 0.9f)));
    }

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
            _grid.lookup(_particles[i].p, Constants::h, [this, i, &density] (size_t j) {
                Vector2f r = _particles[i].p - _particles[j].p;
                if (r.squaredNorm() < sqr(Constants::h)) {
                    density += Constants::m * Kernel::poly6(r);
                }
            });
            _particles[i].density = density;
        }
    }

    void computeForces() {
        for (size_t i = 0; i < _particles.size(); ++i) {
            Vector2f force(0.f);
            _grid.lookup(_particles[i].p, Constants::h, [this, i, &force] (size_t j) {
                const float &density_i = _particles[i].density;
                const float &density_j = _particles[j].density;
                const float &p_i = Constants::k * (density_i - Constants::restDensity);
                const float &p_j = Constants::k * (density_j - Constants::restDensity);
                const Vector2f &v_i = _particles[i].v;
                const Vector2f &v_j = _particles[j].v;
                if (i != j) {
                    Vector2f r = _particles[i].p - _particles[j].p;
                    if (r.squaredNorm() < sqr(Constants::h) && r.squaredNorm() > 0) {
                        // Pressure force
                        //force -= density_i * (p_i / sqr(density_i) + p_j / sqr(density_j)) * Constants::m * Kernel::spikyGrad(r);
                        //force -= 0.5f * (p_i + p_j) * Constants::m / density_j * Kernel::spikyGrad(r);
                        force -= sqr(Constants::m) * (p_i + p_j) / (2.f * density_i * density_j) * Kernel::spikyGrad(r);

                        // Viscosity force
                        force += Constants::viscosity * sqr(Constants::m) * (v_j - v_i) / (density_i * density_j) * Kernel::viscosityLaplace(r);


                    }
                }
            });

//            force -= _particles[i].density * Vector2f(0.f, Constants::g);

            _particles[i].force = force;
            //DBG("particle[%d].force = %s", i, force);
        }
    }

    void computeCollisions(std::function<void(Particle &particle, const Vector2f &n, float d)> handler) {
        for (auto &particle : _particles) {
            if (particle.p.x() < _bounds.min.x()) {
                handler(particle, Vector2f(1.f, 0.f), _bounds.min.x() - particle.p.x());
            }
            if (particle.p.x() > _bounds.max.x()) {
                handler(particle, Vector2f(-1.f, 0.f), particle.p.x() - _bounds.max.x());
            }
            if (particle.p.y() < _bounds.min.y()) {
                handler(particle, Vector2f(0.f, 1.f), _bounds.min.y() - particle.p.y());
            }
            if (particle.p.y() > _bounds.max.y()) {
                handler(particle, Vector2f(0.f, -1.f), particle.p.y() - _bounds.max.y());
            }
        }
    }

    void update(float dt) {
        _t += dt;

        _grid.update(_particles);

        computeDensity();
        computeForces();

        //return;

        Vector2f gd;
        float t = std::fmod(_t * 0.2f, 4.f);
        if (t < 1.f) {
            gd = Vector2f(0.f, -1.f);
        } else if (t < 2.f) {
            gd = Vector2f(1.f, 0.f);
        } else if (t < 3.f) {
            gd = Vector2f(0.f, 1.f);
        } else {
            gd = Vector2f(-1.f, 0.f);
        }
        if (std::fmod(_t, 4.f)) {}

        float invM = 1.f / Constants::m;
        for (auto &particle : _particles) {
            Vector2f a = particle.force * invM;
            a += gd * Constants::g * 0.1;
            //Vector2f a = particle.force / particle.density;
            particle.v += a * dt;
            particle.p += particle.v * dt;
            //particle.p.y() -= dt * 10.f;
        }

        // Collision handling
        computeCollisions([] (Particle &particle, const Vector2f &n, float d) {
            particle.p += n * d;
            //particle.v = particle.v - 2.f * particle.v.dot(n) * n;
            float c = 0.5f;
            particle.v = particle.v - (1 + c) * particle.v.dot(n) * n;

        });

        #if 0
        for (auto &particle : _particles) {
            auto &p = particle.p;
            p.x() = clamp(p.x(), _bounds.min.x(), _bounds.max.x());
            p.y() = clamp(p.y(), _bounds.min.y(), _bounds.max.y());
        }
        #endif
    }


    void voxelizeBox(const Box2f &box) {
        Vector2i min(
            int(std::floor(box.min.x() / Constants::particleSpacing)),
            int(std::floor(box.min.y() / Constants::particleSpacing))
        );
        Vector2i max(
            int(std::floor(box.max.x() / Constants::particleSpacing)),
            int(std::floor(box.max.y() / Constants::particleSpacing))
        );
        for (int y = min.y(); y <= max.y(); ++y) {
            for (int x = min.x(); x <= max.x(); ++x) {
                Vector2f p(x * Constants::particleSpacing, y * Constants::particleSpacing) ;
                _particles.emplace_back(Particle(p));
            }
        }
    }

    const Box2f &bounds() const { return _bounds; }
    const std::vector<Particle> &particles() const { return _particles; }

private:
    Box2f _bounds;
    std::vector<Particle> _particles;
    Grid _grid;

    float _t = 0.f;

};

} // namespace pbs
