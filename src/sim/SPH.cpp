#include "SPH.h"

namespace pbs {

#define HANDLE_BOUNDARIES 1

SPH::SPH(const Scene &scene) {
    _particleRadius = scene.settings.getFloat("particleRadius", _particleRadius);
    _particleRadius2 = sqr(_particleRadius2);
    _particleDiameter = 2.f * _particleRadius;

    _kernelRadius = 4.f * _particleRadius;
    _kernelRadius2 = sqr(_kernelRadius);
    _kernelSupportParticles = int(std::ceil((4.f / 3.f * M_PI * cube(_kernelRadius)) / cube(_particleDiameter)));

    _restDensity = scene.settings.getFloat("restDensity", _restDensity);

    _particleMass = _restDensity / cube(1.f / _particleDiameter);
    _particleMass2 = sqr(_particleMass);

    _gravity = scene.settings.getVector3("gravity", _gravity);

    wcsph.B = _restDensity * sqr(wcsph.cs) / wcsph.gamma;
    wcsph.dt = std::min(0.25f * _kernelRadius / (_particleMass * 9.81f), 0.4f * _kernelRadius / (wcsph.cs * (1.f + 0.6f * wcsph.viscosity)));

    _maxTimestep = 1e-3f;// * 0.1f;


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
    _fluidNormals.resize(_fluidPositions.size());
    _fluidForces.resize(_fluidPositions.size());
    _fluidDensities.resize(_fluidPositions.size());
    _fluidPressures.resize(_fluidPositions.size());

    _boundaryDensities.resize(_boundaryPositions.size());
    _boundaryPressures.resize(_boundaryPositions.size());
    _boundaryActive.resize(_boundaryPositions.size());

    _kernel.init(_kernelRadius);
    _fluidGrid.init(_bounds, _kernelRadius);
    _boundaryGrid.init(_bounds, _kernelRadius);

    // Setup boundary grid
    _boundaryGrid.update(_boundaryPositions, [this] (size_t i, size_t j) {
        std::swap(_boundaryPositions[i], _boundaryPositions[j]);
        std::swap(_boundaryNormals[i], _boundaryNormals[j]);
    });

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

    DBG("# particles = %d", _fluidPositions.size());
    DBG("# boundary particles = %d", _boundaryPositions.size());
}

void SPH::reset() {

}

void SPH::update(float dt) {
    _t += dt;

    {
        ProfileScope profile("Grid Update");
        _fluidGrid.update(_fluidPositions, [this] (size_t i, size_t j) {
            std::swap(_fluidPositions[i], _fluidPositions[j]);
            std::swap(_fluidVelocities[i], _fluidVelocities[j]);
        });
    }

    {
        ProfileScope profile("Activate Boundary");
        activateBoundaryParticles();
    }

    {
        ProfileScope profile("Density Update");
        //computeDensity();
        computeDensitiesAndPressures();
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
        iterate(_fluidPositions.size(), [this, invM, dt] (size_t i) {
            Vector3f a = _fluidForces[i] * invM;
            _fluidVelocities[i] += a * dt;
            _fluidPositions[i] += _fluidVelocities[i] * dt;
        });
    }

    {
        ProfileScope profile("Collision Update");
        // Collision handling
        computeCollisions([this] (size_t i, const Vector3f &n, float d) {
            float c = 0.5f;
            _fluidPositions[i] += n * d;
            _fluidVelocities[i] -= (1 + c) * _fluidVelocities[i].dot(n) * n;
        });
    }
}

void SPH::activateBoundaryParticles() {
    iterate(_boundaryPositions.size(), [this] (size_t i) {
        _boundaryActive[i] = hasNeighbours(_fluidGrid, _fluidPositions, _boundaryPositions[i]);
    });
}

void SPH::computeDensitiesAndPressures() {
#if HANDLE_BOUNDARIES
    iterate(_boundaryPositions.size(), [this] (size_t i) {
        if (!_boundaryActive[i]) {
            return;
        }
        float density = 0.f;
        iterateNeighbours(_fluidGrid, _fluidPositions, _boundaryPositions[i], [this, &density] (size_t j, const Vector3f &r, float r2) {
            density += _kernel.poly6(r2);
        });
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _boundaryPositions[i], [this, &density] (size_t j, const Vector3f &r, float r2) {
            density += _kernel.poly6(r2);
        });
        density *= _particleMass * _kernel.poly6Constant;

        // Tait pressure (WCSPH)
        float t = density / _restDensity;
        float pressure = wcsph.B * ((t*t)*(t*t)*(t*t)*t - 1.f);

        _boundaryDensities[i] = density;
        _boundaryPressures[i] = pressure;
    });
#endif

    iterate(_fluidPositions.size(), [this] (size_t i) {
        float density = 0.f;
        iterateNeighbours(_fluidGrid, _fluidPositions, _fluidPositions[i], [this, &density] (size_t j, const Vector3f &r, float r2) {
            density += _kernel.poly6(r2);
        });
#if HANDLE_BOUNDARIES
        iterateNeighbours(_boundaryGrid, _boundaryPositions, _fluidPositions[i], [this, &density] (size_t j, const Vector3f &r, float r2) {
            density += _kernel.poly6(r2);
        });
#endif
        density *= _particleMass * _kernel.poly6Constant;

        // Tait pressure (WCSPH)
        float t = density / _restDensity;
        float pressure = wcsph.B * ((t*t)*(t*t)*(t*t)*t - 1.f);

        _fluidDensities[i] = density;
        _fluidPressures[i] = pressure;
    });
}

void SPH::computeDensity() {
    iterate(_fluidPositions.size(), [this] (size_t i) {
        float density = 0.f;
        iterateNeighbours(_fluidGrid, _fluidPositions, _fluidPositions[i], [this, &density] (size_t j, const Vector3f &r, float r2) {
            density += _kernel.poly6(r2);
        });
        density *= _particleMass * _kernel.poly6Constant;

        // Tait pressure (WCSPH)
        float t = density / _restDensity;
        float pressure = wcsph.B * ((t*t)*(t*t)*(t*t)*t - 1.f);

        _fluidDensities[i] = density;
        _fluidPressures[i] = pressure;
    });
}

// Compute normals based on [3]
void SPH::computeNormals() {
    iterate(_fluidPositions.size(), [this] (size_t i) {
        Vector3f normal;
        iterateNeighbours(_fluidGrid, _fluidPositions, _fluidPositions[i], [this, &normal] (size_t j, const Vector3f &r, float r2) {
            normal += _kernel.poly6Grad(r, r2) / _fluidDensities[j];
        });
        normal *= _kernelRadius * _particleMass * _kernel.poly6GradConstant;
        _fluidNormals[i] = normal;
    });
}

void SPH::computeForces() {
    iterate(_fluidPositions.size(), [this] (size_t i) {
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
                force -= 2.f * _particleMass2 * (pressure_i / sqr(density_i) + pressure_j / sqr(density_j)) * _kernel.spikyGradConstant * _kernel.spikyGrad(r, rn);
            }
            return true;
        });
#endif

        //const float viscosity = 0.0005f;
        const float viscosity = 0.0001f;
        forceViscosity *= viscosity * _particleMass * _kernel.viscosityLaplaceConstant;

        const float surfaceTension = 2.f;
        forceCohesion *= -surfaceTension * _particleMass2 * _kernel.surfaceTensionConstant;
        forceCurvature *= -surfaceTension * _particleMass;

        force += forceCohesion + forceCurvature + forceViscosity;
        force += _particleMass * _gravity;

        _fluidForces[i] = force;
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



void SPH::buildScene(const Scene &scene) {
    for (const auto &sceneBox : scene.boxes) {
        switch (sceneBox.type) {
        case Scene::Fluid: addFluidParticles(ParticleGenerator::generateVolumeBox(sceneBox.bounds, _particleRadius)); break;
        case Scene::Boundary: addBoundaryParticles(ParticleGenerator::generateBoundaryBox(sceneBox.bounds, _particleRadius)); break;
        }
    }
    for (const auto &sceneSphere : scene.spheres) {
        switch (sceneSphere.type) {
        case Scene::Fluid: addFluidParticles(ParticleGenerator::generateVolumeSphere(sceneSphere.position, sceneSphere.radius, _particleRadius)); break;
        case Scene::Boundary: addBoundaryParticles(ParticleGenerator::generateBoundarySphere(sceneSphere.position, sceneSphere.radius, _particleRadius)); break;
        }
    }
    for (const auto &sceneMesh : scene.meshes) {
        Mesh mesh = ObjReader::load(sceneMesh.filename);
        switch (sceneMesh.type) {
        case Scene::Fluid: addFluidParticles(ParticleGenerator::generateVolumeMesh(mesh, _particleRadius)); break;
        case Scene::Boundary: addBoundaryParticles(ParticleGenerator::generateBoundaryMesh(mesh, _particleRadius)); break;
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


} // namespace pbs
