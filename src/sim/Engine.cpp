#include "Engine.h"
#include "Scene.h"

#include "core/Common.h"
#include "core/DebugMonitor.h"
#include "core/Profiler.h"
#include "core/FileUtils.h"

#include "geometry/ParticleMesher.h"

#include <stb_image_write.h>

namespace pbs {

static inline MatrixXf toMatrix(const std::vector<Vector3f> &data) {
    MatrixXf result;
    result.resize(3, data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        result.col(i) = data[i];
    }
    return std::move(result);
}

Engine::Engine(NVGcontext *ctx, const Vector2i &size, const Vector2i &renderSize) :
    _ctx(ctx),
    _size(size),
    _renderSize(renderSize)
{
    _gridPainter.reset(new GridPainter());
    _boxPainter.reset(new BoxPainter());
    _particlePainter.reset(new ParticleSpherePainter());
    _particleNormalPainter.reset(new ParticleNormalPainter());
    _fluidMeshPainter.reset(new MeshPainter());

    _framebuffer.init(_renderSize, 1);
}

void Engine::loadScene(const filesystem::path &path, const json11::Json &settings) {
    DBG("Loading scene from '%s' ...", path.str());
    _scene = Scene::load(path.str(), settings);
    DBG("%s", _scene.toString());

    _camera.setResolution(_size);
    _camera.setPosition(_scene.camera.position);
    _camera.setTarget(_scene.camera.target);
    _camera.setUp(_scene.camera.up);
    _camera.setFov(_scene.camera.fov);
    _camera.setNear(_scene.camera.near);
    _camera.setFar(_scene.camera.far);

    _sph.reset(new SPH(_scene));
    std::string cachePath = FileUtils::splitExtension(path.str()).first + ".cache";
    _cache.reset(new Cache(cachePath));

    _boundaryMeshPainters.clear();
    for (const auto &mesh : _sph->boundaryMeshes()) {
        _boundaryMeshPainters.emplace_back(new MeshPainter(mesh));
    }
}

void Engine::update(float dt) {
    _sph->update(dt);
}

void Engine::updateStep() {
    _sph->updateStep();
}

float Engine::time() const {
    return _sph->time();
}

void Engine::createFluidMesh(bool anisotropic) {
    ParticleMesher::Parameters params;
    params.particleRadius = _sph->parameters().particleRadius;
    params.particleDiameter = _sph->parameters().particleDiameter;
    params.kernelRadius = _sph->parameters().kernelRadius;
    params.kernelSupportParticles = _sph->parameters().kernelSupportParticles;
    params.particleMass = _sph->parameters().particleMass;
    params.restDensity = _sph->parameters().restDensity;

    MatrixXf positions = toMatrix(_sph->fluidPositions());
    Box3f bounds = _sph->bounds().expanded(_sph->bounds().extents() * 0.05f); // expand bounds by 5% of diagonal
    Vector3f extents = bounds.extents();

    Vector3i cells(
        int(std::ceil(extents.x() / params.particleRadius)),
        int(std::ceil(extents.y() / params.particleRadius)),
        int(std::ceil(extents.z() / params.particleRadius))
    );

    params.isoLevel = anisotropic ? 0.2f : 0.75f;
    _fluidMesh = anisotropic ? ParticleMesher::createMeshAnisotropic(positions, bounds, cells, params)
                             : ParticleMesher::createMeshIsotropic(positions, bounds, cells, params);
    _fluidMeshPainter->setMesh(_fluidMesh);
}

void Engine::clearFluidMesh() {
    _fluidMesh = Mesh();
    _fluidMeshPainter->setMesh(_fluidMesh);
}

void Engine::render() {
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (!_sph) {
        return;
    }

    nanogui::Matrix4f view = _camera.viewMatrix();
    nanogui::Matrix4f proj = _camera.projMatrix();
    nanogui::Matrix4f mv = view;
    nanogui::Matrix4f mvp = proj * view;

    if (_viewOptions.showDomain) {
        //_gridPainter->draw(mvp);
        _boxPainter->draw(mvp, _sph->bounds());
    }
    if (_viewOptions.showFluidParticles) {
        float particleRadius = _sph->parameters().particleRadius * 2.f;
        _particlePainter->draw(mv, proj, toMatrix(_sph->fluidPositions()), nanogui::Color(0.5f, 0.5f, 1.f, 1.f), particleRadius);
    }
    if (_viewOptions.showFluidMesh) {
        _fluidMeshPainter->draw(mvp, nanogui::Color(0.5f, 0.5f, 1.f, 1.f));
    }
    if (_viewOptions.showBoundaryParticles) {
        float particleRadius = _sph->parameters().particleRadius * 2.f;
        _particlePainter->draw(mv, proj, toMatrix(_sph->boundaryPositions()), nanogui::Color(1.f, 0.5f, 0.5f, 1.f), particleRadius);
        //_particleNormalPainter->draw(mvp, toMatrix(_sph->boundaryPositions()), toMatrix(_sph->boundaryNormals()), nanogui::Color(1.f, 1.f, 1.f, 1.f), particleRadius * 2.f);
    }
    if (_viewOptions.showBoundaryMeshes) {
        for (const auto &painter : _boundaryMeshPainters) {
            painter->draw(mvp, nanogui::Color(1.f, 0.5f, 0.5f, 1.f));
        }
    }
    if (_viewOptions.showDebug) {
        renderDebugOverlay();
    }
}

void Engine::setCacheFrame(int frame) {
    DBG("frame = %d", frame);
    readCache(frame, _viewOptions.showFluidParticles, _viewOptions.showFluidMesh);
}

void Engine::setCachePosition(float position) {
    int frame = int(std::floor((_cache->frameCount() - 1) * position));
    setCacheFrame(frame);
}

void Engine::writeCache(int frame, bool particles, bool mesh) {
    _cache->setFrame(frame);
    if (particles) {
        _cache->writeParticles(_sph->fluidPositions());
    }
    if (mesh) {
        _cache->writeMesh(_fluidMesh);
    }
}

void Engine::readCache(int frame, bool particles, bool mesh) {
    _cache->setFrame(frame);
    if (particles) {
        _cache->readParticles(_sph->fluidPositions());
    }
    if (mesh && _cache->readMesh(_fluidMesh)) {
        _fluidMeshPainter->setMesh(_fluidMesh);
    }
}

void Engine::renderToPNG(const filesystem::path &path) {
    _framebuffer.bind();
    _framebuffer.bindRead();
    glViewport(0, 0, _renderSize.x(), _renderSize.y());

    render();

    std::unique_ptr<unsigned char[]> pixels(new unsigned char[_renderSize.prod() * 3]);
    glReadPixels(0, 0, _renderSize.x(), _renderSize.y(), GL_RGB, GL_UNSIGNED_BYTE, pixels.get());
    stbi_write_png(path.str().c_str(), _renderSize.x(), _renderSize.y(), 3, pixels.get() + (3 * _renderSize.x() * (_renderSize.y() - 1)), -3 * _renderSize.x());

    _framebuffer.release();
}

void Engine::renderDebugOverlay() {
    nvgBeginFrame(_ctx, _size[0], _size[1], 1.f);
    nvgResetTransform(_ctx);
    nvgFillColor(_ctx, nanogui::Color(255, 255));
    nvgFontSize(_ctx, 18.f);
    nvgFontBlur(_ctx, 0.f);

    float w = 250.f;
    float t = 180.f;

    float x = _size[0] - w;
    float y = 28.f;

    auto drawTitle = [&] (const std::string &name) {
        nvgFillColor(_ctx, nanogui::Color(150, 150, 255, 200));
        nvgText(_ctx, x, y, name.c_str(), nullptr);
        y += 20.f;
    };

    auto drawProfilerItem = [&] (const std::string &name, double ms) {
        nvgFillColor(_ctx, nanogui::Color(255, 200));
        nvgText(_ctx, x, y, name.c_str(), nullptr);
        nvgText(_ctx, x + t, y, tfm::format("%.1f ms", ms).c_str(), nullptr);
        y += 20.f;
    };

    auto drawDebugMonitorItem = [&] (const DebugMonitor::Item &item) {
        nvgFillColor(_ctx, nanogui::Color(255, 200));
        nvgText(_ctx, x, y, item.name.c_str(), nullptr);
        nvgText(_ctx, x + t, y, item.value.c_str(), nullptr);
        y += 20.f;
    };

    drawTitle("Profiler");
    double totalTime = 0.0;
    for (const auto &item : Profiler::items()) {
        drawProfilerItem(item.name, item.avg);
        totalTime += item.avg;
    }
    drawProfilerItem("Total", totalTime);

    y += 20.f;

    drawTitle("Debug Monitor");
    for (const auto &item : DebugMonitor::items()) {
        drawDebugMonitorItem(item);
    }

    nvgEndFrame(_ctx);
}


} // namespace pbs
