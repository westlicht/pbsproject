#pragma once

#include "Cache.h"
#include "SPH.h"

#include "core/Common.h"
#include "core/Vector.h"
#include "render/Camera.h"
#include "render/Painter.h"
#include "render/GLFramebufferExt.h"

#include <filesystem/path.h>
#include <json11.h>

#include <memory>

namespace pbs {

class Engine {
public:
    struct ViewOptions {
        bool showDomain = true;
        bool showFluidParticles = true;
        bool showFluidMesh = false;
        bool showBoundaryParticles = false;
        bool showBoundaryMeshes = true;
        bool showDebug = true;
        bool showCache = false;
    };

    Engine(NVGcontext *ctx, const Vector2i &size, const Vector2i &renderSize);

    const Camera &camera() const { return _camera; }
          Camera &camera()       { return _camera; }

    const ViewOptions &viewOptions() const { return _viewOptions; }
          ViewOptions &viewOptions()       { return _viewOptions; }

    const Scene &scene() const { return _scene; }
    const SPH &sph() const { return *_sph; }
    
    const Cache &cache() const { return *_cache; }
          Cache &cache()       { return *_cache; }

    void loadScene(const filesystem::path &path, const json11::Json &settings = json11::Json());

    void update(float dt);
    void updateStep();
    float time() const;

    const Mesh &fluidMesh() const { return _fluidMesh; }
    void createFluidMesh(bool anisotropic = false);
    void clearFluidMesh();

    void render();

    void setCacheFrame(int frame);
    void setCachePosition(float position);
    void writeCache(int frame, bool particles, bool mesh);
    void readCache(int frame, bool particles, bool mesh);

    void renderToPNG(const filesystem::path &path);

private:
    void renderDebugOverlay();

    NVGcontext *_ctx;
    Vector2i _size;
    Vector2i _renderSize;

    Camera _camera;
    Scene _scene;
    std::unique_ptr<SPH> _sph;
    std::unique_ptr<Cache> _cache;

    ViewOptions _viewOptions;

    Mesh _fluidMesh;

    std::unique_ptr<GridPainter> _gridPainter;
    std::unique_ptr<BoxPainter> _boxPainter;
    std::unique_ptr<ParticleSpherePainter> _particlePainter;
    std::unique_ptr<ParticleNormalPainter> _particleNormalPainter;
    std::unique_ptr<MeshPainter> _fluidMeshPainter;
    std::vector<std::unique_ptr<MeshPainter>> _boundaryMeshPainters;

    nanogui::GLFramebufferExt _framebuffer;
};

} // namespace pbs
