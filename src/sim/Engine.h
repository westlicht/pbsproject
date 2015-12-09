#pragma once

#include "Cache.h"
#include "SPH.h"

#include "core/Common.h"
#include "core/Vector.h"
#include "render/Camera.h"
#include "render/Painter.h"

#include <filesystem/path.h>
#include <json11.h>

#include <memory>

namespace pbs {

class Engine {
public:
    struct ViewOptions {
        bool showGrid = true;
        bool showParticles = true;
        bool showBoundaryParticles = false;
        bool showMeshes = true;
        bool showDebug = true;
        bool showCache = false;
    };

    Engine(NVGcontext *ctx, const Vector2i &size);

    const Camera &camera() const { return _camera; }
          Camera &camera()       { return _camera; }

    const ViewOptions &viewOptions() const { return _viewOptions; }
          ViewOptions &viewOptions()       { return _viewOptions; }


    void loadScene(const filesystem::path &path, const json11::Json &settings = json11::Json());

    void update(float dt);
    void updateStep();
    float time() const;

    void render();

    Cache &cache() { return *_cache; }
    void setCachePosition(float position);
    void writeCache(int frame);
    void readCache(int frame);

    void savePng(const filesystem::path &path);

private:
    void renderDebugOverlay();

    NVGcontext *_ctx;
    Vector2i _size;

    Camera _camera;
    std::unique_ptr<SPH> _sph;
    std::unique_ptr<Cache> _cache;

    ViewOptions _viewOptions;

    std::unique_ptr<GridPainter> _gridPainter;
    std::unique_ptr<BoxPainter> _boxPainter;
    std::unique_ptr<ParticleSpherePainter> _particlePainter;
    std::unique_ptr<ParticleNormalPainter> _particleNormalPainter;
    std::unique_ptr<MeshPainter> _meshPainter;
    std::vector<std::unique_ptr<MeshPainter>> _boundaryMeshPainters;
};

} // namespace pbs
