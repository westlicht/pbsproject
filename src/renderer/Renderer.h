#include "core/Common.h"
#include "sim/Engine.h"
#include "sim/Scene.h"

#include <nanogui/screen.h>
#include <json11.h>

#include <memory>

namespace pbs {

struct RendererSettings {
    enum RenderMode {
        Particles,
        Mesh,
        SLG,
    };

    std::string filenameScene;
    std::string filenameAnimation;
    int width = 1920;
    int height = 1080;
    int framerate = 30;
    int skipToFrame = -1;
    RenderMode renderMode;

    static std::string renderModeToString(RenderMode renderMode) {
        switch (renderMode) {
        case Particles: return "particles";
        case Mesh:      return "mesh";
        case SLG:       return "slg";
        default:        return "unknown";
        }
    }

    static RenderMode stringToRenderMode(const std::string &renderMode) {
        if (renderMode == "particles") return Particles;
        if (renderMode == "mesh") return Mesh;
        if (renderMode == "slg") return SLG;
        return Particles;
    }
};

// Renderer screen.
class Renderer : public nanogui::Screen {
public:
    Renderer(const RendererSettings &settings);
    ~Renderer();

    bool keyboardEvent(int key, int scancode, int action, int modifiers) override;
    void drawContents() override;

private:
    void initialize();
    void terminate();
    void createVideo();

    void renderSLG();

    void setupEmptyDirectory(const filesystem::path &path);

    RendererSettings _settings;
    Engine _engine;
    Scene _animationScene;
    Scene::Camera _startCamera;
    Scene::Camera _endCamera;

    std::string _basename;
    std::string _animationName;
    int _frameFirst;
    int _frameCount;
    int _frameIndex;

    Timer _timer;
};

} // namespace pbs
