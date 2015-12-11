#include "core/Common.h"
#include "sim/Engine.h"

#include <nanogui/screen.h>
#include <json11.h>

#include <memory>

namespace pbs {

struct SimulatorSettings {
    enum RenderMode {
        Particles,
        Mesh,
    };

    std::string filename;
    std::string tag;
    int width = 1920;
    int height = 1080;
    float duration = 10.f;
    float timescale = 1.f;
    int framerate = 30;
    RenderMode renderMode;
    bool cacheParticles = false;
    bool cacheMesh = false;
    json11::Json sceneSettings;

    static std::string renderModeToString(RenderMode renderMode) {
        switch (renderMode) {
        case Particles: return "particles";
        case Mesh:      return "mesh";
        default:        return "unknown";
        }
    }

    static RenderMode stringToRenderMode(const std::string &renderMode) {
        if (renderMode == "particles") return Particles;
        if (renderMode == "mesh") return Mesh;
        return Particles;
    }
};

// Simulator screen.
class Simulator : public nanogui::Screen {
public:
    Simulator(const SimulatorSettings &settings);
    ~Simulator();

    bool keyboardEvent(int key, int scancode, int action, int modifiers) override;
    void drawContents() override;

private:
    void initialize();
    void terminate();
    void createVideo();

    void setupEmptyDirectory(const filesystem::path &path);

    SimulatorSettings _settings;
    Engine _engine;

    std::string _basename;
    int _frameIndex = 0;
    float _frameTime = 0.f;
    float _frameInterval;

    Timer _timer;
};

} // namespace pbs
