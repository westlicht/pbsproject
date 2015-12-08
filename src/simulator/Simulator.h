#include "core/Common.h"
#include "sim/Engine.h"

#include <nanogui/screen.h>
#include <json11.h>

#include <memory>

namespace pbs {

struct SimulatorSettings {
    std::string filename;
    int width = 1280;
    int height = 720;
    float duration = 10.f;
    float timescale = 1.f;
    float framerate = 30.f;
    json11::Json sceneSettings;
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

    int _frameIndex = 0;
    float _frameTime = 0.f;
    float _frameInterval;

    Timer _timer;
};

} // namespace pbs
