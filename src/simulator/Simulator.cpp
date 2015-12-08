#include "Simulator.h"

#include "core/FileUtils.h"
#include "sim/Scene.h"
#include "sim/SPH.h"

namespace pbs {

Simulator::Simulator(const SimulatorSettings &settings) :
    nanogui::Screen(nanogui::Vector2i(settings.width, settings.height), "Fluid Simulator"),
    _settings(settings),
    _engine(mNVGContext, mSize)
{
    setVisible(true);

    _engine.loadScene(settings.filename, settings.sceneSettings);
    _engine.viewOptions().showDebug = false;

    _frameInterval = _settings.timescale / _settings.framerate;


    FileUtils::createDir("images");

    _timer.reset();
}

Simulator::~Simulator() {
}

bool Simulator::keyboardEvent(int key, int scancode, int action, int modifiers) {
    if (!Screen::keyboardEvent(key, scancode, action, modifiers)) {
        if (action == GLFW_PRESS) {
            switch (key) {
            case GLFW_KEY_ESCAPE:
                setVisible(false);
                break;
            }
        }
    }
    return true;
}

void Simulator::drawContents() {
    _engine.render();

    if (_engine.time() >= _frameTime) {
        _engine.savePng(tfm::format("images/frame%04d.png", _frameIndex));
        _frameTime += _frameInterval;
        ++_frameIndex;
    }

    if (_engine.time() >= _settings.duration) {
        setVisible(false);
    }

    // Show time estimate
    float elapsed = _timer.elapsed() / 1000.f;
    float progress = _engine.time() / _settings.duration;
    float eta = progress != 0.f ? elapsed / progress - elapsed : 0.f;
    std::cout << tfm::format("\rProgress %.1f%% (Elapsed: %s ETA: %s)", progress * 100.f, timeString(elapsed), timeString(eta)) << std::flush;

    _engine.updateStep();
    glfwPostEmptyEvent();
}

} // namespace pbs

