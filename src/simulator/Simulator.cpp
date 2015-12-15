#include "Simulator.h"

#include "core/FileUtils.h"
#include "sim/Scene.h"
#include "sim/SPH.h"

#include <exec-stream.h>
#include <tinydir.h>

namespace pbs {

Simulator::Simulator(const SimulatorSettings &settings) :
    nanogui::Screen(nanogui::Vector2i(1280, 720), "Fluid Simulator"),
    _settings(settings),
    _engine(mNVGContext, mSize, Vector2i(settings.width, settings.height))
{
    setVisible(true);

    _engine.loadScene(settings.filename, settings.sceneSettings);
    _engine.viewOptions().showDebug = false;

    _basename = FileUtils::splitExtension(_settings.filename).first;
    if (!_settings.tag.empty()) {
        _basename += "-" + _settings.tag;
    }

    _frameInterval = _settings.timescale / _settings.framerate;

    setupEmptyDirectory("images");

    _timer.reset();
}

Simulator::~Simulator() {
}

bool Simulator::keyboardEvent(int key, int scancode, int action, int modifiers) {
    if (!Screen::keyboardEvent(key, scancode, action, modifiers)) {
        if (action == GLFW_PRESS) {
            switch (key) {
            case GLFW_KEY_ESCAPE:
                terminate();
                break;
            }
        }
    }
    return true;
}

void Simulator::drawContents() {

    bool writeFrame = _engine.time() >= _frameTime;

    _engine.render();

    if (writeFrame) {
        if (_settings.renderMode == SimulatorSettings::Mesh || _settings.cacheMesh) {
            _engine.createFluidMesh();
        }

        _engine.viewOptions().showDebug = false;
        if (_settings.renderMode == SimulatorSettings::Mesh) {
            _engine.viewOptions().showFluidParticles = false;
            _engine.viewOptions().showFluidMesh = true;

        }

        _engine.renderToPNG(tfm::format("images/frame%04d.png", _frameIndex));
        if (_settings.cacheParticles || _settings.cacheMesh) {
            _engine.writeCache(_frameIndex, _settings.cacheParticles, _settings.cacheMesh);
        }
        _engine.clearFluidMesh();
        _frameTime += _frameInterval;
        ++_frameIndex;

        _engine.viewOptions().showDebug = true;
        _engine.viewOptions().showFluidParticles = true;
        _engine.viewOptions().showFluidMesh = false;
    }

    if (_engine.time() >= _settings.duration) {
        terminate();
    }

    // Show time estimate
    float elapsed = _timer.elapsed();
    float progress = _engine.time() / _settings.duration;
    float eta = progress != 0.f ? elapsed / progress - elapsed : 0.f;
    std::cout << tfm::format("\rProgress %.1f%% (Elapsed: %s ETA: %s)", progress * 100.f, timeString(elapsed), timeString(eta)) << std::flush;

    _engine.updateStep();
    glfwPostEmptyEvent();
}

void Simulator::initialize() {

}

void Simulator::terminate() {
    if (_settings.cacheParticles || _settings.cacheMesh) {
        _engine.cache().setFrameCount(_frameIndex);
        _engine.cache().commit();
    }
    createVideo();
    setVisible(false);
}

void Simulator::createVideo() {
    // FFMPEG candidate locations
    std::vector<std::string> candidates = {
        "/opt/local/bin/ffmpeg",
        "/usr/bin/avconv"
    };

    std::string ffmpeg;
    for (const auto &candidate : candidates) {
        if (filesystem::path(candidate).exists()) {
            ffmpeg = candidate;
            break;
        }
    }
    if (ffmpeg.empty()) {
        std::cerr << "Cannot find FFMPEG to encode video!" << std::endl;
        return;
    }

    std::string filename = _basename + ".mp4";

    std::cout << std::endl;
    std::cout << tfm::format("Encoding video '%s' ...", filename) << std::endl;

    try {
        std::string arguments = tfm::format("-y -framerate %d -i images/frame%%04d.png -pix_fmt yuv420p -vcodec libx264 -r %d -preset slow -crf 10 %s", _settings.framerate, _settings.framerate, filename);
        exec_stream_t es;
        // Set 5 minutes timeout
        es.set_wait_timeout(exec_stream_t::s_out | exec_stream_t::s_err | exec_stream_t::s_child, 300*1000);
        es.start(ffmpeg, arguments);

        std::string s;
        while (std::getline(es.err(), s).good()) {
            std::cout << s << std::endl;
        }

        if (!es.close()) {
            std::cerr << "FFMPEG timeout!" << std::endl;
        }
    } catch (const std::exception &e) {
        std::cerr << "exec-stream error: " << e.what() << "\n";
    }
}

void Simulator::setupEmptyDirectory(const filesystem::path &path) {
    FileUtils::createDir(path.str());
    tinydir_dir dir;
    tinydir_open_sorted(&dir, path.str().c_str());
    for (size_t i = 0; i < dir.n_files; ++i) {
        tinydir_file file;
        tinydir_readfile_n(&dir, &file, i);
        if (file.is_reg) {
            FileUtils::deleteFile((path / file.name).str());
        }
    }
    tinydir_close(&dir);
}

} // namespace pbs

