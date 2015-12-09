#include "Viewer.h"

#include "core/Common.h"
#include "core/Timer.h"
#include "core/DebugMonitor.h"
#include "core/Profiler.h"
#include "geometry/ObjReader.h"
#include "geometry/ParticleMesher.h"
#include "render/Camera.h"
#include "render/Painter.h"
#include "sim/SPH.h"

#include "Config.h"

#include <nanogui/checkbox.h>
#include <nanogui/combobox.h>
#include <nanogui/glutil.h>
#include <nanogui/label.h>
#include <nanogui/layout.h>
#include <nanogui/screen.h>
#include <nanogui/slider.h>
#include <nanogui/textbox.h>
#include <nanogui/window.h>

#include <filesystem/path.h>
#include <stb_image_write.h>
#include <tinydir.h>

#include <memory>
#include <chrono>

namespace pbs {

Viewer::Viewer() :
    Screen(nanogui::Vector2i(1280, 720), "Fluid Simulator"),
    _engine(mNVGContext, mSize)
{
    _sceneNames = findScenes(SCENES_DIR);
    initializeGUI();
    refreshGUI();
}

Viewer::~Viewer() {
}

bool Viewer::mouseMotionEvent(const nanogui::Vector2i &p, const nanogui::Vector2i &rel, int button, int modifiers) {
    if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
        _engine.camera().mouseMotionEvent(p, rel, button, modifiers);
    }
    return true;
}

bool Viewer::mouseButtonEvent(const nanogui::Vector2i &p, int button, bool down, int modifiers) {
    if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
        _engine.camera().mouseButtonEvent(p, button, down, modifiers);
    }
    return true;
}

bool Viewer::scrollEvent(const nanogui::Vector2i &p, const nanogui::Vector2f &rel) {
    if (!Screen::scrollEvent(p, rel)) {
        _engine.camera().scrollEvent(p, rel);
    }
    return true;
}

bool Viewer::keyboardEvent(int key, int scancode, int action, int modifiers) {
    if (!Screen::keyboardEvent(key, scancode, action, modifiers)) {
        if (action == GLFW_PRESS) {
            switch (key) {
            case GLFW_KEY_ESCAPE:
                setVisible(false);
                break;
            case GLFW_KEY_TAB:
                _window->setVisible(!_window->visible());
                break;
            case GLFW_KEY_HOME:
                //_sph->reset();
                break;
            case GLFW_KEY_SPACE:
                _isRunning = !_isRunning;
                break;
            case GLFW_KEY_G:
                _engine.viewOptions().showGrid = !_engine.viewOptions().showGrid;
                break;
            case GLFW_KEY_P:
                _engine.viewOptions().showParticles = !_engine.viewOptions().showParticles;
                break;
            case GLFW_KEY_B:
                _engine.viewOptions().showBoundaryParticles = !_engine.viewOptions().showBoundaryParticles;
                break;
            case GLFW_KEY_M:
                _engine.viewOptions().showMeshes = !_engine.viewOptions().showMeshes;
                break;
            case GLFW_KEY_D:
                _engine.viewOptions().showDebug = !_engine.viewOptions().showDebug;
                break;
            case GLFW_KEY_S:
                _engine.updateStep();
                break;
            }
        }
        refreshGUI();
    }
    return true;
}

void Viewer::drawContents() {
    if (_isAnimation) {
        _engine.update(1.f / _animationFPS);
    } else if (_isRunning) {
        _engine.updateStep();
        glfwPostEmptyEvent();
    }

    _engine.render();
}

void Viewer::initializeGUI() {
    _window = new nanogui::Window(this, "Settings");
    _window->setPosition(Vector2i(15, 15));
    _window->setLayout(new nanogui::GroupLayout());

    new nanogui::Label(_window, "Scene", "sans-bold");

    _sceneComboBox = new nanogui::ComboBox(_window, _sceneNames);
    _sceneComboBox->setCallback([&] (int i) {
        loadScene(filesystem::path(SCENES_DIR) / _sceneNames[i]);
        refreshGUI();
    });
    // Load default scene
    auto it = std::find(_sceneNames.begin(), _sceneNames.end(), "dambreak.json");
    if (it != _sceneNames.end()) {
        loadScene(filesystem::path(SCENES_DIR) / *it);
        _sceneComboBox->setSelectedIndex(std::distance(_sceneNames.begin(), it));
    }

    new nanogui::Label(_window, "Actions", "sans-bold");
    nanogui::Button *createMeshButton = new nanogui::Button(_window, "Create Mesh");
    createMeshButton->setCallback([&] () { createMesh(); refreshGUI(); });
    nanogui::Button *clearMeshButton = new nanogui::Button(_window, "Clear Mesh");
    clearMeshButton->setCallback([&] () { clearMesh(); refreshGUI(); });
    nanogui::Button *renderAnimationButton = new nanogui::Button(_window, "Render Animation");
    renderAnimationButton->setCallback([&] () { renderAnimation(); refreshGUI(); });

    new nanogui::Label(_window, "Display", "sans-bold");
    _showGridCheckBox = new nanogui::CheckBox(_window, "Show Grid & Domain (G)");
    _showGridCheckBox->setCallback([&] (bool b) { _engine.viewOptions().showGrid = b; refreshGUI(); });
    _showParticlesCheckBox = new nanogui::CheckBox(_window, "Show Particles (P)");
    _showParticlesCheckBox->setCallback([&] (bool b) { _engine.viewOptions().showParticles = b; refreshGUI(); });
    _showBoundaryParticlesCheckBox = new nanogui::CheckBox(_window, "Show Boundary Particles (B)");
    _showBoundaryParticlesCheckBox->setCallback([&] (bool b) { _engine.viewOptions().showBoundaryParticles = b; refreshGUI(); });
    _showMeshesCheckBox = new nanogui::CheckBox(_window, "Show Meshes (M)");
    _showMeshesCheckBox->setCallback([&] (bool b) { _engine.viewOptions().showMeshes = b; refreshGUI(); });
    _showDebugCheckBox = new nanogui::CheckBox(_window, "Show Debug (D)");
    _showDebugCheckBox->setCallback([&] (bool b) { _engine.viewOptions().showDebug = b; refreshGUI(); });

    new nanogui::Label(_window, "Options", "sans-bold");
    nanogui::CheckBox *anisotropicMeshCheckBox = new nanogui::CheckBox(_window, "Anisotropic Mesh");
    anisotropicMeshCheckBox->setChecked(_anisotropicMesh);
    anisotropicMeshCheckBox->setCallback([&] (bool b) { _anisotropicMesh = b; refreshGUI(); });

    performLayout(mNVGContext);
    setVisible(true);
}

void Viewer::refreshGUI() {
    _showGridCheckBox->setChecked(_engine.viewOptions().showGrid);
    _showParticlesCheckBox->setChecked(_engine.viewOptions().showParticles);
    _showBoundaryParticlesCheckBox->setChecked(_engine.viewOptions().showBoundaryParticles);
    _showMeshesCheckBox->setChecked(_engine.viewOptions().showMeshes);
    _showDebugCheckBox->setChecked(_engine.viewOptions().showDebug);
}

void Viewer::createMesh() {
#if 0
    ParticleMesher::Parameters params;
    params.particleRadius = _sph->parameters().particleRadius;
    params.particleDiameter = _sph->parameters().particleDiameter;
    params.kernelRadius = _sph->parameters().kernelRadius;
    params.kernelSupportParticles = _sph->parameters().kernelSupportParticles;
    params.particleMass = _sph->parameters().particleMass;
    params.restDensity = _sph->parameters().restDensity;
    params.isoLevel = 0.2f;

    MatrixXf positions = toMatrix(_sph->fluidPositions());
    Box3f bounds = _sph->bounds().expanded(_sph->bounds().extents() * 0.05f); // expand bounds by 5% of diagonal
    Vector3i cells(256);

    Mesh mesh = _anisotropicMesh ?
        ParticleMesher::createMeshAnisotropic(positions, bounds, cells, params) :
        ParticleMesher::createMeshIsotropic(positions, bounds, cells, params);

    _meshPainter->setMesh(mesh);
#endif
}

void Viewer::clearMesh() {
#if 0
    _meshPainter->setMesh(Mesh());
#endif
}

void Viewer::renderAnimation() {
    _window->setVisible(false);
    _isAnimation = true;
    _animationFrame = 0;
}

void Viewer::screenshot(const std::string &filename) {
    std::unique_ptr<unsigned char[]> pixels(new unsigned char[mSize.prod() * 3]);
    glReadPixels(0, 0, mSize.x(), mSize.y(), GL_RGB, GL_UNSIGNED_BYTE, pixels.get());
    stbi_write_png(filename.c_str(), mSize.x(), mSize.y(), 3, pixels.get() + (3 * mSize.x() * (mSize.y() - 1)), -3 * mSize.x());
}

std::vector<std::string> Viewer::findScenes(const std::string &path) {
    std::vector<std::string> scenes;
    tinydir_dir dir;
    tinydir_open_sorted(&dir, path.c_str());
    for (size_t i = 0; i < dir.n_files; ++i) {
        tinydir_file file;
        tinydir_readfile_n(&dir, &file, i);
        if (file.is_reg) {
            scenes.emplace_back(file.name);
        }
    }
    tinydir_close(&dir);        
    return scenes;
}

void Viewer::loadScene(const filesystem::path &path) {
    _engine.loadScene(path);
}

} // namespace pbs