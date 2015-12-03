#include "Painter.h"

#include "core/Common.h"
#include "core/Timer.h"
#include "geometry/ObjReader.h"
#include "geometry/ParticleMesher.h"
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

using namespace nanogui;

// Main screen for 3D viewer.
class Viewer3d : public Screen {
public:
    Viewer3d() : Screen(Vector2i(1200, 800), "PBS Project") {
        _sceneNames = findScenes(SCENES_DIR);
        _sph.reset(new pbs::SPH(pbs::Scene()));
        _viewOrigin = _sph->bounds().center();        
        initializeGUI();
    }

    ~Viewer3d() {
    }

    void framebufferSizeChanged() {
        _arcball.setSize(mSize);
    }

    bool mouseMotionEvent(const Vector2i &p, const Vector2i &rel, int button, int modifiers) override {
        if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
            if (modifiers == 0) {
                _arcball.motion(p);
            } else if (_leftButton && (modifiers & GLFW_MOD_CONTROL)) {
                _viewDistance = pbs::clamp(_viewDistance + rel.y() * (_viewDistance * 0.01f), 0.01f, 100.f);
            } else if (_leftButton && (modifiers & GLFW_MOD_SHIFT)) {
                Matrix4f view = _arcball.matrix();
                Vector3f left = view.block<1,3>(0, 0);
                Vector3f up = view.block<1,3>(1, 0);
                _viewOrigin -= left * rel.x() * (_viewDistance * 0.001f);
                _viewOrigin += up * rel.y() * (_viewDistance * 0.001f);
            }
        }
        return true;
    }

    bool mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers) override {
        if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
            if (button == GLFW_MOUSE_BUTTON_1) {
                _leftButton = down;
            }
            if (button == GLFW_MOUSE_BUTTON_1 && modifiers == 0) {
                _arcball.button(p, down);
            }
        }
        return true;
    }

    bool scrollEvent(const Vector2i &p, const Vector2f &rel) override {
        if (!Screen::scrollEvent(p, rel)) {
            _viewDistance = pbs::clamp(_viewDistance - rel.y() * (_viewDistance * 0.01f), 0.01f, 100.f);
        }
        return true;
    }

    bool keyboardEvent(int key, int scancode, int action, int modifiers) override {
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
                    std::cout << "rewind" << std::endl;
                    break;
                case GLFW_KEY_SPACE:
                    _isRunning = !_isRunning;
                    break;
                case GLFW_KEY_B:
                    _showBounds = !_showBounds;
                    break;
                case GLFW_KEY_G:
                    _showGrid = !_showGrid;
                    break;
                case GLFW_KEY_P:
                    _showParticles = !_showParticles;
                    break;
                case GLFW_KEY_M:
                    _showMeshes = !_showMeshes;
                    break;
                }
            }
            refresh();
        }
        return true;
    }


    void drawContents() override {

        const float dt = _sph->maxTimestep();

        if (_isAnimation) {
            int steps = int((1.f / _animationFPS) / dt);
            for (int i = 0; i < steps; ++i) {
                _sph->update(dt);
            }
            if (_showMeshes) {
                createMesh();
            }
        } else if (_isRunning) {
            _sph->update(dt);
            glfwPostEmptyEvent();
        }

        Matrix4f view, proj, model;
        view = lookAt(Vector3f(0, 0, _viewDistance), Vector3f(0, 0, 0), Vector3f(0, 1, 0));
        const float viewAngle = 30, near = 0.01, far = 100;
        float fH = std::tan(viewAngle / 360.0f * M_PI) * near;
        float fW = fH * (float) mSize.x() / (float) mSize.y();
        proj = frustum(-fW, fW, -fH, fH, near, far);

        model.setIdentity();
        model = translate(model, -_viewOrigin);
        model = _arcball.matrix() * model;

        Matrix4f mv = view * model;
        Matrix4f mvp = proj * view * model;

        if (_showGrid) {
            _gridPainter->draw(mvp);
        }
        if (_showBounds) {
            _boxPainter->draw(mvp, _sph->bounds());
        }
        if (_showParticles) {
            _particlePainter->draw(mv, proj, _sph->positions(), Color(0.5f, 0.5f, 1.f, 1.f));
            _particlePainter->draw(mv, proj, _sph->blockerPositions(), Color(1.f, 0.5f, 0.5f, 1.f));
        }
        if (_showMeshes) {
            _meshPainter->draw(mvp);
        }


        if (_isAnimation) {
            screenshot(tfm::format("frame%04d.png", _animationFrame++));
        }
    }

    void refresh() {
        _showGridCheckBox->setChecked(_showGrid);
        _showBoundsCheckBox->setChecked(_showBounds);
        _showParticlesCheckBox->setChecked(_showParticles);
        _showMeshesCheckBox->setChecked(_showMeshes);
    }

    void initializeGUI() {
        _window = new Window(this, "Settings");
        _window->setPosition(Vector2i(15, 15));
        _window->setLayout(new GroupLayout());

        new Label(_window, "Scene", "sans-bold");

        _sceneComboBox = new ComboBox(_window, _sceneNames);
        _sceneComboBox->setCallback([&] (int i) {
            loadScene(filesystem::path(SCENES_DIR) / _sceneNames[i]);
            refresh();
        });
        // Load default scene
        auto it = std::find(_sceneNames.begin(), _sceneNames.end(), "drop.json");
        if (it != _sceneNames.end()) {
            loadScene(filesystem::path(SCENES_DIR) / *it);
            _sceneComboBox->setSelectedIndex(std::distance(_sceneNames.begin(), it));
        }

        new Label(_window, "Actions", "sans-bold");
        Button *createMeshButton = new Button(_window, "Create Mesh");
        createMeshButton->setCallback([&] () { createMesh(); refresh(); });
        Button *clearMeshButton = new Button(_window, "Clear Mesh");
        clearMeshButton->setCallback([&] () { clearMesh(); refresh(); });
        Button *renderAnimationButton = new Button(_window, "Render Animation");
        renderAnimationButton->setCallback([&] () { renderAnimation(); refresh(); });
        
        _gridPainter.reset(new pbs::GridPainter());
        _boxPainter.reset(new pbs::BoxPainter());
        _particlePainter.reset(new pbs::SphereParticlePainter());
        _meshPainter.reset(new pbs::MeshPainter());

        new Label(_window, "Display", "sans-bold");
        _showGridCheckBox = new CheckBox(_window, "Show Grid (G)");
        _showGridCheckBox->setCallback([&] (bool b) { _showGrid = b; refresh(); });
        _showBoundsCheckBox = new CheckBox(_window, "Show Bounds (B)");
        _showBoundsCheckBox->setCallback([&] (bool b) { _showBounds = b; refresh(); });
        _showParticlesCheckBox = new CheckBox(_window, "Show Particles (P)");
        _showParticlesCheckBox->setCallback([&] (bool b) { _showParticles = b; refresh(); });
        _showMeshesCheckBox = new CheckBox(_window, "Show Meshes (M)");
        _showMeshesCheckBox->setCallback([&] (bool b) { _showMeshes = b; refresh(); });

        new Label(_window, "Options", "sans-bold");
        CheckBox *anisotropicMeshCheckBox = new CheckBox(_window, "Anisotropic Mesh");
        anisotropicMeshCheckBox->setChecked(_anisotropicMesh);
        anisotropicMeshCheckBox->setCallback([&] (bool b) { _anisotropicMesh = b; refresh(); });

        //auto mesh = pbs::ObjReader::load("test.obj");
        //_meshPainter->setMesh(mesh);

        refresh();

        performLayout(mNVGContext);

        setVisible(true);
        framebufferSizeChanged();
    }


private:
    void createMesh() {
        pbs::ParticleMesher::Parameters params;
        params.particleRadius = _sph->parameters().particleRadius;
        params.particleDiameter = _sph->parameters().particleDiameter;
        params.kernelRadius = _sph->parameters().kernelRadius;
        params.kernelSupportParticles = _sph->parameters().kernelSupportParticles;
        params.particleMass = _sph->parameters().particleMass;
        params.restDensity = _sph->parameters().restDensity;
        params.isoLevel = 0.2f;

        pbs::MatrixXf positions = _sph->positions();
        pbs::Box3f bounds = _sph->bounds().expanded(_sph->bounds().extents() * 0.05f); // expand bounds by 5% of diagonal
        pbs::Vector3i cells(256);

        pbs::Mesh mesh = _anisotropicMesh ?
            pbs::ParticleMesher::createMeshAnisotropic(positions, bounds, cells, params) :
            pbs::ParticleMesher::createMeshIsotropic(positions, bounds, cells, params);

        //pbs::ObjWriter::save(mesh, "mc.obj");
        _meshPainter->setMesh(mesh);
    }

    void clearMesh() {
        _meshPainter->setMesh(pbs::Mesh());
    }

    void renderAnimation() {
        _window->setVisible(false);
        _isAnimation = true;
        _animationFrame = 0;
    }

    void screenshot(const std::string &filename) {
        std::unique_ptr<unsigned char[]> pixels(new unsigned char[mSize.prod() * 3]);
        glReadPixels(0, 0, mSize.x(), mSize.y(), GL_RGB, GL_UNSIGNED_BYTE, pixels.get());
        stbi_write_png(filename.c_str(), mSize.x(), mSize.y(), 3, pixels.get() + (3 * mSize.x() * (mSize.y() - 1)), -3 * mSize.x());
    }

    std::vector<std::string> findScenes(const std::string &path) {
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

    void loadScene(const filesystem::path &path) {
        pbs::DBG("Loading scene from '%s' ...", path.str());
        _sph.reset(new pbs::SPH(pbs::Scene::load(path.str())));
        _viewOrigin = _sph->bounds().center();
    }

    Window *_window;
    ComboBox *_sceneComboBox;

    CheckBox *_showGridCheckBox;
    CheckBox *_showBoundsCheckBox;
    CheckBox *_showParticlesCheckBox;
    CheckBox *_showMeshesCheckBox;

    std::vector<std::string> _sceneNames;

    bool _showGrid = true;
    bool _showBounds = true;
    bool _showParticles = true;
    bool _showMeshes = true;
    bool _anisotropicMesh = false;
    bool _isRunning = false;
    bool _leftButton = false;
    Arcball _arcball;

    bool _isAnimation = false;
    float _animationFPS = 60.f;
    int _animationFrame;

    float _viewDistance = 5.f;
    Vector3f _viewOrigin = Vector3f(0.f, 0.f, 0.f);
    std::unique_ptr<pbs::GridPainter> _gridPainter;
    std::unique_ptr<pbs::BoxPainter> _boxPainter;
    std::unique_ptr<pbs::SphereParticlePainter> _particlePainter;
    std::unique_ptr<pbs::MeshPainter> _meshPainter;

    std::unique_ptr<pbs::SPH> _sph;
};
