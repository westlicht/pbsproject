#include "Common.h"
#include "SPH3d.h"
#include "Painter.h"
#include "ObjReader.h"
#include "ParticleMesher.h"
#include "Timer.h"

#include <nanogui/checkbox.h>
#include <nanogui/combobox.h>
#include <nanogui/glutil.h>
#include <nanogui/label.h>
#include <nanogui/layout.h>
#include <nanogui/screen.h>
#include <nanogui/slider.h>
#include <nanogui/textbox.h>
#include <nanogui/window.h>

#include <stb_image_write.h>

#include <memory>
#include <chrono>

using namespace nanogui;

// Main screen for 3D viewer.
class Viewer3d : public Screen {
public:
    Viewer3d(): Screen(Vector2i(1200, 800), "PBS Project") {
        initializeGUI();
        _startTime = std::chrono::high_resolution_clock::now();
        _lastTime = 0;

        _viewOrigin = _sph.bounds().center();
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
                }
            }
        }
        return true;
    }


    void drawContents() override {

        if (_isAnimation) {
            float dt = 0.001f;
            int steps = int((1.f / _animationFPS) / dt);
            pbs::DBG("steps = %d", steps);
            for (int i = 0; i < steps; ++i) {
                pbs::Timer timer;
                _sph.update(dt);
                pbs::DBG("timestep took %s", timer.elapsedString());
            }
            createMesh();
        } else if (_isRunning) {
            pbs::Timer timer;
            _sph.update(0.001f);
            pbs::DBG("timestep took %s", timer.elapsedString());
            #if 0
            double dt = 0.01;
            auto currentTime = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - _startTime).count();
            while (_lastTime + dt < currentTime) {
                std::cout << "timestep"
                _sph.update(dt*0.3);
                _lastTime += dt;
            }
            #endif
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

        Matrix4f mvp = proj * view * model;

        _gridPainter->draw(mvp);
        //_boxPainter->draw(mvp, pbs::Box3f(pbs::Vector3f(-0.1f), pbs::Vector3f(0.1f)));
        _boxPainter->draw(mvp, _sph.bounds());

        if (_showParticles) {
            _particlePainter->draw(mvp, _sph.positions());
        }

        if (_showMeshes) {
            _meshPainter->draw(mvp);
        }


        if (_isAnimation) {
            screenshot(tfm::format("frame%04d.png", _animationFrame++));
        }

#if 0

        auto ctx = mNVGContext;

        auto extents = _sph.bounds().extents();
        auto scale = 0.9f * std::min(width() / extents.x(), height() / extents.y());
        extents *= scale;

        nvgBeginFrame(ctx, mSize[0], mSize[1], mPixelRatio);

        nvgResetTransform(ctx);
        nvgTranslate(ctx, 0.5f * (width() - extents.x()), 0.5f * (height() - extents.y()));

        nvgBeginPath(ctx);
        nvgRect(ctx, 0, 0, extents.x(), extents.y());
        nvgStrokeColor(ctx, Color(255, 255));
        nvgStroke(ctx);

        nvgBeginPath(ctx);
        for (const auto &particle : _sph.particles()) {
            auto p = particle.p * scale;
            float radius = 2.f;
            nvgRect(ctx, p.x() - radius, extents.y() - (p.y() - radius), 2.f * radius, 2.f * radius);
        }
        nvgFillColor(ctx, Color(255, 100));
        nvgFill(ctx);

        nvgEndFrame(ctx);
#endif
    }

    void refresh() {
        _stiffnessSlider->setValue(pbs::rangeToUnit(_sph.settings().stiffness, 0.5f, 10.f));
        _stiffnessTextBox->setValue(tfm::format("%.1f", _sph.settings().stiffness));
        _viscositySlider->setValue(pbs::rangeToUnit(_sph.settings().viscosity, 0.5f, 10.f));
        _viscosityTextBox->setValue(tfm::format("%.1f", _sph.settings().viscosity));
    }

    void initializeGUI() {
        _sceneNames = { "Default" };

        Widget *panel;

        _window = new Window(this, "Settings");
        _window->setPosition(Vector2i(15, 15));
        _window->setLayout(new GroupLayout());

        new Label(_window, "Scene", "sans-bold");

        _sceneComboBox = new ComboBox(_window, _sceneNames);
        _sceneComboBox->setCallback([&] (int i) {
            //loadScene(_sceneNames[i]);
            refresh();
        });

        new Label(_window, "Stiffness", "sans-bold");
        panel = new Widget(_window);
        panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));
        _stiffnessSlider = new Slider(panel);
        _stiffnessSlider->setFixedWidth(100);
        _stiffnessSlider->setCallback([&] (float f) { _sph.settings().stiffness = pbs::unitToRange(f, 0.5f, 10.f); refresh(); });
        _stiffnessTextBox = new TextBox(panel);
        _stiffnessTextBox->setFixedSize(Vector2i(80, 25));

        new Label(_window, "Viscosity", "sans-bold");
        panel = new Widget(_window);
        panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));
        _viscositySlider = new Slider(panel);
        _viscositySlider->setFixedWidth(100);
        _viscositySlider->setCallback([&] (float f) { _sph.settings().viscosity = pbs::unitToRange(f, 0.5f, 10.f); refresh(); });
        _viscosityTextBox = new TextBox(panel);
        _viscosityTextBox->setFixedSize(Vector2i(80, 25));

        new Label(_window, "Actions", "sans-bold");
        Button *createMeshButton = new Button(_window, "Create Mesh");
        createMeshButton->setCallback([&] () { createMesh(); refresh(); });
        Button *clearMeshButton = new Button(_window, "Clear Mesh");
        clearMeshButton->setCallback([&] () { clearMesh(); refresh(); });
        Button *renderAnimationButton = new Button(_window, "Render Animation");
        renderAnimationButton->setCallback([&] () { renderAnimation(); refresh(); });
        
        _gridPainter.reset(new pbs::GridPainter());
        _boxPainter.reset(new pbs::BoxPainter());
        _particlePainter.reset(new pbs::ParticlePainter());
        _meshPainter.reset(new pbs::MeshPainter());

        new Label(_window, "Display", "sans-bold");
        CheckBox *showParticlesCheckBox = new CheckBox(_window, "Show Particles");
        showParticlesCheckBox->setChecked(_showParticles);
        showParticlesCheckBox->setCallback([&] (bool b) { _showParticles = b; refresh(); });
        CheckBox *showMeshesCheckBox = new CheckBox(_window, "Show Meshes");
        showMeshesCheckBox->setChecked(_showMeshes);
        showMeshesCheckBox->setCallback([&] (bool b) { _showMeshes = b; refresh(); });

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
        params.supportParticles = _sph.parameters().supportParticles;
        params.particlesPerUnitVolume = _sph.parameters().particlesPerUnitVolume;
        params.restDensity = _sph.parameters().restDensity;
        params.restSpacing = _sph.parameters().restSpacing;
        params.particleMass = _sph.parameters().particleMass;
        params.h = _sph.parameters().h;
        params.isoLevel = 0.2f;

        pbs::MatrixXf positions = _sph.positions();
        pbs::Box3f bounds = _sph.bounds().expanded(_sph.bounds().extents() * 0.05f); // expand bounds by 5% of diagonal
        pbs::Vector3i cells(256);

        pbs::Mesh mesh = _anisotropicMesh ?
            pbs::ParticleMesher::createMeshAnisotropic(positions, bounds, cells, params) :
            pbs::ParticleMesher::createMeshIsotropic(positions, bounds, cells, params);

        //pbs::ObjWriter::save(mesh, "mc.obj");
        _meshPainter->setMesh(mesh);
    }

    void clearMesh() {
        _meshPainter->setMesh(pbs::Mesh());

        screenshot("test.png");
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

    Window *_window;
    ComboBox *_sceneComboBox;
    Slider *_stiffnessSlider;
    TextBox *_stiffnessTextBox;
    Slider *_viscositySlider;
    TextBox *_viscosityTextBox;

    std::vector<std::string> _sceneNames;

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
    std::unique_ptr<pbs::ParticlePainter> _particlePainter;
    std::unique_ptr<pbs::MeshPainter> _meshPainter;

    pbs::sph3d::SPH _sph;

    std::chrono::high_resolution_clock::time_point _startTime;
    double _lastTime;
};
