#include "core/Common.h"

#include "../SPH.h"

#include <nanogui/combobox.h>
#include <nanogui/glutil.h>
#include <nanogui/label.h>
#include <nanogui/layout.h>
#include <nanogui/screen.h>
#include <nanogui/slider.h>
#include <nanogui/textbox.h>
#include <nanogui/window.h>

#include <memory>
#include <chrono>

using namespace nanogui;

// Main screen for 2D viewer.
class Viewer2d : public Screen {
public:
    Viewer2d(): Screen(Vector2i(1200, 800), "PBS Project") {
        initializeGUI();
        _startTime = std::chrono::high_resolution_clock::now();
        _lastTime = 0;
    }

    ~Viewer2d() {

    }

    bool mouseMotionEvent(const Vector2i &p, const Vector2i &rel, int button, int modifiers) override {
        if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
        }
        return true;
    }

    bool mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers) override {
        if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
        }
        return true;
    }

    void drawContents() override {

        double dt = 0.01;
        auto currentTime = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - _startTime).count();
        while (_lastTime + dt < currentTime) {
            _sph.update(dt*0.3);
            _lastTime += dt;
        }

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

        refresh();

        performLayout(mNVGContext);

        setVisible(true);
    }

private:
    Window *_window;
    ComboBox *_sceneComboBox;
    Slider *_stiffnessSlider;
    TextBox *_stiffnessTextBox;
    Slider *_viscositySlider;
    TextBox *_viscosityTextBox;

    std::vector<std::string> _sceneNames;

    pbs::sph2d::SPH _sph;

    std::chrono::high_resolution_clock::time_point _startTime;
    double _lastTime;
};
