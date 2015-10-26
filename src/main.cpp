#include "Common.h"
#include "SPH.h"

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

#if 0
static inline float unitToRange(float x, float lo, float hi) {
    return lo + nori::clamp(x, 0.f, 1.f) * (hi - lo);
}

static inline float rangeToUnit(float x, float lo, float hi) {
    return nori::clamp((x - lo) / (hi - lo), 0.f, 1.f);
}
#endif

class Viewer : public Screen {
public:
    Viewer(): Screen(Vector2i(1200, 800), "PBS Project") {
        initializeGUI();
        _startTime = std::chrono::high_resolution_clock::now();
        _lastTime = 0;
    }

    ~Viewer() {

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
            _sph.update(dt);
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
        #if 0
        m_focalLengthTextBox->setUnits("mm");
        m_focalLengthTextBox->setValue(tfm::format("%.1f", m_camera->getFocalLength()));
        m_fStopSlider->setValue(rangeToUnit(m_camera->getFStop(), 1.f, 22.f));
        m_fStopTextBox->setValue(tfm::format("%.1f", m_camera->getFStop()));
        m_focusDistanceSlider->setValue(rangeToUnit(m_camera->getFocusDistance(), 0.f, 10.f));
        m_focusDistanceTextBox->setUnits("m");
        m_focusDistanceTextBox->setValue(tfm::format("%.1f", m_camera->getFocusDistance()));
        #endif
    }

    void initializeGUI() {
        #if 0
        Widget *panel;

        m_window = new Window(this, "Lens Settings");
        m_window->setPosition(Vector2i(15, 15));
        m_window->setLayout(new GroupLayout());

        new Label(m_window, "Lens Design", "sans-bold");

        m_lensDesignComboBox = new ComboBox(m_window, m_lensNames);
        m_lensDesignComboBox->setCallback([&] (int i) {
            loadLensSystem(m_lensNames[i]);
            refresh();
        });

        new Label(m_window, "Focal Length", "sans-bold");
        m_focalLengthTextBox = new TextBox(m_window);
        m_focalLengthTextBox->setFixedSize(Vector2i(80, 25));

        new Label(m_window, "F-Stop", "sans-bold");
        panel = new Widget(m_window);
        panel->setLayout(new BoxLayout(BoxLayout::Horizontal, BoxLayout::Middle, 0, 20));
        m_fStopSlider = new Slider(panel);
        m_fStopSlider->setFixedWidth(100);
        m_fStopSlider->setCallback([&] (float f) { m_camera->setFStop(unitToRange(f, 1.f, 22.f)); refresh(); });
        m_fStopTextBox = new TextBox(panel);
        m_fStopTextBox->setFixedSize(Vector2i(80, 25));

        new Label(m_window, "Focus Distance", "sans-bold");
        panel = new Widget(m_window);
        panel->setLayout(new BoxLayout(BoxLayout::Horizontal, BoxLayout::Middle, 0, 20));
        m_focusDistanceSlider = new Slider(panel);
        m_focusDistanceSlider->setFixedWidth(100);
        m_focusDistanceSlider->setCallback([&] (float f) { m_camera->setFocusDistance(unitToRange(f, 0.f, 10.f)); refresh(); });
        m_focusDistanceTextBox = new TextBox(panel);
        m_focusDistanceTextBox->setFixedSize(Vector2i(80, 25));

        new Label(m_window, "Display", "sans-bold");
        panel = new Widget(m_window);
        panel->setLayout(new BoxLayout(BoxLayout::Horizontal, BoxLayout::Middle, 0, 20));
        Slider *scaleSlider = new Slider(panel);
        scaleSlider->setFixedWidth(100);
        scaleSlider->setCallback([&] (float f) { m_scale = f * 10.f + 1.f; });

        refresh();
        #endif

        performLayout(mNVGContext);

        setVisible(true);
    }


private:
    #if 0
    Window *m_window;
    ComboBox *m_lensDesignComboBox;
    TextBox *m_focalLengthTextBox;
    Slider *m_fStopSlider;
    TextBox *m_fStopTextBox;
    Slider *m_focusDistanceSlider;
    TextBox *m_focusDistanceTextBox;
    #endif

    pbs::SPH _sph;

    std::chrono::high_resolution_clock::time_point _startTime;
    double _lastTime;
};

int main(int argc, const char *argv[]) {
    nanogui::init();

    std::unique_ptr<Viewer> screen(new Viewer());
    nanogui::mainloop();
    screen.reset();

    nanogui::shutdown();

    return 0;
}