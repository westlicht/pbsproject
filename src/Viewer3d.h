#include "Common.h"
#include "SPH3d.h"

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

class Viewer3d : public Screen {
public:
    Viewer3d(): Screen(Vector2i(1200, 800), "PBS Project") {
        initializeGUI();
        _startTime = std::chrono::high_resolution_clock::now();
        _lastTime = 0;
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

    void drawContents() override {

        Matrix4f view, proj, model;
        view = lookAt(Vector3f(0, 0, _viewDistance), Vector3f(0, 0, 0), Vector3f(0, 1, 0));
        const float viewAngle = 30, near = 0.01, far = 100;
        float fH = std::tan(viewAngle / 360.0f * M_PI) * near;
        float fW = fH * (float) mSize.x() / (float) mSize.y();
        proj = frustum(-fW, fW, -fH, fH, near, far);

        model.setIdentity();
        model = translate(model, -_viewOrigin);
        //model = translate(model, Vector3f(0.f, -0.5f, 0.0f));
        model = _arcball.matrix() * model;

        Matrix4f mvp = proj * view * model;

        _gridPainter->draw(mvp);
        _boxPainter->draw(mvp, pbs::Box3f(pbs::Vector3f(-0.1f), pbs::Vector3f(0.1f)));

        pbs::Vector3i size(10);
        MatrixXf positions(3, size.prod());
        int index = 0;
        for (int z = 0; z < size.z(); ++z) {
            for (int y = 0; y < size.y(); ++y) {
                for (int x = 0; x < size.x(); ++x) {
                    pbs::Vector3f p((x + 0.5f) / size.x(), (y + 0.5f) / size.y(), (z + 0.5f) / size.z());
                    p = (p - pbs::Vector3f(0.5f)) * 0.1f;
                    positions.col(index++) = p;
                }
            }
        }
        _particlePainter->draw(mvp, positions);

#if 0
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

        _gridPainter.reset(new GridPainter());
        _boxPainter.reset(new BoxPainter());
        _particlePainter.reset(new ParticlePainter());

        refresh();

        performLayout(mNVGContext);

        setVisible(true);
        framebufferSizeChanged();
    }


private:
    // Painter for drawing a grid in XZ plane
    struct GridPainter {
        GLShader shader;
        int vertexCount;

        GridPainter(int size = 10, float spacing = 0.1f)
        {
            shader.init(
                "GridPainter",

                /* Vertex shader */
                "#version 330\n"
                "uniform mat4 mvp;\n"
                "in vec3 position;\n"
                "void main() {\n"
                "    gl_Position = mvp * vec4(position, 1.0);\n"
                "}",

                /* Fragment shader */
                "#version 330\n"
                "out vec4 out_color;\n"
                "void main() {\n"
                "    out_color = vec4(vec3(1.0), 0.4);\n"
                "}"
            );

            MatrixXf positions(3, (4 * (size * 2 + 1)));

            int index = 0;
            for (int i = -size; i <= size; ++i) {
                positions.col(index++) = Vector3f(i, 0.f, -size) * spacing;
                positions.col(index++) = Vector3f(i, 0.f,  size) * spacing;
                positions.col(index++) = Vector3f(-size, 0.f, i) * spacing;
                positions.col(index++) = Vector3f( size, 0.f, i) * spacing;
            }
            vertexCount = index;

            shader.bind();
            shader.uploadAttrib("position", positions);
        }

        void draw(const Matrix4f &mvp) {
            shader.bind();
            shader.setUniform("mvp", mvp);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            shader.drawArray(GL_LINES, 0, vertexCount);
            glDisable(GL_BLEND);
        }
    };

    // Painter for drawing a bounding box
    struct BoxPainter {
        GLShader shader;

        BoxPainter() {
            shader.init(
                "BoxPainter",

                /* Vertex shader */
                "#version 330\n"
                "uniform mat4 mvp;\n"
                "in vec3 position;\n"
                "void main() {\n"
                "    gl_Position = mvp * vec4(position, 1.0);\n"
                "}",

                /* Fragment shader */
                "#version 330\n"
                "out vec4 out_color;\n"
                "void main() {\n"
                "    out_color = vec4(vec3(1.0), 0.4);\n"
                "}"
            );
        }

        void draw(const Matrix4f &mvp, const pbs::Box3f &box) {
            MatrixXf positions(3, 24);
            positions.col(0)  = pbs::Vector3f(box.min.x(), box.min.y(), box.min.z());
            positions.col(1)  = pbs::Vector3f(box.max.x(), box.min.y(), box.min.z());
            positions.col(2)  = pbs::Vector3f(box.min.x(), box.max.y(), box.min.z());
            positions.col(3)  = pbs::Vector3f(box.max.x(), box.max.y(), box.min.z());
            positions.col(4)  = pbs::Vector3f(box.min.x(), box.min.y(), box.min.z());
            positions.col(5)  = pbs::Vector3f(box.min.x(), box.max.y(), box.min.z());
            positions.col(6)  = pbs::Vector3f(box.max.x(), box.min.y(), box.min.z());
            positions.col(7)  = pbs::Vector3f(box.max.x(), box.max.y(), box.min.z());

            positions.col(8)  = pbs::Vector3f(box.min.x(), box.min.y(), box.max.z());
            positions.col(9)  = pbs::Vector3f(box.max.x(), box.min.y(), box.max.z());
            positions.col(10) = pbs::Vector3f(box.min.x(), box.max.y(), box.max.z());
            positions.col(11) = pbs::Vector3f(box.max.x(), box.max.y(), box.max.z());
            positions.col(12) = pbs::Vector3f(box.min.x(), box.min.y(), box.max.z());
            positions.col(13) = pbs::Vector3f(box.min.x(), box.max.y(), box.max.z());
            positions.col(14) = pbs::Vector3f(box.max.x(), box.min.y(), box.max.z());
            positions.col(15) = pbs::Vector3f(box.max.x(), box.max.y(), box.max.z());

            positions.col(16) = pbs::Vector3f(box.min.x(), box.min.y(), box.min.z());
            positions.col(17) = pbs::Vector3f(box.min.x(), box.min.y(), box.max.z());
            positions.col(18) = pbs::Vector3f(box.max.x(), box.min.y(), box.min.z());
            positions.col(19) = pbs::Vector3f(box.max.x(), box.min.y(), box.max.z());
            positions.col(20) = pbs::Vector3f(box.min.x(), box.max.y(), box.min.z());
            positions.col(21) = pbs::Vector3f(box.min.x(), box.max.y(), box.max.z());
            positions.col(22) = pbs::Vector3f(box.max.x(), box.max.y(), box.min.z());
            positions.col(23) = pbs::Vector3f(box.max.x(), box.max.y(), box.max.z());

            shader.bind();
            shader.uploadAttrib("position", positions);
            shader.setUniform("mvp", mvp);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            shader.drawArray(GL_LINES, 0, 24);
            glDisable(GL_BLEND);
        }
    };

    // Painter for drawing particles
    struct ParticlePainter {
        GLShader shader;

        ParticlePainter() {
            shader.init(
                "ParticlePainter",

                /* Vertex shader */
                "#version 330\n"
                "uniform mat4 mvp;\n"
                "in vec3 position;\n"
                "void main() {\n"
                "    gl_Position = mvp * vec4(position, 1.0);\n"
                "}",

                /* Fragment shader */
                "#version 330\n"
                "out vec4 out_color;\n"
                "void main() {\n"
                "    out_color = vec4(1.0);\n"
                "}"
            );
        }

        void draw(const Matrix4f &mvp, const MatrixXf &positions) {
            shader.bind();
            shader.uploadAttrib("position", positions);
            shader.setUniform("mvp", mvp);
            glPointSize(2);
            glEnable(GL_DEPTH_TEST);
            shader.drawArray(GL_POINTS, 0, positions.cols());
        }
    };

    Window *_window;
    ComboBox *_sceneComboBox;
    Slider *_stiffnessSlider;
    TextBox *_stiffnessTextBox;
    Slider *_viscositySlider;
    TextBox *_viscosityTextBox;

    std::vector<std::string> _sceneNames;

    bool _leftButton = false;
    Arcball _arcball;

    float _viewDistance = 5.f;
    Vector3f _viewOrigin = Vector3f(0.f, 0.f, 0.f);
    std::unique_ptr<GridPainter> _gridPainter;
    std::unique_ptr<BoxPainter> _boxPainter;
    std::unique_ptr<ParticlePainter> _particlePainter;

    pbs::sph3d::SPH _sph;

    std::chrono::high_resolution_clock::time_point _startTime;
    double _lastTime;
};
