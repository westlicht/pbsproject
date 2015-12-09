#include "core/Common.h"
#include "sim/Engine.h"
#include "gui/Panel.h"

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

#include <memory>

namespace pbs {

// Viewer screen.
class Viewer : public nanogui::Screen {
public:
    Viewer();
    ~Viewer();

    // Event handlers
    bool mouseMotionEvent(const nanogui::Vector2i &p, const nanogui::Vector2i &rel, int button, int modifiers) override;
    bool mouseButtonEvent(const nanogui::Vector2i &p, int button, bool down, int modifiers) override;
    bool scrollEvent(const nanogui::Vector2i &p, const nanogui::Vector2f &rel) override;
    bool keyboardEvent(int key, int scancode, int action, int modifiers) override;

    void drawContents() override;

private:
    void initializeGUI();
    void refreshGUI();

    void createMesh();
    void clearMesh();
    void renderAnimation();
    void screenshot(const std::string &filename);

    std::vector<std::string> findScenes(const std::string &path);

    void loadScene(const filesystem::path &path);

    nanogui::Window *_window;
    nanogui::ComboBox *_sceneComboBox;

    nanogui::CheckBox *_showGridCheckBox;
    nanogui::CheckBox *_showParticlesCheckBox;
    nanogui::CheckBox *_showBoundaryParticlesCheckBox;
    nanogui::CheckBox *_showMeshesCheckBox;
    nanogui::CheckBox *_showDebugCheckBox;
    nanogui::CheckBox *_showCacheCheckBox;

    nanogui::Panel *_transportPanel;
    nanogui::Slider *_transportSlider;

    std::vector<std::string> _sceneNames;

    bool _anisotropicMesh = false;
    bool _isRunning = false;

    bool _isAnimation = false;
    float _animationFPS = 60.f;
    int _animationFrame;

    Engine _engine;
};

} // namespace pbs
