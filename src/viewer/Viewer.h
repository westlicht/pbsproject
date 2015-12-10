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

struct ViewerSettings {
    std::string filename;
};

// Viewer screen.
class Viewer : public nanogui::Screen {
public:
    Viewer(const ViewerSettings &settings);
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

    void loadScene(const filesystem::path &path);

    nanogui::Window *_window;

    nanogui::CheckBox *_showDomainCheckBox;
    nanogui::CheckBox *_showFluidParticlesCheckBox;
    nanogui::CheckBox *_showFluidMeshCheckBox;
    nanogui::CheckBox *_showBoundaryParticlesCheckBox;
    nanogui::CheckBox *_showBoundaryMeshesCheckBox;
    nanogui::CheckBox *_showDebugCheckBox;
    nanogui::CheckBox *_showCacheCheckBox;

    nanogui::Panel *_transportPanel;
    nanogui::Slider *_transportSlider;

    bool _anisotropicMesh = true;
    bool _isRunning = false;

    Engine _engine;
};

} // namespace pbs
