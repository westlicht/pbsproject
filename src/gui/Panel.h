#pragma once

#include <nanogui/object.h>
#include <nanogui/widget.h>
#include <nanogui/theme.h>

NAMESPACE_BEGIN(nanogui)

class Panel : public Widget {
public:
    Panel(Widget *parent) : Widget(parent) {}

    virtual void draw(NVGcontext *ctx) override {
        nvgSave(ctx);
        nvgBeginPath(ctx);
        nvgRect(ctx, mPos.x(), mPos.y(), mSize.x(), mSize.y());

        nvgFillColor(ctx, mMouseFocus ? mTheme->mWindowFillFocused
                                      : mTheme->mWindowFillUnfocused);
        nvgFill(ctx);
        nvgRestore(ctx);

        Widget::draw(ctx);
    }
};

NAMESPACE_END(nanogui)
