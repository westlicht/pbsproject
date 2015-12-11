#pragma once

#include <nanogui/object.h>
#include <nanogui/glutil.h>

NAMESPACE_BEGIN(nanogui)

class GLFramebufferExt : public GLFramebuffer {
public:
    void bindRead() {
        glBindFramebuffer(GL_READ_FRAMEBUFFER, mFramebuffer);
        glReadBuffer(GL_COLOR_ATTACHMENT0);
    }
};

NAMESPACE_END(nanogui)
