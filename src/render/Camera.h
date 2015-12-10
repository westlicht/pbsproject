#pragma once

#include "core/Common.h"

#include <nanogui/glutil.h>

namespace pbs {

// Camera class piggy-back on Arcball
class Camera {
public:
    void setResolution(const nanogui::Vector2i &resolution) { _resolution = resolution; }
    void setPosition(const nanogui::Vector3f &position) { _position = position; }
    void setTarget(const nanogui::Vector3f &target) { _target = target; }
    void setUp(const nanogui::Vector3f &up) { _up = up; }
    void setFov(float fov) { _fov = fov; }
    void setNear(float near) { _near = near; }
    void setFar(float far) { _far = far; }

    const nanogui::Vector2i &resolution() const { return _resolution; }
    const nanogui::Vector3f &position() const { return _position; }
    const nanogui::Vector3f &target() const { return _target; }
    const nanogui::Vector3f &up() const { return _up; }
    float fov() const { return _fov; }
    float near() const { return _near; }
    float far() const { return _far; }

    nanogui::Matrix4f viewMatrix() const {
        return nanogui::lookAt(_position, _target, _up);
    }

    nanogui::Matrix4f projMatrix() const {
        float fH = std::tan(_fov / 360.0f * M_PI) * _near;
        float fW = fH * float(_resolution.x()) / float(_resolution.y());
        return nanogui::frustum(-fW, fW, -fH, fH, _near, _far);
    }

    bool mouseMotionEvent(const nanogui::Vector2i &p, const nanogui::Vector2i &rel, int button, int modifiers) {
        if (modifiers == 0) {
            _arcball.motion(p);
            if (_leftButton) {
                nanogui::Matrix4f view = _arcball.matrix();
                float viewDistance = (_position - _target).norm();
                nanogui::Vector3f dir = view.block<3,3>(0,0).transpose() * nanogui::Vector3f(0.f, 0.f, 1.f);
                _position = _target + dir * viewDistance;
                _up = view.block<1,3>(1,0);
            }

        } else if (_leftButton && (modifiers & GLFW_MOD_CONTROL)) {
            float viewDistance = (_position - _target).norm();
            viewDistance = clamp(viewDistance + 0.25f * rel.y() * (viewDistance * 0.01f), 0.01f, 10000.f);
            _position = _target + (_position - _target).normalized() * viewDistance;
        } else if (_leftButton && (modifiers & GLFW_MOD_SHIFT)) {
            float viewDistance = (_position - _target).norm();
            nanogui::Matrix4f view = viewMatrix();
            nanogui::Vector3f left = view.block<1,3>(0,0);
            nanogui::Vector3f up = view.block<1,3>(1,0);
            nanogui::Vector3f offset = -left * rel.x() * (viewDistance * 0.001f) + up * rel.y() * (viewDistance * 0.001f);
            _position += offset;
            _target += offset;
        }
        return true;
    }

    bool mouseButtonEvent(const nanogui::Vector2i &p, int button, bool down, int modifiers) {
        if (button == GLFW_MOUSE_BUTTON_1) {
            _leftButton = down;
        }
        if (button == GLFW_MOUSE_BUTTON_1 && modifiers == 0) {
            _arcball.setSize(_resolution);
            _arcball.setState(nanogui::Quaternionf(viewMatrix().block<3,3>(0,0)));
            _arcball.button(p, down);
        }
        return true;
    }

    bool scrollEvent(const nanogui::Vector2i &p, const nanogui::Vector2f &rel) {
        float viewDistance = (_position - _target).norm();
        viewDistance = clamp(viewDistance - rel.y() * (viewDistance * 0.01f), 0.01f, 10000.f);
        _position = _target + (_position - _target).normalized() * viewDistance;
        return true;
    }

private:
    nanogui::Vector2i _resolution;
    nanogui::Vector3f _position;
    nanogui::Vector3f _target;
    nanogui::Vector3f _up;
    float _fov;
    float _near;
    float _far;

    nanogui::Arcball _arcball;
    bool _leftButton = false;
};

} // namespace pbs
