#pragma once

#include "core/Common.h"
#include "core/Vector.h"
#include "core/Box.h"
#include "geometry/Mesh.h"

#include <nanogui/glutil.h>

namespace pbs {

// Painter for drawing a grid in XZ plane.
struct GridPainter {
    nanogui::GLShader shader;
    int vertexCount = 0;

    GridPainter(int size = 10, float spacing = 0.1f) {
        shader.init(
            "GridPainter",
            // Vertex shader
            "#version 330\n"
            "uniform mat4 mvp;\n"
            "in vec3 position;\n"
            "void main() {\n"
            "    gl_Position = mvp * vec4(position, 1.0);\n"
            "}",
            // Fragment shader
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

    void draw(const nanogui::Matrix4f &mvp) {
        shader.bind();
        shader.setUniform("mvp", mvp);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        shader.drawArray(GL_LINES, 0, vertexCount);
        glDisable(GL_BLEND);
    }
};

// Painter for drawing a bounding box.
struct BoxPainter {
    nanogui::GLShader shader;

    BoxPainter() {
        shader.init(
            "BoxPainter",
            // Vertex shader
            "#version 330\n"
            "uniform mat4 mvp;\n"
            "in vec3 position;\n"
            "void main() {\n"
            "    gl_Position = mvp * vec4(position, 1.0);\n"
            "}",
            // Fragment shader
            "#version 330\n"
            "out vec4 out_color;\n"
            "void main() {\n"
            "    out_color = vec4(vec3(1.0), 0.4);\n"
            "}"
        );
    }

    void draw(const nanogui::Matrix4f &mvp, const Box3f &box) {
        MatrixXf positions(3, 24);
        positions.col(0)  = Vector3f(box.min.x(), box.min.y(), box.min.z());
        positions.col(1)  = Vector3f(box.max.x(), box.min.y(), box.min.z());
        positions.col(2)  = Vector3f(box.min.x(), box.max.y(), box.min.z());
        positions.col(3)  = Vector3f(box.max.x(), box.max.y(), box.min.z());
        positions.col(4)  = Vector3f(box.min.x(), box.min.y(), box.min.z());
        positions.col(5)  = Vector3f(box.min.x(), box.max.y(), box.min.z());
        positions.col(6)  = Vector3f(box.max.x(), box.min.y(), box.min.z());
        positions.col(7)  = Vector3f(box.max.x(), box.max.y(), box.min.z());

        positions.col(8)  = Vector3f(box.min.x(), box.min.y(), box.max.z());
        positions.col(9)  = Vector3f(box.max.x(), box.min.y(), box.max.z());
        positions.col(10) = Vector3f(box.min.x(), box.max.y(), box.max.z());
        positions.col(11) = Vector3f(box.max.x(), box.max.y(), box.max.z());
        positions.col(12) = Vector3f(box.min.x(), box.min.y(), box.max.z());
        positions.col(13) = Vector3f(box.min.x(), box.max.y(), box.max.z());
        positions.col(14) = Vector3f(box.max.x(), box.min.y(), box.max.z());
        positions.col(15) = Vector3f(box.max.x(), box.max.y(), box.max.z());

        positions.col(16) = Vector3f(box.min.x(), box.min.y(), box.min.z());
        positions.col(17) = Vector3f(box.min.x(), box.min.y(), box.max.z());
        positions.col(18) = Vector3f(box.max.x(), box.min.y(), box.min.z());
        positions.col(19) = Vector3f(box.max.x(), box.min.y(), box.max.z());
        positions.col(20) = Vector3f(box.min.x(), box.max.y(), box.min.z());
        positions.col(21) = Vector3f(box.min.x(), box.max.y(), box.max.z());
        positions.col(22) = Vector3f(box.max.x(), box.max.y(), box.min.z());
        positions.col(23) = Vector3f(box.max.x(), box.max.y(), box.max.z());

        shader.bind();
        shader.uploadAttrib("position", positions);
        shader.setUniform("mvp", mvp);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        shader.drawArray(GL_LINES, 0, 24);
        glDisable(GL_BLEND);
    }
};

// Painter for drawing particles.
struct ParticlePainter {
    nanogui::GLShader shader;

    ParticlePainter() {
        shader.init(
            "ParticlePainter",
            // Vertex shader
            "#version 330\n"
            "uniform mat4 mvp;\n"
            "in vec3 position;\n"
            "void main() {\n"
            "    gl_Position = mvp * vec4(position, 1.0);\n"
            "}",
            // Fragment shader
            "#version 330\n"
            "out vec4 out_color;\n"
            "void main() {\n"
            "    out_color = vec4(1.0);\n"
            "}"
        );
    }

    void draw(const nanogui::Matrix4f &mvp, const MatrixXf &positions) {
        shader.bind();
        shader.uploadAttrib("position", positions);
        shader.setUniform("mvp", mvp);
        glPointSize(2);
        glEnable(GL_DEPTH_TEST);
        shader.drawArray(GL_POINTS, 0, positions.cols());
    }
};

// Painter for drawing particle spheres.
struct ParticleSpherePainter {
    nanogui::GLShader shader;

    ParticleSpherePainter() {
        shader.init(
            "ParticleSpherePainter",
            // Vertex shader
            "#version 330\n"
            "uniform mat4 mv;\n"
            "in vec3 position;\n"
            "out vec4 vPosition;\n"
            "void main() {\n"
            "    vPosition = mv * vec4(position, 1.0);\n"
            "}",
            // Fragment shader
            "#version 330\n"
            "uniform vec4 color;\n"
            "in vec2 gPosition;\n"
            "out vec4 out_color;\n"
            "void main() {\n"
            "    vec3 n = vec3(2.0 * gPosition, 0.0);\n"
            "    float r2 = dot(n.xy, n.xy);\n"
            "    if (r2 > 1.0) discard;\n"
            "    n.z = 1.0 - sqrt(r2);\n"
            "    vec3 L = normalize(vec3(1.0));\n"
            "    float d = max(0.0, dot(L, n));\n"
            "    out_color = vec4(d * color.xyz, color.w);\n"
            "}",
            // Geometry shader
            "#version 330\n"
            "layout (points) in;\n"
            "layout (triangle_strip) out;\n"
            "layout (max_vertices = 4) out;\n"
            "uniform mat4 proj;\n"
            "uniform float particleRadius;\n"
            "in vec4 vPosition[];\n"
            "out vec2 gPosition;\n"
            "void main() {\n"
            "    vec4 p = vPosition[0];\n"
            "    gPosition = vec2(-1.0, -1.0);\n"
            "    gl_Position = proj * vec4(p.xy + gPosition * particleRadius, p.zw);\n"
            "    EmitVertex();\n"
            "    gPosition = vec2(-1.0, 1.0);\n"
            "    gl_Position = proj * vec4(p.xy + gPosition * particleRadius, p.zw);\n"
            "    EmitVertex();\n"
            "    gPosition = vec2(1.0, -1.0);\n"
            "    gl_Position = proj * vec4(p.xy + gPosition * particleRadius, p.zw);\n"
            "    EmitVertex();\n"
            "    gPosition = vec2(1.0, 1.0);\n"
            "    gl_Position = proj * vec4(p.xy + gPosition * particleRadius, p.zw);\n"
            "    EmitVertex();\n"
            "    EndPrimitive();\n"
            "}"
        );
    }

    void draw(const nanogui::Matrix4f &mv, const nanogui::Matrix4f &proj, const MatrixXf &positions, const nanogui::Color &color, float particleRadius = 0.03f) {
        shader.bind();
        shader.uploadAttrib("position", positions);
        shader.setUniform("mv", mv);
        shader.setUniform("proj", proj);
        shader.setUniform("particleRadius", particleRadius);
        shader.setUniform("color", color);
        glEnable(GL_DEPTH_TEST);
        shader.drawArray(GL_POINTS, 0, positions.cols());
    }
};

// Painter for drawing particle normals.
struct ParticleNormalPainter {
    nanogui::GLShader shader;

    ParticleNormalPainter() {
        shader.init(
            "ParticleNormalPainter",
            // Vertex shader
            "#version 330\n"
            "uniform mat4 mvp;\n"
            "uniform float normalLength;\n"
            "in vec3 position;\n"
            "in vec3 normal;\n"
            "out vec4 vPositionA;\n"
            "out vec4 vPositionB;\n"
            "void main() {\n"
            "    vPositionA = mvp * vec4(position, 1.0);\n"
            "    vPositionB = mvp * vec4(position + normalLength * normal, 1.0);\n"
            "}",
            // Fragment shader
            "#version 330\n"
            "uniform vec4 color;\n"
            "out vec4 out_color;\n"
            "void main() {\n"
            "    out_color = color;\n"
            "}",
            // Geometry shader
            "#version 330\n"
            "layout (points) in;\n"
            "layout (line_strip) out;\n"
            "layout (max_vertices = 2) out;\n"
            "in vec4 vPositionA[];\n"
            "in vec4 vPositionB[];\n"
            "void main() {\n"
            "    gl_Position = vPositionA[0];\n"
            "    EmitVertex();\n"
            "    gl_Position = vPositionB[0];\n"
            "    EmitVertex();\n"
            "    EndPrimitive();\n"
            "}"
        );
    }

    void draw(const nanogui::Matrix4f &mvp, const MatrixXf &positions, const MatrixXf &normals, const nanogui::Color &color, float normalLength = 0.05f) {
        shader.bind();
        shader.uploadAttrib("position", positions);
        shader.uploadAttrib("normal", normals);
        shader.setUniform("mvp", mvp);
        shader.setUniform("normalLength", normalLength);
        shader.setUniform("color", color);
        glEnable(GL_DEPTH_TEST);
        shader.drawArray(GL_POINTS, 0, positions.cols());
    }
};

// Painter for drawing a mesh.
struct MeshPainter {
    nanogui::GLShader shader;
    int triangleCount = 0;

    MeshPainter() {
        shader.init(
            "MeshPainter",
            // Vertex shader
            "#version 330\n"
            "uniform mat4 mvp;\n"
            "in vec3 position;\n"
            "in vec3 normal;\n"
            "out vec3 vNormal;\n"
            "void main() {\n"
            "    gl_Position = mvp * vec4(position, 1.0);\n"
            "    vNormal = normal;\n"
            "}",
            // Fragment shader
            "#version 330\n"
            "in vec3 vNormal;\n"
            "out vec4 out_color;\n"
            "void main() {\n"
            "    float s = mix(abs(dot(vNormal, vec3(0.0, 1.0, 0.0))), 0.2, 0.8);"
            "    out_color = vec4(s);\n"
            "}"
        );
    }

    MeshPainter(const Mesh &mesh) : MeshPainter() {
        setMesh(mesh);
    }

    void setMesh(const Mesh &mesh) {
        shader.bind();
        shader.uploadAttrib("position", mesh.vertices());
        shader.uploadAttrib("normal", mesh.normals());
        shader.uploadIndices(mesh.triangles());
        triangleCount = mesh.triangles().cols();
    }

    void draw(const nanogui::Matrix4f &mvp) {
        shader.bind();
        shader.setUniform("mvp", mvp);
        glEnable(GL_DEPTH_TEST);
        shader.drawIndexed(GL_TRIANGLES, 0, triangleCount);
    }
};

} // namespace pbs
