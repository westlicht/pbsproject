#pragma once

#include "Common.h"
#include "Mesh.h"

#include <nanogui/glutil.h>

namespace pbs {

// Painter for drawing a grid in XZ plane
struct GridPainter {
    GLShader shader;
    int vertexCount = 0;

    GridPainter(int size = 10, float spacing = 0.1f)
    {
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

    void draw(const Matrix4f &mvp, const MatrixXf &positions) {
        shader.bind();
        shader.uploadAttrib("position", positions);
        shader.setUniform("mvp", mvp);
        glPointSize(2);
        glEnable(GL_DEPTH_TEST);
        shader.drawArray(GL_POINTS, 0, positions.cols());
    }
};

// Painter for drawing a mesh.
struct MeshPainter {
    GLShader shader;
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

    void setMesh(const Mesh &mesh) {
        shader.bind();
        shader.uploadAttrib("position", mesh.vertices());
        shader.uploadAttrib("normal", mesh.normals());
        shader.uploadIndices(mesh.triangles());
        triangleCount = mesh.triangles().cols();
    }

    void draw(const Matrix4f &mvp) {
        shader.bind();
        shader.setUniform("mvp", mvp);
        glEnable(GL_DEPTH_TEST);
        shader.drawIndexed(GL_TRIANGLES, 0, triangleCount);
    }
};

} // namespace pbs
