#pragma once

#include "core/Common.h"

#include <ostream>
#include <istream>

namespace pbs {

// Holds a triangle mesh.
// Data is stored in columns of the matrices.
class Mesh {
public:
    Mesh();

    const MatrixXf &vertices() const { return _vertices; }
          MatrixXf &vertices()       { return _vertices; }

    const MatrixXf &normals() const { return _normals; }
          MatrixXf &normals()       { return _normals; }

    const MatrixXf &texcoords() const { return _texcoords; }
          MatrixXf &texcoords()       { return _texcoords; }

    const MatrixXu &triangles() const { return _triangles; }
          MatrixXu &triangles()       { return _triangles; }

    Box3f computeBounds() const;

    void write(std::ostream &os) const;
    void read(std::istream &is);

    static Mesh createBox(const Box3f &box);
    static Mesh createSphere(const Vector3f &position, float radius, int segments = 32);

private:
    void addVertex(const Vector3f &p, const Vector3f &n, const Vector2f &uv);
    void addTriangle(int a, int b, int c);

    MatrixXf _vertices;
    MatrixXf _normals;
    MatrixXf _texcoords;
    MatrixXu _triangles;
};

} // namespace pbs
