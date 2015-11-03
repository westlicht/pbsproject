#pragma once

#include "core/Vector.h"

#include <vector>

namespace pbs {

// Holds a triangle mesh.
// Data is stored in columns of the matrices.
class Mesh {
public:
    const MatrixXf &vertices() const { return _vertices; }
          MatrixXf &vertices()       { return _vertices; }

    const MatrixXf &normals() const { return _normals; }
          MatrixXf &normals()       { return _normals; }

    const MatrixXf &texcoords() const { return _texcoords; }
          MatrixXf &texcoords()       { return _texcoords; }

    const MatrixXu &triangles() const { return _triangles; }
          MatrixXu &triangles()       { return _triangles; }
          
private:
    MatrixXf _vertices;
    MatrixXf _normals;
    MatrixXf _texcoords;
    MatrixXu _triangles;
};

} // namespace pbs
