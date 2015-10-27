#pragma once

#include "Vector.h"

#include <vector>

namespace pbs {

class Mesh {
public:
    struct Triangle {
        uint32_t i1, i2, i3;
    };

    const std::vector<Vector3f> &vertices() const { return _vertices; }
          std::vector<Vector3f> &vertices()       { return _vertices; }

    const std::vector<Vector3f> &normals() const { return _normals; }
          std::vector<Vector3f> &normals()       { return _normals; }

    const std::vector<Vector2f> &texcoords() const { return _texcoords; }
          std::vector<Vector2f> &texcoords()       { return _texcoords; }

    const std::vector<Triangle> &triangles() const { return _triangles; }
          std::vector<Triangle> &triangles()       { return _triangles; }

private:
    std::vector<Vector3f> _vertices;
    std::vector<Vector3f> _normals;
    std::vector<Vector2f> _texcoords;
    std::vector<Triangle> _triangles;
};

} // namespace pbs
