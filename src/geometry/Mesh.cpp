#include "Mesh.h"

#include "core/Vector.h"
#include "core/Box.h"
#include "core/serialize.h"

namespace pbs {

Mesh::Mesh() {
    _vertices.resize(3, 0);
    _normals.resize(3, 0);
    _texcoords.resize(3, 0);
    _triangles.resize(3, 0);
}

Box3f Mesh::computeBounds() const {
    Box3f bounds;
    for (int i = 0; i < _vertices.cols(); ++i) {
        bounds.expandBy(_vertices.col(i));
    }
    return bounds;
}

void Mesh::write(std::ostream &os) const {
    Serialize::writeMatrix(os, _vertices);
    Serialize::writeMatrix(os, _normals);
    Serialize::writeMatrix(os, _texcoords);
    Serialize::writeMatrix(os, _triangles);
}

void Mesh::read(std::istream &is) {
    Serialize::readMatrix(is, _vertices);
    Serialize::readMatrix(is, _normals);
    Serialize::readMatrix(is, _texcoords);
    Serialize::readMatrix(is, _triangles);
}

Mesh Mesh::createBox(const Box3f &box) {
    Vector3f center = box.center();
    Vector3f extents = box.extents() * 0.5;
    Mesh mesh;
    int i = 0;
    for (int f = 0; f < 3; ++f) {
        for (int s = -1; s <= 1; s += 2) {
            Vector3f p(center), n(0.f), u(0.f), v(0.f);
            p[f] += extents[f] * s;
            n[f] = s;
            u[(f + 1) % 3] = extents[(f + 1) % 3];
            v[(f + 2) % 3] = extents[(f + 2) % 3];
            mesh.addVertex(p - u - v, n, Vector2f(0.f, 0.f));
            mesh.addVertex(p + u - v, n, Vector2f(1.f, 0.f));
            mesh.addVertex(p - u + v, n, Vector2f(0.f, 1.f));
            mesh.addVertex(p + u + v, n, Vector2f(1.f, 1.f));
            mesh.addTriangle(i, i + 1, i + 2);
            mesh.addTriangle(i + 1, i + 2, i + 3);
            i += 4;
        }
    }
    return mesh;
}

Mesh Mesh::createSphere(const Vector3f &position, float radius, int segments) {
    Mesh mesh;
    for (int itheta = 0; itheta <= segments / 2; ++itheta) {
        for (int iphi = 0; iphi <= segments; ++iphi) {
            float u = float(iphi) / segments;
            float v = float(itheta) / (segments / 2);
            float phi = u * 2.f * M_PI;
            float theta = v * M_PI;
            float sinTheta = std::sin(theta);
            const Vector3f p(
                std::cos(phi) * sinTheta,
                std::sin(phi) * sinTheta,
                std::cos(theta)
            );
            mesh.addVertex(position + p * radius, p, Vector2f(u, v));
            if (itheta < segments / 2 && iphi < segments) {
                mesh.addTriangle(itheta * (segments + 1) + iphi, itheta * (segments + 1) + iphi + 1, (itheta + 1) * (segments + 1) + iphi);
                mesh.addTriangle(itheta * (segments + 1) + iphi + 1, (itheta + 1) * (segments + 1) + iphi, (itheta + 1) * (segments + 1) + iphi + 1);
            }
        }
    }
    return mesh;
}

void Mesh::addVertex(const Vector3f &p, const Vector3f &n, const Vector2f &uv) {
    _vertices.conservativeResize(Eigen::NoChange, _vertices.cols() + 1);
    _vertices.col(_vertices.cols() - 1) = p;
    _normals.conservativeResize(Eigen::NoChange, _normals.cols() + 1);
    _normals.col(_normals.cols() - 1) = n;
    _texcoords.conservativeResize(Eigen::NoChange, _texcoords.cols() + 1);
    _texcoords.col(_texcoords.cols() - 1) = uv;
}

void Mesh::addTriangle(int a, int b, int c) {
    _triangles.conservativeResize(Eigen::NoChange, _triangles.cols() + 1);
    _triangles.col(_triangles.cols() - 1) = Vector3u(a, b, c);
}

} // namespace pbs
