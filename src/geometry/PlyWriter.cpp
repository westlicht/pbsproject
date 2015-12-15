#include "PlyWriter.h"
#include "Mesh.h"

#include "core/Serialize.h"
#include "core/Vector.h"

#include <fstream>

namespace pbs {

void PlyWriter::save(const Mesh &mesh, const std::string &filename) {

    std::ofstream os(filename);
    if (os.fail()) {
        throw Exception("Unable to open PLY file \"%s\" for writing!", filename);
    }

    bool writeNormals = mesh.normals().cols() > 0;
    bool writeTexcoords = mesh.texcoords().cols() > 0;
    writeTexcoords = false;

    os << "ply" << std::endl;
    os << "format binary_little_endian 1.0" << std::endl;

    os << "element vertex " << mesh.vertices().cols() << std::endl;
    os << "property float x" << std::endl;
    os << "property float y" << std::endl;
    os << "property float z" << std::endl;
    if (writeNormals) {
        os << "property float nx" << std::endl;
        os << "property float ny" << std::endl;
        os << "property float nz" << std::endl;
    }
    if (writeTexcoords) {
        os << "property float s" << std::endl;
        os << "property float t" << std::endl;
    }
    os << "element face " << mesh.triangles().cols() << std::endl;
    os << "property list uchar uint vertex_indices" << std::endl;
    os << "end_header" << std::endl;

    for (int i = 0; i < mesh.vertices().cols(); ++i) {
        Vector3f v = mesh.vertices().col(i);
        Serialize::write(os, v);
        if (writeNormals) {
            Vector3f n = mesh.normals().col(i);
            Serialize::write(os, n);
        }
        if (writeTexcoords) {
            Vector2f tc = mesh.texcoords().col(i);
            Serialize::write(os, tc);
        }
    }

    for (int i = 0; i < mesh.triangles().cols(); ++i) {
        Vector3u t = mesh.triangles().col(i);
        uint8_t n = 3;
        Serialize::write(os, n);
        Serialize::write(os, t);
    }
}

} // namespace pbs
