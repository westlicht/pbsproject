#include "ObjWriter.h"
#include "Mesh.h"

#include <fstream>

namespace pbs {

void ObjWriter::save(const Mesh &mesh, const std::string &filename) {

    std::ofstream os(filename);
    if (os.fail()) {
        throw Exception("Unable to open OBJ file \"%s\" for writing!", filename);
    }

    for (int i = 0; i < mesh.vertices().cols(); ++i) {
        const auto &v = mesh.vertices().col(i);
        os << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    }

    for (int i = 0; i < mesh.normals().cols(); ++i) {
        const auto &n = mesh.normals().col(i);
        os << "vn " << n.x() << " " << n.y() << " " << n.z() << std::endl;
    }

    for (int i = 0; i < mesh.texcoords().cols(); ++i) {
        const auto &tc = mesh.texcoords().col(i);
        os << "vn " << tc.x() << " " << tc.y() << std::endl;
    }

    os << std::endl;

    bool writeNormals = mesh.normals().cols() > 0;
    bool writeTexcoords = mesh.texcoords().cols() > 0;

    for (int i = 0; i < mesh.triangles().cols(); ++i) {
        const auto &t = mesh.triangles().col(i);
        os << "f ";
        for (int j = 0; j < 3; ++j) {
            os << (t[j] + 1) << "/";
            if (writeTexcoords) os << (t[j] + 1);
            os << "/";
            if (writeNormals) os << (t[j] + 1);
            if (j < 2) os << " ";
        }
        os << std::endl;
    }
}

} // namespace pbs
