#include "ObjWriter.h"
#include "Mesh.h"

#include <fstream>

namespace pbs {

void ObjWriter::save(const Mesh &mesh, const std::string &filename) {

    std::ofstream os(filename);
    if (os.fail()) {
        throw Exception("Unable to open OBJ file \"%s\" for writing!", filename);
    }

    for (const auto &v : mesh.vertices()) {
        os << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    }

    for (const auto &n : mesh.normals()) {
        os << "vn " << n.x() << " " << n.y() << " " << n.z() << std::endl;
    }

    for (const auto &tc : mesh.texcoords()) {
        os << "vt " << tc.x() << " " << tc.y() << std::endl;
    }

    os << std::endl;

    bool writeNormals = !mesh.normals().empty();
    bool writeTexcoords = !mesh.texcoords().empty();

    for (const auto &t : mesh.triangles()) {
        uint32_t i[3] = { t.i1 + 1, t.i2 + 1, t.i3 + 1 };
        os << "f ";
        for (int j = 0; j < 3; ++j) {
            os << i[j] << "/";
            if (writeTexcoords) os << i[j];
            os << "/";
            if (writeNormals) os << i[j];
            if (j < 2) os << " ";
        }
        os << std::endl;
    }
}

} // namespace pbs
