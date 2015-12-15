#pragma once

#include <string>

namespace pbs {

class Mesh;

// PLY writer.
class PlyWriter {
public:
    static void save(const Mesh &mesh, const std::string &filename);
};

} // namespace pbs
