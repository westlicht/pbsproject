#pragma once

#include <string>

namespace pbs {

class Mesh;

class ObjWriter {
public:
    static void save(const Mesh &mesh, const std::string &filename);
};

} // namespace pbs
