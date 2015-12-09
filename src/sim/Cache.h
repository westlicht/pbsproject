#pragma once

#include "core/Common.h"

#include "geometry/Mesh.h"

#include <filesystem/path.h>

#include <ostream>
#include <istream>
#include <memory>

namespace pbs {

class Cache {
public:
    Cache(const filesystem::path &dir);

    void clear();
    void commit();
    bool valid() const { return _valid; }

    int frame() const { return _frame; }
    void setFrame(int frame);

    int frameCount() const { return _frameCount; }
    void setFrameCount(int frameCount);

    void writeParticles(const std::vector<Vector3f> &particles);
    bool readParticles(std::vector<Vector3f> &particles);

    void writeMesh(const Mesh &mesh);
    bool readMesh(Mesh &mesh);

private:
    void writeMeta();
    bool readMeta();

    filesystem::path pathFor(const std::string &type, int frame) const;
    std::unique_ptr<std::ostream> ostreamFor(const std::string &type, int frame) const;
    std::unique_ptr<std::istream> istreamFor(const std::string &type, int frame) const;

    filesystem::path _dir;
    bool _valid;
    int _frame;
    int _frameCount;

};

} // namespace pbs
