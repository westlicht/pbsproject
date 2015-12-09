#pragma once

#include "core/Common.h"

#include <filesystem/path.h>

#include <ostream>
#include <istream>
#include <memory>

namespace pbs {

class Cache {
public:
    Cache(const filesystem::path &dir);

    void clear();

    void setFrame(int frame);

    int startFrame() const;
    int endFrame() const;

    void writeParticles(const std::vector<Vector3f> &particles);
    bool readParticles(std::vector<Vector3f> &particles);

private:
    filesystem::path pathFor(const std::string &type, int frame) const;
    std::unique_ptr<std::ostream> ostreamFor(const std::string &type, int frame) const;
    std::unique_ptr<std::istream> istreamFor(const std::string &type, int frame) const;

    filesystem::path _dir;
    int _frame;

};

} // namespace pbs
