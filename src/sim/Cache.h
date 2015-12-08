#pragma once

#include "core/Common.h"

namespace pbs {

class Cache {
public:
    Cache(const std::string &dir);

    void clear();

    int startFrame() const;
    int endFrame() const;

    void writeParticles(int frame, const std::vector<Vector3f> &particles);
    bool readParticles(int frame, std::vector<Vector3f> &particles);



private:
    std::string _dir;
};

} // namespace pbs
