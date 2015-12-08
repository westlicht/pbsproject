#include "Cache.h"

namespace pbs {

Cache::Cache(const std::string &dir) :
    _dir(dir)
{}

void Cache::clear() {

}

int Cache::startFrame() const {
    return 0;
}

int Cache::endFrame() const {
    return 0;
}

void Cache::writeParticles(int frame, const std::vector<Vector3f> &particles) {

}

bool Cache::readParticles(int frame, std::vector<Vector3f> &particles) {
    return false;
}

} // namespace pbs
