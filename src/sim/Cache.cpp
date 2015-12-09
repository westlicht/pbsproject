#include "Cache.h"

#include "core/Serialize.h"
#include "core/Vector.h"

#include <fstream>

namespace pbs {

Cache::Cache(const filesystem::path &dir) :
    _dir(dir),
    _frame(-1)
{
    DBG("cacheDir: %s", dir.str());
}

void Cache::clear() {

}

void Cache::setFrame(int frame) {
    _frame = frame;
}

int Cache::startFrame() const {
    return 0;
}

int Cache::endFrame() const {
    return 0;
}

void Cache::writeParticles(const std::vector<Vector3f> &particles) {
    auto os = ostreamFor("particles", _frame);
    Serialize::writeVector(*os, particles);
}

bool Cache::readParticles(std::vector<Vector3f> &particles) {
    auto is = istreamFor("particles", _frame);
    if (!is) {
        return false;
    }
    Serialize::readVector(*is, particles);
    return true;
}

filesystem::path Cache::pathFor(const std::string &type, int frame) const {
    return _dir / tfm::format("%s-%04d.dat", type, frame);
}

std::unique_ptr<std::ostream> Cache::ostreamFor(const std::string &type, int frame) const {
    return std::move(std::unique_ptr<std::ostream>(new std::ofstream(pathFor(type, frame).str().c_str())));
}

std::unique_ptr<std::istream> Cache::istreamFor(const std::string &type, int frame) const {
    return std::move(std::unique_ptr<std::istream>(new std::ifstream(pathFor(type, frame).str().c_str())));
}

} // namespace pbs
