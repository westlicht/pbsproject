#include "Cache.h"

#include "core/Serialize.h"
#include "core/Vector.h"
#include "core/FileUtils.h"

#include <fstream>

namespace pbs {

Cache::Cache(const filesystem::path &dir) :
    _dir(dir),
    _frame(-1),
    _frameCount(0)
{
    FileUtils::createDir(dir.str());
    _valid = readMeta();
}

void Cache::clear() {
    _valid = false;
}

void Cache::commit() {
    writeMeta();
    _valid = true;
}

void Cache::setFrame(int frame) {
    _frame = frame;
}

void Cache::setFrameCount(int frameCount) {
    _frameCount = frameCount;
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

void Cache::writeMesh(const Mesh &mesh) {
    auto os = ostreamFor("mesh", _frame);
    mesh.write(*os);
}

bool Cache::readMesh(Mesh &mesh) {
    auto is = istreamFor("mesh", _frame);
    if (!is) {
        return false;
    }
    mesh.read(*is);
    return true;
}



void Cache::writeMeta() {
    std::ofstream os((_dir / "metadata").str());
    Serialize::write(os, _frameCount);
}

bool Cache::readMeta() {
    std::ifstream is((_dir / "metadata").str());
    if (!is.good()) {
        return false;
    }
    Serialize::read(is, _frameCount);
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
