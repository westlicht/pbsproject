#pragma once

#include "core/Common.h"
#include "core/Vector.h"
#include "core/Box.h"
#include "core/Properties.h"

#include <string>
#include <vector>

namespace pbs {

class Scene {
public:
    struct World {
        Box3f bounds;
    };
    struct Box {
        Box3f bounds;
    };
    struct Sphere {
        Vector3f position;
        float radius;
    };

    Properties settings;

    World world;
    std::vector<Box> boxes;
    std::vector<Sphere> spheres;

    Scene();

    std::string dump() const;

    static Scene load(const std::string &filename);

private:
};


} // namespace pbs
