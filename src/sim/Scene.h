#pragma once

#include "core/Common.h"
#include "core/Vector.h"
#include "core/Box.h"
#include "core/Properties.h"

#include <filesystem/path.h>

#include <string>
#include <vector>

namespace pbs {

class Scene {
public:
    enum Type {
        Fluid,
        Boundary,
    };

    struct World {
        Box3f bounds;
    };

    struct Shape {
        Type type;
        Shape(const Properties &props);
    };

    struct Box : public Shape {
        Box3f bounds;
        Box(const Properties &props);
    };
    struct Sphere : public Shape {
        Vector3f position;
        float radius;
        Sphere(const Properties &props);
    };
    struct Mesh : public Shape {
        std::string filename;
        Mesh(const Properties &props);
    };

    Properties settings;

    World world;
    std::vector<Box> boxes;
    std::vector<Sphere> spheres;
    std::vector<Mesh> meshes;

    Scene();

    std::string dump() const;

    static Scene load(const std::string &filename);

private:
    static Type typeFromString(const std::string &name);
    static std::string resolvePath(const std::string &path);
    static filesystem::resolver _resolver;

};


} // namespace pbs
