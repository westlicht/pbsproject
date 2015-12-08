#pragma once

#include "core/Common.h"
#include "core/Vector.h"
#include "core/Box.h"
#include "core/Properties.h"

#include <filesystem/path.h>

#include <string>
#include <map>
#include <vector>

namespace pbs {

class Scene {
public:
    enum Type {
        Fluid,
        Boundary,
    };

    struct Camera {
        Vector3f position = Vector3f(0.f, 0.f, 5.f);
        Vector3f target = Vector3f(0.f);
        Vector3f up = Vector3f(0.f, 1.f, 0.f);
        float fov = 30.f;
        float near = 0.1f;
        float far = 100.f;
        Camera() = default;
        Camera(const Properties &props);
        std::string toString() const;
    };

    struct World {
        Box3f bounds = Box3f(Vector3f(-1.f), Vector3f(1.f));
        World() = default;
        World(const Properties &props);
        std::string toString() const;
    };

    struct Shape {
        Type type;
        Shape(const Properties &props);
    };

    struct Box : public Shape {
        Box3f bounds;
        Box(const Properties &props);
        std::string toString() const;
    };

    struct Sphere : public Shape {
        Vector3f position;
        float radius;
        Sphere(const Properties &props);
        std::string toString() const;
    };

    struct Mesh : public Shape {
        std::string filename;
        Mesh(const Properties &props);
        std::string toString() const;
    };

    Properties settings;

    Camera camera;
    World world;
    std::vector<Box> boxes;
    std::vector<Sphere> spheres;
    std::vector<Mesh> meshes;

    static Scene load(const std::string &filename, const json11::Json &settings = json11::Json());

    std::string toString() const;

private:
    static Type typeFromString(const std::string &name);
    static std::string typeToString(Type type);
    static std::string resolvePath(const std::string &path);
    static filesystem::resolver _resolver;

};


} // namespace pbs
