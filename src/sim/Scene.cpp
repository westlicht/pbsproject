#include "Scene.h"

#include "core/Common.h"
#include "core/Vector.h"
#include "core/Box.h"

#include <filesystem/path.h>
#include <filesystem/resolver.h>
#include <tinyformat.h>
#include <json11.h>

#include <fstream>
#include <sstream>

using namespace json11;

namespace pbs {

filesystem::resolver Scene::_resolver;

Scene::Shape::Shape(const Properties &props) {
    type = typeFromString(props.getString("type", "fluid"));
}

Scene::Box::Box(const Properties &props) : Shape(props) {
    bounds = props.getBox3("bounds");
}

Scene::Sphere::Sphere(const Properties &props) : Shape(props) {
    position = props.getVector3("position");
    radius = props.getFloat("radius");
}

Scene::Mesh::Mesh(const Properties &props) : Shape(props) {
    filename = resolvePath(props.getString("filename"));
}

Scene::Scene() {
    world.bounds = Box3f(Vector3f(-1.f), Vector3f(1.f));
}

std::string Scene::dump() const {
    std::stringstream ss;
    ss << world.bounds;

    std::string result;
    result += "Scene[\n";
    result += tfm::format("  settings = %s\n", settings.json().dump());
    result += "  world = [\n";
    result += tfm::format("    bounds = %s\n", world.bounds);
    result += "  ],\n";
    result += "  boxes = [\n";
    for (const auto &box : boxes) {
        result += tfm::format("    Box[bounds=%s],\n", box.bounds);
    }
    result += "  ],\n";
    result += "  spheres = [\n";
    for (const auto &sphere : spheres) {
        result += tfm::format("    Sphere[position=%s, radius=%f],\n", sphere.position, sphere.radius);
    }
    result += "  ],\n";
    result += "  meshes = [\n";
    for (const auto &mesh : meshes) {
        result += tfm::format("    Mesh[filename=%s],\n", mesh.filename);
    }
    result += "  ],\n";
    result += "]";
    return result;
}

Scene Scene::load(const std::string &filename) {
    std::ifstream is(filename);
    std::stringstream ss;
    ss << is.rdbuf();

    std::string err;
    Json jsonRoot = Json::parse(ss.str(), err);
    if (jsonRoot.is_null()) {
        throw Exception("Failed to load scene from '%s' (error: %s)", filename, err);
    }

    _resolver = filesystem::resolver();
    _resolver.prepend(filesystem::path(filename).parent_path());

    Scene scene;
    scene.settings = Properties(jsonRoot["settings"]);

    // Parse scene objects
    Json jsonScene = jsonRoot["scene"];
    if (jsonScene.is_object()) {
        {
            auto jsonWorld = jsonScene["world"];
            Properties props(jsonWorld);
            scene.world.bounds = props.getBox3("bounds");
        }
        for (auto jsonBox : jsonScene["boxes"].array_items()) {
            scene.boxes.emplace_back(Box(Properties(jsonBox)));
        }
        for (auto jsonSphere : jsonScene["spheres"].array_items()) {
            scene.spheres.emplace_back(Sphere(Properties(jsonSphere)));
        }
        for (auto jsonMesh : jsonScene["meshes"].array_items()) {
            scene.meshes.emplace_back(Mesh(Properties(jsonMesh)));
        }
    }

    return scene;
}

Scene::Type Scene::typeFromString(const std::string &name) {
    if (name == "fluid") {
        return Fluid;
    } else if (name == "boundary") {
        return Boundary;
    } else {
        throw Exception("Unknown type name '%s'", name);
    }
}

std::string Scene::resolvePath(const std::string &path) {
    return _resolver.resolve(path).str();
}

} // namespace pbs
