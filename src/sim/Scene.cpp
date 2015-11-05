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

    Scene scene;

    filesystem::resolver resolver;
    resolver.prepend(filesystem::path(filename).parent_path());
    auto resolvePath = [&resolver] (const std::string &path) {
        return resolver.resolve(path).str();
    };

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
            Properties props(jsonBox);
            scene.boxes.emplace_back(Box({
                props.getBox3("bounds")
            }));
        }
        for (auto jsonSphere : jsonScene["spheres"].array_items()) {
            Properties props(jsonSphere);
            scene.spheres.emplace_back(Sphere({
                props.getVector3("position"),
                props.getFloat("radius")
            }));
        }
        for (auto jsonSphere : jsonScene["meshes"].array_items()) {
            Properties props(jsonSphere);
            scene.meshes.emplace_back(Mesh({
                resolvePath(props.getString("filename"))
            }));
        }
    }

    return scene;
}

} // namespace pbs
