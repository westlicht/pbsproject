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

Scene::Camera::Camera(const Properties &props) {
    position = props.getVector3("position", position);
    target = props.getVector3("target", target);
    up = props.getVector3("up", up);
    fov = props.getFloat("fov", fov);
    near = props.getFloat("near", near);
    far = props.getFloat("far", far);
    frame = props.getInteger("frame", frame);
}

std::string Scene::Camera::toString() const {
    return tfm::format(
        "Camera[\n"
        "  position = %s,\n"
        "  target = %s,\n"
        "  up = %s,\n"
        "  fov = %f,\n"
        "  near = %f,\n"
        "  far = %f\n"
        "  frame = %d\n"
        "]",
        position, target, up, fov, near, far, frame
    );
}

Scene::World::World(const Properties &props) {
    bounds = props.getBox3("bounds", bounds);
}

std::string Scene::World::toString() const {
    return tfm::format(
        "World[\n"
        "  bounds = %s\n"
        "]",
        bounds
    );
}

Scene::Shape::Shape(const Properties &props) {
    type = typeFromString(props.getString("type", "fluid"));
}

Scene::Box::Box(const Properties &props) : Shape(props) {
    bounds = props.getBox3("bounds");
}

std::string Scene::Box::toString() const {
    return tfm::format(
        "Box[\n"
        "  type = %s,\n"
        "  bounds = %s\n"
        "]",
        typeToString(type), bounds
    );
}

Scene::Sphere::Sphere(const Properties &props) : Shape(props) {
    position = props.getVector3("position");
    radius = props.getFloat("radius");
}

std::string Scene::Sphere::toString() const {
    return tfm::format(
        "Sphere[\n"
        "  type = %s,\n"
        "  position = %s,\n"
        "  radius = %f\n"
        "]",
        typeToString(type), position, radius
    );
}

Scene::Mesh::Mesh(const Properties &props) : Shape(props) {
    filename = resolvePath(props.getString("filename"));
}

std::string Scene::Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  type = %s,\n"
        "  filename = %s\n"
        "]",
        typeToString(type), filename
    );
}

Scene Scene::load(const std::string &filename, const json11::Json &settings) {
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

    // Patch settings
    auto settingsValues = jsonRoot["settings"].object_items();
    for (auto kv : settings.object_items()) {
        settingsValues[kv.first] = kv.second;
    }
    scene.settings = Properties(json11::Json(settingsValues));

    // Parse scene objects
    Json jsonScene = jsonRoot["scene"];
    if (jsonScene.is_object()) {
        auto jsonCamera = jsonScene["camera"];
        if (jsonCamera.is_object()) {
            Properties props(jsonCamera);
            scene.camera = Camera(jsonCamera);
        }
        auto jsonWorld = jsonScene["world"];
        if (jsonWorld.is_object()) {
            Properties props(jsonWorld);
            scene.world = World(jsonWorld);
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
        for (auto jsonCameraKeyframe : jsonScene["cameraKeyframes"].array_items()) {
            scene.cameraKeyframes.emplace_back(Camera(Properties(jsonCameraKeyframe)));
        }
        // Set default camera
        if (!jsonCamera.is_object()) {
            Vector3f center = scene.world.bounds.center();
            scene.camera.position += center;
            scene.camera.target += center;
        }
    }

    return scene;
}

template<typename T>
static std::string vectorToString(const std::vector<T> &items) {
    std::string result = "[";
    for (const auto &item : items) {
        result += "\n  " + indent(item.toString());
    }
    result += items.empty() ? "]" : "\n]";
    return result;
}

std::string Scene::toString() const {
    return tfm::format(
        "Scene[\n"
        "  settings = %s,\n"
        "  camera = %s,\n"
        "  world = %s,\n"
        "  boxes = %s,\n"
        "  spheres = %s,\n"
        "  meshes = %s\n,"
        "  cameraKeyframes = %s\n"
        "]",
        indent(settings.json().dump()),
        indent(camera.toString()),
        indent(world.toString()),
        indent(vectorToString(boxes)),
        indent(vectorToString(spheres)),
        indent(vectorToString(meshes)),
        indent(vectorToString(cameraKeyframes))
    );
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

std::string Scene::typeToString(Type type) {
    switch (type) {
    case Fluid: return "fluid";
    case Boundary: return "boundary";
    }
    return "unknown";
}

std::string Scene::resolvePath(const std::string &path) {
    return _resolver.resolve(path).str();
}

} // namespace pbs
