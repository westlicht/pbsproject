#include "Properties.h"
#include "Vector.h"
#include "Box.h"

using namespace json11;

namespace pbs {

static bool isVector2(const Json &json) {
    return json.is_array() && json.array_items().size() == 2 && json[0].is_number() && json[1].is_number();
}

static bool isVector3(const Json &json) {
    return json.is_array() && json.array_items().size() == 3 && json[0].is_number() && json[1].is_number() && json[2].is_number();
}

static bool isBox2(const Json &json) {
    return json.is_array() && json.array_items().size() == 2 && isVector2(json[0]) && isVector2(json[1]);
}

static bool isBox3(const Json &json) {
    return json.is_array() && json.array_items().size() == 2 && isVector3(json[0]) && isVector3(json[1]);
}

static Vector2f parseVector2(const Json &json) {
    return Vector2f(json[0].number_value(), json[1].number_value());
}

static Vector3f parseVector3(const Json &json) {
    return Vector3f(json[0].number_value(), json[1].number_value(), json[2].number_value());
}

static Box2f parseBox2(const Json &json) {
    return Box2f(parseVector2(json[0]), parseVector2(json[1]));
}

static Box3f parseBox3(const Json &json) {
    return Box3f(parseVector3(json[0]), parseVector3(json[1]));
}


Properties::Properties(const Json &json) :
    _json(json)
{}

bool Properties::hasObject(const std::string &name) const {
    return _json[name].is_object();
}

Properties Properties::getObject(const std::string &name) const {
    if (!hasObject(name)) throw Exception("No object named '%s'", name);
    return Properties(_json[name]);
}

bool Properties::hasString(const std::string &name) const {
    return _json[name].is_string();
}

std::string Properties::getString(const std::string &name) const {
    if (!hasString(name)) throw Exception("No string property named '%s'", name);
    return _json[name].string_value();
}

std::string Properties::getString(const std::string &name, const std::string &def) const {
    return hasString(name) ? _json[name].string_value() : def;
}

bool Properties::hasBool(const std::string &name) const {
    return _json[name].is_bool();
}

bool Properties::getBool(const std::string &name) const {
    if (!hasBool(name)) throw Exception("No bool property named '%s'", name);
    return _json[name].bool_value();
}

bool Properties::getBool(const std::string &name, bool def) const {
    return hasBool(name) ? _json[name].bool_value() : def;
}

bool Properties::hasFloat(const std::string &name) const {
    return _json[name].is_number();
}

float Properties::getFloat(const std::string &name) const {
    if (!hasFloat(name)) throw Exception("No float property named '%s'", name);
    return _json[name].number_value();
}

float Properties::getFloat(const std::string &name, float def) const {
    return hasFloat(name) ? _json[name].number_value() : def;
}

bool Properties::hasInteger(const std::string &name) const {
    return _json[name].is_number();
}

int Properties::getInteger(const std::string &name) const {
    if (!hasInteger(name)) throw Exception("No integer property named '%s'", name);
    return _json[name].int_value();
}

int Properties::getInteger(const std::string &name, int def) const {
    return hasInteger(name) ? _json[name].int_value() : def;
}

bool Properties::hasVector2(const std::string &name) const {
    return isVector2(_json[name]);
}

Vector2f Properties::getVector2(const std::string &name) const {
    if (!hasVector2(name)) throw Exception("No Vector2 property named '%s'", name);
    return parseVector2(_json[name]);
}

Vector2f Properties::getVector2(const std::string &name, const Vector2f &def) const {
    return hasVector2(name) ? parseVector2(_json[name]) : def;
}

bool Properties::hasVector3(const std::string &name) const {
    return isVector3(_json[name]);
}

Vector3f Properties::getVector3(const std::string &name) const {
    if (!hasVector3(name)) throw Exception("No Vector3 property named '%s'", name);
    return parseVector3(_json[name]);
}

Vector3f Properties::getVector3(const std::string &name, const Vector3f &def) const {
    return hasVector3(name) ? parseVector3(_json[name]) : def;
}

bool Properties::hasBox2(const std::string &name) const {
    return isBox2(_json[name]);
}

Box2f Properties::getBox2(const std::string &name) const {
    if (!hasBox2(name)) throw Exception("No Box2 property named '%s'", name);
    return parseBox2(_json[name]);
}

Box2f Properties::getBox2(const std::string &name, const Box2f &def) const {
    return hasBox2(name) ? parseBox2(_json[name]) : def;
}

bool Properties::hasBox3(const std::string &name) const {
    return isBox3(_json[name]);
}

Box3f Properties::getBox3(const std::string &name) const {
    if (!hasBox3(name)) throw Exception("No Box3 property named '%s'", name);
    return parseBox3(_json[name]);
}

Box3f Properties::getBox3(const std::string &name, const Box3f &def) const {
    return hasBox3(name) ? parseBox3(_json[name]) : def;
}

} // namespace pbs
