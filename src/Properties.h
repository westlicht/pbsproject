#pragma once

#include "Common.h"

#include <json11.h>

#include <string>

namespace pbs {

class Properties {
public:
    Properties() = default;
    Properties(const json11::Json &json);

    const json11::Json &json() const { return _json; }

    bool hasString(const std::string &name) const;
    std::string getString(const std::string &name) const;
    std::string getString(const std::string &name, const std::string &def) const;

    bool hasBool(const std::string &name) const;
    bool getBool(const std::string &name) const;
    bool getBool(const std::string &name, bool def) const;

    bool hasFloat(const std::string &name) const;
    float getFloat(const std::string &name) const;
    float getFloat(const std::string &name, float def) const;

    bool hasInteger(const std::string &name) const;
    int getInteger(const std::string &name) const;
    int getInteger(const std::string &name, int def) const;

    bool hasVector2(const std::string &name) const;
    Vector2f getVector2(const std::string &name) const;
    Vector2f getVector2(const std::string &name, const Vector2f &def) const;

    bool hasVector3(const std::string &name) const;
    Vector3f getVector3(const std::string &name) const;
    Vector3f getVector3(const std::string &name, const Vector3f &def) const;

    bool hasBox2(const std::string &name) const;
    Box2f getBox2(const std::string &name) const;
    Box2f getBox2(const std::string &name, const Box2f &def) const;

    bool hasBox3(const std::string &name) const;
    Box3f getBox3(const std::string &name) const;
    Box3f getBox3(const std::string &name, const Box3f &def) const;

private:
    json11::Json _json;
};

} // namespace pbs
