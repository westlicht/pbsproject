#include "ObjReader.h"
#include "Mesh.h"

#include "core/Vector.h"

#include <fstream>
#include <unordered_map>

namespace pbs {

Mesh ObjReader::load(const std::string &filename) {
    Mesh mesh;

    std::vector<Vector3f> positions;
    std::vector<Vector3f> normals;
    std::vector<Vector2f> texcoords;
    std::vector<uint32_t> indices;
    std::vector<OBJVertex> vertices;
    std::unordered_map<OBJVertex, uint32_t, OBJVertexHash> vertexMap;

    std::ifstream is(filename);
    if (is.fail()) {
        throw Exception("Unable to open OBJ file \"%s\" for reading!", filename);
    }

    std::string line_str;
    while (std::getline(is, line_str)) {
        std::istringstream line(line_str);

        std::string prefix;
        line >> prefix;

        if (prefix == "v") {
            Vector3f p;
            line >> p.x() >> p.y() >> p.z();
            //m_bbox.expandBy(p);
            positions.emplace_back(p);
        } else if (prefix == "vt") {
            Vector2f tc;
            line >> tc.x() >> tc.y();
            texcoords.emplace_back(tc);
        } else if (prefix == "vn") {
            Vector3f n;
            line >> n.x() >> n.y() >> n.z();
            normals.emplace_back(n.normalized());
        } else if (prefix == "f") {
            std::string v1, v2, v3, v4;
            line >> v1 >> v2 >> v3 >> v4;
            OBJVertex verts[6];
            int nVertices = 3;

            verts[0] = OBJVertex(v1);
            verts[1] = OBJVertex(v2);
            verts[2] = OBJVertex(v3);

            if (!v4.empty()) {
                /* This is a quad, split into two triangles */
                verts[3] = OBJVertex(v4);
                verts[4] = verts[0];
                verts[5] = verts[2];
                nVertices = 6;
            }
            /* Convert to an indexed vertex list */
            for (int i=0; i<nVertices; ++i) {
                const OBJVertex &v = verts[i];
                auto it = vertexMap.find(v);
                if (it == vertexMap.end()) {
                    vertexMap[v] = (uint32_t) vertices.size();
                    indices.emplace_back((uint32_t) vertices.size());
                    vertices.emplace_back(v);
                } else {
                    indices.emplace_back(it->second);
                }
            }
        }
    }

    mesh.triangles().resize(3, indices.size() / 3);
    std::memcpy(mesh.triangles().data(), indices.data(), indices.size() * sizeof(uint32_t));

    mesh.vertices().resize(3, vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        mesh.vertices().col(i) = positions[vertices[i].p - 1];
    }

    if (!normals.empty()) {
        mesh.normals().resize(3, vertices.size());
        for (size_t i = 0; i < vertices.size(); ++i) {
            mesh.normals().col(i) = normals[vertices[i].n - 1];
        }
    }

    if (!texcoords.empty()) {
        mesh.texcoords().resize(2, vertices.size());
        for (size_t i = 0; i < vertices.size(); ++i) {
            mesh.texcoords().col(i) = texcoords[vertices[i].uv - 1];
        }
    }

    return std::move(mesh);
}

} // namespace pbs
