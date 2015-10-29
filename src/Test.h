#include "Mesh.h"
#include "VoxelGrid.h"
#include "ObjReader.h"
#include "ObjWriter.h"
#include "MarchingCubes.h"

namespace pbs {

namespace test {

static void obj() {

    Mesh mesh = ObjReader::load("test.obj");
    ObjWriter::save(mesh, "test2.obj");
}

static void marchingCubes() {
    DBG("Marching Cubes ...");

    Vector3i size(100);
    VoxelGridf grid(size);


    for (int z = 0; z < size.z(); ++z) {
        for (int y = 0; y < size.y(); ++y) {
            for (int x = 0; x < size.x(); ++x) {
                Vector3f p((x + 0.5f) / size.x(), (y + 0.5f) / size.y(), (z + 0.5f) / size.z());
                float d = (p - Vector3f(0.5f)).norm() - 0.25f + (std::sin(p.x() * 30.f) * 0.05f + std::cos(p.y() * 50.f) * 0.05f);
                grid(Vector3i(x, y, z)) = d;
            }
        }
    }

    MarchingCubes<float> mc;

    Mesh mesh = mc.generateIsoSurface(grid.data(), 0.f, Box3f(0.f, 1.f), size - Vector3i(1));

    //mc.GenerateSurface(grid.data(), 0.f,  size.x() - 1, size.y() - 1, size.z() - 1, 1.f, 1.f, 1.f);

    ObjWriter::save(mesh, "mc.obj");

}

} // namespace test

} // namespace pbs
