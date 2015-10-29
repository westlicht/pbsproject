#pragma once

#include "Common.h"
#include "Box.h"
#include "Vector.h"

namespace pbs {

class Mesh;

// Marching cubes mesh generation.
template<typename T>
class MarchingCubes {
public:
    MarchingCubes();
    
    Mesh generateIsoSurface(const T *scalarField, T isoLevel, const Box3f &bounds, const Vector3i &cells);

private:
    unsigned int edgeId(unsigned int x, unsigned int y, unsigned int z, unsigned int edgeNo);
    unsigned int vertexId(unsigned int x, unsigned int y, unsigned int z);

    Vector3f calculateIntersection(unsigned int x, unsigned int y, unsigned int z, unsigned int edgeNo);
    Vector3f interpolate(const Vector3f &p1, const Vector3f &p2, T val1, T val2);

    Vector3i _cells;
    Vector3f _origin;
    Vector3f _cellSize;

    const T *_scalarField;
    T _isoLevel;

    static const int _edgeTable[256];
    static const int _triTable[256][16];
};

} // namespace pbs
