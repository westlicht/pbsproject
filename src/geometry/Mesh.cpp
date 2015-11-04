#include "Mesh.h"

#include "core/Box.h"

namespace pbs {

Box3f Mesh::computeBounds() const {
    Box3f bounds;
    for (int i = 0; i < _vertices.cols(); ++i) {
        bounds.expandBy(_vertices.col(i));
    }
    return bounds;
}

} // namespace pbs
