#include "FDTD.hpp"

namespace FDTD {
    FieldSolver::FieldSolver(std::array<std::size_t, 3> shape) : shape(shape) {
        fields.fill(MDVector<floatType, 3>(shape));
    }

    FieldSolver::getFieldAt = FieldSolver::PEC_BC;
}