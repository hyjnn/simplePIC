#include "FDTD.hpp"


namespace FDTD {
    void FieldSolver::compUpdatePEC(std::size_t n, std::size_t i, std::size_t j, std::size_t k) {
        switch (n) {
            case 0:
                if (j == 0 || j == shape[1] || k == 0 || k == shape[2]) return;
                fields[0][i, j, k] += step_ratio/permittivity[i, j, k]*(fields[4][i, j, k] - fields[4][i, j, k-1] + fields[5][i, j, k] - fields[5][i, j-1, k]);
            case 1:
                if (i == 0 || i == shape[0] || k == 0 || k == shape[2]) return;
                fields[1][i, j, k] += step_ratio/permittivity[i, j, k]*(fields[3][i, j, k] - fields[3][i, j, k-1] + fields[5][i, j, k] - fields[5][i-1, j, k]);
            case 2:
                if (i == 0 || i == shape[0] || j == 0 || j == shape[1]) return;
                fields[2][i, j, k] += step_ratio/permittivity[i, j, k]*(fields[3][i, j, k] - fields[3][i, j-1, k] + fields[4][i, j, k] - fields[4][i-1, j, k]);
            case 3:
                fields[3][i, j, k] += step_ratio*(fields[1][i, j, k+1] - fields[1][i, j, k] + fields[2][i, j+1, k] - fields[2][i, j, k]);
            case 4:
                fields[4][i, j, k] += step_ratio*(fields[0][i, j, k+1] - fields[0][i, j, k] + fields[2][i+1, j, k] - fields[2][i, j, k]);
            case 5:
                fields[5][i, j, k] += step_ratio*(fields[0][i, j+1, k] - fields[0][i, j, k] + fields[1][i+1, j, k] - fields[1][i, j, k]);
        }
    }

    FieldSolver::FieldSolver(std::array<std::size_t, 3> shape) : shape(shape), permittivity(MDVector<floatType, 3>(shape)) {
        fields.fill(MDVector<floatType, 3>(shape));
    }
    const std::array<MDVector<floatType, 3>, 6> &FieldSolver::getFields() const {
        return fields;
    }
    std::array<MDVector<floatType, 3>, 6> &FieldSolver::getFields() {
        return fields;
    }
    const MDVector<floatType, 3> &FieldSolver::getPermittivity() const {
        return permittivity;
    }
    MDVector<floatType, 3> &FieldSolver::getPermittivity() {
        return permittivity;
    }
    floatType FieldSolver::getSpaceStep() {
        return space_step;
    }
    floatType FieldSolver::setSpaceStep(floatType new_value) {
        space_step = new_value;
        return space_step;
    }
    floatType FieldSolver::getStepRatio() {
        return step_ratio;
    }

    floatType FieldSolver::setStepRatio(floatType new_value) {
        step_ratio = new_value;
        return step_ratio;
    }

    //I know this looks stupid but I probably will add more BCs later
    void FieldSolver::updateFields() {
        forAllFields(compUpdatePEC);
    }

    void FieldSolver::fieldsToStream(std::ostream &) {
        
    }
}