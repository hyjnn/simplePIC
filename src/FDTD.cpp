#include <array>
#include <exception>
#include <iostream>
#include <format>
#include <fstream>

#include "MDVector.hpp"
#include "FDTD.hpp"

namespace FDTD {
    void FieldSolver::compUpdatePEC(std::size_t n, std::size_t i, std::size_t j, std::size_t k) {
        floatType test;
        switch (n) {
            case 0:
                if (j == 0 || j == shape[1] || k == 0 || k == shape[2]) return;
                fields[0][i, j, k] += step_ratio/permittivity[i, j, k]*(fields[4][i, j, k] - fields[4][i, j, k-1] + fields[5][i, j, k] - fields[5][i, j-1, k]);
                return;
            case 1:
                if (i == 0 || i == shape[0] || k == 0 || k == shape[2]) return;
                fields[1][i, j, k] += step_ratio/permittivity[i, j, k]*(fields[3][i, j, k] - fields[3][i, j, k-1] + fields[5][i, j, k] - fields[5][i-1, j, k]);
                return;
            case 2:
                if (i == 0 || i == shape[0] || j == 0 || j == shape[1]) return;
                fields[2][i, j, k] += step_ratio/permittivity[i, j, k]*(fields[3][i, j, k] - fields[3][i, j-1, k] + fields[4][i, j, k] - fields[4][i-1, j, k]);
                return;
            case 3:
                fields[3][i, j, k] += step_ratio*(fields[1][i, j, k+1] - fields[1][i, j, k] + fields[2][i, j+1, k] - fields[2][i, j, k]);
                return;
            case 4:
                test = step_ratio*(fields[0][i, j, k+1] - fields[0][i, j, k] + fields[2][i+1, j, k] - fields[2][i, j, k]);
                if (test != 0) {
                    std::cout << std::format("Hy[{0}, {1}, {2}] = {3}, Hy[{0}, {1}, {2}-1] = {4}, Hz[{0}, {1}, {2}] = {5}, Hz[{0}, {1}-1, {2}] = {6}\n", i, j, k, fields[0][i, j, k+1], fields[0][i, j, k], fields[2][i+1, j, k], fields[2][i, j, k]) <<
                                 std::format("permittivity[{0}, {1}, {2}] = {3}, i = {0}, j = {1}, k = {2}\n", i, j, k, permittivity[i, j, k]);
                    throw std::runtime_error("wtopa");
                }
                fields[4][i, j, k] += step_ratio*(fields[0][i, j, k+1] - fields[0][i, j, k] + fields[2][i+1, j, k] - fields[2][i, j, k]);
                return;
            case 5:
                fields[5][i, j, k] += step_ratio*(fields[0][i, j+1, k] - fields[0][i, j, k] + fields[1][i+1, j, k] - fields[1][i, j, k]);
                return;
        }
    }

    FieldSolver::FieldSolver(std::array<std::size_t, 3> shape) : shape(shape), permittivity(MDVector<floatType, 3>(shape, 1)) {
        fields.fill(MDVector<floatType, 3>(shape));
    }
    std::array<MDVector<floatType, 3>, 6> &FieldSolver::getFields() {
        return fields;
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

    void FieldSolver::solve(unsigned long long n_steps) {
        for (unsigned long long step = 0; step < n_steps; step++) {
            forAllFields(&FieldSolver::compUpdatePEC);
        }
        time += space_step * step_ratio * n_steps;
    }

    void FieldSolver::exportToFile(std::string filename) {
        std::ofstream outfile(filename, std::ios_base::trunc);
        outfile << std::format("Shape: ({}, {}, {}), xyz step: {}, time step: {}, time: {}\n",
                               shape[0], shape[1], shape[2], space_step, space_step*step_ratio, time);
        outfile << "Ex, Ey, Ez, Hx, Hy, Hz\n";

        for (std::size_t i = 0; i < shape[0]; i++) {
            for (std::size_t j = 0; j < shape[1]; j++) {
                for (std::size_t k = 0; k < shape[2]; k++) {
                    std::string line;
                    for (std::size_t n = 0; n < 5; n++) {
                        line += std::format("{}, ", fields[n][i, j, k]);
                    }
                    line += std::format("{}\n", fields[5][i, j, k]);

                    outfile << line;
                }
            }
        }
    }
}