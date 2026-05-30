#include <array>
#include <exception>
#include <iostream>
#include <format>
#include <fstream>

#include "MDVector.hpp"
#include "FDTD.hpp"

namespace PIC {

    //-----FieldSolver definitions-----

    void FieldSolver::compUpdatePEC(std::size_t n, std::size_t i, std::size_t j) {
        floatType test;
        switch (n) {
            case 0:
                if (j == 0 || j == shape[1] - 1) return;
                fields[0][i, j] += step_ratio/permittivity[i, j]*(fields[5][i, j] - fields[5][i, j-1] - space_step * current[0][i, j]);
                return;
            case 1:
                if (i == 0 || i == shape[0] - 1) return;
                fields[1][i, j] += step_ratio/permittivity[i, j]*(-fields[5][i, j] + fields[5][i-1, j] - space_step * current[1][i, j]);
                return;
            case 2:
                if (i == 0 || j == 0 || i == shape[0] - 1 || j == shape[1] - 1) return;
                fields[2][i, j] += step_ratio/permittivity[i, j]*(-fields[3][i, j] + fields[3][i, j-1] + fields[4][i, j] - fields[4][i-1, j] - space_step * current[2][i, j]);
                return;
            case 3:
                fields[3][i, j] += step_ratio/mu0*(-fields[2][i, j+1] + fields[2][i, j]);
                return;
            case 4:
                fields[4][i, j] += step_ratio/mu0*(fields[2][i+1, j] - fields[2][i, j]);
                return;
            case 5:
                fields[5][i, j] += step_ratio/mu0*(fields[0][i, j+1] - fields[0][i, j] - fields[1][i+1, j] + fields[1][i, j]);
                return;
        }
    }

    FieldSolver::FieldSolver(std::array<std::size_t, 2> shape) : shape(shape), permittivity(shape, 1*eps0) {
        fields.fill(MDVector<floatType, 2>(shape));
        current.fill(MDVector<floatType, 2>(shape));
    }
    std::array<MDVector<floatType, 2>, 6> &FieldSolver::getFields() {
        return fields;
    }
    MDVector<floatType, 2> &FieldSolver::getPermittivity() {
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
        outfile << std::format("Shape: ({}, {}), xyz step: {}, time step: {}, time: {}\n",
                               shape[0], shape[1], space_step, space_step*step_ratio, time);
        outfile << "Ex Ey Ez Hx Hy Hz\n";

        for (std::size_t i = 0; i < shape[0]; i++) {
            for (std::size_t j = 0; j < shape[1]; j++) {
                    std::string line;
                    for (std::size_t n = 0; n < 5; n++) {
                        line += std::format("{} ", fields[n][i, j]);
                    }
                    line += std::format("{}\n", fields[5][i, j]);

                    outfile << line;
            }
        }
    }

    //-----ParticleMover definitions-----

    ParticleMover::ParticleMover(const std::size_t particle_count) : positions({ particle_count, 2 }), velocities({ particle_count, 2 }), charges({ particle_count }) {}

    MDVector<floatType, 2> &ParticleMover::getPositions() {
        return positions;
    }
    MDVector<floatType, 2> &ParticleMover::getVelocities() {
        return velocities;
    }
    MDVector<floatType, 1> &ParticleMover::getCharges() {
        return charges;
    }
    void ParticleMover::move(const std::array<MDVector<floatType, 2>, 6>) {
        for (std::size_t i = 0; i < charges.shape[0]; i++) {
            
        }
    }
}