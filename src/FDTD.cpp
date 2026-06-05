#include <algorithm>
#include <array>
#include <cmath>
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
    std::array<MDVector<floatType, 2>, 3> &FieldSolver::getCurrent() {
        return current;
    }
    MDVector<floatType, 2> &FieldSolver::getPermittivity()
    {
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
    const std::array<std::size_t, 2> &FieldSolver::getShape() const {
        return shape;
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
    MDVector<floatType, 1> &ParticleMover::getMasses() {
        return masses;
    }
    const std::size_t &ParticleMover::getParticleCount() {
        return particle_count;
    }

    void ParticleMover::move(const MDVector<floatType, 2> &fields)
    {
        MDVector<floatType, 1> u_minus({ 2 }), u_prime({ 2 }), u_plus({ 2 });
        floatType A; // helper variable

        for (std::size_t i = 0; i < charges.shape[0]; i++) {
            // Update velocities
            A = time_step * charges[i] / 2 * masses[i];

            u_minus[0] = velocities[i, 0] + A * fields[0, i];
            u_minus[1] = velocities[i, 1] + A * fields[1, i];

            A *= fields[5, i];
            
            u_prime[0] = u_minus[0] + A * u_minus[1];
            u_prime[1] = u_minus[1] - A * u_minus[0];

            u_plus[0] = u_minus[0] + 2 * A * u_prime[1] / (1 + A*A);
            u_plus[0] = u_minus[0] - 2 * A * u_prime[0] / (1 + A*A);

            velocities[i, 0] = u_plus[0] + A * fields[0, i];
            velocities[i, 1] = u_plus[1] + A * fields[1, i];

            // Update positions
            positions[i, 0] += velocities[i, 0] * time_step;
            positions[i, 1] += velocities[i, 1] * time_step;
        }
    }

    //-----SimEngine definitions-----

    std::array<floatType, 2> SimEngine::physToIndex(std::array<floatType, 2> point, std::size_t field_comp) {
        static const struct { floatType x, y; } offsets[6] = { // Lookup table for translating physical to index-like coordinates
            { .5, .0 },
            { .0, .5 },
            { .0, .0 },
            { .0, .5 },
            { .5, .0 },
            { .5, .5 },
        };

        floatType i = point[0] / field_sim.getSpaceStep() - offsets[field_comp].x;
        floatType j = point[1] / field_sim.getSpaceStep() - offsets[field_comp].y;

        return { i, j };
    }

    floatType SimEngine::gatherComponent(std::size_t field_comp, std::size_t particle_num)
    {
        floatType x = particle_sim.getPositions()[particle_num, 0], y = particle_sim.getPositions()[particle_num, 1]; // Physical coordinates of chosen particle
        const auto &fields_all = field_sim.getFields();

        auto [i, j] = physToIndex({ x, y }, field_comp);

        std::size_t imin = std::floor(i), imax = std::ceil(i), jmin = std::floor(j), jmax = std::ceil(j);
        if (i < 0) imin = 0; // It's possible for these to be negative for some of the components.
        if (j < 0) jmin = 0; // The calculation itself doesn't need them in those cases, but we have to prevent a segfault nonetheless.
        // They can exceed the max physical boundary as well, but in this case the padding due to the staggering of the Yee grid saves us.

        MDVector<floatType, 2> fields; // Fields to be splined.
        if (field_comp > 3) {
            fields[0, 0] = fields_all[field_comp][imin, jmin];
            fields[1, 0] = fields_all[field_comp][imax, jmin];
            fields[0, 1] = fields_all[field_comp][imin, jmax];
            fields[1, 1] = fields_all[field_comp][imax, jmax];
        }
        else { // Time averaging for H
            fields[0, 0] = .5 * (fields_all[field_comp][imin, jmin] + prev_fields[field_comp][imin, jmin]);
            fields[1, 0] = .5 * (fields_all[field_comp][imax, jmin] + prev_fields[field_comp][imax, jmin]);
            fields[0, 1] = .5 * (fields_all[field_comp][imin, jmax] + prev_fields[field_comp][imin, jmax]);
            fields[1, 1] = .5 * (fields_all[field_comp][imax, jmax] + prev_fields[field_comp][imax, jmax]);
        }

        std::size_t inearest = i - imin < imax - i ? 0 : 1;
        std::size_t jnearest = j - jmin < jmax - j ? 0 : 1;
        i -= imin; // i and j switch to just their fractional parts
        j -= jmin;

        switch (field_comp) { // Splining logic as described in Vay's paper for spline order n=1, energy conserving.
            case 0:
                [[fallthrough]];
            case 4:
                return (1 - j) * fields[inearest, jmin] + j * fields[inearest, jmax];

            case 1:
                [[fallthrough]];
            case 3:
                return (1 - i) * fields[imin, jnearest] + i * fields[imax, jnearest];

            case 2:
                floatType result = 0;
                result += (1 - j) * ((1 - i) * fields[imin, jmin] + i * fields[imax, jmin]);
                result += j * ((1 - i) * fields[imin, jmax] + i * fields[imax, jmax]);
                return result;
            
            case 5:
                return fields[inearest, jnearest];
        }
    }

    MDVector<floatType, 2> SimEngine::fieldGather() {
        MDVector<floatType, 2> result({ particle_sim.getParticleCount(), 6 });

        for (std::size_t particle_num = 0; particle_num < particle_sim.getParticleCount(); particle_num++) {
            for (std::size_t field_comp = 0; field_comp < 6; field_comp++) {
                result[particle_num, field_comp] = gatherComponent(particle_num, field_comp);
            }
        }

        return result;
    }

    void SimEngine::depositCurrent() {
        auto &currents = field_sim.getCurrent();
        const auto &positions = particle_sim.getPositions();
        const auto &velocities = particle_sim.getVelocities();
        const auto &charges = particle_sim.getCharges();
        const auto &step = field_sim.getSpaceStep();

        for (std::size_t particle_num = 0; particle_num < particle_sim.getParticleCount(); particle_num++) {
            for (std::size_t field_comp = 0; field_comp < 3; field_comp++) {
                floatType x = .5 * (positions[particle_num, 0] + prev_positions[particle_num, 0]);
                floatType y = .5 * (positions[particle_num, 1] + prev_positions[particle_num, 1]);
                auto [i, j] = physToIndex({ x, y }, field_comp);

                std::size_t imin = std::floor(i), imax = std::ceil(i), jmin = std::floor(j), jmax = std::ceil(j);

                if (i > 0) {
                    currents[field_comp][imin, jmax] = (1 - i) * j * charges[particle_num] * velocities[particle_num, field_comp] / step / step / step;
                }
                if (j > 0) {
                    currents[field_comp][imax, jmin] = i * (1 - j) * charges[particle_num] * velocities[particle_num, field_comp] / step / step / step;
                }
                if (i > 0 && j > 0) {
                    currents[field_comp][imin, jmin] = (1 - i) * (1 - j) * charges[particle_num] * velocities[particle_num, field_comp] / step / step / step;
                }
                currents[field_comp][imax, jmax] = i * j * charges[particle_num] * velocities[particle_num, field_comp] / step / step / step;
            }
        }
    }
}