#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <exception>
#include <format>
#include <fstream>
#include <utility>

#include "MDVector.hpp"
#include "FDTD.hpp"

namespace PIC {

    //-----FieldSolver definitions-----

    void FieldSolver::compUpdatePEC(const std::size_t n, const std::size_t i, const std::size_t j) {
        floatType test;
        switch (n) {
            case 0:
                if (j == 0 || j == shape[1] - 1) return;
                fields[0][i, j] += step_ratio/permittivity[i, j]*(fields[5][i, j] - fields[5][i, j-1] - space_step * current[0][i, j]);
                break;
            case 1:
                if (i == 0 || i == shape[0] - 1) return;
                fields[1][i, j] += step_ratio/permittivity[i, j]*(-fields[5][i, j] + fields[5][i-1, j] - space_step * current[1][i, j]);
                break;
            case 2:
                if (i == 0 || j == 0 || i == shape[0] - 1 || j == shape[1] - 1) return;
                fields[2][i, j] += step_ratio/permittivity[i, j]*(-fields[3][i, j] + fields[3][i, j-1] + fields[4][i, j] - fields[4][i-1, j] - space_step * current[2][i, j]);
                break;
            case 3:
                fields[3][i, j] += step_ratio/mu0*(-fields[2][i, j+1] + fields[2][i, j]);
                break;
            case 4:
                fields[4][i, j] += step_ratio/mu0*(fields[2][i+1, j] - fields[2][i, j]);
                break;
            case 5:
                fields[5][i, j] += step_ratio/mu0*(fields[0][i, j+1] - fields[0][i, j] - fields[1][i+1, j] + fields[1][i, j]);
                break;
        }
    }

    FieldSolver::FieldSolver(const std::array<std::size_t, 2> shape) : shape(shape), permittivity(shape, 1*eps0) {
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
    const floatType &FieldSolver::getSpaceStep() const{
        return space_step;
    }
    const floatType &FieldSolver::setSpaceStep(const floatType new_value) {
        space_step = new_value;
        return space_step;
    }
    const floatType &FieldSolver::getStepRatio() const{
        return step_ratio;
    }
    const floatType &FieldSolver::setStepRatio(const floatType new_value) {
        step_ratio = new_value;
        return step_ratio;
    }
    const std::array<std::size_t, 2> &FieldSolver::getShape() const {
        return shape;
    }

    void FieldSolver::init_PEC() {
        for (std::size_t i = 0; i < shape[0]; i++) { // At boundaries perpendicular to y
            fields[0][i, 0] = 0; // Ex=0
            fields[0][i, shape[0] - 1] = 0;
            fields[2][i, 0] = 0; // Ez=0
            fields[2][i, shape[0] - 1] = 0;
            fields[4][i, 0] = 0; // Hy=0
            fields[4][i, shape[0] - 1] = 0;
        }
        for (std::size_t j = 0; j < shape[0]; j++) { // At boundaries perpendicular to x
            fields[1][0, j] = 0; // Ey=0
            fields[1][shape[0] - 1, j] = 0;
            fields[2][0, j] = 0; // Ez=0
            fields[2][shape[0] - 1, j] = 0;
            fields[3][0, j] = 0; // Hx=0
            fields[3][shape[0] - 1, j] = 0;
        }
    }

    void FieldSolver::solve(unsigned long long n_steps)
    {
        for (unsigned long long step = 0; step < n_steps; step++) {
            forAllFields(&FieldSolver::compUpdatePEC);
        }
        time += space_step * step_ratio * n_steps;
    }

    void FieldSolver::exportToFile(const std::string filename) {
        std::ofstream outfile(filename, std::ios_base::trunc);
        outfile << std::format("shape: ({}, {}), xyz step: {}, time step: {}, time: {}\n",
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

    ParticleMover::ParticleMover(const std::size_t particle_count) : positions({ particle_count, 2 }), velocities({ particle_count, 2 }),
                                                                     charges({ particle_count }), masses({ particle_count }) {}

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
    const floatType &ParticleMover::getTimeStep() const {
        return time_step;
    }
    const floatType &ParticleMover::setTimeStep(const floatType val) {
        time_step = val;
        return time_step;
    }

    const std::size_t &ParticleMover::getParticleCount() const {
        return charges.shape[0];
    }

    void ParticleMover::removeOutside(floatType xmax, floatType ymax) {
        for (size_t i = 0; i < getParticleCount(); i++) {
            if (positions[i, 0] < 0 || positions[i, 0] > xmax || positions[i, 1] < 0 || positions[i, 1] > ymax) {
                positions.erase_row(i);
                velocities.erase_row(i);
                charges.erase_row(i);
                masses.erase_row(i);
            }
        }
    }

    void ParticleMover::move()
    {
        for (std::size_t i = 0; i < charges.shape[0]; i++) {
            positions[i, 0] += velocities[i, 0] * time_step;
            positions[i, 1] += velocities[i, 1] * time_step;
        }
    }

    void ParticleMover::kickMove(const MDVector<floatType, 2> &fields)
    {
        MDVector<floatType, 1> u_minus({ 2 }), u_prime({ 2 }), u_plus({ 2 });
        floatType A, B, u_natural; // helper variables

        for (std::size_t i = 0; i < charges.shape[0]; i++) {
            // Update velocities
            A = time_step * charges[i] / (2 * masses[i]);

            u_minus[0] = velocities[i, 0] + A * fields[0, i];
            u_minus[1] = velocities[i, 1] + A * fields[1, i];

            u_natural = u_minus[0]*u_minus[0] / c / c + u_minus[1]*u_minus[1] / c / c; 
            B = A * fields[5, i] * mu0 / std::sqrt(1 + u_natural);
            
            u_prime[0] = u_minus[0] + B * u_minus[1];
            u_prime[1] = u_minus[1] - B * u_minus[0];

            u_plus[0] = u_minus[0] + 2 * B * u_prime[1] / (1 + B*B);
            u_plus[1] = u_minus[1] - 2 * B * u_prime[0] / (1 + B*B);

            velocities[i, 0] = u_plus[0] + A * fields[0, i];
            velocities[i, 1] = u_plus[1] + A * fields[1, i];

            // Update positions
            positions[i, 0] += time_step / std::sqrt(1 / (c*c) + 1 / (velocities[i, 0] * velocities[i, 0]));
            positions[i, 1] += time_step / std::sqrt(1 / (c*c) + 1 / (velocities[i, 1] * velocities[i, 1]));
        }
    }

    //-----SimEngine definitions-----

    std::array<floatType, 2> SimEngine::physToIndex(const std::array<floatType, 2> point, const std::size_t field_comp) {
        static const struct { floatType x, y; } offsets[6] = { // Lookup table for translating physical to index-like coordinates
            { .5, .0 },
            { .0, .5 },
            { .0, .0 },
            { .0, .5 },
            { .5, .0 },
            { .5, .5 },
        };

        floatType i = point[0] / space_step - offsets[field_comp].x;
        floatType j = point[1] / space_step - offsets[field_comp].y;

        return { i, j };
    }

    void SimEngine::syncSteps() {
        field_sim.setSpaceStep(space_step);
        field_sim.setStepRatio(time_step / space_step);
        particle_sim.setTimeStep(time_step);
    }

    void SimEngine::updateTracked() {
        std::size_t last_row;
        for (std::size_t i = 0; i < particle_sim.getCharges().shape[0]; i++) {
            last_row = tracked_positions[i].shape[0];
            tracked_positions[i].append_rows(1);
            tracked_velocities[i].append_rows(1);

            tracked_positions[i][last_row, 0] = particle_sim.getPositions()[i, 0];
            tracked_positions[i][last_row, 1] = particle_sim.getPositions()[i, 1];
            tracked_velocities[i][last_row, 0] = particle_sim.getVelocities()[i, 0];
            tracked_velocities[i][last_row, 1] = particle_sim.getVelocities()[i, 1];
        }
    }

    floatType SimEngine::gatherComponent(const std::size_t field_comp, const std::size_t particle_num)
    {
        floatType x = particle_sim.getPositions()[particle_num, 0], y = particle_sim.getPositions()[particle_num, 1]; // Physical coordinates of chosen particle
        const auto &fields_all = field_sim.getFields();

        auto [i, j] = physToIndex({ x, y }, field_comp);

        std::size_t imin = std::floor(i), imax = imin + 1, jmin = std::floor(j), jmax = jmin + 1;
        if (i < 0) imin = 0; // It's possible for these to be negative for some of the components.
        if (j < 0) jmin = 0; // The calculation itself doesn't need them in those cases, but we have to prevent a segfault nonetheless.
        // They can exceed the max physical boundary as well, but in this case the padding due to the staggering of the Yee grid saves us.

        MDVector<floatType, 2> fields({ 2, 2 }); // Fields to be splined.
        if (field_comp < 3) {
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

        floatType result = 0; // Helper for longer formulas
        switch (field_comp) { // Splining logic as described in Vay's paper for spline order n=1, energy conserving.
            case 0:
                [[fallthrough]];
            case 4:
                return (1 - j) * fields[inearest, 0] + j * fields[inearest, 1];

            case 1:
                [[fallthrough]];
            case 3:
                return (1 - i) * fields[0, jnearest] + i * fields[1, jnearest];

            case 2:
                result += (1 - j) * ((1 - i) * fields[0, 0] + i * fields[1, 0]);
                result += j * ((1 - i) * fields[0, 1] + i * fields[1, 1]);
                return result;
            
            case 5:
                return fields[inearest, jnearest];
            default:
                assert(false && "invalid dimension");
                std::unreachable();
        }
    }

    MDVector<floatType, 2> SimEngine::fieldGather() {
        MDVector<floatType, 2> result({ 6, particle_sim.getParticleCount()});

        for (std::size_t particle_num = 0; particle_num < particle_sim.getParticleCount(); particle_num++) {
            for (std::size_t field_comp = 0; field_comp < 6; field_comp++) {
                result[field_comp, particle_num] = gatherComponent(field_comp, particle_num);
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

        currents.fill(MDVector<floatType, 2>(field_sim.getShape())); // clear previous currents

        for (std::size_t particle_num = 0; particle_num < particle_sim.getParticleCount(); particle_num++) {
            for (std::size_t field_comp = 0; field_comp < 2; field_comp++) {
                floatType x = .5 * (positions[particle_num, 0] + prev_positions[particle_num, 0]);
                floatType y = .5 * (positions[particle_num, 1] + prev_positions[particle_num, 1]);
                auto [i, j] = physToIndex({ x, y }, field_comp);

                std::size_t imin = std::floor(i), imax = imin + 1, jmin = std::floor(j), jmax = jmin + 1;

                floatType test = charges[particle_num] * velocities[particle_num, field_comp] / step / step;

                if (i >= 0) {
                    currents[field_comp][imin, jmax] += (1 - i) * j * test;
                }
                if (j >= 0) {
                    currents[field_comp][imax, jmin] += i * (1 - j) * charges[particle_num] * velocities[particle_num, field_comp] / step / step;
                }
                if (i >= 0 && j >= 0) {
                    currents[field_comp][imin, jmin] += (1 - i) * (1 - j) * charges[particle_num] * velocities[particle_num, field_comp] / step / step;
                }
                currents[field_comp][imax, jmax] += i * j * charges[particle_num] * velocities[particle_num, field_comp] / step / step;
            }
        }
    }

    SimEngine::SimEngine(const std::array<std::size_t, 2> shape, const std::size_t particle_num) : field_sim(shape), particle_sim(particle_num), prev_positions({ particle_num, 2 }),
                                                                                                   tracked_positions(particle_num, MDVector<floatType, 2>({ particle_num, 2 })),
                                                                                                   tracked_velocities(particle_num, MDVector<floatType, 2>({ particle_num, 2 })) {
        prev_fields.fill(MDVector<floatType, 2>(shape));
    }
    SimEngine::SimEngine(std::array<std::size_t, 2> shape, std::size_t particle_num, floatType space_step) : SimEngine(shape, particle_num) {
        this->space_step = space_step;
        this->time_step = space_step / (c*std::sqrt(2));
        syncSteps();
    }
    SimEngine::SimEngine(std::array<std::size_t, 2> shape, std::size_t particle_num, floatType space_step, floatType time_step) : SimEngine(shape, particle_num) {
        this->space_step = space_step;
        this->time_step = time_step;
        syncSteps();
    }

    FieldSolver &SimEngine::getFieldSim() {
        return field_sim;
    }
    ParticleMover &SimEngine::getParticleSim() {
        return particle_sim;
    }

    void SimEngine::exportTracked(std::string filename) {
        std::ofstream outfile(filename, std::ios_base::trunc);
        outfile << std::format("time step: {}, particle count: {}\n", particle_sim.getTimeStep(), particle_sim.getParticleCount());

        for (std::size_t particle = 0; particle < particle_sim.getParticleCount(); particle++) {
            outfile << std::format("particle {}:\n", particle + 1);
            for (std::size_t i = 0; i < tracked_positions[particle].shape[0]; i++) {
                outfile << std::format("{} {} {} {}\n", tracked_positions[particle][i, 0], tracked_positions[particle][i, 1],
                                                        1/std::sqrt(1/(c*c) + 1/(tracked_velocities[particle][i, 0]*tracked_velocities[particle][i, 0])),
                                                        1/std::sqrt(1/(c*c) + 1/(tracked_velocities[particle][i, 1]*tracked_velocities[particle][i, 1])));
            }
        }
    }

    void SimEngine::initialize() {
        prev_positions = particle_sim.getPositions();
        particle_sim.move();
        updateTracked();

        field_sim.init_PEC();
    }

    void SimEngine::run(const unsigned long long n_steps)
    {
        for (std::size_t i = 0; i < n_steps; i++) {
            depositCurrent();
            prev_fields = field_sim.getFields();
            field_sim.solve(1);
            prev_positions = particle_sim.getPositions();
            particle_sim.kickMove(fieldGather());
            particle_sim.removeOutside(space_step * (field_sim.getShape()[0] - 1), space_step * (field_sim.getShape()[1] - 1));
            updateTracked();
            time += time_step;
        }
    }
}