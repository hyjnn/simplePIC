#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <format>
#include <numbers>
#include "FDTD.hpp"
#include "MDVector.hpp"

void fundamentalMode(PIC::FieldSolver &sim) { // Shape is hardcoded to 400 cause I'm lazy and this is just for testing
    auto step = sim.getSpaceStep();
    auto time_step = step * sim.getStepRatio();
    auto length = step * (400 - 1);
    auto frequency = 2 * std::numbers::pi * PIC::c / (2 * length);
    auto& fields = sim.getFields();

    PIC::floatType phase, Ey0, Hz0;
    for (std::size_t i = 0; i < 400; i++) {
        phase = std::numbers::pi * i * step/length;
        Ey0 = -PIC::c * PIC::mu0 * std::sin(phase) * 100;
        Hz0 = std::cos(phase) * std::cos(frequency*time_step) * 100;
        for (std::size_t j = 0; j < 400; j++) {
            fields[1][i, j] = Ey0;
            fields[5][i, j] = Hz0;
        }
    }
}

int main() {
    PIC::SimEngine sim({ 400, 400 }, 1 );
    auto &field_sim = sim.getFieldSim();
    auto &particle_sim = sim.getParticleSim();
    auto &positions = particle_sim.getPositions();
    // fundamentalMode(field_sim);
    auto &Ey = field_sim.getFields()[1].data;
    auto &Hz = field_sim.getFields()[5].data;
    std::ranges::fill(Ey, -40000);
    std::ranges::fill(Hz, 100);

    positions[0, 0] = 200 * field_sim.getSpaceStep();
    positions[0, 1] = 200 * field_sim.getSpaceStep();
    particle_sim.getMasses()[0] = 1;
    particle_sim.getCharges()[0] = 1e10;
    sim.trackParticle(0);

    sim.initialize();

    field_sim.exportToFile(std::format("field_data{}.txt", 0));
    for (int i = 1; i < 5; i++) {
        sim.run(20);
        field_sim.exportToFile(std::format("field_data{}.txt", i));
    }
    sim.exportTracked("particle_data.txt");
}