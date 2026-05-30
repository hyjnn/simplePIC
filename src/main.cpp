#include <array>
#include <cmath>
#include <iostream>
#include <format>
#include <numbers>
#include "FDTD.hpp"
#include "MDVector.hpp"

int main() {
    PIC::FieldSolver sim(std::array<std::size_t, 2>{ 400, 400 });

    auto step = sim.getSpaceStep();
    auto time_step = step * sim.getStepRatio();
    auto length = step * (400 - 1);
    auto frequency = 2 * std::numbers::pi * PIC::c / (2 * length);
    auto& fields = sim.getFields();

    PIC::floatType phase, Ey0, Hz0;
    for (std::size_t i = 0; i < 400; i++) {
        phase = std::numbers::pi * i * step / length;
        Ey0 = -PIC::c*PIC::mu0*std::sin(phase);
        Hz0 = std::cos(phase)*std::cos(frequency * time_step);
        for (std::size_t j = 0; j < 400; j++) {
            fields[1][i, j] = Ey0;
            fields[5][i, j] = Hz0;
        }
    }

    sim.exportToFile(std::format("data{}.txt", 0));
    for (int i = 1; i < 5; i++) {
        sim.solve(20);
        sim.exportToFile(std::format("data{}.txt", i));
    }
}