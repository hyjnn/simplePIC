#include <array>
#include <iostream>
#include <format>
#include "FDTD.hpp"
#include "MDVector.hpp"

int main() {
    FDTD::FieldSolver sim(std::array<std::size_t, 2>{ 400, 400 });

    auto& fields = sim.getFields();
    fields[0][200, 200] = 1;

    sim.exportToFile(std::format("data{}.txt", 0));
    for (int i = 1; i < 5; i++) {
        sim.solve(100);
        sim.exportToFile(std::format("data{}.txt", i));
    }
}