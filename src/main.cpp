#include <array>
#include <iostream>
#include "FDTD.hpp"
#include "MDVector.hpp"

int main() {
    FDTD::FieldSolver sim(std::array<std::size_t, 3>{ 100, 100, 100 });

    sim.solve(1);

    sim.exportToFile("data.txt");
}