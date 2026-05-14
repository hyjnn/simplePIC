#pragma once
#include <array>
#include <vector>

#include "MDVector.hpp"

namespace FDTD {
    using floatType = float; //float type used in calculations

    class FieldSolver {
        //Dimensions of the simulation region.
        std::array<std::size_t, 3> shape;
        //Array storing the Ex, Ey, Ez, Hx, Hy, Hz fields in this order.
        std::array<MDVector<floatType, 3>, 6> fields;

        floatType &PEC_BC(std::size_t, std::size_t, std::size_t, std::size_t);
        

    public:
        FieldSolver(std::array<std::size_t, 3>);
        //Returns n-th EM field component at index (i, j, k).
        floatType &getFieldAt(std::size_t n, std::size_t i, std::size_t j, std::size_t k);
        void updateFields();
    };
}