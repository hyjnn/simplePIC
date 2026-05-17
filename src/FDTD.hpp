#pragma once
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "MDVector.hpp"

namespace FDTD {
    using floatType = float; //float type used in calculations

    class FieldSolver {
        //Dimensions of the simulation region.
        std::array<std::size_t, 3> shape;
        floatType space_step = 1;
        floatType step_ratio = 1/std::sqrt(3.0L); // time_step / space_step
        floatType time = 0;
        //Array storing the Ex, Ey, Ez, Hx, Hy, Hz fields in this order.
        std::array<MDVector<floatType, 3>, 6> fields;
        MDVector<floatType, 3> permittivity; // Relative permittivity

        //Hides scary loop behind function (i.e. calls func for every element in every element in fields)
        template<typename T>
        void forAllFields(T(FieldSolver::*func)(std::size_t, std::size_t, std::size_t, std::size_t));
        //Updates n-th EM field component at index (i, j, k).
        void compUpdatePEC(std::size_t, std::size_t, std::size_t, std::size_t);
        
    public:
        FieldSolver(std::array<std::size_t, 3>); //Constructs a FieldSolver object with specified shape.

        const std::array<MDVector<floatType, 3>, 6> &getFields() const;
        std::array<MDVector<floatType, 3>, 6> &getFields();
        const MDVector<floatType, 3> &getPermittivity() const;
        MDVector<floatType, 3> &getPermittivity();
        floatType getSpaceStep();
        floatType setSpaceStep(floatType);
        floatType getStepRatio();
        floatType setStepRatio(floatType);

        void updateFields();

        void fieldsToStream(std::ostream&);
    };

    template<typename T>
    void FieldSolver::forAllFields(T(FieldSolver::*func)(std::size_t, std::size_t, std::size_t, std::size_t)) {
        for (std::size_t n = 0; n < 6; n++) {
            for (std::size_t i = 0; i < shape[0]; i++) {
                for (std::size_t j = 0; j < shape[1]; j++) {
                    for (std::size_t k = 0; k < shape[2]; k++) {
                        func(n, i, j, k);
                    }
                }
            }
        }
    }
}