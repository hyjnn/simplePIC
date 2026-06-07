/*
    !!! WARNING !!!
    For testing purposes I decicded to make this 2D only for the time being.
    In order for this to make sense with particles inside the sim region,
    we have to only consider the Ex Ey Hz polarization!!

    Things to improve:
    Syncing time_step etc between simulation elements
    FieldSolever shape is not being updated when changed with getFields, this should be redesigned. (e.g. by adding a resize and setters, not allowing the setters to change shape implicitly, or not allowing resizing at all)
    Stop mixing std::arrays and std::vectors of MDVectors - I'll need to add multi-dimensional pushes to MDVector for this to make sense.
    FieldSolver stores currents as a member, meanwhile ParticleMover expects field as an argument to move.
    Make the linalg in ParticleMover::kickMove prettier
*/

#pragma once
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "MDVector.hpp"

namespace PIC {
    using floatType = float; //float type used in calculations
    inline constexpr floatType c = 299792458;
    inline constexpr floatType eps0 = 8.8541878188141414141e-12;
    inline constexpr floatType mu0 = 1.25663706127202020202e-6;

    class FieldSolver {
        //Dimensions of the simulation region.
        std::array<std::size_t, 2> shape;
        floatType space_step = 1;
        floatType step_ratio = 1./(c*std::sqrt(2.)); // time_step / space_step
        floatType time = 0; // Time at which the E field is known, the H field is known half a time step forward.
        // Array storing the Ex, Ey, Ez, Hx, Hy, Hz fields in this order.
        // IMPORTANT: Due to the Yee grid's structure, certain edge elements of the fields lie outside of the simulation region!
        std::array<MDVector<floatType, 2>, 6> fields;
        MDVector<floatType, 2> permittivity; // Absolute permittivity
        std::array<MDVector<floatType, 2>, 3> current;

        //Hides scary loop behind function (i.e. calls func for every element in every element in fields)
        template<typename T>
        void forAllFields(T(FieldSolver::*func)(std::size_t, std::size_t, std::size_t));
        //Updates n-th EM field component at index (i, j, k).
        void compUpdatePEC(std::size_t n, std::size_t i, std::size_t j);
        
    public:
        FieldSolver(std::array<std::size_t, 2>); // Constructs a FieldSolver object with specified shape.

        std::array<MDVector<floatType, 2>, 6> &getFields(); // Currently, changing the shape of fields with this WILL BREAK THE PROGRAM. TO BE CORRECTED!!!!
        std::array<MDVector<floatType, 2>, 3> &getCurrent();
        MDVector<floatType, 2> &getPermittivity();
        floatType getSpaceStep();
        floatType setSpaceStep(floatType);
        floatType getStepRatio();
        floatType setStepRatio(floatType);
        const std::array<std::size_t, 2> &getShape() const;

        void init_PEC(); // Sets relevant boundary components to 0.
        void solve(unsigned long long n_steps);

        void exportToFile(std::string);
    };

    class ParticleMover {
        MDVector<floatType, 2> positions, velocities;
        MDVector<floatType, 1> charges, masses;
        floatType time_step = 1./(c*std::sqrt(2.));

    public:
        ParticleMover(const std::size_t); // Creates object with specified particle count.

        MDVector<floatType, 2>& getPositions();
        MDVector<floatType, 2>& getVelocities();
        MDVector<floatType, 1>& getCharges();
        MDVector<floatType, 1>& getMasses();
        const floatType &getTimeStep();
        const std::size_t &getParticleCount();

        void move(); // Move particles using current velocity
        void kickMove(const MDVector<floatType, 2>&); // Move particles and update velocities in accordance to the supplied fields.
    };

    class SimEngine {
        /*
            This class mixes two coordinate systems - the physical coordinates, which are used by ParticleMover to store the particle coordinates,
            and index-like coordinates, which are really a set of 6 different systems, corresponding to the 6 field components. These come from treating
            the indices of the points on which the given field is known as their coordinates, and extending this linearly to the entire simulation
            region. So, the index-like coordinates are rescaled and shifted physical coordinates.
        */

        FieldSolver field_sim;
        ParticleMover particle_sim;
        floatType time = 0;
        floatType time_step = 1./(c*std::sqrt(2.));

        std::array<MDVector<floatType, 2>, 6> prev_fields;
        MDVector<floatType, 2> prev_positions;
        std::vector<std::size_t> tracked_particles; // Stores indicies of tracked particles.
        std::vector<MDVector<floatType, 2>> tracked_positions; // Stores position history of tracked particles.
        std::vector<MDVector<floatType, 2>> tracked_velocities; // Stores velocity history of tracked particles.

        std::array<floatType, 2> physToIndex(std::array<floatType, 2> point, std::size_t field_comp); // Convert point to index-like coordinates for specified field component.
        void updateTracked();

        floatType gatherComponent(std::size_t field_comp, std::size_t particle_num); // Gathers specified field component onto specified index-like point.
        MDVector<floatType, 2> fieldGather(); // Calculates and returns fields at current particle locations according to the energy conserving scheme as described in Vay.
        void depositCurrent(); // Updates the current distribution in field_sim.

    public:
        SimEngine(std::array<std::size_t, 2> shape, std::size_t particle_num); // Contructs Yee simulation region with given shape and a Boris particle mover with given number of particles.

        FieldSolver &getFieldSim();
        ParticleMover &getParticleSim();

        void trackParticle(std::size_t); // Turns on tracking for particle with chosen index.
        void exportTracked(std::string); // Saves tracked particle positions to file.

        void initialize(); // Should be called before calling run for the first time.
        void run(unsigned long long n_steps); // Requires the current E field to be at time i, B and v at time i + 1/2, and x at i + 1.
    };

    // Applies func to all EM fields' components at all yee grid points, starting with the E field.
    template<typename T>
    void FieldSolver::forAllFields(T(FieldSolver::*func)(std::size_t, std::size_t, std::size_t)) {
        for (std::size_t n = 0; n < 6; n++) {
            std::array<int, 3> skip; // because of the spacial staggering on a Yee grid, we need to ignore some of the elements in fields (they lie outside the sim region)
            if (n < 3) {
                skip = { 0, 0, 0 };
                skip[n] = 1;
            }
            else {
                skip = { 1, 1, 1 };
                skip[n-3] = 0;
            }
            for (std::size_t i = 0; i < shape[0] - skip[0]; i++) {
                for (std::size_t j = 0; j < shape[1] - skip[1]; j++) {
                    (this->*func)(n, i, j);
                }
            }
        }
    }
}