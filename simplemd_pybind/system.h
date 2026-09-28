#ifndef SYSTEM_H
#define SYSTEM_H

#include "Vector.h"
#include "constants.h"
#include <vector>

class System
{
public:
    System(int n);

    int natoms;

    std::vector<PLMD::Vector> positions;

    std::vector<PLMD::Vector> velocities;
    std::vector<PLMD::Vector> forces;

    double potential_energy;
    double kinetic_energy;

    double cell[3];

    std::vector<int> types;
    std::vector<int> atomic_numbers;
    std::vector<double> masses;

    void compute_kinetic_energy();
    PLMD::Vector compute_momentum();
    void zero_forces();
};

#endif