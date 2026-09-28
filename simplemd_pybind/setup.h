#ifndef SETUP_H
#define SETUP_H

#include "system.h"
#include "constants.h"
#include "input.h"
#include "Vector.h"
#include "Random.h"

class Setup
{
public:
    void initialize(System& system, const Input& input, PLMD::Random& random);
private:
    void randomize_velocities(System& system, PLMD::Random& random, const Input& input);
    void remove_center_of_mass_velocity(System& system);
    void rescale_temperature(const Input& input, System& system);
};

#endif