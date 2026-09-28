#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "system.h"
#include "input.h"

class Intergrator
{
private:
    double timestep;
public:
    Intergrator(const Input& input);
    void first_step(System& system);
    void second_step(System& system);
    void update_positions(System& system, const Input& input);
};
#endif
