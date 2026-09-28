#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "system.h"
#include "input.h"

class Integrator
{
private:
    double timestep;
public:
    Integrator(const Input& input);
    void first_step(System& system);
    void second_step(System& system);
    void update_positions(System& system, const Input& input);
};
#endif
