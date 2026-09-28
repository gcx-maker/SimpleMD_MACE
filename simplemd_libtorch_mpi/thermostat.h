#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "system.h"
#include "Random.h"
#include "constants.h"
#include "input.h"

class Thermostat
{
private:
    double temperature;
    double friction;

public:
    Thermostat(const Input& input);
    void apply(System& system, const Input& input, PLMD::Random& random);
};

#endif