#include "thermostat.h"
#include <cmath>

Thermostat::Thermostat(const Input& input)
: temperature(input.temperature), friction(input.friction)
{

}


void Thermostat::apply(System& system, const Input& input, PLMD::Random& random)
{
    double dt = 0.5*input.timestep;
    double c1;
    c1 = exp(-friction*dt);
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        double c2;
        c2 = sqrt((1.0-c1*c1)*Constants::KB*temperature*Constants::EV_AMU_TO_A2_FS2/system.masses[iatom]);
        for (int k=0; k<3; k++)
        {
            system.velocities[iatom][k] = c1*system.velocities[iatom][k]+ c2*random.Gaussian();
        }
    }
}

