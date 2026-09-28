#include "integrator.h"
#include "constants.h"

Integrator::Integrator(const Input& input): timestep(input.timestep)
{

}


void Integrator::first_step(System& system)
{
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        for (int k=0; k<3; k++)
        {
            system.velocities[iatom][k] += 0.5*timestep*system.forces[iatom][k]
            /system.masses[iatom]*Constants::EV_AMU_TO_A2_FS2;
        }
    }
}


void Integrator::update_positions(System& system, const Input& input)
{
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        for (int k=0; k<3; k++)
        {
            system.positions[iatom][k] += system.velocities[iatom][k]*input.timestep;
        }
    }
}


void Integrator::second_step(System& system)
{
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        for (int k=0; k<3; k++)
        {
            system.velocities[iatom][k]+=0.5*timestep*system.forces[iatom][k]
            /system.masses[iatom]*Constants::EV_AMU_TO_A2_FS2;
        }

    }
}
