#include "system.h"

System::System(int n): natoms(n), potential_energy(0.0), kinetic_energy(0.0)
{
    types.resize(n);
    atomic_numbers.resize(n);
    masses.resize(n);

    positions.resize(n);
    velocities.resize(n);
    forces.resize(n);

    for (int i=0; i<3; i++)
    {
        cell[i]=0;
    }
}


void System::compute_kinetic_energy()
{
    kinetic_energy = 0.0;
    for (int iatom=0; iatom<natoms; iatom++)
    {
        for (int i=0; i<3; i++)
        {
            kinetic_energy += 0.5*masses[iatom]
            *velocities[iatom][i]*velocities[iatom][i]
            /Constants::EV_AMU_TO_A2_FS2;
        }
    }
}


PLMD::Vector System::compute_momentum()
{
    PLMD::Vector momentum = {0.0, 0.0, 0.0};
    for (int iatom=0; iatom<natoms; iatom++)
    {
        for (int i=0; i<3; i++)
        {
            momentum[i]+=masses[iatom]*velocities[iatom][i];
        }
    }
    return momentum;
}


void System::zero_forces()
{
    for (int iatom=0; iatom<natoms; iatom++)
    {
        for (int i=0; i<3; i++)
        {
            forces[iatom][i]=0.0;
        }
    }
}
