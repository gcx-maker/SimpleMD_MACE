#include "setup.h"

void Setup::initialize(System& system, const Input& input, PLMD::Random& random)
{
    randomize_velocities(system, random, input);
    remove_center_of_mass_velocity(system);
    rescale_temperature(input, system);
}

void Setup::randomize_velocities(System& system, PLMD::Random& random, const Input& input)
{
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        for (int i=0; i<3; i++)
        {
            system.velocities[iatom][i] = sqrt(Constants::KB * input.temperature 
                * Constants::EV_AMU_TO_A2_FS2 / system.masses[iatom]) * random.Gaussian();
        }
    }
}


void Setup::remove_center_of_mass_velocity(System& system)
{
    PLMD::Vector momentum = system.compute_momentum();
    double total_mass = 0.0;
    for (auto& m:system.masses)
    {
        total_mass += m;
    }
    for (int i=0; i<3; i++)
    {
        momentum[i] /= total_mass;
    }
    for (auto& v:system.velocities)
    {
        for (int i=0; i<3; i++)
        {
            v[i]-=momentum[i];
        }
    }
}

void Setup::rescale_temperature(const Input& input, System& system)
{
    system.compute_kinetic_energy();
    double current_temperature = 2*system.kinetic_energy/((3*system.natoms-3)*Constants::KB);
    double scale = sqrt(input.temperature/current_temperature);
    for (auto& v:system.velocities)
    {
        for (int i=0; i<3; i++)
        {
            v[i]*=scale;
        }
    }
}