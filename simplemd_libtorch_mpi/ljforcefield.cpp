#include "ljforcefield.h"
#include <cmath>
#include <mpi.h>
#include <iostream>
 using namespace std;


LJForceField::LJForceField(const Input& input): ForceField(input.cutoff)
{
    build_lj_table(input);
}


void LJForceField::compute(System& system, const Neighbor& neighbor, const Input& input)
{
    system.potential_energy = 0.0;
    system.zero_forces();

    const auto& point = neighbor.point_list();
    const auto& list = neighbor.neighbor_list();

    const double cutoff2 = cutoff*cutoff;


    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int begin = rank * system.natoms / size;
    int end = (rank+1) * system.natoms / size;

    std::vector<double> local_forces(3*system.natoms, 0.0);
    double local_energy = 0.0;

    for (int iatom=begin; iatom<end; iatom++)
    {
        int type_i = system.types[iatom];
        for (int j=point[iatom]; j<point[iatom+1]; j++)
        {

            int jatom = list[j];

            int type_j = system.types[jatom];

            const LJParameter& lj = lj_table[type_i][type_j];

            double sigma = lj.sigma;
            double epsilon = lj.epsilon;

            PLMD::Vector dr;
            double r2=0.0;
            for (int k=0; k<3; k++)
            {
                dr[k] = system.positions[iatom][k] - system.positions[jatom][k];
                dr[k] -= floor(dr[k]/system.cell[k] + 0.5) * system.cell[k];
                r2 += dr[k]*dr[k];
            }

            if (r2>cutoff2)
            {
                continue;
            }

            double r6 = r2*r2*r2;
            double r12 = r6*r6;
            double s6 = pow(sigma, 6);
            double s12 = s6*s6;

            double engcorrection = 4.0*epsilon*(s12/pow(cutoff2, 6)-s6/pow(cutoff2, 3));

            local_energy += 4.0*epsilon*(s12/r12-s6/r6)-engcorrection;

            double coeff = 24.0*epsilon*(2.0*s12/r12-s6/r6)/r2;

            for (int k=0; k<3; k++)
            {
                local_forces[3*iatom+k] += coeff*dr[k];
                local_forces[3*jatom+k] -= coeff*dr[k];
            }
        }
    }

    double global_energy = 0.0;
    MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std::vector<double> global_forces(
        3 * system.natoms,
        0.0
    );

    MPI_Allreduce(
        local_forces.data(),
        global_forces.data(),
        3 * system.natoms,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    system.potential_energy = global_energy;
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        for (int k=0; k<3; k++)
        {
            system.forces[iatom][k] = global_forces[3*iatom + k];
        }
    }
}

void LJForceField::build_lj_table(const Input& input)
{
    int ntypes = input.atom_type_map.size();
    lj_table.resize(ntypes);
    for (int i=0; i<ntypes; i++)
    {
        lj_table[i].resize(ntypes);
    }
    for (auto& lj:input.ljparameters)
    {
        lj_table[lj.type1][lj.type2] = lj;
        lj_table[lj.type2][lj.type1] = lj;
    }
}