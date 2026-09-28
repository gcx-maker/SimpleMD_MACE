#include "ljforcefield.h"
#include <cmath>

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

    double cutoff2 = cutoff*cutoff;

    for (int iatom=0; iatom<system.natoms; iatom++)
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

            system.potential_energy += 4.0*epsilon*(s12/r12-s6/r6)-engcorrection;

            double coeff = 24.0*epsilon*(2.0*s12/r12-s6/r6)/r2;

            for (int k=0; k<3; k++)
            {
                system.forces[iatom][k] += coeff*dr[k];
                system.forces[jatom][k] -= coeff*dr[k];
            }
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