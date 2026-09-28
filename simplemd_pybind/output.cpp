#include "output.h"


Output::Output()
{
    write_positions_first = true;
    write_statistics_first = true;
    write_statistic_last_time_reopened = 0;
    write_statistics_fp = nullptr;
}


void Output::write_positions(const Input& input, const System& system)
{
    PLMD::Vector pos;
    FILE* fp;
    if (write_positions_first)
    {
        fp = fopen(input.trajfile.c_str(), "w");
        write_positions_first = false;
    }
    else
    {
        fp = fopen(input.trajfile.c_str(), "a");
    }
    fprintf(fp, "%d\n", system.natoms);
    fprintf(fp, "%f %f %f\n", system.cell[0], system.cell[1], system.cell[2]);
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        if (input.wrapatoms)
        {
            for (int k=0; k<3; k++)
            {
                pos[k] = system.positions[iatom][k]-
                floor(system.positions[iatom][k]/system.cell[k])*system.cell[k];
            }
        }
        else
        {
            for (int k=0; k<3; k++)
            {
                pos[k] = system.positions[iatom][k];
            }
        }
        int type = system.types[iatom];
        fprintf(fp, "%s %10.7f %10.7f %10.7f\n", input.type_names[type].c_str(), pos[0], pos[1], pos[2]);
    }
    fclose(fp);
}


void Output::write_final_positions(const Input& input, const System& system)
{
    PLMD::Vector pos;
    FILE *fp;
    fp=fopen(input.outputfile.c_str(), "w");
    fprintf(fp, "%d\n", system.natoms);
    fprintf(fp, "%f %f %f\n", system.cell[0], system.cell[1], system.cell[2]);
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        if (input.wrapatoms) 
        {
            for (int k=0; k<3; k++)
            {
                pos[k] = system.positions[iatom][k]-floor(system.positions[iatom][k]/system.cell[k])* system.cell[k];
            }
        }
        else for (int k=0; k<3; k++) pos[k]=system.positions[iatom][k];
        int type = system.types[iatom];
        fprintf(fp, "%s %10.7f %10.7f %10.7f\n",input.type_names[type].c_str(), pos[0], pos[1], pos[2]);
    }
    fclose(fp);
}

void Output::write_statistics(int step, const Input& input, System& system)
{
    if (write_statistics_first)
    {
        write_statistics_fp = fopen(input.statfile.c_str(), "w");
        fprintf(write_statistics_fp,
            "# Step Time(fs) Temperature(K)"
            "E_pot(eV) E_kin(eV) E_tot(eV)\n");
        write_statistics_first=false;
    }
    if (step-write_statistic_last_time_reopened>100)
    {
        fclose(write_statistics_fp);
        write_statistics_fp=fopen(input.statfile.c_str(), "a");
        write_statistic_last_time_reopened = step;
    }
    system.compute_kinetic_energy();
    fprintf(write_statistics_fp, "%d %f %f %f %f %f\n", 
        step, step*input.timestep, 2*system.kinetic_energy/((3*system.natoms-3)*Constants::KB), 
        system.potential_energy, system.kinetic_energy, system.potential_energy+system.kinetic_energy);
}
