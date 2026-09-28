#include "system.h"
#include "xyzreader.h"
#include "input.h"
#include "setup.h"
#include "constants.h"
#include "neighbor.h"
#include "ljforcefield.h"
#include "integrator.h"
#include "thermostat.h"
#include "Random.h"
#include "output.h"
#include "maceforcefield.h"

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <cstdio> 
#include <cstdlib> 
#include <cstring> 
#include <cmath> 

using namespace std;
using namespace PLMD;

int main(int argc, char* argv[])
{
    py::scoped_interpreter guard{};
    if(argc < 2)
    {
        throw std::runtime_error(
            "Usage: ./main.exe input_file"
        );
    }
    Input input;
    input.read(argv[1]);

    XYZReader xyzreader;
    int n = xyzreader.read_natoms(input);
    fprintf(stdout, "%s %s\n", "Starting configuration           :", input.structurefile.c_str());
    fprintf(stdout, "%s %s\n", "Final configuration              :", input.outputfile.c_str());
    fprintf(stdout, "%s %d\n", "Number of atoms                  :", n);
    fprintf(stdout, "%s %f\n", "Temperature                      :", input.temperature);
    fprintf(stdout, "%s %f\n", "Time step                        :", input.timestep);
    fprintf(stdout, "%s %f\n", "Friction                         :", input.friction);
    fprintf(stdout, "%s %f\n", "Cutoff for forces                :", input.cutoff);
    fprintf(stdout, "%s %f\n", "Skin                             :", input.skin);
    fprintf(stdout, "%s %d\n", "Number of steps                  :", input.nstep);
    fprintf(stdout, "%s %d\n", "Stride for trajectory            :", input.nconfig);
    fprintf(stdout, "%s %s\n", "Trajectory file                  :", input.trajfile.c_str());
    fprintf(stdout, "%s %d\n", "Stride for statistics            :", input.nstat);
    fprintf(stdout, "%s %s\n", "Statistics file                  :", input.statfile.c_str());
    fprintf(stdout, "%s %d\n", "Max average number of neighbors :", input.maxneighbor);
    fprintf(stdout, "%s %d\n", "Seed                             :", input.idum);
    fprintf(stdout, "%s %s\n", "Are atoms wrapped on output      :", (input.wrapatoms?"T":"F"));
    fprintf(stdout, "%s %s\n", "Ensemble                         :", input.ensemble.c_str());
    fprintf(stdout, "%s %s\n", "Potential                        :", input.potential.c_str());

    System system(n);

    xyzreader.read_positions(input, system);

    Random random;
    random.setSeed(input.idum);

    Setup setup;
    setup.initialize(system, input, random);

    system.compute_kinetic_energy();

    cout << "Initial T=" << 2*system.kinetic_energy/((3*system.natoms-3)*Constants::KB) << endl;


    Neighbor neighbor(system, input);
    if (input.potential=="LJ")
    {
        neighbor.build(system);
    }
    Intergrator integrator(input);
    Thermostat thermostat(input);


    std::unique_ptr<ForceField> forcefield;
    if (input.potential=="LJ")
    {
        forcefield = std::make_unique<LJForceField>(input);
    }
    else if (input.potential=="MACE")
    {
        forcefield = std::make_unique<MACEForceField>(input);
    }
    else
    {
        throw std::runtime_error("Unknown potential");
    }
    
    forcefield->compute(system, neighbor, input);



    Output output;

    bool use_thermostat = false;
    if (input.ensemble == "NVT")
    {
        use_thermostat = true;
    }
    else if(input.ensemble=="NVE")
    {
        use_thermostat=false;
    }
    else
    {
        throw std::runtime_error("Unknown ensemble");
    }
    auto start_time = std::chrono::high_resolution_clock::now();  
    for (int step=0; step<=input.nstep; step++)
    {
        if (step%input.nconfig == 0)
        {
            output.write_positions(input, system);
        }
        if (step%input.nstat == 0)
        {
            output.write_statistics(step, input, system);
        }
        if (step==input.nstep)
        {
            break;
        }
        integrator.first_step(system);
        integrator.update_positions(system, input);
        if (input.potential == "LJ")
        {
            if (neighbor.need_update(system))
            {
                neighbor.build(system);
                printf("Neighbor rebuild at %d\n", step);
            }
        }
        forcefield->compute(system, neighbor, input);
        integrator.second_step(system);
        if (use_thermostat)
        {
            thermostat.apply(system,input, random);
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed_time = std::chrono::duration<double>(end_time-start_time).count();
    double steps_per_second = input.nstep / elapsed_time;
    double ns_per_day = steps_per_second *input.timestep *1e-6 *86400.0;
    printf("\n==============================\n");
    printf("MD Performance\n");
    printf("==============================\n");
    printf("Total time       : %.3f s\n", elapsed_time);
    printf("Steps            : %d\n", input.nstep);
    printf("Step time        : %.6f ms\n",elapsed_time/input.nstep*1000.0);
    printf("Steps/sec        : %.2f\n",steps_per_second);
    printf("Simulation speed : %.4f ns/day\n",ns_per_day);
    printf("==============================\n\n");
    
    output.write_final_positions(input, system);
    if (output.write_statistics_fp) fclose(output.write_statistics_fp);
    return 0;
}