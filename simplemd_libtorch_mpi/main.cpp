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


#include <mpi.h>

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <cstdio> 
#include <cstdlib> 
#include <cstring> 
#include <cmath> 
#include <torch/script.h>

using namespace std;
using namespace PLMD;


void broadcast_system(System& system)
{

    MPI_Bcast(&system.positions[0][0],system.natoms*3,MPI_DOUBLE,0,MPI_COMM_WORLD);

    MPI_Bcast(&system.velocities[0][0],system.natoms*3,MPI_DOUBLE,0,MPI_COMM_WORLD);

    MPI_Bcast(system.types.data(),system.natoms,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Bcast(system.masses.data(),system.natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);

    MPI_Bcast(system.cell,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
}



int main(
    int argc,
    char* argv[]
)
{
    MPI_Init(&argc,&argv);

    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);


    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if(argc < 2)
    {
        if(rank==0)
        {
            cerr << "Usage: main.exe input_file" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    Input input;

    input.read(
        argv[1]
    );

    if(rank==0)
    {
        cout
        <<"============================\n"
        <<"MPI MD START\n"
        <<"============================\n";

        cout << "MPI processes : " << size <<endl;
    }

    XYZReader xyzreader;


    int natoms = xyzreader.read_natoms(input);



    if(rank==0)
    {
        fprintf(stdout, "%s %s\n", "Starting configuration           :", input.structurefile.c_str());
        fprintf(stdout, "%s %s\n", "Final configuration              :", input.outputfile.c_str());
        fprintf(stdout, "%s %d\n", "Number of atoms                  :", natoms);
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
        fprintf(stdout, "%s %s\n", "Full_neighbor                    :", (input.full_neighbor?"T":"F"));
        fprintf(stdout, "%s %s\n", "Ensemble                         :", input.ensemble.c_str());
        fprintf(stdout, "%s %s\n", "Potential                        :", input.potential.c_str());
        fprintf(stdout, "%s %s\n", "Model_path                       :", input.model_path.c_str());
        fprintf(stdout, "%s %s\n", "Device                           :", input.device.c_str());
    }

    System system(natoms);
    Random random;

    if(rank==0)
    {
        xyzreader.read_positions(input,system);

        random.setSeed(input.idum);

        Setup setup;

        setup.initialize(system,input,random);
    }


    broadcast_system(system);
    system.compute_kinetic_energy();

    if(rank==0)
    {
        cout << "Initial T = " << 2*system.kinetic_energy/((3*system.natoms-3)*Constants::KB) << endl;
    }

    Neighbor neighbor(system,input);

    neighbor.build(system);


    unique_ptr<ForceField> forcefield;

    if(input.potential=="LJ")
    {
        forcefield =make_unique<LJForceField>(input);
    }
    else if(input.potential=="MACE")
    {
        forcefield =make_unique<MACEForceField>(input,system);
    }
    else
    {
        throw runtime_error("Unknown potential");
    }

    forcefield->compute(system,neighbor,input);

    Integrator integrator(input);

    Thermostat thermostat(input);

    bool use_thermostat=false;

    if(input.ensemble=="NVT")
    {
        use_thermostat=true;
    }
    else if(input.ensemble=="NVE")
    {
        use_thermostat=false;
    }
    else
    {
        throw runtime_error("Unknown ensemble");
    }

    unique_ptr<Output> output;

    if(rank==0)
    {
        output = make_unique<Output>();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    auto start = chrono::high_resolution_clock::now();

    for(int step=0;step<=input.nstep;step++)
    {
        if(rank==0)
        {
            if(step%input.nconfig==0)
            {
                output->write_positions(input,system);
            }

            if(step%input.nstat==0)
            {
                output->write_statistics(step,input,system);
            }
        }

        integrator.first_step(system);

        integrator.update_positions(system,input);

        MPI_Bcast(&system.positions[0][0],system.natoms*3,MPI_DOUBLE,0,MPI_COMM_WORLD);

        if(neighbor.need_update(system))
        {
            neighbor.build(system);
        }

        forcefield->compute(system,neighbor,input);

        integrator.second_step(system);

        if(use_thermostat)
        {
            thermostat.apply(system,input,random);
        }

        MPI_Bcast(&system.velocities[0][0],system.natoms*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    auto end = chrono::high_resolution_clock::now();

    double elapsed = chrono::duration<double>(end-start).count();

    if(rank==0)
    {
        double step_per_second = input.nstep / elapsed;

        double ns_day = step_per_second * input.timestep * 1e-6 * 86400;

        cout << "\n============================\n"
        << "MD Performance\n"
        << "============================\n";

        cout << "Time : " <<elapsed << " s\n";

        cout << "Steps/s : " << step_per_second << endl;

        cout << "Speed : " << ns_day << " ns/day\n";
        
        cout << "============================\n";


        neighbor.report_profile();
        
        forcefield->report_profile();

        output->write_final_positions(input,system);

        if(output->write_statistics_fp)
        {
            fclose(output->write_statistics_fp);
        }
    }
    MPI_Finalize();
    return 0;
}
