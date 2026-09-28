#include "neighbor.h"

#include <cmath>
#include <algorithm>
#include <set>
#include <tuple>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iostream>


static inline double wrap_position(double x, double box)
{
    return x-box*floor(x/box);
}


Neighbor::Neighbor(const System& system, const Input& input)
:cutoff(input.cutoff), skin(input.skin), full_neighbor(input.full_neighbor)
{
    warmup_count = input.warmup_count;
    enable_profile = input.enable_profile;
    build_times.reserve(input.nstep+10);
    point.resize(system.natoms+1);
    list.resize(system.natoms*input.maxneighbor);
    positions0.resize(system.natoms);
}


void Neighbor::build(const System& system)
{
    auto start = std::chrono::high_resolution_clock::now();
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        positions0[iatom] = system.positions[iatom];
    }

    build_cell_list(system);

    compute_list(system);

    generate_csr(system);
    auto end = std::chrono::high_resolution_clock::now();

    if(enable_profile)
    {
        build_calls++;

        if(build_calls > warmup_count)
        {
            build_times.push_back(
                std::chrono::duration<double,std::milli>
                (end-start).count()
            );
        }
    }
}


void Neighbor::build_cell_list(const System& system)
{
    double R = cutoff + skin;
    nx = static_cast<int>(system.cell[0]/R);
    ny = static_cast<int>(system.cell[1]/R);
    nz = static_cast<int>(system.cell[2]/R);

    if (nx<1) nx=1;
    if (ny<1) ny=1;
    if (nz<1) nz=1;

    cellsize_x = system.cell[0]/nx;
    cellsize_y = system.cell[1]/ny;
    cellsize_z = system.cell[2]/nz;

    int ncells = nx*ny*nz;
    cells.clear();
    cells.resize(ncells);
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        double x = wrap_position(system.positions[iatom][0], system.cell[0]);
        double y = wrap_position(system.positions[iatom][1], system.cell[1]);
        double z = wrap_position(system.positions[iatom][2], system.cell[2]);

        int ix = static_cast<int>(x/cellsize_x);
        int iy = static_cast<int>(y/cellsize_y);
        int iz = static_cast<int>(z/cellsize_z);

        if (ix>=nx) ix=nx-1;
        if (iy>=ny) iy=ny-1;
        if (iz>=nz) iz=nz-1;

        int id=ix+nx*(iy+ny*iz);
        cells[id].push_back(iatom);
    }
}


// void Neighbor::compute_list(const System& system)
// {
//     pairs.clear();
//     double R = cutoff + skin;
//     double cutoff2 = R*R;
//     // std::fill(point.begin(), point.end(), 0);

//     for (int iatom=0; iatom<system.natoms; iatom++)
//     {
//         // point[iatom+1] = point[iatom];
//         double x = wrap_position(system.positions[iatom][0], system.cell[0]);
//         double y = wrap_position(system.positions[iatom][1], system.cell[1]);
//         double z = wrap_position(system.positions[iatom][2], system.cell[2]);
        
//         int ix = static_cast<int>(x/cellsize_x);
//         int iy = static_cast<int>(y/cellsize_y);
//         int iz = static_cast<int>(z/cellsize_z);

//         int rx = static_cast<int>(ceil(R/cellsize_x));
//         int ry = static_cast<int>(ceil(R/cellsize_y));
//         int rz = static_cast<int>(ceil(R/cellsize_z));

//         std::vector<int> visited;

//         for (int dx=-rx; dx<=rx; dx++)
//         for (int dy=-ry; dy<=ry; dy++)
//         for (int dz=-rz; dz<=rz; dz++)
//         {
//             int cx = (ix+dx+nx)%nx;
//             int cy = (iy+dy+ny)%ny;
//             int cz = (iz+dz+nz)%nz;

//             int cid = cx+nx*(cy+ny*cz);
//             if (std::find(visited.begin(), visited.end(), cid) != visited.end())
//             {
//                 continue;
//             }
//             visited.push_back(cid);
//             for (int jatom:cells[cid])
//             {
//                 if (jatom<=iatom)
//                     continue;
//                 PLMD::Vector shift;
//                 PLMD::Vector dr;
//                 double r2=0;
//                 for (int k=0; k<3; k++)
//                 {
//                     double delta = system.positions[jatom][k]-system.positions[iatom][k];
//                     int nshift = static_cast<int>(floor(delta/system.cell[k]+0.5));
//                     shift[k] = -nshift*system.cell[k];
//                     dr[k] = delta+shift[k];
//                     // dr[k] = system.positions[iatom][k] - system.positions[jatom][k];
//                     // dr[k]-=floor(dr[k]/system.cell[k]+0.5)*system.cell[k]; 
//                     r2 += dr[k]*dr[k];
//                 }
//                 if (r2>cutoff2)
//                 {
//                     continue;
//                 }
//                 NeighborPair p;
//                 p.sender = iatom;
//                 p.receiver = jatom;
//                 p.shift = shift;
//                 pairs.push_back(p);
//                 // if (point[iatom+1] >= static_cast<int>(list.size()))
//                 // {
//                 //     list.resize(list.size()*2);
//                 // }
//                 // list[point[iatom+1]] = jatom;
//                 // point[iatom+1]++;
//             }
//         }
//     }
// }


void Neighbor::compute_list(const System& system)
{
    pairs.clear();

    double R = cutoff + skin;
    double cutoff2 = R * R;


    int nx_shift = static_cast<int>(
        ceil(R / system.cell[0])
    );

    int ny_shift = static_cast<int>(
        ceil(R / system.cell[1])
    );

    int nz_shift = static_cast<int>(
        ceil(R / system.cell[2])
    );


    std::set<
        std::tuple<int,int,int,int,int>
    > existing;



    for(int iatom=0; iatom<system.natoms; iatom++)
    {

        double x = wrap_position(
            system.positions[iatom][0],
            system.cell[0]
        );

        double y = wrap_position(
            system.positions[iatom][1],
            system.cell[1]
        );

        double z = wrap_position(
            system.positions[iatom][2],
            system.cell[2]
        );


        int ix = static_cast<int>(
            x / cellsize_x
        );

        int iy = static_cast<int>(
            y / cellsize_y
        );

        int iz = static_cast<int>(
            z / cellsize_z
        );



        int rx = static_cast<int>(
            ceil(R / cellsize_x)
        );

        int ry = static_cast<int>(
            ceil(R / cellsize_y)
        );

        int rz = static_cast<int>(
            ceil(R / cellsize_z)
        );



        std::vector<int> visited;



        for(int dx=-rx; dx<=rx; dx++)
        for(int dy=-ry; dy<=ry; dy++)
        for(int dz=-rz; dz<=rz; dz++)
        {


            int cx = (ix+dx+nx)%nx;
            int cy = (iy+dy+ny)%ny;
            int cz = (iz+dz+nz)%nz;


            int cid = cx + nx*(cy+ny*cz);



            if(std::find(
                visited.begin(),
                visited.end(),
                cid
            ) != visited.end())
            {
                continue;
            }

            visited.push_back(cid);



            for(int jatom : cells[cid])
            {


                if(jatom < iatom)
                    continue;



                for(int sx=-nx_shift;
                    sx<=nx_shift;
                    sx++)

                for(int sy=-ny_shift;
                    sy<=ny_shift;
                    sy++)

                for(int sz=-nz_shift;
                    sz<=nz_shift;
                    sz++)
                {


                    if(jatom==iatom &&
                       sx==0 &&
                       sy==0 &&
                       sz==0)
                    {
                        continue;
                    }


                    PLMD::Vector shift;


                    shift[0] =
                        sx * system.cell[0];

                    shift[1] =
                        sy * system.cell[1];

                    shift[2] =
                        sz * system.cell[2];



                    PLMD::Vector dr;

                    double r2 = 0.0;



                    for(int k=0;k<3;k++)
                    {
                        dr[k] =
                            system.positions[jatom][k]
                            +
                            shift[k]
                            -
                            system.positions[iatom][k];

                        r2 += dr[k]*dr[k];
                    }



                    if(r2 > cutoff2)
                        continue;

                    auto key =
                    std::make_tuple(
                        iatom,
                        jatom,
                        sx,
                        sy,
                        sz
                    );


                    if(existing.count(key))
                        continue;


                    existing.insert(key);



                    NeighborPair p;

                    p.sender = iatom;

                    p.receiver = jatom;

                    p.shift = shift;


                    pairs.push_back(p);

                }
            }
        }
    }
}

void Neighbor::generate_csr(const System& system)
{
    point.assign(system.natoms+1, 0);
    for (auto& p:pairs)
    {
        point[p.sender+1]++;
        if (full_neighbor)
        {
            point[p.receiver+1]++;
        }
    }
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        point[iatom+1]+=point[iatom];
    }
    list.resize(point.back());
    shifts.resize(point.back());
    std::vector<int> cursor = point;
    for (auto& p:pairs)
    {
        int i=p.sender;
        int j=p.receiver;
        int index = cursor[i];
        list[index]=j;
        shifts[index]=p.shift;
        cursor[i]++;
        if (full_neighbor)
        {
            int index2 = cursor[j];
            list[index2]=i;
            shifts[index2]=-p.shift;
            cursor[j]++;
        }
    }
}


bool Neighbor::need_update(const System& system)
{
    double threshold = skin*0.5;
    double threshold2 = threshold*threshold;

    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        PLMD::Vector dr;
        double d2=0;
        for (int k=0; k<3; k++)
        {
            dr[k] = system.positions[iatom][k]-positions0[iatom][k];
            dr[k] -= floor(dr[k]/system.cell[k]+0.5)*system.cell[k];
            d2 += dr[k]*dr[k];
        }
        if (d2>threshold2)
        {
            return true;
        }
    }
    return false;
}


const std::vector<int>& Neighbor::point_list() const
{
    return point;
}


const std::vector<int>& Neighbor::neighbor_list() const
{
    return list;
}

const std::vector<PLMD::Vector>& Neighbor::neighbor_shifts() const
{
    return shifts;
}


void Neighbor::report_profile()
{
    if(!enable_profile || build_times.empty())
        return;


    double sum = 0;

    double min_time =
        *std::min_element(
            build_times.begin(),
            build_times.end()
        );

    double max_time =
        *std::max_element(
            build_times.begin(),
            build_times.end()
        );


    for(auto t : build_times)
        sum += t;


    std::cout
    <<"========== Neighbor Build ==========\n"
    <<"calls : "<<build_times.size()<<"\n"
    <<"avg   : "<<sum/build_times.size()<<" ms\n"
    <<"min   : "<<min_time<<" ms\n"
    <<"max   : "<<max_time<<" ms\n"
    <<"====================================\n";


    std::ofstream fout("neighbor_build_times.dat");

    for(auto t : build_times)
    {
        fout << t << "\n";
    }
}