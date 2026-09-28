#include "neighbor.h"

#include <cmath>
#include <algorithm>


static inline double wrap_position(double x, double box)
{
    return x-box*floor(x/box);
}


Neighbor::Neighbor(const System& system, const Input& input)
:cutoff(input.cutoff), skin(input.skin)
{
    point.resize(system.natoms+1);
    list.resize(system.natoms*input.maxneighbor);
    positions0.resize(system.natoms);
}


void Neighbor::build(const System& system)
{
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        positions0[iatom] = system.positions[iatom];
    }

    build_cell_list(system);

    compute_list(system);
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


void Neighbor::compute_list(const System& system)
{
    double R = cutoff + skin;
    double cutoff2 = R*R;
    std::fill(point.begin(), point.end(), 0);

    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        point[iatom+1] = point[iatom];
        double x = wrap_position(system.positions[iatom][0], system.cell[0]);
        double y = wrap_position(system.positions[iatom][1], system.cell[1]);
        double z = wrap_position(system.positions[iatom][2], system.cell[2]);
        
        int ix = static_cast<int>(x/cellsize_x);
        int iy = static_cast<int>(y/cellsize_y);
        int iz = static_cast<int>(z/cellsize_z);

        int rx = static_cast<int>(ceil(R/cellsize_x));
        int ry = static_cast<int>(ceil(R/cellsize_y));
        int rz = static_cast<int>(ceil(R/cellsize_z));

        std::vector<int> visited;

        for (int dx=-rx; dx<=rx; dx++)
        for (int dy=-ry; dy<=ry; dy++)
        for (int dz=-rz; dz<=rz; dz++)
        {
            int cx = (ix+dx+nx)%nx;
            int cy = (iy+dy+ny)%ny;
            int cz = (iz+dz+nz)%nz;

            int cid = cx+nx*(cy+ny*cz);
            if (std::find(visited.begin(), visited.end(), cid) != visited.end())
            {
                continue;
            }
            visited.push_back(cid);
            for (int jatom:cells[cid])
            {
                if (jatom<=iatom)
                    continue;
                PLMD::Vector dr;
                double r2=0;
                for (int k=0; k<3; k++)
                {
                    dr[k] = system.positions[iatom][k] - system.positions[jatom][k];
                    dr[k]-=floor(dr[k]/system.cell[k]+0.5)*system.cell[k]; 
                    r2 += dr[k]*dr[k];
                }
                if (r2>cutoff2)
                {
                    continue;
                }
                if (point[iatom+1] >= static_cast<int>(list.size()))
                {
                    list.resize(list.size()*2);
                }
                list[point[iatom+1]] = jatom;
                point[iatom+1]++;
            }
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