#ifndef NEIGHBOR_H
#define NEIGHBOR_H


#include <vector>
#include "Vector.h"
#include "system.h"
#include "input.h"


class Neighbor
{
private:
    double cutoff;
    double skin;

    int nx;
    int ny;
    int nz;

    double cellsize_x;
    double cellsize_y;
    double cellsize_z;

    std::vector<std::vector<int>> cells;
    std::vector<int> list;
    std::vector<int> point;
    std::vector<PLMD::Vector> positions0;

public:
    Neighbor(const System& system, const Input& input);
    void build(const System&);
    void build_cell_list(const System& system);
    void compute_list(const System& system);
    bool need_update(const System& system);

    const std::vector<int>& point_list() const;
    const std::vector<int>& neighbor_list() const;
};

#endif