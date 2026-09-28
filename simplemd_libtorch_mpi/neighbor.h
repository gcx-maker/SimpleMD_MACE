#ifndef NEIGHBOR_H
#define NEIGHBOR_H


#include <vector>
#include "Vector.h"
#include "system.h"
#include "input.h"

struct NeighborPair
{
    int sender;
    int receiver;
    PLMD::Vector shift;
};

class Neighbor
{
private:
    double cutoff;
    double skin;

    bool full_neighbor;

    int nx;
    int ny;
    int nz;

    double cellsize_x;
    double cellsize_y;
    double cellsize_z;

    std::vector<std::vector<int>> cells;

    std::vector<int> list;
    std::vector<int> point;
    std::vector<PLMD::Vector> shifts;

    std::vector<NeighborPair> pairs;

    std::vector<PLMD::Vector> positions0;

    bool enable_profile=false;
    std::vector<double> build_times;

    int warmup_count;
    int build_calls = 0;

public:
    Neighbor(const System& system, const Input& input);
    void build(const System&);
    void build_cell_list(const System& system);
    void compute_list(const System& system);
    void generate_csr(const System& system);
    bool need_update(const System& system);

    const std::vector<int>& point_list() const;
    const std::vector<int>& neighbor_list() const;
    const std::vector<PLMD::Vector>& neighbor_shifts() const;
    void report_profile();
    double cutoff2 = cutoff * cutoff;
};

#endif