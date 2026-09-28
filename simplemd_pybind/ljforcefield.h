#ifndef LJFORCEFIELD_H
#define LJFORCEFIELD_H

#include "forcefield.h"
#include "input.h"
#include <vector>


class LJForceField : public ForceField
{
private:
    void build_lj_table(const Input& input);
    std::vector<std::vector<LJParameter>> lj_table;
public:
    LJForceField(const Input& input);
    void compute(System& system, const Neighbor& neighbor, const Input& input) override;
};
#endif
