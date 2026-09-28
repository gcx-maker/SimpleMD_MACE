#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include "system.h"
#include "neighbor.h"
#include "input.h"

class ForceField
{
protected:
    double cutoff;
public:
    ForceField(double cutoff):cutoff(cutoff)
    {

    }
    virtual ~ForceField()
    {

    }
    virtual void compute(System& system, const Neighbor& neighbor, const Input& input) = 0;
    virtual void report_profile(){}

};
#endif