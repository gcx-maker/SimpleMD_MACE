#ifndef XYZREADER_H
#define XYZREADER_H

#include <string>
#include <vector>
#include "system.h"
#include "input.h"

class XYZReader
{
public:
    XYZReader();
    int read_natoms(const Input& input);
    void read_positions(const Input& input, System& system);
};

#endif