#ifndef MACEFORCEFIELD_H
#define MACEFORCEFIELD_H

#include <pybind11/embed.h>
#include <pybind11/numpy.h>

#include "system.h"
#include "input.h"
#include "forcefield.h"
#include "neighbor.h"


namespace py = pybind11;


class MACEForceField : public ForceField
{
private:
    py::object module;
    py::object calculate;
    py::array_t<int> atomic_numbers;
    bool atom_number_initialized;
    py::array_t<double> positions;
    bool buffer_initialized;

public:
    MACEForceField(const Input& input);
    void compute(System& system, const Neighbor& neighbor, const Input& input) override;
};

#endif