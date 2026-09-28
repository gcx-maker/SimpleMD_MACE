#include "maceforcefield.h"

#include <iostream>
#include <vector>

MACEForceField::MACEForceField(const Input& input) : ForceField(input.cutoff),
atom_number_initialized(false), buffer_initialized(false)
{
    module = py::module_::import("mace_calculate");
    calculate = module.attr("calculate");
    std::cout << "MACE initialized" << std::endl;
}

void MACEForceField::compute(System& system, const Neighbor& neighbor, const Input& input)
{
    (void)neighbor;
    (void)input;
    if (!buffer_initialized)
    {
        positions = py::array_t<double>({system.natoms, 3});
        buffer_initialized = true;
    }
    if (!atom_number_initialized)
    {
        atomic_numbers = py::array_t<int>({system.natoms}, system.atomic_numbers.data());
        atom_number_initialized = true;
    }

    auto pos_ptr = positions.mutable_data();
    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        pos_ptr[3*iatom] = system.positions[iatom][0];
        pos_ptr[3*iatom+1] = system.positions[iatom][1];
        pos_ptr[3*iatom+2] = system.positions[iatom][2];
    }

    double cell_data[3] = {system.cell[0], system.cell[1], system.cell[2]};
    py::array_t<double> cell({3}, cell_data);

    py::object result = calculate(positions, cell, atomic_numbers);
    py::tuple tuple = result.cast<py::tuple>();

    double energy = tuple[0].cast<double>();
    system.potential_energy = energy;

    py::array_t<double> forces=tuple[1].cast<py::array_t<double>>();
    auto buf = forces.request();

    double* fptr = static_cast<double*>(buf.ptr);

    for (int iatom=0; iatom<system.natoms; iatom++)
    {
        system.forces[iatom][0] = fptr[3*iatom];
        system.forces[iatom][1] = fptr[3*iatom+1];
        system.forces[iatom][2] = fptr[3*iatom+2];
    }
}