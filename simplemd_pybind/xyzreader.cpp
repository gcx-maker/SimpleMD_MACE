#include "xyzreader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>


XYZReader::XYZReader()
{

}


int XYZReader::read_natoms(const Input& input)
{
     std::ifstream file(input.structurefile);
     if (!file)
     {
        throw std::runtime_error(
            "Cannot open structure file:" + input.structurefile
        );
     }
     int natoms;
     file >> natoms;
     if (!file)
     {
        throw std::runtime_error(
            "Failed to read number of atoms"
        );
     }
     return natoms;
}


void XYZReader::read_positions(const Input& input, System& system)
{
    std::ifstream file(input.structurefile);
    if (!file)
    {
        throw std::runtime_error(
            "Cannot open structure file:" + input.structurefile
        );
    }
    std::string line;
    std::getline(file, line);
    std::stringstream ss(line);

    int natoms;
    ss >> natoms;
    
    if (natoms != system.natoms)
    {
        throw std::runtime_error("Atoms number mismatch.");
    }

    std::getline(file,line);
    ss.clear();
    ss.str(line);
    ss >> system.cell[0] >> system.cell[1] >> system.cell[2];

    for (int iatom=0; iatom<natoms; iatom++)
    {
        std::getline(file, line);
        ss.clear();
        ss.str(line);

        std::string atom;
        
        ss >> atom >> system.positions[iatom][0] 
        >> system.positions[iatom][1] >> system.positions[iatom][2];

        if (ss.fail())
        {
            throw std::runtime_error("Invalid xyz coordinate line");
        }

        auto iter = input.atom_type_map.find(atom);
        
        if (iter == input.atom_type_map.end())
        {
            throw std::runtime_error("Unknown atom type:" + atom);
        }

        system.types[iatom] = iter->second.id;
        system.masses[iatom] = iter->second.mass;
        system.atomic_numbers[iatom] = iter->second.atomic_number;
    }
}
