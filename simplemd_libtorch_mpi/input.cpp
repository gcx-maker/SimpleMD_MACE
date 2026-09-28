#include "input.h"
#include <fstream>
#include <sstream>
#include <stdexcept>


Input::Input()
{
    temperature=300.0;
    timestep=0.5;

    friction=0.0;

    cutoff=2.5;
    skin=0.5;

    nstep=1000;
    
    nconfig=10;
    nstat=10;

    warmup_count=0;

    maxneighbor=100;
    idum=-22;

    wrapatoms=false;
    full_neighbor=false;
    enable_profile=false;
    device = "cpu";
}


void Input::read(const std::string& filename)
{
    std::ifstream file(filename);
    
    if (!file)
    {
        throw std::runtime_error("Cannot open input file:" + filename);
    }

    std::string line;

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;

        std::stringstream ss(line);
        std::string key;

        if (!(ss >> key))
            continue;

        if (key[0]=='#')
            continue;

        if (key=="temperature") ss >> temperature;
        else if (key=="timestep") ss >> timestep;
        else if (key=="friction") ss >> friction;
        else if (key=="cutoff") ss >> cutoff;
        else if (key=="skin") ss >> skin;
        else if (key=="nstep") ss >> nstep;
        else if (key=="nconfig")ss >> nconfig >> trajfile;
        else if (key=="nstat") ss >> nstat >> statfile;
        else if (key=="structurefile") ss >> structurefile;
        else if (key=="outputfile") ss >> outputfile;
        else if (key=="model_path") ss >> model_path;
        else if (key=="idum") ss >> idum;
        else if (key=="maxneighbor") ss >> maxneighbor;
        else if (key=="ensemble") ss >> ensemble;
        else if (key=="potential") ss >> potential;
        else if (key=="device") ss >> device;
        else if (key=="warmup_count") ss >> warmup_count;
        else if (key=="wrapatoms")
        {
            std::string value;
            ss >> value;
            if (value=="t" || value=="T") wrapatoms=true;
        }
        else if (key=="enable_profile")
        {
            std::string value;
            ss >> value;
            if (value=="t" || value=="T") enable_profile=true;
        }
        else if (key=="full_neighbor")
        {
            std::string value;
            ss >> value;
            if (value=="t" || value=="T") full_neighbor=true;
        }
        else if (key=="atomtype")
        {
            AtomType atom;
            ss >> atom.id >> atom.name >> atom.mass >> atom.atomic_number;
            atom_type_map[atom.name] = atom;
            type_names.resize(atom.id+1);
            type_names[atom.id] = atom.name;
        }
        else if (key=="pair_coeff")
        {
            LJParameter lj;
            ss >> lj.type1 >> lj.type2 >> lj.sigma >> lj.epsilon;
            ljparameters.push_back(lj);
        }
        else
        {
            throw std::runtime_error("Unknown Keyword: " + key);
        }
    }
}