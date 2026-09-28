#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <vector>
#include <unordered_map>

struct AtomType
{
    int id;
    std::string name;
    double mass;
    int atomic_number;
};


struct LJParameter
{
    int type1 = 0;
    int type2 = 0;

    double sigma = 0.0;
    double epsilon = 0.0;
};


class Input
{
public:
    double temperature;
    double timestep;
    double friction;

    double cutoff;
    double skin;

    int nstep;
    int nconfig;
    int nstat;

    int idum;
    int maxneighbor;

    bool wrapatoms;

    std::string outputfile;
    std::string trajfile;
    std::string statfile;
    std::string structurefile;

    std::string ensemble;
    std::string potential;

    std::unordered_map<std::string, AtomType> atom_type_map;
    std::vector<std::string> type_names;

    std::vector<LJParameter> ljparameters;

    Input();

    void read(const std::string& filename);
};

#endif