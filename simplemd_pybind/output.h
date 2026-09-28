#ifndef OUPUT_H
#define OUPUT_H

#include <vector>
#include <string>
#include "Vector.h"
#include "input.h"
#include "constants.h"
#include "system.h"

class Output
{
public:
    bool write_positions_first;
    bool write_statistics_first;
    int write_statistic_last_time_reopened;
    FILE* write_statistics_fp;

    Output();
    void write_positions(const Input& input, const System& system);
    void write_final_positions(const Input& input, const System& system);
    void write_statistics(int istep, const Input& input, System& system);
};
#endif