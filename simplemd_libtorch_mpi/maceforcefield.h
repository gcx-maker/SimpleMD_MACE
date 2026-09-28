#ifndef MACEFORCEFIELD_H
#define MACEFORCEFIELD_H

#include "system.h"
#include "input.h"
#include "forcefield.h"
#include "neighbor.h"
#include "macedatabuilder.h"

#include <torch/script.h>
#include <torch/torch.h>

#include <vector>
#include <memory>


class MACEForceField : public ForceField
{
private:
    torch::jit::script::Module model;
    torch::Device device;
    std::unique_ptr<MACEDataBuilder> builder;
    bool enable_profile=false;
    std::vector<double> forward_times;
    int warmup_count;
    int forward_calls = 0;

public:
    MACEForceField(const Input& input, const System& system);
    void compute(System& system, const Neighbor& neighbor, const Input& input) override;
    void report_profile();
};

#endif