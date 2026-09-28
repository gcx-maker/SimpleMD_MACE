#include "maceforcefield.h"
#include <chrono>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>


MACEForceField::MACEForceField(const Input& input, const System& system) : ForceField(input.cutoff), device(torch::kCPU)
{
    enable_profile = input.enable_profile;
    forward_times.reserve(input.nstep+10);
    warmup_count = input.warmup_count;
    if (input.device == "cuda")
    {
        if(torch::cuda::is_available())
        {
            device=torch::kCUDA;
        }
        else
        {
            throw std::runtime_error(
                "CUDA requested but unavailable"
            );
        }
    }
    //model = torch::jit::load(input.model_path);
    try
    {
        std::cout << "loading model: "
                << input.model_path
                << std::endl;

        model = torch::jit::load(input.model_path, torch::kCPU);

        std::cout << "model loaded" << std::endl;
    }
    catch (const c10::Error& e)
    {
        std::cerr << "TorchScript load error:\n"
                << e.what()
                << std::endl;
        throw;
    }
    catch (const std::exception& e)
    {
        std::cerr << "std error:\n"
                << e.what()
                << std::endl;
        throw;
    }

    std::vector<int> z_table;
    torch::Tensor atomic_numbers =
        model.attr("atomic_numbers").toTensor().cpu();

    auto ptr = atomic_numbers.data_ptr<int64_t>();

    for (int i = 0; i < atomic_numbers.size(0); i++)
    {
        z_table.push_back(static_cast<int>(ptr[i]));
    }

    builder = std::make_unique<MACEDataBuilder>(z_table, device);
    builder -> initialize(system);

    std::cout << "MACE initialized" << std::endl;
    model.to(device);
    model.to(torch::kFloat32);
    model.eval();
}



void MACEForceField::compute(System& system, const Neighbor& neighbor, const Input& input)
{
    (void)input;
    builder->build(system, neighbor);
    const auto& mace_data = builder->data();
    auto data_dict = c10::impl::GenericDict(c10::StringType::get(), c10::TensorType::get());
    for (auto& item:mace_data)
    {
        data_dict.insert(item.first, item.second);
    }
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(data_dict);
    inputs.push_back(false);
    inputs.push_back(true); 
    inputs.push_back(false);
    inputs.push_back(false);
    inputs.push_back(false);
    inputs.push_back(false);
    inputs.push_back(false);
    inputs.push_back(false);
    inputs.push_back(false);

    double forward_time = 0.0;
    torch::jit::IValue output;

    auto start = std::chrono::high_resolution_clock::now();

    output = model.forward(inputs);

    if(device.is_cuda())
    {
        torch::cuda::synchronize();
    }

    auto end = std::chrono::high_resolution_clock::now();


    forward_time =std::chrono::duration<double,std::milli>(end-start).count();

    if(enable_profile)
    {
        forward_calls++;

        if(forward_calls > warmup_count)
        {
            forward_times.push_back(
                forward_time
            );
        }
    }

    auto output_dict = output.toGenericDict();

    auto energy = output_dict.at("energy").toTensor();
    system.potential_energy = energy.cpu().item<float>();

    auto forces = output_dict.at("forces").toTensor().cpu();
    auto acc = forces.accessor<float,2>();
    for(int iatom=0;iatom<system.natoms;iatom++)
    {
        system.forces[iatom][0]=acc[iatom][0];
        system.forces[iatom][1]=acc[iatom][1];
        system.forces[iatom][2]=acc[iatom][2];
    }
}


void MACEForceField::report_profile()
{
    if(!enable_profile || forward_times.empty())
    {
        return;
    }


    double sum = 0;
    double min_time = *std::min_element(
        forward_times.begin(),
        forward_times.end()
    );

    double max_time = *std::max_element(
        forward_times.begin(),
        forward_times.end()
    );


    for(auto t : forward_times)
        sum += t;


    double avg = sum / forward_times.size();


    double variance = 0;

    for(auto t : forward_times)
    {
        variance += (t - avg) * (t - avg);
    }

    variance /= forward_times.size();

    double stddev = std::sqrt(variance);



    std::cout
    <<"========== MACE Forward ==========\n"
    <<"calls : "<<forward_times.size()<<"\n"
    <<"avg   : "<<avg<<" ms\n"
    <<"min   : "<<min_time<<" ms\n"
    <<"max   : "<<max_time<<" ms\n"
    <<"std   : "<<stddev<<" ms\n"
    <<"==================================\n";



    std::ofstream fout("mace_forward_times.dat");

    for(auto t : forward_times)
    {
        fout << t << "\n";
    }

    fout.close();
}


