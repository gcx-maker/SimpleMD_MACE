// #ifndef MACEDATABUILDER_H
// #define MACEDATABUILDER_H

// #include <torch/torch.h>

// #include <unordered_map>
// #include <string>
// #include <vector>

// #include "system.h"
// #include "neighbor.h"

// class MACEDataBuilder
// {
// public:
//     MACEDataBuilder(const std::vector<int>& z_table, torch::Device device);
//     void initialize(const System& system);
//     void build(const System& system, const Neighbor& neighbor);
//     const std::unordered_map<std::string, torch::Tensor>& data() const;
// private:
//     void build_node_attrs(const System& system);
//     void build_positions(const System& system);
//     void build_edge_index(const Neighbor& neighbor);
//     void build_shifts(const Neighbor& neighbor);
//     void build_batch(const System& system);
//     void build_cell(const System& system);
//     void build_ptr(const System& system);
// private:
//     std::vector<int> z_table;
//     torch::Device device;
//     torch::Dtype dtype;
//     std::unordered_map<std::string, torch::Tensor> mace_data;

//     torch::Tensor node_attrs_tensor;
//     torch::Tensor batch_tensor;
//     torch::Tensor ptr_tensor;

//     torch::Tensor positions_tensor;
//     torch::Tensor positions_cpu_tensor;
//     torch::Tensor edge_index_tensor;
//     torch::Tensor shifts_tensor;
//     torch::Tensor cell_tensor;
// };
// #endif


#pragma once

#include <torch/script.h>
#include <unordered_map>
#include <vector>

#include "system.h"
#include "neighbor.h"


class MACEDataBuilder
{

public:

    MACEDataBuilder(
        const std::vector<int>& z_table,
        torch::Device device
    );


    void initialize(
        const System& system
    );


    void build(
        const System& system,
        const Neighbor& neighbor
    );


    const std::unordered_map<std::string,torch::Tensor>&
    data() const;


public:

    void build_node_attrs(
        const System& system
    );


    void build_positions(
        const System& system
    );


    void build_graph(
        const Neighbor& neighbor,
        const System& system
    );


    void build_cell(
        const System& system
    );


private:


    std::vector<int> z_table;

    torch::Device device;

    torch::Dtype dtype;


    /*
        固定tensor
    */

    torch::Tensor node_attrs_tensor;

    torch::Tensor batch_tensor;

    torch::Tensor ptr_tensor;


    /*
        position buffer
    */

    torch::Tensor positions_cpu_tensor;

    torch::Tensor positions_tensor;



    /*
        graph buffer
    */

    torch::Tensor edge_index_cpu_tensor;

    torch::Tensor edge_index_tensor;


    torch::Tensor shifts_cpu_tensor;

    torch::Tensor shifts_tensor;


    int max_edges;



    /*
        cell
    */

    torch::Tensor cell_tensor;



    std::unordered_map<std::string,torch::Tensor>
    mace_data;

};