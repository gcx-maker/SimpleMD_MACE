// #include "macedatabuilder.h"

// #include <algorithm>
// #include <stdexcept>


// MACEDataBuilder::MACEDataBuilder(const std::vector<int>& z_table, torch::Device device)
// : z_table(z_table), device(device), dtype(torch::kFloat32)
// {

// }


// void MACEDataBuilder::initialize(const System& system)
// {
//     int N = system.natoms;

//     node_attrs_tensor =torch::zeros({N,(int)z_table.size()},
//         torch::TensorOptions().dtype(torch::kFloat32).device(device));


//     batch_tensor =torch::zeros({N},
//         torch::TensorOptions().dtype(torch::kInt64).device(device));


//     ptr_tensor =torch::tensor({0,N},
//         torch::TensorOptions().dtype(torch::kInt64).device(device));


//     positions_tensor =torch::empty({N,3},
//         torch::TensorOptions().dtype(torch::kFloat32).device(device));

//     positions_cpu_tensor =torch::empty({N,3},
//         torch::TensorOptions().dtype(dtype).device(torch::kCPU));

//     cell_tensor =torch::zeros({1,3,3},
//         torch::TensorOptions().dtype(torch::kFloat32).device(device));
//     build_node_attrs(system);
//     build_batch(system);
//     build_ptr(system);
// }


// void MACEDataBuilder::build_node_attrs(const System& system)
// {
//     auto tensor = torch::zeros({(int64_t)system.natoms, (int64_t)z_table.size()}
//     ,torch::TensorOptions().dtype(dtype).device(torch::kCPU));

//     auto accessor = tensor.accessor<float,2>();

//     for (int iatom=0; iatom<system.natoms; iatom++)
//     {
//         int z = system.atomic_numbers[iatom];
//         auto iter = std::find(z_table.begin(), z_table.end(), z);
//         if (iter == z_table.end())
//         {
//             throw std::runtime_error("Atomic number not found in z_table:" + std::to_string(z));
//         }
//         int index = std::distance(z_table.begin(), iter);
//         accessor[iatom][index] = 1.0;
//     }
//     tensor = tensor.to(device);
//     mace_data["node_attrs"] = tensor;
// }


// void MACEDataBuilder::build_batch(const System& system)
// {
//     mace_data["batch"] = torch::zeros({system.natoms},
//         torch::TensorOptions().dtype(torch::kInt64).device(device));
// }

// void MACEDataBuilder::build_positions(const System& system)
// {
//     // auto tensor = torch::empty({system.natoms, 3},
//     //     torch::TensorOptions().dtype(dtype).device(torch::kCPU).requires_grad(true));
    
//     auto accessor = tensor.accessor<float,2>();

//     for (int iatom=0; iatom<system.natoms; iatom++)
//     {
//         accessor[iatom][0] = system.positions[iatom][0];
//         accessor[iatom][1] = system.positions[iatom][1];
//         accessor[iatom][2] = system.positions[iatom][2];
//     }
//     tensor = tensor.to(device);
//     tensor.requires_grad_(true);
//     mace_data["positions"] = tensor;
// }


// void MACEDataBuilder::build_edge_index(const Neighbor& neighbor)
// {
//     const auto& point = neighbor.point_list();
//     const auto& list = neighbor.neighbor_list();

//     int n_edges = list.size();

//     auto tensor = torch::empty({2, n_edges}, 
//         torch::TensorOptions().dtype(torch::kInt64).device(torch::kCPU));
    
//     auto accessor = tensor.accessor<int64_t, 2>();

//     int count = 0;

//     for (int i=0; i<point.size()-1; i++)
//     {
//         int start = point[i];
//         int end = point[i+1];
//         for (int k=start; k<end; k++)
//         {
//             accessor[0][count]=i;
//             accessor[1][count]=list[k];
//             count++;
//         }
//     }
//     tensor = tensor.to(device);
//     mace_data["edge_index"] = tensor;
// }


// void MACEDataBuilder::build_shifts(const Neighbor& neighbor)
// {
//     const auto& shifts = neighbor.neighbor_shifts();

//     int n_edges = shifts.size();


//     auto tensor = torch::zeros({n_edges,3},
//         torch::TensorOptions().dtype(dtype).device(torch::kCPU));


//     auto accessor = tensor.accessor<float,2>();

//     for(int i=0;i<n_edges;i++)
//     {
//         accessor[i][0]=shifts[i][0];
//         accessor[i][1]=shifts[i][1];
//         accessor[i][2]=shifts[i][2];
//     }
//     tensor = tensor.to(device);
//     mace_data["shifts"]=tensor;
// }


// void MACEDataBuilder::build_cell(const System& system)
// {
//     auto tensor = torch::zeros({1, 3, 3}, 
//         torch::TensorOptions().dtype(dtype).device(torch::kCPU));
    
//     auto accessor = tensor.accessor<float, 3>();

//     accessor[0][0][0] = system.cell[0];
//     accessor[0][1][1] = system.cell[1];
//     accessor[0][2][2] = system.cell[2];
    
//     tensor = tensor.to(device);
//     mace_data["cell"] = tensor;
// }


// void MACEDataBuilder::build(const System& system, const Neighbor& neighbor)
// {
//     build_positions(system);
//     build_edge_index(neighbor);
//     build_shifts(neighbor);
//     build_cell(system);
// }


// void MACEDataBuilder::build_ptr(const System& system)
// {
//     auto tensor =torch::tensor({0,system.natoms},
//             torch::TensorOptions().dtype(torch::kInt64).device(device));

//     mace_data["ptr"] = tensor;
// }


// const std::unordered_map<std::string,torch::Tensor>&
// MACEDataBuilder::data() const
// {

//     return mace_data;

// }

#include "macedatabuilder.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>



MACEDataBuilder::MACEDataBuilder(
    const std::vector<int>& z_table,
    torch::Device device
)
:
z_table(z_table),
device(device),
dtype(torch::kFloat32),
max_edges(10000)
{

}




void MACEDataBuilder::initialize(
    const System& system
)
{

    int N = system.natoms;


    /*
        node attrs
    */

    node_attrs_tensor =
        torch::zeros(
            {N,(int)z_table.size()},
            torch::TensorOptions()
            .dtype(dtype)
            .device(device)
        );



    /*
        batch
    */

    batch_tensor =
        torch::zeros(
            {N},
            torch::TensorOptions()
            .dtype(torch::kInt64)
            .device(device)
        );



    /*
        ptr
    */

    ptr_tensor =
        torch::tensor(
            {0,N},
            torch::TensorOptions()
            .dtype(torch::kInt64)
            .device(device)
        );



    /*
        position buffer
    */


    positions_cpu_tensor =
        torch::empty(
            {N,3},
            torch::TensorOptions()
            .dtype(dtype)
            .device(torch::kCPU)
        );


    positions_tensor =
        torch::empty(
            {N,3},
            torch::TensorOptions()
            .dtype(dtype)
            .device(device)
        );


    /*
        graph buffer
    */


    edge_index_cpu_tensor =
        torch::empty(
            {2,max_edges},
            torch::TensorOptions()
            .dtype(torch::kInt64)
            .device(torch::kCPU)
        );


    edge_index_tensor =
        torch::empty(
            {2,max_edges},
            torch::TensorOptions()
            .dtype(torch::kInt64)
            .device(device)
        );



    shifts_cpu_tensor =
        torch::empty(
            {max_edges,3},
            torch::TensorOptions()
            .dtype(dtype)
            .device(torch::kCPU)
        );


    shifts_tensor =
        torch::empty(
            {max_edges,3},
            torch::TensorOptions()
            .dtype(dtype)
            .device(device)
        );



    /*
        cell
    */


    cell_tensor =
        torch::zeros(
            {1,3,3},
            torch::TensorOptions()
            .dtype(dtype)
            .device(device)
        );



    build_node_attrs(system);



    /*
        固定数据
    */

    mace_data["node_attrs"]
        =
        node_attrs_tensor;


    mace_data["batch"]
        =
        batch_tensor;


    mace_data["ptr"]
        =
        ptr_tensor;


    mace_data["positions"]
        =
        positions_tensor;


    mace_data["cell"]
        =
        cell_tensor;


}






void MACEDataBuilder::build_node_attrs(
    const System& system
)
{

    auto cpu_tensor =
    torch::zeros(
        {system.natoms,(int)z_table.size()},
        torch::TensorOptions()
        .dtype(dtype)
        .device(torch::kCPU)
    );


    auto acc =
        cpu_tensor.accessor<float,2>();


    for(int i=0;i<system.natoms;i++)
    {

        int z=system.atomic_numbers[i];

        auto iter =
            std::find(
                z_table.begin(),
                z_table.end(),
                z
            );


        int index =
            std::distance(
                z_table.begin(),
                iter
            );


        acc[i][index]=1.0f;
    }


    node_attrs_tensor.copy_(cpu_tensor);

}






void MACEDataBuilder::build_positions(
    const System& system
)
{

    auto acc =
        positions_cpu_tensor.accessor<float,2>();


    for(int i=0;i<system.natoms;i++)
    {

        acc[i][0]=system.positions[i][0];
        acc[i][1]=system.positions[i][1];
        acc[i][2]=system.positions[i][2];

    }


    positions_tensor.copy_(
        positions_cpu_tensor
    );


    auto pos =
        positions_tensor.detach();


    pos.requires_grad_(true);


    mace_data["positions"]
        =
        pos;

}







void MACEDataBuilder::build_graph(
    const Neighbor& neighbor,
    const System& system
)
{


    const auto& point =
        neighbor.point_list();


    const auto& list =
        neighbor.neighbor_list();


    const auto& shifts =
        neighbor.neighbor_shifts();



    int n_edges =
        list.size();



    if(n_edges > max_edges)
    {

        max_edges =
            n_edges * 2;



        edge_index_cpu_tensor =
            torch::empty(
                {2,max_edges},
                torch::TensorOptions()
                .dtype(torch::kInt64)
                .device(torch::kCPU)
            );


        edge_index_tensor =
            torch::empty(
                {2,max_edges},
                torch::TensorOptions()
                .dtype(torch::kInt64)
                .device(device)
            );



        shifts_cpu_tensor =
            torch::empty(
                {max_edges,3},
                torch::TensorOptions()
                .dtype(dtype)
                .device(torch::kCPU)
            );



        shifts_tensor =
            torch::empty(
                {max_edges,3},
                torch::TensorOptions()
                .dtype(dtype)
                .device(device)
            );

    }



    auto edge_acc =
        edge_index_cpu_tensor.accessor<int64_t,2>();


    auto shift_acc =
        shifts_cpu_tensor.accessor<float,2>();



    int count=0;


    // for(int i=0;i<point.size()-1;i++)
    // {

    //     for(int k=point[i];
    //         k<point[i+1];
    //         k++)
    //     {


    //         edge_acc[0][count]=i;

    //         edge_acc[1][count]=list[k];


    //         shift_acc[count][0]=shifts[k][0];

    //         shift_acc[count][1]=shifts[k][1];

    //         shift_acc[count][2]=shifts[k][2];


    //         count++;

    //     }

    for(int i=0;i<point.size()-1;i++)
    {

        for(int k=point[i];
            k<point[i+1];
            k++)
        {


            int j = list[k];


            PLMD::Vector dr =system.positions[j] - system.positions[i] + shifts[k];

            double r2 = dotProduct(dr,dr);

            if(r2 <= neighbor.cutoff2)
            {

                edge_acc[0][count]=i;
                edge_acc[1][count]=j;

                shift_acc[count][0]=shifts[k][0];
                shift_acc[count][1]=shifts[k][1];
                shift_acc[count][2]=shifts[k][2];

                count++;

            }
        }
    }

    edge_index_tensor
        .slice(1,0,count)
        .copy_(
            edge_index_cpu_tensor
            .slice(1,0,count)
        );



    shifts_tensor
        .slice(0,0,count)
        .copy_(
            shifts_cpu_tensor
            .slice(0,0,count)
        );



    mace_data["edge_index"]
        =
        edge_index_tensor
        .slice(1,0,count)
        .contiguous();



    mace_data["shifts"]
        =
        shifts_tensor
        .slice(0,0,count)
        .contiguous();
}



void MACEDataBuilder::build_cell(
    const System& system
)
{

    auto cpu_cell =
    torch::zeros(
    {1,3,3},
    torch::TensorOptions()
    .dtype(dtype)
    .device(torch::kCPU)
    );


    auto acc =
    cpu_cell.accessor<float,3>();


    acc[0][0][0]=system.cell[0];
    acc[0][1][1]=system.cell[1];
    acc[0][2][2]=system.cell[2];


    cell_tensor.copy_(cpu_cell);

}








void MACEDataBuilder::build(
    const System& system,
    const Neighbor& neighbor
)
{

    build_positions(system);

    build_graph(neighbor, system);

    build_cell(system);

}







const std::unordered_map<std::string,torch::Tensor>&
MACEDataBuilder::data() const
{

    return mace_data;

}