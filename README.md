# SimpleMD_MACE

一个面向机器学习势函数部署的 C++ 分子动力学框架。

本项目探索如何将基于 PyTorch 训练的机器学习势函数（以 MACE 为例）
从 Python 科研环境迁移到 C++ 高性能分子动力学模拟流程中。

项目实现了：
- C++ 分子动力学核心框架
- LJ 经典势函数计算
- MACE 模型 LibTorch 推理
- Pybind11 Python/C++ 混合接口
- MPI 并行分子动力学计算
---

# 项目结构
      SimpleMD_MACE
      |
      ├── simplemd_libtorch_mpi
      │
      │ C++ LibTorch版本
      │
      │ - C++ MD核心
      │ - LJ势函数
      │ - MACE模型推理
      │ - Neighbor List
      │ - MPI通信
      │ - LJ并行计算
      │
      │
      │
      ├── simplemd_pybind11
      │
      │ Python/C++混合版本
      │
      │ - Python调用MACE
      │ - C++负责MD积分

---

# 项目特点

## 1. C++ 分子动力学框架

SimpleMD_MACE 从底层实现分子动力学核心流程：
      
       初始结构
          |
          v
      邻居列表构建
          |
          v
      力计算模块(LJ MACE)
          |
          v
      MD时间积分
          |
          v
     新的原子位置

目前支持：
- NVE 系综
- Langevin 温控
- 周期性边界条件
- Verlet Neighbor List
- Linked-cell neighbor search
- 能量与温度轨迹输出
---

# 2. Lennard-Jones 势函数
项目实现经典 Lennard-Jones 势函数：
支持：
- 多组分相互作用
- 周期边界
- Neighbor List优化
- MPI并行计算
输入示例：

      potential LJ
      atomtype 0 Ar 39.948 18
      atomtype 1 Kr 83.798 36
      
      pair_coeff 0 0 3.405 0.0104
      pair_coeff 0 1 3.5145 0.01207
      pair_coeff 1 1 3.624 0.0140

---

# 3. MACE机器学习势函数支持

SimpleMD_MACE 支持直接在 C++ 中调用 MACE 模型。


# 4. LibTorch C++模型部署

项目采用 TorchScript + LibTorch 实现 MACE 模型部署。

相比传统 Python 调用方式：

优势：

- 无需 Python运行环境
- C++原生集成
- 更容易进行性能优化
---

# 5. Pybind11版本

simplemd_pybind11 提供 Python/C++ 混合调用模式。

计算流程：

      Python
       |
       | MACE inference
       |
      Pybind11
       |
      C++ MD Core
       |
      Integrator

该版本主要用于：

- 快速验证C++ MD框架
- 对比Python和C++性能
- 保持Python生态兼容性

---

# 6. MPI并行版本
`simplemd_libtorch_mpi` 提供MPI并行MD框架。

目前实现：
- MPI进程管理
- 原子数据通信
- LJ势函数并行计算

当前MACE推理仍采用单进程模式。

未来计划：

- MACE多GPU推理
- NCCL通信
- 分布式Neighbor构建
- 大规模AI-MD模拟

---

# 输入文件格式
SimpleMD_MACE通过输入文件控制模拟参数。


## LJ 示例
structurefile ./result/crystal_1.xyz

outputfile ./result/output.xyz

temperature 300

timestep 1.0

friction 0.1

cutoff 8.0

skin 1.0

nstep 2000

nconfig 100 ./result/trajectory.xyz

nstat 100 ./result/energies.dat

ensemble NVE

potential LJ

---

## MACE 示例
structurefile ./result/md_ac_data.xyz

outputfile ./result/output.xyz

temperature 300

timestep 0.5

friction 0.10

cutoff 7.0

skin 1.0

potential MACE

model_path ./IL_liquid_MACE_compiled.model

device cpu

enable_profile T


---

# 性能分析

项目支持C++端性能统计：
enable_profile T

可以统计：
- Neighbor List构建时间
- MACE Forward时间
- 单步MD耗时

---

# 性能测试


## MACE推理测试

测试环境：
GPU:
RTX 3080 Ti

Model:
MACE

Backend:
LibTorch CUDA
System:
208 atoms


      MD Performance
      Total time       : 562.106 s
      Steps            : 20000
      Step time        : 28.105299 ms
      Steps/sec        : 35.58
      Simulation speed : 3.0742 ns/day

      Neighbor Build
      calls : 1454
      avg   : 7.28378 ms
      min   : 6.25228 ms
      max   : 9.2989 ms

      MACE Forward
      calls : 20000
      avg   : 18.711 ms
      min   : 17.3808 ms
      max   : 2032.15 ms
      std   : 22.3867 ms

---

# 开发路线

目前完成：

- [√] C++ MD框架
- [√] LJ势函数
- [√] MACE势函数调用
- [√] LibTorch模型部署
- [√] Neighbor List
- [√] Linked-cell优化
- [√] Pybind11接口
- [√] MPI并行框架

未来计划：
- [ ] CUDA Kernel优化
- [ ] MACE算子优化
- [ ] 多GPU MACE推理
- [ ] 大规模分布式MD
---

