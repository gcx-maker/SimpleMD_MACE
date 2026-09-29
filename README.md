# SimpleMD_MACE

一个基于 C++ 的轻量级分子动力学（Molecular Dynamics, MD）框架，支持传统经典势函数与机器学习势函数（Machine Learning Potential, MLP）的统一调用。

本项目主要探索：

- C++ 分子动力学程序设计
- MACE 机器学习势函数的 C++ 部署
- LibTorch 调用 TorchScript 模型
- Python/C++ 混合计算接口
- MPI 并行分子动力学计算

目标是构建一个面向 AI for Science 的轻量级、高性能分子模拟框架，实现从 Python 科研工作流向 C++ 高性能计算环境的迁移。


---

# 项目结构


---

# 项目特点

## 1. C++ 分子动力学框架

SimpleMD_MACE 从底层实现分子动力学核心流程：
       初始结构
          |
          v
   Neighbor List构建
          |
          v
    力计算模块
          |
  ----------------
  |              |
  v              v
 LJ            MACE
  |              |
  ----------------
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
- 能量与温度输出


---

# 2. Lennard-Jones 势函数


项目实现经典 Lennard-Jones 势函数：

\[
E(r)=4\epsilon
\left[
(\frac{\sigma}{r})^{12}
-
(\frac{\sigma}{r})^6
\right]
\]


支持：

- 双体相互作用
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

---

# 4. LibTorch C++模型部署


项目采用 TorchScript + LibTorch 实现 MACE 模型部署。

相比传统 Python 调用方式：

优势：

- 无需 Python运行环境
- C++原生集成
- 更容易进行性能优化
- 方便进一步接入CUDA/MPI


---

# 5. Pybind11版本


项目同时提供 Python/C++ 混合版本。

用途：

- 快速验证C++实现
- 与Python科研流程结合
- 与ASE等工具进行接口测试


---

# 6. MPI并行版本


`simplemd_libtorch_mpi` 提供MPI并行MD框架。

目前实现：
- MPI进程管理
- 原子数据通信
- LJ势函数并行计算

并行模式：
         MPI World


    +-------------+

    |             |

  Rank0        Rank1

    |             |

Atom set A    Atom set B


    |             |

    +-------------+

          |

    Parallel Force

          |

    MD Update


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


示例：
Neighbor Build

MACE Forward

Total Step Time


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

Performance:
~3 ns/day


---

# 开发路线


目前完成：

- [x] C++ MD框架
- [x] LJ势函数
- [x] MACE势函数调用
- [x] LibTorch模型部署
- [x] Neighbor List
- [x] Linked-cell优化
- [x] Pybind11接口
- [x] MPI并行框架


未来计划：

- [ ] CUDA Kernel优化
- [ ] MACE算子优化
- [ ] 多GPU MACE推理
- [ ] NCCL通信
- [ ] 大规模分布式AI-MD


---

# 项目背景


近年来机器学习势函数（Machine Learning Potential）能够提供接近第一性原理计算精度，同时显著降低计算成本。

然而目前大量机器学习势函数仍依赖Python生态。


构建面向 AI for Science 的高性能分子模拟基础框架。


---

# Author

gcx-maker


Research interests:

- AI for Science
- Machine Learning Potential
- Molecular Dynamics
- High Performance Computing
- C++ Scientific Computing
