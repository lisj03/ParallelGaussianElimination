# ParallelGaussianElimination
The major assignment for the Parallel Computing course at USTB.
# README

### 运行环境

- macbook pro 2021

  芯片：Apple M1 pro

  内存：16GB

  macOS：Sequoia 15.2
- MPICH：4.2.3

### 0-baseline

#### 测试参数

- 向dims数组中添加多个元素，可以实现批处理
- 修改run\_num，可以在多次运行后，计算得到程序的平均单次运行时间

  对于矩阵维度小，运行速度快的情况下，平均时间有效性更好

#### 编译运行

1. 命令行编译（启用c++11）
   ```bash
   g++ -std=c++11 0-baseline.cpp -o baseline
   ```
2. 运行
   ```bash
   ./baseline
   ```

### 1-mpi-block

#### 安装mpi

```bash
brew install mpich
```

#### 代码思路

- 消元过程为块并行
- 回代过程因依赖性强为串行

#### 测试参数

- dim：矩阵维度
- run\_num：运行次数。可以在多次运行后，计算得到程序的平均单次运行时间

  对于矩阵维度小，运行速度快的情况下，平均时间有效性更好

#### 编译运行

1. 编译
   ```c++
   mpicxx -o mpi 1-mpi-block.cpp
   ```
2. 运行
   - 命令行
     ```bash
     mpirun -n 4 ./mpi
     ```
   - shell脚本实现批处理
     ```powershell
     #!/bin/bash  
       
     # 定义要运行的进程数数组  
     np_array=(2 3 4 5 6 7 8)  
       
     # 遍历数组中的每个进程数  
     for np in "${np_array[@]}"  
     do  
       echo "Running with $np processes..."  
       # 使用mpirun运行MPI程序，指定进程数  
       mpirun -n $np ./mpi  
     done   run-mpi.sh
     ```
     ```bash
     ./run-mpi.sh
     ```

#### bug

当dim=2048, 2个进程并行时会报错越界访问的问题

### 2-openmp

#### 安装

```bash
brew install libomp
```

#### 代码思路

- 消元过程为块并行
- 回代过程因依赖性强为串行

#### 测试参数

- 向dims数组中添加多个元素，可以实现批处理
- 修改run\_num，可以在多次运行后，计算得到程序的平均单次运行时间

  对于矩阵维度小，运行速度快的情况下，平均时间有效性更好

#### 编译运行

1. 编译
   ```c++
   clang++ -Xpreprocessor -fopenmp -lomp 2-openmp.cpp -o openmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib
   ```
2. 运行
   ```c++
   ./openmp
   ```
