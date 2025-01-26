#!/bin/bash  
  
# 定义要运行的进程数数组  
np_array=(2 3 4 5 6 7 8)  
  
# 遍历数组中的每个进程数  
for np in "${np_array[@]}"  
do  
  echo "Running with $np processes..."  
  # 使用mpirun运行MPI程序，指定进程数  
  mpirun -n $np ./mpi  
done  
  