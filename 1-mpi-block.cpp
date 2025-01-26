#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <iomanip>
#include "mpi.h"
#include <chrono>

using namespace std;

//打印
void print_matrix(const float *matrix, int dim) {
  std::cout << std::setfill(' '); 
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      std::cout << std::setw(10) << matrix[i * dim + j]; 
    }
    std::cout << '\n';
  }
}

void parallel_back_substitution(float *matrix, float *b, float *x, int dim, int n_rows, int task_id, int num_tasks, int start_row) {
  // 从最后一行开始回代
  for (int row = dim - 1; row >= 0; row--) {
    x[row] = b[row]; 
            
    // 减去已知x值
    for (int col = dim - 1; col > row; col--) {
      x[row] -= matrix[row * dim + col] * x[col];  
    }  
    // 除以对角线元素
    float pivot = matrix[row * dim + row];
    if (fabs(pivot) > 1e-10) {
      x[row] /= pivot;
    }         
  }
}

//计算任务块大小（矩阵行数可以整除进程数）
const int dim = 2048;
const int num_runs = 5; // 运行次数

int main(int argc, char *argv[]) {
  // MPI初始化
  MPI_Init(&argc, &argv);

  // 获得总进程数
  int num_tasks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

  const int n_rows = dim / num_tasks;

  // 获得进程号
  int task_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

  //计算chunk起始下标
  const int start_row = task_id * n_rows;
  const int end_row = start_row + n_rows;

  float *matrix = nullptr;
  float *x_full = nullptr;

  //存储当前进程的chunk
  float m_chunk[dim*n_rows];

  //缓冲区
  float pivot_row[dim];

  float *b = nullptr;
  float b_chunk[n_rows];

  // 0号进程初始化矩阵
  if (task_id == 0) {
    matrix = new float[dim*dim];
    srand(123);
    for (int i = 0; i < dim*dim; i++)
            matrix[i] = rand() / 100.0f;
    cout << "dim=" << dim << endl;
    // cout << "original matrix:" << endl;
    // print_matrix(matrix, dim);

    b = new float[dim];
    for (int i = 0; i < dim; i++)
        b[i] = rand() / 100.0f;

    x_full = new float[dim];
  }

  // 存储非阻塞通信句柄
  std::vector<MPI_Request> requests(num_tasks);

  // 计时变量
  std::chrono::duration<double> total_time(0);
  for (int run = 0; run < num_runs; ++run) {
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    
    // 分散chunk
    MPI_Scatter(matrix, dim * n_rows, MPI_FLOAT, m_chunk,
                dim * n_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    MPI_Scatter(b, n_rows, MPI_FLOAT, b_chunk, n_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //高斯消元
    for (int row = 0; row < end_row; row++) {
        int mapped_rank = row / n_rows;

        if (task_id == mapped_rank) {
            int local_row = row % n_rows;
            float pivot = m_chunk[local_row * dim + row];

            for (int col = row; col < dim; col++) {
                m_chunk[local_row * dim + col] /= pivot;
            }
            b_chunk[local_row] /= pivot;

            // 创建新的请求数组来分别处理矩阵和b向量的发送
            std::vector<MPI_Request> matrix_requests(num_tasks);
            std::vector<MPI_Request> b_requests(num_tasks);

            // 发送更新后的行和对应的b值
            for (int i = mapped_rank + 1; i < num_tasks; i++) {
                MPI_Isend(m_chunk + dim * local_row, dim, MPI_FLOAT, i, 0,
                         MPI_COMM_WORLD, &matrix_requests[i]);
                MPI_Isend(&b_chunk[local_row], 1, MPI_FLOAT, i, 1,
                         MPI_COMM_WORLD, &b_requests[i]);
            }

            for (int elim_row = local_row + 1; elim_row < n_rows; elim_row++) {
                float scale = m_chunk[elim_row * dim + row];
                for (int col = row; col < dim; col++) {
                    m_chunk[elim_row * dim + col] -=
                        m_chunk[local_row * dim + col] * scale;
                }
                b_chunk[elim_row] -= b_chunk[local_row] * scale;
            }

            // 等待所有发送完成
            for (int i = mapped_rank + 1; i < num_tasks; i++) {
                MPI_Wait(&matrix_requests[i], MPI_STATUS_IGNORE);
                MPI_Wait(&b_requests[i], MPI_STATUS_IGNORE);
            }
        } else {
            MPI_Recv(pivot_row, dim, MPI_FLOAT, mapped_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            float pivot_b;
            MPI_Recv(&pivot_b, 1, MPI_FLOAT, mapped_rank, 1, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            for (int elim_row = 0; elim_row < n_rows; elim_row++) {
                float scale = m_chunk[elim_row * dim + row];
                for (int col = row; col < dim; col++) {
                    m_chunk[elim_row * dim + col] -= pivot_row[col] * scale;
                }
                b_chunk[elim_row] -= pivot_b * scale;
            }
        }
    }
    // 收集到0号进程
    MPI_Gather(m_chunk, n_rows * dim, MPI_FLOAT, matrix, n_rows * dim,
               MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(b_chunk, n_rows, MPI_FLOAT, b, n_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if(task_id == 0){
      parallel_back_substitution(matrix, b, x_full, dim, n_rows, task_id, num_tasks, start_row);
    }

    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    total_time += end_time - start_time;
  }

  // 确保在MPI_Finalize之前所有通信都已完成
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Finalize();

  // 计算平均时间
  if (task_id == 0) {
    double average_time = total_time.count() / num_runs;
    cout << "Average execution time over " << num_runs << " runs: " << average_time << " seconds" << endl;
  }
  
  // 释放matrix
  if (task_id == 0 && matrix != nullptr) {
    // cout << "final matrix:" << endl;
    // print_matrix(matrix, dim);

    // cout << "\nSolution vector x:" << endl;
    // for (int i = 0; i < dim; i++) {
    //     cout << setw(10) << x_full[i] << endl;
    // }

    delete[] matrix;
  }

  return 0;
}