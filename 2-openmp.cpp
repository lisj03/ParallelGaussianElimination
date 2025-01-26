#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <omp.h>
#include <chrono>


using namespace std;

// 打印矩阵
void print_matrix(const float *matrix, int dim) {
  std::cout << std::setfill(' '); 
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      std::cout << std::setw(10) << matrix[i * dim + j]; 
    }
    std::cout << '\n';
  }
}

// 回代求解函数
void back_substitution(float *matrix, float *b, float *x, int dim) {
    // 初始化解向量
    for (int i = 0; i < dim; i++) {
        x[i] = 0.0f;
    }

    // 从最后一行开始回代
    for (int row = dim - 1; row >= 0; row--) {
        x[row] = b[row];
        
        // 减去已知x值的影响
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

int main(int argc, char* argv[]) {
    int dims[] = {2048};
    const int num_runs = 5; // 运行次数

    for(int num_threads = 2; num_threads <= 10; num_threads++) {
        cout << "\n=== Testing with " << num_threads << " threads ===" << endl;

        for (int dim : dims) {
            cout << "dim=" << dim << endl;
            float* matrix = new float[dim * dim];
            float* b = new float[dim];  // 添加b向量
            float* x = new float[dim];  // 添加解向量x
            
            // 生成随机矩阵和b向量
            srand(123);
            for (int i = 0; i < dim * dim; i++)
                matrix[i] = rand() / 100.0f;
            for (int i = 0; i < dim; i++)
                b[i] = rand() / 100.0f;

            omp_set_num_threads(num_threads);

            // 计时变量
            std::chrono::duration<double> total_time(0);


            for (int run = 0; run < num_runs; ++run) {
                std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

                #pragma omp parallel
                {
                    for (int i = 0; i < dim; i++) {
                        #pragma omp single
                        {
                            // 取主元并归一化
                            float pivot = matrix[i * dim + i];
                            for (int col = i; col < dim; col++) {
                                matrix[i * dim + col] /= pivot;
                            }
                            b[i] /= pivot;  // 同时更新b向量
                        }

                        // 并行消元
                        #pragma omp for schedule(dynamic)
                        for (int elim_row = i + 1; elim_row < dim; elim_row++) {
                            float scale = matrix[elim_row * dim + i];
                            for (int col = i; col < dim; col++) {
                                matrix[elim_row * dim + col] -= matrix[i * dim + col] * scale;
                            }
                            b[elim_row] -= b[i] * scale;  // 同时更新b向量
                        }
                    }
                }

                // 回代求解（串行执行）
                back_substitution(matrix, b, x, dim);

                std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                total_time += end_time - start_time;
            }

            // 计算平均时间
            double average_time = total_time.count() / num_runs;
            cout << "Average execution time with " << num_threads << " threads over " 
                 << num_runs << " runs: " << average_time << " seconds" << endl;

            // 打印结果（可选，因为每个线程数的结果应该是一样的）
            // cout << "\nSolution vector x:" << endl;
            // for (int i = 0; i < dim; i++) {
            //     cout << setw(10) << x[i] << endl;
            // }

            // 释放内存
            delete[] matrix;
            delete[] b;
            delete[] x;
        }
    }

    return 0;
}