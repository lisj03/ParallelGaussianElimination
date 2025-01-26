#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <chrono>

using namespace std;

//测试参数设置
int dims[] ={8};// 测试的矩阵维度列表
const int num_runs = 1; // 运行次数

// 打印矩阵
void print_matrix(const float *matrix, int dim) {
  std::cout << std::setfill(' '); 
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      std::cout << std::setw(10) << matrix[i * dim + j]; 
    }
    std::cout << '\n';
  }
  cout << endl;
}

// 回代求解
void back_substitution(float *m, float *b, float *x, int dim) {
    for (int i = dim - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < dim; j++) {
            x[i] -= m[i * dim + j] * x[j];
        }
        x[i] /= m[i * dim + i];
    }
}


void gaussian_elimination(float *m, float *b, int dim) {
    for (int row = 0; row < dim; row++) {
        float pivot = m[row * dim + row];
        
        // 处理矩阵行
        for (int col = row; col < dim; col++) {
            m[row * dim + col] /= pivot;
        }
        // 处理b向量
        b[row] /= pivot;

        for (int elim_row = row + 1; elim_row < dim; elim_row++) {
            float scale = m[elim_row * dim + row];
            
            // 消元矩阵行
            for (int col = row; col < dim; col++) {
                m[elim_row * dim + col] -= m[row * dim + col] * scale;
            }
            // 消元b向量
            b[elim_row] -= b[row] * scale;
        }
    }
}

int main() {

  for (int dim : dims) {
    cout << "Testing dimension: " << dim << endl;

    // 生成随机矩阵和b向量
    float* m = new float[dim * dim];
    float* b = new float[dim];  // 方程组右边的常数项
    float* x = new float[dim];  // 解向量
    
    srand(123);
    for (int i = 0; i < dim * dim; i++)
        m[i] = rand() / 100.0f;
    for (int i = 0; i < dim; i++)
        b[i] = rand() / 100.0f;

    // cout << "Original matrix:" << endl;
    // print_matrix(m, dim);
    // cout << "Original b vector:" << endl;
    // for (int i = 0; i < dim; i++)
    //     cout << setw(10) << b[i] << endl;

    std::chrono::duration<double> total_time(0);

    for (int run = 0; run < num_runs; ++run) {
        std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

        // 高斯消元
        gaussian_elimination(m, b, dim);
        // 回代求解
        back_substitution(m, b, x, dim);

        std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
        total_time += end_time - start_time;
    }

    double average_time = total_time.count() / num_runs;
    cout << "Average execution time over " << num_runs << " runs: " << average_time << " seconds" << endl;


    // cout << "final matrix:" << endl;
    // print_matrix(m, dim);
    cout << "\nSolution vector x:" << endl;
    for (int i = 0; i < dim; i++) {
        cout << setw(10) << x[i] << endl;
    }

    // 释放内存
    delete[] m;
    delete[] b;
    delete[] x;
  }
  

  return 0;
}