
/*
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=study_cblas_performance.cpp
*/
// #include <cblas.h>
#include <Accelerate/Accelerate.h>
#include <chrono>
#include <iostream>
#include <vector>

void blasMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec, std::vector<double>& result) {
   cblas_dgemv(CblasRowMajor, CblasNoTrans, matrix.size(), vec.size(), 1.0, &matrix[0][0], vec.size(), &vec[0], 1, 0.0, &result[0], 1);
}

void simpleMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec, std::vector<double>& result) {
   for (size_t i = 0; i < matrix.size(); ++i) {
      result[i] = 0;
      for (size_t j = 0; j < vec.size(); ++j) {
         result[i] += matrix[i][j] * vec[j];
      }
   }
}

int main() {
   std::vector<std::vector<double>> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
   std::vector<double> vec = {1, 2, 3};
   std::vector<double> result(matrix.size());

   auto start = std::chrono::high_resolution_clock::now();
   blasMultiply(matrix, vec, result);
   auto end = std::chrono::high_resolution_clock::now();
   std::cout << "BLAS Result: ";
   for (const auto& i : result) {
      std::cout << i << " ";
   }
   std::cout << "\nTime taken by BLAS: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";

   start = std::chrono::high_resolution_clock::now();
   simpleMultiply(matrix, vec, result);
   end = std::chrono::high_resolution_clock::now();
   std::cout << "Simple Multiply Result: ";
   for (const auto& i : result) {
      std::cout << i << " ";
   }
   std::cout << "\nTime taken by Simple Multiply: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";

   return 0;
}
