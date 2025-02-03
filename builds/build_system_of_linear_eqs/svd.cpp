#include <iostream>
#include <vector>
#include "basic_vectors.hpp"

int main() {

   /*
   ランク欠損の行列の例
   */
   //    std::vector<std::vector<double>> A = {
   //        {1.0, 2.0, 3.0},
   //        {2.0, 3.0, 4.0},
   //        {7.0, 8.0, 9.0}};

   std::vector<std::vector<double>> A = {
       {0.0, 2.0, 3.0},
       {0.0, 3.0, 4.0},
       {0.0, 8.0, 9.0}};

   // Define a vector b (right-hand side of the equation Ax = b)
   std::vector<double> b = {1.0, 1.0, 1.0};

   auto solve_LU = [&](const std::vector<std::vector<double>> &A, const std::vector<double> &b) -> std::vector<double> {
      std::vector<double> x(3);
      lapack_lu lu(A);
      lu.solve(b, x);
      return x;
   };

   auto solve_SVD = [&](const std::vector<std::vector<double>> &A, const std::vector<double> &b) -> std::vector<double> {
      std::vector<double> x(3);
      lapack_svd svd(A);
      svd.solve(b, x);
      return x;
   };

   auto x_svd = solve_SVD(A, b);
   auto x_lu = solve_LU(A, b);

   // Output the solution
   std::cout << "SVD Solution x:" << std::endl;
   for (double val : x_svd)
      std::cout << val << " ";

   std::cout << std::endl;

   std::cout << "LU Solution x:" << std::endl;
   for (double val : x_lu)
      std::cout << val << " ";

   std::cout << std::endl;

   return 0;
}
