#include <iostream>
#include <vector>
#include "basic_vectors.hpp"

int main() {
   // Define a 3x3 matrix
   std::vector<std::vector<double>> A = {
       {1.0, 2.0, 3.0},
       {4.0, 5.0, 6.0},
       {7.0, 8.0, 9.0}};

   // Define a vector b (right-hand side of the equation Ax = b)
   std::vector<double> b = {1.0, 0.0, 0.0};

   // Create an instance of lapack_svd
   lapack_svd svd(A);

   // Vector to hold the solution x
   std::vector<double> x;

   // Solve the linear system using SVD
   svd.solve(b, x);

   // Output the solution
   std::cout << "Solution x:" << std::endl;
   for (double val : x) {
      std::cout << val << " ";
   }
   std::cout << std::endl;

   return 0;
}
