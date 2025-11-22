#include <iostream>
#include <vector>
#include "basic_vectors.hpp"

int main() {
   // Define a 3x3 matrix
   {
      T3Tddd A = {{{1.0, 2.0, 3.0},
                   {2.0, 4.0, 6.0},
                   //  {4.0, 2.0, 6.0},
                   {7.0, 2.0, 9.0}}};

      // Define a vector b (right-hand side of the equation Ax = b)
      Tddd b = {1.0, 1.0, 1.0}, x;
      // Create an instance of lapack_svd
      {
         lapack_svd svd(A, x, b);
         std::cout << "lapack_svd svd(A, x, b)" << std::endl;
         std::cout << "x = " << x << std::endl;
         std::cout << "A.x = " << Dot(A, x) << std::endl;
      }
      {
         lapack_svd svd(x, A, b);
         std::cout << "lapack_svd svd(x, A, b)" << std::endl;
         std::cout << "x = " << x << std::endl;
         std::cout << "x.A = " << Dot(x, A) << std::endl;
      }
      {
         lapack_svd_solve(A, x, b);
         std::cout << "lapack_svd_solve(A, x, b)" << std::endl;
         std::cout << "x = " << x << std::endl;
         std::cout << "A.x = " << Dot(A, x) << std::endl;
      }
   }
   {
      std::vector<std::vector<double>> A = {{1.0, 2.0, 3.0},
                                            {2.0, 4.0, 6.0},
                                            {7.0, 2.0, 9.0}};

      // Define a vector b (right-hand side of the equation Ax = b)
      std::vector<double> b = {1.0, 1.0, 1.0}, x(3);
      // Create an instance of lapack_svd
      {
         lapack_svd svd(A, x, b);
         std::cout << "x = " << x << std::endl;
         std::cout << "A.x = " << Dot(A, x) << std::endl;
      }
      {
         lapack_svd svd(x, A, b);
         std::cout << "x = " << x << std::endl;
         std::cout << "x.A = " << Dot(x, A) << std::endl;
      }
   }
   return 0;
}
