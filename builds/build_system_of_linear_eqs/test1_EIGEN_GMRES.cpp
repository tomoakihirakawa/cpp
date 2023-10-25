
/*DOC_EXTRACT solve_linear_systems0

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test1_EIGEN_GMRES.cpp
make
./test1_EIGEN_GMRES
```

EigenのGMRESを使った結果と比較．

*/

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include "EIGEN_GMRES.hpp"
#include "basic_linear_systems.hpp"

// Define the matrix and vectors
using VV_d = std::vector<std::vector<double>>;
using V_d = std::vector<double>;

const VV_d A = {{4., -1, 0, -1},
                {-1., 4, -1, 0},
                {0., -1, 4, -1},
                {-1., 0, -1, 4}};

const V_d b = {15., 10., 10, 15};
V_d x0(b.size(), 0.);
const V_d solution = {6.875, 5.625, 5.625, 6.875};

// Function to initialize Eigen Matrix from VV_d
Eigen::MatrixXd initEigenMatrix(const VV_d& matrix) {
   Eigen::MatrixXd eigenMatrix(matrix.size(), matrix[0].size());
   for (size_t i = 0; i < matrix.size(); ++i) {
      for (size_t j = 0; j < matrix[0].size(); ++j) {
         eigenMatrix(i, j) = matrix[i][j];
      }
   }
   return eigenMatrix;
}

// Function to initialize Eigen Vector from V_d
Eigen::VectorXd initEigenVector(const V_d& vec) {
   Eigen::VectorXd eigenVector(vec.size());
   for (size_t i = 0; i < vec.size(); ++i) {
      eigenVector(i) = vec[i];
   }
   return eigenVector;
}

int main() {
   // Initialize Eigen matrix and vectors
   Eigen::MatrixXd A_EIGEN = initEigenMatrix(A);
   Eigen::VectorXd b_EIGEN = initEigenVector(b);
   Eigen::VectorXd x_EIGEN(x0.size());

   double tolerance = 1e-12;
   auto x = x0;
   auto width = std::setw(25);
   auto prec = std::setprecision(16);
   for (int max_iter = 1; max_iter <= 3; ++max_iter) {
      // Eigen's GMRES solver
      Eigen::GMRES<Eigen::MatrixXd, Eigen::IdentityPreconditioner> GMRES_EIGEN(A_EIGEN);
      GMRES_EIGEN.setTolerance(tolerance);
      GMRES_EIGEN.setMaxIterations(max_iter);

      // Solve using Eigen's GMRES
      x_EIGEN = GMRES_EIGEN.solveWithGuess(b_EIGEN, x_EIGEN);

      std::cout << Magenta << "Max iterations: " << max_iter << colorOff << std::endl;
      /* -------------------------------------------------------------------------- */
      std::cout << Yellow << width << "\nEigen's GMRES results:\n";
      std::cout << yellow << width << "Number of iterations: " << GMRES_EIGEN.iterations() << colorOff << std::endl;
      std::cout << yellow << width << "Estimated error: " << prec << GMRES_EIGEN.error() << colorOff << std::endl;
      std::cout << yellow << width << "Solution: " << prec << x_EIGEN.transpose() << colorOff << std::endl;
      std::cout << yellow << width << "Actual error: " << prec << (A_EIGEN * x_EIGEN - b_EIGEN).norm() << colorOff << std::endl;

      // Your GMRES solver
      gmres GMRES(A, b, x, 1);
      int count = 1;
      for (int i = 1; i < max_iter; ++i) {
         GMRES.Iterate(A);
         count++;
      }

      x = GMRES.x;
      std::cout << Green << width << "\nCustom GMRES results:\n";
      std::cout << green << width << "Number of iterations: " << count << colorOff << std::endl;
      std::cout << green << width << "Estimated error: " << prec << GMRES.err << colorOff << std::endl;
      std::cout << green << width << "Solution: " << prec << x << colorOff << std::endl;
      std::cout << green << width << "Actual error: " << prec << Norm(Dot(A, x) - b) << colorOff << std::endl;
   }

   return 0;
}
