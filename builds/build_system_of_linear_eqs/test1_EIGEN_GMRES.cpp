/*DOC_EXTRACT solve_linear_systems0

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test1_EIGEN_GMRES.cpp
make
./test1_EIGEN_GMRES
```

EigenのGMRESを使った結果と比較．

*/

#include <iostream>
#include "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense"
#include "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/IterativeLinearSolvers"
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
   for (int max_iter = 1; max_iter <= 6; ++max_iter) {
      // Eigen's GMRES solver
      x = x0;

      // 初期推定を毎回 x0 に戻す（比較のため）
      x_EIGEN = initEigenVector(x0);

      Eigen::GMRES<Eigen::MatrixXd, Eigen::IdentityPreconditioner> GMRES_EIGEN(A_EIGEN);
      GMRES_EIGEN.setTolerance(tolerance);
      GMRES_EIGEN.setMaxIterations(max_iter);

      // Solve using Eigen's GMRES
      x_EIGEN = GMRES_EIGEN.solveWithGuess(b_EIGEN, x_EIGEN);

      std::cout << Magenta << "Max iterations: " << max_iter << colorReset;
      /* -------------------------------------------------------------------------- */
      std::cout << Yellow << width << "\nEigen's GMRES results:\n";
      std::cout << yellow << width << "Number of iterations: " << GMRES_EIGEN.iterations() << colorReset << std::endl;
      std::cout << yellow << width << "Estimated error: " << prec << GMRES_EIGEN.error() << colorReset << std::endl;
      std::cout << yellow << width << "Solution: " << prec << x_EIGEN.transpose() << colorReset << std::endl;
      std::cout << yellow << width << "Actual error: " << prec << (A_EIGEN * x_EIGEN - b_EIGEN).norm() << colorReset << std::endl;

      // 自作 GMRES
      auto A_mult = [](const V_d& v) -> V_d {
         V_d r(A.size(), 0.0);
         for (std::size_t i = 0; i < A.size(); ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < A[i].size(); ++j) sum += A[i][j] * v[j];
            r[i] = sum;
         }
         return r;
      };

      gmres GMRES(A_mult, b, x, max_iter);

      x = GMRES.x;

      const double abs_res_eigen = Norm(Dot(A, std::vector<double>(x_EIGEN.data(), x_EIGEN.data() + x_EIGEN.size())) - b);
      const double rel_res_eigen = abs_res_eigen / Norm(b);
      const double abs_res_custom = Norm(Dot(A, x) - b);
      const double rel_res_custom = abs_res_custom / Norm(b);

      std::cout << Green << width << "\nCustom GMRES results:\n";
      std::cout << green << width << "Estimated error (as implemented): " << prec << GMRES.err << colorReset << std::endl;
      std::cout << green << width << "Abs residual: " << prec << abs_res_custom << colorReset << std::endl;
      std::cout << green << width << "Rel residual: " << prec << rel_res_custom << colorReset << std::endl;
   }

   return 0;
}
