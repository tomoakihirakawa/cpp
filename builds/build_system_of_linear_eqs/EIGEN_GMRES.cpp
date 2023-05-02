#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
//
#include "EIGEN_GMRES.hpp"

int main() {
   // Define a 4x4 matrix
   Eigen::MatrixXd A(4, 4);
   A << 4, -1, 0, -1,
       -1, 4, -1, 0,
       0, -1, 4, -1,
       -1, 0, -1, 4;

   // Define the right-hand side vector
   Eigen::VectorXd b(4);
   b << 15, 10, 10, 15;

   // Initialize the solution vector
   Eigen::VectorXd x(4);

   // Set the tolerance for the solver
   double tolerance = 1e-6;

   for (int max_iter = 1; max_iter <= 10; ++max_iter) {
      // Set up the GMRES solver with a diagonal preconditioner
      Eigen::GMRES<Eigen::MatrixXd, Eigen::IdentityPreconditioner> gmres_solver(A);
      gmres_solver.setTolerance(tolerance);
      gmres_solver.setMaxIterations(max_iter);

      // Solve the linear system
      x = gmres_solver.solveWithGuess(b, x);

      // Output the results
      std::cout << "Max iterations: " << max_iter << std::endl;
      std::cout << "Number of iterations: " << gmres_solver.iterations() << std::endl;
      std::cout << "Estimated error: " << gmres_solver.error() << std::endl;
      std::cout << "Solution: " << std::endl
                << x << std::endl;
   }

   return 0;
}
