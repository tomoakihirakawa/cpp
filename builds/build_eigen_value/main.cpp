#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"

typedef std::vector<double> V_d;
typedef std::vector<std::vector<double>> VV_d;

// int main() {
//    // Define a matrix A
//    VV_d A = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};

//    // Initialize QR object with the matrix A
//    QR qr(A);

//    // The QR decomposition of A is now available in qr.Q and qr.R
//    // Now we will calculate the eigenvalues of A using the QR algorithm
//    // This is a basic implementation and may not converge for some matrices.
//    VV_d Ak = A;
//    for (int i = 0; i < 1; ++i) {
//       std::cout << "A: \n";
//       MatrixForm(Ak, std::setw(10));
//       QR qr_iter(Ak);
//       //   qr_iter.R[0][0] *= -1;
//       //   qr_iter.R[1][1] *= -1;
//       //   qr_iter.Q[0][0] *= -1;
//       //   qr_iter.Q[1][1] *= -1;

//       Ak = Dot(qr_iter.R, qr_iter.Q);  // Assuming the '*' operator multiplies matrices
//       std::cout << "i: " << i << "\n";
//    }
// }

int main() {
   VV_d A = {{6., 5., 0.}, {5., 1., 4.}, {0., 4., 3.}};
   VV_d Ak = A;

   int maxIter = 1;
   double tol = 1e-9;

   auto I = A;
   IdentityMatrix(I);

   for (int i = 0; i < maxIter; ++i) {
      QR qr(Ak);
      Ak = Dot(qr.R, qr.Q);
      std::cout << MatrixForm(qr.Q, std::setw(10)) << std::endl;
      std::cout << MatrixForm(qr.R, std::setw(10)) << std::endl;
      std::cout << MatrixForm(Dot(qr.Q, qr.R), std::setw(10)) << std::endl;

      std::cout << "Iteration: " << i << "\n";
      std::cout << Det(Ak - Dot(Diagonal(Ak), I)) << std::endl;
   }

   // Eigenvalues should be on the diagonal of Ak
   for (int i = 0; i < Ak.size(); ++i) {
      std::cout << "Eigenvalue " << i + 1 << ": " << Ak[i][i] << "\n";
   }

   return 0;
}