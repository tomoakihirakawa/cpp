#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"

typedef std::vector<double> V_d;
typedef std::vector<V_d> VV_d;

int main() {
   // Define a matrix A
   VV_d A = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};

   // Initialize QR object with the matrix A
   QR qr(A);

   // The QR decomposition of A is now available in qr.Q and qr.R
   // Now we will calculate the eigenvalues of A using the QR algorithm
   // This is a basic implementation and may not converge for some matrices.
   VV_d Ak = A;
   for (int i = 0; i < 1; ++i) {
      std::cout << "A: \n";
      MatrixForm(Ak, std::setw(10));
      QR qr_iter(Ak);
      //   qr_iter.R[0][0] *= -1;
      //   qr_iter.R[1][1] *= -1;
      //   qr_iter.Q[0][0] *= -1;
      //   qr_iter.Q[1][1] *= -1;

      Ak = Dot(qr_iter.R, qr_iter.Q);  // Assuming the '*' operator multiplies matrices
      std::cout << "i: " << i << "\n";
   }
}
