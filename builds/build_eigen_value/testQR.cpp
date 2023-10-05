#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"

typedef std::vector<double> V_d;
typedef std::vector<std::vector<double>> VV_d;

/*DOC_EXTRACT 0_1_eigen_value

\insert{QR_decomposition}

*/

// compare with testQR.nb

int main() {

   VV_d A = {{1., 3., 5., 7.},
             {2., 6., 10., 14.},
             {3., 9., 15., 21.},
             {4., 12., 20., 28.},
             {5., 15., 25., 35.}};

   //    A = Transpose(A);

   QR qr(A);

   std::cout << " ------------------- \n";
   std::cout << "A: \n";
   std::cout << MatrixForm(A, 3, 10) << std::endl;
   std::cout << "Q: \n";
   std::cout << MatrixForm(qr.Q, 3, 10) << std::endl;
   std::cout << "R: \n";
   std::cout << MatrixForm(qr.R, 3, 10) << std::endl;
   std::cout << "QR: \n";
   std::cout << std::setw(5) << Dot(qr.Q, qr.R) << std::endl;
   return 0;
}