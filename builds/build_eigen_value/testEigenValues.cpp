#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"

typedef std::vector<double> V_d;
typedef std::vector<std::vector<double>> VV_d;

/*DOC_EXTRACT eigen_value

# 固有値問題

## 固有値の計算

行列$`A`$をQR分解$`A=QR`$し，
$`A_k = Q_k^{-1} A Q_k`$の計算を繰り返すことで，
$`A_k`$の対角成分が$`A`$の固有値に収束することを確認する．

わかりやすいように$`\cdot`$で行列の積を表す．

```math
Q_2^{-1} \cdot (Q_1^{-1} \cdot (Q_0^{-1} \cdot (A = Q_0R_0) \cdot Q_0=Q_1R_1) \cdot Q_1=Q_2 \cdot R_2) \cdot R_2
```

*/

int main() {
   VV_d A = {{6., 5., 0.},
             {5., 1., 4.},
             {0., 4., 3.}};

   V_d eigenvalues = {9.84316, 4.02176, -3.86492};

   VV_d Ak = A;

   int maxIter = 1000;
   double tol = 1e-9;

   auto I = A;
   IdentityMatrix(I);

   for (int i = 0; i < maxIter; ++i) {
      QR qr(Ak);
      Ak = Dot(qr.R, qr.Q);
      //   std::cout << MatrixForm(qr.Q, 5, 10) << std::endl;
      //   std::cout << MatrixForm(qr.R, 5, 10) << std::endl;
      //   std::cout << MatrixForm(Dot(qr.Q, qr.R), 5, 10) << std::endl;

      double total = 0;
      for (auto j = 0; j < Ak.size(); ++j)
         total += std::abs(Det(A - Diagonal(Ak)[j] * I));
      std::cout << "Eigenvalue " << Diagonal(Ak) << ", total: " << total << ", " << i << "\n";
      if (std::abs(total) < tol)
         break;
   }

   std::cout << "True eigenvalue: " << eigenvalues << "\n";

   std::array<std::array<double, 3>, 3> B = {{{6., 5., 0.},
                                              {5., 1., 4.},
                                              {0., 4., 3.}}};

   std::cout << "eigenvalue using function: " << Eigenvalues(B) << "\n";

   // Eigensystem returns pair of eigenvalues and eigenvectors
   std::cout << "eigenvalue using function: " << Eigensystem(B).first << "\n";
   std::cout << "eigenvector using function: " << Eigensystem(B).second << "\n";

   return 0;
}