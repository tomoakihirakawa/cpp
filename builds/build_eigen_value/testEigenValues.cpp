#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"

typedef std::vector<double> V_d;
typedef std::vector<std::vector<double>> VV_d;

/*DOC_EXTRACT 0_0_eigen_value

# 固有値問題

## 固有値の計算

行列$`A`$をQR分解$`A=QR`$し，
$`A_k = Q_k^{-1} A Q_k`$の計算を繰り返すことで，
$`A_k`$の対角成分が$`A`$の固有値に収束することを確認する．

わかりやすいように$`\cdot`$で行列の積を表す．

```math
Q_2^{-1} \cdot (Q_1^{-1} \cdot (Q_0^{-1} \cdot (A = Q_0R_0) \cdot Q_0=Q_1R_1) \cdot Q_1=Q_2 \cdot R_2) \cdot R_2
```

テストケース

```Mathematica
In[23]:= A = Table[i*j + i*i + j*j + i + j, {i, 1, 5}, {j, 1, 5}]
Eigenvalues[N@A]

Out[23]= {{5, 10, 17, 26, 37}, {10, 16, 24, 34, 46}, {17, 24, 33, 44,
  57}, {26, 34, 44, 56, 70}, {37, 46, 57, 70, 85}}

Out[24]= {210.296, -15.5105, 0.214605,
 1.64434*10^-14, -1.77456*10^-15}
```


```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release .. -DSOURCE_FILE=testEigenValues.cpp
make
./testEigenValues
```

*/

int main() {
   // VV_d A = {{6., 5., 0.},
   //           {5., 1., 4.},
   //           {0., 4., 3.}};

   // V_d eigenvalues = {9.84316, 4.02176, -3.86492};

   VV_d A = {{5, 10, 17, 26, 37}, {10, 16, 24, 34, 46}, {17, 24, 33, 44, 57}, {26, 34, 44, 56, 70}, {37, 46, 57, 70, 85}};
   std::array<std::array<double, 5>, 5> B = {{{5, 10, 17, 26, 37}, {10, 16, 24, 34, 46}, {17, 24, 33, 44, 57}, {26, 34, 44, 56, 70}, {37, 46, 57, 70, 85}}};
   V_d eigenvalues = {210.296, -15.5105, 0.214605, 1.64434e-14, -1.77456e-15};

   VV_d Ak = A;

   int maxIter = 1000;
   double tol = 1e-9;

   auto I = A;
   IdentityMatrix(I);

   for (int i = 0; i < maxIter; ++i) {
      QR qr(Ak);
      Ak = Dot(qr.R, qr.Q);
      // std::cout << MatrixForm(qr.Q, 5, 10) << std::endl;
      // std::cout << MatrixForm(qr.R, 5, 10) << std::endl;
      // std::cout << MatrixForm(Dot(qr.Q, qr.R), 5, 10) << std::endl;
      // std::cin.ignore();

      double total = 0;
      for (auto j = 0; j < Ak.size(); ++j)
         total += std::abs(Det(A - Diagonal(Ak)[j] * I));
      std::cout << "Eigenvalue " << Diagonal(Ak) << ", total: " << total << ", " << i << "\n";
      if (std::abs(total) < tol)
         break;
   }

   std::cout << "True eigenvalue: " << eigenvalues << "\n";
   std::cout << "eigenvalue using function: " << Eigenvalues(A) << "\n";

   // Eigensystem returns pair of eigenvalues and eigenvectors
   std::cout << "eigenvalue using function: " << Eigensystem(A).first << "\n";
   std::cout << "eigenvector using function: " << Eigensystem(A).second << "\n";

   return 0;
}