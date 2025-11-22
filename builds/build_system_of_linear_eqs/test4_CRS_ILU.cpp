/*DOC_EXTRACT 0_0_0_CRS

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test4_CRS_ILU.cpp
make
./test4_CRS_ILU
```

*/

#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

using V_d = std::vector<double>;
const VV_d A = {{8., 16., 24., 32}, {2., 7., 12., 17.}, {6., 17., 32., 59.}, {7., 22., 46., 105.}};

int main() {

  std::vector<CRS *> A_CRS(A.size());

#pragma omp parallel for
  for (size_t i = 0; i < A.size(); ++i) {
    A_CRS[i] = new CRS();
    A_CRS[i]->setIndexCRS(i);
    A_CRS[i]->value = 0.;
  }

  /* ------------------------------- 直接LU分解をしてみる ------------------------------ */

  std::vector<std::vector<double>> L(A.size(), V_d(A.size(), 0.));
  std::vector<std::vector<double>> U(A.size(), V_d(A.size(), 0.));

  for (size_t i = 0; i < A.size(); ++i)
    U[i][i] = 1.;
  for (size_t i = 0; i < A.size(); ++i) {
    for (size_t j = 0; j < A.size(); ++j) {
      double sum = 0.;
      if (i >= j) {
        // Lの計算
        for (size_t k = 0; k <= i - 1 && k <= j; ++k)
          sum += L[i][k] * U[k][j];
        L[i][j] = A[i][j] - sum;
      } else {
        // Uの計算
        for (size_t k = 0; k <= i - 1 && k <= j; ++k)
          sum += L[i][k] * U[k][j];
        U[i][j] = (A[i][j] - sum) / L[i][i];
      }
    }
  }

  std::cout << "L:" << std::endl;
  for (const auto &row : L) {
    for (const auto &val : row)
      std::cout << std::setw(10) << val << " ";
    std::cout << std::endl;
  }

  std::cout << "U:" << std::endl;
  for (const auto &row : U) {
    for (const auto &val : row)
      std::cout << std::setw(10) << val << " ";
    std::cout << std::endl;
  }

  std::cout << "Reconstructed A = L * U :" << std::endl;
  double RMSE = 0.;
  std::vector<std::vector<double>> A_reconstructed = Dot(L, U);
  for (size_t row = 0; row < A_reconstructed.size(); ++row) {
    for (size_t col = 0; col < A_reconstructed[row].size(); ++col) {
      double val = A_reconstructed[row][col];
      RMSE += (val - A[row][col]) * (val - A[row][col]);
      std::cout << std::setw(10) << val << " ";
    }
    std::cout << std::endl;
  }

  RMSE = std::sqrt(RMSE / (A.size() * A[0].size()));
  std::cout << "RMSE between A and L*U: " << RMSE << std::endl;

  /* ------------------------------------------------------------ */
  /*                          CRS行列の作成                       */
  /* ------------------------------------------------------------ */

  // A を CRS 化（test3_CRS.cpp と同じ流儀）
  std::cout << "make CRS(A)" << std::endl;
  const double eps = 0.0; // 必要なら疎化の閾値を設定
  for (size_t i = 0; i < A.size(); ++i) {
    A_CRS[i]->setIndexCRS(i);
    A_CRS[i]->value = 0.0;
  }

  for (size_t i = 0; i < A.size(); ++i) {
    for (size_t j = 0; j < A.size(); ++j) {
      const double aij = A[i][j];
      if (std::abs(aij) > eps) {
        A_CRS[i]->set(A_CRS[j], aij);
      }
    }
  }

  // L, U も CRS 化（必要に応じて）
  std::vector<CRS *> L_CRS(A.size()), U_CRS(A.size());
  for (size_t i = 0; i < A.size(); ++i) {
    L_CRS[i] = new CRS();
    L_CRS[i]->setIndexCRS(i);
    L_CRS[i]->value = 0.0;
    U_CRS[i] = new CRS();
    U_CRS[i]->setIndexCRS(i);
    U_CRS[i]->value = 0.0;
  }
}