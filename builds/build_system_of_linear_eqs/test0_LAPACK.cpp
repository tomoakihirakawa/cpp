/*DOC_EXTRACT solve_linear_systems0

## LU分解(LAPACK)

* LU分解は直接的な方法で，あらゆる種類の行列（対称、非対称、正定値、非正定値）に適用できる．
* この方法は反復的な方法よりも計算コストが高くなる可能性があるが，反復法とは異なり，収束性の問題がない．

*/
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"

int main() {

   //    int s = 500;
   //    VV_d A(s, V_d(s));
   //    V_d b(s);
   // #pragma omp parallel
   //    for (auto i = 0; i < s; ++i)
   // #pragma omp single nowait
   //    {
   //       b[i] = RandomReal({-1., 1.});
   //       for (auto j = 0; j < s; ++j)
   //          A[i][j] = RandomReal({-1., 1.});
   //    }
   //    V_d x0(b.size(), 0.);
   /* -------------------------------------------------------------------------- */
   const VV_d A = {{4., -1, 0, -1},
                   {-1., 4, -1, 0},
                   {0., -1, 4, -1},
                   {-1., 0, -1, 4}};
   const V_d b = {15., 10., 10, 15};
   const VV_d inverse = {{0.2916666667, 0.08333333333, 0.04166666667, 0.08333333333},
                         {0.08333333333, 0.2916666667, 0.08333333333, 0.04166666667},
                         {0.04166666667, 0.08333333333, 0.2916666667, 0.08333333333},
                         {0.08333333333, 0.04166666667, 0.08333333333, 0.2916666667}};
   V_d x0(b.size(), 0. /*initial value*/);
   V_d ans = {6.875, 5.625, 5.625, 6.875};
   /* -------------------------------------------------------------------------- */
   //    auto v = diagonal_scaling_vector(A);
   //    for (auto i = 0; i < v.size(); ++i) {
   //       A[i] *= v[i];
   //       b[i] *= v[i];
   //    }

   Timer timer;
   Print("LAPACK LU time:", timer());
   lapack_lu lu(A);
   lu.solve(b, x0);
   Print("LAPACK LU error : ", Norm(Dot(A, x0) - b));
   Print("LAPACK LU time : ", timer());
   Print("LAPACK LU inverse : ", lu.inverse() - inverse);

   lapack_svd svd(A);
   svd.solve(b, x0);
   Print("LAPACK SVD error:", Norm(Dot(A, x0) - b));
   Print("LAPACK SVD time:", timer());
};
