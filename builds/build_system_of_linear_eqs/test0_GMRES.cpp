/*DOC_EXTRACT solve_linear_systems0

# 連立一次方程式の解法

\insert{ArnoldiProcess}

\insert{GMRES}

### テスト

<details>
<summary>HOW TO USE</summary>

![](./WATCHME.gif)

</details>

*/
#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

#define USE_RANDOM_MATRIX
// #define USE_PREDEFINED_MATRIX

int main() {

#if defined(USE_RANDOM_MATRIX)
   int s = 50;
   VV_d A(s, V_d(s));
   V_d b(s);
   #pragma omp parallel
   for (auto i = 0; i < s; ++i)
   #pragma omp single nowait
   {
      b[i] = RandomReal({-1., 1.});
      for (auto j = 0; j < s; ++j)
         A[i][j] = RandomReal({-1., 1.});
   }
   V_d x0(b.size(), 0.);
#elif defined(USE_PREDEFINED_MATRIX)
   const VV_d A = {{4., -1, 0, -1},
                   {-1., 4, -1, 0},
                   {0., -1, 4, -1},
                   {-1., 0, -1, 4}};
   const V_d b = {15., 10., 10, 15};
   V_d x0(b.size(), 0. /*initial value*/);
   V_d ans = {6.875, 5.625, 5.625, 6.875};
#endif
   /* -------------------------------------------------------------------------- */
   //    auto v = diagonal_scaling_vector(A);
   //    for (auto i = 0; i < v.size(); ++i) {
   //       A[i] *= v[i];
   //       b[i] *= v[i];
   //    }

   Timer timer;
   std::cout << "time:" << timer() << std::endl;
   bool finished = false;
   double error;
   int n_max = 50;
   int n_begin = 1;
   auto x0_for_iterate = x0;
   for (auto restart = 0; restart < 5; ++restart) {
      gmres gm_ful(A, b, x0, n_max);
      gmres gm_iterate(A, b, x0_for_iterate, n_begin);
      for (auto i = n_begin + 1; i <= n_max; i++) {
         // gmres gm(A, b, x0, i);
         gm_iterate.Iterate(A);
         
         std::cout << "time:" << timer() << std::endl;
         std::cout << "       Restart : " << restart << std::endl;
         std::cout << "             i : " << i << std::endl;
         std::cout << "estimate error (full) : " << gm_ful.err << std::endl;
         std::cout << "  actual error (full) : " << Norm(Dot(A, gm_ful.x) - b) << std::endl;
         //
         std::cout << "estimate error (iterate) : " << gm_iterate.err << std::endl;
         std::cout << "  actual error (iterate) : " << Norm(Dot(A, gm_iterate.x) - b) << std::endl;
         std::cout << Red << "--------------------------------" << colorOff << std::endl;
      }
      x0_for_iterate = gm_iterate.x;
   }
   std::cout << "time:" << timer() << std::endl;
};
