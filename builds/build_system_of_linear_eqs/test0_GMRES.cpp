/*DOC_EXTRACT solve_linear_systems0

# 連立一次方程式の解法

\insert{ArnoldiProcess}

\insert{GMRES}

* GMRESは反復的な方法で，特に大規模で疎な非対称行列の線形システムを解くのに適している．
* GMRESは一般的に共役勾配法よりも柔軟性があり，非対称行列に対しても使用できる．ただし，反復の回数が増えると計算コストが大きくなる可能性がある．

### テスト

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test0_GMRES.cpp
make
./test0_GMRES
```

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

// #define USE_RANDOM_MATRIX
#define USE_PREDEFINED_MATRIX

int main() {

#if defined(USE_RANDOM_MATRIX)
   int s = 100;
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
   /*
   A = Table[N@RandomInteger[{0, 10}], {i, 1, 10}, {j, 1, 10}]
   b = Table[N@RandomInteger[{0, 10}], {i, 1, 10}]
   Inverse[A].b
   */

   const VV_d A = {{6., 8., 7., 10., 0., 6., 1., 2., 4., 10.},
                   {6., 10., 6., 3., 1., 8., 5., 3., 1., 1.},
                   {9., 10., 8., 1., 7., 9., 5., 0., 9., 5.},
                   {6., 10., 2., 9., 0., 7., 5., 5., 8., 8.},
                   {0., 8., 2., 10., 9., 1., 3., 0., 0., 5.},
                   {3., 3., 6., 9., 0., 9., 10., 7., 9., 3.},
                   {2., 6., 3., 10., 5., 1., 0., 4., 0., 7.},
                   {4., 3., 2., 4., 4., 7., 8., 10., 1., 10.},
                   {7., 1., 3., 0., 6., 9., 7., 1., 5., 0.},
                   {0., 0., 2., 3., 4., 3., 9., 6., 3., 10.}};

   const V_d b = {1., 7., 2., 2., 7., 8., 10., 0., 6., 5.};
   V_d x0(b.size(), 0. /*initial value*/);
   const V_d ans = {6.41725, -0.576444, 1.52755, 1.5933, -0.436665, -7.40106, 4.20872, 0.458307, -0.723366, -1.73434};

   /* ------------------------------------------------------------ */
   /*                          CRS行列の作成                         */
   /* ------------------------------------------------------------ */

   std::cout << "make CRS" << std::endl;

   std::vector<CRS*> A_CRS(A.size(), nullptr);
   for (size_t i = 0; i < A.size(); ++i) {
      A_CRS[i] = new CRS();
      A_CRS[i]->setIndexCRS(i);
      A_CRS[i]->value = x0[i];
   }

   for (size_t i = 0; i < A.size(); ++i)
      for (size_t j = 0; j < A[i].size(); ++j)
         A_CRS[i]->set(A_CRS[j], A[i][j]);

   for (auto& crs : A_CRS)
      crs->setVectorCRS();

   auto A_operator = [&](const V_d& v) -> V_d {
      return Dot(A_CRS, v);
   };

   /* ---------------------------------------------------------- */

#endif

   TimeWatch timer;
   std::cout << "time:" << timer() << std::endl;
   bool finished = false;
   double error;
   int n_max = 15;
   int n_begin = 5;
   auto x0_for_iterate = x0;
   gmres gm_full(A_operator, b, x0, n_max);
   gmres gm_iterate(A_operator, b, x0, n_begin);
   for (auto i = 0; i < 5; ++i) {
      gm_iterate.Iterate(A_operator);
      // std::cout << "restart " << restart << ", " << Green << "terms in an expansion :" << i << colorReset << std::endl;
      std::cout << "     gm_iterate.V.size() : " << gm_iterate.V.size() << std::endl;
      std::cout << "                    time : " << timer() << std::endl;
      std::cout << "estimate error (iterate/full) : " << gm_iterate.err << " / " << gm_full.err << std::endl;
      auto error = Norm(Dot(A, gm_iterate.x) - b);
      std::cout << "  actual error (iterate/full) : " << (error < 1E-10 ? Green : Red) << error << colorReset << " / " << Norm(Dot(A, gm_full.x) - b) << std::endl;
      std::cout << Red << "--------------------------------" << colorReset << std::endl;
   }
   std::cout << "time:" << timer() << std::endl;
};
