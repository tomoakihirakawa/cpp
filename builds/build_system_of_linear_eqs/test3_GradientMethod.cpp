/*DOC_EXTRACT solve_linear_systems0

## 共役勾配法と勾配降下法

### 共役勾配法（Conjugate Gradient, CG）

* 共役勾配法は反復的な方法で，特に大規模で疎な（つまり，ほとんどの要素がゼロである）対称正定値行列の線形システムを解くのに適している．
* この方法の利点は，一般的に反復回数が行列の次元に対して比較的少ないこと．しかし，非対称または非正定値の行列には適用できない．

### 勾配降下法 (Gradient Descent, GD)

* 勾配降下法は最も基本的な最適化アルゴリズムで，線形システムまたは一般的な最適化問題を解くことができる．
* しかし，勾配降下法の収束速度は通常比較的遅く，特に凸でない問題に対しては局所最小値に陥る可能性がある．

*/
#include "lib_measurement.hpp"
//
#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

// #define USE_RANDOM_MATRIX
#define USE_PREDEFINED_MATRIX

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

   // Initialize x
   V_d x = x0;

   // Create an instance of the GradientMethod class
   GradientMethod solver(A);

   // Optionally enable diagonal scaling
   // solver.diagonal_scaling();

   // Measure the time and solve the system using the Conjugate Gradient method
   Timer timer;
   std::cout << "CG method:" << std::endl;
   std::cout << "Start time:" << timer() << std::endl;
   V_d x_CG = solver.solveCG(b, x);
   std::cout << "Error " << Norm(Dot(A, x_CG) - b) << std::endl;
   std::cout << "End time:" << timer() << std::endl;

   // Reset x
   x = x0;

   // Measure the time and solve the system using the Gradient Descent method
   std::cout << "GD method:" << std::endl;
   V_d x_GD = solver.solve(b, x);
   std::cout << "Error " << Norm(Dot(A, x_GD) - b) << std::endl;
   std::cout << "End time:" << timer() << std::endl;

   return 0;
}
