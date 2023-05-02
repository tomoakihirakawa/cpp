#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

int main() {

   int s = 500;
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
   /* -------------------------------------------------------------------------- */
   // const VV_d A = {{4., -1, 0, -1},
   //                 {-1., 4, -1, 0},
   //                 {0., -1, 4, -1},
   //                 {-1., 0, -1, 4}};
   // const V_d b = {15., 10., 10, 15};
   // // 初期値が大事
   // V_d x0(b.size(), 0.);
   // V_d ans = {6.875, 5.625, 5.625, 6.875};
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
   for (auto restart = 0; restart < 1; ++restart) {
      for (auto i = 200; i <= 1000; i += 100) {
         gmres gm(A, b, x0, i);
         //   std::cout << "gm.x = " << gm.x << std::endl;
         std::cout << "time:" << timer() << std::endl;
         error = gm.err;  // Norm(Dot(A, gm.x) - b);
         x0 = gm.x;
         std::cout << "       Restart : " << restart << std::endl;
         std::cout << "             i : " << i << std::endl;
         std::cout << "estimate error : " << gm.err << std::endl;
         std::cout << "   actual error: " << Norm(Dot(A, x0) - b) << std::endl;
         // std::cout << "      Solution : " << x0 << std::endl;
         std::cout << Red << "--------------------------------" << colorOff << std::endl;
      }
   }
   // }
   std::cout << "time:" << timer() << std::endl;
   // std::cout << "gm->y.size() = " << gm.y.size() << std::endl;
   // std::cout << "gm->y = " << gm.y << std::endl;
   // std::cout << "V = " << gm.ap->V << std::endl;
   // std::cout << "H = " << gm.ap->H << std::endl;
   //    MatrixForm(gm.ap->V);
   //    MatrixForm(gm.ap->H);
};
