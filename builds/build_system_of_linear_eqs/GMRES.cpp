#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

int main() {

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

   //    VV_d A = {{0.0247911, 0.161413, 0.625419, 0.465341, 0.794249},
   //              {0.895294, 0.215363, 0.280354, 0.0206005, 0.906597},
   //              {0.457972, 0.76661, 0.590316, 0.535627, 0.00733951},
   //              {0.315392, 0.925959, 0.412796, 0.825637, 0.322538},
   //              {0.572894, 0.0998945, 0.738812, 0.30581, 0.904702}};
   //    V_d b = {1, 2, 3, 4, 5};

   // 初期値が大事
   V_d x0(b.size(), 0.);
   //    auto v = diagonal_scaling_vector(A);
   //    for (auto i = 0; i < v.size(); ++i) {
   //       A[i] *= v[i];
   //       b[i] *= v[i];
   //    }
   Timer timer;
   std::cout << "time:" << timer() << std::endl;
   bool finished = false;
   double error;
   for (auto restart = 0; restart < 3; ++restart) {
      for (auto i = 99; i < 110; i++) {
         gmres gm(A, b, x0, i);
         //   std::cout << "gm.x = " << gm.x << std::endl;
         std::cout << "time:" << timer() << std::endl;
         error = Norm(Dot(A, gm.x) - b);
         Print(i);
         if (error < 1E-10) {
            Print(error, Blue);
            finished = true;
            // break;
         } else {
            Print(error, Green);
            x0 = gm.x;
         }
         std::cout << "restart = " << restart << std::endl;
         std::cout << "i = " << i << std::endl;
         Print(error, Green);
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
