/*DOC_EXTRACT solve_linear_systems1

\insert{compressed_row_storage}

*/

#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

int main() {

   /* -------------------------------------------------------------------------- */
   const VV_d A = {{4., -1, 0, -1},
                   {-1., 4, -1, 0},
                   {0., -1, 4, -1},
                   {-1., 0, -1, 4}};
   const V_d b = {15., 10., 10, 15};
   // 初期値が大事
   V_d x0(b.size(), 0.);
   V_d ans = {6.875, 5.625, 5.625, 6.875};
   /* -------------------------------------------------------------------------- */

   std::unordered_set<CRS*> V_CRS;

   {
      auto CRS = new CRS();
      CRS->setIndexCRS(0);
      CRS->value = 15.;
      V_CRS.emplace(CRS);
   }
   {
      auto CRS = new CRS();
      CRS->setIndexCRS(1);
      CRS->value = 10.;
      V_CRS.emplace(CRS);
   }
   {
      auto CRS = new CRS();
      CRS->setIndexCRS(2);
      CRS->value = 10.;
      V_CRS.emplace(CRS);
   }
   {
      auto CRS = new CRS();
      CRS->setIndexCRS(3);
      CRS->value = 15.;
      V_CRS.emplace(CRS);
   }
   for (const auto& crs0 : V_CRS)
      for (const auto& crs1 : V_CRS) {
         int i = crs0->getIndexCRS();
         int j = crs1->getIndexCRS();
         if (A[i][j] != static_cast<double>(0.))
            crs0->increment(crs1, A[i][j]);
      }

   std::cout << " Dot : " << Dot(V_CRS, V_CRS) << std::endl;
   std::cout << " Dot : " << Dot(V_CRS) << std::endl;
   std::cout << " Dot : " << Dot(A, b) << std::endl;

   Timer timer;
   std::cout << "time:" << timer() << std::endl;
   bool finished = false;
   double error;
   for (auto restart = 0; restart < 1; ++restart) {
      for (auto i = 1; i <= 10; ++i) {
         gmres gm(V_CRS, b, x0, i);
         //   std::cout << "gm.x = " << gm.x << std::endl;
         std::cout << "time:" << timer() << std::endl;
         error = gm.err;  // Norm(Dot(A, gm.x) - b);
         x0 = gm.x;
         std::cout << "       Restart : " << restart << std::endl;
         std::cout << "             i : " << i << std::endl;
         std::cout << "estimate error : " << gm.err << std::endl;
         std::cout << "   actual error: " << Norm(Dot(A, x0) - b) << std::endl;
         std::cout << "      Solution : " << x0 << std::endl;
         std::cout << Red << "--------------------------------" << colorOff << std::endl;
      }
   }
   std::cout << "time:" << timer() << std::endl;
};
