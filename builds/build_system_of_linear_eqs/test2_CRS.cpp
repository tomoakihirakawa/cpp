
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
      auto crs = new CRS();
      crs->setIndexCRS(0);
      crs->value = 15.;
      V_CRS.emplace(crs);
   }
   {
      auto crs = new CRS();
      crs->setIndexCRS(1);
      crs->value = 10.;
      V_CRS.emplace(crs);
   }
   {
      auto crs = new CRS();
      crs->setIndexCRS(2);
      crs->value = 10.;
      V_CRS.emplace(crs);
   }
   {
      auto crs = new CRS();
      crs->setIndexCRS(3);
      crs->value = 15.;
      V_CRS.emplace(crs);
   }
   for (const auto& crs0 : V_CRS)
      for (const auto& crs1 : V_CRS) {
         int i = crs0->getIndexCRS();
         int j = crs1->getIndexCRS();
         if (A[i][j] != static_cast<double>(0.))
            crs0->increment(crs1, A[i][j]);
      }

   // std::cout << " Dot : " << Dot(V_CRS, V_CRS) << std::endl;
   std::cout << " Dot : " << Dot(V_CRS) << std::endl;
   std::cout << " Dot : " << Dot(A, b) << std::endl;
};
