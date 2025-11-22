#include "basic_geometry.hpp"

int main() {
   /* ---------------------------- solid angle test ---------------------------- */
   Tddd org = {0., 0., 0};
   for (auto i = -100; i < 100; i++) {
      Tddd p0 = {-1 + 2., -1., (double)i / 100.};
      Tddd p1 = {1 + 2., -1., (double)i / 100.};
      Tddd p2 = {0. + 2, 1., (double)i / 100.};
      // std::cout << std::setprecision(20) << TetrahedronSolidAngle(org, p0, p1, p2) / M_PI << std::endl;
      auto solidangle0 = TetrahedronSolidAngle_UsingVectorAngle(org, p0, p1, p2);
      auto solidangle1 = TetrahedronSolidAngle(org, p0, p1, p2);
      std::cout << "TetrahedronSolidAngle_UsingVectorAngle " << std::setprecision(10) << Red << solidangle0 << colorOff << std::endl;
      std::cout << "                 TetrahedronSolidAngle " << std::setprecision(10) << Blue << solidangle1 << colorOff << std::endl;
      std::cout << "                                  diff " << Green << (solidangle0 - solidangle1) << colorOff << std::endl;
   }
};
