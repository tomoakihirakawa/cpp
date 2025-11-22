#include "interpolations.hpp"

int main() {

   std::array<std::array<double, 3>, 3> X012 = {{{0, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
   std::array<double, 3> COM = {1. / 3., 1. / 3., 0};
   std::array<double, 3> P012 = {1, 1, 1};

   auto intpP = interpolationTriangleLinear0101(P012);
   auto intpX = interpolationTriangleLinear0101(X012);
   auto n = TriangleNormal(X012);
   double f;
   Tddd force = {0., 0., 0.};
   Tddd torque = {0., 0., 0.};
   for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple) {
      f = intpP(x0, x1) * intpX.J(x0, x1) * w0w1;
      force += f * n;
      torque += Cross(intpX(x0, x1) - COM, f * n);
   }

   std::cout << "force = " << force << std::endl;
   std::cout << "torque = " << torque << std::endl;
};