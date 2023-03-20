/*
シンプルなニュートン法の例
*/
#include "rootFinding.hpp"

double f(const double x) {
   return M_PI + x + sin(x);
}

double dfdx(const double x) {
   return 1 + cos(x);
}

V_d F(const V_d &X) {
   return {sin(X[0]), cos(X[1]), sin(X[2])};
}
VV_d dFdx(const V_d &X) {
   return {{cos(X[0]), 0, 0},
           {0, -sin(X[1]), 0},
           {0, 0, cos(X[2])}};
}

int main() {
   /* ----------------------- 一次元の場合 ----------------------- */
   double xn = 1., xn1 = 0., dx;
   for (auto i = 0; i < 20; i++) {
      dx = -f(xn) / dfdx(xn);
      xn1 = xn + dx;
      xn = xn1;
      std::cout << " i = " << i << ", x = " << xn << ", dx = " << dx << std::endl;
   }
   std::cout << "---------------------------" << std::endl;
   /* ----------------------- 多次元の場合 ----------------------- */
   V_d Xn = {1., 1., 1.};
   V_d Xn1(3), dX(3);
   for (auto i = 0; i < 20; i++) {
      ludcmp lu(dFdx(Xn));
      lu.solve(-F(Xn), dX);
      Xn1 = Xn + dX;
      Xn = Xn1;
      std::cout << " i = " << i << ", Xn = " << Xn << ", dx = " << dX << std::endl;
   }
};
