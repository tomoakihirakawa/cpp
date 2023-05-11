/*DOC_EXTRACT newton
## 準ニュートン法
ニュートン法で使うヤコビアンなどを別のものに置き換えた方法．
*/
#include "minMaxOfFunctions.hpp"

// Himmelblau's function
auto F(double x, double y) { return std::pow(x * x + y - 11, 2) + std::pow(x + y * y - 7, 2); };
auto dFdx(double x, double y) { return 2 * (x * x + y - 11) * 2 * x + 2 * (x + y * y - 7); };
auto dFdy(double x, double y) { return 2 * (x * x + y - 11) + 2 * (x + y * y - 7) * 2 * y; };
auto dFdxy(double x, double y) { return 2 * 1 * 2 * x + 2 * 2 * y; };
auto dFdyx(double x, double y) { return dFdxy(x, y); };
auto dFdxx(double x, double y) { return 2 * 2 * x * 2 * x + 2 * 1; };
auto dFdyy(double x, double y) { return 2 * 1 + 2 * 2 * y * 2 * y; };
auto w = std::setw(20);

int main() {

   // std::array<double, 2> X, dX, gradF, gradF_;
   Vd X, dX, gradF, gradF_;
   X = {5., 5.};  // optimization variables
   dX = {5., 5.};
   BroydenMethod BM(X, X + dX);
   std::cout << w << "iteration" << w << "x" << w << "y" << w << "|F(x,y)|" << w << "|dX|" << std::endl;
   for (auto i = 0; i < 50; i++) {
      auto x = BM.X[0];  // optimization variable
      auto y = BM.X[1];  // optimization variable
      auto dx = BM.dX[0];
      auto dy = BM.dX[1];
      std::cout << w << i << w << x << w << y << w << F(x, y) << w << Norm(BM.dX) << std::endl;
      gradF_ = {dFdx(x - dx, y - dy), dFdy(x - dx, y - dy)};
      gradF = {dFdx(x, y), dFdy(x, y)};
      BM.update(gradF, gradF_, i < 5 ? 0.01 : 1.);
   }
}