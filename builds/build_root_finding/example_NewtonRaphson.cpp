/*DOC_EXTRACT
## ヘッセ行列を利用したニュートン法
**最適か否かを判断するための関数**は１つだけで，**最適化したい変数は複数**である場合でも，
最適化は，ヘッセ行列を利用したニュートン法によって可能である．
この方法で，変数は，関数を根とするのではなく，関数を最大最小（停留点）とする値へと収束する．
*/

#include "rootFinding.hpp"

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
   V_d X = {5., 5.};  // optimization variables
   NewtonRaphson<V_d> nr(X);
   std::cout << w << "iteration" << w << "x" << w << "y" << w << "|F(x,y)|" << w << "|dX|" << std::endl;
   for (auto i = 0; i < 50; i++) {
      double x = nr.X[0];
      double y = nr.X[1];
      V_d dFdX = {dFdx(x, y), dFdy(x, y)};
      VV_d Hessian = {{dFdxx(x, y), dFdxy(x, y)}, {dFdxy(x, y), dFdyy(x, y)}};
      nr.update(dFdX, Hessian);
      std::cout << w << i << w << x << w << y << w << F(x, y) << w << Norm(nr.dX) << std::endl;
   }
}
