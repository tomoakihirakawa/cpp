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
      // std::cout << " i = " << i << ", x = " << xn << ", dx = " << dx << std::endl;
   }
   /* ----------------------- 多変数の場合 ----------------------- */
   {
      /*
      最適化する変数は複数あって，
      最適かどうかを判断する関数は１つの場合，
      ヤコビアンは作れないので，ヤコビアンを使うニュートン法は使えない．
      代わりにヘッシアンを使うニュートン法は使うことができる．
      または，最急降下法（gradient descent）を使うことができる．
       */
      auto F = [](double x, double y) { return std::pow(x * x + y - 11, 2) + std::pow(x + y * y - 7, 2); };
      auto dFdx = [](double x, double y) { return 2 * (x * x + y - 11) * 2 * x + 2 * (x + y * y - 7); };
      auto dFdy = [](double x, double y) { return 2 * (x * x + y - 11) + 2 * (x + y * y - 7) * 2 * y; };
      double x = -1., y = 5.;
      for (auto i = 0; i < 20; i++) {
         std::cout << x << "," << y << "," << F(x, y) << std::endl;
         auto dx = -F(x, y) / dFdx(x, y);
         auto dy = -F(x, y) / dFdy(x, y);
         x += dx;
         y += dy;
         if (Norm(Tdd{x, y}) < 1E-10) {
            std::cout << x << "," << y << "," << F(x, y) << std::endl;
            break;
         }
      }
   }
   /* ----------------------- 停留点に変数を収束させる ----------------------- */
   {
      /*
      最適化する変数は複数あって，
      最適かどうかを判断する関数は１つの場合でも，
      ヘッシアンを使うニュートン法は使うことができる．
      関数の最大最小のような停留点に変数は収束することになる．
      */
      auto F = [](double x, double y) { return std::pow(x * x + y - 11, 2) + std::pow(x + y * y - 7, 2); };
      auto dFdx = [](double x, double y) { return 2 * (x * x + y - 11) * 2 * x + 2 * (x + y * y - 7); };
      auto dFdy = [](double x, double y) { return 2 * (x * x + y - 11) + 2 * (x + y * y - 7) * 2 * y; };
      auto dFdxy = [](double x, double y) { return 2 * 1 * 2 * x + 2 * 2 * y; };
      auto dFdyx = [&](double x, double y) { return dFdxy(x, y); };
      auto dFdxx = [](double x, double y) { return 2 * 2 * x * 2 * x + 2 * 1; };
      auto dFdyy = [](double x, double y) { return 2 * 1 + 2 * 2 * y * 2 * y; };
      double x = -1., y = 5.;
      V_d X = {x, y}, dX(2);
      for (auto i = 0; i < 20; i++) {
         // std::cout << x << "," << y << "," << F(x, y) << std::endl;
         VV_d H = {{dFdxx(x, y), dFdxy(x, y)}, {dFdyx(x, y), dFdyy(x, y)}};
         V_d gradF = {dFdx(x, y), dFdy(x, y)};
         ludcmp lu(H);
         lu.solve(gradF, dX);
         X -= dX;
         x = X[0];
         y = X[1];
         if (Norm(dX) < 1E-10) {
            // std::cout << x << "," << y << "," << F(x, y) << std::endl;
            //
            break;
         }
      }
   }
   /* ----------------------- 多次元の場合 ----------------------- */
   V_d Xn = {1., 1., 1.};
   V_d Xn1(3), dX(3);
   for (auto i = 0; i < 20; i++) {
      ludcmp lu(dFdx(Xn));
      lu.solve(-F(Xn), dX);
      Xn1 = Xn + dX;
      Xn = Xn1;
      // std::cout << " i = " << i << ", Xn = " << Xn << ", dx = " << dX << std::endl;
   }
};
