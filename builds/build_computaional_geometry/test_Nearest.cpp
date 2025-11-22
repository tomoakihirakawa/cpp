#include <array>
#include <iostream>
#include <vector>
#include "basic_geometry.hpp"

int main() {
   // 三角形の頂点を定義
   std::array<std::array<double, 3>, 3> triangle = {
       std::array<double, 3>{0.0, 0.0, 0.4},
       std::array<double, 3>{0.2, 0.0, 0.0},
       std::array<double, 3>{0.0, 0.0, 0.0}};

   // 各テスト点について、三角形上の最寄り点を計算
   std::cout << std::fixed << std::setprecision(5);  // 追加
   std::vector<std::array<std::array<double, 3>, 2>> list;
   std::vector<std::array<double, 3>> start_x;
   int n = 20, m = 30;
   for (auto i = 0; i < n; i++) {
      for (auto j = 0; j < m; j++) {
         auto r = std::sin(i * M_PI / n);
         std::array<double, 3> X = {r * std::cos(j * 2 * M_PI / m), r * std::sin(j * 2 * M_PI / m), -1. + 2. * i / (n - 1)};
         auto Y = Nearest(X, triangle);
         std::cout << X[0] << ", " << X[1] << ", " << X[2] << ", " << Y[0] << ", " << Y[1] << ", " << Y[2] << std::endl;
      }
   }
}
