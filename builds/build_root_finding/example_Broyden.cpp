/*
ニュートン法で使うヤコビアンなどを別のものに置き換えた場合，準ニュートン法と呼ぶ
*/
#include "minMaxOfFunctions.hpp"

int main() {

   tuple_of<10, double> u = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   tuple_of<10, double> v = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

   auto M = TensorProduct(u, v);
   std::cout << M << std::endl;
   IdentityMatrix(M);
   std::cout << M << std::endl;
   BroydenMethod BBB(u, v);
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
   Vd X = {5., 5.};
   VVd J = {{1., 0.}, {0., 1.}};
   Vd dX = {0.01, 0.01};
   /* -------------------------------------------------------------------------- */
   BroydenMethod BM(X, X + dX);
   for (auto i = 0; i < 20; i++) {
      /*Rank One Update*/
      std::cout << BM.X[0] << ", " << BM.X[1] << "," << F(BM.X[0], BM.X[1]) << "," << Norm(BM.dX) << std::endl;
      Vd gradF_ = {dFdx(BM.X[0] - BM.dX[0], BM.X[1] - BM.dX[1]), dFdy(BM.X[0] - BM.dX[0], BM.X[1] - BM.dX[1])};
      Vd gradF = {dFdx(BM.X[0], BM.X[1]), dFdy(BM.X[0], BM.X[1])};

      BM.update(gradF, gradF_, i < 5 ? 0.01 : 1.);

      if (Norm(BM.dX) < 1E-10) {
         std::cout << BM.X[0] << ", " << BM.X[1] << "," << F(BM.X[0], BM.X[1]) << "," << Norm(BM.dX) << std::endl;
         break;
      }
   }
}