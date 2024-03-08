/*DOC_EXTRACT 0_0_0_Newton

# ニュートン法

## ニュートン法

**最適か否かを判断するための関数**，**ゼロまたは最大か最小にしたい関数**を目的関数と呼ぶ．
ヤコビ行列をやヘッセ行列，その両方を使う場合がある．

目的関数の根を見つける場合は，ヤコビ行列を使う．
最適化の問題の多くは，目的関数の最大最小を求めることなので，ヘッセ行列を利用したニュートン法を用いる．

```bash
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_NewtonRaphson_0.cpp
make
./example0_NewtonRaphson_0
```

*/

#include "minMaxOfFunctions.hpp"
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

   for (double x = -5; x <= 5; x++)
      for (double y = -5; y <= 5; y++) {
         std::ofstream ofs_NM("./output_Himmelblau/Himmelblau_Newton" + std::to_string((int)x) + "_" + std::to_string((int)y) + ".csv");
         // ofs_NM << "#x y F(x, y)" << std::endl;
         V_d X = {x, y};  // optimization variables
         NewtonRaphson<V_d> nr(X);
         std::cout << w << "iteration" << w << "x" << w << "y" << w << "|F(x,y)|" << w << "|dX|" << std::endl;
         for (auto i = 0; i < 50; i++) {
            double x = nr.X[0];
            double y = nr.X[1];
            V_d gradF = {dFdx(x, y), dFdy(x, y)};
            VV_d Hessian = {{dFdxx(x, y), dFdxy(x, y)}, {dFdxy(x, y), dFdyy(x, y)}};

            nr.update(gradF, Hessian, 1.);
            ofs_NM << x << "," << y << "," << F(x, y) << std::endl;
            std::cout << w << i << w << x << w << y << w << F(x, y) << w << Norm(nr.dX) << std::endl;
            if (Norm(nr.dX) < 1e-10)
               break;
         }
         ofs_NM.close();
      }
}

// int main() {
//    V_d X = {5., 5.};  // optimization variables
//    GradientMethod gd;
//    gd.initialize(X);
//    std::cout << w << "iteration" << w << "x" << w << "y" << w << "|F(x,y)|" << w << "|dX|" << std::endl;
//    for (auto i = 0; i < 50; i++) {
//       double x = gd.X[0];
//       double y = gd.X[1];
//       V_d dFdX = {dFdx(x, y), dFdy(x, y)};
//       gd.update(dFdX, 0.01);
//       std::cout << w << i << w << x << w << y << w << F(x, y) << w << Norm(gd.dX) << std::endl;
//    }
// }

// int main() {
//    V_d X = {5., 5.};  // optimization variables
//    GradientMethod gd;
//    gd.initialize(X);
//    std::cout << w << "iteration" << w << "x" << w << "y" << w << "|F(x,y)|" << w << "|dX|" << std::endl;
//    for (auto i = 0; i < 50; i++) {
//       double x = gd.X[0];
//       double y = gd.X[1];
//       V_d dFdX = {dFdx(x, y), dFdy(x, y)};
//       auto CD = gd.getFletcherReevesCD(dFdX);
//       gd.update(CD * 0.1);
//       std::cout << w << i << w << x << w << y << w << F(x, y) << w << Norm(gd.dX) << std::endl;
//    }
// }
