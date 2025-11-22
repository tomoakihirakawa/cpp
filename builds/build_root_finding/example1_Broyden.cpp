/*DOC_EXTRACT 0_1_1_quasiNewton

## 準ニュートン法

\insert{Broydens_Method}

ニュートン法で使うヤコビ行列などを別のものに置き換えた方法．

```bash
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_Broyden.cpp
make
./example1_Broyden
```

![./output_Himmelblau/convergence.gif](./output_Himmelblau/convergence.gif)

基本的に，Newton法の方が収束が早いことがわかる．

Newton法では，勾配ベクトルとヘッセ行列を使ったが，Broyden法では，勾配ベクトルのみを使っている．
勾配ベクトルを異なる点で計算して，ヘッセ行列の近似を行う．
勾配ベクトルがゼロになる点を探すのではなく，目的関数ベクトルがゼロになる点を探すこともできるだろう．


*Himmelblau関数*

$$
f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
$$

*1階偏導関数*

$\frac{\partial f}{\partial x} = 4x(x^2 + y - 11) + 2(x + y^2 - 7)$

$\frac{\partial f}{\partial y} = 2(x^2 + y - 11) + 4y(x + y^2 - 7)$

*2階偏導関数（ヘッセ行列の成分）*

$\frac{\partial^2 f}{\partial x^2} = 12x^2 + 4y - 42$

$\frac{\partial^2 f}{\partial y^2} = 12y^2 + 4x - 26$

$\frac{\partial^2 f}{\partial x \partial y} = \frac{\partial^2 f}{\partial y \partial x} = 4x + 4y$

*/

// \label{quasi_newton:broyden}

#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"

// Himmelblau's function
auto F(double x, double y) { return std::pow(x * x + y - 11, 2) + std::pow(x + y * y - 7, 2); };
auto dFdx(double x, double y) { return 2 * (x * x + y - 11) * 2 * x + 2 * (x + y * y - 7); };
auto dFdy(double x, double y) { return 2 * (x * x + y - 11) + 2 * (x + y * y - 7) * 2 * y; };
auto dFdxy(double x, double y) { return 4 * (x + y); };
auto dFdyx(double x, double y) { return dFdxy(x, y); };
auto dFdxx(double x, double y) { return 12 * x * x + 4 * y - 42; };
auto dFdyy(double x, double y) { return 12 * y * y + 4 * x - 26; };
auto w = std::setw(20);

int main() {

   for (double x = -5; x <= 5; x++)
      for (double y = -5; y <= 5; y++) {

         std::array<double, 2> X, dX, gradF, gradF_;
         // Vd X, dX, gradF, gradF_;
         X = {x, y};  // optimization variables
         dX = {0., 1E-10};
         std::cout << w << "iteration" << w << "x" << w << "y" << w << "|F(x,y)|" << w << "|dX|" << std::endl;
         // to visualize the Himmelblau's function
         std::ofstream ofs("./output_Himmelblau/Himmelblau.csv");
         // ofs << "#x y F(x, y)" << std::endl;
         for (auto i = -100; i < 100; i++) {
            for (auto j = -100; j < 100; j++) {
               double x = 5. * i / 100.;
               double y = 5. * j / 100.;
               ofs << x << "," << y << "," << F(x, y) << std::endl;
            }
         }
         ofs.close();
         // visualize the path of the optimization

         //! Newton-Raphson Method
         {
            std::ofstream ofs_NM("./output_Himmelblau/Himmelblau_Newton" + std::to_string((int)x) + "_" + std::to_string((int)y) + ".csv");
            // ofs_NM << "#x y F(x, y)" << std::endl;
            V_d X = {x, y};  // optimization variables
            NewtonRaphson<V_d> nr(X);
            for (auto i = 0; i < 50; i++) {
               double x = nr.X[0];  // optimization variable
               double y = nr.X[1];  // optimization variable
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

         //! Broyden Method
         {
            std::ofstream ofs_BM("./output_Himmelblau/Himmelblau_Broyden" + std::to_string((int)x) + "_" + std::to_string((int)y) + ".csv");
            BroydenMethod<std::array<double, 2>> BM(X, X + dX);
            // ofs_BM << "#x y F(x, y)" << std::endl;
            for (auto i = 0; i < 50; i++) {
               auto x = BM.X[0];  // optimization variable
               auto y = BM.X[1];  // optimization variable
               auto dx = BM.dX[0];
               auto dy = BM.dX[1];
               std::cout << w << i << w << x << w << y << w << F(x, y) << w << Norm(BM.dX) << std::endl;
               ofs_BM << x << "," << y << "," << F(x, y) << std::endl;
               gradF_ = {dFdx(x - dx, y - dy), dFdy(x - dx, y - dy)};
               gradF = {dFdx(x, y), dFdy(x, y)};

               // gradF_ = {dFdx(x, y), dFdy(x, y)};
               // gradF = {dFdx(x + dx, y + dy), dFdy(x + dx, y + dy)};

               BM.updateGoodBroyden(gradF, gradF_, i < 1 ? 0.01 : 1.);
               // BM.update(gradF, gradF_, 1.);

               if (Norm(BM.dX) < 1e-10)
                  break;
            }
            ofs_BM.close();
         }

         //! Broyden Method with BFGS update
         {
            std::ofstream ofs_BM("./output_Himmelblau/Himmelblau_BFGS" + std::to_string((int)x) + "_" + std::to_string((int)y) + ".csv");
            BroydenMethod<std::array<double, 2>> BM(X, X + dX);
            // ofs_BM << "#x y F(x, y)" << std::endl;
            for (auto i = 0; i < 50; i++) {
               auto x = BM.X[0];  // optimization variable
               auto y = BM.X[1];  // optimization variable
               auto dx = BM.dX[0];
               auto dy = BM.dX[1];
               std::cout << w << i << w << x << w << y << w << F(x, y) << w << Norm(BM.dX) << std::endl;
               ofs_BM << x << "," << y << "," << F(x, y) << std::endl;
               gradF_ = {dFdx(x - dx, y - dy), dFdy(x - dx, y - dy)};
               gradF = {dFdx(x, y), dFdy(x, y)};

               // gradF_ = {dFdx(x, y), dFdy(x, y)};
               // gradF = {dFdx(x + dx, y + dy), dFdy(x + dx, y + dy)};

               BM.updateBFGS(gradF, gradF_, i < 1 ? 0.01 : 1.);
               // BM.update(gradF, gradF_, 1.);

               if (Norm(BM.dX) < 1e-10)
                  break;
            }
            ofs_BM.close();
         }

         //! updateInverseBFGS
         {
            std::ofstream ofs_BM("./output_Himmelblau/Himmelblau_InvBFGS" + std::to_string((int)x) + "_" + std::to_string((int)y) + ".csv");
            BroydenMethod<std::array<double, 2>> BM(X, X + dX);
            // ofs_BM << "#x y F(x, y)" << std::endl;
            for (auto i = 0; i < 50; i++) {
               auto x = BM.X[0];  // optimization variable
               auto y = BM.X[1];  // optimization variable
               auto dx = BM.dX[0];
               auto dy = BM.dX[1];
               std::cout << w << i << w << x << w << y << w << F(x, y) << w << Norm(BM.dX) << std::endl;
               ofs_BM << x << "," << y << "," << F(x, y) << std::endl;
               gradF_ = {dFdx(x - dx, y - dy), dFdy(x - dx, y - dy)};
               gradF = {dFdx(x, y), dFdy(x, y)};

               // gradF_ = {dFdx(x, y), dFdy(x, y)};
               // gradF = {dFdx(x + dx, y + dy), dFdy(x + dx, y + dy)};

               BM.updateInverseBFGS(gradF, gradF_, i < 1 ? 0.01 : 1.);
               // BM.update(gradF, gradF_, 1.);

               if (Norm(BM.dX) < 1e-10)
                  break;
            }
            ofs_BM.close();
         }
      }
}