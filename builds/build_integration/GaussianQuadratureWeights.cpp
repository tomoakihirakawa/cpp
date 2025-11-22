#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "basic.hpp"
#include "rootFinding.hpp"

/*DOC_EXTRACT 0_2_0_integration_TrapezoidalRule

ガウス求積法は，ルジャンドル多項式の根を利用して関数の積分を近似する方法である．
`GaussianQuadrature`では，以下のようにして，根の初期値を与え，ニュートンラフソン法を用いて正しい根へと収束させる．

```cpp
const bool is_even = (N % 2 == 0);
const int M = (is_even ? N / 2 : (N + 1) / 2);
int k = 0;
for (int i = 1; i <= M; ++i)
    x[k++] = -std::sin(M_PI * (N + 1 - 2 * i) / (2. * N + 1.));  //! first guess
for (int i = (is_even ? 0 : 1); i < N - M + (is_even ? 0 : 1); ++i)
    x[k++] = -x[M - i - 1];
```

`GaussianQuadrature`のメンバ変数`x`には，ルジャンドル多項式の根が，`w`には，各点でのガウス求積の重みが格納される．

<img src="results_roots_weights.png" width="700">

基本的には台形則による数値積分よりも，ガウス求積法の方が精度が高い．

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=GaussianQuadrature.cpp
make
./GaussianQuadrature
```

##　特異積分


*/

struct GaussianQuadrature {
   std::vector<double> x;
   std::vector<double> first_guess;
   std::vector<double> w;
   const int N;
   GaussianQuadrature(const int N) : N(N), x(N), w(N) { initialize(); };
   GaussianQuadrature(const int N, const double a, const double b) : N(N), x(N), w(N) {
      initialize();
      for (int i = 0; i < N; i++) {
         x[i] = 0.5 * ((b - a) * x[i] + b + a);
         w[i] *= 0.5 * (b - a);
      }
   };
   GaussianQuadrature(const int N, double a, double b, const double beta, const double etaSing) : N(N), x(N), w(N) {
      initialize();
      a = EtaInv(a, beta, etaSing);
      b = EtaInv(b, beta, etaSing);
      double u;
      for (int i = 0; i < N; i++) {
         u = 0.5 * ((b - a) * x[i] + b + a);
         x[i] = Eta(u, beta, etaSing);
         w[i] *= 0.5 * (b - a) * D_Eta(u, beta);
      }
   };

   void initialize() {
      if (N == 1)
         x[0] = 0.;
      else {
         const bool is_even = (N % 2 == 0);
         const int M = (is_even ? N / 2 : (N + 1) / 2);
         int k = 0;
         for (int i = 1; i <= M; ++i)
            x[k++] = -std::sin(M_PI * (N + 1 - 2 * i) / (2. * N + 1.));  //! first guess
         for (int i = (is_even ? 0 : 1); i < N - M + (is_even ? 0 : 1); ++i)
            x[k++] = -x[M - i - 1];
      }
      this->first_guess = x;
      auto i = 0;
      for (auto& y : x) {
         NewtonRaphson<double> NR(y);
         for (auto i = 0; i < 10; i++) {
            NR.update(std::legendre(N, NR.X), D_Legendre(NR.X));
            if (std::abs(NR.dX) < 1e-15)
               break;
         }
         y = NR.X;
         w[i++] = 2 / ((1 - y * y) * std::pow(D_Legendre(y), 2));
      }
   }

   double D_Legendre(const double x, const int derivative = 1) {
      if (N == 0)
         return 0;
      else if (N == 1)
         return 1;
      else
         return std::assoc_legendre(N, derivative, x) * std::pow(1. - x * x, -0.5 * derivative);
   };

   double integrate(const std::function<double(double)>& f) {
      double integral = 0.;
      for (int i = 0; i < N; i++)
         integral = std::fma(w[i], f(x[i]), integral);
      return integral;
   }

   double Eta(const double u, const double beta, const double etaSing) {
      return etaSing + std::copysign(1., u) * std::pow(std::abs(u), beta);
   };

   double EtaInv(const double eta, const double beta, const double etaSing) {
      return std::copysign(1., eta - etaSing) * std::pow(std::abs(eta - etaSing), 1. / beta);
   };

   double D_Eta(const double u, const double beta) {
      if (beta == 1.)
         return 1.;
      else
         return beta * std::pow(std::abs(u), beta - 1.);
      // return std::copysign(1., u) * std::pow(u * u, beta / 2.) * (beta + u) / u;
   };
};

int main() {

   double N = 100;

   std::ofstream ofs_roots("LegendrePolynomials_roots.dat");
   std::ofstream ofs_weight("LegendrePolynomials_weight.dat");

   for (auto order = 1; order < 10; order++) {
      GaussianQuadrature GQW(order);
      for (auto x : GQW.x)
         ofs_roots << std::setprecision(15) << x << " ";
      for (auto w : GQW.w)
         ofs_weight << std::setprecision(15) << w << " ";
      ofs_roots << std::endl;
      ofs_weight << std::endl;
   }
   ofs_roots.close();
   ofs_weight.close();

   {
      const double c = 0.01;
      auto f = [&](double x) { return 1. / (x + c); };
      auto analitical_integral = [&](double a, double b) { return std::log(b + c) - std::log(a + c); };
      std::ofstream ofs("GW_1_x.dat");
      std::ofstream ofs_SING("GW_1_x_SING.dat");
      double a = 0.0001, b = 2.;
      for (int N = 1; N < 50; N++) {
         auto integral = GaussianQuadrature(N, a, b).integrate(f);
         //  auto integral_SING = GaussianQuadrature(N).integrate(f, a, b, 2., 0.);
         auto integral_SING = GaussianQuadrature(N, a, b, 2., 0.).integrate(f);
         std::cout << N << " " << integral << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs << N << " " << integral << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs_SING << N << " " << integral_SING << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral_SING - analitical_integral(a, b))) << std::endl;
      }
      ofs.close();
      ofs_SING.close();
   }

   {
      auto f = [](double x) { return std::sqrt(1. - x * x); };
      auto analitical_integral = [](double a, double b) { return M_PI / 4.; };
      std::ofstream ofs("GW_sqrt1xx.dat");
      std::ofstream ofs_SING("GW_sqrt1xx_SING.dat");
      double a = 0., b = 1.;
      for (int N = 1; N < 50; N++) {
         auto integral = GaussianQuadrature(N, a, b).integrate(f);
         auto GWSING = GaussianQuadrature(N, a, b, 2., 1.);
         double integral_SING = 0;
         for (int i = 0; i < N; i++) {
            integral_SING += GWSING.w[i] * f(GWSING.x[i]);
         }
         std::cout << N << " " << integral << " " << integral_SING << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs << N << " " << integral << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs_SING << N << " " << integral_SING << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral_SING - analitical_integral(a, b))) << std::endl;
      }
      ofs.close();
      ofs_SING.close();
   }

   {
      auto f = [](double x) { return 1. / std::sqrt(1. - x * x); };
      auto analitical_integral = [](double a, double b) { return M_PI / 2.; };
      std::ofstream ofs("GW_1sqrt1xx.dat");
      std::ofstream ofs_SING("GW_1sqrt1xx_SING.dat");
      double a = 0., b = 1.;
      for (int N = 1; N < 50; N++) {
         auto integral = GaussianQuadrature(N, a, b).integrate(f);
         auto GWSING = GaussianQuadrature(N, a, b, 2., 1.);
         double integral_SING = 0;
         for (int i = 0; i < N; i++) {
            integral_SING += GWSING.w[i] * f(GWSING.x[i]);
         }
         std::cout << N << " " << integral << " " << integral_SING << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs << N << " " << integral << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs_SING << N << " " << integral_SING << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral_SING - analitical_integral(a, b))) << std::endl;
      }
      ofs.close();
      ofs_SING.close();
   }

   //! BEMで現れる積分のチェック
   {
      std::ofstream ofs("BEM_integration_accuracy_check.dat");
      auto V = Subdivide<50>(-3, 3);
      GaussianQuadrature GQ6(6, -1, 1.);
      const int MAX_N = 300;
      GaussianQuadrature GQ_MAX(MAX_N, -1, 1.);
      std::array<double, 3> X0 = {std::cos(2 * M_PI / 3.), std::sin(2 * M_PI / 3.), 0};
      std::array<double, 3> X1 = {std::cos(4 * M_PI / 3.), std::sin(4 * M_PI / 3.), 0};
      std::array<double, 3> X2 = {std::cos(6 * M_PI / 3.), std::sin(6 * M_PI / 3.), 0}, A;
      for (const auto& x : V) {
         for (const auto& y : V) {
            for (const auto& z : V) {
               A = {x, y, z};
               std::array<std::array<double, 3>, 3> X012 = {X0, X1, X2};

               auto Integral = [&](const int N = 6) {
                  double integral = 0;
                  if (N == 6) {
                     for (int i = 0; i < N; i++)
                        for (int j = 0; j < N; j++)
                           integral += GQ6.w[i] * GQ6.w[j] * ((1. - GQ6.x[i]) / Norm(Dot(ModTriShape<3>(GQ6.x[i], GQ6.x[j]), X012) - A));
                  } else {
                     std::array<std::array<double, 3>, 3> X012_ = {X0, X1, X2};
                     if (Norm(A - X0) >= Norm(A - X1) && Norm(A - X2) >= Norm(A - X1))
                        X012_ = {X1, X2, X0};
                     else if (Norm(A - X0) >= Norm(A - X2) && Norm(A - X1) >= Norm(A - X2))
                        X012_ = {X2, X0, X1};
                     for (int i = 0; i < N; i++)
                        for (int j = 0; j < N; j++)
                           integral += GQ_MAX.w[i] * GQ_MAX.w[j] * ((1. - GQ_MAX.x[i]) / Norm(Dot(ModTriShape<3>(GQ_MAX.x[i], GQ_MAX.x[j]), X012_) - A));
                  }
                  return integral;
               };

               auto ans_apprx = Integral(MAX_N);
               auto error = std::log10(std::abs(Integral() - ans_apprx));
               // std::cout << "N = " << N << " integral = " << std::setprecision(15) << error << std::endl;
               ofs << x << " " << y << " " << z << " " << std::setprecision(15) << error << std::endl;
            }
         }
      }
   }
   //! compare
   // {
   //    int N = 6;
   //    GaussianQuadrature GQ1(N, -1, 1.);
   //    auto V = GaussianQuadratureWeights(N, -1, 1.);
   //    for (auto i = 0; i < N; i++) {
   //       std::cout << GQ1.x[i] << " " << V[i][0] << " " << GQ1.w[i] << " " << V[i][1] << " diff = " << std::log10(std::abs(GQ1.w[i] - V[i][1])) << " diff = " << std::log10(std::abs(GQ1.x[i] - V[i][0])) << std::endl;
   //    }
   // }
   return 0;
};
