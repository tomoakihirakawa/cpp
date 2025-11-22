#include <cmath>
#include <fstream>
#include <iostream>
/*DOC_EXTRACT 0_0_0_integration_TrapezoidalRule

# 数値積分

## 台形則

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=TrapezoidalRule.cpp
make
./TrapezoidalRule
```

台形則は，関数の積分を台形の面積の和で近似する方法である．

```math
\int_a^b f(x) dx \approx \left(\frac{f(a)+f(b)}{2} + \sum_{i=1}^{N-1} f(a+i\Delta x)\right)\Delta x, \quad \Delta x = \frac{b-a}{N}
```

### 例

* $`\int_0^\pi \sin(x) dx = 2`$

<img src="TrapezoidalRule_sin.png" width="700">

特異性のある積分の例：

* $`\int_\varepsilon^\pi \frac{1}{x} dx = \log{\pi} - \log{\varepsilon}`$

<image src="results_1_x.png" width="700">

* $`\int_0^1 \sqrt{1-x^2} dx = \frac{\pi}{4}`$

<image src="results_1sqrt1xx.png" width="700">

* $`\int_0^1 \frac{1}{\sqrt{1-x^2}} dx = \frac{\pi}{2}`$

<image src="results_sqrt1xx.png" width="700">


*/

template <typename Func, typename arg>
auto TrapezoidalRule(const Func& f, const arg& a, const arg& b, const int N) {
   auto integral = (f(a) + f(b)) / 2.;
   auto dx = (b - a) / N;
   for (int i = 1; i < N; ++i)
      integral += f(a + i * dx);
   integral *= dx;
   return integral;
};

/*
Eta[u_, beta_, etaSing_] := etaSing + Sign[u] Abs[u]^beta;
EtaInv[eta_, beta_, etaSing_] :=
  Sign[eta - etaSing]*Abs[eta - etaSing]^(1/beta);
*/

double Eta(const double u, const double beta, const double etaSing) {
   return etaSing + std::copysign(1., u) * std::pow(std::abs(u), beta);
};

double EtaInv(const double eta, const double beta, const double etaSing) {
   return std::copysign(1., eta - etaSing) * std::pow(std::abs(eta - etaSing), 1. / beta);
};

double D_Eta(const double u, const double beta) {
   return beta * std::pow(std::abs(u), beta - 1.);
   // return std::copysign(1., u) * std::pow(u * u, beta / 2.) * (beta + u) / u;
};

template <typename Func, typename arg>
auto TrapezoidalRule(const Func& f, arg a, arg b, const int N, const double beta, const double etaSing) {
   // auto F = [&](const auto& xi) {
   //    double _11A = 1. / (1. - A);
   //    double dxdxi = std::pow(xi, A * _11A) * _11A;
   //    return f(std::pow(xi, _11A)) * dxdxi;
   // };
   // a = std::pow(a, 1. - A);
   // b = std::pow(b, 1. - A);
   // return TrapezoidalRule(F, a, b, N);
   const double ua = EtaInv(a, beta, etaSing);
   const double ub = EtaInv(b, beta, etaSing);
   auto F = [&](const auto& u) { return f(Eta(u, beta, etaSing)) * D_Eta(u, beta); };
   return TrapezoidalRule(F, ua, ub, N);
};

int main() {

   {
      auto f = [](auto x) { return std::sin(x); };
      auto analitical_integral = [](auto a, auto b) { return -(std::cos(b) - std::cos(a)); };
      std::ofstream ofs("TrapezoidalRule_sin.dat");
      double a = 0, b = M_PI;
      for (int N = 1; N < 100; N++) {
         auto integral = TrapezoidalRule(f, a, b, N);
         std::cout << N << " " << integral << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs << N << " " << integral << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
      }
      ofs.close();
   }

   {
      const double c = 0.01;
      auto f = [&](double x) { return 1. / (x + c); };
      auto analitical_integral = [&](double a, double b) { return std::log(b + c) - std::log(a + c); };
      std::ofstream ofs("TrapezoidalRule_1_x.dat");
      std::ofstream ofs_SING("TrapezoidalRule_1_x_SING.dat");
      double a = 0.0001, b = 2.;
      for (int N = 1; N < 50; N++) {
         auto integral = TrapezoidalRule(f, a, b, N);
         auto integral_SING = TrapezoidalRule(f, a, b, N, 2., 0.);
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
      std::ofstream ofs("TrapezoidalRule_sqrt1xx.dat");
      std::ofstream ofs_SING("TrapezoidalRule_sqrt1xx_SING.dat");
      double a = 0., b = 1.;
      for (int N = 1; N < 50; N++) {
         auto integral = TrapezoidalRule(f, a, b, N);
         auto integral_SING = TrapezoidalRule(f, a, b, N, 2., 1.);
         std::cout << N << " " << integral << " " << integral_SING << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs << N << " " << integral << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
         ofs_SING << N << " " << integral_SING << " " << analitical_integral(a, b) << " " << std::log10(std::abs(integral_SING - analitical_integral(a, b))) << std::endl;
      }
      ofs.close();
      ofs_SING.close();
   }

   return 0;
};
