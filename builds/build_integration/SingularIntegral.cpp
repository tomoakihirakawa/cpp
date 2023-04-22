#include <cmath>
#include <iostream>

template <typename Func, typename arg>
auto TrapezoidalRule(const Func& f, const arg& a, const arg& b, const int N) {
   auto integral = (f(a) + f(b)) / 2.;
   auto dx = (b - a) / N;
   for (int i = 1; i < N; ++i)
      integral += f(a + i * dx);
   integral *= dx;
   return integral;
};

template <typename Func, typename arg>
auto TrapezoidalRule(const Func& f, arg a, arg b, const int N, double A) {
   auto F = [&](const auto& xi) {
      double _11A = 1. / (1. - A);
      double dxdxi = std::pow(xi, A * _11A) * _11A;
      return f(std::pow(xi, _11A)) * dxdxi;
   };
   a = std::pow(a, 1. - A);
   b = std::pow(b, 1. - A);
   return TrapezoidalRule(F, a, b, N);
};

int main() {
   //    auto f = [](auto x) { return 1. / sqrt(x); };
   //    auto analitical_integral = [](auto a, auto b) { return 2 * (std::sqrt(b) - std::sqrt(a)); };

   //    double a = 1E-15, b = 10.;
   //    for (int N = 2; N < 100; N++) {
   //       auto integral = TrapezoidalRule(f, a, b, N, 0.5);
   //       std::cout << N << " " << integral << " " << analitical_integral(a, b) << " " << log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
   //    }
   /* -------------------------------------------------------------------------- */
   auto f = [](auto x) { return 1. / x; };
   auto analitical_integral = [](auto a, auto b) { return std::log(b) - std::log(a); };

   double a = 1E-1, b = 2.;
   for (int N = 2; N < 100; N++) {
      auto integral = TrapezoidalRule(f, a, b, N, 0.99);
      std::cout << N << " " << integral << " " << analitical_integral(a, b) << " " << log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
   }
   /* -------------------------------------------------------------------------- */
   //    auto f = [](auto x) { return 1. / (x * x); };
   //    auto analitical_integral = [](auto a, auto b) { return -1 / b + 1 / a; };

   //    double a = 1E-2, b = 1.;
   //    for (int N = 2; N < 100; N++) {
   //       auto integral = TrapezoidalRule(f, a, b, N, 2.);
   //       std::cout << N << " " << integral << " " << analitical_integral(a, b) << " " << log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
   //    }
   /*
   この方法は，ゼロ付近でA次の特異性を持つ関数の積分の特異性を除去することができる
   */
   return 0;
};
