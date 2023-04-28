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

int main() {

   auto f = [](auto x) { return std::sin(x); };
   auto analitical_integral = [](auto a, auto b) { return -(std::cos(b) - std::cos(a)); };

   double a = 0, b = M_PI;
   for (int N = 1; N < 100; N++) {
      auto integral = TrapezoidalRule(f, a, b, N);
      std::cout << N << " " << integral << " " << analitical_integral(a, b) << " " << log10(std::abs(integral - analitical_integral(a, b))) << std::endl;
   }
   return 0;
};
