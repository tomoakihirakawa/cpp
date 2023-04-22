#include <cmath>
#include <iostream>

// #define polynomial
#define interpolation

#if defined(polynomial)
int main() {

   double N = 100;
   for (int i = -N; i < N; i++) {
      double x = i / N;
      std::cout << x << " ";
      for (auto order = 0; order < 10; order++)
         std::cout << std::legendre(order, x) << " ";
      std::cout << std::endl;
   }

   /*
   1. $ make;./LegendrePolynomials > LegendrePolynomials.dat
   2. $ gnuplot
   3. type:
   file = 'LegendrePolynomials.dat'
   plot file using 1:2 title 'order 0', \
        file using 1:3 title 'order 1', \
        file using 1:4 title 'order 2', \
        file using 1:5 title 'order 3', \
        file using 1:6 title 'order 4', \
        file using 1:7 title 'order 5'
    */
   return 0;
};

#elif defined(interpolation)

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

   double N = 100;

   auto f = [](double x) { return 1. / (1. + 25. * x * x); };

   auto a = [&](int n) {
      return (2. * n + 1.) / 2. * TrapezoidalRule([&](auto x) { return f(x) * std::legendre(n, x); }, -1., 1., 50);
   };

   auto f_approx = [&](auto x, int max_order) {
      double approx = 0;
      for (int n = 0; n < max_order; n++)
         approx += a(n) * std::legendre(n, x);
      return approx;
   };

   for (int i = -N; i < N; i++) {
      double x = 1.1 * i / N;
      std::cout << x << " " << f(x) << " ";
      for (auto order : {10, 20, 30})
         std::cout << f_approx(x, order) << " ";
      std::cout << std::endl;
   }

   /*
   1. $ make;./LegendrePolynomials > LegendrePolynomials.dat
   2. $ gnuplot
   3. type:
   file = 'LegendrePolynomials.dat'
   set xrange [-1.1:1.1]
   set yrange [-0.1:1.1]
   plot file using 1:2 title '1/(1+25*x*x)'  with lines
   replot for [i=3:5] file using 1:i title 'order '.(i-2) with lines
   */
};

#endif