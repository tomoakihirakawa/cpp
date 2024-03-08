#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "rootFinding.hpp"

/*DOC_EXTRACT 0_1_0_integration_LegendrePolynomials

## ルジャンドル多項式，ルジャンドル補間，ガウス・ルジャンドル積分

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=LegendrePolynomials.cpp
make
./LegendrePolynomials
```

Mathematicaの結果とcppの結果を比較:

<img src="results_legendre_polynomial.png" width="700">

<img src="results_D_legendre_polynomial.png" width="700">

`結果のチェック.nb`で結果を確認．

*/

double D_Legendre(int n, double x, const int der_order = 1) {
   if (n == 0) return 0;
   if (n == 1) return 1;
   return std::assoc_legendre(n, der_order, x) * std::pow(1. - x * x, -0.5 * der_order);
};

std::vector<double> first_guess(const int n) {
   const bool is_even = n % 2 == 0;
   std::vector<double> vec;
   for (int i = 1; i <= (is_even ? n / 2 : (n + 1) / 2); i++)
      vec.push_back(std::sin(M_PI * (n + 1 - 2 * i) / (2. * n + 1.)));

   if (n <= 1)
      return vec;

   auto VEC = vec;
   for (auto it = vec.rbegin() + (is_even ? 0 : 1); it != vec.rend(); ++it)
      VEC.push_back(-*it);
   return VEC;
};

#define polynomial
// #define interpolation

#if defined(polynomial)
int main() {

   double N = 100;
   std::ofstream ofs("LegendrePolynomials.dat");
   std::ofstream ofs_D("LegendrePolynomials_D.dat");
   for (int i = -N; i < N; i++) {
      double x = i / N;
      ofs << x << " ";
      ofs_D << x << " ";
      for (auto order = 1; order < 10; order++) {
         ofs << std::legendre(order, x) << " ";
         ofs_D << D_Legendre(order, x) << " ";
      }

      ofs << std::endl;
      ofs_D << std::endl;
   }

   std::ofstream ofs_guess("first_guess.dat");
   for (auto order = 1; order < 10; order++) {
      for (auto v : first_guess(order))
         ofs_guess << v << " ";
      ofs_guess << std::endl;
   }
   ofs_guess.close();

   std::ofstream ofs_roots("LegendrePolynomials_roots.dat");
   std::ofstream ofs_weight("LegendrePolynomials_weight.dat");

   auto weight = [](const double x, const int N) {
      return 2 / ((1 - x * x) * std::pow(D_Legendre(N, x, 1), 2));
   };

   for (auto order = 1; order < 10; order++) {
      for (auto v : first_guess(order)) {
         NewtonRaphson<double> NR(v);
         for (auto i = 0; i < 10; i++) {
            NR.update(std::legendre(order, NR.X), D_Legendre(order, NR.X, 1));
            std::cout << NR.X << " " << std::legendre(order, NR.X) << std::endl;
         }
         ofs_roots << std::setprecision(15) << NR.X << " ";
         ofs_weight << std::setprecision(15) << weight(NR.X, order) << " ";
      }
      ofs_roots << std::endl;
      ofs_weight << std::endl;
   }
   ofs_roots.close();
   ofs_weight.close();
  
   return 0;
};