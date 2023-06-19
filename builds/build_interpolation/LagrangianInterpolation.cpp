#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
#include "basic_arithmetic_array_operations.hpp"

/*DOC_EXTRACT interpolation

## ラグランジュ補間

与えられたデータ点を通る多項式を求める方法の一つにラグランジュ補間がある．

```math
f(x) = \sum_{i=0}^n\dfrac{\prod_{j=0,j\neq i}^n{(x - x_j)}}{\prod_{j=0,j\neq i}^n{(x_i - x_j)}}y_i
```

補間の式の分母は，補間したい横軸$`x`$に依存せず，与えられた横軸データのみに依存するので，
データが与えられたら1度だけ計算しておけばよい．

![](sample_lag.png)

*/

template <typename T>
struct InterpolationLagrange {
   std::vector<double> abscissas;
   std::vector<T> values;
   std::vector<double> denominotor;
   InterpolationLagrange(const std::vector<double> abscissas, const std::vector<T> values)
       : abscissas(abscissas), values(values) {
      this->set();
   };

   void set() {
      denominotor.resize(abscissas.size(), 1.);
      for (auto i = 0; i < abscissas.size(); ++i)
         for (auto j = 0; j < abscissas.size(); ++j)
            if (i != j)
               denominotor[i] *= (abscissas[i] - abscissas[j]);
   };

   T operator()(const double x) {
      T ret;
      ret *= 0.;
      for (auto i = 0; i < abscissas.size(); ++i) {
         T ret_i = this->values[i] / this->denominotor[i];
         for (auto j = 0; j < abscissas.size(); ++j)
            if (i != j)
               ret_i *= (x - abscissas[j]);
         ret += ret_i;
      }
      return ret;
   };
};

int main() {

   auto curve = [](const double x) {
      return sin(x) + x / 5;
   };

   std::vector<double> abscissas;
   std::vector<double> values;
   const int N = 7;

   for (auto& T : Subdivide<N>(0., 10.)) {
      auto t = T + (rand() / (RAND_MAX + 1.0) - 0.5);
      abscissas.push_back(t);
      values.push_back(curve(t));
   }

   InterpolationLagrange<double> intL(abscissas, values);

   {
      auto filename = "lag_data.dat";
      std::ofstream file(filename);
      if (!file.is_open()) {
         std::cerr << "Failed to open the output file." << std::endl;
         return 0;
      }
      file << "# x y index" << std::endl;
      for (auto i = 0; i < abscissas.size(); ++i)
         file << abscissas[i] << " " << values[i] << " " << i << std::endl;
   }

   {
      auto filename = "lag_interpolation.dat";
      std::ofstream file(filename);
      if (!file.is_open()) {
         std::cerr << "Failed to open the output file." << std::endl;
         return 0;
      }
      file << "# x y" << std::endl;
      for (auto t : Subdivide<100>(0., 10.))
         file << t << " " << intL(t) << std::endl;
   }

   {
      auto filename = "lag_exact.dat";
      std::ofstream file(filename);
      if (!file.is_open()) {
         std::cerr << "Failed to open the output file." << std::endl;
         return 0;
      }
      file << "# x y" << std::endl;
      for (auto t : Subdivide<100>(0., 10.))
         file << t << " " << curve(t) << std::endl;
   }

   return 0;
}