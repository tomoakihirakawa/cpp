#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
#include "basic.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "interpolations.hpp"

/*DOC_EXTRACT interpolation

## ラグランジュ補間

与えられたデータ点を通る多項式を求める方法の一つにラグランジュ補間がある．

```math
f(x) = \sum_{i=0}^n\dfrac{\prod_{j=0,j\neq i}^n{(x - x_j)}}{\prod_{j=0,j\neq i,j\neq k}^n{(x_i - x_j)}}y_i
```

微分は，

```math
f(x) = \sum_{i=0}^n\dfrac{\sum_{k=0}^{n}\prod_{j=0,j\neq i}^n{(x - x_j)}}{\prod_{j=0,j\neq i}^n{(x_i - x_j)}}y_i
```

![](sample_lag.png)

*/

int main() {

   auto curve = [](const double x) {
      return sin(x);
   };

   auto D_curve = [](const double x) {
      return cos(x);
   };

   std::vector<double> abscissas;
   std::vector<double> values;
   const int N = 20;

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
      auto filename = "lag_exact_interpolation_derivative.dat";
      std::ofstream file(filename);
      if (!file.is_open()) {
         std::cerr << "Failed to open the output file." << std::endl;
         return 0;
      }
      file << "# x y" << std::endl;
      // for (auto t : Subdivide<100>(-1., 11.))
      //    file << t << " " << curve(t) << " " << D_curve(t) << " " << intL(t) << " " << intL.D(t) << std::endl;
      for (auto t : Subdivide<100>(-1., 11.))
         file << t << " "
              << curve(t) << " "
              << D_curve(t) << " "
              << Dot(intL.N(t), intL.values) << " "
              << Dot(intL.DN(t), intL.values) << std::endl;
   }

   return 0;
}