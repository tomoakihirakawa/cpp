/*DOC_EXTRACT interpolation

## 三角形補間

\insert{interpolation:ModTriShape}

![](sample_tri.png)

*/

// #include "Network.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"

// cmake -DCMAKE_BUID_TYPE=Release ../ -DSOURCE_FILE=TriShape.cpp
// gnuplot
// splot 'triangular_surface.dat' w l

int main() {

   auto filename = "triangular_surface.dat";
   std::ofstream file(filename);
   if (!file.is_open()) {
      std::cerr << "Failed to open the output file." << std::endl;
      return 0;
   }
   file << "# x y z" << std::endl;
   double resolution = 10;
   double step = 1 / resolution;
   std::array<double, 3> a{{0, 0, 0}};
   std::array<double, 3> b{{1, 0, 0}};
   std::array<double, 3> c{{0, 2, -1}};

   //    {
   //       std::array<std::array<double, 3>, 3> M{{a, b, c}};
   //       for (size_t i = 0; i <= resolution; ++i) {
   //          double t0 = step * i;
   //          for (size_t j = 0; j <= resolution; ++j) {
   //             double t1 = step * j;
   //             auto [x, y, z] = Dot(ModTriShape<double, 3>(t0, t1), M);
   //             file << x << " " << y << " " << z << std::endl;
   //          }
   //       }
   //    }
   std::array<double, 3> d = {0.5, 0.0, 0.5};
   std::array<double, 3> e = {0.5, 0.5, 0.5};
   std::array<double, 3> f = {0.0, 0.5, 0.5};

   std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
   for (const auto& t0 : Subdivide<50>(0., 1.)) {
      for (const auto& t1 : Subdivide<50>(0., 1.)) {
         auto [x, y, z] = Dot(ModTriShape<double, 6>(t0, t1, std::array<bool, 3>{1, 1, 1}), M);
         file << x << " " << y << " " << z << std::endl;
      }
      file << std::endl;
   }
   /* -------------------------------------------------------------------------- */
   {
      auto filename = "points.dat";
      std::ofstream file(filename);
      if (!file.is_open()) {
         std::cerr << "Failed to open the output file." << std::endl;
         return 0;
      }
      file << "# x y z index" << std::endl;
      std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
      int index = 0;
      for (auto& [x, y, z] : M)
         file << x << " " << y << " " << z << " " << index++ << std::endl;
   }
}
