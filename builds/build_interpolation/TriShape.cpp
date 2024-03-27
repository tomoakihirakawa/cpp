/*DOC_EXTRACT 1_1_0_interpolation

\insert{interpolation:TriShape}

| 線形補間 | 2次補間 |
| --- | --- |
| <img src="triangle_shape_function_linear.png" width="400"> | <img src="triangle_shape_function_quadratic.png" width="300"> |

\insert{interpolation:ModTriShape}

### 例：補間によって，頂点座標から平面を作成する

<img src="sample_tri.png" width="400">

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=TriShape.cpp
make
./TriShape
```

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

   auto filename_2nd = "triangular_surface_2nd.dat";
   auto filename_2nd_xoo = "triangular_surface_2nd_xoo.dat";
   auto filename_2nd_oxo = "triangular_surface_2nd_oxo.dat";
   auto filename_2nd_oox = "triangular_surface_2nd_oox.dat";
   auto filename_2nd_through_linear = "triangular_surface_2nd_through_linear.dat";
   auto filename_linear = "triangular_surface_linear.dat";
   std::ofstream file_2nd(filename_2nd);
   std::ofstream file_2nd_xoo(filename_2nd_xoo);
   std::ofstream file_2nd_oxo(filename_2nd_oxo);
   std::ofstream file_2nd_oox(filename_2nd_oox);
   std::ofstream file_2nd_through_linear(filename_2nd_through_linear);
   std::ofstream file_linear(filename_linear);

   if (!file_2nd.is_open() || !file_linear.is_open() || !file_2nd_xoo.is_open() || !file_2nd_oxo.is_open() || !file_2nd_oox.is_open()) {
      std::cerr << "Failed to open the output file." << std::endl;
      return 0;
   }

   /* -------------------------------------------------------------------------- */
   const double resolution = 10;
   const double step = 1 / resolution;
   const std::array<double, 3> a{{0, 0, 0}};
   const std::array<double, 3> b{{1, 0, 0}};
   const std::array<double, 3> c{{0, 2, -1}};
   const std::array<double, 3> d = {0.5, 0.0, 0.5};
   const std::array<double, 3> e = {0.5, 0.5, 0.5};
   const std::array<double, 3> f = {0.0, 0.5, 0.5};
   /* -------------------------------------------------------------------------- */
   {
      std::array<std::array<double, 3>, 3> M{{a, b, c}};
      for (size_t i = 0; i <= resolution; ++i) {
         double t0 = step * i;
         for (size_t j = 0; j <= resolution; ++j) {
            double t1 = step * j;
            auto [x, y, z] = Dot(ModTriShape<3>(t0, t1), M);
            file_linear << x << " " << y << " " << z << std::endl;
         }
      }
      file_linear.close();
      std::cout << "Output file: " << filename_linear << std::endl;
   }
   /* -------------------------------------------------------------------------- */
   {
      std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
      for (const auto& t0 : Subdivide<50>(0., 1.)) {
         for (const auto& t1 : Subdivide<50>(0., 1.)) {
            auto [x, y, z] = Dot(ModTriShape<6>(t0, t1, std::array<bool, 3>{1, 1, 1}), M);
            file_2nd << x << " " << y << " " << z << std::endl;
         }
      }
      file_2nd.close();
      std::cout << "Output file: " << filename_2nd << std::endl;
   }
   /* -------------------------------------------------------------------------- */
   {
      std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
      for (const auto& t0 : Subdivide<50>(0., 1.)) {
         for (const auto& t1 : Subdivide<50>(0., 1.)) {
            auto [x, y, z] = Dot(ModTriShape<6>(t0, t1, std::array<bool, 3>{0, 1, 1}), M);
            file_2nd_xoo << x << " " << y << " " << z << std::endl;
         }
      }
      file_2nd_xoo.close();
      std::cout << "Output file: " << filename_2nd_xoo << std::endl;
   }

   /* -------------------------------------------------------------------------- */
   {
      std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
      for (const auto& t0 : Subdivide<50>(0., 1.)) {
         for (const auto& t1 : Subdivide<50>(0., 1.)) {
            auto [x, y, z] = Dot(ModTriShape<6>(t0, t1, std::array<bool, 3>{1, 0, 1}), M);
            file_2nd_oxo << x << " " << y << " " << z << std::endl;
         }
      }
      file_2nd_oxo.close();
      std::cout << "Output file: " << filename_2nd_oxo << std::endl;
   }
   /* -------------------------------------------------------------------------- */
   {
      std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
      for (const auto& t0 : Subdivide<50>(0., 1.)) {
         for (const auto& t1 : Subdivide<50>(0., 1.)) {
            auto [x, y, z] = Dot(ModTriShape<6>(t0, t1, std::array<bool, 3>{1, 1, 0}), M);
            file_2nd_oox << x << " " << y << " " << z << std::endl;
         }
      }
      file_2nd_oox.close();
      std::cout << "Output file: " << filename_2nd_oox << std::endl;
   }
   /* -------------------------------------------------------------------------- */
   {
      std::array<double, 2> t0t1_No4 = {0., 0.5};
      std::array<double, 2> t0t1_No5 = {0.5, 0.};
      std::array<double, 2> t0t1_No3 = {0.5, 0.5};
      std::array<std::array<double, 2>, 3> M_sub{{t0t1_No4, t0t1_No5, t0t1_No3}};
      std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
      for (const auto& t0 : Subdivide<5>(0., 1.)) {
         for (const auto& t1 : Subdivide<5>(0., 1.)) {
            // auto [x, y, z] = Dot(TriShape<6>(Dot(ModTriShape<3>(t0, t1), M_sub), std::array<bool, 3>{1, 1, 1}), M);
            auto [T0, T1] = Dot(ModTriShape<3>(t0, t1), M_sub);
            auto [x, y, z] = Dot(D_TriShape_Quadratic<0, 0>(T0, T1), M);
            auto dxi0 = Dot(D_TriShape_Quadratic<1, 0>(T0, T1), M);
            auto dxi1 = Dot(D_TriShape_Quadratic<0, 1>(T0, T1), M);
            auto dxi01 = Dot(D_TriShape_Quadratic<1, 1>(T0, T1), M);
            file_2nd_through_linear << std::fixed << x << " " << y << " " << z << " " << dxi0 << " " << dxi1 << " " << dxi01 << std::endl;
         }
      }
      file_2nd_through_linear.close();
      std::cout << "Output file: " << filename_2nd_through_linear << std::endl;
   }
   /* -------------------------------------------------------------------------- */
   {
      auto filename = "points.dat";
      std::ofstream file(filename);
      if (!file.is_open()) {
         std::cerr << "Failed to open the output file." << std::endl;
         return 0;
      }
      std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
      int index = 0;
      for (auto& [x, y, z] : M)
         file << x << " " << y << " " << z << std::endl;

      std::cout << "Output file: " << filename << std::endl;
   }
}
