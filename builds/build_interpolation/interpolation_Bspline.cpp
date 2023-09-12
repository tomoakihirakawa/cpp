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

/*DOC_EXTRACT 0_2_0_interpolation

## B-spline補間

与えられたデータ点を通る多項式を求める方法の一つにB-spline補間がある．

### 実行方法

```sh
$ cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=interpolation_Bspline.cpp
$ make
$ ./interpolation_Bspline
$ gnuplot bspline_plot.gnu
```

### コード

\ref{interpolation:Bspline}{Bspline基底関数}を用いて，B-spline補間を行う．

`InterpolationBspline`は，`std::vector<double>`または`std::vector<std::array<double,N>>`を引数に取ることができる．

```cpp
// example for 1D data
std::vector<double> X;
InterpolationBspline intpX(5, abscissas, X);
```

![sample_body_movement_bspline.png](sample_bspline.png)

```cpp
// example for 2D data
std::vector<std::arrray<double,2>> XY;
InterpolationBspline intpXY(5, abscissas, XY);
```

または，クラスを使いまわしたい場合，`set`メンバ関数を用いて，データをセットすることもできる．

```cpp
InterpolationBspline<std::array<double, 2>> intpXY;
intpXY.set(5, abscissas, XY);
```


![sample_body_movement_bspline.png](sample_body_movement_bspline.png)

\insert{RBF}

*/

int main() {

   /* -------------------------------------------------------------------------- */
   {

      std::ifstream file("bspline_sample_data.dat");
      if (!file.is_open()) {
         std::cerr << "Failed to open the input file." << std::endl;
         return 0;
      }

      std::vector<double> abscissas, values;
      std::string line;
      double t, x;
      while (std::getline(file, line)) {
         if (line[0] == '#') continue;
         std::istringstream iss(line);
         iss >> t >> x;
         abscissas.push_back(t);
         values.push_back(x);
      }
      file.close();

      InterpolationBspline intB(5, abscissas, values);

      std::ofstream output("bspline_interpolated_data.dat");
      if (!output.is_open()) {
         std::cerr << "Failed to open the output file." << std::endl;
         return 0;
      }
      output << "# x y z" << std::endl;
      for (auto t : Subdivide<300>(0., 10.)) {
         output << t << " " << intB(t) << " " << intB.D(t) << std::endl;
      }
      output.close();
   }
   /* -------------------------------------------------------------------------- */
   InterpolationBspline<std::array<double, 2>> intpXY;

   for (const std::string& S : {"A", "B", "C"}) {
      std::ifstream file("../build_pybind11/body" + S + ".dat");
      if (!file.is_open()) throw std::runtime_error("Failed to open the input file.");

      std::vector<double> abscissas;
      std::vector<std::array<double, 2>> XY;
      std::string line;
      double t, x, y, z, angle;

      while (std::getline(file, line)) {
         if (line[0] == '#') continue;
         replace(line.begin(), line.end(), ',', ' ');
         std::istringstream iss(line);
         iss >> t >> x >> y >> z >> angle;
         abscissas.push_back(t);
         XY.push_back({x, y});
      }
      file.close();

      intpXY.set(5, abscissas, XY);

      std::ofstream output("bspline_interpolated_body" + S + ".dat");
      if (!output.is_open()) throw std::runtime_error("Failed to open the input file.");
      output << "# t x y" << std::endl;
      for (auto t : Subdivide<300>(0., 4.)) {
         output << t << " " << intpXY(t)[0] << " " << intpXY(t)[1] << std::endl;
      }
      output.close();
   }
   /* -------------------------------------------------------------------------- */
   return 0;
}