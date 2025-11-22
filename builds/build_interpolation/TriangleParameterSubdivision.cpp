#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Network.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"

/*DOC_EXTRACT 1_0_0_interpolation

## 三角形を使った補間

### 三角分割

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=TriangleParameterSubdivision.cpp
make
./TriangleParameterSubdivision
```

* `SymmetricSubdivisionOfTriangle_00_10_01`で三角形を分割
* `SubdivideSquareIntoTriangles` で矩形領域を三角形に分割
* `SymmetricSubdivisionOfTriangle`は`SymmetricSubdivisionOfTriangle_00_10_01`を使って，任意頂点の三角形を分割する

`plot_parametric_subdivision.nb` で描画

| $`(\xi_0,\xi_1)`$, Range: $`(\xi_0,\xi_1)\in[0,1]\times[0,1-\xi_0]`$ | $`(\xi_0,\xi_1)`$, Range: $`(\xi_0,\xi_1)\in[0,1]\times[0,1]`$ | $`(\xi_0,\xi_1(1-\xi_0))`$, Range:$`(\xi_0,\xi_1)\in[0,1]\times[0,1]`$|

| a | b | c |
|:---:|:---:|:---:|
| <img src="output_TriangleParameterSubdivision.gif" width="300"> | <img src="output_SquareParameterSubdivision.gif" width="300"> | <img src="output_SquareParameterSubdivision_into_Triangle.gif" width="300"> |

$`(\xi_0,\xi_1)=(\xi_0,\eta(1-\xi_0))`$とすると，
$`\eta0`$のとき，$(\xi_0,\xi_1)=(\xi_0,0)$，また
$`\eta)=1`$のとき，$(\xi_0,\xi_1)=(\xi_0,1-\xi_0)$
となり，
$`(\xi_0,\eta)\in[0,1]\times[0,1]`$は
$`(\xi_0,\xi_1)=(\xi_0,\eta(1-\xi_0))\in[0,1]\times[0,1-\xi_0]`$に写像される．

このことは，形状関数を使った積分の際に利用できる．

*/

int main() {
   {
      std::ofstream file("./output_TriangleParameterSubdivision");
      for (auto t0t1 : SymmetricSubdivisionOfTriangle_00_10_01(5))  // ３頂点を返す
         file << t0t1[0] << "," << t0t1[1] << "," << t0t1[2] << std::endl;
      file.close();
   }

   {
      std::ofstream file("./output_SubdivideSquareIntoTriangles3D");
      T3Tddd V = {{{0, 0, 0}, {2, -1.5, 0}, {1, 2, 2}}};
      for (auto t0t1 : SymmetricSubdivisionOfTriangle(V, 5))  // ３頂点を返す
         file << t0t1[0] << "," << t0t1[1] << "," << t0t1[2] << std::endl;
      file.close();
   }

   {
      std::ofstream file("./output_SubdivideSquareIntoTriangles2D");
      T3Tdd V = {{{0, 0}, {0.25, 0}, {0, 0.25}}};
      for (auto t0t1 : SymmetricSubdivisionOfTriangle(V, 5))  // ３頂点を返す
         file << t0t1[0] << "," << t0t1[1] << "," << t0t1[2] << std::endl;
      file.close();
   }

   {
      std::ofstream file("./output_SquareParameterSubdivision");
      for (auto t0t1 : SubdivideSquareIntoTriangles(5, 8))  // ３頂点を返す
         file << t0t1[0] << "," << t0t1[1] << "," << t0t1[2] << std::endl;
      file.close();
   }

   std::ofstream file2("./output_TriangleParameterSubdivision_Modified");
   for (auto t0t1 : SubdivideSquareIntoTriangles(5, 8))  // ３頂点を返す
   {
      for (auto [t0, t1] : t0t1) {
         auto t0t1t2 = ModTriShape<3>(t0, t1);
         file2 << t0t1t2[0] << "," << t0t1t2[1] << ",";
      }
      file2 << std::endl;
   }
   file2.close();
};