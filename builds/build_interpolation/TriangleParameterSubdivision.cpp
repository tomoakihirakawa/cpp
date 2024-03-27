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

* `SubdivideTriangleIntoTriangles` で三角形を分割
* `SubdivideSquareIntoTriangles` で矩形領域を三角形に分割

`plot_parametric_subdivision.nb` で描画

<img src="output_TriangleParameterSubdivision.gif" width="400">

<img src="output_SquareParameterSubdivision.gif" width="400">

`ModTriShape`を使うと，(t0,t1)=([0,1],[0,1])領域を(xi0,xi1)=([0,1],[0,1-t0])の三角形に変換できる．

<img src="output_SquareParameterSubdivision_into_Triangle.gif" width="400">

*/

int main() {
   {
      std::ofstream file("./output_TriangleParameterSubdivision");
      for (auto t0t1 : SubdivideTriangleIntoTriangles(5))  // ３頂点を返す
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