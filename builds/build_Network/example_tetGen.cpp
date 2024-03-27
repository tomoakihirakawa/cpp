#include <iostream>
#include "Network.hpp"
#include "tetgen.h"  // Make sure TetGen path is correctly set in your project settings

/*DOC_EXTRACT 1_0_tetGen

# 四面体の生成

## TetGenを使った四面体を生成

[https://wias-berlin.de/software/tetgen](https://wias-berlin.de/software/tetgen)

TetGenを使って四面体を生成し，Networkの四面体へと置き換える．

`tetgenbehavior`は，TetGenのオプションを設定するためのクラスで，`parse_commandline`関数を使ってオプションを設定する．
次のようなオプションがあり意味がある([https://wias-berlin.de/software/tetgen/switches.html](https://wias-berlin.de/software/tetgen/switches.html))：

| オプション | 意味 |
|:---:|:---:|
| `p` | PLC（Piecewise Linear Complex）を四面体化する． 他には，`r`（リージョン）や`y`（境界）などがある． |
| `q` | 最小radius-edge比を指定する．例えば，`q1.4`は最小radius-edge比1.4を指定する． |
| `a` | 最大四面体の体積制約を課す．例えば，`a50.`は最大体積50の四面体の体積制約を課す． |


現在のフォルダに`tetgen1.6.0`を置き，次のコマンドを実行すると，`libtet.a`が生成される．

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./tetgen1.6.0
make
```

これまで使っていたCMakeLists.txt（`./tetgen1.6.0/CMakeLists.txt`ではない）に次の行を追加する．

```cmake
target_link_libraries(${BASE_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/build_Network/libtet.a")
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
```

この`CMakelists.txt`を使って，TetGenを使うプログラムをビルドし，実行する．

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example_tetGen.cpp
make
./example_tetGen
```

<img src="./image.png" width="500px">

*/

int main() {

   /* -------------------------------------------------------------------------- */
   /*                        tetGenを使って四面体を生成する                          */
   /* -------------------------------------------------------------------------- */

   tetgenio in, out;
   tetgenbehavior b;

   // Load a .off file. Replace "socket.off" with your actual file path.
   if (!in.load_off("socket.off")) {
      std::cerr << "Error loading .off file" << std::endl;
      return 1;
   }
   b.parse_commandline("pq2.a50.");

   try {
      tetrahedralize(&b, &in, &out);
   } catch (int e) {
      std::cerr << "TetGen error: " << e << std::endl;
      return 1;
   }

   // From here, you can use the 'out' object which contains the generated tetrahedra.
   std::cout << "Generated " << out.numberoftetrahedra << " tetrahedra." << std::endl;

   /* -------------------------------------------------------------------------- */
   /*                         Networkの四面体へと置き換える                          */
   /* -------------------------------------------------------------------------- */

   auto object = new Network;

   for (int i = 0; i < out.numberofpoints; i++)
      new networkPoint(object, {out.pointlist[i * 3], out.pointlist[i * 3 + 1], out.pointlist[i * 3 + 2]});

   object->setGeometricProperties();
   object->makeBucketPoints(object->getScale() / 50.);

   int count = 0;
   PVDWriter pvd(_HOME_DIR_ + "/output/tetras.pvd");

   /* ----------------------------------- 出力 ----------------------------------- */
   auto output = [&count, &pvd, &object]() {
      auto filename = _HOME_DIR_ + "/output/tetras" + std::to_string(count) + ".vtu";
      std::ofstream ofs(filename);
      vtkUnstructuredGridWrite(ofs, object->getTetras());
      std::cout << "Wrote " << filename << std::endl;
      ofs.close();
      pvd.push(filename, count++);
   };

   /* --------------------------------- 置き換えループ -------------------------------- */
   for (int i = 0; i < out.numberoftetrahedra; i++) {

      /* -------------------------------------------------------------------------- */
      /*                                    置き換え                                  */
      /* -------------------------------------------------------------------------- */

      auto get4Points = [&]() -> std::array<networkPoint*, 4> {
         Tddd X;
         std::array<networkPoint*, 4> points;
         for (int j = 0; j < 4; j++) {
            auto index = out.tetrahedronlist[i * 4 + j];
            X = {out.pointlist[index * 3],
                 out.pointlist[index * 3 + 1],
                 out.pointlist[index * 3 + 2]};
            //! バケツから点を取得
            for (auto& p : object->BucketPoints.getBucket(X)) {
               if (Norm(p->X - X) < 1E-10) {
                  points[j] = p;
                  break;
               }
            }
         }
         return points;
      };

      genTetra(object, get4Points());

      /* -------------------------------------------------------------------------- */
      /*                                     出力                                    */
      /* -------------------------------------------------------------------------- */
      if (i % 1000 == 0) output();
   }
   output();
   pvd.output();

   return 0;
}
