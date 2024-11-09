#include <iostream>
#include "Network.hpp"
#include "tetgen1.6.0/tetgen.h"  // Make sure TetGen path is correctly set in your project settings

/*DOC_EXTRACT 1_0_tetGen

# 四面体の生成

## TetGenを使った四面体を生成

[https://wias-berlin.de/software/tetgen](https://wias-berlin.de/software/tetgen)

TetGenを使って四面体を生成し，Networkの四面体へと置き換え，出力するプログラム．
現在のフォルダに`tetgen1.6.0`を置き（tetgen1.6.0内にCMakelists.txtが保存されている．），次のコマンドを実行すると，`libtet.a`が生成される．
`.a`は，`.o`ファイルをまとめたアーカイブファイルである．

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./tetgen1.6.0
make
```

上のアーカイブを利用するメインのcppプログラムのCMakeLists.txt（`./tetgen1.6.0/CMakeLists.txt`ではない）に次の行を追加する．

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

### `tetgenbehavior`クラス

`tetgenbehavior`クラスは，TetGenのメッシュ生成オプションを設定するために使用され，`parse_commandline`関数を通じてオプションを指定できる．
この関数を使うことで，メッシュ生成の際に必要な条件や制約を細かく調整できる．

例として，pq2.a50.が指定された場合，以下のオプションが適用される．

```cpp
tetgenbehavior b;
b.parse_commandline("pq2.a50.");
```

| オプション | 意味 |
|:---:|:---:|
| p | PLC（Piecewise Linear Complex）を四面体メッシュ化する．その他に，再メッシュ用のrや，境界ポイントの保持を行うyなどのオプションもある．|
| q2 | 最小radius-edge比を2に設定し，品質の高い四面体を生成する．例えば，q1.4なら比率を1.4に設定する．|
| a50. | 四面体の最大体積を50に制限します．例えば，a100.とすると最大体積が100に制限される．|


<figure>
<img src="./image_tetgen_comparison.png" width="600px">
<figcaption>pq2.a50, pq1.a50, pq1.a0.00005の比較</figcaption>
</figure>

*/

int main(int argc, char* argv[]) {

   char *name, *parse_command, *out_file;
   std::cout << "argc: " << argc << std::endl;
   if (argc == 4) {
      for (int i = 0; i < argc; i++)
         std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
      parse_command = argv[1];
      name = argv[2];
      out_file = argv[3];
   } else
      throw std::invalid_argument("Invalid argument");

   /* -------------------------------------------------------------------------- */
   /*                        tetGenを使って四面体を生成する                          */
   /* -------------------------------------------------------------------------- */

   tetgenio in, out;
   tetgenbehavior b;

   // Load a .off file. Replace "socket.off" with your actual file path.

   // if file name incldeu ".off", load_off is used
   // if file name incldeu ".stl", load_stl is used
   std::string filename(name);
   if (filename.contains(".off")) {
      if (!in.load_off(name)) {
         std::cerr << "Error loading .off file" << std::endl;
         return 1;
      }
   } else if (filename.contains(".stl")) {
      if (!in.load_stl(name)) {
         std::cerr << "Error loading .stl file" << std::endl;
         return 1;
      }
   } else if (filename.contains(".obj")) {
      if (!in.load_obj(name)) {
         std::cerr << "Error loading .obj file" << std::endl;
         return 1;
      }
   } else {
      std::cerr << "Error loading file" << std::endl;
      return 1;
   }
   b.parse_commandline(parse_command);

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
   PVDWriter pvd(_HOME_DIR_ + "/output/" + std::string(out_file) + ".pvd");

   /* ----------------------------------- 出力 ----------------------------------- */
   auto output = [&]() {
      auto filename = _HOME_DIR_ + "/output/" + std::string(out_file) + "_" + std::to_string(count) + ".vtu";
      std::ofstream ofs(filename);
      vtkUnstructuredGridWrite(ofs, object->getTetras());
      std::cout << "Wrote " << filename << std::endl;
      ofs.close();
      pvd.push(filename, count++);
   };

   /* --------------------------------- 置き換えループ -------------------------------- */
   for (int i = 0; i < out.numberoftetrahedra; i++) {

      std::array<networkPoint*, 4> points;
      for (int j = 0; j < 4; j++) {
         auto index = out.tetrahedronlist[i * 4 + j];
         Tddd X = {out.pointlist[index * 3], out.pointlist[index * 3 + 1], out.pointlist[index * 3 + 2]};
         //! バケツから点を取得
         for (auto& p : object->BucketPoints.getData(X)) {
            if (Norm(p->X - X) < 1E-10) {
               points[j] = p;
               break;
            }
         }
      }

      genTetra(object, points);

      if (i % 1000 == 0)
         output();
   }
   output();
   pvd.output();

   return 0;
}
