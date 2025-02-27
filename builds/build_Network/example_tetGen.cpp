#include <iostream>
#include "tetgen1.6.0/tetgen.h"  // Make sure TetGen path is correctly set in your project settings
//
#include "Network.hpp"

/*DOC_EXTRACT 1_0_tetGen

# 四面体の生成

## TetGenを使った四面体を生成

* TETLIBRARYを有効にしない
* `terminatetetgen`を修正する
* コンパイラを合わせる

tetgen.h内のTETLIBRARYを有効にすると，コンパイルができなかった．errorによって，abortすることを防ぐために，直接throwするよう以下のように`terminatetetgen`関数を修正した．

```cpp
inline void terminatetetgen(tetgenmesh* m, int x) {
   throw x;
}
```

```cmake
# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 3.5)

project(tetgen)

set(CXX_VERSIONS 13)

# Find C++ compiler
foreach(ver IN LISTS CXX_VERSIONS)
  unset(CXX_COMPILER CACHE)
  string(CONCAT CXX "g++-" ${ver})
  find_program(CXX_COMPILER ${CXX})
  if(CXX_COMPILER)
    message(STATUS "${Green}Found ${CXX}: ${Magenta}${CXX_COMPILER}${ColourReset}")
    set(CMAKE_CXX_COMPILER ${CXX_COMPILER})
    break()
  endif()
endforeach()


# Add an executable to the project using the specified source files.
add_executable(tetgen tetgen.cxx predicates.cxx)

#Add a library to the project using the specified source files.
# In Linux/Unix, it will creates the libtet.a
add_library(tet STATIC tetgen.cxx predicates.cxx)
set_target_properties(tet PROPERTIES PUBLIC_HEADER tetgen.h)

#Set properties on a target.
#We use this here to set -DTETLIBRARY for when compiling the
#library
set_target_properties(tet PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)

if(NOT TETGEN_SKIP_INSTALL)
    include(GNUInstallDirs)
    install(TARGETS tet
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    install(TARGETS tetgen RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
```


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
./example_tetGen pq2.a10. coil.off coil_pq2a10
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
| q2 | 最小circum radius-the min edge比を2に設定する．$`\frac{circum R}{L_{\rm min}} > 2`$
| a50. | 四面体の最大体積を50に制限します．例えば，a100.とすると最大体積が100に制限される．|


<figure>
<img src="./image_tetgen_comparison.png" width="600px">
<figcaption>pq2.a50, pq1.a50, pq1.a0.00005の比較</figcaption>
</figure>

*/

int main(int argc, char* argv[]) {

   std::string filename, out_file, parse_command;
   std::cout << "argc: " << argc << std::endl;
   if (argc == 4) {
      for (int i = 0; i < argc; i++)
         std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
      filename = argv[1];
      parse_command = argv[2];
      out_file = argv[3];
   } else
      throw std::invalid_argument("Invalid argument");

   auto object = new Network(filename);
   object->makeBucketPoints(object->getScale() / 50.);
   object->setGeometricProperties();
   object->tetrahedralize(parse_command);

   /* ----------------------------------- 出力用関数 ----------------------------------- */

   {
      auto filename = "./output/" + std::string(out_file) + ".vtu";
      auto filename_surface = "./output/" + std::string(out_file) + "_surface.vtu";
      std::ofstream ofs(filename);
      std::ofstream ofs_surface(filename_surface);
      std::vector<std::tuple<std::string, std::unordered_map<networkPoint*, std::variant<double, std::array<double, 3>>>>> data;
      std::unordered_map<networkPoint*, std::variant<double, std::array<double, 3>>> data_vector;
      for (auto& p : object->getPoints())
         data_vector[p] = p->X;
      data.push_back({"X", data_vector});

      data_vector.clear();
      for (auto& p : object->getPoints())
         data_vector[p] = (double)p->getFaces().size();
      data.push_back({"face_size", data_vector});

      data_vector.clear();
      for (auto& p : object->getPoints())
         data_vector[p] = (double)p->getNeighbors().size();
      data.push_back({"neighbors_size", data_vector});

      data_vector.clear();
      for (auto& p : object->getPoints())
         data_vector[p] = (double)p->SurfaceQ();
      data.push_back({"any_of_no_tetra", data_vector});

      vtkUnstructuredGridWrite(ofs, object->getTetras(), data);
      std::cout << "Wrote " << filename << std::endl;
      ofs.close();
      vtkPolygonWrite(ofs_surface, object->getSurfaces(), data);
      std::cout << "Wrote " << filename_surface << std::endl;
      ofs_surface.close();
   }
}
