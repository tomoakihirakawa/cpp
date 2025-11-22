# 🐋 `Network` 

節点に隣接する節点や辺や要素の情報を効率的に取得するためには，接続関係を管理し続ける必要がある．`Network`クラスは，接続関係の情報を保持し，その情報をもとに相互にアクセスするための機能を提供する．

* 節点や辺や面の相互アクセス
* メッシュの細分化
* `obj`や`off`オブジェクトデータの読み込みと出力

## ⛵ 点・線・面の接続関係とその整理 

1. `networkFace->Lines`を設定
2. `networkFace->setPoints()`は，`networkFace->Lines`が設定されていることを前提として，`networkFace->Points`と`networkFace->PLPLPL`を設定する．
3. `Network::setGeometricProperties()`は，`f->setGeometricProperties(ToX(f->setPoints()))`を実行している．

## ⛵ 3Dファイルの読み込みと出力 

### 🪼 読み込み `Network` 

[:material-microsoft-visual-studio-code:Network::constructor](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/include/Network.hpp:3587)では，引数として，**OFFファイル**または**OBJファイル**をあたえることができる．
`Load3DFile`クラスを使ってデータを読み込み，`Network`クラスを作成している．

```cpp
auto obj = new Network("./bunny.obj");//ファイルからNetworkオブジェクトを作成
```

### 🪼 出力 `vtkPolygonWrite`，`vtkUnstructuredGridWrite` 

`Network`クラスは，`getFaces`メンバ関数を使って簡単に面の情報を取得できる．

`vtkPolygonWrite`を使うと，`Network`クラスの面の情報を，`vtp`ファイルとして出力できる．
`vtkPolygonWrite`には，`ofstream`と，`std::vector<networkFace*>`や`std::vector<networkLine*>`などを渡し，出力できる．

#### 🪸 面の出力 

面や線や点の情報を出力する場合は，`vtkPolygonWrite`を使う．

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getFaces());
ofs.close();
```

#### 🪸 線の出力 

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getEdges());
ofs.close();
```

#### 🪸 点の出力 

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getPoints());
ofs.close();
```

#### 🪸 四面体の出力 

四面体のような内部構造を持つデータを出力する場合は，`vtkUnstructuredGridWrite`を使う．

```cpp
auto obj = new Network("./input/bunny.off");
obj->tetrahedralize();
auto data1 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
auto data2 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
auto data3 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
for (const auto& p : obj->getPoints()) {
data1[p] = p->X[0];
data2[p] = p->X;
data3[p] = (double)p->Tetras.size();
}
std::vector<std::tuple<std::string, DataMap>> data = {{"x", data1}, {"xyz", data2}, {"tetra_size", data3}};
std::ofstream ofs("./outptut/tetras.vtu");
vtkUnstructuredGridWrite(ofs, obj->getTetras(), data);
ofs.close();
```

[:material-microsoft-visual-studio-code:add_data_default_name](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_load_3d_file.cpp:239)，点に値を付与し，vtpとして出力することもできる．
また，[:material-microsoft-visual-studio-code:add_data_custom_name](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_load_3d_file.cpp:269)を付けることもできる．

#### 🪸 実行方法 

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_load_3d_file.cpp
make
./example0_load_3d_file
```

[:material-microsoft-visual-studio-code:example0_load_3d_file.cpp#L1](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_load_3d_file.cpp:1)

---
### 🪼 `PVDWriter`を使ったpvdファイルの作成方法 

pvdファイルは，ファイルと時間をセットにしてまとめ，paraview上で，3Dファイルのアニメーションを再生するためのファイル．

```cpp
PVDWriter pvd("./bunny_obj.pvd");//出力するpvdファイル名を指定しクラスを作成
pvd.push(filename, time);//`filename`には，`vtp`ファイルなどの3Dファイル名を，`time`には，そのファイルの時間を指定
pvd.output();//最後にpvdファイルを出力
```

| 面のアニメーション | 線のアニメーション |
|:---------------:|:---------------:|
| <img src="sample.gif" width="500px"> | <img src="sample_line.gif" width="500px"> |

💡 QuickTimeで作成したmovファイルをgifに変換するには，次のようにする．

```sh
ffmpeg -i line.mov -filter_complex "[0:v] fps=30, scale=iw*0.5:ih*0.5 [v]" -map "[v]" sample_line.gif
```

[:material-microsoft-visual-studio-code:example0_load_3d_file.cpp#L155](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_load_3d_file.cpp:155)

---
## ⛵ 四面体の操作 

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_manipulation_tetrahedron.cpp
make
./example0_manipulation_tetrahedron
```

* [:material-microsoft-visual-studio-code:make_bucket_tetras](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_manipulation_tetrahedron.cpp:83)，空間分割のためのバケットを作成し，四面体をバケットに登録する例
* [:material-microsoft-visual-studio-code:edge_flip](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_manipulation_tetrahedron.cpp:104)，四面体を持つ表面のエッジフリップテストの例
* [:material-microsoft-visual-studio-code:test_insideQ](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_manipulation_tetrahedron.cpp:139)，座標が四面体の内部か外部かの判定テストの例

[:material-microsoft-visual-studio-code:example0_manipulation_tetrahedron.cpp#L1](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_manipulation_tetrahedron.cpp:1)

## ⛵ 四面体の操作 


```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_manipulation_tetrahedron_simple_simulation.cpp
make
./example0_manipulation_tetrahedron_simple_simulation
```

[:material-microsoft-visual-studio-code:example0_manipulation_tetrahedron_simple_simulation.cpp#L1](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_manipulation_tetrahedron_simple_simulation.cpp:1)

---
## ⛵ ２次補間 

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_quadratic_interpolation.cpp
make
./example0_quadratic_interpolation
```

[:material-microsoft-visual-studio-code:example0_quadratic_interpolation.cpp#L1](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example0_quadratic_interpolation.cpp:1)

---
## ⛵ 階層のある空間分割（木構造） 

シンプルな空間分割クラスを拡張し，木構造による空間分割を試みる．

`has_tree`が`true`の場合，`children`には`Bucket`クラスのポインタが格納される．
`children[i][j][k]`には，上のレベルの`data[i][j][k]`のデータが引き継がれている．
つまり，`children[i][j][k]`は，`data[i][j][k]`のデータをさらに分割したものである．
デフォルトでは，`children[i][j][k]`は内部に８つの`data`を持つ:

`data[0][0][0]`，`data[0][0][1]`，`data[0][1][0]`，`data[0][1][1]`，`data[1][0][0]`，`data[1][0][1]`，`data[1][1][0]`，`data[1][1][1]`．

[buckets_generateTree](not found)は，
バウンディングボックスを範囲と，それを分割する幅を指定する．
分割数を指定するよりも，この方法のように分割幅を指定する方が，自分はわかりやすい．

```cpp
children[i][j][k] = std::make_shared<Buckets<T>>(bounds, this->dL * 0.5 + 1e-10);
```

<img src="example2_tree_faster.gif" width="500px">

レベル０が生成したレベル１のバケットに保存された点を示しており，
白い線は，１階層上のレベル０のバケットの境界を示している．

[:material-microsoft-visual-studio-code:example2_tree.cpp#L2](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example2_tree.cpp:2)

---
# 🐋 四面体の生成 

## ⛵ TetGenを使った四面体を生成 

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
# 🐋 Set  the minimum  required version  of cmake  for a  project. 
cmake_minimum_required(VERSION 3.5)

project(tetgen)

set(CXX_VERSIONS 13)

# 🐋 Find C++ compiler 
foreach(ver IN LISTS CXX_VERSIONS)
unset(CXX_COMPILER CACHE)
string(CONCAT CXX "g++-" ${ver})
find _program(CXX _COMPILER ${CXX})
if(CXX_COMPILER)
message(STATUS "${Green}Found ${CXX}: ${Magenta}${CXX_COMPILER}${ColourReset}")
set(CMAKE _CXX _COMPILER ${CXX_COMPILER})
break()
endif()
endforeach()


# 🐋 Add an executable to the project using the specified source files. 
add_executable(tetgen tetgen.cxx predicates.cxx)

#Add a library to the project using the specified source files.
# 🐋 In Linux/Unix, it will creates the libtet.a 
add_library(tet STATIC tetgen.cxx predicates.cxx)
set_target_properties(tet PROPERTIES PUBLIC_HEADER tetgen.h)

#Set properties on a target.
#We use this here to set -DTETLIBRARY for when compiling the
#library
set_target_properties(tet PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)

if(NOT TETGEN_SKIP_INSTALL)
include(GNUInstallDirs)
install(TARGETS tet
ARCHIVE DESTINATION ${CMAKE _INSTALL _LIBDIR}
PUBLIC _HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(TARGETS tetgen RUNTIME DESTINATION ${CMAKE _INSTALL _BINDIR})
endif()
```


[https://wias-berlin.de/software/tetgen](https://wias-berlin.de/software/tetgen)

TetGenを使って四面体を生成し，Networkの四面体へと置き換え，出力するプログラム．
現在のフォルダに`tetgen1.6.0`を置き（tetgen1.6.0内にCMakelists.txtが保存されている．），次のコマンドを実行すると，`libtet.a`が生成される．
`.a`は，`.o`ファイルをまとめたアーカイブファイルである．

```shell
sh clean
cmake -DCMAKE _BUILD _TYPE=Release ./tetgen1.6.0
make
```

上のアーカイブを利用するメインのcppプログラムのCMakeLists.txt（`./tetgen1.6.0/CMakeLists.txt`ではない）に次の行を追加する．

```cmake
target _link _libraries(${BASE_NAME} "${CMAKE _CURRENT _SOURCE _DIR}/build _Network/libtet.a")
include _directories(${CMAKE_CURRENT_SOURCE_DIR})
```

この`CMakelists.txt`を使って，TetGenを使うプログラムをビルドし，実行する．

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example_tetGen.cpp
make
./example_tetGen pq2.a10. coil.off coil_pq2a10
```

### 🪼 `tetgenbehavior`クラス 

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
| q2 | 最小circum radius-the min edge比を2に設定する．$`\frac{circum R}{L _{\rm min}} > 2`$
| a50. | 四面体の最大体積を50に制限します．例えば，a100.とすると最大体積が100に制限される．|


<figure>
<img src="./image_tetgen_comparison.png" width="600px">
<figcaption>pq2.a50, pq1.a50, pq1.a0.00005の比較</figcaption>
</figure>

[:material-microsoft-visual-studio-code:example_tetGen.cpp#L6](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Network/example_tetGen.cpp:6)

---
