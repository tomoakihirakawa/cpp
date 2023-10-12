# Contents
- [🐋 `Network`](#🐋-`Network`)
    - [⛵ 3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法](#⛵-3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法)
        - [🪼 読み込み `Network`](#🪼-読み込み-`Network`)
        - [🪼 出力 `vtkPolygonWrite`](#🪼-出力-`vtkPolygonWrite`)
- [🐋 空間分割（space_partitioning）](#🐋-空間分割（space_partitioning）)
    - [⛵ 等間隔のシンプルな空間分割](#⛵-等間隔のシンプルな空間分割)
        - [🪼 例](#🪼-例)
    - [⛵ 階層のある空間分割（木構造）](#⛵-階層のある空間分割（木構造）)
- [🐋 CGALを使って四面体を生成する](#🐋-CGALを使って四面体を生成する)
    - [⛵ CGALを使って四面体を生成する](#⛵-CGALを使って四面体を生成する)
    - [⛵ 四面体を生成（制約付き四面分割 constrained tetrahedralization）](#⛵-四面体を生成（制約付き四面分割-constrained-tetrahedralization）)
        - [🪼 `PVDWriter`を使ったpvdファイルの作成方法](#🪼-`PVDWriter`を使ったpvdファイルの作成方法)


---
# 🐋 `Network` 

## ⛵ 3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法 

### 🪼 読み込み `Network` 

[Networkのコンストラクタ](../../include/Network.hpp#L3944)では，拡張子から，
与えられたファイルが，

* OFFファイル
* OBJファイル

かチェクして，`Load3DFile`クラスを使ってデータを読み込み`Network`クラスとして読み込む．

### 🪼 出力 `vtkPolygonWrite` 

`vtkPolygonWrite`には，`ofstream`と，`std::vector<networkFace*>`や`std::vector<networkTetra*>`などを渡し，出力できる．

```cpp
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getFaces());
```

![sample.png](sample.png)

```shell
$ cmake -DCMAKE _BUILD _TYPE=Release ../ -DSOURCE _FILE=load _3d _file.cpp
$ make
$ ./load_3d_file
```

[./example0_load_3d_file.cpp#L4](./example0_load_3d_file.cpp#L4)

---
# 🐋 空間分割（space_partitioning） 

## ⛵ 等間隔のシンプルな空間分割 

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_space_partitioning.cpp
make
./example1_space_partitioning
```

<!-- Key coordinatebounds not found -->

<!-- Key space_partitioning not found -->

### 🪼 例 

この例では，うさぎの３Dモデルを空間分割する．
配列させたバケット内に，うさぎの点または面が含まれるかを判定し，バケットに保存する．

ただ，面は広がりがあるので，複数のバケットに含まれることがある．
面と交わる全バケットを簡単に確実に見つける方法は，現在のところ思いつかない．
なので，今の所は，面を無数の点に分けて，各点を含むバケットに面を保存することで対応している．

![example1_space_partitioning.gif](example1_space_partitioning.gif)

[./example1_space_partitioning.cpp#L6](./example1_space_partitioning.cpp#L6)

---
## ⛵ 階層のある空間分割（木構造） 

シンプルな空間分割クラスを拡張し，木構造による空間分割を試みる．

`has_tree`が`true`の場合，`buckets`には`Bucket`クラスのポインタが格納される．
`buckets[i][j][k]`には，上のレベルの`data[i][j][k]`のデータが引き継がれている．
つまり，`buckets[i][j][k]`は，`data[i][j][k]`のデータをさらに分割したものである．
デフォルトでは，`buckets[i][j][k]`は内部に８つの`data`を持つ:

`data[0][0][0]`，`data[0][0][1]`，`data[0][1][0]`，`data[0][1][1]`，`data[1][0][0]`，`data[1][0][1]`，`data[1][1][0]`，`data[1][1][1]`．

[このツリー生成方法](../../include/lib_spatial_partitioning.hpp#L109)は，
バウンディングボックスを範囲と，それを分割する幅を指定する．
分割数を指定するよりも，この方法のように分割幅を指定する方が，自分はわかりやすい．

```cpp
buckets[i][j][k] = std::make_shared<Buckets<T>>(bounds, this->dL * 0.5 + 1e-10);
```

![example2_tree_faster.gif](example2_tree_faster.gif)

[./example2_tree.cpp#L2](./example2_tree.cpp#L2)

---
# 🐋 CGALを使って四面体を生成する 

## ⛵ CGALを使って四面体を生成する 

```shell
brew install CGAL
```

[./example1_generate_tetra_using_CGAL.cpp#L1](./example1_generate_tetra_using_CGAL.cpp#L1)

---
## ⛵ 四面体を生成（制約付き四面分割 constrained tetrahedralization） 

* PLC: piecewise linear complex
* CDT: constrained Delaunay triangulation

CDTの生成法には，主に２つの方法がある[Schewchuk 2002](Schewchuk 2002)：

* naive gift wrapping algorithm (これはadvancing front algorithmとも呼ばれるものと同じだろう)
* sweep algorithm


[杉原厚吉,計算幾何学](杉原厚吉,計算幾何学)によれば，ドロネー四面体分割以外に，綺麗な四面体分割を作成する方法はほとんど知られていないらしい．
四面体分割は，三角分割の場合のように，最小内角最大性が成り立たたず，スリーバー（sliver）と呼ばれる，外接円が大きくないものの潰れた悪い四面体が作られる可能性がある．
このスリーバーをうまく削除することが重要となる．

[./example2_generate_tetra_constrained2.cpp#L2](./example2_generate_tetra_constrained2.cpp#L2)

---
### 🪼 `PVDWriter`を使ったpvdファイルの作成方法 

pvdファイルは，ファイルと時間をセットにしてまとめ，paraview上で，3Dファイルのアニメーションを再生するためのファイルである．

次のようにして，出力するpvdファイル名を指定しクラスを作成する．

```cpp
PVDWriter pvd("./bunny_obj.pvd");
```

`filename`には，`vtp`ファイルなどの3Dファイル名を，`time`には，そのファイルの時間を指定する：

```cpp
pvd.push(filename, time);
```

最後にpvdファイルを出力する．

```cpp
pvd.output();
```

![sample.gif](sample.gif)

[./example0_load_3d_file.cpp#L53](./example0_load_3d_file.cpp#L53)

---
