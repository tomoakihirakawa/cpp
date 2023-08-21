# Contents

- [🐋`Network`](#🐋`Network`)
    - [⛵️3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法](#⛵️3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法)
        - [🪸読み込み `Network`](#🪸読み込み-`Network`)
        - [🪸出力 `vtkPolygonWrite`](#🪸出力-`vtkPolygonWrite`)
        - [🪸`PVDWriter`を使ったpvdファイルの作成方法](#🪸`PVDWriter`を使ったpvdファイルの作成方法)
    - [⛵️CGALを使って四面体を生成する](#⛵️CGALを使って四面体を生成する)
    - [⛵️四面体を生成（制約付き四面分割 constrained tetrahedralization）](#⛵️四面体を生成（制約付き四面分割-constrained-tetrahedralization）)
        - [🪸Advancing Front Algorithm](#🪸Advancing-Front-Algorithm)
    - [⛵️四面体を生成（制約付き四面分割 constrained tetrahedralization）](#⛵️四面体を生成（制約付き四面分割-constrained-tetrahedralization）)
        - [🪸Advancing Front Algorithm](#🪸Advancing-Front-Algorithm)


---
# 🐋`Network` 

## ⛵️3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法 

### 🪸読み込み `Network` 

[Networkのコンストラクタ](../../include/Network.hpp#L3875)では，拡張子から，
与えられたファイルが，

* OFFファイル
* OBJファイル

かチェクして，`Load3DFile`クラスを使ってデータを読み込み`Network`クラスとして読み込む．

### 🪸出力 `vtkPolygonWrite` 

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


[./expample0_load_3d_file.cpp#L4](./expample0_load_3d_file.cpp#L4)


---
### 🪸`PVDWriter`を使ったpvdファイルの作成方法 

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


[./expample0_load_3d_file.cpp#L53](./expample0_load_3d_file.cpp#L53)


---
## ⛵️CGALを使って四面体を生成する 

```
$ brew install CGAL
```


[./expample1_generate_tetra_using_CGAL.cpp#L2](./expample1_generate_tetra_using_CGAL.cpp#L2)


## ⛵️四面体を生成（制約付き四面分割 constrained tetrahedralization） 

* PLC: piecewise linear complex
* CDT: constrained Delaunay triangulation

CDTの生成法には，主に２つの方法がある[Schewchuk 2002](Schewchuk 2002)：

* naive gift wrapping algorithm (これはadvancing front algorithmとも呼ばれるものと同じだろう)
* sweep algorithm

### 🪸Advancing Front Algorithm


[./expample2_generate_tetra_constrained2.cpp#L2](./expample2_generate_tetra_constrained2.cpp#L2)


## ⛵️四面体を生成（制約付き四面分割 constrained tetrahedralization） 

* PLC: piecewise linear complex
* CDT: constrained Delaunay triangulation

CDTの生成法には，主に２つの方法がある[Schewchuk 2002](Schewchuk 2002)：

* naive gift wrapping algorithm (これはadvancing front algorithmとも呼ばれるものと同じだろう)
* sweep algorithm

### 🪸Advancing Front Algorithm


[./expample3_generate_tetra_constrained.cpp#L2](./expample3_generate_tetra_constrained.cpp#L2)


---
