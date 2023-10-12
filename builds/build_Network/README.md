# Contents
    - [⛵ CGALを使って四面体を生成する](#⛵-CGALを使って四面体を生成する)
    - [⛵ 四面体を生成（制約付き四面分割 constrained tetrahedralization）](#⛵-四面体を生成（制約付き四面分割-constrained-tetrahedralization）)
- [🐋 空間分割（space_partitioning）](#🐋-空間分割（space_partitioning）)
        - [🪼 🪼 CoordinateBounds クラスについて](#🪼-🪼-CoordinateBounds-クラスについて)
            - [🐚 🐚 概要](#🐚-🐚-概要)
            - [🐚 🐚 メンバ変数](#🐚-🐚-メンバ変数)
            - [🐚 🐚 メンバ関数](#🐚-🐚-メンバ関数)
            - [🐚 🐚 オペレーター](#🐚-🐚-オペレーター)
            - [🐚 🐚 有用性](#🐚-🐚-有用性)
    - [⛵ ⛵ Bucket クラスの説明](#⛵-⛵-Bucket-クラスの説明)
        - [🪼 🪼 型エイリアス](#🪼-🪼-型エイリアス)
        - [🪼 🪼 メンバ変数](#🪼-🪼-メンバ変数)
        - [🪼 🪼 メソッド](#🪼-🪼-メソッド)
            - [🐚 🐚 初期化関連](#🐚-🐚-初期化関連)
            - [🐚 🐚 インデックス変換](#🐚-🐚-インデックス変換)
            - [🐚 🐚 データ追加・削除](#🐚-🐚-データ追加・削除)
            - [🐚 🐚 その他](#🐚-🐚-その他)
        - [🪼 🪼 使用例](#🪼-🪼-使用例)
        - [🪼 例](#🪼-例)
- [🐋 木構造による空間分割](#🐋-木構造による空間分割)
- [🐋 `Network`](#🐋-`Network`)
    - [⛵ 3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法](#⛵-3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法)
        - [🪼 読み込み `Network`](#🪼-読み込み-`Network`)
        - [🪼 出力 `vtkPolygonWrite`](#🪼-出力-`vtkPolygonWrite`)
        - [🪼 `PVDWriter`を使ったpvdファイルの作成方法](#🪼-`PVDWriter`を使ったpvdファイルの作成方法)


---
## ⛵ CGALを使って四面体を生成する 

```
$ brew install CGAL
```

[./example1_generate_tetra_using_CGAL.cpp#L2](./example1_generate_tetra_using_CGAL.cpp#L2)

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
# 🐋 空間分割（space_partitioning） 

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_space_partitioning.cpp
make
./example1_space_partitioning
```

### 🪼 🪼 CoordinateBounds クラスについて  

#### 🐚 🐚 概要  
`CoordinateBounds`クラスは，
3次元の座標領域（バウンディングボックス）を管理するためのクラスである．

#### 🐚 🐚 メンバ変数  

| 変数名 | 型 | 説明 |
|--------|----|------|
| bounds | std::array<std::array<double,2>,3> | x, y, zそれぞれの範囲を保持する |
| X      | std::array<double,3>  | 中心座標を保持する |

#### 🐚 🐚 メンバ関数  

| メソッド名     | 戻り値型 | 説明 |
|--------------|----------|------|
| scaledBounds  | std::array<std::array<double,2>,3>    | 指定されたスケールでバウンディングボックスを拡大・縮小する |
| setBounds     | void     | バウンディングボックスを設定する（オーバーロードあり） |
| getXtuple     | const std::array<double,3> & | 中心座標を返す |
| getBounds     | const std::array<std::array<double,2>,3> & | バウンディングボックスの範囲を返す |
| Distance      | Tdd       | 指定座標との最小・最大距離を計算する |
| getVolume     | double    | バウンディングボックスの体積を計算する |
| getScale      | double    | バウンディングボックスのスケールを計算する |
| getCenter     | std::array<double,3>      | バウンディングボックスの中心座標を計算する |
| getVertices   | std::array<std::array<double,3>,8>    | バウンディングボックスの8つの頂点を計算する |

#### 🐚 🐚 オペレーター  

| オペレーター名 | 戻り値型 | 説明 |
|--------------|----------|------|
| ()            | const std::array<std::array<double,2>,3> & | バウンディングボックスの範囲を返す（関数呼び出しオペレータ） |

---

#### 🐚 🐚 有用性  
CoordinateBounds クラスは、3次元空間での領域制限やクエリ処理、衝突判定などに使用できる。簡易的な操作で座標の範囲や距離、体積などを計算できるため、効率的な空間分割やデータ処理が可能である。
[../../include/basic_geometry.hpp#L464](../../include/basic_geometry.hpp#L464)


## ⛵ ⛵ Bucket クラスの説明  

このクラスは，オブジェクトを３次元空間内に配置し，効率的に検索できるようにするための「バケツ（Bucket）」構造を提供します．

⚠️ テンプレート型`T`のオブジェクトは，`getX()`でxyz座標を取得できる必要があります．

### 🪼 🪼 型エイリアス  

| エイリアス   | 説明                                         |
|:------------:|:--------------------------------------------:|
| `sizeType`   | サイズ型（int）                              |
| `ST`         | sizeTypeの別名                               |
| `ST2`        | 2次元サイズ配列（std::array<sizeType, 2>）   |
| `ST3`        | 3次元サイズ配列（std::array<sizeType, 3>）   |
| `ST6`        | 6次元サイズ配列（std::array<sizeType, 6>）   |

### 🪼 🪼 メンバ変数  

| 変数名                | 説明                                              |
|:---------------------:|:-------------------------------------------------:|
| `xbounds`, `ybounds`, `zbounds` | X,Y,Z 座標の境界値．`Tdd` 型                   |
| `xsize`, `ysize`, `zsize` | 各座標のサイズ．`ST` 型                          |
| `bounds`              | 全体のバウンド．`T3Tdd` 型                        |
| `center`              | 空間の中心座標．`Tddd` 型                          |
| `dn`                  | 各座標のサイズ（ST3 型）                           |
| `data`             | バケツ（3D配列）                                  |
| `data_vector`      | バケツのベクター版                                 |
| `data_bool`        | バケツが空であるかどうかを示すbool値の3D配列       |
| `vector_is_set`       | ベクター版が設定されているか                       |
| `all_stored_objects`  | 保存されている全てのオブジェクト                   |
| `map_to_ijk`          | オブジェクトからインデックスへのマッピング         |
| `dL`                  | バケツの１辺の長さ                                 |

### 🪼 🪼 メソッド  

#### 🐚 🐚 初期化関連  

- `initialize(const T3Tdd &boundingboxIN, const double dL_IN)`: バケツを初期化する．

#### 🐚 🐚 インデックス変換  

- `itox(const ST i, const ST j, const ST k) const`: インデックスから座標へ変換．
- `indices(const Tddd &x) const`: 座標からインデックスへ変換．

#### 🐚 🐚 データ追加・削除  

- `add(const Tddd &x, const T p)`: オブジェクトを追加．
- `erase(T const p)`: オブジェクトを削除．

#### 🐚 🐚 その他  

- `none_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const`: 条件に合うオブジェクトがないか確認．

### 🪼 🪼 使用例  

```cpp
// 座標の境界を定義
T3Tdd bounding_box = {{{0.0, 10.0}, {0.0, 10.0}, {0.0, 10.0}}};

// バケツの一辺の長さを定義
double bucket_edge_length = 1.0;

// Bucket インスタンスを初期化
BaseBuckets<MyObject> my_buckets(bounding_box, bucket_edge_length);

// オブジェクトの座標と初期値を定義
Tddd obj1_coordinates = {5.0, 5.0, 5.0};
MyObject obj1 = MyObject(some initializers);

Tddd obj2_coordinates = {6.0, 6.0, 6.0};
MyObject obj2 = MyObject(some initializers);

// オブジェクトを追加
my_buckets.add(obj1_coordinates, obj1);
my_buckets.add(obj2_coordinates, obj2);

// オブジェクトを削除
my_buckets.erase(obj1);
```
[../../include/lib_spatial_partitioning.hpp#L4](../../include/lib_spatial_partitioning.hpp#L4)


---

### 🪼 例 

この例では，バニーのモデルを空間分割し，各バケットに含まれる点群と面を出力している．
面は広がりがあるので，複数のバケットに含まれることがある．
面と交わるバケットに確実に麺を保存するスマートな方法は，現在のところ思いつかない．
今の所は，面を細かい点群に分けて，各点群を含むバケットに面を保存することで対応している．

![anim.gif](anim.gif)

[./example1_space_partitioning.cpp#L6](./example1_space_partitioning.cpp#L6)

# 🐋 木構造による空間分割 

シンプルな空間分割クラスを拡張し，木構造による空間分割を試みる．

[./example2_tree.cpp#L6](./example2_tree.cpp#L6)

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
