# Contents
- [🐋 `Network`](#🐋-`Network`)
    - [⛵ 3Dファイルの読み込みと出力](#⛵-3Dファイルの読み込みと出力)
        - [🪼 読み込み `Network`](#🪼-読み込み-`Network`)
        - [🪼 出力 `vtkPolygonWrite`](#🪼-出力-`vtkPolygonWrite`)
            - [🐚 面の出力](#🐚-面の出力)
            - [🐚 線の出力](#🐚-線の出力)
            - [🐚 実行方法](#🐚-実行方法)
        - [🪼 `PVDWriter`を使ったpvdファイルの作成方法](#🪼-`PVDWriter`を使ったpvdファイルの作成方法)
- [🐋 空間分割（space_partitioning）](#🐋-空間分割（space_partitioning）)
    - [⛵ 等間隔のシンプルな空間分割](#⛵-等間隔のシンプルな空間分割)
        - [🪼 例](#🪼-例)
    - [⛵ ⛵ `Bucket`クラス](#⛵-⛵-`Bucket`クラス)
        - [🪼 🪼 メンバ変数](#🪼-🪼-メンバ変数)
        - [🪼 🪼 メソッド](#🪼-🪼-メソッド)
            - [🐚 🐚 初期化関連](#🐚-🐚-初期化関連)
            - [🐚 🐚 インデックス変換¸](#🐚-🐚-インデックス変換¸)
            - [🐚 🐚 データ追加・削除](#🐚-🐚-データ追加・削除)
            - [🐚 🐚 その他](#🐚-🐚-その他)
    - [⛵ 階層のある空間分割（木構造）](#⛵-階層のある空間分割（木構造）)
    - [⛵ 空間分割の応用例：オブジェクトの接触や交差の判定](#⛵-空間分割の応用例：オブジェクトの接触や交差の判定)
        - [🪼 オブジェクトの接触や交差の判定](#🪼-オブジェクトの接触や交差の判定)
        - [🪼 面同士の接触判定](#🪼-面同士の接触判定)
- [🐋 CGALを使って四面体を生成する](#🐋-CGALを使って四面体を生成する)
    - [⛵ CGALを使って四面体を生成する](#⛵-CGALを使って四面体を生成する)
    - [⛵ 四面体を生成（制約付き四面分割 constrained tetrahedralization）](#⛵-四面体を生成（制約付き四面分割-constrained-tetrahedralization）)


---
# 🐋 `Network` 

数値シミュレーションの多くは，三角形や四面体の３Dメッシュを利用して行われる．
単純に配列にメッシュ情報を格納していては，シミュレーションの質を高めるための工夫を加えることが難しいと思われる．
例えば，ある節点に隣接する節点や辺や要素を取得するのは，配列に格納しているだけでは効率的に行うことができない．
`Network`クラスは，節点や辺や面の接続関係を保持し，接続関係から相互にアクセスできるようにするためのクラスである．
また，メッシュを細分化することもできる．
また，一度作成した`Network`クラスのオブジェクトから，`vtk`ファイルや`obj`ファイルなどを出力することができる．

* データの読み込みと出力
* 節点や辺や面の相互アクセス
* メッシュの細分化

## ⛵ 3Dファイルの読み込みと出力 

### 🪼 読み込み `Network` 

[Networkのコンストラクタ](../../include/Network.hpp#L4038)では，引数として，**OFFファイル**または**OBJファイル**をあたえることができる．
`Load3DFile`クラスを使ってデータを読み込み，`Network`クラスを作成する．

```cpp
auto obj = new Network("./bunny.obj");//ファイルからNetworkオブジェクトを作成
```

### 🪼 出力 `vtkPolygonWrite` 

`Network`クラスは，`getFaces`メンバ関数を使って簡単に面の情報を取得できる．

`vtkPolygonWrite`を使うと，`Network`クラスの面の情報を，`vtp`ファイルとして出力できる．
`vtkPolygonWrite`には，`ofstream`と，`std::vector<networkFace*>`や`std::vector<networkLine*>`などを渡し，出力できる．

#### 🐚 面の出力 

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getFaces());
```

![sample.png](sample.png)

#### 🐚 線の出力 

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getEdges());
```

#### 🐚 実行方法 

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_load_3d_file.cpp
make
./example0_load_3d_file
```

[./example0_load_3d_file.cpp#L1](./example0_load_3d_file.cpp#L1)

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
| ![sample.gif](sample.gif) | ![sample_line.gif](sample_line.gif) |

💡 QuickTimeで作成したmovファイルをgifに変換するには，次のようにする．

```sh
ffmpeg -i line.mov -filter_complex "[0:v] fps=30, scale=iw*0.5:ih*0.5 [v]" -map "[v]" sample_line.gif
```

[./example0_load_3d_file.cpp#L85](./example0_load_3d_file.cpp#L85)

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


## ⛵ ⛵ `Bucket`クラス  

`Bucket`クラスは，オブジェクトを３次元空間内に配置し，効率的に検索できるようにするための「バケツ（Bucket）」構造を提供します．

⚠️ テンプレート型`T`のオブジェクトは，予め`getX()`を使ってxyz座標を取得できるようにしておく必要がある．

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
œa
#### 🐚 🐚 インデックス変換¸  

- `itox(const ST i, const ST j, const ST k) const`: インデックスから座標へ変換．
- `indices(const Tddd &x) const`: 座標からインデックスへ変換．

#### 🐚 🐚 データ追加・削除  

- `add(const Tddd &x, const T p)`: オブジェクトを追加．
- `erase(T const p)`: オブジェクトを削除．

#### 🐚 🐚 その他  

`apply(const Tddd &x, const double d, const std::function<bool(const T &)> &func)`は，バケツの範囲を指定して，その範囲内のオブジェクトに対して関数を適用する．
これと似た関数として，

* `any_of`
* `all_of`
* `none_of`

があり，それぞれ，バケツの範囲内のオブジェクトに対して，関数を適用し，その結果が，それぞれ，`true`，`false`，`false`であれば，`true`を返す．

📝 これらの関数は，`apply`はある点を中心として半径`d`の球状の範囲を指定することができる．これは球状の範囲を指定していることになる．このような範囲指定以外に，直線上の範囲指定や，平面上の範囲指定などもできるようにしたい．
そのためには，バケツのセルと，線分や平面の交差判定を高速に行う関数が必要になる．
ラフに行っても問題ない．
線に関しては細かい分割によってインデックス変換できる．
平面に関しては，平面の方程式を使って，バケツのセルとの交差判定を行う．
[../../include/lib_spatial_partitioning.hpp#L4](../../include/lib_spatial_partitioning.hpp#L4)

[./example1_space_partitioning.cpp#L6](./example1_space_partitioning.cpp#L6)

---
## ⛵ 階層のある空間分割（木構造） 

シンプルな空間分割クラスを拡張し，木構造による空間分割を試みる．

`has_tree`が`true`の場合，`buckets`には`Bucket`クラスのポインタが格納される．
`buckets[i][j][k]`には，上のレベルの`data[i][j][k]`のデータが引き継がれている．
つまり，`buckets[i][j][k]`は，`data[i][j][k]`のデータをさらに分割したものである．
デフォルトでは，`buckets[i][j][k]`は内部に８つの`data`を持つ:

`data[0][0][0]`，`data[0][0][1]`，`data[0][1][0]`，`data[0][1][1]`，`data[1][0][0]`，`data[1][0][1]`，`data[1][1][0]`，`data[1][1][1]`．

[このツリー生成方法](../../include/lib_spatial_partitioning.hpp#L85)は，
バウンディングボックスを範囲と，それを分割する幅を指定する．
分割数を指定するよりも，この方法のように分割幅を指定する方が，自分はわかりやすい．

```cpp
buckets[i][j][k] = std::make_shared<Buckets<T>>(bounds, this->dL * 0.5 + 1e-10);
```

![example2_tree_faster.gif](example2_tree_faster.gif)

レベル０が生成したレベル１のバケットに保存された点を示しており，
白い線は，１階層上のレベル０のバケットの境界を示している．

[./example2_tree.cpp#L2](./example2_tree.cpp#L2)

---
## ⛵ 空間分割の応用例：オブジェクトの接触や交差の判定 

### 🪼 オブジェクトの接触や交差の判定 

`Network`クラスは，`makeBucketPoints`でバケツ`BucketPoints`を準備し，内部に保存している点をバケツに保存する．
同様に，`makeBucketFaces`でバケツを`BucketFaces`を準備し，内部に保存している面をバケツに保存する．

要素の接触や交差の判定には，[`IntersectQ`](../../include/basic_geometry.hpp#L1606)関数を使う．
また，接触判定の高速化のために，空間分割を使う．

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example3_line_face_interaction.cpp
make
./example3_line_face_interaction
```

![./example3/anim.gif](example3/anim_faster.gif)

[./example3_line_face_interaction.cpp#L4](./example3_line_face_interaction.cpp#L4)

---
### 🪼 面同士の接触判定 

[`IntersectQ`](../../include/basic_geometry.hpp#L1606)関数は，交差判定には使えるが，接触判定には使えない．

接触は，ギリギリ交差している状態を指すだろうが，
実際に接触判定を応用する場面では，
交差していなくとも接触していると判定させたい場合が多いだろう．
なので，接触判定条件はより緩く設定されることが多い．

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example3_line_face_interaction.cpp
make
./example3_line_face_interaction
```

![./example3/anim.gif](example3/anim_faster.gif)

[./example4_face2face_contact.cpp#L4](./example4_face2face_contact.cpp#L4)

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
