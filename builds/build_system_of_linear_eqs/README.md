# Contents
- [🐋 🐋 圧縮行格納法 (Compressed Row Storage, CRS)](#🐋-🐋-圧縮行格納法-(Compressed-Row-Storage,-CRS))
    - [⛵ ⛵ 実装方法](#⛵-⛵-実装方法)
    - [⛵ ⛵ CRS構造体の仕様](#⛵-⛵-CRS構造体の仕様)
        - [🪼 🪼 概要](#🪼-🪼-概要)
        - [🪼 🪼 メンバ変数](#🪼-🪼-メンバ変数)
        - [🪼 🪼 メンバ関数](#🪼-🪼-メンバ関数)
        - [🪼 CRSの使用例](#🪼-CRSの使用例)
            - [🪸 CRSは，ある行ベクトルを格納するクラスと考える](#🪸-CRSは，ある行ベクトルを格納するクラスと考える)
            - [🪸 `setIndexCRS`](#🪸-`setIndexCRS`)
            - [🪸 値を格納：`set`と`increment`](#🪸-値を格納：`set`と`increment`)
            - [🪸 `selfDot`](#🪸-`selfDot`)
- [🐋 連立一次方程式の解法](#🐋-連立一次方程式の解法)
    - [⛵ ⛵ Arnoldi過程](#⛵-⛵-Arnoldi過程)
        - [🪼 🪼 基底ベクトルの追加](#🪼-🪼-基底ベクトルの追加)
    - [⛵ ⛵ 一般化最小残差法 (Generalized Minimal Residual Method, GMRES)](#⛵-⛵-一般化最小残差法-(Generalized-Minimal-Residual-Method,-GMRES))
        - [🪼 🪼 GMRESの計算複雑性](#🪼-🪼-GMRESの計算複雑性)
        - [🪼 テスト](#🪼-テスト)
    - [⛵ LU分解(LAPACK)](#⛵-LU分解(LAPACK))
    - [⛵ 共役勾配法と勾配降下法](#⛵-共役勾配法と勾配降下法)
        - [🪼 共役勾配法（Conjugate Gradient, CG）](#🪼-共役勾配法（Conjugate-Gradient,-CG）)
        - [🪼 勾配降下法 (Gradient Descent, GD)](#🪼-勾配降下法-(Gradient-Descent,-GD))


---
# 🐋 🐋 圧縮行格納法 (Compressed Row Storage, CRS)  

CRSは疎行列を表現する一つの手法．
このクラスでは，行列-ベクトル積の高速化と，その行列とベクトルの管理とそれらへのアクセスの容易さを目的とし，
次のような考えで，CRSを実装している．

## ⛵ ⛵ 実装方法  

多くの数値計算で最も高速化したいのは，行列-ベクトル積である．
整数のインデックスを使って，行列-ベクトル積を計算すると，インデックスの管理が大変になる．
例えば，インデックスが消えたり増えたりする場合が大変だ（例えば，要素の節点番号や，粒子法の粒子番号）．

そこで，ポインタをキーとした連装配列として，行列-ベクトル積を計算することが一つの解決策である．
行ベクトルの成分がどの節点や粒子に対応しているかを，ポインタで管理するものである．

また，CRSクラスに，できるだけ`Dot(A,V)`で利用する情報を保存しておくことで，管理とアクセスの容易さを向上させたい．

多くの場合，行ベクトルはある節点や粒子に対して成り立つ方程式を表しているので，
行のインデックスは，その節点や粒子のインデックスと一致する．
強く関連するので，この行ベクトル自体を，それが成り立つ節点や粒子クラスに保存しておくのは，自然なことである．

さらに，そのCRSには方程式`A[i]`だけでなく，節点上または粒子上の`V[i]`も保存しておくことも，自然なことである．

つまり，`Dot(A,V)`の計算を考えて，
CRSには次のような情報を保存しておくことにする．

- 方程式：`A[i]`は，`CRS->column_value`に保存されている．`column_value`は，`std::unordered_map<CRS *, double>`
- 値：`V[i]`は，`CRS->value`に保存されている
- 当然CRSとしてのポインタ
- 行番号(`i`)

特に、Arnoldiプロセスにおける行列-ベクトル積の計算は計算コストが高いです（参照: ArnoldiProcessの行列-ベクトル積）。CRSのDot積を並列化することで、大幅な高速化が可能です（参照: CRSのDot積を並列化）。

## ⛵ ⛵ CRS構造体の仕様  

### 🪼 🪼 概要  

CRS（Compressed Row Storage）構造体は、疎行列の一部を効率的に格納するためのデータ構造であり、高速な線形代数の計算を実現します。

### 🪼 🪼 メンバ変数  

| 変数名 | 型 | 説明 |
|:------:|:--:|:----:|
| `column_value` | `std::unordered_map<CRS *, double>` | 行の非ゼロ要素を格納する連想配列 |
| `value` | `double` | 一般的な値を格納 |
| `diagonal_value` | `double` | 対角要素の値 |
| `tmp_value` | `double` | 一時的な値の格納用 |
| `canUseVector` | `bool` | ベクタが使用可能かのフラグ |
| `value3d` | `std::array<double, 3>` | 3次元空間の値 |
| `__index__` | `std::size_t` | インデックス |

### 🪼 🪼 メンバ関数  

| 関数名 | 引数 | 戻り値 | 説明 |
|:------:|:----:|:------:|:----:|
| `clearColumnValue` | なし | `void` | `column_value`をクリアし、`canUseVector`を`false`に設定する |
| `setIndexCRS` | `std::size_t i` | `void` | インデックス`__index__`を設定する |
| `getIndexCRS` | なし | `std::size_t` | インデックス`__index__`を取得する |
| `at` | `CRS *const p` | `double` | 指定された`p`に対応する`column_value`の値を取得する |
| `contains` | `CRS *const p` | `bool` | 指定された`p`が`column_value`に含まれているかを確認する |
| `increment` | `CRS *const p, const double v` | `void` | 指定された`p`に対する`column_value`の値に`v`を加算、または新規挿入する |
| `setVectorCRS` | なし | `void` | `column_value`を`std::vector`形式に変換し、`canUseVector`を`true`に設定する |
[../../include/basic_linear_systems.hpp#L1140](../../include/basic_linear_systems.hpp#L1140)


### 🪼 CRSの使用例 

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test3_CRS.cpp
make
./test3_CRS
```

ここには，`A`かける`V`をCRSを使って高速に計算する例を示している．

[./test3_CRS.cpp#L1](./test3_CRS.cpp#L1)

---
#### 🪸 CRSは，ある行ベクトルを格納するクラスと考える 

私のプログラムでは，Row-major orderで行列を格納しており，次のように行列を定義している．

```cpp
std::vector<std::vector<double>> Mat; // <- std::vector<ROW VECTOR> Mat;
```

```math
\begin{pmatrix}
\{a _{11} & a _{12} & a _{13} & \cdots & a _{1n}\}&\leftarrow {\text{a ROW VECTOR}} \\
\{a _{21} & a _{22} & a _{23} & \cdots & a _{2n}\}&\\
\{a _{31} & a _{32} & a _{33} & \cdots & a _{3n}\}&\\
\{\vdots & \vdots & \vdots & \ddots & \vdots \}&\\
\{a _{m1} & a _{m2} & a _{m3} & \cdots & a _{mn}\}&
\end{pmatrix}
```

CRSは，このROW VECTORを格納するクラスであり，CRSのベクトルが行列となる．

```cpp
std::vector<CRS*> Mat_CRS(A.size());
```

[./test3_CRS.cpp#L135](./test3_CRS.cpp#L135)

---
#### 🪸 `setIndexCRS` 

CRSは，`CRS->setIndexCRS(i)`のようにして，自身の行番号を保持しておく．
このインデックスは，`std::vector<VRS*>`と掛け算をする相手である`V`の行番号に対等させておく必要がある．

**掛け算`Dot(A,V)`において，CRS（これは行ベクトルと考える）は，自分に保存されている{row index,value}のセットを元に，`V[row index]*value`のようにして足し合わせていく．**

[./test3_CRS.cpp#L162](./test3_CRS.cpp#L162)

---
#### 🪸 値を格納：`set`と`increment` 

値を格納するには，２つの方法があり，`set`と`increment`がある．
このように，インデックスと値を指定するようにしているのは，
値がゼロの場合は，何もせず，インデックスも保存しないようにしているためである．

値を設定する，`set`と`increment`の第一引数は，CRSのポインタである．

[./test3_CRS.cpp#L175](./test3_CRS.cpp#L175)

#### 🪸 `selfDot` 

`selfDot`は，CRSに保存した`A`と`V`を掛け合わせる関数である．

[./test3_CRS.cpp#L197](./test3_CRS.cpp#L197)

---
# 🐋 連立一次方程式の解法 

## ⛵ ⛵ Arnoldi過程  

アーノルディ分解は，Krylov部分空間の生成のために使われる．
GMRESは，Krylov部分空間と呼ばれる線形空間内で反復解を探すアルゴリズム．
この部分空間は，行列Aと初期ベクトルから生成される．
GMRESにとってアーノルディ分解は，ヘッセンベルグ行列を生成するための手段．
得られたヘッセンベルグ行列はQR分解され，直交行列と上三角行列に分解される．

1. 正規化した$`{\bf v} _1`$を与えておく．
2. $`{\bf v} _2 = {\rm Normalize}(\,\,\,\quad\quad\quad\quad\quad A{\bf v} _1 - ((A{\bf v} _1) \cdot {\bf v} _1){\bf v} _1\,\,\qquad\qquad\qquad\qquad\qquad\qquad)`$を計算する．
3. $`{\bf v} _3 = {\rm Normalize}(\quad\quad\quad({\bf w}=A{\bf v} _2 - ((A{\bf v} _2) \cdot {\bf v} _1){\bf v} _1)) - ({\bf w} \cdot {\bf v} _2){\bf v} _2\quad\quad\quad\quad\quad\quad)`$を計算する．
4. $`{\bf v} _4 = {\rm Normalize}(({\bf w}=(({\bf w}=A{\bf v} _3 - ((A{\bf v} _3) \cdot {\bf v} _1){\bf v} _1)) - ({\bf w} \cdot {\bf v} _2){\bf v} _2)) - ({\bf w} \cdot {\bf v} _3){\bf v} _3)`$を計算する．

言い換えると，

1. 正規化した$`{\bf v} _1`$を与えておく．
2. $`{\bf w}=A{\bf v} _1, {\bf v} _2 = {\rm Normalize}({\rm Chop}({\bf w},{\bf v} _1))`$を計算する．
3. $`{\bf w}=A{\bf v} _2, {\bf v} _3 = {\rm Normalize}({\rm Chop}({\rm Chop}({\bf w}, {\bf v} _1), {\bf v} _2))`$を計算する．
4. $`{\bf w}=A{\bf v} _3, {\bf v} _4 = {\rm Normalize}({\rm Chop}({\rm Chop}({\rm Chop}({\bf w}, {\bf v} _1), {\bf v} _2), {\bf v} _3))`$を計算する．

これは，既存のベクトルの成分を削り落として，新しいベクトルを作っているに過ぎない．

💡 ここで最も計算コストがかかるのは，$`{\bf w}=A{\bf v} _i`$の行列-ベクトル積である．

$`A{\bf v} _i`$の直交化の際に，
それに含まれる各基底$`{\bf v} _0,{\bf v} _1,...,{\bf v} _i`$の成分を計算している．
この成分からなる行列が，Hessenberg行列$`H`$である（ほとんど上三角行列のことを指す）．

```math
\begin{align*}
A{\bf v} _1 & = h _{1,1} {\bf v} _1 + h _{2,1} {\bf v} _2\\
A{\bf v} _2 & = h _{1,2} {\bf v} _1 + h _{2,2} {\bf v} _2 + h _{3,2} {\bf v} _3\\
& \dots\\
A{\bf v} _{n} & = h _{1,n} {\bf v} _1 + h _{2,n} {\bf v} _2 + \cdots + h _{n,n+1} {\bf v} _{n+1}
\end{align*}
```

💡 ここで，行数よりも項数が1多いことに注目しよう．

行列を使ってまとめると，

```math
A V _n = V _{n+1} \tilde H _n, \quad V _n = [v _1|v _2|...|v _n],
\quad \tilde H _n = \begin{bmatrix} h _{1,1} & h _{1,2} & \cdots & h _{1,n} & h _{1,n+1} \\ h _{2,1} & h _{2,2} & \cdots & h _{2,n} & h _{2,n+1} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & h _{n,n} & h _{n,n+1} \\ 0 & 0 & \cdots & 0 & h _{n+1,n+1} \end{bmatrix}
```

💡 $`\tilde H _n`$は，Hessenberg行列が１行長くなった行列になっている．これは，前の式において，行数よりも項数が1多いことによる．

これをArnoldi分解という．ここで，$`[v _1|v _2|...|v _n]`$の$`|`$は列ベクトルを連結して行列を形成することを示している．

### 🪼 🪼 基底ベクトルの追加  

基底ベクトルを追加したい場合にどのような操作が必要となるか整理しておこう．
これは，GMRES法の繰り返し計算の中で必要となる．
[../../include/basic_linear_systems.hpp#L1468](../../include/basic_linear_systems.hpp#L1468)


## ⛵ ⛵ 一般化最小残差法 (Generalized Minimal Residual Method, GMRES)  

残差$`\|{\bf b} - A{\bf x}\|`$を最小とするような$`{\bf x}`$を求めたい．
そのような$`{\bf x}`$を，クリロフ部分空間の正規直交基底を用いた，$`{\bf x} _n = V _n {\bf y} _n`$の形で近似し，追い求めていく．
$`n`$はこの表現での展開項数である．$`V _n = \{{\bf v} _1,{\bf v} _2,...,{\bf v} _n\}`$は，アーノルディ過程によって計算する，クリロフ部分空間の正規直交基底である．

1. クリロフ部分空間法の考えから，$`\|{\bf b} - A V _n {\bf y} _n\|`$を最小とするような，$`{\bf y} _n`$を求める問題に書き換える．
2. $`A V _n = V _{n+1} \tilde H _n`$（アーノルディ分解）と書き換える．
3. $`V _{n+1}`$でくくる．
4. QR分解を使って，$`{\bf y} _n`$に関する最小二乗問題を$`{\bf y} _n`$について解く．

```math
\begin{align*}
\|{\bf b} - A{\bf x} _n\| & = \|{\bf b} - A V _n {\bf y} _n\|\\
& = \|{\bf b} - V _{n+1} \tilde H _n {\bf y} _n\|\quad \text{(use Arnoldi decomposition)}\\
& = \|V _{n+1} (\|{\bf b}\| {\bf e} _1 - \tilde H _n {\bf y} _n)\|\\
& = \|\|{\bf b}\| {\bf e} _1 - \tilde H _n {\bf y} _n\|\quad \text{(the dimension has been reduced!)}\\
& = \|\|{\bf b}\| {\bf e} _1 - QR {\bf y} _n\|\quad \text{(use QR decomposition)}\\
\end{align*}
```

ただし，これは理論の話であって，
簡単には，アーノルディ分解で得られたヘッセンベルグ行列をQR分解して，$`{\bf y} _n`$を求める，ということである．

💡 アーノルディ過程が逐次的に計算できるため，展開項数$`n`$を$`n+1`$へと大きくしようとする際に（精度が$`n`$では十分でない場合），GMRESで近似解$`{\bf x} _{n+1}`$を始めから計算しなおす必要はない．$`V _{n+1}`$と$`{\tilde H} _{n+1}`$は，$`V _n`$と$`{\tilde H} _n`$を再利用するようにして計算でき，従って，比較的安く，得られている$`{\bf x} _n`$から$`{\bf x} _{n+1}`$へと更新できる．

### 🪼 🪼 GMRESの計算複雑性  

単純に考えて，$`A`$が$`m \times m`$行列であるとすると，
GMRESの行列ベクトル積の計算量は，$`O(m^2)`$．
もし，クリロフ部分空間の基底を$`n`$個まで展開するとすると，$`O(m^2 n)`$の計算量が必要となる．
さらに，収束するまでの反復回数を$`k`$とすると，$`O(m^2 n k)`$の計算量が必要となる．

LU分解の場合は，$`O(m^3)`$の計算量が必要となる．
従って，$`m`$が大きい場合は，GMRESの方が計算量が少なくて済む．

GMRESと多重極展開法（もし$`m`$が$`m/d`$になったとすると）を組み合わせれば，GMRESは$`O(knm^2/d^2)`$で計算できる．
[../../include/basic_linear_systems.hpp#L1593](../../include/basic_linear_systems.hpp#L1593)


* GMRESは反復的な方法で，特に大規模で疎な非対称行列の線形システムを解くのに適している．
* GMRESは一般的に共役勾配法よりも柔軟性があり，非対称行列に対しても使用できる．ただし，反復の回数が増えると計算コストが大きくなる可能性がある．

### 🪼 テスト 

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test0_GMRES.cpp
make
./test0_GMRES
```

<details>
<summary>HOW TO USE</summary>

![](WATCHME.gif)

</details>

[./test0_GMRES.cpp#L1](./test0_GMRES.cpp#L1)

## ⛵ LU分解(LAPACK) 

* LU分解は直接的な方法で，あらゆる種類の行列（対称、非対称、正定値、非正定値）に適用できる．
* この方法は反復的な方法よりも計算コストが高くなる可能性があるが，反復法とは異なり，収束性の問題がない．

[./test0_LAPACK.cpp#L1](./test0_LAPACK.cpp#L1)

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test1_EIGEN_GMRES.cpp
make
./test1_EIGEN_GMRES
```

EigenのGMRESを使った結果と比較．

[./test1_EIGEN_GMRES.cpp#L2](./test1_EIGEN_GMRES.cpp#L2)

## ⛵ 共役勾配法と勾配降下法 

### 🪼 共役勾配法（Conjugate Gradient, CG） 

* 共役勾配法は反復的な方法で，特に大規模で疎な（つまり，ほとんどの要素がゼロである）対称正定値行列の線形システムを解くのに適している．
* この方法の利点は，一般的に反復回数が行列の次元に対して比較的少ないこと．しかし，非対称または非正定値の行列には適用できない．

### 🪼 勾配降下法 (Gradient Descent, GD) 

* 勾配降下法は最も基本的な最適化アルゴリズムで，線形システムまたは一般的な最適化問題を解くことができる．
* しかし，勾配降下法の収束速度は通常比較的遅く，特に凸でない問題に対しては局所最小値に陥る可能性がある．

[./test3_GradientMethod.cpp#L1](./test3_GradientMethod.cpp#L1)

---
