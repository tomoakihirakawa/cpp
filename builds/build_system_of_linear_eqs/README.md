# Contents
- [🐋 連立一次方程式の解法](#🐋-連立一次方程式の解法)
    - [⛵ ⛵ Arnoldi過程](#⛵-⛵-Arnoldi過程)
    - [⛵ ⛵ 一般化最小残差法/GMRES](#⛵-⛵-一般化最小残差法/GMRES)
        - [🪼 テスト](#🪼-テスト)
    - [⛵ LU分解(LAPACK)](#⛵-LU分解(LAPACK))
    - [⛵ 共役勾配法と勾配降下法](#⛵-共役勾配法と勾配降下法)
        - [🪼 共役勾配法（Conjugate Gradient, CG）](#🪼-共役勾配法（Conjugate-Gradient,-CG）)
        - [🪼 勾配降下法 (Gradient Descent, GD)](#🪼-勾配降下法-(Gradient-Descent,-GD))
    - [⛵ ⛵ Compresed Row Storage (CRS)](#⛵-⛵-Compresed-Row-Storage-(CRS))
- [🐋 🐋 CRS構造体のドキュメント](#🐋-🐋-CRS構造体のドキュメント)
    - [⛵ ⛵ 概要](#⛵-⛵-概要)
    - [⛵ ⛵ メンバ変数](#⛵-⛵-メンバ変数)
    - [⛵ ⛵ メンバ関数](#⛵-⛵-メンバ関数)
    - [⛵ ⛵ テンプレート関数](#⛵-⛵-テンプレート関数)
    - [⛵ ⛵ 使用例](#⛵-⛵-使用例)
    - [⛵ ⛵ 注意](#⛵-⛵-注意)


---
# 🐋 連立一次方程式の解法 

## ⛵ ⛵ Arnoldi過程  

1. 正規化した$`{\bf v} _1`$を与えておく．
2. $`{\bf v} _2 = {\rm Normalize}(\,\,\,\quad\quad\quad\quad\quad A{\bf v} _1 - ((A{\bf v} _1) \cdot {\bf v} _1){\bf v} _1\,\,\qquad\qquad\qquad\qquad\qquad\qquad)`$を計算する．
3. $`{\bf v} _3 = {\rm Normalize}(\quad\quad\quad({\bf w}=A{\bf v} _2 - ((A{\bf v} _2) \cdot {\bf v} _1){\bf v} _1)) - ({\bf w} \cdot {\bf v} _2){\bf v} _2\quad\quad\quad\quad\quad\quad)`$を計算する．
4. $`{\bf v} _4 = {\rm Normalize}(({\bf w}=(({\bf w}=A{\bf v} _3 - ((A{\bf v} _3) \cdot {\bf v} _1){\bf v} _1)) - ({\bf w} \cdot {\bf v} _2){\bf v} _2)) - ({\bf w} \cdot {\bf v} _3){\bf v} _3)`$を計算する．

言い換えると，

1. 正規化した$`{\bf v} _1`$を与えておく．
2. $`{\bf w}=A{\bf v} _1, {\bf v} _2 = {\rm Normalize}({\rm Chop}({\bf w},{\bf v} _1))`$を計算する．
3. $`{\bf w}=A{\bf v} _2, {\bf v} _3 = {\rm Normalize}({\rm Chop}({\rm Chop}({\bf w}, {\bf v} _1), {\bf v} _2))`$を計算する．
4. $`{\bf w}=A{\bf v} _3, {\bf v} _4 = {\rm Normalize}({\rm Chop}({\rm Chop}({\rm Chop}({\bf w}, {\bf v} _1), {\bf v} _2), {\bf v} _3))`$を計算する．

💡 ここで最も計算コストがかかるのは，$`{\bf w}=A{\bf v} _i`$の行列-ベクトル積である．

$`A{\bf v} _i`$の直交化の際に，
それに含まれる各基底$`{\bf v} _0,{\bf v} _1,...,{\bf v} _i`$の成分を計算している．
この成分からなる行列が，Hessenberg行列$`H`$である．

```math
\begin{align*}
A{\bf v} _1 & = h _{1,1} {\bf v} _1 + h _{2,1} {\bf v} _2\\
A{\bf v} _2 & = h _{1,2} {\bf v} _1 + h _{2,2} {\bf v} _2 + h _{3,2} {\bf v} _3\\
& \dots\\
A{\bf v} _{n} & = h _{1,n} {\bf v} _1 + h _{2,n} {\bf v} _2 + \cdots + h _{n,n+1} {\bf v} _{n+1}
\end{align*}
```

行列を使ってまとめると，

```math
A V _n = V _{n+1} \tilde H _n, \quad V _n = [v _1|v _2|...|v _n],
\quad \tilde H _n = \begin{bmatrix} h _{1,1} & h _{1,2} & \cdots & h _{1,n} & h _{1,n+1} \\ h _{2,1} & h _{2,2} & \cdots & h _{2,n} & h _{2,n+1} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & h _{n,n} & h _{n,n+1} \\ 0 & 0 & \cdots & 0 & h _{n+1,n+1} \end{bmatrix}
```

これをArnoldi分解という．ここで，$`[v _1|v _2|...|v _n]`$の$`|`$は列ベクトルを連結して行列を形成することを示している．
[../../include/basic_linear_systems.hpp#L1058](../../include/basic_linear_systems.hpp#L1058)


## ⛵ ⛵ 一般化最小残差法/GMRES  

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

<details>
<summary>なぜアーノルディ分解をするのか？</summary>

* $`A`$は$`m \times m`$とすると
* $`{\bf x}`$と$`{\bf b}`$は，$`m \times 1`$ベクトル（列ベクトル）.
* $`V _n`$は，$`m \times n`$行列で，$`A`$のクリロフ部分空間の基底ベクトルを列に持つ行列．
* $`{\bf y} _n`$は$`n \times 1`$ベクトル．
* $`\tilde H _n`$は$`(n+1) \times n`$行列．

従って，$`n`$が$`m`$よりも大幅に小さい場合，
アーノルディ分解によって作られた問題$`\min\|{\bf b} - V _{n+1}{\tilde H} _n {\bf y} _n\|`$は，
元の問題$`\min\|{\bf b}-A{\bf x}\|`$より計算量が少ない問題となる．

$`A{\bf x} = {\bf b}`$の問題を解くよりも，
$`{\tilde H} _n {\bf y} _n = {\bf b}`$という問題を解く方が計算量が少ない．

</details>

💡 アーノルディ過程が逐次的に計算できるため，展開項数$`n`$を$`n+1`$へと大きくしようとする際に（精度が$`n`$では十分でない場合），GMRESで近似解$`{\bf x} _{n+1}`$を始めから計算しなおす必要はない．$`V _{n+1}`$と$`{\tilde H} _{n+1}`$は，$`V _n`$と$`{\tilde H} _n`$を再利用するようにして計算でき，従って，比較的安く，得られている$`{\bf x} _n`$から$`{\bf x} _{n+1}`$へと更新できる．
[../../include/basic_linear_systems.hpp#L1200](../../include/basic_linear_systems.hpp#L1200)


* GMRESは反復的な方法で，特に大規模で疎な非対称行列の線形システムを解くのに適している．
* GMRESは一般的に共役勾配法よりも柔軟性があり，非対称行列に対しても使用できる．ただし，反復の回数が増えると計算コストが大きくなる可能性がある．

### 🪼 テスト 

<details>
<summary>HOW TO USE</summary>

![](WATCHME.gif)

</details>

[./test0_GMRES.cpp#L1](./test0_GMRES.cpp#L1)

## ⛵ LU分解(LAPACK) 

* LU分解は直接的な方法で，あらゆる種類の行列（対称、非対称、正定値、非正定値）に適用できる．
* この方法は反復的な方法よりも計算コストが高くなる可能性があるが，反復法とは異なり，収束性の問題がない．

[./test0_LAPACK.cpp#L1](./test0_LAPACK.cpp#L1)

EigenのGMRESを使った結果と比較．

[./test1_EIGEN_GMRES.cpp#L6](./test1_EIGEN_GMRES.cpp#L6)

## ⛵ 共役勾配法と勾配降下法 

### 🪼 共役勾配法（Conjugate Gradient, CG） 

* 共役勾配法は反復的な方法で，特に大規模で疎な（つまり，ほとんどの要素がゼロである）対称正定値行列の線形システムを解くのに適している．
* この方法の利点は，一般的に反復回数が行列の次元に対して比較的少ないこと．しかし，非対称または非正定値の行列には適用できない．

### 🪼 勾配降下法 (Gradient Descent, GD) 

* 勾配降下法は最も基本的な最適化アルゴリズムで，線形システムまたは一般的な最適化問題を解くことができる．
* しかし，勾配降下法の収束速度は通常比較的遅く，特に凸でない問題に対しては局所最小値に陥る可能性がある．

[./test3_GradientMethod.cpp#L1](./test3_GradientMethod.cpp#L1)

---
## ⛵ ⛵ Compresed Row Storage (CRS)  

CRSは行列を表現する方法の一つである．
このCRSクラスは，std::unordered_mapを用いて，行列の非ゼロ要素を表現する．
std::unordered_mapのkeyはポインタであり，valueはdoubleである．
CRSクラス自身が，行列の行番号を保存しており，keyであるCRSクラスは行列の列番号を保存している．

[ArnoldiProcessの行列-ベクトル積](../../include/basic_linear_systems.hpp#L1188)は特に計算コストが高い．
[CRSのDot積を並列化](../../include/basic_linear_systems.hpp#L970)すれば，かなり高速化できる．

# 🐋 🐋 CRS構造体のドキュメント  

## ⛵ ⛵ 概要  

CRS（Compressed Row Storage）構造体は、疎行列の一部を効率的に格納するためのデータ構造です。この構造体は疎行列を処理する際に特に有用で、高速な線形代数の計算を可能にします。

## ⛵ ⛵ メンバ変数  

| 変数名 | 型 | 説明 |
|:------:|:--:|:----:|
| `column_value` | `std::unordered_map<CRS *, double>` | 行要素を格納するための連想配列 |
| `value` | `double` | 値を格納する変数 |
| `diagonal_value` | `double` | 対角要素の値 |
| `tmp_value` | `double` | 一時的な値を格納 |
| `canUseVector` | `bool` | ベクタが使用可能かどうか |
| `value3d` | `std::array<double, 3>` | 3次元空間での値 |
| `__index__` | `size_t` | インデックス値 |

## ⛵ ⛵ メンバ関数  

| 関数名 | 引数 | 戻り値 | 説明 |
|:------:|:----:|:------:|:----:|
| `clearColumnValue` | なし | `void` | `column_value`をクリアし、`canUseVector`を`false`に設定 |
| `setIndexCRS` | `size_t i` | `void` | `__index__`を設定 |
| `getIndexCRS` | なし | `size_t` | `__index__`を取得 |
| `at` | `CRS *const p` | `double` | 指定した`p`に対応する`column_value`を取得 |
| `contains` | `CRS *const p` | `bool` | 指定した`p`が`column_value`に含まれているか確認 |
| `increment` | `CRS *const p, const double v` | `void` | `column_value`に値`v`を加算、または新規挿入 |
| `setVectorCRS` | なし | `void` | `column_value`を`std::vector`に変換し、`canUseVector`を`true`に設定 |

## ⛵ ⛵ テンプレート関数  

- 以下のテンプレート関数はDot積を計算します。CRS構造体から派生した型に対して使用できます。

```cpp
double Dot(const std::unordered_map<T *, double> &column_value, const V_d &V);
double Dot(const std::unordered_map<T *, double> &column_value, const std::unordered_set<T *> &V_crs);
V_d Dot(const std::unordered_set<T *> &V_crs);
V_d Dot(const std::unordered_set<T *> &A, const V_d &V);
V_d Dot(const std::vector<T *> &A, const V_d &V);
void DotOutput(const std::unordered_set<T *> &A, const V_d &V, V_d &w);
void DotOutput(const std::vector<T *> &A, const V_d &V, V_d &w);
```

## ⛵ ⛵ 使用例  

```cpp
// CRSオブジェクトの初期化
CRS crs;

// インデックスの設定
crs.setIndexCRS(1);

// 値の追加
crs.increment(other_crs, 5.0);

// ベクタ形式への変換
crs.setVectorCRS();

// Dot積の計算
double result = Dot(crs_map, some_vector);
```

## ⛵ ⛵ 注意  

- `#pragma omp parallel`と`#pragma omp single nowait`はOpenMPを使用して並列処理を行っています。適切なスレッド数を設定して使用してください。

このドキュメントはCRS構造体の詳細と使用方法について説明しています。適切に使用することで、高速な線形代数の計算が可能です。
[../../include/basic_linear_systems.hpp#L748](../../include/basic_linear_systems.hpp#L748)

[./test2_CRS.cpp#L1](./test2_CRS.cpp#L1)

---
