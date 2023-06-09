# Contents

- [🐋連立一次方程式の解法](#🐋連立一次方程式の解法)
    - [⛵️一般化最小残差法(GMRES)](#⛵️一般化最小残差法(GMRES))
    - [⛵️⛵️Arnoldi Process](#⛵️⛵️Arnoldi-Process)
    - [⛵️LU分解(LAPACK)](#⛵️LU分解(LAPACK))
    - [⛵️Compressed Sparse Row (CSR)](#⛵️Compressed-Sparse-Row-(CSR))


---
# 🐋連立一次方程式の解法 

## ⛵️一般化最小残差法(GMRES) 

- ヘッセンベルグ行列$`H`$
- クリロフ部分空間の直交基底$`V`$
- $`H`$をQR分解した行列$`Q`$と$`R`$
- $`g`$は行列$`Q`$の最初の列

ArnoldiProcessによって，$`H`$と$`V`$を求める．このArnoldiProcessクラスの派生クラスとしてGMRESを定義している．

<details>
<summary>HOW TO USE</summary>

![](WATCHME.gif)

</details>

## ⛵️⛵️Arnoldi Process  

ヘッセンベルグ行列$`H[0:k-1]`$は，Aと相似なベクトルであり，同じ固有値を持つ
GMRESで使う場合，$`V0`$にはNormalize(b-A.x0)を与える．
x0は初期値

アーノルディ法は固有値問題の数値解法であり反復解法．
一般的な行列の固有ベクトルと固有値をクリロフ空間の直行基底によって近似する方法計算する方法．

1. 正規化した$`{\bf v} _1`$を与えておく．
2. $`{\bf v} _2 = {\rm Normalize}(\,\,\,\quad\quad\quad\quad\quad A{\bf v} _1 - ((A{\bf v} _1) \cdot {\bf v} _1){\bf v} _1\,\,\qquad\qquad\qquad\qquad\qquad\qquad)`$を計算する．
3. $`{\bf v} _3 = {\rm Normalize}(\quad\quad\quad({\bf w}=A{\bf v} _2 - ((A{\bf v} _2) \cdot {\bf v} _1){\bf v} _1)) - ({\bf w} \cdot {\bf v} _2){\bf v} _2\quad\quad\quad\quad\quad\quad)`$を計算する．
4. $`{\bf v} _4 = {\rm Normalize}(({\bf w}=(({\bf w}=A{\bf v} _3 - ((A{\bf v} _3) \cdot {\bf v} _1){\bf v} _1)) - ({\bf w} \cdot {\bf v} _2){\bf v} _2)) - ({\bf w} \cdot {\bf v} _3){\bf v} _3)`$を計算する．

言い換えると，

1. 正規化した$`{\bf v} _1`$を与えておく．
2. $`{\bf w}=A{\bf v} _1, {\bf v} _2 = {\rm Normalize}({\rm Chop}({\bf w},{\bf v} _1))`$を計算する．
3. $`{\bf w}=A{\bf v} _2, {\bf v} _3 = {\rm Normalize}({\rm Chop}({\rm Chop}({\bf w}, {\bf v} _1), {\bf v} _2))`$を計算する．
4. $`{\bf w}=A{\bf v} _3, {\bf v} _4 = {\rm Normalize}({\rm Chop}({\rm Chop}({\rm Chop}({\bf w}, {\bf v} _1), {\bf v} _2), {\bf v} _3))`$を計算する．

$`A{\bf v} _i`$の直交化の際に，
それに含まれる各基底$`{\bf v} _0,{\bf v} _1,...,{\bf v} _i`$の成分を計算している．
この成分からなる行列が，Hessenberg行列$`H`$である．

$$
\begin{align*}
A{\bf v} _1 & = h _{1,1} {\bf v} _1 + h _{2,1} {\bf v} _2\\
A{\bf v} _2 & = h _{1,2} {\bf v} _1 + h _{2,2} {\bf v} _2 + h _{3,2} {\bf v} _3\\
& \dots\\
A{\bf v} _{n} & = h _{1,n} {\bf v} _1 + h _{2,n} {\bf v} _2 + \cdots + h _{n,n+1} {\bf v} _{n+1}
\end{align*}
$$

行列を使ってまとめると，

```math
A V _n = V _{n+1} \tilde H _n, \quad V _n = [v _1|v _2|...|v _n],
\quad \tilde H _n = \begin{bmatrix} h _{1,1} & h _{1,2} & \cdots & h _{1,n} & h _{1,n+1} \\ h _{2,1} & h _{2,2} & \cdots & h _{2,n} & h _{2,n+1} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & h _{n,n} & h _{n,n+1} \\ 0 & 0 & \cdots & 0 & h _{n+1,n+1} \end{bmatrix}
```

これをArnoldi分解という．ここで，$`[v _1|v _2|...|v _n]`$
の$`|`$は列ベクトルを連結して行列を形成することを示している．

```math
A V _n {\bf y} _n = V _{n+1} \tilde H _n {\bf y} _n
```

```math
{\bf b} - A V _n {\bf y} _n = {\bf b} - V _{n+1} \tilde H _n {\bf y} _n
```

答えとして，$`\|{\bf b} - V _{n+1} \tilde H _n {\bf y} _n\|`$を最小とするような，
$`{\bf y} _n`$を求める．

```math
{\bf x} _n \approx V _n {\bf y} _n
```

[../../include/basic_linear_systems.hpp#L762](../../include/basic_linear_systems.hpp#L762)


[./test0_GMRES.cpp#L1](./test0_GMRES.cpp#L1)


## ⛵️LU分解(LAPACK)


[./test0_LAPACK.cpp#L1](./test0_LAPACK.cpp#L1)


EigenのGMRESを使った結果と比較．


[./test1_EIGEN_GMRES.cpp#L6](./test1_EIGEN_GMRES.cpp#L6)


---
## ⛵️Compressed Sparse Row (CSR) 

CSRは行列を表現する方法の一つである．
このCSRクラスは，std::unordered_mapを用いて，行列の非ゼロ要素を表現する．
std::unordered_mapのkeyはポインタであり，valueはdoubleである．
CSRクラス自身が，行列の行番号を保存しており，keyであるCSRクラスは行列の列番号を保存している．

[ArnoldiProcessの行列-ベクトル積](../../include/basic_linear_systems.hpp#L844)は特に計算コストが高い．
[CSRのDot積を並列化](../../include/basic_linear_systems.hpp#L674)すれば，かなり高速化できる．


[./test2_CSR.cpp#L1](./test2_CSR.cpp#L1)


---
