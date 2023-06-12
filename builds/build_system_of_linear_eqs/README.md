# Contents

- [🐋連立一次方程式の解法](#🐋連立一次方程式の解法)
    - [⛵️⛵️Arnoldi過程](#⛵️⛵️Arnoldi過程)
    - [⛵️⛵️一般化最小残差法/GMRES](#⛵️⛵️一般化最小残差法/GMRES)
        - [🪸テスト](#🪸テスト)
    - [⛵️LU分解(LAPACK)](#⛵️LU分解(LAPACK))
    - [⛵️Compressed Sparse Row (CSR)](#⛵️Compressed-Sparse-Row-(CSR))


---
# 🐋連立一次方程式の解法 

## ⛵️⛵️Arnoldi過程  

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
この成分からなる行列が，Hessenberg行列$H$である．

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

[../../include/basic_linear_systems.hpp#L762](../../include/basic_linear_systems.hpp#L762)



## ⛵️⛵️一般化最小残差法/GMRES  

残差$`\|{\bf b} - A{\bf x}\|`$を最小とするような$`{\bf x}`$を求めたい．
そのような$`{\bf x}`$を，クリロフ部分空間の正規直交基底を用いた，$`{\bf x} _n = V _n {\bf y} _n`$の形で近似解し，追い求めていく．
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
<summary>なぜ，アーノルディ分解をするのか</summary>

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

💡 展開項数$`n`$を$`n+1`$と大きくする際に，始めから計算しなおす必要はない．$`V _{n+1}`$と$`{\tilde H} _{n+1}`$は，$`V _n`$と$`{\tilde H} _n`$を使って計算できる．

[../../include/basic_linear_systems.hpp#L850](../../include/basic_linear_systems.hpp#L850)



### 🪸テスト 

<details>
<summary>HOW TO USE</summary>

![](WATCHME.gif)

</details>


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

[ArnoldiProcessの行列-ベクトル積](../../include/basic_linear_systems.hpp#L840)は特に計算コストが高い．
[CSRのDot積を並列化](../../include/basic_linear_systems.hpp#L674)すれば，かなり高速化できる．


[./test2_CSR.cpp#L1](./test2_CSR.cpp#L1)


---
