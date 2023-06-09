# Contents

- [🐋連立一次方程式の解法](#🐋連立一次方程式の解法)
    - [⛵️一般化最小残差法(GMRES)](#⛵️一般化最小残差法(GMRES))
    - [⛵️⛵️ArnoldiProcess](#⛵️⛵️ArnoldiProcess)
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

## ⛵️⛵️ArnoldiProcess  

ヘッセンベルグ行列$`H[0:k-1]`$は，Aと相似なベクトルであり，同じ固有値を持つ
GMRESで使う場合，$`V0`$にはNormalize(b-A.x0)を与える．
x0は初期値

アーノルディ法は固有値問題の数値解法であり反復解法．
一般的な行列の固有ベクトルと固有値をクリロフ空間の直行基底によって近似する方法計算する方法．

1. 正規化した$`{\bf v} _0`$を与えておく．
2. $`\quad\quad\quad\quad\quad{\bf v} _1 = {\rm Normalize}(A{\bf v} _0 - ((A{\bf v} _0) \cdot {\bf v} _0){\bf v} _0)`$を計算する．
3. $`\quad\quad\quad{\bf v} _2 = {\rm Normalize}((w=A{\bf v} _1 - ((A{\bf v} _1) \cdot {\bf v} _0){\bf v} _0)) - (w \cdot {\bf v} _1){\bf v} _1)`$を計算する．
4. $`{\bf v} _3 = {\rm Normalize}((w=((w=A{\bf v} _2 - ((A{\bf v} _2) \cdot {\bf v} _0){\bf v} _0)) - (w \cdot {\bf v} _1){\bf v} _1)) - (w \cdot {\bf v} _2){\bf v} _2)`$を計算する．

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

[ArnoldiProcessの行列-ベクトル積](../../include/basic_linear_systems.hpp#L798)は特に計算コストが高い．
[CSRのDot積を並列化](../../include/basic_linear_systems.hpp#L674)すれば，かなり高速化できる．


[./test2_CSR.cpp#L1](./test2_CSR.cpp#L1)


---
