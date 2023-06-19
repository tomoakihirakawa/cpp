# Contents

- [🐋補間](#🐋補間)
    - [⛵️三角形補間](#⛵️三角形補間)
    - [⛵️⛵️範囲を修正した三角形形状関数](#⛵️⛵️範囲を修正した三角形形状関数)
    - [⛵️ラグランジュ補間](#⛵️ラグランジュ補間)
    - [⛵️三角形補間](#⛵️三角形補間)
    - [⛵️⛵️範囲を修正した三角形形状関数](#⛵️⛵️範囲を修正した三角形形状関数)


---
# 🐋補間 

## ⛵️三角形補間 

## ⛵️⛵️範囲を修正した三角形形状関数  

普通の三角形形状関数は，$`{\mathbf N}=(N _0,N _1,N _2) = (t _0,t _1,1-t _0-t _1)`$．
これを使った，$`{\rm Dot}({\mathbf N},\{{\mathbf X _0},{\mathbf X _1},{\mathbf X _2}\})`$は，$`t _0,t _1=[0,1]`$で平行四辺形を作る．
$`t _0,t _1=[0,1]`$の範囲で，三角形を形成するように変数変換したいことがある．
そのたびに，変数変換をプログラムするのは面倒なので，予め形状関数自体を変更しておく．
変更した形状関数は，`ModTriShape`にあるように，
3点の場合は，

```math
\begin{align}
N _0 &= t _0 \\
N _1 &= -t _1(t _0-1) \\
N _2 &= (t _0-1)(t _1-1)
\end{align}
```

6点の場合は，

```math
\begin{align}
N _0 &= t _0(2t _0-1) \\
N _1 &= t _1(2t _1-1) \\
N _2 &= (1-t _0-t _1)(2(1-t _0-t _1)-1) \\
N _3 &= 4t _0t _1 \\
N _4 &= 4t _1(1-t _0-t _1) \\
N _5 &= 4t _0(1-t _0-t _1)
\end{align}
```

[../../include/basic_arithmetic_array_operations.hpp#L614](../../include/basic_arithmetic_array_operations.hpp#L614)



![](sample_tri.png)


[./0README.cpp#L1](./0README.cpp#L1)


## ⛵️ラグランジュ補間 

与えられたデータ点を通る多項式を求める方法の一つにラグランジュ補間がある．

```math
f(x) = \sum _{i=0}^n\dfrac{\prod _{j=0,j\neq i}^n{(x - x _j)}}{\prod _{j=0,j\neq i}^n{(x _i - x _j)}}y _i
```

補間の式の分母は，補間したい横軸$`x`$に依存せず，与えられた横軸データのみに依存するので，
データが与えられたら1度だけ計算しておけばよい．

![](sample_lag.png)


[./LagrangianInterpolation.cpp#L10](./LagrangianInterpolation.cpp#L10)


## ⛵️三角形補間 

## ⛵️⛵️範囲を修正した三角形形状関数  

普通の三角形形状関数は，$`{\mathbf N}=(N _0,N _1,N _2) = (t _0,t _1,1-t _0-t _1)`$．
これを使った，$`{\rm Dot}({\mathbf N},\{{\mathbf X _0},{\mathbf X _1},{\mathbf X _2}\})`$は，$`t _0,t _1=[0,1]`$で平行四辺形を作る．
$`t _0,t _1=[0,1]`$の範囲で，三角形を形成するように変数変換したいことがある．
そのたびに，変数変換をプログラムするのは面倒なので，予め形状関数自体を変更しておく．
変更した形状関数は，`ModTriShape`にあるように，
3点の場合は，

```math
\begin{align}
N _0 &= t _0 \\
N _1 &= -t _1(t _0-1) \\
N _2 &= (t _0-1)(t _1-1)
\end{align}
```

6点の場合は，

```math
\begin{align}
N _0 &= t _0(2t _0-1) \\
N _1 &= t _1(2t _1-1) \\
N _2 &= (1-t _0-t _1)(2(1-t _0-t _1)-1) \\
N _3 &= 4t _0t _1 \\
N _4 &= 4t _1(1-t _0-t _1) \\
N _5 &= 4t _0(1-t _0-t _1)
\end{align}
```

[../../include/basic_arithmetic_array_operations.hpp#L614](../../include/basic_arithmetic_array_operations.hpp#L614)



![](sample_tri.png)


[./TriShape.cpp#L1](./TriShape.cpp#L1)


---
