# Contents
- [🐋 Smoothed Particle Hydrodynamics (SPH) ISPH EISPH](#🐋-Smoothed-Particle-Hydrodynamics-(SPH)-ISPH-EISPH)
        - [🪼 CFL条件の設定](#🪼-CFL条件の設定)
    - [⛵ N.S.方程式を解く前の準備](#⛵-N.S.方程式を解く前の準備)
    - [⛵ N.S.方程式を解く前の準備](#⛵-N.S.方程式を解く前の準備)
        - [🪼 `setCorrectionMatrix_gradient`について](#🪼-`setCorrectionMatrix_gradient`について)
        - [🪼 `setCorrectionMatrix_laplacian`について](#🪼-`setCorrectionMatrix_laplacian`について)
            - [🪸 `interp_normal_original`の計算](#🪸-`interp_normal_original`の計算)
            - [🪸 `setCorrectionMatrix`で壁粒子の演算修正用行列を計算](#🪸-`setCorrectionMatrix`で壁粒子の演算修正用行列を計算)
        - [🪼 壁面粒子の抽出と値の計算](#🪼-壁面粒子の抽出と値の計算)
            - [🪸 `isCaptured`が`true`の壁面粒子の流速の計算](#🪸-`isCaptured`が`true`の壁面粒子の流速の計算)
            - [🪸 `isCaptured`の決定](#🪸-`isCaptured`の決定)
        - [🪼 流体の法線方向の計算と水面の判定](#🪼-流体の法線方向の計算と水面の判定)
            - [🪸 流体の法線方向の計算](#🪸-流体の法線方向の計算)
            - [🪸 水面の判定](#🪸-水面の判定)
    - [⛵ 粘性項$`\nabla^2 {\bf u} _i`$の計算（`calcLaplacianU`）](#⛵-粘性項$`\nabla^2-{\bf-u}-_i`$の計算（`calcLaplacianU`）)
    - [⛵ 粘性項$`\nabla^2 {\bf u} _i`$の計算（`calcLaplacianU`）](#⛵-粘性項$`\nabla^2-{\bf-u}-_i`$の計算（`calcLaplacianU`）)
    - [⛵ 粘性項$`\nabla^2 {\bf u} _i`$の計算（`calcLaplacianU`）](#⛵-粘性項$`\nabla^2-{\bf-u}-_i`$の計算（`calcLaplacianU`）)
    - [⛵ ポアソン方程式 $`\nabla ^{n+1} \cdot \left(\frac{1}{\rho ^n} \nabla ^{n} p \right)=b`$](#⛵-ポアソン方程式-$`\nabla-^{n+1}-\cdot-\left(\frac{1}{\rho-^n}-\nabla-^{n}-p-\right)=b`$)
        - [🪼 ポアソン方程式](#🪼-ポアソン方程式)
        - [🪼 右辺，$`b`$，`PoissonRHS`について](#🪼-右辺，$`b`$，`PoissonRHS`について)
        - [🪼 左辺について](#🪼-左辺について)
        - [🪼 水面の計算補助粒子`auxiliaryPoints`](#🪼-水面の計算補助粒子`auxiliaryPoints`)
        - [🪼 次時刻の発散演算，$`\nabla^{n+1} \cdot {\bf b}^n = \sum _j \dfrac{m _j}{\rho _j^{n+1}}({\bf b} _j^n-{\bf b} _i^n)\cdot \nabla W({\bf x} _i^{n+1},{\bf x} _j^{n+1},h)`$](#🪼-次時刻の発散演算，$`\nabla^{n+1}-\cdot-{\bf-b}^n-=-\sum-_j-\dfrac{m-_j}{\rho-_j^{n+1}}({\bf-b}-_j^n-{\bf-b}-_i^n)\cdot-\nabla-W({\bf-x}-_i^{n+1},{\bf-x}-_j^{n+1},h)`$)
        - [🪼 圧力の安定化](#🪼-圧力の安定化)
    - [⛵ ポアソン方程式の解法](#⛵-ポアソン方程式の解法)
    - [⛵ 圧力勾配$`\nabla p^{n+1}`$の計算](#⛵-圧力勾配$`\nabla-p^{n+1}`$の計算)
    - [⛵ 注意点](#⛵-注意点)
    - [⛵ 出力](#⛵-出力)
    - [⛵ 出力（ポリゴン）](#⛵-出力（ポリゴン）)
- [🐋 実行方法](#🐋-実行方法)
- [🐋 Bucketを用いた粒子探索のテスト](#🐋-Bucketを用いた粒子探索のテスト)
- [🐋 テスト](#🐋-テスト)
    - [⛵ 核関数のテスト](#⛵-核関数のテスト)


---
# 🐋 Smoothed Particle Hydrodynamics (SPH) ISPH EISPH 

[README_ABSTRACT.md](./README_ABSTRACT.md)

[README_FOR_STUDENTS.md](./README_FOR_STUDENTS.md)

[./SPH.hpp#L149](./SPH.hpp#L149)

---
### 🪼 CFL条件の設定 

$`\max({\bf u}) \Delta t \leq c _{v} h \cap \max({\bf a}) \Delta t^2 \leq c _{a} h`$
を満たすように，毎時刻$`\Delta t`$を設定する．
$`c _v=0.1,c _a=0.1`$としている．

[./SPH_Functions.hpp#L200](./SPH_Functions.hpp#L200)

---
## ⛵ N.S.方程式を解く前の準備 

壁粒子の法線ベクトル`p->v_to_surface_SPH`を計算する．

[./main.cpp#L270](./main.cpp#L270)

---
## ⛵ N.S.方程式を解く前の準備

[./SPH0_setWall_Freesurface.hpp#L11](./SPH0_setWall_Freesurface.hpp#L11)

---
### 🪼 `setCorrectionMatrix_gradient`について 

[Morikawa et al. (2023)](https://doi.org/10.1016/j.jcpx.2023.100125)で紹介されていた，Randles and Libersky (1996)の勾配演算の精度を改善する行列を計算する．
勾配の演算を修正する行列は，renormalization tensorと呼ばれ，
よく$`i`$番目の粒子に対する修正行列は$`{\bf B} _i`$と書く．
プログラム上では[`grad_corr_M`](./SPH0_setWall_Freesurface.hpp#L68)としている．

```math
{\bf B} _i = \left(\sum _j V _j ({\bf x} _j-{\bf x} _i) \otimes \nabla W _{ij}\right)^{-1}
```

⚠️ `isCaptured`を先に計算しておく必要がある．`isCaptured`が`false`の場合は，`grad_corr_M`は単位行列になる．

[./SPH0_setWall_Freesurface.hpp#L19](./SPH0_setWall_Freesurface.hpp#L19)

### 🪼 `setCorrectionMatrix_laplacian`について 

\cite{Fatehi2011}が提案した，ラプラシアンの演算の精度を改善する行列を計算する．
プログラム上では`laplacian_corr_M`としている．
多くの場合，流速のラプラシアンは次のように計算される．

```math
\nabla^2 {\bf u} _i=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}\cdot{\bf x} _{ij}}
```

これから，流速の勾配を引くことが１段階目の修正である．

```math
\nabla^2 {\bf u} _i=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j - {\bf x} _{ij}\cdot{\nabla \otimes {\bf u} _{ij}})
```

さらに，renomalization tensorを使って，次のように修正する．

```math
\nabla^2 {\bf u} _i=\sum _{j} {\hat{\bf B}} _i:{\bf A} _{ij}({\bf u} _i - {\bf u} _j - {\bf x} _{ij}\cdot{\nabla \otimes {\bf u} _{ij}})
\quad {\bf A} _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\otimes\nabla W _{ij}}{{\bf x} _{ij}\cdot{\bf x} _{ij}}
```

```math
\begin{align}
{\bf M} = \left(\sum _j V _j ({\bf x} _j-{\bf x} _i) \otimes ({\bf x} _j-{\bf x} _i) \nabla^2 W _{ij}\right)^{-1}
\end{align}
```

⚠️ ラプラシアンの修正行列を計算するためには，先に`setCorrectionMatrix_gradient`を計算しておく必要がある．

[./SPH0_setWall_Freesurface.hpp#L259](./SPH0_setWall_Freesurface.hpp#L259)

---
#### 🪸 `interp_normal_original`の計算 

流体粒子と同じ影響半径を使ってしまうと，流体粒子が参照できる範囲ギリギリにある壁粒子の法線方向の値が不正確になる．
そのため，流体粒子の影響半径よりも広い半径を使って，`q->interp_normal_original`の法線方向を計算することが，重要である．
少し大きい半径を`captureRange`としている．

[./SPH0_setWall_Freesurface.hpp#L468](./SPH0_setWall_Freesurface.hpp#L468)

#### 🪸 `setCorrectionMatrix`で壁粒子の演算修正用行列を計算 

`setCorrectionMatrix`で壁粒子の演算修正用行列を計算する．

[./SPH0_setWall_Freesurface.hpp#L534](./SPH0_setWall_Freesurface.hpp#L534)

### 🪼 壁面粒子の抽出と値の計算

[./SPH0_setWall_Freesurface.hpp#L563](./SPH0_setWall_Freesurface.hpp#L563)

#### 🪸 `isCaptured`が`true`の壁面粒子の流速の計算 

次のようにして，鏡写しのように流速を計算する．

```cpp
q->U_SPH = Reflect(q->U_SPH, q->v_to_surface_SPH)
```

[./SPH0_setWall_Freesurface.hpp#L605](./SPH0_setWall_Freesurface.hpp#L605)

---
#### 🪸 `isCaptured`の決定 

法線方向`interp_normal_original`を使って，流体粒子に近くかつ向かい合う方向にある壁粒子を抽出する．
計算に使用する壁粒子を決定し，使用する場合`isCaptured`を`true`にする．

[./SPH0_setWall_Freesurface.hpp#L485](./SPH0_setWall_Freesurface.hpp#L485)

---
### 🪼 流体の法線方向の計算と水面の判定

[./SPH0_setWall_Freesurface.hpp#L698](./SPH0_setWall_Freesurface.hpp#L698)

#### 🪸 流体の法線方向の計算 

✅ [単位法線ベクトル](./SPH0_setWall_Freesurface.hpp#L850): $`{\bf n} _i = {\rm Normalize}\left(-\sum _j {\frac{m _j}{\rho _j} \nabla W _{ij} }\right)`$

単位法線ベクトルは，`interp_normal`としている．

[./SPH0_setWall_Freesurface.hpp#L728](./SPH0_setWall_Freesurface.hpp#L728)

#### 🪸 水面の判定 

水面の判定条件は，少し複雑である．

[./SPH0_setWall_Freesurface.hpp#L1024](./SPH0_setWall_Freesurface.hpp#L1024)

---
## ⛵ 粘性項$`\nabla^2 {\bf u} _i`$の計算（`calcLaplacianU`） 

✅ [流速のラプラシアンの計算方法](./SPH1_lap_div_U3.hpp#L107): $`\nabla^2 {\bf u} _i=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

✅ [流速の発散の計算方法](./SPH1_lap_div_U3.hpp#L97): $`\nabla\cdot{\bf u} _i=\sum _{j}\frac{m _j}{\rho _j}({{\bf u} _j-{\bf u} _i}) \cdot\nabla W _{ij}`$

[./SPH1_lap_div_U.hpp#L7](./SPH1_lap_div_U.hpp#L7)

## ⛵ 粘性項$`\nabla^2 {\bf u} _i`$の計算（`calcLaplacianU`） 

✅ [流速のラプラシアンの計算方法](./SPH1_lap_div_U3.hpp#L107): $`\nabla^2 {\bf u} _i=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

✅ [流速の発散の計算方法](./SPH1_lap_div_U3.hpp#L97): $`\nabla\cdot{\bf u} _i=\sum _{j}\frac{m _j}{\rho _j}({{\bf u} _j-{\bf u} _i}) \cdot\nabla W _{ij}`$

[./SPH1_lap_div_U2.hpp#L7](./SPH1_lap_div_U2.hpp#L7)

## ⛵ 粘性項$`\nabla^2 {\bf u} _i`$の計算（`calcLaplacianU`） 

✅ [流速のラプラシアンの計算方法](./SPH1_lap_div_U3.hpp#L107): $`\nabla^2 {\bf u} _i=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

✅ [流速の発散の計算方法](./SPH1_lap_div_U3.hpp#L97): $`\nabla\cdot{\bf u} _i=\sum _{j}\frac{m _j}{\rho _j}({{\bf u} _j-{\bf u} _i}) \cdot\nabla W _{ij}`$

[./SPH1_lap_div_U3.hpp#L7](./SPH1_lap_div_U3.hpp#L7)

---
## ⛵ ポアソン方程式 $`\nabla ^{n+1} \cdot \left(\frac{1}{\rho ^n} \nabla ^{n} p \right)=b`$ 

### 🪼 ポアソン方程式 

次の時刻の流れ場を発散なし$`\nabla\cdot{\bf u}^{n+1}=0`$としてくれる
$`\frac{D {\bf u}}{D t} =-\frac{1}{\rho} \nabla p +\nu \nabla^2 {\bf u}^n+{\bf g}`$を使って，流速と粒子位置を時間発展させたい．
そのためには，圧力$`p^{n+1}`$を適切に決める必要がある．

$`\frac{D {\bf u}}{D t}`$は．$`\frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t}`$と離散化し条件を考えてみる．

```math
\rho\frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t} =- \nabla p +\mu \nabla^2 {\bf u}^n+\rho{\bf g}
```

次時刻の発散の演算は，次時刻における粒子配置に基づき行われるので，現在の粒子配置に基づく発散演算とは区別すべきである．
現在の微分演算を$`\nabla^{n}`$とし，次時刻の微分演算を$`\nabla^{n+1}`$とする．
$`\nabla^{n+1}`$を上の式に作用させると，

```math
\nabla^{n+1}\cdot {(\rho {\bf u}^{n+1})} = \nabla^{n+1} \cdot{(\rho {\bf u}^n)} - \Delta t \nabla^{n+1} \cdot\left( \nabla p-\mu \nabla^{n2} {\bf u}^n-{\rho\bf g}\right)
```

右辺がゼロとなれば，次時刻の流速の発散がゼロ，$`\nabla^{n+1} \cdot (\rho{\bf u}^{n+1})=0`$になる：

```math
\begin{align*}
&&0 &= \nabla^{n+1} \cdot (\rho{\bf u}^{n}) - \Delta t \nabla^{n+1} \cdot\left(\nabla p-\mu \nabla^{n2} {\bf u}^n-{\rho \bf g}\right)\\
&\rightarrow&\nabla^{n+1} \cdot \nabla p &= \frac{1}{\Delta t}\nabla^{n+1} \cdot (\rho {\bf u}^{n}) + \nabla^{n+1} \cdot\left(\mu \nabla^{n2} {\bf u}^n  + {\rho\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \nabla p &= \nabla^{n+1} \cdot\left(\frac{1}{\Delta t} (\rho {\bf u}^{n}) +\mu \nabla^{n2} {\bf u}^n  + {\rho\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \nabla p & = \nabla^{n+1} \cdot {\bf b}^n= b,\quad  {\bf b}^n=\frac{1}{\Delta t}(\rho {\bf u}^{n}) +\mu \nabla^{n2} {\rho\bf u}^n
\end{align*}
```

重力の発散はゼロなので消した．
ここで，$`\nabla p`$は敢えて，微分演算子や圧力のタイムステップを書かなかった．
$`\nabla p`$は，次時刻の流速の発散をゼロにするためだけの未知ベクトルであって，タイムステップを考える必要はない．

ここでは，演算子が異なると複雑になるので，同じ$`\nabla^{n+1}`$を使うことにする．
これに伴って，次時刻の流速$`{\bf u}^{n+1}`$を計算する際に用いる圧力勾配は，$`\nabla^{n+1} p`$として計算しなければならないことに注意する．

次時刻の微分演算子を使うことにして，$`\nabla p`$を$`\nabla^{n+1} p`$とし，次のようなポアソン方程式を得る．

```math
\nabla^{n+1}\cdot \nabla^{n+1} p =b
```

ただし，水面粒子は$`p=0`$として，上の方程式は使わない．

### 🪼 右辺，$`b`$，`PoissonRHS`について 

この$`b`$を`PoissonRHS`とする．（仮流速は$`{\bf u}^\ast = \frac{\Delta t}{\rho}{\bf b}^n`$と同じ）．
$`{\bf b}^n`$ （[`Poisson_b_vector`](./SPH1_lap_div_U3.hpp#L282)）が計算できるように，$`{\bf u}^n`$と$`\nabla^2 {\bf u}^n`$を計算しておく．

✅ [発散の計算方法](not found): $`b=\nabla\cdot{\bf b}^n=\sum _{j}\frac{m _j}{\rho _j}({\bf b} _j^n-{\bf b} _i^n)\cdot\nabla W _{ij}`$

### 🪼 左辺について 

壁粒子の圧力は時間積分して計算しないので，毎時刻，壁粒子の$`p`$を計算する必要がある．

✅ [ラプラシアンの計算方法](./SPH1_lap_div_U3.hpp#L99): $`\nabla^2 p=\sum _{j}A _{ij}(p _i - p _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

✅ [ラプラシアンの計算方法](not found): $`\nabla^2 p=\sum _{j}A _{ij}(p _i - p _j),\quad A _{ij} = \frac{8 m _j}{(\rho _i+\rho _j)}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

⚠️ 密度$\rho$が粒子に関わらず一定の場合，上の２式は同じになる．しかし，補助粒子の密度は，他の粒子と異なるので，[２つ目のラプラシアンの計算方法](not found)を使うべきだろう．

**ISPH**

- ISPHは作ったポアソン方程式を作成し解くことで圧力を計算する

**EISPH**

1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p`$だけを使って近似）
2. 流体粒子の圧力$`p`$の計算

[EISPHの圧力の設定方法](not found)


$\sum _j A _{ij} (p _i-p _j) = b$において，$p _j^{\rm new} \approx p _j^{\rm old}$とすると，

```math
p _i^{\rm new} = \frac{b + \sum _j A _{ij} p _j^{\rm old}}{\sum _j A _{ij}}
```

となる．

<!---
### 🪼 水面の計算補助粒子`auxiliaryPoints` 

水面においては，流速の発散ゼロ$`\nabla^{n+1} {\bf u}^{n+1}=0`$と$`p^{n+1}=0`$が満たされる必要がある．
水面外部には，粒子がないので，求めた水面圧力は，ゼロであっても，圧力勾配は誤差を含み，$`\nabla^{n+1} {\bf u}^{n+1}=0`$は満足されない．
そこで，[水面の計算補助粒子](not found)を水面外部に追加し，この点を適切計算することで，$`\nabla^{n+1} {\bf u}^{n+1}=0`$が満足されるように工夫する．
--->

[./SPH2_FindPressure.hpp#L7](./SPH2_FindPressure.hpp#L7)

---
### 🪼 次時刻の発散演算，$`\nabla^{n+1} \cdot {\bf b}^n = \sum _j \dfrac{m _j}{\rho _j^{n+1}}({\bf b} _j^n-{\bf b} _i^n)\cdot \nabla W({\bf x} _i^{n+1},{\bf x} _j^{n+1},h)`$ 

$`\nabla^{n+1}`$の計算には，$`\rho^{n+1}`$, $`{\bf x}^{n+1}= {\bf x}^{n} + {\bf u}^{n+1} \Delta t`$が必要である．

* [次時刻の粒子体積](./SPH_Functions.hpp#L310)
* [次時刻の粒子密度](./SPH_Functions.hpp#L299)
* [次時刻の粒子位置](./SPH_Functions.hpp#L315)

[./SPH2_FindPressure.hpp#L104](./SPH2_FindPressure.hpp#L104)

---
* `ROW`は，どの粒子も方程式を保存するかを表す．
* `pO_center`は，圧力の方程式を立てる際の座標を表す（基本的には`ROW`の位置と同じ）．
* `pO`は，影響半径などの情報として使う粒子を表す（基本的には`ROW`と同じ）．

|方程式|目的|
|:---------|---|
| ☑️ [ポアソン方程式](./SPH2_FindPressure.hpp#L237)              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
| ☐ [不透過条件](not found)         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
| ☐ [大気圧条件](not found) | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない． |

各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．

[./SPH2_FindPressure.hpp#L148](./SPH2_FindPressure.hpp#L148)

壁面粒子の圧力の設定方法

ポアソン方程式を解いた場合：壁近傍の粒子が内部方向への圧力を受ける．

壁の法線方向にある流体の圧力を，壁粒子の圧力とした場合（若干の修正をするが）：あまり力を受けない．

[./SPH2_FindPressure.hpp#L327](./SPH2_FindPressure.hpp#L327)

---
### 🪼 圧力の安定化 

$`b = \nabla \cdot {{\bf b}^n} + \alpha \frac{\rho _w - \rho^\ast}{{\Delta t}^2}`$として計算を安定化させる場合がある．
$`\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t`$と近似すると，

```math
\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t,\quad
\frac{D\rho^\ast}{Dt} = - \rho \nabla\cdot{\bf u}^\ast,\quad
\nabla\cdot{\bf u}^\ast = \frac{\Delta t}{\rho} \nabla\cdot{\bf b}^n
```

であることから，$`(\rho _w - \rho^\ast) / {\Delta t^2}`$は，$`\nabla\cdot{\bf b}^n`$となって同じになる．

しかし，実際には，$`\rho^\ast`$は，$\nabla \cdot {{\bf b}^n}`$を使わずに，つまり発散演算を行わずに評価するので，
計算上のようにはまとめることができない．

$`\rho^\ast`$を計算する際に，$`\rho^\ast = \rho _w + \frac{D\rho^\ast}{Dt}\Delta t`$を使った場合，確かに上のようになるが，
実際に粒子を仮位置に移動させその配置から$\rho^*$を計算した場合は，数値計算上のようにまとめることはできない．

`PoissonRHS`,$`b`$の計算方法と同じである場合に限る．
もし，計算方法が異なれば，計算方法の違いによって，安定化の効果も変わってくるだろう．

[./SPH2_FindPressure.hpp#L408](./SPH2_FindPressure.hpp#L408)

---
## ⛵ ポアソン方程式の解法 

ISPHのポアソン方程式を解く場合，[ここではGMRES法](./SPH2_FindPressure.hpp#L658)を使う．

[./SPH2_FindPressure.hpp#L509](./SPH2_FindPressure.hpp#L509)

---
## ⛵ 圧力勾配$`\nabla p^{n+1}`$の計算 

✅ [勾配の計算方法](./SPH3_grad_P.hpp#L166): $`\nabla p _i = \rho _i \sum _{j} m _j (\frac{p _i}{\rho _i^2} + \frac{p _j}{\rho _j^2}) \nabla W _{ij}`$

✅ [勾配の計算方法](./SPH3_grad_P.hpp#L120): $`\nabla p _i = \sum _{j} \frac{m _j}{\rho _i} \left(p _j - p _i\right) \nabla W _{ij}`$

✅ [勾配の計算方法](./SPH3_grad_P.hpp#L131): $`\nabla p _i = \sum _{j} \frac{m _j}{\rho _j} p _j \nabla W _{ij}`$

💡 圧力の方程式を立てる際に，左辺の密度として，流速の発散から見積もった$`\rho^({\rm next})`$を使うことは，
言い換えれば，N.S.方程式の圧力項の計算には，$`\rho^({\rm next})`$を使うと決めたことになる．
なので，圧力勾配$`\nabla p^{n+1}`$の計算にも，$`\rho^({\rm next})`$を使わなければならない．

[./SPH3_grad_P.hpp#L11](./SPH3_grad_P.hpp#L11)

$`\dfrac{D{\bf u}^n}{Dt} = - \frac{1}{\rho} \nabla p^{n+1} + \nu \nabla^2 {\bf u}^n + {\bf g}`$
が計算できた．

[./SPH3_grad_P.hpp#L171](./SPH3_grad_P.hpp#L171)

---
## ⛵ 注意点 

⚠️ 計算がうまく行く設定を知るために，次の箇所をチェックする．

**NEW**

- [壁粒子の速度の決定方法](./SPH0_setWall_Freesurface.hpp#L687)
- [Poissonにおいてどのようにbベクトルを使うか](not found)
- [Poissonにおいてどのようにbベクトルを使うか](not found)
- どのように[壁粒子のb](not found)/[流体粒子のb](./SPH1_lap_div_U3.hpp#L283)を作るか

**壁粒子**

- [壁粒子のラプラシアンの計算方法](./SPH1_lap_div_U3.hpp#L281)
- [圧力の計算方法](./SPH2_FindPressure.hpp#L126)
- [どの位置において方程式を立てるか](./SPH2_FindPressure.hpp#L289)
- [流体として扱う壁粒子を設定するかどうか](./SPH0_setWall_Freesurface.hpp#L455)/[視野角に流体粒子が含まない壁粒子は除外する](not found)
- [壁粒子の圧力をどのように壁面にマッピングするか](not found)
- [壁粒子の法線方向ベクトルの計算方法](./SPH0_setWall_Freesurface.hpp#L850)
- [反射の計算方法](./SPH_Functions.hpp#L460)

**水面粒子**

- [水面粒子の圧力をゼロにするかどうか](not found)
- [補助粒子の設定はどうなっているか](not found)

**その他**

- [密度を更新するかどうか](./SPH_Functions.hpp#L532)
- [圧力の安定化をするかどうか](not found)
- [ルンゲクッタの段数](./from os.py#L145)


壁のwall_as_fluidは繰り返しで計算するのはどうか？

[./SPH_Functions.hpp#L548](./SPH_Functions.hpp#L548)

## ⛵ 出力

[./main.cpp#L417](./main.cpp#L417)

## ⛵ 出力（ポリゴン）

[./main.cpp#L631](./main.cpp#L631)

---
# 🐋 実行方法 

ファイルをダウンロードして，`build_sph`ディレクトリに移動．
⚠️上書きされるので注意．

```sh
git clone https://github.com/tomoakihirakawa/cpp.git
cd ./cpp/builds/build_sph
```

`clean`でCMake関連のファイルを削除して（ゴミがあるかもしれないので），
`cmake`で`Makefile`を生成して，`make`でコンパイルする．

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

次に，入力ファイルを生成．

```sh
python3 input_generator.py
```

例えば，`./input_files/static_pressure_PS0d0125_CSML2d4_RK1`が生成される．
入力ファイルを指定して実行．

```sh
./main ./input_files/static_pressure_PS0d0125_CSML2d4_RK1
```

[./main.cpp#L1](./main.cpp#L1)

---
# 🐋 Bucketを用いた粒子探索のテスト 

Smoothed Particle Hydrodynamics (SPH)では，効率的な近傍粒子探査が必要となる．
このコードでは，Bucketを用いた粒子探索のテストを行う．

結果はVTKファイルに出力される．
* 全ての粒子を表示したものは`all.vtp`
* 中心の粒子を表示したものは`center*.vtp`
* 中心の粒子が探査したセル内にある粒子を表示したものは`inCell*.vtp`
* セル内かつ球内にある粒子を表示したものは`inSphere*.vtp`

- 各セルにある粒子を表示したものは`each_cell*.vtp`
- 各セルの中心位置を表示したものは`each_cell_position*.vtp`

[./test_Buckets.cpp#L2](./test_Buckets.cpp#L2)

---
# 🐋 テスト 

## ⛵ 核関数のテスト 

<!-- Key SPH:kernelFunctions not found -->

プログラムした[3次スプライン関数](not found)と[5次スプライン関数](not found)のテストコード

```sh
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_KernelFunctions.cpp
make
./test_KernelFunctions
```

* 関数の形状を確認．
* 体積積分が1になるかどうかを確認．

| 分割数$N$，体積$`V=(\frac{2r}{N})^3`$   | Sum for 3rd Order | Sum for 5th Order |
| --- | ----------------- | ----------------- |
| 5   | 1.00527           | 0.999206          |
| 10  | 1.00011           | 1.00005           |
| 15  | 0.999972          | 0.999999          |
| 20  | 1                 | 1                 |
| 25  | 1                 | 1                 |

[./test_KernelFunctions.cpp#L1](./test_KernelFunctions.cpp#L1)

---
