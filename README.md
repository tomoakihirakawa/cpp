# Contents

- [🐋 Runge-Kutta Integration of ODE](#🐋-Runge-Kutta-Integration-of-ODE)
- [🐋 Boundary Element Method (BEM-MEL)](<#🐋-Boundary-Element-Method-(BEM-MEL)>)
  - [⛵️ 流速の計算方法](#⛵️-流速の計算方法)
    - [⚓️ 修正流速](#⚓️-修正流速)
  - [⛵️ 境界条件の設定](#⛵️-境界条件の設定)
  - [⛵️ 境界値問題](#⛵️-境界値問題)
    - [⚓️ BIE の離散化](#⚓️-BIEの離散化)
    - [⚓️ 多重節点](#⚓️-多重節点)
- [🐋 準ニュートン法](#🐋-準ニュートン法)
  - [⛵️ ヘッセ行列を利用したニュートン法](#⛵️-ヘッセ行列を利用したニュートン法)
- [🐋 Smoothed Particle Hydrodynamics (SPH) ISPH EISPH](<#🐋-Smoothed-Particle-Hydrodynamics-(SPH)-ISPH-EISPH>)
  - [⛵️ 概要](#⛵️-概要)
    - [⚓️ 前準備](#⚓️-前準備)
    - [⚓️ フラクショナルステップを使って初期値問題を解く](#⚓️-フラクショナルステップを使って初期値問題を解く)
    - [⚓️ 壁面粒子の流速と圧力](#⚓️-壁面粒子の流速と圧力)
    - [⚓️ $`\nabla^2 {\bf u}`$の計算](#⚓️-$`\nabla^2-{\bf-u}`$の計算)
    - [⚓️ `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算](#⚓️-`PoissonRHS`,$`b`$と$`\nabla^2-p^{n+1}`$における$`p^{n+1}`$の係数の計算)
    - [⚓️ 圧力の安定化](#⚓️-圧力の安定化)
    - [⚓️ 圧力勾配$`\nabla p^{n+1}`$の計算 -> $`{D {\bf u}}/{Dt}`$の計算](#⚓️-圧力勾配$`\nabla-p^{n+1}`$の計算-->-$`{D-{\bf-u}}/{Dt}`$の計算)
  - [⛵️ Bucket を用いた粒子探索のテスト](#⛵️-Bucketを用いた粒子探索のテスト)
  - [⛵️ 核関数](#⛵️-核関数)
  - [⛵️ Compressed Sparse Row (CSR)](<#⛵️-Compressed-Sparse-Row-(CSR)>)
  - [⛵️ 一般化最小残差法(GMRES)](<#⛵️-一般化最小残差法(GMRES)>)
  - [⛵️ ArnoldiProcess](#⛵️-ArnoldiProcess)

---

# 🐋 Runge-Kutta Integration of ODE

This C++ program demonstrates the application of various Runge-Kutta methods (first to fourth order) for solving a first-order ordinary differential equation (ODE).
![](builds/build_ODE/runge_kutta/res.png)

[./builds/build_ODE/runge_kutta/main.cpp#L1](./builds/build_ODE/runge_kutta/main.cpp#L1)

---

[![Banner](builds/build_bem/banner.png)](banner.png)

# 🐋 Boundary Element Method (BEM-MEL)

[./builds/build_bem/BEM.hpp#L1](./builds/build_bem/BEM.hpp#L1)

## ⛵️ 流速の計算方法

[./builds/build_bem/BEM_calculateVelocities.hpp#L7](./builds/build_bem/BEM_calculateVelocities.hpp#L7)

### ⚓️ 修正流速

求めた流速から，次の時刻の境界面$`\Omega(t+\Delta t)`$を見積もり，その面上で節点を移動させ歪さを解消する．
修正ベクトルは，$`\Delta t`$で割り，求めた流速$`\nabla \phi`$に足し合わせて，節点を時間発展させる．

ノイマン節点も修正流速を加え時間発展させる．
ただし，ノイマン節点の修正流速に対しては，節点が水槽の角から離れないように，工夫を施している．

[./builds/build_bem/BEM_calculateVelocities.hpp#L354](./builds/build_bem/BEM_calculateVelocities.hpp#L354)

## ⛵️ 境界条件の設定

1. 流体節点が接触する構造物面を保存する
2. 面の境界条件：３節点全てが接触している流体面は Neumann 面，それ以外は Dirichlet 面とする
3. 辺の境界条件：辺を含む２面が Neumann 面なら Neumann 辺，２面が Dirichlet 面なら Dirichlet 面，それ以外は CORNER とする．
4. 点の境界条件：点を含む面全てが Neumann 面なら Neumann 点，面全てが Dirichlet 面なら Dirichlet 点，それ以外は CORNER とする．

[./builds/build_bem/BEM_setBoundaryConditions.hpp#L7](./builds/build_bem/BEM_setBoundaryConditions.hpp#L7)

## ⛵️ 境界値問題

### ⚓️ BIE の離散化

$`\phi`$と$`\phi _n`$に関する BIE は，

$$
\alpha ({\bf{a}})\phi ({\bf{a}}) = \iint _\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
$$

これを線形三角要素と Gauss-Legendre 積分で離散化すると，

$$
\alpha _{i _\circ}(\phi) _{i _\circ}=-\sum\limits _{k _\vartriangle}\sum\limits _{{\xi _1}} {\sum\limits _{{\xi _0}} {\left( {{w _0}{w _1}\left( {\sum\limits _{j=0}^2 {{{\left( {{\phi _n}} \right)} _{k _\vartriangle,j }}{N _{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x} _{i \ _\circ}}} \|}}\left\|\frac{{\partial{\bf{x}}}}{{\partial{\xi _0}}} \times \frac{{\partial{\bf{x}}}}{{\partial{\xi _1}}}\right\|} \right)} }-\sum\limits _{k _\vartriangle}\sum\limits _{{\xi _1}} \sum\limits _{{\xi _0}} {\left( {{w _0}{w _1}\left({\sum\limits _{j =0}^2{{{\left( \phi  \right)} _{k _\vartriangle,j }}{N _{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{{\bf x} _{i _\circ}} - {\bf{x}}\left( \pmb{\xi } \right)}}{{{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x} _{i \ _\circ}}}\|}^3}}} \cdot\left(\frac{{\partial {\bf{x}}}}{{\partial {\xi _0}}}\times\frac{{\partial {\bf{x}}}}{{\partial {\xi _1}}}\right)}\right)}
$$

[./builds/build_bem/BEM_solveBVP.hpp#L226](./builds/build_bem/BEM_solveBVP.hpp#L226)

### ⚓️ 多重節点

このループでは，ある面`integ_f`に隣接する節点{p0,p1,p2}の列,IGIGn[origin(fixed),p0],...に値が追加されていく．
（p0 が多重接点の場合，適切に p0 と同じ位置に別の変数が設定されており，別の面の積分の際に q0 が参照される．）
p0 は，{面,補間添字}で決定することもできる．
{面,補間添字 0}->p0,{面,補間添字 1}->p1,{面,補間添字 2}->p2 というように．

{面 A,補間添字},{面 B,補間添字},{面 C,補間添字}が全て同じ節点 p0 を指していたとする．
普通の節点なら，IGIGn[origin,{p0,nullptr}]を指す．
多重節点なら，IGIGn[origin,{p0,面 A}],IGIGn[origin,{p0,面 B}]を指すようにする．
この操作を言葉で言い換えると，
「n が不連続に変化する点では，その点の隣接面にそれぞれ対して φn を求めるべきである（φ は同じでも）．」
「n が不連続に変化する点では，どの面を積分するかに応じて，参照する φn を区別し切り替える必要がある．」

//@ さて，この段階で p0 が多重節点であるかどうか判断できるだろうか？

{節点，面}-> 列ベクトルのインデックス を決めれるか？

面を区別するかどうかが先にわからないので，face\*のま s まか nullptr とすべきかわからないということ．．．．

PBF_index[{p, Dirichlet, ある要素}]
は存在しないだろう．Dirichlet 節点は，{p, ある要素}からの寄与を，ある面に

[./builds/build_bem/BEM_solveBVP.hpp#L321](./builds/build_bem/BEM_solveBVP.hpp#L321)

IGIGn は 左辺に IG*φn が右辺に IGn*φ が来るように計算しているため，移項する場合，符号を変える必要がある．
$`IG \phi _n = IGn \phi`$

移項前:
$`\begin{bmatrix}IG _0 & IG _1 & IG _2 & IG _3\end{bmatrix} \begin{bmatrix}\phi _{n0} \\ \phi _{n1} \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}IG _{n0} & IG _{n1} & IG _{n2} & IG _{n3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _1 \\ \phi _2 \\ \phi _3\end{bmatrix}`$

移項後:
$`\begin{bmatrix}IG _0 & -IG _{n1} & IG _2 & IG _3\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}IG _{n0} & -IG _1 & IGn _2 & IG _{n3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}`$

多重節点(1 と 3 が多重節点の場合):
$`\begin{bmatrix}0 & 1 & 0 & 0\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}0 & 0 & 0 & 1\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}`$

[./builds/build_bem/BEM_solveBVP.hpp#L383](./builds/build_bem/BEM_solveBVP.hpp#L383)

---

# 🐋 準ニュートン法

ニュートン法で使うヤコビアンなどを別のものに置き換えた方法．

[./builds/build_root_finding/example_Broyden.cpp#L1](./builds/build_root_finding/example_Broyden.cpp#L1)

## ⛵️ ヘッセ行列を利用したニュートン法

**最適か否かを判断するための関数**は１つだけで，**最適化したい変数は複数**である場合でも，
最適化は，ヘッセ行列を利用したニュートン法によって可能である．
この方法で，変数は，関数を根とするのではなく，関数を最大最小（停留点）とする値へと収束する．

[./builds/build_root_finding/example_NewtonRaphson.cpp#L1](./builds/build_root_finding/example_NewtonRaphson.cpp#L1)

---

[![Banner](builds/build_sph/banner.png)](banner.png)

# 🐋 Smoothed Particle Hydrodynamics (SPH) ISPH EISPH

## ⛵️ 概要

### ⚓️ 前準備

1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL 条件を満たすようにタイムステップ間隔 $`\Delta t`$を設定

### ⚓️ フラクショナルステップを使って初期値問題を解く

4. $`\nabla^2 {\bf u}`$の計算
5. `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算
6. 流速の発散から密度 $`{\rho}^\ast`$を計算
7. 次の時刻の圧力 $`p^{n+1}`$を計算
8. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
9. 流体粒子の圧力$`p^{n+1}`$の計算
10. $`\nabla {p^{n+1}}`$が計算でき， $`\frac{D{\bf u}}{D t}=-\frac{1}{\rho}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}`$（粘性率が一定の非圧縮性流れの加速度）を得る．
11. $`\frac{D\bf u}{Dt}`$を使って，流速を更新．流速を使って位置を更新

[./builds/build_sph/SPH.hpp#L211](./builds/build_sph/SPH.hpp#L211)

ISPH を使えば，水面粒子の圧力を簡単にゼロにすることができる．
$`\nabla \cdot {\bf u}^\ast`$は流ればで満たされれば十分であり，壁面表層粒子の圧力を，壁面表層粒子上で$`\nabla \cdot {\bf u}^\ast`$となるように決める必要はない．

[./builds/build_sph/SPH.hpp#L387](./builds/build_sph/SPH.hpp#L387)

### ⚓️ 壁面粒子の流速と圧力

壁粒子の流速を流体粒子の流速に応じて変化させると計算が煩雑になるので，**ここでは**壁面粒子の流速は常にゼロに設定することにした（ゼロで一定というのは不自然ではない）．
一方，壁粒子の圧力がゼロだとするのは不自然で，流体粒子の圧力$`p^{n+1}`$の計算に悪影響を及ぼす．
なので．壁粒子の圧力は各ステップ毎に計算し直す必要がある．

壁面粒子の圧力は，壁面法線方向流速をゼロにするように設定されるべきだろう．

[./builds/build_sph/SPH_Functions.hpp#L216](./builds/build_sph/SPH_Functions.hpp#L216)

### ⚓️ $`\nabla^2 {\bf u}`$の計算

✅ ラプラシアンの計算方法: $`\nabla^2 {\bf u}=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

[./builds/build_sph/SPH_Functions.hpp#L230](./builds/build_sph/SPH_Functions.hpp#L230)

### ⚓️ `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算

次の時刻の流れ場が発散なし$`\nabla\cdot{\bf u}^{n+1}=0`$であることを保証してくれる圧力を使って，
$`\frac{D {\bf u}}{D t} =-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}`$を決定し，時間発展させたい．
そのような圧力を$`p^{n+1}`$と書くことにする．
そのような圧力の条件は，次のようになる．

$$
\begin{align*}
&&\frac{D {\bf u}}{D t} &=-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}\\
&\rightarrow& \frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t} &=-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}\\
&\rightarrow& \nabla \cdot\left(\frac{\rho}{\Delta t} {\bf u}^{n+1}\right) + \nabla^2 p^{n+1} &= \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}^n+\rho {\bf g}\right)\\
&\rightarrow& \nabla^2 p^{n+1} &= b, \quad b = \nabla \cdot {{\bf b}^n} = \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}+\rho {\bf g}\right)
\end{align*}
$$

この$`b`$を`PoissonRHS`とする．（仮流速は$`{\bf u}^\ast = \frac{\Delta t}{\rho}{\bf b}^n`$である．）

✅ 発散の計算方法: $`b=\nabla\cdot{\bf b}^n=\sum _{j}\frac{m _j}{\rho _j}({\bf b} _j^n-{\bf b} _i^n)\cdot\nabla W _{ij}`$

`PoissonRHS`,$`b`$の計算の前に，$`\mu \nabla^2{\bf u}`$を予め計算しておく．

壁粒子の圧力は時間発展させないので，壁粒子の$`p^n`$を計算する必要がある．順で計算する．

1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
2. 流体粒子の圧力$`p^{n+1}`$の計算

✅ ラプラシアンの計算方法: $`\nabla^2 p^{n+1}=\sum _{j}A _{ij}(p _i^{n+1} - p _j^{n+1}),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

[./builds/build_sph/SPH_Functions.hpp#L302](./builds/build_sph/SPH_Functions.hpp#L302)

### ⚓️ 圧力の安定化

$`b = \nabla \cdot {{\bf b}^n} + \alpha \frac{\rho _w - \rho^\ast}{{\Delta t}^2}`$として計算を安定化させる場合がある．
$`\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t`$と近似すると，

$$
\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t,\quad
\frac{D\rho^\ast}{Dt} = - \rho \nabla\cdot{\bf u}^\ast,\quad
\nabla\cdot{\bf u}^\ast = \frac{\Delta t}{\rho} \nabla\cdot{\bf b}^n
$$

であることから，$`(\rho _w - \rho^\ast) / {\Delta t^2}`$は，$`\nabla\cdot{\bf b}^n`$となって同じになる．

しかし，実際には，$`\rho^\ast`$は，$`\nabla \cdot {{\bf b}^n} `$を使わずに，つまり発散演算を行わずに評価するので，
計算上のようにはまとめることができない．

$`\rho^\ast`$を計算する際に，$`\rho^\ast = \rho _w + \frac{D\rho^\ast}{Dt}\Delta t`$を使った場合，確かに上のようになるが，
実際に粒子を仮位置に移動させその配置から$`\rho^\ast`$を計算した場合は，数値計算上のようにまとめることはできない．

`PoissonRHS`,$`b`$の計算方法と同じである場合に限る．
もし，計算方法が異なれば，計算方法の違いによって，安定化の効果も変わってくるだろう．

[./builds/build_sph/SPH_Functions.hpp#L336](./builds/build_sph/SPH_Functions.hpp#L336)

### ⚓️ 圧力勾配$`\nabla p^{n+1}`$の計算 -> $`{D {\bf u}}/{Dt}`$の計算

✅ 勾配の計算方法: $`\nabla p _i = \rho _i \sum _{j} m _j (\frac{p _i}{\rho _i^2} + \frac{p _j}{\rho _j^2}) \nabla W _{ij}`$

✅ 勾配の計算方法: $`\nabla p _i = \sum _{j} \frac{m _j}{\rho _j} p _j \nabla W _{ij}`$

[./builds/build_sph/SPH_Functions.hpp#L450](./builds/build_sph/SPH_Functions.hpp#L450)

## ⛵️ 核関数

3 次スプライン関数と 5 次スプライン関数の実装とテストコード

- 関数の形状を確認．
- 体積積分が 1 になるかどうかを確認．

[./builds/build_sph/test_KernelFunctions.cpp#L1](./builds/build_sph/test_KernelFunctions.cpp#L1)

---

## ⛵️ Bucket を用いた粒子探索のテスト

Smoothed Particle Hydrodynamics (SPH)では，効率的な近傍粒子探査が必要となる．
このコードでは，Bucket を用いた粒子探索のテストを行う．

結果は VTK ファイルに出力される．

- 全ての粒子を表示したものは`all.vtp`
- 中心の粒子を表示したものは`center*.vtp`
- 中心の粒子が探査したセル内にある粒子を表示したものは`inCell*.vtp`
- セル内かつ球内にある粒子を表示したものは`inSphere*.vtp`

* 各セルにある粒子を表示したものは`each_cell*.vtp`
* 各セルの中心位置を表示したものは`each_cell_position*.vtp`

[./builds/build_sph/test_Buckets.cpp#L1](./builds/build_sph/test_Buckets.cpp#L1)

---

## ⛵️ Compressed Sparse Row (CSR)

CSR は行列を表現する方法の一つである．
この CSR クラスは，std::unordered_map を用いて，行列の非ゼロ要素を表現する．
std::unordered_map の key はポインタであり，value は double である．
CSR クラス自身が，行列の行番号を保存しており，key である CSR クラスは行列の列番号を保存している．

[./builds/build_system_of_linear_eqs/CSR.cpp#L1](./builds/build_system_of_linear_eqs/CSR.cpp#L1)

---

## ⛵️ 一般化最小残差法(GMRES)

- ヘッセンベルグ行列$`H`$
- クリロフ部分空間の直交基底$`V`$
- $`H`$を QR 分解した行列$`Q`$と$`R`$
- $`g`$は行列$`Q`$の最初の列

ArnoldiProcess によって，$`H`$と$`V`$を求める．この ArnoldiProcess クラスの派生クラスとして GMRES を定義している．

[./builds/build_system_of_linear_eqs/GMRES.cpp#L1](./builds/build_system_of_linear_eqs/GMRES.cpp#L1)

---

## ⛵️ ArnoldiProcess

ヘッセンベルグ行列$`H[0:k-1]`$は，A と相似なベクトルであり，同じ固有値を持つ
GMRES で使う場合，$`V0`$には Normalize(b-A.x0)を与える．
x0 は初期値

アーノルディ法は固有値問題の数値解法であり反復解法．
一般的な行列の固有ベクトルと固有値をクリロフ空間の直行基底によって近似する方法計算する方法．
https://en.wikipedia.org/wiki/Arnoldi_iteration

[./include/basic_linear_systems.hpp#L678](./include/basic_linear_systems.hpp#L678)

---
