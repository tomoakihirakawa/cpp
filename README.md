# Contents

- [ArnoldiProcess](#ArnoldiProcess)

- [Runge-Kutta Integration of ODE](#Runge-Kutta-Integration-of-ODE)

- [核関数](#核関数)

- [ISPHとEISPH](#ISPHとEISPH)

- [Bucketを用いた粒子探索のテスト](#Bucketを用いた粒子探索のテスト)

- [壁面粒子の流速と圧力](#壁面粒子の流速と圧力)

    - [`PoissonRHS`と $\nabla^2 p^{n+1}$における $p^{n+1}$の係数の計算](#`PoissonRHS`と-$\nabla^2-p^{n+1}$における-$p^{n+1}$の係数の計算)

- [ヘッセ行列を利用したニュートン法](#ヘッセ行列を利用したニュートン法)

- [準ニュートン法](#準ニュートン法)

- [Compressed Sparse Row (CSR)](#Compressed-Sparse-Row-(CSR))

- [一般化最小残差法(GMRES)](#一般化最小残差法(GMRES))



## ArnoldiProcess
ヘッセンベルグ行列$H[0:k-1]$は，Aと相似なベクトルであり，同じ固有値を持つ
   GMRESで使う場合，$V0$にはNormalize(b-A.x0)を与える．
   x0は初期値

   アーノルディ法は固有値問題の数値解法であり反復解法．
   一般的な行列の固有ベクトルと固有値をクリロフ空間の直行基底によって近似する方法計算する方法．
   https://en.wikipedia.org/wiki/Arnoldi_iteration

[./include/basic_linear_systems.hpp#L678](./include/basic_linear_systems.hpp#L678)


 --- 
## Runge-Kutta Integration of ODE
This C++ program demonstrates the application of various Runge-Kutta methods (first to fourth order) for solving a first-order ordinary differential equation (ODE).
![](./builds/build_ODE/runge_kutta/rk.png)

[./builds/build_ODE/runge_kutta/main.cpp#L1](./builds/build_ODE/runge_kutta/main.cpp#L1)


 --- 
## 核関数
3次スプライン関数と5次スプライン関数の実装とテストコード
* 関数の形状を確認．
* 体積積分が1になるかどうかを確認．

[./builds/build_sph/test_KernelFunctions.cpp#L1](./builds/build_sph/test_KernelFunctions.cpp#L1)


 --- 
## ISPHとEISPH
### 前準備
1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL条件を満たすようにタイムステップ間隔 $\Delta t$を設定

### フラクショナルステップを使って初期値問題を解く
4. ${{\bf u}^\ast}$と ${{\bf x}^\ast}$を計算
5. 流速の発散 ${\nabla \cdot {\bf u}^\ast}$の計算

   - Nomeritae et al. (2016)は， ${{\bf u}^\ast}$と ${{\bf x}^\ast}$を使っている
   - Morikawa, D. S., & Asai, M. (2021)， ${{\bf u}^\ast}$は使い， ${{\bf x}^\ast}$は使っていない

6. 流速の発散から密度 ${\rho}^\ast$を計算
7. 次の時刻の圧力 $p^{n+1}$を計算
   - ISPHは， $\nabla^2 {p^{n+1}}=(1-\alpha )\frac{\rho_0}{\Delta t}{\nabla \cdot {\bf u}^\ast}+\alpha \frac{\rho_0-\rho^\ast}{{\Delta t}^2}$を解く
   - EISPHは，陽的に $p^{n+1}$を計算する
8. $\nabla {p^{n+1}}$が計算でき， $\frac{D{\bf u}}{D t}=-\frac{1}{\rho_0}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}$（粘性率が一定の非圧縮性流れの加速度）を得る．
9. $\frac{D\bf u}{Dt}$を使って，流速を更新．流速を使って位置を更新

[./builds/build_sph/SPH.hpp#L214](./builds/build_sph/SPH.hpp#L214)

ISPHを使えば，水面粒子の圧力を簡単にゼロにすることができる．
         $\nabla \cdot {\bf u}^*$は流ればで満たされれば十分であり，壁面表層粒子の圧力を，壁面表層粒子上で$\nabla \cdot {\bf u}^*$となるように決める必要はない．

[./builds/build_sph/SPH.hpp#L387](./builds/build_sph/SPH.hpp#L387)


 --- 
## Bucketを用いた粒子探索のテスト
Smoothed Particle Hydrodynamics (SPH)では，効率的な近傍粒子探査が必要となる．
このコードでは，Bucketを用いた粒子探索のテストを行う．

結果はVTKファイルに出力される．
   * 全ての粒子を表示したものは`all.vtp`
   * 中心の粒子を表示したものは`center*.vtp`
   * 中心の粒子が探査したセル内にある粒子を表示したものは`inCell*.vtp`
   * セル内かつ球内にある粒子を表示したものは`inSphere*.vtp`

   - 各セルにある粒子を表示したものは`each_cell*.vtp`
   - 各セルの中心位置を表示したものは`each_cell_position*.vtp`

[./builds/build_sph/test_Buckets.cpp#L1](./builds/build_sph/test_Buckets.cpp#L1)


 --- 
## 壁面粒子の流速と圧力
壁面粒子の流速は常にゼロとすることは自然なこと．常にゼロとするならば，壁面粒子の流速をマップする方法に悩む必要はない．
一方，壁面粒子の圧力は，各ステップ毎に計算し直す必要がある．

壁面粒子の圧力は，壁面法線方向流速をゼロにするように設定されるべきだろう．

[./builds/build_sph/SPH_Functions.hpp#L215](./builds/build_sph/SPH_Functions.hpp#L215)

### `PoissonRHS`と $\nabla^2 p^{n+1}$における $p^{n+1}$の係数の計算
$$
\begin{align*}
\frac{D {\bf u}}{D t} &=-\frac{1}{\rho} \nabla P+\nu \nabla^2 {\bf u}+{\bf g}\\
\rightarrow \nabla \cdot\left(\frac{\rho}{\Delta t} {\bf u}^{n+1}\right) + \nabla^2 p &= \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}+\rho {\bf g}\right)\\
\rightarrow \nabla^2 p &= b, \quad b = \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}+\rho {\bf g}\right)
\end{align*}
$$

ここの $b$を`PoissonRHS`とする．

**✅ CHECKED:** $\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}$

**✅ CHECKED:** $\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}$

**✅ CHECKED:** $\nabla^2 p^{n+1}=\sum_{j}A_{ij}(p_i^{n+1} - p_j^{n+1}),\quad A_{ij} = \frac{2}{\rho_i}m_j\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}$

[./builds/build_sph/SPH_Functions.hpp#L460](./builds/build_sph/SPH_Functions.hpp#L460)

**✅ CHECKED:** $\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}$

**✅ CHECKED:** $\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}$

[./builds/build_sph/SPH_Functions.hpp#L552](./builds/build_sph/SPH_Functions.hpp#L552)


 --- 
## ヘッセ行列を利用したニュートン法
**最適か否かを判断するための関数**は１つだけで，**最適化したい変数は複数**である場合でも，
最適化は，ヘッセ行列を利用したニュートン法によって可能である．
この方法で，変数は，関数を根とするのではなく，関数を最大最小（停留点）とする値へと収束する．

[./builds/build_root_finding/example_NewtonRaphson.cpp#L1](./builds/build_root_finding/example_NewtonRaphson.cpp#L1)


 --- 
## 準ニュートン法
ニュートン法で使うヤコビアンなどを別のものに置き換えた方法．

[./builds/build_root_finding/example_Broyden.cpp#L1](./builds/build_root_finding/example_Broyden.cpp#L1)


 --- 
## Compressed Sparse Row (CSR)
CSRは行列を表現する方法の一つである．
このCSRクラスは，std::unordered_mapを用いて，行列の非ゼロ要素を表現する．
std::unordered_mapのkeyはポインタであり，valueはdoubleである．
CSRクラス自身が，行列の行番号を保存しており，keyであるCSRクラスは行列の列番号を保存している．

[./builds/build_system_of_linear_eqs/CSR.cpp#L1](./builds/build_system_of_linear_eqs/CSR.cpp#L1)


 --- 
## 一般化最小残差法(GMRES)
- ヘッセンベルグ行列$H$
- クリロフ部分空間の直交基底$V$
- $H$をQR分解した行列$Q$と$R$
- $g$は行列$Q$の最初の列

ArnoldiProcessによって，$H$と$V$を求める．このArnoldiProcessクラスの派生クラスとしてGMRESを定義している．

[./builds/build_system_of_linear_eqs/GMRES.cpp#L1](./builds/build_system_of_linear_eqs/GMRES.cpp#L1)


 --- 
