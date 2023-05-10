# Contents

- [ArnoldiProcess](#ArnoldiProcess)

- [Runge-Kutta Integration of ODE](#Runge-Kutta-Integration-of-ODE)

- [核関数](#核関数)

- [ISPHとEISPH](#ISPHとEISPH)

- [Bucketを用いた粒子探索のテスト](#Bucketを用いた粒子探索のテスト)

- [壁面粒子の流速と圧力](#壁面粒子の流速と圧力)

    - [仮流速の発散$\nabla\cdot{\bf u}^\ast$の計算](#仮流速の発散$\nabla\cdot{\bf-u}^\ast$の計算)

    - [ポアソン方程式を解いて，非圧縮性を満たす圧力を計算する](#ポアソン方程式を解いて，非圧縮性を満たす圧力を計算する)

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

[./builds/build_sph/SPH.hpp#L389](./builds/build_sph/SPH.hpp#L389)


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

### 仮流速の発散$\nabla\cdot{\bf u}^\ast$の計算
後に，次時刻の流れ場が非圧縮性を満たすようにポアソン方程式を立てて圧力$p$を計算する．
ポアソン方程式に，ここで計算する仮流速の発散$\nabla\cdot{\bf u}^\ast$を代入する．

[./builds/build_sph/SPH_Functions.hpp#L459](./builds/build_sph/SPH_Functions.hpp#L459)

### ポアソン方程式を解いて，非圧縮性を満たす圧力を計算する
**💡 NOTE:**
ISPHかEISPHに関わらず，圧力をポアソン方程式から計算する場合は，圧力の $\nabla\cdot{\bf u}^\ast$を利用する．

ISPH
壁粒子の ${p}^n$はわかっておく必要はない．

EISPH
壁粒子の ${p}^n$がわかっておく必要がある．→　壁に鏡写しすることで，壁粒子の ${p}^n$を計算する．

 - [x] $\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}$
 - [x] $\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}$

[./builds/build_sph/SPH_Functions.hpp#L545](./builds/build_sph/SPH_Functions.hpp#L545)

- [x] $\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}$
 - [x] $\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}$

[./builds/build_sph/SPH_Functions.hpp#L621](./builds/build_sph/SPH_Functions.hpp#L621)

- [x] $\nabla^2 p^{n+1} = \frac{2}{\rho_i} \sum_{j} m_j (p_i^{n+1} - p_j^{n+1}) \frac{{{\bf x}_{ij}}\cdot \nabla W_{ij}}{{\bf x}_{ij}}$

[./builds/build_sph/SPH_Functions.hpp#L668](./builds/build_sph/SPH_Functions.hpp#L668)


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
