## ArnoldiProcess
ヘッセンベルグ行列$H[0:k-1]$は，Aと相似なベクトルであり，同じ固有値を持つ
   GMRESで使う場合，$V0$にはNormalize(b-A.x0)を与える．
   x0は初期値

   アーノルディ法は固有値問題の数値解法であり反復解法．
   一般的な行列の固有ベクトルと固有値をクリロフ空間の直行基底によって近似する方法計算する方法．
   https://en.wikipedia.org/wiki/Arnoldi_iteration

[./include/basic_linear_systems.hpp#L677](./include/basic_linear_systems.hpp#L677)


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
