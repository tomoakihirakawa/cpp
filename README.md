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
4x4の行列Aとベクトルbを用いて、Ax=bを解く

[./builds/build_system_of_linear_eqs/CSR.cpp#L1](./builds/build_system_of_linear_eqs/CSR.cpp#L1)


 --- 
## 一般化最小残差 (GMRES)
- ヘッセンベルグ行列$H$
- クリロフ部分空間の直交基底$V$
- $H$をQR分解した行列$Q$と$R$
- $g$は行列$Q$の最初の列

is it ok $Q$?
is it クリ $Q$?
is it クリ$Q$?

- ヘッセンベルグ行列 ![](https://latex.codecogs.com/png.latex?H)
- クリロフ部分空間の直交基底 ![](https://latex.codecogs.com/png.latex?V)
- ![](https://latex.codecogs.com/png.latex?H) をQR分解した行列 ![](https://latex.codecogs.com/png.latex?Q) と ![](https://latex.codecogs.com/png.latex?R)
- ![](https://latex.codecogs.com/png.latex?g) は行列 ![](https://latex.codecogs.com/png.latex?Q) の最初の列

ArnoldiProcessによって，$H$と$V$を求める．このArnoldiProcessクラスの派生クラスとしてGMRESを定義している．

[./builds/build_system_of_linear_eqs/GMRES.cpp#L1](./builds/build_system_of_linear_eqs/GMRES.cpp#L1)


 --- 
