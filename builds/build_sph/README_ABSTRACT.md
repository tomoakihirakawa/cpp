## 概要

### 要素法と粒子法

有限要素法や境界要素法など，
要素を利用する計算手法は，
節点の接続に基づき要素を構成し（補間），
微分方程式を離散化して解く．
基本的には，要素が歪になると計算ができない．
また，上手に要素を再構成するのは大変である．

一方の粒子法は，節点間になんら決まった（要素の様な）パターンを要求せず，再構成という概念がない．
きれいに整列した粒子の方が計算精度は高いが，乱れたとしても計算はできる．

### SPH

粒子法には主に２つの種類がある．
一つは，越塚らによって提案されたMoving Particle Semi-implicit (MPS)法であり，
もう一つは，\cite{Gingold1977}と\cite{Lucy1977}によって提案されたSmoothed Particle Hydrodynamics (SPH)法である．
世界的にはSPH法がよく使われている．

SPHの研究者および産業ユーザーから成る[SPHETIC](https://www.spheric-sph.org/sph-projects-and-codes)というコミュニティがある．
それによるとSPHは，1970年代に天体物理学における非軸対称な現象を研究するために開発され，
その工学への応用は1990年代と2000年代初頭に登場した．
過去二十年で、この手法は多くの応用分野で急速に発展しており、
衝突から破壊，水面波のシミュレーション，流体-構造相互作用に至るまで多岐にわたっている．

### このプログラムの目的

このプログラムは，
ISPHとISPHを簡単化したEISPHを実装したものである．
まずは，不安要素が少ないISPHで安定した計算方法を確立し，
その後，EISPHへと移行する．

### 大まかな計算の流れ

このSPHでは，非圧縮性流体のナビエ・ストークス方程式を解く．

```math
\frac{D\bf u}{Dt} = -\frac{1}{\rho}\nabla {p} + \nu\nabla^2{\bf u} + {\bf g},\quad  \nu=\frac{\mu}{\rho}
```

#### Navier-Stokes方程式を解く前の準備

1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL条件を満たすようにタイムステップ間隔 $`\Delta t`$を設定
4. 水面の判定

#### Navier-Stokes方程式を解く

5. $`\nabla^2 {\bf u}`$の計算
6. `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算
7. 流速の発散から密度 $`{\rho}^\ast`$を計算
8. 次の時刻の圧力 $`p^{n+1}`$を計算
   * 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
   * 流体粒子の圧力$`p^{n+1}`$の計算
9. $`\nabla {p^{n+1}}`$が計算でき， $`\frac{D{\bf u}}{D t}=-\frac{1}{\rho}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}`$（粘性率が一定の非圧縮性流れの加速度）を得る．
10. $`\frac{D\bf u}{Dt}`$を使って，流速を更新．流速を使って位置を更新