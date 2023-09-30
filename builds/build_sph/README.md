# Contents

- [🐋 Smoothed Particle Hydrodynamics (SPH) ISPH EISPH](#🐋-Smoothed-Particle-Hydrodynamics-(SPH)-ISPH-EISPH)
    - [⛵ 概要](#⛵-概要)
        - [🪼 要素法と粒子法](#🪼-要素法と粒子法)
        - [🪼 SPH](#🪼-SPH)
        - [🪼 このプログラムの目的](#🪼-このプログラムの目的)
        - [🪼 大まかな計算の流れ](#🪼-大まかな計算の流れ)
            - [🐚 Navier-Stokes方程式を解く前の準備](#🐚-Navier-Stokes方程式を解く前の準備)
            - [🐚 Navier-Stokes方程式を解く](#🐚-Navier-Stokes方程式を解く)
        - [🪼 CFL条件の設定](#🪼-CFL条件の設定)
    - [⛵ 計算に利用する壁面粒子だけを抽出](#⛵-計算に利用する壁面粒子だけを抽出)
    - [⛵ 壁面粒子の流速の決定](#⛵-壁面粒子の流速の決定)
        - [🪼 勾配演算子の修正をする行列の計算](#🪼-勾配演算子の修正をする行列の計算)
    - [⛵ 流体の法線方向の計算と水面の判定](#⛵-流体の法線方向の計算と水面の判定)
        - [🪼 流体の法線方向の計算](#🪼-流体の法線方向の計算)
        - [🪼 水面の判定](#🪼-水面の判定)
    - [⛵ 粘性項$`\nabla^2 {\bf u} _i`$の計算（`calcLaplacianU`）](#⛵-粘性項$`\nabla^2-{\bf-u}-_i`$の計算（`calcLaplacianU`）)
    - [⛵ ポアソン方程式$`\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) = b`$](#⛵-ポアソン方程式$`\nabla^{n+1}-\cdot-\left(\frac{1}{\rho^n}-\nabla^{n}-p^{n+1}\right)-=-b`$)
        - [🪼 ポアソン方程式](#🪼-ポアソン方程式)
        - [🪼 右辺，$`b`$，`PoissonRHS`について](#🪼-右辺，$`b`$，`PoissonRHS`について)
        - [🪼 左辺について](#🪼-左辺について)
        - [🪼 水面の計算補助粒子`auxiliaryPoints`](#🪼-水面の計算補助粒子`auxiliaryPoints`)
        - [🪼 次時刻の発散演算，$`\nabla^{n+1} \cdot {\bf b}^n = \sum _j \dfrac{m _j}{\rho _j^{n+1}}({\bf b} _j^n-{\bf b} _i^n)\cdot \nabla W({\bf x} _i^{n+1},{\bf x} _j^{n+1},h)`$](#🪼-次時刻の発散演算，$`\nabla^{n+1}-\cdot-{\bf-b}^n-=-\sum-_j-\dfrac{m-_j}{\rho-_j^{n+1}}({\bf-b}-_j^n-{\bf-b}-_i^n)\cdot-\nabla-W({\bf-x}-_i^{n+1},{\bf-x}-_j^{n+1},h)`$)
    - [⛵ ポアソン方程式の解法](#⛵-ポアソン方程式の解法)
    - [⛵ 圧力勾配$`\nabla p^{n+1}`$の計算](#⛵-圧力勾配$`\nabla-p^{n+1}`$の計算)
    - [⛵ 注意点](#⛵-注意点)
- [🐋 実行方法](#🐋-実行方法)
- [🐋 Bucketを用いた粒子探索のテスト](#🐋-Bucketを用いた粒子探索のテスト)
- [🐋 テスト](#🐋-テスト)
    - [⛵ 核関数のテスト](#⛵-核関数のテスト)


---
# 🐋 Smoothed Particle Hydrodynamics (SPH) ISPH EISPH 

## ⛵ 概要 

### 🪼 要素法と粒子法 

有限要素法や境界要素法など，
要素を利用する計算手法は，
節点の接続に基づき要素を構成し（補間），
微分方程式を離散化して解く．
基本的には，要素が歪になると計算ができない．
また，上手に要素を再構成するのは大変である．

一方の粒子法は，節点間になんら決まった（要素の様な）パターンを要求せず，再構成という概念がない．
きれいに整列した粒子の方が計算精度は高いが，乱れたとしても計算はできる．

### 🪼 SPH 

粒子法には主に２つの種類がある．
一つは，越塚らによって提案されたMoving Particle Semi-implicit (MPS)法であり，
もう一つは，[Gingold and Monaghan (1977)](https://academic.oup.com/mnras/article-lookup/doi/10.1093/mnras/181.3.375)と[Lucy (1977)](http://adsabs.harvard.edu/cgi-bin/bib_query?1977AJ.....82.1013L)によって提案されたSmoothed Particle Hydrodynamics (SPH)法である．
世界的にはSPH法がよく使われている．

SPHの研究者および産業ユーザーから成る[SPHETIC](https://www.spheric-sph.org/sph-projects-and-codes)というコミュニティがある．
それによるとSPHは，1970年代に天体物理学における非軸対称な現象を研究するために開発され，
その工学への応用は1990年代と2000年代初頭に登場した．
過去二十年で、この手法は多くの応用分野で急速に発展しており、
衝突から破壊，水面波のシミュレーション，流体-構造相互作用に至るまで多岐にわたっている．

### 🪼 このプログラムの目的 

このプログラムは，
ISPHとISPHを簡単化したEISPHを実装したものである．
まずは，不安要素が少ないISPHで安定した計算方法を確立し，
その後，EISPHへと移行する．

### 🪼 大まかな計算の流れ 

このSPHでは，非圧縮性流体のナビエ・ストークス方程式を解く．

```math
\frac{D\bf u}{Dt} = -\frac{1}{\rho}\nabla {p} + \nu\nabla^2{\bf u} + {\bf g},\quad  \nu=\frac{\mu}{\rho}
```

#### 🐚 Navier-Stokes方程式を解く前の準備 

1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL条件を満たすようにタイムステップ間隔 $`\Delta t`$を設定
4. 水面の判定

#### 🐚 Navier-Stokes方程式を解く 

5. $`\nabla^2 {\bf u}`$の計算
6. `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算
7. 流速の発散から密度 $`{\rho}^\ast`$を計算
8. 次の時刻の圧力 $`p^{n+1}`$を計算
* 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
* 流体粒子の圧力$`p^{n+1}`$の計算
9. $`\nabla {p^{n+1}}`$が計算でき， $`\frac{D{\bf u}}{D t}=-\frac{1}{\rho}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}`$（粘性率が一定の非圧縮性流れの加速度）を得る．
10. $`\frac{D\bf u}{Dt}`$を使って，流速を更新．流速を使って位置を更新


[./SPH.hpp#L210](./SPH.hpp#L210)


---
### 🪼 CFL条件の設定 

$`\max({\bf u}) \Delta t \leq c _{v} h \cap \max({\bf a}) \Delta t^2 \leq c _{a} h`$
を満たすように，毎時刻$`\Delta t`$を設定する．


[./SPH_Functions.hpp#L90](./SPH_Functions.hpp#L90)


---
## ⛵ 計算に利用する壁面粒子だけを抽出 

流体粒子と同じ影響半径を使ってしまうと，流体粒子が参照できる範囲ギリギリにある壁粒子の法線方向の値が不正確になる．
そのため，流体粒子の影響半径よりも広い半径を使って，`q->interp_normal_original`の法線方向を計算することが，重要である．


[./SPH0_setWall_Freesurface.hpp#L273](./SPH0_setWall_Freesurface.hpp#L273)


## ⛵ 壁面粒子の流速の決定 

壁粒子の流速を流体粒子の流速に応じて変化させるとプログラムが煩雑になるので，
**ここでは**壁面粒子の流速は常にゼロに設定することにする．
壁粒子の圧力は，水が圧縮しないように各ステップ毎に計算し直す必要がある．

**フリースリップ条件の設定**

[フリースリップ条件の設定](not found)

| boolian変数 | 意味 |
|:---:|:---:|
|`isSurface` | 水面かどうか |
|`isCaptured` | 計算で用いるかどうか |


[./SPH0_setWall_Freesurface.hpp#L320](./SPH0_setWall_Freesurface.hpp#L320)


---
### 🪼 勾配演算子の修正をする行列の計算 

[Morikawa et al. (2023)](https://doi.org/10.1016/j.jcpx.2023.100125)で紹介されていた，Randles and Libersky (1996)の勾配演算の精度を改善する行列を計算する．
`grad_corr_M`としている．

```math
\begin{align}
{\bf M} = \left(\sum _j V _j ({\bf x} _j-{\bf x} _i) \otimes \nabla W _{ij}\right)^{-1}
\end{align}
```


[./SPH0_setWall_Freesurface.hpp#L11](./SPH0_setWall_Freesurface.hpp#L11)


## ⛵ 流体の法線方向の計算と水面の判定


[./SPH0_setWall_Freesurface.hpp#L375](./SPH0_setWall_Freesurface.hpp#L375)


### 🪼 流体の法線方向の計算 

✅ [単位法線ベクトル](../../builds/build_sph/SPH0_setWall_Freesurface.hpp#L494): $`{\bf n} _i = {\rm Normalize}\left(-\sum _j {\frac{m _j}{\rho _j} \nabla W _{ij} }\right)`$

単位法線ベクトルは，`interp_normal`としている．


[./SPH0_setWall_Freesurface.hpp#L406](./SPH0_setWall_Freesurface.hpp#L406)


### 🪼 水面の判定 

`surface_condition0,1`の両方を満たす場合，水面とする．


[./SPH0_setWall_Freesurface.hpp#L515](./SPH0_setWall_Freesurface.hpp#L515)


---
## ⛵ 粘性項$`\nabla^2 {\bf u} _i`$の計算（`calcLaplacianU`） 

<img src="icons/SELECTED.png" alt="SELECTED" style="width: 1em;"> [流速のラプラシアンの計算方法](../../builds/build_sph/SPH1_lap_div_U.hpp#L44): $`\nabla^2 {\bf u} _i=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

<img src="icons/SELECTED.png" alt="SELECTED" style="width: 1em;"> [流速の発散の計算方法](../../builds/build_sph/SPH1_lap_div_U.hpp#L40): $`\nabla\cdot{\bf u} _i=\sum _{j}\frac{m _j}{\rho _j}({{\bf u} _j-{\bf u} _i}) \cdot\nabla W _{ij}`$


[./SPH1_lap_div_U.hpp#L7](./SPH1_lap_div_U.hpp#L7)


---
## ⛵ ポアソン方程式$`\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) = b`$ 

### 🪼 ポアソン方程式 

次の時刻の流れ場を発散なし$`\nabla\cdot{\bf u}^{n+1}=0`$としてくれる
$`\frac{D {\bf u}}{D t} =-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}`$を使って，流速と粒子位置を時間発展させたい．
そのためには，圧力$`p^{n+1}`$を適切に決める必要がある．

$`\frac{D {\bf u}}{D t}`$は．$`\frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t}`$と離散化し条件を考えてみる．

```math
\frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t} =-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}
```

次時刻の発散の演算は，次時刻における粒子配置に基づき行われるので，現在の粒子配置に基づく発散演算とは区別すべきである．
現在の微分演算を$`\nabla^{n}`$とし，次時刻の微分演算を$`\nabla^{n+1}`$とする．
$`\nabla^{n+1}`$を上の式に作用させると，

```math
\nabla^{n+1}\cdot {\bf u}^{n+1} = \nabla^{n+1} \cdot{\bf u}^{n} - \Delta t \nabla^{n+1} \cdot\left(\frac{1}{\rho} \nabla^{n} p^{n+1}-\nu \nabla^{n2} {\bf u}^n-{\bf g}\right)
```

右辺がゼロとなれば，次時刻の流速の発散がゼロ，$`\nabla^{n+1}{\bf u}^{n+1}=0`$になる：

```math
\begin{align*}
&&0 &= \nabla^{n+1} \cdot{\bf u}^{n} - \Delta t \nabla^{n+1} \cdot\left(\frac{1}{\rho} \nabla^{n} p^{n+1}-\nu \nabla^{n2} {\bf u}^n-{\bf g}\right)\\
&\rightarrow&\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= \frac{1}{\Delta t}\nabla^{n+1} \cdot{\bf u}^{n} + \nabla^{n+1} \cdot\left(\nu^n \nabla^{n2} {\bf u}^n  + {\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= \nabla^{n+1} \cdot\left(\frac{1}{\Delta t}{\bf u}^{n} +\nu^n \nabla^{n2} {\bf u}^n  + {\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= b = \nabla^{n+1} \cdot {\bf b}^n,\quad  {\bf b}^n=\frac{1}{\Delta t}{\bf u}^{n} +\nu^n \nabla^{n2} {\bf u}^n
\end{align*}
```

重力の発散はゼロなので消した．

### 🪼 右辺，$`b`$，`PoissonRHS`について 

この$`b`$を`PoissonRHS`とする．（仮流速は$`{\bf u}^\ast = \frac{\Delta t}{\rho}{\bf b}^n`$と同じ）．
$`{\bf b}^n`$ （[`Poisson_b_vector`](../../builds/build_sph/SPH1_lap_div_U.hpp#L69)）が計算できるように，$`{\bf u}^n`$と$`\nabla^2 {\bf u}^n`$を計算しておく．

✅ [発散の計算方法](../../builds/build_sph/SPH2_FindPressure.hpp#L199): $`b=\nabla\cdot{\bf b}^n=\sum _{j}\frac{m _j}{\rho _j}({\bf b} _j^n-{\bf b} _i^n)\cdot\nabla W _{ij}`$

### 🪼 左辺について 

壁粒子の圧力は時間積分して計算しないので，毎時刻，壁粒子の$`p^{n+1}`$を計算する必要がある．

✅ [ラプラシアンの計算方法](../../builds/build_sph/SPH2_FindPressure.hpp#L193): $`\nabla^2 p^{n+1}=\sum _{j}A _{ij}(p _i^{n+1} - p _j^{n+1}),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

✅ [ラプラシアンの計算方法](../../builds/build_sph/SPH1_lap_div_U.hpp#L48): $`\nabla^2 p^{n+1}=\sum _{j}A _{ij}(p _i^{n+1} - p _j^{n+1}),\quad A _{ij} = \frac{8 m _j\rho _i}{(\rho _i+\rho _j)^2}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

⚠️ 密度$\rho$が粒子に関わらず一定の場合，上の２式は同じになる．しかし，補助粒子の密度は，他の粒子と異なるので，[２つ目のラプラシアンの計算方法](../../builds/build_sph/SPH1_lap_div_U.hpp#L48)を使うべきだろう．

**ISPH**

- ISPHは作ったポアソン方程式を作成し解くことで圧力を計算する

**EISPH**

1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
2. 流体粒子の圧力$`p^{n+1}`$の計算

[EISPHの圧力の設定方法](../../builds/build_sph/SPH2_FindPressure.hpp#L209)


$\sum _j A _{ij} (p _i^{n+1}-p _j^{n+1}) = b$において，$p _j^{n+1} \approx p _j^{n}$とすると，

```math
p _i^{n+1} = \frac{b + \sum _j A _{ij} p _j^{n}}{\sum _j A _{ij}}
```

となる．

### 🪼 水面の計算補助粒子`auxiliaryPoints` 

水面においては，流速の発散ゼロ$`\nabla^{n+1} {\bf u}^{n+1}=0`$と$`p^{n+1}=0`$が満たされる必要がある．
水面外部には，粒子がないので，求めた水面圧力は，ゼロであっても，圧力勾配は誤差を含み，$`\nabla^{n+1} {\bf u}^{n+1}=0`$は満足されない．
そこで，[水面の計算補助粒子](../../include/Network.hpp#L492)を水面外部に追加し，この点を適切計算することで，$`\nabla^{n+1} {\bf u}^{n+1}=0`$が満足されるように工夫する．


[./SPH2_FindPressure.hpp#L7](./SPH2_FindPressure.hpp#L7)


---
### 🪼 次時刻の発散演算，$`\nabla^{n+1} \cdot {\bf b}^n = \sum _j \dfrac{m _j}{\rho _j^{n+1}}({\bf b} _j^n-{\bf b} _i^n)\cdot \nabla W({\bf x} _i^{n+1},{\bf x} _j^{n+1},h)`$ 

$`\nabla^{n+1}`$の計算には，$`\rho^{n+1}`$, $`{\bf x}^{n+1}= {\bf x}^{n} + {\bf u}^{n+1} \Delta t`$が必要である．

* [次時刻の粒子体積](../../builds/build_sph/SPH_Functions.hpp#L183)
* [次時刻の粒子密度](../../builds/build_sph/SPH_Functions.hpp#L165)
* [次時刻の粒子位置](../../builds/build_sph/SPH_Functions.hpp#L186)


[./SPH2_FindPressure.hpp#L89](./SPH2_FindPressure.hpp#L89)


---
各粒子`ROW`が，流体か壁か補助粒子か水面かによって，方程式が異なる．

|方程式|目的|
|:---------|---|
| ☑️ [ポアソン方程式](../../builds/build_sph/SPH2_FindPressure.hpp#L204)              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
| ☐ [不透過条件](../../builds/build_sph/SPH2_FindPressure.hpp#L151)         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
| ☐ [大気圧条件](../../builds/build_sph/SPH2_FindPressure.hpp#L184) | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない． |

各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．


[./SPH2_FindPressure.hpp#L137](./SPH2_FindPressure.hpp#L137)


---
## ⛵ ポアソン方程式の解法 

ISPHのポアソン方程式を解く場合，[ここではGMRES法](../../builds/build_sph/SPH2_FindPressure.hpp#L502)を使う．


[./SPH2_FindPressure.hpp#L405](./SPH2_FindPressure.hpp#L405)


---
## ⛵ 圧力勾配$`\nabla p^{n+1}`$の計算 

✅ [勾配の計算方法](../../builds/build_sph/SPH3_grad_P.hpp#L116): $`\nabla p _i = \rho _i \sum _{j} m _j (\frac{p _i}{\rho _i^2} + \frac{p _j}{\rho _j^2}) \nabla W _{ij}`$

✅ [勾配の計算方法](../../builds/build_sph/SPH3_grad_P.hpp#L64): $`\nabla p _i = \rho _i \sum _{j} m _j \left(p _j - p _i\right) \nabla W _{ij}`$

✅ [勾配の計算方法](../../builds/build_sph/SPH3_grad_P.hpp#L75): $`\nabla p _i = \sum _{j} \frac{m _j}{\rho _j} p _j \nabla W _{ij}`$


[./SPH3_grad_P.hpp#L11](./SPH3_grad_P.hpp#L11)


$`\dfrac{D{\bf u}^n}{Dt} = - \frac{1}{\rho} \nabla p^{n+1} + \nu \nabla^2 {\bf u}^n + {\bf g}`$
が計算できた．


[./SPH3_grad_P.hpp#L121](./SPH3_grad_P.hpp#L121)


---
## ⛵ 注意点 

⚠️ 計算がうまく行く設定を知るために，次の箇所をチェックする．

**NEW**

- [壁粒子の速度の決定方法](../../builds/build_sph/SPH0_setWall_Freesurface.hpp#L357)
- [Poissonにおいてどのようにbベクトルを使うか](../../builds/build_sph/SPH2_FindPressure.hpp#L127)
- [Poissonにおいてどのようにbベクトルを使うか](../../builds/build_sph/SPH2_FindPressure.hpp#L198)
- どのように[壁粒子のb](not found)/[流体粒子のb](../../builds/build_sph/SPH1_lap_div_U.hpp#L91)を作るか

**壁粒子**

- [壁粒子のラプラシアンの計算方法](../../builds/build_sph/SPH1_lap_div_U.hpp#L68)
- [圧力の計算方法](../../builds/build_sph/SPH2_FindPressure.hpp#L106)
- [どの位置において方程式を立てるか](../../builds/build_sph/SPH2_FindPressure.hpp#L278)
- [流体として扱う壁粒子を設定するかどうか](../../builds/build_sph/SPH0_setWall_Freesurface.hpp#L265)/[視野角に流体粒子が含まない壁粒子は除外する](not found)
- [壁粒子の圧力をどのように壁面にマッピングするか](not found)
- [壁粒子の法線方向ベクトルの計算方法](../../builds/build_sph/SPH0_setWall_Freesurface.hpp#L494)
- [反射の計算方法](../../builds/build_sph/SPH_Functions.hpp#L273)

**水面粒子**

- [水面粒子の圧力をゼロにするかどうか](not found)
- [補助粒子の設定はどうなっているか](../../include/Network.hpp#L492)

**その他**

- [密度を更新するかどうか](../../builds/build_sph/SPH_Functions.hpp#L325)
- [圧力の安定化をするかどうか](not found)
- [ルンゲクッタの段数](../../builds/build_sph/from os.py#L145)


壁のwall_as_fluidは繰り返しで計算するのはどうか？


[./SPH_Functions.hpp#L351](./SPH_Functions.hpp#L351)


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



[./main.cpp#L368](./main.cpp#L368)


---
# 🐋 テスト 

## ⛵ 核関数のテスト 

<!-- Key SPH:kernelFunctions not found -->

プログラムした[3次スプライン関数](../../include/kernelFunctions.hpp#L174)と[5次スプライン関数](../../include/kernelFunctions.hpp#L73)のテストコード

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
