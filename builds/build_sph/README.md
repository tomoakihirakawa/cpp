# Contents

- [🐋 Smoothed Particle Hydrodynamics (SPH) ISPH EISPH](#🐋-Smoothed-Particle-Hydrodynamics-(SPH)-ISPH-EISPH)
    - [⛵️ 概要](#⛵️-概要)
        - [⚓️ 前準備](#⚓️-前準備)
        - [⚓️ フラクショナルステップを使って初期値問題を解く](#⚓️-フラクショナルステップを使って初期値問題を解く)
        - [⚓️ 壁面粒子の流速と圧力](#⚓️-壁面粒子の流速と圧力)
        - [⚓️ $`\nabla^2 {\bf u}`$の計算](#⚓️-$`\nabla^2-{\bf-u}`$の計算)
        - [⚓️ `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算](#⚓️-`PoissonRHS`,$`b`$と$`\nabla^2-p^{n+1}`$における$`p^{n+1}`$の係数の計算)
        - [⚓️ 圧力の安定化](#⚓️-圧力の安定化)
        - [⚓️ 圧力勾配$`\nabla p^{n+1}`$の計算 -> $`{D {\bf u}}/{Dt}`$の計算](#⚓️-圧力勾配$`\nabla-p^{n+1}`$の計算-->-$`{D-{\bf-u}}/{Dt}`$の計算)
    - [⛵️ Bucketを用いた粒子探索のテスト](#⛵️-Bucketを用いた粒子探索のテスト)
    - [⛵️ 核関数](#⛵️-核関数)


---
[![Banner](banner.png)](banner.png)

# 🐋 Smoothed Particle Hydrodynamics (SPH) ISPH EISPH

## ⛵️ 概要
### ⚓️ 前準備
1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL条件を満たすようにタイムステップ間隔 $`\Delta t`$を設定

### ⚓️ フラクショナルステップを使って初期値問題を解く

4. $`\nabla^2 {\bf u}`$の計算
5. `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算
6. 流速の発散から密度 $`{\rho}^\ast`$を計算
7. 次の時刻の圧力 $`p^{n+1}`$を計算
1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
2. 流体粒子の圧力$`p^{n+1}`$の計算
8. $`\nabla {p^{n+1}}`$が計算でき， $`\frac{D{\bf u}}{D t}=-\frac{1}{\rho}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}`$（粘性率が一定の非圧縮性流れの加速度）を得る．
9. $`\frac{D\bf u}{Dt}`$を使って，流速を更新．流速を使って位置を更新


[./SPH.hpp#L211](./SPH.hpp#L211)


ISPHを使えば，水面粒子の圧力を簡単にゼロにすることができる．
$`\nabla \cdot {\bf u}^\ast`$は流ればで満たされれば十分であり，壁面表層粒子の圧力を，壁面表層粒子上で$`\nabla \cdot {\bf u}^\ast`$となるように決める必要はない．


[./SPH.hpp#L387](./SPH.hpp#L387)


### ⚓️ 壁面粒子の流速と圧力

壁粒子の流速を流体粒子の流速に応じて変化させると計算が煩雑になるので，**ここでは**壁面粒子の流速は常にゼロに設定することにした（ゼロで一定というのは不自然ではない）．
一方，壁粒子の圧力がゼロだとするのは不自然で，流体粒子の圧力$`p^{n+1}`$の計算に悪影響を及ぼす．
なので．壁粒子の圧力は各ステップ毎に計算し直す必要がある．

壁面粒子の圧力は，壁面法線方向流速をゼロにするように設定されるべきだろう．


[./SPH_Functions.hpp#L216](./SPH_Functions.hpp#L216)


### ⚓️ $`\nabla^2 {\bf u}`$の計算

✅ ラプラシアンの計算方法: $`\nabla^2 {\bf u}=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$


[./SPH_Functions.hpp#L230](./SPH_Functions.hpp#L230)


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


[./SPH_Functions.hpp#L302](./SPH_Functions.hpp#L302)


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


[./SPH_Functions.hpp#L396](./SPH_Functions.hpp#L396)


### ⚓️ 圧力勾配$`\nabla p^{n+1}`$の計算 -> $`{D {\bf u}}/{Dt}`$の計算

✅ 勾配の計算方法: $`\nabla p _i = \rho _i \sum _{j} m _j (\frac{p _i}{\rho _i^2} + \frac{p _j}{\rho _j^2}) \nabla W _{ij}`$

✅ 勾配の計算方法: $`\nabla p _i = \sum _{j} \frac{m _j}{\rho _j} p _j \nabla W _{ij}`$


[./SPH_Functions.hpp#L448](./SPH_Functions.hpp#L448)


## ⛵️ 核関数
3次スプライン関数と5次スプライン関数の実装とテストコード
* 関数の形状を確認．
* 体積積分が1になるかどうかを確認．


[./test_KernelFunctions.cpp#L1](./test_KernelFunctions.cpp#L1)


---
## ⛵️ Bucketを用いた粒子探索のテスト
Smoothed Particle Hydrodynamics (SPH)では，効率的な近傍粒子探査が必要となる．
このコードでは，Bucketを用いた粒子探索のテストを行う．

結果はVTKファイルに出力される．
* 全ての粒子を表示したものは`all.vtp`
* 中心の粒子を表示したものは`center*.vtp`
* 中心の粒子が探査したセル内にある粒子を表示したものは`inCell*.vtp`
* セル内かつ球内にある粒子を表示したものは`inSphere*.vtp`

- 各セルにある粒子を表示したものは`each_cell*.vtp`
- 各セルの中心位置を表示したものは`each_cell_position*.vtp`


[./test_Buckets.cpp#L1](./test_Buckets.cpp#L1)


---