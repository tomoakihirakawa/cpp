# Contents

- [🐋Smoothed Particle Hydrodynamics (SPH) ISPH EISPH](#🐋Smoothed-Particle-Hydrodynamics-(SPH)-ISPH-EISPH)
    - [🪼概要](#🪼概要)
        - [🪸前準備](#🪸前準備)
        - [🪸フラクショナルステップを使って初期値問題を解く](#🪸フラクショナルステップを使って初期値問題を解く)
        - [🪸CFL条件の設定](#🪸CFL条件の設定)
        - [🪸法線方向の計算と水面の判定](#🪸法線方向の計算と水面の判定)
        - [🪸壁面粒子の流速と圧力](#🪸壁面粒子の流速と圧力)
        - [🪸$`\nabla^2 {\bf u} _i`$の計算](#🪸$`\nabla^2-{\bf-u}-_i`$の計算)
        - [🪸圧力の計算　`PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算](#🪸圧力の計算　`PoissonRHS`,$`b`$と$`\nabla^2-p^{n+1}`$における$`p^{n+1}`$の係数の計算)
        - [🪸圧力を決定するための方程式を作成](#🪸圧力を決定するための方程式を作成)
        - [🪸圧力の安定化](#🪸圧力の安定化)
        - [🪸圧力勾配$`\nabla p^{n+1}`$の計算](#🪸圧力勾配$`\nabla-p^{n+1}`$の計算)
    - [🪼注意点](#🪼注意点)
    - [🪼Bucketを用いた粒子探索のテスト](#🪼Bucketを用いた粒子探索のテスト)
    - [🪼核関数](#🪼核関数)


---
[![Banner](banner.png)](banner.png)

# 🐋Smoothed Particle Hydrodynamics (SPH) ISPH EISPH 

## 🪼概要 
### 🪸前準備 
1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL条件を満たすようにタイムステップ間隔 $`\Delta t`$を設定

### 🪸フラクショナルステップを使って初期値問題を解く 

4. 水面の判定
5. $`\nabla^2 {\bf u}`$の計算
6. `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算
7. 流速の発散から密度 $`{\rho}^\ast`$を計算
8. 次の時刻の圧力 $`p^{n+1}`$を計算
1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
2. 流体粒子の圧力$`p^{n+1}`$の計算
9. $`\nabla {p^{n+1}}`$が計算でき， $`\frac{D{\bf u}}{D t}=-\frac{1}{\rho}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}`$（粘性率が一定の非圧縮性流れの加速度）を得る．
10. $`\frac{D\bf u}{Dt}`$を使って，流速を更新．流速を使って位置を更新


[./SPH.hpp#L209](./SPH.hpp#L209)


### 🪸CFL条件の設定 

$`\max({\bf u}) \Delta t \leq c _{v} h \cap \max({\bf a}) \Delta t^2 \leq c _{a} h`$を満たすように，毎時刻$`\Delta t`$を設定する．


[./SPH_Functions.hpp#L22](./SPH_Functions.hpp#L22)


### 🪸法線方向の計算と水面の判定 

✅ 単位法線ベクトル: $`{\bf n} _i = -{\rm Normalize}\left(\sum _j {\frac{m _j}{\rho _j} \nabla W _{ij} }\right)`$


[./SPH_Functions.hpp#L89](./SPH_Functions.hpp#L89)


`surface_condition0,1`の両方を満たす場合，水面とする．


[./SPH_Functions.hpp#L137](./SPH_Functions.hpp#L137)


### 🪸壁面粒子の流速と圧力 

壁粒子の流速を流体粒子の流速に応じて変化させると計算が煩雑になるので，**ここでは**壁面粒子の流速は常にゼロに設定することにした（ゼロで一定というのは不自然ではない）．
一方，壁粒子の圧力がゼロだとするのは不自然で，流体粒子の圧力$`p^{n+1}`$の計算に悪影響を及ぼす．
なので．壁粒子の圧力は各ステップ毎に計算し直す必要がある．

📝 壁面粒子の圧力は，壁面法線方向流速をゼロにするように設定されるべきだろう．


[./SPH_Functions.hpp#L227](./SPH_Functions.hpp#L227)


### 🪸$`\nabla^2 {\bf u} _i`$の計算 

✅ [ラプラシアンの計算方法](../../builds/build_sph/SPH_Functions.hpp#L289): $`\nabla^2 {\bf u} _i=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$

<details>
<summary>見出し部分。ここをクリック。</summary>
<div>
ここが隠れてる部分。
</div>
</details>


[./SPH_Functions.hpp#L241](./SPH_Functions.hpp#L241)


### 🪸圧力の計算　`PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算 

次の時刻の流れ場が発散なし$`\nabla\cdot{\bf u}^{n+1}=0`$であることを保証してくれる圧力を使って，
$`\frac{D {\bf u}}{D t} =-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}`$を決定し，時間発展させたい．
そのような圧力を$`p^{n+1}`$と書くことにする．
そのような圧力の条件は，次のようになる．

$$
\begin{align*}
&&\frac{D {\bf u}}{D t} &=-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}\\
&\rightarrow& \frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t} &=-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}
\end{align*}
$$

非圧縮流体なので，$`\nabla \cdot{\bf u}^{n}`$はゼロであるべきだが，計算誤差が蓄積しゼロからずれてしまう．
そこで次の時刻の$`\nabla \cdot{\bf u}^{n+1}`$をゼロにするように圧力を決定する．

次時刻の発散の演算は，次時刻における粒子配置に基づき行われる．
なので，現在の粒子配置に基づく演算とは区別すべきである．
現在の微分演算を$`\nabla^{n}`$とし，次時刻の微分演算を$`\nabla^{n+1}`$としよう．

$$
\nabla^{n+1}\cdot {\bf u}^{n+1} = \nabla^{n+1} \cdot{\bf u}^{n} - \Delta t \nabla^{n+1} \cdot\left(\frac{1}{\rho} \nabla^{n} p^{n+1}-\nu \nabla^{n2} {\bf u}^n-{\bf g}\right)
$$

次時刻の流速の発散がゼロになるには

$$
\begin{align*}
&&\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= \frac{1}{\Delta t}\nabla^{n+1} \cdot{\bf u}^{n} + \nabla^{n+1} \cdot\left(\nu^n \nabla^{n2} {\bf u}^n  + {\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= \nabla^{n+1} \cdot\left(\frac{1}{\Delta t}{\bf u}^{n} +\nu^n \nabla^{n2} {\bf u}^n  + {\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= b = \nabla^{n+1} \cdot {\bf b}^n,\quad  {\bf b}^n=\frac{1}{\Delta t}{\bf u}^{n} +\nu^n \nabla^{n2} {\bf u}^n
\end{align*}
$$

重力の発散はゼロなので消した．

**右辺について**

この$`b`$を`PoissonRHS`とする．（仮流速は$`{\bf u}^\ast = \frac{\Delta t}{\rho}{\bf b}^n`$と同じ）．`PoissonRHS`,$`b`$の計算の前に，$`\mu \nabla^2{\bf u}`$を予め計算しておく．

✅ [発散の計算方法](not found): $`b=\nabla\cdot{\bf b}^n=\sum _{j}\frac{m _j}{\rho _j}({\bf b} _j^n-{\bf b} _i^n)\cdot\nabla W _{ij}`$

**左辺について**

壁粒子の圧力は時間積分して計算しないので，毎時刻，壁粒子の$`p^n`$を計算する必要がある．

EISPH

1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
2. 流体粒子の圧力$`p^{n+1}`$の計算

ISPH
- ISPHは作ったポアソン方程式を作成し解くことで圧力を計算する

✅ [ラプラシアンの計算方法](../../builds/build_sph/SPH_Functions.hpp#L528): $`\nabla^2 p^{n+1}=\sum _{j}A _{ij}(p _i^{n+1} - p _j^{n+1}),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$


[./SPH_Functions.hpp#L336](./SPH_Functions.hpp#L336)


### 🪸圧力を決定するための方程式を作成 

💡 '次の時刻における流速の発散はゼロになるように'というルールに従えば，次時刻の発散の演算は次時刻の粒子位置において行われるため，今作成するポアソン方程式の発散の演算は，次時刻の粒子位置において行われるべきだ．

各粒子$`A`$に対して，圧力を決定するための方程式を作成する．各粒子$`A`$が，流体か壁か補助粒子か水面かによって，方程式が異なる．

|方程式|目的|
|:---------|---|
| ☑️ [ポアソン方程式](../../builds/build_sph/SPH_Functions.hpp#L515)              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
| ☐ [不透過条件](../../builds/build_sph/SPH_Functions.hpp#L495)         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
| ☐ [大気圧条件](../../builds/build_sph/SPH_Functions.hpp#L503) | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$`h/2`$上なので使わない． |

各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．


[./SPH_Functions.hpp#L477](./SPH_Functions.hpp#L477)


### 🪸圧力の安定化 

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


[./SPH_Functions.hpp#L575](./SPH_Functions.hpp#L575)


### 🪸圧力勾配$`\nabla p^{n+1}`$の計算 

✅ [勾配の計算方法](../../builds/build_sph/SPH_Functions.hpp#L735): $`\nabla p _i = \rho _i \sum _{j} m _j (\frac{p _i}{\rho _i^2} + \frac{p _j}{\rho _j^2}) \nabla W _{ij}`$

✅ [勾配の計算方法](../../builds/build_sph/SPH_Functions.hpp#L737): $`\nabla p _i = \rho _i \sum _{j} m _j \left(p _j - p _i\right) \nabla W _{ij}`$

✅ [勾配の計算方法](../../builds/build_sph/SPH_Functions.hpp#L738): $`\nabla p _i = \sum _{j} \frac{m _j}{\rho _j} p _j \nabla W _{ij}`$


[./SPH_Functions.hpp#L696](./SPH_Functions.hpp#L696)


$`\frac{D{\bf u}^n}{Dt} = - \frac{1}{\rho} \nabla p^{n+1} + \nu \nabla^2 {\bf u}^n + {\bf g}`$が計算できた．


[./SPH_Functions.hpp#L754](./SPH_Functions.hpp#L754)


## 🪼注意点 

⚠️ 計算がうまく行く設定を知るために，次の箇所をチェックする．

- [流体として扱う壁粒子を設定するかどうか](../../builds/build_sph/SPH.hpp#L314)
- [壁粒子の圧力をどのように壁面にマッピングするか](not found)
- [水面粒子の圧力をゼロにするかどうか](not found)
- [密度を更新するかどうか](../../builds/build_sph/SPH_Functions.hpp#L859)
- [圧力の安定化をするかどうか](../../builds/build_sph/SPH_Functions.hpp#L600)
- [ルンゲクッタの段数](../../builds/build_sph/input_generator.py#L143)
- [反射の計算方法](../../builds/build_sph/SPH_Functions.hpp#L799)

壁のwall_as_fluidは繰り返しで計算するのはどうか？


[./SPH_Functions.hpp#L896](./SPH_Functions.hpp#L896)


## 🪼核関数 
3次スプライン関数と5次スプライン関数の実装とテストコード
* 関数の形状を確認．
* 体積積分が1になるかどうかを確認．


[./test_KernelFunctions.cpp#L1](./test_KernelFunctions.cpp#L1)


---
## 🪼Bucketを用いた粒子探索のテスト 
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
