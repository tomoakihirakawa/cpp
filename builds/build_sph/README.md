# Contents

- [🐋Smoothed Particle Hydrodynamics (SPH) ISPH EISPH](#🐋Smoothed-Particle-Hydrodynamics-(SPH)-ISPH-EISPH)
    - [⛵️概要](#⛵️概要)
        - [🪸前準備](#🪸前準備)
        - [🪸フラクショナルステップを使って初期値問題を解く](#🪸フラクショナルステップを使って初期値問題を解く)
        - [🪸CFL条件の設定](#🪸CFL条件の設定)
        - [🪸法線方向の計算と水面の判定](#🪸法線方向の計算と水面の判定)
        - [🪸壁面粒子の流速と圧力](#🪸壁面粒子の流速と圧力)
    - [⛵️$`\nabla^2 {\bf u} _i`$の計算](#⛵️$`\nabla^2-{\bf-u}-_i`$の計算)
        - [🪸高速化のための工夫](#🪸高速化のための工夫)
    - [⛵️ポアソン方程式$`\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) = b`$](#⛵️ポアソン方程式$`\nabla^{n+1}-\cdot-\left(\frac{1}{\rho^n}-\nabla^{n}-p^{n+1}\right)-=-b`$)
        - [🪸ポアソン方程式](#🪸ポアソン方程式)
        - [🪸右辺，$`b`$，`PoissonRHS`について](#🪸右辺，$`b`$，`PoissonRHS`について)
        - [🪸左辺について](#🪸左辺について)
        - [🪸ポアソン方程式の作成のコーディング](#🪸ポアソン方程式の作成のコーディング)
    - [⛵️ポアソン方程式の解法](#⛵️ポアソン方程式の解法)
    - [⛵️圧力勾配$\nabla p^{n+1}$の計算](#⛵️圧力勾配$\nabla-p^{n+1}$の計算)
    - [⛵️注意点](#⛵️注意点)
    - [⛵️Bucketを用いた粒子探索のテスト](#⛵️Bucketを用いた粒子探索のテスト)
    - [⛵️核関数](#⛵️核関数)


---
# 🐋Smoothed Particle Hydrodynamics (SPH) ISPH EISPH 

## ⛵️概要 
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
* 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
* 流体粒子の圧力$`p^{n+1}`$の計算
9. $`\nabla {p^{n+1}}`$が計算でき， $`\frac{D{\bf u}}{D t}=-\frac{1}{\rho}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}`$（粘性率が一定の非圧縮性流れの加速度）を得る．
10. $`\frac{D\bf u}{Dt}`$を使って，流速を更新．流速を使って位置を更新


[./SPH.hpp#L210](./SPH.hpp#L210)


### 🪸CFL条件の設定 

$\max({\bf u}) \Delta t \leq c _{v} h \cap \max({\bf a}) \Delta t^2 \leq c _{a} h$を満たすように，毎時刻$\Delta t$を設定する．


[./SPH_Functions.hpp#L22](./SPH_Functions.hpp#L22)


### 🪸法線方向の計算と水面の判定 

✅ 単位法線ベクトル: ${\bf n} _i = -{\rm Normalize}\left(\sum _j {\frac{m _j}{\rho _j} \nabla W _{ij} }\right)$


[./SPH_Functions.hpp#L89](./SPH_Functions.hpp#L89)


`surface_condition0,1`の両方を満たす場合，水面とする．


[./SPH_Functions.hpp#L137](./SPH_Functions.hpp#L137)


### 🪸壁面粒子の流速と圧力 

壁粒子の流速を流体粒子の流速に応じて変化させるとプログラムが煩雑になるので，
**ここでは**壁面粒子の流速は常にゼロに設定することにする．

壁粒子の圧力は，水が圧縮しないように各ステップ毎に計算し直す必要がある．


[./SPH_Functions.hpp#L227](./SPH_Functions.hpp#L227)


## ⛵️$`\nabla^2 {\bf u} _i`$の計算 

✅ [ラプラシアンの計算方法](../../builds/build_sph/SPH_Functions.hpp#L293): $`\nabla^2 {\bf u} _i=\sum _{j} A _{ij}({\bf u} _i - {\bf u} _j),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$


[./SPH_Functions.hpp#L241](./SPH_Functions.hpp#L241)


### 🪸高速化のための工夫 

何度か行う勾配の計算は，変数は違えど，変数の係数は同じである．
ここで，その係数を`std::unordered_map`で保存しておくことにする．
`A->grad_coeff`と`A->grad_coeff_next`に保存する．


[./SPH_Functions.hpp#L264](./SPH_Functions.hpp#L264)


## ⛵️ポアソン方程式$`\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) = b`$ 

### 🪸ポアソン方程式 

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

次時刻の流速の発散がゼロ，$`\nabla^{n+1}{\bf u}^{n+1}=0`$になるには

```math
\begin{align*}
&&0 &= \nabla^{n+1} \cdot{\bf u}^{n} - \Delta t \nabla^{n+1} \cdot\left(\frac{1}{\rho} \nabla^{n} p^{n+1}-\nu \nabla^{n2} {\bf u}^n-{\bf g}\right)\\
&\rightarrow&\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= \frac{1}{\Delta t}\nabla^{n+1} \cdot{\bf u}^{n} + \nabla^{n+1} \cdot\left(\nu^n \nabla^{n2} {\bf u}^n  + {\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= \nabla^{n+1} \cdot\left(\frac{1}{\Delta t}{\bf u}^{n} +\nu^n \nabla^{n2} {\bf u}^n  + {\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= b = \nabla^{n+1} \cdot {\bf b}^n,\quad  {\bf b}^n=\frac{1}{\Delta t}{\bf u}^{n} +\nu^n \nabla^{n2} {\bf u}^n
\end{align*}
```

重力の発散はゼロなので消した．

### 🪸右辺，$`b`$，`PoissonRHS`について 

この$`b`$を`PoissonRHS`とする．（仮流速は$`{\bf u}^\ast = \frac{\Delta t}{\rho}{\bf b}^n`$と同じ）．
`PoissonRHS`,$`b`$の計算の前に，$`{\bf b}^n`$を予め計算しておく．

✅ [発散の計算方法](../../builds/build_sph/SPH_Functions.hpp#L540): $`b=\nabla\cdot{\bf b}^n=\sum _{j}\frac{m _j}{\rho _j}({\bf b} _j^n-{\bf b} _i^n)\cdot\nabla W _{ij}`$

### 🪸左辺について 

壁粒子の圧力は時間積分して計算しないので，毎時刻，壁粒子の$`p^{n+1}`$を計算する必要がある．

**EISPH**

1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
2. 流体粒子の圧力$`p^{n+1}`$の計算

**ISPH**

- ISPHは作ったポアソン方程式を作成し解くことで圧力を計算する

✅ [ラプラシアンの計算方法](../../builds/build_sph/SPH_Functions.hpp#L552): $`\nabla^2 p^{n+1}=\sum _{j}A _{ij}(p _i^{n+1} - p _j^{n+1}),\quad A _{ij} = \frac{2m _j}{\rho _i}\frac{{{\bf x} _{ij}}\cdot\nabla W _{ij}}{{\bf x} _{ij}^2}`$


[./SPH_Functions.hpp#L340](./SPH_Functions.hpp#L340)


### 🪸ポアソン方程式の作成のコーディング 

各粒子`A`に対して，方程式を作成する．

まずは，[方程式を立てる位置を決める．](../../builds/build_sph/SPH_Functions.hpp#L475)


[./SPH_Functions.hpp#L466](./SPH_Functions.hpp#L466)


各粒子`A`が，流体か壁か補助粒子か水面かによって，方程式が異なる．

|方程式|目的|
|:---------|---|
| ☑️ [ポアソン方程式](../../builds/build_sph/SPH_Functions.hpp#L537)              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
| ☐ [不透過条件](../../builds/build_sph/SPH_Functions.hpp#L518)         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
| ☐ [大気圧条件](../../builds/build_sph/SPH_Functions.hpp#L526) | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない． |

各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．


[./SPH_Functions.hpp#L504](./SPH_Functions.hpp#L504)


## ⛵️ポアソン方程式の解法 

ISPHのポアソン方程式を解く場合，[ここではGMRES法](../../builds/build_sph/SPH_Functions.hpp#L723)を使う．


[./SPH_Functions.hpp#L651](./SPH_Functions.hpp#L651)


## ⛵️圧力勾配$\nabla p^{n+1}$の計算 

✅ [勾配の計算方法](../../builds/build_sph/SPH_Functions.hpp#L783): $\nabla p _i = \rho _i \sum _{j} m _j (\frac{p _i}{\rho _i^2} + \frac{p _j}{\rho _j^2}) \nabla W _{ij}$

✅ [勾配の計算方法](../../builds/build_sph/SPH_Functions.hpp#L785): $\nabla p _i = \rho _i \sum _{j} m _j \left(p _j - p _i\right) \nabla W _{ij}$

✅ [勾配の計算方法](../../builds/build_sph/SPH_Functions.hpp#L786): $\nabla p _i = \sum _{j} \frac{m _j}{\rho _j} p _j \nabla W _{ij}$


[./SPH_Functions.hpp#L740](./SPH_Functions.hpp#L740)


$`\dfrac{D{\bf u}^n}{Dt} = - \frac{1}{\rho} \nabla p^{n+1} + \nu \nabla^2 {\bf u}^n + {\bf g}`$
が計算できた．


[./SPH_Functions.hpp#L802](./SPH_Functions.hpp#L802)


## ⛵️注意点 

⚠️ 計算がうまく行く設定を知るために，次の箇所をチェックする．

- [流体として扱う壁粒子を設定するかどうか](../../builds/build_sph/SPH.hpp#L314)
- [壁粒子の圧力をどのように壁面にマッピングするか](not found)
- [水面粒子の圧力をゼロにするかどうか](not found)
- [密度を更新するかどうか](../../builds/build_sph/SPH_Functions.hpp#L910)
- [圧力の安定化をするかどうか](../../builds/build_sph/SPH_Functions.hpp#L624)
- [ルンゲクッタの段数](../../builds/build_sph/input_generator.py#L143)
- [反射の計算方法](../../builds/build_sph/SPH_Functions.hpp#L850)

壁のwall_as_fluidは繰り返しで計算するのはどうか？


[./SPH_Functions.hpp#L947](./SPH_Functions.hpp#L947)


## ⛵️核関数 

プログラムした[3次スプライン関数](../../include/kernelFunctions.hpp#L122)と[5次スプライン関数](../../include/kernelFunctions.hpp#L73)のテストコード

* 関数の形状を確認．
* 体積積分が1になるかどうかを確認．

| 分割数$N$，体積$V=(\frac{2r}{N})^3$   | Sum for 3rd Order | Sum for 5th Order |
| --- | ----------------- | ----------------- |
| 5   | 1.00527           | 0.999206          |
| 10  | 1.00011           | 1.00005           |
| 15  | 0.999972          | 0.999999          |
| 20  | 1                 | 1                 |
| 25  | 1                 | 1                 |


[./test_KernelFunctions.cpp#L1](./test_KernelFunctions.cpp#L1)


---
## ⛵️Bucketを用いた粒子探索のテスト 
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
