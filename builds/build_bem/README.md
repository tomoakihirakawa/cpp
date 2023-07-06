# Contents

- [🐋BEM-MEL](#🐋BEM-MEL)
    - [⛵️流速の計算方法](#⛵️流速の計算方法)
        - [🪸修正流速（これがないと激しい波の計算は難しい）](#🪸修正流速（これがないと激しい波の計算は難しい）)
        - [🪸エネルギー保存則（計算精度のチェックに利用できる）](#🪸エネルギー保存則（計算精度のチェックに利用できる）)
        - [🪸内部流速の計算方法（使わなくてもいい）](#🪸内部流速の計算方法（使わなくてもいい）)
    - [⛵️境界のタイプを決定する](#⛵️境界のタイプを決定する)
        - [🪸多重節点](#🪸多重節点)
    - [⛵️境界値問題](#⛵️境界値問題)
        - [🪸基礎方程式](#🪸基礎方程式)
        - [🪸境界積分方程式（BIE）](#🪸境界積分方程式（BIE）)
        - [🪸BIEの離散化](#🪸BIEの離散化)
        - [🪸リジッドモードテクニック](#🪸リジッドモードテクニック)
    - [⛵️浮体動揺解析](#⛵️浮体動揺解析)
        - [🪸ノイマン境界面における$`\phi _{nt}`$の求め方](#🪸ノイマン境界面における$`\phi-_{nt}`$の求め方)
    - [⛵️造波装置など](#⛵️造波装置など)
    - [⛵️初期値問題](#⛵️初期値問題)
        - [🪸流速$`\frac{d\bf x}{dt}`$の計算](#🪸流速$`\frac{d\bf-x}{dt}`$の計算)
        - [🪸$`\frac{d\phi}{dt}`$の計算](#🪸$`\frac{d\phi}{dt}`$の計算)
    - [⛵️その他](#⛵️その他)
        - [🪸境界値問題の未知変数](#🪸境界値問題の未知変数)
        - [🪸$`\phi _{nt}`$の計算で必要となる$`{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right)`$について．](#🪸$`\phi-_{nt}`$の計算で必要となる$`{\bf-n}\cdot-\left({\nabla-\phi-\cdot-\nabla\nabla-\phi}\right)`$について．)
- [🐋入力ファイル生成 `input_generator.py`](#🐋入力ファイル生成-`input_generator.py`)
    - [⛵️Usage](#⛵️Usage)
    - [⛵️Customization](#⛵️Customization)
    - [⛵️Output](#⛵️Output)
- [🐋コンパイルと計算の実行方法](#🐋コンパイルと計算の実行方法)
    - [⛵️Prerequisites](#⛵️Prerequisites)
    - [⛵️Building the Code](#⛵️Building-the-Code)
    - [⛵️Running the Simulation](#⛵️Running-the-Simulation)
    - [⛵️Output](#⛵️Output)
        - [🪸計算の流れ](#🪸計算の流れ)


---
# 🐋BEM-MEL 

| 項目 | 詳細|
|---:|:---|
| 要素 | 線形三角要素 |
| 時間発展方法 | 4次のルンゲクッタ |
| 解析領域 | 時間領域 |
| 境界条件 | 水面の境界条件は非線形であるが，非線形のまま解く |


[./BEM.hpp#L1](./BEM.hpp#L1)


## ⛵️流速の計算方法


[./BEM_calculateVelocities.hpp#L7](./BEM_calculateVelocities.hpp#L7)


### 🪸修正流速（これがないと激しい波の計算は難しい） 

求めた流速から，次の時刻の境界面$`\Omega(t+\Delta t)`$を見積もり，その面上で節点を移動させ歪さを解消する．
修正ベクトルは，$`\Delta t`$で割り，求めた流速$`\nabla \phi`$に足し合わせて，節点を時間発展させる．

ノイマン節点も修正流速を加え時間発展させる．
ただし，ノイマン節点の修正流速に対しては，節点が水槽の角から離れないように，工夫を施している．

`calculateVecToSurface`で$`\Omega(t+\Delta t)`$上へのベクトルを計算する．
まず，`vectorTangentialShift2`で接線方向にシフトし，`vectorToNextSurface`で近の$`\Omega(t+\Delta t)`$上へのベクトルを計算する．


[./BEM_calculateVelocities.hpp#L332](./BEM_calculateVelocities.hpp#L332)


### 🪸エネルギー保存則（計算精度のチェックに利用できる） 

流体全体の運動エネルギーは，ラプラス方程式と発散定理を使うと，次のように境界面に沿った積分で表される．

```math
E _K =\frac{\rho}{2} \iint _\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
```

また，流体の位置エネルギーは，次のように表される．

```math
E _P = \frac{\rho}{2} \iint _\Gamma (0,0,g(z - z _0)^2) \cdot {\bf n} d\Gamma
```

<details>

---

<summary>
💡 なぜか？
</summary>

テンソルを使って考えてみると

```math
\begin{align*}
\nabla \cdot (\phi\nabla\phi) &= \frac{\partial\phi}{\partial x _i} \frac{\partial\phi}{\partial x _i} + \phi \frac{\partial^2\phi}{\partial x _i \partial x _i}\\
&= \nabla \phi \cdot \nabla \phi + \phi \nabla^2 \phi\\
&= \nabla \phi \cdot \nabla \phi
\end{align*}
```

よって，

```math
\iiint _\Omega \nabla\phi\cdot\nabla\phi d\Omega = \iiint _\Omega \nabla \cdot (\phi\nabla\phi) d\Omega = \iint _\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
```

---

```math
E _P = \rho g \iiint _\Omega (z - z _0) d\Omega
= \rho g \iiint _\Omega \frac{1}{2} \nabla \cdot (0,0,(z - z _0)^2) d\Omega
= \rho g \iint _\Gamma \frac{1}{2} (0,0,(z - z _0)^2) \cdot {\bf n} d\Gamma
= \frac{1}{2}\rho g \iint _\Gamma (z - z _0)^2 n _z d\Gamma
```

---

</details>


[./BEM_calculateVelocities.hpp#L476](./BEM_calculateVelocities.hpp#L476)


### 🪸内部流速の計算方法（使わなくてもいい） 

[Fochesato2005](https://onlinelibrary.wiley.com/doi/10.1002/fld.838)にあるように，
流体内部の流速$`\nabla \phi`$は，BIEを微分して求めることができる．

```math
u({\bf a}) = \nabla\phi({\bf a}) = \int _{\partial \Omega} \frac{\partial Q}{\partial n} ({\bf x})Q({\bf x}, {\bf a}) - \phi({\bf x}) \frac{\partial Q}{\partial n} ({\bf x}, {\bf a}) d\Gamma
```

```math
Q({\bf x},{\bf a}) = \frac{{\bf r}}{4\pi r^3}, \quad \frac{\partial Q}{\partial n} ({\bf x},{\bf a}) = \frac{1}{4\pi r^3} (3 \mathbf{n} - (\mathbf{r} \cdot \mathbf{n}) \frac{\mathbf{r}}{r^2})
```


[./BEM_calculateVelocities.hpp#L563](./BEM_calculateVelocities.hpp#L563)


## ⛵️境界のタイプを決定する 

0. 流体と物体の衝突を判定し，流体節点が接触する物体面を保存しておく．
[`networkPoint::contact_angle`](../../include/networkPoint.hpp#L171)，
[`networkPoint::isInContact`](../../include/networkPoint.hpp#L178)，
[`networkPoint::addContactFaces`](../../include/networkPoint.hpp#L289)
を使って接触判定を行っている．

つぎに，その情報を使って，境界のタイプを次の順で決める．（物理量を与えるわけではない）

1. 面の境界条件：３節点全てが接触している流体面はNeumann面，それ以外はDirichlet面とする．CORNER面は設定しない．
- Neumann面$`\Gamma^{({\rm N})}`$ : 3点接触流体面
- Dirichlet面$`\Gamma^{({\rm D})}`$ : それ以外の面

2. 辺の境界条件 : 辺を含む２面がNeumann面ならNeumann辺，２面がDirichlet面ならDirichlet辺，それ以外はCORNERとする．
- Neumann辺 : 隣接面2面がNeumann面の辺
- Dirichlet辺 : 隣接面2面がDirichlet面の辺
- CORNER辺 : それ以外の辺（Neumann面とDirichlet面の間にある辺）

3. 点の境界条件：点を含む面全てがNeumann面ならNeumann点，面全てがDirichlet面ならDirichlet点，それ以外はCORNERとする．
- Neumann点 : 隣接面全てがNeumann面である点
- Dirichlet点 : 隣接面全てがDirichlet面である点
- CORNER点 : それ以外の点（Neumann面とDirichlet面の間にある点）

### 🪸多重節点 

💡 面の向き$`\bf n`$がカクッと不連続に変わる節点には，$`\phi`$は同じでも，隣接面にそれぞれ対して異なる$`\phi _n`$を計算できるようにする

💡 $`\bf n`$が不連続に変化する節点まわりの要素は，自分のために用意された$`\phi _n`$を選択し補間に用いなければならない

これを多重節点という．


[./BEM_setBoundaryTypes.hpp#L7](./BEM_setBoundaryTypes.hpp#L7)


## ⛵️境界値問題 

### 🪸基礎方程式 

```math
\begin{align}
\nabla\cdot\nabla \phi& = 0&&\text{in}&&{\bf x} \in \Omega(t),\\
\frac{\partial\phi}{\partial t} +\frac{1}{2}\nabla\phi\cdot\nabla\phi - g z &=0 &&\text{on}&&{\bf x} \in \Gamma^{(\rm D)}(t),\\
\phi _n + {{\bf u} _b}\cdot{{\bf n} _b} &=0&&\text{on}&&{\bf x}\in \Gamma^{(\rm N)}(t),
\end{align}
```

ここで，
$`{\bf x} ={(x,y,z)}`$は空間座標，$`{\bf u} _b`$は物体の流速，
$`{\bf n} _b`$は物体の外向き単位法線ベクトル，
$`\nabla=(\frac{\partial}{\partial x},\frac{\partial}{\partial y},\frac{\partial}{\partial z})`$
である．
また，$`\phi _n`$は境界面上での外向き法線方向の流速を表し，
境界面上の外向き単位法線ベクトル$`\bf n`$を使えば$`\phi _n ={\nabla\phi}\cdot {\bf n}`$で表される．

### 🪸境界積分方程式（BIE） 

**グリーンの定理**

任意の$`\phi`$，$`G`$に対して次が成り立つ（**グリーンの定理**）．

```math
\iiint _\Omega \left(G({\bf x},{\bf a})\nabla^2 \phi({\bf x}) - \phi({\bf x})\nabla^2 G({\bf x},{\bf a})\right)dV
= \iint _\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
```


$`\phi`$がラプラス方程式$`\nabla^2\phi=0`$を満たし，$`G=1/\|{\bf x}-{\bf a}\|`$とすると，
グリーンの定理から$`\phi`$と$`\phi _n`$の関係式，BIEが得られる．

```math
\alpha ({\bf{a}})\phi ({\bf{a}}) = \iint _\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
```

ここで，$`{\bf a}`$は境界面上の位置ベクトルであり，この原点$`{\bf a}`$を固定し$`{\bf x}`$について面積分される．
$`G`$は任意のスカラー関数で$`G=1/\|{\bf x}-{\bf a}\|`$とすることで，グリーンの定理の体積積分が消え，BIEの左辺のように，
原点での立体角$`\alpha\left( {\bf{a}} \right)`$とポテンシャル$`\phi( {\bf{a}})`$の積だけが残る．

この式は，流体内部では，$`\alpha ({\bf{a}})`$は$`1`$とできる．
この式は，$`\bf{a}`$におけるポテンシャル$`\phi ({\bf{a}})`$が，右辺の１重層ポテンシャルと２重層ポテンシャルの和で表されることを示している．
$`G=1/\|{\bf x}-{\bf a}\|`$がラプラス法廷式の基本解であり，$`\phi`$は境界におけるポテンシャルの分布である．


[./BEM_solveBVP.hpp#L7](./BEM_solveBVP.hpp#L7)


### 🪸BIEの離散化 

BIEを線形三角要素とGauss-Legendre積分で離散化すると，

```math
\sum\limits _{k _\vartriangle}\sum\limits _{{\xi _1},{w _1}} {\sum\limits _{{\xi _0},{w _0}} {\left( {{w _0}{w _1}\left( {\sum\limits _{j=0}^2 {{{\left( {{\phi _n}} \right)} _{k _\vartriangle,j }}{N _{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x} _{i _\circ}}} \|}}\left\|\frac{{\partial{\bf{x}}}}{{\partial{\xi _0}}} \times \frac{{\partial{\bf{x}}}}{{\partial{\xi _1}}}\right\|} \right)} }=
```
```math
\alpha _{i _\circ}(\phi) _{i _\circ}-\sum\limits _{k _\vartriangle}\sum\limits _{{\xi _1},{w _1}} \sum\limits _{{\xi _0},{w _0}} {\left( {{w _0}{w _1}\left({\sum\limits _{j =0}^2{{{\left( \phi  \right)} _{k _\vartriangle,j }}{N _{j}}\left( \pmb{\xi } \right)} } \right)\frac{\bf{x}(\pmb{\xi})-{{\bf x} _{i _\circ} }}{{{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x} _{i _\circ}}}\|}^3}}} \cdot\left(\frac{{\partial {\bf{x}}}}{{\partial {\xi _0}}}\times\frac{{\partial {\bf{x}}}}{{\partial {\xi _1}}}\right)}\right)}
```

ここで，$`\phi _{k _\vartriangle,j}`$における$`k _\vartriangle`$は三角形要素の番号，$`j`$は三角形要素の頂点番号．
$`N _j`$は三角形要素の形状関数，$`\pmb{\xi}`$は三角形要素の内部座標，$`w _0,w _1`$はGauss-Legendre積分の重み，$`\alpha _{i _\circ}`$は原点$`i _\circ`$における立体角，$`\phi`$はポテンシャル，$`\phi _n`$は法線方向のポテンシャル，$`\bf{x}`$は空間座標，$`{\bf x} _{i _\circ}`$は原点の空間座標である．

形状関数$`{\pmb N} _j({\pmb \xi}),{\pmb \xi}=(\xi _0,\xi _1)`$は，$`\xi _0,\xi _1`$が$`0`$から$`1`$動くことで，範囲で三角要素全体を動くように定義している．

```math
{\pmb N}({\pmb \xi}) = (N _0({\pmb \xi}),N _1({\pmb \xi}),N _2({\pmb \xi})) = (\xi _0, - \xi _1 (\xi _0 - 1), (\xi _0-1)(\xi _1-1))
```


[./BEM_solveBVP.hpp#L218](./BEM_solveBVP.hpp#L218)


このループでは，BIEの連立一次方程式の係数行列`IGIGn`を作成する作業を行なっている．
`IGIGn`は，ある節点$`i _\circ`$（係数行列の行インデックス）に対する
他の節点$`j _\circ`$（係数行列の列インデックス）の影響度合いのようなものである．
その影響度合いは，他の節点$`j _\circ`$の所属する要素までの距離や向きによって決まることが離散化された式からわかる．

| Variable | Description |
|:--------:|:-----------:|
| `origin` | 原点となる節点$`i _\circ`$ |
| `integ_f` | Element $`k _{\triangle}`$ |
| `t0, t1, ww` | Gaussian points and thier wieghts $`\xi _0, \xi _1, w _0 w _1`$ |
| `p0, p1, p2` | Node of the element $`k _{\triangle}`$ |
| `N012` | Shape function $`\pmb{N} _j`$ |
| `IGIGn` | Coefficient matrices of the left and right sides |
| `nr` | $`\| \pmb{x} - \pmb{x} _{i\circ } \|`$ |
| `tmp` | $`w _0 w _1 \frac{1 - \xi _0}{\| \pmb{x} - \pmb{x} _{i\circ } \|}`$ |
| `cross` | $`\frac{\partial \pmb{x}}{\partial \xi _0} \times \frac{\partial \pmb{x}}{\partial \xi _1}`$ |


[./BEM_solveBVP.hpp#L282](./BEM_solveBVP.hpp#L282)


### 🪸リジッドモードテクニック 

全て$`\phi=1`$とすると，$`\alpha({\bf a}) = -\int\int{\nabla G({\bf x},{\bf a})\cdot{\bf n}({\bf x})dS}`$となり，これを離散化すると，数値積分による評価が難しかった係数行列の対角成分がより精確に計算できる．
これはリジッドモードテクニックと呼ばれている．
$`{\bf x} _{i\circ}`$が$`{\bf x}({\pmb \xi})`$に近い場合，$`G`$は急激に特異的に変化するため，数値積分精度が悪化するが，リジッドモードテクニックによって積分を回避できる．


[./BEM_solveBVP.hpp#L355](./BEM_solveBVP.hpp#L355)


係数行列`IGIGn`は，左辺の$`I _G \phi _n`$，右辺の$`I _{G _n}\phi`$の係数．

```math
(I _G) _{i _\circ,j _\circ} (\phi _n) _{j _\circ} = (I _{Gn}) _{i _\circ,j _\circ}  \phi _{j _\circ}
```

境界条件に応じて，未知変数は$`\phi,\phi _n`$のどちらかに決まる．
未知変数が$`\phi`$の場合（Dirichlet境界条件の場合），
係数行列`IGIGn`中で対応する列を符号変えて入れ替えることで移項したことになる．


移項前:
```math
\begin{bmatrix}I _{G0} & I _{G1} & I _{G2} & I _{G3}\end{bmatrix} \begin{bmatrix}\phi _{n0} \\ \phi _{n1} \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}I _{Gn0} & I _{Gn1} & I _{Gn2} & I _{Gn3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _1 \\ \phi _2 \\ \phi _3\end{bmatrix}
```

移項後:
```math
\begin{bmatrix}I _{G0} & -I _{Gn1} & I _{G2} & I _{G3}\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}I _{Gn0} & -I _{G1} & I _{Gn2} & I _{Gn3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}
```

多重節点(1と3が多重節点の場合):
```math
\begin{bmatrix}0 & 1 & 0 & 0\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}0 & 0 & 0 & 1\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}
```


[./BEM_solveBVP.hpp#L393](./BEM_solveBVP.hpp#L393)


## ⛵️浮体動揺解析 

浮体の重心の運動方程式：

```math
m \frac{d {\boldsymbol U} _{\rm c}}{d t} = \boldsymbol{F} _{\text {ext }}+\boldsymbol{F} _{\text {hydro }}, \quad
\boldsymbol{I} \frac{d {\boldsymbol \Omega} _{\rm c}}{d t} = \boldsymbol{T} _{\text {ext }}+\boldsymbol{T} _{\text {hydro }}
```

$`{\boldsymbol U} _{\rm c}`$は浮体の移動速度．
$`\boldsymbol{F} _{\text {ext }}`$は重力などの外力，$`\boldsymbol{F} _{\text {hydro }}`$は水の力，$`\boldsymbol{T} _{\text {ext }}`$は外力によるトルク，$`\boldsymbol{T} _{\text {hydro }}`$は水の力によるトルク．
浮体が流体から受ける力$`\boldsymbol{F} _{\text {hydro }}`$は，浮体表面の圧力$`p`$を積分することで得られ，
また圧力$`p`$は速度ポテンシャル$`\phi`$を用いて，以下のように書ける．

```math
\boldsymbol{F} _{\text {hydro }}=\int _{S} p\boldsymbol{n}  d S, \quad
p=-\rho\left(\frac{\partial \phi}{\partial t}+\frac{1}{2} (\nabla \phi)^{2}+g z\right)
```

$`\frac{\partial \phi}{\partial t}`$を$`\phi _t`$と書くことにする．この$`\phi _t`$は陽には求められない．
そこで，$`\phi`$と似た方法，BIEを使った方法で$`\phi _t`$を求める．$`\phi`$と$`\phi _n`$の間に成り立つ境界積分方程式と全く同じ式が，$`\phi _t`$と$`\phi _{nt}`$の間にも成り立つ：

```math
\alpha ({\bf{a}})\phi _t ({\bf{a}}) = \iint _\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi _t ({\bf{x}}) - \phi _t ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
```

### 🪸ノイマン境界面における$`\phi _{nt}`$の求め方 

境界面が静止しているかどうかに関わらず，流体と物体との境界では，境界法線方向速度が一致する．
境界面上の位置ベクトルを$`\boldsymbol r`$とする．
表面上のある点の移動速度$`\frac{d\boldsymbol r}{dt}`$と流体粒子の流速$`\nabla \phi`$の間には，次の境界条件が成り立つ．

```math
{\bf n}\cdot\frac{d\boldsymbol r}{dt} =  {\bf n} \cdot \nabla \phi
```

これを微分することで，$`\phi _{nt}`$を$`\phi`$と加速度$`\frac{d{\boldsymbol U} _{\rm c}}{dt}`$と角加速度$`\frac{d{\boldsymbol \Omega} _{\rm c}}{dt}`$を使って表すことができる．
[Wu (1998)](https://www.sciencedirect.com/science/article/pii/S088997469890158X)

```math
\begin{aligned}
&\rightarrow& 0& =\frac{d}{dt}\left({\bf n}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)\right) \\
&\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \frac{d}{dt}\left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)\\
&\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2}-\frac{d}{dt}\nabla \phi\right)\\
&\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2}- {\nabla \phi _t - \nabla \phi \cdot \nabla\nabla \phi}\right)\\
&\rightarrow& \phi _{nt}& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2} - \nabla \phi \cdot \nabla\nabla \phi\right)
\end{aligned}
```

ここの$`\frac{d{\bf n}}{dt}`$と$`\frac{d^2\boldsymbol r}{dt^2}`$は，$`{\boldsymbol U} _{\rm c}`$と$`\boldsymbol \Omega _{\rm c}`$を用いて，

```math
\frac{d^2\boldsymbol r}{dt^2} = \frac{d}{dt}\left({\boldsymbol U} _{\rm c} + \boldsymbol \Omega _{\rm c} \times \boldsymbol r\right),\quad \frac{d{\bf n}}{dt} = {\boldsymbol \Omega} _{\rm c}\times{\bf n}
```

$`\frac{d^2\boldsymbol r}{dt^2}`$を上の式に代入し，$`\phi _{nt}`$を求め，
次にBIEから$`\phi _t`$を求め，次に圧力$p$を求める．
そして，浮体の重さと慣性モーメントを考慮して圧力から求めた$`\frac{d^2\boldsymbol r}{dt^2}`$は，
入力した$`\frac{d^2\boldsymbol r}{dt^2}`$と一致しなければならない．

現状を整理すると，この浮体動揺解析において，知りたい未知変数は，浮体の加速度と角加速度だけ．
しかし，浮体の没水面上にある節点での圧力$`p`$が得られないと，$`\boldsymbol{F} _{\text {hydro }}`$が得られず，運動方程式から浮体加速度が計算できない．
圧力を計算するためには，$`\phi _t`$が必要で，$`\phi _t`$は簡単には得られない，という状況．

物体の加速度は， 節点における$`\{\phi _{nt0},\phi _{nt1},\phi _{nt2},..\} = \Phi _{nt}`$が分かれば求まるが，
逆に$`\phi _{nt}`$は$`\frac{d\boldsymbol U _{\rm c}}{dt}`$が分かれば求まるので

```math
\begin{align*}
&&\frac{d\boldsymbol U _{\rm c}}{dt}& = F\left(\Phi _{nt}\left(\frac{d\boldsymbol U _{\rm c}}{dt}\right)\right)\\
&\rightarrow& Q\left(\frac{d\boldsymbol U _{\rm c}}{dt}\right) &= \frac{d\boldsymbol U _{\rm c}}{dt} - F\left(\Phi _{nt}\left(\frac{d\boldsymbol U _{\rm c}}{dt}\right)\right) =0
\end{align*}
```

のように，ある関数$`Q`$のゼロを探す，根探し問題になる．
$`\phi _{nt}`$は，[ここ](../../builds/build_bem/BEM_solveBVP.hpp#L689)で与えている．


[./BEM_solveBVP.hpp#L572](./BEM_solveBVP.hpp#L572)


```math
\nabla {\bf u} = \nabla \nabla \phi =
\begin{bmatrix} \phi _{xx} & \phi _{xy} & \phi _{xz} \\
\phi _{yx} & \phi _{yy} & \phi _{yz} \\
\phi _{zx} & \phi _{zy} & \phi _{zz}
\end{bmatrix}
```

ヘッセ行列の計算には，要素における変数の勾配の接線成分を計算する[`grad_U_LinearElement`](../../builds/build_bem/BEM_utilities.hpp#L657)を用いる．
節点における変数を$`v`$とすると，$`\nabla v-{\bf n}({\bf n}\cdot\nabla v)`$が計算できる．
要素の法線方向$`{\bf n}`$が$`x`$軸方向$`{(1,0,0)}`$である場合，$`\nabla v - (\frac{\partial}{\partial x},0,0)v`$なので，
$`(0,\frac{\partial v}{\partial y},\frac{\partial v}{\partial z})`$が得られる．


[./BEM_solveBVP.hpp#L654](./BEM_solveBVP.hpp#L654)


## ⛵️造波装置など 

造波板となるobjectに速度を与えることで，造波装置などを模擬することができる．
[強制運動を課す](../../builds/build_bem/main.cpp#L288)


[./BEM_utilities.hpp#L14](./BEM_utilities.hpp#L14)


## ⛵️初期値問題 

節点の位置と速度ポテンシャル$`\phi`$に関する初期値問題を解いて行くことが，シミュレーションである．
言い換えると，節点位置$`\frac{d\bf x}{dt}`$と速度ポテンシャル$`\frac{d\phi}{dt}`$を少しずつ$`\Delta t`$ずつ時間積分することが，シミュレーションである．
ちなみに，$`\frac{d\bf x}{dt}`$や$`\frac{d\phi}{dt}`$を計算するには，境界値問題を解く必要がある．

ある時刻において，境界値問題が解けたら，$`\frac{d\bf x}{dt}`$と$`\frac{d\phi}{dt}`$はどのように計算できるだろうか．

### 🪸流速$`\frac{d\bf x}{dt}`$の計算 

ある三角要素上の接線流速$`\nabla \phi _{\parallel}`$は，線形三角要素補間を使って次のように計算する．

```math
\nabla \phi _{\parallel} = \frac{\bf n}{2A} \times (({\bf x} _2 - {\bf x} _1) \phi _0 +({\bf x} _0 - {\bf x} _2) \phi _1 + ({\bf x} _1 - {\bf x} _0) \phi _2)
```

三角要素上の流速$`\nabla \phi`$は，次のように計算する．

```math
\nabla \phi = \frac{(\phi _n) _0+(\phi _n) _1+(\phi _n) _2}{3} {\bf n} + \nabla \phi _{\parallel}
```

### 🪸$`\frac{d\phi}{dt}`$の計算 

ある流体粒子に乗ってみたときの，速度ポテンシャルの時間変化$`\frac{D \phi}{D t}`$は，次のように計算できる．

```math
\frac{D \phi}{D t} = \frac{\partial \phi}{\partial t} + \nabla \phi \cdot \nabla \phi
```

<details>
<summary>
💡 オイラー的記述
</summary>

$`\phi=\phi(t,{\bf x})`$のように書き表し，位置と空間を独立させ分けて考える方法を，オイラー的記述という．こう書くと，$`\frac{d \phi}{d t}`$は，$`\frac{\partial \phi}{\partial t}`$であり，これは，速度ポテンシャルの純粋な時間変化ではない．純粋な，ある流体粒子の速度ポテンシャルの時間変化を表すためには，位置が時間によって変わると考え，つまり$`\phi=\phi(t,{\bf x}(t))`$と一時的に考えなおし，そして，時間微分する．そうすると$`\frac{d\phi}{dt} = \frac{\partial \phi}{\partial t} + \frac{d\bf x}{dt}\cdot \nabla \phi`$となる．

</details>

ここの$`\frac{\partial \phi}{\partial t}`$の計算は簡単ではない．そこで，ベルヌーイの式（大気圧と接する水面におけるベルヌーイの式は圧力を含まず簡単）を使って，$`\frac{\partial \phi}{\partial t}`$を消去する．


[./BEM_utilities.hpp#L492](./BEM_utilities.hpp#L492)


## ⛵️その他 

### 🪸境界値問題の未知変数 

`isNeumannID_BEM`と`isDirichletID_BEM`は，節点と面の組みが，境界値問題の未知変数かどうかを判定する．
多重節点でない場合は，`{p,nullptr}`が変数のキーとなり，多重節点の場合は，`{p,f}`が変数のキーとなる．


[./BEM_utilities.hpp#L566](./BEM_utilities.hpp#L566)


### 🪸$`\phi _{nt}`$の計算で必要となる$`{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right)`$について． 

$`\nabla`$を，$`(x,y,z)`$の座標系ではなく，
面の法線方向$`{\bf n}`$を$`x`$の代わりにとり，
面に水平な方向を$`t _0,t _1`$とする座標系で考えることにして，$`\nabla^\ast`$と書くことにする．
$`{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right)`$では，$`{\bf n}`$方向成分だけをとる操作をしているので，
新しい座標系でも同じようにすれば，結果は変わらない．

```math
{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right) =  {(1,0,0)}\cdot\left({\nabla^* \phi \cdot \nabla^* \nabla^* \phi}\right).
\quad
\nabla^* \phi = \left(\phi _n, \phi _{t _0}, \phi _{t _1}\right),
\quad \nabla^* \nabla^* \phi =
\begin{bmatrix}
\phi _{nn} & \phi _{nt _0} & \phi _{nt _1} \\
\phi _{t _0n} & \phi _{t _0t _0} & \phi _{t _0t _1} \\
\phi _{t _1n} & \phi _{t _1t _0} & \phi _{t _1t _1}
\end{bmatrix}
```

最後に第１成分だけが残るので，

```math
{(1,0,0)}\cdot\left({\nabla^* \phi \cdot \nabla^* \nabla^* \phi}\right) = \nabla^* \phi \cdot (\phi _{nn}, \phi _{t _0n}, \phi _{t _1n})
```

$`\phi _{nn}`$は，直接計算できないが，ラプラス方程式から$`\phi _{nn}=- \phi _{t _0t _0}- \phi _{t _1t _1}`$となるので，水平方向の勾配の計算から求められる．


[./BEM_utilities.hpp#L624](./BEM_utilities.hpp#L624)


### 🪸計算の流れ 

1. 境界条件の設定
2. 境界値問題（BIE）を解き，$`\phi`$と$`\phi _n`$を求める
3. 三角形の線形補間を使って節点の流速を計算する
4. 次時刻の$`\Omega(t+\Delta t)`$がわかるので，修正流速を計算する
5. 浮体の加速度を計算する．境界値問題（BIE）を解き，$`\phi _t`$と$`\phi _{nt}`$を求め，浮体面上の圧力$`p`$を計算する必要がある
6. 全境界面の節点の位置を更新．ディリクレ境界では$`\phi`$を次時刻の値へ更新


[./main.cpp#L250](./main.cpp#L250)


---
# 🐋入力ファイル生成 `input_generator.py` 

This Python script generates input files for the BEM simulation code. It supports various simulation cases and handles input file generation for each case.

## ⛵️Usage 

1. Make sure the required dependencies are installed.
2. Run the script using the following command:

```
python3 input_generator.py
```

Upon running the script, it will generate input files in JSON format for the specified simulation case. The input files are saved in the `./input_files/` directory.

## ⛵️Customization 

To customize the input file generation for a specific case, follow these steps:

1. Locate the `SimulationCase` variable in the script and set it to the desired case name, e.g., `"Kramer2021"`.
2. Add a new `case` block in the `match SimulationCase:` section to handle the new simulation case.
3. Define the required parameters for the simulation case within the new `case` block, following the examples provided in the script.
4. Update the `inputfiles` variable with the new input objects created for the custom case.

After customizing the script, run it again to generate the input files for the new case.

## ⛵️Output 

The script will generate input files in JSON format for the specified simulation case. The input files will be saved in the `./input_files/` directory. The generated input files can be used to run the BEM simulation.


[./input_generator.py#L1](./input_generator.py#L1)


---
# 🐋コンパイルと計算の実行方法 

This is a C++ implementation of a BEM simulation code. Follow the instructions below to build and run the simulation.

## ⛵️Prerequisites 

- CMake
- LAPACK library
- Python 3 for input generation

## ⛵️Building the Code 

1. Clean the build directory:

```
sh clean
```

2. Configure the build using CMake:

```
cmake -DCMAKE_BUILD_TYPE=Release ../
```

3. Compile the code:

```
make
```

## ⛵️Running the Simulation 

1. Generate input files using the `input_generator.py` script:

```
python3 ./input_generator.py
```

2. Run the simulation with the generated input files:

```
./main ./input_files/Kramer2021_H00d03
```

## ⛵️Output 

The simulation results will be stored in the specified output directory.

[![Banner](sample0.gif)](sample0.gif)

[![Banner](sample1.gif)](sample1.gif)


[./main.cpp#L1](./main.cpp#L1)


---
