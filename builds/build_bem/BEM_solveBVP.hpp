#pragma once

#include "BEM_utilities.hpp"
#include "Network.hpp"
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

## 境界値問題

### 基礎方程式

非粘性渦なし流れを仮定し，ラプラス方程式を満たす速度ポテンシャル$`\phi(t,\bf{x})`$によって流れ場$`\bf{u}(t,\bf{x})=\nabla\phi(t,\bf{x})`$を表す．水面，壁面，浮体表面における境界条件は，

```math
\begin{align}
\nabla\cdot\nabla \phi& = 0&&\text{in}&&{\bf x} \in \Omega(t),\\
\frac{\partial\phi}{\partial t} +\frac{1}{2}\nabla\phi\cdot\nabla\phi + g z &=0 &&\text{on}&&{\bf x} \in \Gamma^{\rm D}(t),\\
\phi_n + {{\bf u}_b}\cdot{{\bf n}_b} &=0&&\text{on}&&{\bf x}\in \Gamma^{\rm N}(t),
\end{align}
```

ここで，
$`{\bf x} ={(x,y,z)}`$は空間座標，$`{\bf u}_b`$は物体の流速，
$`{\bf n}_b`$は物体の外向き単位法線ベクトル，
$`\nabla=(\frac{\partial}{\partial x},\frac{\partial}{\partial y},\frac{\partial}{\partial z})`$
である．
また，$`\phi_n`$は境界面上での外向き法線方向の流速を表し，
境界面上の外向き単位法線ベクトル$`\bf n`$を使えば$`\phi_n ={\nabla\phi}\cdot {\bf n}`$で表される．

### 境界積分方程式（BIE）

**グリーンの定理**

任意の$`\phi`$，$`G`$に対して次が成り立つ（**グリーンの定理**）．

```math
\iiint_\Omega \left(G({\bf x},{\bf a})\nabla^2 \phi({\bf x}) - \phi({\bf x})\nabla^2 G({\bf x},{\bf a})\right)dV
= \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
```


$`\phi`$がラプラス方程式$`\nabla^2\phi=0`$を満たし，$`G=1/\|{\bf x}-{\bf a}\|`$とすると，
グリーンの定理から$`\phi`$と$`\phi_n`$の関係式，BIEが得られる．

```math
\alpha ({\bf{a}})\phi ({\bf{a}}) = \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
```

ここで，$`{\bf a}`$は境界面上の位置ベクトルであり，この原点$`{\bf a}`$を固定し$`{\bf x}`$について面積分される．
$`G`$は任意のスカラー関数で$`G=1/\|{\bf x}-{\bf a}\|`$とすることで，グリーンの定理の体積積分が消え，BIEの左辺のように，
原点での立体角$`\alpha\left( {\bf{a}} \right)`$とポテンシャル$`\phi( {\bf{a}})`$の積だけが残る．

<img src="schematic_BIE.png" width="400px">

この式は，流体内部では，$`\alpha ({\bf{a}})`$は$`1`$とできる．
この式は，$`\bf{a}`$におけるポテンシャル$`\phi ({\bf{a}})`$が，右辺の１重層ポテンシャルと２重層ポテンシャルの和で表されることを示している．
$`G=1/\|{\bf x}-{\bf a}\|`$がラプラス方程式の基本解であり，$`\phi`$は境界におけるポテンシャルの分布である．

*/

// #define solve_equations_on_all_points
// #define solve_equations_on_all_points_rigid_mode
// #define solveBVP_debug

// #define use_CG
// #define use_gmres
#define use_lapack
//
//
struct calculateFluidInteraction {
   const Network *PasObj;
   std::unordered_set<networkFace *> actingFaces;
   Tddd simplified_drag, simplified_drag_torque;
   double area;
   calculateFluidInteraction(const auto &faces /*waterfaces*/, const Network *PasObjIN) : PasObj(PasObjIN) {
      // PasObjと接したfaceの頂点にpressureが設定されている前提
      int count = 0;
      for (const auto &f : faces)
         if (f->Neumann) {
            if (std::ranges::all_of(f->getPoints(), [&](const auto &p) { return std::ranges::any_of(p->getContactFaces(), [&](const auto &F) { return F->getNetwork() == PasObj; }); })) {
               auto [p0, p1, p2] = f->getPoints();
               auto result = this->actingFaces.emplace(f);
               if (result.second)
                  count++;
            }
         }

      // calculate area
      std::array<double, 3> P012;
      std::array<std::array<double, 3>, 3> X012;
      area = 0.;
      for (const auto &f : this->actingFaces) {
         auto [p0, p1, p2] = f->getPoints();
         P012 = {p0->pressure, p1->pressure, p2->pressure};
         X012 = {p0->X, p1->X, p2->X};
         auto intpX = interpolationTriangleLinear0101(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            area += intpX.J(x0, x1) * w0w1;
      }
      std::cout << "接触している面の数:" << count << " 表面積:" << area << std::endl;
   };

   // \label{BEM:surfaceIntegralOfPressure}
   std::array<Tddd, 2> surfaceIntegralOfPressure() {
      Tddd force = {0., 0., 0.}, torque = {0., 0., 0.};
      std::array<double, 3> P012;
      std::array<std::array<double, 3>, 3> X012;
      for (const auto &f : this->actingFaces) {
         auto [p0, p1, p2] = f->getPoints();
         P012 = {p0->pressure, p1->pressure, p2->pressure};
         X012 = {p0->X, p1->X, p2->X};  // auto [pre0, pre1, pre2] = P012;
         // auto [X0, X1, X2] = X012;
         // this->force += (p0 + p1 + p2) / 3. * Cross(X1 - X0, X2 - X0) / 2.;
         auto intpP = interpolationTriangleLinear0101(P012);
         auto intpX = interpolationTriangleLinear0101(X012);
         auto n = TriangleNormal(X012);
         Tddd F_tmp = {0., 0., 0.}, T_tmp = {0., 0., 0.};
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple) {
            auto f = intpP(x0, x1) * intpX.J(x0, x1) * w0w1 * n;
            F_tmp += f;
            T_tmp += Cross(intpX(x0, x1) - this->PasObj->COM, f);
         }
         force += F_tmp;
         torque += T_tmp;
      }
      return {force, torque};
   };

   std::array<Tddd, 2> surfaceIntegralOfVerySimplifiedDrag() {
      this->simplified_drag.fill(0.);
      this->simplified_drag_torque.fill(0.);
      std::array<std::array<double, 3>, 3> X012;
      for (const auto &f : this->actingFaces) {
         auto [p0, p1, p2] = f->getPoints();
         std::array<double, 3> P012 = {p0->pressure, p1->pressure, p2->pressure};
         auto X0 = p0->X;
         auto X1 = p1->X;
         auto X2 = p2->X;
         X012 = {X0, X1, X2};
         const Tddd relative_U0 = p0->U_BEM - PasObj->velocityRigidBody(X0);
         const Tddd relative_U1 = p1->U_BEM - PasObj->velocityRigidBody(X1);
         const Tddd relative_U2 = p2->U_BEM - PasObj->velocityRigidBody(X2);
         // this->force += (p0 + p1 + p2) / 3. * Cross(X1 - X0, X2 - X0) / 2.;
         auto intpRelativeVelocity = interpolationTriangleLinear0101(T3Tddd{relative_U0, relative_U1, relative_U2});
         auto intpX = interpolationTriangleLinear0101(X012);
         const double nu = 10 * 1000 * 1000 * 1.004 * 10E-6;  // m2 /s
         Tddd drag_f;
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple) {
            drag_f = nu * intpRelativeVelocity(x0, x1) * intpX.J(x0, x1) * w0w1;
            this->simplified_drag += drag_f;
            this->simplified_drag_torque += Cross(intpX(x0, x1) - this->PasObj->COM, drag_f);
         }
      }
      return {this->simplified_drag, this->simplified_drag_torque};
   };
};

// b@ -------------------------------------------------------------------------- */
// b@                                   BEM_BVP                                  */
// b@ -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

### BIEの離散化

```math
\alpha ({\bf a})\phi({\bf a})
= \iint_\Gamma {\left({
\frac{1}{\|{\bf x}-{\bf a}\|}
\nabla \phi ({\bf{x}}) + \phi ({\bf{x}})
\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}}
\right) \cdot {\bf{n}}({\bf{x}})dS}
```

面は面上の節点を使って補間され，面積分はこの補間された面上に沿って行われる．
面の法線ベクトル$`{\bf n}=\frac{\frac{{\partial {\bf x}}}{{\partial \xi_0}}\times\frac{{\partial {\bf x}}}{{\partial \xi_1}}}{\left\|\frac{{\partial {\bf x}}}{{\partial \xi_0}}\times\frac{{\partial {\bf x}}}{{\partial \xi_1}}\right\|}`$を代入し，BIEをGauss-Legendre積分で離散化すると，

```math
\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} {\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left( {\sum\limits_{j=0}^2 {{{\left( {{\phi_n}} \right)}_{k_\vartriangle,j }}{N_{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}_{k _\vartriangle}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}\left\|\frac{{\partial{{\bf x}_{k _\vartriangle}}}}{{\partial{\xi_0}}} \times \frac{{\partial{\bf{x}}_{k _\vartriangle}}}{{\partial{\xi_1}}}\right\|} \right)} }=
```

```math
\alpha_{i_\circ}(\phi)_{i_\circ}-\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} \sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{\bf x}_{k \vartriangle}}({\pmb{\xi}})-{{\bf x}_{i_\circ} }}{{{{\| {{{\bf x}_{k_\vartriangle}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}} \cdot\left(\frac{{\partial {\bf{x}}_{k_\vartriangle}}}{{\partial {\xi_0}}}\times\frac{{\partial {\bf{x}}_{k_\vartriangle}}}{{\partial {\xi_1}}}\right)}\right)}
```

離散化では，$`\phi_{i_\circ}`$と$`{\phi_n}_{i_\circ}`$の係数を知りたいので，
$`\phi_{k_\vartriangle}({\pmb{\xi}})`$と$`{\phi_n}_{k_\vartriangle}({\pmb{\xi}})`$と書くのではなく，
$`\phi_{i_\circ}`$と$`{\phi_n}_{i_\circ}`$が見えるように$`\phi_{k_\vartriangle}({\pmb{\xi}})`$と$`{\phi_n}_{k_\vartriangle}({\pmb{\xi}})`$の補間を書いている．

ここで，$`\phi_{k_\vartriangle,j}`$における$`k_\vartriangle`$は三角形要素の番号，$`j`$は三角形要素の頂点番号．
$`N_j`$は三角形要素の形状関数，$`{\pmb{\xi}}`$は三角形要素の内部座標，$`w_0,w_1`$はGauss-Legendre積分の重み，$`\alpha_{i_\circ}`$は原点$`i_\circ`$における立体角，$`\phi`$はポテンシャル，$`\phi_n`$は法線方向のポテンシャル，$`\bf{x}`$は空間座標，$`{\bf x}_{i_\circ}`$は原点の空間座標である．

* $`\phi_{k_\vartriangle}`$は補間で作った関数
* $`\phi_{k_\vartriangle,j}`$は補間を構成する節点$`j`$での値
* $`\phi_{i_\circ}`$はより直接的にある節点$`i_\circ`$での値

NOTE: この段階ではまだ，1.数値積分のパラメタと，2.形状関数のパラメタと元々の面都の対応関係は，指定していない．例えば，やり方によっては$`\xi_1`$のパラメタは，$`\xi_0`$に依存するかもしれない．

補間に使うパラメタを$`{\bf \xi}=(\xi_0, \xi_1)`$として，よく使われる３節点を使う線形補間を使うことにする．
元の面に対応する，線形補間面は，パラメタ上では$`{\xi_0 + \xi_1 = 1}`$を満たす範囲なので，
積分範囲は例えば$`0\leq \xi_0 \leq 1, 0\leq \xi_1 \leq 1-\xi_0`$となる．
しかし，数値積分につかう変数と重みの組み合わせは，コンパイルタイムに決めておき計算を効率化したいので，
この点で，変化する積分範囲は数値積分との相性が悪い．

### 線形三角要素

#### 線形三角要素

<img src="./img/schematic_linear_triangle_element.png" width="400px">

形状関数$`{\pmb N}_j({\pmb \xi}),{\pmb \xi}=(\xi_0,\xi_1)`$は，$`\xi_0,\xi_1`$が$`0`$から$`1`$動くことで，範囲で三角要素全体を動くように定義している．

```math
{\pmb N}({\pmb \xi}) = (N_0({\pmb \xi}),N_1({\pmb \xi}),N_2({\pmb \xi})) = (\xi_0, - \xi_1 (\xi_0 - 1), (\xi_0-1)(\xi_1-1))
```

####  線形三角要素のヤコビアン

線形三角要素のヤコビアンは，$`\|\frac{\partial {\bf{x}}}{\partial {\xi_0}} \times \frac{\partial {\bf{x}}}{\partial {\xi_1}}\|`$である．

```Mathematica
shape[t0_, t1_] := With[{t2 = 1 - t0 - t1, t0m1 = t0 - 1, t1m1 = t1 - 1}, {t0, -t1*t0m1, t0m1*t1m1}];
D0shape[t0_, t1_] = (D[shape[T0, t1], T0] /. T0 -> t0);
D1shape[t0_, t1_] = (D[shape[t0, T1], T1] /. T1 -> t1);
{a, b, c} = {{x0, y0, z0}, {x1, y1, z1}, {x2, y2, z2}}
FullSimplify[Cross[Dot[D[shape[T0, t1], T0], {a, b, c}], Dot[D[shape[t0, T1], T1], {a, b, c}]]]
FullSimplify[Cross[Dot[D[shape[T0, t1], T0], {a, b, c}], Dot[D[shape[t0, T1], T1], {a, b, c}]]/Cross[b - a, c - a]]
```

上の結果は，$`1-\xi_0`$となる．つまり，線形補間の場合，ヤコビアン内の外積は次のように，節点位置を使ってシンプルに計算できる．

```math
\frac{\partial {{\bf x}_{{k _\vartriangle}}}}{\partial {\xi_0}} \times \frac{\partial {{\bf x}_{{k _\vartriangle}}}}{\partial {\xi_1}} = (1-\xi_0) (({{\bf x} _{{k _\vartriangle}_1}}-{{\bf x} _{{k _\vartriangle}_0}})\times({{\bf x} _{{k _\vartriangle}_2}}-{{\bf x} _{{k _\vartriangle} _0}}))
= 2(1-\xi_0)A_{k_\vartriangle}{\bf n}_{k_\vartriangle}
```

これを使えば，BIEは次のように簡単になる．

```math
\sum\limits_{k_\vartriangle}{2A_{k_\vartriangle}}
\sum\limits_{{\xi_1},{w_1}}
{\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left( {\sum\limits_{j=0}^2 {{{\left( {{\phi_n}} \right)}_{k_\vartriangle,j }}{N_{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}
(1-\xi_0)
} \right)} }=
```

```math
\alpha_{i_\circ}(\phi)_{i_\circ}
-\sum\limits_{k_\vartriangle}{2A_{k_\vartriangle}{\bf n}_{k_\vartriangle}}\cdot
\sum\limits_{{\xi_1},{w_1}}
\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{\bf x}_{k _\vartriangle}}({\pmb{\xi}})-{{\bf x}_{i_\circ} }}{{{{\| {{{\bf x}_{k_\vartriangle}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}}} (1-\xi_0)\right)}
```

NOTE: ちなみに，$`\frac{1-\xi_0}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}`$の分子に$`1-\xi_0`$があることで，
関数の特異的な変化を抑えることができる．プログラム上ではこの性質が利用できるように，この分数をまとめて計算している．

*/

std::array<double, 3> weight(double t0, double t1) {

   auto shape = [](double t0, double t1) -> std::array<double, 3> {
      return {t0, t1, 1 - t0 - t1};
   };

   constexpr std::array<std::array<double, 2>, 3> vertex = {{{std::cos(M_PI / 2), std::sin(M_PI / 2)},
                                                             {std::cos(7 * M_PI / 6), std::sin(7 * M_PI / 6)},
                                                             {std::cos(11 * M_PI / 6), std::sin(11 * M_PI / 6)}}};
   auto s = shape(t0, t1);
   auto s_half_half = shape(2. / 3., 2. / 3.);
   auto s_minus_half_half = shape(-1. / 3., 2. / 3.);
   auto s_half_minus_half = shape(2. / 3., -1. / 3.);

   std::array<double, 4> distances;

   distances[0] = Norm(Dot(s, vertex));
   distances[1] = Norm(Dot(s, vertex) - Dot(s_half_half, vertex));
   distances[2] = Norm(Dot(s, vertex) - Dot(s_minus_half_half, vertex));
   distances[3] = Norm(Dot(s, vertex) - Dot(s_half_minus_half, vertex));

   auto weight_func = [&](const double x) {
      const double h = 1 + 1E-10;
      return std::pow(w_Linear(x, h), 1);
      // return std::pow(2. * x, -1.);
      // return std::pow(std::cosh(M_PI * x), -1.);
   };

   // double w0 = weight_func(distances[0]);
   double w1 = weight_func(distances[1]);
   double w2 = weight_func(distances[2]);
   double w3 = weight_func(distances[3]);
   double total = w1 + w2 + w3;
   // w0 /= total;
   w1 /= total;
   w2 /= total;
   w3 /= total;

   return {w1, w2, w3};
}

struct BEM_BVP {
   const bool Neumann = false;
   const bool Dirichlet = true;
   std::vector<Network *> WATERS;
   lapack_lu *lu = nullptr;

   using T_PBF = std::tuple<netP *, netF *>;
   using mapTPBF_Tdd = std::map<T_PBF, Tdd>;
   using mapTPBF_mapTPBF_Tdd = std::map<T_PBF /*タプル*/, mapTPBF_Tdd>;
   using map_P_Vd = std::map<netP *, V_d>;

   //@ 各バケツでのモーメントを次数別に保存する．(ユニーク) p->{k,m,Yn,Y}ベクトル
   using uo_P_uoTiiTdd = std::unordered_map<networkPoint *, std::unordered_map<Tii /*k,m*/, Tdd /*YYn*/>>;
   using V_uo_P_uoTiiTdd = std::vector<uo_P_uoTiiTdd>;
   using VV_uo_P_uoTiiTdd = std::vector<V_uo_P_uoTiiTdd>;
   using VVV_uo_P_uoTiiTdd = std::vector<VV_uo_P_uoTiiTdd>;
   VV_d mat_ukn, mat_kn;
   V_d knowns, b_RHS;
   std::vector<std::vector<std::array<double, 2>>> IGIGn;
   BEM_BVP(std::vector<Network *> WATERS) : WATERS(WATERS) {};
   ~BEM_BVP() {
      if (this->lu)
         delete this->lu;
   };

   /**
   isNeumannID_BEMとisDirichletID_BEMの両方を満たす{p,f}は存在しない．
   */

   void setIGIGn() {

#define use_rigid_mode

      for (const auto &water : WATERS) {
         water->setGeometricProperties();
#pragma omp parallel
         for (const auto &integ_f : water->getSurfaces())
#pragma omp single nowait
            integ_f->setIntegrationInfo();
      }

      this->IGIGn.resize(this->matrix_size, std::vector<std::array<double, 2>>(this->matrix_size));
#pragma omp parallel
      for (auto &IGIGn_i : this->IGIGn)
#pragma omp single nowait
         for (auto &IGIGn_ij : IGIGn_i)
            IGIGn_ij = {0., 0.};

      TimeWatch timer;
      std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;

      /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

      #### 係数行列の作成

      数値シミュレーションでは，境界値問題を$`{\bf A}{\bf x}={\bf b}`$のような線形連立方程式になるよう近似，変形し（離散化），$`{\bf x}`$を求めることが多い．
      BEMでもBIEを離散化してこのような形にする．その際，境界条件に応じて，方程式（$`{\bf A}{\bf x}={\bf b}`$の行）の右辺と左辺が入れ替える必要があるので注意する．
      これは，$`{\bf A}{\bf x}={\bf b}`$の未知変数$`{\bf x}`$と既知変数$`{\bf b}`$がポテンシャル$`\phi`$か法線方向のポテンシャル$`\phi_n`$か，境界条件によって違うからである．
      プログラム上では，係数行列$`\bf A`$やベクトル$`\bf b`$を境界条件に応じて適切に作成すれば，求まる$`\bf x`$が適切なものになる．

      $`\phi`$の係数行列を$`\mathbf{M}`$，$`\phi_n`$の係数行列を$`\mathbf{N}`$，$`\mathbf{\Phi}`$を$`\phi`$のベクトル，$`\mathbf{\Phi_n}`$を$`\phi_n`$のベクトルとして，
      次のような連立一次方程式を得る．

      ```math
      \mathbf{N} \mathbf{\Phi_n} = \mathbf{M} \mathbf{\Phi} \rightarrow {\bf A}{\bf x}={\bf b}
      ```

      このプログラムでは，$`A`$を`IGIGn`，$`b`$を`knowns`としている．

      このループでは，BIEの連立一次方程式の係数行列`IGIGn`を作成する作業を行なっている．
      `IGIGn`は，ある節点$`i_\circ`$（係数行列の行インデックス）に対する
      他の節点$`j_\circ`$（係数行列の列インデックス）の影響度合いのようなものである．
      その影響度合いは，他の節点$`j_\circ`$の所属する要素までの距離や向きによって決まることが離散化された式からわかる．

      | Variable | Description |
      |:--------:|:-----------:|
      | `origin` | 原点となる節点$`i_\circ`$ |
      | `integ_f` | Element $`k_{\triangle}`$ |
      | `t0, t1, ww` | Gaussian points and thier wieghts $`\xi_0, \xi_1, w_0 w_1`$ |
      | `p0, p1, p2` | Node of the element $`k_{\triangle}`$ |
      | `N012` | Shape function $`\pmb{N}_j`$ |
      | `IGIGn` | Coefficient matrices of the left and right sides |
      | `nr` | $`\| \pmb{x} - \pmb{x}_{i\circ } \|`$ |
      | `tmp` | $`w_0 w_1 \frac{1 - \xi_0}{\| \pmb{x} - \pmb{x}_{i\circ } \|}`$ |
      | `cross` | $`\frac{\partial \pmb{x}}{\partial \xi_0} \times \frac{\partial \pmb{x}}{\partial \xi_1}`$ |

      */

      for (const auto &water : WATERS) {
         water->setFacesVector();
         water->setPointsVector();
      }

      const T3Tdd shape_map_center = {{{0., 0.5} /*quad 4 -> linear 0*/, {0.5, 0.} /*quad 5 -> linear 1*/, {0.5, 0.5} /*quad 3 -> linear 2*/}};
      //! これとN3の内積を新なパラメタとして利用すると，(t0,t1)=(1,0)で{0., 0.5}に，(t0,t1)=(0,1)で{0.5, 0.}に，t0=t1=0で{0.5, 0.5}になる．
      const T3Tdd shape_map_l0_face = {{{0.5, 0.} /*0*/, {0., 0.5} /*1*/, {0., 0.} /*2*/}};
      const T3Tdd shape_map_l1_face = {{{0., 0.} /*0*/, {0.5, 0.} /*1*/, {0., 0.5} /*2*/}};
      const T3Tdd shape_map_l2_face = {{{0., 0.5} /*0*/, {0., 0.} /*1*/, {0.5, 0.} /*2*/}};

      constexpr std::array<double, 2> ZEROS2 = {0., 0.};
      constexpr std::array<double, 3> ZEROS3 = {0., 0., 0.};

      // // auto t0_t1_ww_N012_HIGHRESOLUTION = t0_t1_ww_N012_LOWRESOLUTION;

      // for (int i = 0; const auto &[t0, t1, ww] : __array_GW10xGW10__) {
      //    auto t0t1t2 = ModTriShape<3>(t0, t1);
      //    t0_t1_ww_N012_HIGHRESOLUTION.push_back({t0, t1, ww, t0t1t2});
      // }

      /*
      ## 積分の効率化

      `linear_triangle_integration_info`と`pseudo_quadratic_triangle_integration_info`は，
      積分点の位置と重みを事前に計算しておくことで，積分の効率化を図るためのものである．
      `linear_triangle_integration_info`と`pseudo_quadratic_triangle_integration_info`の引数の詳細

      ### `linear_triangle_integration_info`
      | 引数名  | 説明                                                              |
      |---------|-----------------------------------------------------------------|
      | `Tdd`   | 2D パラメータ {[0,1], [0,1]}：積分変数                             |
      | `double`| ガウス重み（積分重み）                                             |
      | `Tddd`  | 3D パラメータ {xi0=[0,1], xi1=[0,1-xi0]}：2Dパラメータと関連        |
      | `Tddd`  | 3D 位置ベクトル（{xi0, xi1, xi2}を使用）                           |
      | `Tddd`  | 外積（dX/dxi0 × dX/dxi1）                                        |
      | `double`| 外積のノルム                                                     |

      ### `pseudo_quadratic_triangle_integration_info`
      | 引数名              | 説明                                                              |
      |---------------------|-----------------------------------------------------------------|
      | `Tdd`               | 2D パラメータ {[0,1], [0,1]}：積分変数                             |
      | `double`            | ガウス重み（積分重み）                                             |
      | `Tddd`              | 3D パラメータ {xi0=[0,1], xi1=[0,1-xi0]}：2Dパラメータと関連        |
      | `std::array<T6d, 4>`| 二次要素の形状関数                                                |
      | `Tddd`              | 3D 位置ベクトル（{xi0, xi1, xi2}を使用）                           |
      | `Tddd`              | 外積（dX/dxi0 × dX/dxi1）                                        |
      | `double`            | 外積のノルム                                                     |
      */

      int count_pseudo_quadratic_element = 0, count_linear_element = 0, total = 0;
      for (const auto &water : WATERS)
         for (const auto &integ_f : water->getSurfaces()) {
            if (integ_f->isLinearElement)
               count_linear_element++;
            else if (integ_f->isPseudoQuadraticElement)
               count_pseudo_quadratic_element++;
            total++;
         }

      std::cout << "線形要素の面の数：" << count_linear_element << " persecent: " << 100. * count_linear_element / total << std::endl;
      std::cout << "擬似二次要素の面の数：" << count_pseudo_quadratic_element << " persecent: " << 100. * count_pseudo_quadratic_element / total << std::endl;

      if (_PSEUDO_QUADRATIC_ELEMENT_)
         std::cout << "擬似二次要素を使ってBIEを離散化" << std::endl;
      else
         std::cout << "線形要素を使ってBIEを離散化" << std::endl;

      for (const auto water : WATERS) {
         double scale = water->getScale();
         auto surfacePoints = water->getSurfacePoints();
         auto surfaces = water->getSurfaces();
#pragma omp parallel
         for (const auto &origin : surfacePoints)
#pragma omp single nowait
         {
            //@ this loop is for the multiple nodes
            for (const auto &[f, index] : origin->f2Index) {
               double origin_ign_rigid_mode = 0.;
               auto &IGIGn_Row = IGIGn[index];
               double nr, ig, ign, ww_nr;
               Tdd ig_ign0, ig_ign1, ig_ign2;
               Tddd R;
               networkPoint *closest_p_to_origin = nullptr;
               std::vector<std::tuple<networkPoint *, networkFace *, double, double>> key_ig_ign;

               //@ for all water faces
               for (const auto &water : WATERS) {
                  // b@ integrate over all faces
                  for (const auto &integ_f : surfaces) {

                     auto [p0, p1, p2] = integ_f->getPoints(origin);
                     closest_p_to_origin = p0;

                     /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

                     WARNING: この`std::vector<std::tuple<networkPoint *, networkFace *, double, double>> key_ig_ign`の`networkFace`は，どの面側から節点を呼び出すかを決めていて，高次補間の場合，積分面と一致しない場合がある．

                     1. fill key_ig_ign
                     2. fill IGIGn_Row

                     */
                     auto dist = Norm((p0->X + p1->X + p2->X) / 3. - origin->X);
                     int how_far = 0;
                     if (dist < scale / 20.)
                        how_far = 1;

                     if (integ_f->isLinearElement) {
                        ig_ign0 = ig_ign1 = ig_ign2 = {0., 0.};
                        if (p0 == origin) {
                           ign = 0.;
                           for (const auto &[t0t1, ww, shape3, X, cross, norm_cross] : integ_f->map_Point_LinearIntegrationInfo_vector[1].at(closest_p_to_origin)) {
                              ig = norm_cross * (ww / (nr = Norm(R = (X - origin->X))));
                              std::get<0>(ig_ign0) += ig * std::get<0>(shape3);
                              std::get<0>(ig_ign1) += ig * std::get<1>(shape3);
                              std::get<0>(ig_ign2) += ig * std::get<2>(shape3);
                           }
                        } else {
                           for (const auto &[t0t1, ww, shape3, X, cross, norm_cross] : integ_f->map_Point_LinearIntegrationInfo_vector[how_far].at(closest_p_to_origin)) {
                              ig = norm_cross * (ww_nr = ww / (nr = Norm(R = (X - origin->X))));
                              ign = Dot(R, cross) * ww_nr / (nr * nr);
                              std::get<0>(ig_ign0) += ig * std::get<0>(shape3);
                              std::get<1>(ig_ign0) -= ign * std::get<0>(shape3);
                              std::get<0>(ig_ign1) += ig * std::get<1>(shape3);
                              std::get<1>(ig_ign1) -= ign * std::get<1>(shape3);
                              std::get<0>(ig_ign2) += ig * std::get<2>(shape3);
                              std::get<1>(ig_ign2) -= ign * std::get<2>(shape3);
                              origin_ign_rigid_mode += ign;
                           }
                        }
                        IGIGn_Row[pf2Index(p0, integ_f)] += ig_ign0;
                        IGIGn_Row[pf2Index(p1, integ_f)] += ig_ign1;
                        IGIGn_Row[pf2Index(p2, integ_f)] += ig_ign2;
                     } else if (integ_f->isPseudoQuadraticElement) {
                        key_ig_ign = integ_f->map_Point_BEM_IGIGn_info_init.at(closest_p_to_origin);
                        for (const auto &[t0t1, ww, Nc_N0_N1_N2, X, cross, norm_cross] : integ_f->map_Point_PseudoQuadraticIntegrationInfo_vector[how_far].at(closest_p_to_origin)) {
                           ig = norm_cross * (ww_nr = ww / (nr = Norm(R = (X - origin->X))));
                           ign = Dot(R, cross) * ww_nr / (nr * nr);
                           for (auto i = 0; i < 6; ++i) {
                              std::get<2>(key_ig_ign[i]) += ig * std::get<0>(Nc_N0_N1_N2)[i];
                              std::get<3>(key_ig_ign[i]) -= ign * std::get<0>(Nc_N0_N1_N2)[i];
                              if (std::get<0>(key_ig_ign[i]) != origin)
                                 origin_ign_rigid_mode += ign * std::get<0>(Nc_N0_N1_N2)[i];

                              std::get<2>(key_ig_ign[i + 6]) += ig * std::get<1>(Nc_N0_N1_N2)[i];
                              std::get<3>(key_ig_ign[i + 6]) -= ign * std::get<1>(Nc_N0_N1_N2)[i];
                              if (std::get<0>(key_ig_ign[i + 6]) != origin)
                                 origin_ign_rigid_mode += ign * std::get<1>(Nc_N0_N1_N2)[i];

                              std::get<2>(key_ig_ign[i + 12]) += ig * std::get<2>(Nc_N0_N1_N2)[i];
                              std::get<3>(key_ig_ign[i + 12]) -= ign * std::get<2>(Nc_N0_N1_N2)[i];
                              if (std::get<0>(key_ig_ign[i + 12]) != origin)
                                 origin_ign_rigid_mode += ign * std::get<2>(Nc_N0_N1_N2)[i];

                              std::get<2>(key_ig_ign[i + 18]) += ig * std::get<3>(Nc_N0_N1_N2)[i];
                              std::get<3>(key_ig_ign[i + 18]) -= ign * std::get<3>(Nc_N0_N1_N2)[i];
                              if (std::get<0>(key_ig_ign[i + 18]) != origin)
                                 origin_ign_rigid_mode += ign * std::get<3>(Nc_N0_N1_N2)[i];
                           }
                        }

                        for (const auto &[p, integ_f, ig, ign] : key_ig_ign)
                           IGIGn_Row[pf2Index(p, integ_f)] += Tdd{ig, ign};  // この面に関する積分において，φまたはφnの寄与
                     }
                  }
               }

               /* -------------------------------------------------------------------------- */
               /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

               ### リジッドモードテクニック（係数行列の対角成分の計算）

               BIEの対角成分の計算で注意が必要なのは，原点$`i_\circ`$の頂点の立体角と，係数の特異性である．

               * 係数行列の対角成分には，立体角$`\alpha`$が含まれており，この計算は面倒である．
               * 係数の計算には，$`\frac{{\mathbf{x}_{k_\vartriangle}({\pmb{\xi}}) - \mathbf{x}_{i_\circ}}}{{\| \mathbf{x}_{k_\vartriangle}({\pmb{\xi}}) - \mathbf{x}_{i_\circ} \|}^3}`$が含まれており，分母が0付近で強い特異性を持つ．

               そこで，素直に幾何学的な観点から立体角を計算するのではなく，BIEの式を使って積分で計算する方法がある．BIEの式に，$`\phi=1`$を代入すると，$`\phi_n`$が消える．結局，対角成分，つまり，原点$`i_\circ`$を頂点上の変数に掛かる係数は，次のようになる．

               ```math
               \sum\limits_{k_\vartriangle} 2 A_{k_\vartriangle} \, \mathbf{n}_{k_\vartriangle} \cdot \sum\limits_{\xi_1, w_1} \sum\limits_{\xi_0, w_0} \left( w_0 w_1 \left( \sum\limits_{j=0}^2 \bar\delta_{(k_\vartriangle, j),i_\circ} N_j({\pmb{\xi}}) \right) \frac{{\mathbf{x}_{k_\vartriangle}({\pmb{\xi}}) - \mathbf{x}_{i_\circ}}}{{\| \mathbf{x}_{k_\vartriangle}({\pmb{\xi}}) - \mathbf{x}_{i_\circ} \|}^3}(1 - \xi_0)\right)
               ```

               $`\bar\delta_{(k_\vartriangle, j),i_\circ}`$は，$`k_\vartriangle`$の$j$番目の頂点が$i_\circ$である場合に0，それ以外は1となる関数である．

               数値計算上は，$`\delta_{(k_\vartriangle, j),i_\circ}`$がゼロの場合は，そもそも係数をインクリメントせず，スキップする．
               これはリジッドモードテクニックと呼ばれていて，分子が小さくなる特異的な計算を省き，立体角の計算もまとめて対角成分を計算することができる方法である．

               ただし，線形要素の場合，原点$`i_\circ`$を頂点とする三角形$`k_{\vartriangle}`$に対する計算，$`{\bf n}_{k_\vartriangle}\cdot ({{\bf x}_{k_\vartriangle}}(\pmb{\xi})-{{\bf x}_{i_\circ}})=0`$となるため，和をとる必要はない．
               よって，そもそも線形要素の場合は，特異的な計算は含まれない．

               */

#if defined(use_rigid_mode)
               std::get<1>(IGIGn_Row[index]) = origin_ign_rigid_mode;
#else
               std::get<1>(IGIGn_Row[index]) += origin->getSolidAngle();
#endif
            }
         }
      }
      std::cout << Green << "離散化にかかった時間" << timer() << colorReset << std::endl;
   };

   /* -------------------------------------------------------------------------- */

   void generateBIEMatrix() {

      TimeWatch timer;
      std::cout << "generateBIEMatrix()" << std::endl;

      /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

      ### 左辺と右辺の入れ替え

      係数行列`IGIGn`は，左辺の$`I_G \phi_n`$，右辺の$`I_{G_n}\phi`$の係数行列を表している．

      ```math
      (I_G)_{i_\circ,j_\circ} (\phi_n)_{j_\circ} = (I_{Gn})_{i_\circ,j_\circ}  \phi_{j_\circ}
      ```

      境界条件に応じて，未知変数は$`\phi,\phi_n`$のどちらかに決まる．
      未知変数が$`\phi`$の場合（Dirichlet境界条件の場合），
      係数行列`IGIGn`中で対応する列を符号変えて入れ替えることで移項したことになる．

      #### ２種類の多重節点

      1. Dirichlet面上であり，かつNeumann面上である多重節点
      2. Dirichlet面上ではなく，完全にNeumann面上にあるが，法線ベクトルが大きく異なる節点

      1の多重節点の場合，BIEの連立一次方程式の係数行列の行を，Dirchlet面上の$`\phi`$とNeumann面上の$`\phi`$の値が一致する，という式に変更する．
      2の場合は，特に変更しない．BIEを解くことで，それぞれの面に対して，$`\phi`$が得られるが，それらの平均値，または重み付け平均値を$`\phi`$として採用する．

      */

      //^ ---------------------------- setIGIGnで計算した結果を基に，係数行列を作成する．---------------------------
      double max_value = 0;
      for (auto i = 0; i < IGIGn.size(); ++i) {
         if (max_value < std::abs(std::get<0>(IGIGn[i][i])))
            max_value = std::abs(std::get<0>(IGIGn[i][i]));
      }

      this->mat_kn.resize(this->matrix_size, V_d(this->matrix_size, 0.));
      this->mat_ukn.resize(this->matrix_size, V_d(this->matrix_size, 0.));
      knowns.resize(this->matrix_size);

      for (const auto water : WATERS) {
         auto surfacePoints = water->getSurfacePoints();
#pragma omp parallel
         for (const auto &a : surfacePoints)
#pragma omp single nowait
            for (const auto &[a_face, i] : a->f2Index) {
               /* -------------------------------------------------------------------------- */
               //^ 左辺と右辺の係数行列を作成する．
               auto &IGIGn_i = IGIGn[i];
               for (const auto water : WATERS)
                  for (const auto &x : surfacePoints)
                     for (const auto &[x_face, j] : x->f2Index) {
                        mat_ukn[i][j] = IGIGn_i[j][0];
                        mat_kn[i][j] = IGIGn_i[j][1];
                        if (isNeumannID_BEM(x, x_face)) {
                           // 未知変数の係数行列は左，既知変数の係数行列は右
                           std::swap(mat_ukn[i][j], mat_kn[i][j]);
                           mat_ukn[i][j] *= -1;
                           mat_kn[i][j] *= -1;
                        }
                     }

               auto isNeumanID = isNeumannID_BEM(a, a_face);
               auto isDirichletID = isDirichletID_BEM(a, a_face);
               /* -------------------------------------------------------------------------- */
               //^ 多重節点の場合，Neumann面上のphiの値は，同じ場所のDirichlet面のphiと一致する
               if (a->CORNER && isNeumanID /*行の変更*/) {
                  std::ranges::fill(mat_ukn[i], 0.);
                  std::ranges::fill(mat_kn[i], 0.);
                  mat_ukn[i][i] = max_value;  // φの系数
                  //! nullptrと指定するのは，全てNeumannのコーナに対して使えないのではないかということでこうした．
                  mat_kn[i][pf2Index(a, nullptr)] = max_value;  // φの系数移行したからマイナス？　いいえ，移項を考慮した上でこれでいい．
               }
               /* -------------------------------------------------------------------------- */
               //^ 既知変数のベクトルを作成する．
               // if (isDirichletID && isNeumanID)
               //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Inconsistent BEM condition: isDirichletID_BEM(P,F) && isNeumannID_BEM(P,F) == true");
               // else if (isDirichletID)
               //    knowns[i] = a->phi_Dirichlet = std::get<0>(a->phiphin);
               // else if (isNeumanID)
               //    knowns[i] = a->phinOnFace.at(a_face);  // はいってない？はいってた．
               // else
               //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Undefined BEM condition.");
               if (isDirichletID) {
                  knowns[i] = a->phiOnFace.at(a_face);
               } else if (isNeumanID) {
                  knowns[i] = a->phinOnFace.at(a_face);
               } else
                  throw std::runtime_error("Error: Boundary type is not defined.");
               /* -------------------------------------------------------------------------- */
            }
      }
      b_RHS = Dot(mat_kn, knowns);

      std::cout << Green << "generateBIEMatrix()" << timer() << colorReset << std::endl;
   };

   /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

   ### 高速多重極展開との関係

   GMRES法は，$`A\cdot x`$の計算を何度も行い，その線形和で解を近似するので，$`A`$をプログラム中で保持せずとも，$`A\cdot x`$を計算することができれば解を求めることができる．
   高速多重極展開は，この$`A\cdot x`$を高速に計算するための手法である．$`A\cdot {\bf x}={\bf b}`$のある行において，具体的な計算を考えてみる．

   ```math
   \begin{align*}
   A\cdot x &= b \\
   \sum\limits_{j=0}^{N-1} A_{i_\circ,j}x_j &= b_{i_\circ} \\
   \end{align*}t
   ```
   　　
   $`\sum\limits_{j=0}^{N-1} A_{{i_\circ},j}x_j = b_{i_\circ}`$は，節点$`{i_\circ}`$を原点節点としてBIEを離散化したものである．

   $`A_{i,j}({\bf a}_i)`$は，$`{\bf a}_i`$に依存しており，$`{\bf a}_i`$が変わると$`A_{i,j}({\bf a}_i)`$も変わる．
   しかし，これをソース点と観測点の関数の積と和の形に変形することできる．
   また，展開中心をソース点付近にとれば，ある変数が小さい場合限っては，その展開は早く収束する．
   ある変数とは具体的には，展開中心からソース点までの距離/展開中心から観測点までの距離である．

   */

   /* -------------------------------------------------------------------------- */

   void isSolutionFinite(const auto &water) const {
      for (const auto &p : water.getSurfacePoints()) {
         if (!isFinite(p->phiphin)) {
            std::stringstream ss;
            ss << p->phiphin;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
         }
         for (const auto &[f, phi] : p->phiOnFace)
            if (!isFinite(phi)) {
               std::stringstream ss;
               for (const auto &[f, phi] : p->phiOnFace)
                  ss << phi << ", ";
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
            }
         for (const auto &[f, phin] : p->phinOnFace)
            if (!isFinite(phin)) {
               std::stringstream ss;
               for (const auto &[f, phin] : p->phiOnFace)
                  ss << phin << ", ";
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
            }
      }
   };

   // b! -------------------------------------------------------------------------- */
   // b!                                    solve                                   */
   // b! -------------------------------------------------------------------------- */

   V_d ans;
   V_d tmp_b_RHS;

   int matrix_size = 0;

   void solve(const Buckets<networkPoint *> &FMM_BucketsPoints, const Buckets<networkFace *> &FMM_BucketsFaces) {

      TimeWatch watch;

      setPhiPhinOnFace(WATERS);
      this->matrix_size = setNodeFaceIndices(WATERS);

#if defined(use_lapack)
      std::cout << Red << "   unknown size : " << this->matrix_size << colorReset << std::endl;
      setIGIGn();
      std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;
      generateBIEMatrix();
      ans.resize(knowns.size(), 0.);
      std::cout << "lapack lu decomposition" << std::endl;
      //! A.x = b
      if (this->lu != nullptr) {
         this->lu->init(mat_ukn);
         this->lu->solve(b_RHS /*既知のベクトル（右辺）*/, ans /*解*/);
      } else
         this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/, ans /*解*/, b_RHS /*既知のベクトル（右辺）*/);

      std::cout << colorReset << "update p->phiphin and p->phinOnFace for Dirichlet boundary" << colorReset << std::endl;
      storePhiPhin(WATERS, ans);

      std::cout << Green << "Elapsed time for solving BIE: " << Red << watch() << colorReset << " s\n";

      for (auto water : WATERS)
         isSolutionFinite(*water);

#elif defined(use_gmres)

      /* -------------------------------------------------------------------------- */
      // std::cout << Red << "   unknown size : " << this->matrix_size << colorReset << std::endl;
      // setIGIGn();
      // std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;
      // generateBIEMatrix();  //! generateBIEMatrix()内部でb_RHSも更新
      // ans.resize(knowns.size(), 0.);
      // std::cout << "lapack lu decomposition" << std::endl;
      /* -------------------------------------------------------------------------- */

      for (auto &net : this->WATERS)
         net->setGeometricProperties();

      for (const auto &water : WATERS)
   #pragma omp parallel
         for (const auto &integ_f : water->getSurfaces())
   #pragma omp single nowait
            integ_f->setIntegrationInfo();

      //@ -------------------------------------------------------------------------- */
      //@                       バケットの作成．極の追加 add                             */
      //@ -------------------------------------------------------------------------- */

      auto obj = WATERS[0];
      auto Center = obj->getCenter();
      auto scale = 1.1 * obj->getScale() / 2.;

      CoordinateBounds bounds(Center[0] - scale, Center[0] + scale, Center[1] - scale, Center[1] + scale, Center[2] - scale, Center[2] + scale);
      Buckets<std::shared_ptr<source4FMM<target4FMM>>> B_poles(bounds, obj->getScale() / 2.);
      std::cout << "バケットの作成．極の追加" << std::endl;
      /* phiOnFaceやphinOnFaceを更新することで，getValuesの結果は自動的に更新される． */
      TimeWatch tw;
      for (auto &F : obj->getSurfaces()) {
         const auto [p0, p1, p2] = F->getPoints();
         const auto closest_p_to_origin = p0;
         //! ソースとして重要なのは，G,Gnを作るための情報のみ. quadpointsには，points_facesが入っている
         if (F->isPseudoQuadraticElement) {

            /* -------------------------------------------------------------------------- */
            /*
            DodecaPointsのvertexは，quadpoint, quadpoint_l0, quadpoint_l1, quadpoint_l2
            それに対応する形状関数は，Nc_N0_N1_N2
            として保存されている．
            */

            auto insert = [](auto &insert_to, const auto &points_faces) {
               int i = 0;
               for (auto &[p, f] : points_faces) {
                  auto key = std::get<1>(pf2ID(p, f));
                  insert_to[i++] = {&p->phiOnFace.at(key), &p->phinOnFace.at(key)};  //! updateされる値へのポインタ かつ キー．
               }
            };
            std::array<std::array<std::array<double *, 2>, 6>, 4> pair_pointer_phiphin;
            insert(pair_pointer_phiphin[0], F->dodecaPoints[0]->quadpoint.points_faces);
            insert(pair_pointer_phiphin[1], F->dodecaPoints[0]->quadpoint_l0.points_faces);
            insert(pair_pointer_phiphin[2], F->dodecaPoints[0]->quadpoint_l1.points_faces);
            insert(pair_pointer_phiphin[3], F->dodecaPoints[0]->quadpoint_l2.points_faces);

            // b% [1] ソース点の配置に関する積分
            for (const auto &[t0t1, quadrature_weight, Nc_N0_N1_N2, X, cross, J_det] : F->map_Point_PseudoQuadraticIntegrationInfo_vector[0].at(closest_p_to_origin)) {

               //^  ---------------------------- update values of sources ---------------------------------------------- */

               auto get_weighted_source_densities = [W = (J_det * quadrature_weight), Nc_N0_N1_N2, pair_pointer_phiphin]() -> std::array<double, 2> {
                  std::array<double, 2> phiphin = {0., 0.};
                  for (int i = 0; i < 4; ++i) {
                     for (int j = 0; j < 6; ++j) {
                        phiphin[0] += W * Nc_N0_N1_N2[i][j] * (*pair_pointer_phiphin[i][j][0]);
                        phiphin[1] += W * Nc_N0_N1_N2[i][j] * (*pair_pointer_phiphin[i][j][1]);
                     }
                  }
                  return phiphin;
               };

               //^  ----------------------------- set up direct integration -------------------------------------------- */
               //^ このソース点の直接積分計算をどのように行うかを指定する関数
               //^ 今のところは，面には1ソース点のみが設置されている．そのため，直接積分は１つの面の積分と考えていい

               auto use_this_soruce_when_set_direct_integration = [&obj, F, points = F->Points](const target4FMM *origin) -> std::vector<std::tuple<std::array<double *, 2>, std::array<double, 2>>> {
                  const bool is_near = (Norm((points[0]->X + points[1]->X + points[2]->X) / 3. - origin->Xtarget) < obj->getScale() / 30.);
                  const auto &array_to_use = is_near ? std::span(__array_GW1xGW1__.begin(), __array_GW1xGW1__.end()) : std::span(__array_GW5xGW5__.begin(), __array_GW5xGW5__.end());
                  //! 観測点(target)を特定の番号に位置するようにする．
                  auto dodecapoint = F->dodecaPoints[points[1] == origin ? 1 : (points[2] == origin ? 2 : 0)];

                  // b% [2] ソース点を直接積分に置き換えたあとの，積分（複数点）
                  std::vector<std::array<double, 2>> ij_WGN_WGnN(24, std::array<double, 2>{0., 0.});
                  std::vector<std::array<double *, 2>> ij_phiphin(24, {nullptr, nullptr});
                  std::array<std::array<double, 6>, 4> Nc_N0_N1_N2;
                  double nr, tmp, WG, WGn_, WGn;
                  std::array<double, 3> R, cross, X;
                  for (const auto &[t0, t1, ww] : array_to_use) {
                     auto [xi0, xi1, xi2] = ModTriShape<3>(t0, t1);
                     cross = dodecapoint->cross(xi0, xi1);
                     Nc_N0_N1_N2 = dodecapoint->N6_new(xi0, xi1);
                     nr = Norm(R = (dodecapoint->X(xi0, xi1) - origin->Xtarget));
                     tmp = (ww * (1. - t0) / nr);
                     WG = (Norm(cross) * tmp);
                     WGn_ = (-Dot(R / nr, cross) * tmp / nr);
                     for (int i = 0; i < 4; ++i) {
                        auto quadpoints = dodecapoint->getQuadPoints(i - 1);
                        for (int j = 0; j < 6; ++j) {
                           WGn = WGn_;
                           auto [p, _] = quadpoints->points_faces[j];
                           ij_WGN_WGnN[6 * i + j] += std::array<double, 2>{WG * Nc_N0_N1_N2[i][j], p == origin ? 0. : WGn * Nc_N0_N1_N2[i][j]};  //! リジッドモードテクニック
                           if (i == 0) {
                              auto [p, key] = pf2ID(quadpoints->points_faces[j]);
                              ij_phiphin[6 * i + j] = {&p->phiOnFace.at(key), &p->phinOnFace.at(key)};
                           }
                        }
                     }
                  }

                  // std::unordered_map<std::array<double *, 2>, std::array<double, 2>> uo_results;
                  // uo_results.reserve(24);
                  // for (int k = 0; k < 24; ++k)
                  //    uo_results[ij_phiphin[k]] += ij_WGN_WGnN[k];
                  // std::vector<std::tuple<std::array<double *, 2>, std::array<double, 2>>> results;
                  // results.reserve(24);
                  // for (const auto &[key, value] : uo_results)
                  //    results.push_back({key, value});

                  std::vector<std::tuple<std::array<double *, 2>, std::array<double, 2>>> phiphin_WGN_WGnN(24, {{nullptr, nullptr}, {0., 0.}});
                  for (int k = 0; k < 24; ++k)
                     phiphin_WGN_WGnN[k] = {ij_phiphin[k], ij_WGN_WGnN[k]};

                  return phiphin_WGN_WGnN;
               };

               //^ -------------------------------------------------------------------------- */
               auto pole = std::make_shared<source4FMM<target4FMM>>(X,
                                                                    F->normal,
                                                                    get_weighted_source_densities,
                                                                    use_this_soruce_when_set_direct_integration);
               B_poles.add(X, pole);
            }

         } else {  //@ ==================================================================================================

            auto key0 = std::get<1>(pf2ID(p0, F));
            auto key1 = std::get<1>(pf2ID(p1, F));
            auto key2 = std::get<1>(pf2ID(p2, F));

            std::array<double *, 2> pair_pointer_phiphin0 = {&p0->phiOnFace.at(key0), &p0->phinOnFace.at(key0)};
            std::array<double *, 2> pair_pointer_phiphin1 = {&p1->phiOnFace.at(key1), &p1->phinOnFace.at(key1)};
            std::array<double *, 2> pair_pointer_phiphin2 = {&p2->phiOnFace.at(key2), &p2->phinOnFace.at(key2)};

            // b% [1] ソース点の配置に関する積分
            for (const auto &[t0t1, quadrature_weight, shape3, X, cross, J_det] : F->map_Point_LinearIntegrationInfo_vector[0].at(closest_p_to_origin)) {

               //!  ---------------------------- update values of sources ---------------------------------------------- */

               auto get_weighted_source_densities = [W = (J_det * quadrature_weight),
                                                     shape3,
                                                     pair_pointer_phiphin0,
                                                     pair_pointer_phiphin1,
                                                     pair_pointer_phiphin2]() -> std::array<double, 2> {
                  return {W * (shape3[0] * (*pair_pointer_phiphin0[0]) + shape3[1] * (*pair_pointer_phiphin1[0]) + shape3[2] * (*pair_pointer_phiphin2[0])),
                          W * (shape3[0] * (*pair_pointer_phiphin0[1]) + shape3[1] * (*pair_pointer_phiphin1[1]) + shape3[2] * (*pair_pointer_phiphin2[1]))};
               };

               //!  ----------------------------- set up direct integration -------------------------------------------- */

               //$ このソース点の直接積分計算をどのように行うかを指定する関数
               auto use_this_soruce_when_set_direct_integration = [&obj,
                                                                   F,
                                                                   points = F->Points,
                                                                   pair_pointer_phiphin0,
                                                                   pair_pointer_phiphin1,
                                                                   pair_pointer_phiphin2](const target4FMM *origin) -> std::vector<std::tuple<std::array<double *, 2>, std::array<double, 2>>> {
                  const bool is_near = (Norm((points[0]->X + points[1]->X + points[2]->X) / 3. - origin->Xtarget) < obj->getScale() / 30.);
                  const auto &array_to_use = is_near ? std::span(__array_GW1xGW1__.begin(), __array_GW1xGW1__.end()) : std::span(__array_GW5xGW5__.begin(), __array_GW5xGW5__.end());

                  const int rotate = points[1] == origin ? 1 : (points[2] == origin ? 2 : 0);
                  const auto [p0, p1, p2] = RotateLeft(points, rotate);  //! poleは面に１つのみ．だからこの積分はOK
                  const std::array<std::array<double, 3>, 3> X012 = {p0->Xtarget, p1->Xtarget, p2->Xtarget};
                  std::array<double, 3> cross = Cross(X012[1] - X012[0], X012[2] - X012[0]), WGN = {0., 0., 0.}, WGnN = {0., 0., 0.}, N012_geometry, R;
                  double norm_cross = Norm(cross);
                  double nr, tmp;
                  // b% [2] ソース点を直接積分に置き換えたあとの，積分
                  for (const auto &[t0, t1, ww] : array_to_use) {
                     N012_geometry = ModTriShape<3>(t0, t1);
                     nr = Norm(R = (Dot(N012_geometry, X012) - origin->Xtarget));
                     tmp = (ww * (1. - t0) / nr);
                     WGN += (norm_cross * tmp) * N012_geometry;
                     WGnN += (-Dot(R / nr, cross) * tmp / nr) * N012_geometry;
                  }
                  if (p0 == origin)
                     std::get<0>(WGnN) = 0.;  //! リジッドモードテクニック
                  std::array<int, 3> idx = {0, 1, 2};
                  if (rotate == 1)
                     idx = {2, 0, 1};
                  else if (rotate == 2)
                     idx = {1, 2, 0};
                  return {{pair_pointer_phiphin0, std::array<double, 2>{WGN[idx[0]], WGnN[idx[0]]}},
                          {pair_pointer_phiphin1, std::array<double, 2>{WGN[idx[1]], WGnN[idx[1]]}},
                          {pair_pointer_phiphin2, std::array<double, 2>{WGN[idx[2]], WGnN[idx[2]]}}};
               };

               auto pole = std::make_shared<source4FMM<target4FMM>>(X, F->normal,
                                                                    get_weighted_source_densities,               //! FMMで使うソースの値を，GMRESの反復計算でアップデートするための関数
                                                                    use_this_soruce_when_set_direct_integration  //! 　このソース点１点に対する直接積分計算をどのように行うかを指定する関数
               );
               B_poles.add(X, pole);
            }
         }
         /* -------------------------------------------------------------------------- */
      }

      std::cout << Magenta << "Add poles" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

      //@ -------------------------------------------------------------------------- */
      //@                                ツリー構造を生成                               */
      //@ -------------------------------------------------------------------------- */

      std::cout << "ツリー構造を生成" << std::endl;
      int max_level = 10;
      B_poles.setLevel(0, max_level);
      B_poles.generateTree([](auto bucket) {
         if (bucket->all_stored_objects_vector.empty())
            return false;
         else
            return bucket->all_stored_objects_vector.size() > 500 && bucket->level < bucket->max_level;
      });
      std::cout << Magenta << "Tree" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

      /* -------------------------------------------------------------------------- */
      if (false) {
         {
            T8Tddd t8tddd = (CoordinateBounds)(obj->bounds);
            std::ofstream ofs("./obj_bounds.vtp");
            vtkPolygonWrite(ofs, t8tddd);
            ofs.close();
         }

         std::vector<T8Tddd> cube_level_deepest;
         B_poles.forEachAtDeepest([&](auto *B) {
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            cube_level_deepest.push_back(t8tddd);
         });
         std::ofstream ofs("./cube_level_deepest_points.vtp");
         vtkPolygonWrite(ofs, cube_level_deepest);
         ofs.close();

         // レベルごとに出力したい場合
         for (int i = 0; i <= B_poles.max_level; ++i) {
            std::vector<T8Tddd> cube_level;
            B_poles.forEachAtLevel(i, [&](auto *B) {
               T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
               cube_level.push_back(t8tddd);
            });
            std::ofstream ofs("./cube_level_points_" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, cube_level);
            ofs.close();
         }
      }
      /* -------------------------------------------------------------------------- */

      //@ -------------------------------------------------------------------------- */
      //@                                  FMM                                       */
      //@ -------------------------------------------------------------------------- */
      auto points = obj->getSurfacePoints();

      initializeFMM(B_poles, points);

      /*
         基本とする形，左辺
         {{G0,G0,G0,G0},{G1,G1,G1,G1}}.{phin,phin,phin,phin} - {{Gn0,Gn0,Gn0,Gn0},{Gn1,Gn1,Gn1,Gn1}}.{phi,phi,phi,phi}
         node1がNeumanの場合，他はDirichletの場合
         左辺，未知変数側：
         {{G0,G0,G0,G0},{G1,G1,G1,G1}}.{phin,0,phin,phin} - {{Gn0,Gn0,Gn0,Gn0},{Gn1,Gn1,Gn1,Gn1}}.{0,phi,0,0}
         右辺，既知変数側：
         - {{G0,G0,G0,G0},{G1,G1,G1,G1}}.{0, phin, 0, 0} + {{Gn0,Gn0,Gn0,Gn0},{Gn1,Gn1,Gn1,Gn1}}.{phi, 0, phi, phi}
      */

      std::vector<std::array<double *, 2>> copy_map;
      copy_map.reserve(points.size() * 2 * 5);
      for (const auto &p : points) {
         p->phiOnFace_copy = p->phiOnFace;
         p->phinOnFace_copy = p->phinOnFace;
         for (const auto &[f, i] : p->f2Index) {
            copy_map.push_back({&p->phiOnFace_copy[f], &p->phiOnFace[f]});
            copy_map.push_back({&p->phinOnFace_copy[f], &p->phinOnFace[f]});
         }
      }

      auto copyPhiPhin = [&](Network *obj) -> void {
         for (auto &[p_copy, value] : copy_map) {
            *p_copy = *value;
            *value = 0.;
         }
      };

      auto restorePhiPhin = [&](Network *obj) -> void {
         for (auto &[p_copy, value] : copy_map)
            *value = *p_copy;
      };

      /* --------------------------------- 対角項を計算 --------------------------------- */

      {
         std::cout << "calculate diagonal coefficient" << std::endl;
         copyPhiPhin(obj);
   #pragma omp parallel
         for (const auto &p : points) {
   #pragma omp single nowait
            for (const auto &[f, i] : p->f2Index) {
               if (isDirichletID_BEM(p, f)) {
                  p->phinOnFace[f] = 0.;
                  p->phiOnFace[f] = 1.;
               } else if (isNeumannID_BEM(p, f)) {
                  p->phinOnFace[f] = 0.;
                  p->phiOnFace[f] = 1.;
               } else
                  throw std::runtime_error("Error: Boundary type is not defined.");
            }
         }
         updateFMM(B_poles);
         TimeWatch tw;
   #pragma omp parallel
         for (const auto &a : points)
   #pragma omp single nowait
         {
            for (const auto &[f, i] : a->f2Index) {
               auto GPhin_GnPhi_near = a->integrate();
               auto GPhin_GnPhi_far = a->integrateFMM();
               auto [GPhin, GnPhi] = GPhin_GnPhi_near + GPhin_GnPhi_far;
               a->diagonal_coefficient = GPhin - GnPhi;
               // std::cout << "diagonal coefficient: " << a->diagonal_coefficient << ", " << a->getSolidAngle() << std::endl;
            }
         }
         std::cout << Magenta << "integration (direct integration is dominant)" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
         restorePhiPhin(obj);
      }

      /* -------------------------------------------------------------------------- */

      std::size_t C = 0;
      for (const auto &p : points)
         C += p->f2Index.size();
      std::vector<double> GMRES_VECTOR(C, 0.);

      /* ---------------------------------------------------------------------------- */

      TimeWatch timewatch;
      int counter = 0;
      auto AdotV_FMM = [&counter, &timewatch, &obj, &B_poles, &GMRES_VECTOR, &points]() -> V_d {
         updateFMM(B_poles);
         auto elp = timewatch();
   #pragma omp parallel
         for (const auto &a : points)
   #pragma omp single nowait
         {
            std::array<double, 2> GPhin_GnPhi;
            for (const auto &[f, i] : a->f2Index) {
               if (a->CORNER && isNeumannID_BEM(a, f)) {
                  GMRES_VECTOR[i] = a->phiOnFace[f];
               } else {
                  GPhin_GnPhi = (a->integrate() + a->integrateFMM()) / a->diagonal_coefficient;
                  GMRES_VECTOR[i] = std::get<0>(GPhin_GnPhi) - (std::get<1>(GPhin_GnPhi) + a->phiOnFace[f]);  //! 基本形
               }
            }
         }
         if (counter++ % 10 == 0)
            std::cout << "Elapsed time : " << Magenta << "updateFMM" << elp << Green << "  integration" << timewatch() << " counter =" << counter << colorReset << std::endl;
         return GMRES_VECTOR;
      };

      /* ---------------------- GMRESで利用する関数を定義する． ------------------------- */

      std::vector<double *> GMRES_input(C);
      for (const auto &p : points) {
         for (const auto &[f, i] : p->f2Index) {
            if (isDirichletID_BEM(p, f)) {
               GMRES_input[i] = &p->phinOnFace[f];
            } else if (isNeumannID_BEM(p, f)) {
               GMRES_input[i] = &p->phiOnFace[f];
            } else
               throw std::runtime_error("Error: Boundary type is not defined.");
         }
      }

      /* ------------------------ GMRES用　return_A_dot_v　--------------------------- */

      TimeWatch watch_;
      auto return_A_dot_v = [&](const V_d &V) -> V_d {
         copyPhiPhin(obj);
         for (int i = 0; i < V.size(); ++i)
            *GMRES_input[i] = V[i];
         AdotV_FMM();  //% 代入した値を使ってFMMでA・Vを計算する
         restorePhiPhin(obj);
         return GMRES_VECTOR;
      };

      /* ---------------------------------- bの計算 ---------------------------------- */

      // 早くする方法は，
      // おそらく
      // ＊　ガウス点の削減（近くは多く，遠くは少なく．がFMMではできない！なにか方法があるはず，数式から考える）
      // ＊　展開工数の削減（遠くは少なく，近くは多く展開？）
      // ＊　プログラムの最適化
      // ＊　展開の収束性の向上（もっと早く収束する展開は存在するのか？）

      //! 既知変数側
      copyPhiPhin(obj);
   #pragma omp parallel
      for (const auto &p : points)
   #pragma omp single nowait
      {
         for (const auto &[f, i] : p->f2Index) {
            if (isDirichletID_BEM(p, f)) {
               // p->phinOnFace[f] = 0.;
               p->phiOnFace[f] = -p->phiOnFace_copy.at(f);
            } else if (isNeumannID_BEM(p, f)) {
               p->phinOnFace[f] = -p->phinOnFace_copy.at(f);
               // p->phiOnFace[f] = 0.;
            } else
               throw std::runtime_error("Error: Boundary type is not defined.");
         }
      }
      //% 代入した値を使ってFMMでA・Vを計算する
      auto b = AdotV_FMM();

      /* ------------------------- Compativility Condition ------------------------ */
   #pragma omp parallel
      for (const auto &p : points)
   #pragma omp single nowait
      {
         if (p->CORNER)
            for (const auto &[f, i] : p->f2Index)
               if (isNeumannID_BEM(p, f))
                  b[i] = p->phiphin[0];
      }
      /* -------------------------------------------------------------------------- */

      //! compatibility condition for p on Gamma^N face and corner
      //! p->CORNER && isNeumannID_BEM(p, f)

      restorePhiPhin(obj);

      int count = 0, n = 0;
      for (const auto &p : obj->getPoints())
         p->b_diff_RHS_FMM = 0.;

      for (const auto &p : obj->getPoints()) {
         for (const auto &[f, i] : p->f2Index) {
            // p->b_RHS_Direct = b_RHS[i];
            p->b_RHS_FMM = b[i];
            // auto diff = std::abs(b[i] - b_RHS[i]);
            // if (p->b_diff_RHS_FMM < diff)
            //    p->b_diff_RHS_FMM = b[i] - b_RHS[i];
         }
      }

      /* --------------------------------- GMRES ------------------------------------- */

      // std::cout << "use gmres" << std::endl;
      std::vector<int> list = {20, 40, 60, 80, 100, 120, 150, 200};
      std::vector<double> error;
      std::vector<double> x0(b.size(), 0.);

      for (const auto &p : obj->getPoints())
         for (const auto &[f, i] : p->f2Index) {
            if (isDirichletID_BEM(p, f))
               x0[i] = p->phinOnFace.at(f);
            else if (isNeumannID_BEM(p, f)) {
               x0[i] = p->phiOnFace.at(f);
            }
         }

      const double torrelance = 1E-13 * obj->getPoints().size();
      int c = 0;
      for (auto gmres_size : list) {
         TimeWatch watch;
         gmres *GMRES = new gmres(return_A_dot_v, b, x0, gmres_size);
         std::cout << "Elapsed time for GMRES: " << watch() << " s" << std::endl;
         std::cout << "gmres size = " << gmres_size << std::endl;
         std::cout << "gmres error = " << GMRES->err << std::endl;
         error.push_back(GMRES->err);
         x0 = GMRES->x;
         //
         auto ANS = [&]() -> V_d {
            copyPhiPhin(obj);
            for (int i = 0; i < x0.size(); ++i)
               *GMRES_input[i] = x0[i];
            auto ret = AdotV_FMM();
            restorePhiPhin(obj);
            return ret;
         };
         double E = Norm(ANS() - b);
         //
         std::cout << "GMRES true error = " << E << std::endl;
         std::cout << "count = " << c++ << std::endl;
         delete GMRES;
         if (E < 1E-10)
            break;
      }

      std::cout << "gmres size list = " << list << std::endl;
      std::cout << "gmres error = " << error << std::endl;

      std::cout << colorReset << "update p->phiphin and p->phinOnFace for Dirichlet boundary" << colorReset << std::endl;
      storePhiPhin(WATERS, x0);

      std::cout << Green << "Elapsed time for solving BIE: " << Red << watch() << colorReset << " s\n";

      for (auto water : WATERS)
         isSolutionFinite(*water);

#endif
   };

   // b! ------------------------------------------------------------------------------ */
   // b!                             solve phi_t and phi_n_t                            */
   // b! ------------------------------------------------------------------------------ */

   /*DOC_EXTRACT 0_4_0_1_FLOATING_BODY_SIMULATION

   ## 浮体動揺解析

   BEM-MELで浮体動揺解析ができるようにするのは簡単ではない．
   浮体に掛かる圧力の計算に必要な$\phi_t$が簡単には求まらないためである．
   これに関しては，\cite{Wu2003}や\cite{Ma2009a}が参考になる．

   ### 浮体の運動方程式

   <img src="schematic_float.png" width="400px" />

   浮体の重心の運動方程式：

   ```math
   m \frac{d {\boldsymbol U}_{\rm c}}{d t} = \boldsymbol{F}_{\text {ext }}+\boldsymbol{F}_{\text {hydro }}, \quad
   \boldsymbol{I} \frac{d {\boldsymbol \Omega}_{\rm c}}{d t} = \boldsymbol{T}_{\text {ext }}+\boldsymbol{T}_{\text {hydro }}
   ```

   NOTE:これらの変数は固定されたグローバル座標系上での変数である．大抵の場合$`m`$は座標系によらないが，$`\boldsymbol{I}`$は，浮体の基本姿勢において定義されたものであり，（浮体の座標系においては変化しないがグローバル座標系においては）浮体の姿勢によって変化する．

   $`{\boldsymbol U}_{\rm c}`$は浮体の移動速度．
   $`\boldsymbol{F}_{\text {ext }}`$は重力などの外力，$`\boldsymbol{F}_{\text {hydro }}`$は水の力，$`\boldsymbol{T}_{\text {ext }}`$は外力によるトルク，$`\boldsymbol{T}_{\text {hydro }}`$は水の力によるトルク．
   浮体が流体から受ける力$`\boldsymbol{F}_{\text {hydro }}`$は，浮体表面の圧力$`p`$を積分することで得られ，
   また圧力$`p`$は速度ポテンシャル$`\phi`$を用いて，以下のように書ける．

   \ref{BEM:surfaceIntegralOfPressure}{圧力積分}と
   \ref{BEM:surfaceIntegralOfTorque}{トルクの積分}：

   ```math
   \boldsymbol{F}_{\text {hydro }}=\iint_{\Gamma_{\rm float}} p\boldsymbol{n}  d S, \quad
   \boldsymbol{T}_{\text {hydro }}=\iint_{\Gamma_{\rm float}} ({\bf x}-{\bf x}_{\rm c})\times (p\boldsymbol{n})  d S, \quad
   p= p({\bf x}) =-\rho\left(\frac{\partial \phi}{\partial t}+\frac{1}{2} \|\nabla \phi\|^{2}+g z\right)
   ```

   $`\frac{\partial \phi}{\partial t}`$を$`\phi_t`$と書くことにする．この$`\phi_t`$は陽には求められない．
   そこで，$`\phi`$と似た方法，BIEを使った方法で$`\phi_t`$を求める．$`\phi`$と$`\phi_n`$の間に成り立つ境界積分方程式と全く同じ式が，$`\phi_t`$と$`\phi_{nt}`$の間にも成り立つ：

   ```math
   \alpha ({\bf{a}})\phi_t ({\bf{a}}) = \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi_t ({\bf{x}}) - \phi_t ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
   \quad\text{on}\quad{\bf x} \in \Gamma(t).
   ```

   */

   /*DOC_EXTRACT 0_4_1_FLOATING_BODY_SIMULATION

   ### 加速度の計算の難しさ

   これは，浮体表面の圧力の計算の困難，もっと言えば$`\phi_t`$の計算の困難に起因する． \cite{Ma2009}によると，$`\phi_t`$の計算方法として以下の４つの方法が提案されている．

   1. 間接的法 (indirect method) ：補助関数を使う方法
   2. モード分解法 (mode-decomposition method)
   3. Dalena and Tanizawa's method
   4. Caoの反復法 (iterative method) \cite{Cao1994}
   5. Maの方法 \cite{Ma2009}

   #### 間接的法，モード分解法

   間接的法，モード分解法は，$`\phi`$に関するBVPと似ているが異なる新たなBVPを解く必要がある．
   この新たなBIEの境界条件は違うが，係数行列は同じ（私は違うと思うのだが）らしい．
   ただ悩ましいので，LU分解のような直接法なら逆行列を保持するので余計な計算が発生しないが，直接法はそもそも遅い．
   そこで反復法を使いたいが，反復法は逆行列を保持しないので，毎回係数行列を計算する必要がある，というジレンマがある．

   #### Dalena and Tanizawa's method

   Dalena and Tanizawa's methodは，$`\phi`$に関するBVPと全く違うBVPを解く必要があり，係数行列も違うので，新たに行列を構成する必要がある．
   \cite{Feng2017}によると，この方法は境界面の局所的な曲率を用いる必要があるため，３次元解析に利用することがとても難しい．

   #### 反復法

   Caoの反復法は，新たなBVPを解く必要がないので，上の問題はないらしい\cite{Ma2009}．
   （これは間違いで，直接法を使った場合はそうだが，反復法（GMRESのような）を使うなら，このCaoの反復法の内部で反復法（GMRESなど）をする必要があり時間がかかり，
   初めの２つに優っているとは言えない．同じ程度の時間がかかる．）

   #### Maの反復法

   | Method |  |
   |:---:|:---:|
   | Indirect method | 浮体１つに対して６つ，新しいBIEを立てる．新たに解く必要があり遅い |
   | Mode-decomposition method | 浮体１つに対して7つ，新しいBIEを立てる．新たに解く必要があり遅い |
   | Dalena and Tanizawa's method | ? |
   | Cao's iterative method | 直接法で解くなら同じBIE係数行列を使えるので，速い．反復法なら，反復法の内部で反復法をするので遅い．|
   | Ma's iterative method | 直接法で解くなら同じBIE係数行列を使えるので，速い．反復法なら，反復法の内部で反復法をするので遅い．|

   ### $`\phi_t`$と$`\phi_{nt}`$に関するBIEの解き方（と$`\phi_{nt}`$の与え方）

   $`\phi_t`$と$`\phi_{nt}`$に関するBIEを解くためには，ディリクレ境界には$`\phi_t`$を，ノイマン境界には$`\phi_{nt}`$を与える．

   #### ディリクレ節点の$`\phi_{nt}`$の与え方(水面：圧力が既知，$`\phi`$が既知)

   このディリクレ境界では，圧力が与えられていないので，このBiEにおいては，ノイマン境界条件を与える．
   ただし，壁が完全に固定されている場合，$`\phi_{nt}`$は0とする．

   #### ディリクレ節点の$`\phi_{t}`$の与え方($`\phi`$を与える造波装置：圧力が未知，$`\phi`$が既知)

   ディリクレ境界では$`\phi_t`$は，圧力が大気圧と決まっているので，ベルヌーイの圧力方程式から$`\phi_t`$を求めることができる．

   #### ノイマン節点での$`\phi_{nt}`$の与え方

   境界面が静止しているかどうかに関わらず，流体と物体との境界では，境界法線方向速度が一致する．
   浮体重心$`{\bf x}_c`$から境界面上の点$\bf x$までの位置ベクトルを$`\boldsymbol r = {\bf x} - {\bf x}_c`$とする．
   表面上のある点の移動速度$`\frac{d\boldsymbol r}{dt}`$と流体粒子の流速$`\nabla \phi`$の間には，次の境界条件が成り立つ．

   ```math
   {\bf n}\cdot\frac{d\boldsymbol r}{dt} =  {\bf n} \cdot \nabla \phi,\quad \frac{d\boldsymbol r}{dt} = \boldsymbol U_{\rm c} + {\boldsymbol \Omega}_{\rm c} \times \boldsymbol r
   ```

   物体上のある点ではこれが常に成り立つ．

   これを微分することで，$`\phi_{nt}`$を$`\phi`$と加速度$`\frac{d{\boldsymbol U}_{\rm c}}{dt}`$と角加速度$`\frac{d{\boldsymbol \Omega}_{\rm c}}{dt}`$を使って表すことができる．
   \cite{Wu1998}

   ```math
   \begin{aligned}
   &\rightarrow& 0& =\frac{d}{dt}\left({\bf n}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)\right) \\
   &\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \frac{d}{dt}\left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)\\
   &\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2}-\left(\frac{\partial}{\partial t}+\frac{d{\boldsymbol r}}{dt}\cdot\nabla\right)\nabla \phi\right)\\
   &\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2}- {\nabla \phi_t - \left(\frac{d\boldsymbol r}{dt} \cdot \nabla\right)\nabla \phi}\right)\\
   &\rightarrow& \phi_{nt}& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2} - \frac{d\boldsymbol r}{dt} \cdot (\nabla\otimes\nabla \phi) \right)
   \end{aligned}
   ```

   ここの$`\frac{d{\bf n}}{dt}`$と$`\frac{d^2\boldsymbol r}{dt^2}`$は，$`{\boldsymbol U}_{\rm c}`$と$`\boldsymbol \Omega_{\rm c}`$を用いて，

   ```math
   \frac{d^2\boldsymbol r}{dt^2}
   = \frac{d}{dt}\left({\boldsymbol U}_{\rm c} + \boldsymbol \Omega_{\rm c} \times \boldsymbol r\right)
   = \frac{d{\boldsymbol U}_{\rm c}}{dt} + \frac{d{\boldsymbol \Omega_{\rm c}}}{dt} \times \boldsymbol r + \boldsymbol \Omega_{\rm c} \times \frac{d\boldsymbol r}{dt}
   ,\quad \frac{d{\bf n}}{dt} = {\boldsymbol \Omega}_{\rm c}\times{\bf n}
   ```

   $`\frac{d \boldsymbol r}{dt}`$は\ref{velocityRigidBody}{`velocityRigidBody`}
   $`\frac{d^2 \boldsymbol r}{dt^2}`$は\ref{accelRigidBody}{`accelRigidBody`}で計算する．

   \ref{BEM:phint_Neumann}{`phin_Neuamnn`}で$`\phi_{nt}`$を計算する．これは\ref{BEM:setPhiPhin_t}{`setPhiPhin_t`}で使っている．

   $`\frac{d^2\boldsymbol r}{dt^2}`$を上の式に代入し，$`\phi_{nt}`$を求め，
   次にBIEから$`\phi_t`$を求め，次に圧力$p$を求める．
   そして，浮体の重さと慣性モーメントを考慮して圧力から求めた$`\frac{d^2\boldsymbol r}{dt^2}`$は，
   入力した$`\frac{d^2\boldsymbol r}{dt^2}`$と一致しなければならない．

   現状を整理すると，この浮体動揺解析において，知りたい未知変数は，浮体の加速度と角加速度だけ．
   しかし，浮体の没水面上にある節点での圧力$`p`$が得られないと，$`\boldsymbol{F}_{\text {hydro }}`$が得られず，運動方程式から浮体加速度が計算できない．
   圧力を計算するためには，$`\phi_t`$が必要で，$`\phi_t`$は簡単には得られない，という状況．

   物体の加速度は， 節点における$`\{\phi_{nt0},\phi_{nt1},\phi_{nt2},..\} = \Phi_{nt}`$が分かれば求まるが，
   逆に$`\phi_{nt}`$は$`\frac{d\boldsymbol U_{\rm c}}{dt}`$と$\frac{d {\boldsymbol \Omega} _{\rm c}}{d t}$が分かれば求まる．また，物体の角加速度に関しても同様である．

   ```math
   m \frac{d\boldsymbol U_{\rm c}}{dt} = \boldsymbol{F} _{\text {ext }}+ F_{\text {hydro}}\left(\Phi_{nt}\left(\frac{d\boldsymbol U_{\rm c}}{dt},\frac{d {\boldsymbol \Omega} _{\rm c}}{d t}\right)\right),\quad
   \boldsymbol{I} \frac{d {\boldsymbol \Omega} _{\rm c}}{d t} = \boldsymbol{T} _{\text {ext }}+\boldsymbol{T} _{\text {hydro }}\left(\Phi_{nt}\left(\frac{d\boldsymbol U_{\rm c}}{dt},\frac{d {\boldsymbol \Omega} _{\rm c}}{d t}\right)\right)
   ```

   これを満たすように，$`\Phi_{nt}`$を求める．これは次のように書き換えて，根探し問題として解く．
   このプログラムでは，\ref{quasi_newton:broyden}{Broyden法}を使って，根探している．

   ```math
   \boldsymbol{0} = m \frac{d\boldsymbol U_{\rm c}}{dt} - \boldsymbol{F} _{\text {ext }} - F_{\text {hydro}}\left(\Phi_{nt}\left(\frac{d\boldsymbol U_{\rm c}}{dt},\frac{d {\boldsymbol \Omega} _{\rm c}}{d t}\right)\right),\quad
   \boldsymbol{0} = \boldsymbol{I} \frac{d {\boldsymbol \Omega} _{\rm c}}{d t} - \boldsymbol{T} _{\text {ext }} - \boldsymbol{T} _{\text {hydro }}\left(\Phi_{nt}\left(\frac{d\boldsymbol U_{\rm c}}{dt},\frac{d {\boldsymbol \Omega} _{\rm c}}{d t} \right)\right)
   ```

   この式を，$`{\boldsymbol Q}\left(\dfrac{d {\boldsymbol U} _{\rm c}}{d t}, \dfrac{d {\boldsymbol \Omega} _{\rm c}}{d t}\right)=(0,0,0,0,0,0)`$
   として，これを満たすような$`\dfrac{d {\boldsymbol U} _{\rm c}}{d t}`$と$`\dfrac{d {\boldsymbol \Omega} _{\rm c}}{d t}`$を求める．
   $`\phi_{nt}`$はこれを満たした$`\dfrac{d {\boldsymbol U} _{\rm c}}{d t}`$と$`\dfrac{d {\boldsymbol \Omega} _{\rm c}}{d t}`$を用いて求める．

   $`\phi_{nt}`$は，\ref{BEM:setphint}{ここ}で与えている．

   この方法は，基本的には\cite{Cao1994}と同じ方法である．

   */

   /*DOC_EXTRACT 0_4_2_FLOATING_BODY_SIMULATION

   ### 補助関数を使った方法

   浮体動揺解析で問題となったのは，圧力の計算に使う$`\phi_t\,{\rm on}\,🚢`$が簡単には求まらないことであったが，
   $`\iint_{\Gamma_{🚢}} \phi_t{\bf n}dS`$と$`\iint_{\Gamma_{🚢}}\phi_{t}({\bf x}-{\bf x}_c)\times{\bf n}dS`$がわかればある場所の圧力はわからないが，
   🚢にかかる力は計算できるのでそれでも問題ない．

   体積積分がゼロとなるように，領域内でラプラス方程式を満たすような$`\varphi`$，
   そして$`\Gamma _{🚢}`$上ではこちらが望む$`\varphi_n`$となり，また$`\Gamma \rm other`$上では$`\varphi=0`$となる
   そんな$`\varphi`$をBIEを使って計算する．この$`\varphi`$を使うと次の式が成り立つ．
   （NOTE：境界上の全ての節点上で$`\varphi`$と$`\varphi_n`$が求まったとする）

   ```math
   \begin{align*}
   0 &= \iint _\Gamma {\left( {\varphi\nabla {\phi_t} ({\bf{x}}) - {\phi_t} ({\bf{x}})\nabla \varphi} \right) \cdot {\bf{n}}({\bf{x}})dS}\\
   \rightarrow 0 &= \iint _{\Gamma _{🚢}+\Gamma _{🌊}+\Gamma _{\rm wall}} \varphi {\phi_{nt}} dS - \iint _{\Gamma _{🚢}+\Gamma _{🌊}+\Gamma _{\rm wall}} {\phi_t} \varphi_n dS\\
   \rightarrow 0 &= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} \varphi {\phi_{nt}} dS - \iint _{\Gamma _{🚢}+\Gamma _{🌊}} {\phi_t} \varphi_n dS\\
   \rightarrow \iint _{\Gamma _{🚢}} {\phi_t} \varphi_n dS &= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} \varphi {\phi_{nt}} dS - \iint _{\Gamma _{🌊}} {\phi_t} \varphi_n dS\\
   \rightarrow \iint_{\Gamma_{🚢}} \phi_t
   \begin{bmatrix}
   \boldsymbol{n} \\
   (\boldsymbol{x} - \boldsymbol{x}_c) \times \boldsymbol{n}
   \end{bmatrix} dS
   &= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} {\boldsymbol{\varphi}_{1-6}} {\phi_{nt}} dS - \iint _{\Gamma _{🌊}} {\phi_t} {\boldsymbol{\varphi}_n}_{1-6} dS\\
   \end{align*}
   ```

   つまり，$`\varphi_n`$を適当に選べば，左辺は知りたかった積分となり，右辺の積分で計算できることになる．

   もし浮体がもう一つあると

   ```math
   \begin{align*}
   \iint_{\Gamma_{🚢}} \phi_t
   \begin{bmatrix}
   \boldsymbol{n} \\
   (\boldsymbol{x} - \boldsymbol{x}_c) \times \boldsymbol{n}
   \end{bmatrix} dS
   & = \iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi}_{1-6}} {\phi_{nt}} dS - \iint _{\Gamma_{🚤}+\Gamma _{🌊}} {\phi_t} {\boldsymbol{\varphi}_n}_{1-6} dS\\
   \rightarrow \iint_{\Gamma_{🚢}} \phi_t
   \begin{bmatrix}
   \boldsymbol{n} \\
   (\boldsymbol{x} - \boldsymbol{x}_c) \times \boldsymbol{n}
   \end{bmatrix} dS
   & = \iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi}_{1-6}} {\phi_{nt}} dS - \iint _{\Gamma _{🌊}} {\phi_t} {\boldsymbol{\varphi}_n}_{1-6} dS
   \end{align*}
   ```

   同じように

   ```math
   \begin{align*}
   \iint_{\Gamma_{🚤}} \phi_t
   \begin{bmatrix}
   \boldsymbol{n} \\
   (\boldsymbol{x} - \boldsymbol{x}_c) \times \boldsymbol{n}
   \end{bmatrix} dS
   & = \iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi}_{7-12}} {\phi_{nt}} dS - \iint _{\Gamma _{🌊}} {\phi_t} {\boldsymbol{\varphi}_n}_{7-12} dS
   \end{align*}
   ```

   $`\iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi}_{1-6}} {\phi_{nt}} dS`$や
   $`\iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi}_{7-12}} {\phi_{nt}} dS`$
   は加速度行列とある既知変数から成る行列の積で表される．こうして，運動方程式の$`\boldsymbol{F}_{\text {hydro }}`$と$`\boldsymbol{T}_{\text {hydro }}`$を加速度によって表すことができ，
   運動方程式は加速度だけに関する連立方程式となる．

   この方法は，\cite{Wu1996}，\cite{Kashiwagi2000}，\cite{Wu2003}で使用されている．
   この方法は，複数の浮体を考えていないが，\cite{Feng2017}はこれを基にして２浮体の場合でも動揺解析を行っている．

   */

   /* ------------------------------------------------------ */

   V_d initializeAcceleration(const std::vector<Network *> &rigidbodies) {
      V_d ACCELS_init;
      for (const auto &net : rigidbodies) {
         // if (net->interp_accel.size() > 3) {
         //    std::cout << Red << "interp_accel" << colorReset << std::endl;
         //    std::ranges::for_each(net->interp_accel(net->RK_Q.get_t()), [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
         // } else
         std::ranges::for_each(net->acceleration, [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
      }
      return ACCELS_init;
   }

   void insertAcceleration(const std::vector<Network *> &rigidbodies, const V_d &BM_X) {
      int i = 0;
      for (const auto &net : rigidbodies) {
         if (net->isFloatingBody) {
            double start_time = 0;
            if (net->inputJSON.at("velocity").size() > 1)
               start_time = std::stod(net->inputJSON.at("velocity")[1]);
            if (simulation_time < start_time)
               std::ranges::for_each(net->acceleration, [&](auto &a_w) { i++; });
            else
               std::ranges::for_each(net->acceleration, [&](auto &a_w) { a_w = BM_X[i++]; });
         } else
            std::ranges::for_each(net->acceleration, [&](auto &a_w) { i++;/*a_w = BM_X[i++];*/ });
      }
   }

   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   V_d Func(V_d ACCELS_IN, const std::vector<Network *> WATERS, const std::vector<Network *> &rigidbodies, V_d &ACCELS_OUT) {
      TimeWatch watch;
      auto ACCELS = ACCELS_IN;
      //* --------------------------------------------------- */
      //*                  加速度 --> phiphin_t                */
      //* --------------------------------------------------- */
      setPhiPhin_t(WATERS);
      knowns.resize(this->matrix_size);
      for (const auto water : WATERS)
#pragma omp parallel
         for (const auto &p : water->getSurfacePoints())
#pragma omp single nowait
            for (const auto &[f, i] : p->f2Index) {
               if (isDirichletID_BEM(p, f))
                  knowns[i] = p->phitOnFace.at(f);
               else if (isNeumannID_BEM(p, f))
                  knowns[i] = p->phintOnFace.at(f);
            }

      std::cout << "knowns" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      ans.resize(knowns.size(), 0);
      tmp_b_RHS.resize(knowns.size(), 0);
#pragma omp parallel for
      for (size_t i = 0; i < mat_kn.size(); ++i)
         tmp_b_RHS[i] = Dot(mat_kn[i], knowns);

      std::cout << "tmp_b_RHS" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      if (this->lu != nullptr)
         this->lu->solve(tmp_b_RHS /*既知のベクトル（右辺）*/, ans /*解*/);

      std::cout << "solved" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      //@ -------------------------------------------------------------------------- */
      //@                    update p->phiphin_t and p->phinOnFace                   */
      //@ -------------------------------------------------------------------------- */

      storePhiPhin_t(WATERS, ans);
      std::cout << Green << "storePhiPhin_t" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

      //* --------------------------------------------------- */
      //*                 phiphin_t --> 圧力                   */
      //* --------------------------------------------------- */

      for (const auto water : WATERS)
         for (const auto &p : water->getSurfacePoints())
            for (const auto &[f, i] : p->f2Index) {
               if (isDirichletID_BEM(p, f))
                  p->pressure = p->pressure_BEM = 0;
               else
                  p->pressure = p->pressure_BEM = -_WATER_DENSITY_ * (std::get<0>(p->phiphin_t) + 0.5 * Dot(p->U_BEM, p->U_BEM) + _GRAVITY_ * p->height());
            }

      std::cout << Green << "pressure" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

      //* --------------------------------------------------- */
      //*              圧力 ---> 力 --> 加速度                  */
      //* --------------------------------------------------- */
      /*DOC_EXTRACT 0_4_0_2_FLOATING_BODY_SIMULATION

      実際の実験では，浮体のある基本的な姿勢における主慣性モーメントが与えられる．$`{\boldsymbol I}`$を主慣性モーメントテンソルとする．

      ```math
      {\boldsymbol I} = \begin{pmatrix}
      I_x & 0 & 0 \\
      0 & I_y & 0 \\
      0 & 0 & I_z
      \end{pmatrix}
      ```

      global座標における浮体の慣性モーメントテンソルを求めるには，次のように考えればいい．

      ```math
      \begin{aligned}
      {\boldsymbol I}\frac{d{\bf \Omega}_{\rm L}}{dt} &= {\bf T}_{\rm L}\\
      {\boldsymbol I}{\rm R}_{g2l} \frac{d{\bf \Omega}_{\rm G}}{dt}& = {\rm R}_{g2l}{\bf T}_{\rm G}\\
      {\rm R}_{g2l}^{-1}{\boldsymbol I}{\rm R}_{g2l} \frac{d{\bf \Omega}_{\rm G}}{dt}& = {\bf T}_{\rm G}\\
      \end{aligned}
      ```

      このことから，global座標における慣性モーメントテンソルは，次のようになる．

      ```math
      {\boldsymbol I}_{\rm G} = {\rm R}_{g2l}^{-1}{\boldsymbol I}{\rm R}_{g2l}
      ```

      $R^{-1}$が$R^{\top}$であることと，$I^{-1}$が対角成分のみの行列であることを利用すれば，次のように書ける．

      ```math
      \begin{align*}
      {\left({{I}_{G}}\right)}_{il}&={\left({{R}^{-{1}}}\right)}_{ij}{I}_{jk}{R}_{kl}={R}_{ji}{I}_{jj}{R}_{jl}={R}_{1i}{R}_{1l}{I}_{x}+{R}_{2i}{R}_{2l}{I}_{y}+{R}_{3i}{R}_{3l}{I}_{z}\\
      {\left({{I}_{G}^{-{1}}}\right)}_{il}&={\left({{R}^{-{1}}}\right)}_{ij}{\left({{I}^{-{1}}}\right)}_{jk}{R}_{kl}={R}_{ji}{\left({{I}^{-{1}}}\right)}_{jj}{R}_{jl}=\frac{{R}_{1i}{R}_{1l}}{{I}_{x}}+\frac{{R}_{2i}{R}_{2l}}{{I}_{y}}+\frac{{R}_{3i}{R}_{3l}}{{I}_{z}}
      \end{align*}
      ```

      この運動方程式から，求めたいのは$`\frac{d{\bf \Omega}_{\rm G}}{dt}`$である．これはとても簡単で，次のように求めることができる．

      ```math
      \frac{d{\bf \Omega}_{\rm G}}{dt} = {\rm R}_{g2l}^{-1}{\boldsymbol I}^{-1}{\rm R}_{g2l} {\bf T}_{\rm G}
      ```

      */
      int i = 0;
      std::vector<networkFace *> all_faces;
      for (auto &water : WATERS) {
         auto surfaces = water->getSurfaces();
         all_faces.insert(all_faces.end(), surfaces.begin(), surfaces.end());
      }
      for (const auto &body : rigidbodies)
         if (body->isFloatingBody) {
            //$ ------------------------------ 係留索から受ける力とトルク ----------------------------- */
            //$ フェアリードの節点が隣の線要素から受けている張力ベクトル --> 浮体が受ける力とトルク
            std::array<double, 3> F_mooring = {0., 0., 0.}, T_mooring = {0., 0., 0.};
            //! simulateはアップデートの際に行なっておく．
            for (auto &mooring_line : body->mooringLines) {
               F_mooring += mooring_line->lastPoint->getForce();
               T_mooring += Cross(mooring_line->lastPoint->X - body->COM, mooring_line->lastPoint->getForce());
            }
            //@ ------------------------------ 浮体が流体力とトルク ----------------------------- */
            auto F_ext = _GRAVITY3_ * body->getMass3D();
            auto tmp = calculateFluidInteraction(all_faces, body);
            auto [F_hydro, T_hydro] = tmp.surfaceIntegralOfPressure();
            auto F = F_ext + F_hydro;
            auto T_GLOBAL = T_hydro;
            F += F_mooring;
            T_GLOBAL += T_mooring;
            body->inputJSON.for_each(
                [&](auto key, auto vec_string) {
                   if (key.contains("spring")) {
                      //% ---------------------- 重心の並進移動によって伸びる線形バネによる係留 --------------- */
                      //% simple spring mooring
                      auto X_k = stod(vec_string);
                      std::array<double, 3> init_fairleader_position = {X_k[0], X_k[1], X_k[2]};
                      std::array<double, 3> current_fairleader_position = body->rigidTransformation(init_fairleader_position);
                      std::array<double, 3> anchor = {X_k[3], X_k[4], X_k[5]};
                      double kx, ky, kz;
                      if (X_k.size() >= 9) {
                         kx = X_k[6];
                         ky = X_k[7];
                         kz = X_k[8];
                      } else if (X_k.size() == 7)
                         kz = ky = kx = X_k[6];

                      double diff = Norm(current_fairleader_position - anchor) - Norm(init_fairleader_position - anchor);
                      std::array<double, 3> f = Tddd{kx, ky, kz} * diff * Normalize(anchor - current_fairleader_position);
                      if (X_k.size() > 9) {
                         double default_tention = X_k[10];
                         f += default_tention * Normalize(anchor - current_fairleader_position);
                      }
                      F += f;
                      //  std::cout << "init_fairleader_position = " << init_fairleader_position << std::endl;
                      //  std::cout << "current_fairleader_position = " << current_fairleader_position << std::endl;
                      //  std::cout << "anchor = " << anchor << std::endl;
                      //  std::cout << "diff = " << diff << std::endl;
                      //  std::cout << "f = " << f << std::endl;
                   } else if (key.contains("linear_cable")) {
                      auto X_k = stod(vec_string);
                      std::array<double, 3> init_fairleader_position = {X_k[0], X_k[1], X_k[2]};
                      std::array<double, 3> current_fairleader_position = body->rigidTransformation(init_fairleader_position);
                      std::array<double, 3> anchor = {X_k[3], X_k[4], X_k[5]};
                      double kx, ky, kz;
                      if (X_k.size() >= 9) {
                         kx = X_k[6];
                         ky = X_k[7];
                         kz = X_k[8];
                      } else if (X_k.size() == 7)
                         kz = ky = kx = X_k[6];
                      double diff = Norm(current_fairleader_position - anchor) - Norm(init_fairleader_position - anchor);
                      Tddd direction = Normalize(anchor - current_fairleader_position);
                      std::array<double, 3> f = Tddd{kx, ky, kz} * (anchor - current_fairleader_position);
                      T_GLOBAL += Cross(current_fairleader_position - body->COM, f);
                      F += f;
                   }
                });

            //^ ---------------------- ダンピング --------------------- */

            if (body->inputJSON.find("damping")) {
               const auto c = stod(body->inputJSON.at("damping"));
               std::array<double, 3> c_xyz = {c[0], c[1], c[2]};
               std::array<double, 3> c_abc = {c[3], c[4], c[5]};
               double start_t = 0, end_t = 1E10;
               if (c.size() == 7)
                  start_t = c[6];
               else if (c.size() == 8) {
                  start_t = c[6];
                  end_t = c[7];
               }
               if (start_t <= simulation_time && simulation_time <= end_t) {
                  F -= c_xyz * body->velocityTranslational();
                  T_GLOBAL -= c_abc * body->velocityRotational();
               }
            }

            //^ ------------------------------------------------------ */

            const auto [mx, my, mz, IG, inv_IG] = body->getInertiaGC();
            double a0 = F[0] / mx;
            double a1 = F[1] / my;
            double a2 = F[2] / mz;
            Tddd A_rotaton = Dot(inv_IG, T_GLOBAL);
            auto [a3, a4, a5] = A_rotaton;
            std::ranges::for_each(T6d{a0, a1, a2, a3, a4, a5}, [&](const auto &a_w) { ACCELS[i++] = a_w; });  // 複数浮体がある場合があるので．
         } else
            i += 6;
      // std::cout << Green << "other" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      ACCELS_OUT = ACCELS;
      return ACCELS - ACCELS_IN;
   };

   /* -------------------------------------------------------------------------- */

   //@ --------------------------------------------------- */
   //@        加速度 --> phiphin_t --> 圧力 --> 加速度        */
   //@ --------------------------------------------------- */
   std::vector<double> solveForPhiPhin_t(const std::vector<Network *> &rigidbodies) {
      std::vector<double> convergence;
      for (auto &water : WATERS)
         water->setGeometricProperties();
      auto ACCELS_init = initializeAcceleration(rigidbodies);

      if (ACCELS_init.empty()) {
         setPhiPhin_t(WATERS);
         return convergence;
      }

      int count = 0;
      double alpha = 1.;

      BroydenMethod BM(ACCELS_init, ACCELS_init);

      V_d ACCELS_OUT = ACCELS_init;

      insertAcceleration(rigidbodies, BM.X - BM.dX);
      auto func_ = Func(BM.X - BM.dX, WATERS, rigidbodies, ACCELS_OUT);
      // #define method_Cao1994
      for (auto j = 0; j < 100; ++j) {
#ifdef method_Cao1994
         V_d ACCELS_IN = ACCELS_OUT;
         std::cout << ACCELS_OUT << std::endl;
         auto func = Func(ACCELS_IN, WATERS, rigidbodies, ACCELS_OUT);
         func = ACCELS_IN - ACCELS_OUT;
         std::cout << ACCELS_OUT << std::endl;
         insertAcceleration(rigidbodies, ACCELS_OUT);
         convergence.push_back(Norm(func));

         auto dX = ACCELS_OUT - ACCELS_IN;
         std::cout << "j = " << j << ", alpha = " << alpha << ", Norm(func) = " << Norm(func) << ", " << Red << "Norm(dX) = " << Norm(dX) << colorReset << std::endl;

         if (Norm(dX) < 1E-10 && Norm(func) < 1E-10 && count++ > 4)
            break;
         else
            count = 0;
#else
         auto func = Func(BM.X, WATERS, rigidbodies, ACCELS_OUT);
         convergence.push_back(Norm(func));
         BM.update(func, func_, j == 0 ? 0.1 : alpha);
         func_ = func;
         insertAcceleration(rigidbodies, BM.X);

         std::cout << "j = " << j << ", alpha = " << alpha << ", Norm(func) = " << Norm(func) << ", " << Red << "Norm(BM.dX) = " << Norm(BM.dX) << colorReset << std::endl;

         if (Norm(BM.dX) < 1E-10 && Norm(func) < 1E-10 && count++ > 4)
            break;
         else if (Norm(BM.dX) < 1E-10 && Norm(func) < 1E-10)
            break;
         else
            count = 0;
#endif
      }
      return convergence;
   };
};
