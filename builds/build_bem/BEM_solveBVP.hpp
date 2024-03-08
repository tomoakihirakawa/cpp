#ifndef BEM_solveBVP_H
#define BEM_solveBVP_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

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
// #define use_gmres 100
#define use_lapack

// #define quad_element
// #define linear_element
// #define liear_and_quad_element
std::unordered_map<std::tuple<netP *, netF *>, int> PBF_index;
struct calculateFluidInteraction {
   const Network *PasObj;
   std::vector<networkFace *> actingFaces;
   Tddd force, torque, simplified_drag, simplified_drag_torque;
   double area;
   T6d acceleration;
   std::vector<std::tuple<std::array<networkPoint *, 3>, Tddd, T3Tddd>> PressureVeticies;
   calculateFluidInteraction(const auto &faces /*waterfaces*/, const Network *PasObjIN)
       : PasObj(PasObjIN), force({0., 0., 0.}), torque({0., 0., 0.}), area(0.), PressureVeticies({}), acceleration({0., 0., 0., 0., 0., 0.}) {
      // PasObjと接したfaceの頂点にpressureが設定されている前提
      int count = 0;
      // set PressureVeticies

      for (const auto &f : faces)
         if (f->Neumann) {
            if (std::ranges::all_of(f->getPoints(), [&](const auto &p) { return std::ranges::any_of(
                                                             
                                                             p->getContactFaces(), [&](const auto &F) { return F->getNetwork() == PasObj; }); })) {
               auto [p0, p1, p2] = f->getPoints();
               this->PressureVeticies.push_back({{p0, p1, p2}, {p0->pressure, p1->pressure, p2->pressure}, ToX(f)});
               this->actingFaces.emplace_back(f);
               count++;
            }
         }

      // calculate area
      for (const auto &[p012, P012, X012] : this->PressureVeticies) {
         auto intpX = interpolationTriangleLinear0101(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            area += intpX.J(x0, x1) * w0w1;
      }
      // std::cout << "接触している面の数:" << count << " 表面積:" << area << std::endl;
   };

   // \label{BEM:surfaceIntegralOfPressure}
   std::array<Tddd, 2> surfaceIntegralOfPressure() {
      this->force.fill(0.);
      this->torque.fill(0.);
      for (const auto &[_, P012, X012] : this->PressureVeticies) {
         // auto [pre0, pre1, pre2] = P012;
         // auto [X0, X1, X2] = X012;
         // this->force += (p0 + p1 + p2) / 3. * Cross(X1 - X0, X2 - X0) / 2.;

         auto intpP = interpolationTriangleLinear0101(P012);
         auto intpX = interpolationTriangleLinear0101(X012);
         auto n = TriangleNormal(X012);
         double f;
         Tddd force = {0., 0., 0.};
         Tddd torque = {0., 0., 0.};
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple) {
            f = intpP(x0, x1) * intpX.J(x0, x1) * w0w1;
            force += f;
            torque += f * Cross(intpX(x0, x1) - this->PasObj->COM, n);
         }
         this->force += force * n;
         this->torque += torque;
      }
      return {this->force, this->torque};
   };

   std::array<Tddd, 2> surfaceIntegralOfVerySimplifiedDrag() {
      this->simplified_drag.fill(0.);
      this->simplified_drag_torque.fill(0.);
      for (const auto &[p012, _, X012] : this->PressureVeticies) {
         auto [p0, p1, p2] = p012;
         auto [X0, X1, X2] = X012;
         const Tddd relative_U0 = p0->U_BEM - PasObj->velocityRigidBody(X0);
         const Tddd relative_U1 = p1->U_BEM - PasObj->velocityRigidBody(X1);
         const Tddd relative_U2 = p2->U_BEM - PasObj->velocityRigidBody(X2);
         // this->force += (p0 + p1 + p2) / 3. * Cross(X1 - X0, X2 - X0) / 2.;
         auto intpRelativeVelocity = interpolationTriangleLinear0101(T3Tddd{relative_U0, relative_U1, relative_U2});
         auto intpX = interpolationTriangleLinear0101(X012);
         const double nu = 1000 * 1000 * 1.004 * 10E-6;  // m2 /s
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

void setPhiPhin(Network &water) {
   /* -------------------------------------------------------------------------- */
   /*                         phinOnFace, phintOnFaceの設定                       */
   /* -------------------------------------------------------------------------- */
   // b! 点
   std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << colorReset << std::endl;

#pragma omp parallel
   for (const auto &p : water.getPoints())
#pragma omp single nowait
   {
      auto setNeumann = [&p](networkFace *const &f) {
         if (isNeumannID_BEM(p, f)) {
            if (f == nullptr) {
               p->phiOnFace.insert({nullptr, 1E+30});
               p->phitOnFace.insert({nullptr, 1E+30});
               p->phinOnFace.insert({nullptr, std::get<1>(p->phiphin) = Dot(uNeumann(p), p->getNormalNeumann_BEM())});
               p->phintOnFace.insert({nullptr, 1E+30});
            } else {
               p->phiOnFace.insert({f, 1E+30});
               p->phitOnFace.insert({f, 1E+30});
               p->phinOnFace.insert({f, Dot(uNeumann(p, f), f->normal)});
               p->phintOnFace.insert({f, 1E+30});
            }
         }
      };

      auto setDirichlet = [&p](networkFace *const &f) {
         if (isDirichletID_BEM(p, f)) {
            p->phiOnFace.insert({nullptr, 1E+30});
            p->phitOnFace.insert({nullptr, 1E+30});
            p->phinOnFace.insert({nullptr, 1E+30});
            p->phintOnFace.insert({nullptr, 1E+30});
         }
      };

      p->phiOnFace.clear();
      p->phitOnFace.clear();
      p->phinOnFace.clear();
      p->phintOnFace.clear();

      // 全ての組{p,f}を調べ，使えるものをチェックする
      setNeumann(nullptr);
      for (const auto &f : p->getFaces())
         setNeumann(f);

      setDirichlet(nullptr);
      for (const auto &f : p->getFaces())
         setDirichlet(f);
   }

   // b! 面
   std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << colorReset << std::endl;
#pragma omp parallel
   for (const auto &f : water.getFaces())
#pragma omp single nowait
   {
      auto [p0, p1, p2] = f->getPoints();
      std::get<0>(f->phiphin) = (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3.;
   }
};

// b@ -------------------------------------------------------------------------- */
// b@                                   BEM_BVP                                  */
// b@ -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

### BIEの離散化

BIEをGauss-Legendre積分で離散化すると，

```math
\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} {\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left( {\sum\limits_{j=0}^2 {{{\left( {{\phi_n}} \right)}_{k_\vartriangle,j }}{N_{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}_{k _\vartriangle}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}\left\|\frac{{\partial{{\bf x}_{k _\vartriangle}}}}{{\partial{\xi_0}}} \times \frac{{\partial{\bf{x}}_{k _\vartriangle}}}{{\partial{\xi_1}}}\right\|} \right)} }=
```

```math
\alpha_{i_\circ}(\phi)_{i_\circ}-\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} \sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{\bf x}_{k \vartriangle}}(\pmb{\xi})-{{\bf x}_{i_\circ} }}{{{{\| {{{\bf x}_{k_\vartriangle}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}} \cdot\left(\frac{{\partial {\bf{x}}_{k_\vartriangle}}}{{\partial {\xi_0}}}\times\frac{{\partial {\bf{x}}_{k_\vartriangle}}}{{\partial {\xi_1}}}\right)}\right)}
```

離散化では，$`\phi_{i_\circ}`$と$`{\phi_n}_{i_\circ}`$の係数を知りたいので，
$`\phi_{k_\vartriangle}(\pmb{\xi})`$と$`{\phi_n}_{k_\vartriangle}(\pmb{\xi})`$と書くのではなく，
$`\phi_{i_\circ}`$と$`{\phi_n}_{i_\circ}`$が見えるように$`\phi_{k_\vartriangle}(\pmb{\xi})`$と$`{\phi_n}_{k_\vartriangle}(\pmb{\xi})`$の補間を書いている．

ここで，$`\phi_{k_\vartriangle,j}`$における$`k_\vartriangle`$は三角形要素の番号，$`j`$は三角形要素の頂点番号．
$`N_j`$は三角形要素の形状関数，$`\pmb{\xi}`$は三角形要素の内部座標，$`w_0,w_1`$はGauss-Legendre積分の重み，$`\alpha_{i_\circ}`$は原点$`i_\circ`$における立体角，$`\phi`$はポテンシャル，$`\phi_n`$は法線方向のポテンシャル，$`\bf{x}`$は空間座標，$`{\bf x}_{i_\circ}`$は原点の空間座標である．

* $`\phi_{k_\vartriangle}`$は補間で作った関数
* $`\phi_{k_\vartriangle,j}`$は補間を構成する節点$`j`$での値
* $`\phi_{i_\circ}`$はより直接的にある節点$`i_\circ`$での値

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
\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{\bf x}_{k _\vartriangle}}(\pmb{\xi})-{{\bf x}_{i_\circ} }}{{{{\| {{{\bf x}_{k_\vartriangle}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}}} (1-\xi_0)\right)}
```

NOTE: ちなみに，$`\frac{1-\xi_0}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}`$の分子に$`1-\xi_0`$があることで，
関数の特異的な変化を抑えることができる．プログラム上ではこの性質が利用できるように，この分数をまとめて計算している．

*/

struct BEM_BVP {
   const bool Neumann = false;
   const bool Dirichlet = true;
   lapack_lu *lu = nullptr;
#if defined(use_lapack)
#else
   ludcmp_parallel *lu;
#endif

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
   V_d knowns;
   std::vector<std::vector<std::array<double, 2>>> IGIGn;
   BEM_BVP(){};
   ~BEM_BVP() {
      if (this->lu) delete this->lu;
   };

   int pf2Index(const networkPoint *p, const networkFace *f) const {
      // if(PBF_index.find(pf2ID(p, f)) == PBF_index.end())
      //    throw std::runtime_error("PBF_indexに存在しない{p,f}が指定されました．");
      return PBF_index.at(pf2ID(p, f));
   };
   /**
   isNeumannID_BEMとisDirichletID_BEMの両方を満たす{p,f}は存在しない．
   */

   void setIGIGn(const std::vector<Network *> &WATERS) {

      for (auto &net : WATERS)
         net->setGeometricProperties();

      IGIGn = std::vector<std::vector<std::array<double, 2>>>(PBF_index.size(), std::vector<std::array<double, 2>>(PBF_index.size(), {0., 0.}));

#define use_rigid_mode
      Timer timer;
      std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;
      std::cout << Magenta << timer() << colorReset << std::endl;
      /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

      #### 係数行列の作成

      実際のプログラムでは，$`{\bf A}{\bf x}={\bf b}`$の形で整理することが多い．
      上のようにBIEは離散化されるが，
      この式を見ても，係数行列$`\bf A`$とベクトル$`\bf b`$を具体的にどのように作成するかわかりにくいかもしれない．

      - $`\phi`$の係数行列を$`\mathbf{M}`$
      - $`\phi_n`$の係数行列を$`\mathbf{N}`$
      - $`\mathbf{\Phi}`$を$`\phi`$のベクトル
      - $`\mathbf{\Phi_n}`$を$`\phi_n`$のベクトル

      として，次のような連立一次方程式を得る．

      ```math
      \mathbf{N} \mathbf{\Phi_n} = \mathbf{M} \mathbf{\Phi}
      ```

      $`{\bf A}{\bf x}={\bf b}`$の形にして，未知変数$`{\bf x}`$を求めるわけだが，
      未知変数が$`\phi`$か$`\phi_n`$かは，境界条件によって決まるので，
      境界条件に応じて，$`{\bf A},{\bf b}`$を間違えずに作成する必要がある．

      ここでは，$`A`$を`IGIGn`，$`b`$を`knowns`としている．

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
      constexpr std::array<double, 2> ZEROS2 = {0., 0.};
      constexpr std::array<double, 3> ZEROS3 = {0., 0., 0.};
      //@ Q：quadratic elementの不安定は，積分点を奇数にするか偶数にするかに依存するか？
      std::vector<std::tuple<double, double, double, Tddd>> t0_t1_ww_N012_LOWRESOLUTION;
      std::vector<std::tuple<double, double, double, Tddd>> t0_t1_ww_N012_HIGHRESOLUTION;
      for (int i = 0; const auto &[t0, t1, ww] : __array_GW6xGW6__)
         t0_t1_ww_N012_LOWRESOLUTION.push_back({t0, t1, ww, ModTriShape<3>(t0, t1)});
      // auto t0_t1_ww_N012_HIGHRESOLUTION = t0_t1_ww_N012_LOWRESOLUTION;
// 
      for (int i = 0;const auto &[t0, t1, ww] : __array_GW13xGW13__)
         t0_t1_ww_N012_HIGHRESOLUTION.push_back({t0,t1,ww,ModTriShape<3>(t0, t1)});

#define USE_LINEAR_ELEMENT
      // #define USE_QUADRATIC_ELEMENT

#if defined(USE_LINEAR_ELEMENT)
      std::cout << "線形要素を使ってBIEを離散化" << std::endl;
#elif defined(USE_QUADRATIC_ELEMENT)
      std::cout << "二次要素を使ってBIEを離散化" << std::endl;
#endif

#pragma omp parallel
      for (const auto &[PBF, index] : PBF_index)
#pragma omp single nowait
      {
         auto [origin, _] = PBF;
         double origin_ign_rigid_mode = 0.;

         double dist_base = 0;  
         for(auto f:origin->getFaces())
               dist_base += f->area;
         dist_base = std::sqrt(dist_base);

         auto &IGIGn_Row = IGIGn[index];
         double nr;
         std::array<double, 3> cross, R;
         std::array<double, 6> N012_potential;
#if defined(USE_LINEAR_ELEMENT)
         std::array<std::tuple<networkPoint *, networkFace *, double, Tddd, double>, 3> ret;
#elif defined(USE_QUADRATIC_ELEMENT)
         std::array<std::tuple<networkPoint *, networkFace *, double, Tddd, double>, 6> ret;
#endif

         // for(const auto &water : WATERS)
         //    for (const auto &integ_f : water->getFacesVector())
         //    {
         //       auto [l0, l1, l2] = integ_f->getLines();
         //       auto fs = l0->getFaces();
         //       if(fs.size() == 2)
         //       {
         //          auto [p0, p1, p2] = integ_f->getPoints(origin);
         //          QuadPoints(l0, fs[0]);
         //          QuadPoints(l1, fs[1]);
         //       }
         //    }

         for (const auto &water : WATERS) {
            double ig;
            Tddd ign, value_for_rigid_mode = {0., 0., 0.};
            for (const auto &integ_f : water->getFacesVector()) {
               value_for_rigid_mode.fill(0.);
               auto [p0, p1, p2] = integ_f->getPoints(origin);
               if(p0!=origin)
               {
                  auto q0 = p0, q1 = p1, q2 = p2;                  
                  if(Norm(p0->X - origin->X) >= Norm(p1->X - origin->X) && Norm(p2->X - origin->X) >= Norm(p1->X - origin->X)){
                     p0 = q1;
                     p1 = q2;
                     p2 = q0;
                  } else if(Norm(p0->X - origin->X) >= Norm(p2->X - origin->X) && Norm(p1->X - origin->X) >= Norm(p2->X - origin->X)){
                     p0 = q2;
                     p1 = q0;
                     p2 = q1;
                  }
               }

               cross = Cross(p1->X - p0->X, p2->X - p0->X);
               const std::array<std::array<double, 3>, 3> X012 = {p0->X, p1->X, p2->X};
               //! 遠い近いの判定基準がないので，とりあえず基準はorginのfaceのsqrt(面積)としておく
               const double dist = 5 * Norm((p0->X+ p1->X+ p2->X)/3. - origin->X);               
#if defined(USE_LINEAR_ELEMENT)
               ret = {{{p0, integ_f, 0., ZEROS3, 0.},
                       {p1, integ_f, 0., ZEROS3, 0.},
                       {p2, integ_f, 0., ZEROS3, 0.}}};
#elif defined(USE_QUADRATIC_ELEMENT)
               constexpr std::array<std::array<double, 2>, 3> M_sub{{{0., 0.5}, {0.5, 0.}, {0.5, 0.5}}};
               auto quadpoint = QuadPoints(origin, integ_f);
               ret = {{{quadpoint.points[0] /*0*/, quadpoint.faces[0], 0., ZEROS3, 0.},
                       {quadpoint.points[1] /*1*/, quadpoint.faces[1], 0., ZEROS3, 0.},
                       {quadpoint.points[2] /*2*/, quadpoint.faces[2], 0., ZEROS3, 0.},
                       {quadpoint.points[3] /*3*/, integ_f, 0., ZEROS3, 0.},
                       {quadpoint.points[4] /*4*/, integ_f, 0., ZEROS3, 0.},
                       {quadpoint.points[5] /*5*/, integ_f, 0., ZEROS3, 0.}}};
#endif

               for (const auto &[t0, t1, ww, N012_geometry] : (dist < dist_base ? t0_t1_ww_N012_HIGHRESOLUTION : t0_t1_ww_N012_LOWRESOLUTION) ) {
                  ig = ww * ((1. - t0)) / (nr = Norm(R = (Dot(N012_geometry, X012) - origin->X)));
                  ign = ig * R / (nr * nr);
#if defined(USE_LINEAR_ELEMENT)
                  for (auto i = 0; i < 3; i++) {
                     FusedMultiplyIncrement(ig, N012_geometry[i], std::get<2>(ret[i]));
                     FusedMultiplyIncrement(-ign, N012_geometry[i], std::get<3>(ret[i]));
                  }
                  value_for_rigid_mode += ign * Total(N012_geometry);
#elif defined(USE_QUADRATIC_ELEMENT)
                  N012_potential = TriShape<6>(Dot(N012_geometry, M_sub), quadpoint.ignore);
                  for (auto i = 0; i < 6; i++) {
                     FusedMultiplyIncrement(ig, N012_potential[i], std::get<2>(ret[i]));
                     FusedMultiplyIncrement(-ign, N012_potential[i], std::get<3>(ret[i]));
                  }
                  value_for_rigid_mode += ign * Total(N012_potential);
#endif
               }

               //! 2Aまたは2Anをかける
               auto Atimes2 = Norm(cross);
               for (auto &[_, __, ig, ign_tmp, ign_result] : ret) {
                  ig *= Atimes2;
                  ign_result = (origin == p0 ? 0. : Dot(ign_tmp, cross));
               }

               //! この面 k triangle
               if (p0 != origin)
                  origin_ign_rigid_mode += Dot(value_for_rigid_mode, cross);

               for (const auto &[p, integ_f, ig, __, ign_result] : ret)
                  IGIGn_Row[pf2Index(p, integ_f)] += Tdd{ig, ign_result};  // この面に関する積分において，φまたはφnの寄与
            }
         }
         /* -------------------------------------------------------------------------- */
         /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

         ### リジッドモードテクニック

         立体角$`\alpha`$はBIEの係数行列の対角成分で正確に計算したいが，正確に計算するプログラムを書くのは面倒である．
         しかし，素直に幾何学的な観点から立体角を計算するのではなく，BIEの式を使って積分で計算することができる．

         立体角はポテンシャルとは関係なく，境界面の形状だけに依存するので，適当に$`\phi=1`$としても立体角は変わらない．
         こうすると，$`\alpha({\bf a}) = -\int\int{\nabla G({\bf x},{\bf a})\cdot{\bf n}({\bf x})dS}`$となり，
         この式は，境界面形状を線形要素を使って離散化すると，次のようになる．

         ```math
         \alpha_{i_\circ}=\sum\limits_{k_\vartriangle}
         \left(
         {2A_{k_\vartriangle}{\bf n}_{k_\vartriangle}}\cdot
         \sum\limits_{{\xi_1},{w_1}}\sum\limits_{{\xi_0},{w_0}}
         {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{\bf x}_{k _\vartriangle}}(\pmb{\xi})-{{\bf x}_{i_\circ} }}{{{{\| {{{\bf x}_{k_\vartriangle}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}}} (1-\xi_0)\right)}
         \right)
         ```

         この式の右辺の一部，
         原点$`i_\circ`$を頂点とする三角形$`k_{\vartriangle}`$に対する和だけを左辺に移項したものは，
         係数行列の対角成分と同じになっている．
         これはリジッドモードテクニックと呼ばれていて，
         分子が小さくなる特異的な計算を省き，立体角の計算もまとめて対角成分を計算することができる方法である．

         ただし，線形要素の場合，原点$`i_\circ`$を頂点とする三角形$`k_{\vartriangle}`$に対する計算，$`{\bf n}_{k_\vartriangle}\cdot ({{\bf x}_{k_\vartriangle}}(\pmb{\xi})-{{\bf x}_{i_\circ}})=0`$となるため，和をとる必要はない．
         よって，そもそも線形要素の場合は，特異的な計算は含まれない．

         */

#if defined(use_rigid_mode)
         std::get<1>(IGIGn_Row[index]) = origin_ign_rigid_mode;
#else
         std::get<1>(IGIGn_Row[index]) += origin->getSolidAngle();
#endif
      }
      std::cout << Green << "離散化にかかった時間" << timer() << colorReset << std::endl;
   };

   void makeMatrix() {
      double max_value = 0;
      for (auto i = 0; i < IGIGn.size(); ++i) {
         if (max_value < std::abs(std::get<0>(IGIGn[i][i])))
            max_value = std::abs(std::get<0>(IGIGn[i][i]));
      }

      mat_kn = mat_ukn = VV_d(PBF_index.size(), V_d(PBF_index.size(), 0.));
#pragma omp parallel
      for (const auto &[i_row, i] : PBF_index)
#pragma omp single nowait
      {
         std::array<double, 2> igign;
         for (const auto &[j_col, j] : PBF_index) {
            igign = IGIGn[i][j];
            // 未知変数の係数行列は左，既知変数の係数行列は右
            if (isNeumannID_BEM(j_col))
               igign = {-std::get<1>(igign), -std::get<0>(igign)};

            /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

            係数行列`IGIGn`は，左辺の$`I_G \phi_n`$，右辺の$`I_{G_n}\phi`$の係数．

            ```math
            (I_G)_{i_\circ,j_\circ} (\phi_n)_{j_\circ} = (I_{Gn})_{i_\circ,j_\circ}  \phi_{j_\circ}
            ```

            境界条件に応じて，未知変数は$`\phi,\phi_n`$のどちらかに決まる．
            未知変数が$`\phi`$の場合（Dirichlet境界条件の場合），
            係数行列`IGIGn`中で対応する列を符号変えて入れ替えることで移項したことになる．


            移項前:
            ```math
            \begin{bmatrix}I_{G0} & I_{G1} & I_{G2} & I_{G3}\end{bmatrix} \begin{bmatrix}\phi _{n0} \\ \phi _{n1} \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}I_{Gn0} & I_{Gn1} & I_{Gn2} & I_{Gn3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _1 \\ \phi _2 \\ \phi _3\end{bmatrix}
            ```

            移項後:
            ```math
            \begin{bmatrix}I_{G0} & -I_{Gn1} & I_{G2} & I_{G3}\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}I_{Gn0} & -I_{G1} & I_{Gn2} & I_{Gn3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}
            ```

            多重節点(1と3が多重節点の場合):
            ```math
            \begin{bmatrix}0 & 1 & 0 & 0\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}0 & 0 & 0 & 1\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}
            ```

            */
            mat_ukn[i][j] = std::get<0>(igign);
            mat_kn[i][j] = std::get<1>(igign);
         }
         auto [a, _] = i_row;
         if (a->CORNER && isNeumannID_BEM(i_row) /*行の変更*/) {
            std::ranges::fill(mat_ukn[i], 0.);
            std::ranges::fill(mat_kn[i], 0.);
            mat_ukn[i][i] = max_value;                    // φの系数
            mat_kn[i][pf2Index(a, nullptr)] = max_value;  // φの系数移行したからマイナス？　いいえ，移項を考慮した上でこれでいい．
         }
      }
   };

   template <typename T1, typename T2, typename T3>
   void storePhiPhinCommon(const Network &water, const V_d &ans, T1 phiphinProperty, T2 phiOnFaceProperty, T3 phinOnFaceProperty) const {
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, f] = PBF;
         if (isDirichletID_BEM(PBF)) {
            (p->*phinOnFaceProperty).at(f) = std::get<1>(p->*phiphinProperty) = ans[i];
            (p->*phiOnFaceProperty).at(f) = std::get<0>(p->*phiphinProperty);
         }
         if (isNeumannID_BEM(PBF))
            (p->*phiOnFaceProperty).at(f) = std::get<0>(p->*phiphinProperty) = ans[i];
      }

      for (const auto &p : water.getPoints())
         if (p->Neumann) {
            double total = 0;
            std::get<0>(p->*phiphinProperty) = 0;
            std::ranges::for_each(p->getFaces(), [&total](const auto &f) { total += f->area; });
            for (const auto &f : p->getFaces()) {
               if ((p->*phiOnFaceProperty).count(f))
                  std::get<0>(p->*phiphinProperty) += (p->*phiOnFaceProperty).at(f) * f->area / total;
               else
                  std::get<0>(p->*phiphinProperty) += (p->*phiOnFaceProperty).at(nullptr) * f->area / total;
            }
         }
   }
   void storePhiPhin(const Network &water, const V_d &ans) const {
      storePhiPhinCommon(water, ans, &networkPoint::phiphin, &networkPoint::phiOnFace, &networkPoint::phinOnFace);
   }

   void storePhiPhin_t(const Network &water, const V_d &ans) const {
      storePhiPhinCommon(water, ans, &networkPoint::phiphin_t, &networkPoint::phitOnFace, &networkPoint::phintOnFace);
   }

   void isSolutionFinite(const auto &water) const {
      for (const auto &p : water.getPoints()) {
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

   void solve(std::vector<Network *> WATERS, const Buckets<networkPoint *> &FMM_BucketsPoints, const Buckets<networkFace *> &FMM_BucketsFaces) {

      for (auto water : WATERS)
         setPhiPhin(*water);

      // std::vector<networkPoint *> points;

      // for (auto water : WATERS)
      //    for (auto &p : water->getPoints())
      //       points.emplace_back(p);

      PBF_index.clear();
      // PBF_index.reserve(3 * points.size());
      int i = 0;

      for (const auto water : WATERS)
         for (const auto &q : water->getPoints())
            for (const auto &id : variableIDs(q))
               if (PBF_index.find(id) == PBF_index.end())
                  PBF_index[id] = i++;

      std::cout << Red << "   unknown size : " << PBF_index.size() << colorReset << std::endl;
      std::cout << Red << "water node size : " << i << colorReset << std::endl;

      setIGIGn(WATERS);

      std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;

      makeMatrix();

      knowns.resize(PBF_index.size());
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, f] = PBF;
         if (isDirichletID_BEM(PBF) && isNeumannID_BEM(PBF))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "isDirichletID_BEM(P,F) && isNeumannID_BEM(P,F)");
         else if (isDirichletID_BEM(PBF)) {
            knowns[i] = p->phi_Dirichlet = std::get<0>(p->phiphin);
         } else if (isNeumannID_BEM(PBF))
            knowns[i] = p->phinOnFace.at(f);  // はいってない？はいってた．
         else
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "cannot be");
      }

      TimeWatch watch;

      ans.resize(knowns.size(), 0.);

      if (this->lu)
         delete this->lu;
#if defined(use_CG)
      this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/);
      std::cout << "The conjugate gradient is used" << std::endl;
      GradientMethod gd1(mat_ukn);
      ans = gd1.solve(ParallelDot(mat_kn, knowns), {}, 1E-1);
      GradientMethod gd2(mat_ukn);
      ans = gd2.solveCG(ParallelDot(mat_kn, knowns), ans);
#elif defined(use_gmres)
      std::cout << "use gmres for BIE" << std::endl;
      auto b = ParallelDot(mat_kn, knowns);
      gmres gm(mat_ukn, b, ans, use_gmres);
      std::cout << Red << "gmres error = " << gm.err << colorReset << std::endl;
      if (!isFinite(gm.err, 1E+5) || !isFinite(gm.x, 1E+5)) {
         std::cout << Red << "gmres failed" << std::endl;
         std::cout << Red << "lapack lu decomposition" << std::endl;
         this->lu = new lapack_lu(mat_ukn, b, ans);
      } else if (simulation_time < 0.001) {
         std::cout << Red << "new lapack lu decomposition" << std::endl;
         this->lu = new lapack_lu(mat_ukn, b, ans);
      } else {
         ans = gm.x;
         for (auto i = 1; i < 5; i++) {
            if (gm.err < 1E-8)
               break;
            gm.Restart(mat_ukn, b, ans, use_gmres);  //\label{SPH:gmres}
            ans = gm.x;
            std::cout << "       gm.err : " << gm.err << std::endl;
            std::cout << " actual error : " << Norm(b - ParallelDot(mat_ukn, ans)) << std::endl;
         }
      }
         // auto err = Norm(Dot(mat_ukn, gm.x) - Dot(mat_kn, knowns));
         // std::cout << err << std::endl;
#elif defined(use_lapack)
      std::cout << "lapack lu decomposition" << std::endl;
      //! A.x = b
      this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/, ans /*解*/, ParallelDot(mat_kn, knowns) /*既知のベクトル（右辺）*/);
#endif

      std::cout << colorReset << "update p->phiphin and p->phinOnFace for Dirichlet boundary" << colorReset << std::endl;

      for (auto water : WATERS)
         storePhiPhin(*water, ans);

      std::cout << Green << "Elapsed time for solving BIE: " << Red << watch() << colorReset << " s\n";

      for (auto water : WATERS)
         isSolutionFinite(*water);
   };

   // b! ------------------------------------------------------------------------------ */
   // b!                             solve phi_t and phi_n_t                            */
   // b! ------------------------------------------------------------------------------ */

   /*DOC_EXTRACT 0_4_0_1_FLOATING_BODY_SIMULATION

   ## 浮体動揺解析

   BEM-MELで浮体動揺解析ができるようにするのは簡単ではない．
   浮体に掛かる圧力の計算に必要な$`\phi_t`$が簡単には求まらないためである．
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

   Dalena and Tanizawa's methodは，$`\phi`$に関するBVPと全く違うBVPを解く必要があり，
   係数行列も違うので，新たに行列を構成する必要がある．
   この行列に関してはあまり研究されていないので，あまり使われない理由だと考えられる．

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

   // \label{BEM:setPhiPhin_t}
   void setPhiPhin_t() const {
#pragma omp parallel
      for (const auto &[PBF, i] : PBF_index)
#pragma omp single nowait
      {
         auto [p, F] = PBF;
         //!!ノイマンの場合はこれでDphiDtは計算できませんよ
         if (isDirichletID_BEM(PBF))
            p->phitOnFace.at(F) = std::get<0>(p->phiphin_t) = p->aphiat(0.);
         else if (isNeumannID_BEM(PBF)) {
            for (auto &[f, phin_t] : p->phintOnFace) {
               // phin_t = std::get<1>(p->phiphin_t) = (f != nullptr) ? phint_Neumann(f) : phint_Neumann(p);  // \label{BEM:setphint}

               phin_t = std::get<1>(p->phiphin_t) = (f != nullptr) ? phint_Neumann(p, f) : phint_Neumann(p);  // \label{BEM:setphint}
            }
         }
      }
   };

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
         } else {
            // if net is not floating, then acceleration is not updated.
            std::ranges::for_each(net->acceleration, [&](auto &a_w) { i++;/*a_w = BM_X[i++];*/ });
         }
      }
   }

   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   V_d Func(const auto &ACCELS_IN, const std::vector<Network *> WATERS, const std::vector<Network *> &rigidbodies) {
      TimeWatch watch;
      auto ACCELS = ACCELS_IN;
      //* --------------------------------------------------- */
      //*                  加速度 --> phiphin_t                */
      //* --------------------------------------------------- */
      setPhiPhin_t();
      knowns.resize(PBF_index.size());
#pragma omp parallel
      for (const auto &[PBF, i] : PBF_index)
#pragma omp single nowait
      {
         auto [p, f] = PBF;
         if (isDirichletID_BEM(PBF))
            knowns[i] = p->phitOnFace.at(f);
         else if (isNeumannID_BEM(PBF))
            knowns[i] = p->phintOnFace.at(f);
      }
      // std::cout << Green << "set knowns" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      ans.resize(knowns.size());
      this->lu->solve(ParallelDot(mat_kn, knowns) /*既知のベクトル（右辺）*/, ans /*解*/);
      // std::cout << Green << "solve by LU" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      //@ -------------------------------------------------------------------------- */
      //@                    update p->phiphin_t and p->phinOnFace                   */
      //@ -------------------------------------------------------------------------- */

      for (auto &water : WATERS)
         storePhiPhin_t(*water, ans);
      // std::cout << Green << "storePhiPhin_t" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

      //* --------------------------------------------------- */
      //*                 phiphin_t --> 圧力                   */
      //* --------------------------------------------------- */

      for (const auto &[PBF, i] : PBF_index) {
         auto [p, f] = PBF;
         if (isDirichletID_BEM(PBF))
            p->pressure = p->pressure_BEM = 0;
         else
            p->pressure = p->pressure_BEM = -_WATER_DENSITY_ * (std::get<0>(p->phiphin_t) + 0.5 * Dot(p->U_BEM, p->U_BEM) + _GRAVITY_ * p->height());
      }

      //* --------------------------------------------------- */
      //*              圧力 ---> 力 --> 加速度                  */
      //* --------------------------------------------------- */
      /*DOC_EXTRACT 0_4_0_2_FLOATING_BODY_SIMULATION

      $`{\bf I}`$は慣性モーメントテンソル（2階のテンソル）．
      実際の実験では，浮体のある基本的な姿勢における主慣性モーメント$`{\bf I}_{\rm principal}={(I_x, I_y, I_z)}`$が与えられる．
      主慣性モーメント$`{\bf I}_{\rm principal}`$から，global座標における浮体の慣性モーメントテンソルを求めるには，次のように考えればいい．

      ```math
      \begin{aligned}
      {\bf I_{\rm principal mat}}\frac{d{\bf \Omega}_{\rm L}}{dt} &= {\bf T}_{\rm L}\\
      {\bf I_{\rm principal mat}} {\rm R}_{g2l} \frac{d{\bf \Omega}_{\rm G}}{dt}& = {\rm R}_{g2l}{\bf T}_{\rm G}\\
      {\rm R}_{g2l}^{-1}{\bf I_{\rm principal mat}} {\rm R}_{g2l} \frac{d{\bf \Omega}_{\rm G}}{dt}& = {\bf T}_{\rm G}\\
      \end{aligned}
      ```
      
      */
      int i = 0;
      std::vector<networkFace *> all_faces;
      for (auto &water : WATERS)
         all_faces.insert(all_faces.end(), water->getFaces().begin(), water->getFaces().end());
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
            //% ---------------------- 重心の並進移動によって伸びる線形バネによる係留 --------------- */
            //% simple spring mooring
            if (body->inputJSON.find("spring")) {
               const auto X_k = stod(body->inputJSON.at("spring"));
               std::array<double, 3> origin = {X_k[0], X_k[1], X_k[2]};
               std::array<double, 3> k = {X_k[3], X_k[4], X_k[5]};
               F += k * (origin - body->COM);
            }
            //^ ---------------------- ダンピング --------------- */
            if (body->inputJSON.find("damping")) {
               const auto c = stod(body->inputJSON.at("damping"));
               std::array<double, 3> c_xyz = {c[0], c[1], c[2]};
               std::array<double, 3> c_abc = {c[3], c[4], c[5]};
               double start_t = 0, end_t = 1E10;
               if(c.size() == 7)
                  start_t = c[6];
               else if(c.size() == 8){
                  start_t = c[6];
                  end_t = c[7];
               }
               if( start_t <= simulation_time && simulation_time <= end_t){
                  F -= c_xyz * body->velocityTranslational();
                  T_GLOBAL -= c_abc * body->velocityRotational();
               }
            }
            //% -------------------------------------------------------------------------- */
            const auto [mx, my, mz, I_accounting_float_attitude] = body->getInertiaGC();
            // auto [Drag_F_hydro, Drag_T_hydro] = tmp.surfaceIntegralOfVerySimplifiedDrag();
            // F += Drag_F_hydro;
            // T += Drag_T_hydro;

            auto [a0, a1, a2] = F / Tddd{mx, my, mz};
            // auto R = body->quaternion.Rv();
            // auto RT = Transpose(R);
            // T3Tddd I_accounting_float_attitude = Dot(R,Dot(T3Tddd{{{Ix, 0., 0.},{0., Iy, 0.},{0., 0., Iz}}}, RT));            
            // I_accounting_float_attitude = I;
            Tddd A_rotaton;
            Solve(I_accounting_float_attitude, A_rotaton, T_GLOBAL);
            auto [a3, a4, a5] = A_rotaton;
            
            std::ranges::for_each(T6d{a0, a1, a2, a3, a4, a5}, [&](const auto &a_w) { ACCELS[i++] = a_w; });  // 複数浮体がある場合があるので．
            // write out details of the body
            // std::cout << Green << "mass = " << body->mass << std::endl;
            // std::cout << Green << "inertia = " << body->getInertiaGC() << std::endl;
         } else
            i += 6;
      // std::cout << Green << "other" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      return ACCELS - ACCELS_IN;
   };

   /* -------------------------------------------------------------------------- */

   //@ --------------------------------------------------- */
   //@        加速度 --> phiphin_t --> 圧力 --> 加速度        */
   //@ --------------------------------------------------- */
   std::vector<double> solveForPhiPhin_t(std::vector<Network *> WATERS, const std::vector<Network *> &rigidbodies) {
      std::vector<double> convergence;
      for (auto &water : WATERS)
         water->setGeometricProperties();
      auto ACCELS_init = initializeAcceleration(rigidbodies);

      if (ACCELS_init.empty()) {
         setPhiPhin_t();
         return convergence;
      }

      int count = 0;
      double alpha = 1.;

      BroydenMethod BM(ACCELS_init, ACCELS_init);

      insertAcceleration(rigidbodies, BM.X - BM.dX);
      auto func_ = Func(BM.X - BM.dX, WATERS, rigidbodies);

      for (auto j = 0; j < 100; ++j) {

         auto func = Func(BM.X, WATERS, rigidbodies);
         convergence.push_back(Norm(func));
         BM.update(func, func_, j == 0 ? 0.1 : alpha);
         func_ = func;
         insertAcceleration(rigidbodies, BM.X);         

         std::cout << "j = " << j << ", alpha = " << alpha << ", Norm(func) = " << Norm(func) << ", " << Red << "Norm(BM.dX) = " << Norm(BM.dX) << colorReset << std::endl;

         if (Norm(BM.dX) < 1E-9 && Norm(func) < 1E-9 && count++ > 4) 
            break;
         else if (Norm(BM.dX) < 1E-12 && Norm(func) < 1E-12)
            break;
         else
            count = 0;

      }
      return convergence;
   };
};

#endif