#ifndef BEM_solveBVP_H
#define BEM_solveBVP_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

/*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

## 境界値問題

### 基礎方程式

```math
\begin{align}
\nabla\cdot\nabla \phi& = 0&&\text{in}&&{\bf x} \in \Omega(t),\\
\frac{\partial\phi}{\partial t} +\frac{1}{2}\nabla\phi\cdot\nabla\phi - g z &=0 &&\text{on}&&{\bf x} \in \Gamma^{(\rm D)}(t),\\
\phi_n + {{\bf u}_b}\cdot{{\bf n}_b} &=0&&\text{on}&&{\bf x}\in \Gamma^{(\rm N)}(t),
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
$`G=1/\|{\bf x}-{\bf a}\|`$がラプラス法廷式の基本解であり，$`\phi`$は境界におけるポテンシャルの分布である．

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
struct calculateFroudeKrylovForce {
   std::vector<networkFace *> actingFaces;
   Tddd force, torque;
   double area;
   T6d acceleration;
   std::vector<std::tuple<Tddd, T3Tddd>> PressureVeticies;

   calculateFroudeKrylovForce(const std::unordered_set<networkFace *> faces /*waterfaces*/, const Network *PasObj)
       : force({0., 0., 0.}), torque({0., 0., 0.}), area(0.), PressureVeticies({}), acceleration({0., 0., 0., 0., 0., 0.}) {
      // PasObjと接したfaceの頂点にpressureが設定されている前提
      int count = 0;
      // set PressureVeticies
      for (const auto &f : faces)
         if (f->Neumann) {
            if (std::ranges::all_of(f->getPoints(), [&](const auto &p) { return std::ranges::any_of(p->getContactFaces(), [&](const auto &F) { return F->getNetwork() == PasObj; }); })) {
               auto [p0, p1, p2] = f->getPoints();
               this->PressureVeticies.push_back({{p0->pressure, p1->pressure, p2->pressure}, ToX(f)});
               this->actingFaces.emplace_back(f);
               count++;
            }
         }

      // calculate area
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpX = interpolationTriangleLinear0101(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            area += intpX.J(x0, x1) * w0w1;
      }
      // std::cout << "接触している面の数:" << count << " 表面積:" << area << std::endl;
   };

   // \label{BEM:surfaceIntegralOfTorque}
   Tddd getFroudeKrylovTorque(const Tddd &COM) {
      this->torque = {0., 0., 0.};
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpP = interpolationTriangleLinear0101(P012);
         auto intpX = interpolationTriangleLinear0101(X012);
         auto n = TriangleNormal(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            this->torque += Cross(intpX(x0, x1) - COM, n * intpP(x0, x1)) * intpX.J(x0, x1) * w0w1;
      }
      return this->torque;
   };

   // \label{BEM:surfaceIntegralOfPressure}
   Tddd surfaceIntegralOfPressure() {
      this->force.fill(0.);
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto [p0, p1, p2] = P012;
         auto [X0, X1, X2] = X012;
         this->force += 1. / 6. * (p0 + p1 + p2) * Cross(X1 - X0, X2 - X0);
      }
      return this->force;
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
\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} {\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left( {\sum\limits_{j=0}^2 {{{\left( {{\phi_n}} \right)}_{k_\vartriangle,j }}{N_{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}\left\|\frac{{\partial{\bf{x}}}}{{\partial{\xi_0}}} \times \frac{{\partial{\bf{x}}}}{{\partial{\xi_1}}}\right\|} \right)} }=
```
```math
\alpha_{i_\circ}(\phi)_{i_\circ}-\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} \sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{\bf{x}(\pmb{\xi})-{{\bf x}_{i_\circ} }}{{{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}} \cdot\left(\frac{{\partial {\bf{x}}}}{{\partial {\xi_0}}}\times\frac{{\partial {\bf{x}}}}{{\partial {\xi_1}}}\right)}\right)}
```

ここで，$`\phi_{k_\vartriangle,j}`$における$`k_\vartriangle`$は三角形要素の番号，$`j`$は三角形要素の頂点番号．
$`N_j`$は三角形要素の形状関数，$`\pmb{\xi}`$は三角形要素の内部座標，$`w_0,w_1`$はGauss-Legendre積分の重み，$`\alpha_{i_\circ}`$は原点$`i_\circ`$における立体角，$`\phi`$はポテンシャル，$`\phi_n`$は法線方向のポテンシャル，$`\bf{x}`$は空間座標，$`{\bf x}_{i_\circ}`$は原点の空間座標である．

#### 線形三角要素

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
\frac{\partial {\bf{x}}}{\partial {\xi_0}} \times \frac{\partial {\bf{x}}}{\partial {\xi_1}} = (1-\xi_0) ((p_1-p_0)\times(p_2-p_0))
```

これを使えば，BIEは次のように簡単になる．

```math
\sum\limits_{k_\vartriangle}{2{\bf n}_{k_\vartriangle}}
\sum\limits_{{\xi_1},{w_1}}
{\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left( {\sum\limits_{j=0}^2 {{{\left( {{\phi_n}} \right)}_{k_\vartriangle,j }}{N_{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}
(1-\xi_0)
} \right)} }=
```
```math
\alpha_{i_\circ}(\phi)_{i_\circ}
-\sum\limits_{k_\vartriangle}{2{\bf n}_{k_\vartriangle}}\cdot
\sum\limits_{{\xi_1},{w_1}}
\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{\bf{x}(\pmb{\xi})-{{\bf x}_{i_\circ} }}{{{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}}
(1-\xi_0)
}\right)}
```

$`(1-\xi_0)`$は，必ず正の値をとるので，絶対値を取る必要はない．
$`{2{\bf n}_{k_\vartriangle}} = ((p_1-p_0)\times(p_2-p_0))`$

NOTE: ちなみに，$`\frac{1-\xi_0}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}`$の分子に$`1-\xi_0`$があることで，
関数の特異的な変化を抑えることができる．プログラム上ではこの性質が利用できるように，この二つをまとめて計算する．

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
      return PBF_index.at(pf2ID(p, f));
   };
   /**
   isNeumannID_BEMとisDirichletID_BEMの両方を満たす{p,f}は存在しない．
   */

   void setIGIGn(Network &water) {
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
      water.setFacesVector();
      water.setPointsVector();
#pragma omp parallel
      for (const auto &[PBF, index] : PBF_index)
#pragma omp single nowait
      {
         auto [origin, _] = PBF;
         double origin_ign_rigid_mode = 0.;
         auto &IGIGn_Row = IGIGn[index];
         double nr, tmp;
         std::array<double, 2> IGIGn, c;
         std::array<double, 3> A, cross, N012, R;
         //
         double sum_area = 0;
         int n = 0;
         for (const auto &f : origin->getFaces()) {
            sum_area += f->area;
            n++;
         }
         double r = std::sqrt(sum_area / M_PI);
         std::array<std::tuple<networkPoint *, networkFace *, std::array<double, 2>>, 3> ret;
         for (const auto &integ_f : water.getFacesVector()) {
            const auto [p0, p1, p2] = integ_f->getPoints(origin);
            ret = {{{p0, integ_f, {0., 0.}}, {p1, integ_f, {0., 0.}}, {p2, integ_f, {0., 0.}}}};
            cross = Cross(p0->X - p2->X, p1->X - p2->X);

            if ((Norm(integ_f->center - origin->X) > 10 * r))
               for (const auto &[t0, t1, ww] : __array_GW5xGW5__) {
                  N012 = ModTriShape<3>(t0, t1);
                  tmp = (1. - t0) / (nr = Norm(R = (std::get<0>(N012) * p0->X + std::get<1>(N012) * p1->X + std::get<2>(N012) * p2->X - origin->X)));
                  tmp *= ww;
                  IGIGn = {tmp, (tmp * (Dot(-R / nr, cross))) / nr};
                  std::get<2>(std::get<0>(ret)) += IGIGn * std::get<0>(N012);  // 補間添字0
                  std::get<2>(std::get<1>(ret)) += IGIGn * std::get<1>(N012);  // 補間添字1
                  std::get<2>(std::get<2>(ret)) += IGIGn * std::get<2>(N012);  // 補間添字2
               }
            else
               for (const auto &[t0, t1, ww] : __array_GW13xGW13__) {
                  N012 = ModTriShape<3>(t0, t1);
                  tmp = (1. - t0) / (nr = Norm(R = (std::get<0>(N012) * p0->X + std::get<1>(N012) * p1->X + std::get<2>(N012) * p2->X - origin->X)));
                  tmp *= ww;
                  IGIGn = {tmp, (tmp * (Dot(-R / nr, cross))) / nr};
                  std::get<2>(std::get<0>(ret)) += IGIGn * std::get<0>(N012);  // 補間添字0
                  std::get<2>(std::get<1>(ret)) += IGIGn * std::get<1>(N012);  // 補間添字1
                  std::get<2>(std::get<2>(ret)) += IGIGn * std::get<2>(N012);  // 補間添字2
               }
            /* -------------------------------------------------------------------------- */
            c = {Norm(cross), origin == p0 ? 0. : 1.};
            for (auto &[_, __, igign] : ret)
               igign *= c;

            for (const auto &[p, which_side_f, igign] : ret) {
               IGIGn_Row[pf2Index(p, which_side_f)] += igign;   // この面に関する積分において，φまたはφnの寄与
               if (p != origin)                                 // for use_rigid_mode
                  origin_ign_rigid_mode -= std::get<1>(igign);  // for use_rigid_mode
            }
         }
         /* -------------------------------------------------------------------------- */
         /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

         ### リジッドモードテクニック

         全て$`\phi=1`$とすると，$`\alpha({\bf a}) = -\int\int{\nabla G({\bf x},{\bf a})\cdot{\bf n}({\bf x})dS}`$となり，これを離散化すると，数値積分による評価が難しかった係数行列の対角成分がより精確に計算できる．
         これはリジッドモードテクニックと呼ばれている．
         $`{\bf x}_{i\circ}`$が$`{\bf x}({\pmb \xi})`$に近い場合，$`G`$は急激に特異的に変化するため，数値積分精度が悪化するが，リジッドモードテクニックによって積分を回避できる．

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

   void solve(Network &water, const Buckets<networkPoint *> &FMM_BucketsPoints, const Buckets<networkFace *> &FMM_BucketsFaces) {

      setPhiPhin(water);

      PBF_index.clear();
      PBF_index.reserve(3 * water.getPoints().size());
      int i = 0;
      for (const auto &q : water.getPoints())
         for (const auto &id : variableIDs(q))
            if (PBF_index.find(id) == PBF_index.end())
               PBF_index[id] = i++;

      std::cout << Red << "   unknown size : " << PBF_index.size() << colorReset << std::endl;
      std::cout << Red << "water node size : " << water.getPoints().size() << colorReset << std::endl;

      setIGIGn(water);

      std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;

      makeMatrix();

      knowns.resize(PBF_index.size());
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, f] = PBF;
         if (isDirichletID_BEM(PBF) && isNeumannID_BEM(PBF))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "isDirichletID_BEM(P,F) && isNeumannID_BEM(P,F)");
         else if (isDirichletID_BEM(PBF))
            knowns[i] = p->phi_Dirichlet = std::get<0>(p->phiphin);
         else if (isNeumannID_BEM(PBF))
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

      storePhiPhin(water, ans);

      std::cout << Green << "Elapsed time for solving BIE: " << Red << watch() << colorReset << " s\n";

      isSolutionFinite(water);
   };

   // b! ------------------------------------------------------------------------------ */
   // b!                             solve phi_t and phi_n_t                            */
   // b! ------------------------------------------------------------------------------ */

   /*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

   ## 浮体動揺解析

   BEM-MELで浮体動揺解析ができるようにするのは簡単ではない．
   浮体に掛かる圧力の計算に必要な$`\phi_t`$が簡単には求まらないためである．
   これに関しては，\cite{Wu2003}が参考になる．

   ### 浮体の運動方程式

   <img src="schematic_float.png" width="400px" />

   浮体の重心の運動方程式：

   ```math
   m \frac{d {\boldsymbol U}_{\rm c}}{d t} = \boldsymbol{F}_{\text {ext }}+\boldsymbol{F}_{\text {hydro }}, \quad
   \boldsymbol{I} \frac{d {\boldsymbol \Omega}_{\rm c}}{d t} = \boldsymbol{T}_{\text {ext }}+\boldsymbol{T}_{\text {hydro }}
   ```

   $`{\boldsymbol U}_{\rm c}`$は浮体の移動速度．
   $`\boldsymbol{F}_{\text {ext }}`$は重力などの外力，$`\boldsymbol{F}_{\text {hydro }}`$は水の力，$`\boldsymbol{T}_{\text {ext }}`$は外力によるトルク，$`\boldsymbol{T}_{\text {hydro }}`$は水の力によるトルク．
   浮体が流体から受ける力$`\boldsymbol{F}_{\text {hydro }}`$は，浮体表面の圧力$`p`$を積分することで得られ，
   また圧力$`p`$は速度ポテンシャル$`\phi`$を用いて，以下のように書ける．

   \ref{BEM:surfaceIntegralOfPressure}{圧力積分}と
   \ref{BEM:surfaceIntegralOfTorque}{トルクの積分}：

   ```math
   \boldsymbol{F}_{\text {hydro }}=\iint_{\Gamma_{\rm float}} p\boldsymbol{n}  d S, \quad
   \boldsymbol{T}_{\text {hydro }}=\iint_{\Gamma_{\rm float}} ({\bf x}-{\bf x}_{\rm c})\times (p\boldsymbol{n})  d S, \quad
   p= p({\bf x}) =-\rho\left(\frac{\partial \phi}{\partial t}+\frac{1}{2} (\nabla \phi)^{2}+g z\right)
   ```

   $`\frac{\partial \phi}{\partial t}`$を$`\phi_t`$と書くことにする．この$`\phi_t`$は陽には求められない．
   そこで，$`\phi`$と似た方法，BIEを使った方法で$`\phi_t`$を求める．$`\phi`$と$`\phi_n`$の間に成り立つ境界積分方程式と全く同じ式が，$`\phi_t`$と$`\phi_{nt}`$の間にも成り立つ：

   ```math
   \alpha ({\bf{a}})\phi_t ({\bf{a}}) = \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi_t ({\bf{x}}) - \phi_t ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
   \quad\text{on}\quad{\bf x} \in \Gamma(t).
   ```

   */

   /*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

   ### $`\phi_t`$と$`\phi_{nt}`$に関するBIEの解き方（と$`\phi_{nt}`$の与え方）

   $`\phi_t`$と$`\phi_{nt}`$に関するBIEを解くためには，ディリクレ境界には$`\phi_t`$を，ノイマン境界には$`\phi_{nt}`$を与える．

   #### ディリクレ節点の$`\phi_{nt}`$の与え方(水面：圧力が既知，$`\phi`$が既知)

   このディリクレ境界では，圧力が与えられていないので，このBiEにおいては，ノイマン境界条件を与える．
   ただし，壁が完全に固定されている場合，$`\phi_{nt}`$は0とする．

   #### ディリクレ節点の$`\phi_{t}`$の与え方($`\phi`$を与える造波装置：圧力が未知，$`\phi`$が既知)

   ディリクレ境界では$`\phi_t`$は，圧力が大気圧と決まっているので，ベルヌーイの圧力方程式から$`\phi_t`$を求めることができる．

   #### ノイマン節点での$`\phi_{nt}`$の与え方

   境界面が静止しているかどうかに関わらず，流体と物体との境界では，境界法線方向速度が一致する．
   境界面上の点の位置ベクトルを$`\boldsymbol r`$とする．
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
   \frac{d^2\boldsymbol r}{dt^2} = \frac{d}{dt}\left({\boldsymbol U}_{\rm c} + \boldsymbol \Omega_{\rm c} \times \boldsymbol r\right),\quad \frac{d{\bf n}}{dt} = {\boldsymbol \Omega}_{\rm c}\times{\bf n}
   ```

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

   */

   /*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

   ```math
   \nabla\otimes{\bf u} = \nabla \otimes \nabla \phi =
   \begin{bmatrix} \phi_{xx} & \phi_{xy} & \phi_{xz} \\
   　　　　　　　　　　\phi_{yx} & \phi_{yy} & \phi_{yz} \\
   　　　　　　　　　　\phi_{zx} & \phi_{zy} & \phi_{zz}
   \end{bmatrix}
   ```

   ヘッセ行列の計算には，要素における変数の勾配の接線成分を計算する\ref{BEM:HessianOfPhi}{`HessianOfPhi`}を用いる．
   節点における変数を$`v`$とすると，$`\nabla v-{\bf n}({\bf n}\cdot\nabla v)`$が計算できる．
   要素の法線方向$`{\bf n}`$が$`x`$軸方向$`{(1,0,0)}`$である場合，$`\nabla v - (\frac{\partial}{\partial x},0,0)v`$なので，
   $`(0,\frac{\partial v}{\partial y},\frac{\partial v}{\partial z})`$が得られる．

   */

   /*DOC_EXTRACT 0_4_1_FLOATING_BODY_SIMULATION

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
   V_d Func(const auto &ACCELS_IN, const Network &water, const std::vector<Network *> &rigidbodies) {
      TimeWatch watch;
      auto ACCELS = ACCELS_IN;

      // {
      //    int i = 0;
      //    for (const auto &net : rigidbodies)
      //       if (calculatePhintQ(net))
      //          std::ranges::for_each(net->acceleration, [&](auto &a_w) { a_w = ACCELS_IN[i++]; });
      // }

      //* --------------------------------------------------- */
      //*                  加速度 --> phiphin_t                */
      //* --------------------------------------------------- */
      setPhiPhin_t();
      // std::cout << Green << "setPhiPhin_t()" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

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
#if defined(use_CG)
      GradientMethod gd(mat_ukn);
      ans = gd.solve(ParallelDot(mat_kn, knowns));
#elif defined(use_gmres)
      std::cout << "use gmres for phiphin_t" << std::endl;
      auto b = ParallelDot(mat_kn, knowns);
      gmres gm(mat_ukn, b, ans, use_gmres);
      std::cout << Red << "gmres for phiphin_t. error = " << gm.err << colorReset << std::endl;
      if (!isFinite(gm.err) || !isFinite(gm.x)) {
         std::cout << Red << "gmres for phiphin_t failed" << std::endl;
         if (this->lu) {
            this->lu->solve(b, ans);
         } else {
            std::cout << Green << "new lapack_lu" << colorReset << std::endl;
            this->lu = new lapack_lu(mat_ukn, b, ans);
         }
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
#elif defined(use_lapack)
      this->lu->solve(ParallelDot(mat_kn, knowns) /*既知のベクトル（右辺）*/, ans /*解*/);
#endif
      // std::cout << Green << "solve by LU" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      //@ -------------------------------------------------------------------------- */
      //@                    update p->phiphin_t and p->phinOnFace                   */
      //@ -------------------------------------------------------------------------- */

      storePhiPhin_t(water, ans);
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
      //*                     圧力 ---> 加速度                  */
      //* --------------------------------------------------- */
      int i = 0;
      for (const auto &net : rigidbodies)
         if (net->isFloatingBody) {
            // std::cout << net->inputJSON.find("velocity") << std::endl;
            // std::cout << net->inputJSON["velocity"][0] << std::endl;
            auto tmp = calculateFroudeKrylovForce(water.getFaces(), net);
            auto [_, __, ___, Ix, Iy, Iz] = net->getInertiaGC();
            /* ------------------------------- 係留策 2023/10/11 ------------------------------- */
            std::array<double, 3> F_mooring = {0., 0., 0.}, T_mooring = {0., 0., 0.};
            net->inputJSON.find("mooring", [&](const auto &values) {
               if (values[0] == "simple_mooring") {
                  std::cout << Yellow << "simple_moorings" << colorReset << std::endl;
                  std::array<double, 3> X_anchor = {stod(values[1]), stod(values[2]), stod(values[3])};            //$ アンカー
                  std::array<double, 3> X_fairlead_initial = {stod(values[4]), stod(values[5]), stod(values[6])};  //$ フェアリード
                  double natural_length = Norm(X_anchor - X_fairlead_initial);
                  std::array<double, 3> X_fairlead = rigidTransformation(net->ICOM, net->COM, net->Q.R(), X_fairlead_initial);
                  double current_length = Norm(X_anchor - X_fairlead);
                  F_mooring = stod(values[7]) * (current_length - natural_length) * Normalize(X_anchor - X_fairlead);
                  T_mooring = Cross(X_fairlead - net->COM, F_mooring);
                  std::cout << Yellow << "F_mooring = " << F_mooring << ", T_mooring = " << T_mooring << colorReset << std::endl;
               }
            });
            /* -------------------------------------------------------------------------- */
            auto F_ext = _GRAVITY3_ * net->getMass3D();
            auto F_hydro = tmp.surfaceIntegralOfPressure();
            auto F = F_hydro + F_ext + F_mooring;
            auto T_hydro = tmp.getFroudeKrylovTorque(net->COM);
            auto T = T_hydro + T_mooring;
            auto [a0, a1, a2] = F / net->getMass3D();
            auto [a3, a4, a5] = T / Tddd{Ix, Iy, Iz};
            std::ranges::for_each(T6d{a0, a1, a2, a3, a4, a5}, [&](const auto &a_w) { ACCELS[i++] = a_w; });  // 複数浮体がある場合があるので．
            // write out details of the body
            // std::cout << Green << "mass = " << net->mass << std::endl;
            // std::cout << Green << "inertia = " << net->getInertiaGC() << std::endl;
         } else
            i += 6;
      // std::cout << Green << "other" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      return ACCELS - ACCELS_IN;
   };

   /* -------------------------------------------------------------------------- */

   //@ --------------------------------------------------- */
   //@        加速度 --> phiphin_t --> 圧力 --> 加速度        */
   //@ --------------------------------------------------- */
   void solveForPhiPhin_t(Network &water, const std::vector<Network *> &rigidbodies) {
      water.setGeometricProperties();
      auto ACCELS_init = initializeAcceleration(rigidbodies);

      if (ACCELS_init.empty()) {
         setPhiPhin_t();
         return;
      }

      int count = 0;
      double alpha = 1.;

      BroydenMethod BM(ACCELS_init, ACCELS_init);

      insertAcceleration(rigidbodies, BM.X - BM.dX);
      auto func_ = Func(BM.X - BM.dX, water, rigidbodies);

      for (auto j = 0; j < 1000; ++j) {

         auto func = Func(BM.X, water, rigidbodies);
         BM.update(func, func_, j == 0 ? 0.1 : alpha);
         func_ = func;
         insertAcceleration(rigidbodies, BM.X);

         std::cout << "j = " << j << ", alpha = " << alpha << ", Norm(func) = " << Norm(func) << ", " << Red << "Norm(BM.dX) = " << Norm(BM.dX) << colorReset << std::endl;

         if (Norm(BM.dX) < 1E-9 && Norm(func) < 1E-9) {
            if (count++ > 4)
               break;
         } else
            count = 0;
         if (Norm(BM.dX) < 1E-12 && Norm(func) < 1E-12)
            break;
      }
   };
};

#endif