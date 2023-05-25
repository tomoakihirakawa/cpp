#ifndef BEM_solveBVP_H
#define BEM_solveBVP_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

// #define solve_equations_on_all_points
// #define solve_equations_on_all_points_rigid_mode
// #define solveBVP_debug

// #define use_CG
// #define use_gmres 20
#define use_lapack

// #define quad_element
// #define linear_element
// #define liear_and_quad_element
std::unordered_map<std::tuple<netP *, netF *>, int> PBF_index;

struct calculateFroudeKrylovForce {
   std::vector<networkFace *> actingFaces;
   Tddd force, force_check, torque;
   double area;
   T6d acceleration;
   std::vector<std::tuple<Tddd, T3Tddd>> PressureVeticies;

   calculateFroudeKrylovForce(const std::unordered_set<networkFace *> faces /*waterfaces*/, const Network *PasObj)
       : force({0., 0., 0.}),
         torque({0., 0., 0.}),
         area(0.),
         PressureVeticies({}),
         acceleration({0., 0., 0., 0., 0., 0.}) {
      // PasObjと接したfaceの頂点にpressureが設定されている前提
      int count = 0;
      for (const auto &f : faces)
         if (f->Neumann) {
            if (std::ranges::all_of(f->getPoints(),
                                    [&](const auto &p) { return std::ranges::any_of(p->getContactFaces(), [&](const auto &F) { return F->getNetwork() == PasObj; }); })) {
               auto [p0, p1, p2] = f->getPoints();
               this->PressureVeticies.push_back({{p0->pressure, p1->pressure, p2->pressure}, ToX(f)});
               this->actingFaces.emplace_back(f);
               count++;
            }
         }
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpX = interpolationTriangleLinear0101(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            area += intpX.J(x0, x1) * w0w1;
      }
      std::cout << "接触している面の数:" << count << std::endl;
      std::cout << "表面積:" << area << std::endl;
   };

   Tddd getFroudeKrylovTorque(const Tddd &COM) {
      /*
      crossの引数の順番に注意
      モーメントの計算が，N=rxP
      */
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

   Tddd surfaceIntegralOfPressure() {
      this->force.fill(0.);
      this->force_check.fill(0.);
      for (const auto &[P012, X012] : this->PressureVeticies) {
         // auto intpP = interpolationTriangleLinear0101(P012);
         // auto intpX = interpolationTriangleLinear0101(X012);
         // auto n = TriangleNormal(X012);
         // for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple) {
         //    this->force += n * intpP(x0, x1) * intpX.J(x0, x1) * w0w1;
         //    force_check += intpP(x0, x1) * intpX.cross(x0, x1) * w0w1;
         // }
         //
         {
            auto [p0, p1, p2] = P012;
            auto [X0, X1, X2] = X012;
            this->force += 1. / 6. * (p0 + p1 + p2) * Cross(X1 - X0, X2 - X0);
         }
      }
      force_check -= force;
      return this->force;
   };
};

bool isNeumannID_BEM(const auto &p, const auto &f) {
   if (p->Neumann || p->CORNER) {
      if (p->isMultipleNode) {
         if (p->MemberQ(f))
            return f->Neumann;
         else
            return false;
      } else
         return (f == nullptr);
   } else
      return false;
};

bool isNeumannID_BEM(const std::tuple<netP *, netF *> &PF) { return isNeumannID_BEM(std::get<0>(PF), std::get<1>(PF)); };

bool isDirichletID_BEM(const auto &p, const auto &f) {
   if (p->Dirichlet || p->CORNER)
      return (f == nullptr);
   else
      return false;
};

bool isDirichletID_BEM(const std::tuple<netP *, netF *> &PF) { return isDirichletID_BEM(std::get<0>(PF), std::get<1>(PF)); };

std::tuple<networkPoint *, networkFace *> pf2ID(const networkPoint *p, const networkFace *f) {
   /**
   NOTE: non-multiple node ID is {p,nullptr}
   NOTE: Iterating over p->getFaces() and p may not get all IDs since p->getFaces() doesn't contain nullptr which is often used for an ID of a non-multiple node.
    */
   if (f == nullptr || !p->isMultipleNode || f->Dirichlet)
      return {const_cast<networkPoint *>(p), nullptr};
   else
      return {const_cast<networkPoint *>(p), const_cast<networkFace *>(f)};
}

std::unordered_set<std::tuple<networkPoint *, networkFace *>> variableIDs(const networkPoint *p) {
   //{p,f}を変換
   // f cannot be nullptr
   //  {p,f} --o--> {p,nullptr}
   //  {p,f} <--x-- {p,nullptr}

   std::unordered_set<std::tuple<networkPoint *, networkFace *>> ret;
   for (const auto &f : p->getFaces())
      ret.emplace(pf2ID(p, f));
   return ret;
};

void setPhiPhin(Network &water) {
   /**
   \phi on Dirichlet nodes have been updated by RK method. \phi_n on Neumann nodes are calculated in this function.
   */
   /* -------------------------------------------------------------------------- */
   /*                         phinOnFace, phintOnFaceの設定                         */
   /* -------------------------------------------------------------------------- */
   // b! 点
   std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << colorOff << std::endl;

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

      setNeumann(nullptr);
      for (const auto &f : p->getFaces())
         setNeumann(f);

      setDirichlet(nullptr);
      for (const auto &f : p->getFaces())
         setDirichlet(f);
   }

   // b! 面
   std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << colorOff << std::endl;
#pragma omp parallel
   for (const auto &f : water.getFaces())
#pragma omp single nowait
   {
      if (f->Neumann) {
         std::get<1>(f->phiphin) = Dot(uNeumann(f), f->normal);
      } else {
         auto [p0, p1, p2] = f->getPoints();
         std::get<0>(f->phiphin) = (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3.;
      }
   }
};

// b@ -------------------------------------------------------------------------- */
// b@                                   BEM_BVP                                  */
// b@ -------------------------------------------------------------------------- */

/*DOC_EXTRACT BEM

## 境界値問題

### BIEの離散化

$\phi$と$\phi_n$に関するBIEは，

$$
\alpha ({\bf{a}})\phi ({\bf{a}}) = \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
$$

これを線形三角要素とGauss-Legendre積分で離散化すると，

$$
\alpha_{i_\circ}(\phi)_{i_\circ}=-\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} {\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left( {\sum\limits_{j=0}^2 {{{\left( {{\phi_n}} \right)}_{k_\vartriangle,j }}{N_{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}\left\|\frac{{\partial{\bf{x}}}}{{\partial{\xi_0}}} \times \frac{{\partial{\bf{x}}}}{{\partial{\xi_1}}}\right\|} \right)} }
$$

$$
-\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} \sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{{\bf x}_{i_\circ}} - {\bf{x}}\left( \pmb{\xi } \right)}}{{{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}} \cdot\left(\frac{{\partial {\bf{x}}}}{{\partial {\xi_0}}}\times\frac{{\partial {\bf{x}}}}{{\partial {\xi_1}}}\right)}\right)}
$$

*/

struct BEM_BVP {
   const bool Neumann = false;
   const bool Dirichlet = true;
#if defined(use_lapack)
   lapack_lu *lu;
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
   BEM_BVP() : lu(nullptr){};
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
      std::cout << Magenta << timer() << colorOff << std::endl;
      /*DOC_EXTRACT BEM

      このループでは，BIEの連立一次方程式の係数行列`IGIGn`を作成する作業を行なっている．
      `IGIGn`は，ある節点$i_\circ$（係数行列の行インデックス）に対する
      他の節点$j_\circ$（係数行列の列インデックス）の影響度合いのようなものである．
      その影響度合いは，他の節点$j_\circ$の所属する要素までの距離や向きによって決まることが離散化された式からわかる．

      | Variable | Description |
      |:--------:|:-----------:|
      | `origin` | 原点となる節点$i_\circ$ |
      | `integ_f` | Element $k_{\triangle}$ |
      | `t0, t1, ww` | Gaussian points and thier wieghts $\xi_0, \xi_1, w_0 w_1$ |
      | `p0, p1, p2` | Node of the element $k_{\triangle}$ |
      | `N012` | Shape function $\pmb{N}_j$ |
      | `IGIGn` | Coefficient matrices of the left and right sides |
      | `nr` | $\| \pmb{x} - \pmb{x}_{i\circ } \|$ |
      | `tmp` | $w_0 w_1 \frac{1 - \xi_0}{\| \pmb{x} - \pmb{x}_{i\circ } \|}$ |
      | `cross` | $\frac{\partial \pmb{x}}{\partial \xi_0} \times \frac{\partial \pmb{x}}{\partial \xi_1}$ |

      */
#pragma omp parallel
      for (const auto &[PBF, index] : PBF_index)
#pragma omp single nowait
      {
         auto [origin, _] = PBF;
         double origin_ign_rigid_mode = 0.;
         auto &IGIGn_Row = IGIGn[index];
         double nr, tmp;
         std::array<double, 2> IGIGn, c;
         std::array<double, 3> X0, X1, X2, A, cross, N012;
         for (const auto &integ_f : water.getFaces()) {
            const auto [p0, p1, p2] = integ_f->getPoints(origin);
            std::array<std::tuple<networkPoint *, networkFace *, std::array<double, 2>>, 3> ret = {{{p0, integ_f, {0., 0.}}, {p1, integ_f, {0., 0.}}, {p2, integ_f, {0., 0.}}}};
            for (const auto &[t0, t1, ww] : __array_GW5xGW5__) {
               N012 = ModTriShape<3>(t0, t1);
               tmp = ww * (1. - t0) / (nr = Norm(N012[0] * p0->X + N012[1] * p1->X + N012[2] * p2->X - origin->X));
               IGIGn = {tmp, tmp / (nr * nr)};
               std::get<2>(std::get<0>(ret)) += IGIGn * std::get<0>(N012);  // 補間添字0
               std::get<2>(std::get<1>(ret)) += IGIGn * std::get<1>(N012);  // 補間添字1
               std::get<2>(std::get<2>(ret)) += IGIGn * std::get<2>(N012);  // 補間添字2
            }
            /* -------------------------------------------------------------------------- */
            cross = Cross(p0->X - p2->X, p1->X - p2->X);
            c = {Norm(cross), Dot(origin->X - p0->X, cross)};
            for (auto &[_, __, igign] : ret)
               igign *= c;

            for (const auto &[p, which_side_f, igign] : ret) {
               IGIGn_Row[pf2Index(p, which_side_f)] += igign;   // この面に関する積分において，φまたはφnの寄与
               if (p != origin)                                 // for use_rigid_mode
                  origin_ign_rigid_mode -= std::get<1>(igign);  // for use_rigid_mode
            }
         }
         /* -------------------------------------------------------------------------- */
#if defined(use_rigid_mode)
         std::get<1>(IGIGn_Row[index]) = origin_ign_rigid_mode;
#else
         /*
         @ ∇^2(1/r)=-1/(4pi)δ(r)
         @ IG*φn=-aφ+IGn*φ
         @ IG*φn=(IGn-a)*φ
         */
         std::get<1>(IGIGn_Row[index]) += origin->getSolidAngle();
#endif
      }
      std::cout << Green << "離散化にかかった時間" << timer() << colorOff << std::endl;
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
            /*DOC_EXTRACT BEM
            IGIGn は 左辺に IG*φn が右辺に IGn*φ が来るように計算しているため，移項する場合，符号を変える必要がある．
            $IG \phi_n = IGn \phi$

            移項前:
            $\begin{bmatrix}IG _0 & IG _1 & IG _2 & IG _3\end{bmatrix} \begin{bmatrix}\phi _{n0} \\ \phi _{n1} \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}IG _{n0} & IG _{n1} & IG _{n2} & IG _{n3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _1 \\ \phi _2 \\ \phi _3\end{bmatrix}$

            移項後:
            $\begin{bmatrix}IG _0 & -IG _{n1} & IG _2 & IG _3\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}IG _{n0} & -IG _1 & IGn _2 & IG _{n3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}$

            多重節点(1と3が多重節点の場合):
            $\begin{bmatrix}0 & 1 & 0 & 0\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}0 & 0 & 0 & 1\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}$

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

   void storePhiPhin(Network &water, const auto &ans) const {
      std::cout << "store ans" << std::endl;
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, f] = PBF;
         if (isDirichletID_BEM(PBF)) {
            p->phinOnFace.at(f) = std::get<1>(p->phiphin) = p->phin_Dirichlet = ans[i];
            p->phiOnFace.at(f) = std::get<0>(p->phiphin);
         }
         if (isNeumannID_BEM(PBF))
            p->phiOnFace.at(f) = std::get<0>(p->phiphin) = ans[i];
      }
      std::cout << "Neumannに限りphiを代入." << std::endl;
      for (const auto &p : water.getPoints())
         if (p->Neumann) {
            double total = 0;
            std::get<0>(p->phiphin) = 0;
            std::ranges::for_each(p->getFaces(), [&total](const auto &f) { total += f->area; });
            for (const auto &f : p->getFaces()) {
               if (p->phiOnFace.count(f))
                  std::get<0>(p->phiphin) += p->phiOnFace.at(f) * f->area / total;
               else
                  std::get<0>(p->phiphin) += p->phiOnFace.at(nullptr) * f->area / total;
            }
         }
   }

   void storePhiPhin_t(const auto &water, const auto &ans) const {
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, f] = PBF;
         if (isDirichletID_BEM(PBF)) {
            p->phintOnFace.at(f) = std::get<1>(p->phiphin_t) = ans[i];
            p->phitOnFace.at(f) = std::get<0>(p->phiphin_t);
         }
         if (isNeumannID_BEM(PBF))
            p->phitOnFace.at(f) = std::get<0>(p->phiphin_t) = ans[i];
      }

      //\phi_tを代入. Neumannに限る
      for (const auto &p : water->getPoints())
         if (p->Neumann) {
            double total = 0;
            std::get<0>(p->phiphin_t) = 0;
            std::ranges::for_each(p->getFaces(), [&total](const auto &f) { total += f->area; });
            for (const auto &f : p->getFaces()) {
               if (p->phitOnFace.count(f))
                  std::get<0>(p->phiphin_t) += p->phitOnFace.at(f) * f->area / total;
               else
                  std::get<0>(p->phiphin_t) += p->phitOnFace.at(nullptr) * f->area / total;
            }
         }
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
      std::cout << Red << "total = " << PBF_index.size() << colorOff << std::endl;

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

      ans.resize(knowns.size());

      if (this->lu)
         delete this->lu;
#if defined(use_CG)
      this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/);
      std::cout << "The conjugate gradient is used" << std::endl;
      GradientMethod gd1(mat_ukn);
      ans = gd1.solve(Dot(mat_kn, knowns), {}, 1E-1);
      GradientMethod gd2(mat_ukn);
      ans = gd2.solveCG(Dot(mat_kn, knowns), ans);
#elif defined(use_gmres)
      std::cout << "gmres for ans" << std::endl;
      gmres gm(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, ans /*解*/, use_gmres);
      ans = gm.x;
      std::cout << "gm.err = " << gm.err << ", isFinite(gm.err) = " << isFinite(gm.err) << std::endl;
      if (real_time < 0.005 || !isFinite(gm.err < 1E-20)) {
         std::cout << "lapack lu decomposition" << std::endl;
         this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, ans /*解*/);
      }
         // auto err = Norm(Dot(mat_ukn, gm.x) - Dot(mat_kn, knowns));
         // std::cout << err << std::endl;
#elif defined(use_lapack)
      std::cout << "lapack lu decomposition" << std::endl;
      this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, ans /*解*/);
#endif

      std::cout << colorOff << "update p->phiphin and p->phinOnFace for Dirichlet boundary" << colorOff << std::endl;

      storePhiPhin(water, ans);

      std::cout << Green << "Elapsed time for solving BIE: " << Red << watch() << colorOff << " s\n";

      isSolutionFinite(water);
   };

   // b! ------------------------------------------------------------------------------ */
   // b!                             solve phi_t and phi_n_t                            */
   // b! ------------------------------------------------------------------------------ */

   /*DOC_EXTRACT BEM

   ## 浮体動揺解析

   浮体の重心の運動方程式：

   $$
   m \frac{d {\boldsymbol U}_{\rm c}}{d t} = \boldsymbol{F}_{\text {ext }}+\boldsymbol{F}_{\text {hydro }}, \quad
   \boldsymbol{I} \frac{d {\boldsymbol \Omega}_{\rm c}}{d t} = \boldsymbol{T}_{\text {ext }}+\boldsymbol{T}_{\text {hydro }}
   $$

   ${\boldsymbol U}_{\rm c}$は浮体の移動速度．
   $\boldsymbol{F}_{\text {ext }}$は重力などの外力，$\boldsymbol{F}_{\text {hydro }}$は水の力，$\boldsymbol{T}_{\text {ext }}$は外力によるトルク，$\boldsymbol{T}_{\text {hydro }}$は水の力によるトルク．
   浮体が流体から受ける力$\boldsymbol{F}_{\text {hydro }}$は，浮体表面の圧力$p$を積分することで得られ，
   また圧力$p$は速度ポテンシャル$\phi$を用いて，以下のように書ける．

   $$
   \boldsymbol{F}_{\text {hydro }}=\int_{S} p\boldsymbol{n}  d S, \quad
   p=-\rho\left(\frac{\partial \phi}{\partial t}+\frac{1}{2} (\nabla \phi)^{2}+g z\right)
   $$

   $\frac{\partial \phi}{\partial t}$を$\phi_t$と書くことにする．この$\phi_t$は陽には求められない．
   そこで，$\phi$と似た方法，BIEを使った方法で$\phi_t$を求める．$\phi$と$\phi_n$の間に成り立つ境界積分方程式と全く同じ式が，$\phi_t$と$\phi_{nt}$の間にも成り立つ：

   $$
   \alpha ({\bf{a}})\phi_t ({\bf{a}}) = \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi_t ({\bf{x}}) - \phi_t ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
   \quad\text{on}\quad{\bf x} \in \Gamma(t).
   $$

   境界面が静止しているかどうかに関わらず，流体と物体との境界では，境界法線方向速度が一致する．
   境界面上の位置ベクトルを$\boldsymbol r$とする．
   表面上のある点の移動速度$\frac{d\boldsymbol r}{dt}$と流体粒子の流速$\nabla \phi$の間には，次の境界条件が成り立つ．

   $$
   {\bf n}\cdot\frac{d\boldsymbol r}{dt} =  {\bf n} \cdot \nabla \phi
   $$

   これを微分することで，$\phi_{nt}$を$\phi$と加速度$\frac{d{\boldsymbol U}_{\rm c}}{dt}$と角加速度$\frac{d{\boldsymbol \Omega}_{\rm c}}{dt}$を使って表すことができる．
   [Wu (1998)](https://www.sciencedirect.com/science/article/pii/S088997469890158X)

   $$
   \begin{aligned}
   &\rightarrow& \frac{d}{dt}({{\bf n}\cdot\frac{d\boldsymbol r}{dt}}) & = \frac{d}{dt}({{\bf n} \cdot \nabla \phi})\\
   &\rightarrow& \frac{d{\bf n}}{dt}\cdot\frac{d\boldsymbol r}{dt} + {\bf n}\cdot\frac{d^2{\boldsymbol r}}{dt^2} & = \frac{d{\bf n}}{dt} \cdot \nabla \phi + {\bf n} \cdot \frac{d}{dt}{\nabla \phi}\\
   &\rightarrow& \frac{d{\bf n}}{dt}\cdot{(\frac{d\boldsymbol r}{dt} - \nabla \phi)} & ={\bf n} \cdot \left(\frac{d}{dt}{\nabla \phi}- \frac{d^2{\boldsymbol r}}{dt^2}\right)\\
   &\rightarrow& \frac{d{\bf n}}{dt}\cdot{(\frac{d\boldsymbol r}{dt} - \nabla \phi)} & ={\bf n} \cdot \left(\phi_t + \nabla \phi\cdot \nabla\nabla \phi - \frac{d^2{\boldsymbol r}}{dt^2}\right)\\
   &\rightarrow& \phi_{nt} &= \frac{d{\bf n}}{dt} \cdot{(\frac{d\boldsymbol r}{dt} - \nabla \phi)} -{\bf n} \cdot \left(\nabla \phi\cdot \nabla\nabla \phi -\frac{d^2{\boldsymbol r}}{dt^2}\right)
   \end{aligned}
   $$

   ここの$\frac{d{\bf n}}{dt}$と$\frac{d^2\boldsymbol r}{dt^2}$は，${\boldsymbol U}_{\rm c}$と$\boldsymbol \Omega_{\rm c}$を用いて，

   $$
   \frac{d^2\boldsymbol r}{dt^2} = \frac{d}{dt}\left({\boldsymbol U}_{\rm c} + \boldsymbol \Omega_{\rm c} \times \boldsymbol r\right),\quad \frac{d{\bf n}}{dt} = {\boldsymbol \Omega}_{\rm c}\times{\bf n}
   $$

   $\frac{d^2\boldsymbol r}{dt^2}$を上の式に代入し，$\phi_{nt}$を求め，
   次にBIEから$\phi_t$を求め，次に圧力$p$を求める．
   そして，浮体の重さと慣性モーメントを考慮して圧力から求めた$\frac{d^2\boldsymbol r}{dt^2}$は，
   入力した$\frac{d^2\boldsymbol r}{dt^2}$と一致しなければならない．

   現状を整理すると，この浮体動揺解析において，知りたい未知変数は，浮体の加速度と角加速度だけ．
   しかし，浮体の没水面上にある節点での圧力$p$が得られないと，$\boldsymbol{F}_{\text {hydro }}$が得られず，運動方程式から浮体加速度が計算できない．
   圧力を計算するためには，$\phi_t$が必要で，$\phi_t$は簡単には得られない，という状況．

   物体の加速度は， 節点における$\{\phi_{nt0},\phi_{nt1},\phi_{nt2},..\} = \Phi_{nt}$が分かれば求まるが，
   逆に$\Phi_{nt}$は$\frac{d\boldsymbol U_{\rm c}}{dt}$が分かれば求まるので

   $$
   \begin{align*}
   &&\frac{d\boldsymbol U_{\rm c}}{dt}& = F\left(\Phi_{nt}\left(\frac{d\boldsymbol U_{\rm c}}{dt}\right)\right)\\
   &\rightarrow& Q\left(\frac{d\boldsymbol U_{\rm c}}{dt}\right) &= \frac{d\boldsymbol U_{\rm c}}{dt} - F\left(\Phi_{nt}\left(\frac{d\boldsymbol U_{\rm c}}{dt}\right)\right) =0
   \end{align*}
   $$

   のように，ある関数$Q$のゼロを探す，根探し問題になる．
   $\phi_{nt}$は，\ref{BEM:setphint}{ここ}で与えている．

   */

   void setPhiPhin_t() const {
#ifdef derivatives_debug
      std::cout << "φtとφntを一部計算👇" << std::endl;
#endif

#ifdef _OPENMP
   #pragma omp parallel
#endif
      for (const auto &[PBF, i] : PBF_index)
#ifdef _OPENMP
   #pragma omp single nowait
#endif
      {
         auto [p, F] = PBF;
         //!!ノイマンの場合はこれでDphiDtは計算できませんよ
         if (isDirichletID_BEM(PBF))
            p->phitOnFace.at(F) = std::get<0>(p->phiphin_t) = p->aphiat(0.);
         if (isNeumannID_BEM(PBF)) {

            // b* p->phintOnFaceは，std::unordered_map<networkFace *, double>
            // b* 節点のphinを保存する．また，多重節点かどうかも，面がnullptrかどうかで判別できる．
            // b* setBoundaryConditionsで決めている．

            /*DOC_EXTRACT BEM

            $$
            \nabla {\bf u} = \nabla \nabla \phi =
            \begin{bmatrix} \phi_{xx} & \phi_{xy} & \phi_{xz} \\
            　　　　　　　　　　\phi_{yx} & \phi_{yy} & \phi_{yz} \\
            　　　　　　　　　　\phi_{zx} & \phi_{zy} & \phi_{zz}
            \end{bmatrix}
            $$

            */

            // \label{BEM:setphint}
            for (auto &[f, phin_t] : p->phintOnFace) {
               auto use_face = (f != nullptr);
               if (use_face) {
                  auto n = f->normal;
                  auto netInContact = NearestContactFace(f)->getNetwork();
                  auto Omega = netInContact->velocityRotational();
                  auto U = uNeumann(p, f);
                  auto A = accelNeumann(p, f);
                  phin_t = Dot(Cross(Omega, n), U - p->U_BEM) + Dot(n, A);
                  auto s0s1s2 = OrthogonalBasis(n);
                  auto [s0, s1, s2] = s0s1s2;
                  auto Hessian = grad_U_LinearElement(f, s0s1s2);
                  phin_t -= std::get<0>(Dot(Tddd{{Dot(p->U_BEM, s0), Dot(p->U_BEM, s1), Dot(p->U_BEM, s2)}}, Hessian));
               } else {
                  auto n = p->getNormalNeumann_BEM();
                  auto netInContact = NearestContactFace(p)->getNetwork();
                  auto Omega = netInContact->velocityRotational();
                  auto U = uNeumann(p);
                  auto A = accelNeumann(p);
                  phin_t = Dot(Cross(Omega, n), U - p->U_BEM) + Dot(n, A);
                  auto s0s1s2 = OrthogonalBasis(n);
                  auto [s0, s1, s2] = s0s1s2;
                  auto Hessian = grad_U_LinearElementNeuamnn(p, s0s1s2);
                  phin_t -= std::get<0>(Dot(Tddd{{Dot(p->U_BEM, s0), Dot(p->U_BEM, s1), Dot(p->U_BEM, s2)}}, Hessian));
               }
               std::get<1>(p->phiphin_t) = phin_t;
            }
         }
      }
   };

   /* ------------------------------------------------------ */
   bool isTarget(Network *net) const {
      // return true;
      if (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating")
         // if (net->inputJSON.find("RigidBody"))
         return true;
      else
         return false;
   };

   V_d initializeAcceleration(const std::vector<Network *> &rigidbodies) {
      V_d ACCELS_init;
      for (const auto &net : rigidbodies) {
         if (isTarget(net)) {
            if (net->inputJSON.at("velocity").size() > 1) {
               std::cout << Red << "net->inputJSON[\" velocity \"][1] = " << net->inputJSON.at("velocity")[1] << colorOff << std::endl;
               double start_time = std::stod(net->inputJSON.at("velocity")[1]);
               std::cout << Red << "start_time = " << start_time << colorOff << std::endl;
               if (real_time < start_time)
                  net->acceleration.fill(0.);
               else
                  std::ranges::for_each(net->acceleration, [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
            } else
               std::ranges::for_each(net->acceleration, [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
         } else
            net->acceleration.fill(0.);
      }
      return ACCELS_init;
   }

   void insertAcceleration(const std::vector<Network *> &rigidbodies, const V_d &BM_X) {
      int i = 0;
      for (const auto &net : rigidbodies) {
         if (isTarget(net)) {
            if (net->inputJSON.at("velocity").size() > 1) {
               std::cout << Red << "net->inputJSON[\" velocity \"][1] = " << net->inputJSON.at("velocity")[1] << colorOff << std::endl;
               double start_time = std::stod(net->inputJSON.at("velocity")[1]);
               std::cout << Red << "start_time = " << start_time << colorOff << std::endl;
               if (real_time < start_time) {
                  net->acceleration.fill(0.);
               } else
                  std::ranges::for_each(net->acceleration, [&](auto &a_w) { a_w = BM_X[i++]; });
            } else
               std::ranges::for_each(net->acceleration, [&](auto &a_w) { a_w = BM_X[i++]; });
         } else
            net->acceleration.fill(0.);
      }
   }

   /* -------------------------------------------------------------------------- */
   V_d Func(const auto &ACCELS_IN, const Network *water, const std::vector<Network *> &rigidbodies) {
      auto ACCELS = ACCELS_IN;
      {
         int i = 0;
         for (const auto &net : rigidbodies)
            if (isTarget(net))
               std::ranges::for_each(net->acceleration, [&](auto &a_w) { a_w = ACCELS_IN[i++]; });
      }

      //* --------------------------------------------------- */
      //*                  加速度 --> phiphin_t                */
      //* --------------------------------------------------- */
      setPhiPhin_t();

      knowns.resize(PBF_index.size());
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, f] = PBF;
         if (isDirichletID_BEM(PBF))
            knowns[i] = p->phitOnFace.at(f);
         if (isNeumannID_BEM(PBF))
            knowns[i] = p->phintOnFace.at(f);  // はいってない？はいってた．
      }

      std::cout << "加速度 --> phiphin_t" << std::endl;
      ans.resize(knowns.size());
#if defined(use_CG)
      GradientMethod gd(mat_ukn);
      ans = gd.solve(Dot(mat_kn, knowns));
#elif defined(use_gmres)
      std::cout << "gmres for phiphin_t" << std::endl;
      gmres gm(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, ans /*解*/, use_gmres);
      ans = gm.x;
      if (!isFinite(gm.err)) {
         std::cout << "gm.err = " << gm.err << std::endl;
         this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, ans /*解*/);
      }
#elif defined(use_lapack)
      this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, ans /*解*/);
#endif
      std::cout << "solved" << std::endl;

      //@ -------------------------------------------------------------------------- */
      //@                    update p->phiphin_t and p->phinOnFace                   */
      //@ -------------------------------------------------------------------------- */

      storePhiPhin_t(water, ans);

      //* --------------------------------------------------- */
      //*                 phiphin_t --> 圧力                   */
      //* --------------------------------------------------- */

      for (const auto &[PBF, i] : PBF_index) {
         auto [p, f] = PBF;
         p->pressure = p->pressure_BEM = -_WATER_DENSITY_ * (std::get<0>(p->phiphin_t) + _GRAVITY_ * p->height() + Dot(p->U_BEM, p->U_BEM) / 2.);
         if (isDirichletID_BEM(PBF))
            p->pressure_BEM = 0;
      }

      //* --------------------------------------------------- */
      //*                     圧力 ---> 加速度                  */
      //* --------------------------------------------------- */
      int i = 0;
      for (const auto &net : rigidbodies)
         if (isTarget(net)) {
            // std::cout << net->inputJSON.find("velocity") << std::endl;
            // std::cout << net->inputJSON["velocity"][0] << std::endl;
            auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
            auto [mx, my, mz, Ix, Iy, Iz] = net->getInertiaGC();
            auto force = tmp.surfaceIntegralOfPressure() + _GRAVITY3_ * net->mass;
            std::cout << "force_check:" << tmp.force_check << std::endl;
            auto torque = tmp.getFroudeKrylovTorque(net->COM);
            auto [a0, a1, a2] = force / Tddd{mx, my, mz};
            auto [a3, a4, a5] = torque / Tddd{Ix, Iy, Iz};
            // net->acceleration = T6d{a0, a1, a2, a3, a4, a5};
            std::ranges::for_each(T6d{a0, a1, a2, a3, a4, a5}, [&](const auto &a_w) { ACCELS[i++] = a_w; });
         }
      return ACCELS - ACCELS_IN;
   };

   //@ --------------------------------------------------- */
   //@        加速度 --> phiphin_t --> 圧力 --> 加速度        */
   //@ --------------------------------------------------- */
   void solveForPhiPhin_t(const Network *water, const std::vector<Network *> &rigidbodies) {

      auto ACCELS_init = initializeAcceleration(rigidbodies);

      if (ACCELS_init.empty()) {
         setPhiPhin_t();
         return;
      }

      auto tmp = ACCELS_init;
      tmp[0] += 1E-10;
      BroydenMethod BM(ACCELS_init, tmp);
      for (auto j = 0; j < 20; ++j) {
         auto func_ = Func(BM.X - BM.dX, water, rigidbodies);
         std::cout << "func_ = " << func_ << std::endl;
         auto func = Func(BM.X, water, rigidbodies);
         std::cout << "func = " << func_ << std::endl;
         BM.update(func, func_, j < 1 ? 1E-10 : 1.);

         insertAcceleration(rigidbodies, BM.X);

         std::cout << "j = " << j << ", " << Red << Norm(func) << colorOff << std::endl;
         std::cout << Red << "func_ = " << func_ << colorOff << std::endl;
         std::cout << Red << "func = " << func << colorOff << std::endl;
         std::cout << Red << "BM.X = " << BM.X << colorOff << std::endl;
         std::cout << Red << "BM.dX = " << BM.dX << colorOff << std::endl;

         if (Norm(func) < 1E-10)
            break;
      }
   };
};

#endif