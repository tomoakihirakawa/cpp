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
   /**
    ## 多重節点
    多重節点という名前は具体性に欠ける．
    普通φnは(節点)にのみ依存する変数だが，nの変化が急なため，不連続性が著しい節点においては，(節点に加え面)にも依存する変数を複数設定する．それらは離散化などで使い分けることになる．
    BIEの離散化における，多重節点扱いについて．
    BIEを数値的に解くために，十分な数の１次方程式を作成する．これは，節点と同じ位置にBIEの原点を取ることで実現できる．
    同じ位置であるにもかかわらず，(節点に加え面)にも依存する変数φnを設定した場合，
    同じ位置であるにもかかわらず，それらを一つ一つを原点として，１次方程式を作成する．
    これらは完全に同じ方程式である．変数の数を節点の数よりも増やしたことによって，方程式の数が増えている．
   */
   std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << colorOff << std::endl;

   for (const auto &p : water.getPoints())
      setIsMultipleNode(p);

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

BIEを離散化すると

$$
\left( \phi_n \right)_{k_\vartriangle,j}
$$


$$
\alpha_{i_{\circ}}\left( \phi \right)_{i_\circ} = \sum_{k_\vartriangle}\sum_{{\xi_1}} \sum_{{\xi_0}} \left( \sum_{j=0}^2 \left( \phi_n \right)_{k_\vartriangle,j} N_{j}\left( \bm{\xi } \right) \frac{w_0 w_1}{\| \mathbf{x}\left( \bm{\xi } \right) - \mathbf{x}_{i_\circ} \|} \left\| \frac{\partial\mathbf{x}}{\partial{\xi_0}} \times \frac{\partial\mathbf{x}}{\partial{\xi_1}} \right\| \right)
$$

$$
\begin{aligned}
\alpha_{i_{\circ}}\left( \phi \right)_{i_\circ} &= \sum_{k_\vartriangle}\sum_{{\xi_1}} \sum_{{\xi_0}} \left( \sum_{j=0}^2 \left( \phi_n \right)_{k_\vartriangle,j} N_{j}\left( \bm{\xi } \right) \frac{w_0 w_1}{\| \mathbf{x}\left( \bm{\xi } \right) - \mathbf{x}_{i_\circ} \|} \left\| \frac{\partial\mathbf{x}}{\partial{\xi_0}} \times \frac{\partial\mathbf{x}}{\partial{\xi_1}} \right\| \right) \\
&- \sum_{k_\vartriangle}\sum_{{\xi_1}} \sum_{{\xi_0}} \left( \sum_{j =0}^2 \left( \phi \right)_{k_\vartriangle,j} N_{j}\left( \bm{\xi } \right) w_0 w_1 \frac{\mathbf{x}_{i_\circ} - \mathbf{x}\left( \bm{\xi } \right)}{\| \mathbf{x}\left( \bm{\xi} \right) - \mathbf{x}_{i_\circ}\|^3} \cdot \left(\frac{\partial \mathbf{x}}{\partial {\xi_0}} \times \frac{\partial \mathbf{x}}{\partial {\xi_1}} \right) \right)
\end{aligned}
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
            /* ------------------------------ 2023/04/03 -------------------------------- */
            /*
            このfor loopでは連立一次方程式の計数行列を作成する作業を行なっている．
            これは，原点をある節点(origin)に固定し，originと各面との位置や向きの関係に依存する，値を各節点に分配する作業である．
            */
            const auto [p0, p1, p2] = integ_f->getPoints(origin);
            std::array<std::tuple<networkPoint *, networkFace *, std::array<double, 2>>, 3> ret = {{{p0, integ_f, {0., 0.}},
                                                                                                    {p1, integ_f, {0., 0.}},
                                                                                                    {p2, integ_f, {0., 0.}}}};
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
               /**
                *
               このループでは，
               ある面integ_fに隣接する節点{p0,p1,p2}の列,IGIGn[origin(fixed),p0],...に値が追加されていく．
               （p0が多重接点の場合，適切にp0と同じ位置に別の変数が設定されており，別の面の積分の際にq0が参照される．）
               p0は，{面,補間添字}で決定することもできる．
               {面,補間添字0}->p0,{面,補間添字1}->p1,{面,補間添字2}->p2というように．
               //@ 多重節点：
               {面A,補間添字},{面B,補間添字},{面C,補間添字}が全て同じ節点p0を指していたとする．
               普通の節点なら，IGIGn[origin,{p0,nullptr}]を指す．
               多重節点なら，IGIGn[origin,{p0,面A}],IGIGn[origin,{p0,面B}]を指すようにする．
               この操作を言葉で言い換えると，
               「nが不連続に変化する点では，その点の隣接面にそれぞれ対してφnを求めるべきである（φは同じでも）．」
               「nが不連続に変化する点では，どの面を積分するかに応じて，参照するφnを区別し切り替える必要がある．」
               //
               //@ さて，この段階でp0が多重節点であるかどうか判断できるだろうか？
               {節点，面}-> 列ベクトルのインデックス を決めれるか？
               //
               面を区別するかどうかが先にわからないので，face*のまsまかnullptrとすべきかわからないということ．．．．
               //
               PBF_index[{p, Dirichlet, ある要素}]
               は存在しないだろう．Dirichlet節点は，{p, ある要素}からの寄与を，ある面に
               */
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
            /*
            IGIGn は 左辺に IG*φn が右辺に IGn*φ が来るように計算しているため，移項する場合，符号を変える必要がある．
            IG*φn = IGn*φ
            移項前:{IG0,IG1,IG2,IG3} . {φn0,φn1,φn2,φn3} = {IGn0,IGn1,IGn2,IGn3} . {φ0,φ1,φ2,φ3}
            移項後:{IG0,-IGn1,IG2,IG3} . {φn0,φ1,φn2,φn3} = {IGn0,-IG1,IGn2,IGn3} . {φ0,φn1,φ2,φ3}
            */
            mat_ukn[i][j] = std::get<0>(igign);
            mat_kn[i][j] = std::get<1>(igign);
         }
         /*
            移項前:{IG0,IG1,IG2,IG3} . {φn0,φn1,φn2,φn3} = {IGn0,IGn1,IGn2,IGn3} . {φ0,φ1,φ2,φ3}
            移項後:{IG0,-IGn1,IG2,IG3} . {φn0,φ1,φn2,φn3} = {IGn0,-IG1,IGn2,IGn3} . {φ0,φn1,φ2,φ3}
            多重節点:{0, 1, 0, 0} . {φn0,φ1,φn2,φn3} = {0, 0, 0, 1} . {φ0,φn1,φ2,φ3}
         */
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
            // /* ∇U=∇∇f={{fxx, fyx, fzx},{fxy, fyy, fzy},{fxz, fyz, fzz}}, ∇∇f=∇∇f^T */
            // // b* p->phintOnFaceは，std::unordered_map<networkFace *, double>
            // // b* 節点のphinを保存する．また，多重節点かどうかも，面がnullptrかどうかで判別できる．
            // // b* setBoundaryConditionsで決めている．
            // for (auto &[f, phin_t] : p->phintOnFace) {
            //    auto use_face = (f != nullptr);
            //    auto n = use_face ? f->normal : p->getNormalNeumann_BEM();
            //    auto netInContact = use_face ? NearestContactFace(f)->getNetwork() : NearestContactFace(p)->getNetwork();
            //    auto w = netInContact->velocityRotational();
            //    auto U = uNeumann(p);
            //    phin_t = Dot(w, U - p->U_BEM) + Dot(n, accelNeumann(p));
            //    auto s0s1s2 = OrthogonalBasis(n);
            //    auto [s0, s1, s2] = s0s1s2;
            //    auto Hessian = use_face ? grad_U_LinearElement(f, s0s1s2) : grad_U_LinearElement(p, s0s1s2);
            //    phin_t -= std::get<0>(Dot(Tddd{{Dot(U, s0), Dot(U, s1), Dot(U, s2)}}, Hessian));
            //    std::get<1>(p->phiphin_t) =  phin_t;
            // }

            /* ∇U=∇∇f={{fxx, fyx, fzx},{fxy, fyy, fzy},{fxz, fyz, fzz}}, ∇∇f=∇∇f^T */
            // b* p->phintOnFaceは，std::unordered_map<networkFace *, double>
            // b* 節点のphinを保存する．また，多重節点かどうかも，面がnullptrかどうかで判別できる．
            // b* setBoundaryConditionsで決めている．
            for (auto &[f, phin_t] : p->phintOnFace) {
               auto use_face = (f != nullptr);
               if (use_face) {
                  auto n = f->normal;
                  auto netInContact = NearestContactFace(f)->getNetwork();
                  auto w = netInContact->velocityRotational();
                  auto U = uNeumann(p, f);
                  auto A = accelNeumann(p, f);
                  phin_t = Dot(w, U - p->U_BEM) + Dot(n, A);
                  auto s0s1s2 = OrthogonalBasis(n);
                  auto [s0, s1, s2] = s0s1s2;
                  auto Hessian = grad_U_LinearElement(f, s0s1s2);
                  phin_t -= std::get<0>(Dot(Tddd{{Dot(U, s0), Dot(U, s1), Dot(U, s2)}}, Hessian));
               } else {
                  auto n = p->getNormalNeumann_BEM();
                  auto netInContact = NearestContactFace(p)->getNetwork();
                  auto w = netInContact->velocityRotational();
                  auto U = uNeumann(p);
                  auto A = accelNeumann(p);
                  phin_t = Dot(w, U - p->U_BEM) + Dot(n, A);
                  auto s0s1s2 = OrthogonalBasis(n);
                  auto [s0, s1, s2] = s0s1s2;
                  auto Hessian = grad_U_LinearElementNeuamnn(p, s0s1s2);
                  phin_t -= std::get<0>(Dot(Tddd{{Dot(U, s0), Dot(U, s1), Dot(U, s2)}}, Hessian));
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
   /* ------------------------------------------------------ */
   V_d Func(V_d ACCELS_IN, const Network *water, const std::vector<Network *> &rigidbodies) {
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

   /* ------------------------------------------------------ */

   V_d initializeAcceleration(const std::vector<Network *> &rigidbodies) {
      V_d ACCELS_init;
      for (const auto &net : rigidbodies) {
         if (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating") {
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
         if (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating") {
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