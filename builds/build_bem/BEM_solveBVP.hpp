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
std::unordered_map<std::tuple<netP *, bool, netF *>, int> PBF_index;

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
            if (all_of(f->getPoints(),
                       [&](const auto &p) { return std::any_of(p->getContactFaces().begin(), p->getContactFaces().end(), [&](const auto &F) { return F->getNetwork() == PasObj; }); })) {
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

struct BEM_BVP {
   // std::unordered_set<networkPoint *> Points;
   // std::unordered_set<networkFace *> Faces;
   const bool Neumann = false;
   const bool Dirichlet = true;
#if defined(use_lapack)
   lapack_lu *lu;
#else
   ludcmp_parallel *lu;
#endif

   using T_PBF = std::tuple<netP *, bool, netF *>;
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
   //% ------------------------------------------------------------------------------ */
   //%                             solve phi_t and phi_n_t                            */
   //% ------------------------------------------------------------------------------ */
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
         auto [p, B, F] = PBF;
         //% ------------------------------------------------------ */
         //%                 ディリクレ境界面上のφtを計算                */
         //% ------------------------------------------------------ */
         //!!ノイマンの場合はこれでDphiDtは計算できませんよ
         if (p->Dirichlet || p->CORNER)
            std::get<0>(p->phiphin_t) = p->aphiat(0.);
         // std::get<0>(p->phiphin_t) = p->DphiDt(0.) - Dot(p->U_BEM, p->U_BEM);
         //% ------------------------------------------------------ */
         //%    ノイマン境界面上の加速度から,ノイマン境界面上のφntを計算     */
         //% ------------------------------------------------------ */
         if (p->Neumann || p->CORNER) {
            /* ∇U=∇∇f={{fxx, fyx, fzx},{fxy, fyy, fzy},{fxz, fyz, fzz}}, ∇∇f=∇∇f^T */
            // b* p->phintOnFaceは，std::unordered_map<networkFace *, double>
            // b* 節点のphinを保存する．また，多重節点かどうかも，面がnullptrかどうかで判別できる．
            // b* setBoundaryConditionsで決めている．
            auto n = p->getNormalNeumann_BEM();
            auto Q = Quaternion();
            for (auto &[f, phin_t] : p->phintOnFace) {
               if (f) {
                  auto netInContact = NearestContactFace(f)->getNetwork();
                  auto w = netInContact->velocityRotational();
                  auto dQdt = Q.d_dt(w);
                  auto U = uNeumann(p, f);
                  phin_t = std::get<1>(p->phiphin_t) = Dot(w, U - p->U_BEM) + Dot(n, accelNeumann(p, f));
                  auto s0s1s2 = OrthogonalBasis(f->normal);
                  auto [s0, s1, s2] = s0s1s2;
                  auto Hessian = grad_U_LinearElement(f, s0s1s2);
                  phin_t -= std::get<0>(Dot(Tddd{{Dot(U, s0), Dot(U, s1), Dot(U, s2)}}, Hessian));
               } else {
                  auto netInContact = NearestContactFace(p)->getNetwork();
                  auto w = netInContact->velocityRotational();
                  auto dQdt = Q.d_dt(w);
                  auto U = uNeumann(p);
                  phin_t = std::get<1>(p->phiphin_t) = Dot(w, U - p->U_BEM) + Dot(n, accelNeumann(p));
                  auto s0s1s2 = OrthogonalBasis(p->getNormal_BEM());
                  auto [s0, s1, s2] = s0s1s2;
                  auto Hessian = grad_U_LinearElement(p, s0s1s2);
                  phin_t -= std::get<0>(Dot(Tddd{{Dot(U, s0), Dot(U, s1), Dot(U, s2)}}, Hessian));
               }
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
   V_d Func(V_d ACCELS_IN, const Network *water, const std::vector<Network *> &rigidbodies) const {
      auto ACCELS = ACCELS_IN;
      {
         int i = 0;
         for (const auto &net : rigidbodies)
            if (isTarget(net))
               for_each(net->acceleration, [&](auto &a_w) { a_w = ACCELS_IN[i++]; });
      }

      //* --------------------------------------------------- */
      //*                  加速度 --> phiphin_t                */
      //* --------------------------------------------------- */
      setPhiPhin_t();
      V_d knowns(PBF_index.size());
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, DorN, f] = PBF;
         if (DorN == Dirichlet)
            knowns[i] = std::get<0>(p->phiphin_t);
         else
            knowns[i] = p->phintOnFace.at(f);  // はいってない？はいってた．
      }

      V_d phiORphin_t(PBF_index.size());
      std::cout << "加速度 --> phiphin_t" << std::endl;

#if defined(use_CG)
      GradientMethod gd(mat_ukn);
      phiORphin_t = gd.solve(Dot(mat_kn, knowns));
#elif defined(use_gmres)
      std::cout << "gmres for phiphin_t" << std::endl;
      gmres gm(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin_t /*解*/, use_gmres);
      phiORphin_t = gm.x;
      if (!isFinite(gm.err)) {
         std::cout << "gm.err = " << gm.err << std::endl;
         this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin_t /*解*/);
      }
      // auto err = Norm(Dot(mat_ukn, gm.x) - Dot(mat_kn, knowns));
      // gmres gm(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin_t /*解*/, 5);
      // auto err = Norm(Dot(mat_ukn, gm.x) - Dot(mat_kn, knowns));
      // std::cout << err << std::endl;
      // /* ------------------------------------------------------ */
      // if (isFinite(err) && err < 0.1)
      // else
      //
#elif defined(use_lapack)
      this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin_t /*解*/);
#endif
      std::cout << "solved" << std::endl;
      //* --------------------------------------------------- */
      //*                 phiphin_t --> 圧力                   */
      //* --------------------------------------------------- */
      // for (const auto &[PBF, i] : PBF_index) {
      //    auto [p, DorN, f] = PBF;
      //    if (DorN == Dirichlet)
      //       std::get<1>(p->phiphin_t) = phiORphin_t[i];
      //    else
      //       std::get<0>(p->phiphin_t) = phiORphin_t[i];
      //    p->pressure = p->pressure_BEM = -_WATER_DENSITY_ * (std::get<0>(p->phiphin_t) + _GRAVITY_ * p->height() + Dot(p->U_BEM, p->U_BEM) / 2.);
      // }
      /* -------------------------------------------------------------------------- */
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, DorN, f] = PBF;
         if (DorN == Dirichlet) {
            // do nothing
         } else if (p->phintOnFace.size() > 1) {
            std::get<0>(p->phiphin_t) = 0;
            std::get<1>(p->phiphin_t) = 0;
         } else
            std::get<0>(p->phiphin_t) = 0;
      }

      for (const auto &[PBF, i] : PBF_index) {
         auto [p, DorN, f] = PBF;
         if (DorN == Dirichlet)
            std::get<1>(p->phiphin_t) = phiORphin_t[i];
         else if (p->phintOnFace.size() > 1) {
            double total = 0;
            for (auto [f, phin_t] : p->phintOnFace)
               total += f->area;
            std::get<0>(p->phiphin_t) += phiORphin_t[i] * f->area / total;
            std::get<1>(p->phiphin_t) += p->phintOnFace.at(f) * f->area / total;
         } else {
            std::get<0>(p->phiphin_t) = phiORphin_t[i];
            std::get<1>(p->phiphin_t) = p->phintOnFace.at(f);
         }
      }
      /* -------------------------------------------------------------------------- */
      // // 足し合わせるので初期化
      // for (const auto &[PBF, i] : PBF_index) {
      //    auto [p, DorN, f] = PBF;
      //    if (DorN == Neumann)
      //       std::get<0>(p->phiphin_t) = 0;
      // }
      // //
      // for (const auto &[PBF, i] : PBF_index) {
      //    auto [p, DorN, f] = PBF;
      //    if (DorN == Dirichlet)
      //       std::get<1>(p->phiphin_t) = phiORphin_t[i];
      //    else if (p->phinOnFace.size() > 1) {
      //       double total = 0;
      //       for (auto [f, phin] : p->phinOnFace)
      //          total += f->area;

      //       std::get<0>(p->phiphin_t) += phiORphin_t[i] * f->area / total;
      //    } else
      //       std::get<0>(p->phiphin_t) = phiORphin_t[i];
      // }

      for (const auto &[PBF, i] : PBF_index) {
         auto [p, DorN, f] = PBF;
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
            for_each(T6d{a0, a1, a2, a3, a4, a5}, [&](const auto &a_w) { ACCELS[i++] = a_w; });
         }
      return ACCELS - ACCELS_IN;
   };

   /* ------------------------------------------------------ */

   void solveForPhiPhin_t(const Network *water, const std::vector<Network *> &rigidbodies) const {
      //@ --------------------------------------------------- */
      //@        加速度 --> phiphin_t --> 圧力 --> 加速度        */
      //@ --------------------------------------------------- */
      V_d ACCELS_init, ACCELS, ACCELS_old, ACCELS_old_old;
      for (const auto &net : rigidbodies) {
         if (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating") {
            if (net->inputJSON.at("velocity").size() > 1) {
               std::cout << Red << "net->inputJSON[\" velocity \"][1] = " << net->inputJSON.at("velocity")[1] << colorOff << std::endl;
               double start_time = std::stod(net->inputJSON.at("velocity")[1]);
               std::cout << Red << "start_time = " << start_time << colorOff << std::endl;
               if (real_time < start_time)
                  net->acceleration.fill(0.);
               else
                  for_each(net->acceleration, [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
            } else
               for_each(net->acceleration, [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
         } else
            net->acceleration.fill(0.);
      }

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
         // for (const auto &net : rigidbodies)
         //    if (isTarget(net))
         //       for_each(net->acceleration, [&](auto &a_w) { a_w = BM.X[i++]; });
         // save --> acceleration
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
                     for_each(net->acceleration, [&](auto &a_w) { a_w = BM.X[i++]; });
               } else
                  for_each(net->acceleration, [&](auto &a_w) { a_w = BM.X[i++]; });
            } else
               net->acceleration.fill(0.);
         }

         std::cout << "j = " << j << ", " << Red << Norm(func) << colorOff << std::endl;
         std::cout << Red << "func_ = " << func_ << colorOff << std::endl;
         std::cout << Red << "func = " << func << colorOff << std::endl;
         std::cout << Red << "BM.X = " << BM.X << colorOff << std::endl;
         std::cout << Red << "BM.dX = " << BM.dX << colorOff << std::endl;

         if (Norm(func) < 1E-10)
            break;
      }
   };

   //% -------------------------------------------------------------------------- */
   //%                             solve phi and phi_n                            */
   //% -------------------------------------------------------------------------- */

   void solve(const Network &water, const Buckets<networkPoint *> &FMM_BucketsPoints, const Buckets<networkFace *> &FMM_BucketsFaces) {
      //* ------------------------------------------------------ */
      //%                     各点で方程式を作る場合                 */
      //* ------------------------------------------------------ */
      std::cout << "各点で方程式を作る場合" << std::endl;
      PBF_index.clear();
      PBF_index.reserve(3 * water.getPoints().size());
      knowns.clear();
      knowns.reserve(3 * water.getPoints().size());
      int i = 0;
      for (const auto &p : water.getPoints()) {
         if (p->Dirichlet || p->CORNER) {
            PBF_index[{p, Dirichlet, nullptr}] = i++;
            knowns.emplace_back(p->phi_Dirichlet = std::get<0>(p->phiphin));
         }
         //! PBF_indexのNeuamnn箇所に関しては，設定済みの　p->phinOnFace　の状態にに任せる
         for (const auto &[f, phin] : p->phinOnFace) {
            PBF_index[{p, Neumann, f}] = i++;
            knowns.emplace_back(phin);
         }
      }

      std::cout << Red << "water.getPoints() = " << i << std::endl;
      IGIGn = std::vector<std::vector<std::array<double, 2>>>(PBF_index.size(), std::vector<std::array<double, 2>>(PBF_index.size(), {0., 0.}));
      mat_kn = mat_ukn = VV_d(PBF_index.size(), V_d(PBF_index.size(), 0.));

      auto PF2index = [&](netP *p, netF *integ_f) {
         auto it = PBF_index.find({p, integ_f->Dirichlet, integ_f});
         if (it != PBF_index.end())
            return it->second;
         else
            return PBF_index[{p, integ_f->Dirichlet, nullptr}];
      };

#define use_rigid_mode
      Timer timer;
      std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;
      std::cout << Magenta << timer() << colorOff << std::endl;
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (const auto &[PBF, index] : PBF_index)
#ifdef _OPENMP
#pragma omp single nowait
#endif
      {
         auto [origin, _, __] = PBF;
         auto it = PBF_index.begin();
         double origin_ign_rigid_mode = 0.;
         auto &IGIGn_Row = IGIGn[index];
         double nr, tmp;
         std::array<double, 2> IGIGn, c;
         std::array<double, 3> X0, X1, X2, A, cross;
         std::array<double, 3> N012;
         for (const auto &integ_f : FMM_BucketsFaces.all_stored_objects) {

            // const auto [p0, p1, p2] = integ_f->getPoints(origin);
            // ret = {{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}};
            // X2 = p2->getXtuple();
            // X0 = p0->getXtuple() - X2;
            // X1 = p1->getXtuple() - X2;
            // A = origin->getXtuple() - X2;
            // cross = Cross(X0, X1);
            // c = {Norm(cross), Dot(A, cross)};

            //
            // for_each(__GW5xGW5__, [&](const auto &GWGW) {
            //    tmp = std::get<2>(GWGW) * (1. - std::get<0>(GWGW)) / (nr = Norm(X0 * std::get<0>(GWGW) + X1 * std::get<1>(GWGW) * (1. - std::get<0>(GWGW)) - A));
            //    std::get<1>(std::get<0>(ret)) += (IGIGn = {tmp, tmp / (nr * nr)}) * std::get<0>(GWGW);
            //    std::get<1>(std::get<1>(ret)) += IGIGn * std::get<1>(GWGW) * (1. - std::get<0>(GWGW));
            //    std::get<1>(std::get<2>(ret)) += IGIGn * (-1. + std::get<0>(GWGW)) * (-1. + std::get<1>(GWGW));
            // });
            //
            /* -------------------------------------------------------------------------- */
            //
            // const auto [p0, p1, p2] = integ_f->getPoints(origin);
            // ret = {{{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}}};
            // X2 = p2->getXtuple();
            // X0 = p0->getXtuple() - X2;
            // X1 = p1->getXtuple() - X2;
            // A = origin->getXtuple() - X2;
            // cross = Cross(X0, X1);
            // c = {Norm(cross), Dot(A, cross)};

            // for_each(__array_GW5xGW5__, [&](const auto &GWGW) {
            //    tmp = std::get<2>(GWGW) * (1. - std::get<0>(GWGW)) / (nr = Norm(X0 * std::get<0>(GWGW) + X1 * std::get<1>(GWGW) * (1. - std::get<0>(GWGW)) - A));
            //    std::get<1>(std::get<0>(ret)) += (IGIGn = {tmp, tmp / (nr * nr)}) * std::get<0>(GWGW);           // 補間添字0
            //    std::get<1>(std::get<1>(ret)) += IGIGn * std::get<1>(GWGW) * (1. - std::get<0>(GWGW));           // 補間添字1
            //    std::get<1>(std::get<2>(ret)) += IGIGn * (-1. + std::get<0>(GWGW)) * (-1. + std::get<1>(GWGW));  // 補間添字2
            // });
            /* -------------------------------------------------------------------------- */
            // // 2023/04/03
            // const auto [p0, p1, p2] = integ_f->getPoints(origin);
            // ret = {{{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}}};
            // X0 = p0->getXtuple();
            // X1 = p1->getXtuple();
            // X2 = p2->getXtuple();
            // std::array<std::array<double, 3>, 3> X012{{ToArray(X0), ToArray(X1), ToArray(X2)}};
            // A = origin->getXtuple();
            // cross = Cross(X0 - X2, X1 - X2);
            // c = {Norm(cross), Dot(A - (X0 + X1 + X2) / 3., cross)};
            // std::array<double, 3> N012;

            // auto add = [&integ_f](const double &t0, const double &t1) {
            //    N012 = integ_f->ModTriShape6(t0, t1);
            // };

            // for_each(__array_GW5xGW5__, [&](const auto &GWGW) {
            //    N012 = ModTriShape<3>(GWGW[0], GWGW[1]);
            //    tmp = GWGW[2] * (1. - std::get<0>(GWGW)) / (nr = Norm(Dot(N012, X012) - ToArray(A)));
            //    IGIGn = {tmp, tmp / (nr * nr)};
            //    //
            //    // std::get<1>(std::get<0>(ret)) += IGIGn * N012[0];  // 補間添字0
            //    // std::get<1>(std::get<1>(ret)) += IGIGn * N012[1];  // 補間添字1
            //    // std::get<1>(std::get<2>(ret)) += IGIGn * N012[2];  // 補間添字2
            //    //
            //    add();
            // });
            // /* -------------------------------------------------------------------------- */
            // 2023/04/03
            // const auto [p0, p1, p2] = integ_f->getPoints(origin);
            // X0 = p0->getXtuple();
            // X1 = p1->getXtuple();
            // X2 = p2->getXtuple();
            // A = origin->getXtuple();
            // cross = Cross(X0 - X2, X1 - X2);
            // c = {Norm(cross), Dot(A - (X0 + X1 + X2) / 3., cross)};
            // std::array<std::array<double, 3>, 3> X012{{ToArray(X0), ToArray(X1), ToArray(X2)}};
            // std::array<double, 3> N012;
            // //
            // // std::array<std::tuple<netP *, Tdd>, 3> ret{{{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}}};
            // //
            // ModTriShape6 shape6(integ_f, origin);
            // std::vector<std::tuple<netP *, Tdd>> ret;
            // for (auto &[p, i] : shape6.p2i)
            //    ret.push_back({p, {0., 0.}});
            // for_each(__array_GW5xGW5__, [&](const auto &GWGW) {
            //    N012 = ModTriShape<3>(GWGW[0], GWGW[1]);
            //    tmp = GWGW[2] * (1. - std::get<0>(GWGW)) / (nr = Norm(Dot(N012, X012) - ToArray(A)));
            //    IGIGn = {tmp, tmp / (nr * nr)};
            //    //
            //    auto shape12 = shape6(GWGW[0], GWGW[1]);
            //    for (auto i = 0; i < ret.size(); ++i)
            //       std::get<1>(ret[i]) *= IGIGn * shape12[i];  // 補間添字
            // });
            // /* -------------------------------------------------------------------------- */
            // // 2023/04/03
            // // const auto [p0, p1, p2] = integ_f->getPoints(origin);
            // const auto [p0, l0, p1, l1, p2, l2] = integ_f->getPointsAndLines(origin);
            // //
            // auto f01 = (*l0)(integ_f);
            // auto p01 = f01->getPointOpposite(l0);
            // std::array<double, 2> p01_w01{{1., 1.}}, p01_w0{{1., 1.}}, p01_w1{{1., 1.}}, p01_w2{{1., 1.}};
            // double n = 2;
            // p01_w01[0] = l0->CORNER ? 0. : 1.;
            // p01_w01 *= 1 / std::pow(Norm(p01->X - 0.5 * (p0->X + p1->X)), n);
            // p01_w0 *= 1 / std::pow(Norm(p0->X - 0.5 * (p0->X + p1->X)), n);
            // p01_w1 *= 1 / std::pow(Norm(p1->X - 0.5 * (p0->X + p1->X)), n);
            // p01_w2 *= 1 / std::pow(Norm(p2->X - 0.5 * (p0->X + p1->X)), n);
            // std::array<double, 2> sum0 = p01_w01 + p01_w0 + p01_w1 + p01_w2;
            // p01_w01 /= sum0;
            // p01_w0 /= sum0;
            // p01_w1 /= sum0;
            // p01_w2 /= sum0;
            // //
            // auto f12 = (*l1)(integ_f);
            // auto p12 = f12->getPointOpposite(l1);
            // std::array<double, 2> p12_w12{{1., 1.}}, p12_w1{{1., 1.}}, p12_w2{{1., 1.}}, p12_w0{{1., 1.}};
            // p12_w12[0] = l1->CORNER ? 0. : 1.;
            // p12_w12 *= 1 / std::pow(Norm(p12->X - 0.5 * (p1->X + p2->X)), n);
            // p12_w1 *= 1 / std::pow(Norm(p1->X - 0.5 * (p1->X + p2->X)), n);
            // p12_w2 *= 1 / std::pow(Norm(p2->X - 0.5 * (p1->X + p2->X)), n);
            // p12_w0 *= 1 / std::pow(Norm(p0->X - 0.5 * (p1->X + p2->X)), n);
            // auto sum1 = p12_w12 + p12_w1 + p12_w2 + p12_w0;
            // p12_w12 /= sum1;
            // p12_w1 /= sum1;
            // p12_w2 /= sum1;
            // p12_w0 /= sum1;
            // //
            // auto f20 = (*l2)(integ_f);
            // auto p20 = f20->getPointOpposite(l2);
            // std::array<double, 2> p20_w20{{1., 1.}}, p20_w2{{1., 1.}}, p20_w0{{1., 1.}}, p20_w1{{1., 1.}};
            // p20_w20[0] = l2->CORNER ? 0. : 1.;
            // p20_w20 *= 1 / std::pow(Norm(p20->X - 0.5 * (p2->X + p0->X)), n);
            // p20_w2 *= 1 / std::pow(Norm(p2->X - 0.5 * (p2->X + p0->X)), n);
            // p20_w0 *= 1 / std::pow(Norm(p0->X - 0.5 * (p2->X + p0->X)), n);
            // p20_w1 *= 1 / std::pow(Norm(p1->X - 0.5 * (p2->X + p0->X)), n);
            // std::array<double, 2> sum2 = p20_w20 + p20_w2 + p20_w0 + p20_w1;
            // p20_w20 /= sum2;
            // p20_w2 /= sum2;
            // p20_w0 /= sum2;
            // p20_w1 /= sum2;
            // //
            // std::array<std::tuple<networkPoint *, networkFace *, std::array<double, 2>>, 6> ret = {{{p0, integ_f, {0., 0.}},
            //                                                                                         {p1, integ_f, {0., 0.}},
            //                                                                                         {p2, integ_f, {0., 0.}},
            //                                                                                         {p01, f01, {0., 0.}},
            //                                                                                         {p12, f12, {0., 0.}},
            //                                                                                         {p20, f20, {0., 0.}}}};
            // //
            // X0 = p0->getXtuple();
            // X1 = p1->getXtuple();
            // X2 = p2->getXtuple();
            // std::array<std::array<double, 3>, 3> X012{{ToArray(X0), ToArray(X1), ToArray(X2)}};
            // A = origin->getXtuple();
            // cross = Cross(X0 - X2, X1 - X2);
            // c = {Norm(cross), Dot(A - (X0 + X1 + X2) / 3., cross)};
            // std::array<double, 3> N012;
            // for (const auto &[t0, t1, ww] : __array_GW5xGW5__) {
            //    N012 = ModTriShape<3>(t0, t1);
            //    tmp = ww * (1. - t0) / (nr = Norm(N012[0] * X0 + N012[1] * X1 + N012[2] * X2 - A));
            //    IGIGn = {tmp, tmp / (nr * nr)};
            //    //
            //    auto N012345 = ModTriShape<6>(t0, t1);
            //    std::get<2>(ret[0]) += IGIGn * N012345[0];  // 補間添字
            //    std::get<2>(ret[1]) += IGIGn * N012345[1];  // 補間添字
            //    std::get<2>(ret[2]) += IGIGn * N012345[2];  // 補間添字
            //    //
            //    std::get<2>(ret[3]) += IGIGn * p01_w01 * N012345[3];  // 補間添字
            //    std::get<2>(ret[0]) += IGIGn * p01_w0 * N012345[3];   // 補間添字
            //    std::get<2>(ret[1]) += IGIGn * p01_w1 * N012345[3];   // 補間添字
            //    std::get<2>(ret[2]) += IGIGn * p01_w2 * N012345[3];   // 補間添字
            //    //
            //    std::get<2>(ret[4]) += IGIGn * p12_w12 * N012345[4];  // 補間添字
            //    std::get<2>(ret[1]) += IGIGn * p12_w1 * N012345[4];   // 補間添字
            //    std::get<2>(ret[2]) += IGIGn * p12_w2 * N012345[4];   // 補間添字
            //    std::get<2>(ret[0]) += IGIGn * p12_w0 * N012345[4];   // 補間添字
            //    //
            //    std::get<2>(ret[5]) += IGIGn * p20_w20 * N012345[5];  // 補間添字
            //    std::get<2>(ret[2]) += IGIGn * p20_w2 * N012345[5];   // 補間添字
            //    std::get<2>(ret[0]) += IGIGn * p20_w0 * N012345[5];   // 補間添字
            //    std::get<2>(ret[1]) += IGIGn * p20_w1 * N012345[5];   // 補間添字
            // }
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
               このループでは，
               ある面integ_fに隣接する節点{p0,p1,p2}の列,IGIGn[origin(fixed),p0],...に値が追加されていく．
               （p0が多重接点の場合，適切にp0と同じ位置に別の変数が設定されており，別の面の積分の際にq0が参照される．）
               //
               p0は，{面,補間添字}で決定することもできる．
               {面,補間添字0}->p0,{面,補間添字1}->p1,{面,補間添字2}->p2というように．
               //
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
               面を区別するかどうかが先にわからないので，face*のままかnullptrとすべきかわからないということ．．．．
               //
               PBF_index[{p, Dirichlet, ある要素}]
               は存在しないだろう．Dirichlet節点は，{p, ある要素}からの寄与を，ある面に
               */
               //
               // auto [p, igign] = p_igign;
               IGIGn_Row[PF2index(p, which_side_f)] += igign;   // この面に関する積分において，φまたはφnの寄与
               if (p != origin)                                 // for use_rigid_mode
                  origin_ign_rigid_mode -= std::get<1>(igign);  // for use_rigid_mode
            }
         }
         /* -------------------------------------------------------------------------- */
         /**
          * # Example Function
          *
          * This is an example function that demonstrates how to use the keywords.
          *
          * NOTE: This is a note.
          * WARNING: This is a warning.
          * TODO: This is a todo item.
          * IMPORTANT: This is an important point.
          * TIP: This is a helpful tip.
          */
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
      //
      std::cout << Green << "離散化にかかった時間" << timer() << colorOff << std::endl;
      /* ------------------------------------------------------ */
      // #define quad_element
#define linear_element
      // #define liear_and_quad_element

      std::cout << "並列化 DONE" << std::endl;
      std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;

      {
         // b@ ------------------------------------------------------ */
         // b@                 系数行列mat_ukn．mat_knの計算             */
         // b@ ------------------------------------------------------ */
#pragma omp parallel
         for (const auto &[PBF, i] : PBF_index)
#pragma omp single nowait
         {
            /*@ mat_ukn, mat_kn, vec_P
              +--+--+--+
            | |--+--+--|
            V |--+--+--|
              +--+--+--+
            */
            auto [p, DorN, f] = PBF;
            // auto &IGIGn_Row = IGIGn[i];
            // p->IGIGn[vec_P[i]]; // PBF_PBF_IGIGn.at(vec_P[i]);
            /* mat_ukn, mat_kn, vec_P
                ====>
              +--+--+--+
              |--+--+--|
              |--+--+--|
              +--+--+--+
            */
            // if (!(p->CORNER && DorN == Neumann /*変更する対象の行*/))  //! OK
            {
               std::array<double, 2> igign;
               for (const auto &[PBF_j, j] : PBF_index) {
                  igign = IGIGn[i][j];
                  // 未知変数の係数行列は左，既知変数の係数行列は右
                  if (std::get<1>(PBF_j) == Neumann)
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
            }
         }

         /* ------------------------------------------------------ */
         double maxpp = 0;
         for (auto i = 0; i < mat_ukn.size(); ++i) {
            // LUするのはmat_uknだけなので，mat_knの最大値を使う必要はない
            if (maxpp < std::abs(mat_ukn[i][i]))
               maxpp = std::abs(mat_ukn[i][i]);
         }
         /* ------------------------------------------------------ */
         // knowns = V_d(PBF_index.size());  //! ok
         //                                  //! Pointsの順番と合わせてとるように注意

#ifdef _OPENMP
#pragma omp parallel
#endif
         for (const auto &[PBF, i] : PBF_index)
#ifdef _OPENMP
#pragma omp single nowait
#endif
         {
            auto [p, DorN, f] = PBF;
            // auto &PBF_IGIGn = IGIGn[i];
            // ディリクレの角点だけが積分した係数を使った方程式を使う．
            // if (DorN == Neumann && p->CORNER /*多重節点の条件*/) {

            /*
               移項前:{IG0,IG1,IG2,IG3} . {φn0,φn1,φn2,φn3} = {IGn0,IGn1,IGn2,IGn3} . {φ0,φ1,φ2,φ3}
               移項後:{IG0,-IGn1,IG2,IG3} . {φn0,φ1,φn2,φn3} = {IGn0,-IG1,IGn2,IGn3} . {φ0,φn1,φ2,φ3}
               多重節点:{0, 1, 0, 0} . {φn0,φ1,φn2,φn3} = {0, 0, 0, 1} . {φ0,φn1,φ2,φ3}
            */

            // b$ ------------------------------------------------------ */
            // b$               mat_unknowns, mat_knownsの計算            */
            // b$ ------------------------------------------------------ */

            if (DorN == Neumann && p->CORNER /*多重節点の条件*/) {
               for (const auto &[PBF_j, j] : PBF_index) {
                  auto [p_, DorN_, f_] = PBF_j;
                  if (p == p_) {
                     if (DorN_ == Neumann && f_ == f /*can be nullptr*/) {
                        mat_ukn[i][j] = maxpp;  // φの系数
                        mat_kn[i][j] = 0;       // φnの系数
                     } else if (DorN_ == Dirichlet && f_ == nullptr /* there is only one in this row*/) {
                        mat_ukn[i][j] = 0;     // φnの系数
                        mat_kn[i][j] = maxpp;  // φの系数移行したからマイナス？　いいえ，移項を考慮した上でこれでいい．
                     } else {
                        mat_ukn[i][j] = 0;
                        mat_kn[i][j] = 0;
                     }
                  } else {
                     mat_ukn[i][j] = 0;
                     mat_kn[i][j] = 0;
                  }
               }
            }
         }

         if (!isFinite(mat_ukn))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mat_ukn is not finite");
         if (!isFinite(mat_kn))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mat_kn is not finite");
      }

      {
         // b% ------------------------------------------------------ */
         // b%                  境界積分方程式を解く                      */
         // b% ------------------------------------------------------ */
         std::cout << "--------------------- 境界積分方程式を解く ---------------------" << std::endl;
         TimeWatch watch;
         /* ------------------------------------------------------ */
         V_d phiORphin(knowns.size(), 0);
         std::cout << "IGIGn.size()= " << IGIGn.size() << std::endl;
         std::cout << "phiORphin.size()= " << phiORphin.size() << std::endl;
         std::cout << "  mat_kn.size() = " << mat_kn.size() << std::endl;
         //* 未知変数の計算
         /* ------------------------------------------------------ */
         // 前処理
         // auto v = diagonal_scaling_vector(mat_ukn);
         // for (auto i = 0; i < v.size(); ++i) {
         //    mat_ukn[i] *= v[i];
         //    mat_kn[i] *= v[i];
         // }
         /* ------------------------------------------------------ */
         if (this->lu)
            delete this->lu;
#if defined(use_CG)
         this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/);
         std::cout << "The conjugate gradient is used" << std::endl;
         GradientMethod gd1(mat_ukn);
         phiORphin = gd1.solve(Dot(mat_kn, knowns), {}, 1E-1);
         GradientMethod gd2(mat_ukn);
         phiORphin = gd2.solveCG(Dot(mat_kn, knowns), phiORphin);
#elif defined(use_gmres)
         std::cout << "gmres for phiORphin" << std::endl;
         gmres gm(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/, use_gmres);
         phiORphin = gm.x;
         std::cout << "gm.err = " << gm.err << ", isFinite(gm.err) = " << isFinite(gm.err) << std::endl;
         if (real_time < 0.005 || !isFinite(gm.err < 1E-20)) {
            std::cout << "lapack lu decomposition" << std::endl;
            this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
         }
         // auto err = Norm(Dot(mat_ukn, gm.x) - Dot(mat_kn, knowns));
         // std::cout << err << std::endl;
#elif defined(use_lapack)
         std::cout << "lapack lu decomposition" << std::endl;
         // this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/);
         // std::cout << "try to solve" << std::endl;
         // this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
         this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/, Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#elif defined(use_lu)
         std::cout << "parallel lu decomposition" << std::endl;
         this->lu = new ludcmp_parallel(mat_ukn /*未知の行列係数（左辺）*/);
         std::cout << "try to solve" << std::endl;
         this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#else
         //* 未知変数の計算
         std::cout << "SVD decomposition" << std::endl;
         SVD svd(mat_ukn);
         svd.solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#endif
         std::cout << Blue << "Elapsed time for solving BIE: " << Red << watch() << colorOff << " s\n";
         // 足し合わせるので初期化
         for (const auto &[PBF, i] : PBF_index) {
            auto [p, DorN, f] = PBF;
            if (DorN == Dirichlet) {
               // do nothing
            } else if (p->phinOnFace.size() > 1) {
               std::get<0>(p->phiphin) = 0;
               std::get<1>(p->phiphin) = 0;
            } else
               std::get<0>(p->phiphin) = 0;
         }

         for (const auto &[PBF, i] : PBF_index) {
            auto [p, DorN, f] = PBF;
            if (DorN == Dirichlet)
               std::get<1>(p->phiphin) = p->phin_Dirichlet = phiORphin[i];
            else if (p->phinOnFace.size() > 1) {
               double total = 0;
               for (auto [f, phin] : p->phinOnFace)
                  total += f->area;
               std::get<0>(p->phiphin) += phiORphin[i] * f->area / total;
               std::get<1>(p->phiphin) += p->phinOnFace.at(f) * f->area / total;
            } else {
               std::get<0>(p->phiphin) = phiORphin[i];
               std::get<1>(p->phiphin) = p->phinOnFace.at(f);
            }
         }

         if (!isFinite(phiORphin))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "phiORphin is not finite");
      }
   };
};

#endif