#ifndef BEM_BVP_H
#define BEM_BVP_H

#include "Network.hpp"
// #define use_CG
// #define use_gmres 20
#define use_lapack
std::unordered_map<std::tuple<netP *, bool, netF *>, int> PBF_index;

struct calculateFroudeKrylovForce {
   std::vector<networkFace *> actingFaces;
   Tddd force, torque;
   double area;
   T6d acceleration;
   std::vector<std::tuple<Tddd, T3Tddd>> PressureVeticies;
   calculateFroudeKrylovForce(const std::unordered_set<networkFace *> faces /*waterfaces*/,
                              const Network *PasObj)
       : force({0., 0., 0.}),
         torque({0., 0., 0.}),
         area(0.),
         PressureVeticies({}),
         acceleration({0., 0., 0., 0., 0., 0.}) {
      // PasObjと接したfaceの頂点にpressureが設定されている前提
      int count = 0;
      for (const auto &f : faces)
         if (f->Neumann) {
            auto [p0, p1, p2] = f->getPoints();
            if (std::any_of(p0->getContactFaces().begin(), p0->getContactFaces().end(), [&](const auto &F) { return F->getNetwork() == PasObj; }) &&
                std::any_of(p1->getContactFaces().begin(), p1->getContactFaces().end(), [&](const auto &F) { return F->getNetwork() == PasObj; }) &&
                std::any_of(p2->getContactFaces().begin(), p2->getContactFaces().end(), [&](const auto &F) { return F->getNetwork() == PasObj; })) {
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
      this->force = {0., 0., 0.};
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpP = interpolationTriangleLinear0101(P012);
         auto intpX = interpolationTriangleLinear0101(X012);
         auto n = TriangleNormal(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            this->force += n * intpP(x0, x1) * intpX.J(x0, x1) * w0w1;
      }
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
   std::vector<std::vector<Tdd>> IGIGn;
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
                  // phin_t = std::get<1>(p->phiphin_t) = Dot(n, Dot(uNeumann(p, f) - grad_phi, dQdt.Rv()) + accelNeumann(p, f) - Dot(grad_phi, grad_U_LinearElement(f)));
                  // phin_t = std::get<1>(p->phiphin_t) = Dot(w, U - p->U_BEM) + Dot(n, accelNeumann(p, f) - Dot(U, Hessian));
                  // phin_t = std::get<1>(p->phiphin_t) = Dot(w, U - p->U_BEM) + Dot(n, accelNeumann(p, f) - Dot(U, Hessian));
                  //
                  phin_t = std::get<1>(p->phiphin_t) = Dot(w, U - p->U_BEM) + Dot(n, accelNeumann(p, f));
                  //
                  auto s0s1s2 = OrthogonalBasis(f->normal);
                  auto [s0, s1, s2] = s0s1s2;
                  auto Hessian = grad_U_LinearElement(f, s0s1s2);
                  phin_t -= std::get<0>(Dot({Dot(U, s0), Dot(U, s1), Dot(U, s2)}, Hessian));
                  //
                  // phin_t = std::get<1>(p->phiphin_t) = Dot(w, uNeumann(p, f) - p->U_BEM) + Dot(n, accelNeumann(p, f) - Dot(p->U_BEM, grad_U_LinearElement(f)));
               } else {
                  auto netInContact = NearestContactFace(p)->getNetwork();
                  auto w = netInContact->velocityRotational();
                  auto dQdt = Q.d_dt(w);
                  // phin_t = std::get<1>(p->phiphin_t) = Dot(n, Dot(uNeumann(p) - p->U_BEM, dQdt.Rv()) + accelNeumann(p) - Dot(p->U_BEM, grad_U_LinearElement(p)));
                  auto U = uNeumann(p);
                  // phin_t = std::get<1>(p->phiphin_t) = Dot(w, U - p->U_BEM) + Dot(n, accelNeumann(p) - Dot(U, Hessian));
                  //
                  phin_t = std::get<1>(p->phiphin_t) = Dot(w, U - p->U_BEM) + Dot(n, accelNeumann(p));
                  //
                  auto s0s1s2 = OrthogonalBasis(p->getNormal_BEM());
                  auto [s0, s1, s2] = s0s1s2;
                  auto Hessian = grad_U_LinearElement(p, s0s1s2);
                  phin_t -= std::get<0>(Dot({Dot(U, s0), Dot(U, s1), Dot(U, s2)}, Hessian));
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

      // read --> acceleration
      for (const auto &net : rigidbodies) {
         if (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating") {
            if (net->inputJSON.at("velocity").size() > 1) {
               std::cout << Red << "net->inputJSON[\" velocity \"][1] = " << net->inputJSON.at("velocity")[1] << colorOff << std::endl;
               double start_time = std::stod(net->inputJSON.at("velocity")[1]);
               std::cout << Red << "start_time = " << start_time << colorOff << std::endl;
               if (real_time < start_time) {
                  for_each(net->acceleration, [&](auto &a_w) { a_w = 0.; });
               } else
                  for_each(net->acceleration, [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
            } else
               for_each(net->acceleration, [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
         } else
            for_each(net->acceleration, [&](auto &a_w) { a_w = 0.; });
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
                     for_each(net->acceleration, [&](auto &a_w) { a_w = 0.; });
                  } else
                     for_each(net->acceleration, [&](auto &a_w) { a_w = BM.X[i++]; });
               } else
                  for_each(net->acceleration, [&](auto &a_w) { a_w = BM.X[i++]; });
            } else
               for_each(net->acceleration, [&](auto &a_w) { a_w = 0.; });
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
      int i = 0;
      for (const auto &p : water.getPoints()) {
         if (p->Dirichlet || p->CORNER)
            PBF_index[{p, Dirichlet, nullptr}] = i++;
         //! PBF_indexのNeuamnn箇所に関しては，設定済みの　p->phinOnFace　の状態にに任せる
         for (const auto &[f, _] : p->phinOnFace)
            PBF_index[{p, Neumann, f}] = i++;
      }

      std::cout << Red << "water.getPoints() = " << i << std::endl;
      IGIGn = std::vector<std::vector<Tdd>>(PBF_index.size(), std::vector<Tdd>(PBF_index.size(), {0., 0.}));
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
         std::array<std::tuple<netP *, Tdd>, 3> ret;
         Tdd IGIGn, c;
         double nr, tmp;
         Tddd X2, X0, X1, A, cross;
         for (const auto &integ_f : FMM_BucketsFaces.all_stored_objects) {
            {
               const auto [p0, p1, p2] = integ_f->getPoints(origin);
               // ret = {{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}};
               ret = {{{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}}};
               X2 = p2->getXtuple();
               X0 = p0->getXtuple() - X2;
               X1 = p1->getXtuple() - X2;
               A = origin->getXtuple() - X2;
               cross = Cross(X0, X1);
               c = {Norm(cross), Dot(A, cross)};
               if (origin == p0 || origin == p1 || origin == p2)
                  std::get<1>(c) = 0;

               // for_each(__GW5xGW5__, [&](const auto &GWGW) {
               //    tmp = std::get<2>(GWGW) * (1. - std::get<0>(GWGW)) / (nr = Norm(X0 * std::get<0>(GWGW) + X1 * std::get<1>(GWGW) * (1. - std::get<0>(GWGW)) - A));
               //    std::get<1>(std::get<0>(ret)) += (IGIGn = {tmp, tmp / (nr * nr)}) * std::get<0>(GWGW);
               //    std::get<1>(std::get<1>(ret)) += IGIGn * std::get<1>(GWGW) * (1. - std::get<0>(GWGW));
               //    std::get<1>(std::get<2>(ret)) += IGIGn * (-1. + std::get<0>(GWGW)) * (-1. + std::get<1>(GWGW));
               // });
               //

               for_each(__array_GW5xGW5__, [&](const auto &GWGW) {
                  tmp = std::get<2>(GWGW) * (1. - std::get<0>(GWGW)) / (nr = Norm(X0 * std::get<0>(GWGW) + X1 * std::get<1>(GWGW) * (1. - std::get<0>(GWGW)) - A));
                  std::get<1>(std::get<0>(ret)) += (IGIGn = {tmp, tmp / (nr * nr)}) * std::get<0>(GWGW);           // 補間添字0
                  std::get<1>(std::get<1>(ret)) += IGIGn * std::get<1>(GWGW) * (1. - std::get<0>(GWGW));           // 補間添字1
                  std::get<1>(std::get<2>(ret)) += IGIGn * (-1. + std::get<0>(GWGW)) * (-1. + std::get<1>(GWGW));  // 補間添字2
               });

               std::get<1>(std::get<0>(ret)) *= c;
               std::get<1>(std::get<1>(ret)) *= c;
               std::get<1>(std::get<2>(ret)) *= c;
            }
            for_each(ret, [&](const auto &p_igign) {
               /*
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
               auto [p, igign] = p_igign;
               IGIGn_Row[PF2index(p, integ_f)] += igign;        // この面に関する積分において，φまたはφnの寄与
               if (p != origin)                                 // for use_rigid_mode
                  origin_ign_rigid_mode -= std::get<1>(igign);  // for use_rigid_mode
            });
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
#ifdef _OPENMP
#pragma omp parallel
#endif
         for (const auto &[PBF, i] : PBF_index)
#ifdef _OPENMP
#pragma omp single nowait
#endif
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
               Tdd igign;
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
         knowns = V_d(PBF_index.size());  //! ok
                                          //! Pointsの順番と合わせてとるように注意

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

            // b$ ------------------------------------------------------ */
            // b$　　　                    knownsの計算　                       */
            // b$ ------------------------------------------------------ */
            /*
             * Dot(mat_ukn,phiORphin) = Dot(mat_kn,knowns)
             * => phiORphin = Dot(mat_ukn^-1, Dot(mat_kn,knowns))
             */

            if (DorN == Dirichlet)
               knowns[i] = p->phi_Dirichlet = std::get<0>(p->phiphin);
            else
               knowns[i] = p->phinOnFace.at(f);  // std::get<1>(p->phiphin);

            if (!isFinite(knowns[i])) {
               std::cout << "p->Dirichlet" << p->Dirichlet << std::endl;
               std::cout << "p->CORNER" << p->CORNER << std::endl;
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "knowns is not finite");
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