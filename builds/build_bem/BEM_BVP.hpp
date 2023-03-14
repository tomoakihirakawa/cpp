#ifndef BEM_BVP_H
#define BEM_BVP_H

#include "Network.hpp"
// #define use_CG
// #define use_gmres 20
#define use_lapack
std::unordered_map<std::tuple<netP *, bool, netF *>, int> PBF_index;

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
   /* ------------------------------------------------------ */
   VV_d mat_ukn, mat_kn;
   V_d knowns;
   std::vector<std::vector<Tdd>> IGIGn;
   /* ------------------------------------------------------ */
   BEM_BVP(){};
   ~BEM_BVP() { delete this->lu; };
   /* ------------------------------------------------------ */
   /* ------------------------------------------------------ */
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
         if (p->Dirichlet || p->CORNER) {
            std::get<0>(p->phiphin_t) = p->DphiDt(0.) - Dot(p->U_BEM, p->U_BEM);
            phin_t_a = std::get<1>(p->phiphin_t_a6) = {0., 0., 0., 0., 0., 0.};  //! 接した物体に加速度に依存しない
         }
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
            for (auto &[f, phin_t_and_phin_t_a6] : p->phintOnFace) {
               auto &[phin_t, phin_t_a6] = phin_t_and_phin_t_a6;
               if (f) {
                  auto netInContact = NearestContactFace(f)->getNetwork();
                  auto w = netInContact->velocityRotational();
                  auto dQdt = Q.d_dt(w);
                  auto n = f->normal;
                  auto [p0, p1, p2] = f->getPoints(p);
                  Tddd phi012 = {std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)};
                  Tddd phin012 = {std::get<1>(p0->phiphin), std::get<1>(p1->phiphin), std::get<1>(p2->phiphin)};
                  Tddd grad_phi = Mean(phin012) * n + gradTangential_LinearElement(phi012, ToX(f));
                  // phin_t = std::get<1>(p->phiphin_t) = Dot(n, Dot(uNeumann(p, f) - grad_phi, dQdt.Rv()) + accelNeumann(p, f) - Dot(grad_phi, grad_U_LinearElement(f)));
                  phin_t = std::get<1>(p->phiphin_t) = Dot(w, uNeumann(p, f) - grad_phi) + Dot(n, accelNeumann(p, f) - Dot(grad_phi, grad_U_LinearElement(f)));
                  std::tie(std::get<0>(phin_t_a6),
                           std::get<1>(phin_t_a6),
                           std::get<2>(phin_t_a6)) = n;  //! 加速度で微分したphin_t_aはここまで簡単になる
                  std::tie(std::get<3>(phin_t_a6),
                           std::get<4>(phin_t_a6),
                           std::get<5>(phin_t_a6)) = Cross(n, p->X - netInContact->COM);  //! 加速度で微分したphin_t_aはここまで簡単になる
               } else {
                  auto netInContact = NearestContactFace(p)->getNetwork();
                  auto w = netInContact->velocityRotational();
                  auto dQdt = Q.d_dt(w);
                  // phin_t = std::get<1>(p->phiphin_t) = Dot(n, Dot(uNeumann(p) - p->U_BEM, dQdt.Rv()) + accelNeumann(p) - Dot(p->U_BEM, grad_U_LinearElement(p)));
                  phin_t = std::get<1>(p->phiphin_t) = Dot(w, uNeumann(p) - p->U_BEM) + Dot(n, accelNeumann(p) - Dot(p->U_BEM, grad_U_LinearElement(p)));
                  std::tie(std::get<0>(phin_t_a6),
                           std::get<1>(phin_t_a6),
                           std::get<2>(phin_t_a6)) = n;  //! 加速度で微分したphin_t_aはここまで簡単になる
                  std::tie(std::get<3>(phin_t_a6),
                           std::get<4>(phin_t_a6),
                           std::get<5>(phin_t_a6)) = Cross(n, p->X - netInContact->COM);  //! 加速度で微分したphin_t_aはここまで簡単になる
               }
            }
         }
      }
   };
   /* ------------------------------------------------------ */
   /* ------------------------------------------------------ */
   void solveForPhiPhin_t(const std::vector<Network *> &rigidbodies) const {
      // 繰り返し計算によって物体の加速度を計算する
      setPhiPhin_t();
      // b* ------------------------------------------------------ */
      // b*                         phiphin_t                      */
      // b* ------------------------------------------------------ */
      V_d knowns(PBF_index.size());
      V_d phiORphin_t(PBF_index.size(), 0);
      // b! ------------------------------------------------ */
      // b!   ddai {φn,φnt} ---BIE---> all ddai{φn,φnt}      */
      // b! ------------------------------------------------ */
      //@ ddai {φn,φnt} = {{0, 0, 0}, {., ., .}}
      for (auto target_body : floating_bodies) {
         double F = 0, dfda = 0;
         // 各成分について解く
         for (auto i = 0; i < 3; ++i) {
            double df = 0., dfda = ;
            for (const auto &[PBF, i] : PBF_index) {
               auto [p, B, face] = PBF;
               if (B == Dirichlet)
                  knowns[i] = 0.;
               else {
                  if (target_body == face->getNetwork()) {
                     //! 微分は値を持つ
                     auto [n0, n1, n2] = n;
                     if (i == 0)
                        knowns[i] = n0;
                     else if (i == 1)
                        knowns[i] = n1;
                     else if (i == 2)
                        knowns[i] = n2;
                  } else
                     knowns[i] = 0.;
               }
            }
            //
            this->lu->solve(Dot(mat_kn, knowns), phiORphin_t);
            //
            for (const auto &[PBF, i] : PBF_index) {
               auto [p, B, F] = PBF;
               knowns[i] = phiORphin_t[i];
            }
            F += f * f;
            dFda += 2 * f * (dfda - 1);
         }
      }
      //% ------------------------------------------------ */
      //%                     dpdaの計算                    */
      //% ------------------------------------------------ */
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, B, F] = PBF;
         if (B == Dirichlet)
            knowns[i] = std::get<0>(p->phiphin_t);
         else
            knowns[i] = std::get<0>(p->phintOnFace.at(F));  //? ここでphintOnFaceを使っているのはなぜか？
      }
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
      std::cout << "lapack for phiphin_t" << std::endl;
      this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin_t /*解*/);
      //
#endif
      std::cout << "solved" << std::endl;
      /* ------------------------------------------------------ */
      for (const auto &[PBF, i] : PBF_index) {
         auto [p, DorN, f] = PBF;
         if (DorN == Dirichlet)
            std::get<1>(p->phiphin_t) = phiORphin_t[i];
         else
            std::get<0>(p->phiphin_t) = phiORphin_t[i];
         p->pressure = p->pressure_BEM = _WATER_DENSITY_ * (-std::get<0>(p->phiphin_t) - _GRAVITY_ * p->height() - Dot(p->U_BEM, p->U_BEM) / 2.);
      }
   };

   void solve(const Network &water, const Buckets<networkPoint *> &FMM_BucketsPoints, const Buckets<networkFace *> &FMM_BucketsFaces) {
      //* ------------------------------------------------------ */
      //%                     各点で方程式を作る場合                 */
      //* ------------------------------------------------------ */
      std::cout << "各点で方程式を作る場合" << std::endl;
      PBF_index.clear();
      PBF_index.reserve(3 * water.getPoints().size());
      int i = 0;
      for (const auto &p : water.getPoints()) {
         if (!p->Neumann)
            PBF_index[{p, Dirichlet, nullptr}] = i++;
         for (const auto &[f, _] : p->phinOnFace)
            PBF_index[{p, Neumann, f}] = i++;
      }

      std::cout << Red << "water.getPoints() = " << i << std::endl;
      IGIGn = std::vector<std::vector<Tdd>>(PBF_index.size(), std::vector<Tdd>(PBF_index.size(), {0., 0.}));
      mat_kn = mat_ukn = VV_d(PBF_index.size(), V_d(PBF_index.size(), 0.));

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
         auto [p0, _, __] = PBF;
         auto it = PBF_index.begin();
         double p0_ign_rigid_mode = 0.;
         auto &IGIGn_Row = IGIGn[index];
         //
         std::tuple<std::tuple<netP *, Tdd>, std::tuple<netP *, Tdd>, std::tuple<netP *, Tdd>> ret;
         Tdd IGIGn, c;  //, p0igign = {0., 0.}, p1igign = {0., 0.}, p2igign = {0., 0.};
         double nr, tmp;
         Tddd X2, X0, X1, A, cross;
         auto origin = p0;
         for (const auto &integ_f : FMM_BucketsFaces.all_stored_objects) {
            {
               const auto [p0, p1, p2] = integ_f->getPoints(origin);
               ret = {{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}};
               X2 = p2->getXtuple();
               X0 = p0->getXtuple() - X2;
               X1 = p1->getXtuple() - X2;
               A = origin->getXtuple() - X2;
               cross = Cross(X0, X1);
               c = {Norm(cross), Dot(A, cross)};
               if (origin == p0 || origin == p1 || origin == p2)
                  std::get<1>(c) = 0;

               for_each(__GW5xGW5__, [&](const auto &GWGW) {
                  tmp = std::get<2>(GWGW) * (1. - std::get<0>(GWGW)) / (nr = Norm(X0 * std::get<0>(GWGW) + X1 * std::get<1>(GWGW) * (1. - std::get<0>(GWGW)) - A));
                  std::get<1>(std::get<0>(ret)) += (IGIGn = {tmp, tmp / (nr * nr)}) * std::get<0>(GWGW);
                  std::get<1>(std::get<1>(ret)) += IGIGn * std::get<1>(GWGW) * (1. - std::get<0>(GWGW));
                  std::get<1>(std::get<2>(ret)) += IGIGn * (-1. + std::get<0>(GWGW)) * (-1. + std::get<1>(GWGW));
               });
               std::get<1>(std::get<0>(ret)) *= c;
               std::get<1>(std::get<1>(ret)) *= c;
               std::get<1>(std::get<2>(ret)) *= c;
            }
            for_each(ret, [&](const auto &p_igign) {
               IGIGn_Row[(it = PBF_index.find({std::get<0>(p_igign), integ_f->Dirichlet, integ_f})) != PBF_index.end()
                             ? it->second
                             : PBF_index[{std::get<0>(p_igign), integ_f->Dirichlet, nullptr}]] += std::get<1>(p_igign);
               /* ------------------------------ */
               if (std::get<0>(p_igign) != p0)                             // for use_rigid_mode
                  p0_ign_rigid_mode -= std::get<1>(std::get<1>(p_igign));  // for use_rigid_mode
            });
         }
         /* -------------------------------------------------------------------------- */
#if defined(use_rigid_mode)
         std::get<1>(IGIGn_Row[index]) = p0_ign_rigid_mode;
#else
         /*
         @ ∇^2(1/r)=-1/(4pi)δ(r)
         @ IG*φn=-aφ+IGn*φ
         @ IG*φn=(IGn-a)*φ
         */
         std::get<1>(IGIGn_Row[index]) += p0->getSolidAngle();
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
            if (!(p->CORNER && DorN == Neumann /*変更する対象の行*/))  //! OK
            {
               Tdd igign;
               for (const auto &[PBF_j, j] : PBF_index) {
                  igign = IGIGn[i][j];
                  /* --------------------------------------- */
                  // 未知変数の係数行列は左，既知変数の係数行列は右
                  if (std::get<1>(PBF_j) == Neumann)
                     igign = {-std::get<1>(igign), -std::get<0>(igign)};
                  //% IGIGn は 左辺に IG*φn が右辺に IGn*φ が来るように計算しているため，移項する場合，符号を変える必要がある．
                  /*
                  IG*φn = IGn*φ

                  {IG0,IG1,IG2,IG3} . {φn0,φn1,φn2,φn3} = {IGn0,IGn1,IGn2,IGn3} . {φ0,φ1,φ2,φ3}
                  {IG0,-IGn1,IG2,IG3} . {φn0,φ1,φn2,φn3} = {IGn0,-IG1,IGn2,IGn3} . {φ0,φn1,φ2,φ3}
                  */
                  /* --------------------------------------- */
                  mat_ukn[i][j] = std::get<0>(igign);
                  mat_kn[i][j] = std::get<1>(igign);
               }
            }
         }

         /* ------------------------------------------------------ */
         double maxpp = 0;
         for (auto i = 0; i < mat_ukn.size(); ++i) {
            // LUするのはmat_uknだけなので，mat_knの最大値を使う必要はない
            //  if (maxpp < std::abs(mat_kn[i][i]))
            //  	maxpp = std::abs(mat_kn[i][i]);
            if (maxpp < std::abs(mat_ukn[i][i]))
               maxpp = std::abs(mat_ukn[i][i]);
         }
         /* ------------------------------------------------------ */
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
            if (p->CORNER && DorN == Neumann /*変更する対象の行*/)  //! OK
            {
               for (const auto &[PBF_j, j] : PBF_index) {
                  auto [p_, DorN_, f_] = PBF_j;
                  if (p == p_ && DorN_ == Neumann && f == f_) {
                     mat_ukn[i][j] = maxpp;  // φの系数
                     mat_kn[i][j] = 0;       // φnの系数
                  } else if (p == p_ && DorN_ == Dirichlet) {
                     mat_ukn[i][j] = 0;     // φnの系数
                     mat_kn[i][j] = maxpp;  // φの系数移行したからマイナス？　いいえ，移項を考慮した上でこれでいい．
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

         // b@ ------------------------------------------------------ */
      }

      {
         std::cout << "knownsの計算" << std::endl;
         // b$ ------------------------------------------------------ */
         // b$　　　                   knownsの計算　                    */
         // b$ ------------------------------------------------------ */
         /**
          * Dot(mat_ukn,phiORphin) = Dot(mat_kn,knowns)
          * => phiORphin = Dot(mat_ukn^-1, Dot(mat_kn,knowns))
          */
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
            auto [p, DorN, F] = PBF;
            if (DorN == Dirichlet)
               knowns[i] = p->phi_Dirichlet = std::get<0>(p->phiphin);
            else
               knowns[i] = p->phinOnFace.at(F);  // std::get<1>(p->phiphin);

            if (!isFinite(knowns[i])) {
               std::cout << "p->Dirichlet" << p->Dirichlet << std::endl;
               std::cout << "p->CORNER" << p->CORNER << std::endl;
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "knowns is not finite");
            }
         }

         // b$ ------------------------------------------------------ */
      }

      {
         // b% ------------------------------------------------------ */
         // b%                  境界積分方程式を解く                      */
         // b% ------------------------------------------------------ */
         std::cout << "--------------------- 境界積分方程式を解く ---------------------" << std::endl;
         /* ------------------------------------------------------ */
         V_d phiORphin(knowns.size(), 0);
         std::cout << "IGIGn.size()= " << IGIGn.size() << std::endl;
         std::cout << "phiORphin.size()= " << phiORphin.size() << std::endl;
         std::cout << "  mat_kn.size() = " << mat_kn.size() << std::endl;
         //* 未知変数の計算
         /* ------------------------------------------------------ */
         // 前処理
         auto v = diagonal_scaling_vector(mat_ukn);
         for (auto i = 0; i < v.size(); ++i) {
            mat_ukn[i] *= v[i];
            mat_kn[i] *= v[i];
         }
         /* ------------------------------------------------------ */
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
         //
         // /* ------------------------------------------------------ */
         // if (isFinite(err) && err < 0.1)
         // else
         /* ------------------------------------------------------ */
#elif defined(use_lapack)
         std::cout << "lapack lu decomposition" << std::endl;
         // this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/);
         // std::cout << "try to solve" << std::endl;
         // this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
         //
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

         for (const auto &[PBF, i] : PBF_index) {
            auto [p, DorN, f] = PBF;
            if (DorN == Dirichlet)
               std::get<1>(p->phiphin) = p->phin_Dirichlet = phiORphin[i];
            else
               std::get<0>(p->phiphin) = phiORphin[i];
         }
         if (!isFinite(phiORphin))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "phiORphin is not finite");
      }
   };
};

#endif