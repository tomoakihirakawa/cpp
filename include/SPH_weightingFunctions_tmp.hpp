#ifndef SPH_weightingFunctions_H
#define SPH_weightingFunctions_H

#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

#define REFLECTION

#define SPP_U 1
#define SPP_pressure_coef -1
#define SPP_pressure 1

#define surface_zero_pressure 0

#define POWER 1.

auto canSetSPP(Network *net, const auto &RigidBodyObject, const auto &p, const Tddd &X, const double &range) {
   {
      const double C = 1.2;
      auto func = [&](const auto &q) { return p != q && Distance(q, X) < C * range; };
      if (net->BucketPoints.any_of(X, C * range, func))
         return false;
   }
   {
      const double C = 1.;
      auto func = [&](const auto &q) { return p != q && Distance(q, X) < C * range; };
      for (const auto &[obj, poly] : RigidBodyObject)
         if (obj->BucketPoints.any_of(X, C * range, func))
            return false;
   }
   return true;
};

// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */

void setNormal_Surface_(auto &net, const std::unordered_set<networkPoint *> &wall_p,
                        const auto &RigidBodyObject,
                        const bool set_surface = true) {
   // b# ------------------------------------------------------ */
   // b#             流体粒子の法線方向の計算，水面の判定              */
   // b# ------------------------------------------------------ */
   // b#  A. Krimi, M. Jandaghian, and A. Shakibaeinia, Water (Switzerland), vol. 12, no. 11, pp. 1–37, 2020.
   DebugPrint("水粒子のオブジェクト外向き法線方向を計算", Green);
#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      // p->interpolated_normal_SPH, q->X - p->Xの方向が完全に一致した際に失敗する
      /* ---------------------- p->interpolated_normal_SPHの計算 --------------------- */
      p->interpolated_normal_SPH_original = {0., 0., 0.};
      p->interpolated_normal_SPH_original_all = {0., 0., 0.};
      net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         if (Between(Distance(p, q), {1E-10, p->radius_SPH})) {
            p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            p->interpolated_normal_SPH_original_all -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         }
      });
      std::vector<Tddd> wall_tangent;
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            if (Distance(p, q) < p->radius_SPH / p->C_SML * 1.1)
               wall_tangent.emplace_back(q->normal_SPH);
            p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            p->interpolated_normal_SPH_original_all -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         });

      p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
      p->interpolated_normal_SPH_all = Normalize(p->interpolated_normal_SPH_original_all);
      // for (const auto &wall_v : wall_tangent) {
      //    p->interpolated_normal_SPH = Normalize(Chop(p->interpolated_normal_SPH, wall_v));
      //    p->interpolated_normal_SPH_all = Normalize(Chop(p->interpolated_normal_SPH_all, wall_v));
      // }
      if (!isFinite(p->interpolated_normal_SPH)) {
         p->interpolated_normal_SPH = {0., 0., 1.};
         p->interpolated_normal_SPH_all = {0., 0., 1.};
      }
      /* ----------------------------------- 検索 ----------------------------------- */
      if (set_surface) {
         p->isSurface = true;
         if (net->BucketPoints.any_of(p->X, (p->radius_SPH / p->C_SML) * 3., [&](const auto &q) {
                {
                   if (Distance(p, q) < (p->radius_SPH / p->C_SML) * 3.) {
                      return p != q && (VectorAngle(p->interpolated_normal_SPH_all, q->X - p->X) < M_PI / 4);
                   } else
                      return false;
                }
                return false;
             }))
            p->isSurface = false;

         if (p->isSurface)
            for (const auto &[obj, poly] : RigidBodyObject)
               if (obj->BucketPoints.any_of(p->X, (p->radius_SPH / p->C_SML) * 3., [&](const auto &q) {
                      {
                         if (Distance(p, q) < (p->radius_SPH / p->C_SML) * 3.) {
                            return p != q && (VectorAngle(p->interpolated_normal_SPH_all, -q->normal_SPH) < M_PI / 180. * 60);
                         } else
                            return false;
                      }
                      return false;
                   }))
                  p->isSurface = false;

         if (surface_zero_pressure && p->isSurface)
            p->p_SPH = 0;
      }
   }

   DebugPrint("壁粒子のオブジェクト外向き法線方向を計算", Green);
   for (const auto &[obj, poly] : RigidBodyObject)
      for (const auto &p : obj->getPoints())
         p->interpolated_normal_SPH_original = {0, 0, 0};
#pragma omp parallel
   for (const auto &p : wall_p) {
#pragma omp single nowait
      {
         p->interpolated_normal_SPH_original = {0., 0., 0.};
         for (const auto &[obj, poly] : RigidBodyObject)
            obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
               if (Between(Distance(p, q), {1E-8, p->radius_SPH}))
                  p->interpolated_normal_SPH_original -= grad_w_Bspline(p->X, q->X, p->radius_SPH);
            });
         p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
         if (!isFinite(p->interpolated_normal_SPH))
            p->interpolated_normal_SPH = {0., 0., 1.};
      }
   }
};
/* -------------------------------------------------------------------------- */
void mapValueOnWall(auto &net,
                    const std::unordered_set<networkPoint *> &wall_p,
                    const auto &RigidBodyObject,
                    const Tiii &indicater = {1, 0, 1}) {
   auto [first, second, do_add] = indicater;
   Timer watch;
   auto initialize = [](networkPoint *PW) {
      PW->div_U_ = PW->volume_ = PW->rho_ = PW->p_SPH_ = PW->total_weight = 0;
      PW->DUDt_SPH_ = PW->U_SPH_ = PW->tmp_U_SPH_ = {0., 0., 0.};
   };

   auto initializeAll = [&]() {
      for (const auto &PW : wall_p)
         initialize(PW);
   };

   double POW, a = 1.;
   bool applyToWall = false;
   bool mirroring = false;
   bool account_volume = false;
   bool do_spp = false;
   auto calc = [&]() {
#pragma omp parallel
      for (const auto &PW : wall_p)
#pragma omp single nowait
      {
         Tddd markerX = PW->X + 2 * PW->normal_SPH,
              accel = {0., 0., 0.},
              mirror_wall = PW->X + PW->normal_SPH,
              mirrored_X,
              n = Normalize(PW->normal_SPH);
         double W;
         int count = 0;
         initialize(PW);
         /* -------------------------------------------------------------------------- */
         auto func = [&](const networkPoint *PF, const Tddd &nearX, bool spp = false) {
            PW->total_N += std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW);
            PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
            PW->p_SPH_ += (spp ? SPP_pressure_coef : 1.) * PF->p_SPH * W;
            PW->DUDt_SPH_ += PF->DUDt_SPH * W;
            PW->rho_ += PF->rho * W;
            PW->U_SPH_ += PF->U_SPH * W;
            PW->tmp_U_SPH_ += PF->tmp_U_SPH * W;
            PW->volume_ += PF->volume * W;
            PW->div_U_ += PF->div_U * W;
            //
            if (mirroring && !applyToWall) {
               mirrored_X = Mirror(nearX, mirror_wall, PW->normal_SPH);
               PW->total_N += std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW);
               PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW));
               PW->p_SPH_ += (spp ? SPP_pressure_coef : 1.) * PF->p_SPH * W;
               PW->DUDt_SPH_ += PF->DUDt_SPH * W;
               PW->rho_ += PF->rho * W;
               PW->U_SPH_ += PF->U_SPH * W;
               PW->tmp_U_SPH_ += PF->tmp_U_SPH * W;
               PW->volume_ += PF->volume * W;
               PW->div_U_ += PF->div_U * W;
            }
            count++;
         };
         net->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &PF) {
            func(PF, PF->X);
            if (do_spp && SPP_pressure && PF->isSurface) {
               auto ps = PF->radius_SPH / PF->C_SML;
               auto nearX = (PF->X + PF->interpolated_normal_SPH * ps);
               if (canSetSPP(net, RigidBodyObject, PF, nearX, ps))
                  func(PF, nearX, true);
            }
         });
         if (applyToWall) {
            for (const auto &[obj, poly] : RigidBodyObject)
               obj->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &qW) {
                  func(qW, qW->X);
               });
            PW->total_weight = 1;
         }
         /* -------------------------------------------------------------------------- */
      }
   };
   //-------------------------
   auto set = [&]() {
      for (const auto &PW : wall_p) {
         {
            if (PW->total_weight > 1E-15) {
               PW->p_SPH = PW->p_SPH_ / PW->total_weight;
               PW->DUDt_SPH = PW->DUDt_SPH_ / PW->total_weight;
               PW->U_SPH = PW->U_SPH_ / PW->total_weight;
               PW->tmp_U_SPH = PW->tmp_U_SPH_ / PW->total_weight;
               PW->div_U = PW->div_U_ / PW->total_weight;
               PW->setDensityVolume(PW->rho_ / PW->total_weight, PW->volume_ / PW->total_weight);
            } else
               initialize(PW);
         }
         std::vector<bool> Bools = {!isFinite(PW->p_SPH), !isFinite(PW->rho), !isFinite(PW->volume), !isFinite(PW->DUDt_SPH), !isFinite(PW->U_SPH), !isFinite(PW->tmp_U_SPH)};
         if (std::any_of(Bools.begin(), Bools.end(), [](bool b) { return b; })) {
            std::stringstream ss;
            ss << Bools;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite" + ss.str());
         }
      }
   };
   //
   //! ------------------------------- first time ------------------------------ */
   if (first) {
      POW = 1.;  // ２などとしてはいけない
      a = 1.;
      initializeAll();
      //
      applyToWall = false;
      mirroring = true;
      account_volume = true;
      do_spp = false;
      calc();
      set();
      setNormal_Surface_(net, wall_p, RigidBodyObject);
      do_spp = true;
      calc();
      set();
      setNormal_Surface_(net, wall_p, RigidBodyObject);
   }
   /* -------------------------------------------------------------------------- */
   //! ------------------------------- second time ------------------------------ */
   if (second) {
      auto set = [&]() {
         for (const auto &PW : wall_p) {
            {
               const double alpha = 0.5;
               if (PW->total_weight > 1E-15) {
                  PW->p_SPH = (1 - alpha) * PW->p_SPH + alpha * PW->p_SPH_ / PW->total_weight;
                  PW->DUDt_SPH = (1 - alpha) * PW->DUDt_SPH + alpha * PW->DUDt_SPH_ / PW->total_weight;
                  PW->U_SPH = (1 - alpha) * PW->U_SPH + alpha * PW->U_SPH_ / PW->total_weight;
                  PW->tmp_U_SPH = (1 - alpha) * PW->tmp_U_SPH + alpha * PW->tmp_U_SPH_ / PW->total_weight;
                  PW->setDensityVolume((1 - alpha) * PW->rho + alpha * PW->rho_ / PW->total_weight, (1 - alpha) * PW->volume + alpha * PW->volume_ / PW->total_weight);
               } else
                  initialize(PW);
            }
            std::vector<bool> Bools = {!isFinite(PW->p_SPH), !isFinite(PW->rho), !isFinite(PW->volume), !isFinite(PW->DUDt_SPH), !isFinite(PW->U_SPH), !isFinite(PW->tmp_U_SPH)};
            if (std::any_of(Bools.begin(), Bools.end(), [](bool b) { return b; })) {
               std::stringstream ss;
               ss << Bools;
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite" + ss.str());
            }
         }
      };
      //
      POW = 1.;
      a = 1.;
      initializeAll();
      applyToWall = true;
      mirroring = false;
      account_volume = true;
      do_spp = true;
      calc();
      set();
      //
      POW = 1.;
      a = 1.;
      initializeAll();
      applyToWall = true;
      mirroring = false;
      account_volume = true;
      do_spp = true;
      calc();
      set();
      //
      POW = 1.;
      a = 1.;
      initializeAll();
      applyToWall = true;
      mirroring = false;
      account_volume = true;
      do_spp = true;
      calc();
      set();
      //
      POW = 1.;
      a = 1.;
      initializeAll();
      applyToWall = true;
      mirroring = false;
      account_volume = true;
      do_spp = true;
      calc();
      set();
   }
   //! ------------------------------- apply condition! ------------------------------ */
   for (const auto &PW : wall_p) {

      if (Norm(PW->normal_SPH) < 1E-12) {
         PW->U_SPH = Chop(PW->U_SPH, PW->normal_SPH);
         PW->tmp_U_SPH = Chop(PW->tmp_U_SPH, PW->normal_SPH);
      }

      // free-slip
      // PW->U_SPH = Reflect(PW->U_SPH, PW->normal_SPH);
      // PW->tmp_U_SPH = Reflect(PW->tmp_U_SPH, PW->normal_SPH);

      // no-slip
      PW->U_SPH *= -1.;
      PW->tmp_U_SPH *= -1.;
      //
      // all-slip
      // 壁面上に粒子を設定した場合，壁上で流速とDUDtをゼロとするのは，妥当な設定
      // PW->DUDt_SPH *= 0.;これはおかしい，ここでのDUDtは粘性と重力のみを考慮したものにすぎない．
      if (Norm(PW->normal_SPH) < 1E-12) {
         PW->U_SPH *= 0.;
         PW->tmp_U_SPH *= 0.;
      }
      //
      if (do_add) {
         Tddd accel = {0., 0., 0.};
         PW->p_SPH = PW->p_SPH + PW->rho * Dot(PW->DUDt_SPH - accel, -2 * PW->normal_SPH);
      }
   }
   DebugPrint(_green, "Elapsed time: ", _red, watch(), "s ", Magenta, "壁面粒子の圧力を計算", (do_add ? ", DUDtを考慮した圧力" : ""));
};

/* -------------------------------------------------------------------------- */

const double damping_factor = 1.;
const double asobi = 0.;

std::unordered_set<networkPoint *> wall_p;

void developByEISPH(Network *net,
                    const auto &RigidBodyObjectIN,
                    double &real_time,
                    const double C_SML,
                    const double particle_spacing,
                    const double max_dt,
                    const int RK_order) {
   try {
      Timer watch;
      std::vector<std::tuple<Network *, Network *>> RigidBodyObject;
      for (const auto &[a, b, _] : RigidBodyObjectIN)
         RigidBodyObject.push_back({a, b});
      //% --------------------------- 平滑化距離の計算 ------------------------------ */
      /*     密度, 平滑化距離      */
      DebugPrint(Green, "固定の平滑化距離の計算: C_SML * particle_spacing = ", C_SML, " * ", particle_spacing, " = ", C_SML * particle_spacing);
      for (const auto &p : net->getPoints()) {
         p->setDensity(_WATER_DENSITY_);
         p->isFreeFalling = false;
         p->isInsideOfBody = false;
         p->p_SPH_SPP = 0;
         p->DUDt_SPH = p->lap_U = {0, 0, 0};
         p->isCaptured = true;
      }
      //@ ------------------------------------------------------- */
      //@                    関連する壁粒子をマーク                   */
      //@ ------------------------------------------------------- */
      wall_p.clear();
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            p->setDensityVolume(0, 0);
            p->isFreeFalling = false;
            p->isCaptured = false;
            p->isSurface = false;
            p->p_SPH = 0;
            p->U_SPH = p->DUDt_SPH = p->lap_U = {0, 0, 0};
         }
      DebugPrint("関連する壁粒子をマーク", Green);
      double C = 1.2;
      for (const auto &[obj, poly] : RigidBodyObject) {
#pragma omp parallel
         for (const auto &p : net->getPoints())
#pragma omp single nowait
         {
            obj->BucketPoints.apply(p->X, p->radius_SPH * C, [&](const auto &q) {
               if (Distance(p, q) < p->radius_SPH * C) {
                  q->isCaptured = true;
                  q->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3.));
               }
            });
         }
      };
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints())
            if (p->isCaptured)
               wall_p.emplace(p);

      DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "関連する壁粒子をマークし，保存");
      //% -------------- バケットの生成, p->radius_SPHの範囲だけ点を取得 --------------- */
      DebugPrint("バケットの生成", Green);
      net->makeBucketPoints(particle_spacing * 0.8);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->makeBucketPoints(particle_spacing * 0.8);
      DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "バケットの生成");
      //% -------------------------------------------------------------------------- */
      //% --------------- CFL条件を満たすようにタイムステップ間隔dtを設定 ----------------- */
      //% -------------------------------------------------------------------------- */
      DebugPrint("近傍粒子探査が終わったら時間ステップを決めることができる", Green);
      double dt = max_dt;
      auto C_CFL_velocity = 0.02;  // dt = C_CFL_velocity*h/Max(U)
      auto C_CFL_accel = 0.1;      // dt = C_CFL_accel*sqrt(h/Max(A))
      for (const auto &p : net->getPoints()) {
         // 速度に関するCFL条件
         auto dt_C_CFL = [&](const auto &q) {
            if (p != q) {
               auto pq = Normalize(p->X - q->X);
               auto distance = Distance(p, q);
               /* ------------------------------------------------ */
               // 相対速度
               double max_dt_vel = C_CFL_velocity * distance / std::abs(Dot(p->U_SPH - q->U_SPH, pq));
               // double max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH - q->U_SPH);
               if (dt > max_dt_vel && isFinite(max_dt_vel))
                  dt = max_dt_vel;
               // 絶対速度
               max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH);
               if (dt > max_dt_vel && isFinite(max_dt_vel))
                  dt = max_dt_vel;
               /* ------------------------------------------------ */
               // 相対速度
               double max_dt_acc = C_CFL_accel * std::sqrt(distance / std::abs(Dot(p->DUDt_SPH - q->DUDt_SPH, pq)));
               // double max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH - q->DUDt_SPH));
               if (dt > max_dt_acc && isFinite(max_dt_acc))
                  dt = max_dt_acc;
               // 絶対速度
               max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH));
               if (dt > max_dt_acc && isFinite(max_dt_acc))
                  dt = max_dt_acc;
            }
         };
         net->BucketPoints.apply(p->X, p->radius_SPH, dt_C_CFL);
         for (const auto &[obj, poly] : RigidBodyObject)
            obj->BucketPoints.apply(p->X, p->radius_SPH, dt_C_CFL);
         double max_dt_vel = C_CFL_velocity * (p->radius_SPH / p->C_SML) / Norm(p->U_SPH);
         if (dt > max_dt_vel && isFinite(max_dt_vel))
            dt = max_dt_vel;
         double max_dt_acc = C_CFL_accel * std::sqrt((p->radius_SPH / p->C_SML) / Norm(p->DUDt_SPH));
         if (dt > max_dt_acc && isFinite(max_dt_acc))
            dt = max_dt_acc;
      }
      std::cout << "dt = " << dt << std::endl;
      //% -------------------------------------------------------------------------- */
      //% -------------------------------------------------------------------------- */
      //@ ----------------------- ルンゲクッタの準備 ------------------- */
      for (const auto &p : net->getPoints()) {
         p->RK_U.initialize(dt, real_time, p->U_SPH, RK_order);
         p->RK_X.initialize(dt, real_time, p->X, RK_order);
         p->RK_P.initialize(dt, real_time, p->p_SPH, RK_order);
      }
      do {
         for (const auto &p : net->getPoints())
            p->setDensity(_WATER_DENSITY_);
         //
         // b@ ======================================================= */
         // b@                  ルンゲクッタを使った時間積分                  */
         // b@ ======================================================= */
         dt = (*net->getPoints().begin())->RK_X.getdt();
         std::cout << "dt = " << dt << std::endl;
         /* -------------------------------------------------------------------------- */
         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, " 流体粒子の法線方向の計算，水面の判定");
         // b! ========================================================================== */
         // b!                          U, DUDt, DPDtを計算                                */
         // b! ========================================================================== */
         {
            /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/
            // ! ------------------------------------------------------ */
            // !                    ∇.∇UとU*を計算                       */
            // ! ------------------------------------------------------ */
            mapValueOnWall(net, wall_p, RigidBodyObject);

#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               Tddd viscous_term = p->lap_U = {0, 0, 0}, rij, Uij, qX;
               double nu_nu;

               auto func = [&](const auto &q) {
                  if (q->isCaptured) {

                     auto nu_lap_U = [&](const Tddd &qX) {
                        rij = p->X - qX;
                        Uij = p->U_SPH - q->U_SPH;
                        nu_nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
                        viscous_term += q->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(p->X, qX, p->radius_SPH) /
                                        ((q->rho + p->rho) * Dot(rij, rij) + std::pow(1E-4 * p->radius_SPH / p->C_SML, 2));
                     };

                     if (p != q)
                        nu_lap_U(q->X);
                     if (SPP_U && q->isSurface) {
                        qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
                        if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-10, p->radius_SPH}))
                           nu_lap_U(qX);
                     }
                  }
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               p->lap_U = p->lap_U_ = viscous_term / (p->mu_SPH / p->rho);
               p->DUDt_SPH = p->DUDt_SPH_ = viscous_term + _GRAVITY3_;  // 後で修正されるDUDt
               p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
               p->tmp_X = p->X + p->tmp_U_SPH * dt;

               // b# -------------------------------------------------------------------------- */
               // b#                           DUDt_SPHとtmp_U_SPHの修正                          */
               // b# -------------------------------------------------------------------------- */

               // #if defined(REFLECTION)
               //                bool isReflect = true;
               //                int count = 0;
               //                while (isReflect && count++ < 10) {
               //                   const auto X = p->tmp_X;
               //                   isReflect = false;
               //                   //
               //                   auto closest = [&]() {
               //                      double distance = 1E+20;
               //                      networkPoint *P = nullptr;
               //                      for (const auto &[obj, poly] : RigidBodyObject) {
               //                         obj->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
               //                            auto tmp = Distance(X, q);
               //                            if (distance > tmp) {
               //                               distance = tmp;
               //                               P = q;
               //                            }
               //                         });
               //                      }
               //                      return P;
               //                   };
               //                   auto p_wall = closest();
               //                   if (p_wall) {
               //                      const auto n = Normalize(p_wall->normal_SPH);
               //                      if (particle_spacing > Distance(p_wall, X))
               //                      // if (Dot(p_wall->X + p_wall->normal_SPH + n * particle_spacing * (1 / 2. + asobi) - X, n) > 0) {
               //                      // if (Dot(p_wall->X + 2. * p_wall->normal_SPH - X, n) > 0)
               //                      {
               //                         //
               //                         // if (Dot(p->tmp_U_SPH, n) < 0) {
               //                         //    // p->DUDt_SPH -= 0.5 * Projection(p->tmp_U_SPH, n) / dt;
               //                         //    // p->DUDt_SPH_ = p->DUDt_SPH;
               //                         //    // p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
               //                         //    p->tmp_X += Projection((p_wall->X + p_wall->normal_SPH + n * particle_spacing / 2.) - p->tmp_X, n);
               //                         //    isReflect = true;
               //                         // }
               //                         if (Dot(p->tmp_U_SPH, n) < 0) {
               //                            p->DUDt_SPH -= (1 + damping_factor) * Projection(p->tmp_U_SPH, n) / dt;
               //                            p->DUDt_SPH_ = p->DUDt_SPH;
               //                            p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
               //                            p->tmp_X = p->X + p->tmp_U_SPH * dt;
               //                            isReflect = true;
               //                         }
               //                      }
               //                   }
               //                };
               // #endif
            }
            // );
            /* -------------------------------------------------------------------------- */
            DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "粘性項の∇.∇Uを計算し，次にU*を計算");
            /* --------------------------------------------------------- */
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
               p->setX(p->tmp_X);
            // p->setX(p->X + p->tmp_U_SPH * dt);  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DebugPrint("U*を使って，仮想的な位置X^*へ粒子を移動", Green);
            // ! ------------------------------------------------------ */
            // !               div(U^*), DρDt=-ρdiv(U)の計算              */
            // ! ------------------------------------------------------ */

            auto loop = [&](const auto &p) {
               p->div_U = p->DrhoDt_SPH = 0.;
               Tddd qp, qX, Uij;
               auto func = [&](const auto &q) {
                  if (q->isCaptured) {
                     Uij = q->tmp_U_SPH - p->tmp_U_SPH;

                     auto div_U = [&](const auto &qX) {
                        p->div_U += q->mass / p->rho * Dot(Uij, grad_w_Bspline(p->X, qX, p->radius_SPH));  // 後藤p.25 (2.89)
                        //-----------
                        // qp = qX - p->X;
                        // p->div_U += q->mass / p->rho * Dot(Uij, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH) / 2.;
                     };

                     if (p != q)
                        div_U(q->X);
                     if (SPP_U && q->isSurface) {
                        qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
                        if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML))
                           if (Between(Distance(p, qX), {1E-10, p->radius_SPH}))
                              div_U(qX);
                     }
                  }
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               if (!isFinite(p->div_U))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "div_U is not a finite");

               p->DrhoDt_SPH = -p->rho * p->div_U;
            };

#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
               loop(p);

            for (const auto &p : net->getPoints()) {
               p->setDensity(p->rho_ = p->rho + p->DrhoDt_SPH * dt);  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               // p->setDensity(_WATER_DENSITY_ * (1. - p->div_U * dt));
            }

            mapValueOnWall(net, wall_p, RigidBodyObject);

            // #pragma omp parallel
            //             for (const auto &[obj, poly] : RigidBodyObject)
            // #pragma omp single nowait
            //                for (const auto &p : obj->getPoints())
            //                   if (Norm(p->normal_SPH) < 1E-12 && p->isCaptured)
            //                      loop(p);
            // ここで壁表層のdivUも計算した，この結果を更新せず使うように

            //
            /* -------------------------------- もう一度 ---------------------------------- */
            //             //
            //             // setWallDensityVolume(net, wall_p, RigidBodyObject);
            //             // setNormal_Surface();
            //             // mapValueOnWall(net, wall_p, RigidBodyObject);
            //             // setWalltmpU_and_U_SPH(net, wall_p, RigidBodyObject);
            //             mapValueOnWall(net, wall_p, RigidBodyObject);

            // #pragma omp parallel
            //             for (const auto &p : net->getPoints())
            // #pragma omp single nowait
            //             {
            //                p->div_U = p->DrhoDt_SPH = 0.;
            //                Tddd qp, qX;
            //                auto func = [&](const auto &q) {
            //                   if (q->isCaptured) {
            //                      if (p != q) {
            //                         // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
            //                         p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
            //                         // qp = q->X - p->X;
            //                         // p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH);
            //                      }
            //                      //
            //                      if (SPP_U && q->isSurface) {
            //                         qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
            //                         // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
            //                         if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML)) {
            //                            // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
            //                            if (Between(Distance(p, qX), {1E-10, p->radius_SPH})) {
            //                               // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, qX, p->radius_SPH));
            //                               p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, qX, p->radius_SPH));
            //                               //
            //                               // qp = qX - p->X;
            //                               // p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH);
            //                            }
            //                         }
            //                      }
            //                   }
            //                };

            //                net->BucketPoints.apply(p->X, p->radius_SPH, func);
            //                for (const auto &[obj, poly] : RigidBodyObject)
            //                   obj->BucketPoints.apply(p->X, p->radius_SPH, func);

            //                if (!isFinite(p->div_U))
            //                   throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "div_U is not a finite");

            //                // p->DrhoDt_SPH = -p->rho * p->div_U;
            //             }

            /* -------------------------------------------------------------------------- */
            // setWallDensityVolume(net, wall_p, RigidBodyObject);
            // setNormal_Surface();
            // mapValueOnWall(net, wall_p, RigidBodyObject);
            // setWalltmpU_and_U_SPH(net, wall_p, RigidBodyObject);
            // setPressureSPP(net, -1);
            // setWallPressure(net, wall_p, RigidBodyObject, 1);
            // setPressureSPP(net, -1);
            DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "div(U^*)を計算");
            // ! ------------------------------------------------------ */
            // !  　　　　　        仮位置における圧力Pの計算                 */
            // ! ------------------------------------------------------ */
            /*
            ここで，疑問
            流体の圧力は，仮の位置に移動する前のもの．壁の圧力も移動前のもの．
            仮の位置に移動した後に，流体の圧力を計算するが，そこで利用する圧力は，移動前の圧力でよいのか？
            */

            std::unordered_map<networkPoint *, NewtonRaphson<double> *> P2Newton;
            for (const auto &p : net->getPoints())
               P2Newton[p] = (new NewtonRaphson(p->p_SPH));

            for (const auto &[obj, poly] : RigidBodyObject)
               for (const auto &p : obj->getPoints())
                  if (Norm(p->normal_SPH) < 1E-12 && p->isCaptured)
                     P2Newton[p] = (new NewtonRaphson(p->p_SPH));

#define Morikawa2019
// #define NewtonMethod
// DebugPrint("仮位置における圧力Pの計算", Magenta);
#if defined(NewtonMethod)
            for (auto i = 0; i < 100; ++i)
#endif
            {
#pragma omp parallel
               for (const auto &[p, newton] : P2Newton)
#pragma omp single nowait
               {
                  int count = 0;
                  Tddd Xij, qX;
                  double Aij, sum_Aij = 0, sum_Aij_Pj = 0, sum_Aij_Pij = 0, sum_Pij = 0;
                  //% -------------------------------------------------------------------------- */
                  auto FUNC = [&](const auto &q) {
                     if (q->isCaptured) {
                        auto func = [&](const auto &qX, bool coef = 1.) {
                           Xij = p->X - qX;
#if defined(Morikawa2019)
                           // Morikawa, D., Senadheera, H., & Asai, M. (2021). Explicit incompressible smoothed particle hydrodynamics in a multi-GPU environment for large-scale simulations. Computational Particle Mechanics, 8(3), 493–510. https://doi.org/10.1007/s40571-020-00347-0
                           Aij = 2. * q->mass / p->rho * Dot(Xij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(Xij, Xij);
#else
                           Aij = 8. * q->mass / std::pow(q->rho + p->rho, 2) * Dot(Xij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(Xij, Xij);
#endif
                           sum_Pij += p->p_SPH - coef * q->p_SPH;
                           sum_Aij += Aij;
                           sum_Aij_Pj += Aij * coef * q->p_SPH;
                           sum_Aij_Pij += Aij * (p->p_SPH - coef * q->p_SPH);
                           count++;
                        };

                        if (p != q)
                           func(q->X);
                        if (SPP_pressure && q->isSurface) {
                           qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
                           if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-8, p->radius_SPH}))
                              func(qX, SPP_pressure_coef);
                        }
                     }
                  };
                  //% -------------------------------------------------------------------------- */
                  net->BucketPoints.apply(p->X, p->radius_SPH, FUNC);
                  for (const auto &[obj, poly] : RigidBodyObject)
                     obj->BucketPoints.apply(p->X, p->radius_SPH, FUNC);

#if defined(Morikawa2019)
                  auto b = _WATER_DENSITY_ * p->div_U / dt;
#else
                  auto b = p->div_U / dt;
#endif

                  if (std::abs(sum_Aij) < 1E-20)
                     p->p_SPH_ = 0;
                  else
                     p->p_SPH_ = (b + sum_Aij_Pj) / sum_Aij;

#if defined(NewtonMethod)
                  auto F = sum_Aij_Pij - b;
                  auto dFdPi = sum_Aij;
                  auto FF = F * F;
                  auto dFFdPi = 2 * F * dFdPi;
                  newton->update(FF, dFFdPi, 2.);
#endif

                  if (!isFinite(p->p_SPH_)) {
                     std::cout << "p->div_U " << p->div_U << std::endl;
                     std::cout << "sum_Aij_Pj " << sum_Aij_Pj << std::endl;
                     std::cout << "sum_Aij " << sum_Aij << std::endl;
                     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");
                  }
               }
#if defined(NewtonMethod)
               double sum = 0;
               for (const auto &[p, newton] : P2Newton) {
                  p->p_SPH = p->p_SPH_ = newton->X;
                  sum += std::abs(newton->dX);
               }
               if (i % 5 == 0)
                  std::cout << Red << "sum " << sum << colorOff << std::endl;

               if (sum < 1E-5) {
                  std::cout << Red << "i=" << i << ", sum " << sum << colorOff << std::endl;
                  break;
               }
#endif
            }

            for (auto &[p, newton] : P2Newton)
               delete newton;

            DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "仮位置における圧力Pの計算");
            // ! ------------------------------------------------------ */
            // !           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
            // ! ------------------------------------------------------ */

            DebugPrint("圧力勾配∇Pを計算 & DU/Dtの計算", Magenta);
            // mapValueOnWall(net, wall_p, RigidBodyObject);

#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               double p_rho2;
               int count = 0;
               Tddd qX;
               p->contact_points_all_SPH = p->contact_points_fluid_SPH = 0;
               p->gradP_SPH = {0., 0., 0.};

               auto func = [&](const auto &q) {
                  if (q->isCaptured) {

                     auto grad_P = [&](const auto &qX) {
                        p_rho2 = q->p_SPH / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                        p->gradP_SPH += p->rho * q->mass * p_rho2 * grad_w_Bspline(p->X, qX, p->radius_SPH);
                     };

                     if (p != q)
                        grad_P(q->X);
                     if (SPP_pressure && q->isSurface) {
                        qX = q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML;
                        if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-8, p->radius_SPH}))
                           grad_P(qX);
                     }
                  }
                  p->contact_points_all_SPH++;
                  if (q->getNetwork()->isFluid)
                     p->contact_points_fluid_SPH++;
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               // p->DUDt_SPH -= p->gradP_SPH / p->rho;
               p->DUDt_SPH -= p->gradP_SPH / _WATER_DENSITY_;
               //
               if (!isFinite(p->DUDt_SPH))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
            }
            DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "圧力勾配∇Pを計算 & DU/Dtの計算");
         }

         //@ -------------------------------------------------------- */
         //@                        粒子の時間発展                      */
         //@ -------------------------------------------------------- */
         DebugPrint("粒子の時間発展", Green);
#pragma omp parallel
         for (const auto &p : net->getPoints())
#pragma omp single nowait
         {
            p->RK_U.push(p->DUDt_SPH);  // 速度
            p->U_SPH = p->RK_U.getX();
            p->RK_X.push(p->U_SPH);  // 位置
            p->setXSingle(p->tmp_X = p->RK_X.getX());
            /* -------------------------------------------------------------------------- */
            int count = 0;
#if defined(REFLECTION)
            auto closest = [&]() {
               double distance = 1E+20;
               networkPoint *P = nullptr;
               for (const auto &[obj, poly] : RigidBodyObject) {
                  obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                     auto tmp = Distance(p->X, q);
                     if (distance > tmp) {
                        distance = tmp;
                        P = q;
                     }
                  });
               }
               return P;
            };
            bool isReflect = true;
            while (isReflect && count++ < 30) {
               // const auto X = p->RK_X.getX(p->U_SPH);
               isReflect = false;
               networkPoint *p_wall;
               if (p_wall = closest())
                  if (particle_spacing > Distance(p_wall->X, p->X)) {
                     auto normal_distance = Norm(Projection(p->X - p_wall->X, p_wall->normal_SPH));
                     if (Dot(p->U_SPH, p_wall->normal_SPH) < 0) {
                        p->DUDt_SPH -= (1. + damping_factor) * Projection(p->U_SPH, p_wall->normal_SPH) / dt;
                        p->RK_U.repush(p->DUDt_SPH);  // 速度
                        p->U_SPH = p->RK_U.getX();
                        p->RK_X.repush(p->U_SPH);  // 位置
                        p->setXSingle(p->tmp_X = p->RK_X.getX());
                        isReflect = true;
                     }
                  }
            };
#endif
         }
         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "粒子の時間発展");
         real_time = (*net->getPoints().begin())->RK_X.gett();
      } while (!((*net->getPoints().begin())->RK_X.finished));

      Print(Yellow, "Elapsed time: ", Red, watch(), "s １タイムステップ終了");

   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

void setDataOmitted(auto &vtp, const auto &Fluid) {
   std::unordered_map<networkPoint *, Tddd> normal_SPH;
   for (const auto &p : Fluid->getPoints())
      normal_SPH[p] = p->normal_SPH;
   vtp.addPointData("normal_SPH", normal_SPH);
   std::unordered_map<networkPoint *, Tddd> interpolated_normal_SPH;
   for (const auto &p : Fluid->getPoints())
      interpolated_normal_SPH[p] = p->interpolated_normal_SPH;
   vtp.addPointData("interpolated_normal_SPH", interpolated_normal_SPH);
   //
   std::unordered_map<networkPoint *, Tddd> interpolated_normal_SPH_original;
   interpolated_normal_SPH_original.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      interpolated_normal_SPH_original[p] = p->interpolated_normal_SPH_original;
   vtp.addPointData("interpolated_normal_SPH_original", interpolated_normal_SPH_original);
   //
   std::unordered_map<networkPoint *, Tddd> U;
   for (const auto &p : Fluid->getPoints())
      U[p] = p->U_SPH;
   vtp.addPointData("U", U);
   std::unordered_map<networkPoint *, double> pressure;
   for (const auto &p : Fluid->getPoints())
      pressure[p] = p->p_SPH;
   vtp.addPointData("pressure", pressure);
};

void setData(auto &vtp, const auto &Fluid) {
   std::unordered_map<networkPoint *, Tddd> U;
   for (const auto &p : Fluid->getPoints())
      U[p] = p->U_SPH;
   vtp.addPointData("U", U);
   //
   std::unordered_map<networkPoint *, double> density;
   density.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      density[p] = p->rho_;
   vtp.addPointData("temporal density", density);
   //
   std::unordered_map<networkPoint *, double> pressure;
   pressure.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      pressure[p] = p->p_SPH;
   vtp.addPointData("pressure", pressure);
   //
   // std::unordered_map<networkPoint *, double> contactpoints;
   // for (const auto &p : Fluid->getPoints())
   //    contactpoints[p] = (double)p->getContactPoints().size();
   // vtp.addPointData("contact points", contactpoints);
   // //
   std::unordered_map<networkPoint *, Tddd> interpolated_normal_SPH;
   interpolated_normal_SPH.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      interpolated_normal_SPH[p] = p->interpolated_normal_SPH;
   vtp.addPointData("interpolated_normal_SPH", interpolated_normal_SPH);
   //
   std::unordered_map<networkPoint *, Tddd> interpolated_normal_SPH_original;
   interpolated_normal_SPH_original.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      interpolated_normal_SPH_original[p] = p->interpolated_normal_SPH_original;
   vtp.addPointData("interpolated_normal_SPH_original", interpolated_normal_SPH_original);
   // //
   std::unordered_map<networkPoint *, Tddd> position;
   position.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      position[p] = p->X;
   vtp.addPointData("position", position);
   // //
   std::unordered_map<networkPoint *, Tddd> tmp_U_SPH;
   for (const auto &p : Fluid->getPoints())
      tmp_U_SPH[p] = p->tmp_U_SPH;
   vtp.addPointData("tmp_U_SPH", tmp_U_SPH);
   // //
   std::unordered_map<networkPoint *, double> div_U;
   for (const auto &p : Fluid->getPoints())
      div_U[p] = p->div_U;
   vtp.addPointData("div_U", div_U);
   // //
   std::unordered_map<networkPoint *, Tddd> gradP_SPH;
   for (const auto &p : Fluid->getPoints())
      gradP_SPH[p] = p->gradP_SPH / p->rho;
   vtp.addPointData("gradP_SPH / rho", gradP_SPH);
   // //
   std::unordered_map<networkPoint *, Tddd> lap_U;
   for (const auto &p : Fluid->getPoints())
      lap_U[p] = p->lap_U;
   vtp.addPointData("lap_U", lap_U);
   // //
   std::unordered_map<networkPoint *, Tddd> DUDt;
   for (const auto &p : Fluid->getPoints())
      DUDt[p] = p->DUDt_SPH;
   vtp.addPointData("DUDt", DUDt);
   // //
   // std::unordered_map<networkPoint *, Tddd> whereToReference;
   // for (const auto &p : Fluid->getPoints())
   //    whereToReference[p] = ToX(p) + 2 * p->normal_SPH - ToX(p);
   // vtp.addPointData("where to reference", whereToReference);
   // //
   std::unordered_map<networkPoint *, double> isSurface;
   isSurface.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      isSurface[p] = p->isSurface;
   vtp.addPointData("isSurface", isSurface);
   // //
   // std::unordered_map<networkPoint *, double> isInsideOfBody;
   // for (const auto &p : Fluid->getPoints())
   //    isInsideOfBody[p] = p->isInsideOfBody;
   // vtp.addPointData("isInsideOfBody", isInsideOfBody);
   // //
   // std::unordered_map<networkPoint *, double> isCaptured;
   // for (const auto &p : Fluid->getPoints())
   //    isCaptured[p] = p->isCaptured ? 1. : -1E+50;
   // vtp.addPointData("isCaptured", isCaptured);
   // //
   // std::unordered_map<networkPoint *, Tddd> repulsive_force_SPH;
   // for (const auto &p : Fluid->getPoints())
   //    repulsive_force_SPH[p] = p->repulsive_force_SPH;
   // vtp.addPointData("repulsive_force_SPH", repulsive_force_SPH);
   // //
   std::unordered_map<networkPoint *, double> contact_points;
   for (const auto &p : Fluid->getPoints())
      contact_points[p] = p->contact_points_fluid_SPH;
   vtp.addPointData("contact_points", contact_points);
   std::unordered_map<networkPoint *, double> contact_points_inradius;
   for (const auto &p : Fluid->getPoints())
      contact_points_inradius[p] = p->contact_points_all_SPH;
   vtp.addPointData("contact_points_inradius", contact_points_inradius);
};

#endif
//