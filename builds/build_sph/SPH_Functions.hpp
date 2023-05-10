#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"
/* -------------------------------------------------------------------------- */

double dt_CFL(const double dt_IN, const auto &net, const auto &RigidBodyObject) {
   double dt = dt_IN;
   const auto C_CFL_velocity = 0.02;  // dt = C_CFL_velocity*h/Max(U)
   const auto C_CFL_accel = 0.1;      // dt = C_CFL_accel*sqrt(h/Max(A))
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
   return dt;
}
/* -------------------------------------------------------------------------- */

Tddd SPP_X(const networkPoint *p, const double c = 1.0) {
   // return p - (p->COM_SPH - p);
   // return p->X - 2 * p->COM_SPH;
   return p->X + c * p->interpolated_normal_SPH * p->radius_SPH / p->C_SML;
};

Tddd SPP_X_tmp(const networkPoint *p, const double c = 1.0) {
   return p->tmp_X + c * p->interpolated_normal_SPH * p->radius_SPH / p->C_SML;
};

auto canSetSPP(const auto &target_nets, const auto &p) {

#if defined(USE_SPP_Fluid)
   auto X = SPP_X(p);
   auto range = p->radius_SPH / p->C_SML;
   if (!p->isSurface)
      return false;

   // auto d = Distance(p, X);

   // auto func = [&](const auto &q) { return p != q && Distance(q, X) < d; };

   // if (net->BucketPoints.any_of(X, d, func))
   //    return false;
   // for (const auto &[obj, poly] : RigidBodyObject)
   //    if (obj->BucketPoints.any_of(X, d, func))
   //       return false;

   const double C = 1.2;
   for (const auto &net : target_nets)
      if (net->BucketPoints.any_of(X, C * range,
                                   [&](const auto &q) { return p != q && Distance(q, X) < C * range; }))
         return false;

   // {
   //    const double C = 1.2;
   //    auto func = [&](const auto &q) { return p != q && Distance(q, X) < C * range; };
   //    if (net->BucketPoints.any_of(X, C * range, func))
   //       return false;
   // }
   // {
   //    const double C = 1.2;
   //    auto func = [&](const auto &q) { return p != q && Distance(q, X) < C * range; };
   //    for (const auto &[obj, poly] : RigidBodyObject)
   //       if (obj->BucketPoints.any_of(X, C * range, func))
   //          return false;
   // }
   return true;
#else
   return false;
#endif
};
/* -------------------------------------------------------------------------- */
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
      p->COM_SPH.fill(0.);
      p->interpolated_normal_SPH_original.fill(0.);
      p->interpolated_normal_SPH_original_all.fill(0.);
      double total_vol = 0, w;

      net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
         p->COM_SPH += (q->X - p->X) * w;
         total_vol += w;
         // if (Between(Distance(p, q), {1E-10, p->radius_SPH}))
         {
            p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            p->interpolated_normal_SPH_original_all -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         }
      });

      std::vector<Tddd> wall_tangent;
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
            p->COM_SPH += (q->X - p->X) * w;
            total_vol += w;
            if (Distance(p, q) < p->radius_SPH / p->C_SML * 1.5)
               wall_tangent.emplace_back(q->normal_SPH);
            p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            p->interpolated_normal_SPH_original_all -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         });

      p->COM_SPH /= total_vol;
      p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
      p->interpolated_normal_SPH_all = Normalize(p->interpolated_normal_SPH_original_all);

      for (const auto &wall_v : wall_tangent) {
         p->interpolated_normal_SPH = Normalize(Chop(p->interpolated_normal_SPH, wall_v));
         p->interpolated_normal_SPH_all = Normalize(Chop(p->interpolated_normal_SPH_all, wall_v));
      }

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
                      return p != q && (VectorAngle(p->interpolated_normal_SPH_all, q->X - p->X) < std::numbers::pi / 4);
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
                            return p != q && (VectorAngle(p->interpolated_normal_SPH_all, -q->normal_SPH) < std::numbers::pi / 180. * 60);
                         } else
                            return false;
                      }
                      return false;
                   }))
                  p->isSurface = false;
#ifdef surface_zero_pressure
         if (p->isSurface)
            p->p_SPH = 0;
#endif
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
               // if (Between(Distance(p, q), {1E-8, p->radius_SPH}))
               p->interpolated_normal_SPH_original -= grad_w_Bspline(p->X, q->X, p->radius_SPH);
            });
         p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
         if (!isFinite(p->interpolated_normal_SPH))
            p->interpolated_normal_SPH = {0., 0., 1.};
      }
   }
};

/*DOC_EXTRACT
## 壁面粒子の流速と圧力
壁面粒子の流速は常にゼロとすることは自然なこと．常にゼロとするならば，壁面粒子の流速をマップする方法に悩む必要はない．
一方，壁面粒子の圧力は，各ステップ毎に計算し直す必要がある．

壁面粒子の圧力は，壁面法線方向流速をゼロにするように設定されるべきだろう．


*/

void mapValueOnWall(auto &net,
                    const std::unordered_set<networkPoint *> &wall_p,
                    const auto &RigidBodyObject,
                    const Tiii &indicater = {1, 0, 1}) {
   setPressureSPP(net, RigidBodyObject);
   auto [first, second, do_add] = indicater;
   Timer watch;
   auto initialize = [](networkPoint *PW) {
      PW->div_U_ = PW->volume_ = PW->rho_ = PW->p_SPH_ = PW->total_weight = PW->total_weight_ = 0;
      PW->U_SPH_ = PW->tmp_U_SPH_ = PW->grad_div_U_ = {0., 0., 0.};
      PW->DUDt_SPH_ = PW->gradP_SPH_ = {0., 0., 0.};
      PW->tmp_X = PW->X;
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
   // b$ ---------------------------------- calc ---------------------------------- */
   auto calc = [&]() {
#pragma omp parallel
      for (const auto &PW : wall_p)
#pragma omp single nowait
      {
         const Tddd markerX = PW->X + 2 * PW->normal_SPH;
         double W;
         int count = 0;
         // initialize(PW);
         /* -------------------------------------------------------------------------- */
         auto func = [&](const networkPoint *PF, const Tddd &nearX, const bool spp = false) {
#ifdef USE_SPP_Wall
            double cU = spp ? SPP_U_coef_of_Wall : 1.;
            double cp = spp ? SPP_p_coef_of_Wall : 1.;
            double cRho = spp ? SPP_rho_coef_of_Wall : 1.;
#else
            double cU = 1.;
            double cp = 1.;
            double cRho = 1.;
#endif
            // if (PF->getNetwork()->isFluid)
            if (PF->isCaptured) {
               W = w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a);
               PW->total_N += W;
               PW->total_weight += W;
               PW->total_weight_ += W;
               PW->U_SPH_ += cU * PF->U_SPH * W;
               PW->tmp_U_SPH_ += cU * PF->tmp_U_SPH * W;
               PW->grad_div_U_ += PF->grad_div_U * W;
               PW->p_SPH_ += cp * (spp ? PF->p_SPH : PF->p_SPH) * W;
               PW->div_U_ += cU * PF->div_U * W;
               PW->rho_ += cRho * PF->rho * W;
               PW->volume_ += PF->volume * W;
               PW->DUDt_SPH_ += PF->DUDt_SPH * W;
               PW->gradP_SPH_ += PF->gradP_SPH * W;
            }
            count++;
         };
         /* -------------------------------------------------------------------------- */
         net->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &PF) {
            func(PF, PF->X);
#ifdef USE_SPP_Wall
            if (do_spp && PF->isSurface) {
               // if (canSetSPP(net, RigidBodyObject, PF))
               func(PF, SPP_X(PF), true);
               // func(PF, SPP_X(PF, 2.), true);
            }
#endif
         });
         for (const auto &[obj, _] : RigidBodyObject)
            obj->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &PF) { func(PF, PF->X); });
      }
   };
   // b% ----------------------------------- set ---------------------------------- */
   auto set = [&]() {
      for (const auto &PW : wall_p) {
         {
            PW->p_SPH = PW->p_SPH_ / PW->total_weight;
            PW->U_SPH = PW->U_SPH_ / PW->total_weight;
            PW->tmp_U_SPH = PW->tmp_U_SPH_ / PW->total_weight;
            PW->grad_div_U = PW->grad_div_U_ / PW->total_weight;
            PW->div_U = PW->div_U_ / PW->total_weight;
            PW->DUDt_SPH = PW->DUDt_SPH_ / PW->total_weight;
         }
         std::vector<bool> Bools = {!isFinite(PW->p_SPH), !isFinite(PW->rho), !isFinite(PW->volume), !isFinite(PW->U_SPH), !isFinite(PW->tmp_U_SPH)};
         if (std::ranges::any_of(Bools, [](bool b) { return b; })) {
            std::stringstream ss;
            ss << Bools << std::endl;
            ss << "PW->total_weight = " << PW->total_weight << std::endl;
            ss << "PW->p_SPH = " << PW->p_SPH << std::endl;
            ss << "PW->isFluid = " << PW->isFluid << std::endl;
            ss << "PW->rho = " << PW->rho << std::endl;
            ss << "PW->volume = " << PW->volume << std::endl;
            // ss << "PW->lap_tmpU = " << PW->lap_tmpU << std::endl;
            ss << "PW->U_SPH = " << PW->U_SPH << std::endl;
            ss << "PW->tmp_U_SPH = " << PW->tmp_U_SPH << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite" + ss.str());
         }
      }
   };
   // b% -------------------------------------------------------------------------- */
   POW = 1.;  // ２などとしてはいけない
   a = 1.;
   initializeAll();
   applyToWall = false;
   mirroring = false;
   account_volume = true;
   do_spp = false;
   // for (auto i = 0; i < 2; i++) {
   calc();
   set();
   // }
   setNormal_Surface_(net, wall_p, RigidBodyObject);
   //! ------------------------------- apply condition! ------------------------------ */
   for (const auto &PW : wall_p) {

      // free-slip
      // PW->U_SPH = Reflect(PW->U_SPH, PW->normal_SPH);
      // PW->tmp_U_SPH = Reflect(PW->tmp_U_SPH, PW->normal_SPH);

      // no-slip
      // PW->U_SPH *= -1.;
      // PW->tmp_U_SPH *= -1.;

      PW->U_SPH *= 0.;
      PW->tmp_U_SPH *= 0.;

      // if (Norm(PW->normal_SPH) < 1E-12) {
      //    PW->U_SPH *= 0;
      //    PW->tmp_U_SPH *= 0;
      // }

      // ゼロとした方が，悪い影響を受けないのではないだろうか？
      // PW->lap_tmpU_ = Projection(PW->lap_tmpU_, PW->normal_SPH);
      // PW->U_SPH = Projection(PW->U_SPH, PW->normal_SPH);
      // PW->tmp_U_SPH = Projection(PW->tmp_U_SPH, PW->normal_SPH);

      //
      // all-slip
      // 壁面上に粒子を設定した場合，壁上で流速とDUDtをゼロとするのは，妥当な設定
      // PW->lap_tmpU *= 0.;これはおかしい，ここでのDUDtは粘性と重力のみを考慮したものにすぎない．
      // PW->U_SPH *= 0.;
      // PW->tmp_U_SPH *= 0.;

      if (do_add) {
         Tddd accel = {0., 0., 0.};
         PW->p_SPH = PW->p_SPH + _WATER_DENSITY_ * Dot(PW->ViscousAndGravityForce - accel, -2 * PW->normal_SPH);
         // PW->p_SPH = PW->p_SPH + _WATER_DENSITY_ * Dot(PW->tmp_ViscousAndGravityForce - accel, -2 * PW->normal_SPH);
         // PW->p_SPH = PW->p_SPH + _WATER_DENSITY_ * Dot(PW->mu_SPH / PW->rho * PW->lap_tmpU - PW->grad_div_U /*/ dtは入れ込み済み*/ + _GRAVITY3_ - accel, -2 * PW->normal_SPH);
         // PW->p_SPH = PW->p_SPH + _WATER_DENSITY_ * Dot(PW->DUDt_SPH - accel, -2 * PW->normal_SPH);
         // auto nu = PW->mu_SPH / PW->rho;
         // PW->p_SPH = PW->p_SPH + _WATER_DENSITY_ * Dot(nu * PW->lap_tmpU + _GRAVITY3_ - accel, -2 * PW->normal_SPH);
      }
   }
   DebugPrint(_green, "Elapsed time: ", _red, watch(), "s ", Magenta, "壁面粒子の圧力を計算", (do_add ? ", DUDtを考慮した圧力" : ""));
};

/* -------------------------------------------------------------------------- */

#define Morikawa2019
#define new_method

// b$ ------------------------------------------------------ */
// b$                    ∇.∇UとU*を計算                       */
// b$ ------------------------------------------------------ */

auto Lap_U(const auto &points, const std::unordered_set<Network *> &target_nets, const double dt) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      // ここで初期化することは問題ない．この点pに対して再度計算することはないので，途中で初期化される心配はない．
      A->checked_points_in_radius_SPH = A->checked_points_in_radius_of_fluid_SPH = A->checked_points_SPH = 0;
      A->lap_U.fill(0.);
      //$ ------------------------------------------ */
      auto add = [&](const auto &B, const auto &qX, const double coef = 1.) {
         const auto rij = qX - A->X;
         if (Between(Norm(rij), {1E-10, A->radius_SPH})) {

#if defined(Morikawa2019)
            const auto Uij = A->U_SPH - coef * B->U_SPH;
            A->lap_U += 2 * B->mass / A->rho * Dot_grad_w_Bspline_Dot(A->X, qX, A->radius_SPH) * Uij;
#elif defined(Nomeritae2016)
            const auto Uij = coef * B->U_SPH - A->U_SPH;
            const auto nu_nu = B->mu_SPH / B->rho + A->mu_SPH / A->rho;
            A->lap_U += 1 / (A->mu_SPH / A->rho) * B->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(A->X, qX, A->radius_SPH) /
                        ((B->rho + A->rho) * Dot(rij, rij));
#endif
            // just counting
            A->checked_points_in_radius_SPH++;
            if (B->getNetwork()->isFluid || B->isFluid)
               A->checked_points_in_radius_of_fluid_SPH++;
         }
         A->checked_points_SPH++;
      };
      //$ ------------------------------------------ */
      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               add(B, B->X);
#ifdef USE_SPP_Fluid
               if (B->isSurface && canSetSPP(target_nets, B)) add(B, SPP_X(B), SPP_U_coef);
#endif
            }
         });
      //$ ------------------------------------------ */
      A->DUDt_SPH_ = A->lap_U * (A->mu_SPH / A->rho) + _GRAVITY3_;  // 後で修正されるDUDt
      A->ViscousAndGravityForce = A->DUDt_SPH = A->DUDt_SPH_;
      A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
      A->tmp_X = A->X + A->tmp_U_SPH * dt;
   }
};
// b$ ------------------------------------------------------ */
// b$                     rho^*  の計算                       */
// b$ ------------------------------------------------------ */

void setTmpDensity(const std::unordered_set<networkPoint *> &points,
                   const double dt) {
#if defined(Morikawa2019)
      /* use rho_ calculated at div_tmpU */
#elif defined(Nomeritae2016)
   for (const auto &p : points)
      p->rho_ = p->rho + (p->DrhoDt_SPH = -p->rho * p->div_tmpU) * dt;
#elif defined(Barcarolo2013)
   throw std::runtime_error("not implemented");
#endif
}

// b% -------------------------------------------------------------------------- */

/*DOC_EXTRACT
### `PoissonRHS`と $\nabla^2 p^{n+1}$における $p^{n+1}$の係数の計算

$$
\begin{align*}
\frac{D {\bf u}}{D t} &=-\frac{1}{\rho} \nabla P+\nu \nabla^2 {\bf u}+{\bf g}\\
\rightarrow \nabla \cdot\left(\frac{\rho}{\Delta t} {\bf u}^{n+1}\right) + \nabla^2 p &= \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}+\rho {\bf g}\right)\\
\rightarrow \nabla^2 p &= b, \quad b = \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}+\rho {\bf g}\right)
\end{align*}
$$

ここの $b$を`PoissonRHS`とする．

CHECKED: $\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}$

CHECKED: $\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}$

CHECKED: $\nabla^2 p^{n+1} = \sum_{j}A_{ij}(p_i^{n+1} - p_j^{n+1})$
$A_{ij}=\frac{2}{{\rho}_i}m_j$
$\frac{{\bf x}_{ij}\cdot \nabla W_{ij}}{{\bf x}_{ij}^2}$
*/

void PoissonEquation(const std::unordered_set<networkPoint *> &points,
                     const std::unordered_set<Network *> &target_nets,
                     const double dt,
                     const bool isWall = false) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
      A->PoissonRHS = 0;
      std::array<double, 3> B_VALUE;
      const auto A_VALUE = (A->rho / dt * A->U_SPH) + (A->mu_SPH * A->lap_U) + (A->rho * _GRAVITY3_);
      //
      A->column_value.clear();
      auto markerX = A->X;
      double total_weight = 0, P = 0;
      if (isWall)
         markerX += (1. - 1E-12) * A->normal_SPH;
      //
      //% ----------------- PoissonRHS ------------------------- */
      auto add = [&](const auto &B, const auto &qX, const double coef = 1.) {
         B_VALUE = (B->rho / dt * B->U_SPH) + (B->mu_SPH * B->lap_U) + (B->rho * _GRAVITY3_);
         A->PoissonRHS += B->volume * Dot(B_VALUE - A_VALUE, grad_w_Bspline(markerX, qX, A->radius_SPH));
#if defined(Morikawa2019)
         Aij = 2. * B->mass / _WATER_DENSITY_ * Dot_grad_w_Bspline_Dot(markerX, qX, A->radius_SPH);
#elif defined(Nomeritae2016)
         Aij = 2. * B->mass / std::pow((B->rho_ + A->rho_) / 2., 2.) * Dot_grad_w_Bspline_Dot(markerX, qX, A->radius_SPH) * A->rho_;
#elif defined(Barcarolo2013)
         Aij = 2. * B->mass / _WATER_DENSITY_ * Dot_grad_w_Bspline_Dot(markerX, qX, A->radius_SPH);
#endif
         sum_Aij += Aij;
         sum_Aij_Pj += Aij * B->p_SPH;

         A->increment(B, -Aij);
         A->increment(A, Aij);

         // for mapping to wall
         total_weight += B->volume * w_Bspline(Norm(markerX - qX), A->radius_SPH);
         P += B->p_SPH * B->volume * w_Bspline(Norm(markerX - qX), A->radius_SPH);
      };

      //% ------------------------------------------------------- */

      for (const auto &net : target_nets)
         net->BucketPoints.apply(markerX, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               add(B, B->X);
#ifdef USE_SPP_Fluid
               if (B->isSurface && canSetSPP(target_nets, B)) add(B, SPP_X(B), SPP_p_coef);
#endif
            }
         });

      A->p_SPH_ = (A->PoissonRHS + sum_Aij_Pj) / sum_Aij;
      if (isWall) {
         if (total_weight > 1E-10)
            A->p_SPH_ = P / total_weight;
         else
            A->p_SPH_ = 0;
      }
   };
};

void setPressure(const std::unordered_set<networkPoint *> &points) {
   for (const auto &p : points)
      p->p_SPH = p->p_SPH_;
}

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */
/*DOC_EXTRACT

CHECKED: $\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}$

CHECKED: $\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}$

*/
void gradP(const std::unordered_set<networkPoint *> &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->gradP_SPH.fill(0.);
      //% ------------------------------------------ */
      auto func = [&](const auto &B, const auto &qX, const double coef = 1.) {
#if defined(Morikawa2019)
         A->gradP_SPH += A->rho * B->mass * (coef * B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, qX, A->radius_SPH);
#elif defined(Nomeritae2016)
         A->gradP_SPH += B->volume * B->p_SPH * grad_w_Bspline(A->X, qX, A->radius_SPH);
#elif defined(Barcarolo2013)
         A->gradP_SPH += A->rho * B->mass * (coef * B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, qX, A->radius_SPH);
#endif
      };
      //% ------------------------------------------ */
      auto FUNC = [&](const auto &B) {
         if (B->isCaptured) {
            func(B, B->X);
#ifdef USE_SPP_Fluid
            if (B->isSurface && canSetSPP(target_nets, B)) func(B, SPP_X(B), SPP_p_coef);
#endif
         }
      };

      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, FUNC);

#if defined(Morikawa2019)
      A->DUDt_SPH -= A->gradP_SPH / _WATER_DENSITY_;
#elif defined(Nomeritae2016)
      A->DUDt_SPH -= A->gradP_SPH / A->rho;
#endif

      if (!isFinite(A->DUDt_SPH))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
   }
}

//@ -------------------------------------------------------- */
//@                        粒子の時間発展                      */
//@ -------------------------------------------------------- */

void updateParticles(const auto &points, const auto &RigidBodyObject, const double &particle_spacing, const double dt) {
   DebugPrint("粒子の時間発展", Green);
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
   {
      // テスト
      auto U = p->U_SPH;
      auto X_last = p->X;
      p->RK_U.push(p->DUDt_SPH);  // 速度
      p->U_SPH = p->RK_U.getX();  // * 0.5 + U * 0.5;
      p->RK_X.push(p->U_SPH);     // 位置
      p->setXSingle(p->tmp_X = p->RK_X.getX());
      // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
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
      bool isReflected = true;
      while (isReflected && count++ < 30) {
         // const auto X = p->RK_X.getX(p->U_SPH);
         isReflected = false;
         networkPoint *closest_wall_point;
         if (closest_wall_point = closest()) {
            auto ovre_run = ((1. - asobi) * particle_spacing - Distance(closest_wall_point->X, p->X)) / 2.;
            if (ovre_run > 0.) {
               auto normal_distance = Norm(Projection(p->X - closest_wall_point->X, closest_wall_point->normal_SPH));
               if (Dot(p->U_SPH, closest_wall_point->normal_SPH) < 0) {
                  p->DUDt_SPH -= (1. + reflection_factor) * Projection(p->U_SPH, closest_wall_point->normal_SPH) / dt;
                  p->RK_U.repush(p->DUDt_SPH);  // 速度
                  p->U_SPH = p->RK_U.getX();    //* 0.5 + U * 0.5;
                  p->RK_X.repush(p->U_SPH);     // 位置
                  p->setXSingle(p->tmp_X = p->RK_X.getX());
                  //

                  // p->DUDt_SPH += (ovre_run * closest_wall_point->normal_SPH) / dt / dt;
                  // p->RK_U.repush(p->DUDt_SPH);  // 速度
                  // p->U_SPH = p->RK_U.getX();    //* 0.5 + U * 0.5;
                  // p->RK_X.repush(p->U_SPH);     // 位置
                  // p->setXSingle(p->tmp_X = p->RK_X.getX());

                  isReflected = true;
               }
            }
         }
      };
#endif
   }
}

#endif