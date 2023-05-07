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

// b$ ------------------------------------------------------ */
// b$                    ∇.∇UとU*を計算                       */
// b$ ------------------------------------------------------ */

auto Lap_U(const netPp A, const std::unordered_set<Network *> &target_nets) {
   // ここで初期化することは問題ない．この点pに対して再度計算することはないので，途中で初期化される心配はない．
   A->checked_points_in_radius_SPH = A->checked_points_in_radius_of_fluid_SPH = A->checked_points_SPH = 0;
   A->lap_U_ = A->lap_U = {0., 0., 0.};
   //$ ------------------------------------------ */
   auto func = [&](const auto &B, const auto &qX, const double coef = 1.) {
      const auto rij = qX - A->X;
      if (Between(Norm(rij), {1E-10, A->radius_SPH})) {
         const auto Uij = coef * B->U_SPH - A->U_SPH;
         const auto nu_nu = B->mu_SPH / B->rho + A->mu_SPH / A->rho;
         A->lap_U_ += 1 / (A->mu_SPH / A->rho) * B->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(A->X, qX, A->radius_SPH) /
                      ((B->rho + A->rho) * Dot(rij, rij));

         A->checked_points_in_radius_SPH++;
         if (B->getNetwork()->isFluid || B->isFluid)
            A->checked_points_in_radius_of_fluid_SPH++;
      }
      A->checked_points_SPH++;
   };
   //$ ------------------------------------------ */
   //$ ------------------------------------------ */
   // auto func = [&](const auto &B, const Tddd &qX) {
   //    const auto rij = A->X - qX;
   //    const auto Uij = A->U_SPH - coef *B->U_SPH;
   //    A->lap_U_ += 2 * B->mass / A->rho * Dot(rij, grad_w_Bspline(A->X, qX, A->radius_SPH)) / Dot(rij, rij) * Uij;
   // };
   //$ ------------------------------------------ */
   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
         if (B->isCaptured) {
            if (A != B)
               func(B, B->X);
#ifdef USE_SPP_Fluid
            if (B->isSurface) {
               if (canSetSPP(target_nets, B)) {
                  func(B, SPP_X(B), SPP_U_coef);
                  // func(B, SPP_X(B, 2.), SPP_U_coef);
               }
            }
#endif
         }
      });
};

auto Lap_U(const auto &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
      Lap_U(A, target_nets);
};

auto setLap_U(const auto &points, const double dt) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->lap_U = A->lap_U_;
      A->ViscousAndGravityForce = A->DUDt_SPH = A->DUDt_SPH_ = A->lap_U_ * (A->mu_SPH / A->rho) + _GRAVITY3_;  // 後で修正されるDUDt
      A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
      A->tmp_X = A->X + A->tmp_U_SPH * dt;
   }
};

// b$ ------------------------------------------------------ */
// b$      div(U^*),　DρDt=-ρdiv(U),　ラプラシアン∇.∇U* の計算      */
// b$ ------------------------------------------------------ */

void div_tmpU(const netPp A, const std::unordered_set<Network *> &target_nets) {
   A->div_tmpU = A->div_U = 0;  // this is div(U^*) not div(U^n)
   A->lap_tmpU = A->grad_div_U = {0., 0., 0.};
   double nu_nu;
   //@ ------------------------------------------ */
   // to use both A and spp
   auto func = [&](const auto &B, const auto &qX, const double coef = 1.) {
      if (Between(Distance(A, qX), {1E-10, A->radius_SPH})) {
         // 後藤p.25 (2.89)
         auto Uij = coef * B->U_SPH - A->U_SPH;
         A->div_U += B->mass / A->rho * Dot(Uij, grad_w_Bspline(A->X, qX, A->radius_SPH));
         //
         Uij = coef * B->tmp_U_SPH - A->tmp_U_SPH;
         A->div_tmpU += B->mass / A->rho * Dot(Uij, grad_w_Bspline(A->X, qX, A->radius_SPH));
         //
         auto rij = qX - A->X;
         nu_nu = B->mu_SPH / B->rho + A->mu_SPH / A->rho;
         A->lap_tmpU += 1 / (A->mu_SPH / A->rho) * B->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(A->X, qX, A->radius_SPH) /
                        ((B->rho + A->rho) * Dot(rij, rij));
      }
   };
   //@ ------------------------------------------ */
   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
         if (B->isCaptured) {
            func(B, B->X);
#ifdef USE_SPP_Fluid
            if (B->isSurface) {
               if (canSetSPP(target_nets, B)) {
                  func(B, SPP_X(B), SPP_U_coef);
                  // func(B, SPP_X(B, 2.), SPP_U_coef);
               }
            }
#endif
         }
      });
};

void div_tmpU(const std::unordered_set<networkPoint *> &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
      div_tmpU(p, target_nets);
};

void div_tmpU(const auto &pointsA, const auto &pointsB, const std::unordered_set<Network *> &target_nets) {
   div_tmpU(pointsA, target_nets);
   div_tmpU(pointsB, target_nets);
};

auto set_tmpLap_U(const auto &points) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
   {
      p->tmp_ViscousAndGravityForce = p->lap_tmpU * (p->mu_SPH / p->rho) + _GRAVITY3_;  // 後で修正されるDUDt
   }
};

// b% -------------------------------------------------------------------------- */
// b% -------------------------------------------------------------------------- */
// #define Morikawa2019
// #define NewtonMethod

void calculateTemporalPressure(const netPp A, const std::unordered_set<Network *> &target_nets, const double dt) {
   double sum_Aij = 0, Aij, sum_Aij_Pj = 0, sum_Aij_Pij = 0, sum_Pij = 0, pressure_SPH = 0, qP;
   Tddd Xij;
   A->p_SPH_ = 0;
   //% ------------------------------------------ */
   auto func = [&](const auto &B, const auto &qX, const double coef = 1.) {
      Xij = A->X - qX;
      qP = coef * B->p_SPH;
      if (Between(Distance(A, qX), {1E-10, A->radius_SPH})) {
         auto share = Dot(Xij, grad_w_Bspline(A->X, qX, A->radius_SPH)) / Dot(Xij, Xij);
#if defined(Morikawa2019)
         // Morikawa, D., Senadheera, H., & Asai, M. (2021). Explicit incompressible smoothed particle hydrodynamics in a multi-GPU environment for large-scale simulations. Computational Particle Mechanics, 8(3), 493–510. https://doi.org/10.1007/s40571-020-00347-0
         Aij = 2. * B->mass / A->rho * share;
#else
         // Nomeritae
         // Shao and Lo
         // auto Aij = B->mass * 8. / std::pow(B->rho + A->rho, 2) * share;, which is the same as
         Aij = 2. * B->mass / std::pow((B->rho + A->rho) / 2., 2) * share;
#endif
         sum_Aij += Aij;
         sum_Pij += A->p_SPH - qP;
         sum_Aij_Pj += Aij * qP;
         sum_Aij_Pij += Aij * (A->p_SPH - qP);
      }
   };
   //% ------------------------------------------ */
   auto FUNC = [&](const auto &B) {
      if (B->isCaptured) {
         func(B, B->X);
#ifdef USE_SPP_Fluid
         if (B->isSurface)
            if (canSetSPP(target_nets, B)) {
               func(B, SPP_X(B), SPP_p_coef);
               // func(B, SPP_X(B, 2.), SPP_p_coef);
            }
#endif
      }
   };
   //% -------------------------------------------------------------------------- */
   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH, FUNC);
      // for (const auto &[obj, poly] : RigidBodyObject)
      //    obj->BucketPoints.apply(A->X, A->radius_SPH, FUNC);

#if defined(Morikawa2019)
   auto b = _WATER_DENSITY_ * A->div_tmpU / dt;
#else
   auto b = A->div_tmpU / dt;
#endif
   A->p_SPH_ = (b + sum_Aij_Pj) / sum_Aij;

#if defined(NewtonMethod)
   A->div_U_error = std::abs(sum_Aij_Pij - b);
   auto F = sum_Aij_Pij - b;
   auto dFdPi = sum_Aij;

   #ifdef surface_zero_pressure
   if (A->isSurface) {
      F += A->p_SPH;
      sum_Aij += 1;
   }
   #endif

   auto FF = F * F / 2.;
   auto dFFdPi = F * dFdPi;
   A->NR_pressure.update(FF, dFFdPi, 2.);

#endif
};

void calculateTemporalPressure(const std::unordered_set<networkPoint *> &points,
                               const std::unordered_set<Network *> &target_nets,
                               const double dt) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
      calculateTemporalPressure(p, target_nets, dt);
};

void setPressure(const std::unordered_set<networkPoint *> &points) {
   for (const auto &p : points)
      p->p_SPH = p->p_SPH_;
}

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */

void gradP(const netPp A, const std::unordered_set<Network *> &target_nets) {
   // A->checked_points_in_radius_SPH = A->checked_points_in_radius_of_fluid_SPH = 0;
   A->gradP_SPH.fill(0.);
   //% ------------------------------------------ */
   auto func = [&](const auto &B, const auto &qX, const double coef = 1.) {
      if (Between(Distance(A, qX), {1E-10, A->radius_SPH})) {
         // double p_rho2 = coef * B->p_SPH / std::pow(B->rho, 2) + A->p_SPH / std::pow(A->rho, 2);
         // auto qP = (spp ? SPP_p_coef * B->p_SPH_SPP : B->p_SPH);
         // A->gradP_SPH += (qP + A->p_SPH) * B->mass / _WATER_DENSITY_ * grad_w_Bspline(A->X, qX, A->radius_SPH);
         A->gradP_SPH += A->rho * (coef * B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * B->mass * grad_w_Bspline(A->X, qX, A->radius_SPH);
      }
   };
   //% ------------------------------------------ */
   auto FUNC = [&](const auto &B) {
      if (B->isCaptured) {
         func(B, B->X);
#ifdef USE_SPP_Fluid
         if (B->isSurface) {
            if (canSetSPP(target_nets, B)) {
               func(B, SPP_X(B), SPP_p_coef);
               // func(B, SPP_X(B, 2.), SPP_p_coef);
            }
         }
#endif
      }
   };

   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH, FUNC);

   // A->DUDt_SPH -= A->gradP_SPH / A->rho;
   A->DUDt_SPH -= A->gradP_SPH / _WATER_DENSITY_;
   //
   if (!isFinite(A->DUDt_SPH))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
};

void gradP(const std::unordered_set<networkPoint *> &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
      gradP(p, target_nets);
}
// b$ ------------------------------------------------------ */
// b$                    ∇²Pを計算                           */
// b$ ------------------------------------------------------ */

auto Lap_P(const netPp A, const std::unordered_set<Network *> &target_nets) {
   A->column_value.clear();
   auto func = [&](const auto &B, const auto &qX, const double coef = 1.) {
      const auto rij = qX - A->X;
      if (Between(Distance(A, qX), {1E-10, A->radius_SPH})) {
         double v = 2 * B->mass * Dot(rij, grad_w_Bspline(A->X, qX, A->radius_SPH));
         v /= B->rho * Dot(rij, rij);
         A->increment(B, v);
         A->increment(A, -v);
      }
   };

   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
         if (B->isCaptured) {
            if (A != B)
               func(B, B->X);
#ifdef USE_SPP_Fluid
            if (B->isSurface) {
               if (canSetSPP(target_nets, B)) {
                  func(B, SPP_X(B), SPP_p_coef);
                  // func(B, SPP_X(B, 2.), SPP_p_coef);
               }
            }
#endif
         }
      });
};

auto Lap_P(const auto &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
      Lap_P(A, target_nets);
};

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