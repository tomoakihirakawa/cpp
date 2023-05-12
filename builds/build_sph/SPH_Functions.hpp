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

/*DOC_EXTRACT SPH
### 壁面粒子の流速と圧力

壁粒子の流速を流体粒子の流速に応じて変化させると計算が煩雑になるので，**ここでは**壁面粒子の流速は常にゼロに設定することにした（ゼロで一定というのは不自然ではない）．
一方，壁粒子の圧力がゼロだとするのは不自然で，流体粒子の圧力$p^{n+1}$の計算に悪影響を及ぼす．
なので．壁粒子の圧力は各ステップ毎に計算し直す必要がある．

壁面粒子の圧力は，壁面法線方向流速をゼロにするように設定されるべきだろう．

*/

#define Morikawa2019
#define new_method

/*DOC_EXTRACT SPH
### $\nabla^2 {\bf u}$の計算

ラプラシアンの計算方法：

CHECKED: $\nabla^2 {\bf u}=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}$

*/

// b$ ------------------------------------------------------ */
// b$                    ∇.∇UとU*を計算                       */
// b$ ------------------------------------------------------ */

auto calcLaplacianU(const auto &points, const std::unordered_set<Network *> &target_nets, const double dt) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->checked_points_in_radius_SPH = A->checked_points_in_radius_of_fluid_SPH = A->checked_points_SPH = 0;
      A->lap_U.fill(0.);
      //$ ------------------------------------------ */
      auto add = [&](const auto &B, const auto &qX, const double coef = 1.) {
         const auto rij = qX - A->X;
         if (Between(Norm(rij), {1E-12, A->radius_SPH})) {
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

void setTmpDensity(const std::unordered_set<networkPoint *> &points, const double dt) {
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

/*DOC_EXTRACT SPH
### `PoissonRHS`,$b$と$\nabla^2 p^{n+1}$における$p^{n+1}$の係数の計算

$$
\begin{align*}
&&\frac{D {\bf u}}{D t} &=-\frac{1}{\rho} \nabla P+\nu \nabla^2 {\bf u}+{\bf g}\\
&\rightarrow& \frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t} &=-\frac{1}{\rho} \nabla P+\nu \nabla^2 {\bf u}+{\bf g}\\
&\rightarrow& \nabla \cdot\left(\frac{\rho}{\Delta t} {\bf u}^{n+1}\right) + \nabla^2 p &= \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}+\rho {\bf g}\right)\\
&\rightarrow& \nabla^2 p &= b, \quad b = \nabla \cdot {{\bf b}^n} = \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}+\rho {\bf g}\right)
\end{align*}
$$

このように${{\bf b}^n}$を定義し，また，ここの$b$を`PoissonRHS`とする．
仮流速は${\bf u}^* = \frac{\Delta t}{\rho}{\bf b}^n$である．

発散の計算方法：

CHECKED: $\nabla\cdot{\bf u}=\sum_{j}\frac{m_j}{\rho_j} \frac{{\bf x}_{ij}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}$

`PoissonRHS`,$b$の計算の前に，$\mu \nabla^2{\bf u}$を予め計算しておく．
今の所，次の順で計算すること．

1. 壁粒子の圧力の計算（流体粒子の現在の圧力$p^n$だけを使って近似）
2. 流体粒子の圧力$p^{n+1}$の計算

ラプラシアンの計算方法：

CHECKED: $`\nabla^2 p^{n+1}=\sum_{j}A_{ij}(p_i^{n+1} - p_j^{n+1}),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$
*/

/*DOC_EXTRACT SPH
### 圧力の安定化

$b =(1-\alpha) \nabla \cdot {{\bf b}^n} + \alpha \frac{\rho - \rho^*}{{\Delta t}^2}$として計算を安定化させる場合がある．

$$
\begin{equation}
\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t,\quad
\frac{D\rho^\ast}{Dt} = - \rho \nabla\cdot{\bf u}^\ast,\quad
\nabla\cdot{\bf u}^\ast = \frac{\Delta t}{\rho} \nabla\cdot{\bf b}^n
\end{equation}
$$

であることから，$(\rho - \rho^*) / {\Delta t^2} = -\nabla\cdot{\bf b}^n$なので，
この安定化は何もしておらず，本来の$b$の計算方法$b =\nabla \cdot {{\bf b}^n} $と同じように見える．

$\rho^*$を計算する際に，$\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t$を使った場合，確かに上のようになるが，
実際に粒子を仮位置に移動させその配置から$\rho^*$を計算した場合は，数値計算上のようにまとめることはできない．

`PoissonRHS`,$b$の計算方法と同じである場合に限る．
もし，計算方法が異なれば，計算方法の違いによって，安定化の効果も変わってくるだろう．

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
         markerX += 1. * A->normal_SPH;
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
      A->div_tmpU = A->PoissonRHS * dt / A->rho;
      A->DrhoDt_SPH = -A->rho * A->div_tmpU;
      A->rho_ = A->rho += A->DrhoDt_SPH * dt;
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

#if defined(Morikawa2019)
      const double alpha = 0.1 * dt;
      // A->PoissonRHS += alpha * (A->rho - A->rho_) / (dt * dt);
      A->PoissonRHS *= 0.5;
#endif

      A->p_SPH_ = (A->PoissonRHS + sum_Aij_Pj) / sum_Aij;
      if (isWall) {
         if (total_weight > 1E-13)
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
/*DOC_EXTRACT SPH

### 圧力勾配$\nabla p^{n+1}$の計算 -> ${D {\bf u}}/{Dt}$の計算

勾配の計算方法：

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