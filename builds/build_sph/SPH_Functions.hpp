#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

networkPoint *getClosestExcludeRigidBodyInlcudeFirstLayer(networkPoint *p, auto &target_nets) {
   double distance = 1E+20;
   networkPoint *P = nullptr;
   for (const auto &obj : target_nets) {
      obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         if (q->isFluid || q->isFirstWallLayer) {
            auto tmp = Distance(p, q);
            if (distance > tmp) {
               distance = tmp;
               P = q;
            }
         }
      });
   }
   return P;
};

networkPoint *getClosestParticle(networkPoint *p, const Network *obj) {
   double distance = 1E+20;
   networkPoint *P = nullptr;
   obj->BucketPoints.apply(p->X, p->radius_SPH * 1.5, [&](const auto &q) {
      auto tmp = Distance(p, q);
      if (distance > tmp) {
         distance = tmp;
         P = q;
      }
   });
   return P;
};

networkPoint *getClosestParticle(networkPoint *p, const std::vector<Network *> objs) {
   double distance = 1E+20;
   networkPoint *P = nullptr;
   for (const auto &obj : objs) {
      auto q = getClosestParticle(p, obj);
      if (q != nullptr) {
         auto tmp = Distance(p, q);
         if (distance > tmp) {
            distance = tmp;
            P = q;
         }
      }
   }
   return P;
};

networkPoint *getClosestFluid(networkPoint *p, auto &target_nets) {
   // std::cout << p << " getClosestFluid" << std::endl;
   double distance = 1E+20;
   networkPoint *P = nullptr;
   for (const auto &obj : target_nets) {
      if (obj->isFluid)
         obj->BucketPoints.apply(p->X, p->radius_SPH * 1.5, [&](const auto &q) {
            if (q->isFluid) {
               auto tmp = Distance(p, q);
               if (distance > tmp) {
                  distance = tmp;
                  P = q;
               }
            }
         });
   }
   return P;
};

/*DOC_EXTRACT SPH

### CFL条件の設定

$`\max({\bf u}) \Delta t \leq c_{v} h \cap \max({\bf a}) \Delta t^2 \leq c_{a} h`$
を満たすように，毎時刻$`\Delta t`$を設定する．

*/

double dt_CFL(const double dt_IN, const auto &net, const auto &RigidBodyObject) {
   double dt = dt_IN;
   const auto C_CFL_velocity = 0.05;  // dt = C_CFL_velocity*h/Max(U)
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

#define Morikawa2019

/* -------------------------------------------------------------------------- */
Tddd aux_position(const networkPoint *p, const double &c) {
   // auto c = p->radius_SPH / p->C_SML;
   return p->X + c * Normalize(p->interpolated_normal_SPH);
   // return p->X - (p->COM_SPH - p->X);
};

Tddd aux_position_next(const networkPoint *p) {
   auto q = p->surfacePoint;
   auto c = q->radius_SPH / q->C_SML;
#if defined(USE_RungeKutta)
   return q->RK_X.getX(q->U_SPH) + c * Normalize(q->interpolated_normal_SPH_next);
#elif defined(USE_LeapFrog)
   return q->LPFG_X.get_x(q->U_SPH) + c * Normalize(q->interpolated_normal_SPH_next);
#endif
};

// \label{SPH:rho_next}
double rho_next_(auto p) {
   if (p->isAuxiliary)
      p = p->surfacePoint;
   if (p->getNetwork()->isRigidBody)
      return _WATER_DENSITY_;
   else {
#if defined(USE_RungeKutta)
      return p->RK_rho.getX(p->DrhoDt_SPH);
#elif defined(USE_LeapFrog)
      return p->LPFG_rho.get_x(p->DrhoDt_SPH);
#endif
   }
};

// \label{SPH:volume_next}
double V_next_(const auto &p) {
   return p->mass / rho_next(p);
};

// \label{SPH:position_next}
std::array<double, 3> X_next_(const auto &p) {
   if (p->isAuxiliary)
      return p->X + p->surfacePoint->LPFG_X.get_x(p->surfacePoint->U_SPH) - p->surfacePoint->X;
   else if (p->getNetwork()->isRigidBody)
      return p->X;
   else
#if defined(USE_RungeKutta)
      return p->RK_X.getX(p->U_SPH);
#elif defined(USE_LeapFrog)
      return p->LPFG_X.get_x(p->U_SPH);
         // return p->X + p->U_SPH * p->LPFG_X.get_dt();
#endif
};

/* -------------------------------------------------------------------------- */

// \label{SPH:rho_next}
double rho_next(auto p) {
   // return rho_next_(p);
   // return rho_next_(p);
   return _WATER_DENSITY_;
};

// \label{SPH:volume_next}
double V_next(const auto &p) {
   return p->mass / rho_next(p);
};

// \label{SPH:position_next}
std::array<double, 3> X_next(const auto &p) {
   return p->X;
   // return X_next_(p);
};

/* -------------------------------------------------------------------------- */

#include "SPH_setWall_Freesurface.hpp"

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT SPH

## $`\nabla^2 {\bf u}_i`$の計算

CHECKED: \ref{SPH:lapU}{ラプラシアンの計算方法}: $`\nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

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
      A->div_U = 0.;
      A->lap_U.fill(0.);
      A->b_vector.fill(0.);
      //
      A->grad_coeff.clear();
      A->grad_coeff_next.clear();
      //$ ------------------------------------------ */
      /*DOC_EXTRACT SPH

      ### 高速化のための工夫

      何度か行う勾配の計算は，変数は違えど，変数の係数は同じである．
      ここで，その係数を`std::unordered_map`で保存しておくことにする．
      `A->grad_coeff`と`A->grad_coeff_next`に保存する．

      NOTE: `A->grad_coeff`と`A->grad_coeff_next`は，自身もキーとして含む．使う時に注意する．

      */

      A->density_based_on_positions = 0;

      auto add_lap_U = [&](const auto &B) {
         const auto Uij = A->U_SPH - B->U_SPH;
         A->density_based_on_positions += A->rho * A->volume * w_Bspline(Norm(A->X - B->X), A->radius_SPH);
         A->div_U += B->volume * Dot(B->U_SPH - A->U_SPH, grad_w_Bspline(A->X, B->X, A->radius_SPH));
         A->lap_U += 2 * B->mass / A->rho * Uij * Dot_grad_w_Bspline_Dot(A->X, B->X, A->radius_SPH);  //\label{SPH:lapU}

         // just counting
         if (Between(Distance(A, B), {1E-12, A->radius_SPH})) {
            A->checked_points_in_radius_SPH++;
            if (B->getNetwork()->isFluid || B->isFluid)
               A->checked_points_in_radius_of_fluid_SPH++;

            // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP1}
            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->mass / A->rho * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP2}
            // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP3}
         }
         A->checked_points_SPH++;
      };
      // sum 計算
      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               // 全ての壁粒子の流速はゼロなのだから，isCapturedされていないものを含めても問題ない．
               // ただこの後の圧力の計算においては，isCapturedされていないものは含めない．圧力方程式をうまく立てれないから．
               add_lap_U(B);
#if defined(USE_SHARED_AUX)
               if (B->isSurface && B == A)
                  for (const auto &AUX : B->auxiliaryPoints)
                     add_lap_U(AUX);
#endif
            }
         });
#if defined(USE_SIMPLE_SINGLE_AUX)
      if (A->isSurface)
         for (const auto &AUX : A->auxiliaryPoints)
            PoissonEquation(AUX);
#endif
      //$ ------------------------------------------ */
      //\label{SPH:lapU_for_wall}
      if (A->getNetwork()->isRigidBody) {
         A->DUDt_SPH_.fill(0.);
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH.fill(0.);
         A->tmp_U_SPH.fill(0.);
         A->tmp_X = A->X;
         A->DrhoDt_SPH = 0;
         A->b_vector.fill(0.);  // + _GRAVITY3_;
      } else {
         A->DUDt_SPH_ = A->DUDt_SPH;
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH = nu * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt
         A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
         A->tmp_X = A->X + A->tmp_U_SPH * dt;
         A->DrhoDt_SPH = -A->rho * A->div_U;
         A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
      }
      //$ ------------------------------------------ */
      // \label{SPH:Poisson_b_vector}

      // auto add_b_vector = [&](const auto &B) {
      //    auto w = B->volume * w_Bspline(Norm(A->X - B->X), A->radius_SPH);
      //    A->b_vector += w * (B->U_SPH / dt + B->mu_SPH / B->rho * B->lap_U);  // + (A->rho * _GRAVITY3_);
      // };
      // // sum 計算
      // for (const auto &net : target_nets)
      //    net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
      //       if (B->isCaptured) {
      //          add_b_vector(B);
      //       }
      //    });

      if (A->vec_time_SPH.size() > 10) {

#if defined(USE_RungeKutta)
         double current_time = A->RK_X.get_t();
         double next_time = current_time + A->RK_X.get_dt();
#elif defined(USE_LeapFrog)
         double current_time = A->LPFG_X.get_t();
         double next_time = current_time + dt;
#endif
         std::vector<double> times = {next_time, current_time};
         std::array<double, 3> U1, U2, U3, U4;
         U1 = A->U_SPH;
         if (*(A->vec_time_SPH.rbegin()) == current_time) {
            times.push_back(*(A->vec_time_SPH.rbegin() + 1));
            U2 = *(A->vec_U_SPH.rbegin() + 1);
            times.push_back(*(A->vec_time_SPH.rbegin() + 2));
            U3 = *(A->vec_U_SPH.rbegin() + 2);
            // times.push_back(*(A->vec_time_SPH.rbegin() + 3));
            // U4 = *(A->vec_U_SPH.rbegin() + 3);
         } else {
            times.push_back(*(A->vec_time_SPH.rbegin() + 0));
            U2 = *(A->vec_U_SPH.rbegin() + 0);
            times.push_back(*(A->vec_time_SPH.rbegin() + 1));
            U3 = *(A->vec_U_SPH.rbegin() + 1);
            // times.push_back(*(A->vec_time_SPH.rbegin() + 2));
            // U4 = *(A->vec_U_SPH.rbegin() + 2);
         }
         InterpolationLagrange<double> lag(times);
         auto D = lag.DN(current_time);
         A->b_vector = -(D[1] * U1 + D[2] * U2 + D[3] * U3 /* + D[4] * U4*/) + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
      }
   }
};

#include "SPH_FindPressure.hpp"

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */

/*DOC_EXTRACT SPH

## 圧力勾配$`\nabla p^{n+1}`$の計算

CHECKED: \ref{SPH:gradP1}{勾配の計算方法}: $`\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}`$

CHECKED: \ref{SPH:gradP2}{勾配の計算方法}: $`\nabla p_i = \rho_i \sum_{j} m_j \left(p_j - p_i\right) \nabla W_{ij}`$

CHECKED: \ref{SPH:gradP3}{勾配の計算方法}: $`\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}`$

*/

void gradP(const std::unordered_set<networkPoint *> &points, const std::unordered_set<Network *> &target_nets) {

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->gradP_SPH.fill(0.);
      auto add_gradP_SPH = [&](const auto &B, const double &coef = 1.) {
         A->gradP_SPH += coef * A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP1}0.2647

         // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->mass / A->rho * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP2}

         //\label{SPH:gradP2}は，Aij = 1*..を使うと，0.132
         //\label{SPH:gradP2}は，Aij = 3*..を使うと，0.221
         //\label{SPH:gradP2}は，Aij = 2*..を使うと，0.21198

         // A->gradP_SPH += (B->p_SPH - A->p_SPH) * V_next(B) * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP2}
         // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP3}0.34
         //\label{SPH:gradP3}は，Aij = 2*..を使うと，0.34
         //\label{SPH:gradP3}は，Aij = 2*..を使うと，Lagrangeを使うと，0.46
         //\label{SPH:gradP3}は，Aij = 3*..を使うと，0.49
         //\label{SPH:gradP3}は，Aij = 3*..を使うと，Lagrangeを使うと，0.4は超えたが，綺麗な結果ではなく，つぶれた．
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，_nextを使うと，0.24
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，X_next以外の_nextを使うと，0.277>
         //\label{SPH:gradP3}は，Aij = 2*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.34
         //\label{SPH:gradP3}は，Aij = 2.5*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.45
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.49
         //\label{SPH:gradP3}は，Aij = 4*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.457
      };

      for (const auto &net : target_nets) {
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               add_gradP_SPH(B);
#if defined(USE_SHARED_AUX)
               if (B->isSurface && B == A)
                  for (const auto &AUX : B->auxiliaryPoints)
                     add_gradP_SPH(AUX);
#endif
            }
         });
      }
#if defined(USE_SIMPLE_SINGLE_AUX)
      if (A->isSurface)
         for (const auto &AUX : A->auxiliaryPoints)
            add_gradP_SPH(AUX);
#endif

      /*DOC_EXTRACT SPH

      $`\dfrac{D{\bf u}^n}{Dt} = - \frac{1}{\rho} \nabla p^{n+1} + \nu \nabla^2 {\bf u}^n + {\bf g}`$
      が計算できた．

      */

      A->DUDt_SPH -= A->gradP_SPH / A->rho;

      if (!isFinite(A->DUDt_SPH))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
   }
}

//@ -------------------------------------------------------- */
//@                        粒子の時間発展                      */
//@ -------------------------------------------------------- */
#define REFLECTION
void updateParticles(const auto &points,
                     const std::unordered_set<Network *> &target_nets,
                     const auto &RigidBodyObject,
                     const double &particle_spacing,
                     const double dt) {
   DebugPrint("粒子の時間発展", Green);
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
   {
      // テスト
      auto U = p->U_SPH;
      auto X_last = p->X;
#if defined(USE_RungeKutta)
      auto getX_next = [&](const auto &p) { return p->RK_X.getX(p->U_SPH); };
      if (p->RK_X.steps == 1) {
         p->RK_U.push(p->DUDt_SPH);  // 速度
         p->U_SPH = p->RK_U.getX();
         p->RK_X.push(p->U_SPH);  // 位置
         p->setXSingle(p->tmp_X = p->RK_X.getX());
      } else {
         p->RK_X.push(p->U_SPH);  // 位置
         p->setXSingle(p->tmp_X = p->RK_X.getX());
         p->RK_U.push(p->DUDt_SPH);  // 速度
         p->U_SPH = p->RK_U.getX();
      }
         // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
#elif defined(USE_LeapFrog)
      auto getX_next = [&](const auto &p) { return X_next(p); };
      p->DUDt_modify_SPH.fill(0.);
      p->LPFG_X.push(p->DUDt_SPH);  // 速度
      p->U_SPH = p->LPFG_X.get_v();
      p->setXSingle(p->tmp_X = p->LPFG_X.get_x());
         // auto getX = [&](const auto &p) { return p->X; };
#endif

#if defined(REFLECTION)
      int count = 0;
      //\label{SPH:reflection}
      const double reflection_factor = .8;
      const double asobi = 0.05;

      auto closest = [&]() {
         double distance = 1E+20;
         networkPoint *P = nullptr;
         for (const auto &[obj, _] : RigidBodyObject) {
            obj->BucketPoints.apply(X_next(p), p->radius_SPH, [&](const auto &q) {
               auto tmp = Distance(X_next(p), q);
               if (distance > tmp) {
                  distance = tmp;
                  P = q;
               }
            });
         }
         return P;
      };

      bool isReflected = true;
      while (isReflected && count++ < 10) {
         isReflected = false;
         for (const auto &[obj, _] : RigidBodyObject)
            obj->BucketPoints.apply(getX_next(p), p->radius_SPH, [&](const auto &Pwall) {
               auto dist = Distance(Pwall->X, getX_next(p));
               if (dist < (1 - asobi) * particle_spacing) {
                  auto alpha = (particle_spacing - dist) / particle_spacing;
                  auto n = Pwall->interpolated_normal_SPH;
                  if (Dot(p->U_SPH, n) < 0) {
                     auto tmp = -(1. + reflection_factor) * Projection(p->U_SPH, n) / dt;
                     p->DUDt_modify_SPH += tmp;
                     p->DUDt_SPH += tmp;
   #if defined(USE_RungeKutta)
                     p->RK_U.repush(p->DUDt_SPH);  // 速度
                     if (p->RK_X.steps == 1) {
                        p->U_SPH = p->RK_U.getX();
                        p->RK_X.repush(p->U_SPH);  // 位置
                        p->setXSingle(p->tmp_X = p->RK_X.getX());
                     } else {
                        p->RK_U.repush(p->DUDt_SPH);  // 速度
                        p->U_SPH = p->RK_U.getX();
                     }
                     isReflected = true;
   #elif defined(USE_LeapFrog)
                     p->LPFG_X.repush(p->DUDt_SPH);  // 速度
                     p->U_SPH = p->LPFG_X.get_v();
                     p->setXSingle(p->tmp_X = p->LPFG_X.get_x());
                     isReflected = true;
   #endif
                  }
               }
            });
      };

#endif
   }

   // \label{SPH:update_density}
   for (const auto &A : points) {
#if defined(USE_RungeKutta)
      A->DrhoDt_SPH = -A->rho * A->div_U;
      A->RK_rho.push(A->DrhoDt_SPH);  // 密度
                                      // A->setDensity(A->RK_rho.get_x());
      A->setDensity(_WATER_DENSITY_);
#elif defined(USE_LeapFrog)
      A->DrhoDt_SPH = -A->rho * A->div_U;
      A->LPFG_rho.push(A->DrhoDt_SPH);
      // A->setDensity(A->LPFG_rho.get_x());
      A->setDensity(_WATER_DENSITY_);
#endif
   }
}

/*DOC_EXTRACT SPH

## 注意点

WARNING: 計算がうまく行く設定を知るために，次の箇所をチェックする．

**壁粒子**

- \ref{SPH:lapU_for_wall}{壁粒子のラプラシアンの計算方法}
- \ref{SPH:setPoissonEquation}{圧力の計算方法}
   - \ref{SPH:whereToMakeTheEquation}{どの位置において方程式を立てるか}
- \ref{SPH:capture_condition_1st}{流体として扱う壁粒子を設定するかどうか}/\ref{SPH:capture_condition_2nd}{視野角に流体粒子が含まない壁粒子は除外する}
- \ref{SPH:map_fluid_pressure_to_wall}{壁粒子の圧力をどのように壁面にマッピングするか}
- \ref{SPH:interpolated_normal_SPH}{壁粒子の法線方向ベクトルの計算方法}
- \ref{SPH:reflection}{反射の計算方法}

**水面粒子**

- \ref{SPH:water_surface_pressure}{水面粒子の圧力をゼロにするかどうか}
- \ref{SPH:auxiliaryPoints}{補助粒子の設定はどうなっているか}

**その他**

- \ref{SPH:update_density}{密度を更新するかどうか}
- \ref{SPH:pressure_stabilization}{圧力の安定化をするかどうか}
- \ref{SPH:RK_order}{ルンゲクッタの段数}


壁のwall_as_fluidは繰り返しで計算するのはどうか？

*/

#endif