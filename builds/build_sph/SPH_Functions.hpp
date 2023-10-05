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

networkPoint *getClosestFluid(networkPoint *p, const auto &target_nets) {
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

networkPoint *getClosest(const Tddd X, const double range, const auto &target_nets, const auto &condition) {
   double distance = 1E+20;
   networkPoint *P = nullptr;
   for (const auto &obj : target_nets) {
      if (obj->isFluid)
         obj->BucketPoints.apply(X, range, [&](const auto &q) {
            if (condition(q)) {
               auto tmp = Distance(X, q);
               if (distance > tmp) {
                  distance = tmp;
                  P = q;
               }
            }
         });
   }
   return P;
};

/*DOC_EXTRACT 0_1_0_SPH

### CFL条件の設定

$`\max({\bf u}) \Delta t \leq c_{v} h \cap \max({\bf a}) \Delta t^2 \leq c_{a} h`$
を満たすように，毎時刻$`\Delta t`$を設定する．

*/

double dt_CFL(const double dt_IN, const auto &net, const auto &RigidBodyObject) {
   double dt = dt_IN;
   const auto C_CFL_velocity = 0.05;  // dt = C_CFL_velocity*h/Max(U)
   const auto C_CFL_accel = 0.05;     // dt = C_CFL_accel*sqrt(h/Max(A))
   for (const auto &p : net->getPoints()) {
      // 速度に関するCFL条件
      auto dt_C_CFL = [&](const auto &q) {
         if (p != q) {
            auto pq = Normalize(p->X - q->X);
            auto distance = Distance(p, q);
            /* ------------------------------------------------ */
            // 相対速度
            // double max_dt_vel = C_CFL_velocity * distance / std::abs(Dot(p->U_SPH - q->U_SPH, pq));
            double max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH - q->U_SPH);
            if (dt > max_dt_vel && isFinite(max_dt_vel))
               dt = max_dt_vel;
            // 絶対速度
            // max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH);
            if (dt > max_dt_vel && isFinite(max_dt_vel))
               dt = max_dt_vel;
            /* ------------------------------------------------ */
            // 相対速度
            // double max_dt_acc = C_CFL_accel * std::sqrt(distance / std::abs(Dot(p->DUDt_SPH - q->DUDt_SPH, pq)));
            double max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH - q->DUDt_SPH));
            if (dt > max_dt_acc && isFinite(max_dt_acc))
               dt = max_dt_acc;
            // 絶対速度
            // max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH));
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
   return p->X + c * Normalize(p->interp_normal);
   // return p->X - (p->COM_SPH - p->X);
};

Tddd aux_position_next(const networkPoint *p) {
   auto q = p->surfacePoint;
   auto c = q->radius_SPH / q->C_SML;
#if defined(USE_RungeKutta)
   return q->RK_X.getX(q->U_SPH) + c * Normalize(q->interp_normal_next);
#elif defined(USE_LeapFrog)
   return q->LPFG_X.get_x(q->U_SPH) + c * Normalize(q->interp_normal_next);
#endif
};

/* -------------------------------------------------------------------------- */
// \label{SPH:rho_next}
double rho_next(auto p) {
   // return _WATER_DENSITY_;
   /* -------------------------------------------------------------------------- */
   //    if (p->isAuxiliary)
   //       return rho_next(p->surfacePoint);
   //    else if (p->getNetwork()->isRigidBody)
   //       return _WATER_DENSITY_;
   //    else {
   // #if defined(USE_RungeKutta)
   return p->RK_rho.getX(p->DrhoDt_SPH);
   //@ これを使った方が安定するようだ
   // #elif defined(USE_LeapFrog)
   //       return p->rho + p->DrhoDt_SPH + p->RK_rho.get_dt();
   // #endif
   //    }
};

// \label{SPH:volume_next}
double V_next(const auto &p) { return p->mass / rho_next(p); };

// \label{SPH:position_next}
std::array<double, 3> X_next(const auto &p) {
   if (p->getNetwork()->isRigidBody)
      return p->X;
   else {
#if defined(USE_RungeKutta)
      // return p->RK_X.getX(p->RK_U.getX(p->DUDt_SPH));
      //
      // Tddd U;
      // if (p->interp_U_Bspline != nullptr) {
      //    U = (*p->interp_U_Bspline)(p->RK_X.getNextTime());
      // } else
      // if (p->interp_U_lag != nullptr) {
      //    U = (*p->interp_U_lag)(p->RK_X.getNextTime());
      // } else
      // U = p->U_SPH;

      // if (p->interp_U_lag != nullptr) {
      //    U = (*p->interp_U_lag)(p->RK_X.getNextTime());
      // } else
      //    U = p->U_SPH;
      return p->RK_X.getX(p->U_SPH);
         // return p->RK_X.getX(p->U_SPH);
#elif defined(USE_LeapFrog)
      return p->X + p->U_SPH * p->LPFG_X.get_dt() / 2.;
#endif
   }
};

/* -------------------------------------------------------------------------- */
#include "SPH0_setWall_Freesurface.hpp"
//
#include "SPH1_lap_div_U.hpp"
//
#include "SPH2_FindPressure.hpp"
//
#include "SPH3_grad_p.hpp"
//

void set_interp_U(networkPoint *A) {
   int N = 6;
   if (A->vec_U_SPH.size() > N + 4 && A->vec_time_SPH.size() > N + 4) {
      // std::cout << "set_interp_U" << std::endl;
      // std::cout << "A->vec_U_SPH = " << A->vec_U_SPH << std::endl;
      // std::cout << "A->vec_time_SPH = " << A->vec_time_SPH << std::endl;
#if defined(USE_RungeKutta)
      double current_time = A->RK_X.get_t();
#elif defined(USE_LeapFrog)
      double current_time = A->LPFG_X.get_t();
#endif
      std::vector<double> times(N);
      std::vector<std::array<double, 3>> U(N);
      for (int i = 0; i < N; i++) {
         times[i] = *(A->vec_time_SPH.rbegin() + i);
         U[i] = *(A->vec_U_SPH.rbegin() + i);
      }

      if (A->interp_U_lag != nullptr)
         delete A->interp_U_lag;
      A->interp_U_lag = new InterpolationLagrange<std::array<double, 3>>(times, U);
      // std::cout << "A->interp_U_lag has been set" << std::endl;

      if (A->interp_U_Bspline != nullptr)
         delete A->interp_U_Bspline;
      A->interp_U_Bspline = new InterpolationBspline(3, times, U);
      // std::cout << "A->interp_U_Bspline has been set" << std::endl;
   }
};

//@ -------------------------------------------------------- */
//@                        粒子の時間発展                      */
//@ -------------------------------------------------------- */
void set_nearest_wall_p_next(networkPoint *p, auto &RigidBodyObject) {
   double nearest_dist = 1E+20, d;
   networkPoint *P = nullptr;
   std::array<double, 3> p_to_q;
   for (const auto &[obj, _] : RigidBodyObject) {
      obj->BucketPoints.apply(X_next(p), p->radius_SPH, [&](const auto &q) {
         p_to_q = q->X - X_next(p);
         d = Norm(p_to_q);
         if (nearest_dist > d) {
            nearest_dist = d;
            P = q;
         }
      });
   }
   p->nearest_wall_p_next = P;
};

void set_nearest_wall_p_next(auto &points, auto &RigidBodyObject) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
      set_nearest_wall_p_next(p, RigidBodyObject);
};

void set_nearest_wall_p(networkPoint *p, auto &RigidBodyObject) {
   double nearest_dist = 1E+20, d;
   networkPoint *P = nullptr;
   std::array<double, 3> p_to_q;
   for (const auto &[obj, _] : RigidBodyObject) {
      obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         p_to_q = q->X - p->X;
         d = Norm(p_to_q);
         if (nearest_dist > d) {
            nearest_dist = d;
            P = q;
         }
      });
   }
   p->nearest_wall_p = P;
};

void set_nearest_wall_p(auto &points, auto &RigidBodyObject) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
      set_nearest_wall_p(p, RigidBodyObject);
};

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
      auto U = p->U_SPH;
      auto X_last = p->X;
#if defined(USE_RungeKutta)
      auto X_next = [&](const auto &p) {
         // if (p->RK_X.steps == 1)
         //    return p->RK_X.getX();
         // else
         return p->RK_X.getX(p->U_SPH);
      };
      // if (p->RK_X.steps == 1) {
      //    p->RK_U.push(p->DUDt_SPH);  // 速度
      //    p->U_SPH = p->RK_U.getX();
      //    p->RK_X.push(p->U_SPH);  // 位置
      //    p->setXSingle(p->tmp_X = p->RK_X.getX());
      // } else
      {
         p->RK_X.push(p->U_SPH);  // 位置
         p->setXSingle(p->tmp_X = p->RK_X.getX());
         p->RK_U.push(p->DUDt_SPH);  // 速度
         p->U_SPH = p->RK_U.getX();
      }
         // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
#elif defined(USE_LeapFrog)
      auto X_next = [&](const auto &p) { return p->X; };
      p->DUDt_modify_SPH.fill(0.);
      p->LPFG_X.push(p->DUDt_SPH);  // 速度
      p->U_SPH = p->LPFG_X.get_v();
      p->setXSingle(p->tmp_X = p->LPFG_X.get_x());
         // auto getX = [&](const auto &p) { return p->X; };
#endif

#if defined(REFLECTION)

      auto closest = [&]() {
         double distance = 1E+20;
         networkPoint *P = nullptr;
         for (const auto &[obj, _] : RigidBodyObject) {
            obj->BucketPoints.apply(X_next(p), p->radius_SPH, [&](const auto &q) {
               auto p_to_q = q->X - X_next(p);
               auto tmp = Norm(p_to_q);
               if (distance > tmp) {
                  distance = tmp;
                  P = q;
               }
            });
         }
         return P;
      };

      int count = 0;
      //\label{SPH:reflection}
      const double reflection_factor = .1;
      // const double reflection_factor = 1.;
      const double c = 0.;

      bool isReflected = true;
      p->DUDt_modify_SPH.fill(0.);
      while (isReflected && count++ < 2) {
         isReflected = false;
         auto closest_p = closest();
         if (closest_p != nullptr) {
            // for (const auto &[obj, _] : RigidBodyObject)
            //    obj->BucketPoints.apply(X_next(p), p->radius_SPH, [&](const auto &Pwall) {
            // auto dist = Distance(X_next(closest_p), X_next(p));
            auto d_ps = particle_spacing;
            auto d0 = (1 - c) * particle_spacing;
            auto n = Normalize(closest_p->interp_normal_original);
            auto v_f2w = X_next(closest_p) - X_next(p);
            auto n_d_f2w = Norm(Projection(v_f2w, n));
            if (Norm(v_f2w) < d0) {
               auto ratio = (d0 - n_d_f2w) / d0;
               if (Dot(p->U_SPH, n) < 0) {
                  auto tmp = Norm(d_ps - n_d_f2w) * n / dt;
                  p->DUDt_modify_SPH += tmp;
                  p->DUDt_SPH += tmp;
                     // p->DUDt_SPH -= Projection(tmp, p->DUDt_SPH);
   #if defined(USE_RungeKutta)
                  p->RK_U.repush(p->DUDt_SPH);  // 速度
                  // if (p->RK_X.steps == 1) {
                  //    p->U_SPH = p->RK_U.getX();
                  //    p->RK_X.repush(p->U_SPH);  // 位置
                  //    p->setXSingle(p->tmp_X = p->RK_X.getX());
                  // } else
                  {
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
         }
         // });
      };

#endif
   }

   // \label{SPH:update_density}
   for (const auto &A : points) {
#if defined(USE_RungeKutta)
      A->DrhoDt_SPH = -A->rho * A->div_U;
      A->RK_rho.push(A->DrhoDt_SPH);  // 密度

      // if (A->RK_rho.finished)
      //    A->setDensity(_WATER_DENSITY_);
      // else
      // A->setDensity(A->RK_rho.get_x());
      // A->setDensity((A->RK_rho.get_x() + _WATER_DENSITY_) / 2.);
      // A->setDensity(A->RK_rho.get_x());
      A->setDensity(_WATER_DENSITY_);
#elif defined(USE_LeapFrog)
      A->DrhoDt_SPH = -A->rho * A->div_U;
      A->LPFG_rho.push(A->DrhoDt_SPH);
      // リープフロッグの密度更新はこれでほんとにいいのか？
      // divの補正＆ぁプラしなんの補正
      //     A->setDensity(_WATER_DENSITY_);
      // if (A->LPFG_rho.finished)
      A->setDensity(_WATER_DENSITY_);
         // else
         // A->setDensity(A->LPFG_rho.get_x());
#endif

      set_interp_U(A);
   }

   set_nearest_wall_p_next(points, RigidBodyObject);
   set_nearest_wall_p(points, RigidBodyObject);
}

/*DOC_EXTRACT 0_3_0_SPH

## 注意点

WARNING: 計算がうまく行く設定を知るために，次の箇所をチェックする．

**NEW**

- \ref{SPH:wall_particle_velocity}{壁粒子の速度の決定方法}
- \ref{SPH:how_to_use_b_vector_in_Poisson0}{Poissonにおいてどのようにbベクトルを使うか}
- \ref{SPH:how_to_use_b_vector_in_Poisson1}{Poissonにおいてどのようにbベクトルを使うか}
- どのように\ref{SPH:how_to_set_wall_b_vector}{壁粒子のb}/\ref{SPH:how_to_set_fluid_b_vector}{流体粒子のb}を作るか

**壁粒子**

- \ref{SPH:lapU_for_wall}{壁粒子のラプラシアンの計算方法}
- \ref{SPH:setPoissonEquation}{圧力の計算方法}
   - \ref{SPH:whereToMakeTheEquation}{どの位置において方程式を立てるか}
- \ref{SPH:capture_condition_1st}{流体として扱う壁粒子を設定するかどうか}/\ref{SPH:capture_condition_2nd}{視野角に流体粒子が含まない壁粒子は除外する}
- \ref{SPH:map_fluid_pressure_to_wall}{壁粒子の圧力をどのように壁面にマッピングするか}
- \ref{SPH:interp_normal}{壁粒子の法線方向ベクトルの計算方法}
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