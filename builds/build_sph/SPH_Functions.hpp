#ifndef SPH_Functions_H
#define SPH_Functions_H
// old
#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/* -------------------------------------------------------------------------- */

bool canInteract(const networkPoint *A, const networkPoint *B) {
   if (!B->isCaptured)
      return false;
   /* ----------------------------------- */
   // if (A->isAuxiliary && B->isAuxiliary)
   //    return false;
   // else if (A->isAuxiliary && !B->isAuxiliary)
   //    return true;
   // else if (!A->isAuxiliary && B->isAuxiliary) {
   //    if (B->surfacePoint == A)
   //       return true;
   //    else
   //       return false;
   // } else if (!A->isAuxiliary && !B->isAuxiliary)
   //    return true;

   // return false;
   /* ----------------------------------- */
   // if (A->isAuxiliary && B->isAuxiliary)
   //    return false;
   // else if (A->isAuxiliary || B->isAuxiliary) {
   //    if (A->isAuxiliary && A->surfacePoint == B)
   //       return true;
   //    if (B->isAuxiliary && B->surfacePoint == A)
   //       return true;
   // }
   /* ----------------------------------- */
   // const double c = .3;
   // if (B->isAuxiliary) {
   //    if (Distance(A, B) < c * A->SML())
   //       return true;
   //    else
   //       return false;
   // }
   // if (B->auxPoint != nullptr) {
   //    if (B->auxPoint == A)
   //       return false;
   //    if (Distance(A, B) > c * A->SML())
   //       return true;
   //    else
   //       return false;
   // }

   return true;
};

/* -------------------------------------------------------------------------- */
void deleteAuxiliaryPoints(const auto net) {
   std::cout << "delete auxiliary points" << std::endl;
   auto points = net->getPoints();
   for (auto p : points) {
      p->auxPoint = nullptr;
      if (p->isAuxiliary) {
         p->surfacePoint->auxPoint = nullptr;
         delete p;
      }
   }

   for (auto p : net->getPoints())
      p->isAuxiliary = false;

   net->remakeBucketPoints();
   std::cout << "delete auxiliary points done" << std::endl;
}

Tddd aux_position(const networkPoint *p) {
   auto sp = p->surfacePoint;
#ifdef SET_AUX_AT_PARTICLE_SPACING
   return sp->X + sp->particle_spacing * Normalize(Dot(sp->inv_grad_corr_M, sp->interp_normal_original));
#elif defined(SET_AUX_AT_MASS_CENTER)
   return sp->X - sp->vec2COM;
#else
   return sp->X + sp->particle_spacing * Normalize(Dot(sp->inv_grad_corr_M, sp->interp_normal_original));
#endif
};

void setAuxiliaryPoints(const auto net, const std::function<Tddd(networkPoint *)> &position = aux_position) {
   Print("setAuxiliaryPoints");
   auto points = net->getPoints();
   for (const auto &p : points) {
      if (p->hasAuxiliary()) {

         auto q = new networkPoint(net, p->X);
         q->surfacePoint = p;
         p->auxPoint = q;
         q->setX(position(q));
         //
         q->grad_corr_M = p->grad_corr_M;
         q->grad_corr_M_next = p->grad_corr_M_next;
         q->inv_grad_corr_M = p->inv_grad_corr_M_next;
         q->inv_grad_corr_M_next = p->inv_grad_corr_M_next;
         //
         // q->grad_corr_M_rigid = p->grad_corr_M_next_rigid;
         // q->grad_corr_M_next_rigid = p->grad_corr_M_next_rigid;
         // q->inv_grad_corr_M_rigid = p->inv_grad_corr_M_next_rigid;
         // q->inv_grad_corr_M_next_rigid = p->inv_grad_corr_M_next_rigid;
         q->laplacian_corr_M = p->laplacian_corr_M;
         q->laplacian_corr_M_next = p->laplacian_corr_M_next;
         //
         q->isSurface = p->isSurface;
         q->isSurface_next = p->isSurface_next;
         q->isNeumannSurface = p->isNeumannSurface;
         q->isAuxiliary = true;
         //
         q->b_vector = p->b_vector;
         q->U_SPH = p->U_SPH;  // - 2 * Chop(p->U_SPH, p->interp_normal);
         //
         q->intp_density = p->intp_density;
         q->intp_density_next = p->intp_density_next;
         //
         // q->U_SPH.fill(0.);
         q->v_to_surface_SPH = p->v_to_surface_SPH;
         q->interp_normal = p->interp_normal;
         q->interp_normal_next = p->interp_normal_next;
         q->intp_normal_Eigen = p->intp_normal_Eigen;
         q->interp_normal_original = p->interp_normal_original;
         q->interp_normal_original_next = p->interp_normal_original_next;
         q->intp_density = p->intp_density;
         //
         q->div_U = p->div_U;
         q->DUDt_SPH = p->DUDt_SPH;
         q->lap_U = p->lap_U;
         q->p_SPH = p->p_SPH;
         q->rho = p->rho;
         q->setDensityVolume(_WATER_DENSITY_, p->volume);
         q->particle_spacing = p->particle_spacing;
         q->C_SML_next = p->C_SML_next;
         q->C_SML = p->C_SML;
         q->isFluid = p->isFluid;
         q->isFirstWallLayer = false;
         q->isCaptured = true;
         //
         auto dt = p->RK_X.getdt();
         // q->RK_U.initialize(dt, simulation_time, q->U_SPH, 1);
         // q->RK_X.initialize(dt, simulation_time, q->X, 1);
         // q->RK_P.initialize(dt, simulation_time, q->p_SPH, 1);
         // q->RK_rho.initialize(dt, simulation_time, q->rho, 1);
         q->RK_U = p->RK_U;
         q->RK_U.Xinit = p->RK_U.Xinit;
         q->RK_U.t_init = p->RK_U.t_init;
         q->RK_U.dt_fixed = p->RK_U.dt_fixed;
         q->RK_U.dt = p->RK_U.dt;
         q->RK_U.steps = p->RK_U.steps;
         q->RK_U.current_step = p->RK_U.current_step;
         q->RK_U._dX = p->RK_U._dX;
         q->RK_U.dX = p->RK_U.dX;
         //
         q->RK_X = p->RK_X;
         // q->RK_X.Xinit = p->RK_X.Xinit;
         q->RK_X.Xinit = q->X;
         q->RK_X.t_init = p->RK_X.t_init;
         q->RK_X.dt_fixed = p->RK_X.dt_fixed;
         q->RK_X.dt = p->RK_X.dt;
         q->RK_X.steps = p->RK_X.steps;
         q->RK_X.current_step = p->RK_X.current_step;
         q->RK_X._dX = p->RK_X._dX;
         q->RK_X.dX = p->RK_X.dX;
         //
         q->RK_rho = p->RK_rho;
         q->RK_rho.Xinit = p->RK_rho.Xinit;
         q->RK_rho.t_init = p->RK_rho.t_init;
         q->RK_rho.dt_fixed = p->RK_rho.dt_fixed;
         q->RK_rho.dt = p->RK_rho.dt;
         q->RK_rho.steps = p->RK_rho.steps;
         q->RK_rho.current_step = p->RK_rho.current_step;
         q->RK_rho._dX = p->RK_rho._dX;
         q->RK_rho.dX = p->RK_rho.dX;
      }
   }
   Print("setAuxiliaryPoints done");
   Print("remakeBucketPoints");
   net->remakeBucketPoints();
   Print("remakeBucketPoints done");
}

/* -------------------------------------------------------------------------- */

void setSML(const auto &target_nets) {
   DebugPrint("setSML", Yellow);
   /* -------------------------------- C_SMLの調整 -------------------------------- */
   const double C_SML_max = 3.;
   const double C_SML_min = 2.;
   const double C_SML_min_rigid = 2.;
   for (const auto &NET : target_nets)
      if (NET->isFluid) {
         {
#pragma omp parallel
            for (const auto &p : NET->getPoints())
#pragma omp single nowait
            {
               p->isNearSurface = false;
               if (p->isCaptured) {
                  double closest_d = 1E+20, closest_d_next = 1E+20, d;
                  networkPoint *closest_q = nullptr, *closest_q_next = nullptr;
                  for (const auto &net : target_nets) {
                     net->BucketPoints.apply(p->X, 3.5 * p->particle_spacing, [&](const auto &q) {
                        if ((q->isSurface || q->isNeumannSurface) && (d = Norm(q->X - p->X)) < closest_d) {
                           closest_d = d;
                           closest_q = q;
                           p->C_SML = closest_d / p->particle_spacing + 0.5;
                        }
                        if ((q->isSurface || q->isNeumannSurface) && (d = Norm(X_next(q) - X_next(p))) < closest_d_next) {
                           closest_d_next = d;
                           closest_q_next = q;
                           p->C_SML_next = closest_d_next / p->particle_spacing + 0.5;
                        }
                        if (q->isSurface && Distance(p, q) < 2. * p->particle_spacing)
                           p->isNearSurface = true;
                     });
                  }
                  if (closest_q != nullptr) {
                     if (closest_q->isSurface)
                        p->C_SML = std::clamp(p->C_SML, C_SML_min, C_SML_max);
                     else if (closest_q->isNeumannSurface)
                        p->C_SML = std::clamp(p->C_SML, C_SML_min, C_SML_max);
                  } else
                     p->C_SML = std::clamp(p->C_SML, C_SML_max, C_SML_max);

                  if (closest_q_next != nullptr) {
                     if (closest_q_next->isSurface)
                        p->C_SML_next = std::clamp(p->C_SML_next, C_SML_min, C_SML_max);
                     else if (closest_q_next->isNeumannSurface)
                        p->C_SML_next = std::clamp(p->C_SML_next, C_SML_min, C_SML_max);
                  } else
                     p->C_SML_next = std::clamp(p->C_SML_next, C_SML_max, C_SML_max);
               }
            }
         }
      } else {
#pragma omp parallel
         for (const auto &p : NET->getPoints())
#pragma omp single nowait
            if (p->isCaptured) {
               auto markerX = p->X + 2 * p->v_to_surface_SPH, markerX_next = X_next(p) + 2 * p->v_to_surface_SPH;
               double closest_d = 1E+20, closest_d_next = 1E+20, d;
               networkPoint *closest_q = nullptr, *closest_q_next = nullptr;
               for (const auto &net : target_nets) {
                  net->BucketPoints.apply(markerX, 3.5 * p->particle_spacing, [&](const auto &q) {
                     if ((q->isSurface || q->isNeumannSurface) && (d = Norm(q->X - markerX)) < closest_d) {
                        closest_d = d;
                        closest_q = q;
                        p->C_SML = closest_d / p->particle_spacing + 0.5;
                     }
                     if ((q->isSurface || q->isNeumannSurface) && (d = Norm(X_next(q) - markerX_next)) < closest_d_next) {
                        closest_d_next = d;
                        closest_q_next = q;
                        p->C_SML_next = closest_d_next / p->particle_spacing + 0.5;
                     }
                  });
               }
               if (closest_q != nullptr) {
                  if (closest_q->isSurface)
                     p->C_SML = std::clamp(p->C_SML, C_SML_min_rigid, C_SML_max);
                  else if (closest_q->isNeumannSurface)
                     p->C_SML = std::clamp(p->C_SML, C_SML_min_rigid, C_SML_max);
               } else
                  p->C_SML = std::clamp(p->C_SML, C_SML_min_rigid, C_SML_max);
               if (closest_q_next != nullptr) {
                  if (closest_q_next->isSurface)
                     p->C_SML_next = std::clamp(p->C_SML_next, C_SML_min_rigid, C_SML_max);
                  else if (closest_q_next->isNeumannSurface)
                     p->C_SML_next = std::clamp(p->C_SML_next, C_SML_min_rigid, C_SML_max);
               } else
                  p->C_SML_next = std::clamp(p->C_SML_next, C_SML_max, C_SML_max);
            }
      }
   DebugPrint("setSML done", Yellow);
};

/*DOC_EXTRACT 0_1_0_SPH

### CFL条件の設定

$`\max({\bf u}) \Delta t \leq c_{v} h \cap \max({\bf a}) \Delta t^2 \leq c_{a} h`$
を満たすように，毎時刻$`\Delta t`$を設定する．
$`c_v=0.1,c_a=0.1`$としている．

*/

double dt_CFL(const double dt_IN, const auto &net, const auto &RigidBodyObject) {
   double dt = dt_IN;
   const auto C_CFL_velocity = 0.1;  // dt = C_CFL_velocity*h/Max(U)
   const auto C_CFL_accel = 0.1;     // dt = C_CFL_accel*sqrt(h/Max(A))
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
            max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH);
            if (dt > max_dt_vel && isFinite(max_dt_vel))
               dt = max_dt_vel;
            /* ------------------------------------------------ */
            // 相対速度
            // double max_dt_acc = C_CFL_accel * std::sqrt(distance / std::abs(Dot(p->DUDt_SPH - q->DUDt_SPH, pq)));
            double max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH - q->DUDt_SPH));
            if (dt > max_dt_acc && isFinite(max_dt_acc))
               dt = max_dt_acc;
            // 絶対速度
            max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH));
            if (dt > max_dt_acc && isFinite(max_dt_acc))
               dt = max_dt_acc;
         }
      };
      net->BucketPoints.apply(p->X, p->SML(), dt_C_CFL);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->SML(), dt_C_CFL);
      double max_dt_vel = C_CFL_velocity * (p->particle_spacing) / Norm(p->U_SPH);
      if (dt > max_dt_vel && isFinite(max_dt_vel))
         dt = max_dt_vel;
      double max_dt_acc = C_CFL_accel * std::sqrt((p->particle_spacing) / Norm(p->DUDt_SPH));
      if (dt > max_dt_acc && isFinite(max_dt_acc))
         dt = max_dt_acc;
   }
   return dt;
}

/* -------------------------------------------------------------------------- */
Tddd U_next(const networkPoint *p) {
   try {
      if (p->isAuxiliary) {
         //! aux pointの動きはsurface pointに合わせる
         return p->surfacePoint->RK_U.getX(p->surfacePoint->DUDt_SPH);
      } else
         return p->RK_U.getX(p->DUDt_SPH);
      // return p->RK_U.getX();
   } catch (...) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "U_next");
   }
}

Tddd X_next(const networkPoint *p) {
   try {
      //! aux pointの位置は，忘れずに正しくRK_Xにセットすること
      if (p->isAuxiliary) {
         auto sp = p->surfacePoint;
#ifdef SET_AUX_AT_PARTICLE_SPACING
         return X_next(sp) + sp->particle_spacing * Normalize(Dot(sp->inv_grad_corr_M, sp->interp_normal_original));
#elif defined(SET_AUX_AT_MASS_CENTER)
         return X_next(sp) - sp->vec2COM_next;
#else
         return p->RK_X.getX(U_next(p));
#endif
      } else if (!p->isFluid)
         return MOVE_WALL_PARTICLE ? p->RK_X.getX(U_next(p)) : p->X;
      else
         return p->RK_X.getX(U_next(p));
   } catch (...) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "X_next");
   }
}

// \label{SPH:rho_next}
double rho_next(auto p) {
   // if (p->getNetwork()->isRigidBody)
   // if (p->getNetwork()->isRigidBody && !p->isFirstWallLayer)
   return _WATER_DENSITY_;
   /* -------------------------------------------------------------------------- */
   //    if (p->isAuxiliary)
   //       return rho_next(p->surfacePoint);
   //    else if (p->getNetwork()->isRigidBody)
   //       return _WATER_DENSITY_;
   //    else {
   // #if defined(USE_RungeKutta)
   // return p->RK_rho.getX(-p->rho * p->div_U);

   //@ これを使った方が安定するようだ
   // #elif defined(USE_LeapFrog)
   // return p->rho + p->DrhoDt_SPH + p->RK_rho.get_dt();
   // #endif
   //    }
};

// \label{SPH:volume_next}
double V_next(const auto &p) { return p->mass / rho_next(p); };

// \label{SPH:position_next}
/* -------------------------------------------------------------------------- */

std::tuple<networkPoint *, Tddd> closest(const networkPoint *p, const auto &RigidBodyObject) {
   double closest_d = 1E+20, d;
   Tddd p_to_q, closest_p_to_q;
   networkPoint *closest_q = nullptr;
   // for (const auto &[obj, _] : RigidBodyObject) {
   for (const auto &obj : RigidBodyObject) {
      obj->BucketPoints.apply(p->X, p->SML(), [&](const auto &q) {
         p_to_q = q->X - p->X;
         if (closest_d > (d = Norm(p_to_q))) {
            closest_d = d;
            closest_p_to_q = p_to_q;
            closest_q = q;
         }
      });
   }
   return {closest_q, closest_p_to_q};
};

std::tuple<networkPoint *, Tddd> closest_next(const networkPoint *p, const auto &RigidBodyObject) {
   double closest_d = 1E+20, d;
   Tddd p_to_q, closest_p_to_q;
   networkPoint *closest_q = nullptr;
   // for (const auto &[obj, _] : RigidBodyObject) {
   for (const auto &obj : RigidBodyObject) {
      obj->BucketPoints.apply(X_next(p), 1.1 * p->SML(), [&](const auto &q) {
         p_to_q = X_next(q) - X_next(p);
         if (closest_d > (d = Norm(p_to_q))) {
            closest_d = d;
            closest_p_to_q = p_to_q;
            closest_q = q;
         }
      });
   }
   return {closest_q, closest_p_to_q};
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
//@ -------------------------------------------------------- */
//@                        粒子の時間発展                      */
//@ -------------------------------------------------------- */

#define REFLECTION

void updateParticles(const auto &points,
                     const std::unordered_set<Network *> &target_nets,
                     // const std::vector<std::tuple<Network *, Network *>> &RigidBodyObject,
                     const auto &RigidBodyObject,
                     const double &particle_spacing) {
   try {
      DebugPrint("粒子の時間発展", Green);

#pragma omp parallel
      for (const auto &p : points)
#pragma omp single nowait
      {
         auto U = p->U_SPH;
         auto X_last = p->X;
#if defined(USE_RungeKutta)
         // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
         p->RK_U.push(p->DUDt_SPH);  // 速度
         p->U_SPH = p->RK_U.getX();
         p->RK_X.push(p->U_SPH);  // 位置
         p->setXSingle(p->RK_X.getX());
#elif defined(USE_LeapFrog)
         auto X_next = [&](const auto &p) { return p->X; };
         p->DUDt_modify_SPH.fill(0.);
         p->LPFG_X.push(p->DUDt_SPH);  // 速度
         p->U_SPH = p->LPFG_X.get_v();
         p->setXSingle(p->tmp_X = p->LPFG_X.get_x());
            // auto getX = [&](const auto &p) { return p->X; };
#endif

#if defined(REFLECTION)
         int count = 0;
         //\label{SPH:reflection}
         const double reflection_factor = .1;
         const double c = 0.;
         bool isReflected = true;
         p->DUDt_modify_SPH.fill(0.);
         /* --------------------------------------------------------- */
         // using apply of RigidBodyObject calculate reflection
         //       for (const auto &[rigid, _] : RigidBodyObject) {
         //          rigid->BucketPoints.apply(X_next(p), 1.5 * particle_spacing, [&](const auto &Pwall) {
         //             auto v2wall_next = X_next(Pwall) - X_next(p);
         //             auto v2wall = Pwall->X - p->X;
         //             // auto dist = Distance(X_next(Pwall), X_next(p));
         //             auto dist_next = Norm(v2wall_next);
         //             auto dist = Norm(v2wall);
         //             auto a = 0.5;
         //             // auto mid_dist = Distance(X_next(Pwall), ((1. - a) * X_next(p) + a * p->X));
         //             auto n = Normalize(Pwall->interp_normal_original);
         //             // if (Dot(p->U_SPH, n) < 0 && dist < particle_spacing) {
         //             if (dist < 1.1 * particle_spacing || dist_next < 1.1 * particle_spacing) {
         //                auto tmp = -5E-3 * Projection(p->U_SPH, n) / p->RK_X.get_dt();
         //                p->DUDt_modify_SPH += tmp;
         //                p->DUDt_SPH += tmp;
         // #if defined(USE_RungeKutta)
         //                p->RK_U.repush(p->DUDt_SPH);
         //                p->U_SPH = p->RK_U.getX();
         //                isReflected = true;
         // #elif defined(USE_LeapFrog)
         //                          p->LPFG_X.repush(p->DUDt_SPH);  // 速度
         //                          p->U_SPH = p->LPFG_X.get_v();
         //                          p->setXSingle(p->tmp_X = p->LPFG_X.get_x());
         //                          isReflected = true;
         // #endif
         //             }
         //          });
         //       }
         /* --------------------------------------------------------- */
         while (isReflected && count++ < 1) {
            isReflected = false;
            auto d_ps = particle_spacing;
            auto d0 = (1 - c) * particle_spacing;
            for (const auto &[closest_p, v_f2w] : {/*closest(p, RigidBodyObject), */ closest_next(p, RigidBodyObject)}) {
               if (closest_p != nullptr) {
                  auto n = Normalize(closest_p->interp_normal_original);
                  auto n_d_f2w = Norm(Projection(v_f2w, n));
                  auto ratio = (d0 - n_d_f2w) / d0;
                  if (Norm(v_f2w) < 1. * d0 && Norm(closest_p->X - p->X) < 1. * d0) {
                     // auto ratio = (d0 - n_d_f2w) / d0;
                     if (Dot(p->U_SPH, n) < 0) {
                        auto tmp = -0.05 * Projection(p->U_SPH, n) / p->RK_X.get_dt();
                        tmp += 0.5 * _GRAVITY_ * n;
                        // auto tmp = -0.02 * Projection(p->U_SPH, n) / p->RK_X.get_dt();
                        // auto tmp = -0.01 * Projection(p->U_SPH, n) / p->RK_X.get_dt();
                        p->DUDt_modify_SPH += tmp;
                        p->DUDt_SPH += tmp;
   #if defined(USE_RungeKutta)
                        // p->RK_U.repush(p->DUDt_SPH);  // 速度
                        // p->U_SPH = p->RK_U.getX();
                        p->RK_U.repush(p->DUDt_SPH);  // 速度
                        p->U_SPH = p->RK_U.getX();
                        p->RK_X.repush(p->U_SPH);  // 位置
                        p->setXSingle(p->RK_X.getX());
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
            }
         };
            /* ------------------------------------------------------- */
#endif
      }

      // \label{SPH:update_density}
      for (const auto &A : points) {
#if defined(USE_RungeKutta)
         A->DrhoDt_SPH = -A->rho * A->div_U;
         A->RK_rho.push(A->DrhoDt_SPH);  // 密度

         // if (A->RK_rho.finished)
         A->setDensity(_WATER_DENSITY_);
            // else
            // A->setDensity(A->RK_rho.get_x());
            // A->setDensity((A->RK_rho.get_x() + _WATER_DENSITY_) / 2.);
            // A->setDensity(A->RK_rho.get_x());
            // A->setDensity(_WATER_DENSITY_);
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
      }
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in updateParticles");
   };
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