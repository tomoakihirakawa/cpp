#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/* -------------------------------------------------------------------------- */

void showParticleInfo(const networkPoint* A) {
   // show details
   std::cout << "A->DUDt_SPH = " << A->DUDt_SPH << std::endl;
   std::cout << "A->DUDt_modify_SPH = " << A->DUDt_modify_SPH << std::endl;
   std::cout << "A->mu_SPH = " << A->mu_SPH << std::endl;
   std::cout << "A->rho = " << A->rho << std::endl;
   std::cout << "A->lap_U = " << A->lap_U << std::endl;
   std::cout << "A->U_XSPH = " << A->U_XSPH << std::endl;
   std::cout << "A->U_SPH = " << A->U_SPH << std::endl;
   std::cout << "A->b_vector = " << A->b_vector << std::endl;
   std::cout << "A->DrhoDt_SPH = " << A->DrhoDt_SPH << std::endl;
   std::cout << "A->div_U = " << A->div_U << std::endl;
   std::cout << "A->laplacian_corr_M = " << A->laplacian_corr_M << std::endl;
   std::cout << "A->laplacian_corr_M_next = " << A->laplacian_corr_M_next << std::endl;
   std::cout << "A->grad_corr_M = " << A->grad_corr_M << std::endl;
   std::cout << "A->grad_corr_M_next = " << A->grad_corr_M_next << std::endl;
   std::cout << "A->inv_grad_corr_M = " << A->inv_grad_corr_M << std::endl;
   std::cout << "A->inv_grad_corr_M_next = " << A->inv_grad_corr_M_next << std::endl;
   std::cout << "A->nabla_otimes_U = " << A->nabla_otimes_U << std::endl;
   std::cout << "A->C_SML = " << A->C_SML << std::endl;
   std::cout << "A->C_SML_next = " << A->C_SML_next << std::endl;
   std::cout << "A->isFluid = " << A->isFluid << std::endl;
   std::cout << "A->isSurface = " << A->isSurface << std::endl;
   std::cout << "A->isSurface_next = " << A->isSurface_next << std::endl;
   std::cout << "A->isNeumannSurface = " << A->isNeumannSurface << std::endl;
   std::cout << "A->isAuxiliary = " << A->isAuxiliary << std::endl;
   std::cout << "A->isCaptured = " << A->isCaptured << std::endl;
   std::cout << "A->isNearSurface = " << A->isNearSurface << std::endl;
   std::cout << "A->particle_spacing = " << A->particle_spacing << std::endl;
   std::cout << "A->volume = " << A->volume << std::endl;
   std::cout << "A->mass = " << A->mass << std::endl;
   std::cout << "A->Eigenvalues_of_M1 = " << A->Eigenvalues_of_M1 << std::endl;
   std::cout << "A->Eigenvectors_of_M1 = " << A->Eigenvectors_of_M1 << std::endl;
}

bool canInteract(const networkPoint* A, const networkPoint* B) {
   return B->isCaptured;
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

/* -------------------------------- C_SMLの調整 -------------------------------- */

void setSML(const auto& target_nets, const std::array<double, 2> C_SML_min_max = { 2.7, 2.7 }) {
   DebugPrint("setSML", Yellow);
   const auto [C_SML_min, C_SML_max] = C_SML_min_max;
   const auto [C_SML_min_rigid, C_SML_max_rigid] = C_SML_min_max;

   //! 初期化
   for (const auto& NET : target_nets)
      for (const auto& p : NET->getPoints())
      {
         p->min_particle_spacing = p->particle_spacing;
         p->min_particle_spacing_next = p->particle_spacing;
      }

   for (const auto& NET : target_nets)
      if (NET->isFluid) {
#pragma omp parallel
         for (const auto& p : NET->getPoints())
#pragma omp single nowait
         {
            p->isNearSurface = false;
            if (p->isCaptured) {
               double closest_d = 1E+20, closest_d_next = 1E+20, d;
               networkPoint* closest_q = nullptr, * closest_q_next = nullptr;
               for (const auto& net : target_nets) {
                  net->BucketPoints.apply(p->X, 3.5 * p->particle_spacing, [&](const auto& q) {
                     if ((q->isSurface /* || q->isNeumannSurface*/) && (d = Norm(q->X - p->X)) < closest_d) {
                        closest_d = d;
                        closest_q = q;
                        p->C_SML = closest_d / p->particle_spacing + 1;
                     }
                     if ((q->isSurface /* || q->isNeumannSurface*/) && (d = Norm(X_next(q) - X_next(p))) < closest_d_next) {
                        closest_d_next = d;
                        closest_q_next = q;
                        p->C_SML_next = closest_d_next / p->particle_spacing + 1;
                     }
                     if (q->isSurface && Distance(p, q) < 1.5 * p->particle_spacing)
                        p->isNearSurface = true;

                     if (p != q) {
                        p->min_particle_spacing = std::min(p->min_particle_spacing, Norm(q->X - p->X));
                        p->min_particle_spacing_next = std::min(p->min_particle_spacing_next, Norm(X_next(q) - X_next(p)));
                     }
                     });
               }
               if (closest_q != nullptr) {
                  if (closest_q->isSurface)
                     p->C_SML = std::clamp(p->C_SML, C_SML_min, C_SML_max);
                  // else if (closest_q->isNeumannSurface)
                  //    p->C_SML = std::clamp(p->C_SML, C_SML_min, C_SML_max);
               }
               else
                  p->C_SML = std::clamp(p->C_SML, C_SML_max, C_SML_max);

               if (closest_q_next != nullptr) {
                  if (closest_q_next->isSurface)
                     p->C_SML_next = std::clamp(p->C_SML_next, C_SML_min, C_SML_max);
                  // else if (closest_q_next->isNeumannSurface)
                  //    p->C_SML_next = std::clamp(p->C_SML_next, C_SML_min, C_SML_max);
               }
               else
                  p->C_SML_next = std::clamp(p->C_SML_next, C_SML_max, C_SML_max);
            }
         }
      }
      else {
#pragma omp parallel
         for (const auto& p : NET->getPoints())
#pragma omp single nowait
            if (p->isCaptured) {
               auto markerX = p->X + 2 * p->v_to_surface_SPH, markerX_next = X_next(p) + 2 * p->v_to_surface_SPH;
               double closest_d = 1E+20, closest_d_next = 1E+20, d;
               networkPoint* closest_q = nullptr, * closest_q_next = nullptr;
               for (const auto& net : target_nets) {
                  net->BucketPoints.apply(markerX, 3.5 * p->particle_spacing, [&](const auto& q) {
                     if ((q->isSurface /* || q->isNeumannSurface*/) && (d = Norm(q->X - markerX)) < closest_d) {
                        closest_d = d;
                        closest_q = q;
                        p->C_SML = closest_d / p->particle_spacing + 1;
                     }
                     if ((q->isSurface /* || q->isNeumannSurface*/) && (d = Norm(X_next(q) - markerX_next)) < closest_d_next) {
                        closest_d_next = d;
                        closest_q_next = q;
                        p->C_SML_next = closest_d_next / p->particle_spacing + 1;
                     }
                     });
               }
               if (closest_q != nullptr) {
                  if (closest_q->isSurface)
                     p->C_SML = std::clamp(p->C_SML, C_SML_min_rigid, C_SML_max_rigid);
                  // else if (closest_q->isNeumannSurface)
                  //    p->C_SML = std::clamp(p->C_SML, C_SML_min_rigid, C_SML_max);
               }
               else
                  p->C_SML = std::clamp(p->C_SML, C_SML_min_rigid, C_SML_max_rigid);
               if (closest_q_next != nullptr) {
                  if (closest_q_next->isSurface)
                     p->C_SML_next = std::clamp(p->C_SML_next, C_SML_min_rigid, C_SML_max_rigid);
                  // else if (closest_q_next->isNeumannSurface)
                  //    p->C_SML_next = std::clamp(p->C_SML_next, C_SML_min_rigid, C_SML_max);
               }
               else
                  p->C_SML_next = std::clamp(p->C_SML_next, C_SML_max_rigid, C_SML_max_rigid);
            }
      }

   //! fluidのpoints_in_SMLを更新
   for (const auto& NET : target_nets)
      if (NET->isFluid) {
#pragma omp parallel
         for (const auto& p : NET->getPoints())
#pragma omp single nowait
         {
            p->points_in_SML.clear();
            for (const auto& net : target_nets)
               net->BucketPoints.apply(p->X, p->SML_grad(), [&](const auto& q) {
               if (Norm(q->X - p->X) < p->SML_grad() && p != q && canInteract(p, q))
                  p->points_in_SML.emplace_back(q);
                  });
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

double dt_CFL(const double dt_IN, const auto& net, const auto& RigidBodyObject) {
   // double dt = dt_IN;
   const auto C_CFL_velocity = 0.1;           // dt = C_CFL_velocity*h/Max(U)
   const auto C_CFL_accel = 0.1;              // dt = C_CFL_accel*sqrt(h/Max(A))
   const auto C_CFL_viscosity = 0.1;          // dt = C_CFL_viscosity*h^2/Max(Viscosity)
   const auto C_CFL_velocity_relative = 0.1;  // dt = C_CFL_velocity*h/Max(U)
   const auto C_CFL_accel_relative = 0.1;     // dt = C_CFL_accel*sqrt(h/Max(A))j
   // 音速
#pragma omp parallel
   for (const auto& p : net->getPoints())
#pragma omp single nowait
   {
      p->dt_CFL = dt_IN;  // 速度に関するCFL条件
      auto dt_C_CFL = [&](const auto& q) {
         for (double distance : V_d{ Distance(p, q), p->particle_spacing })
            if (p != q) {
               auto pq = Normalize(p->X - q->X);
               double max_dt_vel, max_dt_acc;
               /* ------------------------------------------------ */
               // 音速
               // max_dt_vel = C_CFL_velocity * distance / Cs;
               // if (p->dt_CFL > max_dt_vel && isFinite(max_dt_vel))
               //    p->dt_CFL = max_dt_vel;
               /* ------------------------------------------------ */
               // 相対速度
               max_dt_vel = C_CFL_velocity_relative * distance / Norm(p->U_SPH - q->U_SPH);
               if (p->dt_CFL > max_dt_vel && isFinite(max_dt_vel))
                  p->dt_CFL = max_dt_vel;
               // 絶対速度
               max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH);
               if (p->dt_CFL > max_dt_vel && isFinite(max_dt_vel))
                  p->dt_CFL = max_dt_vel;
               /* ------------------------------------------------ */
               // 相対速度
               max_dt_acc = C_CFL_accel_relative * std::sqrt(distance / Norm(p->DUDt_SPH - q->DUDt_SPH));
               if (p->dt_CFL > max_dt_acc && isFinite(max_dt_acc))
                  p->dt_CFL = max_dt_acc;
               // double max_dt_acc = C_CFL_accel_relative * std::sqrt(distance / (Norm(p->DUDt_update_SPH) + Norm(q->DUDt_update_SPH)));
               // if (p->dt_CFL > max_dt_acc && isFinite(max_dt_acc))
               //    p->dt_CFL = max_dt_acc;
               // 絶対速度
               max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH));
               if (p->dt_CFL > max_dt_acc && isFinite(max_dt_acc))
                  p->dt_CFL = max_dt_acc;
               /* -------------------------------------------------------------------------- */
               double max_dt_visc = distance * distance / std::abs(8. * p->mu_SPH / _WATER_DENSITY_);
               if (p->dt_CFL > max_dt_visc && isFinite(max_dt_visc))
                  p->dt_CFL = max_dt_visc;
            }
         };

      net->BucketPoints.apply(p->X, p->SML_grad(), dt_C_CFL);
      for (const auto& [obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->SML_grad(), dt_C_CFL);

      // double max_dt_vel = C_CFL_velocity * (p->particle_spacing) / Norm(p->U_SPH);
      // if (dt > max_dt_vel && isFinite(max_dt_vel))
      //    dt = max_dt_vel;
      // double max_dt_acc = C_CFL_accel * std::sqrt((p->particle_spacing) / Norm(p->DUDt_SPH));
      // if (dt > max_dt_acc && isFinite(max_dt_acc))
      //    dt = max_dt_acc;
   }

   // find smallest dt
   double dt = dt_IN;
   for (const auto& p : net->getPoints())
      if (p->dt_CFL < dt_IN)
         dt = p->dt_CFL;

   return dt;
}

#include "SPH_Auxiliary.hpp"

/* -------------------------------------------------------------------------- */
Tddd U_next(const networkPoint* p) {
   return p->RK_U.getX(p->DUDt_SPH);
}

Tddd X_next(const networkPoint* p) {
   if (!p->isFluid)
      return p->X;
   else
      return p->RK_X.getX(U_next(p));
}

// \label{SPH:rho_next}
double rho_next(const networkPoint* p) {
   // if (!p->isFluid)
   //    return _WATER_DENSITY_;
   // else
   // {
      // return p->RK_rho.getX(-p->rho * p->div_U);
      // return 0.01 * p->RK_rho.getX(p->DrhoDt_SPH) + 0.99 * _WATER_DENSITY_;

   auto tmp = -p->rho * p->div_U;
   return p->RK_rho.getX(tmp / std::sqrt(1. + 0.01 * std::abs(tmp)));

   // if (tmp > 0)
   //    return p->RK_rho.getX(std::log(tmp + 1.));
   // else
   //    return p->RK_rho.getX(-std::log(-tmp + 1.));

   // const double a = 0.6;
   // return a * p->RK_rho.getX(-p->rho * p->div_U_next) + (1. - a) * _WATER_DENSITY_;
// }
};

// \label{SPH:volume_next}
double V_next(const networkPoint* p) {
   return p->mass / rho_next(p);
};

// \label{SPH:position_next}
/* -------------------------------------------------------------------------- */

std::vector<networkPoint*> getNearPoints(const networkPoint* p, const auto& RigidBodyObject, double range = -1.) {
   double closest_d = 1E+20, d;
   Tddd p_to_q, closest_p_to_q;
   networkPoint* closest_q = nullptr;
   if (range < 0.)
      range = p->SML();
   std::vector<networkPoint*> nearPoints;
   // for (const auto &[obj, _] : RigidBodyObject) {
   for (const auto& obj : RigidBodyObject) {
      obj->BucketPoints.apply(p->X, range, [&](const auto& q) {
         p_to_q = q->X - p->X;
         if (Norm(p_to_q) < range)
            nearPoints.emplace_back(q);
         });
   }
   return nearPoints;
};

std::tuple<networkPoint*, Tddd> findClosest(const networkPoint* p, const std::vector<networkPoint*>& nearPoints) {
   double closest_d = 1E+20, d;
   Tddd p_to_q, closest_p_to_q;
   networkPoint* closest_q = nullptr;
   for (const auto& q : nearPoints) {
      p_to_q = q->X - p->X;
      if (closest_d > (d = Norm(p_to_q))) {
         closest_d = d;
         closest_p_to_q = p_to_q;
         closest_q = q;
      }
   };
   return { closest_q, closest_p_to_q };
};

std::tuple<networkPoint*, Tddd> closest(const networkPoint* p, const auto& RigidBodyObject, double range = -1.) {
   double closest_d = 1E+20, d;
   Tddd p_to_q, closest_p_to_q;
   networkPoint* closest_q = nullptr;
   if (range < 0.)
      range = p->SML();
   // for (const auto &[obj, _] : RigidBodyObject) {
   for (const auto& obj : RigidBodyObject) {
      obj->BucketPoints.apply(p->X, range, [&](const auto& q) {
         p_to_q = q->X - p->X;
         if (closest_d > (d = Norm(p_to_q))) {
            closest_d = d;
            closest_p_to_q = p_to_q;
            closest_q = q;
         }
         });
   }
   return { closest_q, closest_p_to_q };
};

std::tuple<networkPoint*, Tddd> closest_next(const networkPoint* p, const auto& RigidBodyObject) {
   double closest_d = 1E+20, d;
   Tddd p_to_q, closest_p_to_q;
   networkPoint* closest_q = nullptr;
   // for (const auto &[obj, _] : RigidBodyObject) {
   for (const auto& obj : RigidBodyObject) {
      obj->BucketPoints.apply(X_next(p), 1.1 * p->SML(), [&](const auto& q) {
         p_to_q = X_next(q) - X_next(p);
         if (closest_d > (d = Norm(p_to_q))) {
            closest_d = d;
            closest_p_to_q = p_to_q;
            closest_q = q;
         }
         });
   }
   return { closest_q, closest_p_to_q };
};

/* -------------------------------------------------------------------------- */
#include "SPH0_setWall_Freesurface.hpp"
//
#include "SPH1_lap_div_U.hpp"
//
#include "SPH2_FindPressure.hpp"
//
#include "SPH3_grad_P.hpp"
//
//@ -------------------------------------------------------- */
//@                        粒子の時間発展                      */
//@ -------------------------------------------------------- */

//! calculate p->nabla_otimes_U using the velocity at next time step
//! nabla_otimes_U is TensorProduct of velocity gradient
void calculate_nabla_otimes_U_next(const auto& points, const auto& target_nets) {
   try {
      DebugPrint("calculate_nabla_otimes_U_next", Green);
#pragma omp parallel
      for (const auto& p : points)
#pragma omp single nowait
      {
         T3Tddd nabla_otimes_U;
         // p->U_XSPH = U_next(p);
         p->U_XSPH.fill(0.);
         Fill(nabla_otimes_U, 0.);
         const double c_xsph = 0.01;

         auto add = [&](const auto& q) {
            auto grad = grad_w_Bspline(p, q);
            auto pU = p->RK_U.getX(p->DUDt_SPH);
            auto qU = q->RK_U.getX(q->DUDt_SPH);
            if (p->isFluid && q->isFluid && !q->isAuxiliary && q != p)
               p->U_XSPH -= c_xsph * (U_next(p) - U_next(q)) * V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), p->SML_next());  //\label{SPH:U_XSPH}
            nabla_otimes_U += q->volume * TensorProduct(grad, qU - pU);
            };

         for (const auto& net : target_nets)
            net->BucketPoints.apply(p->X, p->SML_grad(), add);

         p->nabla_otimes_U = nabla_otimes_U;
      }
   }
   catch (...) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in calculate_nabla_otimes_U_next");
   }
}

#define REFLECTION
// #define USE_ALE
void updateParticles(const auto& net,
   const std::unordered_set<Network*>& target_nets,
   // const std::vector<std::tuple<Network *, Network *>> &RigidBodyObject,
   const auto& RigidBodyObject,
   const double& particle_spacing) {
   try {
      DebugPrint("粒子の時間発展", Green);
      auto points = net->getPoints();
#pragma omp parallel
      for (const auto& p : points)
#pragma omp single nowait
      {
         auto U = p->U_SPH;
         auto X_last = p->X;
         const double dt = p->RK_X.getdt();
         p->RK_U.push(p->DUDt_SPH);
         p->U_SPH = p->RK_U.getX();  // + p->U_XSPH;
         auto U_SPH_original = p->U_SPH;
         p->RK_X.push(p->U_SPH);  // 位置の更新はU_XSPHを使う
         p->setXSingle(p->RK_X.getX());
         p->DUDt_update_SPH = p->DUDt_SPH;
#if defined(REFLECTION)
         int count = 0;
         //\label{SPH:reflection}
         const double reflection_factor = .1;
         const double c = 0.;
         bool isReflected = true;
         p->DUDt_modify_SPH.fill(0.);
         // var_M1が大きすぎる場合は計算を修正しよう．
         static const double N = 10000.;
         static const double c_reflection = 1. / N;
         const double d_ps = particle_spacing;
         const double d0 = (1 - c) * particle_spacing;
         auto vector_points = getNearPoints(p, RigidBodyObject, 1.5 * particle_spacing);
         while (isReflected && count++ < N) {  // 厳しく
            isReflected = false;
            for (const auto& [closest_p, f2w] : { findClosest(p, vector_points) /* , closest_next(p, RigidBodyObject)*/ }) {
               if (closest_p != nullptr) {
                  auto dist_f2w = Norm(f2w);
                  // auto n = closest_p->interp_normal;
                  auto n = closest_p->normal_SPH;
                  // auto n = Normalize(closest_p->interp_normal_original);
                  auto normal_dist_f2w = Norm(Projection(f2w, n));  //! nomal distance from the fluid particle to the wall particle
                  auto ratio = 1. - normal_dist_f2w / d0;
                  // if (Norm(f2w) < 1. * d0 && Norm(closest_p->X - p->X) < 1. * d0) {
                  if (dist_f2w < particle_spacing) {
                     if (normal_dist_f2w < 0.8 * particle_spacing && Dot(p->U_SPH, n) < 0) {
                        p->DUDt_modify_SPH -= c_reflection * Projection(U_SPH_original, n) / dt;
                        p->DUDt_update_SPH = p->DUDt_SPH + p->DUDt_modify_SPH;
                        p->RK_U.repush(p->DUDt_update_SPH);  // 速度
                        p->U_SPH = p->RK_U.getX();           // + p->U_XSPH;
                        p->RK_X.repush(p->U_SPH);            // 位置
                        p->setXSingle(p->RK_X.getX());
                        isReflected = true;
                     }
                  }
               }
            }
         };
         /* ------------------------------------------------------- */
#endif
      }

      // #pragma omp parallel
      //       for (const auto& p : points)
      // #pragma omp single nowait
      //       {
      //          double SML = p->SML();
      //          double total_volume = 0.;

      //          p->intp_density = 0.;
      //          net->BucketPoints.apply(p->X, 1.1 * SML, [&](const auto& q) {
      //             double w = q->volume * w_Bspline(Norm(p->X - q->X), SML);
      //             p->intp_density += q->density * w;
      //             total_volume += w;
      //             });
      //          if (total_volume > 0.)
      //             p->intp_density /= total_volume;
      //          else
      //             p->intp_density = _WATER_DENSITY_;

      //          p->intp_density_next = p->intp_density;

      //       }

            /* -------------------------------------------------------------------------- */
            // % div_U の計算
            //       setCorrectionMatrix(std::unordered_set<Network*>{net});
            // #pragma omp parallel
            //       for (const auto& A : points)
            // #pragma omp single nowait
            //       {
            //          A->div_U = 0.;
            //          for (const auto& net : target_nets) {
            //             net->BucketPoints.apply(A->X, A->SML_grad(), [&](const auto& B) {
            //                if (canInteract(A, B))
            //                   FusedMultiplyIncrement(-B->volume, Dot(A->U_SPH - B->U_SPH, grad_w_Bspline(A, B)), A->div_U);
            //                });
            //          }
            //          A->DrhoDt_SPH = -A->rho * A->div_U;
            //       }

            //% 密度の更新
            // \label{SPH:update_density}
            // const double a = 0.5;
      for (const auto& A : points) {
         double SML = A->SML();
         double total_volume = 0.;

         A->intp_density = 0.;
         for (const auto& net : target_nets)
            net->BucketPoints.apply(A->X, 1.1 * SML, [&](const auto& q) {
            double w = q->volume * w_Bspline(Norm(A->X - q->X), SML);
            A->intp_density += q->density * w;
            total_volume += w;
               });
         A->intp_density;
         A->intp_density_next = A->intp_density;
      }

      for (const auto& A : points) {
         const double a = 0.9;
         A->setDensity(0.01 * A->intp_density + 0.99 * (a * rho_next(A) + (1. - a) * _WATER_DENSITY_));
         A->RK_rho.push(A->DrhoDt_SPH);
      }
      /* -------------------------------------------------------------------------- */

   }
   catch (std::exception& e) {
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