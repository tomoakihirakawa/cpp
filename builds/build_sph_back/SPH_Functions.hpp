#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

void setSML(const auto &target_nets) {
   /* -------------------------------- C_SMLの調整 -------------------------------- */
   double C_SML_max = 3.;
   // double C_SML_min = 1.866;
   double C_SML_min = 2.;
   for (const auto &NET : target_nets)
      if (NET->isFluid) {
         {
#pragma omp parallel
            for (const auto &p : NET->getPoints())
#pragma omp single nowait
               if (p->isCaptured) {
                  double closest_d = 1E+20, closest_d_next = 1E+20, d;
                  networkPoint *closest_q = nullptr, *closest_q_next = nullptr;
                  for (const auto &net : target_nets) {
                     net->BucketPoints.apply(p->X, 3.5 * p->particle_spacing, [&](const auto &q) {
                        if ((q->isSurface || q->isNeumannSurface) && (d = Norm(q->X - p->X)) < closest_d) {
                           closest_d = d;
                           closest_q = q;
                           p->C_SML = closest_d / p->particle_spacing;
                        }
                        if ((q->isSurface || q->isNeumannSurface) && (d = Norm(X_next(q) - X_next(p))) < closest_d_next) {
                           closest_d_next = d;
                           closest_q_next = q;
                           p->C_SML_next = closest_d_next / p->particle_spacing;
                        }
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
                        p->C_SML = closest_d / p->particle_spacing;
                     }
                     if ((q->isSurface || q->isNeumannSurface) && (d = Norm(X_next(q) - markerX_next)) < closest_d_next) {
                        closest_d_next = d;
                        closest_q_next = q;
                        p->C_SML_next = closest_d_next / p->particle_spacing;
                     }
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
};

networkPoint *getClosestExcludeRigidBodyInlcudeFirstLayer(networkPoint *p, auto &target_nets) {
   double distance = 1E+20;
   networkPoint *P = nullptr;
   for (const auto &obj : target_nets) {
      obj->BucketPoints.apply(p->X, p->SML(), [&](const auto &q) {
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
   obj->BucketPoints.apply(p->X, p->SML() * 1.5, [&](const auto &q) {
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
         obj->BucketPoints.apply(p->X, p->SML() * 1.5, [&](const auto &q) {
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

Tddd aux_position(const networkPoint *p, const double &c) {
   // auto c = p->particle_spaincing;
   return p->X + c * Normalize(p->interp_normal);
   // return p->X - (p->COM_SPH - p->X);
};

Tddd aux_position_next(const networkPoint *p) {
   auto q = p->surfacePoint;
   auto c = q->SML() / q->C_SML;
#if defined(USE_RungeKutta)
   return q->RK_X.getX(q->U_SPH) + c * Normalize(q->interp_normal_next);
#elif defined(USE_LeapFrog)
   return q->LPFG_X.get_x(q->U_SPH) + c * Normalize(q->interp_normal_next);
#endif
};

/* -------------------------------------------------------------------------- */
// \label{SPH:rho_next}
double rho_next(auto p) {
   // if (p->getNetwork()->isRigidBody)
   return _WATER_DENSITY_;
   /* -------------------------------------------------------------------------- */
   //    if (p->isAuxiliary)
   //       return rho_next(p->surfacePoint);
   //    else if (p->getNetwork()->isRigidBody)
   //       return _WATER_DENSITY_;
   //    else {
   // #if defined(USE_RungeKutta)
   // return p->RK_rho.getX(p->DrhoDt_SPH);
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
      // return p->RK_X.getX(p->U_SPH);

      // auto du = (p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_) * p->RK_X.get_dt();
      // return p->RK_X.getX(p->U_SPH + du);
      // return p->RK_X.getX(p->RK_U.getX(p->DUDt_SPH));
      return p->RK_X.getX(p->U_SPH);
#elif defined(USE_LeapFrog)
      return p->X + p->U_SPH * p->LPFG_X.get_dt() / 2.;
#endif
   }
};

/* -------------------------------------------------------------------------- */
double w_Bspline(const networkPoint *p, const networkPoint *q, const double r) {
   return w_Bspline(Norm(p->X - q->X), r);
};

std::array<double, 3> grad_w_Bspline_helper(const networkPoint *p, const Tddd &pX, const Tddd &qX, const double r) {
   // if (p->isFluid && p->isNotSurfaceButNearSurface)
   //    return grad_w_Bspline(pX, qX, r, p->inv_grad_corr_M_mirror);
   // else
   return grad_w_Bspline(pX, qX, r, p->inv_grad_corr_M);
};

std::array<double, 3> grad_w_Bspline(const networkPoint *p, const networkPoint *q, const double r) {
   return grad_w_Bspline_helper(p, p->X, q->X, r);
};

std::array<double, 3> grad_w_Bspline(const networkPoint *p, const networkPoint *q) {
   return grad_w_Bspline_helper(p, p->X, q->X, p->SML());
};

std::array<double, 3> grad_w_Bspline(const std::tuple<networkPoint *, Tddd> &p_X, const networkPoint *q) {
   return grad_w_Bspline_helper(std::get<0>(p_X), std::get<1>(p_X), q->X, std::get<0>(p_X)->SML());
};

std::array<double, 3> grad_w_Bspline(const networkPoint *p, const Tddd &X) {
   return grad_w_Bspline_helper(p, p->X, X, p->SML());
};

double Dot_grad_w_Bspline_helper(const networkPoint *p, const Tddd &pX, const Tddd &qX, const double r) {
   // if (p->isFluid && p->isNotSurfaceButNearSurface)
   //    return Dot_grad_w_Bspline(pX, qX, r, p->inv_grad_corr_M_mirror);
   // else
   return Dot_grad_w_Bspline(pX, qX, r, p->inv_grad_corr_M);
};

double Dot_grad_w_Bspline(const networkPoint *p, const networkPoint *q, const double r) {
   return Dot_grad_w_Bspline_helper(p, p->X, q->X, r);
};

double Dot_grad_w_Bspline(const networkPoint *p, const networkPoint *q) {
   return Dot_grad_w_Bspline_helper(p, p->X, q->X, p->SML());
};

double Dot_grad_w_Bspline(const std::tuple<networkPoint *, Tddd> &p_X, const networkPoint *q) {
   return Dot_grad_w_Bspline_helper(std::get<0>(p_X), std::get<1>(p_X), q->X, std::get<0>(p_X)->SML());
};

double Dot_grad_w_Bspline(const std::tuple<networkPoint *, Tddd> &p_X, const Tddd &Y) {
   return Dot_grad_w_Bspline_helper(std::get<0>(p_X), std::get<1>(p_X), Y, std::get<0>(p_X)->SML());
};

/* -------------------------------------------------------------------------- */
double w_Bspline_next(const networkPoint *p, const networkPoint *q, const double r) {
   return w_Bspline(Norm(X_next(p) - X_next(q)), r);
};

std::array<double, 3> grad_w_Bspline_next_helper(const networkPoint *p, const Tddd &pX, const Tddd &qX, const double r) {
   // if (p->isFluid && p->isNotSurfaceButNearSurface)
   //    return grad_w_Bspline(pX, qX, r, p->inv_grad_corr_M_next_mirror);
   // else
   return grad_w_Bspline(pX, qX, r, p->inv_grad_corr_M_next);
};

std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
   return grad_w_Bspline_next_helper(p, X_next(p), X_next(q), p->SML_next());
};

std::array<double, 3> grad_w_Bspline_next(const std::tuple<networkPoint *, Tddd> &p_X, const networkPoint *q) {
   return grad_w_Bspline_next_helper(std::get<0>(p_X), std::get<1>(p_X), X_next(q), std::get<0>(p_X)->SML_next());
};

std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const Tddd &X) {
   return grad_w_Bspline_next_helper(p, X_next(p), X, p->SML_next());
};

double Dot_grad_w_Bspline_next_helper(const networkPoint *p, const Tddd &pX, const Tddd &qX, const double r) {
   // if (p->isFluid && p->isNotSurfaceButNearSurface)
   //    return Dot_grad_w_Bspline(pX, qX, r, p->inv_grad_corr_M_next_mirror);
   // else
   return Dot_grad_w_Bspline(pX, qX, r, p->inv_grad_corr_M_next);
};

double Dot_grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
   return Dot_grad_w_Bspline_next_helper(p, X_next(p), X_next(q), p->SML_next());
};

double Dot_grad_w_Bspline_next(const std::tuple<networkPoint *, Tddd> &p_X, const networkPoint *q) {
   return Dot_grad_w_Bspline_next_helper(std::get<0>(p_X), std::get<1>(p_X), X_next(q), std::get<0>(p_X)->SML_next());
};

double Dot_grad_w_Bspline_next(const std::tuple<networkPoint *, Tddd> &p_X, const Tddd &Y) {
   return Dot_grad_w_Bspline_next_helper(std::get<0>(p_X), std::get<1>(p_X), Y, std::get<0>(p_X)->SML_next());
};

/* -------------------------------------------------------------------------- */
// std::array<double, 3> grad_w_Bspline_next(const std::tuple<networkPoint *, Tddd> &p_X, const Tddd &Y) {
//    auto [p, X] = p_X;
//    if (p->isFluid && p->isNotSurfaceButNearSurface)
//       return grad_w_Bspline(X, Y, p->SML_next(), p->inv_grad_corr_M_next_mirror);
//    else
//       return grad_w_Bspline(X, Y, p->SML_next(), p->inv_grad_corr_M_next);
// };

// std::array<double, 3> grad_w_Bspline_next(const std::tuple<networkPoint *, Tddd> &p_X, const networkPoint *q) {
//    auto [p, X] = p_X;
//    if (p->isFluid && p->isNotSurfaceButNearSurface)
//       return grad_w_Bspline(X, X_next(q), p->SML_next(), p->inv_grad_corr_M_next_mirror);
//    else
//       return grad_w_Bspline(X, X_next(q), p->SML_next(), p->inv_grad_corr_M_next);
// };

// std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
//    if (p->isFluid && p->isNotSurfaceButNearSurface)
//       return grad_w_Bspline(X_next(p), X_next(q), p->SML_next(), p->inv_grad_corr_M_next_mirror);
//    else
//       return grad_w_Bspline(X_next(p), X_next(q), p->SML_next(), p->inv_grad_corr_M_next);
// };

// std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const Tddd &X) {
//    if (p->isFluid && p->isNotSurfaceButNearSurface)
//       return grad_w_Bspline(X_next(p), X, p->SML_next(), p->inv_grad_corr_M_next_mirror);
//    else
//       return grad_w_Bspline(X_next(p), X, p->SML_next(), p->inv_grad_corr_M_next);
// };

// std::array<double, 3> grad_w_Bspline_next_mirror(const networkPoint *p, const Tddd &X) {
//    return grad_w_Bspline(X_next(p), X, p->SML_next(), p->inv_grad_corr_M_next_mirror);
// };

// double Dot_grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
//    if (p->isFluid && p->isNotSurfaceButNearSurface)
//       return Dot_grad_w_Bspline(X_next(p), X_next(q), p->SML_next(), p->inv_grad_corr_M_next_mirror);
//    else
//       return Dot_grad_w_Bspline(X_next(p), X_next(q), p->SML_next(), p->inv_grad_corr_M_next);
//    // return Dot_grad_w_Bspline_Dot(X_next(p), X_next(q), p->SML_next());
// };

// double Dot_grad_w_Bspline_next(const std::tuple<networkPoint *, Tddd> &p_X, const networkPoint *q) {
//    auto [p, X] = p_X;
//    if (p->isFluid && p->isNotSurfaceButNearSurface)
//       return Dot_grad_w_Bspline(X, X_next(q), p->SML_next(), p->inv_grad_corr_M_next_mirror);
//    else
//       return Dot_grad_w_Bspline(X, X_next(q), p->SML_next(), p->inv_grad_corr_M_next);
// };

// double Dot_grad_w_Bspline_next(const std::tuple<networkPoint *, Tddd> &p_X, const Tddd &Y) {
//    auto [p, X] = p_X;
//    if (p->isFluid && p->isNotSurfaceButNearSurface)
//       return Dot_grad_w_Bspline(X, Y, p->SML_next(), p->inv_grad_corr_M_next_mirror);
//    else
//       return Dot_grad_w_Bspline(X, Y, p->SML_next(), p->inv_grad_corr_M_next);
// };

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
      obj->BucketPoints.apply(X_next(p), p->SML(), [&](const auto &q) {
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
      obj->BucketPoints.apply(p->X, p->SML(), [&](const auto &q) {
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

std::tuple<networkPoint *, Tddd> closest(const networkPoint *p, const auto &RigidBodyObject) {
   double closest_d = 1E+20, d;
   Tddd p_to_q, closest_p_to_q;
   networkPoint *closest_q = nullptr;
   for (const auto &[obj, _] : RigidBodyObject) {
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
   for (const auto &[obj, _] : RigidBodyObject) {
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

#define REFLECTION
void updateParticles(const auto &points,
                     const std::unordered_set<Network *> &target_nets,
                     const std::vector<std::tuple<Network *, Network *>> &RigidBodyObject,
                     const double &particle_spacing,
                     const double dt) {
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
         p->RK_X.push(p->U_SPH);  // 位置
         p->setXSingle(p->tmp_X = p->RK_X.getX());
         p->RK_U.push(p->DUDt_SPH);  // 速度
         p->U_SPH = p->RK_U.getX();
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
            for (const auto &[closest_p, v_f2w] : {closest(p, RigidBodyObject) /*, closest_next(p, RigidBodyObject)*/}) {
               if (closest_p != nullptr) {
                  auto n = Normalize(closest_p->interp_normal_original);
                  auto n_d_f2w = Norm(Projection(v_f2w, n));
                  auto ratio = (d0 - n_d_f2w) / d0;
                  if (Norm(v_f2w) < 1. * d0 && Norm(closest_p->X - p->X) < d0) {
                     // auto ratio = (d0 - n_d_f2w) / d0;
                     if (Dot(p->U_SPH, n) < 0) {
                        auto tmp = -0.1 * ratio * Projection(p->U_SPH, n) / p->RK_X.get_dt();
                        p->DUDt_modify_SPH += tmp;
                        p->DUDt_SPH += tmp;
   #if defined(USE_RungeKutta)
                        p->RK_U.repush(p->DUDt_SPH);  // 速度
                        p->U_SPH = p->RK_U.getX();
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