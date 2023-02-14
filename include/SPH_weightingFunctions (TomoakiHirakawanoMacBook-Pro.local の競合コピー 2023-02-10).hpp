#ifndef SPH_weightingFunctions_H
#define SPH_weightingFunctions_H

#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

#define REFLECTION

#define SPP_U 1
#define SPP_U_coef 1
#define SPP_U_coef_of_Wall 1
#define SPP_p 1
#define SPP_p_coef -1
#define SPP_p_coef_of_Wall -1
#define SPP_rho_coef_of_Wall 0.5

#define surface_zero_pressure 0

#define POWER 1.

const double damping_factor = 0.8;
const double asobi = 0.;

/* -------------------------------------------------------------------------- */
Tddd SPP_X(const networkPoint *p) {
   // return p - (p->COM_SPH - p);
   // return p->X - 2 * p->COM_SPH;
   return p->X + p->interpolated_normal_SPH * p->radius_SPH / p->C_SML;
};

auto canSetSPP(Network *net, const auto &RigidBodyObject, const auto &p) {
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
void setPressureForInitialGuess(Network *net) {

#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      p->p_SPH_ = p->total_weight = 0;
      double w = 0, own_w = 0;
      auto func = [&](const auto &q) {
         p->total_weight += (w = q->volume * std::pow(w_Bspline(Norm(p->X - q->X), p->radius_SPH), POWER));
         p->p_SPH_ += q->p_SPH * w;
      };
      net->BucketPoints.apply(p->X, p->radius_SPH, func);
   }
#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      p->p_SPH = p->p_SPH_ / p->total_weight;
   }
};

/* -------------------------------------------------------------------------- */

void setPressureSPP(Network *net, const auto &RigidBodyObject) {

#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      // if (p->isSurface && surface_zero_pressure)
      //    p->p_SPH = 0;
      // /* -------------------------------------------------------------------------- */
      // p->p_SPH_SPP = p->total_weight = 0;
      // double w = 0, own_w = 0;
      // auto func = [&](const auto &q) {
      //    p->total_weight += (w = q->volume * std::pow(w_Bspline(Norm(p->X - q->X), p->radius_SPH), POWER));

      //    if (q == p)
      //       own_w = w;
      //    else
      //       p->p_SPH_SPP += q->p_SPH * w;
      // };
      // net->BucketPoints.apply(p->X, p->radius_SPH, func);
      // for (const auto &[obj, poly] : RigidBodyObject)
      //    obj->BucketPoints.apply(p->X, p->radius_SPH, func);
      // p->p_SPH_SPP /= (1 - own_w);
      // p->p_SPH_SPP /= (p->total_weight - own_w);
      // /* -------------------------------------------------------------------------- */
      p->p_SPH_SPP = p->p_SPH;
      /* -------------------------------------------------------------------------- */
      // p->p_SPH_SPP = p->total_weight = 0;
      // double w = 0, own_w = 0;
      // auto func = [&](const auto &p_fluid) {
      //    p->total_weight += (w = p_fluid->volume * std::pow(w_Bspline(Norm(p->X - p_fluid->X), p->radius_SPH), POWER));
      //    p->p_SPH_SPP += p_fluid->p_SPH * w;
      // };
      // net->BucketPoints.apply(p->X, p->radius_SPH, func);
      // for (const auto &[obj, poly] : RigidBodyObject)
      //    obj->BucketPoints.apply(p->X, p->radius_SPH, func);
      /* -------------------------------------------------------------------------- */
      // p->p_SPH_SPP = p->total_weight = 0;
      // double w = 0;
      // net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &p_fluid) {
      //    w = p_fluid->volume * w_Bspline(Norm(p->X - p_fluid->X), p->radius_SPH * p->C_SML);
      //    p->p_SPH_SPP += coef * p->p_SPH * w;
      //    if (p_fluid->isSurface) {
      //       p->p_SPH_SPP += coef * p->p_SPH * w;
      //    }
      // });
      // p->p_SPH_SPP /= (w_Bspline(0, p->radius_SPH * p->C_SML) + w_Bspline(q->radius_SPH / q->C_SML, p->radius_SPH * p->C_SML));
   }
};

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
      p->COM_SPH = {0., 0., 0.};
      p->interpolated_normal_SPH_original = {0., 0., 0.};
      p->interpolated_normal_SPH_original_all = {0., 0., 0.};
      double total_vol = 0, w;
      net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
         p->COM_SPH += (q->X - p->X) * w;
         total_vol += w;
         if (Between(Distance(p, q), {1E-10, p->radius_SPH})) {
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
            if (Distance(p, q) < p->radius_SPH / p->C_SML * 1.1)
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

void mapValueOnWall_first(auto &net, const std::unordered_set<networkPoint *> &wall_p,
                          const auto &RigidBodyObject) {
   Timer watch;
   double POW, a = 1.;
   bool applyToWall = false;
   bool mirroring = false;
   bool account_volume = true;
   bool do_spp = false;
   auto calc = [&]() {
#pragma omp parallel
      for (const auto &PW : wall_p)
#pragma omp single nowait
      {
         PW->total_weight = 0;
         PW->DUDt_SPH_ = PW->gradP_SPH_ = {0., 0., 0.};
         Tddd markerX = PW->X + 2 * PW->normal_SPH,
              accel = {0., 0., 0.},
              mirror_wall = PW->X + PW->normal_SPH,
              mirrored_X,
              n = Normalize(PW->normal_SPH);
         double W;
         int count = 0;
         /* -------------------------------------------------------------------------- */
         auto func = [&](const networkPoint *PF, const Tddd &nearX, const bool spp = false) {
            if (PF->getNetwork()->isFluid) {
               PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->DUDt_SPH_ += PF->DUDt_SPH * W;
               PW->gradP_SPH_ += PF->gradP_SPH * W;

            } else if (PF->isCaptured_) {
               PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->gradP_SPH_ += PF->gradP_SPH * W;
            }
         };
         /* -------------------------------------------------------------------------- */
         net->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &PF) {
            func(PF, PF->X);
            if (do_spp && SPP_p && PF->isSurface) {
               if (canSetSPP(net, RigidBodyObject, PF))
                  func(PF, SPP_X(PF), true);
            }
         });
         for (const auto &[obj, poly] : RigidBodyObject)
            obj->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &qW) {
               if (PW->isCaptured_)
                  func(qW, qW->X);
            });
      }
   };

   for (const auto &PW : wall_p)
      if (PW->total_weight > 1E-15) {
         PW->DUDt_SPH = (PW->DUDt_SPH_) / PW->total_weight;
         PW->gradP_SPH = (PW->gradP_SPH_) / PW->total_weight;
      }
};

/* -------------------------------------------------------------------------- */

void mapValueOnWall(auto &net,
                    const std::unordered_set<networkPoint *> &wall_p,
                    const auto &RigidBodyObject,
                    const Tiii &indicater = {1, 0, 1}) {
   setPressureSPP(net, RigidBodyObject);
   auto [first, second, do_add] = indicater;
   Timer watch;
   auto initialize = [](networkPoint *PW) {
      PW->div_U_ = PW->volume_ = PW->rho_ = PW->p_SPH_ = PW->total_weight = PW->total_weight_ = 0;
      PW->U_SPH_ = PW->tmp_U_SPH_ = PW->grad_div_U_ = PW->lap_tmpU_ = {0., 0., 0.};
      PW->DUDt_SPH_ = PW->gradP_SPH_ = {0., 0., 0.};
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
         auto func = [&](const networkPoint *PF, const Tddd &nearX, const bool spp = false) {
            double cU = spp ? SPP_U_coef_of_Wall : 1.;
            double cp = spp ? SPP_p_coef_of_Wall : 1.;
            double cRho = spp ? SPP_rho_coef_of_Wall : 1.;
            if (PF->getNetwork()->isFluid) {
               PW->total_N += std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW);
               PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->total_weight_ += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->lap_tmpU_ += cU * PF->lap_tmpU * W;
               PW->U_SPH_ += cU * PF->U_SPH * W;
               PW->tmp_U_SPH_ += cU * PF->tmp_U_SPH * W;
               PW->grad_div_U_ += PF->grad_div_U * W;
               PW->p_SPH_ += cp * (spp ? PF->p_SPH_SPP : PF->p_SPH) * W;
               PW->div_U_ += cU * PF->div_U * W;
               PW->rho_ += cRho * PF->rho * W;
               PW->volume_ += PF->volume * W;
               // PW->DUDt_SPH_ += PF->DUDt_SPH * W;
               // PW->gradP_SPH_ += PF->gradP_SPH * W;
            } else if (PF->isCaptured_) {
               // 壁粒子でもこの壁粒子の流速はゼロと確定している．
               // 他の壁粒子は，ミラリングによって補う．
               PW->total_weight_ += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->lap_tmpU_ += PF->lap_tmpU * W;
               PW->U_SPH_ += 0. * PF->U_SPH * W;
               PW->tmp_U_SPH_ += 0. * PF->tmp_U_SPH * W;
               PW->grad_div_U_ += PF->grad_div_U * W;
               // PW->rho_ += PF->rho * W;
               // PW->volume_ += PW->volume * W;
            }
            if (PF->getNetwork()->isFluid)
               if (mirroring && !applyToWall) {
                  mirrored_X = Mirror(nearX, mirror_wall, PW->normal_SPH);
                  PW->total_N += std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW);
                  PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW));
                  PW->total_weight_ += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW));
                  PW->p_SPH_ += cp * (spp ? PF->p_SPH_SPP : PF->p_SPH) * W;
                  PW->lap_tmpU_ += cU * PF->lap_tmpU * W;
                  PW->grad_div_U_ += PF->grad_div_U * W;
                  PW->U_SPH_ += cU * PF->U_SPH * W;
                  PW->tmp_U_SPH_ += cU * PF->tmp_U_SPH * W;
                  PW->div_U_ += cU * PF->div_U * W;
                  PW->rho_ += cRho * PF->rho * W;
                  PW->volume_ += PF->volume * W;
                  // PW->DUDt_SPH_ += PF->DUDt_SPH * W;
                  // PW->gradP_SPH_ += PF->gradP_SPH * W;
               }
            count++;
         };
         /* -------------------------------------------------------------------------- */
         auto funcAtPreX = [&](const networkPoint *PF, const Tddd &nearX, const bool spp = false) {
            if (PF->getNetwork()->isFluid) {
               PW->total_weight__ += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->DUDt_SPH_ += PF->DUDt_SPH * W;
               PW->gradP_SPH_ += PF->gradP_SPH * W;
               if (mirroring && !applyToWall) {
                  mirrored_X = Mirror(nearX, mirror_wall, PW->normal_SPH);
                  PW->total_weight__ += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW));
                  PW->DUDt_SPH_ += PF->DUDt_SPH * W;
                  PW->gradP_SPH_ += PF->gradP_SPH * W;
               }
            }

            count++;
         };
         /* -------------------------------------------------------------------------- */
         net->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &PF) {
            func(PF, PF->X);
            funcAtPreX(PF, PF->X);
            if (do_spp && SPP_p && PF->isSurface) {
               if (canSetSPP(net, RigidBodyObject, PF))
                  func(PF, SPP_X(PF), true);
            }
         });
         if (applyToWall) {
            for (const auto &[obj, poly] : RigidBodyObject)
               obj->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &qW) {
                  func(qW, qW->X);
               });
            PW->total_weight = 1;
         } else {
            for (const auto &[obj, poly] : RigidBodyObject)
               obj->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &qW) {
                  if (PW->isCaptured_)
                     func(qW, qW->X);
               });
         }
      }
   };

   auto set = [&]() {
      for (const auto &PW : wall_p) {
         {
            if (PW->total_weight > 1E-15) {
               PW->p_SPH = PW->p_SPH_ / PW->total_weight;
               PW->lap_tmpU = PW->lap_tmpU_ / PW->total_weight_;
               PW->U_SPH = PW->U_SPH_ / PW->total_weight_;
               PW->tmp_U_SPH = PW->tmp_U_SPH_ / PW->total_weight_;
               PW->grad_div_U = PW->grad_div_U_ / PW->total_weight_;
               // 以下二つは，流体と同様に変更する
               // if (PW->isCaptured && !PW->isCaptured_) {
               PW->setDensityVolume(PW->rho_ / PW->total_weight, PW->volume_ / PW->total_weight);
               PW->div_U = PW->div_U_ / PW->total_weight;
               //
               PW->DUDt_SPH = PW->DUDt_SPH_ / PW->total_weight__;
               PW->gradP_SPH = PW->gradP_SPH_ / PW->total_weight__;
               // }
            } else
               initialize(PW);
         }
         std::vector<bool> Bools = {!isFinite(PW->p_SPH), !isFinite(PW->rho), !isFinite(PW->volume), !isFinite(PW->lap_tmpU), !isFinite(PW->U_SPH), !isFinite(PW->tmp_U_SPH)};
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
   //! ------------------------------- second time ------------------------------ */
   if (second) {
      auto set = [&]() {
         for (const auto &PW : wall_p) {
            {
               const double alpha = 0.5;
               if (PW->total_weight > 1E-15) {
                  PW->p_SPH = (1 - alpha) * PW->p_SPH + alpha * PW->p_SPH_ / PW->total_weight;
                  PW->lap_tmpU = (1 - alpha) * PW->lap_tmpU + alpha * PW->lap_tmpU_ / PW->total_weight;
                  PW->U_SPH = (1 - alpha) * PW->U_SPH + alpha * PW->U_SPH_ / PW->total_weight;
                  PW->tmp_U_SPH = (1 - alpha) * PW->tmp_U_SPH + alpha * PW->tmp_U_SPH_ / PW->total_weight;
                  PW->setDensityVolume((1 - alpha) * PW->rho + alpha * PW->rho_ / PW->total_weight, (1 - alpha) * PW->volume + alpha * PW->volume_ / PW->total_weight);
               } else
                  initialize(PW);
            }
            std::vector<bool> Bools = {!isFinite(PW->p_SPH), !isFinite(PW->rho), !isFinite(PW->volume), !isFinite(PW->lap_tmpU), !isFinite(PW->U_SPH), !isFinite(PW->tmp_U_SPH)};
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

      // free-slip
      // PW->U_SPH = Reflect(PW->U_SPH, PW->normal_SPH);
      // PW->tmp_U_SPH = Reflect(PW->tmp_U_SPH, PW->normal_SPH);

      // no-slip
      // PW->U_SPH *= -1.;
      // PW->tmp_U_SPH *= -1.;

      // if (PW->isCaptured_) {
      PW->U_SPH *= 0;
      PW->tmp_U_SPH *= 0;
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
         // PW->p_SPH = PW->p_SPH + _WATER_DENSITY_ * Dot(PW->mu_SPH / PW->rho * PW->lap_tmpU - PW->grad_div_U /*/ dtは入れ込み済み*/ + _GRAVITY3_ - accel, -2 * PW->normal_SPH);
         PW->p_SPH = PW->p_SPH + _WATER_DENSITY_ * Dot(PW->DUDt_SPH - accel, -2 * PW->normal_SPH);
      }
   }
   DebugPrint(_green, "Elapsed time: ", _red, watch(), "s ", Magenta, "壁面粒子の圧力を計算", (do_add ? ", DUDtを考慮した圧力" : ""));
};

/* -------------------------------------------------------------------------- */
std::unordered_set<networkPoint *> wall_p_surface;
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
         p->pre_X = p->X;
      }
      //@ ------------------------------------------------------- */
      //@                    関連する壁粒子をマーク                   */
      //@ ------------------------------------------------------- */
      wall_p.clear();
      wall_p_surface.clear();
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            p->setDensityVolume(0, 0);
            p->isFreeFalling = false;
            p->isCaptured = p->isCaptured_ = false;
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
               if (Norm(q->normal_SPH) < 1E-12 && Distance(p, q) < p->radius_SPH / p->C_SML) {
                  q->isCaptured_ = true;
                  std::cout << "What" << std::endl;
                  std::cout << "p->normal_SPH = " << q->normal_SPH << std::endl;
                  std::cout << "Distance(p, q) = " << Distance(p, q) << std::endl;
                  std::cin.ignore();
               }
            });
         }
      };
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints())
            if (p->isCaptured) {
               wall_p.emplace(p);
               if (p->isCaptured_) {
                  wall_p_surface.emplace(p);
               }
            }
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
      auto C_CFL_velocity = 0.05;  // dt = C_CFL_velocity*h/Max(U)
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
            p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3));
         //
         // b@ ======================================================= */
         // b@                  ルンゲクッタを使った時間積分                  */
         // b@ ======================================================= */
         dt = (*net->getPoints().begin())->RK_X.getdt();
         std::cout << "dt = " << dt << std::endl;
         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, " 流体粒子の法線方向の計算，水面の判定");
         // b! ========================================================================== */
         // b!                          U, DUDt, DPDtを計算                                */
         // b! ========================================================================== */
         {
            /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/
            // # ------------------------------------------------------ */
            // #                    ∇.∇UとU*を計算                       */
            // # ------------------------------------------------------ */
            mapValueOnWall(net, wall_p, RigidBodyObject);
            mapValueOnWall_first(net, wall_p, RigidBodyObject);

#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               p->lap_U = p->lap_U_ = {0, 0, 0};
               auto FUNC = [&](const auto &q) {
                  if (q->isCaptured) {
                     //@ ------------------------------------------ */
                     auto func = [&](const Tddd &qX) {
                        auto rij = p->X - qX;
                        auto Uij = p->U_SPH - q->U_SPH;
                        auto nu_nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
                        p->lap_U_ += 1 / (p->mu_SPH / p->rho) * q->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(p->X, qX, p->radius_SPH) /
                                     ((q->rho + p->rho) * Dot(rij, rij) + std::pow(1E-4 * p->radius_SPH / p->C_SML, 2));
                     };
                     //@ ------------------------------------------ */
                     //@ ------------------------------------------ */
                     // auto func = [&](const Tddd &qX) {
                     //    auto rij = p->X - qX;
                     //    auto Uij = p->U_SPH - q->U_SPH;
                     //    p->lap_U_ += 2 * q->mass / p->rho * Dot(rij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(rij, rij) * Uij;
                     // };
                     //@ ------------------------------------------ */
                     if (p != q)
                        func(q->X);
                     if (SPP_U && q->isSurface) {
                        if (canSetSPP(net, RigidBodyObject, q))
                           func(SPP_X(q));
                     }
                  }
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, FUNC);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, FUNC);

               p->lap_U = p->lap_U_;
               p->DUDt_SPH = p->DUDt_SPH_ = p->lap_U_ * (p->mu_SPH / p->rho) + _GRAVITY3_;  // 後で修正されるDUDt
               p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
               p->tmp_X = (p->pre_X = p->X) + p->tmp_U_SPH * dt;
            }

#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               p->setX(p->tmp_X);
            }

            /* -------------------------------------------------------------------------- */
            DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "粘性項の∇.∇Uを計算し，次にU*を計算");
            /* --------------------------------------------------------- */

            mapValueOnWall(net, wall_p, RigidBodyObject);
            DebugPrint("U*を使って，仮想的な位置X^*へ粒子を移動", Green);
            // ! ------------------------------------------------------ */
            // !               div(U^*), DρDt=-ρdiv(U)の計算              */
            // ! ------------------------------------------------------ */

            auto loop = [&](const auto &p) {
               p->div_U = 0;  // this is div(U^*) not div(U^n)
               double own_w = 0;
               auto FUNC = [&](const auto &q) {
                  if (q->isCaptured) {
                     auto Uij = q->tmp_U_SPH - p->tmp_U_SPH;
                     //@ ------------------------------------------ */
                     auto func = [&](const auto &qX) {
                        if (Distance(p, qX) > 1E-12)
                           p->div_U += q->mass / p->rho * Dot(Uij, grad_w_Bspline(p->X, qX, p->radius_SPH));  // 後藤p.25 (2.89)
                     };
                     //@ ------------------------------------------ */
                     func(q->X);
                     if (SPP_U && q->isSurface) {
                        if (canSetSPP(net, RigidBodyObject, q))
                           func(SPP_X(q));
                     }
                  }
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, FUNC);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, FUNC);

               if (!isFinite(p->div_U))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "div_U is not a finite");
            };

#pragma omp parallel
            for (const auto &p : Join(net->getPoints(), wall_p_surface))
#pragma omp single nowait
               loop(p);

            for (const auto &p : Join(net->getPoints(), wall_p_surface)) {
               p->DrhoDt_SPH = -p->rho * p->div_U;
               p->setDensity(p->rho_ = p->rho + p->DrhoDt_SPH * dt);
            }

            /*
            nablaは同じでないといけないので，
            div(U) = div(grad(P^n+1))
            のdivは同じように計算されなければならないだろう:同じ密度，体積を使う．
            */

#pragma omp parallel
            for (const auto &p : Join(net->getPoints(), wall_p_surface))
#pragma omp single nowait
               loop(p);

            mapValueOnWall(net, wall_p, RigidBodyObject);
            DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "div(U^*)を計算");

            // # ------------------------------------------------------ */
            // #                        ∇.∇U*                           */
            // # ------------------------------------------------------ */

#pragma omp parallel
            for (const auto &p : Join(net->getPoints(), wall_p_surface))
#pragma omp single nowait
            {
               p->grad_div_U = p->lap_tmpU = {0, 0, 0};
               auto FUNC = [&](const auto &q) {
                  if (q->isCaptured) {
                     //@ ------------------------------------------ */
                     auto func = [&](const Tddd &qX) {
                        auto rij = p->X - qX;
                        auto Uij = p->tmp_U_SPH - q->tmp_U_SPH;
                        auto nu_nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
                        p->lap_tmpU += 1 / (p->mu_SPH / p->rho) * q->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(p->X, qX, p->radius_SPH) /
                                       ((q->rho + p->rho) * Dot(rij, rij) + std::pow(1E-4 * p->radius_SPH / p->C_SML, 2));
                        //
                        p->grad_div_U += q->div_U * q->volume * grad_w_Bspline(p->X, qX, p->radius_SPH) / dt;
                     };
                     //@ ------------------------------------------ */
                     //@ ------------------------------------------ */
                     // auto func = [&](const Tddd &qX) {
                     //    auto rij = p->X - qX;
                     //    auto Uij = p->tmp_U_SPH - q->tmp_U_SPH;
                     //    p->lap_tmpU += 2 * q->mass / p->rho * Dot(rij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(rij, rij) * Uij;
                     // };
                     //@ ------------------------------------------ */
                     if (p != q)
                        func(q->X);
                     if (SPP_U && q->isSurface) {
                        if (canSetSPP(net, RigidBodyObject, q))
                           func(SPP_X(q));
                     }
                  }
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, FUNC);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, FUNC);
            }

            // ! ------------------------------------------------------ */
            // !  　　　　　        仮位置における圧力Pの計算                 */
            // ! ------------------------------------------------------ */

            // setPressureForInitialGuess(net);
            // setPressureSPP(net, RigidBodyObject);

            mapValueOnWall(net, wall_p, RigidBodyObject);

            std::vector<std::tuple<networkPoint *, NewtonRaphson<double> *, double, double>> P2Newton;

            for (const auto &p : net->getPoints()) {
               if (real_time < 0.0001) {
                  p->dp_SPH = 0;
                  p->dp_SPH_ = 0;
               }
               double a = 0.8;
               p->dp_SPH_ = (1 - a) * p->dp_SPH_ + a * p->dp_SPH;
               P2Newton.push_back({p, new NewtonRaphson(p->p_SPH), 0., 2.});
               p->dp_SPH = p->p_SPH;
            }
            // for (const auto &[obj, poly] : RigidBodyObject)
            //    for (const auto &p : obj->getPoints())
            //       if (p->isCaptured_)
            //          P2Newton.push_back({p, new NewtonRaphson(p->p_SPH), 0., 2.});

// #define Morikawa2019
#define NewtonMethod
            DebugPrint("仮位置における圧力Pの計算", Magenta);
            int loopnum = 5;
#if defined(NewtonMethod)
            for (auto i = 0; i < loopnum; ++i)
#endif
            {

#pragma omp parallel
               for (auto &[p, newton, v, l] : P2Newton)
#pragma omp single nowait
               {
                  double sum_Aij = 0, sum_Aij_Pj = 0, sum_Aij_Pij = 0, sum_Pij = 0, pressure_SPH = 0;
                  p->p_SPH_ = 0;
                  p->gradP_SPH = p->gradP_SPH_ = {0., 0., 0.};
                  //% -------------------------------------------------------------------------- */
                  auto FUNC = [&](const auto &q) {
                     if (q->isCaptured) {
                        //@ ------------------------------------------ */
                        auto func = [&](const auto &qX, const bool spp = false) {
                           auto Xij = p->X - qX;
                           auto qP = (spp ? SPP_p_coef * q->p_SPH_SPP : q->p_SPH);
                           if (Distance(p, qX) > 1E-12) {
                              auto share = Dot(Xij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(Xij, Xij);
#if defined(Morikawa2019)
                              // Morikawa, D., Senadheera, H., & Asai, M. (2021). Explicit incompressible smoothed particle hydrodynamics in a multi-GPU environment for large-scale simulations. Computational Particle Mechanics, 8(3), 493–510. https://doi.org/10.1007/s40571-020-00347-0
                              auto Aij = 2. * q->mass / p->rho * share;
#else
                              // Nomeritae
                              // Shao and Lo
                              // auto Aij = q->mass * 8. / std::pow(q->rho + p->rho, 2) * share;, which is the same as
                              auto Aij = 2. * q->mass / std::pow((q->rho + p->rho) / 2., 2) * share;
#endif

                              sum_Aij += Aij;
                              sum_Pij += p->p_SPH - qP;
                              sum_Aij_Pj += Aij * qP;
                              sum_Aij_Pij += Aij * (p->p_SPH - qP);

                              // WALL
                              p->gradP_SPH += p->rho * (qP / (q->rho * q->rho) + p->p_SPH / (p->rho * p->rho)) * q->mass * grad_w_Bspline(p->X, qX, p->radius_SPH);
                              p->gradP_SPH_ += p->rho * (1. / (p->rho * p->rho)) * q->mass * grad_w_Bspline(p->X, qX, p->radius_SPH);
                           }
                        };
                        //
                        //@ ------------------------------------------ */
                        func(q->X);
                        if (SPP_p && q->isSurface) {
                           if (canSetSPP(net, RigidBodyObject, q))
                              func(SPP_X(q), true);
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
                  p->p_SPH_ = (b + sum_Aij_Pj) / sum_Aij;

#if defined(NewtonMethod)
                  p->div_U_error = std::abs(sum_Aij_Pij - b);
                  // if (!p->getNetwork()->isFluid) {
                  //    auto n = Normalize(p->normal_SPH);
                  //    auto F = Dot(p->U_SPH + dt * (p->DUDt_SPH - p->gradP_SPH / _WATER_DENSITY_), n);
                  //    auto dFdPi = Dot(-dt * p->gradP_SPH_ / _WATER_DENSITY_, n);
                  //    auto FF = F * F / 2.;
                  //    auto dFFdPi = F * dFdPi;
                  //    if (FF > 1E-5 && dFFdPi > 1E-5)
                  //       newton->update(FF, dFFdPi, l);
                  // } else
                  {
                     auto F = v = sum_Aij_Pij - b;
                     auto dFdPi = sum_Aij;
                     auto FF = F * F / 2.;
                     auto dFFdPi = F * dFdPi;
                     newton->update(FF, dFFdPi, l);
                  }
#endif
                  if (!isFinite(p->p_SPH_) || !isFinite(v) || !isFinite(sum_Aij_Pij)) {
                     std::cout << "p->p_SPH_ = " << p->p_SPH_ << std::endl;
                     std::cout << "p->div_U = " << p->div_U << std::endl;
                     std::cout << "b = " << b << std::endl;
                     std::cout << "sum_Aij_Pij = " << sum_Aij_Pj << std::endl;
                     std::cout << "sum_Aij_Pj = " << sum_Aij_Pj << std::endl;
                     std::cout << "sum_Aij = " << sum_Aij << std::endl;
                     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");
                  }
               }
#if defined(NewtonMethod)
               double sum = 0;
               for (const auto &[p, newton, v, _] : P2Newton) {

                  if (i == 0)
                     p->p_SPH = p->p_SPH_;
                  else
                     p->p_SPH = p->p_SPH_ = newton->X;

                  // if (p->isSurface && surface_zero_pressure)
                  //    p->p_SPH = p->p_SPH_ = 0;

                  // if (p->isSurface && SPP_p){
                  //    p->p_SPH_SPP = newtonSPP->X;
                  // }
                  sum += std::abs(v);
               }
               if (i == 0)
                  std::cout << "EISPH" << Red << "sum/N " << sum / (P2Newton.size()) << colorOff << std::endl;
               else if (i == loopnum - 1)
                  std::cout << " 最後" << Red << "sum/N " << sum / (P2Newton.size()) << colorOff << std::endl;
#endif
            }

            for (auto &[p, newton, _, __] : P2Newton) {
               delete newton;
            }

            for (const auto &p : net->getPoints()) {
               p->dp_SPH = p->p_SPH - p->dp_SPH;
               p->DPDt_SPH = (p->p_SPH - p->RK_P.getX()) / dt;
               // p->p_SPH = p->p_SPH_;
               p->RK_P.push(p->DPDt_SPH);  // 圧力
               // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
               p->p_SPH = p->p_SPH_;
            }

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
               p->contact_points_all_SPH = p->contact_points_fluid_SPH = 0;
               p->gradP_SPH = {0., 0., 0.};
               auto FUNC = [&](const auto &q) {
                  if (q->isCaptured) {
                     //@ ------------------------------------------ */
                     auto func = [&](const auto &qX, const bool spp = false) {
                        if (Distance(p, qX) > 1E-12) {
                           // double p_rho2 = coef * q->p_SPH / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                           // auto qP = (spp ? SPP_p_coef * q->p_SPH_SPP : q->p_SPH);
                           // p->gradP_SPH += (qP + p->p_SPH) * q->mass / _WATER_DENSITY_ * grad_w_Bspline(p->X, qX, p->radius_SPH);

                           auto qP = (spp ? SPP_p_coef * q->p_SPH_SPP : q->p_SPH);
                           p->gradP_SPH += p->rho * (qP / (q->rho * q->rho) + p->p_SPH / (p->rho * p->rho)) * q->mass * grad_w_Bspline(p->X, qX, p->radius_SPH);
                        }
                     };
                     //@ ------------------------------------------ */
                     func(q->X);
                     if (SPP_p && q->isSurface) {
                        if (canSetSPP(net, RigidBodyObject, q))
                           func(SPP_X(q), true);
                     }
                  }
                  p->contact_points_all_SPH++;
                  if (q->getNetwork()->isFluid)
                     p->contact_points_fluid_SPH++;
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, FUNC);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, FUNC);

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
            bool isReflect = true;
            while (isReflect && count++ < 30) {
               // const auto X = p->RK_X.getX(p->U_SPH);
               isReflect = false;
               networkPoint *p_wall;
               if (p_wall = closest())
                  if (particle_spacing * (1. - asobi) > Distance(p_wall->X, p->X)) {
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
   std::unordered_map<networkPoint *, double> isWallSurface;
   for (const auto &p : Fluid->getPoints())
      isWallSurface[p] = p->isCaptured_;
   vtp.addPointData("isWallSurface", isWallSurface);
   std::unordered_map<networkPoint *, double> pressure;
   for (const auto &p : Fluid->getPoints())
      pressure[p] = p->p_SPH;
   vtp.addPointData("pressure", pressure);
   std::unordered_map<networkPoint *, double> density;
   density.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      density[p] = p->rho;
   vtp.addPointData("density", density);
};

void setData(auto &vtp, const auto &Fluid) {
   std::unordered_map<networkPoint *, Tddd> U;
   for (const auto &p : Fluid->getPoints())
      U[p] = p->U_SPH;
   vtp.addPointData("U", U);
   //
   std::unordered_map<networkPoint *, double> div_U_error;
   div_U_error.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      div_U_error[p] = p->div_U_error;
   vtp.addPointData("div U error", div_U_error);
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