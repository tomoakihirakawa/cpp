#ifndef SPH_weightingFunctions_H
   #define SPH_weightingFunctions_H

   #include "kernelFunctions.hpp"
   #include "vtkWriter.hpp"

   #define REFLECTION

   // #define surface_zero_pressure

   //$ ------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   #define USE_SPP_Fluid
   #if defined(USE_SPP_Fluid)
      #define SPP_U_coef 1.
      #define SPP_p_coef -1.
   #endif
   /* -------------------------------------------------------------------------- */
   #define USE_SPP_Wall
   #if defined(USE_SPP_Wall)
      #define SPP_U_coef_of_Wall 1.
      #define SPP_p_coef_of_Wall -1.
      #define SPP_DUDt_coef_of_Wall 1.
      #define SPP_rho_coef_of_Wall 1.
   #endif
/* -------------------------------------------------------------------------- */

   #define POWER 1.

const double reflection_factor = .5;
const double asobi = 0.;

   #include "SPH_Functions.hpp"
/* -------------------------------------------------------------------------- */

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
      p->p_SPH_SPP = p->total_weight = 0;
      double w = 0, own_w = 0;
      auto func = [&](const auto &q) {
         p->total_weight += (w = q->volume * std::pow(w_Bspline(Norm(p->X - q->X), p->radius_SPH), POWER));

         if (q == p)
            own_w = w;
         else
            p->p_SPH_SPP += q->p_SPH * w;
      };
      net->BucketPoints.apply(p->X, p->radius_SPH, func);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->radius_SPH, func);
      p->p_SPH_SPP /= (1 - own_w);
      // p->p_SPH_SPP /= (p->total_weight - own_w);
      // /* -------------------------------------------------------------------------- */
      // p->p_SPH_SPP = p->p_SPH;
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
      p->COM_SPH.fill(0.);
      p->interpolated_normal_SPH_original.fill(0.);
      p->interpolated_normal_SPH_original_all.fill(0.);
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

void setWallDUDt(auto &net, const std::unordered_set<networkPoint *> &wall_p, const auto &RigidBodyObject) {
   Timer watch;
   const double POW = 1., a = 1.;
   const bool mirroring = true;
   const bool account_volume = true;
   #pragma omp parallel
   for (const auto &PW : wall_p)
   #pragma omp single nowait
   {
      PW->total_weight = 0;
      PW->tmp_ViscousAndGravityForce_ = PW->ViscousAndGravityForce_ = PW->DUDt_SPH_ = PW->gradP_SPH_ = {0., 0., 0.};
      const Tddd markerX = PW->X + 2 * PW->normal_SPH;
      const Tddd accel = {0., 0., 0.};
      const Tddd mirror_wall = PW->X + PW->normal_SPH;
      const Tddd n = Normalize(PW->normal_SPH);
      double W;
      int count = 0;
      /* -------------------------------------------------------------------------- */
      auto func = [&](const networkPoint *PF, const Tddd &nearX, const double coef = 1.) {
         if (PF->getNetwork()->isFluid || PF->isFluid) {
            // if (PF->getNetwork()->isFluid) {
            PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
            PW->DUDt_SPH_ += coef * PF->DUDt_SPH * W;
            PW->ViscousAndGravityForce_ += coef * PF->ViscousAndGravityForce * W;
            PW->tmp_ViscousAndGravityForce_ += coef * PF->tmp_ViscousAndGravityForce * W;
            if (mirroring) {
               auto mirrored_X = Mirror(nearX, mirror_wall, n);
               PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW));
               PW->DUDt_SPH_ += coef * PF->DUDt_SPH * W;
               PW->ViscousAndGravityForce_ += coef * PF->ViscousAndGravityForce * W;
               PW->tmp_ViscousAndGravityForce_ += coef * PF->tmp_ViscousAndGravityForce * W;
            }
         }
      };
      /* -------------------------------------------------------------------------- */
      net->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &PF) {
         func(PF, PF->X);
   #ifdef USE_SPP_Wall
         if (PF->isSurface) {
            // if (canSetSPP(net, RigidBodyObject, PF))
            func(PF, SPP_X(PF), SPP_DUDt_coef_of_Wall);
            func(PF, SPP_X(PF, 2), SPP_DUDt_coef_of_Wall);
         }
   #endif
      });
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &qW) {
            func(qW, qW->X);
         });
      // for (const auto &[obj, poly] : RigidBodyObject)
      //    obj->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &qW) {
      //       func(qW, qW->X);
      //    });
   }

   for (const auto &PW : wall_p) {
      if (isFinite(1. / PW->total_weight, 1E+15) && !PW->isFluid) {
         PW->DUDt_SPH = PW->DUDt_SPH_ / PW->total_weight;
         PW->ViscousAndGravityForce = PW->ViscousAndGravityForce_ / PW->total_weight;
         PW->tmp_ViscousAndGravityForce = PW->tmp_ViscousAndGravityForce_ / PW->total_weight;
      }
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
         const Tddd markerX = PW->X + 2 * PW->normal_SPH;
         const Tddd accel = {0., 0., 0.};
         const Tddd mirror_wall = PW->X + PW->normal_SPH;
         const Tddd n = Normalize(PW->normal_SPH);
         double W;
         int count = 0;
         initialize(PW);
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
            if (PF->getNetwork()->isFluid) {
               PW->total_N += std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW);
               PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->total_weight_ += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->lap_tmpU_ += cU * PF->lap_tmpU * W;
               PW->U_SPH_ += cU * PF->U_SPH * W;
               PW->tmp_U_SPH_ += cU * PF->tmp_U_SPH * W;
               PW->grad_div_U_ += PF->grad_div_U * W;
               PW->p_SPH_ += cp * (spp ? PF->p_SPH : PF->p_SPH) * W;
               PW->div_U_ += cU * PF->div_U * W;
               PW->rho_ += cRho * PF->rho * W;
               PW->volume_ += PF->volume * W;
               PW->DUDt_SPH_ += PF->DUDt_SPH * W;
               PW->gradP_SPH_ += PF->gradP_SPH * W;
               if (mirroring && !applyToWall) {
                  const auto mirrored_X = Mirror(nearX, mirror_wall, PW->normal_SPH);
                  PW->total_N += std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW);
                  PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW));
                  PW->total_weight_ += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(mirrored_X - markerX), PW->radius_SPH * a), POW));
                  PW->p_SPH_ += cp * (spp ? PF->p_SPH : PF->p_SPH) * W;
                  PW->lap_tmpU_ += cU * PF->lap_tmpU * W;
                  PW->grad_div_U_ += PF->grad_div_U * W;
                  PW->U_SPH_ += cU * PF->U_SPH * W;
                  PW->tmp_U_SPH_ += cU * PF->tmp_U_SPH * W;
                  PW->div_U_ += cU * PF->div_U * W;
                  PW->rho_ += cRho * PF->rho * W;
                  PW->volume_ += PF->volume * W;
                  PW->DUDt_SPH_ += PF->DUDt_SPH * W;
                  PW->gradP_SPH_ += PF->gradP_SPH * W;
               }
            } else if (PF->isFluid) {
               // 壁粒子でもこの壁粒子の流速はゼロと確定している．
               // 他の壁粒子は，ミラリングによって補う．
               PW->total_N += std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW);
               // PW->total_weight += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               PW->total_weight_ += (W = (account_volume ? PF->volume : 1.) * std::pow(w_Bspline(Norm(nearX - markerX), PW->radius_SPH * a), POW));
               // PW->lap_tmpU_ += cU * PF->lap_tmpU * W;
               PW->U_SPH_ += PF->U_SPH * W;
               PW->tmp_U_SPH_ += PF->tmp_U_SPH * W;
               // PW->grad_div_U_ += PF->grad_div_U * W;
               // PW->p_SPH_ += cp * (spp ? PF->p_SPH_SPP : PF->p_SPH) * W;
               // PW->div_U_ += cU * PF->div_U * W;
               PW->rho_ += cRho * PF->rho * W;
               PW->volume_ += PF->volume * W;
               // PW->DUDt_SPH_ += 0 * PF->DUDt_SPH * W;
               // PW->gradP_SPH_ += PF->gradP_SPH * W;
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
               func(PF, SPP_X(PF, 2.), true);
            }
   #endif
         });
         for (const auto &[obj, poly] : RigidBodyObject)
            obj->BucketPoints.apply(markerX, PW->radius_SPH, [&](const auto &qW) {
               func(qW, qW->X);
            });
      }
   };
   //-------------------------
   auto set = [&]() {
      for (const auto &PW : wall_p) {
         {
            if (PW->total_weight > 1E-15) {
               PW->p_SPH = PW->p_SPH_ / PW->total_weight;
               PW->lap_tmpU = PW->lap_tmpU_ / PW->total_weight_;
               if (!PW->isFluid) {
                  /*流体として扱う壁粒子は，流速をmapによって更新しない*/
                  PW->U_SPH = PW->U_SPH_ / PW->total_weight_;
                  PW->tmp_U_SPH = PW->tmp_U_SPH_ / PW->total_weight_;
                  PW->grad_div_U = PW->grad_div_U_ / PW->total_weight_;
                  PW->div_U = PW->div_U_ / PW->total_weight;
                  PW->DUDt_SPH = (PW->DUDt_SPH_) / PW->total_weight;
               }
            } else
               initialize(PW);
         }
         std::vector<bool> Bools = {!isFinite(PW->p_SPH), !isFinite(PW->rho), !isFinite(PW->volume), !isFinite(PW->lap_tmpU), !isFinite(PW->U_SPH), !isFinite(PW->tmp_U_SPH)};
         if (std::any_of(Bools.begin(), Bools.end(), [](bool b) { return b; })) {
            std::stringstream ss;
            ss << Bools << std::endl;
            ss << "PW->p_SPH = " << PW->p_SPH << std::endl;
            ss << "PW->isFluid = " << PW->isFluid << std::endl;
            ss << "PW->rho = " << PW->rho << std::endl;
            ss << "PW->volume = " << PW->volume << std::endl;
            ss << "PW->lap_tmpU = " << PW->lap_tmpU << std::endl;
            ss << "PW->U_SPH = " << PW->U_SPH << std::endl;
            ss << "PW->tmp_U_SPH = " << PW->tmp_U_SPH << std::endl;
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
      PW->U_SPH *= -1.;
      PW->tmp_U_SPH *= -1.;

      // PW->U_SPH *= 0.;
      // PW->tmp_U_SPH *= 0.;

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
std::unordered_set<networkPoint *> wall_as_fluid;
std::unordered_set<networkPoint *> wall_p_surface;
std::unordered_set<networkPoint *> wall_p;

void test_Bucket(const auto &water, const auto &nets, const std::string &output_directory, const auto &r) {
   std::filesystem::create_directory(output_directory);

   {
      vtkPolygonWriter<networkPoint *> vtp;
      for (const auto &net : nets)
         for (const auto &p : net->getPoints())
            vtp.add(p);
      std::ofstream ofs(output_directory + "all.vtp");
      vtp.write(ofs);
      ofs.close();
   }

   int i = 0;
   for (const auto &p : water->getPoints()) {
      {
         vtkPolygonWriter<networkPoint *> vtp;
         vtp.add(p);
         std::ofstream ofs(output_directory + "center" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
         ofs.close();
      }
      {
         std::unordered_map<networkPoint *, double> Bspline3;
         std::unordered_map<networkPoint *, double> Bspline5;
         double sum3 = 0, sum5 = 0;
         vtkPolygonWriter<networkPoint *> vtp;
         int j = 0;
         for (const auto &net : nets)
            net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
               vtp.add(q);
               sum3 += Bspline3[q] = w_Bspline3(Norm(q->X - p->X), p->radius_SPH) * p->volume;
               sum5 += Bspline5[q] = w_Bspline5(Norm(q->X - p->X), p->radius_SPH) * p->volume;
               j++;
            });
         vtp.addPointData("w_Bspline3", Bspline3);
         vtp.addPointData("w_Bspline5", Bspline5);

         std::cout << "number of points in cell: " << j << std::endl;
         std::ofstream ofs(output_directory + "inCell" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
         ofs.close();
      }
      {
         std::unordered_map<networkPoint *, double> Bspline3;
         std::unordered_map<networkPoint *, double> Bspline5;
         double sum3 = 0, sum5 = 0;
         vtkPolygonWriter<networkPoint *> vtp;
         int j = 0;
         for (const auto &net : nets)
            net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
               if (Distance(p->X, q->X) < p->radius_SPH) {
                  vtp.add(q);
                  sum3 += Bspline3[q] = w_Bspline3(Norm(q->X - p->X), p->radius_SPH) * p->volume;
                  sum5 += Bspline5[q] = w_Bspline5(Norm(q->X - p->X), p->radius_SPH) * p->volume;
                  j++;
               }
            });
         vtp.addPointData("w_Bspline3", Bspline3);
         vtp.addPointData("w_Bspline5", Bspline5);

         std::cout << "        points: " << j << std::endl;
         std::cout << "sum B-spline 3: " << sum3 << std::endl;
         std::cout << "sum B-spline 5: " << sum5 << std::endl;
         std::ofstream ofs(output_directory + "inSphere" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
         ofs.close();
      }
      i++;
   }

   // check if each cell properly stores objects
   {
      int I = 0, J = 0;
      for (auto i = 0; i < water->BucketPoints.xsize; ++i)
         for (auto j = 0; j < water->BucketPoints.ysize; ++j)
            for (auto k = 0; k < water->BucketPoints.zsize; ++k) {
               {
                  vtkPolygonWriter<networkPoint *> vtp;
                  for (const auto &net : nets)
                     vtp.add(net->BucketPoints.buckets[i][j][k]);
                  std::ofstream ofs(output_directory + "each_cell" + std::to_string(I++) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
               {
                  vtkPolygonWriter<Tddd> vtp;
                  for (const auto &net : nets)
                     vtp.add(net->BucketPoints.itox(i, j, k));
                  std::ofstream ofs(output_directory + "each_cell_position" + std::to_string(J++) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
            }
   }
};

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
      std::unordered_set<Network *> net_RigidBody;
      // b# -------------------------------------------------------------------------- */
      for (const auto &[a, b, _] : RigidBodyObjectIN) {
         RigidBodyObject.push_back({a, b});
         net_RigidBody.emplace(a);
      }
      // b# -------------- バケットの生成, p->radius_SPHの範囲だけ点を取得 --------------- */
      DebugPrint("バケットの生成", Green);
      net->makeBucketPoints(particle_spacing * 0.8);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->makeBucketPoints(particle_spacing * 0.8);
      DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "バケットの生成");

      // test_Bucket(net, Append(net_RigidBody, net), "./test_SPH_Bucket/", particle_spacing);

      // void setDataOmitted(auto &vtp, const auto &Fluid) {
      //    std::unordered_map<networkPoint *, Tddd> normal_SPH;
      //    for (const auto &p : Fluid->getPoints())
      //       normal_SPH[p] = p->normal_SPH;
      //    vtp.addPointData("normal_SPH", normal_SPH);}

      // b# --------------------------- 平滑化距離の計算 ------------------------------ */
      /*     密度, 平滑化距離      */
      DebugPrint(Green, "固定の平滑化距離の計算: C_SML * particle_spacing = ", C_SML, " * ", particle_spacing, " = ", C_SML * particle_spacing);
      for (const auto &p : net->getPoints()) {
         // p->C_SML = C_SML * std::pow(_WATER_DENSITY_ / p->rho, 3);  // C_SMLを密度の応じて変えてみる．
         p->setDensity(_WATER_DENSITY_);
         p->isFreeFalling = false;
         p->isInsideOfBody = false;
         p->p_SPH_SPP = 0;
         p->DUDt_SPH = p->lap_U = {0, 0, 0};
         p->isCaptured = true;
      }
      // b# ------------------------------------------------------- */
      // b#                    関連する壁粒子をマーク                   */
      // b# ------------------------------------------------------- */

      wall_as_fluid.clear();
      wall_p_surface.clear();
      wall_p.clear();

      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            p->setDensityVolume(0, 0);
            p->isFluid = false;
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
               if (Norm(q->normal_SPH) < 1E-12 && Distance(p, q) < p->radius_SPH / p->C_SML * 1.1) {
                  q->isCaptured_ = true;
                  q->isFluid = true;
               }
            });
         }
      };
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints())
            if (p->isCaptured) {
               wall_p.emplace(p);
               if (p->isFluid)
                  wall_as_fluid.emplace(p);
            }

      DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "関連する壁粒子をマークし，保存");
      // b# -------------------------------------------------------------------------- */
      // b# --------------- CFL条件を満たすようにタイムステップ間隔dtを設定 ----------------- */
      // b# -------------------------------------------------------------------------- */
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
      // b# -------------------------------------------------------------------------- */
      // b# -------------------------------------------------------------------------- */
      // b# ----------------------- ルンゲクッタの準備 ------------------- */
      for (const auto &p : net->getPoints()) {
         p->RK_U.initialize(dt, real_time, p->U_SPH, RK_order);
         p->RK_X.initialize(dt, real_time, p->X, RK_order);
         p->RK_P.initialize(dt, real_time, p->p_SPH, RK_order);
      }
      // b# ======================================================= */
      // b#                  ルンゲクッタを使った時間積分                */
      // b# ======================================================= */
      /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/
      do {
         for (const auto &p : net->getPoints())
            p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3));
         dt = (*net->getPoints().begin())->RK_X.getdt();
         std::cout << "dt = " << dt << std::endl;
         mapValueOnWall(net, wall_p, RigidBodyObject);

         // b$ -------------------------------------------------------------------------- */

         Lap_U(net->getPoints(), Append(net_RigidBody, net));
         setLap_U(net->getPoints(), dt);

         Lap_U(wall_as_fluid, Append(net_RigidBody, net));
         setLap_U(wall_as_fluid, dt);

         // b$ -------------------------------------------------------------------------- */

         setWallDUDt(net, wall_p, RigidBodyObject);

         // b$ ---------------------------------- (1) 位置の更新 for (2) 密度の更新--------------------------------- */

   #pragma omp parallel
         for (const auto &p : net->getPoints())
   #pragma omp single nowait
            p->setX(p->tmp_X);

         mapValueOnWall(net, wall_p, RigidBodyObject);

         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "粘性項の∇.∇Uを計算し，次にU*を計算");
         DebugPrint("U*を使って，仮想的な位置X^*へ粒子を移動", Green);

         div_tmpU(net->getPoints(), Append(net_RigidBody, net));
         div_tmpU(wall_as_fluid, Append(net_RigidBody, net));

         set_tmpLap_U(net->getPoints());
         set_tmpLap_U(wall_as_fluid);
         // b$ ---------------------------------- (2) 密度の更新 --------------------------------- */
         for (const auto &p : net->getPoints()) {
            p->setDensity(p->rho_ = p->rho + (p->DrhoDt_SPH = -p->rho * p->div_U) * dt);
         }
         for (const auto &p : wall_as_fluid) {
            p->setDensity(p->rho_ = p->rho + (p->DrhoDt_SPH = -p->rho * p->div_U) * dt);
         }

         mapValueOnWall(net, wall_p, RigidBodyObject);

         /*
         nablaは同じでないといけないので，
         div(U) = div(grad(P^n+1))
         のdivは同じように計算されなければならないだろう:同じ密度，体積を使う．
         */

         // b$ ---------------------- 仮位置における div(U*) and laplacian(U*) ----------------------*/

         div_tmpU(net->getPoints(), Append(net_RigidBody, net));
         div_tmpU(wall_as_fluid, Append(net_RigidBody, net));

         set_tmpLap_U(net->getPoints());
         set_tmpLap_U(wall_as_fluid);

         // setWallDUDt(net, wall_p, RigidBodyObject);  // 壁の粘性項を計算する

         // b$ -------------------------------------------------------------------------- */

         // mapValueOnWall(net, wall_p, RigidBodyObject, {1, 0, 0});  // ここではDUDtなど考慮せずに壁の圧力を蹴っているす．圧力勾配を計算する時のみDUDtを考慮するのがいいかもしれない．
         // b% ------------------------------------------------------ */
         // b%  　　　　　             仮位置における圧力Pの計算                   */
         // b% ------------------------------------------------------ */
   #ifdef surface_zero_pressure
         for (const auto &p : net->getPoints())
            if (p->isSurface)
               p->p_SPH = 0;
   #endif

         mapValueOnWall(net, wall_p, RigidBodyObject);

         for (const auto &p : net->getPoints())
            p->NR_pressure.initialize(p->p_SPH);

         for (const auto &p : wall_as_fluid)
            p->NR_pressure.initialize(p->p_SPH);

         /*
         もし水面の圧力をgrad(p)を計算する直前でゼロとしてしまうと，
         水面は落ちてきてしまい，計算が終わってしまう．
         圧力を計算する直前だけ水面の圧力をゼロとすると，
         水面の圧力はゼロとはならないが，圧力勾配の計算結果で水面が落ち込むことはない．
         */
         //
         DebugPrint("仮位置における圧力Pの計算", Magenta);
         //
   // #define ISPH
   #if defined(ISPH)

         DebugPrint("activate");
         V_d b(net->getPoints().size()), x0(net->getPoints().size());
         size_t i = 0;
         for (const auto &p : net->getPoints()) {
            p->setIndexCSR(i);
            b[i] = p->value = p->div_tmpU;
            x0[i] = p->p_SPH;
            i++;
         }

         DebugPrint("Lap_P");
         Lap_P(net->getPoints(), Append(net_RigidBody, net));
         for (const auto &p : net->getPoints())
            if (p->isSurface) {
               p->column_value.clear();
               p->increment(p, 1.);
               p->value = 0.;
            }
         gmres gm(net->getPoints(), b, x0, 100);
         std::cout << "gm.err : " << gm.err << std::endl;
         for (const auto &p : net->getPoints()) {
            p->p_SPH = gm.x[p->getIndexCSR()];
         }
   #else
         // #define NewtonMethod
         int loopnum = 1;
      #if defined(NewtonMethod)
         for (auto i = 0; i < loopnum; ++i)
      #endif
         {

      #pragma omp parallel
            for (const auto &p : net->getPoints())
      #pragma omp single nowait
               nextPressure(p, Append(net_RigidBody, net), dt);

      #pragma omp parallel
            for (const auto &p : wall_as_fluid)
      #pragma omp single nowait
               nextPressure(p, Append(net_RigidBody, net), dt);

      #if defined(NewtonMethod)
            double sum = 0;
            for (const auto &p : net->getPoints()) {
               if (i == 0)
                  p->p_SPH = p->p_SPH_;
               else
                  p->p_SPH = p->p_SPH_ = p->NR_pressure.X;
               sum += p->div_U_error;
            }
            for (const auto &p : wall_as_fluid) {
               if (i == 0)
                  p->p_SPH = p->p_SPH_;
               else
                  p->p_SPH = p->p_SPH_ = p->NR_pressure.X;
               sum += p->div_U_error;
            }
            if (i == 0)
               std::cout << "EISPH" << Red << "sum/N " << sum / (net->getPoints().size() + wall_as_fluid.size()) << colorOff << std::endl;
            else if (i == loopnum - 1)
               std::cout << " 最後" << Red << "sum/N " << sum / (net->getPoints().size() + wall_as_fluid.size()) << colorOff << std::endl;
      #endif
         }

         for (const auto &p : net->getPoints()) {
            p->dp_SPH = p->p_SPH - p->dp_SPH;
            p->DPDt_SPH = (p->p_SPH - p->RK_P.getX()) / dt;
            // p->p_SPH = p->p_SPH_;
            p->RK_P.push(p->DPDt_SPH);  // 圧力
            // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
            p->p_SPH = p->p_SPH_;
         }

            /* -------------------------------------------------------------------------- */
            // if (real_time > 0.01) {
            //    DebugPrint("activate");
            //    V_d b, x0;
            //    size_t i = 0;
            //    std::unordered_set<networkPoint *> points;
            //    points.reserve(net->getPoints().size() + wall_as_fluid.size());
            //    b.reserve(net->getPoints().size() + wall_as_fluid.size());
            //    x0.reserve(net->getPoints().size() + wall_as_fluid.size());
            //    for (const auto &p : net->getPoints()) {
            //       if (p->isCaptured) {
            //          points.emplace(p);
            //          p->setIndexCSR(i++);
            //          p->exclude(true);
            //          b.emplace_back(p->value = p->rho_ * p->div_tmpU / dt);
            //          x0.emplace_back(p->p_SPH);
            //       }
            //    }
            //    for (const auto &p : wall_as_fluid) {
            //       if (p->isCaptured) {
            //          points.emplace(p);
            //          p->setIndexCSR(i++);
            //          p->exclude(true);
            //          b.emplace_back(p->value = p->rho_ * p->div_tmpU / dt);
            //          x0.emplace_back(p->p_SPH);
            //       }
            //    }

            //    DebugPrint("Lap_P");
            //    Lap_P(points, Append(net_RigidBody, net));

            //    for (const auto &p : net->getPoints())
            //       if (p->isSurface) {
            //          p->column_value.clear();
            //          p->increment(p, 1.);
            //          p->value = 0.;
            //          b[p->getIndexCSR()] = 0.;
            //       }

            //    gmres gm(points, b, x0, 200);
            //    std::cout << "gm.err : " << gm.err << std::endl;

            //    for (const auto &p : points) {
            //       p->p_SPH = gm.x[p->getIndexCSR()];
            //    }
            // }
            /* -------------------------------------------------------------------------- */

   #endif
         //% ---------------------------- 計算した圧力をマップするかどうか ---------------------------- */

         // mapValueOnWall(net, wall_p, RigidBodyObject);

         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "仮位置における圧力Pの計算");

         // b% ------------------------------------------------------ */
         // b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
         // b% ------------------------------------------------------ */

         DebugPrint("圧力勾配∇Pを計算 & DU/Dtの計算", Magenta);

   #pragma omp parallel
         for (const auto &p : net->getPoints())
   #pragma omp single nowait
            gradP(p, Append(net_RigidBody, net));

         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "圧力勾配∇Pを計算 & DU/Dtの計算");

         //@ -------------------------------------------------------- */
         //@                        粒子の時間発展                      */
         //@ -------------------------------------------------------- */

         DebugPrint("粒子の時間発展", Green);
   #pragma omp parallel
         for (const auto &p : net->getPoints())
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
   std::unordered_map<networkPoint *, double> isFluid;
   for (const auto &p : Fluid->getPoints())
      isFluid[p] = p->isFluid;
   vtp.addPointData("isFluid", isFluid);
   std::unordered_map<networkPoint *, double> pressure;
   for (const auto &p : Fluid->getPoints())
      pressure[p] = p->p_SPH;
   vtp.addPointData("pressure", pressure);
   //
   std::unordered_map<networkPoint *, double> density;
   density.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      density[p] = p->rho;
   vtp.addPointData("density", density);
   //
   std::unordered_map<networkPoint *, Tddd> tmp_ViscousAndGravityForce;
   for (const auto &p : Fluid->getPoints())
      tmp_ViscousAndGravityForce[p] = Projection(p->tmp_ViscousAndGravityForce, p->normal_SPH);
   vtp.addPointData("Projectioned　tmp_ViscousAndGravityForce", tmp_ViscousAndGravityForce);
   //
   std::unordered_map<networkPoint *, Tddd> ViscousAndGravityForce;
   for (const auto &p : Fluid->getPoints())
      ViscousAndGravityForce[p] = Projection(p->ViscousAndGravityForce, p->normal_SPH);
   vtp.addPointData("Projectioned　ViscousAndGravityForce", ViscousAndGravityForce);
   //
   std::unordered_map<networkPoint *, Tddd> DUDt;
   for (const auto &p : Fluid->getPoints())
      DUDt[p] = p->DUDt_SPH;
   vtp.addPointData("DUDt", DUDt);
};

void setData(auto &vtp, const auto &Fluid, const Tddd &X = {1E+50, 1E+50, 1E+50}) {

   std::unordered_map<networkPoint *, Tddd> ViscousAndGravityForce;
   for (const auto &p : Fluid->getPoints())
      ViscousAndGravityForce[p] = p->ViscousAndGravityForce;
   vtp.addPointData("ViscousAndGravityForce", ViscousAndGravityForce);

   std::unordered_map<networkPoint *, Tddd> tmp_ViscousAndGravityForce;
   for (const auto &p : Fluid->getPoints())
      tmp_ViscousAndGravityForce[p] = p->tmp_ViscousAndGravityForce;
   vtp.addPointData("tmp_ViscousAndGravityForce", tmp_ViscousAndGravityForce);
   //

   if (isFinite(X)) {
      std::unordered_map<networkPoint *, double> W;
      for (const auto &p : Fluid->getPoints())
         W[p] = w_Bspline(Norm(X - p->X), p->radius_SPH);
      vtp.addPointData("W", W);
   }
   //
   std::unordered_map<networkPoint *, Tddd> U;
   for (const auto &p : Fluid->getPoints())
      U[p] = p->U_SPH;
   vtp.addPointData("U", U);

   std::unordered_map<networkPoint *, double> C_SML;
   for (const auto &p : Fluid->getPoints())
      C_SML[p] = p->C_SML;
   vtp.addPointData("C_SML", C_SML);

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