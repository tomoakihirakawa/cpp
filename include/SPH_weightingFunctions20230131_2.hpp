#ifndef SPH_weightingFunctions_H
#define SPH_weightingFunctions_H

#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

#define REFLECTION

#define SPP_U 1
#define SPP_pressure_coef -1
#define SPP_pressure 1

#define surface_zero_pressure 0

#define POWER 1.
#define _CSML_ 1.
// #define normalize_threshold 10.

auto canSetSPP(Network *net, const auto &RigidBodyObject, const auto &p, const Tddd &X, const double &range) {
   auto func = [&](const auto &q) { return p != q && Distance(q, X) < 1.2 * range; };
   if (net->BucketPoints.any_of(X, 1.5 * range, func))
      return false;
   for (const auto &[obj, poly] : RigidBodyObject)
      if (obj->BucketPoints.any_of(X, 1.5 * range, func))
         return false;
   return true;
};
// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */

void setWallDensityVolume(auto &net, const std::unordered_set<networkPoint *> &wall_p, const auto &RigidBodyObject) {
   Timer watch;
   for (const auto &[obj, poly] : RigidBodyObject)
      for (const auto &p : obj->getPoints())
         p->setDensityVolume(_WATER_DENSITY_, 0.);

#pragma omp parallel
   for (const auto &p_wall : wall_p)
#pragma omp single nowait
   {
      Tddd Xi = p_wall->X + 2 * p_wall->normal_SPH, mirror_wall = p_wall->X + p_wall->normal_SPH, mirrored_X, n = Normalize(p_wall->normal_SPH);
      double weight;
      int count = 0;
      p_wall->rho_ = p_wall->volume_ = p_wall->total_weight = 0;

      if (net->BucketPoints.any_of(Xi, p_wall->radius_SPH / p_wall->C_SML * 2, [&](const auto &p_fluid) { return Distance(Xi, p_fluid) < 2. * p_wall->radius_SPH / p_wall->C_SML; })) {
         auto func = [&](const networkPoint *p_fluid, const Tddd &whereToSet) {
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
            p_wall->rho_ += p_fluid->rho * weight;
            p_wall->volume_ += p_fluid->volume * weight;
            //
            mirrored_X = Mirror(whereToSet, mirror_wall, p_wall->normal_SPH);
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
            p_wall->rho_ += p_fluid->rho * weight;
            p_wall->volume_ += p_fluid->volume * weight;
            count++;
         };

         net->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &p_fluid) {
            func(p_fluid, p_fluid->X);
            // if (SPP_pressure && p_fluid->isSurface) {
            //    auto ps = p_fluid->radius_SPH / p_fluid->C_SML;
            //    auto whereToSet = (p_fluid->X + p_fluid->interpolated_normal_SPH * ps);
            //    if (canSetSPP(net, RigidBodyObject, p_fluid, whereToSet, ps))
            //       func(p_fluid, whereToSet);
            // }
         });
         // for (const auto &[obj, poly] : RigidBodyObject) {
         //    obj->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &q) {
         //       auto vol = std::pow(q->radius_SPH / q->C_SML, 3);
         //       p_wall->total_weight += (weight = vol * std::pow(w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH * _CSML_), POWER));
         //       p_wall->rho_ += _WATER_DENSITY_ * weight;
         //       p_wall->volume_ += vol * weight;
         //    });
         // }
      }
   }

   for (const auto &p_wall : wall_p) {
      {
         if (p_wall->total_weight < 1E-15)
            p_wall->setDensityVolume(_WATER_DENSITY_, 1E-15);
         else
            p_wall->setDensityVolume(p_wall->rho_ / p_wall->total_weight, p_wall->volume_ / p_wall->total_weight);
      }
      if (!isFinite(p_wall->rho_) || !isFinite(p_wall->volume_))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");
   }
};
/* -------------------------------------------------------------------------- */
void setWallDUDt(auto &net, const std::unordered_set<networkPoint *> &wall_p,
                 const auto &RigidBodyObject,
                 const int do_add = 0) {
   Timer watch;
   for (const auto &p : wall_p)
      p->p_SPH_ = p->p_SPH = p->total_weight = 0;
#pragma omp parallel
   for (const auto &p_wall : wall_p)
#pragma omp single nowait
   {
      p_wall->total_weight = 0;
      p_wall->DUDt_SPH_ = {0., 0., 0.};
      Tddd Xi = p_wall->X + 2 * p_wall->normal_SPH, accel, mirror_wall = p_wall->X + p_wall->normal_SPH, mirrored_X, n = Normalize(p_wall->normal_SPH);
      double weight, nu, add, value;
      int count = 0;

      auto func = [&](const networkPoint *p_fluid, const Tddd &whereToSet, const bool spp = false) {
         accel = {0., 0., 0.};
         mirrored_X = Mirror(whereToSet, mirror_wall, p_wall->normal_SPH);
         if (spp) {
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
            // p_wall->p_SPH_ += SPP_pressure_coef * p_fluid->p_SPH_SPP * weight;
            p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
            // p_wall->p_SPH_ += SPP_pressure_coef * p_fluid->p_SPH_SPP * weight;
            p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
         } else {
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
            // p_wall->p_SPH_ += p_fluid->p_SPH * weight;
            p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
            // p_wall->p_SPH_ += p_fluid->p_SPH * weight;
            p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
            // 傾きをつけるためには片方はaddなし
         }
         //
         count++;
      };

      net->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &p_fluid) {
         func(p_fluid, p_fluid->X);
         if (SPP_pressure && p_fluid->isSurface) {
            auto ps = p_fluid->radius_SPH / p_fluid->C_SML;
            auto whereToSet = (p_fluid->X + p_fluid->interpolated_normal_SPH * ps);
            if (canSetSPP(net, RigidBodyObject, p_fluid, whereToSet, ps))
               func(p_fluid, whereToSet, true);
         }
      });
   }
   //
   for (const auto &p_wall : wall_p) {
      {
         if (p_wall->total_weight < 1E-20)
            p_wall->DUDt_SPH_ = {0., 0., 0.};
         else {

            p_wall->DUDt_SPH_ /= p_wall->total_weight;
         }
      }
      if (!isFinite(p_wall->DUDt_SPH_))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");
   }
};

/* -------------------------------------------------------------------------- */

void setWallPressure(auto &net, const std::unordered_set<networkPoint *> &wall_p,
                     const auto &RigidBodyObject,
                     const int do_add = 0) {

   setPressureSPP(net);

   Timer watch;
   for (const auto &p : wall_p)
      p->p_SPH_ = p->p_SPH = p->total_weight = 0;
#pragma omp parallel
   for (const auto &p_wall : wall_p)
#pragma omp single nowait
   {
      p_wall->rho_ = p_wall->p_SPH_ = p_wall->total_weight = 0;
      // p_wall->DUDt_SPH_ = {0., 0., 0.};
      Tddd Xi = p_wall->X + 2 * p_wall->normal_SPH, accel, mirror_wall = p_wall->X + p_wall->normal_SPH, mirrored_X, n = Normalize(p_wall->normal_SPH);
      double weight, nu, add, value;
      int count = 0;

      auto func = [&](const networkPoint *p_fluid, const Tddd &whereToSet, const bool spp = false) {
         accel = {0., 0., 0.};
         mirrored_X = Mirror(whereToSet, mirror_wall, p_wall->normal_SPH);
         if (spp) {
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
            p_wall->p_SPH_ += SPP_pressure_coef * p_fluid->p_SPH_SPP * weight;
            // p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
            p_wall->rho_ += p_fluid->rho * weight;

            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
            p_wall->p_SPH_ += SPP_pressure_coef * p_fluid->p_SPH_SPP * weight;
            // p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
            p_wall->rho_ += p_fluid->rho * weight;
         } else {
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
            p_wall->p_SPH_ += p_fluid->p_SPH * weight;
            // p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
            p_wall->rho_ += p_fluid->rho * weight;
            p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
            p_wall->p_SPH_ += p_fluid->p_SPH * weight;
            // p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
            p_wall->rho_ += p_fluid->rho * weight;
            // 傾きをつけるためには片方はaddなし
         }
         //
         count++;
      };

      net->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &p_fluid) {
         func(p_fluid, p_fluid->X);
         if (SPP_pressure && p_fluid->isSurface) {
            auto ps = p_fluid->radius_SPH / p_fluid->C_SML;
            auto whereToSet = (p_fluid->X + p_fluid->interpolated_normal_SPH * ps);
            if (canSetSPP(net, RigidBodyObject, p_fluid, whereToSet, ps))
               func(p_fluid, whereToSet, true);
         }
      });
   }
   //
   for (const auto &p_wall : wall_p) {
      {
         if (p_wall->total_weight < 1E-20) {
            p_wall->DUDt_SPH_ = {0., 0., 0.};
            p_wall->p_SPH = p_wall->p_SPH_ = 0.;
         } else {
            // p_wall->DUDt_SPH_ /= p_wall->total_weight;

            {
               p_wall->p_SPH_ /= p_wall->total_weight;
               p_wall->rho_ /= p_wall->total_weight;
            }
            Tddd accel = {0., 0., 0.};
            auto mirrored_X = Mirror(p_wall->X, p_wall->X + p_wall->normal_SPH, p_wall->normal_SPH);
            if (do_add)
               p_wall->p_SPH = p_wall->p_SPH_ + p_wall->rho_ * Dot(p_wall->DUDt_SPH_ - accel, p_wall->X - mirrored_X);
            else
               p_wall->p_SPH = p_wall->p_SPH_;
         }
      }
      if (!isFinite(p_wall->p_SPH))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");
   }
};

// void setWallPressure(auto &net, const std::unordered_set<networkPoint *> &wall_p,
//                      const auto &RigidBodyObject,
//                      const int do_add = 0) {

//    setPressureSPP(net);

//    Timer watch;
//    for (const auto &p : wall_p)
//       p->p_SPH_ = p->p_SPH = p->total_weight = 0;
// #pragma omp parallel
//    for (const auto &p_wall : wall_p)
// #pragma omp single nowait
//    {
//       p_wall->p_SPH_ = p_wall->total_weight = p_wall->rho_ = 0;
//       p_wall->DUDt_SPH_ = {0., 0., 0.};
//       Tddd Xi = p_wall->X + 2 * p_wall->normal_SPH, accel, mirror_wall = p_wall->X + p_wall->normal_SPH, mirrored_X, n = Normalize(p_wall->normal_SPH);
//       double weight, nu, add, value;
//       int count = 0;

//       auto func = [&](const networkPoint *p_fluid, const Tddd &whereToSet, const bool spp = false) {
//          accel = {0., 0., 0.};
//          mirrored_X = Mirror(whereToSet, mirror_wall, p_wall->normal_SPH);
//          // add = do_add ? (p_fluid->rho * Dot(p_fluid->DUDt_SPH_ - accel, -2 * p_wall->normal_SPH)) : 0;

//          if (do_add > 0)
//             add = do_add ? ((do_add == 2 ? p_fluid->rho : p_fluid->rho) * Dot((do_add == 2 ? p_fluid->DUDt_SPH : p_fluid->DUDt_SPH_) - accel, mirrored_X - whereToSet)) : 0;
//          if (do_add == 0)
//             add = 0;
//          // if (spp) {
//          //    p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
//          //    p_wall->p_SPH_ += SPP_pressure_coef * p_fluid->p_SPH_SPP * weight;
//          //    p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
//          //    p_wall->rho_ += p_fluid->rho * weight;
//          //    p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
//          //    p_wall->p_SPH_ += SPP_pressure_coef * p_fluid->p_SPH_SPP * weight;
//          //    p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
//          //    p_wall->rho_ += p_fluid->rho * weight;
//          // } else {
//          //    p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
//          //    p_wall->p_SPH_ += p_fluid->p_SPH * weight;
//          //    p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
//          //    p_wall->rho_ += p_fluid->rho * weight;
//          //    p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
//          //    p_wall->p_SPH_ += p_fluid->p_SPH * weight;
//          //    p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
//          //    p_wall->rho_ += p_fluid->rho * weight;
//          //    // 傾きをつけるためには片方はaddなし
//          // }
//          /* -------------------------------------------------------------------------- */
//          if (spp) {
//             p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
//             p_wall->p_SPH_ += SPP_pressure_coef * (p_fluid->p_SPH_SPP + add) * weight;
//             p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
//             p_wall->rho_ += p_fluid->rho * weight;
//             p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
//             p_wall->p_SPH_ += SPP_pressure_coef * p_fluid->p_SPH_SPP * weight;
//             p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
//             p_wall->rho_ += p_fluid->rho * weight;
//          } else {
//             p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
//             p_wall->p_SPH_ += (p_fluid->p_SPH + add) * weight;
//             p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
//             p_wall->rho_ += p_fluid->rho * weight;
//             p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
//             p_wall->p_SPH_ += p_fluid->p_SPH * weight;
//             p_wall->DUDt_SPH_ += p_fluid->DUDt_SPH_ * weight;
//             p_wall->rho_ += p_fluid->rho * weight;
//             // 傾きをつけるためには片方はaddなし
//          }
//          //
//          count++;
//       };

//       net->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &p_fluid) {
//          func(p_fluid, p_fluid->X);
//          if (SPP_pressure && p_fluid->isSurface) {
//             auto ps = p_fluid->radius_SPH / p_fluid->C_SML;
//             auto whereToSet = (p_fluid->X + p_fluid->interpolated_normal_SPH * ps);
//             if (canSetSPP(net, RigidBodyObject, p_fluid, whereToSet, ps))
//                func(p_fluid, whereToSet, true);
//          }
//       });
//    }
//    //
//    for (const auto &p_wall : wall_p) {
//       {
//          if (p_wall->total_weight < 1E-15) {
//             p_wall->DUDt_SPH_ = {0., 0., 0.};
//             p_wall->p_SPH = p_wall->p_SPH_ = 0.;
//          } else {
//             p_wall->DUDt_SPH_ /= p_wall->total_weight;
//             p_wall->p_SPH_ /= p_wall->total_weight;
//             p_wall->rho_ /= p_wall->total_weight;
//             Tddd accel = {0., 0., 0.};
//             // if (do_add == 1)
//             //    p_wall->p_SPH = p_wall->p_SPH_ + _WATER_DENSITY_ * Dot(p_wall->DUDt_SPH_ - accel, -2 * p_wall->normal_SPH);
//             // else if (do_add == 2)
//             //    p_wall->p_SPH = p_wall->p_SPH_ + _WATER_DENSITY_ * Dot(_GRAVITY3_ - accel, -2 * p_wall->normal_SPH);
//             // else
//             p_wall->p_SPH = p_wall->p_SPH_;
//          }
//       }
//       if (!isFinite(p_wall->p_SPH))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");
//    }
// };
/* -------------------------------------------------------------------------- */
void setWalltmpU_and_U_SPH(auto &net, const std::unordered_set<networkPoint *> &wall_p,
                           const auto &RigidBodyObject) {
   Timer watch;
   for (const auto &p : wall_p)
      p->U_SPH = p->U_SPH_ = p->tmp_U_SPH = p->tmp_U_SPH_ = {0., 0., 0.};

#pragma omp parallel
   for (const auto &p_wall : wall_p)
#pragma omp single nowait
   {
      p_wall->U_SPH_ = p_wall->tmp_U_SPH_ = {0., 0., 0.};
      p_wall->total_weight = 0;
      Tddd Xi = p_wall->X + 2 * p_wall->normal_SPH, mirrored_X, mirror_wall = p_wall->X + p_wall->normal_SPH;
      double weight;
      int count = 0;

      auto func = [&](const networkPoint *p_fluid, const Tddd &whereToSet) {
         p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(whereToSet - Xi), p_wall->radius_SPH * _CSML_), POWER));
         p_wall->U_SPH_ += p_fluid->U_SPH * weight;
         p_wall->tmp_U_SPH_ += p_fluid->tmp_U_SPH * weight;
         //
         mirrored_X = Mirror(whereToSet, mirror_wall, p_wall->normal_SPH);
         p_wall->total_weight += (weight = p_fluid->volume * std::pow(w_Bspline(Norm(mirrored_X - Xi), p_wall->radius_SPH * _CSML_), POWER));
         p_wall->U_SPH_ += Reflect(p_fluid->U_SPH, p_wall->normal_SPH) * weight;
         p_wall->tmp_U_SPH_ += Reflect(p_fluid->tmp_U_SPH, p_wall->normal_SPH) * weight;
         //
         count++;
      };

      net->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &p_fluid) {
         func(p_fluid, p_fluid->X);
         if (SPP_U && p_fluid->isSurface) {
            auto ps = p_fluid->radius_SPH / p_fluid->C_SML;
            auto whereToSet = (p_fluid->X + p_fluid->interpolated_normal_SPH * ps);
            if (canSetSPP(net, RigidBodyObject, p_fluid, whereToSet, ps)) {
               func(p_fluid, whereToSet);
            }
         }
      });
      /*周囲の壁粒子の流速は反転などせずにそのまま読み込む*/
      if (!count) {
         p_wall->U_SPH_ = {0., 0., 0.};
         p_wall->tmp_U_SPH_ = {0., 0., 0.};
      }
   }

   for (const auto &p_wall : wall_p) {
      if (p_wall->total_weight < 1E-20) {
         p_wall->U_SPH_ = {0., 0., 0.};
         p_wall->tmp_U_SPH_ = {0., 0., 0.};
      } else {

         p_wall->U_SPH = p_wall->U_SPH_ / p_wall->total_weight;
         p_wall->tmp_U_SPH = p_wall->tmp_U_SPH_ / p_wall->total_weight;
         // free-slip
         // p_wall->U_SPH = Reflect(p_wall->U_SPH, p_wall->normal_SPH);
         // p_wall->tmp_U_SPH = Reflect(p_wall->tmp_U_SPH, p_wall->normal_SPH);

         // no-slip
         p_wall->U_SPH *= -1.;
         p_wall->tmp_U_SPH *= -1.;
      }
      if (!isFinite(p_wall->U_SPH) || !isFinite(p_wall->tmp_U_SPH))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");
   }
};
/* -------------------------------------------------------------------------- */
void setPressureSPP(Network *net) {
#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      p->p_SPH_SPP = p->total_weight = 0;
      double w = 0;
      net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &p_fluid) {
         // p->total_weight += (w = w_Bspline(Norm(p->X - p_fluid->X), p->radius_SPH));
         p->total_weight += (w = p_fluid->volume * std::pow(w_Bspline(Norm(p->X - p_fluid->X), p->radius_SPH * _CSML_), POWER));
         p->p_SPH_SPP += p->p_SPH * w;
      });
      p->p_SPH_SPP /= p->total_weight;
   }
};
/* -------------------------------------------------------------------------- */
const double damping_factor = 1.;
const double asobi = 0.;

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
      std::cout << Green << "固定の平滑化距離の計算: C_SML * particle_spacing = " << C_SML << " * " << particle_spacing << " = " << C_SML * particle_spacing << colorOff << std::endl;
      for (const auto &p : net->getPoints()) {
         p->setDensity(_WATER_DENSITY_);
         p->isFreeFalling = false;
         p->isInsideOfBody = false;
         p->p_SPH_SPP = 0;
         p->DUDt_SPH = p->lap_U = {0, 0, 0};
      }
      //@ ------------------------------------------------------- */
      //@                    関連する壁粒子をマーク                   */
      //@ ------------------------------------------------------- */
      wall_p.clear();
      DebugPrint("関連する壁粒子をマーク", Green);
      double C = 1.2;
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            p->isCaptured = false;
            p->setDensityVolume(_WATER_DENSITY_, 1E-20);
            p->isFreeFalling = false;
            p->isSurface = false;
            p->p_SPH = 1E+50;
            p->U_SPH = p->DUDt_SPH = p->lap_U = {0, 0, 0};
         }

      for (const auto &[obj, poly] : RigidBodyObject) {
#pragma omp parallel
         for (const auto &p : net->getPoints())
#pragma omp single nowait
         {
            obj->BucketPoints.apply(p->X, p->radius_SPH * C, [&](const auto &q) {
               if (Distance(p, q) < p->radius_SPH * C) {
                  q->isCaptured = true;
                  q->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3.));
                  q->p_SPH = 0;
                  q->U_SPH = {0, 0, 0};
               }
            });
         }
      };

      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints())
            if (p->isCaptured)
               wall_p.emplace(p);

      std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "関連する壁粒子をマークし，保存" << colorOff << std::endl;
      //% -------------- バケットの生成, p->radius_SPHの範囲だけ点を取得 --------------- */
      DebugPrint("バケットの生成", Green);
      net->makeBucketPoints(particle_spacing);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->makeBucketPoints(particle_spacing);
      std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "バケットの生成" << colorOff << std::endl;
      //% -------------------------------------------------------------------------- */
      //% --------------- CFL条件を満たすようにタイムステップ間隔dtを設定 ----------------- */
      //% -------------------------------------------------------------------------- */
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
               //  double max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH);
               //  if (dt > max_dt_vel && isFinite(max_dt_vel))
               //     dt = max_dt_vel;
               /* ------------------------------------------------ */
               // 相対速度
               double max_dt_acc = C_CFL_accel * std::sqrt(distance / std::abs(Dot(p->DUDt_SPH - q->DUDt_SPH, pq)));
               // double max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH - q->DUDt_SPH));
               if (dt > max_dt_acc && isFinite(max_dt_acc))
                  dt = max_dt_acc;
               // 絶対速度
               // double max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH));
               // if (dt > max_dt_acc && isFinite(max_dt_acc))
               //    dt = max_dt_acc;
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
         //% --------------------------- 平滑化距離の計算 ------------------------------ */
         /*     密度, 平滑化距離      */
         std::cout << Green << "固定の平滑化距離の計算: C_SML * particle_spacing = " << C_SML << " * " << particle_spacing << " = " << C_SML * particle_spacing << colorOff << std::endl;
         for (const auto &p : net->getPoints()) {
            p->setDensity(_WATER_DENSITY_);
            p->isFreeFalling = false;
            p->isInsideOfBody = false;
         }

         // b@ ======================================================= */
         // b@                  ルンゲクッタを使った時間積分                  */
         // b@ ======================================================= */
         dt = (*net->getPoints().begin())->RK_X.getdt();
         std::cout << "dt = " << dt << std::endl;

         auto setNormal_Surface = [&](const bool set_surface = true) {
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
               p->interpolated_normal_SPH_original = {0., 0., 0.};
               p->interpolated_normal_SPH_original_all = {0., 0., 0.};
               net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                  if (Between(Distance(p, q), {1E-10, p->radius_SPH})) {
                     p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     p->interpolated_normal_SPH_original_all -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                  }
               });
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                     p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     p->interpolated_normal_SPH_original_all -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                  });

               p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
               p->interpolated_normal_SPH_all = Normalize(p->interpolated_normal_SPH_original_all);
               if (!isFinite(p->interpolated_normal_SPH)) {
                  p->interpolated_normal_SPH = {0., 0., 1.};
                  p->interpolated_normal_SPH_all = {0., 0., 1.};
               }
               /* ----------------------------------- 検索 ----------------------------------- */
               if (set_surface) {
                  p->isSurface = true;
                  if (net->BucketPoints.any_of(p->X, (p->radius_SPH / p->C_SML) * 3., [&](const auto &q) {
                         if (Distance(p, q) < (p->radius_SPH / p->C_SML) * 3.) {
                            return p != q && (VectorAngle(p->interpolated_normal_SPH_all, q->X - p->X) < M_PI / 4);
                         } else
                            return false;
                      }))
                     p->isSurface = false;

                  if (p->isSurface)
                     for (const auto &[obj, poly] : RigidBodyObject)
                        if (obj->BucketPoints.any_of(p->X, (p->radius_SPH / p->C_SML) * 2.5, [&](const auto &q) {
                               if (Distance(p, q) < (p->radius_SPH / p->C_SML) * 2.5) {
                                  return p != q && (VectorAngle(p->interpolated_normal_SPH_all, -q->normal_SPH) < M_PI / 180. * 60);
                               } else
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
         //
         /* -------------------------------------------------------------------------- */
         std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << " 流体粒子の法線方向の計算，水面の判定" << colorOff << std::endl;
         // b! ========================================================================== */
         // b!                          U, DUDt, DPDtを計算                                */
         // b! ========================================================================== */
         {
            /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/
            // ! ------------------------------------------------------ */
            // !                    ∇.∇UとU*を計算                       */
            // ! ------------------------------------------------------ */
            setNormal_Surface();
            setWallDensityVolume(net, wall_p, RigidBodyObject);
            setNormal_Surface();
            DebugPrint("粘性項の∇.∇Uを計算し，次にU*を計算", Magenta);
            setWalltmpU_and_U_SPH(net, wall_p, RigidBodyObject);
            /* -------------------------------------------------------------------------- */
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               Tddd viscous_term = p->lap_U = {0, 0, 0}, rij, qX;
               double nu;
               auto funcOffunc = [&](const networkPoint *q, const Tddd &qX) {
                  rij = p->X - qX;
                  nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
                  viscous_term += q->mass * 8 * nu * Dot(p->U_SPH - q->U_SPH, rij) * grad_w_Bspline(p->X, qX, p->radius_SPH) / ((q->rho + p->rho) * Dot(rij, rij) + std::pow(1E-4 * p->radius_SPH, 2));
               };

               auto func = [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && !q->isCaptured)
                     return;
                  if (p != q)
                     funcOffunc(q, q->X);
                  if (SPP_U && q->isSurface) {
                     qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
                     if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-10, p->radius_SPH}))
                        funcOffunc(q, qX);
                  }
               };
               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);
               p->lap_U = p->lap_U_ = viscous_term / (p->mu_SPH / p->rho);
               p->DUDt_SPH = p->DUDt_SPH_ = viscous_term + _GRAVITY3_;  // 後で修正されるDUDt
               p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
               p->tmp_X = p->X + p->tmp_U_SPH * dt;
               // p->tmp_X = p->X + p->U_SPH * dt + p->DUDt_SPH * dt * dt / 2.;
               // b# -------------------------------------------------------------------------- */
               // b#                           DUDt_SPHとtmp_U_SPHの修正                          */
               // b# -------------------------------------------------------------------------- */

#if defined(REFLECTION)
               bool isReflect = true;
               while (isReflect) {
                  const auto X = p->tmp_X;
                  isReflect = false;
                  //
                  auto closest = [&]() {
                     double distance = 1E+20;
                     networkPoint *P = nullptr;
                     for (const auto &[obj, poly] : RigidBodyObject) {
                        obj->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
                           auto tmp = Distance(X, q);
                           if (distance > tmp) {
                              distance = tmp;
                              P = q;
                           }
                        });
                     }
                     return P;
                  };
                  auto p_wall = closest();
                  if (p_wall) {
                     const auto n = Normalize(p_wall->normal_SPH);
                     if (particle_spacing > Distance(p_wall, X))
                        if (Dot(p_wall->X + p_wall->normal_SPH + n * particle_spacing * (1 / 2. + asobi) - X, n) > 0) {
                           if (Dot(p->tmp_U_SPH, n) < 0) {
                              p->DUDt_SPH -= (1 + damping_factor) * Projection(p->tmp_U_SPH, n) / dt;
                              p->DUDt_SPH_ = p->DUDt_SPH;
                              p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
                              p->tmp_X = p->X + p->tmp_U_SPH * dt;
                              isReflect = true;
                           }
                        }
                  }
               };
#endif
            }
            /* -------------------------------------------------------------------------- */
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粘性項の∇.∇Uを計算し，次にU*を計算" << colorOff << std::endl;
            /* --------------------------------------------------------- */
            DebugPrint("U*を使って，仮想的な位置X^*へ粒子を移動", Green);

#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
               p->setX(p->tmp_X);
            // p->setX(p->X + p->tmp_U_SPH * dt);  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            /* -------------------------------------------------------------------------- */
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粘性項の∇.∇Uを計算し，次にU*を計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !               div(U^*), DρDt=-ρdiv(U)の計算              */
            // ! ------------------------------------------------------ */
            setWallDensityVolume(net, wall_p, RigidBodyObject);
            setNormal_Surface();
            setWalltmpU_and_U_SPH(net, wall_p, RigidBodyObject);
            DebugPrint("div(U^*)を計算", Magenta);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               p->div_U = p->DrhoDt_SPH = 0.;
               Tddd qp, qX;
               auto func = [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && !q->isCaptured)
                     return;
                  //
                  if (p != q) {
                     // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
                     p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
                     // qp = q->X - p->X;
                     // p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH) / 2.;
                  }
                  //
                  if (SPP_U && q->isSurface) {
                     qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
                     // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
                     if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML)) {
                        // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
                        if (Between(Distance(p, qX), {1E-10, p->radius_SPH})) {
                           // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, qX, p->radius_SPH));
                           p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, qX, p->radius_SPH));
                           //
                           // qp = qX - p->X;
                           // p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH) / 2.;
                        }
                     }
                  }
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               if (!isFinite(p->div_U))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "div_U is not a finite");

               p->DrhoDt_SPH = -p->rho * p->div_U;
            }

            for (const auto &p : net->getPoints()) {
               p->setDensity(p->rho_ = p->rho + p->DrhoDt_SPH * dt);  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               // p->setDensity(_WATER_DENSITY_ * (1. - p->div_U * dt));
            }
            /* -------------------------------- もう一度 ---------------------------------- */
            setWallDensityVolume(net, wall_p, RigidBodyObject);
            setWalltmpU_and_U_SPH(net, wall_p, RigidBodyObject);

#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               p->div_U = p->DrhoDt_SPH = 0.;
               Tddd qp, qX;
               auto func = [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && !q->isCaptured)
                     return;
                  //
                  if (p != q) {
                     // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
                     p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
                     // qp = q->X - p->X;
                     // p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH);
                  }
                  //
                  if (SPP_U && q->isSurface) {
                     qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
                     // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
                     if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML)) {
                        // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
                        if (Between(Distance(p, qX), {1E-10, p->radius_SPH})) {
                           // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, qX, p->radius_SPH));
                           p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, qX, p->radius_SPH));
                           //
                           // qp = qX - p->X;
                           // p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH);
                        }
                     }
                  }
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               if (!isFinite(p->div_U))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "div_U is not a finite");

               p->DrhoDt_SPH = -p->rho * p->div_U;
            }

            // #pragma omp parallel
            //             for (const auto &p : net->getPoints())
            // #pragma omp single nowait
            //             {
            //                Tddd viscous_term = {0, 0, 0};

            //                auto funcOffunc = [&](const networkPoint *q, const Tddd &qX) {
            //                   auto rij = p->X - qX;
            //                   // auto Uij = p->U_SPH - q->U_SPH;
            //                   auto nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
            //                   auto tmp = q->mass * 8 * nu * Dot(p->tmp_U_SPH - q->tmp_U_SPH, rij) * grad_w_Bspline(p->X, qX, p->radius_SPH);
            //                   tmp /= (q->rho + p->rho) * Dot(rij, rij) + std::pow(1E-4 * p->radius_SPH, 2);
            //                   viscous_term += tmp;
            //                };

            //                auto func = [&](const auto &q) {
            //                   if (q->getNetwork()->isRigidBody && !q->isCaptured)
            //                      return;
            //                   if (p != q)
            //                      funcOffunc(q, q->X);
            //                   if (SPP_U && q->isSurface) {
            //                      // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
            //                      auto qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
            //                      // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
            //                      if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-10, p->radius_SPH}))
            //                         funcOffunc(q, qX);
            //                   }
            //                };
            //                net->BucketPoints.apply(p->X, p->radius_SPH, func);
            //                for (const auto &[obj, poly] : RigidBodyObject)
            //                   obj->BucketPoints.apply(p->X, p->radius_SPH, func);
            //                p->lap_U_ = viscous_term / (p->mu_SPH / p->rho);
            //                p->DUDt_SPH_ = viscous_term + _GRAVITY3_;  // 後で修正されるDUDt
            //                p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH_ * dt;
            //                p->tmp_X = p->X + p->tmp_U_SPH * dt;

            // #if defined(REFLECTION)
            //                bool isReflect = true;
            //                while (isReflect) {
            //                   // const auto X = p->X + p->tmp_U_SPH * dt;
            //                   const auto X = p->tmp_X;
            //                   isReflect = false;
            //                   //
            //                   auto closest = [&]() {
            //                      double distance = 1E+20;
            //                      networkPoint *P = nullptr;
            //                      for (const auto &[obj, poly] : RigidBodyObject) {
            //                         obj->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
            //                            auto tmp = Distance(X, q);
            //                            if (distance > tmp) {
            //                               distance = tmp;
            //                               P = q;
            //                            }
            //                         });
            //                      }
            //                      return P;
            //                   };
            //                   auto p_wall = closest();
            //                   const double damping_factor = 1.;
            //                   if (p_wall) {
            //                      const auto n = Normalize(p_wall->normal_SPH);
            //                      if (particle_spacing > Distance(p_wall, X))
            //                         if (Dot(p_wall->X + p_wall->normal_SPH + n * (particle_spacing * 0.5) - X, n) > 0) {
            //                            if (Dot(p->tmp_U_SPH, n) < 0) {
            //                               p->DUDt_SPH_ -= (1 + damping_factor) * Projection(p->tmp_U_SPH, n) / dt;
            //                               p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH_ * dt;
            //                               p->tmp_X = p->X + p->U_SPH * dt;
            //                               // p->tmp_X = Mirror(p->tmp_X, p_wall->X + p_wall->normal_SPH + n * particle_spacing / 2, p_wall->normal_SPH);
            //                               // p->tmp_X = p->X + p->U_SPH * dt + p->DUDt_SPH * dt * dt / 2.;
            //                               isReflect = true;
            //                            }
            //                         }
            //                   }
            //                };
            // #endif
            //             }
            /* -------------------------------------------------------------------------- */
            setWallDensityVolume(net, wall_p, RigidBodyObject);
            setWalltmpU_and_U_SPH(net, wall_p, RigidBodyObject);
            setWallPressure(net, wall_p, RigidBodyObject, 0);
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "div(U^*)を計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !  　　　　　        仮位置における圧力Pの計算                 */
            // ! ------------------------------------------------------ */
            /*
            ここで，疑問
            流体の圧力は，仮の位置に移動する前のもの．壁の圧力も移動前のもの．
            仮の位置に移動した後に，流体の圧力を計算するが，そこで利用する圧力は，移動前の圧力でよいのか？
            */

            // #define Morikawa2019

            DebugPrint("仮位置における圧力Pの計算", Magenta);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               int count = 0;
               Tddd Xij, qX;
               double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
               auto func = [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && !q->isCaptured)
                     return;
                  if (p != q) {
                     Xij = p->X - q->X;
#if defined(Morikawa2019)
                     // Morikawa, D., Senadheera, H., & Asai, M. (2021). Explicit incompressible smoothed particle hydrodynamics in a multi-GPU environment for large-scale simulations. Computational Particle Mechanics, 8(3), 493–510. https://doi.org/10.1007/s40571-020-00347-0
                     Aij = 2. * q->mass / p->rho * Dot(Xij, grad_w_Bspline(p->X, q->X, p->radius_SPH)) / Dot(Xij, Xij);
                     // Aij = 2. * q->volume * Dot(Xij, grad_w_Bspline(p->X, q->X, p->radius_SPH)) / Dot(Xij, Xij);
#else
                     Aij = q->mass * 8. / std::pow(q->rho + p->rho, 2) * Dot(Xij, grad_w_Bspline(p->X, q->X, p->radius_SPH)) / Dot(Xij, Xij);
#endif
                     sum_Aij += Aij;
                     sum_Aij_Pj += Aij * q->p_SPH;
                     count++;
                  }
                  if (SPP_pressure && q->isSurface) {
                     // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
                     qX = (q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML);
                     // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
                     if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-8, p->radius_SPH})) {
                        Xij = p->X - qX;
#if defined(Morikawa2019)
                        Aij = 2. * q->mass / p->rho * Dot(Xij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(Xij, Xij);
                        // Aij = 2. * q->volume * Dot(Xij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(Xij, Xij);
#else
                        Aij = q->mass * 8. / std::pow(q->rho + p->rho, 2) * Dot(Xij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(Xij, Xij);
#endif
                        sum_Aij += Aij;
                        sum_Aij_Pj += Aij * SPP_pressure_coef * q->p_SPH_SPP;
                     }
                  }
                  //% -------------------------------------------------------------------------- */
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               if (std::abs(sum_Aij) < 1E-20)
                  p->p_SPH_ = 0;
               else {
#if defined(Morikawa2019)
                  p->p_SPH_ = (_WATER_DENSITY_ * p->div_U / dt + sum_Aij_Pj) / sum_Aij;
#else
                  p->p_SPH_ = (p->div_U / dt + sum_Aij_Pj) / sum_Aij;
#endif
               }

               if (!isFinite(p->p_SPH_)) {
                  std::cout << "p->div_U " << p->div_U << std::endl;
                  std::cout << "sum_Aij_Pj " << sum_Aij_Pj << std::endl;
                  std::cout << "sum_Aij " << sum_Aij << std::endl;
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");
               }
            }

            for (const auto &p : net->getPoints()) {
               p->DPDt_SPH = (p->p_SPH_ - p->p_SPH) / dt;
               // p->p_SPH = p->p_SPH_;
               p->RK_P.push(p->DPDt_SPH);  // 圧力
               // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
               p->p_SPH = p->p_SPH_;
            }

            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "仮位置における圧力Pの計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
            // ! ------------------------------------------------------ */

            DebugPrint("圧力勾配∇Pを計算 & DU/Dtの計算", Magenta);
            setWallDensityVolume(net, wall_p, RigidBodyObject);
            setWalltmpU_and_U_SPH(net, wall_p, RigidBodyObject);
            setWallDUDt(net, wall_p, RigidBodyObject);
            setWallPressure(net, wall_p, RigidBodyObject, 1);

#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               double p_rho2;
               int count = 0;
               Tddd qX;
               p->contact_points_all_SPH = p->contact_points_fluid_SPH = 0;
               p->gradP_SPH = {0., 0., 0.};

               auto func = [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && !q->isCaptured && Distance(q, p) >= p->radius_SPH)
                     return;
                  //
                  p->contact_points_all_SPH++;
                  if (q->getNetwork()->isFluid)
                     p->contact_points_fluid_SPH++;

                  if (p != q) {
                     p_rho2 = q->p_SPH / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                     p->gradP_SPH += p->rho * q->mass * p_rho2 * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     count++;
                  }
                  if (SPP_pressure && q->isSurface) {
                     qX = q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML;
                     if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-8, p->radius_SPH})) {
                        p_rho2 = SPP_pressure_coef * q->p_SPH_SPP / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                        p->gradP_SPH += p->rho * q->mass * p_rho2 * grad_w_Bspline(p->X, qX, p->radius_SPH);
                        count++;
                     }
                  }
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               p->DUDt_SPH -= p->gradP_SPH / p->rho;
               //
               if (!isFinite(p->DUDt_SPH))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
            }
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "圧力勾配∇Pを計算 & DU/Dtの計算" << colorOff << std::endl;
            /* -------------------------------------------------------------------------- */
            setPressureSPP(net);
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
            /* -------------------------------------------------------------------------- */

#if defined(REFLECTION)
            bool isReflect = true;
            while (isReflect) {
               // const auto X = p->RK_X.getX(p->U_SPH);
               const auto X = p->X;
               isReflect = false;

               auto closest = [&]() {
                  double distance = 1E+20;
                  networkPoint *P = nullptr;
                  for (const auto &[obj, poly] : RigidBodyObject) {
                     obj->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
                        auto tmp = Distance(X, q);
                        if (p != q && distance > tmp) {
                           distance = tmp;
                           P = q;
                        }
                     });
                  }
                  return P;
               };
               auto p_wall = closest();
               if (p_wall) {
                  auto n = Normalize(p_wall->normal_SPH);
                  if (particle_spacing > Distance(p_wall, X))
                     if (Dot(p_wall->X + p_wall->normal_SPH + n * particle_spacing * (1 / 2. + asobi) - X, n) > 0) {
                        if (Dot(p->U_SPH, n) < 0) {
                           p->DUDt_SPH -= (1 + damping_factor) * Projection(p->U_SPH, n) / dt;
                           //
                           p->RK_U.repush(p->DUDt_SPH);  // 速度
                           p->U_SPH = p->RK_U.getX();
                           p->RK_X.repush(p->U_SPH);  // 位置
                           p->setXSingle(p->tmp_X = p->RK_X.getX());
                           //
                           isReflect = true;
                        }
                     }
               }
            };
#endif
         }
         std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粒子の時間発展" << colorOff << std::endl;
         real_time = (*net->getPoints().begin())->RK_X.gett();
      } while (!((*net->getPoints().begin())->RK_X.finished));

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
   std::unordered_map<networkPoint *, double> pressure;
   for (const auto &p : Fluid->getPoints())
      pressure[p] = p->p_SPH;
   vtp.addPointData("pressure", pressure);
};

void setData(auto &vtp, const auto &Fluid) {
   std::unordered_map<networkPoint *, Tddd> U;
   for (const auto &p : Fluid->getPoints())
      U[p] = p->U_SPH;
   vtp.addPointData("U", U);
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