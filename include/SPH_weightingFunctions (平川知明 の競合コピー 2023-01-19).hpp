#ifndef SPH_weightingFunctions_H
#define SPH_weightingFunctions_H

#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

Tddd Reflect(const Tddd &v, const Tddd &n) {
   return v - 2. * n * Dot(v, n);
};

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
   std::unordered_map<networkPoint *, double> contact_points_SPH;
   for (const auto &p : Fluid->getPoints())
      contact_points_SPH[p] = p->contact_points_SPH;
   vtp.addPointData("contact_points_SPH", contact_points_SPH);
};

#define REFLECTION
#define OurWall_LoopNumber 1

#define SPP_setRigidBodyObject_Pressure 1
#define SPP_setRigidBodyObject_tmpU_and_U_SPH 1
#define SPP_lap_U 1
#define SPP_div_U_and_DrhoDt 1
#define SPP_pressure 1
#define SPP_grad_P 1

#define SPP_lap_U_self 0
#define SPP_div_U_and_DrhoDt_self 0
#define SPP_pressure_self 0
#define SPP_grad_P_self 0

#define surface_zero_pressure 0

const double power = 1.;

auto canSetSPP(Network *net, const auto &RigidBodyObject, const auto &p, const Tddd &X, const double &range) {
   auto func = [&](const auto &q) { return p != q && Distance(q, X) < 1.1 * range; };
   if (net->BucketPoints.any_of(X, range, func))
      return false;
   for (const auto &[obj, poly] : RigidBodyObject)
      if (obj->BucketPoints.any_of(X, range, func))
         return false;
   return true;
};
// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */
void setRigidBodyObject_DensityVolume(auto &net, const std::unordered_set<networkPoint *> &wall_p, const auto &RigidBodyObject) {
   Timer watch;
   for (const auto &[obj, poly] : RigidBodyObject)
      for (const auto &p : obj->getPoints()) {
         p->setDensityVolume(_WATER_DENSITY_, 0.);
      }
   const double a = 0.5;
   double W = 0;
   for (auto i = 0; i < OurWall_LoopNumber; ++i) {
#pragma omp parallel
      for (const auto &p_wall : wall_p)
#pragma omp single nowait
      {
         p_wall->rho_ = p_wall->volume_ = p_wall->total_weight = 0;
         Tddd Xi = p_wall->X + 2 * p_wall->normal_SPH, accel;
         double weight;
         int count = 0;
         net->BucketPoints.apply(Xi, 1.1 * p_wall->radius_SPH, [&](const auto &q) {
            if (Dot(p_wall->normal_SPH, q->X - p_wall->X) < 0.)
               return;

            weight = q->volume * std::pow(w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH), power);
            p_wall->total_weight += weight;
            p_wall->rho_ += q->rho * weight;
            p_wall->volume_ += q->volume * weight;

            count++;
            if (SPP_setRigidBodyObject_Pressure && q->isSurface)
               for (auto I = 1; I <= 1; I++) {
                  auto qX = (q->X + q->interpolated_normal_SPH * I * q->radius_SPH / q->C_SML);
                  if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML)) {
                     weight = q->volume * std::pow(w_Bspline(Norm(qX - Xi), p_wall->radius_SPH), power);
                     p_wall->total_weight += weight;
                     p_wall->rho_ += q->rho * weight;
                     p_wall->volume_ += q->volume * weight;
                     count++;
                     /* ----------------------------------------------------- */
                     auto X = qX - 2 * (Dot(qX - p_wall->X, Normalize(p_wall->normal_SPH)) * Normalize(p_wall->normal_SPH) - p_wall->normal_SPH);
                     weight = q->volume * std::pow(w_Bspline(Norm(X - Xi), p_wall->radius_SPH), power);
                     p_wall->total_weight += weight;
                     p_wall->rho_ += q->rho * weight;
                     p_wall->volume_ += q->volume * weight;
                     count++;
                     /* ----------------------------------------------------- */
                  }
               }
            /* ----------------------------------------------------- */
            auto X = q->X - 2 * (Dot(q->X - p_wall->X, Normalize(p_wall->normal_SPH)) * Normalize(p_wall->normal_SPH) - p_wall->normal_SPH);
            weight = q->volume * std::pow(w_Bspline(Norm(X - Xi), p_wall->radius_SPH), power);
            p_wall->total_weight += weight;
            p_wall->rho_ += q->rho * weight;
            p_wall->volume_ += q->volume * weight;
            /* ----------------------------------------------------- */
         });
      }

      for (const auto &p_wall : wall_p) {
         {
            if (p_wall->total_weight < 1E-15)
               p_wall->setDensityVolume(_WATER_DENSITY_, 1E-15);
            else
               // p_wall->setDensityVolume(p_wall->rho_, p_wall->volume_);
               p_wall->setDensityVolume(p_wall->rho_ / p_wall->total_weight,
                                        p_wall->volume_ / p_wall->total_weight);
         }
      }
   }
   /* -------------------------------------------------------------------------- */
   std::stringstream ss;
#if defined(Morikawa2021)
   ss << Blue << "setRigidBodyObject_Pressure (Morikawa2021) Elapsed time: " << Red << watch();
#elif defined(Morikawa2021_IDW)
   ss << Blue << "setRigidBodyObject_Pressure (Morikawa2021 with IDW) Elapsed time: " << Red << watch();
#else
   ss << Blue << "setRigidBodyObject_Pressure (another method) Elapsed time: " << Red << watch();
#endif
   DebugPrint(ss.str());
};
/* -------------------------------------------------------------------------- */
void setRigidBodyObject_Pressure(auto &net, const std::unordered_set<networkPoint *> &wall_p,
                                 const auto &RigidBodyObject,
                                 const bool do_add) {
   setRigidBodyObject_DensityVolume(net, wall_p, RigidBodyObject);
   Timer watch;
   for (const auto &p : wall_p)
      p->p_SPH__ = p->p_SPH_ = p->p_SPH = p->total_weight = 0;

   const double a = 0.5;
   double W = 0;

   for (auto i = 0; i < OurWall_LoopNumber; ++i) {
#pragma omp parallel
      for (const auto &p_wall : wall_p)
#pragma omp single nowait
      {
         p_wall->p_SPH_ = p_wall->total_weight = 0;
         Tddd Xi = p_wall->X + 2 * p_wall->normal_SPH, accel;
         double weight, nu, add;
         auto n = Normalize(p_wall->normal_SPH);
         int count = 0;
         net->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &q) {
            if (Dot(n, q->X - p_wall->X) < 0.)
               return;

            weight = q->volume * std::pow(w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH), power);
            p_wall->total_weight += weight;
            accel = {0., 0., 0.};
            nu = q->mu_SPH / q->rho;
            // add = q->rho * Dot(nu * q->lap_U + _GRAVITY3_ - accel, -2 * p_wall->normal_SPH);
            add = _WATER_DENSITY_ * Dot(q->DUDt_SPH - accel, -2 * p_wall->normal_SPH);
            if (!do_add)
               add = 0;
            p_wall->p_SPH_ += (q->p_SPH + add) * weight;
            count++;
            if (SPP_setRigidBodyObject_Pressure && q->isSurface)
               for (auto I = 1; I <= 1; I++) {
                  auto qX = (q->X + q->interpolated_normal_SPH * I * q->radius_SPH / q->C_SML);
                  // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
                  if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML)) {
                     weight = q->volume * std::pow(w_Bspline(Norm(qX - Xi), p_wall->radius_SPH), power);
                     p_wall->total_weight += weight;
                     accel = {0., 0., 0.};
                     //
                     nu = q->mu_SPH / q->rho;
                     // add = q->rho * Dot(nu * q->lap_U + _GRAVITY3_ - accel, -2 * p_wall->normal_SPH);
                     add = _WATER_DENSITY_ * Dot(q->DUDt_SPH - accel, -2 * p_wall->normal_SPH);
                     if (!do_add)
                        add = 0;
                     p_wall->p_SPH_ += -(q->p_SPH + add) * weight;
                     count++;
                     /* ----------------------------------------------------- */
                     auto X = qX - 2 * (Dot(qX - p_wall->X, Normalize(p_wall->normal_SPH)) * Normalize(p_wall->normal_SPH) - p_wall->normal_SPH);
                     weight = q->volume * std::pow(w_Bspline(Norm(X - Xi), p_wall->radius_SPH), power);
                     p_wall->total_weight += weight;
                     accel = {0., 0., 0.};
                     //
                     nu = q->mu_SPH / q->rho;
                     // add = q->rho * Dot(nu * q->lap_U + _GRAVITY3_ - accel, -2 * p_wall->normal_SPH);
                     add = _WATER_DENSITY_ * Dot(q->DUDt_SPH - accel, -2 * p_wall->normal_SPH);
                     if (!do_add)
                        add = 0;
                     p_wall->p_SPH_ += -(q->p_SPH + add) * weight;
                     count++;
                     /* ----------------------------------------------------- */
                  }
               }
            /* ----------------------------------------------------- */
            auto X = q->X - 2 * (Dot(q->X - p_wall->X, Normalize(p_wall->normal_SPH)) * Normalize(p_wall->normal_SPH) - p_wall->normal_SPH);
            weight = q->volume * std::pow(w_Bspline(Norm(X - Xi), p_wall->radius_SPH), power);
            p_wall->total_weight += weight;
            add = _WATER_DENSITY_ * Dot(q->DUDt_SPH - accel, -2 * p_wall->normal_SPH);
            if (!do_add)
               add = 0;
            p_wall->p_SPH_ += (q->p_SPH + add) * weight;
            /* ----------------------------------------------------- */
         });
      }
      //
      for (const auto &p_wall : wall_p) {
         {
            if (p_wall->total_weight < 1E-15)
               p_wall->p_SPH = p_wall->p_SPH_;
            else
               // p_wall->p_SPH = p_wall->p_SPH_;
               p_wall->p_SPH = p_wall->p_SPH_ / p_wall->total_weight;
         }
      }
   }
   /* -------------------------------------------------------------------------- */
   std::stringstream ss;
#if defined(Morikawa2021)
   ss << Blue << "setRigidBodyObject_Pressure (Morikawa2021) Elapsed time: " << Red << watch();
#elif defined(Morikawa2021_IDW)
   ss << Blue << "setRigidBodyObject_Pressure (Morikawa2021 with IDW) Elapsed time: " << Red << watch();
#else
   ss << Blue << "setRigidBodyObject_Pressure (another method) Elapsed time: " << Red << watch();
#endif
   DebugPrint(ss.str());
};
/* -------------------------------------------------------------------------- */
void setRigidBodyObject_tmpU_and_U_SPH(auto &net, const std::unordered_set<networkPoint *> &wall_p, const auto &RigidBodyObject) {
   setRigidBodyObject_DensityVolume(net, wall_p, RigidBodyObject);
   Timer watch;
   for (const auto &p : wall_p)
      p->U_SPH = p->U_SPH_ = p->tmp_U_SPH = p->tmp_U_SPH_ = {0., 0., 0.};

   const double a = 0.5;
   double W = 0;
   for (auto i = 0; i < OurWall_LoopNumber; ++i) {
#pragma omp parallel
      for (const auto &p_wall : wall_p)
#pragma omp single nowait
      // if (Norm(p_wall->normal_SPH) <= 1.1 * (OurWall_LoopNumber - i) / OurWall_LoopNumber * p_wall->radius_SPH)
      {
         p_wall->U_SPH_ = p_wall->tmp_U_SPH_ = {0., 0., 0.};
         p_wall->total_weight = 0;
         // p_wall->rho_ = 0.;
         auto Xi = p_wall->X + 2 * p_wall->normal_SPH;
         double weight, l = Norm(p_wall->normal_SPH);
         Tddd U, n = Normalize(p_wall->normal_SPH);
         int count = 0;
         auto func = [&](const auto &q) {
            if (Dot(n, q->X - p_wall->X) < 0.)
               return;

            // free-slip
            // if (Norm(q->X - Xi) <= p_wall->radius_SPH) {
            //    weight = q->volume * w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH);
            //    p_wall->total_weight += weight;
            //    p_wall->U_SPH_ += (q->U_SPH - 2 * Dot(q->U_SPH, n) * n) * weight;
            //    p_wall->tmp_U_SPH_ += (q->tmp_U_SPH - 2 * Dot(q->tmp_U_SPH, n) * n) * weight;
            //    count++;
            //    if (SPP_setRigidBodyObject_tmpU_and_U_SPH && q->isSurface)
            //       for (auto I = 1; I <= 1; I++) {
            //          auto qX = (q->X + q->interpolated_normal_SPH * I * q->radius_SPH / q->C_SML);
            //          // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
            //          if (canSetSPP(net, RigidBodyObject, qX, q->radius_SPH / q->C_SML)) {
            //             weight = q->volume * w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH);
            //             p_wall->total_weight += weight;
            //             p_wall->U_SPH_ += (q->U_SPH - 2 * Dot(q->U_SPH, n) * n) * weight;
            //             p_wall->tmp_U_SPH_ += (q->tmp_U_SPH - 2 * Dot(q->tmp_U_SPH, n) * n) * weight;
            //          }
            //       }
            // }

            // no-slip
            if (Norm(q->X - Xi) <= p_wall->radius_SPH) {
               weight = q->volume * std::pow(w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH), power);
               p_wall->total_weight += weight;
               p_wall->U_SPH_ += q->U_SPH * weight;
               p_wall->tmp_U_SPH_ += q->tmp_U_SPH * weight;
               count++;
               if (SPP_setRigidBodyObject_tmpU_and_U_SPH && q->isSurface)
                  for (auto I = 1; I <= 1; I++) {
                     auto qX = (q->X + q->interpolated_normal_SPH * I * q->radius_SPH / q->C_SML);
                     // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
                     if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML)) {
                        weight = q->volume * std::pow(w_Bspline(Norm(qX - Xi), p_wall->radius_SPH), power);
                        p_wall->total_weight += weight;
                        p_wall->U_SPH_ += q->U_SPH * weight;
                        p_wall->tmp_U_SPH_ += q->tmp_U_SPH * weight;
                        count++;
                        /* ----------------------------------------------------- */
                        auto X = qX - 2 * (Dot(qX - p_wall->X, Normalize(p_wall->normal_SPH)) * Normalize(p_wall->normal_SPH) - p_wall->normal_SPH);
                        weight = q->volume * std::pow(w_Bspline(Norm(X - Xi), p_wall->radius_SPH), power);
                        p_wall->total_weight += weight;
                        p_wall->U_SPH_ += -q->U_SPH * weight;
                        p_wall->tmp_U_SPH_ += -q->tmp_U_SPH * weight;
                        /* ----------------------------------------------------- */
                     }
                  }
            }

            /* ----------------------------------------------------- */
            auto X = q->X - 2 * (Dot(q->X - p_wall->X, Normalize(p_wall->normal_SPH)) * Normalize(p_wall->normal_SPH) - p_wall->normal_SPH);
            weight = q->volume * std::pow(w_Bspline(Norm(X - Xi), p_wall->radius_SPH), power);
            p_wall->total_weight += weight;
            p_wall->U_SPH_ += -q->U_SPH * weight;
            p_wall->tmp_U_SPH_ += -q->tmp_U_SPH * weight;
            /* ----------------------------------------------------- */
         };
         net->BucketPoints.apply(Xi, p_wall->radius_SPH, func);
         /*周囲の壁粒子の流速は反転などせずにそのまま読み込む*/
         if (!count) {
            p_wall->U_SPH_ = {0., 0., 0.};
            p_wall->tmp_U_SPH_ = {0., 0., 0.};
         }
      }

      for (const auto &p_wall : wall_p) {
         if (p_wall->total_weight < 1E-15) {
            p_wall->U_SPH = -p_wall->U_SPH_;
            p_wall->tmp_U_SPH = -p_wall->tmp_U_SPH_;
         } else {
            // p_wall->U_SPH = p_wall->U_SPH_;
            // p_wall->tmp_U_SPH = p_wall->tmp_U_SPH_;
            p_wall->U_SPH = -p_wall->U_SPH_ / p_wall->total_weight;
            p_wall->tmp_U_SPH = -p_wall->tmp_U_SPH_ / p_wall->total_weight;
         }
      }
   }

   std::stringstream ss;
#if defined(Morikawa2021)
   ss << Blue << "setRigidBodyObject_U_SPH (Morikawa2021) Elapsed time: " << Red << watch();
#elif defined(Morikawa2021_IDW)
   ss << Blue << "setRigidBodyObject_U_SPH (Morikawa2021 with IDW) Elapsed time: " << Red << watch();
#else
   ss << Blue << "setRigidBodyObject_U_SPH (another method) Elapsed time: " << Red << watch();
#endif
   DebugPrint(ss.str());
};
/* -------------------------------------------------------------------------- */
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
      }
      /* -------------------------------------------------------------------------- */
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            // p->setDensity(_WATER_DENSITY_);
            // p->p_SPH = 0;
            // p->U_SPH = {0, 0, 0};
            p->setDensityVolume(_WATER_DENSITY_, 1E-20);
            p->isFreeFalling = false;
            p->p_SPH = 1E+50;
            p->U_SPH = {0, 0, 0};
         }
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
               if (dt > max_dt_vel && isFinite(max_dt_vel))
                  dt = max_dt_vel;
               // 絶対速度
               //  double max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH);
               //  if (dt > max_dt_vel && isFinite(max_dt_vel))
               //     dt = max_dt_vel;
               /* ------------------------------------------------ */
               // 相対速度
               double max_dt_acc = C_CFL_accel * std::sqrt(distance / std::abs(Dot(p->DUDt_SPH - q->DUDt_SPH, pq)));
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
      //@ ------------------------------------------------------- */
      //@                    関連する壁粒子をマーク                   */
      //@ ------------------------------------------------------- */
      DebugPrint("関連する壁粒子をマーク", Green);
      std::unordered_set<networkPoint *> wall_p;
      double C = 1.;
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints())
            p->isCaptured = false;
      for (const auto &[obj, poly] : RigidBodyObject) {
#pragma omp parallel
         for (const auto &p : net->getPoints())
#pragma omp single nowait
         {
            obj->BucketPoints.apply(p->X, p->radius_SPH * C, [&](const auto &q) {
               if (Distance(p, q) > p->radius_SPH * 1.1 || Dot(q->normal_SPH, p->X - q->X) < 0.) {
                  return;
               } else {
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
      //@ ----------------------- ルンゲクッタの準備 ------------------- */
      // int RK_order = 2;

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
            if (set_surface) {
               DebugPrint("水粒子のオブジェクト外向き法線方向を計算", Green);
               for (const auto &p : net->getPoints())
                  p->isSurface = false;
            }
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               // p->interpolated_normal_SPH, q->X - p->Xの方向が完全に一致した際に失敗する
               /* ---------------------- p->interpolated_normal_SPHの計算 --------------------- */
               p->interpolated_normal_SPH_original = {0., 0., 0.};
               p->interpolated_normal_SPH_original_all = {0., 0., 0.};
               net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                  if (Between(Distance(p, q), {1E-8, p->radius_SPH})) {
                     p->interpolated_normal_SPH_original -= grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     p->interpolated_normal_SPH_original_all -= grad_w_Bspline(p->X, q->X, p->radius_SPH);
                  }
               });
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                     p->interpolated_normal_SPH_original -= 0.8 * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     p->interpolated_normal_SPH_original_all -= grad_w_Bspline(p->X, q->X, p->radius_SPH);
                  });

               p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
               p->interpolated_normal_SPH_all = Normalize(p->interpolated_normal_SPH_original_all);
               if (!isFinite(p->interpolated_normal_SPH)) {
                  p->interpolated_normal_SPH = {0., 0., 1.};
                  p->interpolated_normal_SPH_all = {0., 0., 1.};
               }
               // if (p->RK_U.current_step == 0) {
               /* ----------------------------------- 検索 ----------------------------------- */
               // double r = (p->radius_SPH / C_SML) * 2.5;  // チェック範囲 : particle spacing * 2.5
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
                        if (obj->BucketPoints.any_of(p->X, (p->radius_SPH / p->C_SML) * 3., [&](const auto &q) {
                               if (Distance(p, q) < (p->radius_SPH / p->C_SML) * 3.) {
                                  return p != q && (VectorAngle(p->interpolated_normal_SPH_all, q->X - p->X) < M_PI / 5);
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
                     obj->BucketPoints.apply(p->X, p->radius_SPH * C, [&](const auto &q) {
                        if (Between(Distance(p, q), {1E-8, p->radius_SPH}))
                           p->interpolated_normal_SPH_original -= grad_w_Bspline(p->X, q->X, p->radius_SPH * C);
                     });
                  p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
                  if (!isFinite(p->interpolated_normal_SPH))
                     p->interpolated_normal_SPH = {0., 0., 1.};
               }
            }
         };
         setNormal_Surface();
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
            DebugPrint("粘性項の∇.∇Uを計算し，次にU*を計算", Magenta);
            setRigidBodyObject_tmpU_and_U_SPH(net, wall_p, RigidBodyObject);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               Tddd viscous_term = p->lap_U = {0, 0, 0};
               auto func = [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && Dot(q->normal_SPH, p->X - q->X) < 0.)
                     return;
                  if (Between(Distance(p, q), {1E-8, p->radius_SPH})) {
                     auto rij = p->X - q->X;
                     auto Uij = p->U_SPH - q->U_SPH;
                     auto nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
                     //
                     auto tmp = q->mass * 8 * nu * Dot(p->U_SPH - q->U_SPH, rij) * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     tmp /= (q->rho + p->rho) * Dot(rij, rij) + std::pow(1E-4 * p->radius_SPH, 2);
                     viscous_term += tmp;
                     //
                     // auto tmp = 2 / p->rho * q->mass * Dot(rij, grad_w_Bspline(p->X, q->X, p->radius_SPH)) / Dot(rij, rij) * Uij;
                     // viscous_term += p->mu_SPH / p->rho * tmp;
                     //
                     // auto tmp = p->rho * 8. * q->mass / std::pow(q->rho + p->rho, 2) * Uij * Dot(rij, grad_w_Bspline(p->X, q->X, p->radius_SPH)) / Dot(rij, rij);  // div(grad(u))/rho
                     // viscous_term += p->mu_SPH / p->rho * tmp;
                  }
                  if ((SPP_lap_U && q->isSurface) || (SPP_lap_U_self && q->isSurface && p == q))
                     for (auto I = 1; I <= 1; ++I) {
                        // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
                        auto qX = (q->X + q->interpolated_normal_SPH * I * q->radius_SPH / q->C_SML);
                        // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
                        if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-8, p->radius_SPH})) {
                           auto rij = p->X - qX;
                           auto Uij = p->U_SPH - q->U_SPH;
                           auto nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
                           //
                           auto tmp = q->mass * 8 * nu * Dot(p->U_SPH - q->U_SPH, rij) * grad_w_Bspline(p->X, qX, p->radius_SPH);
                           tmp /= (q->rho + p->rho) * Dot(rij, rij) + std::pow(1E-4 * p->radius_SPH, 2);
                           viscous_term += tmp;
                           //
                           // auto tmp = 2 / p->rho * q->mass * Dot(rij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(rij, rij) * Uij;
                           // viscous_term += p->mu_SPH / p->rho * tmp;
                        }
                     }
               };
               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);
               p->lap_U = viscous_term / (p->mu_SPH / p->rho);
               p->DUDt_SPH = viscous_term + _GRAVITY3_;  // 後で修正されるDUDt
               p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
               p->tmp_X = p->X + p->tmp_U_SPH * dt;
// b# -------------------------------------------------------------------------- */
// b#                           DUDt_SPHとtmp_U_SPHの修正                          */
// b# -------------------------------------------------------------------------- */
#if defined(REFLECTION)
               bool isReflect = true;
               while (isReflect) {
                  isReflect = false;
                  auto closest = [&]() {
                     // auto X = p->X + p->tmp_U_SPH * dt;
                     auto X = p->tmp_X;
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
                  auto q = closest();
                  if (q) {
                     // auto n = q->interpolated_normal_SPH;
                     auto n = Normalize(q->normal_SPH);
                     auto Ureflect = -2. * Dot(p->tmp_U_SPH, n) * n;
                     // auto X = p->X + p->tmp_U_SPH * dt;
                     auto X = p->tmp_X;
                     // if (Distance(q, X) < std::sqrt(2.) * particle_spacing && std::abs(Dot(n, q->X - X)) < particle_spacing && Dot(p->tmp_U_SPH, n) < 0) {
                     if (Distance(q, X) < particle_spacing && Dot(p->tmp_U_SPH, n) < 0) {
                        p->DUDt_SPH += Ureflect / dt;
                        p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
                        p->tmp_X += 2 * Dot((q->X + q->normal_SPH + n * particle_spacing / 2.) - p->tmp_X, n);
                        isReflect = true;
                     }
                  }
               };
#endif
               /* -------------------------------------------------------------------------- */
            }
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粘性項の∇.∇Uを計算し，次にU*を計算" << colorOff << std::endl;
            /* --------------------------------------------------------- */
            DebugPrint("U*を使って，仮想的な位置X^*へ粒子を移動", Green);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
               p->setX(p->tmp_X);
            // p->setX(p->X + p->tmp_U_SPH * dt);
            /* -------------------------------------------------------------------------- */
            setNormal_Surface();
            /* -------------------------------------------------------------------------- */
            setRigidBodyObject_tmpU_and_U_SPH(net, wall_p, RigidBodyObject);
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粘性項の∇.∇Uを計算し，次にU*を計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !               div(U^*), DρDt=-ρdiv(U)の計算              */
            // ! ------------------------------------------------------ */
            DebugPrint("div(U^*)を計算", Magenta);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               p->div_U = p->DrhoDt_SPH = 0.;
               Tddd qp;
               auto func = [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && Dot(q->normal_SPH, p->X - q->X) < 0.)
                     return;

                  if (Between(Distance(p, q), {1E-8, p->radius_SPH})) {
                     p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
                     //
                     // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
                     //
                     // qp = q->X - p->X;
                     // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH);
                  }
                  if ((SPP_div_U_and_DrhoDt && q->isSurface) || (SPP_div_U_and_DrhoDt_self && q == p && q->isSurface))
                     for (auto I = 1; I <= 1; ++I) {
                        auto qX = (q->X + q->interpolated_normal_SPH * I * q->radius_SPH / q->C_SML);
                        // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
                        if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-8, p->radius_SPH})) {
                           // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
                           p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, qX, p->radius_SPH));
                           //
                           // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, qX, p->radius_SPH));
                           //
                           // qp = qX - p->X;
                           // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH);
                        }
                     }
               };
               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               if (!isFinite(p->div_U))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "div_U is not a finite");

               p->DrhoDt_SPH = -p->rho * p->div_U;
               // p->setDensity(p->rho + p->DrhoDt_SPH * dt);
               p->setDensity(_WATER_DENSITY_ * (1. - p->div_U * dt));
            }
            /* -------------------------------------------------------------------------- */
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "div(U^*)を計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !  　　　　　        仮位置における圧力Pの計算                 */
            // ! ------------------------------------------------------ */

#define Morikawa2019

            DebugPrint("流体粒子をもとに壁粒子の圧力の計算", Magenta);
            setRigidBodyObject_Pressure(net, wall_p, RigidBodyObject, false);  // 内部では，次回利用する壁粒子の数を少し広くしている(p->isCapture)．
            DebugPrint("仮位置における圧力Pの計算", Magenta);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               int count = 0;
               Tddd Xij;
               double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
               auto func = [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && Dot(q->normal_SPH, p->X - q->X) < 0.)
                     return;

                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  if (Between(Distance(p, q), {1E-8, p->radius_SPH})) {
                     Xij = p->X - q->X;
#if defined(Morikawa2019)
                     // Morikawa, D., Senadheera, H., & Asai, M. (2021). Explicit incompressible smoothed particle hydrodynamics in a multi-GPU environment for large-scale simulations. Computational Particle Mechanics, 8(3), 493–510. https://doi.org/10.1007/s40571-020-00347-0
                     Aij = 2. * q->mass / p->rho * Dot(Xij, grad_w_Bspline(p->X, q->X, p->radius_SPH)) / Dot(Xij, Xij);
#else
                     Aij = q->mass * 8. / std::pow(q->rho + p->rho, 2) * Dot(Xij, grad_w_Bspline(p->X, q->X, p->radius_SPH)) / Dot(Xij, Xij);
#endif
                     sum_Aij += Aij;
                     sum_Aij_Pj += Aij * q->p_SPH;
                     count++;
                  }
                  if ((SPP_pressure && q->isSurface) || (SPP_pressure_self && q->isSurface && q == p))
                     for (auto I = 1; I <= 1; ++I) {
                        // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
                        auto qX = (q->X + q->interpolated_normal_SPH * I * q->radius_SPH / q->C_SML);
                        // auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
                        if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-8, p->radius_SPH})) {
                           Xij = p->X - qX;
#if defined(Morikawa2019)
                           Aij = 2. * q->mass / p->rho * Dot(Xij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(Xij, Xij);
#else
                           Aij = q->mass * 8. / std::pow(q->rho + p->rho, 2) * Dot(Xij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(Xij, Xij);
#endif
                           sum_Aij += Aij;
                           // sum_Aij_Pj += Aij * 0.;
                           sum_Aij_Pj += Aij * (-q->p_SPH);
                        }
                     }
                  //% -------------------------------------------------------------------------- */
               };

               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);

               if (count) {
#if defined(Morikawa2019)
                  p->p_SPH_ = (_WATER_DENSITY_ * p->div_U / dt + sum_Aij_Pj) / sum_Aij;
#else
                  p->p_SPH_ = (p->div_U / dt + sum_Aij_Pj) / sum_Aij;
#endif
               } else
                  p->p_SPH_ = 0;

               if (!isFinite(p->p_SPH_)) {
                  p->p_SPH_ = 0;
                  // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "pressure_SPH_ is not a finite");
               }
            }

            for (const auto &p : net->getPoints()) {
               p->DPDt_SPH = (p->p_SPH_ - p->p_SPH) / dt;
               // p->p_SPH = p->p_SPH_;
               p->RK_P.push(p->DPDt_SPH);  // 圧力
               // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
               p->p_SPH = p->p_SPH_;  // これをいれてうまく行ったことはない．
            }

            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "仮位置における圧力Pの計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
            // ! ------------------------------------------------------ */

            DebugPrint("圧力勾配∇Pを計算 & DU/Dtの計算", Magenta);
            setRigidBodyObject_Pressure(net, wall_p, RigidBodyObject, true);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               double p_rho2;
               int count = 0;
               p->contact_points_SPH = 0;
               p->gradP_SPH = {0., 0., 0.};
               net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                  if (q->getNetwork()->isRigidBody && Dot(q->normal_SPH, p->X - q->X) < 0.)
                     return;

                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  if (Between(Distance(p, q), {1E-8, p->radius_SPH})) {
                     p->contact_points_SPH++;
                     p_rho2 = q->p_SPH / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                     p->gradP_SPH += p->rho * q->mass * p_rho2 * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     count++;
                  }
                  if ((SPP_grad_P && q->isSurface) || (SPP_grad_P_self && q->isSurface && q == p)) {
                     // これは自身をも含むことがあるので，自身のために特別に作る必要はない．
                     auto qX = q->X + q->interpolated_normal_SPH * q->radius_SPH / q->C_SML;
                     // auto qX = (q->X + q->interpolated_normal_SPH * std::pow(q->volume, 1 / 3.));
                     if (canSetSPP(net, RigidBodyObject, q, qX, q->radius_SPH / q->C_SML) && Between(Distance(p, qX), {1E-8, p->radius_SPH})) {
                        // p_rho2 = 0. / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                        p_rho2 = -q->p_SPH / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                        p->gradP_SPH += p->rho * q->mass * p_rho2 * grad_w_Bspline(p->X, qX, p->radius_SPH);
                        count++;
                     }
                  }
               });
               //
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                     if (q->getNetwork()->isRigidBody && Dot(q->normal_SPH, p->X - q->X) < 0.)
                        return;

                     if (Between(Distance(p, q), {1E-8, p->radius_SPH})) {
                        p_rho2 = q->p_SPH / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                        p->gradP_SPH += p->rho * q->mass * p_rho2 * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     }
                  });
               p->DUDt_SPH -= p->gradP_SPH / p->rho;
               //
               if (!isFinite(p->DUDt_SPH))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
            }
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "圧力勾配∇Pを計算 & DU/Dtの計算" << colorOff << std::endl;

            /* -------------------------------------------------------------------------- */
            // @ ------------------------------------------------------ */
            // @                        圧力の計算                        */
            // @ ------------------------------------------------------ */
            //             DebugPrint("圧力の計算", Magenta);
            // #pragma omp parallel
            //             for (const auto &p : net->getPoints())
            // #pragma omp single nowait
            //             {
            //                p->p_SPH_ = 0.;
            //                Tddd qp;
            //                auto func = [&](const auto &q) {
            //                   qp = q->X - p->X;
            //                   p->p_SPH_ += q->volume * q->p_SPH * w_Bspline(Norm(qp), p->radius_SPH);
            //                   if (SPP_div_U_and_DrhoDt && q->isSurface)
            //                      for (auto I = 1; I <= 1; ++I) {
            //                         // auto qX = (q->X + q->interpolated_normal_SPH * I * q->radius_SPH / q->C_SML);
            //                         auto qX = (q->X + q->interpolated_normal_SPH * I * std::pow(q->volume, 1 / 3.));
            //                         qp = qX - p->X;
            //                         p->p_SPH_ += q->volume * q->p_SPH * w_Bspline(Norm(qp), p->radius_SPH);
            //                      }
            //                };

            //                net->BucketPoints.apply(p->X, p->radius_SPH, func);
            //                for (const auto &[obj, poly] : RigidBodyObject)
            //                   obj->BucketPoints.apply(p->X, p->radius_SPH, func);

            //                if (!isFinite(p->p_SPH_))
            //                   throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p_SPH_ is not a finite");
            //             }
            //             for (const auto &p : net->getPoints()) {
            //                p->p_SPH = p->p_SPH_;
            //             }
         }
         //@ -------------------------------------------------------- */
         //@                        粒子の時間発展                      */
         //@ -------------------------------------------------------- */
         DebugPrint("粒子の時間発展", Green);
#pragma omp parallel
         for (const auto &p : net->getPoints())
#pragma omp single nowait
         {
            p->RK_X.push(p->U_SPH);  // 位置
            p->setXSingle(p->RK_X.getX());
            p->RK_U.push(p->DUDt_SPH);  // 速度
            p->U_SPH = p->RK_U.getX();
/* -------------------------------------------------------------------------- */
// チェック　
#if defined(REFLECTION)
            // // これだけでは，破綻することがある
            // bool isReflect = true;
            // while (isReflect) {
            //    isReflect = false;
            //    auto closest = [&]() {
            //       auto X = p->RK_X.getX(p->U_SPH);
            //       double distance = 1E+20;
            //       networkPoint *P = nullptr;
            //       for (const auto &[obj, poly] : RigidBodyObject) {
            //          obj->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
            //             auto tmp = Distance(X, q);
            //             if (p != q && distance > tmp) {
            //                distance = tmp;
            //                P = q;
            //             }
            //          });
            //       }
            //       return P;
            //    };
            //    auto q = closest();
            //    if (q) {
            //       // auto n = q->interpolated_normal_SPH;

            //       // auto n = q->interpolated_normal_SPH;
            //       auto n = Normalize(q->normal_SPH);

            //       auto Ureflect = -2. * Dot(p->U_SPH, n) * n;
            //       auto X = p->RK_X.getX(p->U_SPH);
            //       // if (Distance(q, X) < std::sqrt(2.) * particle_spacing && std::abs(Dot(n, q->X - X)) < particle_spacing && Dot(p->U_SPH, n) < 0) {
            //       if (Distance(q, X) < particle_spacing && Dot(p->U_SPH, n) < 0) {
            //          p->DUDt_SPH += Ureflect / dt;
            //          p->RK_U.repush(p->DUDt_SPH);
            //          p->U_SPH = p->RK_U.getX();
            //          isReflect = true;
            //       }
            //    }
            // };
#endif
            /* -------------------------------------------------------------------------- */
            // p->RK_P.push(p->DPDt_SPH);  // 圧力
            // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
            // p->setDensity(_WATER_DENSITY_);
            // if (!isFinite(ToX(p)) && !isFinite(p->U_SPH))
            //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "ToX(p) is not a finite");
         }
         //
         // for (const auto &p : net->getPoints())
         //    p->setDensity(_WATER_DENSITY_);
         // for (const auto &[obj, poly] : RigidBodyObject)
         //    for (const auto &p : obj->getPoints())
         //       p->setDensity(_WATER_DENSITY_);

         std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粒子の時間発展" << colorOff << std::endl;

         // for (const auto &obj : RigidBodyObject) {
         //    vtkPolygonWriter<networkPoint *> vtp;
         //    vtp.add(ToVector(obj->getPoints()));
         //    setData(vtp, obj);
         //    std::ofstream ofs("./output/" + obj->getName() + ".vtp");
         //    vtp.write(ofs);
         //    ofs.close();
         //    std::cin.ignore();
         // }
         real_time = (*net->getPoints().begin())->RK_X.gett();
      } while (!((*net->getPoints().begin())->RK_X.finished));
      // for (const auto &[obj, poly] : RigidBodyObject) {
      //    auto points = net->getPoints();
      //    for (const auto &p : points) {
      //       auto [inside, cell, _] = poly->isInside_MethodOctree(p->X);
      //       if (inside) {
      //          std::cout << Blue << "delete" << std::endl;
      //          delete p;
      //       }
      //    }
      // }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

#endif