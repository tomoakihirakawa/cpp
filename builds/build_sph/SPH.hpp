#ifndef SPH_weightingFunctions_H
#define SPH_weightingFunctions_H

#define USE_LeapFrog

#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

/* -------------------------------------------------------------------------- */

// #define USE_SPP_Fluid
#if defined(USE_SPP_Fluid)
   #define SPP_U_coef 1.
   #define SPP_p_coef -1.
#endif

/* -------------------------------------------------------------------------- */

// #define USE_SPP_Wall
#if defined(USE_SPP_Wall)
   #define SPP_U_coef_of_Wall 1.
   #define SPP_p_coef_of_Wall -1.
   #define SPP_DUDt_coef_of_Wall 1.
   #define SPP_rho_coef_of_Wall 1.
#endif

/* -------------------------------------------------------------------------- */

#define POWER 1.

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
   int K = 0;
   for (const auto &water : nets) {
      int I = 0, J = 0;
      for (auto i = 0; i < water->BucketPoints.xsize; ++i)
         for (auto j = 0; j < water->BucketPoints.ysize; ++j)
            for (auto k = 0; k < water->BucketPoints.zsize; ++k) {
               {
                  vtkPolygonWriter<networkPoint *> vtp;
                  for (const auto &p : water->BucketPoints.buckets[i][j][k])
                     vtp.add(p);
                  std::ofstream ofs(output_directory + "each_cell_" + std::to_string(K) + "_" + std::to_string(I++) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
               {
                  vtkPolygonWriter<Tddd> vtp;
                  vtp.add(water->BucketPoints.itox(i, j, k));
                  std::ofstream ofs(output_directory + "each_cell_position_" + std::to_string(K) + "_" + std::to_string(J++) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
            }
      K++;
   }
};

/*DOC_EXTRACT SPH
[![Banner](banner.png)](banner.png)

# Smoothed Particle Hydrodynamics (SPH) ISPH EISPH

## 概要
### 前準備
1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL条件を満たすようにタイムステップ間隔 $\Delta t$を設定

### フラクショナルステップを使って初期値問題を解く

4. 水面の判定
5. $\nabla^2 {\bf u}$の計算
6. `PoissonRHS`,$b$と$\nabla^2 p^{n+1}$における$p^{n+1}$の係数の計算
7. 流速の発散から密度 ${\rho}^\ast$を計算
8. 次の時刻の圧力 $p^{n+1}$を計算
   1. 壁粒子の圧力の計算（流体粒子の現在の圧力$p^n$だけを使って近似）
   2. 流体粒子の圧力$p^{n+1}$の計算
9. $\nabla {p^{n+1}}$が計算でき， $\frac{D{\bf u}}{D t}=-\frac{1}{\rho}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}$（粘性率が一定の非圧縮性流れの加速度）を得る．
10. $\frac{D\bf u}{Dt}$を使って，流速を更新．流速を使って位置を更新

*/

std::unordered_set<networkPoint *> wall_as_fluid;
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
      std::unordered_set<Network *> net_RigidBody;
      // b# -------------------------------------------------------------------------- */
      for (const auto &[a, b, _] : RigidBodyObjectIN) {
         RigidBodyObject.push_back({a, b});
         net_RigidBody.emplace(a);
      }
      // b# -------------- バケットの生成, p->radius_SPHの範囲だけ点を取得 --------------- */
      net->setGeometricProperties();
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->setGeometricProperties();
      //
      DebugPrint("バケットの生成", Green);
      auto bucket_spacing = particle_spacing * 1.;
      net->makeBucketPoints(bucket_spacing);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->makeBucketPoints(bucket_spacing);
      DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "バケットの生成");

      // test_Bucket(net, Append(net_RigidBody, net), "./test_SPH_Bucket/", bucket_spacing);

      // void setDataOmitted(auto &vtp, const auto &Fluid) {
      //    std::unordered_map<networkPoint *, Tddd> normal_SPH;
      //    for (const auto &p : Fluid->getPoints())
      //       normal_SPH[p] = p->normal_SPH;
      //    vtp.addPointData("normal_SPH", normal_SPH);}

      DebugPrint(Green, "固定の平滑化距離の計算: C_SML * particle_spacing = ", C_SML, " * ", particle_spacing, " = ", C_SML * particle_spacing);
      for (const auto &p : net->getPoints()) {
         // p->C_SML = C_SML * std::pow(_WATER_DENSITY_ / p->rho, 3);  // C_SMLを密度の応じて変えてみる．
         // p->setDensity(_WATER_DENSITY_);
         p->isFreeFalling = false;
         p->isInsideOfBody = false;
         p->p_SPH_SPP = 0;
         p->DUDt_SPH = p->lap_U = {0, 0, 0};
         p->isCaptured = true;
         p->tmp_X = p->X;
      }
      /* ---------------------------- 流れの計算に関与する壁粒子を保存 ---------------------------- */
      wall_as_fluid.clear();
      wall_p.clear();

      // 初期化
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            p->setDensityVolume(0, 0);
            p->isFluid = false;
            p->isFreeFalling = false;
            p->isCaptured = p->isCaptured_ = false;
            p->isSurface = false;
            p->p_SPH = 0;
            p->U_SPH = p->DUDt_SPH = p->lap_U = {0, 0, 0};
            p->tmp_X = p->X;
         }
      DebugPrint("関連する壁粒子をマーク", Green);
      // capture wall particles
#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
      {
         const double captureRange = p->radius_SPH;
         const double captureRange_wall_as_fluid = p->radius_SPH;
         //
         for (const auto &[obj, poly] : RigidBodyObject) {
            obj->BucketPoints.apply(p->X, captureRange, [&](const auto &q) {
               if (Distance(p, q) < captureRange) {
                  q->isCaptured = true;
                  q->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3.));
                  if (Distance(p, q) < captureRange_wall_as_fluid)  //\label{SPH:select_wall_as_fluid}
                     q->isFluid = true;
               }
            });
         }
      };
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            if (p->isCaptured)
               wall_p.emplace(p);
            if (p->isFluid)
               wall_as_fluid.emplace(p);
         }
      DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "関連する壁粒子をマークし，保存");

      // b# ----------- CFL条件を満たすようにタイムステップ間隔dtを設定 ------------ */

      DebugPrint("近傍粒子探査が終わったら時間ステップを決めることができる", Green);
      double dt = dt_CFL(max_dt, net, RigidBodyObject);
      std::cout << "dt = " << dt << std::endl;

      // b# --------------------------- ルンゲクッタの準備 ------------------- */
      for (const auto &p : net->getPoints()) {
#if defined(USE_RungeKutta)
// p->RK_U.initialize(dt, real_time, p->U_SPH, RK_order);
// p->RK_X.initialize(dt, real_time, p->X, RK_order);
// p->RK_P.initialize(dt, real_time, p->p_SPH, RK_order);
// p->RK_rho.initialize(dt, real_time, p->rho, RK_order);
//
#elif defined(USE_LeapFrog)
         p->LPFG_X.initialize(dt, real_time, p->X, p->U_SPH);
         p->LPFG_rho.initialize(dt, real_time, p->rho, p->DrhoDt_SPH);
#endif
      }
      // b# ======================================================= */
      // b#                  ルンゲクッタを使った時間積分                 */
      // b# ======================================================= */
      /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/
      bool finished = false;
      do {
         setNormal_Surface(net, wall_p, RigidBodyObject);
         // for (const auto &p : net->getPoints())
         //    p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3));
#if defined(USE_RungeKutta)
         dt = (*net->getPoints().begin())->RK_X.getdt();
#elif defined(USE_LeapFrog)
            // dt is fixed
#endif
         const auto DT = dt;

         std::cout << "DT = " << DT << std::endl;
         for (const auto &p : wall_p) {
            p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3));
            p->tmp_U_SPH.fill(0.);
            p->U_SPH.fill(0.);
            p->DUDt_SPH.fill(0.);
            p->DUDt_SPH_.fill(0.);
            p->lap_U.fill(0.);
            p->tmp_X = p->X;
            p->rho_ = p->rho;
         }

         //@ ∇.∇UとU*を計算
         DebugPrint("∇.∇UとU*を計算");
         calcLaplacianU(net->getPoints(), Append(net_RigidBody, net), DT);
         calcLaplacianU(net->surfaceNet->getPoints(), Append(net_RigidBody, net), DT);
         calcLaplacianU(wall_p, Append(net_RigidBody, net), DT);
         // mapValueOnWall(net, wall_p, RigidBodyObject);

         //@ 圧力 p^n+1の計算
         DebugPrint("圧力 p^n+1の計算", Magenta);

         auto getX = [&](const auto p) {
            if (p->LPFG_X.is_first)
               return p->X + dt * p->U_SPH + 0.5 * dt * dt * p->DUDt_SPH_;
            else
               return p->X;
         };

         PoissonEquation(wall_p, {net}, DT, particle_spacing, getX, true);
         setPressure(wall_p);

         PoissonEquation(net->getPoints(), Append(net_RigidBody, net), DT, particle_spacing, getX);

         // debug of surfaceNet
         std::cout << "net->surfaceNet->getPoints().size() = " << net->surfaceNet->getPoints().size() << std::endl;
         PoissonEquation(net->surfaceNet->getPoints(), Append(net_RigidBody, net), DT, particle_spacing, getX);
         // PoissonEquation(wall_p, Append(net_RigidBody, net), DT);
         setPressure(net->getPoints());
         // setPressure(wall_p);
         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "圧力勾配∇Pを計算 & DU/Dtの計算");
         solvePoisson(net->getPoints(), wall_p, Append(net_RigidBody, net));

         //@ 圧力勾配 grad(P)の計算 -> DU/Dtの計算
         DebugPrint("圧力勾配∇Pを計算 & DU/Dtの計算", Magenta);
         gradP(net->getPoints(), Append(net_RigidBody, net));

         //@ 粒子の時間発展
         DebugPrint("粒子の時間発展", Green);
         updateParticles(net->getPoints(), Append(net_RigidBody, net), RigidBodyObject, particle_spacing, DT);
         net->setGeometricProperties();

         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "粒子の時間発展");

#if defined(USE_RungeKutta)
         real_time = (*net->getPoints().begin())->RK_X.get_t();
         finished = (*net->getPoints().begin())->RK_X.finished;
#elif defined(USE_LeapFrog)
         real_time = (*net->getPoints().begin())->LPFG_X.get_t();
         finished = (*net->getPoints().begin())->LPFG_X.finished;
#endif
      } while (!finished);

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
   //
   std::unordered_map<networkPoint *, Tddd> gradP_SPH;
   for (const auto &p : Fluid->getPoints())
      gradP_SPH[p] = p->gradP_SPH / p->rho;
   vtp.addPointData("gradP_SPH / rho", gradP_SPH);
   // //
   std::unordered_map<networkPoint *, Tddd> lap_U;
   for (const auto &p : Fluid->getPoints())
      lap_U[p] = p->mu_SPH / p->rho * p->lap_U;
   vtp.addPointData("nu*lapU", lap_U);
   // //
   std::unordered_map<networkPoint *, Tddd> lap_U__GRAVITY3_;
   for (const auto &p : Fluid->getPoints())
      lap_U__GRAVITY3_[p] = p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
   vtp.addPointData("nu*lapU + g", lap_U__GRAVITY3_);
   // //
   std::unordered_map<networkPoint *, Tddd> dudt;
   for (const auto &p : Fluid->getPoints())
      dudt[p] = -p->gradP_SPH / p->rho + p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
   vtp.addPointData("- grad(U)/rho + nu*lapU + g", dudt);
   //
   std::unordered_map<networkPoint *, double> div_U;
   div_U.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      div_U[p] = p->div_U;
   vtp.addPointData("div U", div_U);
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

   std::unordered_map<networkPoint *, double> len_column_value;
   for (const auto &p : Fluid->getPoints())
      len_column_value[p] = p->column_value.size();
   vtp.addPointData("len_column_value", len_column_value);

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
   std::unordered_map<networkPoint *, double> rho;
   rho.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      rho[p] = p->rho;
   vtp.addPointData("rho", rho);
   //
   std::unordered_map<networkPoint *, double> volume;
   volume.reserve(Fluid->getPoints().size());
   for (const auto &p : Fluid->getPoints())
      volume[p] = p->volume;
   vtp.addPointData("volume", volume);
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
      lap_U[p] = p->mu_SPH / p->rho * p->lap_U;
   vtp.addPointData("nu*lapU", lap_U);
   // //
   std::unordered_map<networkPoint *, Tddd> lap_U__GRAVITY3_;
   for (const auto &p : Fluid->getPoints())
      lap_U__GRAVITY3_[p] = p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
   vtp.addPointData("nu*lapU + g", lap_U__GRAVITY3_);
   // //
   std::unordered_map<networkPoint *, Tddd> dudt;
   for (const auto &p : Fluid->getPoints())
      dudt[p] = -p->gradP_SPH / p->rho + p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
   vtp.addPointData("- grad(U)/rho + nu*lapU + g", dudt);
   // //
   std::unordered_map<networkPoint *, Tddd> DUDt;
   for (const auto &p : Fluid->getPoints())
      DUDt[p] = p->DUDt_SPH;
   vtp.addPointData("DUDt", DUDt);
   // //
   std::unordered_map<networkPoint *, Tddd> bucket_index;
   for (const auto &p : Fluid->getPoints()) {
      std::array<int, 3> ijk_int = p->getNetwork()->BucketPoints.map_to_ijk.at(p);
      std::array<double, 3> ijk_double = {static_cast<double>(ijk_int[0]), static_cast<double>(ijk_int[1]), static_cast<double>(ijk_int[2])};
      bucket_index[p] = ijk_double;
   }
   vtp.addPointData("bucket_index", bucket_index);
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
   std::unordered_map<networkPoint *, double> checked_points_in_radius_of_fluid_SPH;
   for (const auto &p : Fluid->getPoints())
      checked_points_in_radius_of_fluid_SPH[p] = p->checked_points_in_radius_of_fluid_SPH;
   vtp.addPointData("checked_points_in_radius_of_fluid_SPH", checked_points_in_radius_of_fluid_SPH);
   //
   std::unordered_map<networkPoint *, double> radius;
   for (const auto &p : Fluid->getPoints())
      radius[p] = p->radius_SPH;
   vtp.addPointData("radius", radius);
   //
   std::unordered_map<networkPoint *, double> checked_points_in_radius_SPH;
   for (const auto &p : Fluid->getPoints())
      checked_points_in_radius_SPH[p] = p->checked_points_in_radius_SPH;
   vtp.addPointData("checked_points_in_radius_SPH", checked_points_in_radius_SPH);
   //
   std::unordered_map<networkPoint *, double> checked_points_SPH;
   for (const auto &p : Fluid->getPoints())
      checked_points_SPH[p] = p->checked_points_SPH;
   vtp.addPointData("checked_points_SPH", checked_points_SPH);
};

#endif