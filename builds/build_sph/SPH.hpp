#ifndef SPH_weightingFunctions_H
#define SPH_weightingFunctions_H

#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

#define REFLECTION

// #define surface_zero_pressure

/* -------------------------------------------------------------------------- */

#define USE_SPP_Fluid
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

/*DOC_EXTRACT
## ISPHとEISPHの計算過程

1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL条件を満たすようにタイムステップ間隔 $dt$を設定

4. ${{\bf u}^\ast}$と ${{\bf x}^\ast}$を計算
5. 流速の発散 ${\nabla \cdot {\bf u}^\ast}$の計算

   - Nomeritae et al. (2016)は， ${{\bf u}^\ast}$と ${{\bf x}^\ast}$を使っている
   - Morikawa, D. S., & Asai, M. (2021)， ${{\bf u}^\ast}$は使い， ${{\bf x}^\ast}$は使っていない

6. 流速の発散から密度 ${\rho}^\ast$を計算
7. 次の時刻の圧力 $p^{n+1}$を計算
   - ISPHは， $\nabla^2 {p^{n+1}}=(1-\alpha )\frac{\rho_0}{\Delta t}{\nabla \cdot {\bf u}^\ast}+\alpha \frac{\rho_0-\rho^\ast}{{\Delta t}^2}$を解く
   - EISPHは，陽的に $p^{n+1}$を計算する
8. $\nabla {p^{n+1}}$が計算でき， $\frac{D{\bf u}}{D t}=-\frac{1}{\rho_0}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}$（粘性率が一定の非圧縮性流れの加速度）を得る．
9. $\frac{D\bf u}{Dt}$を使って，流速を更新．流速を使って位置を更新

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
         p->setDensity(_WATER_DENSITY_);
         p->isFreeFalling = false;
         p->isInsideOfBody = false;
         p->p_SPH_SPP = 0;
         p->DUDt_SPH = p->lap_U = {0, 0, 0};
         p->isCaptured = true;
      }
      /* ---------------------------- 流れの計算に関与する壁粒子を保存 ---------------------------- */
      wall_as_fluid.clear();
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
#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
      {
         for (const auto &[obj, poly] : RigidBodyObject) {
            obj->BucketPoints.apply(p->X, p->radius_SPH * C, [&](const auto &q) {
               if (Distance(p, q) < p->radius_SPH * C) {
                  q->isCaptured = true;
                  q->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3.));
                  if (Distance(p, q) < p->radius_SPH / p->C_SML * 1.8) {
                     q->isFluid = true;
                  }
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

      // b# --------------- CFL条件を満たすようにタイムステップ間隔dtを設定 ----------------- */

      DebugPrint("近傍粒子探査が終わったら時間ステップを決めることができる", Green);
      double dt = dt_CFL(max_dt, net, RigidBodyObject);
      std::cout << "dt = " << dt << std::endl;

      // b# --------------------------- ルンゲクッタの準備 ------------------- */

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
#ifdef surface_zero_pressure
         for (const auto &p : net->getPoints())
            if (p->isSurface)
               p->p_SPH = 0;
#endif

         for (const auto &p : net->getPoints())
            p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3));
         dt = (*net->getPoints().begin())->RK_X.getdt();
         std::cout << "dt = " << dt << std::endl;
         mapValueOnWall(net, wall_p, RigidBodyObject);

         //@ ∇.∇UとU*を計算
         Lap_U(net->getPoints(), Append(net_RigidBody, net));
         setLap_U(net->getPoints(), dt);
         Lap_U(wall_as_fluid, Append(net_RigidBody, net));
         setLap_U(wall_as_fluid, dt);

         mapValueOnWall(net, wall_p, RigidBodyObject);

         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "粘性項の∇.∇Uを計算し，次にU*を計算");

         //% ---------------------------  仮の位置に移動して仮密度の更新 ------------------------*/
         for (const auto &p : net->getPoints()) p->setX(p->tmp_X);

         //@ 発散の計算
         div_tmpU(net->getPoints(), Append(net_RigidBody, net));
         div_tmpU(wall_as_fluid, Append(net_RigidBody, net));

         for (const auto &p : net->getPoints())
            p->rho_ = p->rho + (p->DrhoDt_SPH = -p->rho * p->div_tmpU) * dt;
         // for (const auto &p : wall_as_fluid)
         //    p->rho_ = p->rho + (p->DrhoDt_SPH = -p->rho * p->div_tmpU) * dt;

         //@ --------------------------------  元の位置に移動 --------------------------------*/
         for (const auto &p : net->getPoints()) p->setX(p->pre_X);

         //@ 圧力 p^n+1の計算
         mapValueOnWall(net, wall_p, RigidBodyObject);
         DebugPrint("仮位置における圧力Pの計算", Magenta);

         calculateTemporalPressure(net->getPoints(), Append(net_RigidBody, net), dt);
         // calculateTemporalPressure(wall_as_fluid, Append(net_RigidBody, net), dt);
         setPressure(net->getPoints());
         // setPressure(wall_as_fluid);
         mapValueOnWall(net, wall_p, RigidBodyObject);
// #define ISPH
#ifdef ISPH
         /*DOC_EXTRACT
         ISPHを使えば，水面粒子の圧力を簡単にゼロにすることができる．
         $\nabla \cdot {\bf u}^*$は流ればで満たされれば十分であり，壁面表層粒子の圧力を，壁面表層粒子上で$\nabla \cdot {\bf u}^*$となるように決める必要はない．
          */
         if (real_time > 0.0) {
            DebugPrint("activate");
            V_d b, x0;
            size_t i = 0;
            std::unordered_set<networkPoint *> points;
            points.reserve(net->getPoints().size() + wall_as_fluid.size());
            b.reserve(net->getPoints().size() + wall_as_fluid.size());
            x0.reserve(net->getPoints().size() + wall_as_fluid.size());
            const double alpha = 0.01;
            for (const auto &p : net->getPoints()) {
               if (p->isCaptured) {
                  points.emplace(p);
                  p->setIndexCSR(i++);
                  p->exclude(true);
                  b.emplace_back(p->value = (1 - alpha) * p->rho * p->div_tmpU / dt + alpha * (p->rho - p->rho_) / dt / dt);
                  x0.emplace_back(p->p_SPH);
               }
            }
            for (const auto &p : wall_as_fluid) {
               if (p->isCaptured && p->isFluid) {
                  points.emplace(p);
                  p->setIndexCSR(i++);
                  p->exclude(true);
                  b.emplace_back(p->value = (1 - alpha) * p->rho * p->div_tmpU / dt + alpha * (p->rho - p->rho_) / dt / dt);
                  x0.emplace_back(p->p_SPH);
               }
            }

            DebugPrint("Lap_P");
            Lap_P(net->getPoints(), Append(net_RigidBody, net));
            Lap_P_for_Wall(wall_as_fluid, Append(net_RigidBody, net), dt);

            for (const auto &p : net->getPoints())
               if (p->isSurface) {
                  p->column_value.clear();
                  p->increment(p, 1.);
                  p->value = 0.;
                  b[p->getIndexCSR()] = 0.;
               }

            for (const auto &p : wall_as_fluid) {
               if (p->isCaptured && p->isFluid)
                  b[p->getIndexCSR()] = p->value;
               else
                  p->p_SPH = 0;
            }

            std::cout << "gmres" << std::endl;
            gmres gm(points, b, x0, 100);
            std::cout << "gm.err : " << gm.err << std::endl;

            for (const auto &p : points) {
               p->p_SPH = gm.x[p->getIndexCSR()];
            }
         }
#endif

         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "仮位置における圧力Pの計算");

         //@ 圧力勾配 grad(P)の計算 -> DU/Dtの計算
         DebugPrint("圧力勾配∇Pを計算 & DU/Dtの計算", Magenta);

         gradP(net->getPoints(), Append(net_RigidBody, net));
         gradP(wall_as_fluid, Append(net_RigidBody, net));

         DebugPrint(Green, "Elapsed time: ", Red, watch(), "s ", Magenta, "圧力勾配∇Pを計算 & DU/Dtの計算");

         //@ -------------------------------------------------------- */
         //@                        粒子の時間発展                      */
         //@ -------------------------------------------------------- */

         DebugPrint("粒子の時間発展", Green);
         updateParticles(net->getPoints(), RigidBodyObject, particle_spacing, dt);

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