// old
#ifndef SPH_H
#define SPH_H

#define USE_RungeKutta
// #define USE_LeapFrog

// ラプラシナアンの修正をするかどうか
#define USE_Laplacian_correction

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
            net->BucketPoints.apply(p->X, p->SML(), [&](const auto &q) {
               vtp.add(q);
               sum3 += Bspline3[q] = w_Bspline3(Norm(q->X - p->X), p->SML()) * p->volume;
               sum5 += Bspline5[q] = w_Bspline5(Norm(q->X - p->X), p->SML()) * p->volume;
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
            net->BucketPoints.apply(p->X, p->SML(), [&](const auto &q) {
               if (Distance(p->X, q->X) < p->SML()) {
                  vtp.add(q);
                  sum3 += Bspline3[q] = w_Bspline3(Norm(q->X - p->X), p->SML()) * p->volume;
                  sum5 += Bspline5[q] = w_Bspline5(Norm(q->X - p->X), p->SML()) * p->volume;
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

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_0_0_SPH

# Smoothed Particle Hydrodynamics (SPH) ISPH EISPH

## 概要

### 要素法と粒子法

有限要素法や境界要素法など，
要素を利用する計算手法は，
節点の接続に基づき要素を構成し（補間），
微分方程式を離散化して解く．
基本的には，要素が歪になると計算ができない．
また，上手に要素を再構成するのは大変である．

一方の粒子法は，節点間になんら決まった（要素の様な）パターンを要求せず，再構成という概念がない．
きれいに整列した粒子の方が計算精度は高いが，乱れたとしても計算はできる．

### SPH

粒子法には主に２つの種類がある．
一つは，越塚らによって提案されたMoving Particle Semi-implicit (MPS)法であり，
もう一つは，\cite{Gingold1977}と\cite{Lucy1977}によって提案されたSmoothed Particle Hydrodynamics (SPH)法である．
世界的にはSPH法がよく使われている．

SPHの研究者および産業ユーザーから成る[SPHETIC](https://www.spheric-sph.org/sph-projects-and-codes)というコミュニティがある．
それによるとSPHは，1970年代に天体物理学における非軸対称な現象を研究するために開発され，
その工学への応用は1990年代と2000年代初頭に登場した．
過去二十年で、この手法は多くの応用分野で急速に発展しており、
衝突から破壊，水面波のシミュレーション，流体-構造相互作用に至るまで多岐にわたっている．

### このプログラムの目的

このプログラムは，
ISPHとISPHを簡単化したEISPHを実装したものである．
まずは，不安要素が少ないISPHで安定した計算方法を確立し，
その後，EISPHへと移行する．

### 大まかな計算の流れ

このSPHでは，非圧縮性流体のナビエ・ストークス方程式を解く．

```math
\frac{D\bf u}{Dt} = -\frac{1}{\rho}\nabla {p} + \nu\nabla^2{\bf u} + {\bf g},\quad  \nu=\frac{\mu}{\rho}
```

#### Navier-Stokes方程式を解く前の準備

1. バケットの生成
2. 流れの計算に関与する壁粒子を保存
3. CFL条件を満たすようにタイムステップ間隔 $`\Delta t`$を設定
4. 水面の判定

#### Navier-Stokes方程式を解く

5. $`\nabla^2 {\bf u}`$の計算
6. `PoissonRHS`,$`b`$と$`\nabla^2 p^{n+1}`$における$`p^{n+1}`$の係数の計算
7. 流速の発散から密度 $`{\rho}^\ast`$を計算
8. 次の時刻の圧力 $`p^{n+1}`$を計算
   * 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
   * 流体粒子の圧力$`p^{n+1}`$の計算
9. $`\nabla {p^{n+1}}`$が計算でき， $`\frac{D{\bf u}}{D t}=-\frac{1}{\rho}\nabla {p^{n+1}} + \frac{1}{\nu}\nabla^2{\bf u} + {\bf g}`$（粘性率が一定の非圧縮性流れの加速度）を得る．
10. $`\frac{D\bf u}{Dt}`$を使って，流速を更新．流速を使って位置を更新


*/

/* -------------------------------------------------------------------------- */

std::unordered_set<networkPoint *> wall_p;

void developByEISPH(Network *net,
                    const auto &RigidBodyObjectIN,
                    double &simulation_time,
                    const double C_SML,
                    const double particle_spacing,
                    const double max_dt,
                    const int RK_order) {
   try {
      std::vector<std::tuple<Network *, Network *>> RigidBodyObject;
      std::unordered_set<Network *> net_RigidBody;
      std::vector<Network *> all_net;
      all_net.push_back(net);

      // b# -------------------------------------------------------------------------- */

      for (const auto &[a, b, _] : RigidBodyObjectIN) {
         RigidBodyObject.push_back({a, b});
         net_RigidBody.emplace(a);
         all_net.push_back(a);
      }

      // b# -------------- バケットの生成, p->SML()の範囲だけ点を取得 --------------- */
      TimeWatch watch;

      const auto bucket_spacing = particle_spacing * 2.;
#pragma omp parallel
      for (const auto &obj : all_net)
#pragma omp single nowait
         obj->makeBucketPoints(bucket_spacing);
      Print(Green, "バケットの生成", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");

      for (const auto &p : net->getPoints()) {
         // p->C_SML = C_SML * std::pow(_WATER_DENSITY_ / p->rho, 3);  // C_SMLを密度の応じて変えてみる．
         // p->setDensity(_WATER_DENSITY_);
         p->isFreeFalling = false;
         p->isInsideOfBody = false;
         p->p_SPH_SPP = 0;
         p->DUDt_SPH = p->lap_U = {0, 0, 0};
         p->isCaptured = true;
      }

      // CFL条件を満たすようにタイムステップ間隔dtを設定
      DebugPrint("近傍粒子探査が終わったら時間ステップを決めることができる", Green);
      double dt = dt_CFL(max_dt, net, RigidBodyObject);
      std::cout << "dt = " << dt << std::endl;

      deleteAuxiliaryPoints(net);
      // ルンゲクッタの準備
      for (auto &N : Append(net_RigidBody, net))
         for (const auto &p : N->getPoints()) {
#if defined(USE_RungeKutta)
            p->RK_U.initialize(dt, simulation_time, p->U_SPH, RK_order);
            p->RK_X.initialize(dt, simulation_time, p->X, RK_order);
            p->RK_P.initialize(dt, simulation_time, p->p_SPH, RK_order);
            p->RK_rho.initialize(dt, simulation_time, p->rho, RK_order);
#elif defined(USE_LeapFrog)
            p->LPFG_X.initialize(dt, simulation_time, p->X, p->U_SPH);
            p->LPFG_rho.initialize(dt, simulation_time, p->rho, p->DrhoDt_SPH);
#endif
         }

      // ルンゲクッタを使った時間積分
      /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/

      bool finished = false;

      do {
         //% 壁粒子と水面粒子の設定
         deleteAuxiliaryPoints(net);
         setSML(Append(net_RigidBody, net));
         setCorrectionMatrix(Append(net_RigidBody, net));
         setWall(net, RigidBodyObject, particle_spacing, wall_p);
         Print(Green, "setWall", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         setFreeSurface(net, RigidBodyObject);
         Print(Green, "setFreeSurface", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* -------------------------------------------------------------------------- */
         //% 粘性項の計算
         calcLaplacianU(net->getPoints(), Append(net_RigidBody, net));
         calcLaplacianU(wall_p, Append(net_RigidBody, net));
         Print(Green, "calcLaplacianU", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* ------------------------------- 次時刻の計算が可能に ------------------------------- */
         setCorrectionMatrix(Append(net_RigidBody, net));
         setWall(net, RigidBodyObject, particle_spacing, wall_p);
         setFreeSurface(net, RigidBodyObject);
         // setCorrectionMatrix(Append(net_RigidBody, net));
         // setAuxiliaryPoints(net);
         Print(Green, "recalculation", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* -------------------------------------------------------------------------- */
         //% Poisson方程式の離散化
         setPoissonEquation(wall_p, Append(net_RigidBody, net), particle_spacing);
         setPoissonEquation(net->getPoints(), Append(net_RigidBody, net), particle_spacing);
         Print(Green, "setPoissonEquation", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* -------------------------------------------------------------------------- */
         //% 圧力 p^n+1の計算
         solvePoisson(net->getPoints(), wall_p);
         Print(Green, "solvePoisson", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* -------------------------------------------------------------------------- */
         //% 圧力勾配 grad(P)の計算 -> DU/Dtの計算
         gradP(net->getPoints(), Append(net_RigidBody, net));
         Print(Green, "calculate grad(P) and increment DU/Dt", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* -------------------------------------------------------------------------- */
         //% 粒子の時間発展
         updateParticles(net->getPoints(), Append(net_RigidBody, net), RigidBodyObject, particle_spacing);
         Print(Green, "updateParticles", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         net->setGeometricProperties();
#if defined(USE_RungeKutta)
         simulation_time = (*net->getPoints().begin())->RK_X.get_t();
         finished = (*net->getPoints().begin())->RK_X.finished;
#elif defined(USE_LeapFrog)
         simulation_time = (*net->getPoints().begin())->LPFG_X.get_t();
         finished = (*net->getPoints().begin())->LPFG_X.finished;
#endif

      } while (!finished);

      Print(Green, "1タイムステップ終了", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");

   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in developByEISPH");
   };
};

/* -------------------------------------------------------------------------- */

void setDataOmitted(auto &vtp, const auto &Fluid) {
   std::unordered_map<networkPoint *, double> uo_double;
   std::unordered_map<networkPoint *, Tddd> uo_3d;
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->Eigenvalues_of_M;
   vtp.addPointData("Eigenvalues_of_M", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->Eigenvalues_of_M_rigid;
   // vtp.addPointData("Eigenvalues_of_M_rigid", uo_3d);

   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->C_SML;
   vtp.addPointData("C_SML", uo_double);

   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->var_Eigenvalues_of_M;
   vtp.addPointData("var_Eigenvalues_of_M", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->min_Eigenvalues_of_M;
   vtp.addPointData("min_Eigenvalues_of_M", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->Eigenvalues_of_M1;
   vtp.addPointData("Eigenvalues_of_M1", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->var_Eigenvalues_of_M1;
   vtp.addPointData("var_Eigenvalues_of_M1", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->min_Eigenvalues_of_M1;
   vtp.addPointData("min_Eigenvalues_of_M1", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = (p->var_Eigenvalues_of_M > 0.125);
   vtp.addPointData("isSurface_var_Eigenvalues_of_M", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = (p->var_Eigenvalues_of_M1 > 0.15);
   vtp.addPointData("isSurface_var_Eigenvalues_of_M1", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->b_vector;
   vtp.addPointData("b_vector", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->v_to_surface_SPH;
   vtp.addPointData("v_to_surface_SPH", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = Normalize(p->intp_normal_Eigen);
   vtp.addPointData("intp_normal_Eigen", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->interp_normal;
   vtp.addPointData("interp_normal", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->interp_normal_original;
   vtp.addPointData("interp_normal_original", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_Min_gradM;
   // vtp.addPointData("grad_Min_gradM", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = (double)(p->column_value.find(p) != p->column_value.end());
   vtp.addPointData("contains diagonal", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->DrhoDt_SPH;
   vtp.addPointData("DrhoDt", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->U_SPH;
   vtp.addPointData("U", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isFluid;
   vtp.addPointData("isFluid", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isCaptured ? 1 : -1E+50;
   vtp.addPointData("isCaptured", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isFirstWallLayer;
   vtp.addPointData("isFirstWallLayer", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isAuxiliary;
   vtp.addPointData("isAuxiliary", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isSurface;
   vtp.addPointData("isSurface", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isNeumannSurface;
   vtp.addPointData("isNeumannSurface", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->p_SPH;
   vtp.addPointData("pressure", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->rho;
   vtp.addPointData("density", uo_double);
   for (auto i = 0; i < 6; i++) {
      for (const auto &p : Fluid->getPoints()) {
         if (p->vector_to_polygon.size() > i)
            uo_3d[p] = p->vector_to_polygon[i];
         else
            uo_3d[p] = {0, 0, 0};
      }
      vtp.addPointData("vector_to_polygon" + std::to_string(i), uo_3d);
   }
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->vector_to_polygon.size();
   vtp.addPointData("vector_to_polygon size", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->getContactFaces().size();
   vtp.addPointData("ContactFaces size", uo_double);

   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->pressure_equation_index;
   vtp.addPointData("pressure_equation_index", uo_double);

   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->intp_density;
   vtp.addPointData("intp_density", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->volume;
   vtp.addPointData("volume", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = Projection(p->tmp_ViscousAndGravityForce, p->v_to_surface_SPH);
   vtp.addPointData("Projectioned　tmp_ViscousAndGravityForce", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = Projection(p->ViscousAndGravityForce, p->v_to_surface_SPH);
   vtp.addPointData("Projectioned　ViscousAndGravityForce", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->DUDt_modify_SPH;
   vtp.addPointData("DUDt_modify_SPH", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_double[p] = p->div_U;
   vtp.addPointData("div U", uo_double);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = (p->nearest_wall_p_next != nullptr) ? p->nearest_wall_p_next->X - p->X : Tddd{1E+50, 1E+50, 1E+50};
   vtp.addPointData("vector to nearest wall point next", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = (p->nearest_wall_p != nullptr) ? p->nearest_wall_p->X - p->X : Tddd{1E+50, 1E+50, 1E+50};
   vtp.addPointData("vector to nearest wall point", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->DUDt_SPH;
   vtp.addPointData("NS: DUDt", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->mu_SPH / p->rho * p->lap_U;
   vtp.addPointData("NS: nu*lapU", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
   vtp.addPointData("NS: nu*lapU + g", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = -p->gradP_SPH / p->rho + p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
   vtp.addPointData("NS: - gradP/rho + nu*lapU + g", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = -p->gradP_SPH / p->rho;
   vtp.addPointData("NS: - gradP/rho", uo_3d);
};

#endif