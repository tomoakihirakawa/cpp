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
         p->total_weight += (w = q->volume * std::pow(w_Bspline(Norm(p->X - q->X), p->SML()), POWER));
         p->p_SPH_ += q->p_SPH * w;
      };
      net->BucketPoints.apply(p->X, p->SML(), func);
   }
#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      p->p_SPH = p->p_SPH_ / p->total_weight;
   }
};

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
         p->total_weight += (w = q->volume * std::pow(w_Bspline(Norm(p->X - q->X), p->SML()), POWER));

         if (q == p)
            own_w = w;
         else
            p->p_SPH_SPP += q->p_SPH * w;
      };
      net->BucketPoints.apply(p->X, p->SML(), func);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->SML(), func);
      p->p_SPH_SPP /= (1 - own_w);
      // p->p_SPH_SPP /= (p->total_weight - own_w);
      // /* -------------------------------------------------------------------------- */
      // p->p_SPH_SPP = p->p_SPH;
      /* -------------------------------------------------------------------------- */
      // p->p_SPH_SPP = p->total_weight = 0;
      // double w = 0, own_w = 0;
      // auto func = [&](const auto &p_fluid) {
      //    p->total_weight += (w = p_fluid->volume * std::pow(w_Bspline(Norm(p->X - p_fluid->X), p->SML()), POWER));
      //    p->p_SPH_SPP += p_fluid->p_SPH * w;
      // };
      // net->BucketPoints.apply(p->X, p->SML(), func);
      // for (const auto &[obj, poly] : RigidBodyObject)
      //    obj->BucketPoints.apply(p->X, p->SML(), func);
      /* -------------------------------------------------------------------------- */
      // p->p_SPH_SPP = p->total_weight = 0;
      // double w = 0;
      // net->BucketPoints.apply(p->X, p->SML(), [&](const auto &p_fluid) {
      //    w = p_fluid->volume * w_Bspline(Norm(p->X - p_fluid->X), p->SML() * p->C_SML);
      //    p->p_SPH_SPP += coef * p->p_SPH * w;
      //    if (p_fluid->isSurface) {
      //       p->p_SPH_SPP += coef * p->p_SPH * w;
      //    }
      // });
      // p->p_SPH_SPP /= (w_Bspline(0, p->SML() * p->C_SML) + w_Bspline(q->SML() / q->C_SML, p->SML() * p->C_SML));
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

std::unordered_set<networkPoint *> wall_p;

void developByEISPH(Network *net,
                    const auto &RigidBodyObjectIN,
                    double &simulation_time,
                    const double C_SML,
                    const double particle_spacing,
                    const double max_dt,
                    const int RK_order) {
   try {
      TimeWatch watch;
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
      DebugPrint("バケットの生成", Green);
      const auto bucket_spacing = particle_spacing * 2.;
#pragma omp parallel
      for (const auto &obj : all_net)
#pragma omp single nowait
         obj->makeBucketPoints(bucket_spacing);
      std::cout << Green << "バケットの生成" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";
      DebugPrint(Green, "固定の平滑化距離の計算: C_SML * particle_spacing = ", C_SML, " * ", particle_spacing, " = ", C_SML * particle_spacing);

      // 初期化
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

      // CFL条件を満たすようにタイムステップ間隔dtを設定
      DebugPrint("近傍粒子探査が終わったら時間ステップを決めることができる", Green);
      double dt = dt_CFL(max_dt, net, RigidBodyObject);
      std::cout << "dt = " << dt << std::endl;

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
         setSML(Append(net_RigidBody, net));
         // {
         //    auto points = net->getPoints();
         //    for (auto p : points)
         //       if (p->isAuxiliary)
         //          delete p;

         //    {
         //       auto points = net->getPoints();
         //       for (auto p : points) {
         //          p->isAuxiliary = false;
         //       }
         //    }
         // }
         /* -------------------------------------------------------------------------- */
         // 流れの計算に関与する壁粒子を保存
         // setVectorToPolygon(net, RigidBodyObject, C_SML * particle_spacing);
         setWall(net, RigidBodyObject, particle_spacing, wall_p);
         std::cout << Green << "setWall" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";
         setFreeSurface(net, RigidBodyObject);
         std::cout << Green << "setFreeSurface" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";
         // for (const auto &p : net->getPoints())
         //    p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3));
         /* -------------------------------------------------------------------------- */

         if (false) {
            auto points = net->getPoints();
            for (const auto &p : points) {
               if (p->hasAuxiliary()) {
                  auto q = new networkPoint(net, p->X);
                  q->surfacePoint = p;
                  p->auxPoint = q;
                  //
                  q->grad_corr_M = p->grad_corr_M;
                  q->grad_corr_M_next = p->grad_corr_M_next;
                  q->inv_grad_corr_M = p->inv_grad_corr_M_next;
                  q->inv_grad_corr_M_next = p->inv_grad_corr_M_next;
                  //
                  q->grad_corr_M_rigid = p->grad_corr_M_next_rigid;
                  q->grad_corr_M_next_rigid = p->grad_corr_M_next_rigid;
                  q->inv_grad_corr_M_rigid = p->inv_grad_corr_M_next_rigid;
                  q->inv_grad_corr_M_next_rigid = p->inv_grad_corr_M_next_rigid;
                  //
                  q->isSurface = p->isSurface;
                  q->isSurface_next = p->isSurface_next;
                  q->isNeumannSurface = p->isNeumannSurface;
                  q->isAuxiliary = true;
                  //
                  q->b_vector = p->b_vector;
                  q->U_SPH = p->U_SPH;
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
                  q->isFluid = true;
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
                  q->RK_X.Xinit = p->RK_X.Xinit;
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
            net->remakeBucketPoints();
            // #pragma omp parallel
            //             for (const auto &A : net->getPoints())
            // #pragma omp single nowait
            //                if (A->isCaptured) {
            //                   setCorrectionMatrix(A, all_net);
            //                }
         }
         /* -------------------------------------------------------------------------- */

#if defined(USE_RungeKutta)
         dt = (*net->getPoints().begin())->RK_X.getdt();
         const auto DT = dt;
#elif defined(USE_LeapFrog)
         const auto DT = dt / 2.;
         delta_t = DT;
#endif

         std::cout << "DT = " << DT << std::endl;

         //@ ∇.∇UとU*を計算
         calcLaplacianU(net->getPoints(), Append(net_RigidBody, net));
         calcLaplacianU(wall_p, Append(net_RigidBody, net));
         std::cout << Green << "∇.∇UとU*を計算" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";
         setPoissonEquation(wall_p, Append(net_RigidBody, net), particle_spacing);
         setPoissonEquation(net->getPoints(), Append(net_RigidBody, net), particle_spacing);
         // debug of surfaceNet
         // std::cout << "net->surfaceNet->getPoints().size() = " << net->surfaceNet->getPoints().size() << std::endl;
         // setPoissonEquation(wall_p, Append(net_RigidBody, net), DT);
         // setPressure(net->getPoints());
         // setPressure(wall_p);
         std::cout << Green << "ポアソン方程式を立てる" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

         //@ 圧力 p^n+1の計算
         // if (simulation_time < 0.01)
         solvePoisson(net->getPoints(), wall_p);
         std::cout << Green << "圧力 p^n+1の計算" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

         // set free surface pressure using EISPH way
         // setPoissonEquation(net->surfaceNet->getPoints(), Append(net_RigidBody, net), DT, particle_spacing);
         // checkIfPoissonEqIsSatisfied(net, Append(net_RigidBody, net));

         //@ 圧力勾配 grad(P)の計算 -> DU/Dtの計算
         gradP(net->getPoints(), Append(net_RigidBody, net));
         std::cout << Green << "圧力勾配∇Pを計算 & DU/Dtの計算" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

         //@ 粒子の時間発展
         std::cout << Green << "粒子の時間発展" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";
         updateParticles(net->getPoints(), Append(net_RigidBody, net), RigidBodyObject, particle_spacing, DT);
         net->setGeometricProperties();

         std::cout << Green << "粒子の時間発展" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

#if defined(USE_RungeKutta)
         simulation_time = (*net->getPoints().begin())->RK_X.get_t();
         finished = (*net->getPoints().begin())->RK_X.finished;
#elif defined(USE_LeapFrog)
         simulation_time = (*net->getPoints().begin())->LPFG_X.get_t();
         finished = (*net->getPoints().begin())->LPFG_X.finished;
#endif

      } while (!finished);

      std::cout << Green << "１タイムステップ終了" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in developByEISPH");
   };
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

void setDataOmitted(auto &vtp, const auto &Fluid) {
   std::unordered_map<networkPoint *, double> uo_double;
   std::unordered_map<networkPoint *, Tddd> uo_3d;
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->Eigenvalues_of_M;
   vtp.addPointData("Eigenvalues_of_M", uo_3d);
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->Eigenvalues_of_M_rigid;
   vtp.addPointData("Eigenvalues_of_M_rigid", uo_3d);

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
   for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_Min_gradM;
   vtp.addPointData("grad_Min_gradM", uo_3d);
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

// void setData(auto &vtp, const auto &Fluid, const Tddd &X = {1E+50, 1E+50, 1E+50}) {

//    std::unordered_map<networkPoint *, Tddd> ViscousAndGravityForce;
//    for (const auto &p : Fluid->getPoints())
//       ViscousAndGravityForce[p] = p->ViscousAndGravityForce;
//    vtp.addPointData("ViscousAndGravityForce", ViscousAndGravityForce);

//    std::unordered_map<networkPoint *, Tddd> tmp_ViscousAndGravityForce;
//    for (const auto &p : Fluid->getPoints())
//       tmp_ViscousAndGravityForce[p] = p->tmp_ViscousAndGravityForce;
//    vtp.addPointData("tmp_ViscousAndGravityForce", tmp_ViscousAndGravityForce);

//    if (isFinite(X)) {
//       std::unordered_map<networkPoint *, double> W;
//       for (const auto &p : Fluid->getPoints())
//          W[p] = w_Bspline(Norm(X - p->X), p->SML());
//       vtp.addPointData("W", W);
//    }

//    std::unordered_map<networkPoint *, Tddd> U;
//    for (const auto &p : Fluid->getPoints())
//       U[p] = p->U_SPH;
//    vtp.addPointData("U", U);

//    std::unordered_map<networkPoint *, double> len_column_value;
//    for (const auto &p : Fluid->getPoints())
//       len_column_value[p] = p->column_value.size();
//    vtp.addPointData("len_column_value", len_column_value);

//    std::unordered_map<networkPoint *, double> C_SML;
//    for (const auto &p : Fluid->getPoints())
//       C_SML[p] = p->C_SML;
//    vtp.addPointData("C_SML", C_SML);

//    std::unordered_map<networkPoint *, double> div_U_error;
//    div_U_error.reserve(Fluid->getPoints().size());
//    for (const auto &p : Fluid->getPoints())
//       div_U_error[p] = p->div_U_error;
//    vtp.addPointData("div U error", div_U_error);
//    //
//    std::unordered_map<networkPoint *, double> rho;
//    rho.reserve(Fluid->getPoints().size());
//    for (const auto &p : Fluid->getPoints())
//       rho[p] = p->rho;
//    vtp.addPointData("rho", rho);
//    //
//    std::unordered_map<networkPoint *, double> volume;
//    volume.reserve(Fluid->getPoints().size());
//    for (const auto &p : Fluid->getPoints())
//       volume[p] = p->volume;
//    vtp.addPointData("volume", volume);
//    //
//    std::unordered_map<networkPoint *, double> pressure;
//    pressure.reserve(Fluid->getPoints().size());
//    for (const auto &p : Fluid->getPoints())
//       pressure[p] = p->p_SPH;
//    vtp.addPointData("pressure", pressure);
//    //
//    // std::unordered_map<networkPoint *, double> contactpoints;
//    // for (const auto &p : Fluid->getPoints())
//    //    contactpoints[p] = (double)p->getContactPoints().size();
//    // vtp.addPointData("contact points", contactpoints);
//    // //
//    std::unordered_map<networkPoint *, Tddd> interp_normal;
//    interp_normal.reserve(Fluid->getPoints().size());
//    for (const auto &p : Fluid->getPoints())
//       interp_normal[p] = p->interp_normal;
//    vtp.addPointData("interp_normal", interp_normal);
//    //
//    std::unordered_map<networkPoint *, Tddd> interp_normal_original;
//    interp_normal_original.reserve(Fluid->getPoints().size());
//    for (const auto &p : Fluid->getPoints())
//       interp_normal_original[p] = p->interp_normal_original;
//    vtp.addPointData("interp_normal_original", interp_normal_original);
//    // //
//    std::unordered_map<networkPoint *, Tddd> position;
//    position.reserve(Fluid->getPoints().size());
//    for (const auto &p : Fluid->getPoints())
//       position[p] = p->X;
//    vtp.addPointData("position", position);
//    // //
//    std::unordered_map<networkPoint *, Tddd> tmp_U_SPH;
//    for (const auto &p : Fluid->getPoints())
//       tmp_U_SPH[p] = p->tmp_U_SPH;
//    vtp.addPointData("tmp_U_SPH", tmp_U_SPH);
//    // //
//    std::unordered_map<networkPoint *, double> div_U;
//    for (const auto &p : Fluid->getPoints())
//       div_U[p] = p->div_U;
//    vtp.addPointData("div_U", div_U);
//    // //
//    std::unordered_map<networkPoint *, Tddd> gradP_SPH;
//    for (const auto &p : Fluid->getPoints())
//       gradP_SPH[p] = p->gradP_SPH / p->rho;
//    vtp.addPointData("gradP_SPH / rho", gradP_SPH);
//    // //
//    std::unordered_map<networkPoint *, Tddd> lap_U;
//    for (const auto &p : Fluid->getPoints())
//       lap_U[p] = p->mu_SPH / p->rho * p->lap_U;
//    vtp.addPointData("nu*lapU", lap_U);
//    // //
//    std::unordered_map<networkPoint *, Tddd> lap_U__GRAVITY3_;
//    for (const auto &p : Fluid->getPoints())
//       lap_U__GRAVITY3_[p] = p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
//    vtp.addPointData("nu*lapU + g", lap_U__GRAVITY3_);
//    // //
//    std::unordered_map<networkPoint *, Tddd> dudt;
//    for (const auto &p : Fluid->getPoints())
//       dudt[p] = -p->gradP_SPH / p->rho + p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
//    vtp.addPointData("- grad(U)/rho + nu*lapU + g", dudt);
//    // //
//    std::unordered_map<networkPoint *, Tddd> DUDt;
//    for (const auto &p : Fluid->getPoints())
//       DUDt[p] = p->DUDt_SPH;
//    vtp.addPointData("DUDt", DUDt);
//    // //
//    std::unordered_map<networkPoint *, Tddd> bucket_index;
//    for (const auto &p : Fluid->getPoints()) {
//       std::array<int, 3> ijk_int = p->getNetwork()->BucketPoints.map_to_ijk.at(p);
//       std::array<double, 3> ijk_double = {static_cast<double>(ijk_int[0]), static_cast<double>(ijk_int[1]), static_cast<double>(ijk_int[2])};
//       bucket_index[p] = ijk_double;
//    }
//    vtp.addPointData("bucket_index", bucket_index);
//    // //
//    // std::unordered_map<networkPoint *, Tddd> whereToReference;
//    // for (const auto &p : Fluid->getPoints())
//    //    whereToReference[p] = ToX(p) + 2 * p->v_to_surface_SPH - ToX(p);
//    // vtp.addPointData("where to reference", whereToReference);
//    // //
//    std::unordered_map<networkPoint *, double> isSurface;
//    isSurface.reserve(Fluid->getPoints().size());
//    for (const auto &p : Fluid->getPoints())
//       isSurface[p] = p->isSurface;
//    vtp.addPointData("isSurface", isSurface);
//    // //
//    // std::unordered_map<networkPoint *, double> isInsideOfBody;
//    // for (const auto &p : Fluid->getPoints())
//    //    isInsideOfBody[p] = p->isInsideOfBody;
//    // vtp.addPointData("isInsideOfBody", isInsideOfBody);
//    // //
//    std::unordered_map<networkPoint *, double> isCaptured;
//    for (const auto &p : Fluid->getPoints())
//       isCaptured[p] = p->isCaptured ? 1. : -1E+50;
//    vtp.addPointData("isCaptured", isCaptured);
//    // //
//    // std::unordered_map<networkPoint *, Tddd> repulsive_force_SPH;
//    // for (const auto &p : Fluid->getPoints())
//    //    repulsive_force_SPH[p] = p->repulsive_force_SPH;
//    // vtp.addPointData("repulsive_force_SPH", repulsive_force_SPH);
//    // //
//    std::unordered_map<networkPoint *, double> checked_points_in_radius_of_fluid_SPH;
//    for (const auto &p : Fluid->getPoints())
//       checked_points_in_radius_of_fluid_SPH[p] = p->checked_points_in_radius_of_fluid_SPH;
//    vtp.addPointData("checked_points_in_radius_of_fluid_SPH", checked_points_in_radius_of_fluid_SPH);
//    //
//    std::unordered_map<networkPoint *, double> radius;
//    for (const auto &p : Fluid->getPoints())
//       radius[p] = p->SML();
//    vtp.addPointData("radius", radius);
//    //
//    std::unordered_map<networkPoint *, double> checked_points_in_SML();
//    for (const auto &p : Fluid->getPoints())
//       checked_points_in_SML()[p] = p->checked_points_in_SML();
//    vtp.addPointData("checked_points_in_SML()", checked_points_in_SML());
//    //
//    std::unordered_map<networkPoint *, double> checked_points_SPH;
//    for (const auto &p : Fluid->getPoints())
//       checked_points_SPH[p] = p->checked_points_SPH;
//    vtp.addPointData("checked_points_SPH", checked_points_SPH);
// };

#endif