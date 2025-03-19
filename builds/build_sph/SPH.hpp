#ifndef SPH_H
#define SPH_H

#define USE_ISPH
// #define USE_ESPH

#define USE_LAPLACIAN_CORRECTION
#define USE_GRAD_CORRECTION

#define USE_SYMMETRIC_FORM_FOR_PRESSURE_GRADIENT
// #define USE_SUBTRACTIVE_FORM_FOR_PRESSURE_GRADIENT

// #define SET_AUX_AT_PARTICLE_SPACING
// #define SET_AUX_AT_MASS_CENTER
// #define USE_AUX_PREESURE

#define USE_RungeKutta
// #define USE_LeapFrog

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

void test_Bucket(const auto& water, const auto& nets, const std::string& output_directory, const auto& r) {
   std::filesystem::create_directory(output_directory);

   {
      vtkPolygonWriter<networkPoint*> vtp;
      for (const auto& net : nets)
         for (const auto& p : net->getPoints())
            vtp.add(p);
      std::ofstream ofs(output_directory + "all.vtp");
      vtp.write(ofs);
      ofs.close();
   }

   int i = 0;
   for (const auto& p : water->getPoints()) {
      {
         vtkPolygonWriter<networkPoint*> vtp;
         vtp.add(p);
         std::ofstream ofs(output_directory + "center" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
         ofs.close();
      }
      {
         std::unordered_map<networkPoint*, double> Bspline3;
         std::unordered_map<networkPoint*, double> Bspline5;
         double sum3 = 0, sum5 = 0;
         vtkPolygonWriter<networkPoint*> vtp;
         int j = 0;
         for (const auto& net : nets)
            net->BucketPoints.apply(p->X, p->SML(), [&](const auto& q) {
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
         std::unordered_map<networkPoint*, double> Bspline3;
         std::unordered_map<networkPoint*, double> Bspline5;
         double sum3 = 0, sum5 = 0;
         vtkPolygonWriter<networkPoint*> vtp;
         int j = 0;
         for (const auto& net : nets)
            net->BucketPoints.apply(p->X, p->SML(), [&](const auto& q) {
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
   for (const auto& water : nets) {
      int I = 0, J = 0;
      for (auto i = 0; i < water->BucketPoints.xsize; ++i)
         for (auto j = 0; j < water->BucketPoints.ysize; ++j)
            for (auto k = 0; k < water->BucketPoints.zsize; ++k) {
               {
                  vtkPolygonWriter<networkPoint*> vtp;
                  for (const auto& p : water->BucketPoints.buckets[i][j][k])
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

[README_ABSTRACT.md](./README_ABSTRACT.md)

[README_FOR_STUDENTS.md](./README_FOR_STUDENTS.md)

*/

/* -------------------------------------------------------------------------- */

std::unordered_set<networkPoint*> wall_p;

void developByEISPH(Network* net,
   const auto& RigidBodyObjectIN,
   double& simulation_time,
   const double CSML,
   const double particle_spacing,
   const double max_dt,
   const int RK_order) {
   try {
      std::vector<std::tuple<Network*, Network*>> RigidBodyObject;
      std::unordered_set<Network*> net_RigidBody;
      std::vector<Network*> all_net;
      all_net.push_back(net);

      // b# -------------------------------------------------------------------------- */

      for (const auto& [a, b, _] : RigidBodyObjectIN) {
         RigidBodyObject.push_back({ a, b });
         net_RigidBody.emplace(a);
         all_net.push_back(a);
      }

      // b# -------------- バケットの生成, p->SML()の範囲だけ点を取得 --------------- */
      TimeWatch watch;

      const double bucket_spacing = particle_spacing;
#pragma omp parallel
      for (const auto& obj : all_net)
#pragma omp single nowait
         obj->makeBucketPoints(0.8 * bucket_spacing);
      // obj->makeBucketPoints(bucket_spacing);
      Print(Green, "バケットの生成", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");

      for (const auto& p : net->getPoints()) {
         p->isFreeFalling = false;
         p->isInsideOfBody = false;
         p->p_SPH_SPP = 0;
         p->DUDt_SPH = p->lap_U = { 0, 0, 0 };
         p->isCaptured = true;
      }

      // CFL条件を満たすようにタイムステップ間隔dtを設定
      DebugPrint("近傍粒子探査が終わったら時間ステップを決めることができる", Green);
      double dt = dt_CFL(max_dt, net, RigidBodyObject);
      std::cout << "dt = " << dt << std::endl;

      if (dt < 1E-10)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "dt is too small");

      // deleteAuxiliaryPoints(net);
      // ルンゲクッタの準備
      for (auto& N : Append(net_RigidBody, net))
         for (const auto& p : N->getPoints()) {
            p->RK_U.initialize(dt, simulation_time, p->U_SPH, RK_order);
            p->RK_X.initialize(dt, simulation_time, p->X, RK_order);
            p->RK_P.initialize(dt, simulation_time, p->p_SPH, RK_order);
            p->RK_rho.initialize(dt, simulation_time, p->rho, RK_order);
         }

      // ルンゲクッタを使った時間積分
      /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/

      bool finished = false;

      do {
         watch();
         //% 壁粒子と水面粒子の設定
         // deleteAuxiliaryPoints(net);
         setSML(Append(net_RigidBody, net), { CSML, CSML });  // 2.7は大きすぎる 2.55は小さすぎる
         captureWallParticle(net, RigidBodyObject, particle_spacing, wall_p);
         Print(Green, "captureWallParticle", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         setCorrectionMatrix(Append(net_RigidBody, net));
         Print(Green, "setCorrectionMatrix", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         setWall(net, RigidBodyObject, particle_spacing);
         Print(Green, "setWall", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         setFreeSurface(net, net_RigidBody);
         Print(Green, "setFreeSurface", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* -------------------------------------------------------------------------- */
         //% 粘性項の計算
         calcLaplacianU(Join(net->getPoints(), wall_p), Append(net_RigidBody, net));
         Print(Green, "calcLaplacianU", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         setCorrectionMatrix(Append(net_RigidBody, net));
         Print(green, "setCorrectionMatrix", blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* ------------------------------- 次時刻の計算が可能に ------------------------ */
         Print(green, "recalculation", blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         setWall(net, RigidBodyObject, particle_spacing);
         Print(green, "setWall", blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         setFreeSurface(net, net_RigidBody);
         Print(green, "setFreeSurface", blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         /* -------------------------------------------------------------------------- */
         //% Poisson方程式の離散化
         setPoissonEquation(Join(net->getPoints(), wall_p), Append(net_RigidBody, net), particle_spacing);
         Print(Green, "setPoissonEquation", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         DebugPrint(Red, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
         /* -------------------------------------------------------------------------- */
         //% 圧力 p^n+1の計算
#if defined(USE_ISPH)
         solvePoisson(Join(net->getPoints(), wall_p));
         Print(Green, "ISPH: solvePoisson", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
#elif defined(USE_ESPH)
         if (simulation_time < 1E-10)
            solvePoisson(Join(net->getPoints(), wall_p));
         Print(Yellow, "EISPH: this does not solve Poisson equation", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
#endif
         /* -------------------------------------------------------------------------- */
         //% 圧力勾配 grad(P)の計算 -> DU/Dtの計算
         gradP(net->getPoints(), Append(net_RigidBody, net));
         Print(Green, "calculate grad(P) and increment DU/Dt", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         DebugPrint(Red, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
         //! NEW -------------------------------------------------------------------------- */
         // calculate_nabla_otimes_U_next(Join(net->getPoints(), wall_p), Append(net_RigidBody, net));
         /* -------------------------------------------------------------------------- */
         //% 粒子の時間発展
         updateParticles(net, Append(net_RigidBody, net), net_RigidBody, particle_spacing);
         Print(Green, "updateParticles", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");
         net->setGeometricProperties();

         finished = (*net->getPoints().begin())->RK_X.finished;

      } while (!finished);
      simulation_time = (*net->getPoints().begin())->RK_X.get_t();

      Print(Green, "1タイムステップ終了", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s");

   }
   catch (std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in developByEISPH");
   };
};

/* -------------------------------------------------------------------------- */

void setDataOmitted(auto& vtp, const auto& Fluid) {
   std::unordered_map<networkPoint*, double> uo_double;
   std::unordered_map<networkPoint*, Tddd> uo_3d;
   uo_double.reserve(Fluid->getPoints().size());
   uo_3d.reserve(Fluid->getPoints().size());

   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->Eigenvalues_of_M;
   vtp.addPointData("Eigenvalues_of_M", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->Eigenvalues_of_M_rigid;
   // vtp.addPointData("Eigenvalues_of_M_rigid", uo_3d);
   //
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->checked_points_in_radius_SPH;
   vtp.addPointData("checked_points_in_radius_SPH", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->checked_points_in_radius_of_fluid_SPH;
   vtp.addPointData("checked_points_in_radius_of_fluid_SPH", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->checked_points_SPH;
   vtp.addPointData("checked_points_SPH", uo_double);
   //
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->SML();
   vtp.addPointData("SML()", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->var_Eigenvalues_of_M;
   vtp.addPointData("var_Eigenvalues_of_M", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->min_Eigenvalues_of_M;
   vtp.addPointData("min_Eigenvalues_of_M", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->Eigenvalues_of_M1;
   vtp.addPointData("Eigenvalues_of_M1", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->var_Eigenvalues_of_M1;
   vtp.addPointData("var_Eigenvalues_of_M1", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = p->min_Eigenvalues_of_M1;
   // vtp.addPointData("min_Eigenvalues_of_M1", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = (p->var_Eigenvalues_of_M > 0.125);
   // vtp.addPointData("isSurface_var_Eigenvalues_of_M", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = (p->var_Eigenvalues_of_M1 > 0.15);
   // vtp.addPointData("isSurface_var_Eigenvalues_of_M1", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->b_vector;
   vtp.addPointData("b_vector", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->v_to_surface_SPH;
   vtp.addPointData("v_to_surface_SPH", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->Integral_gradP_W_SPH;
   vtp.addPointData("Integral_gradP_W_SPH", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = Normalize(p->intp_normal_Eigen);
   // vtp.addPointData("intp_normal_Eigen", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->interp_normal;
   vtp.addPointData("interp_normal", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->min_Eigenvector_of_M;
   vtp.addPointData("min_Eigenvector_of_M", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->vec2COM;
   // vtp.addPointData("vec2COM", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->vec2COM_next;
   // vtp.addPointData("vec2COM_next", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_Min_gradM;
   // vtp.addPointData("grad_Min_gradM", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = (double)(p->column_value.find(p) != p->column_value.end());
   // vtp.addPointData("contains diagonal", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = p->DrhoDt_SPH;
   // vtp.addPointData("DrhoDt", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->U_SPH;
   vtp.addPointData("U", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->isFluid;
   vtp.addPointData("isFluid", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isCaptured ? 1 : -1E+50;
   // vtp.addPointData("isCaptured", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->isFirstWallLayer;
   vtp.addPointData("isFirstWallLayer", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isAuxiliary;
   // vtp.addPointData("isAuxiliary", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->isSurface;
   vtp.addPointData("isSurface", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = p->isSurface_next;
   // vtp.addPointData("isSurface_next", uo_double);

   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->isNearSurface;
   vtp.addPointData("isNearSurface", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->isNeumannSurface;
   vtp.addPointData("isNeumannSurface", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->p_SPH;
   vtp.addPointData("pressure", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->rho;
   vtp.addPointData("density", uo_double);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_corr_M[0];
   // vtp.addPointData("grad_corr_M0", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_corr_M[1];
   // vtp.addPointData("grad_corr_M1", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_corr_M[2];
   // vtp.addPointData("grad_corr_M2", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_corr_M_next[0];
   // vtp.addPointData("grad_corr_M_next0", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_corr_M_next[1];
   // vtp.addPointData("grad_corr_M_next1", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->grad_corr_M_next[2];
   // vtp.addPointData("grad_corr_M_next2", uo_3d);

   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->inv_grad_corr_M[0];
   vtp.addPointData("inv_grad_corr_M0", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->inv_grad_corr_M[1];
   vtp.addPointData("inv_grad_corr_M1", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->inv_grad_corr_M[2];
   vtp.addPointData("inv_grad_corr_M2", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->inv_grad_corr_M_next[0];
   // vtp.addPointData("inv_grad_corr_M_next0", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->inv_grad_corr_M_next[1];
   // vtp.addPointData("inv_grad_corr_M_next1", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->inv_grad_corr_M_next[2];
   // vtp.addPointData("inv_grad_corr_M_next2", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M[0];
   // vtp.addPointData("laplacian_corr_M0", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M[1];
   // vtp.addPointData("laplacian_corr_M1", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M[2];
   // vtp.addPointData("laplacian_corr_M2", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M_next[0];
   // vtp.addPointData("laplacian_corr_M_next0", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M_next[1];
   // vtp.addPointData("laplacian_corr_M_next1", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M_next[2];
   // vtp.addPointData("laplacian_corr_M_next2", uo_3d);

   // for (auto i = 0; i < 6; i++) {
   //    for (const auto &p : Fluid->getPoints()) {
   //       if (p->vector_to_polygon.size() > i)
   //          uo_3d[p] = p->vector_to_polygon[i];
   //       else
   //          uo_3d[p] = {0, 0, 0};
   //    }
   //    vtp.addPointData("vector_to_polygon" + std::to_string(i), uo_3d);
   // }
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = p->vector_to_polygon.size();
   // vtp.addPointData("vector_to_polygon size", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = p->getContactFaces().size();
   // vtp.addPointData("ContactFaces size", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->pressure_equation_index;
   vtp.addPointData("pressure_equation_index", uo_double);
   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->intp_density;
   vtp.addPointData("intp_density", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_double[p] = p->volume;
   // vtp.addPointData("volume", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = Projection(p->tmp_ViscousAndGravityForce, p->v_to_surface_SPH);
   // vtp.addPointData("Projectioned　tmp_ViscousAndGravityForce", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = Projection(p->ViscousAndGravityForce, p->v_to_surface_SPH);
   // vtp.addPointData("Projectioned　ViscousAndGravityForce", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->DUDt_modify_SPH;
   vtp.addPointData("DUDt_modify_SPH", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->DUDt_modify_SPH_2;
   vtp.addPointData("DUDt_modify_SPH_2", uo_3d);

   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->U_XSPH;
   vtp.addPointData("U_XSPH", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = Diagonal(p->laplacian_corr_M_next);
   // vtp.addPointData("laplacian_corr_M_next", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M_next[0];
   // vtp.addPointData("laplacian_corr_M_next_x", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M_next[1];
   // vtp.addPointData("laplacian_corr_M_next_y", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M_next[2];
   // vtp.addPointData("laplacian_corr_M_next_z", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = Diagonal(p->laplacian_corr_M);
   // vtp.addPointData("laplacian_corr_M", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M[0];
   // vtp.addPointData("laplacian_corr_M_x", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M[1];
   // vtp.addPointData("laplacian_corr_M_y", uo_3d);

   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = p->laplacian_corr_M[2];
   // vtp.addPointData("laplacian_corr_M_z", uo_3d);

   for (const auto& p : Fluid->getPoints()) uo_double[p] = p->div_U;
   vtp.addPointData("div U", uo_double);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = (p->nearest_wall_p_next != nullptr) ? p->nearest_wall_p_next->X - p->X : Tddd{1E+50, 1E+50, 1E+50};
   // vtp.addPointData("vector to nearest wall point next", uo_3d);
   // for (const auto &p : Fluid->getPoints()) uo_3d[p] = (p->nearest_wall_p != nullptr) ? p->nearest_wall_p->X - p->X : Tddd{1E+50, 1E+50, 1E+50};
   // vtp.addPointData("vector to nearest wall point", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->DUDt_SPH;
   vtp.addPointData("NS: DUDt", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->mu_SPH / p->rho * p->lap_U;
   vtp.addPointData("NS: nu*lapU", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
   vtp.addPointData("NS: nu*lapU + g", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = -p->gradP_SPH / p->rho + p->mu_SPH / p->rho * p->lap_U + _GRAVITY3_;
   vtp.addPointData("NS: - gradP/rho + nu*lapU + g", uo_3d);
   for (const auto& p : Fluid->getPoints()) uo_3d[p] = -p->gradP_SPH / p->rho;
   vtp.addPointData("NS: - gradP/rho", uo_3d);
};

#endif