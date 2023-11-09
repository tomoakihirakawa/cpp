#ifndef SPH_setWall_Freesurface_H
#define SPH_setWall_Freesurface_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT 0_1_1_preparation_wall_and_freesurface

## N.S.方程式を解く前の準備

*/

/* -------------------------------------------------------------------------- */

void setCorrectionMatrix(networkPoint *A, const auto &all_nets) {

   // if (A->isAuxiliary)
   //    return;

   DebugPrintLevel(2, "setCorrectionMatrix start", Blue);

   /*DOC_EXTRACT 0_1_1_set_free_surface

   ### `setCorrectionMatrix`について

   \cite{Morikawa2023}で紹介されていた，\cite{Randles1996}の勾配演算の精度を改善する行列を計算する．
   `grad_corr_M`としている．

   ```math
   \begin{align}
   {\bf M} = \left(\sum_j V_j ({\bf x}_j-{\bf x}_i) \otimes \nabla W_{ij}\right)^{-1}
   \end{align}
   ```

   WARNING: `isCaptured`を先に計算しておく必要がある．`isCaptured`が`false`の場合は，`grad_corr_M`は単位行列になる．

   */

   // キャプチャされなかった場合
   Fill(A->grad_U, 0.);
   IdentityMatrix(A->grad_corr_M);
   IdentityMatrix(A->grad_corr_M_next);
   IdentityMatrix(A->inv_grad_corr_M);
   IdentityMatrix(A->inv_grad_corr_M_next);
   //
   // IdentityMatrix(A->inv_grad_corr_M);
   // IdentityMatrix(A->inv_grad_corr_M_next);
   //
   // IdentityMatrix(A->grad_corr_M_rigid);
   // IdentityMatrix(A->grad_corr_M_next_rigid);
   // IdentityMatrix(A->inv_grad_corr_M_rigid);
   // IdentityMatrix(A->inv_grad_corr_M_next_rigid);
   //
   // IdentityMatrix(A->Mat1);
   // IdentityMatrix(A->Mat2);
   // IdentityMatrix(A->Mat3);
   // IdentityMatrix(A->Mat_B);

   // キャプチャーされた粒子のみを対象にする．
   if (!A->isCaptured)
      return;

   // 初期化
   Fill(A->v_reeDW, 0.);
   Fill(A->v_eeDW, 0.);
   Fill(A->v_rrDW, 0.);
   Fill(A->v_reeDW_next, 0.);
   Fill(A->v_eeDW_next, 0.);
   Fill(A->v_rrDW_next, 0.);
   Fill(A->grad_U, 0.);
   Fill(A->grad_corr_M, 0.);
   Fill(A->grad_corr_M_next, 0.);
   // Fill(A->grad_corr_M_mirror, 0.);
   // Fill(A->grad_corr_M_next_mirror, 0.);

   //@ 流体内部の位置における剛体の演算修正用行列として使おう
   // Fill(A->grad_corr_M_rigid, 0.);
   // Fill(A->grad_corr_M_next_rigid, 0.);

   // A->Mat1.fill({0., 0., 0.});
   // A->Mat2.fill({0., 0., 0.});
   // A->Mat3.fill({0., 0., 0.});
   // A->Mat_B.fill({0., 0., 0.});

   Tddd Xij, Uij;

   DebugPrintLevel(2, "setCorrectionMatrix: initiation done. start add", Blue);

   auto add = [&](const auto &B) {
      auto grad_w = grad_w_Bspline(A->X, B->X, A->SML());  //! DO NOT USE grad_w_Bspline_next(A, B) because it will be applied correction matrix
      if (B->isCaptured) {
         A->grad_corr_M += B->volume * TensorProduct(grad_w, B->X - A->X);
         A->grad_corr_M_next += V_next(B) * TensorProduct(grad_w_Bspline(X_next(A), X_next(B), A->SML_next()), X_next(B) - X_next(A));  //! DO NOT USE grad_w_Bspline_next(A, B) because it will be applied correction matrix
         //! mirrror
         // A->grad_corr_M_mirror += B->volume * TensorProduct(grad_w, B->X - A->X);
         // A->grad_corr_M_next_mirror += V_next(B) * TensorProduct(grad_w_Bspline(X_next(A), X_next(B), A->SML_next()), X_next(B) - X_next(A));  //! DO NOT USE grad_w_Bspline_next(A, B) because it will be applied correction matrix
         // if (A->isFluid && !A->isSurface && B->isSurface) {
         //    auto BX = Mirror(A->X, B->X, B->X - A->X);
         //    auto grad_w = grad_w_Bspline(A->X, BX, A->SML());  //! DO NOT USE grad_w_Bspline_next(A, B) because it will be applied correction matrix
         //    A->grad_corr_M_mirror += B->volume * TensorProduct(grad_w, BX - A->X);
         //    auto BX_next = Mirror(X_next(A), X_next(B), X_next(B) - X_next(A));
         //    auto grad_w_next = grad_w_Bspline(X_next(A), BX_next, A->SML_next());  //! DO NOT USE grad_w_Bspline_next(A, B) because it will be applied correction matrix
         //    A->grad_corr_M_next_mirror += V_next(B) * TensorProduct(grad_w_next, BX_next - X_next(A));
         // }
         Uij = A->U_SPH - B->U_SPH;
         A->grad_U += B->volume * TensorProduct(-Uij, grad_w);
         //! ラプラシアンの修正
         {
            auto DW = grad_w_Bspline(A, B);
            auto e = Normalize(A->X - B->X);
            auto r = A->X - B->X;
            A->v_reeDW += B->volume * TensorProduct(TensorProduct(TensorProduct(r, e), e), DW);
            A->v_eeDW += B->volume * TensorProduct(TensorProduct(e, e), DW);
            A->v_rrDW += B->volume * TensorProduct(TensorProduct(r, r), DW);
         }
         {
            auto DW = grad_w_Bspline_next(A, B);
            auto e = Normalize(X_next(A) - X_next(B));
            auto r = X_next(A) - X_next(B);
            A->v_reeDW_next += B->volume * TensorProduct(TensorProduct(TensorProduct(r, e), e), DW);
            A->v_eeDW_next += B->volume * TensorProduct(TensorProduct(e, e), DW);
            A->v_rrDW_next += B->volume * TensorProduct(TensorProduct(r, r), DW);
         }
      }
      // if (A->getNetwork()->isRigidBody) {
      //    //@ 流体内部の位置における剛体の演算修正用行列として使おう
      //    auto n = A->v_to_surface_SPH;
      //    A->grad_corr_M_rigid += B->volume * TensorProduct(grad_w_Bspline(A->X + 2. * n, B->X, A->SML()), B->X - (A->X + 2. * n));
      //    A->grad_corr_M_next_rigid += V_next(B) * TensorProduct(grad_w_Bspline(X_next(A) + 2. * n, X_next(B), A->SML_next()), X_next(B) - (X_next(A) + 2. * n));
      // }
   };

   std::ranges::for_each(all_nets, [&](const auto &NET) {
      NET->BucketPoints.apply(A->X, 1.2 * A->SML(), [&](const auto &B) {
         if (canInteract(A, B))
            add(B);
      });
   });

   DebugPrintLevel(2, "setCorrectionMatrix: add done", Blue);

   A->Eigenvalues_of_M = {0., 0., 0.};

   const double small = 1E-10;
   if (Norm(A->grad_corr_M[0]) > small && Norm(A->grad_corr_M[1]) > small && Norm(A->grad_corr_M[2]) > small)
      A->inv_grad_corr_M = Inverse(A->grad_corr_M);

   if (Norm(A->grad_corr_M_next[0]) > small && Norm(A->grad_corr_M_next[1]) > small && Norm(A->grad_corr_M_next[2]) > small)
      A->inv_grad_corr_M_next = Inverse(A->grad_corr_M_next);
   // mirror
   // A->inv_grad_corr_M_mirror = Inverse(A->grad_corr_M_mirror);
   // A->inv_grad_corr_M_next_mirror = Inverse(A->grad_corr_M_next_mirror);
   //
   // A->inv_grad_corr_M_rigid = Inverse(A->grad_corr_M_rigid);
   // A->inv_grad_corr_M_next_rigid = Inverse(A->grad_corr_M_next_rigid);

   {
      auto [lambdas, vectors] = Eigensystem(A->grad_corr_M, 1E-10, 100.);
      A->Eigenvalues_of_M = lambdas;
      A->Eigenvectors_of_M = vectors;
   }

   // {
   //    auto [lambdas, vectors] = Eigensystem(A->grad_corr_M_rigid, 1E-10, 100.);
   //    A->Eigenvalues_of_M_rigid = lambdas;
   //    A->Eigenvectors_of_M_rigid = vectors;
   // }

   {
      auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M, 1E-10, 100.);
      A->Eigenvalues_of_M1 = lambdas;
      A->Eigenvectors_of_M1 = vectors;
   }

   {
      auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M_next, 1E-10, 100.);
      A->Eigenvalues_of_M1_next = lambdas;
      A->Eigenvectors_of_M1_next = vectors;
   }
   //
   A->var_Eigenvalues_of_M = 0;
   for (const auto &v : A->Eigenvalues_of_M)
      A->var_Eigenvalues_of_M += std::pow(1. - std::abs(v), 2.);
   for (const auto &v : A->Eigenvalues_of_M1)
      A->var_Eigenvalues_of_M1 += std::pow(1. - std::abs(v), 2.);
   //
   A->var_Eigenvalues_of_M = std::sqrt(A->var_Eigenvalues_of_M) / 3.;
   A->var_Eigenvalues_of_M1 = std::sqrt(A->var_Eigenvalues_of_M1) / 3.;
   //
   A->var_Eigenvalues_of_M1_next = 0.;
   for (const auto &v : A->Eigenvalues_of_M1_next)
      A->var_Eigenvalues_of_M1_next += std::pow(1. - std::abs(v), 2.);
   A->var_Eigenvalues_of_M1_next = std::sqrt(A->var_Eigenvalues_of_M1_next) / 3.;
   //
   A->min_Eigenvalues_of_M = Min(A->Eigenvalues_of_M);
   A->min_Eigenvalues_of_M1 = Min(A->Eigenvalues_of_M1);

   /* -------------------------------------------------------------------------- */
   DebugPrintLevel(2, "setCorrectionMatrix: Eigenvalues done. start calc laplacian_corr_M", Blue);

   // Convert B into a 9x9 matrix
   /*
   B:T=-I
   3x3 : 3x3x3x3
   */
   auto toMatrix = [](const std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> &B) {
      std::vector<std::vector<double>> M(9, std::vector<double>(9));
      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
               for (int l = 0; l < 3; ++l)
                  M[i * 3 + k][j * 3 + l] = B[i][j][k][l];
      return M;
   };
   const std::vector<double> I = {-1., 0., 0., 0., -1., 0., 0., 0., -1.};
   std::vector<double> x(9), x_next(9);
   auto M = toMatrix(A->v_reeDW + Dot(Dot(A->v_eeDW, A->inv_grad_corr_M), A->v_rrDW));
   lapack_lu(M, x, I);
   //
   auto M_next = toMatrix(A->v_reeDW_next + Dot(Dot(A->v_eeDW_next, A->inv_grad_corr_M_next), A->v_rrDW_next));
   lapack_lu(M_next, x_next, I);

   for (auto i = 0; i < 3; ++i)
      for (auto j = 0; j < 3; ++j) {
         A->laplacian_corr_M[i][j] = x[i * 3 + j];
         A->laplacian_corr_M_next[i][j] = x_next[i * 3 + j];
      }

   DebugPrintLevel(2, "setCorrectionMatrix done", Blue);
};

void setCorrectionMatrix(const auto &all_nets) {
   for (const auto &net : all_nets)
#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
         if (p->isCaptured)
            setCorrectionMatrix(p, all_nets);
}
//! -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_1_2_setWall
### 壁面粒子の抽出と値の計算
*/

void setWall(const auto &net, const auto &RigidBodyObject, const auto &particle_spacing, auto &wall_p) {
   DebugPrint("setWall", Cyan);
   try {
      // 初期化
      std::vector<Network *> all_nets = {net};
      for (const auto &[obj, poly] : RigidBodyObject)
         all_nets.push_back(obj);

      for (const auto &p : net->getPoints()) {
         p->isFreeFalling = false;
         p->isFirstWallLayer = false;
      }

      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3.));
            p->interp_normal_original.fill(0.);
            // p->grad_Min_gradM.fill(0.);
            p->isFluid = false;
            p->isFreeFalling = false;
            p->isCaptured = p->isCaptured_ = false;
            p->isSurface = false;
            p->p_SPH = 0;
            p->U_SPH = p->DUDt_SPH = p->lap_U = {0, 0, 0};
            p->isFirstWallLayer = false;
            p->isChecked = false;
            p->intp_density = 0.;
         }

      DebugPrint("isCapturedなどを計算", Cyan);

#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
      {
         const double captureRange = 3.5 * p->particle_spacing;  //\label{SPH:capture_condition_1st}
         // const double captureRange_wall_as_fluid = p->SML();
         for (const auto &[obj, poly] : RigidBodyObject) {
            obj->BucketPoints.apply(p->X, captureRange, [&](const auto &q) {
               // p->grad_Min_gradM -= p->volume * (Min(p->Eigenvalues_of_M_rigid) - Min(q->Eigenvalues_of_M_rigid)) * grad_w_Bspline(p->X, q->X, captureRange);  //! DO NOT USE grad_w_Bspline(p, q)

               if (Distance(p, q) < captureRange && !q->isChecked) {

                  q->isChecked = true;

                  /*DOC_EXTRACT 0_1_2_setWall

                  #### `interp_normal_original`の計算

                  流体粒子と同じ影響半径を使ってしまうと，流体粒子が参照できる範囲ギリギリにある壁粒子の法線方向の値が不正確になる．
                  そのため，流体粒子の影響半径よりも広い半径を使って，`q->interp_normal_original`の法線方向を計算することが，重要である．
                  少し大きい半径を`captureRange`としている．

                  */

                  for (const auto &[obj, poly] : RigidBodyObject)
                     obj->BucketPoints.apply(q->X, captureRange, [&](const auto &Q) {
                        q->interp_normal_original -= Q->rho * Q->volume * grad_w_Bspline(q->X, Q->X, captureRange);  //! DO NOT USE grad_w_Bspline(p, q)
                     });
                  const auto N = q->interp_normal_original;
                  q->interp_normal = Normalize(q->interp_normal_original);
                  if (!isFinite(q->interp_normal))
                     q->interp_normal = {0., 0., 1.};

                  /*DOC_EXTRACT 0_1_2_setWallb

                  #### `isCaptured`の決定

                  法線方向`interp_normal_original`を使って，流体粒子に近くかつ向かい合う方向にある壁粒子を抽出する．
                  計算に使用する壁粒子を決定し，使用する場合`isCaptured`を`true`にする．

                  */
                  q->isCaptured = false;
                  q->isFirstWallLayer = false;
                  q->isNeumannSurface = false;

                  auto R = captureRange;
                  q->isCaptured = net->BucketPoints.any_of(q->X, R,
                                                           [&](const auto &Q) {
                                                              bool canSeeNear = Distance(q, Q) < 1.1 * Q->SML() && q != Q && isFlat(N, Q->X - q->X, M_PI / 5);
                                                              bool veryClose = Distance(q, Q) < std::sqrt(1.5) * particle_spacing && q != Q;

                                                              bool canSeeNear_next = Distance(X_next(q), X_next(Q)) < 1.1 * Q->SML() && q != Q && isFlat(N, X_next(Q) - X_next(q), M_PI / 5);
                                                              bool veryClose_next = Distance(X_next(q), X_next(Q)) < std::sqrt(1.5) * particle_spacing && q != Q;

                                                              return canSeeNear || veryClose || canSeeNear_next || veryClose_next;
                                                           });

                  /* -------------------------------------------------------------------------- */
                  if (q->isCaptured) {
                     R = particle_spacing * 1.8;
                     q->isFirstWallLayer = net->BucketPoints.any_of(q->X, R,
                                                                    [&](const auto &Q) {
                                                                       bool canSeeNear = Distance(q, Q) < R && q != Q && isFlat(N, Q->X - q->X, M_PI / 6);
                                                                       bool veryClose = Distance(q, Q) < 1.25 * particle_spacing && q != Q;
                                                                       return canSeeNear || veryClose;
                                                                    });
                     /* -------------------------------------------------------------------------- */
                     R = particle_spacing * 1.75;
                     q->isNeumannSurface = net->BucketPoints.any_of(q->X, R,
                                                                    [&](const auto &Q) {
                                                                       bool canSeeNear = Distance(q, Q) < R && q != Q && isFlat(N, Q->X - q->X, M_PI / 5);
                                                                       bool veryClose = Distance(q, Q) < 1.25 * particle_spacing && q != Q;
                                                                       return canSeeNear || veryClose;
                                                                    });
                  }

                  if (Norm(q->interp_normal_original) > 7500) {
                     q->v_to_surface_SPH = q->normal_SPH = q->particle_spacing * 0.5 * Normalize(q->interp_normal_original);
                  } else if (Norm(q->interp_normal_original) > 900) {
                     q->v_to_surface_SPH = q->normal_SPH = q->particle_spacing * 1.5 * Normalize(q->interp_normal_original);
                  } else {
                     q->v_to_surface_SPH = q->normal_SPH = q->particle_spacing * 2.5 * Normalize(q->interp_normal_original);
                  }

                  /* -------------------------------------------------------------------------- */

                  /*DOC_EXTRACT 0_1_2_setWall

                  #### `isCaptured`が`true`の壁面粒子の流速の計算

                  次のようにして，鏡写しのように流速を計算する．

                  ```cpp
                  q->U_SPH = Reflect(q->U_SPH, q->v_to_surface_SPH)
                  ```

                  */

                  if (q->isCaptured) {
                     q->U_SPH.fill(0.);
                     Tddd b_vector, DUDt_SPH;
                     q->intp_density = 0.;
                     q->marker_X = q->X + 2 * q->v_to_surface_SPH;
                     Tddd marker_X_next = X_next(q) + 2 * q->v_to_surface_SPH;
                     double total_w = 0, r = q->SML(), w;
                     net->BucketPoints.apply(q->marker_X, 1.1 * r, [&](const auto &Q) {
                        q->U_SPH += Q->U_SPH * (w = Q->volume * w_Bspline(Norm(q->marker_X - Q->X), r));
                        b_vector += Q->b_vector * w;
                        DUDt_SPH += Q->DUDt_SPH * w;
                        q->intp_density += Q->rho * Q->volume * w_Bspline(Norm(q->marker_X - Q->X), r);
                        total_w += w;
                     });
                     if (total_w < 1E-5) {
                        q->U_SPH.fill(0.);
                        // q->intp_density = _WATER_DENSITY_;
                     } else {
                        q->U_SPH /= total_w;
                        // q->b_vector = b_vector / total_w;
                        q->DUDt_SPH /= total_w;
                        q->marker_U = q->U_SPH;  //! そのままの値
                        // q->intp_density /= total_w;
                     }
                     // q->U_SPH = -q->U_SPH;
                     // q->U_SPH.fill(0.);
                     q->U_SPH = Reflect(q->U_SPH, q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
                     // q->U_SPH = Chop(Reflect(q->U_SPH, q->v_to_surface_SPH), q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
                     q->RK_U.initialize(q->RK_U.get_dt(), q->RK_U.t_init, q->U_SPH, q->RK_U.steps);
                     // if (q->isFirstWallLayer)
                     //    q->U_SPH = Chop(q->U_SPH, q->v_to_surface_SPH);
                     // q->U_SPH = Projection(Reflect(q->U_SPH, q->v_to_surface_SPH), q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
                  }
               }
            });
         }
         // p->grad_Min_gradM = -Normalize(p->grad_Min_gradM);
      };

      DebugPrint("setCorrectionMatrix", Cyan);

      /*DOC_EXTRACT 0_1_2_setWall

      #### `setCorrectionMatrix`で壁粒子の演算修正用行列を計算

      `setCorrectionMatrix`で壁粒子の演算修正用行列を計算する．

      */

      wall_p.clear();
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints())
            if (p->isCaptured)
               wall_p.emplace(p);

      for (const auto &[obj, poly] : RigidBodyObject)
         // #pragma omp parallel
         for (const auto &p : obj->getPoints())
            // #pragma omp single nowait
            setCorrectionMatrix(p, all_nets);

   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in setWall");
   };

   DebugPrint("setWall done", Cyan);
};

/*DOC_EXTRACT 0_1_3_set_free_surface

### 流体の法線方向の計算と水面の判定

*/

void setFreeSurface(auto &net, const auto &RigidBodyObject) {
   try {
      DebugPrint("水粒子のオブジェクト外向き法線方向を計算", Green);
      // refference: A. Krimi, M. Jandaghian, and A. Shakibaeinia, Water (Switzerland), vol. 12, no. 11, pp. 1–37, 2020.

      std::vector<Network *> all_nets = {net};
      for (const auto &[obj, poly] : RigidBodyObject)
         all_nets.push_back(obj);

      auto water = net;

      //! -------------------------------------------------------------------------- */

      // for (const auto &net : all_nets)
#pragma omp parallel
      for (const auto &A : net->getPoints())
#pragma omp single nowait
         if (A->isCaptured) {
            /* -------------------------------------------------------------------------- */
            /*DOC_EXTRACT 0_1_3_set_free_surface

            #### `setCorrectionMatrix`で流体粒子の演算修正用行列を計算

            `setCorrectionMatrix`で壁粒子の演算修正用行列を計算する．

            */

            setCorrectionMatrix(A, all_nets);

            /* -------------------------------------------------------------------------- */

            /*DOC_EXTRACT 0_1_3_set_free_surface

            #### 流体の法線方向の計算

            CHECKED \ref{SPH:interp_normal}{単位法線ベクトル}: $`{\bf n}_i = {\rm Normalize}\left(-\sum_j {\frac{m_j}{\rho_j} \nabla W_{ij} }\right)`$

            単位法線ベクトルは，`interp_normal`としている．

            */

            auto p = A;
            // 初期化
            p->COM_SPH.fill(0.);
            p->totalMass_SPH = 0.;
            p->intp_density = 0.;
            p->intp_density_next = 0.;
            p->interp_normal_original.fill(0.);
            p->interp_normal_water.fill(0.);
            p->interp_normal_original_next.fill(0.);
            p->interp_normal_water_next.fill(0.);
            //
            p->intp_normal_Eigen.fill(0.);
            //
            p->interp_normal_rigid.fill(0.);
            p->interp_normal_rigid_next.fill(0.);
            //
            // p->interpolated_skewness.fill(0.);
            double w, total_weight = 0, total_weight_next = 0;
            // std::vector<networkPoint *> samples;

            net->BucketPoints.apply(p->X, p->SML(), [&](const auto &q) {
               // w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());

               double w;
               if (Distance(p, q) < p->SML()) {
                  p->COM_SPH += q->mass * q->X;
                  p->totalMass_SPH += q->mass;
               }
               w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());
               total_weight += w;
               p->intp_density += q->rho * w;
               w = V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), p->SML_next());
               total_weight_next += w;
               p->intp_density_next += rho_next(q) * w;

               p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p, q);                 // #OK to use grad_w_Bspline(p, q) because p is surface particle
               p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline_next(p, q);  // #OK to use grad_w_Bspline_next(p, q) because p is surface particle
               p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p, q);     // #OK to use grad_w_Bspline(p, q) because p is surface particle

               auto A = p;
               auto B = q;

               // 歪度
               // p->sample_SPH += 1;
               // p->mean_SPH += q->X;
               // samples.push_back(q->X);
            });
            p->interp_normal_water = p->interp_normal_original;
            p->interp_normal_water_next = p->interp_normal_original_next;
            std::vector<Tddd> near_wall_particle, near_wall_particle_next;

            const double a = 1.;
            for (const auto &[obj, poly] : RigidBodyObject)
               obj->BucketPoints.apply(p->X, 1.2 * p->SML(), [&](const auto &q) {
                  if (canInteract(p, q)) {
                     // w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());
                     if (Distance(p, q) < p->SML() || Distance(X_next(p), X_next(q)) < p->SML()) {
                        w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());
                        total_weight += w;
                        p->intp_density += q->rho * w;
                        w = V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), p->SML());
                        total_weight_next += w;
                        p->intp_density_next += rho_next(q) * w;
                        p->COM_SPH += q->mass * q->X;
                        p->totalMass_SPH += q->mass;
                     }
                     if (Distance(p, q) < p->particle_spacing * 1.5)
                        near_wall_particle.emplace_back(q->v_to_surface_SPH);
                     p->interp_normal_original -= a * q->rho * q->volume * grad_w_Bspline(p, q);  // #OK to use grad_w_Bspline(p, q) because p is surface particle
                     p->interp_normal_rigid -= q->rho * q->volume * grad_w_Bspline(p, q);         // #OK to use grad_w_Bspline(p, q) because p is surface particle
                     if (Distance(X_next(p), X_next(q)) < p->particle_spacing * 1.5)
                        near_wall_particle_next.emplace_back(q->v_to_surface_SPH);
                     p->interp_normal_original_next -= a * rho_next(q) * V_next(q) * grad_w_Bspline_next(p, q);  // #OK to use grad_w_Bspline_next(p, q) because p is surface particle
                     p->interp_normal_rigid_next -= rho_next(q) * V_next(q) * grad_w_Bspline_next(p, q);         // #OK to use grad_w_Bspline_next(p, q) because p is surface particle
                     p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p, q);         // #OK to use grad_w_Bspline(p, q) because p is surface particle
                  }
                  // 歪度
                  // p->sample_SPH += 1;
                  // p->mean_SPH += q->X;
                  // samples.push_back(q->X);
               });

            // p->mean_SPH /= p->sample_SPH;
            // p->skewness_SPH.fill(0.);
            // for (auto &s : samples)
            //    p->skewness_SPH += (s->X - p->mean_SPH) * (s->X - p->mean_SPH) * (s->X - p->mean_SPH);

            p->COM_SPH /= p->totalMass_SPH;

            // \label{SPH:interp_normal}
            p->interp_normal = Normalize(p->interp_normal_original);  //\label{SPH:interp_normal}
            p->interp_normal_original_choped = p->interp_normal_original;
            p->interp_normal_original_next_choped = p->interp_normal_original_next;

            for (const auto &n : near_wall_particle) {
               p->interp_normal_original_choped = Chop(p->interp_normal_original_choped, n);
               p->interp_normal_original_next_choped = Chop(p->interp_normal_original_next_choped, n);
               p->interp_normal = Normalize(Chop(p->interp_normal, n));
               if (!isFinite(p->interp_normal))
                  p->interp_normal = {0., 0., 1.};
            }

            p->interp_normal_next = Normalize(p->interp_normal_original_next);
            for (const auto &n : near_wall_particle_next) {
               p->interp_normal_next = Normalize(Chop(p->interp_normal_next, n));
               if (!isFinite(p->interp_normal_next))
                  p->interp_normal_next = {0., 0., 1.};
            }
            /* -------------------------- Surface condition ---------------------------- */
            // if (A->isFluid)
            //    A->isSurface = A->var_Eigenvalues_of_M1 > 0.3;
            // 上の判定は余計に水面を判定してしまうので，ここで内部のものを水面粒子から除外する．
            // if (p->isSurface &&
            //     p->var_Eigenvalues_of_M1 < 0.3 /*誤診ではないか？*/)
            // {
            const auto r = A->particle_spacing * 2.;
            // auto surface_condition0 = [&](const auto &q) {
            //    return q->isCaptured &&
            //           Distance(p, q) < radius &&
            //           p != q &&
            //           isFlat(p->interp_normal_original, q->X - p->X, M_PI / 5);
            // }
            A->isSurface = (A->var_Eigenvalues_of_M1 > 0.2) &&
                           std::ranges::none_of(all_nets,
                                                [&](const auto &net) {
                                                   return net->BucketPoints.any_of(A->X, A->particle_spacing * 2.,
                                                                                   [&](const auto &q) {
                                                                                      return q->isCaptured &&
                                                                                             Distance(A, q) < r && A != q &&
                                                                                             isFlat(A->interp_normal_original, q->X - A->X, M_PI / 6);
                                                                                   });
                                                });

            A->isSurface_next = (A->var_Eigenvalues_of_M1_next > 0.2) &&
                                std::ranges::none_of(all_nets, [&](const auto &net) {
                                   return net->BucketPoints.any_of(X_next(A), A->particle_spacing * 2.,
                                                                   [&](const auto &q) { return q->isCaptured && Distance(X_next(A), X_next(q)) < r &&
                                                                                               A != q &&
                                                                                               isFlat(A->interp_normal_original_next, q->X - A->X, M_PI / 6); });
                                });

            // p->isSurface = net->BucketPoints.none_of(p->X, radius, surface_condition0) &&
            //                std::ranges::none_of(RigidBodyObject, [&](const auto &net) {
            //                   return std::get<0>(net)->BucketPoints.any_of(p->X, radius, surface_condition0);
            //                });
            // }
            // if (p->intp_density_next < 0.5 * _WATER_DENSITY_)
            //    p->isSurface = true;

            /* -------------------------------------------------------------------------- */

            /*DOC_EXTRACT 0_1_3_set_free_surface

            #### 水面の判定

            水面の判定条件は，少し複雑である．

            */

            /* ------------------------ Neumann surface condition ----------------------- */

            const auto r_short = (p->particle_spacing) * 2.;

            auto surface = net->BucketPoints.none_of(p->X, p->SML(), [&](const auto &q) {
               return Distance(p, q) < r_short && p != q && (isFlat(p->interp_normal_water, q->X - p->X, M_PI / 5));
            });

            if (surface)
               p->isNeumannSurface = std::ranges::any_of(RigidBodyObject, [&](const auto &net) {
                  return std::get<0>(net)->BucketPoints.any_of(p->X, p->SML(), [&](const auto &q) {
                     return Distance(p, q) < r_short && p != q && (isFlat(p->interp_normal_water, q->X - p->X, M_PI / 5));
                  });
               });
            else
               p->isNeumannSurface = false;
         }

      // generate auxiliary particles
      // if (false)

      Print("水粒子のオブジェクト外向き法線方向を計算 done", Green);
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in setFreeSurface");
   };
};

#endif