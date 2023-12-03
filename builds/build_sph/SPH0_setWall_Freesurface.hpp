#ifndef SPH_setWall_Freesurface_H
#define SPH_setWall_Freesurface_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

//! -------------------------------------------------------------------------- */

void setCorrectionMatrix(const auto &all_nets) {

   /*DOC_EXTRACT 0_1_1_preparation_wall_and_freesurface

   ## N.S.方程式を解く前の準備

   */

   /* -------------------------------------------------------------------------- */

   /*DOC_EXTRACT 0_1_1_set_free_surface

   ### `setCorrectionMatrix_gradient`について

   \cite{Morikawa2023}で紹介されていた，\cite{Randles1996}の勾配演算の精度を改善する行列を計算する．
   `grad_corr_M`としている．

   ```math
   \begin{align}
   {\bf M} = \left(\sum_j V_j ({\bf x}_j-{\bf x}_i) \otimes \nabla W_{ij}\right)^{-1}
   \end{align}
   ```
   WARNING: `isCaptured`を先に計算しておく必要がある．`isCaptured`が`false`の場合は，`grad_corr_M`は単位行列になる．

   */
   /* -------------------------------------------------------------------------- */
   /*                             勾配の補正行列の計算                               */
   /* -------------------------------------------------------------------------- */

   auto setCorrectionMatrix_gradient = [&all_nets](networkPoint *A) {
      // if (A->isAuxiliary)
      //    return;

      DebugPrintLevel(2, "setCorrectionMatrix start", Blue);

      // キャプチャされなかった場合
      Fill(A->grad_U, 0.);
      Fill(A->grad_corr_M, 0.);
      Fill(A->grad_corr_M_next, 0.);
      IdentityMatrix(A->inv_grad_corr_M);
      IdentityMatrix(A->inv_grad_corr_M_next);
      // キャプチャーされた粒子のみを対象にする．
      if (!A->isCaptured)
         return;

      //@ 流体内部の位置における剛体の演算修正用行列として使おう

      Tddd Xij, Uij;

      DebugPrintLevel(2, "勾配の補正行列の計算", Blue);

      double total_w = 0, total_w_next = 0;
      auto add_for_grad_correction = [&](const auto &B) {
         if (A != B) {
            total_w += w_Bspline(Norm(B->X - A->X), A->SML());
            total_w_next += w_Bspline(Norm(X_next(A) - X_next(B)), A->SML_next());
            //! DO NOT USE grad_w_Bspline_next(A, B) because it will be applied correction matrix
            A->grad_corr_M += TensorProduct(B->volume * grad_w_Bspline(A->X, B->X, A->SML()), B->X - A->X);
            A->grad_corr_M_next += TensorProduct(V_next(B) * grad_w_Bspline(X_next(A), X_next(B), A->SML_next()), X_next(B) - X_next(A));
         }
      };

      std::ranges::for_each(all_nets, [&](const auto &NET) {
         NET->BucketPoints.apply(A->X, 1.2 * A->SML(), [&](const auto &B) {
            if (canInteract(A, B))
               add_for_grad_correction(B);
         });
      });

      DebugPrintLevel(2, "勾配の補正行列の計算 done", Blue);
      A->Eigenvalues_of_M.fill(0.);

      try {
         if (std::ranges::none_of(A->grad_corr_M, [](const auto &v) { return Norm(v) == 0.; })) {
            A->inv_grad_corr_M = Inverse(A->grad_corr_M);
         }

         if (std::ranges::none_of(A->grad_corr_M_next, [](const auto &v) { return Norm(v) == 0.; })) {
            A->inv_grad_corr_M_next = Inverse(A->grad_corr_M_next);
         }

      } catch (std::exception &e) {
         std::cout << "error in setCorrectionMatrix" << std::endl;
         std::cout << "A->grad_corr_M" << std::endl;
         std::cout << A->grad_corr_M << std::endl;
         std::cout << "A->grad_corr_M_next" << std::endl;
         std::cout << A->grad_corr_M_next << std::endl;
         std::cout << "A->inv_grad_corr_M" << std::endl;
         std::cout << A->inv_grad_corr_M << std::endl;
         std::cout << "A->inv_grad_corr_M_next" << std::endl;
         std::cout << A->inv_grad_corr_M_next << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in setCorrectionMatrix");
      }

      {
         auto [lambdas, vectors] = Eigensystem(A->grad_corr_M, 1E-10, 100);
         A->Eigenvalues_of_M = lambdas;
         A->Eigenvectors_of_M = vectors;
      }
      {
         auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M, 1E-10, 100);
         A->Eigenvalues_of_M1 = lambdas;
         A->Eigenvectors_of_M1 = vectors;
      }

      {
         auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M_next, 1E-10, 100);
         A->Eigenvalues_of_M1_next = lambdas;
         A->Eigenvectors_of_M1_next = vectors;
      }
      /* -------------------------------------------------------------------------- */
      A->var_Eigenvalues_of_M = 0.;
      A->var_Eigenvalues_of_M1 = 0;

      for (const auto &v : A->Eigenvalues_of_M)
         A->var_Eigenvalues_of_M += std::pow(1. - std::abs(v), 2.);
      for (const auto &v : A->Eigenvalues_of_M1)
         A->var_Eigenvalues_of_M1 += std::pow(1. - std::abs(v), 2.);
      //
      A->var_Eigenvalues_of_M = std::sqrt(A->var_Eigenvalues_of_M) / 3.;
      A->var_Eigenvalues_of_M1 = std::sqrt(A->var_Eigenvalues_of_M1) / 3.;
      /* -------------------------------------------------------------------------- */
      A->var_Eigenvalues_of_M_next = 0.;
      A->var_Eigenvalues_of_M1_next = 0.;

      for (const auto &v : A->Eigenvalues_of_M)
         A->var_Eigenvalues_of_M_next += std::pow(1. - std::abs(v), 2.);
      for (const auto &v : A->Eigenvalues_of_M1)
         A->var_Eigenvalues_of_M1_next += std::pow(1. - std::abs(v), 2.);
      //
      A->var_Eigenvalues_of_M_next = std::sqrt(A->var_Eigenvalues_of_M_next) / 3.;
      A->var_Eigenvalues_of_M1_next = std::sqrt(A->var_Eigenvalues_of_M1_next) / 3.;
      /* -------------------------------------------------------------------------- */
      A->min_Eigenvalues_of_M = Min(A->Eigenvalues_of_M);
      A->min_Eigenvalues_of_M1 = Min(A->Eigenvalues_of_M1);
      A->max_Eigenvalues_of_M = Max(A->Eigenvalues_of_M);
      A->max_Eigenvalues_of_M1 = Max(A->Eigenvalues_of_M1);
   };

   /*DOC_EXTRACT 0_1_1_set_free_surface

   ### `setCorrectionMatrix_laplacian`について

   ラプラシアンの演算の精度を改善する行列を計算する．
   `laplacian_corr_M`としている．

   ```math
   \begin{align}
   {\bf M} = \left(\sum_j V_j ({\bf x}_j-{\bf x}_i) \otimes ({\bf x}_j-{\bf x}_i) \nabla^2 W_{ij}\right)^{-1}
   \end{align}
   ```

   WARNING: ラプラシアンの修正行列を計算するためには，先に`setCorrectionMatrix_gradient`を計算しておく必要がある．

   */
   /* -------------------------------------------------------------------------- */
   /*                             ラプラシアンの補正行列の計算                        */
   /* -------------------------------------------------------------------------- */

   auto setCorrectionMatrix_laplacian = [&all_nets](networkPoint *A) {
      if (!A->isCaptured)
         return;

      DebugPrintLevel(2, "ラプラシアンの補正行列の計算", Blue);

      Fill(A->v_reeDW, 0.);
      Fill(A->v_eeDW, 0.);
      Fill(A->v_rrDW, 0.);
      Fill(A->v_reeDW_next, 0.);
      Fill(A->v_eeDW_next, 0.);
      Fill(A->v_rrDW_next, 0.);

      Tddd DW, r, e;
      auto add_for_laplacian_correction = [&](const auto &B) {
         if (A != B) {
            if (Distance(A, B) <= A->SML()) {
               DW = B->volume * grad_w_Bspline(A->X, B->X, A->SML());  //! DO NOT USE grad_w_Bspline(A, B) because it will be applied correction matrix
               r = A->X - B->X;
               e = Normalize(r);
               A->v_reeDW += TensorProduct(TensorProduct(TensorProduct(r, e), e), DW);
               A->v_eeDW += TensorProduct(TensorProduct(e, e), DW);
               A->v_rrDW += TensorProduct(TensorProduct(r, r), DW);
            }
            if (Distance(X_next(A), X_next(B)) <= A->SML_next()) {
               DW = V_next(B) * grad_w_Bspline(X_next(A), X_next(B), A->SML_next());  //! DO NOT USE grad_w_Bspline_next(A, B) because it will be applied correction matrix
               r = X_next(A) - X_next(B);
               e = Normalize(r);
               A->v_reeDW_next += TensorProduct(TensorProduct(TensorProduct(r, e), e), DW);
               A->v_eeDW_next += TensorProduct(TensorProduct(e, e), DW);
               A->v_rrDW_next += TensorProduct(TensorProduct(r, r), DW);
            }
         }
      };

      std::ranges::for_each(all_nets, [&](const auto &NET) {
         NET->BucketPoints.apply(A->X, 1.2 * A->SML(), [&](const auto &B) {
            if (canInteract(A, B)) add_for_laplacian_correction(B);
         });
      });

      // Convert B into a 9x9 matrix
      /*
      B:T=-I
      3x3 : 3x3x3x3
      */
      auto toMatrix = [](const std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> &B) {
         std::vector<std::vector<double>> M(9, std::vector<double>(9));
         int i = 0, j = 0, k = 0, l = 0;
         for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
               for (k = 0; k < 3; ++k)
                  for (l = 0; l < 3; ++l)
                     M[i * 3 + k][j * 3 + l] = B[i][j][k][l];
         return M;
      };
      const std::vector<double> I = {-1., 0., 0., 0., -1., 0., 0., 0., -1.};
      std::vector<double> x(9);
      auto M = toMatrix(A->v_reeDW + Dot(Dot(A->v_eeDW, A->inv_grad_corr_M), A->v_rrDW));
      lapack_lu(M, x, I);
      //
      std::vector<double> x_next(9);
      M = toMatrix(A->v_reeDW_next + Dot(Dot(A->v_eeDW_next, A->inv_grad_corr_M_next), A->v_rrDW_next));
      lapack_lu(M, x_next, I);

      for (auto i = 0; i < 3; ++i)
         for (auto j = 0; j < 3; ++j) {
            A->laplacian_corr_M[i][j] = x[i * 3 + j];
            A->laplacian_corr_M_next[i][j] = x_next[i * 3 + j];
         }

      DebugPrintLevel(2, "ラプラシアンの補正行列の計算 done", Blue);
   };

   //$ -------------------------------------------------------------------------- */
   //$ -------------------------------------------------------------------------- */

   for (const auto &net : all_nets)
#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
         if (p->isCaptured)
            setCorrectionMatrix_gradient(p);

   for (const auto &net : all_nets)
#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
         if (p->isCaptured)
            setCorrectionMatrix_laplacian(p);
   //$ -------------------------------------------------------------------------- */
   //$ -------------------------------------------------------------------------- */
}

//! -------------------------------------------------------------------------- */

void captureWallParticle(const auto &net, const auto &RigidBodyObject, const auto &particle_spacing, auto &wall_p) {
   DebugPrint("captureWallParticle", Cyan);
   try {
      // 初期化
      std::vector<Network *> all_nets = {net}, rigidbody_nets;
      for (const auto &[obj, poly] : RigidBodyObject) {
         all_nets.push_back(obj);
         rigidbody_nets.push_back(obj);
      }

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
                                                              if (Q->isAuxiliary)
                                                                 return false;
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

                  // if (Norm(q->interp_normal_original) > 7500) {
                  //    q->v_to_surface_SPH = q->normal_SPH = q->particle_spacing * 0.5 * Normalize(q->interp_normal_original);
                  // } else if (Norm(q->interp_normal_original) > 900) {
                  //    q->v_to_surface_SPH = q->normal_SPH = q->particle_spacing * 1.5 * Normalize(q->interp_normal_original);
                  // } else {
                  //    q->v_to_surface_SPH = q->normal_SPH = q->particle_spacing * 2.5 * Normalize(q->interp_normal_original);
                  // }
               }
            });
         }
         // p->grad_Min_gradM = -Normalize(p->grad_Min_gradM);
      };

      DebugPrint("setCorrectionMatrix for wall particles", Cyan);

      /*DOC_EXTRACT 0_1_2_setWall

      #### `setCorrectionMatrix`で壁粒子の演算修正用行列を計算

      `setCorrectionMatrix`で壁粒子の演算修正用行列を計算する．

      */

      wall_p.clear();
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints())
            if (p->isCaptured)
               wall_p.emplace(p);

      setCorrectionMatrix(rigidbody_nets);

#pragma omp parallel
      for (const auto &p : wall_p)
#pragma omp single nowait
      {
         p->interp_normal = Normalize(Dot(p->inv_grad_corr_M, p->interp_normal_original));
         if (!isFinite(p->interp_normal))
            p->interp_normal = {0., 0., 1.};
      }
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in setWall");
   };

   DebugPrint("captureWallParticle done", Cyan);
};

/*DOC_EXTRACT 0_1_2_setWall
### 壁面粒子の抽出と値の計算
*/

void setWall(const auto &net, const auto &RigidBodyObject, const auto &particle_spacing) {
   DebugPrint("setWall", Cyan);
   try {
      // 初期化
      std::vector<Network *> all_nets = {net}, rigidbody_nets;
      for (const auto &[obj, poly] : RigidBodyObject) {
         all_nets.push_back(obj);
         rigidbody_nets.push_back(obj);
      }

      // for (const auto &p : net->getPoints()) {
      //    p->isFreeFalling = false;
      //    p->isFirstWallLayer = false;
      // }

      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3.));
            // p->interp_normal_original.fill(0.);
            // p->grad_Min_gradM.fill(0.);
            p->isFluid = false;
            p->isFreeFalling = false;
            // p->isCaptured = p->isCaptured_ = false;
            p->isSurface = false;
            p->p_SPH = 0;
            p->U_SPH = p->DUDt_SPH = p->lap_U = {0, 0, 0};
            // p->isFirstWallLayer = false;
            p->isChecked = false;
            p->intp_density = 0.;
         }

#pragma omp parallel
      for (const auto &[obj, poly] : RigidBodyObject)
#pragma omp single nowait
         for (const auto &q : obj->getPoints()) {
            if (q->isCaptured) {

               /*DOC_EXTRACT 0_1_2_setWall

               #### `isCaptured`が`true`の壁面粒子の流速の計算

               次のようにして，鏡写しのように流速を計算する．

               ```cpp
               q->U_SPH = Reflect(q->U_SPH, q->v_to_surface_SPH)
               ```

               */
               q->U_SPH.fill(0.);
               std::array<double, 3> DUDT{};
               double div_U = 0;
               q->intp_density = 0.;
               //! ちょっとした修正11/28
               if (q->isFirstWallLayer)
                  q->marker_X = q->X + q->v_to_surface_SPH;
               else
                  q->marker_X = q->X + 2. * q->v_to_surface_SPH;
               double total_vol_w = 0, r = q->SML(), vol_w;
               net->BucketPoints.apply(q->marker_X, 1.1 * r, [&](const auto &Q) {
                  if (Norm(q->marker_X - Q->X) > r)
                     return;
                  vol_w = Q->volume * w_Bspline(Norm(q->marker_X - Q->X), r);
                  q->U_SPH += Q->U_SPH * vol_w;
                  div_U += Q->div_U * vol_w;
                  DUDT += Q->DUDt_SPH * vol_w;
                  q->intp_density += Q->rho * vol_w;
                  total_vol_w += vol_w;
               });

               if (total_vol_w != 0.) {
                  DUDT /= total_vol_w;
                  div_U /= total_vol_w;
                  q->U_SPH /= total_vol_w;
                  q->marker_U = q->U_SPH;  //! そのままの値
               }
               // q->U_SPH = -q->U_SPH;
               // q->U_SPH.fill(0.);
               q->U_SPH = Reflect(q->U_SPH, q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
               // q->DUDt_SPH = Reflect(q->DUDt_SPH, q->v_to_surface_SPH);
               q->div_U = div_U;
               // q->U_SPH = Chop(Reflect(q->U_SPH, q->v_to_surface_SPH), q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
               q->RK_U.initialize(q->RK_U.get_dt(), q->RK_U.t_init, q->U_SPH, q->RK_U.steps);
               // if (q->isFirstWallLayer)
               //    q->U_SPH = Chop(q->U_SPH, q->v_to_surface_SPH);
               // q->U_SPH = Projection(Reflect(q->U_SPH, q->v_to_surface_SPH), q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
            }
         };

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
      for (const auto &obj : RigidBodyObject)
         all_nets.push_back(obj);

      auto water = net;

      // for (const auto &net : all_nets)
#pragma omp parallel
      for (const auto &A : net->getPoints())
#pragma omp single nowait
         if (A->isCaptured) {

            /*DOC_EXTRACT 0_1_3_set_free_surface

            #### 流体の法線方向の計算

            CHECKED \ref{SPH:interp_normal}{単位法線ベクトル}: $`{\bf n}_i = {\rm Normalize}\left(-\sum_j {\frac{m_j}{\rho_j} \nabla W_{ij} }\right)`$

            単位法線ベクトルは，`interp_normal`としている．

            */

            auto p = A;
            // 初期化
            p->COM_SPH.fill(0.);
            p->vec2COM.fill(0.);
            p->vec2COM_next.fill(0.);
            double total_w_vec2COM = 0;
            double total_w_vec2COM_next = 0;
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
            double w, total_volume_w = 0, total_volume_w_next = 0;
            // std::vector<networkPoint *> samples;

            net->BucketPoints.apply(p->X, 1.2 * p->SML(), [&](const auto &q) {
               // w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());
               if (Distance(p, q) < p->SML()) {
                  p->COM_SPH += q->mass * q->X;
                  p->totalMass_SPH += q->mass;

                  if (p != q) {
                     p->vec2COM += (q->X - p->X);
                     p->vec2COM_next += (X_next(q) - X_next(p));
                     total_w_vec2COM += 1;
                     total_w_vec2COM_next += 1;
                  }
               }
               double w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());
               total_volume_w += w;
               p->intp_density += q->rho * w;
               w = V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), p->SML_next());
               total_volume_w_next += w;
               p->intp_density_next += rho_next(q) * w;

               p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML());  //! DO NOT USE grad_w_Bspline(p, q) because this is original
               p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
               p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p->X, q->X, p->SML());

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
            for (const auto &obj : RigidBodyObject)
               obj->BucketPoints.apply(p->X, 1.2 * p->SML(), [&](const auto &q) {
                  if (canInteract(p, q)) {
                     // w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());
                     if (Distance(p, q) < p->SML() || Distance(X_next(p), X_next(q)) < p->SML()) {
                        w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());
                        total_volume_w += w;
                        p->intp_density += q->rho * w;
                        w = V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), p->SML());
                        total_volume_w_next += w;
                        p->intp_density_next += rho_next(q) * w;
                        p->COM_SPH += q->mass * q->X;
                        p->totalMass_SPH += q->mass;
                        if (p != q) {
                           p->vec2COM += (q->X - p->X);
                           p->vec2COM_next += (X_next(q) - X_next(p));
                           total_w_vec2COM += 1;
                           total_w_vec2COM_next += 1;
                        }
                     }
                     if (Distance(p, q) < p->particle_spacing * 1.5)
                        near_wall_particle.emplace_back(q->v_to_surface_SPH);
                     p->interp_normal_original -= a * q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML());
                     p->interp_normal_rigid -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML());
                     if (Distance(X_next(p), X_next(q)) < p->particle_spacing * 1.5)
                        near_wall_particle_next.emplace_back(q->v_to_surface_SPH);
                     p->interp_normal_original_next -= a * rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
                     p->interp_normal_rigid_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
                     p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p->X, q->X, p->SML());
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
            if (total_w_vec2COM != 0)
               p->vec2COM /= total_w_vec2COM;
            if (total_w_vec2COM_next != 0)
               p->vec2COM_next /= total_w_vec2COM_next;

            // \label{SPH:interp_normal}
            p->interp_normal = Normalize(Dot(p->inv_grad_corr_M, p->interp_normal_original));  //\label{SPH:interp_normal}
            p->interp_normal_next = Normalize(Dot(p->inv_grad_corr_M_next, p->interp_normal_original_next));
            p->interp_normal_original_choped = p->interp_normal_original;
            p->interp_normal_original_next_choped = p->interp_normal_original_next;

            // for (const auto &n : near_wall_particle) {
            //    p->interp_normal_original_choped = Chop(p->interp_normal_original_choped, n);
            //    p->interp_normal_original_next_choped = Chop(p->interp_normal_original_next_choped, n);
            //    auto tmp = Normalize(Chop(p->interp_normal, n));
            //    if (isFlat(p->interp_normal, tmp, 20 * M_PI / 180.))
            //       p->interp_normal = tmp;
            //    if (!isFinite(p->interp_normal))
            //       p->interp_normal = {0., 0., 1.};
            // }

            // for (const auto &n : near_wall_particle_next) {
            //    p->interp_normal_next = Normalize(Chop(p->interp_normal_next, n));
            //    if (!isFinite(p->interp_normal_next))
            //       p->interp_normal_next = {0., 0., 1.};
            // }
            /* -------------------------- Surface condition ---------------------------- */
            // if (A->isFluid)
            //    A->isSurface = A->var_Eigenvalues_of_M1 > 0.3;
            // 上の判定は余計に水面を判定してしまうので，ここで内部のものを水面粒子から除外する．
            // if (p->isSurface &&
            //     p->var_Eigenvalues_of_M1 < 0.3 /*誤診ではないか？*/)
            // {
            // auto surface_condition0 = [&](const auto &q) {
            //    return q->isCaptured &&
            //           Distance(p, q) < radius &&
            //           p != q &&
            //           isFlat(p->interp_normal_original, q->X - p->X, M_PI / 5);
            // }

            /*
            var_Eigenvalues_of_M1 > 0.3なんかは，まともに微分計算ができないと思われる．
            しかし，流体内部の粒子にもこの条件が当てはまる場合がある．
            */

            // 流体粒子に対してはr=2.5*ps，\theta=30度
            const auto surface_check_r_for_fluid = A->particle_spacing * 2.5;
            // 壁粒子に対してはr=2*ps，\theta=30度．理由壁は長く伸びていることがあるため，長距離の粒子を認識しないようにするため．
            const auto surface_check_r_for_wall = A->particle_spacing * 2.;
            // 以上に該当する粒子があった場合は，水面粒子として認識しない．
            A->isSurface = A->var_Eigenvalues_of_M1 > 0.3 ||
                           (A->intp_density < _WATER_DENSITY_ * 1.05 && A->var_Eigenvalues_of_M1 > 0.1 &&
                            std::ranges::none_of(all_nets, [&](const auto &net) {
                               return net->BucketPoints.any_of(A->X, A->particle_spacing * 3.,
                                                               [&](const auto &q) {
                                                                  if (q->isFluid)
                                                                     return q->isCaptured &&
                                                                            Distance(A, q) < surface_check_r_for_fluid && A != q &&
                                                                            isFlat(Dot(p->inv_grad_corr_M, p->interp_normal_original), q->X - A->X, M_PI / 6);
                                                                  else
                                                                     return q->isCaptured &&
                                                                            Distance(A, q) < surface_check_r_for_wall && A != q &&
                                                                            isFlat(Dot(p->inv_grad_corr_M, p->interp_normal_original), q->X - A->X, M_PI / 6);
                                                               });
                            }));
            if (A->intp_density > _WATER_DENSITY_ * 1.02)
               A->isSurface = false;
            if (A->var_Eigenvalues_of_M1 > 0.5)
               A->isSurface = true;

            A->isSurface_next = A->var_Eigenvalues_of_M1_next > 0.3 ||
                                (A->intp_density_next < _WATER_DENSITY_ * 1.05 && A->var_Eigenvalues_of_M1_next > 0.1 &&
                                 std::ranges::none_of(all_nets, [&](const auto &net) {
                                    return net->BucketPoints.any_of(A->X, A->particle_spacing * 3.,
                                                                    [&](const auto &q) {
                                                                       if (q->isFluid)
                                                                          return q->isCaptured &&
                                                                                 Distance(X_next(A), X_next(q)) < surface_check_r_for_fluid && A != q &&
                                                                                 isFlat(Dot(p->inv_grad_corr_M_next, p->interp_normal_original_next), X_next(q) - X_next(A), M_PI / 6);
                                                                       else
                                                                          return q->isCaptured &&
                                                                                 Distance(X_next(A), X_next(q)) < surface_check_r_for_wall && A != q &&
                                                                                 isFlat(Dot(p->inv_grad_corr_M_next, p->interp_normal_original_next), X_next(q) - X_next(A), M_PI / 6);
                                                                    });
                                 }));

            if (A->intp_density_next > _WATER_DENSITY_ * 1.02)
               A->isSurface_next = false;
            if (A->var_Eigenvalues_of_M1_next > 0.5)
               A->isSurface_next = true;

            // A->isSurface = std::ranges::none_of(all_nets, [&](const auto &net) {
            //    return net->BucketPoints.any_of(A->X, A->particle_spacing * 2.,
            //                                    [&](const auto &q) {
            //                                       return q->isCaptured &&
            //                                              Distance(A, q) < r && A != q &&
            //                                              isFlat(A->interp_normal_original, q->X - A->X, M_PI / 6);
            //                                    });
            // });

            // A->isSurface_next = std::ranges::none_of(all_nets, [&](const auto &net) {
            //    return net->BucketPoints.any_of(A->X, A->particle_spacing * 2.,
            //                                    [&](const auto &q) {
            //                                       return q->isCaptured &&
            //                                              Distance(X_next(A), X_next(q)) < r && A != q &&
            //                                              isFlat(A->interp_normal_original, X_next(q) - X_next(A), M_PI / 6);
            //                                    });
            // });

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
                  return net->BucketPoints.any_of(p->X, p->SML(), [&](const auto &q) {
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

void modify_interp_normal_original(auto &net, const auto &RigidBodyObject) {
   try {
      DebugPrint("interp_normal_originalを修正", Green);
#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
         if (p->isCaptured) {
            p->interp_normal_original.fill(0.);
            p->interp_normal_original_next.fill(0.);
            net->BucketPoints.apply(p->X, 1.2 * p->SML(), [&](const auto &q) {
               p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML());  //! DO NOT USE grad_w_Bspline(p, q) because this is original
               p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
            });
            for (const auto &obj : RigidBodyObject)
               obj->BucketPoints.apply(p->X, 1.2 * p->SML(), [&](const auto &q) {
                  if (canInteract(p, q)) {
                     p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML());
                     p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
                  }
               });
         }

      Print("interp_normal_originalを修正 done", Green);
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in setFreeSurface");
   };
};
//! -------------------------------------------------------------------------- */

//! ここでSPH_kernel_helper_functions.hppを読み込んでいるのは，上の関数でhelperを使わないようにするため
#include "SPH_kernel_helper_functions.hpp"

#endif