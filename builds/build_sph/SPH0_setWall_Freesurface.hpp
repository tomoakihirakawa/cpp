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
   勾配の演算を修正する行列は，renormalization tensorと呼ばれ，
   よく$`i`$番目の粒子に対する修正行列は$`{\bf B}_i`$と書く．
   プログラム上では\ref{SPH:grad_corr_M}{`grad_corr_M`}としている．

   ```math
   {\bf B}_i = \left(\sum_j V_j ({\bf x}_j-{\bf x}_i) \otimes \nabla W_{ij}\right)^{-1}
   ```

   WARNING: `isCaptured`を先に計算しておく必要がある．`isCaptured`が`false`の場合は，`grad_corr_M`は単位行列になる．

   */
   /* -------------------------------------------------------------------------- */
   /*                             勾配の補正行列の計算                               */
   /* -------------------------------------------------------------------------- */

   auto setCorrectionMatrix_gradient = [&all_nets](networkPoint *A) {
      DebugPrintLevel(2, "setCorrectionMatrix start", Blue);
      A->checked_points_in_radius_SPH = 0;
      A->checked_points_in_radius_of_fluid_SPH = 0;
      A->checked_points_in_radius_SPH_next = 0;
      A->checked_points_in_radius_of_fluid_SPH_next = 0;
      A->checked_points_SPH = 0;
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
            A->grad_corr_M += TensorProduct(B->volume * grad_w_Bspline(A->X, B->X, A->SML_grad()), B->X - A->X);  //\label{SPH:grad_corr_M}
            A->grad_corr_M_next += TensorProduct(V_next(B) * grad_w_Bspline(X_next(A), X_next(B), A->SML_grad_next()), X_next(B) - X_next(A));
            //
            if (Distance(A, B) <= A->SML_grad()) {
               A->checked_points_in_radius_SPH++;
               if (B->getNetwork()->isFluid || B->isFluid)
                  A->checked_points_in_radius_of_fluid_SPH++;
            }
            if (Distance(X_next(A), X_next(B)) <= A->SML_grad_next()) {
               A->checked_points_in_radius_SPH_next++;
               if (B->getNetwork()->isFluid || B->isFluid)
                  A->checked_points_in_radius_of_fluid_SPH_next++;
            }
            A->checked_points_SPH++;
         }
      };

      if (A->isFluid) {
         for (const auto &B : A->points_in_SML)
            add_for_grad_correction(B);
      } else {
         for (const auto &NET : all_nets) {
            NET->BucketPoints.apply(A->X, A->SML_grad(), [&](const auto &B) {
               if (canInteract(A, B))
                  add_for_grad_correction(B);
            });
         };
      }

      DebugPrintLevel(2, "勾配の補正行列の計算 done", Blue);
      A->Eigenvalues_of_M.fill(0.);
      A->is_grad_corr_M_singular = false;
      A->is_grad_corr_M_next_singular = false;
      try {
         // if (std::abs(Det(A->grad_corr_M))) {
         //    A->inv_grad_corr_M = Inverse(A->grad_corr_M);
         //    if (!isFinite(A->inv_grad_corr_M, 1e+2)) {
         //       IdentityMatrix(A->inv_grad_corr_M);
         //       A->is_grad_corr_M_singular = true;
         //    }
         // } else {
         //    IdentityMatrix(A->inv_grad_corr_M);
         //    A->is_grad_corr_M_singular = true;
         // }

         // if (std::abs(Det(A->grad_corr_M_next)) > 1E-13) {
         //    A->inv_grad_corr_M_next = Inverse(A->grad_corr_M_next);
         //    //! gradがOKでも，inv_gradが大きな値になる場合がある．//経験的に100を超えると，計算が不安定になる．
         //    if (!isFinite(A->inv_grad_corr_M_next, 1e+2)) {
         //       IdentityMatrix(A->inv_grad_corr_M_next);
         //       A->is_grad_corr_M_next_singular = true;
         //    }
         // } else {
         //    IdentityMatrix(A->inv_grad_corr_M_next);
         //    A->is_grad_corr_M_next_singular = true;
         // }

         auto clamp = [](std::array<std::array<double, 3>, 3> M, const double max_value) {
            double max = 0.;
            for (auto i = 0; i < 3; ++i)
               for (auto j = 0; j < 3; ++j)
                  max = std::max(max, std::abs(M[i][j]));

            if (max > max_value) {
               for (auto i = 0; i < 3; ++i)
                  for (auto j = 0; j < 3; ++j)
                     M[i][j] = M[i][j] * max_value / max;
               // IdentityMatrix(M);
            }
            return M;
         };

         auto modify = [](auto &grad_corr_M) {
            double too_small = 1E-2;
            if (Norm(grad_corr_M[0]) < too_small)
               grad_corr_M[0] = {1., 0., 0.};
            if (Norm(grad_corr_M[1]) < too_small)
               grad_corr_M[1] = {0., 1., 0.};
            if (Norm(grad_corr_M[2]) < too_small)
               grad_corr_M[2] = {0., 0., 1.};
         };

         modify(A->grad_corr_M);
         modify(A->grad_corr_M_next);

         A->inv_grad_corr_M = clamp(Inverse(A->grad_corr_M), 100.);
         A->inv_grad_corr_M_next = clamp(Inverse(A->grad_corr_M_next), 100.);

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

      A->var_Eigenvalues_of_M = 0.;
      A->var_Eigenvalues_of_M1 = 0.;
      A->var_Eigenvalues_of_M_next = 0.;
      A->var_Eigenvalues_of_M1_next = 0.;

      //! current
      if (!A->is_grad_corr_M_singular) {
         //$ grad_corr_M
         {
            auto [lambdas, vectors] = Eigensystem(A->grad_corr_M, 1E-5, 100);
            A->Eigenvalues_of_M = lambdas;
            A->Eigenvectors_of_M = vectors;

            for (const auto &v : A->Eigenvalues_of_M)
               A->var_Eigenvalues_of_M += std::pow(1. - std::abs(v), 2.);

            A->var_Eigenvalues_of_M = std::sqrt(A->var_Eigenvalues_of_M) / 3.;
         }

         //$ inv_grad_corr_M
         {
            auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M, 1E-5, 100);
            A->Eigenvalues_of_M1 = lambdas;
            A->Eigenvectors_of_M1 = vectors;

            for (const auto &v : A->Eigenvalues_of_M1)
               A->var_Eigenvalues_of_M1 += std::pow(1. - std::abs(v), 2.);
            A->var_Eigenvalues_of_M1 = std::sqrt(A->var_Eigenvalues_of_M1) / 3.;
         }
      } else {
         A->var_Eigenvalues_of_M = 10.;
         A->var_Eigenvalues_of_M1 = 10.;
      }

      //! next
      if (!A->is_grad_corr_M_next_singular) {
         //$ grad_corr_M_next
         {
            auto [lambdas, vectors] = Eigensystem(A->grad_corr_M_next, 1E-5, 100);
            A->Eigenvalues_of_M_next = lambdas;
            A->Eigenvectors_of_M_next = vectors;

            for (const auto &v : A->Eigenvalues_of_M_next)
               A->var_Eigenvalues_of_M_next += std::pow(1. - std::abs(v), 2.);
            A->var_Eigenvalues_of_M_next = std::sqrt(A->var_Eigenvalues_of_M_next) / 3.;
         }

         //$ inv_grad_corr_M_next
         {
            auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M_next, 1E-5, 100);
            A->Eigenvalues_of_M1_next = lambdas;
            A->Eigenvectors_of_M1_next = vectors;

            for (const auto &v : A->Eigenvalues_of_M1_next)
               A->var_Eigenvalues_of_M1_next += std::pow(1. - std::abs(v), 2.);
            A->var_Eigenvalues_of_M1_next = std::sqrt(A->var_Eigenvalues_of_M1_next) / 3.;
         }
      } else {
         A->var_Eigenvalues_of_M_next = 1.;
         A->var_Eigenvalues_of_M1_next = 1.;
      }

      A->var_Eigenvalues_of_M = std::clamp(A->var_Eigenvalues_of_M, 0., 1.);
      A->var_Eigenvalues_of_M_next = std::clamp(A->var_Eigenvalues_of_M_next, 0., 1.);
      /* -------------------------------------------------------------------------- */
      A->min_Eigenvalues_of_M = Min(A->Eigenvalues_of_M);
      A->min_Eigenvalues_of_M1 = Min(A->Eigenvalues_of_M1);
      A->max_Eigenvalues_of_M = Max(A->Eigenvalues_of_M);
      A->max_Eigenvalues_of_M1 = Max(A->Eigenvalues_of_M1);
      /* -------------------------------------------------------------------------- */

      if (!isFinite(A->var_Eigenvalues_of_M1) || A->var_Eigenvalues_of_M1 > 1E+2 || !isFinite(A->inv_grad_corr_M, 10.)) {
         IdentityMatrix(A->inv_grad_corr_M);
         A->is_grad_corr_M_singular = true;
      }

      if (!isFinite(A->var_Eigenvalues_of_M1_next) || A->var_Eigenvalues_of_M1_next > 1E+2 || !isFinite(A->inv_grad_corr_M_next, 10.)) {
         IdentityMatrix(A->inv_grad_corr_M_next);
         A->is_grad_corr_M_next_singular = true;
      }
   };

   for (const auto &net : all_nets) {
#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
         if (p->isCaptured)
            setCorrectionMatrix_gradient(p);
   }
   //
   /*DOC_EXTRACT 0_1_1_set_free_surface

   ### `setCorrectionMatrix_laplacian`について

   \cite{Fatehi2011}が提案した，ラプラシアンの演算の精度を改善する行列を計算する．
   プログラム上では`laplacian_corr_M`としている．
   多くの場合，流速のラプラシアンは次のように計算される．

   ```math
   \nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}\cdot{\bf x}_{ij}}
   ```

   これから，流速の勾配を引くことが１段階目の修正である．

   ```math
   \nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j - {\bf x}_{ij}\cdot{\nabla \otimes {\bf u}_{ij}})
   ```

   さらに，renomalization tensorを使って，次のように修正する．

   ```math
   \nabla^2 {\bf u}_i=\sum_{j} {\hat{\bf B}}_i:{\bf A}_{ij}({\bf u}_i - {\bf u}_j - {\bf x}_{ij}\cdot{\nabla \otimes {\bf u}_{ij}})
   \quad {\bf A}_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\otimes\nabla W_{ij}}{{\bf x}_{ij}\cdot{\bf x}_{ij}}
   ```

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

#ifdef USE_LAPLACIAN_CORRECTION
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
      IdentityMatrix(A->laplacian_corr_M);
      IdentityMatrix(A->laplacian_corr_M_next);

      // treatMatrix(A->inv_grad_corr_M, 1E-13, false);
      // treatMatrix(A->inv_grad_corr_M_next, 1E-13, false);

      Tddd DW, r, e;
      auto add_for_laplacian_correction = [&](const auto &B) {
         if (A != B) {
            if (Distance(A, B) <= A->SML_grad()) {
               DW = B->volume * grad_w_Bspline(A->X, B->X, A->SML_grad());  //! DO NOT USE grad_w_Bspline(A, B) because it will be applied correction matrix
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
         NET->BucketPoints.apply(A->X, 1.1 * A->SML_grad(), [&](const auto &B) {
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
      auto M = toMatrix(A->v_reeDW) + toMatrix(Dot(Dot(A->v_eeDW, A->inv_grad_corr_M), A->v_rrDW));
      lapack_lu(M, x, I);

      std::vector<double> x_next(9);
      auto M_next = toMatrix(A->v_reeDW_next) + toMatrix(Dot(Dot(A->v_eeDW_next, A->inv_grad_corr_M_next), A->v_rrDW_next));
      lapack_lu(M_next, x_next, I);

      for (auto i = 0; i < 3; ++i)
         for (auto j = 0; j < 3; ++j) {
            A->laplacian_corr_M[i][j] = x[i * 3 + j];
            A->laplacian_corr_M_next[i][j] = x_next[i * 3 + j];
         }

      // for (auto i = 0; i < 3; ++i)
      //    for (auto j = 0; j < 3; ++j) {
      //       if (i == j) {
      //          const double max = 1.1;
      //          auto v = x[i * 3 + j];
      //          auto v_abs = std::abs(v);
      //          A->laplacian_corr_M[i][j] = v_abs < max ? v : max * v / v_abs;
      //          v = x_next[i * 3 + j];
      //          v_abs = std::abs(v);
      //          A->laplacian_corr_M_next[i][j] = v_abs < max ? v : max * v / v_abs;
      //       }
      //    }

      //! 諦める
      /* -------------------------------------------------------------------------- */

      // treatMatrix(A->laplacian_corr_M_next, 1E-13, false);
      // treatMatrix(A->laplacian_corr_M, 1E-13, false);

      // if (!isFinite(A->laplacian_corr_M, 1e+13))
      // if (A->isSurface)
      //    IdentityMatrix(A->laplacian_corr_M);

      // if (!isFinite(A->laplacian_corr_M_next, 1e+13))
      // if (A->isSurface)
      //    IdentityMatrix(A->laplacian_corr_M_next);

      //
      if (!isFinite(A->laplacian_corr_M) || A->is_grad_corr_M_singular)
         IdentityMatrix(A->laplacian_corr_M);

      if (!isFinite(A->laplacian_corr_M_next) || A->is_grad_corr_M_next_singular)
         IdentityMatrix(A->laplacian_corr_M_next);

      /* -------------------------------------------------------------------------- */

      DebugPrintLevel(2, "ラプラシアンの補正行列の計算 done", Blue);
   };

   for (const auto &net : all_nets)
   #pragma omp parallel
      for (const auto &p : net->getPoints())
   #pragma omp single nowait
         if (p->isCaptured)
            setCorrectionMatrix_laplacian(p);
#endif
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

      const double captureRange = 3.3 * particle_spacing;  //\label{SPH:capture_condition_1st}

#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
      {
         // const double captureRange_wall_as_fluid = p->SML_grad();
         for (const auto &[obj, poly] : RigidBodyObject) {
            obj->BucketPoints.apply(p->X, captureRange, [&](const auto &q) {
               // p->grad_Min_gradM -= p->volume * (Min(p->Eigenvalues_of_M_rigid) - Min(q->Eigenvalues_of_M_rigid)) * grad_w_Bspline(p->X, q->X, captureRange);  //! DO NOT USE grad_w_Bspline(p, q)
               if (!q->isChecked && Distance(p, q) < captureRange) {
                  q->isChecked = true;

                  /*DOC_EXTRACT 0_1_2_setWall

                  #### `interp_normal_original`の計算

                  流体粒子と同じ影響半径を使ってしまうと，流体粒子が参照できる範囲ギリギリにある壁粒子の法線方向の値が不正確になる．
                  そのため，流体粒子の影響半径よりも広い半径を使って，`q->interp_normal_original`の法線方向を計算することが，重要である．
                  少し大きい半径を`captureRange`としている．

                  */

                  Tddd N = {0., 0., 0.};
                  for (const auto &[obj, poly] : RigidBodyObject)
                     obj->BucketPoints.apply(q->X, captureRange, [&](const auto &Q) {
                        N -= Q->rho * Q->volume * grad_w_Bspline(q->X, Q->X, captureRange);  //! DO NOT USE grad_w_Bspline(p, q)
                     });
                  q->interp_normal_original = N;

                  /*DOC_EXTRACT 0_1_2_setWallb

                  #### `isCaptured`の決定

                  法線方向`interp_normal_original`を使って，流体粒子に近くかつ向かい合う方向にある壁粒子を抽出する．
                  計算に使用する壁粒子を決定し，使用する場合`isCaptured`を`true`にする．

                  */
                  q->isCaptured = false;
                  q->isFirstWallLayer = false;
                  q->isNeumannSurface = false;

                  auto R = captureRange;
                  q->isCaptured = net->BucketPoints.any_of(q->X, R, [&](const auto &Q) {
                     if (Q->isAuxiliary)
                        return false;
                     bool canSeeNear = Distance(q, Q) < Q->SML() && q != Q && isFlat(N, Q->X - q->X, M_PI / 5);
                     bool veryClose = Distance(q, Q) < std::sqrt(1.5) * particle_spacing && q != Q;

                     bool canSeeNear_next = Distance(X_next(q), X_next(Q)) < Q->SML() && q != Q && isFlat(N, X_next(Q) - X_next(q), M_PI / 5);
                     bool veryClose_next = Distance(X_next(q), X_next(Q)) < std::sqrt(1.5) * particle_spacing && q != Q;

                     return canSeeNear || veryClose || canSeeNear_next || veryClose_next;
                  });

                  if (q->isCaptured) {
                     const double c = 1.5;
                     R = particle_spacing * 1.8;
                     q->isFirstWallLayer = net->BucketPoints.any_of(q->X, R, [&](const auto &Q) {
                        bool canSeeNear = Distance(q, Q) < R && q != Q && isFlat(N, Q->X - q->X, M_PI / 6);
                        bool veryClose = Distance(q, Q) < c * particle_spacing && q != Q;
                        return canSeeNear || veryClose;
                     });
                     /* -------------------------------------------------------------------------- */
                     R = particle_spacing * 1.75;
                     q->isNeumannSurface = net->BucketPoints.any_of(q->X, R, [&](const auto &Q) {
                        bool canSeeNear = Distance(q, Q) < R && q != Q && isFlat(N, Q->X - q->X, M_PI / 5);
                        bool veryClose = Distance(q, Q) < c * particle_spacing && q != Q;
                        return canSeeNear || veryClose;
                     });
                  }
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

#pragma omp parallel
      for (const auto &p : wall_p)
#pragma omp single nowait
      {
         p->interp_normal = Normalize(p->interp_normal_original);
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

      for (const auto &[obj, poly] : RigidBodyObject)
#pragma omp parallel
         for (const auto &q : obj->getPoints())
#pragma omp single nowait
         {
            if (q->isCaptured) {

               /*DOC_EXTRACT 0_1_2_setWall

               #### `isCaptured`が`true`の壁面粒子の流速の計算

               次のようにして，鏡写しのように流速を計算する．

               ```cpp
               q->U_SPH = Reflect(q->U_SPH, q->v_to_surface_SPH)
               ```

               */
               q->U_SPH.fill(0.);
               std::array<double, 3> DUDT;
               DUDT.fill(0.);
               double div_U = 0, div_U_next = 0;
               q->intp_density = 0.;
               double density = 0., density_next = 0.;
               double volume = 0., volume_next = 0.;
               //! ちょっとした修正11/28 -> 効果なし
               // if (q->isFirstWallLayer)
               //    q->marker_X = q->X + q->v_to_surface_SPH;
               // else
               q->marker_X = q->X + 2. * q->v_to_surface_SPH;
               double total_vol_w = 0, r = q->SML(), vol_w, vol_w_next, total_vol_w_next = 0;
               int count = 0;
               net->BucketPoints.apply(q->marker_X, 1.1 * r, [&](const auto &Q) {
                  if (Norm(q->marker_X - Q->X) > r)
                     return;
                  vol_w = Q->volume * w_Bspline(Norm(q->marker_X - Q->X), r);
                  vol_w_next = Q->volume * w_Bspline(Norm(q->marker_X - X_next(Q)), r);
                  q->U_SPH += Q->U_SPH * vol_w;
                  density += Q->rho * vol_w;
                  density_next += Q->rho * vol_w_next;
                  div_U += Q->div_U * vol_w;
                  div_U_next += Q->div_U_next * vol_w_next;
                  volume += vol_w;
                  volume_next += vol_w_next;
                  DUDT += Q->DUDt_SPH * vol_w;
                  q->intp_density += Q->rho * vol_w;
                  total_vol_w += vol_w;
                  total_vol_w_next += vol_w_next;
                  count += 1;
               });

               if (total_vol_w >= 0.01 /*should be close to unity*/) {
                  DUDT /= total_vol_w;
                  div_U /= total_vol_w;
                  q->U_SPH /= total_vol_w;
                  q->marker_U = q->U_SPH;  //! そのままの値
                  density /= total_vol_w;
                  volume /= total_vol_w;
                  if (density < 1E-10)
                     density = _WATER_DENSITY_;
               }
               if (total_vol_w_next >= 0.01) {
                  div_U_next /= total_vol_w_next;
                  density_next /= total_vol_w_next;
                  volume_next /= total_vol_w_next;
                  if (density_next < 1E-10)
                     density_next = _WATER_DENSITY_;
               }
               // q->U_SPH = -q->U_SPH;
               // q->U_SPH.fill(0.);
               q->U_SPH = Reflect(q->U_SPH, q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
               // q->DUDt_SPH = Reflect(q->DUDt_SPH, q->v_to_surface_SPH);
               q->div_U = div_U;
               q->div_U_next = div_U_next;
               // q->U_SPH = Chop(Reflect(q->U_SPH, q->v_to_surface_SPH), q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
               q->RK_U.initialize(q->RK_U.dt, q->RK_U.t_init, q->U_SPH, q->RK_U.steps);
               if (!isFinite(density))
                  density = _WATER_DENSITY_;
               if (!isFinite(density_next))
                  density_next = _WATER_DENSITY_;

               density = 0.99 * density + 0.01 * _WATER_DENSITY_;
               density_next = 0.99 * density_next + 0.01 * _WATER_DENSITY_;
               q->RK_rho.initialize(q->RK_rho.dt, q->RK_rho.t_init, density, q->RK_rho.steps);
               q->setDensityVolume(density, std::pow(particle_spacing, 3.));
               //
               // q->setDensity(_WATER_DENSITY_);
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

      for (const auto &A : net->getPoints()) {
         A->isSurface = false;
         A->isSurface_tmp = false;
         A->isSurface_next = false;
         A->isSurface_next_tmp = false;
      }

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

            const double SML = p->SML_grad();
            const double SML_next = p->SML_grad_next();

            net->BucketPoints.apply(p->X, 1.1 * SML, [&](const auto &q) {
               // w = q->volume * w_Bspline(Norm(p->X - q->X), SML);
               if (Distance(p, q) < SML) {
                  p->COM_SPH += q->mass * q->X;
                  p->totalMass_SPH += q->mass;

                  if (p != q) {
                     p->vec2COM += (q->X - p->X);
                     p->vec2COM_next += (X_next(q) - X_next(p));
                     total_w_vec2COM += 1;
                     total_w_vec2COM_next += 1;
                  }
               }
               double w = q->volume * w_Bspline(Norm(p->X - q->X), SML);
               total_volume_w += w;
               p->intp_density += q->rho * w;
               w = V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), SML_next);
               total_volume_w_next += w;
               p->intp_density_next += rho_next(q) * w;

               p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML_grad());  //! DO NOT USE grad_w_Bspline(p, q) because this is original
               p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_grad_next());
               p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p->X, q->X, p->SML_grad());

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
               obj->BucketPoints.apply(p->X, 1.1 * SML, [&](const auto &q) {
                  if (canInteract(p, q)) {
                     // w = q->volume * w_Bspline(Norm(p->X - q->X), SML);
                     if (Distance(p, q) < SML || Distance(X_next(p), X_next(q)) < SML) {
                        w = q->volume * w_Bspline(Norm(p->X - q->X), SML);
                        total_volume_w += w;
                        p->intp_density += q->rho * w;
                        w = V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), SML_next);
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
                     p->interp_normal_original -= a * q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML_grad());
                     p->interp_normal_rigid -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML_grad());
                     if (Distance(X_next(p), X_next(q)) < p->particle_spacing * 1.5)
                        near_wall_particle_next.emplace_back(q->v_to_surface_SPH);
                     p->interp_normal_original_next -= a * rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_grad_next());
                     p->interp_normal_rigid_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_grad_next());
                     p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p->X, q->X, p->SML_grad());
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

            //! 流体粒子に対してはr=2.5*ps，\theta=30度
            const double ratio_flud = 2.5;
            // 2.2OK
            // 2.4OUT
            const double min_ratio_flud = 2.25;
            const double surface_check_r_for_fluid_ = A->particle_spacing * ratio_flud;
            const double surface_check_r_for_fluid_surface = A->particle_spacing * min_ratio_flud;
            //! 壁粒子に対してはr=2*ps，\theta=30度．理由壁は長く伸びていることがあるため，長距離の粒子を認識しないようにするため．
            //! surface_check_r_for_wall_は最低の値
            //! 規則正しく並んでいる場合を考えると，水面粒子の探査範囲は，大体，std::sqrt(3.)かstd::sqrt(2.)あたりが妥当な値だろう．
            const double surface_check_r_for_wall_ = A->particle_spacing * 2.15;  // std::sqrt(3.);
            const double ratio = surface_check_r_for_wall_ / surface_check_r_for_fluid_;
            const double surface_check_agnle_for_fluid = M_PI / 5.;  // * (1. - A->var_Eigenvalues_of_M1);
            const double surface_check_angle_for_wall = M_PI / 4.;   // * (1. - A->var_Eigenvalues_of_M1);
            // if (A->var_Eigenvalues_of_M1 > 0.)
            // {
            //! 分布が均等でないものは水面になりやすくする
            double surface_check_r_for_fluid = surface_check_r_for_fluid_;  // * std::clamp(1. - 0.5 * A->var_Eigenvalues_of_M, min_ratio_flud / ratio_flud, 1.);
            double surface_check_r_for_wall = surface_check_r_for_wall_;    // * std::clamp(1. - A->var_Eigenvalues_of_M, ratio, 1.);
            // }

            // 以上に該当する粒子があった場合は，水面粒子として認識しない．
            auto n_current = Dot(p->inv_grad_corr_M, p->interp_normal_original);
            A->isSurface_tmp = A->is_grad_corr_M_singular ||
                               A->checked_points_in_radius_SPH < 3 ||
                               A->checked_points_in_radius_of_fluid_SPH < 3 ||
                               // A->var_Eigenvalues_of_M > 0.35 ||
                               (A->var_Eigenvalues_of_M > 0.1 && std::ranges::none_of(all_nets, [&](const auto &net) { return net->BucketPoints.any_of(A->X, A->particle_spacing * 3.,
                                                                                                                                                       [&](const auto &q) {
                                                                                                                                                          if (q->isCaptured && A != q) {
                                                                                                                                                             auto radius = q->isFluid ? surface_check_r_for_fluid : surface_check_r_for_wall;
                                                                                                                                                             auto angle = q->isFluid ? surface_check_agnle_for_fluid : surface_check_angle_for_wall;
                                                                                                                                                             auto c = A->X;  // + radius * n_current;
                                                                                                                                                             // Sphere sphere(c, radius);
                                                                                                                                                             // return sphere.isInside(q->X);
                                                                                                                                                             auto n_of_q = Dot(q->inv_grad_corr_M, q->interp_normal_original);
                                                                                                                                                             return Distance(c, q->X) < radius &&
                                                                                                                                                                    isFlat(n_current, q->X - A->X, angle) &&
                                                                                                                                                                    !q->is_grad_corr_M_singular &&
                                                                                                                                                                    (/*この文は，水面をキープするためのもので，これによって，水面らしい面どうしの衝突はギリギリまで水面としてキープされるはず．*/
                                                                                                                                                                     (q->isFluid &&
                                                                                                                                                                      (A->var_Eigenvalues_of_M > 0.8 || Norm(n_current) > 2E+4) &&
                                                                                                                                                                      (q->var_Eigenvalues_of_M > 0.8 || Norm(n_of_q) > 2E+4)
                                                                                                                                                                      /*確実な法線ベクトルを意味する*/)
                                                                                                                                                                         ? Dot(n_current, Dot(q->inv_grad_corr_M, q->interp_normal_original)) > 0.
                                                                                                                                                                         : true);
                                                                                                                                                          } else
                                                                                                                                                             return false;
                                                                                                                                                          // return q->isCaptured &&
                                                                                                                                                          //        Distance(A, q) < surface_check_r_for_fluid && A != q &&
                                                                                                                                                          //        isFlat(Dot(p->inv_grad_corr_M, p->interp_normal_original), q->X - A->X, M_PI / 5);
                                                                                                                                                       }); }));
            // if (A->intp_density > _WATER_DENSITY_ * 1.02)
            //    A->isSurface = false;
            // if (A->var_Eigenvalues_of_M1 > 0.6)
            //    A->isSurface = true;

            // if (A->var_Eigenvalues_of_M1_next > 0.3)
            // {
            //! 分布が均等でないものは水面になりやすくする
            surface_check_r_for_fluid = surface_check_r_for_fluid_;  // * std::clamp(1. - 0.5 * A->var_Eigenvalues_of_M_next, min_ratio_flud / ratio_flud, 1.);
            surface_check_r_for_wall = surface_check_r_for_wall_;    // * std::clamp(1. - A->var_Eigenvalues_of_M_next, ratio, 1.);
            // }
            // surfaceではないが，var_Eigenvalues_of_M1_nextが大きいものはEISPHで計算する，
            auto n_next = Dot(p->inv_grad_corr_M_next, p->interp_normal_original_next);
            A->isSurface_next = A->is_grad_corr_M_next_singular ||
                                A->checked_points_in_radius_SPH_next < 3 ||
                                A->checked_points_in_radius_of_fluid_SPH_next < 3 ||
                                //   A->var_Eigenvalues_of_M_next > 0.35 ||
                                (A->var_Eigenvalues_of_M_next > 0.1 && std::ranges::none_of(all_nets, [&](const auto &net) { return net->BucketPoints.any_of(A->X, A->particle_spacing * 3.,
                                                                                                                                                             [&](const auto &q) {
                                                                                                                                                                if (q->isCaptured && A != q) {
                                                                                                                                                                   auto radius = q->isFluid ? surface_check_r_for_fluid : surface_check_r_for_wall;
                                                                                                                                                                   auto angle = q->isFluid ? surface_check_agnle_for_fluid : surface_check_angle_for_wall;
                                                                                                                                                                   auto c = X_next(A);  // + radius * n_next;
                                                                                                                                                                   //   Sphere sphere(c, radius);
                                                                                                                                                                   //   return sphere.isInside(X_next(q));
                                                                                                                                                                   auto n_of_q = Dot(q->inv_grad_corr_M_next, q->interp_normal_original_next);
                                                                                                                                                                   return Distance(c, X_next(q)) < radius &&
                                                                                                                                                                          isFlat(n_next, X_next(q) - X_next(A), angle) &&
                                                                                                                                                                          !q->is_grad_corr_M_next_singular &&
                                                                                                                                                                          (/*この文は，水面をキープするためのもので，これによって，水面らしい面どうしの衝突はギリギリまで水面としてキープされるはず．*/
                                                                                                                                                                           (q->isFluid &&
                                                                                                                                                                            (A->var_Eigenvalues_of_M_next > 0.8 || Norm(n_next) > 2E+4) &&
                                                                                                                                                                            (q->var_Eigenvalues_of_M_next > 0.8 || Norm(n_of_q) > 2E+4 /*確実な法線ベクトルを意味する*/))
                                                                                                                                                                               ? Dot(n_next, n_of_q) > 0.
                                                                                                                                                                               : true);
                                                                                                                                                                } else
                                                                                                                                                                   return false;
                                                                                                                                                                //   return q->isCaptured &&
                                                                                                                                                                //          Distance(X_next(A), X_next(q)) < range && A != q &&
                                                                                                                                                                //          isFlat(Dot(p->inv_grad_corr_M_next, p->interp_normal_original_next), X_next(q) - X_next(A), M_PI / 5);
                                                                                                                                                             }); }));

            if (A->isSurface_tmp)
               A->isSurface_count += 1;
            else
               A->isSurface_count = 0;
            A->isSurface = (A->isSurface_count >= (A->isSurface_count_init + 1));
            A->isSurface_last_tmp = A->isSurface_tmp;
            A->isSurface_next_tmp = A->isSurface_next;

            // if (A->intp_density_next > _WATER_DENSITY_ * 1.02)
            //    A->isSurface_next = false;
            // if (A->var_Eigenvalues_of_M1_next > 0.5)
            //    A->isSurface_next = true;

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

            auto surface = net->BucketPoints.none_of(p->X, SML, [&](const auto &q) {
               return Distance(p, q) < r_short && p != q && (isFlat(p->interp_normal_water, q->X - p->X, M_PI / 5));
            });

            if (surface)
               p->isNeumannSurface = std::ranges::any_of(RigidBodyObject, [&](const auto &net) {
                  return net->BucketPoints.any_of(p->X, SML, [&](const auto &q) {
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
            net->BucketPoints.apply(p->X, p->SML_grad(), [&](const auto &q) {
               p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML_grad());  //! DO NOT USE grad_w_Bspline(p, q) because this is original
               p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_grad_next());
            });
            for (const auto &obj : RigidBodyObject)
               obj->BucketPoints.apply(p->X, p->SML_grad(), [&](const auto &q) {
                  if (canInteract(p, q)) {
                     p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->SML_grad());
                     p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->SML_grad_next());
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