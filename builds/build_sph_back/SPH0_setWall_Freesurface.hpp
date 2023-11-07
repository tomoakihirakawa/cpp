#ifndef SPH_setWall_Freesurface_H
#define SPH_setWall_Freesurface_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT 0_1_1_preparation_wall_and_freesurface

## N.S.方程式を解く前の準備

*/

/* -------------------------------------------------------------------------- */

void setCorrectionMatrix(networkPoint *A, const std::vector<Network *> &all_nets) {

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
   //
   // キャプチャされなかった場合
   A->grad_U.fill({0., 0., 0.});
   IdentityMatrix(A->grad_corr_M);
   IdentityMatrix(A->grad_corr_M_next);
   IdentityMatrix(A->inv_grad_corr_M);
   IdentityMatrix(A->inv_grad_corr_M_next);
   //
   IdentityMatrix(A->inv_grad_corr_M);
   IdentityMatrix(A->inv_grad_corr_M_next);
   //
   IdentityMatrix(A->grad_corr_M_rigid);
   IdentityMatrix(A->grad_corr_M_next_rigid);
   IdentityMatrix(A->inv_grad_corr_M_rigid);
   IdentityMatrix(A->inv_grad_corr_M_next_rigid);
   //
   // IdentityMatrix(A->Mat1);
   // IdentityMatrix(A->Mat2);
   // IdentityMatrix(A->Mat3);
   // IdentityMatrix(A->Mat_B);

   // キャプチャーされた粒子のみを対象にする．
   if (!A->isCaptured)
      return;

   // 初期化
   A->grad_U.fill({0., 0., 0.});
   A->grad_corr_M.fill({0., 0., 0.});
   A->grad_corr_M_next.fill({0., 0., 0.});
   //
   A->grad_corr_M_mirror.fill({0., 0., 0.});
   A->grad_corr_M_next_mirror.fill({0., 0., 0.});
   //
   A->grad_corr_M_rigid.fill({0., 0., 0.});
   A->grad_corr_M_next_rigid.fill({0., 0., 0.});

   // A->Mat1.fill({0., 0., 0.});
   // A->Mat2.fill({0., 0., 0.});
   // A->Mat3.fill({0., 0., 0.});
   // A->Mat_B.fill({0., 0., 0.});

   Tddd Xij, Uij;
   auto add = [&](const auto &B) {
      auto grad_w = grad_w_Bspline(A->X, B->X, A->SML());  //! DO NOT USE grad_w_Bspline_next(A, B)
      if (B->isCaptured) {
         A->grad_corr_M += B->volume * TensorProduct(B->X - A->X, grad_w);
         A->grad_corr_M_next += V_next(B) * TensorProduct(X_next(B) - X_next(A), grad_w_Bspline(X_next(A), X_next(B), A->SML_next()));  //! DO NOT USE grad_w_Bspline_next(A, B)
         //! mirrror
         A->grad_corr_M_mirror += B->volume * TensorProduct(B->X - A->X, grad_w);
         A->grad_corr_M_next_mirror += V_next(B) * TensorProduct(X_next(B) - X_next(A), grad_w_Bspline(X_next(A), X_next(B), A->SML_next()));  //! DO NOT USE grad_w_Bspline_next(A, B)
         if (A->isFluid && !A->isSurface && B->isSurface) {
            auto BX = Mirror(A->X, B->X, B->X - A->X);
            auto grad_w = grad_w_Bspline(A->X, BX, A->SML());  //! DO NOT USE grad_w_Bspline_next(A, B)
            A->grad_corr_M_mirror += B->volume * TensorProduct(BX - A->X, grad_w);
            auto BX_next = Mirror(X_next(A), X_next(B), X_next(B) - X_next(A));
            auto grad_w_next = grad_w_Bspline(X_next(A), BX_next, A->SML_next());  //! DO NOT USE grad_w_Bspline_next(A, B)
            A->grad_corr_M_next_mirror += V_next(B) * TensorProduct(BX_next - X_next(A), grad_w_next);
         }
         Uij = A->U_SPH - B->U_SPH;
         A->grad_U += B->volume * TensorProduct(-Uij, grad_w);
      }
      if (B->getNetwork()->isRigidBody) {
         A->grad_corr_M_rigid += B->volume * TensorProduct(B->X - A->X, grad_w);
         A->grad_corr_M_next_rigid += V_next(B) * TensorProduct(X_next(B) - X_next(A), grad_w_Bspline(X_next(A), X_next(B), A->SML_next()));  //! DO NOT USE grad_w_Bspline_next(A, B)
      }
   };

   std::ranges::for_each(all_nets, [&](const auto &NET) {
      NET->BucketPoints.apply(A->X, 1.2 * A->SML(), [&](const auto &B) { add(B); });
   });

   A->Eigenvalues_of_M = {0., 0., 0.};
   A->inv_grad_corr_M = Inverse(A->grad_corr_M);
   A->inv_grad_corr_M_next = Inverse(A->grad_corr_M_next);
   // mirror
   A->inv_grad_corr_M_mirror = Inverse(A->grad_corr_M_mirror);
   A->inv_grad_corr_M_next_mirror = Inverse(A->grad_corr_M_next_mirror);
   //
   A->inv_grad_corr_M_rigid = Inverse(A->grad_corr_M_rigid);
   A->inv_grad_corr_M_next_rigid = Inverse(A->grad_corr_M_next_rigid);

   {
      auto [lambdas, vectors] = Eigensystem(A->grad_corr_M, 1E-10, 100.);
      A->Eigenvalues_of_M = lambdas;
      A->Eigenvectors_of_M = vectors;
   }

   {
      auto [lambdas, vectors] = Eigensystem(A->grad_corr_M_rigid, 1E-10, 100.);
      A->Eigenvalues_of_M_rigid = lambdas;
      A->Eigenvectors_of_M_rigid = vectors;
   }

   {
      auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M_next, 1E-10, 100.);
      A->Eigenvalues_of_M1 = lambdas;
      A->Eigenvectors_of_M1 = vectors;
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
   A->min_Eigenvalues_of_M = Min(A->Eigenvalues_of_M);
   A->min_Eigenvalues_of_M1 = Min(A->Eigenvalues_of_M1);
};

//! -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_1_2_setWall
### 壁面粒子の抽出と値の計算
*/

void setWall(const auto &net, const auto &RigidBodyObject, const auto &particle_spacing, auto &wall_p) {
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
            p->grad_Min_gradM.fill(0.);
            p->isFluid = false;
            p->isFreeFalling = false;
            p->isCaptured = p->isCaptured_ = false;
            p->isSurface = false;
            p->p_SPH = 0;
            p->U_SPH = p->DUDt_SPH = p->lap_U = {0, 0, 0};
            p->tmp_X = p->X;
            p->isFirstWallLayer = false;
            p->isChecked = false;
            p->intp_density = 0.;
         }

      DebugPrint("isFirstWallLayerを決める", Green);

#pragma omp parallel
      for (const auto &p : net->getPoints())
#pragma omp single nowait
      {
         const double captureRange = 1.42 * p->SML();  //\label{SPH:capture_condition_1st}
         // const double captureRange_wall_as_fluid = p->SML();
         for (const auto &[obj, poly] : RigidBodyObject) {
            obj->BucketPoints.apply(p->X, captureRange, [&](const auto &q) {
               p->grad_Min_gradM -= p->volume * (Min(p->Eigenvalues_of_M_rigid) - Min(q->Eigenvalues_of_M_rigid)) * grad_w_Bspline(p, q, captureRange);
               if (Distance(p, q) < captureRange && !q->isChecked) {

                  q->isChecked = true;

                  /*DOC_EXTRACT 0_1_2_setWall

                  #### `interp_normal_original`の計算

                  流体粒子と同じ影響半径を使ってしまうと，流体粒子が参照できる範囲ギリギリにある壁粒子の法線方向の値が不正確になる．
                  そのため，流体粒子の影響半径よりも広い半径を使って，`q->interp_normal_original`の法線方向を計算することが，重要である．
                  少し大きい半径を`captureRange`としている．

                  */

                  for (const auto &[obj, poly] : RigidBodyObject)
                     obj->BucketPoints.apply(q->X, captureRange,
                                             [&](const auto &Q) {
                                                q->interp_normal_original -= Q->rho * Q->volume * grad_w_Bspline(q, Q, captureRange);
                                             });
                  const auto N = q->interp_normal_original;
                  q->interp_normal = Normalize(q->interp_normal_original);
                  if (!isFinite(q->interp_normal))
                     q->interp_normal = {0., 0., 1.};

                  /*DOC_EXTRACT 0_1_2_setWall

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
                                                              bool canSeeNear = Distance(q, Q) < R && q != Q && isFlat(N, Q->X - q->X, M_PI / 8);
                                                              bool veryClose = Distance(q, Q) < std::sqrt(1.5) * particle_spacing && q != Q;

                                                              bool canSeeNear_next = Distance(X_next(q), X_next(Q)) < R && q != Q && isFlat(N, X_next(Q) - X_next(q), M_PI / 8);
                                                              bool veryClose_next = Distance(X_next(q), X_next(Q)) < std::sqrt(1.5) * particle_spacing && q != Q;

                                                              return canSeeNear || veryClose || canSeeNear_next || veryClose_next;
                                                           });

                  // ここを修正した
                  // if (q->isCaptured) {
                  //    auto X = q->X + 2 * q->normal_SPH;
                  //    q->isCaptured = net->BucketPoints.any_of(X, particle_spacing * 1.25, [&](const auto &Q) {
                  //       return Distance(Q, X) < particle_spacing * 1.25 && q != Q;
                  //    });
                  // }

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
                  if (q->isFirstWallLayer) {
                     q->v_to_surface_SPH = q->normal_SPH = q->particle_spacing * 0.5 * Normalize(q->interp_normal_original);
                  } else if (q->isCaptured) {
                     q->v_to_surface_SPH = q->normal_SPH = q->particle_spacing * 1.5 * Normalize(q->interp_normal_original);
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
                     q->intp_density = 0.;
                     auto markerX = q->X + 2 * q->v_to_surface_SPH;
                     double total_w = 0, r = q->SML(), w;
                     net->BucketPoints.apply(markerX, 1.1 * r, [&](const auto &Q) {
                        q->U_SPH += Q->U_SPH * (w = Q->volume * w_Bspline(Norm(markerX - Q->X), r));
                        q->intp_density += Q->rho * Q->volume * w_Bspline(Norm(markerX - Q->X), r);

                        total_w += w;
                     });
                     if (total_w == 0.) {
                        q->U_SPH.fill(0.);
                        // q->intp_density = _WATER_DENSITY_;
                     } else {
                        q->U_SPH /= total_w;
                        // q->intp_density /= total_w;
                     }
                     // q->U_SPH = -q->U_SPH;
                     // q->U_SPH.fill(0.);
                     q->U_SPH = Reflect(q->U_SPH, q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
                     // if (q->isFirstWallLayer)
                     //    q->U_SPH = Chop(q->U_SPH, q->v_to_surface_SPH);
                     // q->U_SPH = Projection(Reflect(q->U_SPH, q->v_to_surface_SPH), q->v_to_surface_SPH);  //\label{SPH:wall_particle_velocity}
                  }
               }
            });
         }
         p->grad_Min_gradM = -Normalize(p->grad_Min_gradM);
      };

      // DebugPrint("壁粒子のオブジェクト外向き法線方向を計算", Green);

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
#pragma omp parallel
         for (const auto &p : obj->getPoints())
#pragma omp single nowait
            setCorrectionMatrix(p, all_nets);

   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in setWall");
   };
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

            A->isAuxiliary = false;

            /* -------------------------------------------------------------------------- */
            /*DOC_EXTRACT 0_1_3_set_free_surface

            #### `setCorrectionMatrix`で流体粒子の演算修正用行列を計算

            `setCorrectionMatrix`で壁粒子の演算修正用行列を計算する．

            */

            setCorrectionMatrix(A, all_nets);

            // if (A->isFluid)
            //    A->isSurface = A->var_Eigenvalues_of_M1 > 0.3;
            /* -------------------------- Surface condition ---------------------------- */
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
            A->isSurface = std::ranges::none_of(all_nets,
                                                [&](const auto &net) {
                                                   return net->BucketPoints.any_of(A->X, A->particle_spacing * 2.,
                                                                                   [&](const auto &q) {
                                                                                      return q->isCaptured &&
                                                                                             Distance(A, q) < r && A != q &&
                                                                                             isFlat(A->interp_normal_original, q->X - A->X, M_PI / 5);
                                                                                   });
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

            // p->isNotSurfaceButNearSurface = false;

            net->BucketPoints.apply(p->X, p->SML(), [&](const auto &q) {
               // w = q->volume * w_Bspline(Norm(p->X - q->X), p->SML());

               // if (p->isFluid && !p->isSurface && q->isSurface && Distance(p, q) < p->SML())
               //    p->isNotSurfaceButNearSurface = true;

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

               p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p, q);
               p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline_next(p, q);
               p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p, q);

               auto A = p;
               auto B = q;
               // if (A->isNotSurfaceButNearSurface && !A->isSurface && B->isSurface) {
               //    auto BX = Mirror(A->X, B->X, B->X - A->X);
               //    w = q->volume * w_Bspline(Norm(p->X - BX), p->SML());
               //    total_weight += w;
               //    p->intp_density += q->rho * w;
               //    //

               //    auto BX_next = Mirror(X_next(A), X_next(B), X_next(B) - X_next(A));
               //    w = V_next(q) * w_Bspline(Norm(X_next(p) - BX_next), p->SML_next());
               //    total_weight_next += w;
               //    p->intp_density_next += rho_next(q) * w;
               //    // auto BX_next = Mirror(X_next(A), X_next(B), X_next(B) - X_next(A));
               //    // auto grad_w_next = grad_w_Bspline(X_next(A), BX_next, A->SML());
               //    // A->grad_corr_M_next_mirror += V_next(B) * TensorProduct(BX_next - X_next(A), grad_w_next);
               // }

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
                  if (q->isCaptured) {
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
                     p->interp_normal_original -= a * q->rho * q->volume * grad_w_Bspline(p, q);
                     p->interp_normal_rigid -= q->rho * q->volume * grad_w_Bspline(p, q);
                     if (Distance(X_next(p), X_next(q)) < p->particle_spacing * 1.5)
                        near_wall_particle_next.emplace_back(q->v_to_surface_SPH);
                     p->interp_normal_original_next -= a * rho_next(q) * V_next(q) * grad_w_Bspline_next(p, q);
                     p->interp_normal_rigid_next -= rho_next(q) * V_next(q) * grad_w_Bspline_next(p, q);
                     p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p, q);
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