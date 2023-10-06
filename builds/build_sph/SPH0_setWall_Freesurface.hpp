#ifndef SPH_setWall_Freesurface_H
#define SPH_setWall_Freesurface_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/* -------------------------------------------------------------------------- */

void setCorrectionMatrix(networkPoint *A, const std::vector<Network *> &all_nets) {

   /*DOC_EXTRACT 0_1_1_set_free_surface

   ### 勾配演算子の修正をする行列の計算

   \cite{Morikawa2023}で紹介されていた，\cite{Randles1996}の勾配演算の精度を改善する行列を計算する．
   `grad_corr_M`としている．

   ```math
   \begin{align}
   {\bf M} = \left(\sum_j V_j ({\bf x}_j-{\bf x}_i) \otimes \nabla W_{ij}\right)^{-1}
   \end{align}
   ```

   */

   // キャプチャされなかった場合
   A->grad_U.fill({0., 0., 0.});
   IdentityMatrix(A->grad_corr_M);
   IdentityMatrix(A->grad_corr_M_next);
   IdentityMatrix(A->inv_grad_corr_M);
   IdentityMatrix(A->inv_grad_corr_M_next);
   IdentityMatrix(A->Mat1);
   IdentityMatrix(A->Mat2);
   IdentityMatrix(A->Mat3);
   IdentityMatrix(A->Mat_B);

   // キャプチャーされた粒子のみを対象にする．
   if (!A->isCaptured)
      return;

   // 初期化
   A->grad_U.fill({0., 0., 0.});
   A->grad_corr_M.fill({0., 0., 0.});
   A->grad_corr_M_next.fill({0., 0., 0.});
   A->Mat1.fill({0., 0., 0.});
   A->Mat2.fill({0., 0., 0.});
   A->Mat3.fill({0., 0., 0.});
   A->Mat_B.fill({0., 0., 0.});

   Tddd Xij, nz_Xij;
   auto add = [&](const auto &B) {
      auto grad_w = grad_w_Bspline(A->X, B->X, A->radius_SPH);
      A->grad_corr_M += B->volume * TensorProduct(B->X - A->X, grad_w);
      A->grad_corr_M_next += V_next(B) * TensorProduct(X_next(B) - X_next(A), grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH));

      const auto Uij = A->U_SPH - B->U_SPH;
      A->grad_U += B->volume * TensorProduct(-Uij, grad_w);

      Xij = A->X - B->X;
      nz_Xij = Normalize(Xij);
      A->Mat1 += B->volume * TensorProduct(Xij * nz_Xij * nz_Xij, grad_w);
      A->Mat2 += B->volume * TensorProduct(nz_Xij * nz_Xij, grad_w);
      A->Mat3 += B->volume * TensorProduct(Xij * Xij, grad_w);
   };

   for (const auto &NET : all_nets)
      NET->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
         if (B->isCaptured) add(B);
      });

   A->Eigenvalues_of_M = {0., 0., 0.};
   A->inv_grad_corr_M = Inverse(A->grad_corr_M);
   A->inv_grad_corr_M_next = Inverse(A->grad_corr_M_next);

   {
      auto [lambdas, vectors] = Eigensystem(A->grad_corr_M, 1E-4, 100.);
      A->Eigenvalues_of_M = lambdas;
      A->Eigenvectors_of_M = vectors;
   }
   {
      auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M, 1E-4, 100.);
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
   //
   if (A->isFluid)
      A->isSurface = A->var_Eigenvalues_of_M1 > 0.175 && A->min_Eigenvalues_of_M1 > 0.9;
   //
   T3Tddd I;
   IdentityMatrix(I);
   A->Mat_B = Dot(-I, Inverse(A->Mat1 + Dot(Dot(A->Mat2, A->inv_grad_corr_M), A->Mat3)));

   /* --------------------------------------------------------- */
   if (false) {
      auto add_to_map = [&](const auto &p, const auto &c) {
         if (A->map_p_grad.find(p) == p->map_p_grad.end())
            A->map_p_grad[p] = c;
         else
            A->map_p_grad[p] += c;
      };

      auto add = [&](const auto &B) {
         auto c = B->volume * Dot(A->X - B->X, Dot(grad_w_Bspline(A->X, B->X, A->radius_SPH), A->inv_grad_corr_M));
         add_to_map(B, -c);
         add_to_map(A, c);
         // 符号についての疑問が残る．逆になっているように思える．
      };

      A->map_p_grad.clear();
      for (const auto &net : all_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) add(B);
         });

      // copy    std::unordered_map<networkPoint *, double> map_p_grad to  std::vector<std::tuple<networkPoint *, double>> vector_p_grad;
      A->vector_p_grad.clear();
      for (const auto &[p, v] : A->map_p_grad)
         A->vector_p_grad.push_back(std::tuple<networkPoint *, double>{p, v});
   }
};

/* -------------------------------------------------------------------------- */

void setVectorToPolygon(const auto &water, const auto &RigidBodyObject, const double particle_spacing) {

   for (const auto &p : water->getPoints())
      p->clearContactFaces();

   //!!! 衝突の判定がよくエラーが出る箇所
   for (const auto &[_, body] : RigidBodyObject) {
#pragma omp parallel
      for (const auto &p : water->getPoints())
#pragma omp single nowait
      {
         //! ここも重要：点と面の衝突をどのようにすれば矛盾なく判定できるか．
         // \label{BEM:detection_range}
         p->radius = particle_spacing;
         p->addContactFaces(body->getBucketFaces(), false);
      }
   }

   for (const auto &[_, body] : RigidBodyObject) {
#pragma omp parallel
      for (const auto &p : water->getPoints())
#pragma omp single nowait
      {
         std::unordered_map<networkLine *, int> shared_lines;
         std::unordered_map<networkPoint *, int> shared_points;
         if (p->getContactFaces().size() > 1)
            for (const auto &f : p->getContactFaces()) {
               // count
               for (const auto &l : f->getLines())
                  if (shared_lines.find(l) != shared_lines.end())
                     shared_lines[l]++;
                  else
                     shared_lines[l] = 1;

               // count
               for (const auto &p : f->getPoints())
                  if (shared_points.find(p) != shared_points.end())
                     shared_points[p]++;
                  else
                     shared_points[p] = 1;
            }

         // erase elements in shared_lines and shared_points if the number is 1
         for (auto it = shared_lines.begin(); it != shared_lines.end();) {
            if (it->first->isAdjacentFacesFlat(M_PI / 180))
               it = shared_lines.erase(it);
            else
               ++it;
         }

         for (auto it = shared_points.begin(); it != shared_points.end();) {
            if (it->second <= 2 || Between(it->first->getSolidAngle(), {2 * M_PI * 0.8, 2 * M_PI * 1.2}))
               it = shared_points.erase(it);
            else
               ++it;
         }

         p->vector_to_polygon.clear();
         p->vector_to_polygon_next.clear();
         // store vectors to polygon into p->vector_to_polygon using Nearest function

         // to faces
         Tddd vec;
         for (const auto &f : p->getContactFaces()) {
            vec = Nearest(p->X, ToX(f)) - p->X;
            if (Norm(vec) < p->radius && isFlat(vec, -f->normal, M_PI / 180. * 60)) {
               p->vector_to_polygon.push_back(vec);
               p->vector_to_polygon_next.push_back(Nearest(X_next(p), ToX(f)) - X_next(p));
            }
         }
         // to lines
         for (const auto &[l, _] : shared_lines) {
            vec = Nearest(p->X, ToX(l)) - p->X;
            if (Norm(vec) < p->radius && isFlat(vec, -l->getNormal(), M_PI / 180. * 60)) {
               if (std::ranges::none_of(p->vector_to_polygon, [&](const auto &v) { return isFlat(v, vec, M_PI / 180. * 20); })) {
                  p->vector_to_polygon.push_back(vec);
                  p->vector_to_polygon_next.push_back(Nearest(X_next(p), ToX(l)) - X_next(p));
               }
            }
         }

         // to points
         for (const auto &[q, _] : shared_points) {
            vec = q->X - p->X;
            if (Norm(vec) < p->radius && isFlat(vec, -p->getNormal_BEM(), M_PI / 180. * 20))
               if (std::ranges::none_of(p->vector_to_polygon, [&](const auto &v) { return isFlat(v, vec, M_PI / 180. * 20); })) {
                  p->vector_to_polygon.push_back(vec);
                  p->vector_to_polygon_next.push_back(q->X - X_next(p));
               }
         }
      }
   }
};

//! -------------------------------------------------------------------------- */

void setWall(const auto &net, const auto &RigidBodyObject, const auto &particle_spacing, auto &wall_p) {

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
         p->isFluid = false;
         p->isFreeFalling = false;
         p->isCaptured = p->isCaptured_ = false;
         p->isSurface = false;
         p->p_SPH = 0;
         p->U_SPH = p->DUDt_SPH = p->lap_U = {0, 0, 0};
         p->tmp_X = p->X;
         p->isFirstWallLayer = false;
         p->isChecked = false;
      }

   DebugPrint("isFirstWallLayerを決める", Green);

#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      const double captureRange = 1.42 * p->radius_SPH;  //\label{SPH:capture_condition_1st}
      // const double captureRange_wall_as_fluid = p->radius_SPH;
      for (const auto &[obj, poly] : RigidBodyObject) {
         obj->BucketPoints.apply(p->X, captureRange, [&](const auto &q) {
            if (Distance(p, q) < captureRange && !q->isChecked) {

               q->isChecked = true;

               /*DOC_EXTRACT 0_1_0_set_normal_of_wall

               ## 壁面粒子の抽出と値の計算

               ### `interp_normal_original`の計算

               流体粒子と同じ影響半径を使ってしまうと，流体粒子が参照できる範囲ギリギリにある壁粒子の法線方向の値が不正確になる．
               そのため，流体粒子の影響半径よりも広い半径を使って，`q->interp_normal_original`の法線方向を計算することが，重要である．
               少し大きい半径を`captureRange`としている．

               */

               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(q->X, captureRange,
                                          [&](const auto &Q) {
                                             q->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(q->X, Q->X, captureRange);
                                          });
               const auto N = q->interp_normal_original;
               q->interp_normal = Normalize(q->interp_normal_original);
               if (!isFinite(q->interp_normal))
                  q->interp_normal = {0., 0., 1.};

               /*DOC_EXTRACT 0_1_0_set_normal_of_wall

               ### `isCaptured`の決定

               法線方向`interp_normal_original`を使って，流体粒子に近くかつ向かい合う方向にある壁粒子を抽出する．
               計算に使用する壁粒子を決定し，使用する場合`isCaptured`を`true`にする．

               */
               q->isCaptured = false;
               q->isFirstWallLayer = false;
               q->isNeumannSurface = false;

               auto R = captureRange;
               q->isCaptured = net->BucketPoints.any_of(q->X, R,
                                                        [&](const auto &Q) {
                                                           bool canSeeNear = Distance(q, Q) < R && q != Q && isFlat(N, Q->X - q->X, M_PI / 6);
                                                           bool veryClose = Distance(q, Q) < 1.25 * particle_spacing && q != Q;
                                                           return canSeeNear || veryClose;
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
               /* -------------------------------------------------------------------------- */

               /*DOC_EXTRACT 0_1_0_set_normal_of_wall

               ### `isCaptured`が`true`の壁面粒子の流速の計算

               次のようにして，鏡写しのように流速を計算する．

               ```cpp
               q->U_SPH = Reflect(q->U_SPH, q->normal_SPH)
               ```

               */

               if (q->isCaptured) {
                  q->U_SPH.fill(0.);
                  auto X = q->X + 2 * q->normal_SPH;
                  double total_w = 0, r = q->radius_SPH, w;
                  net->BucketPoints.apply(X, 1.1 * r, [&](const auto &Q) {
                     q->U_SPH += Q->U_SPH * (w = Q->volume * w_Bspline(Norm(X - Q->X), r));
                     total_w += w;
                  });

                  if (total_w == 0.)
                     q->U_SPH.fill(0.);
                  else
                     q->U_SPH /= total_w;
                  // q->U_SPH = -q->U_SPH;
                  // q->U_SPH.fill(0.);
                  q->U_SPH = Reflect(q->U_SPH, q->normal_SPH);  //\label{SPH:wall_particle_velocity}
                  // if (q->isFirstWallLayer)
                  //    q->U_SPH = Chop(q->U_SPH, q->normal_SPH);
                  // q->U_SPH = Projection(Reflect(q->U_SPH, q->normal_SPH), q->normal_SPH);  //\label{SPH:wall_particle_velocity}
               }
            }
         });
      }
   };

   // DebugPrint("壁粒子のオブジェクト外向き法線方向を計算", Green);

   /*DOC_EXTRACT 0_1_0_set_normal_of_wall

   ### `setCorrectionMatrix`で壁粒子の演算修正用行列を計算

   `setCorrectionMatrix`で壁粒子の演算修正用行列を計算する．

   */

   wall_p.clear();
   for (const auto &[obj, poly] : RigidBodyObject)
      for (const auto &p : obj->getPoints()) {
         if (p->isCaptured) {
            wall_p.emplace(p);
         }
         setCorrectionMatrix(p, all_nets);
      }
};

/*DOC_EXTRACT 0_1_1_set_free_surface

## 流体の法線方向の計算と水面の判定

*/

void setFreeSurface(auto &net, const auto &RigidBodyObject) {

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

         setCorrectionMatrix(A, all_nets);

         /* -------------------------------------------------------------------------- */

         /*DOC_EXTRACT 0_1_1_set_free_surface

         ### 流体の法線方向の計算

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
         double w;
         // std::vector<networkPoint *> samples;
         net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            // w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
            double w;
            if (Distance(p, q) < p->radius_SPH) {
               p->COM_SPH += q->mass * q->X;
               p->totalMass_SPH += q->mass;
            }
            p->intp_density += q->rho * q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
            p->intp_density_next += rho_next(q) * V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), p->radius_SPH);
            p->interp_normal_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            p->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
            p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            //
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
            obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
               if (q->isCaptured) {
                  // w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
                  if (Distance(p, q) < p->radius_SPH) {
                     p->intp_density += q->rho * q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
                     p->intp_density_next += rho_next(q) * V_next(q) * w_Bspline(Norm(X_next(p) - X_next(q)), p->radius_SPH);
                     p->COM_SPH += q->mass * q->X;
                     p->totalMass_SPH += q->mass;
                  }
                  if (Distance(p, q) < p->radius_SPH / p->C_SML * 1.5)
                     near_wall_particle.emplace_back(q->normal_SPH);
                  p->interp_normal_original -= a * q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                  p->interp_normal_rigid -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                  //
                  if (Distance(X_next(p), X_next(q)) < p->radius_SPH / p->C_SML * 1.5)
                     near_wall_particle_next.emplace_back(q->normal_SPH);
                  p->interp_normal_original_next -= a * rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
                  p->interp_normal_rigid_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
                  //
                  p->intp_normal_Eigen += q->var_Eigenvalues_of_M * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
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

         /*DOC_EXTRACT 0_1_1_set_free_surface

         ### 水面の判定

         `surface_condition0,1`の両方を満たす場合，水面とする．

         */

         /* ------------------------ Neumann surface condition ----------------------- */

         const auto r_short = (p->radius_SPH / p->C_SML) * 2.;

         auto surface = net->BucketPoints.none_of(p->X, p->radius_SPH, [&](const auto &q) {
            return Distance(p, q) < r_short && p != q && (isFlat(p->interp_normal_water, q->X - p->X, M_PI / 5));
         });

         if (surface)
            p->isNeumannSurface = std::ranges::any_of(RigidBodyObject, [&](const auto &net) {
               return std::get<0>(net)->BucketPoints.any_of(p->X, p->radius_SPH, [&](const auto &q) {
                  return Distance(p, q) < r_short && p != q && (isFlat(p->interp_normal_water, q->X - p->X, M_PI / 5));
               });
            });
         else
            p->isNeumannSurface = false;

         /* -------------------------- Surface condition ---------------------------- */
         // 上の判定は余計に水面を判定してしまうので，ここで内部のものを水面粒子から除外する．
         if (p->isSurface &&
             p->var_Eigenvalues_of_M1 < 0.4 /*誤診ではないか？*/) {
            const auto radius = p->radius_SPH;  //(p->radius_SPH / p->C_SML) * 2.;
            auto surface_condition0 = [&](const auto &q) {
               return q->isCaptured &&
                      Distance(p, q) < radius &&
                      p != q &&
                      isFlat(p->interp_normal_original_choped, q->X - p->X, M_PI / 5);
            };
            p->isSurface = net->BucketPoints.none_of(p->X, radius, surface_condition0) &&
                           std::ranges::none_of(RigidBodyObject, [&](const auto &net) {
                              return std::get<0>(net)->BucketPoints.any_of(p->X, radius, surface_condition0);
                           });
         }
         if (p->intp_density_next < 0.5 * _WATER_DENSITY_)
            p->isSurface = true;
      }

   DebugPrint("水粒子のオブジェクト外向き法線方向を計算 done", Green);
};

#endif