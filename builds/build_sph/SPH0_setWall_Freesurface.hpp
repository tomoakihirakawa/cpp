#ifndef SPH_setWall_Freesurface_H
#define SPH_setWall_Freesurface_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

//! -------------------------------------------------------------------------- */

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
               if (std::ranges::none_of(p->vector_to_polygon, [&](const auto &v) { return VectorAngle(v, vec) < M_PI / 180. * 20; })) {
                  p->vector_to_polygon.push_back(vec);
                  p->vector_to_polygon_next.push_back(Nearest(X_next(p), ToX(l)) - X_next(p));
               }
            }
         }

         // to points
         for (const auto &[q, _] : shared_points) {
            vec = q->X - p->X;
            if (Norm(vec) < p->radius && VectorAngle(vec, -p->getNormal_BEM()) < M_PI / 180. * 20)
               if (std::ranges::none_of(p->vector_to_polygon, [&](const auto &v) { return VectorAngle(v, vec) < M_PI / 180. * 20; })) {
                  p->vector_to_polygon.push_back(vec);
                  p->vector_to_polygon_next.push_back(q->X - X_next(p));
               }
         }
      }
   }
};

//! -------------------------------------------------------------------------- */

std::array<double, 3> optimizeFunction(const std::array<double, 3> &X0,
                                       const double &vol0,
                                       const auto &p,
                                       const double radius,
                                       const auto &objects) {
   std::array<double, 3> normal0 = {0., 0., 0.};
   double rho0 = 0.;
   for (const auto &obj : objects) {
      obj->BucketPoints.apply(p->X, radius, [&](const auto &q) {
         if (Distance(p, q) < radius) {
            auto mass = q->rho * q->volume;
            normal0 -= mass * grad_w_Bspline(q->X, p->X, radius);
            rho0 += mass * w_Bspline(Norm(q->X - p->X), radius);
         }
      });
   }

   auto X = X0;
   auto vol = vol0;

   // volとXを更新する
   auto normal = normal0 - _WATER_DENSITY_ * vol * grad_w_Bspline(X, p->X, radius);
   auto rho = rho0 + _WATER_DENSITY_ * vol * w_Bspline(Norm(X - p->X), radius);
   // std::cout << "zero is better, initial normal0 = " << normal0 / _WATER_DENSITY_ << std::endl;
   // std::cout << "1 is better, initial rho0 = " << (rho0 - _WATER_DENSITY_) / _WATER_DENSITY_ << std::endl;
   // std::cout << "zero is better, initial normal = " << normal / _WATER_DENSITY_ << std::endl;
   // std::cout << "1 is better, initial rho = " << (rho - _WATER_DENSITY_) / _WATER_DENSITY_ << std::endl;
   // auto tmp = p->radius_SPH / p->C_SML + std::pow(vol, 1. / 3.);

   return normal / _WATER_DENSITY_ + (rho - _WATER_DENSITY_) / _WATER_DENSITY_;  // + (Norm(X - p->X) - tmp) / tmp;
};

/*DOC_EXTRACT SPH

## 壁面粒子の流速と圧力

壁粒子の流速を流体粒子の流速に応じて変化させるとプログラムが煩雑になるので，**ここでは**壁面粒子の流速は常にゼロに設定することにする．
壁粒子の圧力は，水が圧縮しないように各ステップ毎に計算し直す必要がある．

**フリースリップ条件の設定**

\ref{SPH:freeslip}{フリースリップ条件の設定}

*/

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
         p->setDensityVolume(0, 0);
         p->isAir = p->isFluid = false;
         p->isFreeFalling = false;
         p->isCaptured = p->isCaptured_ = false;
         p->isSurface = false;
         p->p_SPH = 0;
         p->U_SPH = p->DUDt_SPH = p->lap_U = {0, 0, 0};
         p->tmp_X = p->X;
         p->isFirstWallLayer = false;
      }

   DebugPrint("isFirstWallLayerを決める", Green);

#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      // 安定かを入れよう．
      // ここでも結構変わる
      const double captureRange = 1.1 * p->radius_SPH;  //\label{SPH:capture_condition_1st}
      // const double captureRange_wall_as_fluid = p->radius_SPH;
      for (const auto &[obj, poly] : RigidBodyObject) {
         obj->BucketPoints.apply(p->X, captureRange, [&](const auto &q) {
            if (Distance(p, q) < captureRange) {

               q->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3.));
               q->interpolated_normal_SPH_original.fill(0.);

               // captureされた点の法線方向を計算
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(q->X, q->radius_SPH, [&](const auto &Q) {
                     q->interpolated_normal_SPH_original -= q->rho * q->volume * grad_w_Bspline(q->X, Q->X, Q->radius_SPH);
                  });

               q->interpolated_normal_SPH = Normalize(q->interpolated_normal_SPH_original);
               if (!isFinite(q->interpolated_normal_SPH))
                  q->interpolated_normal_SPH = {0., 0., 1.};

               // captureされた点のうち流体粒子に近い点は第１層目としてマーク
               auto firstWallLayerCondition = [&](const auto &Q) {
                  return Distance(q, Q) < particle_spacing * 1.8 && q != Q && (VectorAngle(q->interpolated_normal_SPH, Q->X - q->X) < M_PI / 5);
               };
               q->isFirstWallLayer = net->BucketPoints.any_of(q->X, q->radius_SPH, firstWallLayerCondition);
               /* -------------------------------------------------------------------------- */
               const auto r_short = (q->radius_SPH / q->C_SML) * 2.;
               q->isNeumannSurface = net->BucketPoints.any_of(q->X, q->radius_SPH, [&](const auto &b) {
                  return Distance(q, b) < r_short && q != b && (VectorAngle(q->interpolated_normal_SPH_original, b->X - q->X) < M_PI / 7);
               });

               /* -------------------------------------------------------------------------- */

               // capture
               auto nearWallCondition = [&](const auto &Q) {
                  return Distance(q, Q) < captureRange && q != Q && (VectorAngle(q->interpolated_normal_SPH, Q->X - q->X) < M_PI / 7);
               };
               q->isCaptured = net->BucketPoints.any_of(q->X, captureRange, nearWallCondition);
            }
         });
      }
   };

   DebugPrint("壁粒子のオブジェクト外向き法線方向を計算", Green);
   wall_p.clear();
   for (const auto &[obj, poly] : RigidBodyObject)
      for (const auto &p : obj->getPoints())
         if (p->isCaptured) {
            wall_p.emplace(p);
         }

   // \label{SPH:freeslip}
   for (const auto &p : wall_p)
   // if (p->isFirstWallLayer)
   {
      p->U_SPH.fill(0.);
      double total_w = 0;
      auto X = p->X + 2 * p->normal_SPH;
      double r = p->radius_SPH;
      net->BucketPoints.apply(X, 1.1 * r, [&](const auto &q) {
         auto w = q->volume * w_Bspline(Norm(X - q->X), r);
         p->U_SPH += q->U_SPH * w;
         total_w += w;
      });
      if (total_w == 0.)
         p->U_SPH.fill(0.);
      else
         p->U_SPH /= total_w;
      // p->U_SPH = -p->U_SPH;
      // p->U_SPH.fill(0.);
      p->U_SPH = Reflect(p->U_SPH, p->normal_SPH);  //\label{SPH:wall_particle_velocity}
      // if (p->isFirstWallLayer)
      //    p->U_SPH = Chop(p->U_SPH, p->normal_SPH);
      // p->U_SPH = Projection(Reflect(p->U_SPH, p->normal_SPH), p->normal_SPH);  //\label{SPH:wall_particle_velocity}
   }
};

/*DOC_EXTRACT SPH

## 法線方向の計算と水面の判定

*/
// #ifndef USE_AIR_PARTICLE
// #define USE_MIRROR_PARTICLE
// #define USE_ONE_AUXP
// #define USE_ALL_AUXP
// #endif
//
void setFreeSurface(auto &net, const auto &RigidBodyObject) {

   DebugPrint("水粒子のオブジェクト外向き法線方向を計算", Green);
   // refference: A. Krimi, M. Jandaghian, and A. Shakibaeinia, Water (Switzerland), vol. 12, no. 11, pp. 1–37, 2020.

   /*DOC_EXTRACT SPH

   ### 法線方向の計算

   CHECKED \ref{SPH:interpolated_normal_SPH}{単位法線ベクトル}: $`{\bf n}_i = {\rm Normalize}\left(-\sum_j {\frac{m_j}{\rho_j} \nabla W_{ij} }\right)`$

   単位法線ベクトルは，`interpolated_normal_SPH`としている．

   */

   std::vector<Network *> all_nets = {net};
   for (const auto &[obj, poly] : RigidBodyObject)
      all_nets.push_back(obj);

   auto netPoints = ToVector(net->getPoints());

   auto water = net;

#pragma omp parallel
   for (const auto &p : netPoints)
#pragma omp single nowait
   {
      // 初期化
      p->COM_SPH.fill(0.);
      p->totalMass_SPH = 0.;
      p->intp_density = 0.;
      p->intp_density_next = 0.;
      p->interpolated_normal_SPH_original.fill(0.);
      p->interpolated_normal_SPH_water.fill(0.);
      p->interpolated_normal_SPH_original_next.fill(0.);
      p->interpolated_normal_SPH_water_next.fill(0.);
      //
      p->interpolated_normal_SPH_rigid.fill(0.);
      p->interpolated_normal_SPH_rigid_next.fill(0.);
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
         p->interpolated_normal_SPH_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         p->interpolated_normal_SPH_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
         //
         // 歪度
         // p->sample_SPH += 1;
         // p->mean_SPH += q->X;
         // samples.push_back(q->X);
      });
      p->interpolated_normal_SPH_water = p->interpolated_normal_SPH_original;
      p->interpolated_normal_SPH_water_next = p->interpolated_normal_SPH_original_next;
      std::vector<Tddd> near_wall_particle, near_wall_particle_next;
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
               p->interpolated_normal_SPH_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
               p->interpolated_normal_SPH_rigid -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);

               if (Distance(X_next(p), X_next(q)) < p->radius_SPH / p->C_SML * 1.5)
                  near_wall_particle_next.emplace_back(q->normal_SPH);
               p->interpolated_normal_SPH_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
               p->interpolated_normal_SPH_rigid_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
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

      // \label{SPH:interpolated_normal_SPH}
      p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);  //\label{SPH:interpolated_normal_SPH}
      p->interpolated_normal_SPH_original_choped = p->interpolated_normal_SPH_original;
      p->interpolated_normal_SPH_original_next_choped = p->interpolated_normal_SPH_original_next;

      for (const auto &n : near_wall_particle) {
         p->interpolated_normal_SPH_original_choped = Chop(p->interpolated_normal_SPH_original_choped, n);
         p->interpolated_normal_SPH_original_next_choped = Chop(p->interpolated_normal_SPH_original_next_choped, n);
         p->interpolated_normal_SPH = Normalize(Chop(p->interpolated_normal_SPH, n));
         if (!isFinite(p->interpolated_normal_SPH))
            p->interpolated_normal_SPH = {0., 0., 1.};
      }

      p->interpolated_normal_SPH_next = Normalize(p->interpolated_normal_SPH_original_next);
      for (const auto &n : near_wall_particle_next) {
         p->interpolated_normal_SPH_next = Normalize(Chop(p->interpolated_normal_SPH_next, n));
         if (!isFinite(p->interpolated_normal_SPH_next))
            p->interpolated_normal_SPH_next = {0., 0., 1.};
      }
   }

#pragma omp parallel
   for (const auto &p : netPoints)
#pragma omp single nowait
   {
      /*DOC_EXTRACT SPH

      ### 水面の判定

      `surface_condition0,1`の両方を満たす場合，水面とする．

      */

      /* ------------------------ Neumann surface condition ----------------------- */

      const auto r_short = (p->radius_SPH / p->C_SML) * 2.;

      auto surface = net->BucketPoints.none_of(p->X, p->radius_SPH, [&](const auto &q) {
         return Distance(p, q) < r_short && p != q && (VectorAngle(p->interpolated_normal_SPH_water, q->X - p->X) < M_PI / 5);
      });

      if (surface)
         p->isNeumannSurface = std::ranges::any_of(RigidBodyObject, [&](const auto &net) {
            return std::get<0>(net)->BucketPoints.any_of(p->X, p->radius_SPH, [&](const auto &q) {
               return Distance(p, q) < r_short && p != q && (VectorAngle(p->interpolated_normal_SPH_water, q->X - p->X) < M_PI / 5);
            });
         });
      else
         p->isNeumannSurface = false;

      /* -------------------------- Surface condition ---------------------------- */

      // p->isSurface = true;

      // const auto radius = (p->radius_SPH / p->C_SML) * 2.;

      // auto surface_condition0 = [&](const auto &q) {
      //    return Distance(p, q) < radius && p != q && (VectorAngle(p->interpolated_normal_SPH_original, q->X - p->X) < M_PI / 6);
      // };

      // p->isSurface = net->BucketPoints.none_of(p->X, radius, surface_condition0);

      // auto surface_condition1 = [&](const auto &q) {
      //    // 0.8よりも小さいものは，条件を厳しくする
      //    // 0.8よりも大きいものは，条件を緩める
      //    return ((p->intp_density / _WATER_DENSITY_ < 0.8) && Distance(p, q) < radius && p != q && (VectorAngle(p->interpolated_normal_SPH_original, -q->normal_SPH) < M_PI / 3)) ||
      //           ((p->intp_density / _WATER_DENSITY_ > 0.8) && Distance(p, q) < radius * 1.2 && p != q && (VectorAngle(p->interpolated_normal_SPH_original, -q->normal_SPH) < M_PI / 2));
      // };

      // for (const auto &[obj, poly] : RigidBodyObject)
      //    if (p->isSurface) {
      //       p->isSurface = obj->BucketPoints.none_of(p->X, radius, surface_condition1);
      //    }

      // if (p->intp_density / _WATER_DENSITY_ < 0.89)
      //    p->isSurface = true;
   }

   /*DOC_EXTRACT SPH

   ## 水面補助粒子の作成

   */

   //! -------------------------------------------------------------------------- */
#pragma omp parallel
   for (const auto &A : netPoints)
#pragma omp single nowait
   {

      A->grad_corr_M.fill({0., 0., 0.});
      A->grad_corr_M_next.fill({0., 0., 0.});
      A->grad_U.fill({0., 0., 0.});
      A->Mat1.fill({0., 0., 0.});
      A->Mat2.fill({0., 0., 0.});
      A->Mat3.fill({0., 0., 0.});

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

      for (const auto &net : all_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured)
               add(B);
         });

      A->Eigenvalues_of_M = {0., 0., 0.};
      A->inv_grad_corr_M = Inverse(A->grad_corr_M);
      //
      A->inv_grad_corr_M_next = Inverse(A->grad_corr_M_next);
      //
      {
         auto [lambdas, vectors] = Eigensystem(A->grad_corr_M, 1E-4, 300.);
         A->Eigenvalues_of_M = lambdas;
         A->Eigenvectors_of_M = vectors;
      }
      {
         auto [lambdas, vectors] = Eigensystem(A->inv_grad_corr_M, 1E-4, 300.);
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
      A->isSurface = A->var_Eigenvalues_of_M > 0.125;
      //
      T3Tddd I;
      IdentityMatrix(I);
      A->Mat_B = Dot(-I, Inverse(A->Mat1 + Dot(Dot(A->Mat2, A->inv_grad_corr_M), A->Mat3)));

      /* --------------------------------------------------------- */

      {
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
   }

   //! -------------------------------------------------------------------------- */

#if defined(USE_ONE_AUXP) || defined(USE_ALL_AUXP)
   DebugPrint("水面ネットワークの初期化", Green);
   if (net->surfaceNet == nullptr)
      net->surfaceNet = new Network();

   DebugPrint("水面補助粒子の作成", Green);
   for (const auto &p : net->getPoints()) {
      p->isAuxiliary = false;
      if (p->isSurface) {

         int num = 2, count = 0;
         // if (Norm(p->interpolated_normal_SPH_rigid) > 1E-10)
         //    num = 2;
         // else
         //    num = 1;

         for (auto &auxp : p->auxiliaryPoints) {
            if (count++ < num) {
               if (auxp == nullptr) {
                  auxp = new networkPoint(net->surfaceNet, p->X);
                  auxp->isAuxiliary = true;
               }
            } else if (auxp != nullptr) {
               delete auxp;
               auxp = nullptr;
            }
         }
      } else
         for (auto &auxp : p->auxiliaryPoints) {
            if (auxp != nullptr) {
               delete auxp;
               auxp = nullptr;
            }
         }
   }

   // DebugPrint("水面ネットワークの初期化", Green);
   // if (net->surfaceNet != nullptr)
   //    delete net->surfaceNet;
   // net->surfaceNet = new Network();

   // DebugPrint("水面補助粒子の作成", Green);
   // for (const auto &p : net->getPoints()) {
   //    p->isAuxiliary = false;
   //    p->auxiliaryPoints.fill(nullptr);
   //    if (p->isSurface)
   //       for (auto &auxp : p->auxiliaryPoints)
   //          auxp = new networkPoint(net->surfaceNet, {0., 0., 0.});
   // }

   DebugPrint("水面補助粒子の作成", Green);

   #pragma omp parallel
   for (const auto &p : net->getPoints())
   #pragma omp single nowait
   {
      double d = 0;
      if (p->isSurface) {
         int count = 0;
         for (auto &auxp : p->auxiliaryPoints)
            if (auxp != nullptr)
               count++;
         for (auto i = 0; i < p->auxiliaryPoints.size(); i++) {
            auto &auxp = p->auxiliaryPoints[i];
            if (auxp != nullptr) {
               auto radius_SPH = p->radius_SPH;
               auto C_SML = p->C_SML;
               auxp->radius_SPH = radius_SPH;
               auxp->C_SML = C_SML;
               auxp->surfacePoint = p;
               auxp->isAuxiliary = true;
               auxp->isCaptured = true;
               auxp->isSurface = false;
               auxp->p_SPH = p->p_SPH;
               auxp->U_SPH = p->U_SPH;
               auxp->b_vector = p->b_vector;
               auxp->DUDt_SPH = p->DUDt_SPH;
               auxp->volume = p->volume;
               auxp->setXSingle(p->X);                                                                                        // 初期値
               auxp->setXSingle(p->X + (i + 1) * p->radius_SPH / p->C_SML * Normalize(p->interpolated_normal_SPH_original));  // 初期値
               auxp->setDensityVolume(_WATER_DENSITY_, p->volume);
               //
               auxp->volume_next = auxp->volume;
               auxp->X_next = X_next(p) + (i + 1) * p->radius_SPH / p->C_SML * Normalize(p->interpolated_normal_SPH_original_next);
               auxp->mass_next = _WATER_DENSITY_ * auxp->volume;
   #if defined(USE_RungeKutta)
               auxp->RK_U = p->RK_U;
               auxp->RK_X = p->RK_X;
               auxp->RK_P = p->RK_P;
               auxp->RK_rho = p->RK_rho;
   #elif defined(USE_LeapFrog)
               auxp->LPFG_X = p->LPFG_X;
               auxp->LPFG_rho = p->LPFG_rho;
   #endif
               auto opt_func = [&](const auto &p,
                                   const std::array<double, 3> &p_X,
                                   const auto &grad_to_minimize_more,
                                   double density,
                                   const std::array<double, 3> &q_X1,
                                   const std::array<double, 3> &q_X2,
                                   const double &vol1,
                                   const double &vol2) {
                  auto rho = _WATER_DENSITY_;
                  auto f = density + rho * vol1 * w_Bspline(Norm(p_X - q_X1), p->radius_SPH) + rho * vol2 * w_Bspline(Norm(p_X - q_X2), p->radius_SPH) - rho;
                  auto F2 = grad_to_minimize_more - rho * vol1 * grad_w_Bspline(p_X, q_X1, p->radius_SPH) - rho * vol2 * grad_w_Bspline(p_X, q_X2, p->radius_SPH);
                  return Dot(F2, F2) + f * f;
               };

               auto get_X_V = [&](const auto &center,
                                  const auto &first_dir,
                                  const auto &grad_to_minimize,
                                  double intp_density) {
                  const double rho = _WATER_DENSITY_;
                  double min_f = 1E+20, best_vol, best_rho, vol = p->volume, R, f, opt_value, best_vol1, best_vol2, vol1, vol2;
                  int N = 100;
                  std::array<double, 3> X1, X2, F, best_X1, best_X2, unit_normal1 = Normalize(first_dir), unit_normal2, grad_to_minimize_more;
                  std::array<double, 3> C1, C2;
                  // const double a = 10.;
                  auto V_init = vol;
                  for (auto k = 1; k < 50; ++k) {
                     vol2 = vol1 = vol = 5. * V_init * (double)k / (double)N;
                     // if (k == 0)
                     // C1 = center + 0.5 * p->radius_SPH * unit_normal1;
                     // else
                     //    C1 = best_X1;
                     for (auto i = 1; i < N; ++i) {
                        X1 = center + p->radius_SPH * (0.2 + 0.7 * (double)i / (double)N) * unit_normal1;
                        // X1 = C1 + std::pow(a, -k) * 0.5 * p->radius_SPH * (-1. + 2. * (double)i / (double)N) * unit_normal1;
                        grad_to_minimize_more = grad_to_minimize - rho * vol1 * grad_w_Bspline(center, X1, p->radius_SPH);
                        unit_normal2 = Normalize(grad_to_minimize_more);
                        // C2 = center + 0.5 * p->radius_SPH * unit_normal2;
                        for (auto j = 1; j < N; ++j) {
                           X2 = center + p->radius_SPH * (0.3 + 0.7 * (double)j / (double)N) * unit_normal2;
                           // X2 = C2 + std::pow(a, -k) * 0.5 * p->radius_SPH * (-1. + 2. * (double)j / (double)N) * unit_normal2;
                           if ((opt_value = opt_func(p, center, grad_to_minimize_more, intp_density, X1, X2, vol, vol)) < min_f) {
                              min_f = opt_value;
                              best_vol2 = best_vol1 = best_vol = vol;
                              best_rho = rho;
                              best_X1 = X1;
                              best_X2 = X2;
                           }
                        }
                     }
                  }

                  // N = 300;
                  // for (auto i = 1; i < N; ++i) {
                  //    vol1 = 5. * best_vol * (double)i / (double)N;
                  //    grad_to_minimize_more = grad_to_minimize - rho * vol1 * grad_w_Bspline(center, best_X1, p->radius_SPH);
                  //    for (auto j = 1; j < N; ++j) {
                  //       vol2 = 5. * best_vol * (double)j / (double)N;
                  //       if ((opt_value = opt_func(p, center, grad_to_minimize_more, intp_density, best_X1, best_X2, vol1, vol2)) < min_f) {
                  //          min_f = opt_value;
                  //          best_vol1 = vol1;
                  //          best_vol2 = vol2;
                  //       }
                  //    }
                  // }
                  return std::tuple<Tddd, Tddd, double, double, double>{best_X1, best_X2, best_vol1, best_vol2, best_rho};
               };

               // {
               //    auto [X1, X2, V1, V2, RHO] = get_X_V(p->X,
               //                                         p->interpolated_normal_SPH_original_choped /*base direction*/,
               //                                         p->interpolated_normal_SPH_original /*minimize this by adding*/,
               //                                         p->intp_density);
               //    p->auxiliaryPoints[0]->setDensityVolume(RHO, V1);
               //    p->auxiliaryPoints[1]->setDensityVolume(RHO, V2);
               //    p->auxiliaryPoints[0]->setXSingle(X1);  // 初期値
               //    p->auxiliaryPoints[1]->setXSingle(X2);  // 初期値
               // }
               // {
               //    auto [X1, X2, V1, V2, RHO] = get_X_V(X_next(p),
               //                                         p->interpolated_normal_SPH_original_next_choped,
               //                                         p->interpolated_normal_SPH_original_next,
               //                                         p->intp_density_next);
               //    p->auxiliaryPoints[0]->volume_next = V1;
               //    p->auxiliaryPoints[1]->volume_next = V2;
               //    p->auxiliaryPoints[0]->X_next = X1;
               //    p->auxiliaryPoints[1]->X_next = X2;
               //    p->auxiliaryPoints[0]->mass_next = RHO * V1;
               //    p->auxiliaryPoints[1]->mass_next = RHO * V2;
               // }
            }
         }
      }
   }

      // #pragma omp parallel
      //    for (const auto &p : net->getPoints())
      // #pragma omp single nowait
      //    {
      //       if (p->isSurface) {
      //          for (auto &AUX : p->auxiliaryPoints) {
      //             if (AUX != nullptr)
      //                AUX->NR_double.initialize(AUX->volume);
      //          }
      //       }
      //    }
      //
      //    for (auto i = 0; i <= 100; ++i) {
      // #pragma omp parallel
      //       for (const auto &p : net->getPoints())
      // #pragma omp single nowait
      //       {
      //          if (p->isSurface) {
      //             for (auto &AUX : p->auxiliaryPoints) {
      //                auto F1 = p->interpolated_normal_SPH_original;
      //                auto F2 = p->intp_density - _WATER_DENSITY_;
      //                auto dF1dx = -AUX->rho * grad_w_Bspline(p->X, AUX->X, p->radius_SPH);
      //                auto dF2dx = AUX->rho * w_Bspline(Norm(p->X - AUX->X), p->radius_SPH);
      //                for (const auto &net : all_nets)
      //                   net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
      //                      if (q->isSurface)
      //                         for (const auto &aux : q->auxiliaryPoints)
      //                            if (aux != nullptr) {
      //                               F1 -= aux->rho * aux->volume * grad_w_Bspline(p->X, aux->X, p->radius_SPH);
      //                               F2 += aux->rho * aux->volume * w_Bspline(Norm(p->X - aux->X), p->radius_SPH);
      //                            }
      //                   });
      //                auto F = Dot(F1, F1) / 2. + F2 * F2 / 2.;
      //                auto dFdx = Dot(dF1dx, F1) + dF2dx * F2;
      //                AUX->NR_double.update(F, dFdx, 0.001);
      //             }
      //          }
      //       }
      // #pragma omp parallel
      //       for (const auto &p : net->getPoints())
      // #pragma omp single nowait
      //          if (p->isSurface) {
      //             for (auto &AUX : p->auxiliaryPoints) {
      //                AUX->volume_next = AUX->volume = AUX->NR_double.X;
      //                AUX->setDensityVolume(_WATER_DENSITY_, AUX->volume);
      //                AUX->mass_next = _WATER_DENSITY_ * AUX->volume_next;
      //             }
      //          }
      //    }

      /* -------------------------------------------------------------------------- */

   #pragma omp parallel
   for (const auto &p : net->getPoints())
   #pragma omp single nowait
   {
      // 初期化
      p->totalMass_SPH = 0.;
      p->interpolated_normal_SPH_original_modified.fill(0.);
      //
      networkPoint *closest_surface_point = nullptr;
      double min_distance = 1e10, distance;
      //
      for (const auto &net : all_nets)
         net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            // w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
            if (Distance(p, q) < p->radius_SPH) {
               auto w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
               p->totalMass_SPH += q->mass;
            }
            p->interpolated_normal_SPH_original_modified -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            if (q->isSurface) {
   #if defined(USE_ONE_AUXP)
               if ((distance = Distance(p->X, q->X)) < min_distance) {
                  min_distance = distance;
                  closest_surface_point = q;
               }
   #elif defined(USE_ALL_AUXP)
            for (const auto &AUX : q->auxiliaryPoints)
               if (AUX != nullptr) p->interpolated_normal_SPH_original_modified -= AUX->rho * AUX->volume * grad_w_Bspline(p->X, AUX->X, p->radius_SPH);
   #endif
            }
         });

   #if defined(USE_ONE_AUXP)
      // if (p->isSurface)
      if (closest_surface_point != nullptr)
         for (const auto &AUX : closest_surface_point->auxiliaryPoints)
            if (AUX != nullptr) {
               auto q = AUX;
               p->interpolated_normal_SPH_original_modified -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            }
   #endif
   }
   net->surfaceNet->setGeometricProperties();

#endif

   DebugPrint("水粒子のオブジェクト外向き法線方向を計算 done", Green);
};

#endif