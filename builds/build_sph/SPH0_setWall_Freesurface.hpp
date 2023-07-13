#ifndef SPH_setWall_Freesurface_H
#define SPH_setWall_Freesurface_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

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
         p->isFluid = false;
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
      // ここでも結構変わる
      const double captureRange = p->radius_SPH;  //\label{SPH:capture_condition_1st}

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
                  return Distance(q, Q) < particle_spacing * 1.5 && q != Q && (VectorAngle(q->interpolated_normal_SPH, Q->X - q->X) < M_PI / 8);
               };
               q->isFirstWallLayer = net->BucketPoints.any_of(q->X, q->radius_SPH, firstWallLayerCondition);

               // capture
               auto nearWallCondition = [&](const auto &Q) {
                  return Distance(q, Q) < captureRange && q != Q;  // && (VectorAngle(q->interpolated_normal_SPH, Q->X - q->X) < M_PI / 6);
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
      if (p->isFirstWallLayer) {
         p->U_SPH.fill(0.);
         double total_w = 0;
         auto X = p->X + 2 * p->normal_SPH;
         net->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
            if (Distance(X, q) < p->radius_SPH) {
               auto w = q->volume * w_Bspline(Norm(X - q->X), p->radius_SPH);
               p->U_SPH += q->U_SPH * w;
               total_w += w;
            }
         });
         if (total_w == 0.)
            p->U_SPH.fill(0.);
         else
            p->U_SPH /= total_w;
         // p->U_SPH = Reflect(p->U_SPH, p->normal_SPH);  //\label{SPH:wall_particle_velocity}
         p->U_SPH = Projection(Reflect(p->U_SPH, p->normal_SPH), p->normal_SPH);  //\label{SPH:wall_particle_velocity}
      }
};

std::tuple<double, std::array<double, 3>> interpDensityAndGrad(const auto &p, const auto &all_nets) {
   // 法線方向ではない．単なる密度と密度勾配の補間
   std::array<double, 3> grad_rho = {0., 0., 0.};
   double rho = 0.;
   std::vector<Tddd> near_wall_particle;
   for (const auto &net : all_nets)
      net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         if (Distance(p, q) < p->radius_SPH) {
            grad_rho += q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            // grad_rho += p->rho * q->mass * (q->rho / (q->rho * q->rho) + p->rho / (p->rho * p->rho)) * grad_w_Bspline(p->X, q->X, p->radius_SPH);  //\label{SPH:gradP1}
            // grad_rho += (q->rho - p->rho) * q->mass / p->rho * grad_w_Bspline(p->X, q->X, p->radius_SPH);  //\label{SPH:gradP2}
            rho += q->rho * q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
         }

         if (Distance(p, q) < p->radius_SPH / p->C_SML)
            near_wall_particle.emplace_back(q->normal_SPH);
         // このようにすることで，自身の補助点を無視して欠損のある方向を指す．
         // if (p != q && q->isSurface)
         //    for (const auto &Q : q->auxiliaryPoints)
         //       if (Q != nullptr) {
         //          grad_rho += q->rho * Q->volume * grad_w_Bspline(p->X, Q->X, p->radius_SPH);
         //          rho += q->rho * Q->volume * w_Bspline(Norm(p->X - Q->X), p->radius_SPH);
         //       }
      });
   // for (const auto &n : near_wall_particle)
   //    grad_rho = Chop(grad_rho, n);

   return {rho, grad_rho};
};

/*DOC_EXTRACT SPH

## 法線方向の計算と水面の判定

*/

// #define USE_SIMPLE_SINGLE_AUX
#define USE_SHARED_AUX

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

#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      // 初期化
      p->COM_SPH.fill(0.);
      p->totalMass_SPH = 0.;
      p->intp_density = 0.;
      p->interpolated_normal_SPH_original.fill(0.);
      p->interpolated_normal_SPH_original_next.fill(0.);
      double w;

      net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         // w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
         if (Distance(p, q) < p->radius_SPH) {
            auto w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
            p->intp_density += q->rho * w;
            p->COM_SPH += q->mass * q->X;
            p->totalMass_SPH += q->mass;
         }
         p->interpolated_normal_SPH_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         p->interpolated_normal_SPH_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
      });

      std::vector<Tddd> near_wall_particle, near_wall_particle_next;
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            if (q->isCaptured) {
               // w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
               if (Distance(p, q) < p->radius_SPH) {
                  auto w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
                  p->intp_density += q->rho * w;
                  p->COM_SPH += q->mass * q->X;
                  p->totalMass_SPH += q->mass;
               }

               if (Distance(p, q) < p->radius_SPH / p->C_SML * 1.5)
                  near_wall_particle.emplace_back(q->normal_SPH);
               p->interpolated_normal_SPH_original -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);

               if (Distance(X_next(p), X_next(q)) < p->radius_SPH / p->C_SML * 1.5)
                  near_wall_particle_next.emplace_back(q->normal_SPH);
               p->interpolated_normal_SPH_original_next -= q->rho * V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
            }
         });

      p->COM_SPH /= p->totalMass_SPH;

      // \label{SPH:interpolated_normal_SPH}
      p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);  //\label{SPH:interpolated_normal_SPH}
      p->interpolated_normal_SPH_original_choped = p->interpolated_normal_SPH_original;

      for (const auto &n : near_wall_particle) {
         p->interpolated_normal_SPH_original_choped = Chop(p->interpolated_normal_SPH_original_choped, n);
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
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      /*DOC_EXTRACT SPH

      ### 水面の判定

      `surface_condition0,1`の両方を満たす場合，水面とする．

      */

      p->isSurface = true;
      const auto radius = (p->radius_SPH / p->C_SML) * 3.;

      auto surface_condition0 = [&](const auto &q) {
         return Distance(p, q) < radius && p != q && (VectorAngle(p->interpolated_normal_SPH, q->X - p->X) < std::numbers::pi / 4);
      };

      auto surface_condition1 = [&](const auto &q) {
         return Distance(p, q) < radius && p != q && (VectorAngle(p->interpolated_normal_SPH, -q->normal_SPH) < std::numbers::pi / 180. * 45);
      };

      if (net->BucketPoints.any_of(p->X, radius, surface_condition0))
         p->isSurface = false;

      if (p->isSurface)
         for (const auto &[obj, poly] : RigidBodyObject)
            if (obj->BucketPoints.any_of(p->X, radius, surface_condition1))
               p->isSurface = false;
   }

   /*DOC_EXTRACT SPH

   ## 水面補助粒子の作成

   */

   DebugPrint("水面ネットワークの初期化", Green);
   if (net->surfaceNet == nullptr)
      net->surfaceNet = new Network();

   DebugPrint("水面補助粒子の作成", Green);
   for (const auto &p : net->getPoints()) {
      p->isAuxiliary = false;
      if (p->isSurface) {
         for (auto &auxp : p->auxiliaryPoints) {
            if (auxp == nullptr) {
               auxp = new networkPoint(net->surfaceNet, {0., 0., 0.});
               auxp->isAuxiliary = true;
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
         for (auto &auxp : p->auxiliaryPoints) {
            auto radius_SPH = p->radius_SPH;
            auto C_SML = p->C_SML;
            d += radius_SPH / C_SML;
            auxp->setXSingle(aux_position(p, d));
            auxp->radius_SPH = radius_SPH;
            auxp->C_SML = C_SML;
            auxp->surfacePoint = p;
            auxp->isAuxiliary = true;
            auxp->isCaptured = true;
            auxp->isSurface = false;
            auxp->p_SPH = p->p_SPH;
            auxp->U_SPH = p->U_SPH;
            auxp->b_vector.fill(0.);
            auxp->DUDt_SPH.fill(0.);
            auxp->volume = d * d * d;
            auxp->setDensityVolume(_WATER_DENSITY_, p->volume);
#if defined(USE_RungeKutta)
            auxp->RK_U = p->RK_U;
            auxp->RK_X = p->RK_X;
            auxp->RK_P = p->RK_P;
            auxp->RK_rho = p->RK_rho;
#elif defined(USE_LeapFrog)
            auxp->LPFG_X = p->LPFG_X;
            auxp->LPFG_rho = p->LPFG_rho;
#endif
            /* -------------------------------------------------------------------------- */
            auto VEC = -(p->COM_SPH - p->X);
            auxp->setXSingle(p->X + VEC);  // 初期値
            // auxp->setXSingle((p->X - p->radius_SPH / 2. * Normalize(p->COM_SPH - p->X)));  // 初期値
            // 重心位置を使う場合
            auxp->volume = std::pow(p->radius_SPH / p->C_SML, 3);
            auxp->setDensityVolume(_WATER_DENSITY_, auxp->volume);
            auto r0 = std::pow(p->volume, 1. / 3.) / 2.;
            // optimizeFunction(auxp->X, auxp->volume, p, auxp->radius_SPH, auxp->C_SML);
            // auxp->volume = (_WATER_DENSITY_ - p->intp_density) / (_WATER_DENSITY_ * w_Bspline(Norm(auxp->X - p->X), p->radius_SPH));
            auxp->setDensityVolume(_WATER_DENSITY_, auxp->volume);
            double min_f = 1E+20, best_vol, best_dist, best_rho;
            double vol, rho, dist, R;
            std::array<double, 3> X;
            // auto unit_normal = Normalize(p->interpolated_normal_SPH);
            // auto n_vec = (p->interpolated_normal_SPH_original + p->interpolated_normal_SPH_original_choped) / 2.;
            auto n_vec = p->interpolated_normal_SPH_original;
            auto unit_normal = Normalize(n_vec);
            int N = 1000;
            for (auto i = 1; i < N; ++i) {
               R = p->radius_SPH * i / (double)N;
               dist = R + r0;
               X = p->X + dist * unit_normal;
               for (auto j = 1; j < N; ++j) {
                  rho = _WATER_DENSITY_;
                  vol = std::pow(2. * R, 3.) * (2 * j / (double)N);
                  //
                  // auto [intp_density, intp_grad_density] = interpDensityAndGrad(p, all_nets);
                  // auto f = intp_density + auxp->rho * vol * w_Bspline(dist, p->radius_SPH) - p->rho;
                  // auto F = intp_grad_density + auxp->rho * vol * grad_w_Bspline(p->X, X, p->radius_SPH);
                  //
                  // auto [intp_density, intp_grad_density] = interpDensityAndGrad(p, all_nets);
                  auto f = p->intp_density + rho * vol * w_Bspline(dist, p->radius_SPH) - p->rho;
                  auto F = n_vec - rho * vol * grad_w_Bspline(p->X, X, p->radius_SPH);
                  // auto F = p->interpolated_normal_SPH_original_choped - auxp->rho * vol * grad_w_Bspline(p->X, X, p->radius_SPH);
                  //
                  f /= _WATER_DENSITY_;
                  F /= _WATER_DENSITY_;
                  auto opt_func = 0;
                  opt_func += f * f;
                  opt_func += Dot(F, F);
                  if (opt_func < min_f) {
                     min_f = opt_func;
                     best_vol = vol;
                     best_dist = dist;
                     best_rho = rho;
                     // std::cout << "auxp = " << auxp
                     //           << ", Dot(F, F) = " << Dot(F, F)
                     //           << ", f*f = " << f * f
                     //           << ", opt_func = " << opt_func
                     //           << ", dist = " << dist
                     //           << ", r0 = " << r0 << std::endl;
                  }
               }
            }

            auxp->setDensityVolume(best_rho, best_vol);
            auxp->setXSingle(p->X + best_dist * unit_normal);  // 初期値
         }
      }
   }

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
            if (q->isSurface) {
               if ((distance = Distance(p->X, q->X)) < min_distance) {
                  min_distance = distance;
                  closest_surface_point = q;
               }
            }

            p->interpolated_normal_SPH_original_modified -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         });

      if (closest_surface_point != nullptr)
         for (const auto &AUX : closest_surface_point->auxiliaryPoints) {
            auto q = AUX;
            p->interpolated_normal_SPH_original_modified -= q->rho * q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         }
   }

   net->surfaceNet->setGeometricProperties();
};

#endif