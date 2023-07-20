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
                  return Distance(q, Q) < particle_spacing * 1.8 && q != Q && (VectorAngle(q->interpolated_normal_SPH, Q->X - q->X) < M_PI / 6);
               };
               q->isFirstWallLayer = net->BucketPoints.any_of(q->X, q->radius_SPH, firstWallLayerCondition);

               // capture
               auto nearWallCondition = [&](const auto &Q) {
                  return Distance(q, Q) < captureRange && q != Q /* && (VectorAngle(q->interpolated_normal_SPH, Q->X - q->X) < M_PI / 6)*/;
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
      net->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
         {
            auto w = q->volume * w_Bspline(Norm(X - q->X), p->radius_SPH);
            p->U_SPH += q->U_SPH * w;
            total_w += w;
         }
      });
      if (total_w == 0.)
         p->U_SPH.fill(0.);
      else
         p->U_SPH /= total_w;
      // p->U_SPH = -p->U_SPH;
      p->U_SPH = Reflect(p->U_SPH, p->normal_SPH);  //\label{SPH:wall_particle_velocity}
      // p->U_SPH = Projection(Reflect(p->U_SPH, p->normal_SPH), p->normal_SPH);  //\label{SPH:wall_particle_velocity}
   }
};

/*DOC_EXTRACT SPH

## 法線方向の計算と水面の判定

*/

#define USE_ONE_AUXP
// #define USE_ALL_AUXP

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

      std::vector<networkPoint *> samples;
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

         int num = 2, count = 0;
         // if (Norm(p->interpolated_normal_SPH_rigid) > 1E-10)
         //    num = 2;
         // else
         //    num = 1;

         for (auto &auxp : p->auxiliaryPoints) {
            if (count++ < num) {
               if (auxp == nullptr) {
                  auxp = new networkPoint(net->surfaceNet, {0., 0., 0.});
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
               auxp->volume = std::pow(p->radius_SPH / p->C_SML, 3);
               auxp->setXSingle(p->X + (i + 1) * p->radius_SPH / p->C_SML * Normalize(p->interpolated_normal_SPH_original));  // 初期値
               auxp->setDensityVolume(_WATER_DENSITY_, auxp->volume);
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

               {
                  auto [X1, X2, V1, V2, RHO] = get_X_V(p->X,
                                                       p->interpolated_normal_SPH_original_choped /*base direction*/,
                                                       p->interpolated_normal_SPH_original /*minimize this by adding*/,
                                                       p->intp_density);
                  p->auxiliaryPoints[0]->setDensityVolume(RHO, V1);
                  p->auxiliaryPoints[1]->setDensityVolume(RHO, V2);
                  p->auxiliaryPoints[0]->setXSingle(X1);  // 初期値
                  p->auxiliaryPoints[1]->setXSingle(X2);  // 初期値
               }
               {
                  auto [X1, X2, V1, V2, RHO] = get_X_V(X_next(p),
                                                       p->interpolated_normal_SPH_original_next_choped,
                                                       p->interpolated_normal_SPH_original_next,
                                                       p->intp_density_next);
                  p->auxiliaryPoints[0]->volume_next = V1;
                  p->auxiliaryPoints[1]->volume_next = V2;
                  p->auxiliaryPoints[0]->X_next = X1;
                  p->auxiliaryPoints[1]->X_next = X2;
                  p->auxiliaryPoints[0]->mass_next = RHO * V1;
                  p->auxiliaryPoints[1]->mass_next = RHO * V2;
               }
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

   // #pragma omp parallel
   //    for (const auto &p : net->getPoints())
   // #pragma omp single nowait
   //    {
   //       if (p->isSurface) {
   //          for (auto &AUX : p->auxiliaryPoints) {
   //             if (AUX != nullptr)
   //                AUX->NR_double.initialize(AUX->volume_next);
   //          }
   //       }
   //    }

   //    for (auto i = 0; i <= 100; ++i) {
   // #pragma omp parallel
   //       for (const auto &p : net->getPoints())
   // #pragma omp single nowait
   //       {
   //          if (p->isSurface) {
   //             for (auto &AUX : p->auxiliaryPoints) {
   //                auto F1 = p->interpolated_normal_SPH_original_next;
   //                auto F2 = p->intp_density_next;
   //                auto dF1dx = -rho_next(AUX) * grad_w_Bspline(X_next(p), X_next(AUX), p->radius_SPH);
   //                auto dF2dx = rho_next(AUX) * w_Bspline(Norm(X_next(p) - X_next(AUX)), p->radius_SPH);
   //                for (const auto &net : all_nets)
   //                   net->BucketPoints.apply(p->X, p->radius_SPH * 1.1, [&](const auto &q) {
   //                      if (q->isSurface)
   //                         for (const auto &aux : q->auxiliaryPoints)
   //                            if (aux != nullptr) {
   //                               F1 -= rho_next(aux) * V_next(aux) * grad_w_Bspline(X_next(p), X_next(aux), p->radius_SPH);
   //                               F2 += rho_next(aux) * V_next(aux) * w_Bspline(Norm(X_next(p) - X_next(aux)), p->radius_SPH);
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
   //                AUX->volume_next = AUX->NR_double.X;
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
};

#endif