#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"

template <typename T>
struct InterpolationLagrange {
   std::vector<double> abscissas;
   std::vector<T> values;
   std::vector<double> denominotor;

   InterpolationLagrange(const std::vector<double> abscissas)
       : abscissas(abscissas){};

   InterpolationLagrange(const std::vector<double> abscissas, const std::vector<T> values)
       : abscissas(abscissas), values(values) {
      if (abscissas.size() != values.size()) {
         throw std::invalid_argument("Size of abscissas and values vectors must be the same");
      }
      this->set();
   };

   void set() {
      denominotor.resize(abscissas.size(), 1.);
      for (auto i = 0; i < abscissas.size(); ++i)
         for (auto j = 0; j < abscissas.size(); ++j)
            if (i != j)
               denominotor[i] *= (abscissas[i] - abscissas[j]);
   };

   T operator()(const double x) {
      T ret, N = 1;
      ret *= 0.;
      for (auto i = 0; i < abscissas.size(); ++i) {
         for (auto j = 0; j < abscissas.size(); ++j) {
            if (i != j)
               N *= (x - this->abscissas[j]) / (this->abscissas[i] - this->abscissas[j]);
         }
         ret += N * this->values[i];
         N = 1;
      }
      return ret;
   };

   T D(const double x) {
      T ret;
      ret *= 0.;
      for (auto i = 0; i < abscissas.size(); ++i) {
         for (auto j = 0; j < abscissas.size(); ++j) {
            if (i != j) {
               T temp = this->values[i] / (this->abscissas[i] - this->abscissas[j]);
               for (auto k = 0; k < abscissas.size(); ++k) {
                  if (k != i && k != j) {
                     temp *= (x - this->abscissas[k]) / (this->abscissas[i] - this->abscissas[k]);
                  }
               }
               ret += temp;
            }
         }
      }
      return ret;
   }

   std::vector<T> N(const double x) {
      std::vector<T> ret(abscissas.size(), 1.);
      for (auto i = 0; i < abscissas.size(); ++i) {
         for (auto j = 0; j < abscissas.size(); ++j)
            if (i != j)
               ret[i] *= (x - this->abscissas[j]) / (this->abscissas[i] - this->abscissas[j]);
         // ret[i] *= this->values[i];
      }
      return ret;
   };

   std::vector<T> DN(const double x) {
      std::vector<T> ret(abscissas.size(), 0.);
      for (auto i = 0; i < abscissas.size(); ++i) {
         for (auto j = 0; j < abscissas.size(); ++j) {
            if (i != j) {
               T temp = 1. / (this->abscissas[i] - this->abscissas[j]);
               for (auto k = 0; k < abscissas.size(); ++k) {
                  if (k != i && k != j) {
                     temp *= (x - this->abscissas[k]) / (this->abscissas[i] - this->abscissas[k]);
                  }
               }
               ret[i] += temp;
            }
         }
      }
      return ret;
   };
};

networkPoint *getClosestExcludeRigidBody(networkPoint *p, auto &target_nets) {
   double distance = 1E+20;
   networkPoint *P = nullptr;
   for (const auto &obj : target_nets)
      if (!obj->isRigidBody) {
         obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            auto tmp = Distance(p, q);
            if (distance > tmp) {
               distance = tmp;
               P = q;
            }
         });
      }
   return P;
};

/*DOC_EXTRACT SPH

### CFL条件の設定

$`\max({\bf u}) \Delta t \leq c_{v} h \cap \max({\bf a}) \Delta t^2 \leq c_{a} h`$
を満たすように，毎時刻$`\Delta t`$を設定する．

*/

double dt_CFL(const double dt_IN, const auto &net, const auto &RigidBodyObject) {
   double dt = dt_IN;
   const auto C_CFL_velocity = 0.02;  // dt = C_CFL_velocity*h/Max(U)
   const auto C_CFL_accel = 0.1;      // dt = C_CFL_accel*sqrt(h/Max(A))
   for (const auto &p : net->getPoints()) {
      // 速度に関するCFL条件
      auto dt_C_CFL = [&](const auto &q) {
         if (p != q) {
            auto pq = Normalize(p->X - q->X);
            auto distance = Distance(p, q);
            /* ------------------------------------------------ */
            // 相対速度
            double max_dt_vel = C_CFL_velocity * distance / std::abs(Dot(p->U_SPH - q->U_SPH, pq));
            // double max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH - q->U_SPH);
            if (dt > max_dt_vel && isFinite(max_dt_vel))
               dt = max_dt_vel;
            // 絶対速度
            max_dt_vel = C_CFL_velocity * distance / Norm(p->U_SPH);
            if (dt > max_dt_vel && isFinite(max_dt_vel))
               dt = max_dt_vel;
            /* ------------------------------------------------ */
            // 相対速度
            double max_dt_acc = C_CFL_accel * std::sqrt(distance / std::abs(Dot(p->DUDt_SPH - q->DUDt_SPH, pq)));
            // double max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH - q->DUDt_SPH));
            if (dt > max_dt_acc && isFinite(max_dt_acc))
               dt = max_dt_acc;
            // 絶対速度
            max_dt_acc = C_CFL_accel * std::sqrt(distance / Norm(p->DUDt_SPH));
            if (dt > max_dt_acc && isFinite(max_dt_acc))
               dt = max_dt_acc;
         }
      };
      net->BucketPoints.apply(p->X, p->radius_SPH, dt_C_CFL);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->radius_SPH, dt_C_CFL);
      double max_dt_vel = C_CFL_velocity * (p->radius_SPH / p->C_SML) / Norm(p->U_SPH);
      if (dt > max_dt_vel && isFinite(max_dt_vel))
         dt = max_dt_vel;
      double max_dt_acc = C_CFL_accel * std::sqrt((p->radius_SPH / p->C_SML) / Norm(p->DUDt_SPH));
      if (dt > max_dt_acc && isFinite(max_dt_acc))
         dt = max_dt_acc;
   }
   return dt;
}

#define Morikawa2019

/* -------------------------------------------------------------------------- */
Tddd aux_position(const networkPoint *p) {
   auto c = p->radius_SPH / p->C_SML;
   return p->X + c * Normalize(p->interpolated_normal_SPH);
};

Tddd aux_position_next(const networkPoint *p) {
   auto q = p->surfacePoint;
   auto c = q->radius_SPH / q->C_SML;
#if defined(USE_RungeKutta)
   return q->RK_X.getX(q->U_SPH) + c * Normalize(q->interpolated_normal_SPH_next);
#elif defined(USE_LeapFrog)
   return q->LPFG_X.get_x(q->U_SPH) + c * Normalize(q->interpolated_normal_SPH_next);
#endif
};

// \label{SPH:rho_next}
double rho_next_(auto p) {
   if (p->isAuxiliary)
      p = p->surfacePoint;
   if (p->getNetwork()->isRigidBody)
      return _WATER_DENSITY_;
   else {
#if defined(USE_RungeKutta)
      return p->RK_rho.getX(p->DrhoDt_SPH);
#elif defined(USE_LeapFrog)
      return p->LPFG_rho.get_x(p->DrhoDt_SPH);
#endif
   }
};

// \label{SPH:volume_next}
double V_next_(const auto &p) {
   return p->mass / rho_next(p);
};

// \label{SPH:position_next}
std::array<double, 3> X_next_(const auto &p) {
   if (p->isAuxiliary)
      return aux_position_next(p);
   else if (p->getNetwork()->isRigidBody)
      return p->X;
   else
#if defined(USE_RungeKutta)
      return p->RK_X.getX(p->U_SPH);
#elif defined(USE_LeapFrog)
      return p->LPFG_X.get_x(p->U_SPH);
         // return p->X + p->U_SPH * p->LPFG_X.get_dt();
#endif
};

/* -------------------------------------------------------------------------- */

// \label{SPH:rho_next}
double rho_next(auto p) {
   // return rho_next_(p);
   // return rho_next_(p);
   return _WATER_DENSITY_;
};

// \label{SPH:volume_next}
double V_next(const auto &p) {
   return p->mass / rho_next(p);
};

// \label{SPH:position_next}
std::array<double, 3> X_next(const auto &p) {
   return p->X;
   // return X_next_(p);
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT SPH

## 壁面粒子の流速と圧力

壁粒子の流速を流体粒子の流速に応じて変化させるとプログラムが煩雑になるので，**ここでは**壁面粒子の流速は常にゼロに設定することにする．
壁粒子の圧力は，水が圧縮しないように各ステップ毎に計算し直す必要がある．

*/

void setWall(const auto &net, const auto &RigidBodyObject, const auto &particle_spacing, auto &wall_p) {
   // wall_as_fluid.clear();
   wall_p.clear();

   // 初期化
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
      }
   DebugPrint("関連する壁粒子をマーク", Green);
// capture wall particles
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
               q->isCaptured = true;
               q->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3.));

               p->interpolated_normal_SPH_original.fill(0.);
               // if (Distance(p, q) < captureRange_wall_as_fluid)  //\label{SPH:select_wall_as_fluid}
               //    q->isFluid = true;
            }
         });
      }
   };

   DebugPrint("壁粒子のオブジェクト外向き法線方向を計算", Green);
   for (const auto &[obj, poly] : RigidBodyObject) {
#pragma omp parallel
      for (const auto &p : obj->getPoints())
#pragma omp single nowait
         if (p->isCaptured) {
            {
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                     p->interpolated_normal_SPH_original -= grad_w_Bspline(p->X, q->X, p->radius_SPH);
                  });
               p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
               if (!isFinite(p->interpolated_normal_SPH))
                  p->interpolated_normal_SPH = {0., 0., 1.};
            }

            // 不必要な壁粒子を除外．
            //\label{SPH:capture_condition_2nd}
            auto captureCondition = [&](const auto &q) {
               return Distance(p, q) < q->radius_SPH && p != q && (VectorAngle(p->interpolated_normal_SPH, q->X - p->X) < std::numbers::pi / 4);
            };

            if (net->BucketPoints.any_of(p->X, p->radius_SPH, captureCondition))
               p->isCaptured = true;
            else
               p->isCaptured = false;
         }
   }
   for (const auto &[obj, poly] : RigidBodyObject)
      for (const auto &p : obj->getPoints()) {
         if (p->isCaptured)
            wall_p.emplace(p);
         // if (p->isFluid)
         //    wall_as_fluid.emplace(p);
      }
};

/*DOC_EXTRACT SPH

## 法線方向の計算と水面の判定

*/

void setFreeSurface(auto &net, const auto &RigidBodyObject) {

   DebugPrint("水粒子のオブジェクト外向き法線方向を計算", Green);
// refference: A. Krimi, M. Jandaghian, and A. Shakibaeinia, Water (Switzerland), vol. 12, no. 11, pp. 1–37, 2020.

/*DOC_EXTRACT SPH

### 法線方向の計算

CHECKED \ref{SPH:interpolated_normal_SPH}{単位法線ベクトル}: $`{\bf n}_i = {\rm Normalize}\left(-\sum_j {\frac{m_j}{\rho_j} \nabla W_{ij} }\right)`$

単位法線ベクトルは，`interpolated_normal_SPH`としている．

*/
#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      // 初期化
      p->COM_SPH.fill(0.);
      p->interpolated_normal_SPH_original.fill(0.);
      p->interpolated_normal_SPH_original_next.fill(0.);
      double total_vol = 0, w;

      net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
         p->COM_SPH += (q->X - p->X) * w;
         total_vol += w;
         p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         p->interpolated_normal_SPH_original_next -= V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
      });

      std::vector<Tddd> near_wall_particle, near_wall_particle_next;
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            if (q->isCaptured) {
               w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
               p->COM_SPH += (q->X - p->X) * w;
               total_vol += w;

               if (Distance(p, q) < p->radius_SPH / p->C_SML * 1.5)
                  near_wall_particle.emplace_back(q->normal_SPH);
               p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);

               if (Distance(X_next(p), X_next(q)) < p->radius_SPH / p->C_SML * 1.5)
                  near_wall_particle_next.emplace_back(q->normal_SPH);
               p->interpolated_normal_SPH_original_next -= V_next(q) * grad_w_Bspline(X_next(p), X_next(q), p->radius_SPH);
            }
         });

      p->COM_SPH /= total_vol;

      // \label{SPH:interpolated_normal_SPH}
      p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);  //\label{SPH:interpolated_normal_SPH}
      for (const auto &n : near_wall_particle) {
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
         return Distance(p, q) < radius && p != q && (VectorAngle(p->interpolated_normal_SPH, -q->normal_SPH) < std::numbers::pi / 180. * 60);
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
   if (net->surfaceNet != nullptr)
      delete net->surfaceNet;

   net->surfaceNet = new Network();

   DebugPrint("水面補助粒子の作成", Green);
   for (const auto &p : net->getPoints()) {
      p->isAuxiliary = false;
      p->auxiliaryPoints.fill(nullptr);
      double d = 0;
      if (p->isSurface) {
         for (auto &auxp : p->auxiliaryPoints) {

            // double distance = 1E+20;
            // net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            //    if (q != p) {
            //       auto tmp = Distance(p->X, q);
            //       if (distance > tmp) {
            //          distance = tmp;
            //       }
            //    }
            // });
            // d = distance;
            // d += p->radius_SPH / p->C_SML;

            auxp = new networkPoint(net->surfaceNet, aux_position(p));
            auxp->radius_SPH = p->radius_SPH;
            auxp->surfacePoint = p;
            auxp->isAuxiliary = true;
            auxp->isSurface = false;
            auxp->p_SPH = p->p_SPH;
            auxp->U_SPH = p->U_SPH;
            auxp->setDensityVolume(p->rho, p->volume);
#if defined(USE_RungeKutta)
            auxp->RK_U = p->RK_U;
            auxp->RK_X = p->RK_X;
            auxp->RK_P = p->RK_P;
            auxp->RK_rho = p->RK_rho;
#elif defined(USE_LeapFrog)
            auxp->LPFG_X = p->LPFG_X;
            auxp->LPFG_rho = p->LPFG_rho;
#endif
         }
      }
   }
   net->surfaceNet->setGeometricProperties();
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT SPH

## $`\nabla^2 {\bf u}_i`$の計算

CHECKED: \ref{SPH:lapU}{ラプラシアンの計算方法}: $`\nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

*/

// b$ ------------------------------------------------------ */
// b$                    ∇.∇UとU*を計算                       */
// b$ ------------------------------------------------------ */

auto calcLaplacianU(const auto &points, const std::unordered_set<Network *> &target_nets, const double dt) {

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->checked_points_in_radius_SPH = A->checked_points_in_radius_of_fluid_SPH = A->checked_points_SPH = 0;
      A->div_U = 0.;
      A->lap_U.fill(0.);
      A->b_vector.fill(0.);
      //
      A->grad_coeff.clear();
      A->grad_coeff_next.clear();
      //$ ------------------------------------------ */
      /*DOC_EXTRACT SPH

      ### 高速化のための工夫

      何度か行う勾配の計算は，変数は違えど，変数の係数は同じである．
      ここで，その係数を`std::unordered_map`で保存しておくことにする．
      `A->grad_coeff`と`A->grad_coeff_next`に保存する．

      NOTE: `A->grad_coeff`と`A->grad_coeff_next`は，自身もキーとして含む．使う時に注意する．

      */
      //$ ------------------------------------------ */
      auto add_to_unmap = [&](const auto &key, const Tddd coef) {
         auto it = A->grad_coeff.find(key);
         if (it != A->grad_coeff.end())
            it->second += coef;
         else
            A->grad_coeff.emplace_hint(it, key, coef);
      };
      auto add_to_unmap_next = [&](const auto &key, const Tddd coef) {
         auto it = A->grad_coeff_next.find(key);
         if (it != A->grad_coeff_next.end())
            it->second += coef;
         else
            A->grad_coeff_next.emplace_hint(it, key, coef);
      };
      auto add_lap_U = [&](const auto &B) {
         const auto Uij = A->U_SPH - B->U_SPH;
         A->div_U += B->volume * Dot(B->U_SPH - A->U_SPH, grad_w_Bspline(A->X, B->X, A->radius_SPH));
         A->lap_U += 2 * B->mass / A->rho * Uij * Dot_grad_w_Bspline_Dot(A->X, B->X, A->radius_SPH);  //\label{SPH:lapU}

         // just counting
         if (Between(Distance(A, B), {1E-12, A->radius_SPH})) {
            A->checked_points_in_radius_SPH++;
            if (B->getNetwork()->isFluid || B->isFluid)
               A->checked_points_in_radius_of_fluid_SPH++;

            // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP1}
            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->mass / A->rho * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP2}
            // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP3}

            {
               auto coef = B->mass / A->rho * grad_w_Bspline(A->X, B->X, A->radius_SPH);
               add_to_unmap(A, -coef);
               add_to_unmap(B, coef);
            }
            {
               auto coef = B->mass / A->rho * grad_w_Bspline(A->X + dt * A->U_SPH, B->X + dt * B->U_SPH, A->radius_SPH);
               add_to_unmap_next(A, -coef);
               add_to_unmap_next(B, coef);
            }
         }
         A->checked_points_SPH++;
      };
      // sum 計算
      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            // if (B->isCaptured)
            {
               // 全ての壁粒子の流速はゼロなのだから，isCapturedされていないものを含めても問題ない．
               // ただこの後の圧力の計算においては，isCapturedされていないものは含めない．圧力方程式をうまく立てれないから．
               add_lap_U(B);
            }
         });
      //$ ------------------------------------------ */
      //\label{SPH:lapU_for_wall}
      if (A->getNetwork()->isRigidBody) {
         A->DUDt_SPH_ *= 0;
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH *= 0;
         A->tmp_U_SPH *= 0;
         A->tmp_X = A->X;
         A->DrhoDt_SPH *= 0;
      } else {
         A->DUDt_SPH_ = A->DUDt_SPH;
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH = nu * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt
         A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
         A->tmp_X = A->X + A->tmp_U_SPH * dt;
         A->DrhoDt_SPH = -A->rho * A->div_U;
      }
      // 予めb_vectorを計算するようにした．
      //     どっちがいいのか
      //         しかし，これはsumを取る必要がない．
      //             だたこれで，任意の場所でb_vectorを計算しやすくなる．
      //                 Laplacianが水面補助粒子を考慮していないこと，
      //                     補助粒子が密度変化を考慮できないこと
      //                         など
      //                             gradが

      //$ ------------------------------------------ */
      // \label{SPH:Poisson_b_vector}

      // auto add_b_vector = [&](const auto &B) {
      //    auto w = B->volume * w_Bspline(Norm(A->X - B->X), A->radius_SPH);
      //    A->b_vector += w * (B->U_SPH / dt + B->mu_SPH / B->rho * B->lap_U);  // + (A->rho * _GRAVITY3_);
      // };
      // // sum 計算
      // for (const auto &net : target_nets)
      //    net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
      //       if (B->isCaptured) {
      //          add_b_vector(B);
      //       }
      //    });

      A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;

      if (A->vec_time_SPH.size() > 10) {

#if defined(USE_RungeKutta)
         double current_time = A->RK_X.getTime();
         double next_time = A->RK_X.getNextTime();
#elif defined(USE_LeapFrog)
         double current_time = A->LPFG_X.get_t();
         double next_time = A->LPFG_X.get_t() + dt;
#endif
         std::vector<double> time3 = {next_time, current_time};
         std::array<double, 3> U1, U2, U3;
         U1 = A->U_SPH;
         if (*(A->vec_time_SPH.rbegin()) == current_time) {
            time3.push_back(*(A->vec_time_SPH.rbegin() + 1));
            U2 = *(A->vec_U_SPH.rbegin() + 1);
            time3.push_back(*(A->vec_time_SPH.rbegin() + 2));
            U3 = *(A->vec_U_SPH.rbegin() + 2);
         } else {
            time3.push_back(*(A->vec_time_SPH.rbegin() + 0));
            U2 = *(A->vec_U_SPH.rbegin() + 0);
            time3.push_back(*(A->vec_time_SPH.rbegin() + 1));
            U3 = *(A->vec_U_SPH.rbegin() + 1);
         }
         InterpolationLagrange<double> lag(time3);
         auto D = lag.DN(current_time);
         A->b_vector = -(D[1] * U1 + D[2] * U2 + D[3] * U3) + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
      }

      //
      A->b_vector3[2] = A->b_vector3[1];
      A->b_vector3[1] = A->b_vector3[0];
      A->b_vector3[0] = A->b_vector;
   }
};

// b% -------------------------------------------------------------------------- */

/*DOC_EXTRACT SPH

## ポアソン方程式$`\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) = b`$

### ポアソン方程式

次の時刻の流れ場を発散なし$`\nabla\cdot{\bf u}^{n+1}=0`$としてくれる
$`\frac{D {\bf u}}{D t} =-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}`$を使って，流速と粒子位置を時間発展させたい．
そのためには，圧力$`p^{n+1}`$を適切に決める必要がある．

$`\frac{D {\bf u}}{D t}`$は．$`\frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t}`$と離散化し条件を考えてみる．

```math
\frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t} =-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}
```

次時刻の発散の演算は，次時刻における粒子配置に基づき行われるので，現在の粒子配置に基づく発散演算とは区別すべきである．
現在の微分演算を$`\nabla^{n}`$とし，次時刻の微分演算を$`\nabla^{n+1}`$とする．
$`\nabla^{n+1}`$を上の式に作用させると，

```math
\nabla^{n+1}\cdot {\bf u}^{n+1} = \nabla^{n+1} \cdot{\bf u}^{n} - \Delta t \nabla^{n+1} \cdot\left(\frac{1}{\rho} \nabla^{n} p^{n+1}-\nu \nabla^{n2} {\bf u}^n-{\bf g}\right)
```

右辺がゼロとなれば，次時刻の流速の発散がゼロ，$`\nabla^{n+1}{\bf u}^{n+1}=0`$になる：

```math
\begin{align*}
&&0 &= \nabla^{n+1} \cdot{\bf u}^{n} - \Delta t \nabla^{n+1} \cdot\left(\frac{1}{\rho} \nabla^{n} p^{n+1}-\nu \nabla^{n2} {\bf u}^n-{\bf g}\right)\\
&\rightarrow&\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= \frac{1}{\Delta t}\nabla^{n+1} \cdot{\bf u}^{n} + \nabla^{n+1} \cdot\left(\nu^n \nabla^{n2} {\bf u}^n  + {\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= \nabla^{n+1} \cdot\left(\frac{1}{\Delta t}{\bf u}^{n} +\nu^n \nabla^{n2} {\bf u}^n  + {\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) &= b = \nabla^{n+1} \cdot {\bf b}^n,\quad  {\bf b}^n=\frac{1}{\Delta t}{\bf u}^{n} +\nu^n \nabla^{n2} {\bf u}^n
\end{align*}
```

重力の発散はゼロなので消した．

### 右辺，$`b`$，`PoissonRHS`について

この$`b`$を`PoissonRHS`とする．（仮流速は$`{\bf u}^* = \frac{\Delta t}{\rho}{\bf b}^n`$と同じ）．
$`{\bf b}^n`$ （\ref{SPH:Poisson_b_vector}{`Poisson_b_vector`}）が計算できるように，$`{\bf u}^n`$と$`\nabla^2 {\bf u}^n`$を計算しておく．

CHECKED: \ref{SPH:div_b_vector}{発散の計算方法}: $`b=\nabla\cdot{\bf b}^n=\sum_{j}\frac{m_j}{\rho_j}({\bf b}_j^n-{\bf b}_i^n)\cdot\nabla W_{ij}`$

### 左辺について

壁粒子の圧力は時間積分して計算しないので，毎時刻，壁粒子の$`p^{n+1}`$を計算する必要がある．

**EISPH**

   1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
   2. 流体粒子の圧力$`p^{n+1}`$の計算

**ISPH**

   - ISPHは作ったポアソン方程式を作成し解くことで圧力を計算する

CHECKED: \ref{SPH:lapP}{ラプラシアンの計算方法}: $`\nabla^2 p^{n+1}=\sum_{j}A_{ij}(p_i^{n+1} - p_j^{n+1}),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

### 水面の計算補助粒子`auxiliaryPoints`

水面においては，流速の発散ゼロ$`\nabla^{n+1} {\bf u}^{n+1}=0`$と$`p^{n+1}=0`$が満たされる必要がある．
水面外部には，粒子がないので，求めた水面圧力は，ゼロであっても，圧力勾配は誤差を含み，$`\nabla^{n+1} {\bf u}^{n+1}=0`$は満足されない．
そこで，\ref{SPH:auxiliaryPoints}{水面の計算補助粒子}を水面外部に追加し，この点を適切計算することで，$`\nabla^{n+1} {\bf u}^{n+1}=0`$が満足されるように工夫する．

*/

/*DOC_EXTRACT SPH

### 次時刻の発散演算，$`\nabla^{n+1} \cdot {\bf b}^n = \sum_j \dfrac{m_j}{\rho_j^{n+1}}({\bf b}_j^n-{\bf b}_i^n)\cdot \nabla W({\bf x}_i^{n+1},{\bf x}_j^{n+1},h)`$

$`\nabla^{n+1}`$の計算には，$`\rho^{n+1}`$, $`{\bf x}^{n+1}= {\bf x}^{n} + {\bf u}^{n+1} \Delta t`$が必要である．

* \ref{SPH:volume_next}{次時刻の粒子体積}
* \ref{SPH:rho_next}{次時刻の粒子密度}
* \ref{SPH:position_next}{次時刻の粒子位置}

*/

// \label{SPH:setPoissonEquation}
void setPoissonEquation(const std::unordered_set<networkPoint *> &points,
                        const std::unordered_set<Network *> &target_nets,
                        const double dt, const double &particle_spacing) {

   // これは現在の粒子位置で計算するが，この微分は次時刻の粒子位置で計算する
   // auto Poisson_b_vector = [&](const networkPoint *A, const double dt) {
   //    return A->RK_U.get_U0_for_SPH() / dt + A->mu_SPH / A->rho * A->lap_U;  // + (A->rho * _GRAVITY3_);
   // };

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
      A->PoissonRHS = 0;
      A->column_value.clear();
      Tddd origin_x, origin_b;
      auto origin = A;
      auto b_vector = [&](const auto &B) {
         return B->b_vector;
         // Crank-Nicolson second order
         // auto tmp = 3 * B->b_vector3[0] - B->b_vector3[1];
         // return tmp / 2;
         // Crank-Nicolson third order
         // auto tmp = 23 * B->b_vector3[0] - 16 * B->b_vector3[1] + 5 * B->b_vector3[2];
         // return tmp / 12;
      };

      /*DOC_EXTRACT SPH

      ### ポアソン方程式の作成のコーディング

      各粒子`A`に対して，方程式を作成する．

      まずは，\ref{SPH:whereToMakeTheEquation}{方程式を立てる位置を決める．}

      */

      // \label{SPH:whereToMakeTheEquation}
      if (A->isAuxiliary) {
         origin = A->surfacePoint;
         origin_x = X_next(origin);
         origin_b = b_vector(origin);
      } else if (A->getNetwork()->isRigidBody) {
         // origin = getClosestExcludeRigidBody(A, target_nets);
         // origin_x = X_next(origin);
         // origin_x = X_next(A);  // + 1.0001 * A->normal_SPH;
         // origin_b = b_vector(A);

         // origin_b.fill(0.);
         origin_x = X_next(A);  // + A->radius_SPH / 2. * Normalize(A->normal_SPH);
         origin_b = b_vector(origin);
         // sum 計算
         for (const auto &net : target_nets)
            net->BucketPoints.apply(origin_x, A->radius_SPH, [&](const auto &B) {
               if (B->isCaptured) {
                  auto w = B->volume * w_Bspline(Norm(origin_x - B->X), A->radius_SPH);
                  origin_b += B->b_vector * w;
               }
            });

         // origin_b = A->b_vector;
         //
         // origin_x = X_next(A) + A->normal_SPH + RandomReal({-1., 1.}) * 1E-10;
         // auto origin = getClosestExcludeRigidBody(A, target_nets);
         // origin_b *= 0;  // Poisson_b_vector(origin, dt);

         // origin_x = X_next(origin);
         // origin_b = Poisson_b_vector(origin, dt);
         // origin_b = Poisson_b(origin_x, origin->radius_SPH, dt, target_nets);
      } else {
         origin = A;
         origin_x = X_next(A);
         origin_b = b_vector(A);
      }

      double total_weight = 0, P_wall = 0, dP;
      A->density_based_on_positions = 0;

      /*DOC_EXTRACT SPH

      各粒子`A`が，流体か壁か補助粒子か水面かによって，方程式が異なる．

      |方程式|目的|
      |:---------|---|
      | IMPLEMENTED  \ref{SPH:PoissonEquation}{ポアソン方程式}              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
      | NOTIMPLEMENTED  \ref{SPH:ImpermeableCondition}{不透過条件}         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
      | NOTIMPLEMENTED  \ref{SPH:AtmosphericPressureCondition}{大気圧条件} | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない． |

      各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．

      */

      auto ImpermeableCondition = [&](const auto &B /*column id*/) {  // \label{SPH:ImpermeableCondition}
         A->PoissonRHS -= (V_next(B) * Dot(b_vector(B), Normalize(origin->normal_SPH)) * w_Bspline(Norm(origin_x - B->X), origin->radius_SPH));
         auto coeff = V_next(B) * Dot(grad_w_Bspline(origin_x, B->X, origin->radius_SPH), Normalize(origin->normal_SPH));  // こっちはOKだろう．
         A->increment(B, coeff);
      };

      auto AtmosphericPressureCondition = [&](const auto &p) {  //  \label{SPH:AtmosphericPressureCondition}
         // pの圧力を完全にゼロにする条件
         A->PoissonRHS = 0;
         A->increment(p, 1.);
      };

      auto PoissonEquation = [&](const auto &B /*column id*/) {  // \label{SPH:PoissonEquation}
         if (!B->isAuxiliary) {
            A->PoissonRHS += V_next(B) * Dot(b_vector(B) - origin_b, grad_w_Bspline(origin_x, X_next(B), origin->radius_SPH));  // \label{SPH:div_b_vector}
            A->density_based_on_positions += B->volume * w_Bspline(Norm(origin_x - X_next(B)), origin->radius_SPH);
         }
         Aij = 2. * B->mass / rho_next(origin) * Dot_grad_w_Bspline_Dot(origin_x, X_next(B), origin->radius_SPH);  //\label{SPH:lapP}
         // for ISPH
         A->increment(origin, Aij / origin->rho);
         A->increment(B, -Aij / origin->rho);
         // for EISPH
         sum_Aij_Pj += Aij * B->p_SPH;
         sum_Aij += Aij;
      };

      /*DOC_EXTRACT SPH

      ### ポアソン方程式の作成

      */

      if (A->isAuxiliary) {
         AtmosphericPressureCondition(A->surfacePoint);
      } else if (A->getNetwork()->isRigidBody) {
         for (const auto &net : target_nets) {
            net->BucketPoints.apply(origin_x, A->radius_SPH * 1.1, [&](const auto &B) {
               if (B->isCaptured) {
                  PoissonEquation(B);
                  if (B->isSurface)
                     for (const auto &AUX : B->auxiliaryPoints)
                        PoissonEquation(AUX);
                  // for mapping to wall
                  total_weight += B->volume * w_Bspline(Norm(origin_x - X_next(B)), A->radius_SPH);
                  dP = Dot(X_next(A) - origin_x, B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
                  P_wall += (B->p_SPH + dP) * B->volume * w_Bspline(Norm(origin_x - X_next(B)), A->radius_SPH);
               }
            });
         }
      } else {
         for (const auto &net : target_nets) {
            net->BucketPoints.apply(origin_x, A->radius_SPH * 1.1, [&](const auto &B) {
               if (B->isCaptured) {
                  PoissonEquation(B);
                  if (B->isSurface)
                     for (const auto &AUX : B->auxiliaryPoints)
                        PoissonEquation(AUX);
                  // for mapping to wall
                  total_weight += B->volume * w_Bspline(Norm(origin_x - X_next(B)), A->radius_SPH);
                  dP = Dot(X_next(A) - origin_x, B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
                  P_wall += (B->p_SPH + dP) * B->volume * w_Bspline(Norm(origin_x - X_next(B)), A->radius_SPH);
               }
            });
         }
      }
      /* -------------------------------------------------------------------------- */
#if defined(Morikawa2019)
         /* SPH
         ### 圧力の安定化

         $`b = \nabla \cdot {{\bf b}^n} + \alpha \frac{\rho_w - \rho^*}{{\Delta t}^2}`$として計算を安定化させる場合がある．
         $`\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t`$と近似すると，

         ```math
         \rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t,\quad
         \frac{D\rho^\ast}{Dt} = - \rho \nabla\cdot{\bf u}^\ast,\quad
         \nabla\cdot{\bf u}^\ast = \frac{\Delta t}{\rho} \nabla\cdot{\bf b}^n
         ```

         であることから，$`(\rho_w - \rho^*) / {\Delta t^2}$は，$\nabla\cdot{\bf b}^n`$となって同じになる．

         しかし，実際には，$`\rho^*$は，$\nabla \cdot {{\bf b}^n}`$を使わずに，つまり発散演算を行わずに評価するので，
         計算上のようにはまとめることができない．

         $`\rho^*`$を計算する際に，$`\rho^\ast = \rho_w + \frac{D\rho^\ast}{Dt}\Delta t`$を使った場合，確かに上のようになるが，
         実際に粒子を仮位置に移動させその配置から$\rho^*$を計算した場合は，数値計算上のようにまとめることはできない．

         `PoissonRHS`,$b$の計算方法と同じである場合に限る．
         もし，計算方法が異なれば，計算方法の違いによって，安定化の効果も変わってくるだろう．

         */
         // if (A->isFluid) {
         //    // \label{SPH:pressure_stabilization}
         //    const double alpha = 0.1 * dt;
         //    // A->PoissonRHS += alpha * (_WATER_DENSITY_ - A->density_based_on_positions) / (dt * dt);
         //    A->PoissonRHS += alpha * (_WATER_DENSITY_ - A->rho) / (dt * dt);
         // }
#endif
      //% ------------------------------------------------------- */
      // A->div_tmpU = A->PoissonRHS * dt / A->rho;
      // A->DrhoDt_SPH = -A->rho * A->div_tmpU;
      // A->rho_ = A->rho + A->DrhoDt_SPH * dt;

      A->p_SPH_ = (A->PoissonRHS + sum_Aij_Pj) / sum_Aij;

      if (A->getNetwork()->isRigidBody) {
         if (total_weight > 0.001)
            A->p_SPH_ = P_wall / total_weight;
         else
            A->p_SPH_ = 0;
      }
   };
};

void setPressure(const std::unordered_set<networkPoint *> &points) {
   for (const auto &p : points)
      p->p_SPH = p->p_SPH_;
}

/*DOC_EXTRACT SPH

## ポアソン方程式の解法

ISPHのポアソン方程式を解く場合，\ref{SPH:gmres}{ここではGMRES法}を使う．

*/

#define USE_LAPACK

void solvePoisson(const std::unordered_set<networkPoint *> &fluid_particle,
                  const std::unordered_set<networkPoint *> &wall_as_fluid,
                  const std::unordered_set<Network *> &target_nets) {

   std::unordered_set<networkPoint *> points;
   points.reserve(fluid_particle.size() + wall_as_fluid.size() + 1000);

   for (const auto &p : fluid_particle) {
      points.emplace(p);
      if (p->isSurface)
         for (const auto &AUX : p->auxiliaryPoints)
            points.emplace(AUX);
   }

   for (const auto &p : wall_as_fluid)
      points.emplace(p);

   for (auto i = 0; const auto &p : points)
      p->setIndexCSR(i++);

   V_d b(points.size()), x0(points.size(), 0);

   for (const auto &p : points) {
      b[p->getIndexCSR()] = p->PoissonRHS;
      x0[p->getIndexCSR()] = p->p_SPH;
   }

   /* ------------------ preconditioning using diagonal value ------------------ */
   for (const auto &p : points) {
      double max = 0;
      // find max
      for (const auto &[_, v] : p->column_value)
         if (std::abs(v) > max)
            max = std::abs(v);
      // normalize
      b[p->getIndexCSR()] /= max;
      for (auto &[_, v] : p->column_value)
         v /= max;
   }
#if defined(USE_GMRES)
   for (auto i = 1; i < 5; i++) {
      gmres gm(points, b, x0, 100);  //\label{SPH:gmres}
      x0 = gm.x;
      std::cout << " gm.err : " << gm.err << std::endl;
   }
   gmres gm(points, b, x0, 100);
   std::cout << " gm.err : " << gm.err << std::endl;

   for (const auto &p : points)
      x0[p->getIndexCSR()] = p->p_SPH = gm.x[p->getIndexCSR()];

   std::cout << " gm.err : " << gm.err << std::endl;
#elif defined(USE_LAPACK)
   VV_d A(b.size(), V_d(b.size(), 0.));
   for (const auto &p : points) {
      auto i = p->getIndexCSR();
      for (const auto &[q, v] : p->column_value) {
         auto j = q->getIndexCSR();
         A[i][j] = v;
      }
   }
   lapack_lu lu(A);
   lu.solve(b, x0);
   for (const auto &p : points)
      p->p_SPH = x0[p->getIndexCSR()];
#elif defined(USE_LAPACK_SVD)
   VV_d A(b.size(), V_d(b.size(), 0.));
   for (const auto &p : points) {
      auto i = p->getIndexCSR();
      for (const auto &[q, v] : p->column_value) {
         auto j = q->getIndexCSR();
         A[i][j] = v;
      }
   }
   lapack_svd svd(A);
   svd.solve(b, x0);
   for (const auto &p : points)
      p->p_SPH = x0[p->getIndexCSR()];
#endif

   auto error = Norm(b - Dot(points, x0));
   std::cout << "actual error : " << error << std::endl;

   if (!isFinite(error))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error is not a finite");

   for (const auto &p : points)
      p->column_value.clear();
};

/* -------------------------------------------------------------------------- */

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */

/*DOC_EXTRACT SPH

## 圧力勾配$`\nabla p^{n+1}`$の計算

CHECKED: \ref{SPH:gradP1}{勾配の計算方法}: $`\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}`$

CHECKED: \ref{SPH:gradP2}{勾配の計算方法}: $`\nabla p_i = \rho_i \sum_{j} m_j \left(p_j - p_i\right) \nabla W_{ij}`$

CHECKED: \ref{SPH:gradP3}{勾配の計算方法}: $`\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}`$

*/

void gradP(const std::unordered_set<networkPoint *> &points, const std::unordered_set<Network *> &target_nets) {

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->gradP_SPH.fill(0.);

      auto add_gradP_SPH = [&](const auto &B) {
         // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP1}0.2647

         // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->mass / A->rho * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP2}
         //\label{SPH:gradP2}は，Aij = 1*..を使うと，0.132
         //\label{SPH:gradP2}は，Aij = 3*..を使うと，0.221
         //\label{SPH:gradP2}は，Aij = 2*..を使うと，0.21198

         //   A->gradP_SPH += (B->p_SPH - A->p_SPH) * V_next(B) * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP2}
         A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP3}0.34
         //\label{SPH:gradP3}は，Aij = 2*..を使うと，0.34
         //\label{SPH:gradP3}は，Aij = 2*..を使うと，Lagrangeを使うと，0.46
         //\label{SPH:gradP3}は，Aij = 3*..を使うと，0.49
         //\label{SPH:gradP3}は，Aij = 3*..を使うと，Lagrangeを使うと，0.4は超えたが，綺麗な結果ではなく，つぶれた．
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，_nextを使うと，0.24
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，X_next以外の_nextを使うと，0.277>
         //\label{SPH:gradP3}は，Aij = 2*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.34
         //\label{SPH:gradP3}は，Aij = 2.5*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.45
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.49
         //\label{SPH:gradP3}は，Aij = 4*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.457
      };

      for (const auto &net : target_nets) {
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               add_gradP_SPH(B);
               if (B->isSurface)
                  for (const auto &AUX : B->auxiliaryPoints)
                     add_gradP_SPH(AUX);
            }
         });
      }

      /*DOC_EXTRACT SPH

      $`\dfrac{D{\bf u}^n}{Dt} = - \frac{1}{\rho} \nabla p^{n+1} + \nu \nabla^2 {\bf u}^n + {\bf g}`$
      が計算できた．

      */

      A->DUDt_SPH -= A->gradP_SPH / A->rho;

      if (!isFinite(A->DUDt_SPH))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
   }
}

//@ -------------------------------------------------------- */
//@                        粒子の時間発展                      */
//@ -------------------------------------------------------- */
#define REFLECTION
void updateParticles(const auto &points,
                     const std::unordered_set<Network *> &target_nets,
                     const auto &RigidBodyObject,
                     const double &particle_spacing,
                     const double dt) {
   DebugPrint("粒子の時間発展", Green);
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
   {
      // テスト
      auto U = p->U_SPH;
      auto X_last = p->X;
#if defined(USE_RungeKutta)
      p->RK_X.push(p->U_SPH);  // 位置
      p->setXSingle(p->tmp_X = p->RK_X.getX());
      //
      p->RK_U.push(p->DUDt_SPH);  // 速度
      p->U_SPH = p->RK_U.getX();
         // auto getX = [&](const auto &p) { return p->RK_X.getX(p->U_SPH); };
         // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
#elif defined(USE_LeapFrog)
      p->LPFG_X.push(p->DUDt_SPH);  // 速度
      p->U_SPH = p->LPFG_X.get_v();
      p->setXSingle(p->tmp_X = p->LPFG_X.get_x());
         // auto getX = [&](const auto &p) { return p->X; };
#endif

#if defined(REFLECTION)
      int count = 0;
      //\label{SPH:reflection}
      const double reflection_factor = .5;
      const double asobi = 0.;

      auto closest = [&]() {
         double distance = 1E+20;
         networkPoint *P = nullptr;
         for (const auto &[obj, _] : RigidBodyObject) {
            obj->BucketPoints.apply(X_next(p), p->radius_SPH, [&](const auto &q) {
               auto tmp = Distance(X_next(p), q);
               if (distance > tmp) {
                  distance = tmp;
                  P = q;
               }
            });
         }
         return P;
      };

      bool isReflected = true;
      while (isReflected && count++ < 5) {
         isReflected = false;
         networkPoint *closest_wall_point;
         if (closest_wall_point = closest()) {
            auto modify_position = particle_spacing * Normalize(X_next(p) - closest_wall_point->X) + closest_wall_point->X;
            auto ovre_run = ((1. - asobi) * particle_spacing - Distance(closest_wall_point->X, X_next(p))) / 2.;
            if (Distance(closest_wall_point->X, X_next(p)) < particle_spacing)
               if (ovre_run > 0.) {
                  auto normal_distance = Norm(Projection(X_next(p) - closest_wall_point->X, closest_wall_point->normal_SPH));
                  if (Dot(p->U_SPH, closest_wall_point->normal_SPH) < 0) {
   #if defined(USE_RungeKutta)
                     p->DUDt_SPH -= (1. + reflection_factor) * Projection(p->U_SPH, closest_wall_point->normal_SPH) / dt;
                     p->RK_U.repush(p->DUDt_SPH);  // 速度
                     p->U_SPH = p->RK_U.getX();
                     //
                     // p->RK_X.repush(p->U_SPH);  // 位置
                     // p->setXSingle(p->tmp_X = p->RK_X.X_next());
                     isReflected = true;
   #elif defined(USE_LeapFrog)
                     p->DUDt_SPH -= (1. + reflection_factor) * Projection(p->U_SPH, closest_wall_point->normal_SPH) / dt;
                     p->LPFG_X.repush(p->DUDt_SPH);  // 速度
                     p->U_SPH = p->LPFG_X.get_v();
                     p->setXSingle(p->tmp_X = modify_position);
                     isReflected = true;
   #endif
                     /* -------------------------------------------------------------------------- */
                     // p->DUDt_SPH -= (1. + reflection_factor) * Projection(p->U_SPH, closest_wall_point->normal_SPH) / dt;
                     // p->RK_U.repush(p->DUDt_SPH);  // 速度
                     // p->U_SPH = p->RK_U.X_next();
                     // // p->RK_X.repush(p->U_SPH);  // 位置
                     // // p->setXSingle(p->tmp_X = p->RK_X.X_next());
                     // isReflected = true;
                     /* -------------------------------------------------------------------------- */
                  }
               }
         }
      };

#endif
   }

   // \label{SPH:update_density}
   for (const auto &A : points) {
#if defined(USE_RungeKutta)
      A->DrhoDt_SPH = -A->rho * A->div_U;
      A->RK_rho.push(A->DrhoDt_SPH);  // 密度
      A->setDensity(A->RK_rho.get_x());
#elif defined(USE_LeapFrog)
      A->DrhoDt_SPH = -A->rho * A->div_U;
      A->LPFG_rho.push(A->DrhoDt_SPH);
      // A->setDensity(A->LPFG_rho.get_x());
      A->setDensity(_WATER_DENSITY_);
#endif
   }
}

/*DOC_EXTRACT SPH

## 注意点

WARNING: 計算がうまく行く設定を知るために，次の箇所をチェックする．

**壁粒子**

- \ref{SPH:lapU_for_wall}{壁粒子のラプラシアンの計算方法}
- \ref{SPH:setPoissonEquation}{圧力の計算方法}
   - \ref{SPH:whereToMakeTheEquation}{どの位置において方程式を立てるか}
- \ref{SPH:capture_condition_1st}{流体として扱う壁粒子を設定するかどうか}/\ref{SPH:capture_condition_2nd}{視野角に流体粒子が含まない壁粒子は除外する}
- \ref{SPH:map_fluid_pressure_to_wall}{壁粒子の圧力をどのように壁面にマッピングするか}
- \ref{SPH:interpolated_normal_SPH}{壁粒子の法線方向ベクトルの計算方法}
- \ref{SPH:reflection}{反射の計算方法}

**水面粒子**

- \ref{SPH:water_surface_pressure}{水面粒子の圧力をゼロにするかどうか}
- \ref{SPH:auxiliaryPoints}{補助粒子の設定はどうなっているか}

**その他**

- \ref{SPH:update_density}{密度を更新するかどうか}
- \ref{SPH:pressure_stabilization}{圧力の安定化をするかどうか}
- \ref{SPH:RK_order}{ルンゲクッタの段数}


壁のwall_as_fluidは繰り返しで計算するのはどうか？

*/

#endif