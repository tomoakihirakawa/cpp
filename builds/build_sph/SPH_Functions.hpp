#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"

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
#define new_method

Tddd X_SPP(const networkPoint *p, const double c = 1.) {
   return p->X + c * Normalize(p->interpolated_normal_SPH);
};

/*DOC_EXTRACT SPH

## 法線方向の計算と水面の判定

*/

void setNormal_Surface(auto &net, const std::unordered_set<networkPoint *> &wall_p, const auto &RigidBodyObject) {

   DebugPrint("水粒子のオブジェクト外向き法線方向を計算", Green);
// refference: A. Krimi, M. Jandaghian, and A. Shakibaeinia, Water (Switzerland), vol. 12, no. 11, pp. 1–37, 2020.
#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      // 初期化
      p->COM_SPH.fill(0.);
      p->interpolated_normal_SPH_original.fill(0.);
      p->interpolated_normal_SPH_original_all.fill(0.);
      double total_vol = 0, w;

      /*DOC_EXTRACT SPH

      ### 法線方向の計算

      CHECKED \ref{SPH:interpolated_normal_SPH}{単位法線ベクトル}: ${\bf n}_i = {\rm Normalize}\left(-\sum_j {\frac{m_j}{\rho_j} \nabla W_{ij} }\right)$

      単位法線ベクトルは，`interpolated_normal_SPH`としている．

      */

      net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
         w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
         p->COM_SPH += (q->X - p->X) * w;
         total_vol += w;
         p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         p->interpolated_normal_SPH_original_all -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
      });

      std::vector<Tddd> normal_of_near_wall_particle;
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
            w = q->volume * w_Bspline(Norm(p->X - q->X), p->radius_SPH);
            p->COM_SPH += (q->X - p->X) * w;
            total_vol += w;
            if (Distance(p, q) < p->radius_SPH / p->C_SML * 1.5)
               normal_of_near_wall_particle.emplace_back(q->normal_SPH);
            p->interpolated_normal_SPH_original -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
            p->interpolated_normal_SPH_original_all -= q->volume * grad_w_Bspline(p->X, q->X, p->radius_SPH);
         });

      p->COM_SPH /= total_vol;
      p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);  //\label{SPH:interpolated_normal_SPH}
      p->interpolated_normal_SPH_all = Normalize(p->interpolated_normal_SPH_original_all);

      for (const auto &n : normal_of_near_wall_particle) {
         p->interpolated_normal_SPH = Normalize(Chop(p->interpolated_normal_SPH, n));
         p->interpolated_normal_SPH_all = Normalize(Chop(p->interpolated_normal_SPH_all, n));
      }

      if (!isFinite(p->interpolated_normal_SPH)) {
         p->interpolated_normal_SPH = {0., 0., 1.};
         p->interpolated_normal_SPH_all = {0., 0., 1.};
      }

      /*DOC_EXTRACT SPH

      ### 水面の判定

      `surface_condition0,1`の両方を満たす場合，水面とする．

      */

      p->isSurface = true;
      const auto radius = (p->radius_SPH / p->C_SML) * 3.;

      auto surface_condition0 = [&](const auto &q) {
         return Distance(p, q) < radius && p != q && (VectorAngle(p->interpolated_normal_SPH_all, q->X - p->X) < std::numbers::pi / 4);
      };

      auto surface_condition1 = [&](const auto &q) {
         return Distance(p, q) < radius && p != q && (VectorAngle(p->interpolated_normal_SPH_all, -q->normal_SPH) < std::numbers::pi / 180. * 60);
      };

      if (net->BucketPoints.any_of(p->X, radius, surface_condition0))
         p->isSurface = false;

      if (p->isSurface)
         for (const auto &[obj, poly] : RigidBodyObject)
            if (obj->BucketPoints.any_of(p->X, radius, surface_condition1))
               p->isSurface = false;
   }

   /* -------------------------------------------------------------------------- */
   // 水面ネットワークの初期化
   for (const auto &p : net->getPoints())
      p->auxiliaryPoints.fill(nullptr);

   DebugPrint("水面ネットワークの初期化", Green);
   if (net->surfaceNet != nullptr)
      delete net->surfaceNet;

   net->surfaceNet = new Network();

   DebugPrint("水面補助粒子の作成", Green);
   for (const auto &p : net->getPoints()) {
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
            d += p->radius_SPH / p->C_SML;

            auxp = new networkPoint(net->surfaceNet, X_SPP(p, d));
            auxp->radius_SPH = p->radius_SPH;
            auxp->surfacePoint = p;
            auxp->isAuxiliary = true;
            auxp->isSurface = false;
            auxp->p_SPH = p->p_SPH;
            auxp->U_SPH = p->U_SPH;
            auxp->setDensityVolume(p->rho, p->volume);
         }
         p->isAuxiliary = false;
      }
   }
   net->surfaceNet->setGeometricProperties();
   /* -------------------------------------------------------------------------- */

   DebugPrint("壁粒子のオブジェクト外向き法線方向を計算", Green);
   for (const auto &[obj, poly] : RigidBodyObject)
      for (const auto &p : obj->getPoints())
         p->interpolated_normal_SPH_original.fill(0.);
#pragma omp parallel
   for (const auto &p : wall_p) {
#pragma omp single nowait
      {
         p->interpolated_normal_SPH_original.fill(0.);
         for (const auto &[obj, poly] : RigidBodyObject)
            obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
               p->interpolated_normal_SPH_original -= grad_w_Bspline(p->X, q->X, p->radius_SPH);
            });
         p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
         if (!isFinite(p->interpolated_normal_SPH))
            p->interpolated_normal_SPH = {0., 0., 1.};
      }
   }
};

/*DOC_EXTRACT SPH

### 壁面粒子の流速と圧力

壁粒子の流速を流体粒子の流速に応じて変化させるとプログラムが煩雑になるので，
**ここでは**壁面粒子の流速は常にゼロに設定することにする．

壁粒子の圧力は，水が圧縮しないように各ステップ毎に計算し直す必要がある．

*/

#define Morikawa2019
#define new_method

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
      //$ ------------------------------------------ */
      auto add_lap_U = [&](const auto &B) {
         if (!B->isAuxiliary) {

            const auto Uij = A->U_SPH - B->U_SPH;
            A->div_U += B->volume * Dot(B->U_SPH - A->U_SPH, grad_w_Bspline(A->X, B->X, A->radius_SPH));
            A->lap_U += 2 * B->mass / A->rho * Uij * Dot_grad_w_Bspline_Dot(A->X, B->X, A->radius_SPH);  //\label{SPH:lapU}

            // just counting
            if (Between(Distance(A, B), {1E-12, A->radius_SPH})) {
               A->checked_points_in_radius_SPH++;
               if (B->getNetwork()->isFluid || B->isFluid)
                  A->checked_points_in_radius_of_fluid_SPH++;

               // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(getX(A), getX(B), A->radius_SPH);  //\label{SPH:gradP1}
               // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->mass / A->rho * grad_w_Bspline(getX(A), getX(B), A->radius_SPH);  //\label{SPH:gradP2}
               // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(getX(A), getX(B), A->radius_SPH);  //\label{SPH:gradP3}

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
         }
      };
      //$ ------------------------------------------ */
      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               add_lap_U(B);
               // if (B->isSurface)
               //    for (const auto &AUX : B->auxiliaryPoints)
               //       add_lap_U(AUX);
            }
         });
      //$ ------------------------------------------ */
      A->DUDt_SPH_ = A->DUDt_SPH;
      double nu = A->mu_SPH / A->rho;
      A->DUDt_SPH = nu * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt
      A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
      A->tmp_X = A->X + A->tmp_U_SPH * dt;
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

次時刻の流速の発散がゼロ，$`\nabla^{n+1}{\bf u}^{n+1}=0`$になるには

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

void PoissonEquation(const std::unordered_set<networkPoint *> &points,
                     const std::unordered_set<Network *> &target_nets,
                     const double dt, const double &particle_spacing,
                     const std::function<Tddd(const networkPoint *)> &getX) {

   auto V_next = [&](const auto &p) {
      if (p->isAuxiliary || p->getNetwork()->isRigidBody)
         return p->mass / p->rho;
      else {
#if defined(USE_RungeKutta)
         return p->mass / p->RK_rho.getX(-p->rho * p->div_U);
#else
         return p->mass / p->rho;
#endif
      }
   };

   auto rho_next = [&](const auto &p) {
      if (p->isAuxiliary || p->getNetwork()->isRigidBody)
         return p->rho;
      else {
#if defined(USE_RungeKutta)
         return p->RK_rho.getX(-p->rho * p->div_U);
#else
         return p->mass / p->rho;
#endif
      }
   };

   // これは現在の粒子位置で計算するが，この微分は次時刻の粒子位置で計算する
   // auto Poisson_b_vector = [&](const networkPoint *A, const double dt) {
   //    return A->RK_U.get_U0_for_SPH() / dt + A->mu_SPH / A->rho * A->lap_U;  // + (A->rho * _GRAVITY3_);
   // };

   // \label{SPH:Poisson_b_vector}
   auto Poisson_b_vector = [&](const networkPoint *A, const double dt) {
      return A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + (A->rho * _GRAVITY3_);//
   };

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
      A->PoissonRHS = 0;
      A->column_value.clear();
      Tddd origin_x, origin_b;

      /*DOC_EXTRACT SPH

      ### ポアソン方程式の作成のコーディング

      各粒子`A`に対して，方程式を作成する．

      まずは，\ref{SPH:whereToMakeTheEquation}{方程式を立てる位置を決める．}

      */

      // \label{SPH:whereToMakeTheEquation}
      if (A->isAuxiliary) {
         origin_x = getX(A->surfacePoint);
         origin_b = Poisson_b_vector(A->surfacePoint, dt);
      } else if (A->getNetwork()->isRigidBody) {
         origin_x = getX(A);
         origin_b = Poisson_b_vector(A, dt);
      } else {
         origin_x = getX(A);
         origin_b = Poisson_b_vector(A, dt);
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

      // \label{SPH:ImpermeableCondition}
      auto ImpermeableCondition = [&](const auto &B /*column id*/) {
         A->PoissonRHS -= (V_next(B) * Dot(Poisson_b_vector(B, dt), Normalize(A->normal_SPH)) * w_Bspline(Norm(origin_x - B->X), A->radius_SPH));
         auto coeff = V_next(B) * Dot(grad_w_Bspline(origin_x, B->X, A->radius_SPH), Normalize(A->normal_SPH));  // こっちはOKだろう．
         A->increment(B, coeff);
      };

      // 壁付近の水面との違いが出るので修正し，完全にゼロとする．
      //  \label{SPH:AtmosphericPressureCondition}
      // auto AtmosphericPressureCondition = [&]() {
      //    A->PoissonRHS = 0;
      //    A->increment(A, 1.);
      // };

      auto AtmosphericPressureCondition = [&](const auto &B /*column id*/) {
         A->PoissonRHS = 0;
         A->increment(B, V_next(B) * w_Bspline(Norm(origin_x - getX(B)), A->radius_SPH));
      };

      // \label{SPH:PoissonEquation}
      auto PoissonEquation = [&](const auto &B /*column id*/) {
         if (!B->isAuxiliary) {
            A->PoissonRHS += V_next(B) * Dot(Poisson_b_vector(B, dt) - origin_b, grad_w_Bspline(origin_x, getX(B), A->radius_SPH));  // \label{SPH:div_b_vector}
            A->density_based_on_positions += B->volume * w_Bspline(Norm(origin_x - getX(B)), A->radius_SPH);
         }

         Aij = 2. * B->mass / rho_next(A) * Dot_grad_w_Bspline_Dot(origin_x, getX(B), A->radius_SPH);  //\label{SPH:lapP}

         // for ISPH
         A->increment(A, Aij / A->rho);
         A->increment(B, -Aij / A->rho);

         // for EISPH
         sum_Aij_Pj += Aij * B->p_SPH;
         sum_Aij += Aij;
      };

      if (A->isAuxiliary) {
         A->PoissonRHS = 0;
         A->increment(A->surfacePoint, 1.);
      } else
         for (const auto &net : target_nets) {
            net->BucketPoints.apply(origin_x, A->radius_SPH * 1.1, [&](const auto &B) {
               if (B->isCaptured) {
                  PoissonEquation(B);
                  if (B->isSurface)
                     for (const auto &AUX : B->auxiliaryPoints)
                        PoissonEquation(AUX);

                  // for mapping to wall
                  total_weight += B->volume * w_Bspline(Norm(origin_x - getX(B)), A->radius_SPH);
                  dP = Dot(getX(A) - origin_x, B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
                  P_wall += (B->p_SPH + dP) * B->volume * w_Bspline(Norm(origin_x - getX(B)), A->radius_SPH);
               }
            });
         }

         /* -------------------------------------------------------------------------- */
#if defined(Morikawa2019)
            /* SPH
            ### 圧力の安定化

            $b = \nabla \cdot {{\bf b}^n} + \alpha \frac{\rho_w - \rho^*}{{\Delta t}^2}$として計算を安定化させる場合がある．
            $\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t$と近似すると，

            $$
            \rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t,\quad
            \frac{D\rho^\ast}{Dt} = - \rho \nabla\cdot{\bf u}^\ast,\quad
            \nabla\cdot{\bf u}^\ast = \frac{\Delta t}{\rho} \nabla\cdot{\bf b}^n
            $$

            であることから，$(\rho_w - \rho^*) / {\Delta t^2}$は，$\nabla\cdot{\bf b}^n$となって同じになる．

            しかし，実際には，$\rho^*$は，$\nabla \cdot {{\bf b}^n} $を使わずに，つまり発散演算を行わずに評価するので，
            計算上のようにはまとめることができない．

            $\rho^*$を計算する際に，$\rho^\ast = \rho_w + \frac{D\rho^\ast}{Dt}\Delta t$を使った場合，確かに上のようになるが，
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

void solvePoisson(const std::unordered_set<networkPoint *> &fluid_particle,
                  const std::unordered_set<networkPoint *> &wall_as_fluid,
                  const std::unordered_set<Network *> &target_nets) {

   size_t i = 0;
   std::unordered_set<networkPoint *> points;
   points.reserve(fluid_particle.size() + wall_as_fluid.size() + 1000);

   for (const auto &A : fluid_particle) {
      points.emplace(A);

      if (A->isSurface)
         for (const auto &AUX : A->auxiliaryPoints) {
            points.emplace(AUX);
            // respect AUX's pressure
            double c = 1.;
            AUX->PoissonRHS *= c;
            for (auto &[_, v] : AUX->column_value)
               v *= c;
         }
   }

   for (const auto &p : wall_as_fluid)
      points.emplace(p);

   for (const auto &p : points)
      p->setIndexCSR(i++);

   V_d b(points.size()), x0(points.size(), 0);

   for (const auto &p : points)
      b[p->getIndexCSR()] = p->PoissonRHS;

   // // store diagonal value
   // for (const auto &p : points) {
   //    // find max
   //    double max = 0;
   //    for (const auto &[_, v] : p->column_value)
   //       if (std::abs(v) > max)
   //          max = std::abs(v);
   //    p->diagonal_value = max;
   // }

   // // preconditioning using diagonal value
   // for (const auto &p : points) {
   //    b[p->getIndexCSR()] /= p->diagonal_value;
   //    for (auto &[_, v] : p->column_value)
   //       v /= p->diagonal_value;
   // }

   int N = 200;
   DebugPrint("gmres iteration ", N, Green);
   gmres gm(points, b, x0, N);  //\label{SPH:gmres}
   std::cout << " gm.err : " << gm.err << std::endl;
   gm.Restart(points, b, gm.x, N);
   std::cout << " gm.err : " << gm.err << std::endl;
   gm.Restart(points, b, gm.x, N);
   std::cout << " gm.err : " << gm.err << std::endl;
   gm.Restart(points, b, gm.x, N);
   std::cout << " gm.err : " << gm.err << std::endl;
   gm.Restart(points, b, gm.x, N);
   std::cout << " gm.err : " << gm.err << std::endl;

   // for (auto j = 0; j < 5; ++j) {
   //    std::cout << "j = " << j << std::endl;
   //    for (auto i = 0; i < 100; ++i) {
   //       std::cout << "i = " << i << std::endl;
   //       gm.Iterate(points);
   //       if (std::abs(gm.err) < 1)
   //          break;
   //       std::cout << "i, j = " << i << ", " << j << " gm.err : " << gm.err << std::endl;
   //    }
   //    x0 = gm.x;
   //    if (std::abs(gm.err) < 1)
   //       break;
   //    else
   //       gm.Restart(points, b, x0, N);

   //    std::cout << "j = " << j << " gm.err : " << gm.err << std::endl;
   // }

   for (const auto &p : points)
      x0[p->getIndexCSR()] = p->p_SPH = gm.x[p->getIndexCSR()];

   std::cout << " gm.err : " << gm.err << std::endl;
   std::cout << "actual error : " << Norm(b - Dot(points, x0)) << std::endl;

   // for (const auto &p : points)
   //    p->column_value.clear();
};

/* -------------------------------------------------------------------------- */

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */
/*DOC_EXTRACT SPH

## 圧力勾配$\nabla p^{n+1}$の計算

CHECKED: \ref{SPH:gradP1}{勾配の計算方法}: $\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}$

CHECKED: \ref{SPH:gradP2}{勾配の計算方法}: $\nabla p_i = \rho_i \sum_{j} m_j \left(p_j - p_i\right) \nabla W_{ij}$

CHECKED: \ref{SPH:gradP3}{勾配の計算方法}: $\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}$

*/

void gradP(const std::unordered_set<networkPoint *> &points,
           const std::unordered_set<Network *> &target_nets,
           const std::function<Tddd(const networkPoint *)> &getX) {

   auto V_next = [&](const auto &p) {
      if (p->isAuxiliary) {
         return p->mass / p->rho;
      } else if (p->getNetwork()->isRigidBody)
         return p->mass / p->rho;
      else
#if defined(USE_RungeKutta)
         return p->mass / p->RK_rho.getX(-p->rho * p->div_U);
#else
         return p->mass / p->rho;
#endif
   };

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->gradP_SPH.fill(0.);

      auto add_gradP_SPH = [&](const auto &B) {
         // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP1}
         // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->mass / A->rho * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP2}
         A->gradP_SPH += (B->p_SPH - A->p_SPH) * V_next(B) * grad_w_Bspline(getX(A), getX(B), A->radius_SPH);  //\label{SPH:gradP2}
                                                                                                               // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP3}
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
      auto getX = [&](const auto &p) { return p->RK_X.getX(p->U_SPH); };
         // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
#elif defined(USE_LeapFrog)
      p->LPFG_X.push(p->DUDt_SPH);  // 速度
      p->U_SPH = p->LPFG_X.get_v();
      p->setXSingle(p->tmp_X = p->LPFG_X.get_x());
      auto getX = [&](const auto &p) { return p->X; };
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
            obj->BucketPoints.apply(getX(p), p->radius_SPH, [&](const auto &q) {
               auto tmp = Distance(getX(p), q);
               if (distance > tmp) {
                  distance = tmp;
                  P = q;
               }
            });
         }
         return P;
      };

      bool isReflected = true;
      while (isReflected && count++ < 1) {
         isReflected = false;
         networkPoint *closest_wall_point;
         if (closest_wall_point = closest()) {
            auto modify_position = particle_spacing * Normalize(getX(p) - closest_wall_point->X) + closest_wall_point->X;
            auto ovre_run = ((1. - asobi) * particle_spacing - Distance(closest_wall_point->X, getX(p))) / 2.;
            if (Distance(closest_wall_point->X, getX(p)) < particle_spacing)
               if (ovre_run > 0.) {
                  auto normal_distance = Norm(Projection(getX(p) - closest_wall_point->X, closest_wall_point->normal_SPH));
                  if (Dot(p->U_SPH, closest_wall_point->normal_SPH) < 0) {
   #if defined(USE_RungeKutta)
                     p->DUDt_SPH -= (1. + reflection_factor) * Projection(p->U_SPH, closest_wall_point->normal_SPH) / dt;
                     p->RK_U.repush(p->DUDt_SPH);  // 速度
                     p->U_SPH = p->RK_U.getX();
                     //
                     // p->RK_X.repush(p->U_SPH);  // 位置
                     // p->setXSingle(p->tmp_X = p->RK_X.getX());
                     isReflected = true;
   #elif defined(USE_LeapFrog)
                     p->DUDt_SPH -= (1. + reflection_factor) * Projection(p->U_SPH, closest_wall_point->normal_SPH) / dt;
                     p->LPFG_X.repush(p->DUDt_SPH);  // 速度
                     p->U_SPH = p->LPFG_X.get_v();
                     // p->setXSingle(p->tmp_X = modify_position);
                     isReflected = true;
   #endif
                     /* -------------------------------------------------------------------------- */
                     // p->DUDt_SPH -= (1. + reflection_factor) * Projection(p->U_SPH, closest_wall_point->normal_SPH) / dt;
                     // p->RK_U.repush(p->DUDt_SPH);  // 速度
                     // p->U_SPH = p->RK_U.getX();
                     // // p->RK_X.repush(p->U_SPH);  // 位置
                     // // p->setXSingle(p->tmp_X = p->RK_X.getX());
                     // isReflected = true;
                     /* -------------------------------------------------------------------------- */
                  }
               }
         }
      };

#endif
   }

   // \label{SPH:update_density}

   //    for (const auto &A : points) {
   // #if defined(USE_RungeKutta)
   //       A->DrhoDt_SPH = -A->rho * A->div_U;
   //       A->RK_rho.push(A->DrhoDt_SPH);  // 密度
   //       A->setDensity(A->RK_rho.getX());
   // #elif defined(USE_LeapFrog)
   //       A->DrhoDt_SPH = -A->rho * A->div_U;
   //       A->LPFG_rho.push(A->DrhoDt_SPH);
   //       A->setDensity(A->rho + A->DrhoDt_SPH * dt);
   // #endif
   //    }
}

/*DOC_EXTRACT SPH

## 注意点

WARNING: 計算がうまく行く設定を知るために，次の箇所をチェックする．

- \ref{SPH:select_wall_as_fluid}{流体として扱う壁粒子を設定するかどうか}
- \ref{SPH:map_fluid_pressure_to_wall}{壁粒子の圧力をどのように壁面にマッピングするか}
- \ref{SPH:water_surface_pressure}{水面粒子の圧力をゼロにするかどうか}
- \ref{SPH:update_density}{密度を更新するかどうか}
- \ref{SPH:pressure_stabilization}{圧力の安定化をするかどうか}
- \ref{SPH:RK_order}{ルンゲクッタの段数}
- \ref{SPH:reflection}{反射の計算方法}

壁のwall_as_fluidは繰り返しで計算するのはどうか？

*/

#endif