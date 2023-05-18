#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"

/*DOC_EXTRACT SPH

### CFL条件の設定

$\max({\bf u}) \Delta t \leq c_{v} h \cap \max({\bf a}) \Delta t^2 \leq c_{a} h$

を満たすように，毎時刻$\Delta t$を設定する．

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
   return p->X + c * p->interpolated_normal_SPH * p->radius_SPH / p->C_SML;
};

void setNormal_Surface(auto &net, const std::unordered_set<networkPoint *> &wall_p, const auto &RigidBodyObject) {

   DebugPrint("水粒子のオブジェクト外向き法線方向を計算", Green);
// refference: A. Krimi, M. Jandaghian, and A. Shakibaeinia, Water (Switzerland), vol. 12, no. 11, pp. 1–37, 2020.
#pragma omp parallel
   for (const auto &p : net->getPoints())
#pragma omp single nowait
   {
      /*DOC_EXTRACT SPH

      ### 法線方向の計算と水面の判定

      CHECKED 単位法線ベクトル: ${\bf n}_i = -{\rm Normalize}\left(\sum_j {\frac{m_j}{\rho_j} \nabla W_{ij} }\right)$

      */

      // 初期化
      p->COM_SPH.fill(0.);
      p->interpolated_normal_SPH_original.fill(0.);
      p->interpolated_normal_SPH_original_all.fill(0.);
      double total_vol = 0, w;

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
      p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH_original);
      p->interpolated_normal_SPH_all = Normalize(p->interpolated_normal_SPH_original_all);

      for (const auto &n : normal_of_near_wall_particle) {
         p->interpolated_normal_SPH = Normalize(Chop(p->interpolated_normal_SPH, n));
         p->interpolated_normal_SPH_all = Normalize(Chop(p->interpolated_normal_SPH_all, n));
      }

      if (!isFinite(p->interpolated_normal_SPH)) {
         p->interpolated_normal_SPH = {0., 0., 1.};
         p->interpolated_normal_SPH_all = {0., 0., 1.};
      }

      // 水面の判定
      /*DOC_EXTRACT SPH

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

   DebugPrint("水面粒子の作成", Green);
   for (const auto &p : net->getPoints()) {
      double d = 0;
      if (p->isSurface) {
         for (auto &auxp : p->auxiliaryPoints) {
            d += 1.;
            auxp = new networkPoint(net->surfaceNet, X_SPP(p, d));
            auxp->surfacePoint = p;
            auxp->isAuxiliary = true;
            auxp->p_SPH = p->p_SPH;
            auxp->U_SPH = p->U_SPH;
            auxp->rho = p->rho;
            auxp->volume = p->volume;
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

壁粒子の流速を流体粒子の流速に応じて変化させると計算が煩雑になるので，**ここでは**壁面粒子の流速は常にゼロに設定することにした（ゼロで一定というのは不自然ではない）．
一方，壁粒子の圧力がゼロだとするのは不自然で，流体粒子の圧力$p^{n+1}$の計算に悪影響を及ぼす．
なので．壁粒子の圧力は各ステップ毎に計算し直す必要がある．

壁面粒子の圧力は，壁面法線方向流速をゼロにするように設定されるべきだろう．

*/

#define Morikawa2019
#define new_method

/*DOC_EXTRACT SPH
### $\nabla^2 {\bf u}_i$の計算

CHECKED: ラプラシアンの計算方法: $\nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}$

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
      A->lap_U.fill(0.);
      //$ ------------------------------------------ */
      auto add = [&](const auto &B, const auto &qX, const double c = 1.) {
         const auto rij = qX - A->X;
#if defined(Morikawa2019)
         const auto Uij = A->U_SPH - c * B->U_SPH;
         A->lap_U += 2 * B->mass / A->rho * Uij * Dot_grad_w_Bspline_Dot(A->X, qX, A->radius_SPH);
#elif defined(Nomeritae2016)
         const auto Uij = c * B->U_SPH - A->U_SPH;
         const auto nu_nu = B->mu_SPH / B->rho + A->mu_SPH / A->rho;
         A->lap_U += 1 / (A->mu_SPH / A->rho) * B->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(A->X, qX, A->radius_SPH) /
                     ((B->rho + A->rho) * Dot(rij, rij));
#endif

         // 水面
         // if (B->isSurface)
         //    for (const auto &AUX : B->auxiliaryPoints)
         //       A->lap_U += 2 * B->mass / A->rho * Dot_grad_w_Bspline_Dot(A->X, AUX->X, A->radius_SPH) * Uij;

         // just counting
         if (Between(Norm(rij), {1E-12, A->radius_SPH})) {
            A->checked_points_in_radius_SPH++;
            if (B->getNetwork()->isFluid || B->isFluid)
               A->checked_points_in_radius_of_fluid_SPH++;
         }
         A->checked_points_SPH++;
      };
      //$ ------------------------------------------ */
      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) add(B, B->X);
         });
      //$ ------------------------------------------ */
      A->DUDt_SPH_ = A->DUDt_SPH;
      A->ViscousAndGravityForce = A->DUDt_SPH = A->lap_U * (A->mu_SPH / A->rho) + _GRAVITY3_;  // 後で修正されるDUDt
      A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
      A->tmp_X = A->X + A->tmp_U_SPH * dt;
   }
};
// b$ ------------------------------------------------------ */
// b$                     rho^*  の計算                       */
// b$ ------------------------------------------------------ */

void setTmpDensity(const std::unordered_set<networkPoint *> &points, const double dt) {
#if defined(Morikawa2019)
      /* use rho_ calculated at div_tmpU */
#elif defined(Nomeritae2016)
   for (const auto &p : points)
      p->rho_ = p->rho + (p->DrhoDt_SPH = -p->rho * p->div_tmpU) * dt;
#elif defined(Barcarolo2013)
   throw std::runtime_error("not implemented");
#endif
}

// b% -------------------------------------------------------------------------- */

/*DOC_EXTRACT SPH

### `PoissonRHS`,$b$と$\nabla^2 p^{n+1}$における$p^{n+1}$の係数の計算

次の時刻の流れ場が発散なし$\nabla\cdot{\bf u}^{n+1}=0$であることを保証してくれる圧力を使って，
$\frac{D {\bf u}}{D t} =-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}$を決定し，時間発展させたい．
そのような圧力を$p^{n+1}$と書くことにする．
そのような圧力の条件は，次のようになる．

$$
\begin{align*}
&&\frac{D {\bf u}}{D t} &=-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}\\
&\rightarrow& \frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t} &=-\frac{1}{\rho} \nabla p^{n+1}+\nu \nabla^2 {\bf u}^n+{\bf g}\\
&\rightarrow& \nabla \cdot\left(\frac{\rho}{\Delta t} {\bf u}^{n+1}\right) + \nabla^2 p^{n+1} &= \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}^n+\rho {\bf g}\right)\\
&\rightarrow& \nabla^2 p^{n+1} &= b, \quad b = \nabla \cdot {{\bf b}^n} = \nabla \cdot \left(\frac{\rho}{\Delta t} {\bf u}^n+\mu \nabla^2 {\bf u}+\rho {\bf g}\right)
\end{align*}
$$

この$b$を`PoissonRHS`とする．（仮流速は${\bf u}^* = \frac{\Delta t}{\rho}{\bf b}^n$である．）

CHECKED: 発散の計算方法: $b=\nabla\cdot{\bf b}^n=\sum_{j}\frac{m_j}{\rho_j}({\bf b}_j^n-{\bf b}_i^n)\cdot\nabla W_{ij}$

`PoissonRHS`,$b$の計算の前に，$\mu \nabla^2{\bf u}$を予め計算しておく．

壁粒子の圧力は時間発展させないので，壁粒子の$p^n$を計算する必要がある．順で計算する．

1. 壁粒子の圧力の計算（流体粒子の現在の圧力$p^n$だけを使って近似）
2. 流体粒子の圧力$p^{n+1}$の計算

CHECKED: ラプラシアンの計算方法: $`\nabla^2 p^{n+1}=\sum_{j}A_{ij}(p_i^{n+1} - p_j^{n+1}),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

*/

void PoissonEquation(const std::unordered_set<networkPoint *> &points,
                     const std::unordered_set<Network *> &target_nets,
                     const double dt, const double &particle_spacing,
                     const std::function<Tddd(const networkPoint *)> &getX,
                     const bool isWall = false) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
      A->PoissonRHS = 0;
      auto b = [&](const auto &A) { return (A->rho / dt * A->U_SPH) + (A->mu_SPH * A->lap_U) + (A->rho * _GRAVITY3_); };
      // auto b_ = [&](const auto &A) { return (A->mu_SPH * A->lap_U) + (A->rho * _GRAVITY3_); };
      auto b_at_X = b(A);
      auto a = A;
      //
      A->column_value.clear();
      auto markerX = getX(A);
      double total_weight = 0, P_wall = 0, dP;

      if (isWall) {
         // markerX += 1. * A->normal_SPH;  //\label{SPH:map_fluid_pressure_to_wall}
         markerX = A->X + 1.8 * A->normal_SPH;  //\label{SPH:map_fluid_pressure_to_wall}
         // markerX = 0.1 * A->normal_SPH;  //\label{SPH:map_fluid_pressure_to_wall}
      }
      if (A->isAuxiliary) {
         markerX = A->surfacePoint->X;  // + 0.00001 * particle_spacing * Normalize(A->surfacePoint->X - A->X);  //\label{SPH:map_fluid_pressure_to_wall}
         a = A->surfacePoint;
         b_at_X = b(a);
      }
      if (isWall || A->isAuxiliary) {
         b_at_X.fill(0.);
         auto process = [&](const auto &B /*use this property, p, U, rho, lap_U*/, const auto &qX /*use this position*/) {
            b_at_X += b(B) * B->volume * w_Bspline(Norm(qX - markerX), B->radius_SPH);
         };
         for (const auto &net : target_nets)
            net->BucketPoints.apply(markerX, a->radius_SPH, [&](const auto &B) {
               if (B->isCaptured) {
                  process(B, getX(B));
                  if (B->isSurface)
                     for (const auto &AUX : B->auxiliaryPoints)
                        process(B, AUX->X);
               }
            });
      }

      A->density_based_on_positions = 0;
      //% ----------------- PoissonRHS ------------------------- */
      // \label{SPH:PoissonEquation}
      auto PoissonEquation = [&](const auto &ID /*column id*/, const auto &B /*use this property, p, U, rho, lap_U*/, const auto &qX /*use this position*/) {
         A->PoissonRHS += B->volume * Dot(b(B) - b_at_X, grad_w_Bspline(markerX, qX, a->radius_SPH));
         // A->PoissonRHS += B->volume * Dot(b(B), grad_w_Bspline(markerX, qX, a->radius_SPH));
         A->density_based_on_positions += B->volume * w_Bspline(Norm(markerX - qX), a->radius_SPH);
#if defined(Morikawa2019)
         Aij = 2. * B->mass / a->rho * Dot_grad_w_Bspline_Dot(markerX, qX, a->radius_SPH);
#elif defined(Nomeritae2016)
         Aij = 2. * B->mass / std::pow((B->rho_ + A->rho_) / 2., 2.) * Dot_grad_w_Bspline_Dot(markerX, qX, A->radius_SPH) * A->rho_;
#elif defined(Barcarolo2013)
         Aij = 2. * B->mass / A->rho * Dot_grad_w_Bspline_Dot(markerX, qX, A->radius_SPH);
#endif
         sum_Aij_Pj += Aij * B->p_SPH;
         sum_Aij += Aij;
         A->increment(ID, -Aij);
         A->increment(A, Aij);
      };

      // \label{SPH:AuxiliaryEquation}
      auto AuxiliaryEquation = [&](const auto &ID /*column id*/, const auto &B /*use this property, p, U, rho, lap_U*/, const auto &qX /*use this position*/) {
         A->PoissonRHS += B->volume * Dot(/*b(B)*/ -b_at_X, grad_w_Bspline(markerX, qX, a->radius_SPH));
         // A->PoissonRHS += B->volume * Dot(b(B), grad_w_Bspline(markerX, qX, a->radius_SPH));
         A->density_based_on_positions += B->volume * w_Bspline(Norm(markerX - qX), a->radius_SPH);
#if defined(Morikawa2019)
         Aij = 2. * B->mass / a->rho * Dot_grad_w_Bspline_Dot(markerX, qX, a->radius_SPH);
#elif defined(Nomeritae2016)
         Aij = 2. * B->mass / std::pow((B->rho_ + A->rho_) / 2., 2.) * Dot_grad_w_Bspline_Dot(markerX, qX, A->radius_SPH) * A->rho_;
#elif defined(Barcarolo2013)
         Aij = 2. * B->mass / A->rho * Dot_grad_w_Bspline_Dot(markerX, qX, A->radius_SPH);
#endif
         sum_Aij_Pj += Aij * B->p_SPH;
         sum_Aij += Aij;
         A->increment(ID, -Aij);
         A->increment(A, Aij);
      };

      // \label{SPH:ImpermeableCondition}
      const double c = 1.;
      auto ImpermeableCondition = [&](const auto &ID /*column id*/, const auto &B /*use this property, p, U, rho, lap_U*/, const auto &qX /*use this position*/) {
         A->PoissonRHS -= c * (B->volume * Dot(b(B), Normalize(A->normal_SPH)) * w_Bspline(Norm(markerX - qX), a->radius_SPH));
         auto coeff = B->volume * Dot(grad_w_Bspline(markerX, qX, a->radius_SPH), Normalize(A->normal_SPH));  // こっちはOKだろう．
         A->increment(ID, c * coeff);
      };

      // \label{SPH:AtmosphericPressureCondition}
      auto AtmosphericPressureCondition = [&](const auto &ID /*column id*/, const auto &B /*use this property, p, U, rho, lap_U*/, const auto &qX /*use this position*/) {
         // A->PoissonRHS -= c * (B->volume * Dot(b(B), Normalize(A->normal_SPH)) * w_Bspline(Norm(markerX - qX), a->radius_SPH));
         A->PoissonRHS = 0;
         A->increment(ID, B->volume * w_Bspline(Norm(markerX - qX), a->radius_SPH));
      };

      /*DOC_EXTRACT SPH

      ### 圧力を決定するための方程式を作成

      各粒子$A$に対して，圧力を決定するための方程式を作成する．各粒子$A$が，流体か壁か補助粒子か水面かによって，方程式が異なる．

         - [x] \ref{SPH:PoissonEquation}{ポアソン方程式}　次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する．
         - [x] \ref{SPH:AuxiliaryEquation}{補助方程式}　水面上部に粒子を補い，水面での圧力が大気圧になるように圧力を決定する．
         - [] \ref{SPH:ImpermeableCondition}{不透過条件} この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない．
         - [] \ref{SPH:AtmosphericPressureCondition}{大気圧条件}　この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない．

      各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．

      */
      auto add = [&](const auto &B, const auto &qX) {
         if (A->isAuxiliary) {
            AtmosphericPressureCondition(B, B, qX);
            if (B->isSurface)
               for (const auto &AUX : B->auxiliaryPoints)
                  AtmosphericPressureCondition(AUX, B, AUX->X);
         } else {
            PoissonEquation(B, B, qX);
            if (B->isSurface)
               for (const auto &AUX : B->auxiliaryPoints)
                  AuxiliaryEquation(AUX, B, AUX->X);
         }
         // for mapping to wall
         total_weight += B->volume * w_Bspline(Norm(markerX - qX), A->radius_SPH);
         dP = Dot(getX(A) - markerX, B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
         P_wall += (B->p_SPH + dP) * B->volume * w_Bspline(Norm(markerX - qX), A->radius_SPH);
      };
      //% ------------------------------------------------------- */
      for (const auto &net : target_nets) {
         net->BucketPoints.apply(markerX, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured)
               add(B, getX(B));
         });
      }
      /* -------------------------------------------------------------------------- */
#if defined(Morikawa2019)
         /*DOC_EXTRACT SPH
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
      A->div_tmpU = A->PoissonRHS * dt / A->rho;
      A->DrhoDt_SPH = -A->rho * A->div_tmpU;
      A->rho_ = A->rho + A->DrhoDt_SPH * dt;

      A->p_SPH_ = (A->PoissonRHS + sum_Aij_Pj) / sum_Aij;

      if (isWall) {
         if (total_weight > 0.001)
            A->p_SPH_ = P_wall / total_weight;
         else
            A->p_SPH_ = 0;
      }
   };
};

/* -------------------------------------------------------------------------- */

void setPressure(const std::unordered_set<networkPoint *> &points) {
   for (const auto &p : points)
      p->p_SPH = p->p_SPH_;
}

/* -------------------------------------------------------------------------- */
/*DOC_EXTRACT SPH

### ISPH

TODO: ISPHの解がもとまらないのはなぜか？

- \ref{SPH:map_fluid_pressure_to_wall}{壁粒子の圧力を計算する位置には留意する}

*/
void solvePoisson(const std::unordered_set<networkPoint *> &fluid_particle,
                  const std::unordered_set<networkPoint *> &wall_as_fluid,
                  const std::unordered_set<Network *> &target_nets) {
   size_t i = 0;
   std::unordered_set<networkPoint *> points;
   points.reserve(fluid_particle.size() + wall_as_fluid.size() + 1000);

   DebugPrint("fluid_particle", Green);
   for (const auto &A : fluid_particle) {
      points.emplace(A);

      if (A->isSurface)
         for (const auto &AUX : A->auxiliaryPoints)
            points.emplace(AUX);
   }

   DebugPrint("wall_as_fluid", Green);
   for (const auto &p : wall_as_fluid)
      points.emplace(p);

   for (const auto &p : points)
      p->setIndexCSR(i++);

   // for (const auto &p : fluid_particle)
   //    if (p->isSurface) {
   //       auto max = 0;
   //       for (const auto &[_, v] : p->column_value)
   //          if (std::abs(v) > max)
   //             max = std::abs(v);
   //       p->column_value.clear();
   //       p->increment(p, max);
   //       p->PoissonRHS = 0.;
   //       p->p_SPH = 0;
   //    }

   //
   // for (const auto &p : fluid_particle)
   //    if (p->isSurface) {
   //       auto a0 = p->auxiliaryPoints[0];
   //       a0->column_value.clear();
   //       // std::cout << "max : " << max << std::endl;
   //       // a0->increment(p, -max);
   //       a0->increment(a0, 10.);
   //       p->PoissonRHS = 0.;
   //       p->p_SPH = 0;
   //    }

   V_d b, x0;
   b.resize(points.size());
   x0.resize(points.size());
   int count = 0;
   for (const auto &p : points) {
      if (p->column_value.empty())
         count++;
   }

   DebugPrint("set b and x0", Green);
   for (const auto &p : points) {
      b[p->getIndexCSR()] = p->PoissonRHS;
      x0[p->getIndexCSR()] = 0;  // p->p_SPH;
   }

   DebugPrint("gmres", Green);
   gmres gm(points, b, x0, 100);
   std::cout << "gm.err : " << gm.err << std::endl;

   for (const auto &p : points) {
      p->p_SPH = gm.x[p->getIndexCSR()];
      p->column_value.clear();
   }
};

/* -------------------------------------------------------------------------- */

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */
/*DOC_EXTRACT SPH

### 圧力勾配$\nabla p^{n+1}$の計算 -> ${D {\bf u}}/{Dt}$の計算

CHECKED: 勾配の計算方法: $\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}$

CHECKED: 勾配の計算方法: $\nabla p_i = \rho_i \sum_{j} m_j \left(p_j - p_i\right) \nabla W_{ij}$

CHECKED: 勾配の計算方法: $\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}$

*/
void gradP(const std::unordered_set<networkPoint *> &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->gradP_SPH.fill(0.);
      //% ------------------------------------------ */
      auto func = [&](const auto &B, const auto &qX, const double c = 1.) {
#if defined(Morikawa2019)
         // A->gradP_SPH += A->rho * B->mass * (c * B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, qX, A->radius_SPH);
         A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(A->X, qX, A->radius_SPH);
            // A->gradP_SPH += B->mass / A->rho * (B->p_SPH - A->p_SPH) * grad_w_Bspline(A->X, qX, A->radius_SPH);
#elif defined(Nomeritae2016)
         A->gradP_SPH += B->volume * c * B->p_SPH * grad_w_Bspline(A->X, qX, A->radius_SPH);
#elif defined(Barcarolo2013)
         A->gradP_SPH += A->rho * B->mass * (c * B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, qX, A->radius_SPH);
#endif

         if (B->isSurface)
            for (const auto &AUX : B->auxiliaryPoints) {
               // A->gradP_SPH += A->rho * B->mass * (c * AUX->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, AUX->X, A->radius_SPH);
               A->gradP_SPH += AUX->p_SPH * B->mass / B->rho * grad_w_Bspline(A->X, AUX->X, A->radius_SPH);
               // A->gradP_SPH += B->mass / A->rho * (AUX->p_SPH - A->p_SPH) * grad_w_Bspline(A->X, AUX->X, A->radius_SPH);
            }
      };

      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured)
               func(B, B->X);
         });
         //% ------------------------------------------ */
#if defined(Morikawa2019)
      A->DUDt_SPH -= A->gradP_SPH / A->rho;
#elif defined(Nomeritae2016)
      A->DUDt_SPH -= A->gradP_SPH / A->rho;
#endif

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
      p->RK_U.push(p->DUDt_SPH);  // 速度
      p->U_SPH = p->RK_U.getX();

      p->RK_X.push(p->U_SPH);  // 位置
      p->setXSingle(p->tmp_X = p->RK_X.getX());
      auto getX = [&](const auto &p) { return p->X; };
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
      const double reflection_factor = 1.;
      const double asobi = 0.;
      auto closest = [&]() {
         double distance = 1E+20;
         networkPoint *P = nullptr;
         for (const auto &[obj, poly] : RigidBodyObject) {
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
      while (isReflected && count++ < 30) {
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
                     p->RK_X.repush(p->U_SPH);  // 位置
                     p->setXSingle(p->tmp_X = p->RK_X.getX());
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
   // // // \label{SPH:update_density}
   // #pragma omp parallel
   //    for (const auto &A : points)
   // #pragma omp single nowait
   //    {
   //       A->div_U = 0.;
   //       A->rho_ = 0.;
   //       //% ----------------- div_U ------------------------- */
   //       auto add = [&](const auto &B, const auto &qX, const double coef = 1.) {
   //          A->div_U += B->volume * Dot(B->U_SPH - A->U_SPH, grad_w_Bspline(A->X, qX, A->radius_SPH));
   //          A->rho_ += B->volume * w_Bspline(Norm(A->X - qX), A->radius_SPH);
   //       };
   //       for (const auto &net : target_nets)
   //          net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
   //             if (B->isCaptured) {
   //                add(B, B->X);
   // #ifdef USE_SPP_Fluid
   //                if (B->isSurface && canSetSPP(target_nets, B)) add(B, SPP_X(B), SPP_p_coef);
   // #endif
   //             }
   //          });
   //       //% ------------------------------------------------ */
   //    };
   //    for (const auto &A : points) {
   //       A->DrhoDt_SPH = -A->rho * A->div_U;
   //       A->RK_rho.push(A->DrhoDt_SPH);  // 密度
   //       A->setDensity(A->RK_rho.getX());

   //       // A->setDensity(A->rho_);
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