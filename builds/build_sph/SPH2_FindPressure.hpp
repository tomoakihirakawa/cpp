#ifndef SPH_FindPressure_H
#define SPH_FindPressure_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT 0_2_1_SPH

## ポアソン方程式 $`\nabla^{n+1} \cdot \left(\frac{1}{\rho^n} \nabla^{n} p^{n+1}\right) = b`$

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

CHECKED: \ref{SPH:lapP1}{ラプラシアンの計算方法}: $`\nabla^2 p^{n+1}=\sum_{j}A_{ij}(p_i^{n+1} - p_j^{n+1}),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

CHECKED: \ref{SPH:lapP2}{ラプラシアンの計算方法}: $`\nabla^2 p^{n+1}=\sum_{j}A_{ij}(p_i^{n+1} - p_j^{n+1}),\quad A_{ij} = \frac{8 m_j\rho_i}{(\rho_i+\rho_j)^2}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

WARNING: 密度$\rho$が粒子に関わらず一定の場合，上の２式は同じになる．しかし，補助粒子の密度は，他の粒子と異なるので，\ref{SPH:lapP2}{２つ目のラプラシアンの計算方法}を使うべきだろう．

**ISPH**

   - ISPHは作ったポアソン方程式を作成し解くことで圧力を計算する

**EISPH**

   1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p^n`$だけを使って近似）
   2. 流体粒子の圧力$`p^{n+1}`$の計算

   \ref{SPH:EISPH_wall_pressure}{EISPHの圧力の設定方法}


$\sum_j A_{ij} (p_i^{n+1}-p_j^{n+1}) = b$において，$p_j^{n+1} \approx p_j^{n}$とすると，

```math
 p_i^{n+1} = \frac{b + \sum_j A_{ij} p_j^{n}}{\sum_j A_{ij}}
```

となる．

### 水面の計算補助粒子`auxiliaryPoints`

水面においては，流速の発散ゼロ$`\nabla^{n+1} {\bf u}^{n+1}=0`$と$`p^{n+1}=0`$が満たされる必要がある．
水面外部には，粒子がないので，求めた水面圧力は，ゼロであっても，圧力勾配は誤差を含み，$`\nabla^{n+1} {\bf u}^{n+1}=0`$は満足されない．
そこで，\ref{SPH:auxiliaryPoints}{水面の計算補助粒子}を水面外部に追加し，この点を適切計算することで，$`\nabla^{n+1} {\bf u}^{n+1}=0`$が満足されるように工夫する．

*/

/*DOC_EXTRACT 0_2_2_SPH

### 次時刻の発散演算，$`\nabla^{n+1} \cdot {\bf b}^n = \sum_j \dfrac{m_j}{\rho_j^{n+1}}({\bf b}_j^n-{\bf b}_i^n)\cdot \nabla W({\bf x}_i^{n+1},{\bf x}_j^{n+1},h)`$

$`\nabla^{n+1}`$の計算には，$`\rho^{n+1}`$, $`{\bf x}^{n+1}= {\bf x}^{n} + {\bf u}^{n+1} \Delta t`$が必要である．

* \ref{SPH:volume_next}{次時刻の粒子体積}
* \ref{SPH:rho_next}{次時刻の粒子密度}
* \ref{SPH:position_next}{次時刻の粒子位置}

*/

void setPressure(const std::unordered_set<networkPoint *> &points) {
   for (const auto &p : points)
      p->p_SPH = p->p_SPH_;
}

// \label{SPH:setPoissonEquation}
void setPoissonEquation(const std::unordered_set<networkPoint *> &points,
                        const std::unordered_set<Network *> &target_nets,
                        const double &particle_spacing) {
#pragma omp parallel
   for (const auto &ROW : points)
#pragma omp single nowait
   {

#if defined(USE_RungeKutta)
      const double dt = ROW->RK_X.get_dt();
#elif defined(USE_LeapFrog)
      const double dt = ROW->LPFG_X.get_dt();
#endif

      double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
      ROW->PoissonRHS = 0;
      ROW->clearColumnValue();
      Tddd pO_x, pO_b;
      auto pO = ROW;

      // \label{SPH:how_to_use_b_vector_in_Poisson0}
      auto b_vector = [&](const auto &B) {
         if (B->isAuxiliary)
            return B->surfacePoint->b_vector;
         else
            return B->b_vector;
      };

      double total_weight = 0, dP = 0., pressure_for_wall = 0.;

      /*DOC_EXTRACT 0_2_3_SPH

      各粒子`ROW`が，流体か壁か補助粒子か水面かによって，方程式が異なる．

      |方程式|目的|
      |:---------|---|
      | IMPLEMENTED  \ref{SPH:PoissonEquation}{ポアソン方程式}              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
      | NOTIMPLEMENTED  \ref{SPH:ImpermeableCondition}{不透過条件}         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
      | NOTIMPLEMENTED  \ref{SPH:AtmosphericPressureCondition}{大気圧条件} | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない． |

      各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．

      */

      auto ImpermeableCondition = [&](const auto &B /*column id*/) {  // \label{SPH:ImpermeableCondition}
         auto A = pO;

         // auto grad = Dot(grad_w_Bspline(A->X, B->X, A->radius_SPH), A->inv_grad_corr_M);
         // auto grad = grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);

         // auto grad = grad_w_Bspline(A->X, B->X, A->radius_SPH);
         auto grad = Dot(grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH), A->inv_grad_corr_M_next);
         auto cB = -B->mass * (1. / (rho_next(B) * rho_next(B))) * grad;
         auto cA = -B->mass * (1. / (rho_next(A) * rho_next(A))) * grad;

         auto n = Normalize(A->interp_normal_original_choped);
         ROW->increment(B, Dot(n, cB));  // [Newton]
         ROW->increment(A, Dot(n, cA));  // [Newton]
         // next force by pressure gradient

         // ROW->increment(B, Dot(Normalize(A->interp_normal_water), cB));
         // ROW->increment(A, Dot(Normalize(A->interp_normal_water), cA));

         // auto dt = ROW->LPFG_X.get_dt();
         // ROW->PoissonRHS = Dot(A->interp_normal_water, -A->DUDt_SPH);

         // auto nu_lap_U_g = A->DUDt_SPH;
         // auto nu_lap_U_g = 0.;
         ROW->PoissonRHS = 0;  // Dot(n, -nu_lap_U_g);
         // auto nu_lap_U_g = B->DUDt_SPH;
         // auto w = B->volume * w_Bspline(Norm(X_next(A) - X_next(B)), A->radius_SPH);
         // ROW->PoissonRHS += Dot(A->interp_normal_water, -nu_lap_U_g) * w;

         if (!isFinite(cB) || !isFinite(cA) || !isFinite(ROW->PoissonRHS))
            throw std::runtime_error("cB or cA is not finite.");
      };

      auto AtmosphericPressureCondition = [&](const auto &p) {  //  \label{SPH:AtmosphericPressureCondition}
         // pの圧力を完全にゼロにする条件
         ROW->PoissonRHS = 0;
         ROW->clearColumnValue();
         ROW->increment(p, 1.);
      };

      auto PoissonEquationWithX = [&](const auto &B /*column id*/, const Tddd &X, const bool isMirror = false) {
         if (Distance(pO_x, X) < pO->radius_SPH) {
            Aij = 2. * V_next(B) * Dot_grad_w_Bspline_Dot(pO_x, X, pO->radius_SPH);  //\label{SPH:lapP1}
            Aij /= pO->rho;
            ROW->increment(pO, Aij);
            ROW->increment(B, -Aij);
            for (const auto &[p, v] : pO->vector_p_grad) ROW->increment(p, -Aij * v);
            // \label{SPH:how_to_use_b_vector_in_Poisson1}
            ROW->PoissonRHS += V_next(B) * Dot(b_vector(B) - pO_b, grad_w_Bspline(pO_x, X, pO->radius_SPH));  // \label{SPH:div_b_vector}
            sum_Aij_Pj += Aij * B->p_SPH;
            sum_Aij += Aij;
         }
      };
      // \label{SPH:PoissonEquation}
      auto PoissonEquation = [&](const auto &B /*column id*/) {
         PoissonEquationWithX(B, X_next(B));
      };

      //\label{SPH:EISPH_wall_pressure}
      auto EISPH_wall_pressure = [&]() {
         total_weight = 0;
         pressure_for_wall = 0;
         double dP = 0, w = 0;
         auto markerX = ROW->X + 2 * ROW->normal_SPH;
         auto n = Normalize(ROW->normal_SPH);
         for (const auto &net : target_nets)
            if (!net->isRigidBody) {
               net->BucketPoints.apply(markerX, ROW->radius_SPH, [&](const auto &B) {
                  if (Distance(markerX, B->X) < ROW->radius_SPH) {
                     w = B->volume * w_Bspline(Norm(B->X - markerX), ROW->radius_SPH);
                     dP = Dot(Projection(ROW->X - B->X, n), B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
                     total_weight += w;
                     pressure_for_wall += (B->p_SPH + dP) * w;
                  }
               });
            }
      };

      auto ISPH_wall_pressure = [&](const auto &B /*column id*/, const double &coef = 1.) {
         auto X_next = [](const auto& p){return p->X;};
         auto V_next = [](const auto& p){return p->volume;};
         auto normal = ROW->normal_SPH;
         auto markerX = ROW->X + 2. * ROW->normal_SPH;
         auto r = ROW->radius_SPH;
         if (Distance(markerX, X_next(B)) < r) {
            auto w = V_next(B) * w_Bspline(Norm(X_next(B) - markerX), r);
            total_weight += w;
            ROW->increment(B, w);
            auto dir = Projection(X_next(ROW) - X_next(B), Normalize(normal));
            auto dP = Dot(dir, B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
            ROW->PoissonRHS -= dP * w;
         } };

      auto ISPH_wall_pressure_next = [&](const auto &B /*column id*/, const double &coef = 1.) {
         auto markerX = X_next(ROW) + 2. * ROW->normal_SPH;
         auto r = ROW->radius_SPH;
         if (Distance(markerX, X_next(B)) < r) {
            auto w = V_next(B) * w_Bspline(Norm(X_next(B) - markerX), r);
            total_weight += w;
            ROW->increment(B, w);
            auto dir = X_next(ROW) - X_next(B);
            // auto dir = Projection(X_next(ROW) - X_next(B), Normalize(ROW->normal_SPH));
            auto dP = Dot(dir, B->mu_SPH * B->lap_U + rho_next(B) * _GRAVITY3_);
            ROW->PoissonRHS -= dP * w;
         }
      };

      bool find_any_in_dist = false;
      double the_dist = ROW->radius_SPH / 2.5;

      auto addPoissonEquation = [&](const auto &PoissonEquation, const auto &pO_x, const std::unordered_set<Network *> NETS) {
         // find closest surface point
         networkPoint *closest_surface_point = nullptr;
         double min_distance = 1e10, distance;
         for (const auto &net : NETS) {
            net->BucketPoints.apply(pO_x, pO->radius_SPH, [&](const auto &B) {
               if (pO->isAuxiliary && pO->surfacePoint == B) return;
               if (B->isCaptured) {
                  PoissonEquation(B);
                  if (Distance(B, pO) < the_dist) find_any_in_dist = true;
               }
            });
         }
      };

      int type;

      // \label{SPH:whereToMakeTheEquation}
      if (ROW->isSurface) {
         type = 0;
         pO = ROW;
         pO_x = X_next(pO);
         AtmosphericPressureCondition(pO);
         // b% EISPH
         ROW->p_SPH = ROW->p_EISPH = 0;
         pO->p_SPH = pO->p_EISPH = 0;
      } else if (ROW->getNetwork()->isRigidBody) {
         type = 2;
         // b% EISPH
         EISPH_wall_pressure();
         if (ROW->getNetwork()->isRigidBody) {
            if (total_weight == 0.)
               ROW->p_SPH = ROW->p_EISPH = 0;
            else
               ROW->p_SPH = ROW->p_EISPH = pressure_for_wall / total_weight;
         }
         {
            type = 6;
            // b$ ISPH like EISPH
            pO = ROW;
            pO_x = pO->X + 2. * pO->normal_SPH;
            total_weight = 0;
            pressure_for_wall = 0;

            std::unordered_set<Network *> nets;
            for (const auto &n : target_nets)
               if (n->isFluid)
                  nets.emplace(n);

            addPoissonEquation(ISPH_wall_pressure_next, pO_x, nets);
            // addPoissonEquation(ISPH_wall_pressure, pO_x, nets);

            if (total_weight != 0.) {
               for (auto &[_, v] : ROW->column_value)
                  v /= total_weight;
               ROW->PoissonRHS /= total_weight;
            }

            ROW->increment(ROW, -1.);
         }
      } else {
         type = 7;
         // b$ ISPH
         pO = ROW;
         pO_x = X_next(pO);
         pO_b = b_vector(pO);

         addPoissonEquation(PoissonEquation, pO_x, target_nets);

         // const double alpha = 0.1 * dt;
         // if (ROW->isFluid)
         //    ROW->PoissonRHS += alpha * (_WATER_DENSITY_ - ROW->rho) / (dt * dt);

         // b% EISPH
         ROW->p_SPH = ROW->p_EISPH = (ROW->PoissonRHS + sum_Aij_Pj) / sum_Aij;
      }

      /* -------------------------------------------------------------------------- */

      if (!isFinite(pO->PoissonRHS)) {
         // check type of particle
         std::cout << "pO->vector_to_polygon_next.size() : " << pO->vector_to_polygon_next.size() << std::endl;
         std::cout << "pO->isFluid : " << pO->isFluid << std::endl;
         std::cout << "pO->isAuxiliary : " << pO->isAuxiliary << std::endl;
         std::cout << "pO->isSurface : " << pO->isSurface << std::endl;
         std::cout << "pO->getNetwork()->isRigidBody : " << pO->getNetwork()->isRigidBody << std::endl;
         std::cout << "pO->isFirstWallLayer : " << pO->isFirstWallLayer << std::endl;
         std::cout << "pO->PoissonRHS : " << pO->PoissonRHS << std::endl;
         std::cout << "pO->p_SPH : " << pO->p_SPH << std::endl;
         std::cout << "type : " << type << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "pO->PoissonRHS is not a finite");
      }

      /* -------------------------------------------------------------------------- */

#if defined(Morikawa2019)
         /* DOC_EXTRACT SPH

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
#endif
      //% ------------------------------------------------------- */
      // ROW->div_tmpU = ROW->PoissonRHS * dt / ROW->rho;
      // ROW->DrhoDt_SPH = -ROW->rho * ROW->div_tmpU;
      // ROW->rho_ = ROW->rho + ROW->DrhoDt_SPH * dt;

      if (ROW->column_value.empty()) {
         // show detail of this particle
         std::cout << "sum_Aij_Pj : " << sum_Aij_Pj << std::endl;
         std::cout << "sum_Aij : " << sum_Aij << std::endl;
         std::cout << "total_weight : " << total_weight << std::endl;
         std::cout << "ROW->PoissonRHS : " << ROW->PoissonRHS << std::endl;
         std::cout << "ROW->p_SPH_ : " << ROW->p_SPH_ << std::endl;
         std::cout << "ROW->isFluid : " << ROW->isFluid << std::endl;
         std::cout << "ROW->isAuxiliary :" << ROW->isAuxiliary << std::endl;
         std::cout << "ROW->isSurface :" << ROW->isSurface << std::endl;
         std::cout << "ROW->getNetwork()->isRigidBody :" << ROW->getNetwork()->isRigidBody << std::endl;
         std::cout << "ROW->isFirstWallLayer :" << ROW->isFirstWallLayer << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "empty column_value");
      }
   };
};

/*DOC_EXTRACT 0_2_4_SPH

## ポアソン方程式の解法

ISPHのポアソン方程式を解く場合，\ref{SPH:gmres}{ここではGMRES法}を使う．

*/

// #define USE_LAPACK
#define USE_GMRES

#if defined(USE_MIRROR_PARTICLE)
void solvePoisson(const std::unordered_set<networkPoint *> &fluid_particle)
#else
void solvePoisson(const std::unordered_set<networkPoint *> &fluid_particle,
                  const std::unordered_set<networkPoint *> &wall_as_fluid)

#endif
{
   std::vector<networkPoint *> points;

#if defined(USE_MIRROR_PARTICLE)
   points.reserve(fluid_particle.size() + 1000);
#else
   points.reserve(fluid_particle.size() + wall_as_fluid.size() + 1000);
#endif
   // 解く粒子の集合を保存

   for (const auto &p : fluid_particle)
      points.emplace_back(p);

#ifndef USE_MIRROR_PARTICLE
   for (const auto &p : wall_as_fluid)
      // if (p->isFirstWallLayer)
      points.emplace_back(p);
#endif

#if defined(USE_ONE_AUXP) || defined(USE_ALL_AUXP)
   for (const auto &p : fluid_particle)
      if (p->isSurface)
         for (const auto &AUX : p->auxiliaryPoints)
            if (AUX != nullptr)
               points.emplace_back(AUX);
#endif

   /* -------------------------------------------------------------------------- */

   for (auto i = 0; const auto &p : points)
      p->setIndexCSR(i++);

   V_d b(points.size()), x0(points.size(), 0);

   for (const auto &p : points) {
      int index = p->getIndexCSR();
      if (!isFinite(p->PoissonRHS))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->PoissonRHS is not a finite");
      if (!isFinite(p->p_SPH))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->p_SPH is not a finite");
      b[index] = p->PoissonRHS;
      x0[index] = p->p_SPH;
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

      p->setVectorCSR();
   }

#if defined(USE_GMRES)

   int size = 50;
   {
      gmres gm(points, b, x0, size);  //\label{SPH:gmres}
      std::cout << " gm.err : " << gm.err << std::endl;
      if (isFinite(gm.err))
         x0 = gm.x;
      else
         size = 100;
   }

   //
   gmres gm(points, b, x0, size);  //\label{SPH:gmres}
   x0 = gm.x;
   double error;
   for (auto i = 1; i < 3; i++) {
      if (gm.err < 1E-5)
         break;
      gm.Restart(points, b, x0, size);  //\label{SPH:gmres}
      x0 = gm.x;
      std::cout << "       gm.err : " << (error = gm.err) << std::endl;
      std::cout << " actual error : " << Norm(b - Dot(points, x0)) << std::endl;
   }

   for (const auto &p : points)
      p->p_SPH = gm.x[p->getIndexCSR()];

#elif defined(USE_LAPACK)
   VV_d A(b.size(), V_d(b.size(), 0.));
   for (const auto &p : points) {
      auto i = p->getIndexCSR();
      if (p->column_value.empty())
         throw std::runtime_error("empty column_value");
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

   if (!isFinite(error))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error is not a finite");
};

/* --------------- check if the Poisson equation is satisfied. -------------- */

void checkIfPoissonEqIsSatisfied(const auto &net, const auto &target_nets) {
   auto b_vector = [&](const auto &B) {
      if (B->isAuxiliary)
         return B->surfacePoint->b_vector;
      else
         return B->b_vector;
   };

   {
      double error_tot = 0;
      for (const auto &A : net->getPoints())
      // if (A->isSurface)
      {
         double eq = 0, Aij = 0, error = 0, Aij_tot = 0;
#if defined(USE_MIRROR_PARTICLE)
         net->BucketPoints.apply(X_next(A), A->radius_SPH, [&](const auto &B) {
            // 綺麗に配列させ設置した壁粒子の圧力を，内部の流体粒子の圧力から決める際に，理想的な値を設定することが困難と思われる．
            eq += V_next(B) * Dot(b_vector(B) - b_vector(A), grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH));
            Aij = 2. * V_next(B) * Dot_grad_w_Bspline_Dot(X_next(A), X_next(B), A->radius_SPH);
            Aij /= A->rho;
            eq -= -B->p_SPH * Aij;
            Aij_tot += Aij;
            for (const auto &v : B->vector_to_polygon_next) {
               auto X = X_next(B) + 2. * v;
               eq += V_next(B) * Dot(b_vector(B) - b_vector(A), grad_w_Bspline(X_next(A), X, A->radius_SPH));
               Aij = 2. * V_next(B) * Dot_grad_w_Bspline_Dot(X_next(A), X, A->radius_SPH);
               Aij /= A->rho;
               eq -= -B->p_SPH * Aij;
               Aij_tot += Aij;
            }
         });
#else
         for (const auto &net : target_nets) {
            net->BucketPoints.apply(X_next(A), A->radius_SPH, [&](const auto &B) {
               eq += V_next(B) * Dot(b_vector(B) - b_vector(A), grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH));
               Aij = 2. * V_next(B) * Dot_grad_w_Bspline_Dot(X_next(A), X_next(B), A->radius_SPH);
               Aij /= A->rho;
               eq -= -B->p_SPH * Aij;
               Aij_tot += Aij;
            });
         }
#endif
         A->p_SPH_ = eq / Aij_tot;
      };
      // これを満たしたとしても，制度が良くないため，非圧縮性を満たすことができない．
      // この方程式自体の精度を改善する必要がある．

      // for (const auto &A : net->getPoints())
      //    if (!A->isSurface)
      //       A->p_SPH = A->p_SPH_;
   }
   // {
   //    double error_tot = 0;
   //    for (const auto &A : net->getPoints())
   //       if (A->isSurface) {
   //          double eq = 0, Aij = 0, error = 0;
   //          for (const auto &net : _ALL_NET_) {
   //             net->BucketPoints.apply(X_next(A), A->radius_SPH, [&](const auto &B) {
   //                eq += V_next(B) * Dot(b_vector(B) - b_vector(A), grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH));
   //                Aij = 2. * V_next(B) * Dot_grad_w_Bspline_Dot(X_next(A), X_next(B), A->radius_SPH);
   //                Aij /= A->rho;
   //                if (A->isSurface)
   //                   eq -= A->auxiliaryPoints[0]->p_SPH * Aij;
   //                else
   //                   eq -= A->p_SPH * Aij;
   //                eq -= -B->p_SPH * Aij;
   //             });
   //          }
   //          error_tot += std::abs(eq);
   //       };
   //    std::cout << Red << "error_tot : " << error_tot << std::endl;
   // }
};

#endif