#ifndef SPH_FindPressure_H
#define SPH_FindPressure_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

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
   for (const auto &ROW : points)
#pragma omp single nowait
   {
      double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
      ROW->PoissonRHS = 0;
      ROW->column_value.clear();
      Tddd pO_x, pO_b;
      auto pO = ROW;
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

      各粒子`ROW`に対して，方程式を作成する．

      まずは，\ref{SPH:whereToMakeTheEquation}{方程式を立てる位置を決める．}

      */

      double total_weight = 0, dP = 0., P_wall = 0.;

      /*DOC_EXTRACT SPH

      各粒子`ROW`が，流体か壁か補助粒子か水面かによって，方程式が異なる．

      |方程式|目的|
      |:---------|---|
      | IMPLEMENTED  \ref{SPH:PoissonEquation}{ポアソン方程式}              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
      | NOTIMPLEMENTED  \ref{SPH:ImpermeableCondition}{不透過条件}         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
      | NOTIMPLEMENTED  \ref{SPH:AtmosphericPressureCondition}{大気圧条件} | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない． |

      各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．

      */

      auto ImpermeableCondition = [&](const auto &B /*column id*/) {  // \label{SPH:ImpermeableCondition}
         // ROW->PoissonRHS -= (V_next(B) * Dot(b_vector(B), Normalize(pO->interpolated_normal_SPH)) * w_Bspline(Norm(pO_x - B->X), pO->radius_SPH));
         // auto coeff = V_next(B) * Dot(grad_w_Bspline(pO_x, B->X, pO->radius_SPH), Normalize(pO->interpolated_normal_SPH));  // こっちはOKだろう．
         // ROW->increment(B, coeff);

         // auto coeff = /*B->p_SPH **/ B->mass / B->rho * grad_w_Bspline(pO->X, B->X, pO->radius_SPH);  //\label{SPH:gradP3}0.34
         // ROW->increment(B, Dot(coeff,pO->interpolated_normal_SPH));

         // auto n = Normalize(pO_x - ROW->X);
         // // auto w = B->mass / rho_next(pO) * w_Bspline(Norm(pO_x - B->X), pO->radius_SPH);
         // ROW->PoissonRHS = -Dot(pO->DUDt_SPH,n);
         // auto coeff = B->mass * grad_w_Bspline(pO_x, B->X, pO->radius_SPH);  //\label{SPH:gradP1}0.2647
         // ROW->increment(B, Dot(coeff, n)/ (B->rho * B->rho));
         // ROW->increment(pO, Dot(coeff, n)/ (pO->rho * pO->rho));

         auto n = Normalize(pO->interpolated_normal_SPH);
         auto w = B->mass / rho_next(pO) * w_Bspline(Norm(pO_x - B->X), pO->radius_SPH);
         ROW->PoissonRHS -= Dot(ROW->DUDt_SPH, n);
         auto coeff = B->mass * grad_w_Bspline(pO_x, B->X, ROW->radius_SPH);  //\label{SPH:gradP1}0.2647
         ROW->increment(B, Dot(coeff, n) / (B->rho * B->rho));
         ROW->increment(pO, Dot(coeff, n) / (pO->rho * pO->rho));
      };

      auto AtmosphericPressureCondition = [&](const auto &p) {  //  \label{SPH:AtmosphericPressureCondition}
         // pの圧力を完全にゼロにする条件
         ROW->PoissonRHS = 0;
         ROW->increment(p, 1.);
      };

      auto EquivalentPressure = [&](const auto &p) {
         // pの圧力を完全にゼロにする条件
         ROW->PoissonRHS = 0;
         ROW->column_value = {{pO, 1.}, {ROW, -1.}};
      };

      // \label{SPH:PoissonEquation}
      auto PoissonEquation = [&](const auto &B /*column id*/, const double &coef = 1.) {
         if (Distance(pO, B) < pO->radius_SPH) {
            // if (!B->isAuxiliary)
            {
               ROW->PoissonRHS += V_next(B) * Dot(b_vector(B) - pO_b, grad_w_Bspline(pO_x, X_next(B), pO->radius_SPH)) * coef;  // \label{SPH:div_b_vector}
               // \label{SPH:pressure_stabilization}
               // if (pO->isFluid && !pO->isAuxiliary) {
               //    const double alpha = 0.1;
               //    ROW->PoissonRHS += alpha * (_WATER_DENSITY_ - pO->density_based_on_positions) / (dt * dt);
               // }
            }
            Aij = 2. * B->mass / rho_next(pO) * Dot_grad_w_Bspline_Dot(pO_x, X_next(B), pO->radius_SPH) * coef;  //\label{SPH:lapP}
            // for ISPH
            ROW->increment(pO, Aij / pO->rho);
            ROW->increment(B, -Aij / pO->rho);
            // for EISPH
            sum_Aij_Pj += Aij * B->p_SPH;
            sum_Aij += Aij;
         }
      };
      auto addPoissonEquation = [&]() {
         for (const auto &net : target_nets) {
            net->BucketPoints.apply(pO_x, pO->radius_SPH * 1.5, [&](const auto &B) {
               if (B->isCaptured) {
                  PoissonEquation(B);
#if defined(USE_SHARED_AUX)
                  if (B->isSurface) {
                     if (pO == B) {
                        for (const auto &AUX : B->auxiliaryPoints)
                           PoissonEquation(AUX);
                     }
                     //  else if (!pO->isSurface) {
                     //    for (const auto &AUX : B->auxiliaryPoints) {
                     //       PoissonEquation(AUX, AUX->volume * w_Bspline(Norm(AUX->X - AUX->X), AUX->radius_SPH));
                     //       net->BucketPoints.apply(AUX->X, AUX->radius_SPH, [&](const auto &C) {
                     //          PoissonEquation(C, C->volume * w_Bspline(Norm(C->X - AUX->X), AUX->radius_SPH));
                     //       });
                     //    }
                     // }
                  }
#endif
                  // for mapping to wall
                  total_weight += B->volume * w_Bspline(Norm(pO_x - X_next(B)), ROW->radius_SPH);
                  dP = Dot(X_next(ROW) - pO_x, B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
                  P_wall += (B->p_SPH + dP) * B->volume * w_Bspline(Norm(pO_x - X_next(B)), ROW->radius_SPH);
               }
            });
         }
#if defined(USE_SIMPLE_SINGLE_AUX)
         if (pO->isSurface)
            for (const auto &AUX : pO->auxiliaryPoints)
               PoissonEquation(AUX);
#endif
      };
      /*DOC_EXTRACT SPH

      ### ポアソン方程式の作成

      */

      // \label{SPH:whereToMakeTheEquation}
      if (ROW->isAuxiliary) {
         pO = ROW->surfacePoint;  // Aが安定する．
         // pO = ROW;  // Aが安定する．
         pO_x = X_next(pO);
         pO_b = b_vector(pO);
         AtmosphericPressureCondition(pO);
      } else if (ROW->getNetwork()->isRigidBody) {
         auto tmp = getClosestFluid(ROW, target_nets);
         if (tmp == nullptr)
            tmp = getClosestExcludeRigidBodyInlcudeFirstLayer(ROW, target_nets);
         if (tmp != nullptr)
            pO = tmp;
         pO_x = X_next(pO);
         pO_b = b_vector(pO);
         addPoissonEquation();
      } else {
         pO = ROW;
         pO_x = X_next(ROW);
         pO_b = b_vector(ROW);
         addPoissonEquation();
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
#endif
      //% ------------------------------------------------------- */
      // ROW->div_tmpU = ROW->PoissonRHS * dt / ROW->rho;
      // ROW->DrhoDt_SPH = -ROW->rho * ROW->div_tmpU;
      // ROW->rho_ = ROW->rho + ROW->DrhoDt_SPH * dt;

      ROW->p_SPH_ = (ROW->PoissonRHS + sum_Aij_Pj) / sum_Aij;

      if (ROW->getNetwork()->isRigidBody) {
         if (total_weight > 0.001)
            ROW->p_SPH_ = P_wall / total_weight;
         else
            ROW->p_SPH_ = 0;
      }

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

void setPressure(const std::unordered_set<networkPoint *> &points) {
   for (const auto &p : points)
      p->p_SPH = p->p_SPH_;
}

/*DOC_EXTRACT SPH

## ポアソン方程式の解法

ISPHのポアソン方程式を解く場合，\ref{SPH:gmres}{ここではGMRES法}を使う．

*/

// #define USE_LAPACK
#define USE_GMRES

void solvePoisson(const std::unordered_set<networkPoint *> &fluid_particle,
                  const std::unordered_set<networkPoint *> &wall_as_fluid,
                  const std::unordered_set<Network *> &target_nets) {

   std::unordered_set<networkPoint *> points;
   points.reserve(fluid_particle.size() + wall_as_fluid.size() + 1000);

   // 解く粒子の集合を保存

   for (const auto &p : fluid_particle) {
      points.emplace(p);
      if (p->isSurface)
         for (const auto &AUX : p->auxiliaryPoints)
            points.emplace(AUX);
   }

   for (const auto &p : wall_as_fluid)
      // if (p->isFirstWallLayer)
      points.emplace(p);

   /* -------------------------------------------------------------------------- */

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

      p->setVectorCSR();
   }
#if defined(USE_GMRES)

   gmres gm(points, b, x0, 100);  //\label{SPH:gmres}
   for (auto i = 1; i < 10; i++) {
      x0 = gm.x;
      std::cout << " gm.err : " << gm.err << std::endl;
      auto error = Norm(b - Dot(points, x0));
      std::cout << "actual error : " << error << std::endl;
      if (gm.err < 1E-5)
         break;
      gm.Restart(points, b, x0, 100);  //\label{SPH:gmres}
   }
   // gmres gm(ToVector(points), b, x0, 100);
   // std::cout << " gm.err : " << gm.err << std::endl;

   for (const auto &p : points)
      p->p_SPH = gm.x[p->getIndexCSR()];

   std::cout << " gm.err : " << gm.err << std::endl;
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

   auto error = Norm(b - Dot(points, x0));
   std::cout << "actual error : " << error << std::endl;

   if (!isFinite(error))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error is not a finite");
};

#endif