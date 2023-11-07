#ifndef SPH_FindPressure_H
#define SPH_FindPressure_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT 0_2_1_set_pressure_eq

## ポアソン方程式 $`\nabla ^{n+1} \cdot \left(\frac{1}{\rho ^n} \nabla ^{n} p \right)=b`$

### ポアソン方程式

次の時刻の流れ場を発散なし$`\nabla\cdot{\bf u}^{n+1}=0`$としてくれる
$`\frac{D {\bf u}}{D t} =-\frac{1}{\rho} \nabla p +\nu \nabla^2 {\bf u}^n+{\bf g}`$を使って，流速と粒子位置を時間発展させたい．
そのためには，圧力$`p^{n+1}`$を適切に決める必要がある．

$`\frac{D {\bf u}}{D t}`$は．$`\frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t}`$と離散化し条件を考えてみる．

```math
\rho\frac{{\bf u}^{n+1} - {\bf u}^{n}}{\Delta t} =- \nabla p +\mu \nabla^2 {\bf u}^n+\rho{\bf g}
```

次時刻の発散の演算は，次時刻における粒子配置に基づき行われるので，現在の粒子配置に基づく発散演算とは区別すべきである．
現在の微分演算を$`\nabla^{n}`$とし，次時刻の微分演算を$`\nabla^{n+1}`$とする．
$`\nabla^{n+1}`$を上の式に作用させると，

```math
\nabla^{n+1}\cdot {(\rho {\bf u}^{n+1})} = \nabla^{n+1} \cdot{(\rho {\bf u}^n)} - \Delta t \nabla^{n+1} \cdot\left( \nabla p-\mu \nabla^{n2} {\bf u}^n-{\rho\bf g}\right)
```

右辺がゼロとなれば，次時刻の流速の発散がゼロ，$`\nabla^{n+1} \cdot (\rho{\bf u}^{n+1})=0`$になる：

```math
\begin{align*}
&&0 &= \nabla^{n+1} \cdot (\rho{\bf u}^{n}) - \Delta t \nabla^{n+1} \cdot\left(\nabla p-\mu \nabla^{n2} {\bf u}^n-{\rho \bf g}\right)\\
&\rightarrow&\nabla^{n+1} \cdot \nabla p &= \frac{1}{\Delta t}\nabla^{n+1} \cdot (\rho {\bf u}^{n}) + \nabla^{n+1} \cdot\left(\mu \nabla^{n2} {\bf u}^n  + {\rho\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \nabla p &= \nabla^{n+1} \cdot\left(\frac{1}{\Delta t} (\rho {\bf u}^{n}) +\mu \nabla^{n2} {\bf u}^n  + {\rho\bf g}\right)\\
&\rightarrow& \nabla^{n+1} \cdot \nabla p & = \nabla^{n+1} \cdot {\bf b}^n= b,\quad  {\bf b}^n=\frac{1}{\Delta t}(\rho {\bf u}^{n}) +\mu \nabla^{n2} {\rho\bf u}^n
\end{align*}
```

重力の発散はゼロなので消した．
ここで，$`\nabla p`$は敢えて，微分演算子や圧力のタイムステップを書かなかった．
$`\nabla p`$は，次時刻の流速の発散をゼロにするためだけの未知ベクトルであって，タイムステップを考える必要はない．

ここでは，演算子が異なると複雑になるので，同じ$`\nabla^{n+1}`$を使うことにする．
これに伴って，次時刻の流速$`{\bf u}^{n+1}`$を計算する際に用いる圧力勾配は，$`\nabla^{n+1} p`$として計算しなければならないことに注意する．

次時刻の微分演算子を使うことにして，$`\nabla p`$を$`\nabla^{n+1} p`$とし，次のようなポアソン方程式を得る．

```math
\nabla^{n+1}\cdot \nabla^{n+1} p =b
```

ただし，水面粒子は$`p=0`$として，上の方程式は使わない．

### 右辺，$`b`$，`PoissonRHS`について

この$`b`$を`PoissonRHS`とする．（仮流速は$`{\bf u}^* = \frac{\Delta t}{\rho}{\bf b}^n`$と同じ）．
$`{\bf b}^n`$ （\ref{SPH:Poisson_b_vector}{`Poisson_b_vector`}）が計算できるように，$`{\bf u}^n`$と$`\nabla^2 {\bf u}^n`$を計算しておく．

CHECKED: \ref{SPH:div_b_vector}{発散の計算方法}: $`b=\nabla\cdot{\bf b}^n=\sum_{j}\frac{m_j}{\rho_j}({\bf b}_j^n-{\bf b}_i^n)\cdot\nabla W_{ij}`$

### 左辺について

壁粒子の圧力は時間積分して計算しないので，毎時刻，壁粒子の$`p`$を計算する必要がある．

CHECKED: \ref{SPH:lapP1}{ラプラシアンの計算方法}: $`\nabla^2 p=\sum_{j}A_{ij}(p_i - p_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

CHECKED: \ref{SPH:lapP2}{ラプラシアンの計算方法}: $`\nabla^2 p=\sum_{j}A_{ij}(p_i - p_j),\quad A_{ij} = \frac{8 m_j}{(\rho_i+\rho_j)}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

WARNING: 密度$\rho$が粒子に関わらず一定の場合，上の２式は同じになる．しかし，補助粒子の密度は，他の粒子と異なるので，\ref{SPH:lapP2}{２つ目のラプラシアンの計算方法}を使うべきだろう．

**ISPH**

   - ISPHは作ったポアソン方程式を作成し解くことで圧力を計算する

**EISPH**

   1. 壁粒子の圧力の計算（流体粒子の現在の圧力$`p`$だけを使って近似）
   2. 流体粒子の圧力$`p`$の計算

   \ref{SPH:EISPH_wall_pressure}{EISPHの圧力の設定方法}


$\sum_j A_{ij} (p_i-p_j) = b$において，$p_j^{\rm new} \approx p_j^{\rm old}$とすると，

```math
 p_i^{\rm new} = \frac{b + \sum_j A_{ij} p_j^{\rm old}}{\sum_j A_{ij}}
```

となる．

<!---
### 水面の計算補助粒子`auxiliaryPoints`

水面においては，流速の発散ゼロ$`\nabla^{n+1} {\bf u}^{n+1}=0`$と$`p^{n+1}=0`$が満たされる必要がある．
水面外部には，粒子がないので，求めた水面圧力は，ゼロであっても，圧力勾配は誤差を含み，$`\nabla^{n+1} {\bf u}^{n+1}=0`$は満足されない．
そこで，\ref{SPH:auxiliaryPoints}{水面の計算補助粒子}を水面外部に追加し，この点を適切計算することで，$`\nabla^{n+1} {\bf u}^{n+1}=0`$が満足されるように工夫する．
--->

*/

/*DOC_EXTRACT 0_2_2_SPH

### 次時刻の発散演算，$`\nabla^{n+1} \cdot {\bf b}^n = \sum_j \dfrac{m_j}{\rho_j^{n+1}}({\bf b}_j^n-{\bf b}_i^n)\cdot \nabla W({\bf x}_i^{n+1},{\bf x}_j^{n+1},h)`$

$`\nabla^{n+1}`$の計算には，$`\rho^{n+1}`$, $`{\bf x}^{n+1}= {\bf x}^{n} + {\bf u}^{n+1} \Delta t`$が必要である．

* \ref{SPH:volume_next}{次時刻の粒子体積}
* \ref{SPH:rho_next}{次時刻の粒子密度}
* \ref{SPH:position_next}{次時刻の粒子位置}

*/

#define USE_NEXT_POSITION true

void setPressure(const std::unordered_set<networkPoint *> &points) {
   for (const auto &p : points)
      p->p_SPH = p->p_SPH_;
}

/* -------------------------------------------------------------------------- */
/*                     setPoissonEquation A x = b                             */
/* -------------------------------------------------------------------------- */

// \label{SPH:setPoissonEquation}
void setPoissonEquation(const std::unordered_set<networkPoint *> &points,
                        const std::unordered_set<Network *> &target_nets,
                        const double &particle_spacing) {
   try {
      auto net = (*points.begin())->getNetwork();
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
         double total_weight = 0, dP = 0., pressure_for_wall = 0.;

         /*DOC_EXTRACT 0_2_2_set_pressure_eq

         各粒子`ROW`が，流体か壁か補助粒子か水面かによって，方程式が異なる．

         |方程式|目的|
         |:---------|---|
         | IMPLEMENTED  \ref{SPH:PoissonEquation}{ポアソン方程式}              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
         | NOTIMPLEMENTED  \ref{SPH:ImpermeableCondition}{不透過条件}         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
         | NOTIMPLEMENTED  \ref{SPH:AtmosphericPressureCondition}{大気圧条件} | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない． |

         各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．

         */

         auto applyOverPoints = [&](const auto &equation, const auto &pO_x, const std::unordered_set<Network *> NETS) {
            const auto r = USE_NEXT_POSITION ? ROW->SML_next() : ROW->SML();
            for (const auto &net : NETS) {
               net->BucketPoints.apply(pO_x, 1.2 * r, [&](const auto &B) {
                  if (B->isCaptured) {
                     // isSurface = true is just a reference point
                     if (pO->isSurface)
                        return;  // simple equation
                     else {
                        if (pO->isAuxiliary && !B->isSurface)
                           equation(B);
                        else if (!pO->isAuxiliary && !B->isAuxiliary)
                           equation(B);
                     }
                  }
               });
            }
         };

         /* ------------------- 壁粒子の圧力の方程式 ------------------- */

         auto ImpermeableCondition = [&](const auto &B /*column id*/) {  // \label{SPH:ImpermeableCondition}
            auto A = pO;
            auto c = rho_next(A) * B->mass;
            auto grad = grad_w_Bspline(X_next(A) + A->normal_SPH, X_next(B), A->SML_next());
            grad = Dot(grad, A->inv_grad_corr_M_next);
            auto cB = (c / std::pow(rho_next(B), 2)) * grad;
            auto cA = (c / std::pow(rho_next(A), 2)) * grad;
            //
            auto n = Normalize(A->interp_normal);
            ROW->increment(B, Dot(n, cB));  // [Newton]
            ROW->increment(A, Dot(n, cA));  // [Newton]
            ROW->PoissonRHS = 0;            // Dot(n, -nu_lap_U_g);
         };

         auto EISPH_wall_pressure = [&](const auto &B /*column id*/, const double &coef = 1.) {
            const auto markerX = pO_x;  // X_next(ROW) + 2. * ROW->v_to_surface_SPH;
            const auto r = USE_NEXT_POSITION ? ROW->SML_next() : ROW->SML();
            const auto BX = USE_NEXT_POSITION ? X_next(B) : B->X;
            const auto ROWX = USE_NEXT_POSITION ? X_next(ROW) : ROW->X;
            if (Distance(markerX, BX) < r) {
               auto w = B->volume * w_Bspline(Norm(BX - markerX), r);
               total_weight += w;
               auto dP = Dot(ROWX - BX, B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
               ROW->p_EISPH += (B->p_SPH + dP) * w;
               // auto dir = Projection(X_next(ROW) - BX, Normalize(ROW->v_to_surface_SPH));
            }
         };

         auto ISPH_wall_pressure = [&](const auto &B /*column id*/, const double &coef = 1.) {
            auto markerX = pO_x;  // X_next(ROW) + 2. * ROW->v_to_surface_SPH;
            const auto r = USE_NEXT_POSITION ? ROW->SML_next() : ROW->SML();
            const auto BX = USE_NEXT_POSITION ? X_next(B) : B->X;
            const auto ROWX = USE_NEXT_POSITION ? X_next(ROW) : ROW->X;
            if (Distance(markerX, BX) < r) {
               auto w = B->volume * w_Bspline(Norm(BX - markerX), r);
               total_weight += w;
               auto dP = Dot(ROWX - BX, B->mu_SPH * B->lap_U + B->rho * _GRAVITY3_);
               ROW->increment(B, w);
               ROW->PoissonRHS -= dP * w;  // auto dir = Projection(X_next(ROW) - BX, Normalize(ROW->v_to_surface_SPH));
            }
         };

         /* ------------------ 流体粒子の圧力の方程式 ------------------ */
         // \label{SPH:PoissonEquation}
         auto PoissonEquation = [&](const auto &B /*column id*/) {
            // 上の行の意味は，普通の流体粒子の圧力の計算には，補助粒子は関係ない．
            // しかし，補助粒子の圧力の計算には，普通の流体粒子は関係ある．

            const auto BX = USE_NEXT_POSITION ? X_next(B) : B->X;
            const auto r = USE_NEXT_POSITION ? ROW->SML_next() : ROW->SML();
            if (Distance(pO_x, BX) < r) {
               //\label{SPH:lapP1}
               // Aij = 2. * B->volume * Dot_grad_w_Bspline_Dot(pO_x, BX, r, USE_NEXT_POSITION ? pO->inv_grad_corr_M_next : pO->inv_grad_corr_M);
               Aij = 2. * B->volume * Dot_grad_w_Bspline(pO_x, BX, r, USE_NEXT_POSITION ? pO->inv_grad_corr_M_next : pO->inv_grad_corr_M);

               // //! 修正
               // auto rij = (pO_x - BX);
               // applyOverPoints([&](const auto &Q) {
               //    const auto QX = USE_NEXT_POSITION ? X_next(Q) : Q->X;
               //    auto rij_dot_p = Aij * Q->volume * Dot(rij, grad_w_Bspline(pO_x, QX, r, USE_NEXT_POSITION ? pO->inv_grad_corr_M_next : pO->inv_grad_corr_M));
               //    ROW->increment(pO, rij_dot_p);
               //    ROW->increment(Q, -rij_dot_p);
               // },
               //                 pO_x, target_nets);
               // //!
               ROW->increment(pO, Aij);
               ROW->increment(B, -Aij);
               //
               ROW->PoissonRHS += B->volume * Dot(B->b_vector - pO->b_vector, grad_w_Bspline(pO_x, BX, r, USE_NEXT_POSITION ? pO->inv_grad_corr_M_next : pO->inv_grad_corr_M));  // \label{SPH:div_b_vector}

               // for (auto i = 0; i < 3; ++i) {
               //    auto c = ROW->rho * B->mass;
               //    auto grad = USE_NEXT_POSITION ? grad_w_Bspline_next(ROW, B) : grad_w_Bspline(ROW, B);
               //    auto tmp = B->b_vector[i] * (c / std::pow(B->rho, 2)) * grad;
               //    tmp += ROW->b_vector[i] * (c / std::pow(ROW->rho, 2)) * grad;
               //    ROW->PoissonRHS += tmp[i];
               // }

               // \label{SPH:how_to_use_b_vector_in_Poisson1}
               // ROW->PoissonRHS += B->volume * Dot(B->b_vector - pO->b_vector, grad_w_Bspline(pO_x, X, r));  // \label{SPH:div_b_vector}
               // ROW->PoissonRHS += B->volume * Dot(B->b_vector - pO->b_vector, Dot(grad_w_Bspline(pO_x, X, r), pO->inv_grad_corr_M_next));  // \label{SPH:div_b_vector}
               // ROW->PoissonRHS += B->volume * Dot(B->b_vector - pO->b_vector, grad_w_Bspline(pO_x, X, r));  // \label{SPH:div_b_vector}
               sum_Aij_Pj += Aij * B->p_SPH;
               sum_Aij += Aij;
               // if (pO->isAuxiliary) {
               //    std::cout << "auxiliary found\n is empty now? " << ROW->column_value.empty() << ", Aij = " << Aij << std::endl;
               // }

               //%鏡写しにする．
               // if (ROW->isNotSurfaceButNearSurface && B->isSurface) {
               //    auto X = Mirror(pO_x, BX, pO_x - BX);
               //    auto Aij = 2. * B->volume * Dot_grad_w_Bspline(pO_x, X, r, USE_NEXT_POSITION ? pO->inv_grad_corr_M_next : pO->inv_grad_corr_M);
               //    ROW->increment(pO, Aij);
               //    ROW->increment(B, -Aij);
               //    // 鏡などで差はゼロになる．
               //    //  ROW->PoissonRHS += B->volume * Dot(pO->b_vector - pO->b_vector, grad_w_Bspline(pO_x, X, r, USE_NEXT_POSITION ? pO->inv_grad_corr_M_next : pO->inv_grad_corr_M));
               //    sum_Aij_Pj += Aij * B->p_SPH;
               //    sum_Aij += Aij;
               // }
               //%
            }
            // PoissonEquationWithX(B, B->X);
         };

         /* ---------------------------------------------------------- */
         //
         int type;

         std::unordered_set<Network *> fluid_nets;
         for (const auto &n : target_nets)
            if (n->isFluid)
               fluid_nets.emplace(n);
         // \label{SPH:whereToMakeTheEquation}
         if (ROW->isSurface) {
            ROW->pressure_equation_index = type = 0;
            // b% EISPH
            ROW->p_SPH = ROW->p_EISPH = 0;
            // b@ ISPH
            ROW->clearColumnValue();
            ROW->CRS::set(ROW, 1.);
            ROW->PoissonRHS = 0;
         } else if (ROW->getNetwork()->isRigidBody /* && !ROW->isFirstWallLayer*/) {
            ROW->pressure_equation_index = type = 2;
            ROW->p_SPH = ROW->p_EISPH = 0;
            pO = ROW;
            pO_x = (USE_NEXT_POSITION ? X_next(pO) : pO->X) + 2. * pO->v_to_surface_SPH;
            //  b% EISPH (initial guess) EISPHは方程式を立てる必要がなく，直接圧力を計算する
            total_weight = 0;
            applyOverPoints(EISPH_wall_pressure, pO_x, fluid_nets);
            if (total_weight == 0.)
               ROW->p_SPH = ROW->p_EISPH = 0;
            else
               ROW->p_SPH = ROW->p_EISPH / total_weight;

            // b@ ISPH like EISPH．　ISPHは方程式を立てる必要がある
            total_weight = 0;
            applyOverPoints(ISPH_wall_pressure, pO_x, fluid_nets);

            if (total_weight != 0.) {
               for (auto &[_, v] : ROW->column_value)
                  v /= total_weight;
               ROW->PoissonRHS /= total_weight;
            }
            ROW->increment(ROW, -1.);
         } else {
            // b@ ISPH
            ROW->pressure_equation_index = type = 3;
            pO = ROW;
            pO_x = USE_NEXT_POSITION ? X_next(pO) : pO->X;
            applyOverPoints(PoissonEquation, pO_x, target_nets);
            //! 安定化
            /*DOC_EXTRACT 0_2_3_set_pressure_eq

            ### 圧力の安定化

            $`b = \nabla \cdot {{\bf b}^n} + \alpha \frac{\rho_w - \rho^*}{{\Delta t}^2}`$として計算を安定化させる場合がある．
            $`\rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t`$と近似すると，

            ```math
            \rho^\ast = \rho + \frac{D\rho^\ast}{Dt}\Delta t,\quad
            \frac{D\rho^\ast}{Dt} = - \rho \nabla\cdot{\bf u}^\ast,\quad
            \nabla\cdot{\bf u}^\ast = \frac{\Delta t}{\rho} \nabla\cdot{\bf b}^n
            ```

            であることから，$`(\rho_w - \rho^*) / {\Delta t^2}`$は，$`\nabla\cdot{\bf b}^n`$となって同じになる．

            しかし，実際には，$`\rho^*`$は，$\nabla \cdot {{\bf b}^n}`$を使わずに，つまり発散演算を行わずに評価するので，
            計算上のようにはまとめることができない．

            $`\rho^*`$を計算する際に，$`\rho^\ast = \rho_w + \frac{D\rho^\ast}{Dt}\Delta t`$を使った場合，確かに上のようになるが，
            実際に粒子を仮位置に移動させその配置から$\rho^*$を計算した場合は，数値計算上のようにまとめることはできない．

            `PoissonRHS`,$`b`$の計算方法と同じである場合に限る．
            もし，計算方法が異なれば，計算方法の違いによって，安定化の効果も変わってくるだろう．

            */
            // if (!ROW->isSurface) {
            const auto dt = pO->RK_X.get_dt();
            // const double alpha = ROW->SML() > 2.0 ? dt : 0.;
            const double alpha = dt * 0.01;
            // ROW->PoissonRHS += alpha * (ROW->intp_density - ROW->intp_density_next) / std::pow(dt, 2);
            auto DrhoDt = (ROW->intp_density_next - ROW->intp_density) / dt;
            // ROW->PoissonRHS -= alpha * DrhoDt / dt;
            // }
            // b% EISPH
            ROW->p_SPH = ROW->p_EISPH = (ROW->PoissonRHS + sum_Aij_Pj) / sum_Aij;
         }

         /* -------------------------------------------------------------------------- */

         if (!isFinite(pO->PoissonRHS)) {
            // check type of particle

            bool find = false;
            applyOverPoints([&](const auto &p) {if (p->isCaptured && Distance(p->X, pO->X) < pO->SML()) find = true; }, ROW->X, fluid_nets);
            std::cout << "find : " << find << std::endl;

            std::cout << "pO->vector_to_polygon_next.size() : " << pO->vector_to_polygon_next.size() << std::endl;
            std::cout << "pO->isFluid : " << pO->isFluid << std::endl;
            std::cout << "pO->isAuxiliary : " << pO->isAuxiliary << std::endl;
            std::cout << "pO->isSurface : " << pO->isSurface << std::endl;
            std::cout << "pO->getNetwork()->isRigidBody : " << pO->getNetwork()->isRigidBody << std::endl;
            std::cout << "pO->isFirstWallLayer : " << pO->isFirstWallLayer << std::endl;
            std::cout << "pO->PoissonRHS : " << pO->PoissonRHS << std::endl;
            std::cout << "pO->p_SPH : " << pO->p_SPH << std::endl;
            std::cout << "pO->X : " << pO->X << std::endl;
            std::cout << "X_next(pO) :" << X_next(pO) << std::endl;
            std::cout << "type : " << type << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "pO->PoissonRHS is not a finite");
         }

         /* -------------------------------------------------------------------------- */

         // ROW->div_tmpU = ROW->PoissonRHS * dt / ROW->rho;
         // ROW->DrhoDt_SPH = -ROW->rho * ROW->div_tmpU;
         // ROW->rho_ = ROW->rho + ROW->DrhoDt_SPH * dt;

         if (ROW->column_value.empty()) {

            bool find = false;
            applyOverPoints([&](const auto &p) {if (p->isCaptured && Distance(p->X, pO->X) < pO->SML()) find = true; }, ROW->X, fluid_nets);
            std::cout << "find : " << find << std::endl;

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
            std::cout << "ROW->X :" << ROW->X << std::endl;
            std::cout << "X_next(ROW) :" << X_next(ROW) << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "empty column_value");
         }
      };
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in setPoissonEquation");
   };
};

/* -------------------------------------------------------------------------- */
/*                                solvePoisson                                */
/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_2_4_set_pressure_eq

## ポアソン方程式の解法

ISPHのポアソン方程式を解く場合，\ref{SPH:gmres}{ここではGMRES法}を使う．

*/

// #define USE_LAPACK
#define USE_GMRES

#if defined(USE_GMRES)
gmres<std::vector<networkPoint *>> *GMRES = nullptr;
#endif

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
      p->setIndexCRS(i++);

   V_d b(points.size()), x0(points.size(), 0);

   for (const auto &p : points) {
      int index = p->getIndexCRS();
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
      b[p->getIndexCRS()] /= max;
      for (auto &[_, v] : p->column_value)
         v /= max;

      p->setVectorCRS();
   }

#if defined(USE_GMRES)

   int size = 70;
   if (GMRES == nullptr)
      GMRES = new gmres(points, b, x0, size);  //\label{SPH:gmres}
   else
      GMRES->Restart(points, b, x0, size);  //\label{SPH:gmres}
   // gmres gm(points, b, x0, size);  //\label{SPH:gmres}

   x0 = GMRES->x;
   double torr = 1E-10;
   double error = GMRES->err;
   std::cout << Red << "       GMRES->err : " << (error = GMRES->err) << std::endl;
   std::cout << red << " actual error : " << Norm(b_minus_A_dot_V(b, points, x0)) << std::endl;
   if (error > torr)
      for (auto i = 1; i < 5; i++) {
         std::cout << "Restart : " << i << std::endl;
         GMRES->Restart(points, b, x0, size);  //\label{SPH:gmres}
         x0 = GMRES->x;
         std::cout << Red << "       GMRES->err : " << (error = GMRES->err) << std::endl;
         std::cout << red << " actual error : " << Norm(b_minus_A_dot_V(b, points, x0)) << std::endl;
         if (GMRES->err < torr)
            break;
      }

   for (const auto &p : points)
      p->p_SPH = GMRES->x[p->getIndexCRS()];

#elif defined(USE_LAPACK)
   VV_d A(b.size(), V_d(b.size(), 0.));
   for (const auto &p : points) {
      auto i = p->getIndexCRS();
      if (p->column_value.empty())
         throw std::runtime_error("empty column_value");
      for (const auto &[q, v] : p->column_value) {
         auto j = q->getIndexCRS();
         A[i][j] = v;
      }
   }
   lapack_lu lu(A);
   lu.solve(b, x0);
   for (const auto &p : points)
      p->p_SPH = x0[p->getIndexCRS()];
#elif defined(USE_LAPACK_SVD)
   VV_d A(b.size(), V_d(b.size(), 0.));
   for (const auto &p : points) {
      auto i = p->getIndexCRS();
      for (const auto &[q, v] : p->column_value) {
         auto j = q->getIndexCRS();
         A[i][j] = v;
      }
   }
   lapack_svd svd(A);
   svd.solve(b, x0);
   for (const auto &p : points)
      p->p_SPH = x0[p->getIndexCRS()];
#endif

   if (!isFinite(error))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error is not a finite");
};

#endif