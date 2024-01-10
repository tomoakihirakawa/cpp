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

// #define USE_NEXT_POSITION true

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
      std::unordered_set<Network *> fluid_nets;
      for (const auto &n : target_nets)
         if (n->isFluid)
            fluid_nets.emplace(n);

#pragma omp parallel
      for (const auto &ROW : points)
#pragma omp single nowait
      {

#if defined(USE_RungeKutta)
         const double dt = ROW->RK_X.get_dt();
#elif defined(USE_LeapFrog)
         const double dt = ROW->LPFG_X.get_dt();
#endif

         /*DOC_EXTRACT 0_2_2_set_pressure_eq

         * `ROW`は，どの粒子も方程式を保存するかを表す．
         * `pO_center`は，圧力の方程式を立てる際の座標を表す（基本的には`ROW`の位置と同じ）．
         * `pO`は，影響半径などの情報として使う粒子を表す（基本的には`ROW`と同じ）．

         |方程式|目的|
         |:---------|---|
         | IMPLEMENTED  \ref{SPH:PoissonEquation}{ポアソン方程式}              | 次時刻の流速の発散をゼロにする（非圧縮性を満たす）ように圧力を決定する． |
         | NOTIMPLEMENTED  \ref{SPH:ImpermeableCondition}{不透過条件}         | この式は圧力勾配がそれ以外の力を打ち消すように圧力を決定する．壁面付近の圧力が滑らかにならないため使わない． |
         | NOTIMPLEMENTED  \ref{SPH:AtmosphericPressureCondition}{大気圧条件} | この式は水面粒子の圧力をゼロに固定する．圧力がゼロであるべき場所は水面から$h/2$上なので使わない． |

         各方程式は，`equation(列番号を指定する粒子ポインタ, 計算に使われる物性値を持つ粒子ポインタ, 方程式を立てる位置)`の形で使用する．

         */

         double sum_Aij = 0, sum_Aij_Pj = 0;
         ROW->PoissonRHS = 0;
         ROW->clearColumnValue();
         auto pO = ROW;
         Tddd pO_center = X_next(pO);
         Tddd pO_center_mirror = pO_center;
         double total_volume_w = 0, total_w = 0;
         T3Tddd M;  // correection matrix

         /* -------------------------------------------------------------------------- */

         const double r = pO->SML_next();
         auto applyOverPoints = [&r, &pO, &pO_center](const auto &equation, const std::unordered_set<Network *> NETS) {
            for (const auto &net : NETS) {
               net->BucketPoints.apply(pO_center, 1.2 * r, [&](const auto &B) {
                  if (canInteract(pO, B))
                     equation(B);
               });
            }
         };

         //! -------------------------------------------------------------------------- */

         auto average_surrounding_pressure = [&ROW, &pO, &pO_center, &total_volume_w, &total_w](const auto &B /*column id*/) {
            const auto r = pO->SML();
            // const auto r = 2. * pO->particle_spacing;
            const auto BX = X_next(B);
            if (Distance(pO_center, BX) < r)
               // if (B->isFluid && pO->surfacePoint != B) {
               if (!B->isSurface && !B->isAuxiliary) {
                  auto w = w_Bspline(Norm(BX - pO_center), r);
                  auto vol_w = V_next(B) * w;
                  total_volume_w += vol_w;
                  total_w += w;
                  ROW->increment(B, vol_w);
                  ROW->PoissonRHS = 0.;
               }
         };

         //% ------------------- 壁粒子の圧力の方程式 ------------------- */

         auto EISPH_wall_pressure = [&ROW, &pO, &pO_center, &pO_center_mirror, &total_volume_w, &total_w](const auto &B /*column id*/) {
            const auto r = pO->SML_next();
            // const auto r = 1.9999999999 * pO->particle_spacing;
            const auto BX = X_next(B);
            if (Distance(pO_center, BX) <= r) {
               auto w = w_Bspline(Norm(BX - pO_center), r);
               auto vol_w = V_next(B) * w;
               total_volume_w += vol_w;
               total_w += w;
               auto dP = Dot(pO_center_mirror - BX, B->mu_SPH * B->lap_U + rho_next(B) * _GRAVITY3_);
               // ROW->p_EISPH += (B->p_SPH_last + dP) * vol_w;  // auto dir = Projection(X_next(ROW) - BX, Normalize(ROW->v_to_surface_SPH));
               FusedMultiplyIncrement(vol_w, B->p_SPH_last + dP, ROW->p_EISPH);
            }
         };

         auto ISPH_wall_pressure = [&ROW, &pO, &pO_center, &pO_center_mirror, &total_volume_w, &total_w](const auto &B /*column id*/) {
            const auto r = pO->SML_next();
            // const auto r = 1.9999999999 * pO->particle_spacing;
            const auto BX = X_next(B);
            if (Distance(pO_center, BX) <= r && !B->isAuxiliary) {
               auto w = w_Bspline(Norm(BX - pO_center), r);
               auto vol_w = V_next(B) * w;
               total_volume_w += vol_w;
               total_w += w;
               auto dP = Dot(pO_center_mirror - BX, B->mu_SPH * B->lap_U + rho_next(B) * _GRAVITY3_);
               ROW->increment(B, vol_w);
               // ROW->PoissonRHS -= dP * vol_w;  // auto dir = Projection(X_next(ROW) - BX, Normalize(ROW->v_to_surface_SPH));
               FusedMultiplyIncrement(vol_w, -dP, ROW->PoissonRHS);
            }
         };

         //% ------------------ 流体粒子の圧力の方程式 ------------------ */
         // \label{SPH:PoissonEquation}

         std::vector<std::tuple<networkPoint *, std::array<double, 3>>> near_particles_and_V_GradW;

         auto make_near_particles_and_V_GradW = [&pO, &pO_center, &near_particles_and_V_GradW](const auto &Q) {
            if (Distance(pO_center, X_next(Q)) <= pO->SML_next() && pO != Q)
               near_particles_and_V_GradW.emplace_back(Q, V_next(Q) * grad_w_Bspline_next(pO, pO_center, Q));
         };

         auto PoissonEquation = [&ROW, &pO, &sum_Aij_Pj, &sum_Aij, &pO_center, &applyOverPoints, &target_nets, &near_particles_and_V_GradW](const auto &B /*column id*/) {
            if (pO != B) {
               const auto BX = X_next(B);
               const auto r = pO->SML_next();
               if (Distance(pO_center, BX) <= r) {
                  // const double Aij = 2. * V_next(B) * Dot_grad_w_Bspline_next(pO, pO_center, B);  //\label{SPH:lapP1}
                  double Aij;
                  // if (pO->isNearSurface) {
                  // Aij = 2. * V_next(B) * Dot_grad_w_Bspline(pO_center, X_next(B), pO->SML_next());  //! 水面付近は修正によってエラーが起きやすいのでそのままにする．
                  // } else
                  Aij = 2. * V_next(B) * Dot_grad_w_Bspline_next(pO, pO_center, B);  //! 内部は比較的安定しているようなので修正する
                  // Aij /= 2.;
                  // var_Eifen_が大きすぎる内部の点は，計算がこんなんだろうから，EISPHの近似で計算する
                  /* -------------------------------------------------------------------------- */
                  //! 修正
                  const auto DelX = (pO_center - BX);
                  double c, total_c = 0;  //! 毎回すると遅いので
                  if (!pO->is_grad_corr_M_next_singular)
                     for (const auto &[Q, vol_grad] : near_particles_and_V_GradW) {
                        total_c += (c = Dot(Aij * DelX, vol_grad));
                        ROW->increment(Q, -c);
                     }
                  ROW->increment(pO, Aij + total_c);
                  ROW->increment(B, -Aij);
                  /* --------------------------------------------------------------------------- */

                  // ROW->PoissonRHS += V_next(B) * Dot(B->b_vector - pO->b_vector, grad_w_Bspline_next(pO, pO_center, B));
                  FusedMultiplyIncrement(V_next(B), Dot(B->b_vector - pO->b_vector, grad_w_Bspline_next(pO, pO_center, B)), ROW->PoissonRHS);
                  // ROW->PoissonRHS += V_next(B) * Dot(rho_next(B) * U_next(B) / B->RK_U.get_dt() - rho_next(pO) * U_next(pO) / pO->RK_U.get_dt(),  grad_w_Bspline_next(pO, pO_center, B));

                  //% for EISPH
                  // sum_Aij_Pj += Aij * B->p_SPH_last;
                  FusedMultiplyIncrement(Aij, B->p_SPH_last, sum_Aij_Pj);
                  sum_Aij += Aij;
               }
               // PoissonEquationWithX(B, B->X);
            }
         };

         // b$ -------------------------------------------------------------------------- */
         // b$                         粒子の種類によって，圧力の方程式を変える．                    */
         // b$ -------------------------------------------------------------------------- */

         // \label{SPH:whereToMakeTheEquation}
         // if (ROW->isSurface || ROW->isAuxiliary) {
         if (ROW->isAuxiliary) {
            // if auxiliary just average value of surrounding fluid particles
            // ROW->pressure_equation_index = 1;
            /* -------------------------------------------------------------------------- */
            // pO = ROW->surfacePoint;
            // pO_center = X_next(ROW->surfacePoint);
            /* -------------------------------------------------------------------------- */
            applyOverPoints(average_surrounding_pressure, target_nets);
            if (total_volume_w != 0.)
               for (auto &[_, v] : ROW->column_value)
                  v /= total_volume_w;
            ROW->increment(ROW, -1.);
            /* -------------------------------------------------------------------------- */
            // pO = ROW->surfacePoint;
            // pO_center = X_next(pO);
            // ROW->pressure_equation_index = 3;
            // applyOverPoints(PoissonEquation, target_nets);
            // ROW->p_SPH = ROW->p_EISPH = (ROW->PoissonRHS + sum_Aij_Pj) / sum_Aij;
            /* -------------------------------------------------------------------------- */
            // ROW->CRS::set(ROW, 1.);
            // ROW->PoissonRHS = 0;
         } else if (ROW->isSurface_next) {
            // if ((ROW->isSurface && !ROW->isAuxiliary && ROW->auxPoint != nullptr) || (ROW->isSurface && ROW->isAuxiliary && ROW->surfacePoint != nullptr)) {
            // if ((ROW->isSurface && !ROW->isAuxiliary && ROW->auxPoint != nullptr) || (ROW->isSurface && ROW->isAuxiliary && ROW->surfacePoint != nullptr)) {
            // if ((ROW->isSurface && ROW->isAuxiliary && ROW->surfacePoint != nullptr)) {
            ROW->pressure_equation_index = 0;
            // b% EISPH
            ROW->p_SPH = ROW->p_EISPH = 0;
            // b@ ISPH
            ROW->clearColumnValue();
            ROW->CRS::set(ROW, 1.);
            ROW->PoissonRHS = 0;
         } else if (ROW->getNetwork()->isRigidBody /*&& !ROW->isFirstWallLayer*/) {
            //! 壁面の圧力はPoissonを解かない方が計算が安定するようだ
            //! 　流体に近く内部の壁面粒子だけにPoisson方程式を解かせたとしても安定しない
            // } else if (false) {
            /*DOC_EXTRACT 0_2_2_set_pressure_eq

            壁面粒子の圧力の設定方法

            ポアソン方程式を解いた場合：壁近傍の粒子が内部方向への圧力を受ける．

            壁の法線方向にある流体の圧力を，壁粒子の圧力とした場合（若干の修正をするが）：あまり力を受けない．

            */
            ROW->pressure_equation_index = 2;
            ROW->p_SPH = ROW->p_EISPH = 0;
            pO = ROW;
            pO_center_mirror = X_next(pO);
            //! pO_center_mirror: これは壁粒子自身の壁内の位置座標
            //! pO_center : よくマーカーと言われる流体内の位置座標
            //! ちょっとした修正11/28 -> 効果なし
            // if (ROW->isFirstWallLayer)
            //    pO_center = pO_center_mirror + pO->v_to_surface_SPH;
            // else
            pO_center = pO_center_mirror + 2. * pO->v_to_surface_SPH;

            //  b% EISPH (initial guess) EISPHは方程式を立てる必要がなく，直接圧力を計算する
            total_volume_w = 0;
            applyOverPoints(EISPH_wall_pressure, fluid_nets);
            if (total_volume_w != 0.)
               ROW->p_SPH = ROW->p_EISPH / total_volume_w;

            // b@ ISPH like EISPH．　ISPHは方程式を立てる必要がある
            total_volume_w = 0;
            applyOverPoints(ISPH_wall_pressure, fluid_nets);
            if (total_volume_w != 0.) {
               for (auto &[_, v] : ROW->column_value)
                  v /= total_volume_w;
               ROW->PoissonRHS /= total_volume_w;
            }
            ROW->increment(ROW, -1.);
         } else {
            // b@ ISPH
            ROW->pressure_equation_index = 3;
            pO = ROW;
            // if (ROW->isAuxiliary)
            //    pO = ROW->surfacePoint;

            pO_center = X_next(pO);
            applyOverPoints(make_near_particles_and_V_GradW, target_nets);
            applyOverPoints(PoissonEquation, target_nets);
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
            // const auto dt = pO->RK_X.get_dt();
            // const double alpha = ROW->SML() > 2.0 ? dt : 0.;
            // const double gamma = dt * 0.001;  // 逆効果2023/12/26
            // ROW->PoissonRHS = (1. - gamma) * ROW->PoissonRHS + gamma * (_WATER_DENSITY_ - ROW->intp_density) / std::pow(dt, 2);
            // auto DrhoDt = (ROW->intp_density_next - ROW->intp_density) / dt;
            // auto DrhoDt = /;
            // }
            // b% EISPH
            ROW->p_SPH = ROW->p_EISPH = (ROW->PoissonRHS + sum_Aij_Pj) / sum_Aij;
         }

         if (!isFinite(pO->PoissonRHS) || ROW->column_value.empty()) {
            ROW->clearColumnValue();
            ROW->CRS::set(ROW, 1.);
            ROW->PoissonRHS = 0;
         }

         /* -------------------------------------------------------------------------- */

         if (!isFinite(pO->PoissonRHS)) {
            // check type of particle
            std::cout << "pO->pressure_equation_index : " << ROW->pressure_equation_index << std::endl;
            std::cout << "pO->column_value.size() : " << pO->column_value.size() << std::endl;
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
            std::cout << "total_w : " << total_w << std::endl;
            std::cout << "pO->b_vector : " << pO->b_vector << std::endl;
            std::cout << "pO->lap_U : " << pO->lap_U << std::endl;
            std::cout << "pO->DUDt_SPH : " << pO->DUDt_SPH << std::endl;
            std::cout << "pO->U_SPH : " << pO->U_SPH << std::endl;
            std::cout << "pO->U_XSPH : " << pO->U_XSPH << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "pO->PoissonRHS is not a finite");
         }

         /* -------------------------------------------------------------------------- */

         // ROW->div_tmpU = ROW->PoissonRHS * dt / ROW->rho;
         // ROW->DrhoDt_SPH = -ROW->rho * ROW->div_tmpU;
         // ROW->rho_ = ROW->rho + ROW->DrhoDt_SPH * dt;

         if (ROW->column_value.empty()) {
            // show detail of this particle
            std::cout << "ROW->pressure_equation_index : " << ROW->pressure_equation_index << std::endl;
            std::cout << "sum_Aij_Pj : " << sum_Aij_Pj << std::endl;
            std::cout << "sum_Aij : " << sum_Aij << std::endl;
            std::cout << "total_volume_w : " << total_volume_w << std::endl;
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

void solvePoisson(const std::unordered_set<networkPoint *> &all_particle) {
   DebugPrint(Blue, "solvePoisson", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
   std::vector<networkPoint *> points;
   points.reserve(all_particle.size());

   for (const auto &p : all_particle)
      points.emplace_back(p);

   for (auto i = 0; const auto &p : points)
      p->setIndexCRS(i++);

   V_d b(points.size()), x0(points.size(), 0);

   for (const auto &p : points) {
      int index = p->getIndexCRS();
      if (p->column_value.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "empty column_value");
      if (!isFinite(p->PoissonRHS))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->PoissonRHS is not a finite");
      if (!isFinite(p->p_SPH)) {
         // print out information of p
         std::cout << "p->column_value.size() : " << p->column_value.size() << std::endl;
         std::cout << "p->vector_to_polygon_next.size() : " << p->vector_to_polygon_next.size() << std::endl;
         std::cout << "p->isFluid : " << p->isFluid << std::endl;
         std::cout << "p->isAuxiliary : " << p->isAuxiliary << std::endl;
         std::cout << "p->isSurface : " << p->isSurface << std::endl;
         std::cout << "p->getNetwork()->isRigidBody : " << p->getNetwork()->isRigidBody << std::endl;
         std::cout << "p->isFirstWallLayer : " << p->isFirstWallLayer << std::endl;
         std::cout << "p->PoissonRHS : " << p->PoissonRHS << std::endl;
         std::cout << "p->p_SPH : " << p->p_SPH << std::endl;
         std::cout << "p->X : " << p->X << std::endl;
         std::cout << "X_next(p) :" << X_next(p) << std::endl;
         std::cout << "p->b_vector : " << p->b_vector << std::endl;
         // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->p_SPH is not a finite");
         p->p_SPH = 0.;
      }
      b[index] = p->PoissonRHS;
      x0[index] = p->p_SPH;
   }

   DebugPrint(Blue, "preconditioning using diagonal value", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);

   /* ------------------ preconditioning using diagonal value ------------------ */

#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
   {
      double max = 0, value = 0;
      for (const auto &[_, v] : p->column_value)
         if (std::abs(v) > max) {
            max = std::abs(v);
            value = 1. / v;
         }

      // normalize
      b[p->getIndexCRS()] *= value;
      for (auto &[_, v] : p->column_value)
         v *= value;

      p->setVectorCRS();
   }

   /* -------------------------------------------------------------------------- */

#if defined(USE_GMRES)
   DebugPrint(Blue, "solve Poisson equation", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
   int size = 100;
   if (GMRES == nullptr) {
      GMRES = new gmres(points, b, x0, size);  //\label{SPH:gmres}
   } else {
      GMRES->Restart(points, b, x0, size);  //\label{SPH:gmres}
   }
   // gmres gm(points, b, x0, size);  //\label{SPH:gmres}

   DebugPrint(Blue, "solve Poisson equation", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
   x0 = GMRES->x;
   double torr = 1E-12;
   double error = GMRES->err;
   std::cout << Red << "       GMRES->err : " << GMRES->err << std::endl;
   std::cout << red << " actual error : " << (error = Norm(b_minus_A_dot_V(b, points, x0))) << std::endl;
   if (GMRES->err > torr)
      for (auto i = 1; i < 10; i++) {
         std::cout << "Restart : " << i << std::endl;
         GMRES->Restart(points, b, x0, size);  //\label{SPH:gmres}
         x0 = GMRES->x;
         std::cout << Red << "       GMRES->err : " << GMRES->err << std::endl;
         std::cout << red << " actual error : " << (error = Norm(b_minus_A_dot_V(b, points, x0))) << std::endl;
         if (GMRES->err < torr)
            break;
      }

   DebugPrint(Blue, "set pressure", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);

   for (const auto &p : points)
      p->p_SPH = p->p_SPH_last = GMRES->x[p->getIndexCRS()];

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

   DebugPrint(Blue, "solvePoisson end", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
};

#endif