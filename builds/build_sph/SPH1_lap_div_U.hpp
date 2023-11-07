#ifndef SPH_lap_div_U_H
#define SPH_lap_div_U_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT 0_2_0_lap_div_U

## 粘性項$`\nabla^2 {\bf u}_i`$の計算（`calcLaplacianU`）

SELECTED: \ref{SPH:lapU}{流速のラプラシアンの計算方法}: $`\nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

SELECTED: \ref{SPH:divU}{流速の発散の計算方法}: $`\nabla\cdot{\bf u}_i=\sum_{j}\frac{m_j}{\rho_j}({{\bf u}_j-{\bf u}_i}) \cdot\nabla W_{ij}`$

*/

auto calcLaplacianU(const auto &points, const std::unordered_set<Network *> &target_nets) {
   try {
#pragma omp parallel
      for (const auto &A : points)
#pragma omp single nowait
      {
#if defined(USE_RungeKutta)
         const double dt = A->RK_X.get_dt();
#elif defined(USE_LeapFrog)
         const double dt = A->LPFG_X.get_dt();
#endif

         A->checked_points_in_radius_SPH = A->checked_points_in_radius_of_fluid_SPH = A->checked_points_SPH = 0;
         A->div_U = 0.;
         A->div_U_next = 0.;
         A->lap_U.fill(0.);
         A->lap_U_next.fill(0.);
         A->convection_term.fill(0.);
         A->b_vector.fill(0.);
         A->U_XSPH.fill(0.);
         // A->grad_U.fill({0., 0., 0.});
         double Aij, Aij_next;
         Tddd Uij;
         //
         auto applyOverPoints = [&](const auto &equation) {
            const auto r = A->SML();
            for (const auto &net : target_nets) {
               net->BucketPoints.apply(A->X, 1.1 * r, [&](const auto &B) {
                  if (canInteract(A, B))
                     equation(B);
               });
            }
         };

         const double c_xsph = 0.001;
         double total_w = 0;
         auto add = [&](const auto &B) {
            Uij = A->U_SPH - B->U_SPH;
            auto w = B->volume * w_Bspline(Norm(A->X - B->X), A->SML());

            //
            if (A->isFluid && B->isFluid) {
               A->U_XSPH += c_xsph * (-Uij) * w;  //\label{SPH:U_XSPH}
               total_w += w;
            }
            //
            A->div_U += B->volume * Dot(-Uij, grad_w_Bspline(A, B));  //\label{SPH:divU}
            Aij = 2. * B->volume * Dot_grad_w_Bspline(A, B);          //\label{SPH:lapP1}

            // #if defined(USE_Laplacian_correction)
            //! 修正
            const auto DelX = (A->X - B->X);
            applyOverPoints([&](const auto &Q) {
               //
               for (int i = 0; i < 3; i++)
                  A->lap_U[i] -= Aij * Q->volume * Dot(DelX, grad_w_Bspline(A, Q)) * (Q->U_SPH[i] - A->U_SPH[i]);
               //
               // A->lap_U -= Aij * Q->volume * DelX * grad_w_Bspline(A, Q) * (Q->U_SPH - A->U_SPH);
               // A->lap_U_next -= Aij * V_next(Q) * (Q->U_SPH - A->U_SPH) * Dot(DelX, grad_w_Bspline(A, Q));
            });

            // const auto DelX = (A->X - B->X);
            // auto c = Dot(DelX, B->volume * grad_w_Bspline(A, B));
            // Aij *= (1. + c);

            //!
            // #endif
            A->lap_U += Aij * Uij;  //\label{SPH:lapU}
            // A->lap_U_next += Aij_next * Uij;  // Uij_nextはわからないので，Uijで代用
            // just counting
            if (Between(Distance(A, B), {1E-13, A->SML()})) {
               A->checked_points_in_radius_SPH++;
               if (B->getNetwork()->isFluid || B->isFluid)
                  A->checked_points_in_radius_of_fluid_SPH++;
            }
            A->checked_points_SPH++;
         };

         applyOverPoints(add);
         // if (total_w > 1E-10)
         //    A->U_XSPH /= total_w;

         auto add_next = [&](const auto &B) {
            Uij = U_next(A) - U_next(B);
            A->div_U_next += V_next(B) * Dot(-Uij, grad_w_Bspline_next(A, B));  //\label{SPH:divU}
            Aij = 2. * B->volume * Dot_grad_w_Bspline_next(A, B);               //\label{SPH:lapP1}
            A->lap_U_next += Aij * Uij;                                         //\label{SPH:lapU}
         };

         applyOverPoints(add_next);

         // A->lap_U = Dot(A->Mat_B, A->lap_U);
         //$ ------------------------------------------ */
         // \label{SPH:lapU_for_wall}
         // \label{SPH:Poisson_b_vector}
         if (A->getNetwork()->isRigidBody) {
            // A->DUDt_SPH_.fill(0.);
            A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt
            A->DUDt_SPH += A->U_XSPH / dt;
            A->rho = _WATER_DENSITY_;
            double nu = A->mu_SPH / A->rho;
            // A->DUDt_SPH.fill(0.);
            // A->tmp_U_SPH.fill(0.);
            A->tmp_X = A->X;
            A->DrhoDt_SPH = 0;
            A->b_vector.fill(0.);
            // if (A->isNeumannSurface)

            A->DrhoDt_SPH = -A->rho * A->div_U;

            auto rho = A->rho + A->DrhoDt_SPH * dt;

            // A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);                            // 最も自然な結果を返す
            A->b_vector = A->rho * A->U_SPH / dt + A->mu_SPH * A->lap_U + A->rho * _GRAVITY3_;  // 最も自然な結果を返す
            A->b_vector += A->rho * A->U_XSPH / dt;

         } else {
            A->DUDt_SPH_ = A->DUDt_SPH;
            A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt
            A->DUDt_SPH += A->U_XSPH / dt;

            A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
            A->tmp_X = A->X + A->tmp_U_SPH * dt;
            A->DrhoDt_SPH = -A->rho * A->div_U;

            // \label{SPH:how_to_set_fluid_b_vector}
            // A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
            A->b_vector = A->rho * A->U_SPH / dt + A->mu_SPH * A->lap_U_next + A->rho * _GRAVITY3_;
            A->b_vector += A->rho * A->U_XSPH / dt;
         }
      }

   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in calcLaplacianU");
   };
};

#endif