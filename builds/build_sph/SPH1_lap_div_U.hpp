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
         A->lap_U.fill(0.);
         A->b_vector.fill(0.);
         A->U_XSPH.fill(0.);
         // A->grad_U.fill({0., 0., 0.});
         double Aij;
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

         double total_w = 0, w;
         auto add = [&](const auto &B) {
            Uij = A->U_SPH - B->U_SPH;
            w = B->volume * w_Bspline(Norm(A->X - B->X), A->SML());

            if (A->isFluid && B->isFluid) {
               A->U_XSPH += c_xsph * (-Uij) * w;  //\label{SPH:U_XSPH}
               total_w += w;
            }

            if (A->isFluid && B->isFluid) {
               if (A->isSurface || B->isSurface) {
                  // A->lap_U += _GRAVITY_ * B->volume * w_Bspline(Norm(A->X - B->X), A->particle_spacing) * Normalize(A->X - B->X);
                  A->lap_U += _GRAVITY_ * B->volume * w_Bspline(Norm(A->X - B->X), (1 - 1E-10) * A->particle_spacing) * Normalize(A->X - B->X) / (A->mu_SPH / A->rho);
                  // auto X_online = A->X + A->particle_spacing * Normalize(B->X - A->X);
                  // auto pro_Uij = Projection(Uij, A->X - X_online);
                  // A->lap_U += 0.1 * (-pro_Uij) * w_Bspline(Norm(B->X - X_online), A->particle_spacing);
                  // A->U_XSPH += c_xsph * (-Uij) * w;  //\label{SPH:U_XSPH}
                  // total_w += w;
               }
            }
            if ((A->isFluid && !B->isFluid) || (!A->isFluid && B->isFluid)) {
               A->lap_U += _GRAVITY_ * B->volume * w_Bspline(Norm(A->X - B->X), (1 - 1E-10) * A->particle_spacing) * Normalize(A->X - B->X) / (A->mu_SPH / A->rho);
            }

            A->div_U += B->volume * Dot(-Uij, grad_w_Bspline(A, B));  //\label{SPH:divU}
            Aij = 2. * B->volume * Dot_grad_w_Bspline(A, B);          //\label{SPH:lapP1}

            //! 修正
            // const auto DelX = (A->X - B->X);
            // applyOverPoints([&](const auto &Q) {
            //    A->lap_U -= Aij * Q->volume * Dot(DelX, grad_w_Bspline(A, Q)) * (Q->U_SPH - A->U_SPH);
            // });

            A->lap_U += Aij * Uij;  //\label{SPH:lapU}

            if (Between(Distance(A, B), {1E-13, A->SML()})) {
               A->checked_points_in_radius_SPH++;
               if (B->getNetwork()->isFluid || B->isFluid)
                  A->checked_points_in_radius_of_fluid_SPH++;
            }
            A->checked_points_SPH++;
         };

         applyOverPoints(add);

         if (total_w > too_small_total_w)
            A->U_XSPH /= total_w;

         /* ----------------------- calculate DUDt and b_vector ---------------------- */

         // \label{SPH:lapU_for_wall}
         // \label{SPH:Poisson_b_vector}
         // \label{SPH:how_to_set_fluid_b_vector}
         if (A->getNetwork()->isRigidBody) {
            A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + _GRAVITY3_;  //! 後で修正されるDUDt
            A->DUDt_SPH += A->U_XSPH / dt;
            //
            A->DrhoDt_SPH = -A->rho * A->div_U;
            A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
         } else {
            A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + _GRAVITY3_;  //! 後で修正されるDUDt
            A->DUDt_SPH += A->U_XSPH / dt;
            //
            A->DrhoDt_SPH = -A->rho * A->div_U;
            A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
         }
      }

   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in calcLaplacianU");
   };
};

#endif