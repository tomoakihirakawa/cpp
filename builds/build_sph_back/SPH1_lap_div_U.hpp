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
         A->lap_U_next.fill(0.);
         A->convection_term.fill(0.);
         A->b_vector.fill(0.);

         // A->grad_U.fill({0., 0., 0.});
         double Aij, Aij_next;
         //
         auto applyOverPoints = [&](const auto &equation, const auto &pO_x, const std::unordered_set<Network *> NETS) {
            const auto r = A->SML();
            for (const auto &net : NETS) {
               net->BucketPoints.apply(pO_x, 1.1 * r, [&](const auto &B) {
                  if (B->isCaptured)
                     equation(B);
               });
            }
         };
         //
         auto add = [&](const auto &B) {
            const auto Uij = A->U_SPH - B->U_SPH;
            A->div_U += B->volume * Dot(-Uij, grad_w_Bspline(A, B));    //\label{SPH:divU}
            Aij = 2. * B->volume * Dot_grad_w_Bspline(A, B);            //\label{SPH:lapP1}
            Aij_next = 2. * V_next(B) * Dot_grad_w_Bspline_next(A, B);  //\label{SPH:lapP1}

#if defined(USE_Laplacian_correction)
            //! 修正
            // Aij *= -1;
            const auto rij = (A->X - B->X);
            applyOverPoints([&](const auto &Q) {
               A->lap_U -= Aij * Dot(rij, Q->volume * (Q->U_SPH - A->U_SPH) * grad_w_Bspline(A, Q));
               A->lap_U_next -= Aij_next * Dot(rij, V_next(Q) * (Q->U_SPH - A->U_SPH) * grad_w_Bspline_next(A, Q));
            },
                            A->X, target_nets);
//!
#endif
            A->lap_U += Aij * Uij;            //\label{SPH:lapU}
            A->lap_U_next += Aij_next * Uij;  // Uij_nextはわからないので，Uijで代用
            // just counting
            if (Between(Distance(A, B), {1E-13, A->SML()})) {
               A->checked_points_in_radius_SPH++;
               if (B->getNetwork()->isFluid || B->isFluid)
                  A->checked_points_in_radius_of_fluid_SPH++;
            }
            A->checked_points_SPH++;
         };

         for (const auto &net : target_nets)
            net->BucketPoints.apply(A->X, 1.1 * A->SML(), [&](const auto &B) {
               if (B->isCaptured) add(B);
            });

         // A->lap_U = Dot(A->Mat_B, A->lap_U);
         //$ ------------------------------------------ */
         // \label{SPH:lapU_for_wall}
         // \label{SPH:Poisson_b_vector}
         if (A->getNetwork()->isRigidBody) {
            A->DUDt_SPH_.fill(0.);
            A->rho = _WATER_DENSITY_;
            double nu = A->mu_SPH / A->rho;
            A->DUDt_SPH.fill(0.);
            A->tmp_U_SPH.fill(0.);
            A->tmp_X = A->X;
            A->DrhoDt_SPH = 0;
            A->b_vector.fill(0.);
            // if (A->isNeumannSurface)

            A->DrhoDt_SPH = -A->rho * A->div_U;
            A->b_vector = A->rho * A->U_SPH / dt + A->mu_SPH * A->lap_U + A->rho * _GRAVITY3_;  // 最も自然な結果を返す
            // auto rho = A->intp_density;
            // auto rho = A->rho;
            // A->b_vector = rho * A->U_SPH / dt + A->mu_SPH * A->lap_U + rho * _GRAVITY3_;  // 最も自然な結果を返す

         } else {
            A->DUDt_SPH_ = A->DUDt_SPH;
            A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt

            A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
            A->tmp_X = A->X + A->tmp_U_SPH * dt;
            A->DrhoDt_SPH = -A->rho * A->div_U;
            // \label{SPH:how_to_set_fluid_b_vector}

            A->b_vector = A->rho * A->U_SPH / dt + A->mu_SPH * A->lap_U + A->rho * _GRAVITY3_;
            // auto rho = A->intp_density;
            // auto rho = A->rho;
            // A->b_vector = rho * A->U_SPH / dt + A->mu_SPH * A->lap_U + rho * _GRAVITY3_;

            //          if (A->vec_time_SPH.size() > 10) {
            // #if defined(USE_RungeKutta)
            //             double current_time = A->RK_X.get_t();
            //             double next_time = current_time + A->RK_X.get_dt();
            // #elif defined(USE_LeapFrog)
            //             double current_time = A->LPFG_X.get_t();
            //             double next_time = current_time + dt;
            // #endif

            //             std::vector<double> times = {next_time, current_time};
            //             std::array<double, 3> U1, U2, U3;
            //             U1 = A->U_SPH;
            //             if (*(A->vec_time_SPH.rbegin()) == current_time) {
            //                times.push_back(*(A->vec_time_SPH.rbegin() + 1));
            //                U2 = *(A->vec_U_SPH.rbegin() + 1);
            //                // times.push_back(*(A->vec_time_SPH.rbegin() + 2));
            //                // U3 = *(A->vec_U_SPH.rbegin() + 2);
            //             } else {
            //                times.push_back(*(A->vec_time_SPH.rbegin() + 0));
            //                U2 = *(A->vec_U_SPH.rbegin() + 0);
            //                // times.push_back(*(A->vec_time_SPH.rbegin() + 1));
            //                // U3 = *(A->vec_U_SPH.rbegin() + 1);
            //             }
            //             // size of times is 4
            //             InterpolationLagrange<double> lag(times);
            //             auto D = lag.DN(current_time);
            //             A->b_vector = -(D[1] * U1 + D[2] * U2) + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
            //          }
         }
      }

   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in calcLaplacianU");
   };
};

#endif