#ifndef SPH_lap_div_U_H
#define SPH_lap_div_U_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT SPH

## (1) $`\nabla^2 {\bf u}_i`$の計算（`calcLaplacianU`）

CHECKED: \ref{SPH:lapU}{流速のラプラシアンの計算方法}: $`\nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

CHECKED: \ref{SPH:divU}{流速の発散の計算方法}: $`\nabla\cdot{\bf u}_i=\sum_{j}\frac{m_j}{\rho_j}({{\bf u}_j-{\bf u}_i}) \cdot\nabla W_{ij}`$

*/

auto calcLaplacianU(const auto &points, const std::unordered_set<Network *> &target_nets, const double dt) {

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->checked_points_in_radius_SPH = A->checked_points_in_radius_of_fluid_SPH = A->checked_points_SPH = 0;
      A->div_U = 0.;
      A->lap_U.fill(0.);
      A->b_vector.fill(0.);

      A->grad_coeff.clear();
      A->grad_coeff_next.clear();

      A->density_based_on_positions = 0;

      auto add = [&](const auto &B) {
         const auto Uij = A->U_SPH - B->U_SPH;
         A->density_based_on_positions += B->rho * B->volume * w_Bspline(Norm(A->X - B->X), A->radius_SPH);
         A->div_U += B->volume * Dot(B->U_SPH - A->U_SPH, grad_w_Bspline(A->X, B->X, A->radius_SPH));  //\label{SPH:divU}
         A->lap_U += 2 * B->mass / A->rho * Uij * Dot_grad_w_Bspline_Dot(A->X, B->X, A->radius_SPH);   //\label{SPH:lapU}

         // just counting
         if (Between(Distance(A, B), {1E-12, A->radius_SPH})) {
            A->checked_points_in_radius_SPH++;
            if (B->getNetwork()->isFluid || B->isFluid)
               A->checked_points_in_radius_of_fluid_SPH++;
            // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP1}
            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->mass / A->rho * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP2}
            // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP3}
         }
         A->checked_points_SPH++;
      };

      // sum 計算
      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               add(B);
#if defined(USE_SHARED_AUX)
               if (B->isSurface && B == A)
                  for (const auto &AUX : B->auxiliaryPoints)
                     add(AUX);
#endif
            }
         });
#if defined(USE_SIMPLE_SINGLE_AUX)
      if (A->isSurface)
         for (const auto &AUX : A->auxiliaryPoints)
            PoissonEquation(AUX);
#endif
      //$ ------------------------------------------ */
      // \label{SPH:lapU_for_wall}
      // \label{SPH:Poisson_b_vector}
      if (A->getNetwork()->isRigidBody) {
         A->DUDt_SPH_.fill(0.);
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH.fill(0.);
         A->tmp_U_SPH.fill(0.);
         A->tmp_X = A->X;
         A->DrhoDt_SPH = 0.;
         // \label{SPH:how_to_set_wall_b_vector}
         // A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
         A->b_vector = A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
         //
         // A->b_vector.fill(0.);
      } else {
         A->DUDt_SPH_ = A->DUDt_SPH;
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH = nu * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt
         A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
         A->tmp_X = A->X + A->tmp_U_SPH * dt;
         A->DrhoDt_SPH = -A->rho * A->div_U;

         // \label{SPH:how_to_set_fluid_b_vector}
         A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;

         if (A->vec_time_SPH.size() > 10) {
#if defined(USE_RungeKutta)
            double current_time = A->RK_X.get_t();
            double next_time = current_time + A->RK_X.get_dt();
#elif defined(USE_LeapFrog)
            double current_time = A->LPFG_X.get_t();
            double next_time = current_time + dt;
#endif
            std::vector<double> times = {next_time, current_time};
            std::array<double, 3> U1, U2, U3, U4;
            U1 = A->U_SPH;
            if (*(A->vec_time_SPH.rbegin()) == current_time) {
               times.push_back(*(A->vec_time_SPH.rbegin() + 1));
               U2 = *(A->vec_U_SPH.rbegin() + 1);
               times.push_back(*(A->vec_time_SPH.rbegin() + 2));
               U3 = *(A->vec_U_SPH.rbegin() + 2);
               // times.push_back(*(A->vec_time_SPH.rbegin() + 3));
               // U4 = *(A->vec_U_SPH.rbegin() + 3);
            } else {
               times.push_back(*(A->vec_time_SPH.rbegin() + 0));
               U2 = *(A->vec_U_SPH.rbegin() + 0);
               times.push_back(*(A->vec_time_SPH.rbegin() + 1));
               U3 = *(A->vec_U_SPH.rbegin() + 1);
               // times.push_back(*(A->vec_time_SPH.rbegin() + 2));
               // U4 = *(A->vec_U_SPH.rbegin() + 2);
            }
            InterpolationLagrange<double> lag(times);
            auto D = lag.DN(current_time);
            A->b_vector = -(D[1] * U1 + D[2] * U2 + D[3] * U3 /* + D[4] * U4*/) + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
         }
      }

      //$ ------------------------------------------ */

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
   }
};

#endif