#ifndef SPH_lap_div_U_H
#define SPH_lap_div_U_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT 0_2_0_SPH

## 粘性項$`\nabla^2 {\bf u}_i`$の計算（`calcLaplacianU`）

SELECTED: \ref{SPH:lapU}{流速のラプラシアンの計算方法}: $`\nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

SELECTED: \ref{SPH:divU}{流速の発散の計算方法}: $`\nabla\cdot{\bf u}_i=\sum_{j}\frac{m_j}{\rho_j}({{\bf u}_j-{\bf u}_i}) \cdot\nabla W_{ij}`$

*/

auto calcLaplacianU(const auto &points, const std::unordered_set<Network *> &target_nets) {

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

      // A->grad_U.fill({0., 0., 0.});
      double Aij;
      //
      auto add = [&](const auto &B) {
         const auto Uij = A->U_SPH - B->U_SPH;
         // A->div_U += B->volume * Dot(-Uij, grad_w_Bspline(A->X, B->X, A->radius_SPH));  //\label{SPH:divU}
         A->div_U += B->volume * Dot(-Uij, Dot(grad_w_Bspline(A->X, B->X, A->radius_SPH), A->inv_grad_corr_M));  //\label{SPH:divU}
         // A->grad_U += B->volume * TensorProduct(-Uij, Dot(grad_w_Bspline(A->X, B->X, A->radius_SPH), A->inv_grad_corr_M));  //\label{SPH:gradU}
         Aij = 2. * B->volume * Dot_grad_w_Bspline_Dot(A->X, B->X, A->radius_SPH);  //\label{SPH:lapP1}
         // Aij = 2. * B->volume * Dot_grad_w_Bspline_Dot_Modified(A->X, B->X, A->radius_SPH, A->inv_grad_corr_M);  //\label{SPH:lapP1}
         A->lap_U += Aij * Uij;  //\label{SPH:lapU}
         // A->lap_U -= Dot(A->X - B->X, A->grad_U);
         for (const auto &[p, v] : A->vector_p_grad)
            A->lap_U -= Aij * v * p->U_SPH;
         // A->lap_U += -8. * B->mass / (A->rho + B->rho) * Dot(Uij, A->X - B->X) / Dot(A->X - B->X, A->X - B->X) * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:lapP2}
         // just counting
         if (Between(Distance(A, B), {1E-13, A->radius_SPH})) {
            A->checked_points_in_radius_SPH++;
            if (B->getNetwork()->isFluid || B->isFluid)
               A->checked_points_in_radius_of_fluid_SPH++;
            // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP1}
            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->mass / A->rho * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP2}
            // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP3}
         }
         A->checked_points_SPH++;
      };

      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
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
         if (A->isNeumannSurface) {
            A->DrhoDt_SPH = -A->rho * A->div_U;
            A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;最も自然な結果を返す
         }
      } else {
         A->DUDt_SPH_ = A->DUDt_SPH;
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH = nu * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt

         A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
         A->tmp_X = A->X + A->tmp_U_SPH * dt;
         A->DrhoDt_SPH = -A->rho * A->div_U;
         // \label{SPH:how_to_set_fluid_b_vector}
         A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;

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
};

#endif