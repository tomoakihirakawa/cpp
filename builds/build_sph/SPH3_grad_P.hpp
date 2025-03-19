#ifndef SPH_grad_P_H
#define SPH_grad_P_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */

/*DOC_EXTRACT 0_2_5_SPH

## 圧力勾配$`\nabla p^{n+1}`$の計算

CHECKED: \ref{SPH:gradP1}{勾配の計算方法}: $`\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}`$

CHECKED: \ref{SPH:gradP2}{勾配の計算方法}: $`\nabla p_i = \sum_{j} \frac{m_j}{\rho_i} \left(p_j - p_i\right) \nabla W_{ij}`$

CHECKED: \ref{SPH:gradP3}{勾配の計算方法}: $`\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}`$

NOTE:圧力の方程式を立てる際に，左辺の密度として，流速の発散から見積もった$`\rho^({\rm next})`$を使うことは，
言い換えれば，N.S.方程式の圧力項の計算には，$`\rho^({\rm next})`$を使うと決めたことになる．
なので，圧力勾配$`\nabla p^{n+1}`$の計算にも，$`\rho^({\rm next})`$を使わなければならない．

*/

void gradP(const std::unordered_set<networkPoint*>& points, const std::unordered_set<Network*>& target_nets) {
   DebugPrint("setWall", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__, " ", Cyan);

#pragma omp parallel
   for (const auto& A : points)
#pragma omp single nowait
   {
      A->gradP_W_SPH.fill(0.);//\label{SPH:gradP1}
      A->p_SPH_ = 0.;
      A->Integral_gradP_W_SPH.fill(0.);
      auto add_gradP_W_SPH = [&](const auto& B) {
         auto V = (B->p_SPH) * w_Bspline(0., A->SML_next());
         A->gradP_W_SPH += V_next(B) * V * grad_w_Bspline_next(A, B);
         A->p_SPH_ += V_next(B) * B->p_SPH * w_Bspline(Norm(X_next(A) - X_next(B)), A->SML_next());
         };

      const double r = A->SML_grad_next();
      for (const auto& net : target_nets) {
         net->BucketPoints.apply(A->X, 1.1 * r, [&](const auto& B) {
            if (canInteract(A, B))
               add_gradP_W_SPH(B);
            });
      }
   }

#pragma omp parallel
   for (const auto& A : points)
#pragma omp single nowait
   {
      A->gradP_SPH.fill(0.);//\label{SPH:gradP1}
      auto add_gradP_SPH = [&](const auto& B) {
         // A->gradP_SPH += V_next(B) * ((B->p_SPH + A->p_SPH) * grad_w_Bspline_next(A, B));

         // A->gradP_SPH += ((V_next(B) * B->p_SPH + V_next(A) * A->p_SPH) * grad_w_Bspline_next(A, B));

         A->gradP_SPH += rho_next(B) * B->mass * (B->p_SPH / std::pow(rho_next(B), 2) + A->p_SPH / std::pow(rho_next(A), 2)) * grad_w_Bspline_next(A, B);
         A->Integral_gradP_W_SPH += V_next(B) * B->gradP_W_SPH * w_Bspline(Norm(X_next(A) - X_next(B)), A->SML_next());
         };

      for (const auto& net : target_nets) {
         net->BucketPoints.apply(A->X, 1.1 * A->SML_grad_next(), [&](const auto& B) {
            if (canInteract(A, B))
               add_gradP_SPH(B);
            });
      }

      if (!isFinite(A->DUDt_SPH)) {
         std::cout << "A->DUDt_SPH = " << A->DUDt_SPH << std::endl;
         std::cout << "A->mu_SPH = " << A->mu_SPH << std::endl;
         std::cout << "A->rho = " << A->rho << std::endl;
         std::cout << "A->lap_U = " << A->lap_U << std::endl;
         std::cout << "A->U_XSPH = " << A->U_XSPH << std::endl;
         std::cout << "A->U_SPH = " << A->U_SPH << std::endl;
         std::cout << "A->b_vector = " << A->b_vector << std::endl;
         std::cout << "A->DrhoDt_SPH = " << A->DrhoDt_SPH << std::endl;
         std::cout << "A->div_U = " << A->div_U << std::endl;
         std::cout << "A->gradP_SPH = " << A->gradP_SPH << std::endl;
         std::cout << "A->laplacian_corr_M = " << A->laplacian_corr_M << std::endl;
         std::cout << "A->grad_corr_M = " << A->grad_corr_M << std::endl;
         std::cout << "A->grad_corr_M_next = " << A->grad_corr_M_next << std::endl;
         std::cout << "A->inv_grad_corr_M = " << A->inv_grad_corr_M << std::endl;
         std::cout << "A->inv_grad_corr_M_next = " << A->inv_grad_corr_M_next << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
      }
   }

   for (const auto& A : points) {
      if (A->isAuxiliary)
         A->DUDt_SPH -= Chop(A->gradP_SPH / A->rho, A->surfacePoint->interp_normal_original_next);
      else
         A->DUDt_SPH -= A->gradP_SPH / A->rho;
   }

   DebugPrint("setWall done", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__, " ", Cyan);
}

#endif