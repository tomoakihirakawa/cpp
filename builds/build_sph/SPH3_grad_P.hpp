#ifndef SPH_grad_P_H
#define SPH_grad_P_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */

/*DOC_EXTRACT SPH

## 圧力勾配$`\nabla p^{n+1}`$の計算

CHECKED: \ref{SPH:gradP1}{勾配の計算方法}: $`\nabla p_i = \rho_i \sum_{j} m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}) \nabla W_{ij}`$

CHECKED: \ref{SPH:gradP2}{勾配の計算方法}: $`\nabla p_i = \rho_i \sum_{j} m_j \left(p_j - p_i\right) \nabla W_{ij}`$

CHECKED: \ref{SPH:gradP3}{勾配の計算方法}: $`\nabla p_i = \sum_{j} \frac{m_j}{\rho_j} p_j \nabla W_{ij}`$

*/

void gradP(const std::unordered_set<networkPoint *> &points, const std::unordered_set<Network *> &target_nets) {

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->gradP_SPH.fill(0.);
      auto add_gradP_SPH = [&](const auto &B, const double &coef = 1.) {
         A->gradP_SPH += coef * A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP1}

         // A->gradP_SPH = Dot(A->gradP_SPH, A->inv_grad_corr_M);

         // A->gradP_SPH += coef * A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * Dot(grad_w_Bspline(A->X, B->X, A->radius_SPH), A->inv_grad_corr_M);  //\label{SPH:gradP1}

         // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->volume * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP2}

         //\label{SPH:gradP2}は，Aij = 1*..を使うと，0.132
         //\label{SPH:gradP2}は，Aij = 3*..を使うと，0.221
         //\label{SPH:gradP2}は，Aij = 2*..を使うと，0.21198

         // A->gradP_SPH += (B->p_SPH - A->p_SPH) * V_next(B) * grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH);  //\label{SPH:gradP2}
         // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(A->X, B->X, A->radius_SPH);  //\label{SPH:gradP3}0.34
         //\label{SPH:gradP3}は，Aij = 2*..を使うと，0.34
         //\label{SPH:gradP3}は，Aij = 2*..を使うと，Lagrangeを使うと，0.46
         //\label{SPH:gradP3}は，Aij = 3*..を使うと，0.49
         //\label{SPH:gradP3}は，Aij = 3*..を使うと，Lagrangeを使うと，0.4は超えたが，綺麗な結果ではなく，つぶれた．
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，_nextを使うと，0.24
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，X_next以外の_nextを使うと，0.277>
         //\label{SPH:gradP3}は，Aij = 2*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.34
         //\label{SPH:gradP3}は，Aij = 2.5*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.45
         //\label{SPH:gradP3}は，Aij = 3*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.49
         //\label{SPH:gradP3}は，Aij = 4*.. さらに，X_next以外の_nextを使う．しかし，実際は密度は一定とすると，0.457
      };

      networkPoint *closest_surface_point = nullptr;
      double min_distance = 1e10, distance;

      for (const auto &net : target_nets) {

         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               add_gradP_SPH(B);

#if defined(USE_ONE_AUXP)
               if (B->isSurface) {
                  if ((distance = Distance(A->X, X_next(B))) < min_distance) {
                     min_distance = distance;
                     closest_surface_point = B;
                  }
               }
#elif defined(USE_ALL_AUXP)
            if (B->isSurface)
               for (const auto &AUX : B->auxiliaryPoints) 
               if (AUX != nullptr)
                  add_gradP_SPH(AUX);
#endif
            }
         });
      }
//
#if defined(USE_ONE_AUXP)
      // if (A->isSurface)
      if (closest_surface_point != nullptr)
         for (const auto &AUX : closest_surface_point->auxiliaryPoints)
            if (AUX != nullptr)
               add_gradP_SPH(AUX);
#endif

      /*DOC_EXTRACT SPH

      $`\dfrac{D{\bf u}^n}{Dt} = - \frac{1}{\rho} \nabla p^{n+1} + \nu \nabla^2 {\bf u}^n + {\bf g}`$
      が計算できた．

      */
      // if (!A->isSurface)
      //    A->gradP_SPH = Dot(A->gradP_SPH, A->inv_grad_corr_M);
      //
      A->DUDt_SPH -= A->gradP_SPH / A->rho;

      if (!isFinite(A->DUDt_SPH))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
   }
}

#endif