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

void gradP(const std::unordered_set<networkPoint *> &points, const std::unordered_set<Network *> &target_nets) {
   try {
#pragma omp parallel
      for (const auto &A : points)
#pragma omp single nowait
      {
         A->gradP_SPH.fill(0.);
         const auto AX = USE_NEXT_POSITION ? X_next(A) : A->X;
         const double r = USE_NEXT_POSITION ? A->SML_next() : A->SML();
         // A->grad_corr_M.fill({0., 0., 0.});
         auto add_gradP_SPH_Where = [&](const auto &B) {
            const auto BX = USE_NEXT_POSITION ? X_next(B) : B->X;

            // if (A->isSurface || (!A->isSurface )) {
            //\label{SPH:gradP1}

            auto c = A->rho * B->mass;

            // if (B->isSurface && !B->isNeumannSurface) {
            //    const double a = 0.1;
            //    if (B->isSurface)
            //       c *= (1 - a);
            //    else if (B->isAuxiliary)
            //       c *= a;
            // }

            // if ((A->isAuxiliary && B->isSurface) || (A->isSurface && B->isAuxiliary))
            //    return;  // c=0

            // else if (!A->isSurface) {
            //    const double a = 0.1;
            //    if (B->isSurface)
            //       c *= (1 - a);
            //    else if (B->isAuxiliary)
            //       c *= a;
            // }

            auto grad = USE_NEXT_POSITION ? grad_w_Bspline_next(A, B) : grad_w_Bspline(A, B);
            A->gradP_SPH += B->p_SPH * (c / std::pow(B->rho, 2)) * grad;
            A->gradP_SPH += A->p_SPH * (c / std::pow(A->rho, 2)) * grad;
            //
            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->volume * grad;  //\label{SPH:gradP2}

            //%鏡写しにした分
            // if (A->isNotSurfaceButNearSurface && B->isSurface) {
            //    auto X = Mirror(AX, BX, BX - AX);
            //    auto c = A->rho * A->mass;
            //    auto grad = USE_NEXT_POSITION ? grad_w_Bspline_next(A, X) : grad_w_Bspline(A, X);
            //    A->gradP_SPH += A->p_SPH * (c / std::pow(A->rho, 2)) * grad;
            //    A->gradP_SPH += A->p_SPH * (c / std::pow(A->rho, 2)) * grad;
            //    // auto X = Mirror(pO_x, BX, pO_x - BX);
            //    // auto Aij = 2. * B->volume * Dot_grad_w_Bspline_Dot(pO_x, BX, r, USE_NEXT_POSITION ? pO->inv_grad_corr_M_next : pO->inv_grad_corr_M);
            //    // ROW->increment(pO, Aij);
            //    // ROW->increment(B, -Aij);
            //    // ROW->PoissonRHS += B->volume * Dot(B->b_vector - pO->b_vector, grad_w_Bspline(pO_x, X, r, USE_NEXT_POSITION ? pO->inv_grad_corr_M_next : pO->inv_grad_corr_M));
            //    // sum_Aij_Pj += Aij * B->p_SPH;
            //    // sum_Aij += Aij;
            // }
            //%

            // A->gradP_SPH += B->mass / A->rho * (B->p_SPH + A->p_SPH) * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP1}
            // A->grad_corr_M += B->volume * TensorProduct(X - A->X, grad_w_Bspline(A->X, X, A->SML()));

            // if (A->isSurface || A->isNeumannSurface)
            //    A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH_ / (A->rho * A->rho)) * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP1}
            // else

            // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH_ / (A->rho * A->rho)) * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP1}

            // if (A->isSurface)
            //    A->gradP_SPH += (B->p_SPH - A->p_SPH_) * B->volume * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP2}
            // else

            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->volume * grad_w_Bspline(X_next(A), X_next(B), A->SML());  //\label{SPH:gradP2}

            //
            // }
            // A->gradP_SPH = Dot(A->gradP_SPH, A->inv_grad_corr_M);

            // A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * Dot(grad_w_Bspline(A->X, X, A->SML()), A->inv_grad_corr_M);  //\label{SPH:gradP1}

            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->volume * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP2}

            //\label{SPH:gradP2}は，Aij = 1*..を使うと，0.132
            //\label{SPH:gradP2}は，Aij = 3*..を使うと，0.221
            //\label{SPH:gradP2}は，Aij = 2*..を使うと，0.21198

            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * V_next(B) * grad_w_Bspline(X_next(A), X_next(B), A->SML());  //\label{SPH:gradP2}
            // A->gradP_SPH += B->p_SPH * B->mass / B->rho * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP3}0.34
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

         auto add_gradP_SPH = [&](const auto &B) {
            return add_gradP_SPH_Where(B);
         };

         // networkPoint *closest_surface_point = nullptr;
         // double min_distance = 1e10, distance;

         // if (A->isSurface) {
         //    auto X = A->X - (A->COM_SPH - A->X);
         //    A->grad_corr_M += A->volume * TensorProduct(X - A->X, grad_w_Bspline(A->X, X, A->SML()));
         //    A->gradP_SPH += A->rho * A->mass * (A->p_SPH / (A->rho * A->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP1}
         // }

         for (const auto &net : target_nets) {
            net->BucketPoints.apply(A->X, 1.2 * r, [&](const auto &B) {
               if (B->isCaptured) {
                  add_gradP_SPH(B);
                  // if (B->isSurface) {
                  //    if ((distance = Distance(A->X, X_next(B))) < min_distance) {
                  //       min_distance = distance;
                  //       closest_surface_point = B;
                  //    }
                  // }
               }
            });
         }

         // if (closest_surface_point != nullptr) {
         //    auto B = closest_surface_point;
         //    auto X = A->X - (A->COM_SPH - A->X);
         //    A->grad_corr_M += B->volume * TensorProduct(X - A->X, grad_w_Bspline(A->X, X, A->SML()));
         //    A->gradP_SPH += A->rho * B->mass * (B->p_SPH / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP1}
         // }

         // A->inv_grad_corr_M = Inverse(A->grad_corr_M);

         /*DOC_EXTRACT 0_2_5_SPH

         $`\dfrac{D{\bf u}^n}{Dt} = - \frac{1}{\rho} \nabla p^{n+1} + \nu \nabla^2 {\bf u}^n + {\bf g}`$
         が計算できた．

         */

         // if (A->isSurface)
         //    A->gradP_SPH = (A->gradP_SPH + Dot(A->gradP_SPH, A->inv_grad_corr_M)) / 2.;
         // else
         // if (!(A->isSurface && !A->isNeumannSurface))
         // if (A->isNeumannSurface)
         //    A->DUDt_SPH = Chop(A->DUDt_SPH, A->interp_normal_water);

         if (!isFinite(A->DUDt_SPH))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
      }

      for (const auto &A : points) {
         // if (!A->isSurface && !A->isAuxiliary)
         A->DUDt_SPH -= A->gradP_SPH / A->rho;
         // else
         //    A->DUDt_SPH.fill(0.);
         // else if (A->isAuxiliary) {
         //    {
         //       //
         //       const double a = 0.1;
         //       auto n = A->surfacePoint->interp_normal_original_next;
         //       A->surfacePoint->DUDt_SPH -= (1. - a) * Projection(A->surfacePoint->gradP_SPH / A->surfacePoint->rho, n);
         //       A->surfacePoint->DUDt_SPH -= a * Projection(A->gradP_SPH / A->rho, n);  //! 中央が噴射する
         //       //
         //       // auto n = A->surfacePoint->interp_normal_original_next;
         //       // n += A->surfacePoint->interp_normal_original;
         //       // n /= 2.;
         //       // Tddd v = Chop(A->gradP_SPH / A->rho, n) + Projection(A->surfacePoint->gradP_SPH / A->surfacePoint->rho, n);
         //       // A->surfacePoint->DUDt_SPH -= optimumVector(std::vector<Tddd>{v, A->surfacePoint->gradP_SPH / A->surfacePoint->rho}, {0., 0., 0.}, 1E-12);
         //    }
         // }
      }
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in gradP");
   };
}

#endif