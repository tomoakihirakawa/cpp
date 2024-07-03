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
   DebugPrint("setWall", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__, " ", Cyan);
   try {
#pragma omp parallel
      for (const auto &A : points)
#pragma omp single nowait
      {
         A->gradP_SPH.fill(0.);
         // A->grad_corr_M.fill({0., 0., 0.});
         auto add_gradP_SPH = [&](const auto &B) {
            //\label{SPH:gradP1}

            auto c = rho_next(B) * B->mass;

            // const double a = 0.5;
            // if (A->isSurface) {
            //    if (B->isAuxiliary) {
            //       c *= a;
            //    }
            // } else if (!A->isSurface) {
            //    if (B->isSurface && !B->isAuxiliary) {
            //       c *= (1 - a);
            //    } else if (B->isSurface && B->isAuxiliary) {
            //       c *= a;
            //    }
            // }

            // if ((A->isAuxiliary && B->isSurface) || (A->isSurface && B->isAuxiliary))
            //    return;  // c=0

            // if (B->isSurface) {
            //    c *= 0.5;
            // } else if (B->isSurface) {
            //    c *= 0.5;
            // }

            // auto grad = grad_w_Bspline_next(A, B);

            auto grad = grad_w_Bspline_next(A, B);

            // #ifdef USE_SYMMETRIC_FORM_FOR_PRESSURE_GRADIENT
            //! ダムブレイクの計算で壁面付近の粒子の挙動はこっちの方が良かった．
            //! この方法は運動量を保存するらしい．
            if (A->var_Eigenvalues_of_M_next < 0.03)
               FusedMultiplyIncrement(V_next(B), (B->p_SPH - A->p_SPH) * grad, A->gradP_SPH);
            else
               FusedMultiplyIncrement(B->p_SPH * (c / std::pow(rho_next(B), 2)) + A->p_SPH * (c / std::pow(rho_next(A), 2)), grad, A->gradP_SPH);
            //
            // #elif defined(USE_SUBTRACTIVE_FORM_FOR_PRESSURE_GRADIENT)
            //             //! 静水の計算でこっちの方が良かった．
            //             FusedMultiplyIncrement(V_next(B), (B->p_SPH - A->p_SPH) * grad, A->gradP_SPH);
            // #endif
            // #ifdef USE_SYMMETRIC_FORM_FOR_PRESSURE
            //             //! これまでの方法
            //             // A->gradP_SPH += B->p_SPH * (c / std::pow(rho_next(B), 2)) * grad;
            //             // A->gradP_SPH += A->p_SPH * (c / std::pow(rho_next(A), 2)) * grad;
            //             FusedMultiplyIncrement(B->p_SPH * (c / std::pow(rho_next(B), 2)), grad, A->gradP_SPH);
            //             FusedMultiplyIncrement(A->p_SPH * (c / std::pow(rho_next(A), 2)), grad, A->gradP_SPH);
            // #elif defined(USE_RANDLES_LIBERSKY_FORM_FOR_PRESSURE)
            //             FusedMultiplyIncrement(V_next(B), (B->p_SPH - A->p_SPH) * grad, A->gradP_SPH);
            // #endif

            // A->gradP_SPH += V_next(B) * B->p_SPH * grad_w_Bspline_next(A, B);

            // A->gradP_SPH += (B->p_SPH - A->p_SPH) * B->volume * grad;  //\label{SPH:gradP2}

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

         // auto add_gradP_SPH = [&](const auto &B) {
         //    return add_gradP_SPH_Where(B);
         // };

         // networkPoint *closest_surface_point = nullptr;
         // double min_distance = 1e10, distance;

         // if (A->isSurface) {
         //    auto X = A->X - (A->COM_SPH - A->X);
         //    A->grad_corr_M += A->volume * TensorProduct(X - A->X, grad_w_Bspline(A->X, X, A->SML()));
         //    A->gradP_SPH += A->rho * A->mass * (A->p_SPH / (A->rho * A->rho) + A->p_SPH / (A->rho * A->rho)) * grad_w_Bspline(A->X, X, A->SML());  //\label{SPH:gradP1}
         // }

         const double r = A->SML_next();
         for (const auto &net : target_nets) {
            net->BucketPoints.apply(A->X, 1.1 * r, [&](const auto &B) {
               if (canInteract(A, B) && A != B) {
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

      // for (const auto &A : points) {
      //    if (A->isAuxiliary)
      //       A->surfacePoint->gradP_SPH = A->gradP_SPH;
      // }

      for (const auto &A : points) {
         if (A->isAuxiliary)
            A->DUDt_SPH -= Chop(A->gradP_SPH / A->rho, A->surfacePoint->interp_normal_original_next);
         else
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

      // for (const auto &A : points)
      //    if (A->isAuxiliary) {
      //       A->surfacePoint->DUDt_SPH = A->DUDt_SPH;
      //    }

   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in gradP");
   };
   DebugPrint("setWall done", __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__, " ", Cyan);
}

#endif