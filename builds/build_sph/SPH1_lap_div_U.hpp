#ifndef SPH_lap_div_U_H
#define SPH_lap_div_U_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT 0_2_0_lap_div_U

## 粘性項$`\nabla^2 {\bf u}_i`$の計算（`calcLaplacianU`）

SELECTED: \ref{SPH:lapU}{流速のラプラシアンの計算方法}: $`\nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

SELECTED: \ref{SPH:divU}{流速の発散の計算方法}: $`\nabla\cdot{\bf u}_i=\sum_{j}\frac{m_j}{\rho_j}({{\bf u}_j-{\bf u}_i}) \cdot\nabla W_{ij}`$

*/

#define USE_PRE_CALC_tensorproduct_grad_Uij

auto calcLaplacianU(const auto &points, const std::unordered_set<Network *> &target_nets) {

   std::vector<Network *> RigidBodyObject, FluidObject;
   for (const auto &net : target_nets)
      if (net->isRigidBody)
         RigidBodyObject.push_back(net);
   for (const auto &net : target_nets)
      if (net->isFluid)
         FluidObject.push_back(net);

   try {
      // b! ----------------------------- 基本的な粘性項と重力項の計算 ----------------------------- */
      {
#pragma omp parallel
         for (const auto &A : points)
#pragma omp single nowait
         {
            const double dt = A->RK_X.getTimeAtNextStep() - A->RK_X.getTimeAtCurrentStep();

            A->div_U = 0.;
            A->lap_U.fill(0.);
            A->b_vector.fill(0.);
            A->U_XSPH.fill(0.);
            A->U_next.fill(0.);
            A->DUDt_modify_SPH.fill(0.);
            A->DUDt_SPH.fill(0.);
            A->DrhoDt_SPH = 0.;
            A->DUDt_modify_SPH_2.fill(0.);

            //! ラプラシアンの修正行列の一部を計算する
            Fill(A->tensorproduct_grad_Uij, 0.);
            Fill(A->tensorproduct_grad_Uij_next, 0.);

            const auto r = A->SML();
            const auto center = A->X;

            if (A->isFluid) {
               for (const auto &Q : A->points_in_SML) {
                  A->tensorproduct_grad_Uij += TensorProduct(Q->volume * grad_w_Bspline(A, Q), Q->U_SPH - A->U_SPH);
                  A->tensorproduct_grad_Uij_next += TensorProduct(V_next(Q) * grad_w_Bspline_next(A, Q), U_next(Q) - U_next(A));
               }
            } else
               for (const auto &net : target_nets) {
                  net->BucketPoints.applyVec1D(center, 1.1 * r, [&](const auto &vecQ) {
                     for (const auto &Q : vecQ)
                        if (canInteract(A, Q) && A != Q) {
                           A->tensorproduct_grad_Uij += TensorProduct(Q->volume * grad_w_Bspline(A, Q), Q->U_SPH - A->U_SPH);
                           A->tensorproduct_grad_Uij_next += TensorProduct(V_next(Q) * grad_w_Bspline_next(A, Q), U_next(Q) - U_next(A));
                        }
                  });
               }

            double c_xsph = 0.04;

            // if (A->isFluid) {
            //    // c_xsph = std::clamp(0., 0.03, std::pow(A->var_Eigenvalues_of_M, 2));
            //    c_xsph = std::clamp(0.03, 0.06, 0.1 * A->var_Eigenvalues_of_M);
            // }
            const double csml_factor = 1.;
            double total_w = 0, w, vol_w, total_w_U_next = 0;
            Tddd Uij = {0., 0., 0.}, DelX, vec_A2B, A_no_oikoshi_velocity;
            auto add = [&](const std::vector<networkPoint *> &vecB) {
               for (const auto &B : vecB)
                  if (canInteract(A, B) && A != B) {

                     vec_A2B = B->X - A->X;
                     A_no_oikoshi_velocity = A->RK_U.getX() - (B->isFluid ? B->RK_U.getX() : std::array<double, 3>{0., 0., 0.});

                     Uij = A->U_SPH - B->U_SPH;
                     if (A->isFluid && B->isFluid && !B->isAuxiliary && A != B && !B->is_grad_corr_M_singular) {
                        // if (((A->isFluid && B->isFluid) || (A->isFluid && B->isFirstWallLayer)) && !B->isAuxiliary && A != B) {
                        vol_w = B->volume * w_Bspline(Norm(A->X - B->X), A->SML());
                        // if (B->isSurface)
                        //    A->U_XSPH += 2 * c_xsph * (-Uij) * vol_w;  //\label{SPH:U_XSPH}
                        // else
                        A->U_XSPH -= c_xsph * Uij * vol_w;  //\label{SPH:U_XSPH}
                        total_w += vol_w;
                     }

                     // A->div_U += B->volume * Dot(-Uij, grad_w_Bspline(A, B));       //\label{SPH:divU}
                     FusedMultiplyIncrement(-B->volume, Dot(Uij, grad_w_Bspline(A, B)), A->div_U);
                     const double Aij = 2. * B->volume * Dot_grad_w_Bspline(A, B);  //\label{SPH:lapP1}
                     // const double Aij = 2. * B->volume * Dot_grad_w_Bspline(A->X, B->X, A->SML());

                     //! 修正
                     DelX = (A->X - B->X);
                     if (!A->is_grad_corr_M_singular)
                        FusedMultiplyIncrement(-Aij, Dot(DelX, A->tensorproduct_grad_Uij), A->lap_U);

                     //\label{SPH:lapU}
                     FusedMultiplyIncrement(Aij, Uij, A->lap_U);

                     /* -------------------------------------------------------------------------- */
                     // if (Dot(A_no_oikoshi_velocity, vec_A2B) > 0 /*まずは，被るかどうかのチェック*/) {
                     //    double h = A->particle_spacing;
                     //    auto tmp = (0.00001 * Projection(A_no_oikoshi_velocity, vec_A2B) / dt + 0.000001 * _GRAVITY_ * Normalize(vec_A2B)) * w_Bspline(Norm(vec_A2B), h) / w_Bspline(0., h);
                     //    if (B->isSurface || !B->isFluid)
                     //       tmp *= 5.;
                     //    A->DUDt_modify_SPH_2 -= tmp;
                     // }
                  }
            };

            // applyOverPoints(add);
            if (A->isFluid)
               add(A->points_in_SML);
            else
               for (const auto &net : target_nets) {
                  net->BucketPoints.applyVec1D(center, 1.1 * r, [&](const auto &vecQ) {
                     add(vecQ);
                  });
               }

            // if (A->isSurface && A->isFluid) {
            //    // A->DUDt_modify_SPH = 0.0001 * Dot(A->inv_grad_corr_M, A->interp_normal_original) / dt;
            //    A->DUDt_modify_SPH_2 = Chop(A->DUDt_modify_SPH_2, Dot(A->inv_grad_corr_M, A->interp_normal_original));
            // }

            if (total_w_U_next > too_small_total_w)
               A->U_next /= total_w_U_next;
         }

         for (const auto &A : points) {
            /* ----------------------- calculate DUDt and b_vector ---------------------- */
            // \label{SPH:lapU_for_wall}
            // \label{SPH:Poisson_b_vector}
            // \label{SPH:how_to_set_fluid_b_vector}
            const double dt = A->RK_X.getTimeAtNextStep() - A->RK_X.getTimeAtCurrentStep();
            if (A->getNetwork()->isRigidBody) {
               A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt;
               A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
               A->DUDt_SPH += _GRAVITY3_;
               A->DrhoDt_SPH = -A->rho * A->div_U;
            } else {
               A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt + A->DUDt_modify_SPH_2;
               A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
               A->DUDt_SPH += _GRAVITY3_;
               A->DrhoDt_SPH = -A->rho * A->div_U;
            }
            if (!isFinite(A->DUDt_SPH)) {
               // show details
               std::cout << "A->DUDt_SPH = " << A->DUDt_SPH << std::endl;
               std::cout << "A->mu_SPH = " << A->mu_SPH << std::endl;
               std::cout << "A->rho = " << A->rho << std::endl;
               std::cout << "A->tensorproduct_grad_Uij = " << A->tensorproduct_grad_Uij << std::endl;
               std::cout << "A->lap_U = " << A->lap_U << std::endl;
               std::cout << "A->U_XSPH = " << A->U_XSPH << std::endl;
               std::cout << "dt = " << dt << std::endl;
               std::cout << "A->U_SPH = " << A->U_SPH << std::endl;
               std::cout << "A->b_vector = " << A->b_vector << std::endl;
               std::cout << "A->DrhoDt_SPH = " << A->DrhoDt_SPH << std::endl;
               std::cout << "A->div_U = " << A->div_U << std::endl;
               std::cout << "A->laplacian_corr_M = " << A->laplacian_corr_M << std::endl;
               std::cout << "A->laplacian_corr_M_next = " << A->laplacian_corr_M_next << std::endl;
               std::cout << "A->grad_corr_M = " << A->grad_corr_M << std::endl;
               std::cout << "A->grad_corr_M_next = " << A->grad_corr_M_next << std::endl;
               std::cout << "A->inv_grad_corr_M = " << A->inv_grad_corr_M << std::endl;
               std::cout << "A->inv_grad_corr_M_next = " << A->inv_grad_corr_M_next << std::endl;
               // std::cout << "total_w = " << total_w << std::endl;
               // std::cout << "total_w_U_next = " << total_w_U_next << std::endl;
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
            }
         }
      }

      //% div_U_nextの計算
      /* -------------------------------------------------------------------------- */

      // b% ------------------------------------ さらなる修正 ------------------------------------ */

      // #pragma omp parallel
      //       for (const auto &A : points)
      // #pragma omp single nowait
      //       {
      //          A->DUDt_modify_SPH.fill(0.);
      //          A->DUDt_modify_SPH_2.fill(0.);
      //          if (A->isFluid) {
      //             const auto r = A->SML();
      //             double c = 0.001;
      //             A->DUDt_modify_SPH_2.fill(0.);
      //             const double dt = A->RK_X.get_dt();
      //             auto add = [&](const auto &B) {
      //                // if (A->isFluid && B->isFluid && A != B) {
      //                if (A != B) {
      //                   // const auto vec_A2B = X_next(B) - X_next(A);
      //                   // const auto A_no_oikoshi_velocity = U_next(A) - U_next(B);
      //                   const auto vec_A2B = B->X - A->X;
      //                   const auto A_no_oikoshi_velocity = A->RK_U.getX() - B->RK_U.getX();
      //                   //! Aの速度はBの速度に対してどうなっているか．Aは追い越そうとしているか？
      //                   //! Dot(vec_A2B, relative_velocity) > 0追い越し
      //                   if (Dot(A_no_oikoshi_velocity, vec_A2B) > 0 /*まずは，被るかどうかのチェック*/) {
      //                      double h = A->particle_spacing;
      //                      // auto chopped_vector = -Chop(vec_A2B, A_no_oikoshi_velocity);  // 重なりがなくなる，離れる方向なのでマイナス
      //                      // auto repulsive_dir = Normalize(chopped_vector);
      //                      // auto overlapped_distance = Norm(chopped_vector);
      //                      // local_DUDt_modify_SPH_2 += Norm(A_no_oikoshi_velocity) * repulsive_dir * w_Bspline(Norm(overlapped_distance), h) / max_w;
      //                      // local_DUDt_modify_SPH_2 -= c * Projection(A_no_oikoshi_velocity, vec_A2B) * w_Bspline(Norm(vec_A2B), h) / max_w;
      //                      auto tmp = (0.0001 * Projection(A_no_oikoshi_velocity, vec_A2B) / dt + 0.0001 * _GRAVITY_ * Normalize(vec_A2B)) * w_Bspline(Norm(vec_A2B), h) / w_Bspline(0., h);
      //                      if (B->isSurface)
      //                         tmp *= 10.;
      //                      A->DUDt_modify_SPH_2 -= tmp;
      //                      // local_DUDt_modify_SPH_2 -= c * vec_A2B * w_Bspline(Norm(vec_A2B), h) / max_w;
      //                   }
      //                }
      //             };

      //             for (const auto &net : target_nets) {
      //                if (A->isFluid) {
      //                   for (const auto &B : A->points_in_SML)
      //                      add(B);
      //                } else {
      //                   net->BucketPoints.apply(A->X, 1.1 * r, [&](const auto &B) {
      //                      if (canInteract(A, B))
      //                         add(B);
      //                   });
      //                }
      //             }

      //             // A->DUDt_modify_SPH_2 = local_DUDt_modify_SPH_2;
      //             // double max_scale = 0.001 * Norm(A->DUDt_SPH);
      //             // double original_scale = Norm(A->DUDt_modify_SPH_2);
      //             // if (original_scale > 0.)
      //             //    A->DUDt_modify_SPH_2 *= std::min(original_scale, max_scale) / original_scale;
      //          }
      //       }

      // #pragma omp parallel
      //       for (const auto &A : points)
      // #pragma omp single nowait
      //       {
      //          const double dt = A->RK_X.get_dt();

      //          if (A->getNetwork()->isRigidBody) {
      //             A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt + A->DUDt_modify_SPH_2;
      //             A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
      //             A->DUDt_SPH += _GRAVITY3_;
      //             A->DrhoDt_SPH = -A->rho * A->div_U;
      //          } else {
      //             A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt + A->DUDt_modify_SPH_2;
      //             A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
      //             A->DUDt_SPH += _GRAVITY3_;
      //             A->DrhoDt_SPH = -A->rho * A->div_U;
      //          }
      //       }

      // #pragma omp parallel
      //       for (const auto &A : points)
      // #pragma omp single nowait
      //       {
      //          /* ----------------------- calculate DUDt and b_vector ---------------------- */
      //          // \label{SPH:lapU_for_wall}
      //          // \label{SPH:Poisson_b_vector}
      //          // \label{SPH:how_to_set_fluid_b_vector}
      //          const double dt = A->RK_X.get_dt();
      //          if (A->getNetwork()->isRigidBody) {
      //             A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt + A->DUDt_modify_SPH_2;
      //             A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
      //             A->DUDt_SPH += _GRAVITY3_;
      //             A->DrhoDt_SPH = -A->rho * A->div_U;
      //          } else {
      //             A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt + A->DUDt_modify_SPH_2;
      //             A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
      //             A->DUDt_SPH += _GRAVITY3_;
      //             A->DrhoDt_SPH = -A->rho * A->div_U;
      //          }
      //       }

      /* -------------------------------------------------------------------------- */

#pragma omp parallel
      for (const auto &A : points)
#pragma omp single nowait
      {
         const double dt = A->RK_X.getTimeAtNextStep() - A->RK_X.getTimeAtCurrentStep();
         A->div_U_next = 0.;
         const auto U_A = U_next(A);
         if (A->isFluid) {
            // add(A->points_in_SML);
            for (const auto &B : A->points_in_SML)
               FusedMultiplyIncrement(-V_next(B), Dot(U_A - U_next(B), grad_w_Bspline_next(A, B)), A->div_U_next);
         } else
            for (const auto &net : target_nets) {
               net->BucketPoints.apply(X_next(A), 1.1 * A->SML(), [&](const auto &B) {
                  if (canInteract(A, B))
                     FusedMultiplyIncrement(-V_next(B), Dot(U_A - U_next(B), grad_w_Bspline_next(A, B)), A->div_U_next);
               });
            }
         A->DrhoDt_SPH_next = -A->rho * A->div_U_next;
      }

   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in calcLaplacianU");
   };
};

#endif