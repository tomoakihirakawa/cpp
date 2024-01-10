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
#if defined(USE_RungeKutta)
            const double dt = A->RK_X.get_dt();
#elif defined(USE_LeapFrog)
            const double dt = A->LPFG_X.get_dt();
#endif
            A->div_U = 0.;
            A->lap_U.fill(0.);
            A->b_vector.fill(0.);
            A->U_XSPH.fill(0.);
            A->U_next.fill(0.);
            A->DUDt_modify_SPH.fill(0.);
            A->DUDt_SPH.fill(0.);
            A->DrhoDt_SPH = 0.;

            //! ラプラシアンの修正行列の一部を計算する
            Fill(A->tensorproduct_grad_Uij, 0.);
            Fill(A->tensorproduct_grad_Uij_next, 0.);
            // A->grad_U.fill({0., 0., 0.});
            const auto r = A->SML();
            const auto center = A->X;
            auto applyOverPoints = [&](const auto &equation) {
               for (const auto &net : target_nets) {
                  net->BucketPoints.apply(center, 1.2 * r, [&](const auto &B) {
                     if (canInteract(A, B))
                        equation(B);
                  });
               }
            };

            applyOverPoints([&](const auto &Q) {
               if (canInteract(A, Q) && A != Q) {
                  A->tensorproduct_grad_Uij += TensorProduct(Q->volume * grad_w_Bspline(A, Q), Q->U_SPH - A->U_SPH);
                  A->tensorproduct_grad_Uij_next += TensorProduct(V_next(Q) * grad_w_Bspline_next(A, Q), U_next(Q) - U_next(A));

                  // A->tensorproduct_grad_Uij += TensorProduct(Q->volume * grad_w_Bspline(A->X, Q->X, A->SML()), Q->U_SPH - A->U_SPH);
                  // A->tensorproduct_grad_Uij_next += TensorProduct(V_next(Q) * grad_w_Bspline(X_next(A), X_next(Q), A->SML_next()), U_next(Q) - U_next(A));
               }
            });

            const double c_xsph = 0.02;  // 少しは必要のようだ
            const double csml_factor = 1.;
            double total_w = 0, w, vol_w;
            double total_w_U_next = 0;
            auto add = [&](const auto &B) {
               if (A->isAuxiliary) {
                  // if (!B->isAuxiliary && !B->isSurface) {
                  if (!B->isAuxiliary) {
                     w = B->volume * w_Bspline(Norm(A->X - B->X), A->SML());
                     total_w_U_next += w;
                     A->U_next += B->U_SPH * w;
                  }
               }

               const auto Uij = A->U_SPH - B->U_SPH;

               if (A->isFluid && B->isFluid && !B->isAuxiliary && A != B) {
                  // if (((A->isFluid && B->isFluid) || (A->isFluid && B->isFirstWallLayer)) && !B->isAuxiliary && A != B) {
                  vol_w = B->volume * w_Bspline(Norm(A->X - B->X), A->particle_spacing * 2.1);
                  // if (B->isSurface)
                  //    A->U_XSPH += 2 * c_xsph * (-Uij) * vol_w;  //\label{SPH:U_XSPH}
                  // else
                  A->U_XSPH += c_xsph * (-Uij) * vol_w;  //\label{SPH:U_XSPH}
                  total_w += vol_w;
               }

               // A->div_U += B->volume * Dot(-Uij, grad_w_Bspline(A, B));       //\label{SPH:divU}
               FusedMultiplyIncrement(B->volume, -Dot(Uij, grad_w_Bspline(A, B)), A->div_U);
               const double Aij = 2. * B->volume * Dot_grad_w_Bspline(A, B);  //\label{SPH:lapP1}
               // const double Aij = 2. * B->volume * Dot_grad_w_Bspline(A->X, B->X, A->SML());

               //! 修正
               const auto DelX = (A->X - B->X);
#ifdef USE_PRE_CALC_tensorproduct_grad_Uij
               // A->lap_U -= Aij * Dot(DelX, A->tensorproduct_grad_Uij);
               if (!A->is_grad_corr_M_singular)
                  FusedMultiplyIncrement(-Aij, Dot(DelX, A->tensorproduct_grad_Uij), A->lap_U);
#else
               applyOverPoints([&](const auto &Q) {
                  if (A != Q)
                     A->lap_U -= Aij * Q->volume * Dot(DelX, grad_w_Bspline(A, Q)) * (Q->U_SPH - A->U_SPH);
                  // A->lap_U -= Aij * Q->volume * Dot(DelX, grad_w_Bspline(A->X, Q->X, A->SML())) * (Q->U_SPH - A->U_SPH);
               });
#endif
               // A->lap_U += Aij * Uij;  //\label{SPH:lapU}
               FusedMultiplyIncrement(Aij, Uij, A->lap_U);
            };

            applyOverPoints(add);

            if (total_w_U_next > too_small_total_w)
               A->U_next /= total_w_U_next;
            // if (total_w != 0.)
            //    A->U_XSPH /= total_w;
         }

#pragma omp parallel
         for (const auto &A : points)
#pragma omp single nowait
         {
            /* ----------------------- calculate DUDt and b_vector ---------------------- */
            // \label{SPH:lapU_for_wall}
            // \label{SPH:Poisson_b_vector}
            // \label{SPH:how_to_set_fluid_b_vector}
            const double dt = A->RK_X.get_dt();
            if (A->getNetwork()->isRigidBody) {
               A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt;  // + A->DUDt_modify_SPH;
               A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
               A->DUDt_SPH += _GRAVITY3_;
               A->DrhoDt_SPH = -A->rho * A->div_U;
            } else {
               A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt;  // + A->DUDt_modify_SPH;
               A->b_vector = A->rho * (A->U_SPH / dt + A->DUDt_SPH);
               A->DUDt_SPH += _GRAVITY3_;
               A->DrhoDt_SPH = -A->rho * A->div_U;
            }
            if (!isFinite(A->DUDt_SPH)) {
               // show details
               std::cout << "A->DUDt_SPH = " << A->DUDt_SPH << std::endl;
               std::cout << "A->mu_SPH = " << A->mu_SPH << std::endl;
               std::cout << "A->rho = " << A->rho << std::endl;
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
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
            }
         }
      }
      //
      // b% ------------------------------------ さらなる修正 ------------------------------------ */

      // #pragma omp parallel
      //       for (const auto &A : points)
      // #pragma omp single nowait
      //       {
      //          A->DUDt_modify_SPH.fill(0.);
      //          A->DUDt_modify_SPH_2.fill(0.);
      //          if (A->isFluid) {
      // #if defined(USE_RungeKutta)
      //             const double dt = A->RK_X.get_dt();
      // #elif defined(USE_LeapFrog)
      //             const double dt = A->LPFG_X.get_dt();
      // #endif

      //             const auto r = A->SML();
      //             const auto center = A->X;

      //             /* -------------------------------------------------------------------------- */
      //             A->interp_normal_original_next.fill(0.);
      //             Tddd interp_normal_original_next_mod = {0., 0., 0.};
      //             for (const auto &net : target_nets)
      //                net->BucketPoints.apply(A->X, 1.2 * A->SML_next(), [&](const auto &q) {
      //                   A->interp_normal_original_next -= rho_next(q) * V_next(q) * grad_w_Bspline(X_next(A), X_next(q), A->SML_next());
      //                   interp_normal_original_next_mod -= rho_next(q) * V_next(q) * grad_w_Bspline_next(A, q);
      //                });
      //             /* -------------------------------------------------------------------------- */

      //             // 修正１

      //             // const double typical_velocity = A->particle_spacing / dt;
      //             // const double typical_acceleration = typical_velocity / dt;
      //             // const auto V_org = A->interp_normal_original_next / _WATER_DENSITY_;  // non-dimensional
      //             // if (A->isSurface) {
      //             //    const auto V_mod = interp_normal_original_next_mod / _WATER_DENSITY_;
      //             //    const auto [Vn, Vs] = DecomposeVector(V_org, V_mod);
      //             //    A->DUDt_modify_SPH = dt * 0.001 * Vs * typical_acceleration;
      //             // } else {
      //             //    A->DUDt_modify_SPH = dt * 0.001 * V_org * typical_acceleration;
      //             // }

      //             // 修正２
      //             auto applyOverPoints = [&](const auto &equation) {
      //                for (const auto &net : target_nets) {
      //                   net->BucketPoints.apply(center, 1.2 * r, [&](const auto &B) {
      //                      if (canInteract(A, B))
      //                         equation(B);
      //                   });
      //                }
      //             };
      //             // 次はこれで試す
      //             auto add = [&](const auto &B) {
      //                if (A == B || A->isAuxiliary || B->isAuxiliary || !A->isFluid)
      //                   return;
      //                if (A->isFluid && B->isFluid) {
      //                   const auto vec_A2B = X_next(B) - X_next(A);
      //                   // const auto vec_A2B = B->X - A->X;
      //                   const double ratio = Norm(vec_A2B) / A->particle_spacing;
      //                   const auto relative_velocity = U_next(A) - U_next(B);  //! この意味はAがBに近づく速度
      //                   if (ratio < 0.9 && Dot(relative_velocity, vec_A2B) < 0) {
      //                      // if (A->isSurface) {
      //                      // auto [vec_A2B_normal, vec_A2B_tangential] = DecomposeVector(vec_A2B, Normalize(Dot(A->inv_grad_corr_M, A->interp_normal_original)));
      //                      //    A->DUDt_modify_SPH -= 0.0001 * Projection(relative_velocity, vec_A2B_tangential) / dt;
      //                      // } else
      //                      // {
      //                      // A->DUDt_modify_SPH -= 0.1 * (1. - ratio) * Projection(relative_velocity, vec_A2B) / dt;
      //                      A->DUDt_modify_SPH_2 -= 0.01 * Projection(relative_velocity, vec_A2B) / dt;
      //                      // }
      //                   }
      //                }
      //             };

      //             applyOverPoints(add);

      //             // suppress the change
      //             double max_scale = 0.1 * Norm(A->DUDt_SPH);
      //             double original_scale = Norm(A->DUDt_modify_SPH_2);
      //             if (original_scale > 0.)
      //                A->DUDt_modify_SPH_2 *= std::min(original_scale, max_scale) / original_scale;
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

      // b$ --------------------------------- 壁面での反射 --------------------------------- */
      //! 壁の反射を先に取り込めば計算精度が上がるかもしれない
      //       {
      // #pragma omp parallel
      //          for (const auto &A : points)
      // #pragma omp single nowait
      //             if (A->isFluid) {
      // #if defined(USE_RungeKutta)
      //                const double dt = A->RK_X.get_dt();
      // #elif defined(USE_LeapFrog)
      //                const double dt = A->LPFG_X.get_dt();
      // #endif
      //                A->DUDt_modify_SPH.fill(0.);
      //                auto p = A;
      //                auto d0 = A->particle_spacing;
      //                for (const auto &[closest_p, f2w] : {/*closest(p, RigidBodyObject) ,*/ closest_next(p, RigidBodyObject)}) {
      //                   if (closest_p != nullptr) {
      //                      auto dist_f2w = Norm(f2w);
      //                      auto n = closest_p->interp_normal;
      //                      auto normal_dist_f2w = Norm(Projection(f2w, n));  //! nomal distance from the fluid particle to the wall particle
      //                      auto ratio = (d0 - normal_dist_f2w) / d0;
      //                      if (dist_f2w < std::sqrt(2.) * A->particle_spacing && normal_dist_f2w < A->particle_spacing) {
      //                         if (Dot(U_next(p), n) < 0) {
      //                            auto tmp = -0.1 * Projection(U_next(p), n) / p->RK_X.get_dt();
      //                            p->DUDt_modify_SPH += tmp;
      //                         }
      //                      }
      //                   }
      //                }
      //             }

      // #pragma omp parallel
      //          for (const auto &A : points)
      // #pragma omp single nowait
      //          {
      // #if defined(USE_RungeKutta)
      //             const double dt = A->RK_X.get_dt();
      // #elif defined(USE_LeapFrog)
      //             const double dt = A->LPFG_X.get_dt();
      // #endif
      //             if (A->getNetwork()->isRigidBody) {
      //                A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt + A->DUDt_modify_SPH;
      //                A->b_vector = FusedMultiplyAdd(A->rho / dt, (A->U_SPH + A->U_XSPH), A->mu_SPH * A->lap_U);
      //                A->DUDt_SPH += _GRAVITY3_;
      //                A->DrhoDt_SPH = -A->rho * A->div_U;
      //             } else {
      //                A->DUDt_SPH = A->mu_SPH / A->rho * A->lap_U + A->U_XSPH / dt + A->DUDt_modify_SPH;
      //                A->b_vector = FusedMultiplyAdd(A->rho / dt, (A->U_SPH + A->U_XSPH), A->mu_SPH * A->lap_U);
      //                A->DUDt_SPH += _GRAVITY3_;
      //                A->DrhoDt_SPH = -A->rho * A->div_U;
      //             }
      //          }
      //          // setCorrectionMatrix(target_nets);
      //       }
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in calcLaplacianU");
   };
};

#endif