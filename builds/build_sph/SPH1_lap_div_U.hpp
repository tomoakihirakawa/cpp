#ifndef SPH_lap_div_U_H
#define SPH_lap_div_U_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/*DOC_EXTRACT SPH

## (1) $`\nabla^2 {\bf u}_i`$の計算（`calcLaplacianU`）

CHECKED: \ref{SPH:lapU}{流速のラプラシアンの計算方法}: $`\nabla^2 {\bf u}_i=\sum_{j} A_{ij}({\bf u}_i - {\bf u}_j),\quad A_{ij} = \frac{2m_j}{\rho_i}\frac{{{\bf x}_{ij}}\cdot\nabla W_{ij}}{{\bf x}_{ij}^2}`$

CHECKED: \ref{SPH:divU}{流速の発散の計算方法}: $`\nabla\cdot{\bf u}_i=\sum_{j}\frac{m_j}{\rho_j}({{\bf u}_j-{\bf u}_i}) \cdot\nabla W_{ij}`$

*/

auto calcLaplacianU(const auto &points, const std::unordered_set<Network *> &target_nets) {

   //! -------------------------------------------------------------------------- */

#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {

      A->grad_corr_M.fill({0., 0., 0.});
      A->grad_corr_M_next.fill({0., 0., 0.});
      A->grad_U.fill({0., 0., 0.});
      A->Mat1.fill({0., 0., 0.});
      A->Mat2.fill({0., 0., 0.});
      A->Mat3.fill({0., 0., 0.});

      Tddd Xij, nz_Xij;
      auto add = [&](const auto &B) {
         auto grad_w = grad_w_Bspline(A->X, B->X, A->radius_SPH);
         A->grad_corr_M += B->volume * TensorProduct(B->X - A->X, grad_w);
         A->grad_corr_M_next += V_next(B) * TensorProduct(X_next(B) - X_next(A), grad_w_Bspline(X_next(A), X_next(B), A->radius_SPH));

         const auto Uij = A->U_SPH - B->U_SPH;
         A->grad_U += B->volume * TensorProduct(-Uij, grad_w);

         Xij = A->X - B->X;
         nz_Xij = Normalize(Xij);
         A->Mat1 += B->volume * TensorProduct(Xij * nz_Xij * nz_Xij, grad_w);
         A->Mat2 += B->volume * TensorProduct(nz_Xij * nz_Xij, grad_w);
         A->Mat3 += B->volume * TensorProduct(Xij * Xij, grad_w);

#if defined(USE_MIRROR_PARTICLE)
         for (auto i = 0; i < B->vector_to_polygon.size(); ++i) {
            auto X = B->X + 2. * B->vector_to_polygon[i];
            auto Xnext = B->X + 2. * B->vector_to_polygon_next[i];
            A->grad_corr_M += B->volume * TensorProduct(X - A->X, grad_w_Bspline(A->X, X, A->radius_SPH));
            A->grad_corr_M_next += V_next(B) * TensorProduct(Xnext - X_next(A), grad_w_Bspline(X_next(A), Xnext, A->radius_SPH));
         }
#endif
      };

      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) add(B);
         });

      A->inv_grad_corr_M = Inverse(A->grad_corr_M);
      A->inv_grad_corr_M_next = Inverse(A->grad_corr_M_next);
      T3Tddd I;
      IdentityMatrix(I);
      A->Mat_B = Dot(-I, Inverse(A->Mat1 + Dot(Dot(A->Mat2, A->inv_grad_corr_M), A->Mat3)));

      /* --------------------------------------------------------- */

      {
         auto add_to_map = [&](const auto &p, const auto &c) {
            if (A->map_p_grad.find(p) == p->map_p_grad.end())
               A->map_p_grad[p] = c;
            else
               A->map_p_grad[p] += c;
         };

         auto add = [&](const auto &B) {
            auto c = B->volume * Dot(A->X - B->X, Dot(grad_w_Bspline(A->X, B->X, A->radius_SPH), A->inv_grad_corr_M));
            add_to_map(B, -c);
            add_to_map(A, c);
            // 符号についての疑問が残る．逆になっているように思える．
         };

         A->map_p_grad.clear();
         for (const auto &net : target_nets)
            net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
               if (B->isCaptured) add(B);
            });

         // copy    std::unordered_map<networkPoint *, double> map_p_grad to  std::vector<std::tuple<networkPoint *, double>> vector_p_grad;
         A->vector_p_grad.clear();
         for (const auto &[p, v] : A->map_p_grad)
            A->vector_p_grad.push_back(std::tuple<networkPoint *, double>{p, v});
      }
   }

   //! -------------------------------------------------------------------------- */

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
      A->grad_coeff.clear();
      A->grad_coeff_next.clear();
      // A->grad_U.fill({0., 0., 0.});
      double Aij;
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

#if defined(USE_MIRROR_PARTICLE)
         // 綺麗に配列させ設置した壁粒子の圧力を，内部の流体粒子の圧力から決める際に，理想的な値を設定することが困難と思われる．
         for (auto i = 0; i < B->vector_to_polygon.size(); ++i) {
            auto v = B->vector_to_polygon[i];
            auto v_next = B->vector_to_polygon_next[i];
            auto X = B->X + 2. * v;
            auto Xnext = B->X + 2. * v_next;
            auto U = Reflect(B->U_SPH, v);

            const auto Uij = A->U_SPH - U;
            A->density_based_on_positions += B->rho * B->volume * w_Bspline(Norm(A->X - X), A->radius_SPH);
            // A->div_U += B->volume * Dot(-Uij, grad_w_Bspline(A->X, X, A->radius_SPH));         //\label{SPH:divU}
            // A->lap_U += 2 * B->volume * Uij * Dot_grad_w_Bspline_Dot(A->X, X, A->radius_SPH);  //\label{SPH:lapU}

            A->div_U += B->volume * Dot(-Uij, Dot(grad_w_Bspline(A->X, X, A->radius_SPH), A->inv_grad_corr_M));  //\label{SPH:divU}
            A->lap_U += 2 * B->volume * Uij * Dot_grad_w_Bspline_Dot(A->X, X, A->radius_SPH);                    //\label{SPH:lapU}

            // A->lap_U += -8. * B->mass / (A->rho + B->rho) * Dot(Uij, A->X - X) / Dot(A->X - X, A->X - X) * grad_w_Bspline(A->X, X, A->radius_SPH);  //\label{SPH:lapP2}
         }
#endif

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
      networkPoint *closest_surface_point = nullptr;
      double min_distance = 1e10, distance;

      for (const auto &net : target_nets)
         net->BucketPoints.apply(A->X, A->radius_SPH, [&](const auto &B) {
            if (B->isCaptured) {
               add(B);

               if (B->isSurface) {
#if defined(USE_ONE_AUXP)
                  if ((distance = Distance(A, B)) < min_distance) {
                     min_distance = distance;
                     closest_surface_point = B;
                  }
#elif defined(USE_ALL_AUXP)
                   for (const auto &AUX : B->auxiliaryPoints){
                      if (AUX != nullptr)
                        add(AUX);
                   }
#endif
               }
            }
         });
//
#if defined(USE_ONE_AUXP)
      // if (A->isSurface)
      if (closest_surface_point != nullptr)
         for (const auto &AUX : closest_surface_point->auxiliaryPoints)
            if (AUX != nullptr)
               add(AUX);
#endif

#if defined(USE_SIMPLE_SINGLE_AUX)
      if (A->isSurface)
         for (const auto &AUX : A->auxiliaryPoints)
            if (AUX != nullptr)
               PoissonEquation(AUX);
#endif

      // A->lap_U = Dot(A->Mat_B, A->lap_U);
      //$ ------------------------------------------ */
      // \label{SPH:lapU_for_wall}
      // \label{SPH:Poisson_b_vector}
      if (A->isAir) {
         A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;最も自然な結果を返す
         A->DrhoDt_SPH = 0;
      } else if (A->getNetwork()->isRigidBody) {
         A->DUDt_SPH_.fill(0.);
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH.fill(0.);
         A->tmp_U_SPH.fill(0.);
         A->tmp_X = A->X;
         A->DrhoDt_SPH = 0.;
         // \label{SPH:how_to_set_wall_b_vector}
         A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;最も自然な結果を返す

         // if (A->isFirstWallLayer)
         // A->b_vector = A->U_SPH / dt + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
         // else {
         // A->lap_U.fill(0.);
         // A->b_vector.fill(0.);
         // }
      } else {
         A->DUDt_SPH_ = A->DUDt_SPH;
         double nu = A->mu_SPH / A->rho;
         A->DUDt_SPH = nu * A->lap_U + _GRAVITY3_;  // 後で修正されるDUDt

         // for (const auto &v : A->vector_to_polygon) {
         //    auto half_ps = A->radius_SPH / A->C_SML / 2;
         //    if (Norm(v) < half_ps) {
         //       if (Dot(A->U_SPH, v) > 0) {
         //          // auto alpha = (half_ps - Norm(v)) / half_ps;
         //          auto alpha = half_ps / Norm(v) - 1;
         //          A->DUDt_SPH += Normalize(v) * _GRAVITY3_ * alpha;
         //          // U = Chop(U, v) + Projection(U, v) * (2. * Norm(v) / half_ps - 1.);
         //          // U += Projection(U, v) * (half_ps / Norm(v));
         //       }
         //    }
         // }

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
         //             std::array<double, 3> U1, U2, U3, U4;
         //             U1 = A->U_SPH;
         //             if (*(A->vec_time_SPH.rbegin()) == current_time) {
         //                times.push_back(*(A->vec_time_SPH.rbegin() + 1));
         //                U2 = *(A->vec_U_SPH.rbegin() + 1);
         //                times.push_back(*(A->vec_time_SPH.rbegin() + 2));
         //                U3 = *(A->vec_U_SPH.rbegin() + 2);
         //                times.push_back(*(A->vec_time_SPH.rbegin() + 3));
         //                U4 = *(A->vec_U_SPH.rbegin() + 3);
         //             } else {
         //                times.push_back(*(A->vec_time_SPH.rbegin() + 0));
         //                U2 = *(A->vec_U_SPH.rbegin() + 0);
         //                times.push_back(*(A->vec_time_SPH.rbegin() + 1));
         //                U3 = *(A->vec_U_SPH.rbegin() + 1);
         //                times.push_back(*(A->vec_time_SPH.rbegin() + 2));
         //                U4 = *(A->vec_U_SPH.rbegin() + 2);
         //             }
         //             InterpolationLagrange<double> lag(times);
         //             auto D = lag.DN(current_time);
         //             A->b_vector = -(D[1] * U1 + D[2] * U2 + D[3] * U3 + D[4] * U4) + A->mu_SPH / A->rho * A->lap_U;  // + _GRAVITY3_;
         //          }
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