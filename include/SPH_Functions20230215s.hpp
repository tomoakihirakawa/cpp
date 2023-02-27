#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"

/* -------------------------------------------------------------------------- */

Tddd SPP_X(const networkPoint *p) {
   // return p - (p->COM_SPH - p);
   // return p->X - 2 * p->COM_SPH;
   return p->X + p->interpolated_normal_SPH * p->radius_SPH / p->C_SML;
};

auto canSetSPP(const auto &target_nets, const auto &p) {

#if defined(USE_SPP_Fluid)
   auto X = SPP_X(p);
   auto range = p->radius_SPH / p->C_SML;
   if (!p->isSurface)
      return false;

   // auto d = Distance(p, X);

   // auto func = [&](const auto &q) { return p != q && Distance(q, X) < d; };

   // if (net->BucketPoints.any_of(X, d, func))
   //    return false;
   // for (const auto &[obj, poly] : RigidBodyObject)
   //    if (obj->BucketPoints.any_of(X, d, func))
   //       return false;

   const double C = 1.2;
   for (const auto &net : target_nets)
      if (net->BucketPoints.any_of(X, C * range,
                                   [&](const auto &q) { return p != q && Distance(q, X) < C * range; }))
         return false;

   // {
   //    const double C = 1.2;
   //    auto func = [&](const auto &q) { return p != q && Distance(q, X) < C * range; };
   //    if (net->BucketPoints.any_of(X, C * range, func))
   //       return false;
   // }
   // {
   //    const double C = 1.2;
   //    auto func = [&](const auto &q) { return p != q && Distance(q, X) < C * range; };
   //    for (const auto &[obj, poly] : RigidBodyObject)
   //       if (obj->BucketPoints.any_of(X, C * range, func))
   //          return false;
   // }
   return true;
#else
   return false;
#endif
};

// b$ ------------------------------------------------------ */
// b$                    ∇.∇UとU*を計算                       */
// b$ ------------------------------------------------------ */

auto Lap_U(const netPp A, const std::unordered_set<Network *> &target_nets) {
   // ここで初期化することは問題ない．この点pに対して再度計算することはないので，途中で初期化される心配はない．
   A->lap_U_ = A->lap_U = {0., 0., 0.};
   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH,
                              [&](const auto &B) {
                                 if (B->isCaptured) {
                                    //$ ------------------------------------------ */
                                    auto func = [&](const Tddd &qX, const double coef = 1.) {
                                       auto rij = qX - A->X;
                                       auto Uij = coef * B->U_SPH - A->U_SPH;
                                       auto nu_nu = B->mu_SPH / B->rho + A->mu_SPH / A->rho;
                                       A->lap_U_ += 1 / (A->mu_SPH / A->rho) * B->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(A->X, qX, A->radius_SPH) /
                                                    ((B->rho + A->rho) * Dot(rij, rij) + std::pow(1E-4 * A->radius_SPH / A->C_SML, 2));
                                    };
                                    //$ ------------------------------------------ */
                                    //$ ------------------------------------------ */
                                    // auto func = [&](const Tddd &qX) {
                                    //    auto rij = A->X - qX;
                                    //    auto Uij = A->U_SPH - coef *B->U_SPH;
                                    //    A->lap_U_ += 2 * B->mass / A->rho * Dot(rij, grad_w_Bspline(A->X, qX, A->radius_SPH)) / Dot(rij, rij) * Uij;
                                    // };
                                    //$ ------------------------------------------ */
                                    if (A != B)
                                       func(B->X);
#ifdef USE_SPP_Fluid
                                    if (B->isSurface) {
                                       // if (canSetSPP(target_nets, B))
                                       func(SPP_X(B), SPP_U_coef);
                                    }
#endif
                                 }
                              });
};

auto Lap_U(const auto &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
      Lap_U(A, target_nets);
};

auto setLap_U(const auto &points, const double dt) {
#pragma omp parallel
   for (const auto &A : points)
#pragma omp single nowait
   {
      A->lap_U = A->lap_U_;
      A->ViscousAndGravityForce = A->DUDt_SPH = A->DUDt_SPH_ = A->lap_U_ * (A->mu_SPH / A->rho) + _GRAVITY3_;  // 後で修正されるDUDt
      A->tmp_U_SPH = A->U_SPH + A->DUDt_SPH * dt;
      A->tmp_X = A->X + A->tmp_U_SPH * dt;
   }
};

// b$ ------------------------------------------------------ */
// b$      div(U^*),　DρDt=-ρdiv(U),　ラプラシアン∇.∇U* の計算   */
// b$ ------------------------------------------------------ */

void div_tmpU(const netPp A, const std::unordered_set<Network *> &target_nets) {
   A->div_tmpU = A->div_U = 0;  // this is div(U^*) not div(U^n)
   A->lap_tmpU = A->grad_div_U = {0., 0., 0.};
   Tddd Uij, rij;
   double nu_nu;
   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH,
                              [&](const auto &B) {
                                 if (B->isCaptured) {
                                    //@ ------------------------------------------ */
                                    // to use both A and spp
                                    auto func = [&](const auto &qX, const double coef = 1.) {
                                       if (Distance(A, qX) > 1E-12) {
                                          // 後藤p.25 (2.89)
                                          Uij = coef * B->U_SPH - A->U_SPH;
                                          A->div_U += B->mass / A->rho * Dot(Uij, grad_w_Bspline(A->X, qX, A->radius_SPH));
                                          //
                                          Uij = coef * B->tmp_U_SPH - A->tmp_U_SPH;
                                          A->div_tmpU += B->mass / A->rho * Dot(Uij, grad_w_Bspline(A->X, qX, A->radius_SPH));
                                          //
                                          rij = qX - A->X;
                                          nu_nu = B->mu_SPH / B->rho + A->mu_SPH / A->rho;
                                          A->lap_tmpU += 1 / (A->mu_SPH / A->rho) * B->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(A->X, qX, A->radius_SPH) /
                                                         ((B->rho + A->rho) * Dot(rij, rij) + std::pow(1E-4 * A->radius_SPH / A->C_SML, 2));
                                       }
                                    };
                                    //@ ------------------------------------------ */
                                    func(B->X);
#ifdef USE_SPP_Fluid
                                    if (B->isSurface) {
                                       // if (canSetSPP(target_nets, B))
                                       func(SPP_X(B), SPP_U_coef);
                                    }
#endif
                                 }
                              });
};

void div_tmpU(const auto &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
      div_tmpU(p, target_nets);
};

void div_tmpU(const auto &pointsA, const auto &pointsB, const std::unordered_set<Network *> &target_nets) {
   div_tmpU(pointsA, target_nets);
   div_tmpU(pointsB, target_nets);
};

auto set_tmpLap_U(const auto &points) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
   {
      p->tmp_ViscousAndGravityForce = p->lap_tmpU * (p->mu_SPH / p->rho) + _GRAVITY3_;  // 後で修正されるDUDt
   }
};

// b% -------------------------------------------------------------------------- */
// b% -------------------------------------------------------------------------- */
// #define Morikawa2019
#define NewtonMethod

void nextPressure(const netPp A, const std::unordered_set<Network *> &target_nets, const double dt) {
   double sum_Aij = 0, Aij, sum_Aij_Pj = 0, sum_Aij_Pij = 0, sum_Pij = 0, pressure_SPH = 0, qP;
   Tddd Xij;
   A->p_SPH_ = 0;
   A->gradP_SPH = A->gradP_SPH_ = {0., 0., 0.};
   //% -------------------------------------------------------------------------- */
   auto FUNC = [&](const auto &B) {
      if (B->isCaptured) {
         //% ------------------------------------------ */
         auto func = [&](const auto &qX, const double coef = 1.) {
            Xij = A->X - qX;
            qP = coef * B->p_SPH;
            if (Distance(A, qX) > 1E-12) {
               auto share = Dot(Xij, grad_w_Bspline(A->X, qX, A->radius_SPH)) / Dot(Xij, Xij);
#if defined(Morikawa2019)
               // Morikawa, D., Senadheera, H., & Asai, M. (2021). Explicit incompressible smoothed particle hydrodynamics in a multi-GPU environment for large-scale simulations. Computational Particle Mechanics, 8(3), 493–510. https://doi.org/10.1007/s40571-020-00347-0
               Aij = 2. * B->mass / A->rho * share;
#else
               // Nomeritae
               // Shao and Lo
               // auto Aij = B->mass * 8. / std::pow(B->rho + A->rho, 2) * share;, which is the same as
               Aij = 2. * B->mass / std::pow((B->rho + A->rho) / 2., 2) * share;
#endif

               sum_Aij += Aij;
               sum_Pij += A->p_SPH - qP;
               sum_Aij_Pj += Aij * qP;
               sum_Aij_Pij += Aij * (A->p_SPH - qP);

               // WALL
               A->gradP_SPH += A->rho * (qP / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * B->mass * grad_w_Bspline(A->X, qX, A->radius_SPH);
               A->gradP_SPH_ += A->rho * (1. / (A->rho * A->rho)) * B->mass * grad_w_Bspline(A->X, qX, A->radius_SPH);
            }
         };
         //% ------------------------------------------ */
         func(B->X);
#ifdef USE_SPP_Fluid
         if (B->isSurface)
            // if (canSetSPP(target_nets, B))
            func(SPP_X(B), SPP_p_coef);
#endif
      }
   };
   //% -------------------------------------------------------------------------- */
   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH, FUNC);
      // for (const auto &[obj, poly] : RigidBodyObject)
      //    obj->BucketPoints.apply(A->X, A->radius_SPH, FUNC);

#if defined(Morikawa2019)
   auto b = _WATER_DENSITY_ * A->div_tmpU / dt;
#else
   auto b = A->div_tmpU / dt;
#endif
   A->p_SPH_ = (b + sum_Aij_Pj) / sum_Aij;

#if defined(NewtonMethod)
   A->div_U_error = std::abs(sum_Aij_Pij - b);
   auto F = sum_Aij_Pij - b;
   auto dFdPi = sum_Aij;

#ifdef surface_zero_pressure
   if (A->isSurface) {
      F += A->p_SPH;
      sum_Aij += 1;
   }
#endif

   auto FF = F * F / 2.;
   auto dFFdPi = F * dFdPi;
   A->NR_pressure.update(FF, dFFdPi, 2.);

#endif
   if (!isFinite(A->p_SPH_) || !isFinite(sum_Aij_Pij)) {
      std::cout << "A->p_SPH_ = " << A->p_SPH_ << std::endl;
      std::cout << "A->div_tmpU = " << A->div_tmpU << std::endl;
      std::cout << "b = " << b << std::endl;
      std::cout << "sum_Aij_Pij = " << sum_Aij_Pj << std::endl;
      std::cout << "sum_Aij_Pj = " << sum_Aij_Pj << std::endl;
      std::cout << "sum_Aij = " << sum_Aij << std::endl;
      // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not a finite");

      A->p_SPH = A->p_SPH_ = 0;
   }
};

// b% ------------------------------------------------------ */
// b%           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
// b% ------------------------------------------------------ */

void gradP(const netPp A, const std::unordered_set<Network *> &target_nets) {
   A->contact_points_all_SPH = A->contact_points_fluid_SPH = 0;
   A->gradP_SPH = {0., 0., 0.};
   auto FUNC = [&](const auto &B) {
      if (B->isCaptured) {
         //% ------------------------------------------ */
         auto func = [&](const auto &qX, const double coef = 1.) {
            if (Distance(A, qX) > 1E-12) {
               // double p_rho2 = coef * B->p_SPH / std::pow(B->rho, 2) + A->p_SPH / std::pow(A->rho, 2);
               // auto qP = (spp ? SPP_p_coef * B->p_SPH_SPP : B->p_SPH);
               // A->gradP_SPH += (qP + A->p_SPH) * B->mass / _WATER_DENSITY_ * grad_w_Bspline(A->X, qX, A->radius_SPH);
               auto qP = coef * B->p_SPH;
               A->gradP_SPH += A->rho * (qP / (B->rho * B->rho) + A->p_SPH / (A->rho * A->rho)) * B->mass * grad_w_Bspline(A->X, qX, A->radius_SPH);
            }
         };
         //% ------------------------------------------ */
         func(B->X);
#ifdef USE_SPP_Fluid
         if (B->isSurface) {
            // if (canSetSPP(target_nets, B))
            func(SPP_X(B), SPP_p_coef);
         }
#endif
      }
      A->contact_points_all_SPH++;
      if (B->getNetwork()->isFluid || B->isFluid)
         A->contact_points_fluid_SPH++;
   };

   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH, FUNC);

   // A->DUDt_SPH -= A->gradP_SPH / A->rho;
   A->DUDt_SPH -= A->gradP_SPH / _WATER_DENSITY_;
   //
   if (!isFinite(A->DUDt_SPH))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
};

#endif