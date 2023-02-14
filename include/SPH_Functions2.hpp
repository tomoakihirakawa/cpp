#ifndef SPH_Functions_H
#define SPH_Functions_H

#include "Network.hpp"

// b$ ------------------------------------------------------ */
// b$                    ∇.∇UとU*を計算                       */
// b$ ------------------------------------------------------ */

auto init_Lap_U(const auto &points) {
   for (const auto &p : points)
      p->lap_U_ = p->lap_U = {0., 0., 0.};
};

auto Lap_U(const netPp p, const std::unordered_set<Network *> &target_nets) {
   for (const auto &net : target_nets)
      net->BucketPoints.apply(p->X, p->radius_SPH,
                              [&](const auto &q) {
                                 if (q->isCaptured) {
                                    //$ ------------------------------------------ */
                                    auto func = [&](const Tddd &qX, const double coef = 1.) {
                                       auto rij = qX - p->X;
                                       auto Uij = coef * q->U_SPH - p->U_SPH;
                                       auto nu_nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
                                       p->lap_U_ += 1 / (p->mu_SPH / p->rho) * q->mass * 8 * nu_nu * Dot(Uij, rij) * grad_w_Bspline(p->X, qX, p->radius_SPH) /
                                                    ((q->rho + p->rho) * Dot(rij, rij) + std::pow(1E-4 * p->radius_SPH / p->C_SML, 2));
                                    };
                                    //$ ------------------------------------------ */
                                    //$ ------------------------------------------ */
                                    // auto func = [&](const Tddd &qX) {
                                    //    auto rij = p->X - qX;
                                    //    auto Uij = p->U_SPH - coef *q->U_SPH;
                                    //    p->lap_U_ += 2 * q->mass / p->rho * Dot(rij, grad_w_Bspline(p->X, qX, p->radius_SPH)) / Dot(rij, rij) * Uij;
                                    // };
                                    //$ ------------------------------------------ */
                                    if (p != q)
                                       func(q->X);
#ifdef USE_SPP_Fluid
                                    if (q->isSurface) {
                                       if (canSetSPP(net, RigidBodyObject, q))
                                          func(SPP_X(q), SPP_U_coef);
                                    }
#endif
                                 }
                              });
};

auto Lap_U(const auto &points, const std::unordered_set<Network *> &target_nets) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
      Lap_U(p, target_nets);
};

auto setLap_U(const auto &points, const double dt) {
#pragma omp parallel
   for (const auto &p : points)
#pragma omp single nowait
   {
      p->lap_U = p->lap_U_;
      p->ViscousAndGravityForce = p->DUDt_SPH = p->DUDt_SPH_ = p->lap_U_ * (p->mu_SPH / p->rho) + _GRAVITY3_;  // 後で修正されるDUDt
      p->tmp_U_SPH = p->U_SPH + p->DUDt_SPH * dt;
      p->tmp_X = p->X + p->tmp_U_SPH * dt;
   }
};

// b$ ------------------------------------------------------ */
// b$      div(U^*),　DρDt=-ρdiv(U),　ラプラシアン∇.∇U* の計算   */
// b$ ------------------------------------------------------ */

void init_div_tmpU(const auto &points) {
   for (const auto &A : points) {
      A->div_tmpU = A->div_U = 0;  // this is div(U^*) not div(U^n)
      A->lap_tmpU = A->grad_div_U = {0., 0., 0.};
   }
};

void div_tmpU(const netPp A, const std::unordered_set<Network *> &target_nets) {
   Tddd Uij, rij;
   double nu_nu;
   for (const auto &net : target_nets)
      net->BucketPoints.apply(A->X, A->radius_SPH,
                              [&](const auto &B) {
                                 if (B->isCaptured) {
                                    //@ ------------------------------------------ */
                                    // to use both p and spp
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
                                       if (canSetSPP(net, RigidBodyObject, B))
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

#endif