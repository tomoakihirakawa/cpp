#ifndef SPH_Auxiliary_H
#define SPH_Auxiliary_H

#include "Network.hpp"
#include "minMaxOfFunctions.hpp"

/* -------------------------------------------------------------------------- */
void deleteAuxiliaryPoints(const auto net) {
   std::cout << "delete auxiliary points" << std::endl;
   auto points = net->getPoints();
   for (auto p : points) {
      p->auxPoint = nullptr;
      if (p->isAuxiliary) {
         p->surfacePoint->auxPoint = nullptr;
         delete p;
      }
   }

   for (auto p : net->getPoints())
      p->isAuxiliary = false;

   net->remakeBucketPoints();
   std::cout << "delete auxiliary points done" << std::endl;
}

Tddd aux_position(const networkPoint *p) {
   auto sp = p->surfacePoint;
#ifdef SET_AUX_AT_PARTICLE_SPACING
   return sp->X + sp->particle_spacing * Normalize(Dot(sp->inv_grad_corr_M, sp->interp_normal_original));
#elif defined(SET_AUX_AT_MASS_CENTER)
   return sp->X - sp->vec2COM;
#else
   return sp->X + sp->particle_spacing * Normalize(Dot(sp->inv_grad_corr_M, sp->interp_normal_original));
#endif
};

// void setAuxiliaryPoints(const auto net, const std::function<Tddd(networkPoint *)> &position = aux_position) {
//    Print("setAuxiliaryPoints");
//    auto points = net->getPoints();
//    for (const auto &p : points) {
//       if (p->hasAuxiliary()) {

//          auto q = new networkPoint(net, p->X);
//          q->surfacePoint = p;
//          p->auxPoint = q;
//          q->setX(position(q));
//          //
//          q->grad_corr_M = p->grad_corr_M;
//          q->grad_corr_M_next = p->grad_corr_M_next;
//          q->inv_grad_corr_M = p->inv_grad_corr_M_next;
//          q->inv_grad_corr_M_next = p->inv_grad_corr_M_next;
//          //
//          // q->grad_corr_M_rigid = p->grad_corr_M_next_rigid;
//          // q->grad_corr_M_next_rigid = p->grad_corr_M_next_rigid;
//          // q->inv_grad_corr_M_rigid = p->inv_grad_corr_M_next_rigid;
//          // q->inv_grad_corr_M_next_rigid = p->inv_grad_corr_M_next_rigid;
//          q->laplacian_corr_M = p->laplacian_corr_M;
//          q->laplacian_corr_M_next = p->laplacian_corr_M_next;
//          //
//          q->isSurface = p->isSurface;
//          q->isSurface_next = p->isSurface_next;
//          q->isNeumannSurface = p->isNeumannSurface;
//          q->isAuxiliary = true;
//          //
//          q->b_vector = p->b_vector;
//          q->U_SPH = p->U_SPH;  // - 2 * Chop(p->U_SPH, p->interp_normal);
//          //
//          q->intp_density = p->intp_density;
//          q->intp_density_next = p->intp_density_next;
//          //
//          // q->U_SPH.fill(0.);
//          q->v_to_surface_SPH = p->v_to_surface_SPH;
//          q->interp_normal = p->interp_normal;
//          q->interp_normal_next = p->interp_normal_next;
//          q->intp_normal_Eigen = p->intp_normal_Eigen;
//          q->interp_normal_original = p->interp_normal_original;
//          q->interp_normal_original_next = p->interp_normal_original_next;
//          q->intp_density = p->intp_density;
//          //
//          q->div_U = p->div_U;
//          q->DUDt_SPH = p->DUDt_SPH;
//          q->lap_U = p->lap_U;
//          q->p_SPH = p->p_SPH;
//          q->rho = p->rho;
//          q->setDensityVolume(_WATER_DENSITY_, p->volume);
//          q->particle_spacing = p->particle_spacing;
//          q->C_SML_next = p->C_SML_next;
//          q->C_SML = p->C_SML;
//          q->isFluid = p->isFluid;
//          q->isFirstWallLayer = false;
//          q->isCaptured = true;
//          //
//          auto dt = p->RK_X.getdt();
//          // q->RK_U.initialize(dt, simulation_time, q->U_SPH, 1);
//          // q->RK_X.initialize(dt, simulation_time, q->X, 1);
//          // q->RK_P.initialize(dt, simulation_time, q->p_SPH, 1);
//          // q->RK_rho.initialize(dt, simulation_time, q->rho, 1);
//          q->RK_U = p->RK_U;
//          q->RK_U.Xinit = p->RK_U.Xinit;
//          q->RK_U.t_init = p->RK_U.t_init;
//          q->RK_U.dt_fixed = p->RK_U.dt_fixed;
//          q->RK_U.dt = p->RK_U.dt;
//          q->RK_U.steps = p->RK_U.steps;
//          q->RK_U.current_step = p->RK_U.current_step;
//          q->RK_U._dX = p->RK_U._dX;
//          q->RK_U.dX = p->RK_U.dX;
//          //
//          q->RK_X = p->RK_X;
//          // q->RK_X.Xinit = p->RK_X.Xinit;
//          q->RK_X.Xinit = q->X;
//          q->RK_X.t_init = p->RK_X.t_init;
//          q->RK_X.dt_fixed = p->RK_X.dt_fixed;
//          q->RK_X.dt = p->RK_X.dt;
//          q->RK_X.steps = p->RK_X.steps;
//          q->RK_X.current_step = p->RK_X.current_step;
//          q->RK_X._dX = p->RK_X._dX;
//          q->RK_X.dX = p->RK_X.dX;
//          //
//          q->RK_rho = p->RK_rho;
//          q->RK_rho.Xinit = p->RK_rho.Xinit;
//          q->RK_rho.t_init = p->RK_rho.t_init;
//          q->RK_rho.dt_fixed = p->RK_rho.dt_fixed;
//          q->RK_rho.dt = p->RK_rho.dt;
//          q->RK_rho.steps = p->RK_rho.steps;
//          q->RK_rho.current_step = p->RK_rho.current_step;
//          q->RK_rho._dX = p->RK_rho._dX;
//          q->RK_rho.dX = p->RK_rho.dX;
//       }
//    }
//    Print("setAuxiliaryPoints done");
//    Print("remakeBucketPoints");
//    net->remakeBucketPoints();
//    Print("remakeBucketPoints done");
// }
#endif