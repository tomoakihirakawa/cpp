#ifndef SPH_KERNEL_HELPER_FUNCTIONS_HPP
#define SPH_KERNEL_HELPER_FUNCTIONS_HPP

// # -------------------------------------------------------------------------- */

std::array<double, 3> grad_w_Bspline(const networkPoint *p, const networkPoint *q) {
#ifdef USE_GRAD_CORRECTION
   return Dot(p->inv_grad_corr_M, grad_w_Bspline(p->X, q->X, p->SML()));
#else
   return grad_w_Bspline(p->X, q->X, p->SML());
#endif
}

// # -------------------------------------------------------------------------- */

double Dot_grad_w_Bspline(const networkPoint *p, const networkPoint *q) {
#ifdef USE_LAPLACIAN_CORRECTION
   return Dot(p->inv_grad_corr_M, grad_w_Bspline(p->X, q->X, p->SML()));
#else
   return Dot_grad_w_Bspline(p->X, q->X, p->SML());
#endif
}

//! -------------------------------------------------------------------------- */

std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
#ifdef USE_GRAD_CORRECTION
   return Dot(p->inv_grad_corr_M_next, grad_w_Bspline(X_next(p), X_next(q), p->SML_next()));
#else
   return grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
#endif
}

std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const Tddd &X, const networkPoint *q) {
#ifdef USE_GRAD_CORRECTION
   return Dot(p->inv_grad_corr_M_next, grad_w_Bspline(X, X_next(q), p->SML_next()));
#else
   return grad_w_Bspline(X, X_next(q), p->SML_next());
#endif
}

//! -------------------------------------------------------------------------- */

double Dot_grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
#ifdef USE_LAPLACIAN_CORRECTION
   return Dot_grad_w_Bspline(X_next(p), X_next(q), p->SML_next(), p->laplacian_corr_M_next);
#else
   return Dot_grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
#endif
}

double Dot_grad_w_Bspline_next(const networkPoint *p, const Tddd &X, const networkPoint *q) {
#ifdef USE_LAPLACIAN_CORRECTION
   return Dot_grad_w_Bspline(X, X_next(q), p->SML_next(), p->laplacian_corr_M_next);
#else
   return Dot_grad_w_Bspline(X, X_next(q), p->SML_next());
#endif
}

#endif