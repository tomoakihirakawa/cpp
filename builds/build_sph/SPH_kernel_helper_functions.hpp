#ifndef SPH_KERNEL_HELPER_FUNCTIONS_HPP
#define SPH_KERNEL_HELPER_FUNCTIONS_HPP

// # -------------------------------------------------------------------------- */
std::array<double, 3> grad_w_Bspline_helper(const networkPoint *p, const Tddd &pX, const Tddd &qX) {
#ifdef USE_GRAD_LAPLACIAN_CORRECTION
   return grad_w_Bspline(pX, qX, p->SML(), p->inv_grad_corr_M);
#else
   return grad_w_Bspline(pX, qX, p->SML());
#endif
}

std::array<double, 3> grad_w_Bspline(const networkPoint *p, const networkPoint *q) {
   return grad_w_Bspline_helper(p, p->X, q->X);
}

std::array<double, 3> grad_w_Bspline(const networkPoint *p, const Tddd &qX) {
   return grad_w_Bspline_helper(p, p->X, qX);
}

std::array<double, 3> grad_w_Bspline(const networkPoint *p, Tddd &pX, const networkPoint *q) {
   return grad_w_Bspline_helper(p, pX, q->X);
}

// # -------------------------------------------------------------------------- */

double Dot_grad_w_Bspline_helper(const networkPoint *p, const Tddd &pX, const Tddd &qX) {
#ifdef USE_GRAD_LAPLACIAN_CORRECTION
   return Dot_grad_w_Bspline(pX, qX, p->SML(), p->laplacian_corr_M);
#else
   return Dot_grad_w_Bspline(pX, qX, p->SML());
#endif
   // return Dot_grad_w_Bspline(pX, qX, p->SML());
}

double Dot_grad_w_Bspline(const networkPoint *p, const networkPoint *q) {
   return Dot_grad_w_Bspline_helper(p, p->X, q->X);
}

double Dot_grad_w_Bspline(const networkPoint *p, Tddd &X, const networkPoint *q) {
   return Dot_grad_w_Bspline_helper(p, X, q->X);
}

double Dot_grad_w_Bspline(const networkPoint *p, Tddd &X, const Tddd &qX) {
   return Dot_grad_w_Bspline_helper(p, X, qX);
}

//! -------------------------------------------------------------------------- */

std::array<double, 3> grad_w_Bspline_next_helper(const networkPoint *p, const Tddd &pX, const Tddd &qX) {
#ifdef USE_GRAD_LAPLACIAN_CORRECTION
   return grad_w_Bspline(pX, qX, p->SML_next(), p->inv_grad_corr_M_next);
#else
   return grad_w_Bspline(pX, qX, p->SML_next());
#endif
}

std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
   return grad_w_Bspline_next_helper(p, X_next(p), X_next(q));
}

std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const Tddd &X, const networkPoint *q) {
   return grad_w_Bspline_next_helper(p, X, X_next(q));
}

std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const Tddd &X) {
   return grad_w_Bspline_next_helper(p, X_next(p), X);
}

//! -------------------------------------------------------------------------- */

double Dot_grad_w_Bspline_next_helper(const networkPoint *p, const Tddd &pX, const Tddd &qX) {
#ifdef USE_GRAD_LAPLACIAN_CORRECTION
   return Dot_grad_w_Bspline(pX, qX, p->SML_next(), p->laplacian_corr_M_next);
#else
   return Dot_grad_w_Bspline(pX, qX, p->SML_next());
#endif
   // return Dot_grad_w_Bspline(pX, qX, p->SML_next());
}

double Dot_grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
   return Dot_grad_w_Bspline_next_helper(p, X_next(p), X_next(q));
}

double Dot_grad_w_Bspline_next(const networkPoint *p, Tddd &X, const networkPoint *q) {
   return Dot_grad_w_Bspline_next_helper(p, X, X_next(q));
}

double Dot_grad_w_Bspline_next(const networkPoint *p, Tddd &X, const Tddd &qX) {
   return Dot_grad_w_Bspline_next_helper(p, X, qX);
}

#endif