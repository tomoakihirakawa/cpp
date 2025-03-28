#ifndef SPH_KERNEL_HELPER_FUNCTIONS_HPP
#define SPH_KERNEL_HELPER_FUNCTIONS_HPP

// # -------------------------------------------------------------------------- */

std::array<double, 3> grad_w_Bspline(const networkPoint* p, const networkPoint* q) {
#ifdef USE_GRAD_CORRECTION
   return Dot(p->inv_grad_corr_M, grad_w_Bspline(p->X, q->X, p->SML_grad()));
#else
   return grad_w_Bspline(p->X, q->X, p->SML_grad());
#endif
}

// # -------------------------------------------------------------------------- */

double Dot_grad_w_Bspline(const networkPoint* p, const networkPoint* q) {
#ifdef USE_LAPLACIAN_CORRECTION
   return Dot_grad_w_Bspline(p->X, q->X, p->SML_grad(), p->laplacian_corr_M);
#else
   return Dot_grad_w_Bspline(p->X, q->X, p->SML_grad());
   // const std::array<double, 3> Xij = p->X - q->X;
      // const double r = Norm(Xij);
      // return (r / p->SML() > 1. || r < 1E-13) ? 0. : Dot(Xij / (r * r), grad_w_Bspline(p, q));
#endif
}

//! -------------------------------------------------------------------------- */

std::array<double, 3> grad_w_Bspline_next(const networkPoint* p, const networkPoint* q) {
#ifdef USE_GRAD_CORRECTION
   return Dot(p->inv_grad_corr_M_next, grad_w_Bspline(X_next(p), X_next(q), p->SML_grad_next()));
#else
   return grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
#endif
}

std::array<double, 3> grad_w_Bspline_next(const networkPoint* p, const Tddd& X, const networkPoint* q) {
#ifdef USE_GRAD_CORRECTION
   return Dot(p->inv_grad_corr_M_next, grad_w_Bspline(X, X_next(q), p->SML_grad_next()));
#else
   return grad_w_Bspline(X, X_next(q), p->SML_next());
#endif
}

std::array<double, 3> grad_w_Bspline_for_dot_grad(const networkPoint* p, const networkPoint* q) {
#ifdef USE_GRAD_CORRECTION
   return Dot(p->inv_grad_corr_M_next, grad_w_Bspline_for_dot_grad(X_next(p), X_next(q), p->SML_grad_next()));
#else
   return grad_w_Bspline_for_dot_grad(X_next(p), X_next(q), p->SML_next());
#endif
}


std::array<double, 3> grad_w_Bspline_for_dot_grad(const networkPoint* p, const Tddd& X, const networkPoint* q) {
#ifdef USE_GRAD_CORRECTION
   return Dot(p->inv_grad_corr_M_next, grad_w_Bspline_for_dot_grad(X, X_next(q), p->SML_grad_next()));
#else
   return grad_w_Bspline_for_dot_grad(X, X_next(q), p->SML_next());
#endif
}

//! -------------------------------------------------------------------------- */

double Dot_grad_w_Bspline_next(const networkPoint* p, const networkPoint* q) {
#ifdef USE_LAPLACIAN_CORRECTION
   return Dot_grad_w_Bspline(X_next(p), X_next(q), p->SML_grad_next(), p->laplacian_corr_M_next);
#else
   return Dot_grad_w_Bspline(X_next(p), X_next(q), p->SML_grad_next());
   // const std::array<double, 3> Xij = X_next(p) - X_next(q);
      // const double r = Norm(Xij);
      // return (r / p->SML_grad_next() > 1. || r < 1E-13) ? 0. : Dot(Xij / (r * r), grad_w_Bspline_next(p, q));
#endif
}

// std::array<double, 3> grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
// #ifdef USE_GRAD_CORRECTION
//    return Dot(p->inv_grad_corr_M_next, grad_w_Bspline(X_next(p), X_next(q), p->SML_next()));
// #else
//    return grad_w_Bspline(X_next(p), X_next(q), p->SML_next());
// #endif
// }

// # -------------------------------------------------------------------------- */

// double Dot_grad_w_Bspline_next(const networkPoint *p, const networkPoint *q) {
//    //! これはラプラシアンの計算に使われる
//    const std::array<double, 3> Xij = X_next(p) - X_next(q);
//    const double h = p->SML_grad_next();
//    const double r = Norm(Xij);
//    if (r / h > 1. || r < 1E-13)
//       return 0.;
//    else
//       // return Total(Total(p->inv_grad_corr_M_next * TensorProduct(Xij / (r * r), grad_w_Bspline(X_next(p), X_next(q), h))));
//       return Dot(Xij / (r * r), grad_w_Bspline(X_next(p), X_next(q), h));
//    // return Dot(Xij / (r * r), grad_w_Bspline_next(p, q));
// }

double Dot_grad_w_Bspline_next(const networkPoint* p, const Tddd& X, const networkPoint* q) {
#ifdef USE_LAPLACIAN_CORRECTION
   return Dot_grad_w_Bspline(X, X_next(q), p->SML_grad_next(), p->laplacian_corr_M_next);
#else
   return Dot_grad_w_Bspline(X, X_next(q), p->SML_grad_next());
   // const std::array<double, 3> Xij = X_next(p) - X_next(q);
      // const double r = Norm(Xij);
      // return (r / p->SML_grad_next() > 1. || r < 1E-13) ? 0. : Dot(Xij / (r * r), grad_w_Bspline_next(p, X, q));
#endif
}

#endif