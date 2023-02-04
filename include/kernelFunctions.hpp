#ifndef kernelFunctions_H
#define kernelFunctions_H

#include "basic.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

//* ------------------------------------------------------ */
//*                  Kernel functions 核関数                */
//* ------------------------------------------------------ */

//! ------------------ Multiquadric 多重二乗 ----------------- */
double kernel_MQ(const V_d &x, const V_d &a, const double e) { return sqrt(pow(e * Norm(x - a), 2.) + 1.); };
V_d grad_kernel_MQ(const V_d &x, const V_d &a, const double e) { return (x - a) * (e * e / sqrt(pow(e * Norm(x - a), 2.) + 1.)); };
double laplacian_kernel_MQ(const V_d &x, const V_d &a /*補間点と考える*/, const double e) {
   double exyz = pow(e * Norm(x - a), 2.);
   return e * e * (3. + 2. * exyz) / pow(1. + exyz, 1.5);
};
/* ------------------------------------------------------ */
double kernel_MQ(const Tdd &x, const Tdd &a, const double e) { return sqrt(pow(e * Norm(x - a), 2.) + 1.); };
Tdd grad_kernel_MQ(const Tdd &x, const Tdd &a, const double e) { return (x - a) * (e * e / sqrt(pow(e * Norm(x - a), 2.) + 1.)); };
double laplacian_kernel_MQ(const Tdd &x, const Tdd &a, const double e) {
   double exyz = pow(e * Norm(x - a), 2.);
   return e * e * (3. + 2. * exyz) / pow(1. + exyz, 1.5);
};
double kernel_MQ(const Tddd &x, const Tddd &a, const double e) { return sqrt(pow(e * Norm(x - a), 2.) + 1.); };
Tddd grad_kernel_MQ(const Tddd &x, const Tddd &a, const double e) { return (x - a) * (e * e / sqrt(pow(e * Norm(x - a), 2.) + 1.)); };
double laplacian_kernel_MQ(const Tddd &x, const Tddd &a, const double e) {
   double exyz = pow(e * Norm(x - a), 2.);
   return e * e * (3. + 2. * exyz) / pow(1. + exyz, 1.5);
};
/* ------------------------------------------------------ */
double kernel_TPS(const Tdd &x, const Tdd &a, const double e) {
   double r = Norm(x - a);
   if (r < 1E-15)
      return 0.;
   return pow(r, 2.) * log(r * e);
};
Tdd grad_kernel_TPS(const Tdd &x, const Tdd &a, const double e) {
   double r = Norm(x - a);
   if (r < 1E-15)
      return {0., 0.};
   return -(a - x) * (1. + 2. * log(e * r));
};
double laplacian_kernel_TPS(const Tdd &x, const Tdd &a, const double e) {
   double r = Norm(x - a);
   if (r < 1E-15)
      return 0.;
   return 5. + 6. * log(r * e);
};
double kernel_TPS(const Tddd &x, const Tddd &a, const double e) {
   double r = Norm(x - a);
   if (r < 1E-15)
      return 0.;
   return pow(r, 2.) * log(r * e);
};
Tddd grad_kernel_TPS(const Tddd &x, const Tddd &a, const double e) {
   double r = Norm(x - a);
   if (r < 1E-15)
      return {0., 0., 0.};
   return -(a - x) * (1. + 2. * log(e * r));
};
double laplacian_kernel_TPS(const Tddd &x, const Tddd &a, const double e) {
   double r = Norm(x - a);
   if (r < 1E-15)
      return 0.;
   return pow(r, 2.) * log(r * e);
};
//! --------------------------------------------------------- */
//! ------------------ 3rd order Bspline functions ----------------- */
// smlは2hとしている
V_d grad_q_by_first_kernel_Bspline3(const V_d &x, const V_d &a, double e) {
   //! grad by x (the first argument)
   return (x - a) / Norm(x - a) / e;
};
V_d grad_q_by_second_kernel_Bspline3(const V_d &x, const V_d &a, double e) {
   //! grad by x (the first argument)
   return -(x - a) / Norm(x - a) / e;
};
double kernel_Bspline3(const double r, double h) {
   h /= 2.;
   double q = r / h;
   if (q <= 2.) {
      if (1 < q)
         return (1. / (M_PI * h * h * h)) * std::pow(2. - q, 3.) / 4.;
      else
         return (1. / (M_PI * h * h * h)) * (1. - 3. / 2. * q * q + 3. / 4. * q * q * q);
   } else
      return 0.;
};

double D_kernel_Bspline3(const double r, double h) {
   h /= 2.;
   double q = r / h;
   if (q <= 2.) {
      if (1 < q)
         return (1. / (M_PI * h * h * h)) * std::pow(2. - q, 2.) * (-3. / 4.);
      else
         return (1. / (M_PI * h * h * h)) * (-3. * q + 9. / 4. * q * q);
   } else
      return 0.;
};
double D2_kernel_Bspline3(const double r, double h) {
   h /= 2.;
   double q = r / h;
   if (q <= 2.) {
      if (1 < q)
         return (1. / (M_PI * h * h * h)) * (3 * (2 - q)) / 2.;
      else
         return (1. / (M_PI * h * h * h)) * (-3 + (9 * q) / 2.);
   } else
      return 0.;
};
//! ------------------ 5th order Bspline functions ----------------- */
double kernel_Bspline5(const double r, double h) {
   h /= 3.;
   double q = r / h;
   if (q <= 3.) {
      double alpha = 1 / (M_PI * h * h * h * 120.);
      if (2. < q)
         return alpha * std::pow(3. - q, 5);
      else if (1. <= q)
         return alpha * (std::pow(3. - q, 5) - 6 * std::pow(2. - q, 5));
      else
         return alpha * (std::pow(3. - q, 5) - 6 * std::pow(2. - q, 5) + 15 * std::pow(1. - q, 5));
   } else
      return 0.;
};
double D_kernel_Bspline5(const double r, double h) {
   h /= 3;
   double q = r / h;
   if (q <= 3.) {
      double alpha = 1 / (M_PI * h * h * h * 120.);
      if (2. < q)
         return -5. * alpha * std::pow(3. - q, 4);
      else if (1. <= q)
         return -5. * alpha * (std::pow(3. - q, 4) - 6 * std::pow(2. - q, 4));
      else
         return -5. * alpha * (std::pow(3. - q, 4) - 6 * std::pow(2. - q, 4) + 15 * std::pow(1. - q, 4));
   } else
      return 0.;
};

/* ------------------------------------------------------ */
double kernel_Bspline3(const Tddd &xi, const Tddd &xj, const double h) {
   return kernel_Bspline3(Norm(xi - xj), h);
};
Tddd grad_kernel_Bspline3(const Tddd &xi, const Tddd &xj, const double h) {
   double r = Norm(xi - xj);
   return -(xi - xj) / (r * h / 2.) * D_kernel_Bspline3(r, h);
};
//
double kernel_Bspline5(const Tddd &xi, const Tddd &xj, const double h) { return kernel_Bspline5(Norm(xi - xj), h); };
Tddd grad_kernel_Bspline5(const Tddd &xi, const Tddd &xj, const double h) {
   double r = Norm(xi - xj);
   return -(xi - xj) / (r * h / 3.) * D_kernel_Bspline5(r, h);
};
//! --------------------------------- ５次スプライン -------------------------------- */
double w_Bspline5(double s, const double &h) {
   if ((s /= h) < 1.) {
      if (s < 0.333333333333333333)
         return (std::pow(1 - s, 5) - 6. * std::pow(0.6666666666666666 - s, 5) + 15. * std::pow(0.333333333333333333 - s, 5)) * 2187. / (40. * M_PI * h * h * h);
      else if (s < 0.6666666666666666)
         return (std::pow(1 - s, 5) - 6. * std::pow(0.6666666666666666 - s, 5)) * 2187. / (40. * M_PI * h * h * h);
      else
         return (std::pow(1 - s, 5)) * 2187. / (40. * M_PI * h * h * h);
   } else
      return 0.;
};
Tddd grad_w_Bspline5(const Tddd &xi, const Tddd &xj, const double h) {
   double r = Norm(xi - xj);
   double s = r / h;
   if (s < 1.) {
      if (s < 0.333333333333333333)
         return (xi - xj) / (r * h) * (-5 * std::pow(1. - s, 4) + 30. * std::pow(0.6666666666666666 - s, 4) - 75. * std::pow(0.333333333333333333 - s, 4)) * 2187. / (40. * M_PI * h * h * h);
      else if (s < 0.6666666666666666)
         return (xi - xj) / (r * h) * (-5 * std::pow(1. - s, 4) + 30. * std::pow(0.6666666666666666 - s, 4)) * 2187. / (40. * M_PI * h * h * h);
      else
         return (xi - xj) / (r * h) * (-5 * std::pow(1. - s, 4)) * 2187. / (40. * M_PI * h * h * h);
   } else
      return {0., 0., 0.};
};
//! --------------------------------- 3次スプライン -------------------------------- */
double w_Bspline3(const double &r, double h) {
   h /= 2;
   double q = r / h;
   if (q < 2.) {
      if (q < 1.)
         return (1. - 3. / 2. * q * q * (1. - q / 2.)) / (M_PI * h * h * h);
      else
         return std::pow(2. - q, 3.) / (M_PI * h * h * h) / 4.;
   } else
      return 0.;
};
Tddd grad_w_Bspline3(const Tddd &xi, const Tddd &xj, double h) {
   h /= 2;
   double r = Norm(xi - xj);
   double q = r / h;
   Tddd dqdx = (xi - xj) / (r * h);
   if (q < 2.) {
      if (q < 1.)
         return (-3. / 2. * 2 * q * dqdx * (1. - q / 2.) - 3. / 2. * q * q * (-dqdx / 2.)) / (M_PI * h * h * h);
      else
         return -3 * std::pow(2. - q, 2.) * dqdx / (M_PI * h * h * h) / 4.;
   } else
      return {0., 0., 0.};
};
/* -------------------------------------------------------------------------- */
//
auto &w_Bspline = w_Bspline5;
auto &grad_w_Bspline = grad_w_Bspline5;
// auto &w_Bspline = w_Bspline3;
// auto &grad_w_Bspline = grad_w_Bspline3;

#endif