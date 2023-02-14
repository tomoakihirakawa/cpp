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
//! --------------------------------- ５次スプライン -------------------------------- */
double w_Bspline5(double q, const double &h) {
   constexpr double a = 2187. / (40. * M_PI);
   if ((q /= h) > 1.)
      return 0;
   else if (q < 0.333333333333333333)
      return (std::pow(1 - q, 5) - 6. * std::pow(0.6666666666666666 - q, 5) + 15. * std::pow(0.333333333333333333 - q, 5)) * a / (h * h * h);
   else if (q < 0.6666666666666666)
      return (std::pow(1 - q, 5) - 6. * std::pow(0.6666666666666666 - q, 5)) * a / (h * h * h);
   else
      return (std::pow(1 - q, 5)) * a / (h * h * h);
};
Tddd grad_w_Bspline5(const Tddd &xi, const Tddd &xj, const double h) {
   constexpr double a = 2187. / (40. * M_PI);
   double r = Norm(xi - xj);
   double q = r / h;
   if (q > 1.)
      return {0., 0., 0.};
   else if (q < 0.333333333333333333)
      return (xi - xj) * (-5 * std::pow(1. - q, 4) + 30. * std::pow(0.6666666666666666 - q, 4) - 75. * std::pow(0.333333333333333333 - q, 4)) * a / (h * h * h * r * h);
   else if (q < 0.6666666666666666)
      return (xi - xj) * (-5 * std::pow(1. - q, 4) + 30. * std::pow(0.6666666666666666 - q, 4)) * a / (h * h * h * r * h);
   else
      return (xi - xj) * (-5 * std::pow(1. - q, 4)) * a / (h * h * h * r * h);
};
//! --------------------------------- 3次スプライン -------------------------------- */
double w_Bspline3(const double &r, const double &h) {
   double q = r / h;
   if (q > 1.)
      return 0.;
   else if (q < 0.5)
      return (1. - 6 * q * q + 6 * q * q * q) / (M_PI * h * h * h);
   else
      return 2. * std::pow(1 - q, 3) / (M_PI * h * h * h);
};
Tddd grad_w_Bspline3(const Tddd &xi, const Tddd &xj, const double h) {
   // q = r / (3*h)なので，
   // grad(q) = r / (3*h) = Normalize(xi - xj) / (3 * h)
   auto r = Norm(xi - xj);
   double q = r / h;
   /*
   dqdx = dqdx(q=r/h) = (xi - xj) / (r * h)
   dWdq =
   */
   if (q > 1.)
      return {0., 0., 0.};
   else if (q < 0.5)
      return (xi - xj) * 6 * (-2 * q + 3 * q * q) / (M_PI * h * h * h * r * h);
   else
      return -(xi - xj) * 6. * std::pow(1 - q, 2) / (M_PI * h * h * h * r * h);
};
/* -------------------------------------------------------------------------- */
auto &w_Bspline = w_Bspline5;
auto &grad_w_Bspline = grad_w_Bspline5;

// auto &w_Bspline = w_Bspline3;
// auto &grad_w_Bspline = grad_w_Bspline3;

#endif