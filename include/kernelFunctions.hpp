#ifndef kernelFunctions_H
#define kernelFunctions_H

#include "basic_vectors.hpp"

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
// \label{SPH:w_Bspline5}
double w_Bspline5(double q, const double &h) {
   constexpr double a = 2187. / (40. * std::numbers::pi);
   constexpr double one_third = 1.0 / 3.0;
   constexpr double two_thirds = 2.0 / 3.0;

   if ((q /= h) > 1.)
      return 0;
   else if (q < one_third)
      return (std::pow(1 - q, 5) - 6. * std::pow(two_thirds - q, 5) + 15. * std::pow(one_third - q, 5)) * a / (h * h * h);
   else if (q < two_thirds)
      return (std::pow(1 - q, 5) - 6. * std::pow(two_thirds - q, 5)) * a / (h * h * h);
   else
      return (std::pow(1 - q, 5)) * a / (h * h * h);
};

Tddd grad_w_Bspline5(const Tddd &xi, const Tddd &xj, const double h) {
   constexpr double a = 2187. / (40. * std::numbers::pi);
   constexpr double one_third = 1.0 / 3.0;
   constexpr double two_thirds = 2.0 / 3.0;

   const double r = Norm(xi - xj);
   const double q = r / h;

   if (q > 1.)
      return {0., 0., 0.};
   else if (q < one_third) {
      const double dinom = (h * h * h * r * h);
      if (dinom == 0.0)
         return {0., 0., 0.};
      else
         return (xi - xj) * (-5 * std::pow(1. - q, 4) + 30. * std::pow(two_thirds - q, 4) - 75. * std::pow(one_third - q, 4)) * a / dinom;
   } else if (q < two_thirds)
      return (xi - xj) * (-5 * std::pow(1. - q, 4) + 30. * std::pow(two_thirds - q, 4)) * a / (h * h * h * r * h);
   else
      return (xi - xj) * (-5 * std::pow(1. - q, 4)) * a / (h * h * h * r * h);
};
double Dot_grad_w_Bspline5_Dot(const Tddd &xi, const Tddd &xj, const double h) {
   const Tddd Xij = xi - xj;
   const double r = Norm(Xij);
   const double q = r / h;
   if (q > 1. || r < 1E-14)
      return 0.;
   else
      return Dot(Xij / (r * r), grad_w_Bspline5(xi, xj, h));
};

//! --------------------------------- 3次スプライン -------------------------------- */

// \label{SPH:w_Bspline3}
#include <numbers>
double w_Bspline3(const double &r, const double &h) {
   const double q = r / h;
   if (q > 1.)
      return 0.;
   else if (q < 0.5)
      return 8. * (1. - 6. * q * q + 6. * q * q * q) / (std::numbers::pi * h * h * h);
   else
      return 8. * 2. * std::pow(1. - q, 3) / (std::numbers::pi * h * h * h);
};

Tddd grad_w_Bspline3(const Tddd &xi, const Tddd &xj, const double h) {
   const auto r = Norm(xi - xj);
   const double q = r / h;
   const Tddd dqdr = (xi - xj) / (r * h);
   const double dinom = std::numbers::pi * h * h * h * h * r;
   if (q > 1. || dinom < 1E-14)
      return {0., 0., 0.};
   else if (q < 0.5)
      return (xi - xj) * (-96. + 144. * q) * q / dinom;
   else
      return -48. * (xi - xj) * std::pow(1. - q, 2) / dinom;
};

double Dot_grad_w_Bspline3_Dot(const Tddd &xi, const Tddd &xj, const double h) {
   const Tddd Xij = xi - xj;
   const double r = Norm(Xij);
   const double q = r / h;
   if (q > 1. || r < 1E-14)
      return 0.;
   else
      return Dot(Xij / (r * r), grad_w_Bspline3(xi, xj, h));
};

// auto &w_Bspline = w_Bspline5;
// auto &grad_w_Bspline = grad_w_Bspline5;
// auto &Dot_grad_w_Bspline_Dot = Dot_grad_w_Bspline5_Dot;

auto &w_Bspline = w_Bspline3;
auto &grad_w_Bspline = grad_w_Bspline3;
auto &Dot_grad_w_Bspline_Dot = Dot_grad_w_Bspline3_Dot;

#endif