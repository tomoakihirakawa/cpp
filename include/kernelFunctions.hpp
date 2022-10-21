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
   // h<=1だけあたいを持つように変更
   h /= 3.;
   double q = r / h;
   if (q <= 3.) {
      if (2. < q)
         return std::pow(3. - q, 5) / (M_PI * h * h * h * 120.);
      else if (1. <= q)
         return (std::pow(3. - q, 5) - 6 * std::pow(2. - q, 5)) / (M_PI * h * h * h * 120.);
      else
         return (std::pow(3. - q, 5) - 6 * std::pow(2. - q, 5) + 15 * std::pow(1. - q, 5)) / (M_PI * h * h * h * 120.);
   } else
      return 0.;
};
double D_kernel_Bspline5(const double r, double h) {
   // h<=1だけあたいを持つように変更
   h /= 3;
   double q = r / h;
   if (q <= 3.) {
      if (2. < q)
         return ((-5) * std::pow(3. - q, 4)) / (M_PI * h * h * h * 120.);
      else if (1. <= q)
         return ((-5) * std::pow(3. - q, 4) + 30 * std::pow(2. - q, 4)) / (M_PI * h * h * h * 120.);
      else
         return ((-5) * std::pow(3. - q, 4) + 30 * std::pow(2. - q, 4) - 75 * std::pow(1 - q, 4)) / (M_PI * h * h * h * 120.);
   } else
      return 0.;
};
/* ------------------------------------------------------ */
double kernel_Bspline3(const Tddd &xi, const Tddd &xj, const double h) { return kernel_Bspline3(Norm(xi - xj), h); };
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

#endif