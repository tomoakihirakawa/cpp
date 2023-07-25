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
/*DOC_EXTRACT SPH:kernelFunctions

## 5次スプライン関数

３次元のシミュレーションで用いられる，5次スプライン関数：

```math
W_{5}(q,h) = \frac{2187}{40\pi h^3}
\begin{cases}
(1-q)^5 - 6(\frac{2}{3}-q)^5 + 15(\frac{1}{3}-q)^5 & (0 \leq q < \frac{1}{3}) \\
(1-q)^5 - 6(\frac{2}{3}-q)^5 & (\frac{1}{3} \leq q < \frac{2}{3}) \\
(1-q)^5 & (\frac{2}{3} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
```

その勾配は，

```math
\nabla W_{5}(q,h) = \frac{2187}{40\pi h^4}
\begin{cases}
-5(1-q)^4 + 30(\frac{2}{3}-q)^4 - 75(\frac{1}{3}-q)^4 & (0 \leq q < \frac{1}{3}) \\
-5(1-q)^4 + 30(\frac{2}{3}-q)^4 & (\frac{1}{3} \leq q < \frac{2}{3}) \\
-5(1-q)^4 & (\frac{2}{3} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
```

ヘッシアンは，

```math
\nabla^2 W_{5}(q,h) = \frac{2187}{40\pi h^5}
\begin{cases}
20(1-q)^3 - 120(\frac{2}{3}-q)^3 + 300(\frac{1}{3}-q)^3 & (0 \leq q < \frac{1}{3}) \\
20(1-q)^3 - 120(\frac{2}{3}-q)^3 & (\frac{1}{3} \leq q < \frac{2}{3}) \\
20(1-q)^3 & (\frac{2}{3} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
```

*/
double w_Bspline5(double q, const double &h) {
   constexpr double a = 2187. / (40. * M_PI);
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
   constexpr double a = 2187. / (40. * M_PI);
   constexpr double one_third = 1.0 / 3.0;
   constexpr double two_thirds = 2.0 / 3.0;

   const double r = Norm(xi - xj);
   const double q = r / h;
   const Tddd grad_q = (xi - xj) / (r * h);
   const double dinom = h * h * h;
   const double c = a / dinom;

   if (q > 1. || r * h * h * h * h == 0.0)
      return {0., 0., 0.};
   else if (q < one_third) {
      if (r * h * h * h * h == 0.0)
         return {0., 0., 0.};
      else
         return grad_q * -(5 * std::pow(1. - q, 4) - 30. * std::pow(two_thirds - q, 4) + 75. * std::pow(one_third - q, 4)) * c;
   } else if (q < two_thirds)
      return grad_q * -(5 * std::pow(1. - q, 4) - 30. * std::pow(two_thirds - q, 4)) * c;
   else
      return grad_q * -(5 * std::pow(1. - q, 4)) * c;
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
/*DOC_EXTRACT SPH:kernelFunctions

## 3次スプライン関数

３次元のシミュレーションで用いられる，3次スプライン関数：

```math
W_{3}(q,h) = \frac{8}{\pi h^3}
\begin{cases}
(1-6q^2+6q^3) & (0 \leq q < \frac{1}{2}) \\
2(1-q)^3 & (\frac{1}{2} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
```

その勾配は，

```math
\nabla W_{3}(q,h) =
-{\nabla q}\frac{48}{\pi h^4}
\begin{cases}
-2q+3q^2 & (0 \leq q < \frac{1}{2}) \\
(1-q)^2 & (\frac{1}{2} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
,\quad
\nabla q = \frac{\nabla r}{h} = \frac{{\bf x}_i - {\bf x}_j}{r h}
```

ヘッシアンは，

```math
\nabla\otimes\nabla W_{3}(q,h) =
-\frac{48}{\pi h^4}
\begin{cases}
-2\nabla \otimes{(q\nabla q)}+3\nabla \otimes{(q^2\nabla q)} & (0 \leq q < \frac{1}{2}) \\
\nabla \otimes ((1-q)^2\nabla q) & (\frac{1}{2} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
,\quad
\nabla q = \frac{\nabla r}{h} = \frac{{\bf x}_i - {\bf x}_j}{r h}\\
= -\frac{48}{\pi h^4}
\begin{cases}
-2\nabla \otimes{(\frac{r}{h}\frac{{\bf x}_i - {\bf x}_j}{r})}+3\nabla \otimes{(\frac{r^2}{h^2}\frac{{\bf x}_i - {\bf x}_j}{r})} & (0 \leq q < \frac{1}{2}) \\
\nabla \otimes{(\frac{r^2}{h^2}\frac{{\bf x}_i - {\bf x}_j}{r})} & (\frac{1}{2} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
```


*/

#include <numbers>
double w_Bspline3(const double &r, const double &h) {
   const double q = r / h;
   if (q > 1.)
      return 0.;
   else if (q < 0.5)
      return 8. * (1. - 6. * q * q + 6. * q * q * q) / (M_PI * h * h * h);
   else
      return 8. * 2. * std::pow(1. - q, 3) / (M_PI * h * h * h);
};

double ddr_w_Bspline3(const double &r, const double &h) {
   const double q = r / h;
   const double dqdr = 1. / h;
   if (q > 1.)
      return 0.;
   else if (q < 0.5)
      return 8. * (-12. * q + 18. * q * q) * dqdr / (M_PI * h * h * h);
   else
      return 8. * 2. * 3 * std::pow(1. - q, 2) * (-dqdr) / (M_PI * h * h * h);
};

double ddr2_w_Bspline3(const double &r, const double &h) {
   const double q = r / h;
   const double dqdr = 1. / h;
   if (q > 1.)
      return 0.;
   else if (q < 0.5)
      return 8. * ((-12. + 36. * q) * dqdr) / (M_PI * h * h * h);
   else
      return 8. * 2. * 3 * 2 * std::pow(1. - q, 1) * (-dqdr) * (-dqdr) / (M_PI * h * h * h);
};

std::array<double, 3> grad_w_Bspline3(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
   const auto r = Norm(xi - xj);
   const double q = r / h;
   const std::array<double, 3> dqdr = (xi - xj) / (r * h);
   const double dinom = M_PI * h * h * h * h * r;
   if (q > 1. || dinom < 1E-13)
      return {0., 0., 0.};
   else if (q < 0.5)
      return (xi - xj) * (-96. + 144. * q) * q / dinom;
   else
      return -48. * (xi - xj) * std::pow(1. - q, 2) / dinom;
   //
   // const auto r = Norm(xi - xj);
   // const double q = r / h;
   // const std::array<double, 3> grad_q = (xi - xj) / (r * h);
   // const double dinom = h * h * h;
   // const double c = 48. / M_PI / dinom;
   // if (q > 1. || dinom < 1E-14)
   //    return {0., 0., 0.};
   // else if (q < 0.5)
   //    return -c * (-2 * q + 3 * q * q) * grad_q;
   // else
   //    return -c * std::pow(1. - q, 2) * grad_q;
};

// std::array<std::array<double, 3>, 3> H_grad_w_Bspline3(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
//    const auto r = Norm(xi - xj);
//    const auto grad_r = (xi - xj) / Norm(xi - xj);
//    const double q = r / h;
//    const std::array<double, 3> grad_q = (xi - xj) / (r * h);
//    std::array<std::array<double, 3>, 3> I = {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}};
//    const auto H_q = TensorProduct((xi - xj) / h, grad_r / (r * r)) + (I / h) / r;
//    const double dinom = h * h * h;
//    const double c = 8. / M_PI / dinom;
//    if (q > 1. || dinom < 1E-14)
//       return {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
//    else if (q < 0.5)
//       return c * (-6. * 2 * q + 6. * 3 * q * q) * grad_q;
//    else
//       return c * 2. * 3 * std::pow(1. - q, 2) * (-grad_q);
// };
//

double Dot_grad_w_Bspline3_Dot(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
   const std::array<double, 3> Xij = xi - xj;
   const double r = Norm(Xij);
   const double q = r / h;
   if (q > 1. || r < 1E-13)
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
auto &ddr_w_Bspline = ddr_w_Bspline3;

#endif