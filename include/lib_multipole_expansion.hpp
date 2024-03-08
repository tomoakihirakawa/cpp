#ifndef lib_multipole_expansion_H
#define lib_multipole_expansion_H

#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include "basic_arithmetic_array_operations.hpp"

double G(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) { return 1 / Norm(X_near - X_far); }

std::array<double, 3> gradG(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) {
   return -(X_near - X_far) / std::pow(Norm(X_near - X_far), 3.);
}

/*DOC_EXTRACT Multipole_Expansion

# 多重極展開

## Green関数の多重極展開

次のGreen関数を考える．

```math
G({\bf x},{\bf a}) = \frac{1}{\|{\bf x}-{\bf a}\|},
\quad \nabla G({\bf x},{\bf a}) = -\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}
```

グリーン関数は，球面調和関数を使って近似できる．
近似を$`G_{\rm apx}({\bf x},{\bf a},{\bf c})`$とする．

```math
G_{\rm apx}(n, {\bf x},{\bf a},{\bf c}) = \sum_{k=0}^n \sum_{m=-k}^k \left( \frac{r_{\rm near}}{r_{\rm far}} \right)^k \frac{1}{r_{\rm far}} Y(k, -m, a_{\rm near}, b_{\rm near}) Y(k, m, a_{\rm far}, b_{\rm far})=
{\bf Y}^\ast({\bf x},{\bf c})\cdot{\bf Y}({\bf a},{\bf c})
```

```math
{\bf Y}^\ast({\bf x},{\bf c}) = r_{\rm near}^k Y(k, -m, a_{\rm near},b_{\rm near}), \quad {\bf Y}({\bf a},{\bf c}) = r_{\rm far}^{-k-1} Y(k, m, a_{\rm far}, b_{\rm far})
```

ここで，$`(r_{\rm near},a_{\rm near},b_{\rm near})`$は，球面座標系に$`{\bf x}-{\bf c}`$を変換したものであり，
$`(r_{\rm far},a_{\rm far},b_{\rm far})`$は，球面座標系に$`{\bf a}-{\bf c}`$を変換したもの．$`Y(k, m, a, b)`$は球面調和関数：

```math
Y(k, m, a, b) = \sqrt{\frac{(k - |m|)!}{(k + |m|)!}} P_k^{|m|}(\cos(a)) e^{i mb}
```

$`P_k^m(x)`$はルジャンドル陪関数：

```math
P_k^m(x) = \frac{(-1)^m}{2^k k!} (1-x^2)^{m/2} \frac{d^{k+m}}{dx^{k+m}}(x^2-1)^k
```

### 球面座標系への変換

$`{\bf x}=(x,y,z)`$から球面座標$`(r,a,b)`$への変換は次のように行う．

```math
r = \|{\bf x}\|, \quad a = \arctan \frac{\sqrt{x^2 + y^2}}{z}, \quad b = \arctan \frac{y}{x}
```

$`r_\parallel=\sqrt{x^2+y^2}`$とする．$`\frac{\partial}{\partial t}(\arctan(f(t))) = \frac{f'(t)}{1 + f(t)^2}`$なので，
$`(r,a,b)`$の$`(x,y,z)`$に関する勾配は次のようになる．

```math
\nabla r = \frac{\bf x}{r},\quad
\nabla a = \frac{1}{r^2r_\parallel} \left(xz,yz,-r_\parallel^2\right),\quad
\nabla b = \frac{1}{r_\parallel^2} \left(-y,x,0\right)
```

*/

// \label{ToSphericalCoordinates}
std::array<double, 3> ToSphericalCoordinates(const std::array<double, 3>& X) {
   return {Norm(X),
           std::atan2(std::hypot(std::get<0>(X), std::get<1>(X)), std::get<2>(X)),
           std::atan2(std::get<1>(X), std::get<0>(X))};
};

std::array<std::array<double, 3>, 3> gradSphericalCoordinates(const std::array<double, 3>& X) {
   auto [x, y, z] = X;
   auto r = Norm(X);
   // auto R = std::sqrt(x * x + y * y);
   // use hypothesis
   auto R = std::hypot(x, y);
   return {std::array<double, 3>{x / r, y / r, z / r},
           std::array<double, 3>{x * z / (r * r * R), y * z / (r * r * R), -R / (r * r)},
           std::array<double, 3>{-y / (R * R), x / (R * R), 0.}};
};

// Compute the factorial of a_near given number
constexpr int factorial(const int n) {
   if (n < 0)
      throw std::runtime_error("factorial of negative number");
   if (n == 0 || n == 1)
      return 1;
   return n * factorial(n - 1);
}

double P(const double k, const double m, const double x) { return std::pow(-1., m) * std::assoc_legendre(k, m, x); };

// Compute the spherical harmonic function sph(k, m, a_far, b_far)
constexpr std::complex<double> Y(const int k, const int m, const double a_far, const double b_far) {
   if (k < 0 || std::abs(m) > k) return std::complex<double>(0.0, 0.0);
   const double assocLegendre = std::sqrt(static_cast<double>(factorial(k - std::abs(m))) / static_cast<double>(factorial(k + std::abs(m)))) * P(k, std::abs(m), std::cos(a_far));
   return std::polar(assocLegendre, m * b_far);
}

double Gapx(unsigned p,
            std::array<double, 3> X_near,
            std::array<double, 3> X_far,
            const std::array<double, 3>& center) {
   X_near -= center;
   X_far -= center;
   auto [r_near, a_near, b_near] = ToSphericalCoordinates(X_near);
   auto [r_far, a_far, b_far] = ToSphericalCoordinates(X_far);
   double inv_r_far = 1. / r_far;
   double r_near_inv_r = r_near / r_far;
   double c;
   std::complex<double> accum = 0;
   int k, m;
   for (k = 0; k <= p; ++k) {
      c = std::pow(r_near_inv_r, k) * inv_r_far;
      for (m = -k; m <= k; ++m)
         accum += c * Y(k, -m, a_near, b_near) * Y(k, m, a_far, b_far);
   }
   return accum.real();
}

// std::array<std::complex<double>, 2> Gapx_separated(unsigned p,
//                                                    std::array<double, 3> X_near,
//                                                    std::array<double, 3> X_far,
//                                                    const std::array<double, 3>& center) {
//    X_near -= center;
//    X_far -= center;
//    auto [r_near, a_near, b_near] = ToSphericalCoordinates(X_near);
//    auto [r_far, a_far, b_far] = ToSphericalCoordinates(X_far);
//    double inv_r_far = 1. / r_far;
//    double r_near_inv_r = r_near / r_far;
//    double c;
//    std::complex<double> accum_near = 0, accum_far = 0;
//    int k, m;
//    for (k = 0; k <= p; ++k) {
//       auto c_far = std::pow(1 / r_far, k + 1);
//       for (m = -k; m <= k; ++m) {
//          accum_near += r_near * Y(k, -m, a_near, b_near);
//          accum_far += c_far * Y(k, m, a_far, b_far);
//       }
//    }
//    return {accum_near, accum_far};
// }

// template <std::size_t max_n>
// struct multipole_expansion {
//    std::array<double, 3> center;
//    std::vector<std::array<std::array<double, 3>, 3>> vertex_position;
//    std::vector<std::array<double, 3>> vertex_value;

//    std::array<std::complex<double>, (max_n + 1) * (max_n + 1)> M;

//    multipole_expansion(std::vector<std::array<std::array<double, 3>, 3>> vertex_position,
//                        std::vector<std::array<double, 3>> vertex_value,
//                        const std::array<double, 3>& center)
//        : center(center), vec_X_near(vec_X_near) {
//       M.fill(0);

//       for (const auto& X012 : vec_X_near) {
//          auto G = [](const std::array<double, 3>& X, const std::array<double, 3>& a) {
//             return 1. / (nr = Norm(R = (Dot(N012_geometry, X012) - a)));
//          };
//          for (const auto& [t0, t1, ww, N012_geometry] : t0_t1_ww_N012_LOWRESOLUTION) {
//             ig = ww * ((1. - t0)) * G;
//             ign = ig * R / (nr * nr);
//             for (auto i = 0; i < 3; i++) {
//                FusedMultiplyIncrement(ig, N012_geometry[i], std::get<2>(ret[i]));
//                FusedMultiplyIncrement(-ign, N012_geometry[i], std::get<3>(ret[i]));
//             }
//             value_for_rigid_mode += ign * Total(N012_geometry);
//          }
//       }
//    }

//    void addPoles(const std::array<double, 3>& X_near) {
//       auto [r_near, a_near, b_near] = ToSphericalCoordinates(X_near - this->center);
//       int m;
//       double c;
//       for (auto n = 0; n <= this->max_n; ++n) {
//          c = std::pow(r_near, n);
//          for (m = -n; m <= n; ++m)
//             M[(n * n + 1 /*一つ小さいnのサイズ*/) + n + m] += c * Y(n, -m, a_near, b_near);
//       }
//    }

//    const std::complex<double>& getM(const int n, const int m) { return M[(n * n + 1 /*一つ小さいnのサイズ*/) + n + m]; };
// };

double dPdx(const double k, const double m, const double x) {
   double tmp = 1. / std::sqrt(1 - x * x);
   return std::pow(-1., m) * tmp * (m * x * tmp * std::assoc_legendre(k, m, x) + std::assoc_legendre(k, m + 1, x));
}

constexpr std::complex<double> dYda(const int k, const int m, const double a_far, const double b_far) {
   if (k < 0 || std::abs(m) > k) {
      return std::complex<double>(0.0, 0.0);
   }
   const double assocLegendre = std::sqrt(static_cast<double>(factorial(k - std::abs(m))) / static_cast<double>(factorial(k + std::abs(m)))) * dPdx(k, std::abs(m), std::cos(a_far)) * (-std::sin(a_far));
   return std::polar(assocLegendre, m * b_far);
}

constexpr std::complex<double> dYdb(const int k, const int m, const double a_far, const double b_far) {
   if (k < 0 || std::abs(m) > k) {
      return std::complex<double>(0.0, 0.0);
   }
   const double assocLegendre = std::sqrt(static_cast<double>(factorial(k - std::abs(m))) / static_cast<double>(factorial(k + std::abs(m)))) * P(k, std::abs(m), std::cos(a_far));
   return std::polar(assocLegendre * m, m * b_far) * std::complex<double>(0.0, 1.0);
}

std::array<double, 3> gradGapx(unsigned p,
                               const std::array<double, 3>& X_near_IN,
                               const std::array<double, 3>& X_far_IN,
                               const std::array<double, 3>& center) {
   auto [r_near, a_near, b_near] = ToSphericalCoordinates(X_near_IN - center);
   auto [r_far, a_far, b_far] = ToSphericalCoordinates(X_far_IN - center);
   auto inv_r_far = 1. / r_far;
   auto r_near_inv_r = r_near * inv_r_far;
   std::array<double, 3> grad_sphe_Gapx = {0., 0., 0.};
   std::complex<double> Y_far;
   for (int k = 0; k <= p; ++k)
      for (int m = -k; m <= k; ++m) {
         Y_far = Y(k, m, a_far, b_far);
         std::get<0>(grad_sphe_Gapx) += (k * std::pow(r_near_inv_r, k - 1) * std::pow(inv_r_far, 2.) * Y(k, -m, a_near, b_near) * Y_far).real();
         std::get<1>(grad_sphe_Gapx) += (std::pow(r_near_inv_r, k) * inv_r_far * dYda(k, -m, a_near, b_near) * Y_far).real();
         std::get<2>(grad_sphe_Gapx) += (std::pow(r_near_inv_r, k) * inv_r_far * dYdb(k, -m, a_near, b_near) * Y_far).real();
      }
   return Dot(grad_sphe_Gapx, gradSphericalCoordinates(X_far_IN - center));
}

// struct poleSummation {

//    std::complex<double> G_near = 0;
//    std::complex<double> G_far = 0;

//    poleSummation(){};

//    nearPoles(const std::vector<std::array<double, 3>>& X_near, const std::array<double, 3>& center) {
//       for (const auto& X : X_near)
//          G_near += Gapx_separated(10, X, Tddd{0., 0., 0.}, center);
//    };

//    farPoles(const std::vector<std::array<double, 3>>& X_far, const std::array<double, 3>& center) {
//       for (const auto& X : X_far)
//          G_near += Gapx_separated(10, Tddd{0., 0., 0.}, X, center);
//    };
// };

#endif