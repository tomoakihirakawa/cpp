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

// Compute the factorial of a_near given number
// static const std::vector<uint64_t> factorial_list = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000};
static const std::vector<double> factorial_list_as_double = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000};

constexpr double factorial(const int n) {
   if (n < 0)
      throw std::runtime_error("factorial of negative number");
   else if (n >= factorial_list_as_double.size())
      return n * factorial(n - 1);
   else
      return factorial_list_as_double[n];
}

std::array<double, 3> ToSphericalCoordinates(const double x, const double y, const double z) {
   return {std::hypot(x, y, z), std::atan2(std::hypot(x, y), z), std::atan2(y, x)};
};

// \label{ToSphericalCoordinates}
std::array<double, 3> ToSphericalCoordinates(const std::array<double, 3>& X) {
   return ToSphericalCoordinates(std::get<0>(X), std::get<1>(X), std::get<2>(X));
};

std::array<std::array<double, 3>, 3> gradSphericalCoordinates(const std::array<double, 3>& X) {
   const auto [x, y, z] = X;
   const double r = Norm(X);
   const double R = std::hypot(x, y);
   const double R2 = (x * x + y * y);
   const double x_r = x / r;
   const double y_r = y / r;
   const double z_r = z / r;
   return {std::array<double, 3>{x_r, y_r, z_r},
           std::array<double, 3>{x_r * z_r / R, y_r * z_r / R, -R / (r * r)},
           std::array<double, 3>{-y / R2, x / R2, 0.}};
};

/*

grad R = spherical_grad_R \cdot JacobianCartesian2Spherical

*/

std::array<std::array<double, 3>, 3> JacobianCartesian2Spherical(const std::array<double, 3>& X) {
   auto [x, y, z] = X;
   const double R = std::hypot(x, y);
   const double R2 = R * R;
   const double r = std::hypot(x, y, z);
   const double x_r = x / r;
   const double y_r = y / r;
   const double z_r = z / r;
   const std::array<double, 3> grad_r = std::array<double, 3>{x_r, y_r, z_r};
   const std::array<double, 3> grad_theta = std::array<double, 3>{x_r * z_r / R, y_r * z_r / R, -R / (r * r)};
   const std::array<double, 3> grad_phi = std::array<double, 3>{-y / R2, x / R2, 0.};
   return {grad_r, grad_theta, grad_phi};
};

/* -------------------------------------------------------------------------- */

template <typename T>
void complex_fused_multiply_increment(const std::complex<T>& a, const std::complex<T>& b, std::complex<T>& c) {
   c = std::complex<T>(std::fma(a.real(), b.real(), std::fma(-a.imag(), b.imag(), c.real())),
                       std::fma(a.real(), b.imag(), std::fma(a.imag(), b.real(), c.imag())));
}

std::complex<double> Dot(const std::array<double, 3>& normal, const std::array<std::complex<double>, 3>& abc) {
   return std::get<0>(abc) * std::get<0>(normal) + std::get<1>(abc) * std::get<1>(normal) + std::get<2>(abc) * std::get<2>(normal);
};

std::complex<double> Dot(const std::array<std::complex<double>, 3>& abc, const std::array<double, 3>& normal) {
   return std::get<0>(abc) * std::get<0>(normal) + std::get<1>(abc) * std::get<1>(normal) + std::get<2>(abc) * std::get<2>(normal);
};

std::array<std::complex<double>, 3> Dot(const std::array<std::complex<double>, 3>& abc, const std::array<std::array<double, 3>, 3>& A) {
   std::array<std::complex<double>, 3> ret;
   std::get<0>(ret) = abc[0] * A[0][0] + abc[1] * A[1][0] + abc[2] * A[2][0];
   std::get<1>(ret) = abc[0] * A[0][1] + abc[1] * A[1][1] + abc[2] * A[2][1];
   std::get<2>(ret) = abc[0] * A[0][2] + abc[1] * A[1][2] + abc[2] * A[2][2];
   return ret;
};

/*

\cite{Liu_2009}p.70には，体球調和関数を使った高速多重極展開のアルゴリズムが書かれている．
以下は，その定義通りの体球調和関数を計算する関数である．

*/

constexpr std::complex<double> complex_zero(0., 0.);
constexpr std::array<std::complex<double>, 3> complex_zero3D{complex_zero, complex_zero, complex_zero};

//! (-1)^mを含まないルジャンドル陪関数
double D_assoc_legendre(int n, int m, const double x) {
   if (std::abs(x) > 1.0)
      throw std::invalid_argument("Argument out of bounds for Legendre polynomials.");
   n = std::abs(n);
   m = std::abs(m);
   return (-m * x / (1. - x * x) * std::assoc_legendre(n, m, x) + std::assoc_legendre(n, m + 1, x) / std::sqrt(1. - x * x));
}

std::complex<double> SolidHarmonicR(int n, int m, const std::array<double, 3>& X) {
   n = std::abs(n);
   m = std::abs(m);
   if (n + m < 0)
      return std::complex<double>(0.0, 0.0);
   const auto [rho, theta, phi] = ToSphericalCoordinates(X);
   return std::polar((std::pow(rho, n) / factorial(n + m)) * std::assoc_legendre(n, m, std::cos(theta)), m * phi);
}

std::array<std::complex<double>, 3> Grad_SolidHarmonicR(int n, int m, const std::array<double, 3>& X) {
   n = std::abs(n);
   m = std::abs(m);
   if (n + m < 0)
      return complex_zero3D;
   const auto [rho, theta, phi] = ToSphericalCoordinates(X);
   const double rho_n = std::pow(rho, n);
   const double TMP = rho_n / factorial(n + m);
   const double tmp = TMP * std::assoc_legendre(n, m, std::cos(theta));
   const std::complex<double> dRdrho = std::polar(n * (tmp / rho), m * phi);
   const std::complex<double> dRdtheta = std::polar(TMP * D_assoc_legendre(n, m, std::cos(theta)) * (-std::sin(theta)), m * phi);
   const std::complex<double> dRdphi = std::polar(m * tmp, m * phi) * std::complex<double>(0, m);

   return Dot(std::array<std::complex<double>, 3>{dRdrho, dRdtheta, dRdphi}, JacobianCartesian2Spherical(X));
}

std::complex<double> SolidHarmonicS(int n, int m, const std::array<double, 3>& X) {
   n = std::abs(n);
   m = std::abs(m);
   if (n - m < 0)
      return std::complex<double>(0.0, 0.0);

   const auto [rho, theta, phi] = ToSphericalCoordinates(X);
   return std::polar((factorial(n - m) / std::pow(rho, n + 1.)) * std::assoc_legendre(n, m, std::cos(theta)), m * phi);
}

/* -------------------------------------------------------------------------- */

#include "basic_exception.hpp"

template <int N>
struct ExpCoeffs {
   std::array<double, 3> X;
   std::array<std::array<std::complex<double>, N + 1 /*0~N*/>, 2 * N> coeffs;
   std::array<std::array<std::complex<double>, N + 1 /*0~N*/>, 2 * N> coeffs_;

   std::complex<double> zero = {0., 0.};
   const std::complex<double>& get_coeffs(int n, int m) const {
      if (!indexedQ(n, m)) {
         // std::cout << "n = " << n << ", m = " << m << ", n+m = " << n + m << std::endl;
         // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__);
         return zero;
      }
      return coeffs[n][n + m];
   }
   const std::complex<double>& get_coeffs_(int n, int m) const {
      if (!indexedQ(n, m)) {
         // std::cout << "n = " << n << ", m = " << m << ", n+m = " << n + m << std::endl;
         // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__);
         return zero;
      }
      return coeffs_[n][n + m];
   }

   bool indexedQ(const int n, const int m) const { return n <= N && -n <= m && m <= n; }

   ExpCoeffs(const std::array<double, 3>& XIN) : X(XIN) {
      for (auto& cc : coeffs) cc.fill(0);
      for (auto& cc : coeffs_) cc.fill(0);
   }

   void increment(const Tddd& XIN, const std::array<double, 2> weights, const Tddd& normal) {
      const Tddd R = XIN - this->X;
      auto set_coeffs = [&](int n, int m) -> std::array<std::complex<double>, 2> {
         return {SolidHarmonicR(n, m, R) * weights[0],
                 Dot(normal, Grad_SolidHarmonicR(n, m, R)) * weights[1]};
      };
      this->set(set_coeffs);
   };

   void set(const auto& func) {
      std::array<std::complex<double>, 2> cc;
      for (int n = 0; n <= N; ++n) {
         cc = func(n, 0);
         coeffs[n][n] += std::get<0>(cc);
         coeffs_[n][n] += std::get<1>(cc);
         for (int m = 1; m <= n; ++m) {
            cc = func(n, m);
            coeffs[n][n + m] += std::get<0>(cc);
            coeffs_[n][n + m] += std::get<1>(cc);
            //
            cc = func(n, -m);
            coeffs[n][n - m] += std::get<0>(cc);
            coeffs_[n][n - m] += std::get<1>(cc);
         }
      }
   }

   std::complex<double> IG_using_M(const std::array<double, 3>& O) {
      std::array<double, 3> Xc2O = O - this->X;
      std::complex<double> ret = 0;
      for (int n = 0; n <= N; ++n) {
         complex_fused_multiply_increment(std::conj(SolidHarmonicS(n, 0, Xc2O)), this->get_coeffs(n, 0), ret);
         for (int m = 1; m <= n; ++m) {
            complex_fused_multiply_increment(std::conj(SolidHarmonicS(n, m, Xc2O)), this->get_coeffs(n, m), ret);
            complex_fused_multiply_increment(std::conj(SolidHarmonicS(n, -m, Xc2O)), this->get_coeffs(n, -m), ret);
         }
      }

      return ret;
   }

   std::complex<double> IG_using_L(const std::array<double, 3>& O) {
      std::array<double, 3> Xc2O = O - this->X;
      std::complex<double> ret = 0;
      for (int n = 0; n <= N; ++n) {
         complex_fused_multiply_increment(SolidHarmonicR(n, 0, Xc2O), this->get_coeffs(n, 0), ret);
         for (int m = 1; m <= n; ++m) {
            complex_fused_multiply_increment(SolidHarmonicR(n, m, Xc2O), this->get_coeffs(n, m), ret);
            complex_fused_multiply_increment(SolidHarmonicR(n, -m, Xc2O), this->get_coeffs(n, -m), ret);
         }
      }
      return ret;
   }

   std::complex<double> IGn_using_M(const std::array<double, 3>& O) {
      std::array<double, 3> Xc2O = O - this->X;
      std::complex<double> ret = 0;
      for (int n = 0; n <= N; ++n) {
         complex_fused_multiply_increment(std::conj(SolidHarmonicS(n, 0, Xc2O)), this->get_coeffs_(n, 0), ret);
         for (int m = 1; m <= n; ++m) {
            complex_fused_multiply_increment(std::conj(SolidHarmonicS(n, m, Xc2O)), this->get_coeffs_(n, m), ret);
            complex_fused_multiply_increment(std::conj(SolidHarmonicS(n, -m, Xc2O)), this->get_coeffs_(n, -m), ret);
         }
      }
      return ret;
   }

   std::complex<double> IGn_using_L(const std::array<double, 3>& O) {
      std::array<double, 3> Xc2O = O - this->X;
      std::complex<double> ret = 0;
      for (int n = 0; n <= N; ++n) {
         complex_fused_multiply_increment(SolidHarmonicR(n, 0, Xc2O), this->get_coeffs_(n, 0), ret);
         for (int m = 1; m <= n; ++m) {
            complex_fused_multiply_increment(SolidHarmonicR(n, m, Xc2O), this->get_coeffs_(n, m), ret);
            complex_fused_multiply_increment(SolidHarmonicR(n, -m, Xc2O), this->get_coeffs_(n, -m), ret);
         }
      }
      return ret;
   }
};

//! Returns the shifted ExpCoeffs of M to the new center X
template <int N>
ExpCoeffs<N> M2M(const ExpCoeffs<N>& M, const std::array<double, 3>& XIN) {
   ExpCoeffs<N> M_shifted(XIN);
   std::array<double, 3> r = M.X - XIN;
   std::complex<double> R;

   auto set_coeffs = [&](int n, int m) -> std::array<std::complex<double>, 2> {
      std::array<std::complex<double>, 2> ret = {0, 0};
      int m_;
      for (int n_ = 0; n_ <= n; ++n_) {
         m_ = 0;
         R = SolidHarmonicR(n_, m_, r);
         complex_fused_multiply_increment(R, M.get_coeffs(n - n_, m - m_), std::get<0>(ret));
         complex_fused_multiply_increment(R, M.get_coeffs_(n - n_, m - m_), std::get<1>(ret));
         for (m_ = 1; m_ <= n_; ++m_) {
            R = SolidHarmonicR(n_, m_, r);
            complex_fused_multiply_increment(R, M.get_coeffs(n - n_, m - m_), std::get<0>(ret));
            complex_fused_multiply_increment(R, M.get_coeffs_(n - n_, m - m_), std::get<1>(ret));
            //
            R = SolidHarmonicR(n_, -m_, r);
            complex_fused_multiply_increment(R, M.get_coeffs(n - n_, m + m_), std::get<0>(ret));
            complex_fused_multiply_increment(R, M.get_coeffs_(n - n_, m + m_), std::get<1>(ret));
         }
      }
      return ret;
   };

   M_shifted.set(set_coeffs);

   return M_shifted;
};

//! Returns the LocalExp of M at the new center X
template <int N>
ExpCoeffs<N> M2L(const ExpCoeffs<N>& M, const std::array<double, 3>& XIN) {
   ExpCoeffs<N> LocalExp(XIN);
   std::complex<double> S_;
   std::array<double, 3> r = XIN - M.X;

   auto set_coeffs = [&](int n, int m) -> std::array<std::complex<double>, 2> {
      std::array<std::complex<double>, 2> ret = {0, 0};
      int m_;
      double sign = (n % 2 == 0) ? 1.0 : -1.0;
      for (int n_ = 0; n_ <= N; ++n_) {
         m_ = 0;
         S_ = sign * std::conj(SolidHarmonicS(n + n_, m + m_, r));
         complex_fused_multiply_increment(S_, M.get_coeffs(n_, m_), std::get<0>(ret));
         complex_fused_multiply_increment(S_, M.get_coeffs_(n_, m_), std::get<1>(ret));
         for (m_ = 1; m_ <= n_; ++m_) {
            S_ = sign * std::conj(SolidHarmonicS(n + n_, m + m_, r));
            complex_fused_multiply_increment(S_, M.get_coeffs(n_, m_), std::get<0>(ret));
            complex_fused_multiply_increment(S_, M.get_coeffs_(n_, m_), std::get<1>(ret));
            //
            S_ = sign * std::conj(SolidHarmonicS(n + n_, m + -m_, r));
            complex_fused_multiply_increment(S_, M.get_coeffs(n_, -m_), std::get<0>(ret));
            complex_fused_multiply_increment(S_, M.get_coeffs_(n_, -m_), std::get<1>(ret));
         }
      }
      return ret;
   };

   LocalExp.set(set_coeffs);

   return LocalExp;
};

//! Returns the shifted LocalExp of L to the new center X
template <int N>
ExpCoeffs<N> L2L(const ExpCoeffs<N>& L, const std::array<double, 3>& XIN) {
   ExpCoeffs<N> LocalExp(XIN);
   std::complex<double> R;
   std::array<double, 3> r = XIN - L.X;

   auto set_coeffs = [&](int n, int m) -> std::array<std::complex<double>, 2> {
      int m_;
      std::array<std::complex<double>, 2> ret = {0, 0};
      for (int n_ = n; n_ <= N; ++n_) {
         m_ = 0;
         R = SolidHarmonicR(n_ - n, m_ - m, r);
         complex_fused_multiply_increment(R, L.get_coeffs(n_, m_), std::get<0>(ret));
         complex_fused_multiply_increment(R, L.get_coeffs_(n_, m_), std::get<1>(ret));
         for (m_ = 1; m_ <= n_; ++m_) {
            R = SolidHarmonicR(n_ - n, m_ - m, r);
            complex_fused_multiply_increment(R, L.get_coeffs(n_, m_), std::get<0>(ret));
            complex_fused_multiply_increment(R, L.get_coeffs_(n_, m_), std::get<1>(ret));

            R = SolidHarmonicR(n_ - n, -m_ - m, r);
            complex_fused_multiply_increment(R, L.get_coeffs(n_, -m_), std::get<0>(ret));
            complex_fused_multiply_increment(R, L.get_coeffs_(n_, -m_), std::get<1>(ret));
         }
      }
      return ret;
   };

   LocalExp.set(set_coeffs);

   return LocalExp;
}

/* -------------------------------------------------------------------------- */

//! (-1)^mを含むルジャンドル陪関数
double dPdx(const double k, const double m, const double x) {
   double tmp = 1. / std::sqrt(1 - x * x);
   return std::pow(-1., m) * tmp * (m * x * tmp * std::assoc_legendre(k, m, x) + std::assoc_legendre(k, m + 1, x));
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

// template <int max_n>
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

// double dPdx(const double k, const double m, const double x) {
//    double tmp = 1. / std::sqrt(1 - x * x);
//    return std::pow(-1., m) * tmp * (m * x * tmp * std::assoc_legendre(k, m, x) + std::assoc_legendre(k, m + 1, x));
// }

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
//
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