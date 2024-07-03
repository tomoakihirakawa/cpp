#ifndef lib_multipole_expansion_H
#define lib_multipole_expansion_H

#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
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
static const std::array<double, 21> factorial_list_as_double = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880., 3.6288e6, 3.99168e7, 4.790016e8, 6.2270208e9, 8.71782912e10, 1.307674368e12,
                                                                2.0922789888e13, 3.55687428096e14, 6.402373705728e15, 1.21645100408832e17, 2.43290200817664e18};

constexpr double factorial(const int n) {
   if (n < 0)
      throw std::runtime_error("factorial of negative number");
   else if (n < factorial_list_as_double.size())
      return factorial_list_as_double[n];
   else
      return n * factorial(n - 1);
}

// struct SphericalCoordinates {

//    const std::array<double, 3> X;
//    const double r2D;
//    const double rho, theta, phi;
//    const double div_rho, div_r2D;

//    SphericalCoordinates(const std::array<double, 3>& X)
//        : X(X),
//          r2D(std::hypot(std::get<0>(X), std::get<1>(X))),
//          rho(std::hypot(std::get<0>(X), std::get<1>(X), std::get<2>(X))),
//          div_r2D(1. / r2D),
//          div_rho(1. / rho),
//          theta(std::atan2(std::hypot(std::get<0>(X), std::get<1>(X)), std::get<2>(X))),
//          phi(std::atan2(std::get<1>(X), std::get<0>(X))) {}

//    std::array<std::array<double, 3>, 3> JacobianCartesian2Spherical() const {
//       const double x_r = std::get<0>(X) * div_rho;
//       const double y_r = std::get<1>(X) * div_rho;
//       const double z_r = std::get<2>(X) * div_rho;
//       const std::array<double, 3> grad_r = {x_r, y_r, z_r};
//       const std::array<double, 3> grad_theta = {x_r * z_r * div_r2D,
//                                                 y_r * z_r * div_r2D,
//                                                 -r2D * div_rho * div_rho};
//       const std::array<double, 3> grad_phi = {-std::get<1>(X) * div_r2D * div_r2D, std::get<0>(X) * div_r2D * div_r2D, 0.};
//       return {grad_r, grad_theta, grad_phi};
//    }
// };

/*

\cite{Liu_2009}p.70には，体球調和関数を使った高速多重極展開のアルゴリズムが書かれている．
以下は，その定義通りの体球調和関数を計算する関数である．

*/

constexpr std::complex<double> complex_zero(0., 0.);
constexpr std::array<std::complex<double>, 3> complex_zero3D{complex_zero, complex_zero, complex_zero};

double AssocLegendre(const int n, int m, const double x) {
   if (n < std::abs(m))
      return 0;
   double a = 1.0;
   if (m < 0) {
      m = -m;  //! m is now positive
      if (m % 2 == 0)
         a = factorial(n - m) / factorial(n + m);
      else
         a = -factorial(n - m) / factorial(n + m);
   }

   if (x < 0) {
      if ((n + m) % 2 == 0)
         return a * std::assoc_legendre(n, m, -x);
      else
         return -a * std::assoc_legendre(n, m, -x);
   } else
      return a * std::assoc_legendre(n, m, x);
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

std::array<std::complex<double>, 3> operator*(const std::complex<double>& a, const std::array<double, 3>& A) {
   return {a * A[0], a * A[1], a * A[2]};
};

std::array<std::complex<double>, 3> operator+(const std::array<std::complex<double>, 3>& a, const std::array<double, 3>& A) {
   return {a[0] + A[0], a[1] + A[1], a[2] + A[2]};
};

std::array<std::complex<double>, 3> operator+(const std::array<double, 3>& A, const std::array<std::complex<double>, 3>& a) {
   return {a[0] + A[0], a[1] + A[1], a[2] + A[2]};
};

// 球面調和関数
// std::complex<double> sph_harmonics(int l, int m, double theta, double phi) {
//    const double s = std::sqrt(4.0 * M_PI / (2.0 * l + 1.));
//    if (m < 0) {
//       m = -m;
//       return ((m % 2 == 0) ? 1.0 : -1.0) * std::polar(std::sph_legendre(l, m, theta) * s, -m * phi);
//    } else
//       return std::polar(std::sph_legendre(l, m, theta) * s, m * phi);
// }

// CHRISTOPHE FOCHESATO AND FRE´DE´RIC DIAS
// std::complex<double> sph_harmonics(const int l, const int m, const double theta, const double phi) {
//    auto error = [&]() {
//       std::stringstream ss;
//       ss << "m = " << m << " is out of range for l = " << l;
//       throw std::runtime_error(ss.str());
//    };
//    const double s = std::sqrt(2.0 * M_PI / (l + 0.5));
//    if (m < 0) {
//       if (l < 0 || -m > l) error();
//       return std::polar(std::sph_legendre(l, -m, theta) * s, m * phi);
//    } else {
//       if (l < 0 || m > l) error();
//       return std::polar(std::sph_legendre(l, m, theta) * s, m * phi);
//    }
// }

struct SphericalCoordinates {
   const std::array<double, 3> X;
   const double r2D;
   const double rho, theta, phi;
   const double div_rho, div_r2D;

   SphericalCoordinates(const std::array<double, 3>& X)
       : X(X),
         r2D(std::hypot(std::get<0>(X), std::get<1>(X))),
         rho(std::hypot(std::get<0>(X), std::get<1>(X), std::get<2>(X))),
         div_r2D(1. / r2D),
         div_rho(1. / rho),
         theta(std::atan2(r2D, std::get<2>(X))),
         phi(std::atan2(std::get<1>(X), std::get<0>(X))) {}

   std::complex<double> sph_harmonics(const int l, const int m) const {
      auto error = [&]() {
         std::stringstream ss;
         ss << "m = " << m << " is out of range for l = " << l;
         throw std::runtime_error(ss.str());
      };
      const double s = std::sqrt(2.0 * M_PI / (l + 0.5));
      if (m < 0) {
         if (l < 0 || -m > l) error();
         return std::polar(std::sph_legendre(l, -m, theta) * s, m * phi);
      } else {
         if (l < 0 || m > l) error();
         return std::polar(std::sph_legendre(l, m, theta) * s, m * phi);
      }
   }

   std::complex<double> conj_sph_harmonics(const int l, const int m) const {
      auto error = [&]() {
         std::stringstream ss;
         ss << "m = " << m << " is out of range for l = " << l;
         throw std::runtime_error(ss.str());
      };
      const double s = std::sqrt(2.0 * M_PI / (l + 0.5));
      if (m < 0) {
         if (l < 0 || -m > l) error();
         return std::polar(std::sph_legendre(l, -m, theta) * s, -m * phi);
      } else {
         if (l < 0 || m > l) error();
         return std::polar(std::sph_legendre(l, m, theta) * s, -m * phi);
      }
   }

   std::complex<double> SolidHarmonicR(const int n, const int m) const {
      return n < std::abs(m) ? std::complex<double>(0., 0.) : (std::pow(rho, n) * sph_harmonics(n, -m));
   }

   std::complex<double> conj_SolidHarmonicR(const int n, const int m) const {
      return n < std::abs(m) ? std::complex<double>(0., 0.) : (std::pow(rho, n) * conj_sph_harmonics(n, -m));
   }

   std::complex<double> SolidHarmonicS(const int n, const int m) const {
      return n < std::abs(m) ? std::complex<double>(0.0, 0.0) : (std::pow(rho, -n - 1.) * sph_harmonics(n, m));
   }

   std::complex<double> conj_SolidHarmonicS(const int n, const int m) const {
      return n < std::abs(m) ? std::complex<double>(0.0, 0.0) : (std::pow(rho, -n - 1.) * conj_sph_harmonics(n, m));
   }

   std::complex<double> sph_harmonics_divided_by_sph_harmonics(const int l, const int m, const int m2) const {
      int M = m < 0 ? -m : m;
      int M2 = m2 < 0 ? -m2 : m2;
      if (M == M2)
         return std::polar(1., (m - m2) * phi /*ここはこのままでOK*/);
      else
         return std::polar(std::sph_legendre(l, M, theta) / std::sph_legendre(l, M2, theta), (m - m2) * phi /*ここはこのままでOK*/);
   }

   std::array<std::complex<double>, 3> Grad_SolidHarmonicR(const int n, const int m) const {
      double c = std::sqrt((m < 0) ? (n + m) * (n - m + 1.) : (n - m) * (n + m + 1.));
      return std::pow(rho, n - 2) * sph_harmonics(n, -m) * ((double)n * X

                                                            + ((double)m / std::tan(theta) + sph_harmonics_divided_by_sph_harmonics(n, -m - 1, -m) * std::polar(c, phi)) * std::array<double, 3>{{std::get<0>(X) * std::get<2>(X) * div_r2D, std::get<1>(X) * std::get<2>(X) * div_r2D, -r2D}}

                                                            + std::pow(rho * div_r2D, 2) * std::complex<double>(m) * std::array<double, 3>{{-std::get<1>(X), std::get<0>(X), 0.}});
   }
};

std::array<double, 3> ToSphericalCoordinates(const double x, const double y, const double z) {
   return {std::hypot(x, y, z), std::atan2(std::hypot(x, y), z), std::atan2(y, x)};
};

// \label{ToSphericalCoordinates}
std::array<double, 3> ToSphericalCoordinates(const std::array<double, 3>& X) {
   return ToSphericalCoordinates(std::get<0>(X), std::get<1>(X), std::get<2>(X));
};

template <typename T>
void complex_fused_multiply_increment(const std::complex<T>& a, const std::complex<T>& b, std::complex<T>& c) {
   // c = std::complex<T>(std::fma(a.real(), b.real(), std::fma(-a.imag(), b.imag(), c.real())),
   //                     std::fma(a.real(), b.imag(), std::fma(a.imag(), b.real(), c.imag())));
   c.real(c.real() + a.real() * b.real() - a.imag() * b.imag());
   c.imag(c.imag() + a.real() * b.imag() + a.imag() * b.real());
}
//

// //! (-1)^mを含まないルジャンドル陪関数
// double D_assoc_legendre(int n, int m, const double x) {
//    return (-m * x / (1. - x * x) * AssocLegendre(n, m, x) + AssocLegendre(n, m + 1, x) / std::sqrt(1. - x * x));
// }

// std::array<std::complex<double>, 3> Grad_SolidHarmonicR(int n, int m, const SphericalCoordinates& sph) {
//    return sph.Grad_SolidHarmonicR(n, m);
// }

/* -------------------------------------------------------------------------- */

#include "basic_exception.hpp"

struct pole4FMM {
   Tddd X;
   Tdd weights;
   Tddd normal;
   pole4FMM(const Tddd& X, const Tdd& weights, const Tddd& normal) : X(X), weights(weights), normal(normal) {}
};

template <int N>
struct ExpCoeffs {

   /*

   N = 0 [ 0,  ,  ,  ,  ] n=0   total 1 = (N+1)^2             center : 0, in 1D 0
   N = 1 [-1, 0, 1,  ,  ] n=1   total 1 + 3 = 4 = (N+1)^2     center : 1, in 1D 2
   N = 2 [-2,-1, 0, 1, 2] n=2=N total 4 + 5 = 9 = (N+1)^2     center : 2, in 1D 6
         <-- 2*N+1=5 -->
   */

   std::array<double, 3> X;
   std::array<std::complex<double>, (N + 1) * (N + 1)> coeffs;
   std::array<std::complex<double>, (N + 1) * (N + 1)> coeffs_;
   std::array<std::array<int, 2>, (N + 1) * (N + 1)> index_map;

   const std::complex<double> zero = {0., 0.};

   int index(int n, int m) const { return n * (n + 1) + m; }

   const std::complex<double>& get_coeffs(const int n, const int m) const { return coeffs[index(n, m)]; }

   const std::complex<double>& get_coeffs_(const int n, const int m) const { return coeffs_[index(n, m)]; }

   void initialize(const std::array<double, 3>& XIN) {
      this->X = XIN;
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      int n, m;
      for (n = 0; n <= N; ++n)
         for (m = -n; m <= n; ++m)
            index_map[index(n, m)] = {n, m};
   }

   // Default constructor
   ExpCoeffs() : X{0.0, 0.0, 0.0} {
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      int n, m;
      for (n = 0; n <= N; ++n)
         for (m = -n; m <= n; ++m)
            index_map[index(n, m)] = {n, m};
   }

   ExpCoeffs(const std::array<double, 3>& XIN) : X(XIN) {
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      int n, m;
      for (n = 0; n <= N; ++n)
         for (m = -n; m <= n; ++m)
            index_map[index(n, m)] = {n, m};
   }

   void increment(const Tddd& XIN, const std::array<double, 2> weights, const Tddd& normal) {
      SphericalCoordinates R(XIN - this->X);
      auto set_coeffs = [&](int n, int m) -> std::array<std::complex<double>, 2> {
         return {R.SolidHarmonicR(n, m) * std::get<0>(weights), Dot(normal, R.Grad_SolidHarmonicR(n, m)) * std::get<1>(weights)};
      };
      this->increment_coeffs(set_coeffs);
   }

   template <typename T>
   void increment(const T& poles) {
      for (const auto& pole : poles) {
         SphericalCoordinates R(pole->X - this->X);
         if constexpr (std::is_pointer_v<typename T::value_type> || std::is_same_v<typename T::value_type, std::shared_ptr<pole4FMM>>) {
            this->increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
               return {R.SolidHarmonicR(n, m) * std::get<0>(pole->weights),
                       Dot(pole->normal, R.Grad_SolidHarmonicR(n, m)) * std::get<1>(pole->weights)};
            });
         } else {
            this->increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
               return {R.SolidHarmonicR(n, m) * std::get<0>(pole.weights),
                       Dot(pole.normal, R.Grad_SolidHarmonicR(n, m)) * std::get<1>(pole.weights)};
            });
         }
      }
   }

   void increment_coeffs(const auto& func) {
      std::array<std::complex<double>, 2> cc;
      int ind = 0;
      for (const auto& [n, m] : this->index_map) {
         cc = func(n, m);
         coeffs[ind] += std::get<0>(cc);
         coeffs_[ind] += std::get<1>(cc);
         ind++;
      }
   }

   std::complex<double> IG_using_M(const std::array<double, 3>& O) const {
      SphericalCoordinates Xc2O(O - this->X);
      std::complex<double> ret = 0;
      for (const auto& [n, m] : this->index_map)
         complex_fused_multiply_increment(Xc2O.conj_SolidHarmonicS(n, m), this->get_coeffs(n, m), ret);
      return ret;
   }

   std::complex<double> IG_using_L(const std::array<double, 3>& O) const {
      SphericalCoordinates Xc2O(O - this->X);
      std::complex<double> ret = 0;
      for (const auto& [n, m] : this->index_map)
         complex_fused_multiply_increment(Xc2O.SolidHarmonicR(n, m), this->get_coeffs(n, m), ret);
      return ret;
   }

   std::complex<double> IGn_using_M(const std::array<double, 3>& O) const {
      SphericalCoordinates Xc2O(O - this->X);
      std::complex<double> ret = 0;
      for (const auto& [n, m] : this->index_map)
         complex_fused_multiply_increment(Xc2O.conj_SolidHarmonicS(n, m), this->get_coeffs_(n, m), ret);
      return ret;
   }

   std::complex<double> IGn_using_L(const std::array<double, 3>& O) const {
      SphericalCoordinates Xc2O(O - this->X);
      std::complex<double> ret = 0;
      for (const auto& [n, m] : this->index_map)
         complex_fused_multiply_increment(Xc2O.SolidHarmonicR(n, m), this->get_coeffs_(n, m), ret);
      return ret;
   }
};

template <int N>
void M2M(const ExpCoeffs<N>& M, ExpCoeffs<N>& M_shifted) {
   try {
      SphericalCoordinates r(M.X - M_shifted.X);
      std::complex<double> R;
      M_shifted.increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
         std::array<std::complex<double>, 2> ret = {0, 0};
         for (int n_ = 0; n_ <= n; ++n_)
            for (int m_ = -n_; m_ <= n_; ++m_)
               if (std::abs(m_) <= n_ && 0 <= M.index(n - n_, m - m_) && M.index(n - n_, m - m_) < (N + 1) * (N + 1)) {
                  R = r.SolidHarmonicR(n_, m_);
                  complex_fused_multiply_increment(R, M.get_coeffs(n - n_, m - m_), std::get<0>(ret));
                  complex_fused_multiply_increment(R, M.get_coeffs_(n - n_, m - m_), std::get<1>(ret));
               }
         return ret;
      });
   } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
      std::stringstream ss;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   }
}

template <int N>
void M2L(const ExpCoeffs<N>& M, ExpCoeffs<N>& LocalExp) {
   try {
      std::complex<double> S_;
      SphericalCoordinates r(LocalExp.X - M.X);
      std::array<std::complex<double>, 2> ret = {0, 0};
      int ind = 0;
      double sign;
      LocalExp.increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
         ret = {0, 0};
         ind = 0;
         sign = (n % 2 == 0) ? 1.0 : -1.0;
         for (const auto& [n_, m_] : M.index_map) {
            if (std::abs(m_ + m) <= n_ + n) {
               S_ = sign * r.conj_SolidHarmonicS(n + n_, m + m_);
               complex_fused_multiply_increment(S_, M.coeffs[ind], std::get<0>(ret));
               complex_fused_multiply_increment(S_, M.coeffs_[ind], std::get<1>(ret));
            }
            ind++;
         }
         return ret;
      });
   } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
      std::stringstream ss;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   }
}

template <int N>
void L2L(const ExpCoeffs<N>& L, ExpCoeffs<N>& LocalExp) {
   try {
      std::complex<double> R;
      SphericalCoordinates r(LocalExp.X - L.X);
      LocalExp.increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
         int m_;
         std::array<std::complex<double>, 2> ret = {0, 0};
         for (int n_ = n; n_ <= N; ++n_)
            for (m_ = -n_; m_ <= n_; ++m_)
               if (std::abs(m_ - m) <= n_ - n && 0 <= L.index(n_, m_) && L.index(n_, m_) < (N + 1) * (N + 1)) {

                  R = r.SolidHarmonicR(n_ - n, m_ - m);
                  complex_fused_multiply_increment(R, L.get_coeffs(n_, m_), std::get<0>(ret));
                  complex_fused_multiply_increment(R, L.get_coeffs_(n_, m_), std::get<1>(ret));
               }
         return ret;
      });
   } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
      std::stringstream ss;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   }
}

/* -------------------------------------------------------------------------- */

double Gapx(unsigned p,
            std::array<double, 3> X_near,
            std::array<double, 3> X_far,
            const std::array<double, 3>& center) {
   SphericalCoordinates sph_near(X_near - center);
   SphericalCoordinates sph_far(X_far - center);

   std::complex<double> accum = 0;
   int k, m;
   for (k = 0; k <= p; ++k)
      for (m = -k; m <= k; ++m)
         accum += sph_near.SolidHarmonicR(k, m) * sph_far.SolidHarmonicS(k, m);
   return accum.real();
}

std::array<double, 3> gradGapx(unsigned p,
                               const std::array<double, 3>& X_near_IN,
                               const std::array<double, 3>& X_far_IN,
                               const std::array<double, 3>& center) {
   SphericalCoordinates sph_near(X_near_IN - center);
   SphericalCoordinates sph_far(X_far_IN - center);
   int k, m;

   std::array<std::complex<double>, 3> grad = {0., 0., 0.};

   for (k = 0; k <= p; ++k)
      for (m = -k; m <= k; ++m)
         grad += sph_near.Grad_SolidHarmonicR(k, m) * sph_far.SolidHarmonicS(k, m);

   return {grad[0].real(), grad[1].real(), grad[2].real()};
}

#endif