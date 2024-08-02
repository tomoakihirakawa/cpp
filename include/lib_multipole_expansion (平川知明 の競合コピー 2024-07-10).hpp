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
#include "lib_multipole_expansion_constants.hpp"

double G(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) { return 1 / Norm(X_near - X_far); }

std::array<double, 3> gradG(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) {
   return -(X_near - X_far) / std::pow(Norm(X_near - X_far), 3.);
}

/*_DOC_EXTRACT Multipole_Expansion

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

constexpr double
factorial(const int n) {
   if (n < 0)
      throw std::runtime_error("factorial of negative number");
   else if (n < factorial_list_as_double.size())
      return factorial_list_as_double[n];
   else
      return n * factorial(n - 1);
}

/* -------------------------------------------------------------------------- */

/*

\cite{Liu_2009}p.70には，体球調和関数を使った高速多重極展開のアルゴリズムが書かれている．
以下は，その定義通りの体球調和関数を計算する関数である．

*/

constexpr std::complex<double> complex_zero(0., 0.);
constexpr std::array<std::complex<double>, 3> complex_zero3D{complex_zero, complex_zero, complex_zero};

std::complex<double> Dot(const std::array<double, 3>& normal, const std::array<std::complex<double>, 3>& abc) {
   return std::get<0>(abc) * std::get<0>(normal) + std::get<1>(abc) * std::get<1>(normal) + std::get<2>(abc) * std::get<2>(normal);
};

std::complex<double> Dot(const std::array<std::complex<double>, 3>& abc, const std::array<double, 3>& normal) {
   return Dot(normal, abc);
};

std::array<std::complex<double>, 3> operator+(std::array<std::complex<double>, 3> a, const std::array<double, 3>& A) {
   std::get<0>(a) += std::get<0>(A);
   std::get<1>(a) += std::get<1>(A);
   std::get<2>(a) += std::get<2>(A);
   return a;
};

std::array<std::complex<double>, 3> operator+(const std::array<double, 3>& A, const std::array<std::complex<double>, 3>& a) {
   return a + A;
};

struct SphericalCoordinates {
   const std::array<double, 3> X;
   const double r2D;
   const double rho, theta, phi;
   const double div_rho, div_r2D;
   const double x_r, y_r, z_r;

   SphericalCoordinates(const std::array<double, 3>& X)
       : X(X),
         r2D(std::hypot(std::get<0>(X), std::get<1>(X))),
         rho(std::hypot(std::get<0>(X), std::get<1>(X), std::get<2>(X))),
         div_r2D(1. / r2D),
         div_rho(1. / rho),
         theta(std::atan2(r2D, std::get<2>(X))),
         phi(std::atan2(std::get<1>(X), std::get<0>(X))),
         x_r(std::get<0>(X) / rho),
         y_r(std::get<1>(X) / rho),
         z_r(std::get<2>(X) / rho) {}

   std::complex<double> sph_harmonics(const int l, const int m) const {
      auto error = [&]() {
         std::stringstream ss;
         ss << "m = " << m << " is out of range for l = " << l;
         throw std::runtime_error(ss.str());
      };
      double s = std::sqrt(4.0 * M_PI / (double)(2 * l + 1)) * std::pow(-1.0, m);
      if (std::abs(m) <= l) {
         if (m < 0) {
            if (l < 0) error();
            return std::polar(std::sph_legendre(l, -m, theta) * s, m * phi);
         } else {
            if (l < 0) error();
            return std::polar(std::sph_legendre(l, m, theta) * s, m * phi);
         }
      } else {
         return 0.0;
      }
   }

   std::complex<double> SolidHarmonicR(const int n, const int m) const { return std::pow(rho, n) * sph_harmonics(n, -m); }

   std::complex<double> SolidHarmonicS(const int n, const int m) const { return std::pow(rho, -n - 1) * sph_harmonics(n, m); }

   std::array<std::complex<double>, 2> SolidHarmonicR_Grad_SolidHarmonicR_normal(const int n, const int m, const std::array<double, 3>& normal) const {

      auto Rnm = std::pow(rho, n) * sph_harmonics(n, -m);
      //
      auto Rnm1 = std::pow(rho, n) * sph_harmonics(n, -m - 1);
      const std::array<double, 3> grad_r = {x_r, y_r, z_r};
      const std::array<double, 3> grad_theta = {x_r * z_r * div_r2D, y_r * z_r * div_r2D, -r2D * div_rho * div_rho};
      const std::array<double, 3> grad_phi = {-std::get<1>(X) * div_r2D * div_r2D, std::get<0>(X) * div_r2D * div_r2D, 0.};
      double c = std::sqrt((n - std::abs(m)) * (n + std::abs(m) + 1));
      auto grad_Rnm1_dot_normal = (double)n / rho * Rnm * Dot(normal, grad_r);
      // auto B = ((double)m * (std::get<2>(X) / r2D) /*/ std::tan(theta)*/ * Rnm + std::polar(c, phi) * SolidHarmonicR(n, m + 1)) * Dot(normal, grad_theta);
      grad_Rnm1_dot_normal += ((double)m * (std::get<2>(X) / r2D) * Rnm + std::polar(c, phi) * Rnm1) * Dot(normal, grad_theta);
      grad_Rnm1_dot_normal += -std::complex<double>(1.) * (double)m * Rnm * Dot(normal, grad_phi);
      //
      return {Rnm, grad_Rnm1_dot_normal};
   }
};

template <typename T>
void complex_fused_multiply_increment(const std::complex<T>& a, const std::complex<T>& b, std::complex<T>& c) {
   c.real(c.real() + a.real() * b.real() - a.imag() * b.imag());
   c.imag(c.imag() + a.real() * b.imag() + a.imag() * b.real());
}

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

   /* -------------------------------------------------------------------------- */

   std::array<double, 3> X;
   std::array<std::complex<double>, (N + 1) * (N + 1)> coeffs;
   std::array<std::complex<double>, (N + 1) * (N + 1)> coeffs_;
   // std::array<std::array<int, 2>, (N + 1) * (N + 1)> index_map;
   std::array<std::array<int, 2>, (N + 1) * (N + 1)> nm_set;

   std::vector<std::tuple<int, int, int, int, std::complex<double>>> jknm_M2L;

   const std::complex<double> zero = {0., 0.};

   int index(int n, int m) const { return n * (n + 1) + m; }

   const std::complex<double>& get_coeffs(const int n, const int m) const { return coeffs[index(n, m)]; }

   const std::complex<double>& get_coeffs_(const int n, const int m) const { return coeffs_[index(n, m)]; }

   void set_nm_set() {
      int ind = 0;
      for (int n = 0; n <= N; ++n)
         for (int m = -n; m <= n; ++m)
            nm_set[ind++] = {n, m};

      jknm_M2L.clear();
      for (int j = 0; j <= N; ++j)
         for (int k = -j; k <= j; ++k)
            for (int n = 0; n <= N; ++n)
               for (int m = -n; m <= n; ++m) {
                  auto AAA = AAA_M2L_FMM[n][m + N_AAA_M2L_FMM][j][k + N_AAA_M2L_FMM];
                  if (AAA.real() != 0.0 || AAA.imag() != 0.0)
                     jknm_M2L.push_back({j, k, n, m, AAA});
               }
   }

   void initialize(const std::array<double, 3>& XIN) {
      this->X = XIN;
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      set_nm_set();
   }

   // Default constructor
   ExpCoeffs() : X{0.0, 0.0, 0.0} {
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      set_nm_set();
   }

   ExpCoeffs(const std::array<double, 3>& XIN) : X(XIN) {
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      set_nm_set();
   }

   // void increment(const Tddd& XIN, const std::array<double, 2> weights, const Tddd& normal) {
   //    SphericalCoordinates R(XIN - this->X);
   //    auto set_coeffs = [&](int n, int m) -> std::array<std::complex<double>, 2> {
   //       return {R.SolidHarmonicR(n, m) * std::get<0>(weights), Dot(normal, R.Grad_SolidHarmonicR(n, m)) * std::get<1>(weights)};
   //    };
   //    this->increment_coeffs(set_coeffs);
   // }

   template <typename T>
   void increment_moments(const T& poles) {
      for (const auto& pole : poles) {
         // SphericalCoordinates R(this->X - pole->X);
         SphericalCoordinates R(pole->X - this->X);
         if constexpr (std::is_pointer_v<typename T::value_type> || std::is_same_v<typename T::value_type, std::shared_ptr<pole4FMM>>) {
            this->increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
               // return {R.SolidHarmonicR(n, m) * std::get<0>(pole->weights),
               //         R.Grad_SolidHarmonicR_normal(n, m, pole->normal) * std::get<1>(pole->weights)};
               return R.SolidHarmonicR_Grad_SolidHarmonicR_normal(n, m, pole->normal) * pole->weights;
            });
         } else {
            this->increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
               // return {R.SolidHarmonicR(n, m) * std::get<0>(pole.weights),
               //         R.Grad_SolidHarmonicR_normal(n, m, pole.normal) * std::get<1>(pole.weights)};
               return R.SolidHarmonicR_Grad_SolidHarmonicR_normal(n, m, pole->normal) * pole.weights;
            });
         }
      }
   }

   void increment_coeffs(const auto& func) {
      std::array<std::complex<double>, 2> cc;
      int ind = 0;
      for (const auto& [n, m] : this->nm_set) {
         cc = func(n, m);
         coeffs[ind] += std::get<0>(cc);
         coeffs_[ind] += std::get<1>(cc);
         ind++;
      }
   }

   std::array<std::complex<double>, 2> IGIGn_using_L(const std::array<double, 3>& a) const {
      // SphericalCoordinates Xc2O(this->X - a);
      SphericalCoordinates P(a - this->X);
      std::array<std::complex<double>, 2> ret = {0, 0};
      for (const auto& [n, m] : this->nm_set) {
         // auto R = P.SolidHarmonicR(n, -m);
         auto R = std::pow(P.rho, n) * P.sph_harmonics(n, m);
         complex_fused_multiply_increment(R, this->get_coeffs(n, m) /*L*/, ret[0]);
         complex_fused_multiply_increment(R, this->get_coeffs_(n, m) /*L*/, ret[1]);
      }
      return ret;
   }

   /* -------------------------------------------------------------------------- */

   void M2M(const ExpCoeffs<N>& M) {
      // SphericalCoordinates r(M_shifted.X - M.X);
      SphericalCoordinates rab(M.X - this->X);
      std::complex<double> R, c = {1., 0.}, AAA;
      this->increment_coeffs([&](int j, int k) -> std::array<std::complex<double>, 2> {
         std::array<std::complex<double>, 2> ret = {0, 0};
         for (int n = 0; n <= j; ++n)
            for (int m = -n; m <= n; ++m) {
               AAA = AAA_M2M_FMM[n][m + N_AAA_M2M_FMM][j][k + N_AAA_M2M_FMM];
               if (AAA.real() != 0.0 || AAA.imag() != 0.0) {
                  R = AAA * std::pow(rab.rho, n) * rab.sph_harmonics(n, -m);
                  complex_fused_multiply_increment(R, M.get_coeffs(j - n, k - m), std::get<0>(ret));
                  complex_fused_multiply_increment(R, M.get_coeffs_(j - n, k - m), std::get<1>(ret));
               }
            }
         return ret;
      });
   }

   void M2L(const ExpCoeffs<N>& M) {
      std::complex<double> S_, AAA;
      // SphericalCoordinates r(LocalExp.X - M.X);
      SphericalCoordinates rab(M.X - this->X);
      this->increment_coeffs([&](int j, int k) -> std::array<std::complex<double>, 2> {
         std::array<std::complex<double>, 2> ret = {0, 0};
         for (int n = 0; n <= N; ++n)
            for (int m = -n; m <= n; ++m) {
               AAA = AAA_M2L_FMM[n][m + N_AAA_M2L_FMM][j][k + N_AAA_M2L_FMM];
               if (AAA.real() != 0.0 || AAA.imag() != 0.0) {
                  S_ = AAA * rab.sph_harmonics(j + n, m - k) / std::pow(rab.rho, j + n + 1);
                  complex_fused_multiply_increment(S_, M.get_coeffs(n, m), std::get<0>(ret));
                  complex_fused_multiply_increment(S_, M.get_coeffs_(n, m), std::get<1>(ret));
               }
            }
         return ret;
      });
   }

   void L2L(const ExpCoeffs<N>& L) {
      std::complex<double> R, AAA;
      // SphericalCoordinates r(LocalExp.X - L.X);
      SphericalCoordinates rab(L.X - this->X);
      // this->increment_coeffs([&](int j, int k) -> std::array<std::complex<double>, 2> {
      //    std::array<std::complex<double>, 2> ret = {0, 0};
      //    for (int n = j; n <= N; ++n)
      //       for (int m = -n; m <= n; ++m) {
      //          AAA = AAA_L2L_FMM[n][m + N_AAA_L2L_FMM][j][k + N_AAA_L2L_FMM];
      //          if (AAA.real() != 0.0 || AAA.imag() != 0.0) {
      //             R = AAA * rab.sph_harmonics(n - j, m - k) * std::pow(rab.rho, n - j);
      //             complex_fused_multiply_increment(R, L.get_coeffs(n, m), std::get<0>(ret));
      //             complex_fused_multiply_increment(R, L.get_coeffs_(n, m), std::get<1>(ret));
      //          }
      //       }
      //    return ret;
      // });

      for (const auto& [j, k, n, m, AAA] : this->jknm_M2L) {
         if (AAA.real() != 0.0 || AAA.imag() != 0.0) {
            R = AAA * rab.sph_harmonics(n - j, m - k) * std::pow(rab.rho, n - j);
            complex_fused_multiply_increment(R, L.get_coeffs(n, m), this->get_coeffs(j, k));
            complex_fused_multiply_increment(R, L.get_coeffs_(n, m), this->get_coeffs_(j, k));
         }
      }
   }
};

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

// std::array<double, 3> gradGapx(unsigned p,
//                                const std::array<double, 3>& X_near_IN,
//                                const std::array<double, 3>& X_far_IN,
//                                const std::array<double, 3>& center) {
//    SphericalCoordinates sph_near(X_near_IN - center);
//    SphericalCoordinates sph_far(X_far_IN - center);
//    int k, m;

//    std::array<std::complex<double>, 3> grad = {0., 0., 0.};

//    for (k = 0; k <= p; ++k)
//       for (m = -k; m <= k; ++m)
//          grad += sph_near.Grad_SolidHarmonicR(k, m) * sph_far.SolidHarmonicS(k, m);

//    return {grad[0].real(), grad[1].real(), grad[2].real()};
// }

std::array<double, 3> ToSphericalCoordinates(const double x, const double y, const double z) {
   return {std::hypot(x, y, z), std::atan2(std::hypot(x, y), z), std::atan2(y, x)};
};

// \label{ToSphericalCoordinates}
std::array<double, 3> ToSphericalCoordinates(const std::array<double, 3>& X) {
   return ToSphericalCoordinates(std::get<0>(X), std::get<1>(X), std::get<2>(X));
};

#endif