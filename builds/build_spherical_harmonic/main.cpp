#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include "basic_arithmetic_array_operations.hpp"

// Compute the factorial of a_near given number
constexpr int factorial(const int n) {
   if (n < 0)
      throw std::runtime_error("factorial of negative number");
   if (n == 0 || n == 1)
      return 1;
   return n * factorial(n - 1);
}

// Compute the spherical harmonic function sph(k, m, a_far, b_far)
constexpr std::complex<double> Y(const int k, const int m, const double a_far, const double b_far) {
   if (k < 0 || abs(m) > k) {
      return std::complex<double>(0.0, 0.0);
   }
   const double assocLegendre = std::sqrt(static_cast<double>(factorial(k - abs(m)) / factorial(k + abs(m)))) * std::pow(-1, m) * std::assoc_legendre(k, abs(m), std::cos(a_far));
   const double realPart = assocLegendre * std::cos(m * b_far);
   const double imagPart = -assocLegendre * std::sin(m * b_far);
   return std::complex<double>(realPart, imagPart);
}

std::complex<double> G(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) { return 1 / Norm(X_near - X_far); }

/*DOC_EXTRACT BEM

## 多重極展開(Multipole Expansion)

Green関数を次のようにする．

$$
G({\bf x},{\bf a}) = \frac{1}{\|{\bf x}-{\bf a}\|}
$$

近似解 $G_{approx}({\bf x- \bf c},{\bf a - \bf c})$ を以下の式で定義する：

$$
G_{approx}(n, {\bf x- \bf c},{\bf a - \bf c}) \approx \sum_{k=0}^{n} \sum_{m=-k}^{k} \left( \frac{r_{near}}{r_{far}} \right)^k \frac{1}{r_{far}} Y(k, -m, a_{near}, b_{near}) Y(k, m, a_{far}, b_{far})
$$

$$
Y(k, m, a, b) = \sqrt{\frac{(k - |m|)!}{(k + |m|)!}}(-1)^m P_k^{|m|}(\cos(a)) \left(\cos(m b) - i \sin(m b)\right)
$$

ここで，

- $Y(k, m, a, b)$ は球面調和関数
- $r_{near}$ と $r_{far}$ はベクトル ${\bf x - c}$ と ${\bf a - c}$ のノルム
- $a_{near}$, $b_{near}$, $a_{far}$, $b_{far}$ はベクトル ${\bf x - c}$ と ${\bf a - c}$ の球面座標

${\bf c}=(x,y,0)$を変化させてプロットした結果：

| | **n=3** | **n=6** | **n=9** |
|:----:|:---:|:---:|:---:|
| **$`{\bf x} = (0,0,0),{\bf a} = (5,5,5)`$** | ![n3_A_5_5_5](./output_n3_A_5_5_5.png) | ![n6_A_5_5_5](./output_n6_A_5_5_5.png) | ![n9_A_5_5_5](./output_n9_A_5_5_5.png) |
| **$`{\bf x} = (0,0,0),{\bf a} = (10,10,10)`$** | ![n3_A_10_10_10](./output_n3_A_10_10_10.png) | ![n6_A_10_10_10](./output_n6_A_10_10_10.png) | ![n9_A_10_10_10](./output_n9_A_10_10_10.png) |

この結果からわかるように，Green関数の実際の値は，${\bf c}$によって変わらないが，$G_{approx}$の値は${\bf c}$によって変化し，
${\bf c}$が${\bf x}$に近いところでは，$G_{approx}$の値は$G$の値に近づく．

$a_{near},b_{near}$は，より小さければ精度が良く，
また，$a_{far},b_{far}$は，より大きければ精度が良くなる．

*/

std::array<double, 3> ToSphericalCoordinates(const std::array<double, 3>& X) {
   return {Norm(X),
           std::atan2(std::sqrt(std::get<0>(X) * std::get<0>(X) + std::get<1>(X) * std::get<1>(X)), std::get<2>(X)),
           std::atan2(std::get<1>(X), std::get<0>(X))};
};

std::complex<double> G_approx(unsigned p,
                              const std::array<double, 3>& X_near,
                              const std::array<double, 3>& X_far) {
   auto [r_near, a_near, b_near] = ToSphericalCoordinates(X_near);
   auto [r_far, a_far, b_far] = ToSphericalCoordinates(X_far);

   auto r_far_inv_r = 1. / r_far;
   auto r_near_inv_r = r_near * r_far_inv_r;
   std::complex<double> accum = 0;
   for (int k = 0; k <= p; ++k)
      for (int m = -k; m <= k; ++m)
         accum += std::pow(r_near_inv_r, k) * r_far_inv_r * Y(k, -m, a_near, b_near) * Y(k, m, a_far, b_far);
   return accum;
}

// int main() {
//    std::array<double, 3> X_near = {1, 1, 0};
//    std::array<double, 3> X_far = {100, 100, 100};
//    std::cout << "  G  ,   G_approx " << std::endl;
//    for (int i = 1; i < 40; i++) {
//       std::cout << G(X_near, X_far) << " , " << G_approx(i, X_near, X_far) << "\n";
//    }
// }

int main() {
   std::array<double, 3> A = {5, 5, 5};
   std::array<double, 3> X = {0, 0, 0};

   for (double x = -20.0; x <= 20.0; x += .1) {
      for (double y = -20.0; y <= 20.0; y += .1) {
         double z = 0;
         std::array<double, 3> center = {x, y, z};
         auto error = std::log(std::abs(G(X - center, A - center) - G_approx(9, X - center, A - center)));
         std::cout << x << " " << y << " " << 0 << " " << error << "\n";
      }
      std::cout << "\n";
   }
}
