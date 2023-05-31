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


| | A(5,5,5) | A(10,10,10) |
|:----:|:---:|:---:|
| **n=3** | ![n3_A_5_5_5](output_n3_A_5_5_5.png)  | ![n3_A_10_10_10](output_n3_A_10_10_10.png) |
| **n=6** | ![n6_A_5_5_5](output_n6_A_5_5_5.png)  | ![n6_A_10_10_10](output_n6_A_10_10_10.png) |
| **n=9** | ![n9_A_5_5_5](output_n9_A_5_5_5.png)  | ![n9_A_10_10_10](output_n9_A_10_10_10.png) |


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

   for (double x = -20.0; x <= 20.0; x += .5) {
      for (double y = -20.0; y <= 20.0; y += .5) {
         double z = 0;
         std::array<double, 3> center = {x, y, z};
         auto error = std::log(std::abs(G(X - center, A - center) - G_approx(9, X - center, A - center)));
         std::cout << x << " " << y << " " << 0 << " " << error << "\n";
      }
      std::cout << "\n";
   }
}
