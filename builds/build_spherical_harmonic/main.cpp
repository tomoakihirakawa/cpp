#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include "basic_arithmetic_array_operations.hpp"

// Compute the factorial of a given number
constexpr int factorial(const int n) {
   if (n < 0)
      throw std::runtime_error("factorial of negative number");
   if (n == 0 || n == 1)
      return 1;
   return n * factorial(n - 1);
}

// Compute the spherical harmonic function sph(k, m, theta, phi)
constexpr std::complex<double> Y(const int k, const int m, const double theta, const double phi) {
   if (k < 0 || abs(m) > k) {
      return std::complex<double>(0.0, 0.0);
   }
   const double assocLegendre = std::sqrt(static_cast<double>(factorial(k - abs(m)) / factorial(k + abs(m)))) * std::pow(-1, m) * std::assoc_legendre(k, abs(m), std::cos(theta));
   const double realPart = assocLegendre * std::cos(m * phi);
   const double imagPart = -assocLegendre * std::sin(m * phi);
   return std::complex<double>(realPart, imagPart);
}

std::complex<double> G(const std::array<double, 3>& x0y0z0, const std::array<double, 3>& x1y1z1) { return 1 / Norm(x0y0z0 - x1y1z1); }

std::complex<double> G_approx(unsigned p,
                              const std::array<double, 3>& x0y0z0,
                              const std::array<double, 3>& x1y1z1) {
   auto [x0, y0, z0] = x0y0z0;
   auto [x1, y1, z1] = x1y1z1;
   //
   auto r0 = Norm(x0y0z0);
   auto a = std::atan2(std::sqrt(x0 * x0 + y0 * y0), z0);
   auto b = std::atan2(y0, x0);
   //
   auto r1 = Norm(x1y1z1);
   auto theta = std::atan2(std::sqrt(x1 * x1 + y1 * y1), z1);
   auto phi = std::atan2(y1, x1);

   auto inv_r = 1. / r1;
   auto r0_inv_r = r0 * inv_r;
   std::complex<double> accum = 0;
   for (int k = 0; k <= p; ++k)
      for (int m = -k; m <= k; ++m)
         accum += std::pow(r0_inv_r, k) * inv_r * Y(k, -m, a, b) * Y(k, m, theta, phi);
   return accum;
}

int main() {
   std::array<double, 3> x0y0z0 = {1, 0, 0};
   std::array<double, 3> x1y1z1 = {100, 100, 100};
   std::cout << "  G  ,   G_approx " << std::endl;
   for (int i = 1; i < 100; i++) {
      std::cout << G(x0y0z0, x1y1z1) << " , " << G_approx(i, x0y0z0, x1y1z1) << "\n";
   }
}