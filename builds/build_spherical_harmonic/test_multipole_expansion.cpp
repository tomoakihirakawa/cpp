#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include "basic_arithmetic_array_operations.hpp"

// Compute the factorial of a given number
int factorial(int n) {
   if (n == 0 || n == 1)
      return 1;
   return n * factorial(n - 1);
}

// Compute the spherical harmonic function sph(k, m, theta, phi)
std::complex<double> Y(int k, int m, double theta, double phi) {
   if (k < 0 || abs(m) > k) {
      return std::complex<double>(0.0, 0.0);
   }

   double assocLegendre = std::sqrt(static_cast<double>(factorial(k - abs(m))) / (factorial(k + abs(m)))) * std::pow(-1., m) * std::assoc_legendre(k, abs(m), std::cos(theta));
   double realPart = assocLegendre * std::cos(m * phi);
   double imagPart = -assocLegendre * std::sin(m * phi);

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
   auto a = std::atan2(sqrt(x0 * x0 + y0 * y0), z0);
   auto b = std::atan2(y0, x0);
   //
   auto r1 = Norm(x1y1z1);
   auto theta = std::atan2(sqrt(x1 * x1 + y1 * y1), z1);
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

   std::array<double, 3> A = {0, 0, 0};
   std::array<double, 3> V0 = {0, 0, 0};
   std::array<double, 3> V1 = {100, 100, 100};
   std::array<double, 3> x0y0z0 = V0 + A;
   std::array<double, 3> x1y1z1 = V1 + A;

   // Open a file to output the data
   std::ofstream dataFile("fmm_error_data.dat");

   // Write headers for the data file
   dataFile << "# Order Error" << std::endl;

   // Calculate the exact potential
   auto exact_potential = G(x0y0z0, x1y1z1);

   for (int i = 1; i < 100; i++) {
      // Calculate the approximated potential
      auto approx_potential = G_approx(i, x0y0z0, x1y1z1);

      // Calculate the error between the exact potential and the approximated potential
      double error = std::abs(exact_potential - approx_potential);
      std::cout << "Order: " << i << " Error: " << error << std::endl;

      // Write the data to the file
      dataFile << i << " " << error << std::endl;
   }

   // Close the data file
   dataFile.close();

   // Create a Gnuplot script file
   std::ofstream scriptFile("fmm_error_plot.gp");
   scriptFile << "set xlabel 'Expansion Order'\n";
   scriptFile << "set ylabel 'Error'\n";
   scriptFile << "set logscale y\n";
   scriptFile << "set title 'FMM Potential Approximation Error'\n";
   scriptFile << "plot 'fmm_error_data.dat' using 1:2 with linespoints title 'Error'\n";
   scriptFile.close();

   return 0;
}
