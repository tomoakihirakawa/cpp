#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include "lib_multipole_expansion.hpp"

int main() {
   // Test the sph_harmonics_ function
   int l = 2;
   int m = 1;
   std::array<double, 3> X = {1., 2., 3.};
   SphericalCoordinates sph(X);
   //! real part
   int max = 15;
   for (int j = 0; j < 2; ++j)
      for (int k = -j; k <= j; ++k) {
         std::cout << "{";
         for (int n = 0; n <= max; ++n)
            for (int m = -n; m <= n; ++m) {
               auto sph_nm = sph.sph_harmonics_(j + n, m - k) * AAA_M2L_FMM[j][k + N_AAA_M2L_FMM][n][m + N_AAA_M2L_FMM];
               std::cout << "{" << n << "," << m << "," << sph_nm.real() << ((n == m && m == max) ? "}" : "},");
            }
         std::cout << "}" << std::endl;
         std::cout << std::endl;
      }

   //! imaginary part
   // std::cout << "{";
   // for (int j = 1; j < max; ++j)
   //    for (int k = -max; k < max; ++k) {
   //       auto sph_nm = sph.sph_harmonics_(j + n, m - k) * AAA_M2L_FMM[j][k + N_AAA_M2L_FMM][l][m + N_AAA_M2L_FMM];
   //       std::cout << "{" << l << "," << m << "," << sph_lm.imag() << ((m == l && m == max - 1) ? "}" : "},");
   //    }
   // std::cout << "}" << std::endl;

   return 0;
}
