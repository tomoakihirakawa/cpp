
/*

cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=study_fma_performance.cpp

*/

#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
   double a = 1.1111111111111111e10;
   double b = 1.9e-10;
   double c = -1.1111111111111111e10;
   double traditional_result = 0.0;
   double fma_result = 0.0;

   for (int i = 0; i < 100000; ++i) {
      traditional_result += (a * b) + c;            // Traditional computation
      fma_result = std::fma(a, b, fma_result + c);  // FMA computation
   }
   std::cout << std::setw(25) << "Traditional Result: " << std::setprecision(17) << traditional_result << "\n";
   std::cout << std::setw(25) << "FMA Result: " << std::setprecision(17) << fma_result << "\n";
}