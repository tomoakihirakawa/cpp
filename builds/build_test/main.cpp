#include "basic.hpp"

int main() {

   Tddd n = {1E-20, 1E-20, 1E-20};
   // Tddd n = {0, 0, 0};
   // n += (double)1;
   std::cout << "(int)0 = " << std::setprecision(20) << (int)0 << std::endl;
   std::cout << "(double)0 = " << std::setprecision(20) << (double)0 << std::endl;
   std::cout << "n = " << std::setprecision(20) << n << std::endl;
   std::cout << "Total(n) = " << std::setprecision(20) << Total(n) << std::endl;
   std::cout << "Norm(n) = " << std::setprecision(20) << Norm(n) << std::endl;
   std::cout << "Normalize(n) = " << std::setprecision(20) << Normalize(n) << std::endl;
   std::cout << std::setprecision(20) << 1. - Norm(Normalize(n)) << std::endl;
}