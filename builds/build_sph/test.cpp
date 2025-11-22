#include <iostream>
#include <limits>

int main() {
   std::cout << "SIZE_MAX: " << SIZE_MAX << std::endl;
   std::cout << "static_cast<size_t>(-1): " << static_cast<size_t>(-1) << std::endl;
   std::cout << "static_cast<size_t>(-1): " << static_cast<size_t>(1) << std::endl;
   std::cout << "static_cast<size_t>(-1): " << static_cast<int>(1) << std::endl;
   return 0;
}