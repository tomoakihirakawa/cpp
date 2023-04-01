#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <type_traits>
// #include "basic.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "lib_measurement.hpp"

/*
The main issue with this code is that it's trying to deduce the lambda function types using template type arguments, which is not possible.
*/

int main() {
   // Test cross product
   std::array<double, 3> vec0 = {1.0, 2.0, 3.0};
   std::array<double, 3> vec1 = {1.0, 1.0, 1.0};
   std::array<double, 4> vec2 = {4.0, 5.0, 6.0, 6.0};
   std::array<std::array<double, 3>, 4> M{{{2.0, 2.0, 2.0},
                                           {3.0, 3.0, 3.0},
                                           {4.0, 4.0, 4.0},
                                           {4.0, 4.0, 4.0}}};
   std::cout << Dot(Cross(vec1, vec0), vec1) << std::endl;
   std::cout << Dot(M, vec1) << std::endl;
   std::cout << Dot(vec2, M) << std::endl;

   return 0;
}
