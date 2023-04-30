#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <type_traits>
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "lib_measurement.hpp"

int main() {

   std::array<double, 2> arr = {3, 4};
   std::array<std::array<double, 2>, 2> M0 = {{{1, 2}, {3, 4}}};
   std::array<std::array<double, 2>, 2> M1 = {{{5, 6}, {7, 8}}};
   std::cout << RotateLeft(arr, 2) << std::endl;
   std::cout << RotateRight(arr, 3) << std::endl;
   std::cout << Join(arr, arr) << std::endl;
   std::cout << Norm(arr) << std::endl;
   std::cout << Normalize(arr) << std::endl;
   std::cout << Dot(Normalize(arr), Normalize(arr)) << std::endl;
   std::cout << arr / (Norm(arr)) << std::endl;
   std::cout << Total(Normalize(arr)) << std::endl;
   std::cout << Dot(M0, M1) << std::endl;
   std::cout << Dot(M0, arr) << std::endl;
   std::cout << Dot(M1, arr) << std::endl;

   std::array<std::array<double, 3>, 3> M = {{{51, 6, 7}, {8, 9, 10}, {11, 12, 13}}};
   std::cout << Inverse(M) << std::endl;
}