// TypeTraitsExample1.cpp
#include <array>
#include <iostream>
#include <tuple>
#include <type_traits>

template <typename T>
struct is_double_array : std::false_type {};

template <std::size_t N>
struct is_double_array<std::array<double, N>> : std::true_type {};

int main() {
   int x = 5;
   double y = 10.0;
   int arr[3] = {1, 2, 3};
   std::array<int, 3> arr_std_int = {1, 2, 3};
   std::array<double, 3> arr_std_double = {1.0, 2.0, 3.0};
   std::tuple<int, double, char> tup = std::make_tuple(1, 2.0, 'a');

   std::cout << "Is x of integral type? " << std::boolalpha << std::is_integral_v<decltype(x)> << std::endl;
   std::cout << "Is y of integral type? " << std::boolalpha << std::is_integral_v<decltype(y)> << std::endl;
   std::cout << "Is arr a C-style array? " << std::boolalpha << std::is_array_v<decltype(arr)> << std::endl;
   std::cout << "Is arr_std_int an std::array of ints? " << std::boolalpha << std::is_same_v<decltype(arr_std_int), std::array<int, 3>> << std::endl;
   std::cout << "Is arr_std_double an std::array of doubles (regardless of size)? " << std::boolalpha << is_double_array<decltype(arr_std_double)>::value << std::endl;
   std::cout << "Is tup an std::tuple? " << std::boolalpha << std::is_same_v<decltype(tup), std::tuple<int, double, char>> << std::endl;
   std::cout << "Is the type of x the same as y? " << std::boolalpha << std::is_same_v<decltype(x), decltype(y)> << std::endl;

   return 0;
}
