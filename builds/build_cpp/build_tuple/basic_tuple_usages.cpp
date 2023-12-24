#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include "basic_vectors.hpp"

template <typename T>
struct test {
   using Tnd = tuple_of<std::tuple_size<T>::value, double>;
   Tnd X;
   tuple_of<std::tuple_size<T>::value, Tnd> J;
   test(const T& Xin) : X(Xin){};
};

int main() {
   /*
   tuple's Helper classes
   ------------------------
   std::tuple_size
   */
   std::tuple<double, double, double> u = {1., 1., 1.};
   constexpr auto n = std::tuple_size<decltype(u)>::value;  // タプルの長さを取得
   std::cout << "n=" << n << std::endl;
   /* -------------------------------------------------------------------------- */
   std::tuple<decltype(u), decltype(u), decltype(u)> M = {{1., 1., 1.},
                                                          {1., 1., 1.},
                                                          {1., 1., 1.}};
   constexpr auto m = std::tuple_size<decltype(M)>::value;  // タプルの長さを取得
   std::cout << "m=" << m << std::endl;
   /* -------------------------------------------------------------------------- */
   /*
    クラスに与えられたタプルの長さに応じて，同じ長さの2階のタプルを自動で準備するクラス．
   */
   {
      tuple_of<3, double> u;
      constexpr auto n = std::tuple_size<decltype(u)>::value;  // タプルの長さを取得
      std::cout << "n=" << n << std::endl;
   }
   {
      tuple_tuple_of<3, double> u;
      constexpr auto n = std::tuple_size<decltype(u)>::value;  // タプルの長さを取得
      std::cout << "n=" << n << std::endl;
      auto e = std::get<0>(u);
      constexpr auto n0 = std::tuple_size<decltype(e)>::value;  // タプルの長さを取得
      std::cout << "n0=" << n0 << std::endl;
   }
   {
      tuple_of<3, double> u;
      test t(u);
   }
   {
      tuple_of<3, double> u;
      std::cout << Inverse(tuple_of<3, tuple_of<3, double>>{{1, 2, 3}, {3, 4, 3}, {6, 4, 2}}) << std::endl;
   }
}