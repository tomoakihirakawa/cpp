// #include "Network.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"

template <typename T, size_t N>
std::array<T, N> TriShape(T t0, T t1) {
   static_assert(N == 3 || N == 6, "Unsupported shape function size. Only 3 or 6 are supported.");
   auto t2 = 1 - t0 - t1;
   if constexpr (N == 3) {
      return {t0, t1, t2};
   } else if constexpr (N == 6) {
      return {t0 * (2 * t0 - 1),
              t1 * (2 * t1 - 1),
              t2 * (2 * t2 - 1),
              4 * t0 * t1,
              4 * t1 * t2,
              4 * t0 * t2};
   }
}

template <size_t N>
std::array<double, N> TriShape(double t0, double t1) { return TriShape<double, N>(t0, t1); }

template <typename T, size_t N>
std::array<T, N> ModTriShape(T t0, T t1) {
   // b! どのポインターにどれだけの係数を加えるかを得られるようにする
   // b!　これは，線形補間であっても利用できる
   //
   static_assert(N == 3 || N == 6, "Unsupported shape function size. Only 3 or 6 are supported.");
   auto t2 = 1 - t0 - t1;
   auto t0m1 = t0 - 1;
   auto t1m1 = t1 - 1;
   if constexpr (N == 3) {
      return {t0, -t1 * t0m1, t0m1 * t1m1};
   } else if constexpr (N == 6) {
      //   return {t0 * (-1 + 2 * t0),
      //           t0m1 * t1 * (1 + 2 * t0m1 * t1),
      //           t1m1 * t0m1 * (1 + 2 * t1m1 * t0 - 2 * t1),
      //           -4 * t0m1 * t0 * t1,
      //           -4 * std::pow(t0m1, 2) * t1m1 * t1,
      //           4 * t1m1 * t0m1 * t0};

      std::array<T, N> deflt{{t0 * (-1 + 2 * t0),
                              t0m1 * t1 * (1 + 2 * t0m1 * t1),
                              t1m1 * t0m1 * (1 + 2 * t1m1 * t0 - 2 * t1),
                              -4 * t0m1 * t0 * t1,
                              -4 * std::pow(t0m1, 2) * t1m1 * t1,
                              4 * t1m1 * t0m1 * t0}};

      auto [M0_0, M0_1, M0_2] = deflt[0] * TriShape<3>(1., 1.);

      auto [M1_0, M1_1, M1_2] = deflt[1] * TriShape<3>(1., 1.);

      auto [M2_0, M2_1, M2_2] = deflt[2] * TriShape<3>(1., 1.);

      deflt[0] = 0;
      deflt[3] += M0_1;
      deflt[4] += M0_2;
      deflt[5] += M0_0;

      //   deflt[1] = 0;
      //   deflt[3] += M1_0;
      //   deflt[4] += M1_1;
      //   deflt[5] += M1_2;

      //   deflt[2] = 0;
      //   deflt[4] += M2_0;
      //   deflt[5] += M2_1;
      //   deflt[3] += M2_2;

      return deflt;
   }
}

template <size_t N>
std::array<double, N> ModTriShape(double t0, double t1) { return ModTriShape<double, N>(t0, t1); }

int main() {

   auto filename = "triangular_surface.dat";
   std::ofstream file(filename);
   if (!file.is_open()) {
      std::cerr << "Failed to open the output file." << std::endl;
      return 0;
   }
   file << "# x y z" << std::endl;
   double resolution = 10;
   double step = 1 / resolution;
   std::array<double, 3> a{{0, 0, 0}};
   std::array<double, 3> b{{1, 0, 0}};
   std::array<double, 3> c{{0, 1, 1}};
   //    {
   //       std::array<std::array<double, 3>, 3> M{{a, b, c}};
   //       for (size_t i = 0; i <= resolution; ++i) {
   //          double t0 = step * i;
   //          for (size_t j = 0; j <= resolution; ++j) {
   //             double t1 = step * j;
   //             auto [x, y, z] = Dot(ModTriShape<double, 3>(t0, t1), M);
   //             file << x << " " << y << " " << z << std::endl;
   //          }
   //       }
   //    }
   {
      std::array<double, 3> d = {0.5, 0.0, 0.5};
      std::array<double, 3> e = {0.5, 0.5, 0.5};
      std::array<double, 3> f = {0.0, 0.5, 0.5};
      std::array<std::array<double, 3>, 6> M{{a, b, c, d, e, f}};
      for (size_t i = 0; i <= resolution; ++i) {
         double t0 = step * i;
         for (size_t j = 0; j <= resolution; ++j) {
            double t1 = step * j;
            auto [x, y, z] = Dot(ModTriShape<6>(t0, t1), M);
            file << x << " " << y << " " << z << std::endl;
         }
      }
   }
}
