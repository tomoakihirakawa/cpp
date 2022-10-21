#ifndef basic_arithmetic_vector_operations_H
#define basic_arithmetic_vector_operations_H

#include <algorithm>  //transformなど
#include <numeric>    //数値のシーケンスの処理に特化したアルゴリズム
#include <tuple>
#include <type_traits>
#include <vector>
#include "basic_alias.hpp"
/* ------------------------------------------------------ */

auto Power(const auto T, const auto N) { return std::pow(T, N); };
auto Sqrt(const auto T) { return std::sqrt(T); };
auto Tanh(const auto T) { return std::tanh(T); };
auto Cosh(const auto T) { return std::cosh(T); };
auto Sin(const auto T) { return std::sin(T); };
auto Cos(const auto T) { return std::cos(T); };
auto Sech(const auto T) { return 1 / std::cosh(T); };
auto Csch(const auto T) { return 1 / std::sinh(T); };
auto ArcTan(const auto T) { return std::atan(T); };
auto ArcTan(const auto Y, const auto X) { return std::atan2(Y, X); };

/* ------------------------------------------------------ */
/*                    タプルダブル7                         */
/* ------------------------------------------------------ */
T7d operator-(T7d v) {
   std::get<0>(v) *= -1;
   std::get<1>(v) *= -1;
   std::get<2>(v) *= -1;
   std::get<3>(v) *= -1;
   std::get<4>(v) *= -1;
   std::get<5>(v) *= -1;
   std::get<6>(v) *= -1;
   return v;
};
/* ------------------------------------------------------ */
T7d &operator-=(T7d &v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   std::get<2>(v) -= u;
   std::get<3>(v) -= u;
   std::get<4>(v) -= u;
   std::get<5>(v) -= u;
   std::get<6>(v) -= u;
   return v;
};
T7d &operator+=(T7d &v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   std::get<4>(v) += u;
   std::get<5>(v) += u;
   std::get<6>(v) += u;
   return v;
};
T7d &operator*=(T7d &v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   std::get<4>(v) *= u;
   std::get<5>(v) *= u;
   std::get<6>(v) *= u;
   return v;
};
T7d &operator/=(T7d &v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   std::get<2>(v) /= u;
   std::get<3>(v) /= u;
   std::get<4>(v) /= u;
   std::get<5>(v) /= u;
   std::get<6>(v) /= u;
   return v;
}; /* ------------------------------------------------------ */
T7d &operator-=(T7d &v, const T7d &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   std::get<3>(v) -= std::get<3>(u);
   std::get<4>(v) -= std::get<4>(u);
   std::get<5>(v) -= std::get<5>(u);
   std::get<6>(v) -= std::get<6>(u);
   return v;
};
T7d &operator+=(T7d &v, const T7d &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   std::get<3>(v) += std::get<3>(u);
   std::get<4>(v) += std::get<4>(u);
   std::get<5>(v) += std::get<5>(u);
   std::get<6>(v) += std::get<6>(u);
   return v;
};
T7d &operator*=(T7d &v, const T7d &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   std::get<3>(v) *= std::get<3>(u);
   std::get<4>(v) *= std::get<4>(u);
   std::get<5>(v) *= std::get<5>(u);
   std::get<6>(v) *= std::get<6>(u);
   return v;
};
T7d &operator/=(T7d &v, const T7d &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   std::get<3>(v) /= std::get<3>(u);
   std::get<4>(v) /= std::get<4>(u);
   std::get<5>(v) /= std::get<5>(u);
   std::get<6>(v) /= std::get<6>(u);
   return v;
};
/* ------------------------------------------------------ */
T7d operator-(T7d v, const T7d &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   std::get<3>(v) -= std::get<3>(u);
   std::get<4>(v) -= std::get<4>(u);
   std::get<5>(v) -= std::get<5>(u);
   std::get<6>(v) -= std::get<6>(u);
   return v;
};
T7d operator+(T7d v, const T7d &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   std::get<3>(v) += std::get<3>(u);
   std::get<4>(v) += std::get<4>(u);
   std::get<5>(v) += std::get<5>(u);
   std::get<6>(v) += std::get<6>(u);
   return v;
};
T7d operator*(T7d v, const T7d &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   std::get<3>(v) *= std::get<3>(u);
   std::get<4>(v) *= std::get<4>(u);
   std::get<5>(v) *= std::get<5>(u);
   std::get<6>(v) *= std::get<6>(u);
   return v;
};
T7d operator/(T7d v, const T7d &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   std::get<3>(v) /= std::get<3>(u);
   std::get<4>(v) /= std::get<4>(u);
   std::get<5>(v) /= std::get<5>(u);
   std::get<6>(v) /= std::get<6>(u);
   return v;
};
/* ------------------------------------------------------ */
T7d operator-(T7d v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   std::get<2>(v) -= u;
   std::get<3>(v) -= u;
   std::get<4>(v) -= u;
   std::get<5>(v) -= u;
   std::get<6>(v) -= u;
   return v;
};
T7d operator+(T7d v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   std::get<4>(v) += u;
   std::get<5>(v) += u;
   std::get<6>(v) += u;
   return v;
};
T7d operator*(T7d v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   std::get<4>(v) *= u;
   std::get<5>(v) *= u;
   std::get<6>(v) *= u;
   return v;
};
T7d operator/(T7d v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   std::get<2>(v) /= u;
   std::get<3>(v) /= u;
   std::get<4>(v) /= u;
   std::get<5>(v) /= u;
   std::get<6>(v) /= u;
   return v;
};
/* ------------------------------------------------------ */
T7d operator-(const double u, T7d v) {
   std::get<0>(v) = u - std::get<0>(v);
   std::get<1>(v) = u - std::get<1>(v);
   std::get<2>(v) = u - std::get<2>(v);
   std::get<3>(v) = u - std::get<3>(v);
   std::get<4>(v) = u - std::get<4>(v);
   std::get<5>(v) = u - std::get<5>(v);
   std::get<6>(v) = u - std::get<6>(v);
   return v;
};
T7d operator+(const double u, T7d v) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   std::get<4>(v) += u;
   std::get<5>(v) += u;
   std::get<6>(v) += u;
   return v;
};
T7d operator*(const double u, T7d v) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   std::get<4>(v) *= u;
   std::get<5>(v) *= u;
   std::get<6>(v) *= u;
   return v;
};
T7d operator/(const double u, T7d v) {
   std::get<0>(v) = u / std::get<0>(v);
   std::get<1>(v) = u / std::get<1>(v);
   std::get<2>(v) = u / std::get<2>(v);
   std::get<3>(v) = u / std::get<3>(v);
   std::get<4>(v) = u / std::get<4>(v);
   std::get<5>(v) = u / std::get<5>(v);
   std::get<6>(v) = u / std::get<6>(v);
   return v;
};

/* ------------------------------------------------------ */
/*                    タプルダブル6                         */
/* ------------------------------------------------------ */
T6d operator-(T6d v) {
   std::get<0>(v) *= -1;
   std::get<1>(v) *= -1;
   std::get<2>(v) *= -1;
   std::get<3>(v) *= -1;
   std::get<4>(v) *= -1;
   std::get<5>(v) *= -1;
   return v;
};
/* ------------------------------------------------------ */
T6d &operator-=(T6d &v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   std::get<2>(v) -= u;
   std::get<3>(v) -= u;
   std::get<4>(v) -= u;
   std::get<5>(v) -= u;
   return v;
};
T6d &operator+=(T6d &v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   std::get<4>(v) += u;
   std::get<5>(v) += u;
   return v;
};
T6d &operator*=(T6d &v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   std::get<4>(v) *= u;
   std::get<5>(v) *= u;
   return v;
};
T6d &operator/=(T6d &v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   std::get<2>(v) /= u;
   std::get<3>(v) /= u;
   std::get<4>(v) /= u;
   std::get<5>(v) /= u;
   return v;
}; /* ------------------------------------------------------ */
T6d &operator-=(T6d &v, const T6d &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   std::get<3>(v) -= std::get<3>(u);
   std::get<4>(v) -= std::get<4>(u);
   std::get<5>(v) -= std::get<5>(u);
   return v;
};
T6d &operator+=(T6d &v, const T6d &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   std::get<3>(v) += std::get<3>(u);
   std::get<4>(v) += std::get<4>(u);
   std::get<5>(v) += std::get<5>(u);
   return v;
};
T6d &operator*=(T6d &v, const T6d &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   std::get<3>(v) *= std::get<3>(u);
   std::get<4>(v) *= std::get<4>(u);
   std::get<5>(v) *= std::get<5>(u);
   return v;
};
T6d &operator/=(T6d &v, const T6d &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   std::get<3>(v) /= std::get<3>(u);
   std::get<4>(v) /= std::get<4>(u);
   std::get<5>(v) /= std::get<5>(u);
   return v;
};
/* ------------------------------------------------------ */
T6d operator-(T6d v, const T6d &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   std::get<3>(v) -= std::get<3>(u);
   std::get<4>(v) -= std::get<4>(u);
   std::get<5>(v) -= std::get<5>(u);
   return v;
};
T6d operator+(T6d v, const T6d &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   std::get<3>(v) += std::get<3>(u);
   std::get<4>(v) += std::get<4>(u);
   std::get<5>(v) += std::get<5>(u);
   return v;
};
T6d operator*(T6d v, const T6d &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   std::get<3>(v) *= std::get<3>(u);
   std::get<4>(v) *= std::get<4>(u);
   std::get<5>(v) *= std::get<5>(u);
   return v;
};
T6d operator/(T6d v, const T6d &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   std::get<3>(v) /= std::get<3>(u);
   std::get<4>(v) /= std::get<4>(u);
   std::get<5>(v) /= std::get<5>(u);
   return v;
};
/* ------------------------------------------------------ */
T6d operator-(T6d v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   std::get<2>(v) -= u;
   std::get<3>(v) -= u;
   std::get<4>(v) -= u;
   std::get<5>(v) -= u;
   return v;
};
T6d operator+(T6d v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   std::get<4>(v) += u;
   std::get<5>(v) += u;
   return v;
};
T6d operator*(T6d v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   std::get<4>(v) *= u;
   std::get<5>(v) *= u;
   return v;
};
T6d operator/(T6d v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   std::get<2>(v) /= u;
   std::get<3>(v) /= u;
   std::get<4>(v) /= u;
   std::get<5>(v) /= u;
   return v;
};
/* ------------------------------------------------------ */
T6d operator-(const double u, T6d v) {
   std::get<0>(v) = u - std::get<0>(v);
   std::get<1>(v) = u - std::get<1>(v);
   std::get<2>(v) = u - std::get<2>(v);
   std::get<3>(v) = u - std::get<3>(v);
   std::get<4>(v) = u - std::get<4>(v);
   std::get<5>(v) = u - std::get<5>(v);
   return v;
};
T6d operator+(const double u, T6d v) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   std::get<4>(v) += u;
   std::get<5>(v) += u;
   return v;
};
T6d operator*(const double u, T6d v) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   std::get<4>(v) *= u;
   std::get<5>(v) *= u;
   return v;
};
T6d operator/(const double u, T6d v) {
   std::get<0>(v) = u / std::get<0>(v);
   std::get<1>(v) = u / std::get<1>(v);
   std::get<2>(v) = u / std::get<2>(v);
   std::get<3>(v) = u / std::get<3>(v);
   std::get<4>(v) = u / std::get<4>(v);
   std::get<5>(v) = u / std::get<5>(v);
   return v;
};
/* ------------------------------------------------------ */
/*                    タプルダブル4                         */
/* ------------------------------------------------------ */
T4d operator-(T4d v) {
   std::get<0>(v) = -std::get<0>(v);
   std::get<1>(v) = -std::get<1>(v);
   std::get<2>(v) = -std::get<2>(v);
   std::get<3>(v) = -std::get<3>(v);
   return v;
};
T4d operator-(T4d v, const T4d &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   std::get<3>(v) -= std::get<3>(u);
   return v;
};
T4d operator+(T4d v, const T4d &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   std::get<3>(v) += std::get<3>(u);
   return v;
};
T4d operator*(T4d v, const T4d &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   std::get<3>(v) *= std::get<3>(u);
   return v;
};
T4d operator/(T4d v, const T4d &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   std::get<3>(v) /= std::get<3>(u);
   return v;
};
/* ------------------------------------------------------ */
T4d operator-(T4d v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   std::get<2>(v) -= u;
   std::get<3>(v) -= u;
   return v;
};
T4d operator+(T4d v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   return v;
};
T4d operator+(const double u, T4d v) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   return v;
};
T4d operator*(T4d v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   return v;
};
T4d operator*(const double u, T4d v) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   return v;
};
T4d operator/(const double u, T4d v) {
   std::get<0>(v) = u / std::get<0>(v);
   std::get<1>(v) = u / std::get<1>(v);
   std::get<2>(v) = u / std::get<2>(v);
   std::get<3>(v) = u / std::get<3>(v);
   return v;
};
T4d operator/(T4d v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   std::get<2>(v) /= u;
   std::get<3>(v) /= u;
   return v;
};
/* ------------------------------------------------------ */
T4d &operator-=(T4d &v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   std::get<2>(v) -= u;
   std::get<3>(v) -= u;
   return v;
};
T4d &operator+=(T4d &v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   std::get<3>(v) += u;
   return v;
};
T4d &operator*=(T4d &v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   std::get<3>(v) *= u;
   return v;
};
T4d &operator/=(T4d &v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   std::get<2>(v) /= u;
   std::get<3>(v) /= u;
   return v;
}; /* ------------------------------------------------------ */
T4d &operator-=(T4d &v, const T4d &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   std::get<3>(v) -= std::get<3>(u);
   return v;
};
T4d &operator+=(T4d &v, const T4d &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   std::get<3>(v) += std::get<3>(u);
   return v;
};
T4d &operator*=(T4d &v, const T4d &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   std::get<3>(v) *= std::get<3>(u);
   return v;
};
T4d &operator/=(T4d &v, const T4d &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   std::get<3>(v) /= std::get<3>(u);
   return v;
};
/* ------------------------------------------------------ */
/*                    タプルダブル3                         */
/* ------------------------------------------------------ */
Tddd operator-(Tddd v) {
   std::get<0>(v) = -std::get<0>(v);
   std::get<1>(v) = -std::get<1>(v);
   std::get<2>(v) = -std::get<2>(v);
   return v;
};
/* ------------------------------------------------------ */
Tddd &operator-=(Tddd &v, const Tddd &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   return v;
};
Tddd &operator-=(Tddd &v, const double d) {
   std::get<0>(v) -= d;
   std::get<1>(v) -= d;
   std::get<2>(v) -= d;
   return v;
};
Tddd &operator+=(Tddd &v, const Tddd &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   return v;
};
Tddd &operator+=(Tddd &v, const double d) {
   std::get<0>(v) += d;
   std::get<1>(v) += d;
   std::get<2>(v) += d;
   return v;
};
Tddd &operator*=(Tddd &v, const Tddd &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   return v;
};
Tddd &operator*=(Tddd &v, const double d) {
   std::get<0>(v) *= d;
   std::get<1>(v) *= d;
   std::get<2>(v) *= d;
   return v;
};
Tddd &operator/=(Tddd &v, const Tddd &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   return v;
};
Tddd &operator/=(Tddd &v, const double d) {
   std::get<0>(v) /= d;
   std::get<1>(v) /= d;
   std::get<2>(v) /= d;
   return v;
};
/* ------------------------------------------------------ */
Tddd operator-(Tddd v, const Tddd &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   return v;
};
Tddd operator+(Tddd v, const Tddd &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   return v;
};
Tddd operator*(Tddd v, const Tddd &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   return v;
};
Tddd operator/(Tddd v, const Tddd &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   return v;
};
/* ------------------------------------------------------ */
Tddd operator-(Tddd v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   std::get<2>(v) -= u;
   return v;
};
Tddd operator+(Tddd v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   return v;
};
Tddd operator*(Tddd v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   return v;
};
Tddd operator/(Tddd v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   std::get<2>(v) /= u;
   return v;
};
/* ------------------------------------------------------ */
std::vector<Tddd> operator-(std::vector<Tddd> V, const double u) {
   for (auto &v : V)
      v -= u;
   return V;
};
std::vector<Tddd> operator+(std::vector<Tddd> V, const double u) {
   for (auto &v : V)
      v += u;
   return V;
};
std::vector<Tddd> operator*(std::vector<Tddd> V, const double u) {
   for (auto &v : V)
      v *= u;
   return V;
};
std::vector<Tddd> operator/(std::vector<Tddd> V, const double u) {
   for (auto &v : V)
      v /= u;
   return V;
};
/* ------------------------------------------------------ */
Tddd operator-(const double u, Tddd v) {
   std::get<0>(v) = u - std::get<0>(v);
   std::get<1>(v) = u - std::get<1>(v);
   std::get<2>(v) = u - std::get<2>(v);
   return v;
};
Tddd operator+(const double u, Tddd v) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   return v;
};
Tddd operator*(const double u, Tddd v) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   return v;
};
Tddd operator/(const double u, Tddd v) {
   std::get<0>(v) = u / std::get<0>(v);
   std::get<1>(v) = u / std::get<1>(v);
   std::get<2>(v) = u / std::get<2>(v);
   return v;
};
/* ------------------------------------------------------ */
/*                    タプルダブル2                         */
/* ------------------------------------------------------ */
Tdd operator-(Tdd v) {
   std::get<0>(v) *= -1.;
   std::get<1>(v) *= -1.;
   return v;
};
/* ------------------------------------------------------ */
Tdd &operator-=(Tdd &v, const Tdd &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   return v;
};
Tdd &operator-=(Tdd &v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   return v;
};
Tdd &operator+=(Tdd &v, const Tdd &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   return v;
};
Tdd &operator+=(Tdd &v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   return v;
};
Tdd &operator*=(Tdd &v, const Tdd &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   return v;
};
Tdd &operator*=(Tdd &v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   return v;
};
Tdd &operator/=(Tdd &v, const Tdd &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   return v;
};
Tdd &operator/=(Tdd &v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   return v;
};
/* ------------------------------------------------------ */
Tdd operator-(Tdd v, const Tdd &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   return v;
};
Tdd operator+(Tdd v, const Tdd &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   return v;
};
Tdd operator*(Tdd v, const Tdd &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   return v;
};
Tdd operator/(Tdd v, const Tdd &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   return v;
};
/* ------------------------------------------------------ */
Tdd operator-(Tdd v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   return v;
};
Tdd operator+(Tdd v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   return v;
};
Tdd operator*(Tdd v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   return v;
};
Tdd operator/(Tdd v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   return v;
};
/* ------------------------------------------------------ */
Tdd operator-(const double u, Tdd v) {
   std::get<0>(v) = u - std::get<0>(v);
   std::get<1>(v) = u - std::get<1>(v);
   return v;
};
Tdd operator+(const double u, Tdd v) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   return v;
};
Tdd operator*(const double u, Tdd v) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   return v;
};
Tdd operator/(const double u, Tdd v) {
   std::get<0>(v) = u / std::get<0>(v);
   std::get<1>(v) = u / std::get<1>(v);
   return v;
};
/* ------------------------------------------------------ */
/*                     行列T3Tddd                          */
/* ------------------------------------------------------ */
// Tddd operator-(Tddd v)
// {
// 	std::get<0>(v) = -std::get<0>(v);
// 	std::get<1>(v) = -std::get<1>(v);
// 	std::get<2>(v) = -std::get<2>(v);
// 	return v;
// };
// /* ------------------------------------------------------ */
T3Tddd &operator-=(T3Tddd &v, const T3Tddd &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   return v;
};
T3Tddd &operator-=(T3Tddd &v, const double d) {
   std::get<0>(v) -= d;
   std::get<1>(v) -= d;
   std::get<2>(v) -= d;
   return v;
};
T3Tddd &operator+=(T3Tddd &v, const T3Tddd &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   return v;
};
T3Tddd &operator+=(T3Tddd &v, const double d) {
   std::get<0>(v) += d;
   std::get<1>(v) += d;
   std::get<2>(v) += d;
   return v;
};
T3Tddd &operator*=(T3Tddd &v, const T3Tddd &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   return v;
};
T3Tddd &operator*=(T3Tddd &v, const double d) {
   std::get<0>(v) *= d;
   std::get<1>(v) *= d;
   std::get<2>(v) *= d;
   return v;
};
T3Tddd &operator/=(T3Tddd &v, const T3Tddd &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   return v;
};
T3Tddd &operator/=(T3Tddd &v, const double d) {
   std::get<0>(v) /= d;
   std::get<1>(v) /= d;
   std::get<2>(v) /= d;
   return v;
};
/* ------------------------------------------------------ */
T3Tddd operator-(T3Tddd v, const T3Tddd &u) {
   std::get<0>(v) -= std::get<0>(u);
   std::get<1>(v) -= std::get<1>(u);
   std::get<2>(v) -= std::get<2>(u);
   return v;
};
T3Tddd operator+(T3Tddd v, const T3Tddd &u) {
   std::get<0>(v) += std::get<0>(u);
   std::get<1>(v) += std::get<1>(u);
   std::get<2>(v) += std::get<2>(u);
   return v;
};
T3Tddd operator*(T3Tddd v, const T3Tddd &u) {
   std::get<0>(v) *= std::get<0>(u);
   std::get<1>(v) *= std::get<1>(u);
   std::get<2>(v) *= std::get<2>(u);
   return v;
};
T3Tddd operator/(T3Tddd v, const T3Tddd &u) {
   std::get<0>(v) /= std::get<0>(u);
   std::get<1>(v) /= std::get<1>(u);
   std::get<2>(v) /= std::get<2>(u);
   return v;
};
/* ------------------------------------------------------ */
T3Tddd operator-(T3Tddd v, const double u) {
   std::get<0>(v) -= u;
   std::get<1>(v) -= u;
   std::get<2>(v) -= u;
   return v;
};
T3Tddd operator+(T3Tddd v, const double u) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   return v;
};
T3Tddd operator*(T3Tddd v, const double u) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   return v;
};
T3Tddd operator/(T3Tddd v, const double u) {
   std::get<0>(v) /= u;
   std::get<1>(v) /= u;
   std::get<2>(v) /= u;
   return v;
};
/* ------------------------------------------------------ */
T3Tddd operator-(const double u, T3Tddd v) {
   std::get<0>(v) = u - std::get<0>(v);
   std::get<1>(v) = u - std::get<1>(v);
   std::get<2>(v) = u - std::get<2>(v);
   return v;
};
T3Tddd operator+(const double u, T3Tddd v) {
   std::get<0>(v) += u;
   std::get<1>(v) += u;
   std::get<2>(v) += u;
   return v;
};
T3Tddd operator*(const double u, T3Tddd v) {
   std::get<0>(v) *= u;
   std::get<1>(v) *= u;
   std::get<2>(v) *= u;
   return v;
};
T3Tddd operator/(const double u, T3Tddd v) {
   std::get<0>(v) = u / std::get<0>(v);
   std::get<1>(v) = u / std::get<1>(v);
   std::get<2>(v) = u / std::get<2>(v);
   return v;
};
/* ------------------------------------------------------ */
std::vector<Tddd> operator-(std::vector<Tddd> V, const Tddd &u) {
   for (auto &v : V)
      v -= u;
   return V;
};
std::vector<Tddd> operator+(std::vector<Tddd> V, const Tddd &u) {
   for (auto &v : V)
      v += u;
   return V;
};
std::vector<Tddd> operator*(std::vector<Tddd> V, const Tddd &u) {
   for (auto &v : V)
      v *= u;
   return V;
};
std::vector<Tddd> operator/(std::vector<Tddd> V, const Tddd &u) {
   for (auto &v : V)
      v /= u;
   return V;
};
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(std::vector<T> ret) {
   std::transform(ret.begin(), ret.end(), ret.begin(), [](const T tmp) { return -tmp; });
   return ret;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> ret) {
   std::transform(ret.begin(), ret.end(), ret.begin(), [](const std::vector<T> &tmp) { return -tmp; });
   return ret;
};
///////////////////////////////////////////////////
// vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator*(std::vector<T> v, const T din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const T c) { return c * din; });
   return v;
};
//! matrix * scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator*(std::vector<std::vector<T>> v, const T din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &c) { return c * din; });
   return v;
};
//@ scaler * vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator*(const T din, std::vector<T> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const T c) { return c * din; });
   return v;
};
//! scaler * matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator*(const T din, std::vector<std::vector<T>> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &c) { return c * din; });
   return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<std::vector<T>>> operator*(std::vector<std::vector<std::vector<T>>> v, const T &din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<std::vector<T>> &c) { return c * din; });
   return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<std::vector<T>>> operator*(const T &din, std::vector<std::vector<std::vector<T>>> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<std::vector<T>> &c) { return c * din; });
   return v;
};
//@ vector * vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator*(std::vector<T> v, const std::vector<T> &w) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const T a, const T b) { return a * b; });
   return v;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator*=(std::vector<T> &v, const T &w) { return (v = v * w); };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator*=(std::vector<T> &v, const std::vector<T> &w) { return (v = v * w); };
// //vector x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>> operator*(const std::vector<T>& v, std::vector<std::vector<T>> w){
//   std::transform(w.begin(),w.end(),v.cbegin(),w.begin(),[](const T& a, const T& b){ return a*b; });
//   return w;
// };
// //matrix x vector
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator*(std::vector<std::vector<T>> w, const std::vector<T>& v){
//   std::transform(w.begin(),w.end(),v.cbegin(),w.begin(),[](const T& a, const T& b){ return a*b; });
//   return w;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>>& operator*=(std::vector<std::vector<T>>& v, const std::vector<T>& w){
//   return (v = v * w);
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>>& operator*=(std::vector<std::vector<T>>& v, const std::vector<std::vector<T>>& w){
//   return (v = v * w);
// };
/////////////////////////////////////////////////
//@vector / scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator/(std::vector<T> v, const T din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const T &tmp) { return tmp / din; });
   return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator/(std::vector<std::vector<T>> v, const T din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp) { return tmp / din; });
   return v;
};
// vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator/(const T din, std::vector<T> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const T &tmp) { return din / tmp; });
   return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator/(const T din, std::vector<std::vector<T>> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp) { return din / tmp; });
   return v;
};
// vector x vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator/(std::vector<T> v, const std::vector<T> &w) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](T a, T b) { return a / b; });
   return v;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator/=(std::vector<T> &v, const T &w) { return (v = v / w); };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator/=(std::vector<T> &v, const std::vector<T> &w) { return (v = v / w); };
// //vector x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator/(const std::vector<T>& w, const std::vector< std::vector<T> >& v){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return b/a; });
//   return ret;
// };
// //matrix x vector
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator/(const std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return a/b; });
//   return ret;
// };
// //matrix x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator/(const std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, std::vector<T> b){ return a/b; });
//   return ret;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> >& operator/=(std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   v = v / w;
//   return v;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> >& operator/=(std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   v = v / w;
//   return v;
// };
/////////////////////////////////////////////////
//@ vector - scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(std::vector<T> v, const T din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const T tmp) { return tmp - din; });
   return v;
};
//! matrix - scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> v, const T &din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp) { return tmp - din; });
   return v;
};
//@ scaler - vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(const T din, std::vector<T> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const T tmp) { return din - tmp; });
   return v;
};
//! scaler - matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(const T din, std::vector<std::vector<T>> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp) { return din - tmp; });
   return v;
};
//@ vector - vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(std::vector<T> v, const std::vector<T> &w) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const T a, const T b) { return a - b; });
   return v;
};
//! matrix - matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> v, const std::vector<std::vector<T>> &w) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const std::vector<T> &b) { return a - b; });
   return v;
};
//%matrix - vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> v, const std::vector<T> &w) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b) { return a - b; });
   return v;
};
//%vector - matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(const std::vector<T> &w, std::vector<std::vector<T>> v) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b) { return b - a; });  //逆にする
   return v;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator-=(std::vector<T> &v, const T &w) { return (v = v - w); };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator-=(std::vector<T> &v, const std::vector<T> &w) { return (v = v - w); };
//! matrix -= matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator-=(std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   return (v = v - w);
};
//%vector -= matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator-=(std::vector<T> &v, const std::vector<std::vector<T>> &w) {
   return (v = v - w);
};
//%matrix -= vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator-=(std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   return (v = v - w);
};
// //vector x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator-(const std::vector<T>& w, const std::vector< std::vector<T> >& v){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return b-a; });
//   return ret;
// };
// //matrix x vector
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator-(const std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return a-b; });
//   return ret;
// };
// //matrix x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator-(const std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, std::vector<T> b){ return a-b; });
//   return ret;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> >& operator-=(std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   v = v - w;
//   return v;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> >& operator-=(std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   v = v - w;
//   return v;
// };
/////////////////////////////////////////////////
//@vector + scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator+(std::vector<T> v, const T &din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const T &tmp) { return tmp + din; });
   return v;
};
//! matrix + scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(std::vector<std::vector<T>> v, const T &din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp) { return tmp + din; });
   return v;
};
//@scaler + vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator+(const T &din, std::vector<T> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const T &tmp) { return tmp + din; });
   return v;
};
//! scaler + matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(const T &din, std::vector<std::vector<T>> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp) { return tmp + din; });
   return v;
};
//@vector + vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator+(std::vector<T> v, const std::vector<T> &w) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const T &a, const T &b) { return a + b; });
   return v;
};
//! matrix + matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(std::vector<std::vector<T>> v, const std::vector<std::vector<T>> &w) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const std::vector<T> &b) { return a + b; });
   return v;
};
//%matrix + vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(std::vector<std::vector<T>> v, const std::vector<T> &w) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b) { return a + b; });
   return v;
};
//%vector + matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(const std::vector<T> &w, std::vector<std::vector<T>> v) {
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b) { return a + b; });
   return v;
};
//@vector += scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator+=(std::vector<T> &v, const T &w) {
   return (v = v + w);
};
//@vector += vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator+=(std::vector<T> &v, const std::vector<T> &w) {
   return (v = v + w);
};
//! matrix += matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator+=(std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   return (v = v + w);
};
//%vector += matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator+=(std::vector<T> &v, const std::vector<std::vector<T>> &w) {
   return (v = v + w);
};
//%matrix += vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator+=(std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   return (v = v + w);
};
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

namespace std {
namespace {

// Code from boost
// Reciprocal of the golden ratio helps spread entropy
//     and handles duplicates.
// See Mike Seymour in magic-numbers-in-boosthash-combine:
//     https://stackoverflow.com/questions/4948780

template <class T>
inline void hash_combine(std::size_t &seed, T const &v) {
   seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Recursive template code derived from Matthieu M.
template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
struct HashValueImpl {
   static void apply(size_t &seed, Tuple const &tuple) {
      HashValueImpl<Tuple, Index - 1>::apply(seed, tuple);
      hash_combine(seed, get<Index>(tuple));
   }
};

template <class Tuple>
struct HashValueImpl<Tuple, 0> {
   static void apply(size_t &seed, Tuple const &tuple) {
      hash_combine(seed, get<0>(tuple));
   }
};
}  // namespace

template <typename... TT>
struct hash<std::tuple<TT...>> {
   size_t
   operator()(std::tuple<TT...> const &tt) const {
      size_t seed = 0;
      HashValueImpl<std::tuple<TT...>>::apply(seed, tt);
      return seed;
   }
};

template <typename T1, typename T2, typename T3>
struct hash<std::tuple<T1, T2, T3>> {
   size_t
   operator()(const std::tuple<T1, T2, T3> &tt) const {
      size_t seed = 0;
      HashValueImpl<std::tuple<T1, T2, T3>>::apply(seed, tt);
      return seed;
   }
};

template <typename T1, typename T2>
struct hash<std::pair<T1, T2>> {
   size_t
   operator()(const std::pair<T1, T2> &tt) const {
      size_t seed = 0;
      HashValueImpl<std::pair<T1, T2>>::apply(seed, tt);
      return seed;
   }
};

}  // namespace std

#include <unordered_map>

std::unordered_map<std::tuple<int, int>, std::tuple<double, double>> operator*(const double din, const std::unordered_map<std::tuple<int, int>, std::tuple<double, double>> &v) {
   auto ret = v;
   for (auto &[ii, dd] : ret)
      dd *= din;
   return ret;
};

std::unordered_map<std::tuple<int, int>, std::tuple<double, double>> operator*(const std::unordered_map<std::tuple<int, int>, std::tuple<double, double>> &v,
                                                                               const double din) {
   auto ret = v;
   for (auto &[ii, dd] : ret)
      dd *= din;
   return ret;
};

std::unordered_map<std::tuple<int, int>, std::tuple<double, double>> operator/(const double din, const std::unordered_map<std::tuple<int, int>, std::tuple<double, double>> &v) {
   auto ret = v;
   double m = 1. / din;
   for (auto &[ii, dd] : ret)
      dd *= m;
   return ret;
};

std::unordered_map<std::tuple<int, int>, std::tuple<double, double>> operator/(const std::unordered_map<std::tuple<int, int>, std::tuple<double, double>> &v,
                                                                               const double din) {
   auto ret = v;
   double m = 1. / din;
   for (auto &[ii, dd] : ret)
      dd *= m;
   return ret;
};

std::unordered_map<Tii, Tdd> &operator+=(std::unordered_map<Tii, Tdd> &ii_dd, const std::unordered_map<Tii, Tdd> &jj_dd) {
   std::unordered_map<Tii, Tdd>::iterator it;
   for (auto &[jj, dd] : jj_dd)
      if ((it = ii_dd.find(jj)) != ii_dd.end())
         it->second += dd;
      else
         ii_dd[jj] = dd;
   return ii_dd;
};

   // std::unordered_map<Tii, Tddd> &operator+=(std::unordered_map<Tii, Tddd> &ii_dd, const std::unordered_map<Tii, Tddd> &jj_dd)
   // {
   // 	std::unordered_map<Tii, Tddd>::const_iterator it;
   // 	for (auto &[ii, dd] : ii_dd)
   // 		if ((it = jj_dd.find(ii)) != jj_dd.end())
   // 			dd += it->second;
   // 	return ii_dd;
   // };

#endif