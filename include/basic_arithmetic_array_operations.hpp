#ifndef basic_arithmetic_array_operations_H
#define basic_arithmetic_array_operations_H

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <unordered_set>
#include "basic_alias.hpp"

/*
   intended behavior:
   {{a, b, c}, {a, b, c}, {a, b, c}} + {a, b, c}
   = {a, b, c} + {{a, b, c}, {a, b, c}, {a, b, c}}
   = {{2 a, a + b, a + c}, {a + b, 2 b, b + c}, {a + c, b + c, 2 c}}

   {{a, b, c}, {a, b, c}, {a, b, c}} - {a, b, c}
   = {{0, -a + b, -a + c}, {a - b, 0, -b + c}, {a - c, b - c, 0}}

   {a, b, c} - {{a, b, c}, {a, b, c}, {a, b, c}}
   = {{0, a - b, a - c}, {-a + b, 0, b - c}, {-a + c, -b + c, 0}}

   {{a, b, c}, {a, b, c}, {a, b, c}}*{a, b, c}
   = {a, b, c}*{{a, b, c}, {a, b, c}, {a, b, c}}
   = {{a^2,a b,a c},{a b,b^2,b c},{a c,b c,c^2}}

   {{a, b, c}, {a, b, c}, {a, b, c}}/{a, b, c}
   = {a, b, c}/{{a, b, c}, {a, b, c}, {a, b, c}}
   = {{1,b/a,c/a},{a/b,1,c/b},{a/c,b/c,1}}

   {{a, b, c}, {a, b, c}, {a, b, c}}/{a, b, c}
   = {a, b, c}/{{a, b, c}, {a, b, c}, {a, b, c}}
   = {{1,b/a,c/a},{a/b,1,c/b},{a/c,b/c,1}}

   {a, b, c}/{{a, b, c}, {a, b, c}, {a, b, c}}
   = {{1,a/b,a/c},{b/a,1,b/c},{c/a,c/b,1}}
*/

template <typename T>
struct is_std_array : std::false_type {};

template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};
/* -------------------------------------------------------------------------- */
template <size_t N, typename T>
constexpr std::array<T, N> operator-(std::array<T, N> arr) noexcept {
   for (auto& a : arr)
      a = -a;
   return arr;
}
/* -------------------------------- addition -------------------------------- */
// 1d array + value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>&> operator+=(std::array<T, N>& arr /*ref*/, const TT& d) noexcept {
   for (auto& a : arr) a += d;
   return arr;
}
// 1d array + 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N>& operator+=(std::array<T, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) += (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}

template <typename T>
constexpr std::array<T, 3>& operator+=(std::array<T, 3>& arr /*ref*/, const std::array<T, 3>& ARR) noexcept {
   std::get<0>(arr) += std::get<0>(ARR);
   std::get<1>(arr) += std::get<1>(ARR);
   std::get<2>(arr) += std::get<2>(ARR);
   return arr;
}

// 2d array + 1d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator+=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) += (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}
// 2d array + 2d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator+=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) += (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}

template <typename T>
constexpr std::array<std::array<T, 3>, 3>& operator+=(std::array<std::array<T, 3>, 3>& arr /*ref*/, const std::array<std::array<T, 3>, 3>& ARR) noexcept {
   std::get<0>(std::get<0>(arr)) += std::get<0>(std::get<0>(ARR));
   std::get<0>(std::get<1>(arr)) += std::get<0>(std::get<1>(ARR));
   std::get<0>(std::get<2>(arr)) += std::get<0>(std::get<2>(ARR));
   std::get<1>(std::get<0>(arr)) += std::get<1>(std::get<0>(ARR));
   std::get<1>(std::get<1>(arr)) += std::get<1>(std::get<1>(ARR));
   std::get<1>(std::get<2>(arr)) += std::get<1>(std::get<2>(ARR));
   std::get<2>(std::get<0>(arr)) += std::get<2>(std::get<0>(ARR));
   std::get<2>(std::get<1>(arr)) += std::get<2>(std::get<1>(ARR));
   std::get<2>(std::get<2>(arr)) += std::get<2>(std::get<2>(ARR));
   return arr;
}

/* -------------------------------------------------------------------------- */
// 1d array + value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator+(std::array<T, N> arr /*copy*/, const TT& d) noexcept { return arr += d; }
// value + 1d array
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator+(const TT& d, std::array<T, N> arr /*copy*/) noexcept { return arr + d; }
// 1d array + 1d array
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator+(std::array<T, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr += ARR; }
// 2d array + 1d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator+(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept {
   return arr += ARR;
}
// 1d array + 2d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N> operator+(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept { return arr += ARR; }
/* ------------------------------- subtraction ------------------------------ */
// 1d array -= value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>&> operator-=(std::array<T, N>& arr /*ref*/, const TT& d) noexcept {
   for (auto& a : arr) a -= d;
   return arr;
}
// 1d array -= 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N>& operator-=(std::array<T, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) -= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}

template <typename T>
constexpr std::array<T, 3>& operator-=(std::array<T, 3>& arr /*ref*/, const std::array<T, 3>& ARR) noexcept {
   std::get<0>(arr) -= std::get<0>(ARR);
   std::get<1>(arr) -= std::get<1>(ARR);
   std::get<2>(arr) -= std::get<2>(ARR);
   return arr;
}

// 2d array -= 1d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator-=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) -= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}
// 2d array -= 2d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator-=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) -= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}

template <typename T>
constexpr std::array<std::array<T, 3>, 3>& operator-=(std::array<std::array<T, 3>, 3>& arr /*ref*/, const std::array<std::array<T, 3>, 3>& ARR) noexcept {
   std::get<0>(std::get<0>(arr)) -= std::get<0>(std::get<0>(ARR));
   std::get<0>(std::get<1>(arr)) -= std::get<0>(std::get<1>(ARR));
   std::get<0>(std::get<2>(arr)) -= std::get<0>(std::get<2>(ARR));
   std::get<1>(std::get<0>(arr)) -= std::get<1>(std::get<0>(ARR));
   std::get<1>(std::get<1>(arr)) -= std::get<1>(std::get<1>(ARR));
   std::get<1>(std::get<2>(arr)) -= std::get<1>(std::get<2>(ARR));
   std::get<2>(std::get<0>(arr)) -= std::get<2>(std::get<0>(ARR));
   std::get<2>(std::get<1>(arr)) -= std::get<2>(std::get<1>(ARR));
   std::get<2>(std::get<2>(arr)) -= std::get<2>(std::get<2>(ARR));
   return arr;
}
/* -------------------------------------------------------------------------- */
// 1d array - value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator-(std::array<T, N> arr /*copy*/, const TT& d) noexcept {
   for (auto& a : arr) a -= d;
   return arr;
}
// value - 1d array
template <size_t N, typename T>
constexpr std::enable_if_t<!is_std_array<T>::value, std::array<T, N>> operator-(const T& d, std::array<T, N> arr /*copy*/) noexcept {
   for (auto& a : arr)
      a = d - a;
   return arr;
}
// 1d array - 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N> operator-(std::array<T, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr -= ARR; }

// 2d array - 1d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator-(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr -= ARR; }

// 1d array - 2d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator-(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept {
   for (auto i = 0; auto& a : arr)
      a = ARR[i++] - a;
   return arr;
}

template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<T>::value && !is_std_array<TT>::value, std::array<std::array<TT, M>, N>> operator-(const std::array<std::array<TT, M>, N>& ARR, const std::array<std::array<T, M>, N>& arr) noexcept {
   std::array<std::array<TT, M>, N> result{};
   for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
         result[i][j] = ARR[i][j] - arr[i][j];
      }
   }
   return result;
}

template <size_t N0, std::size_t N1, std::size_t N2, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<T>::value && !is_std_array<TT>::value, std::array<std::array<std::array<TT, N2>, N1>, N0>> operator-(const std::array<std::array<std::array<TT, N2>, N1>, N0>& ARR, const std::array<std::array<std::array<T, N2>, N1>, N0>& arr) noexcept {
   std::array<std::array<std::array<TT, N2>, N1>, N0> result{};
   for (size_t i = 0; i < N0; ++i) {
      for (size_t j = 0; j < N1; ++j) {
         for (size_t k = 0; k < N2; ++k) {
            result[i][j][k] = ARR[i][j][k] - arr[i][j][k];
         }
      }
   }
   return result;
}

/* ----------------------------- multiplication ----------------------------- */
// 1d array *= value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>&> operator*=(std::array<T, N>& arr /*ref*/, const TT& d) noexcept {
   for (auto& a : arr) a *= d;
   return arr;
}
// 1d array *= 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N>& operator*=(std::array<T, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) *= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}
// 2d array *= 1d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator*=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) *= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}
// 2d array *= 2d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator*=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) *= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}

// 1d array * value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator*(std::array<T, N> arr /*copy*/, const TT& d) noexcept {
   for (auto& a : arr) a *= d;
   return arr;
}
// value * 1d array
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator*(const TT& d, std::array<T, N> arr /*copy*/) noexcept {
   for (auto& a : arr) a *= d;
   return arr;
}
// 1d array * 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N> operator*(std::array<T, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr *= ARR; }
// 2d array * 1d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator*(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr *= ARR; }
// 1d array * 2d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N> operator*(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept {
   for (auto i = 0; auto& a : arr) a *= ARR[i++];
   return arr;
}

// 2d array * 2d array
template <size_t N, std::size_t M>
constexpr std::array<std::array<double, M>, N> operator*(const std::array<std::array<double, M>, N>& ARR, std::array<std::array<double, M>, N> arr /*copy*/) noexcept {
   for (auto i = 0; i < N; ++i) {
      for (auto j = 0; j < M; ++j) {
         arr[i][j] *= ARR[i][j];
      }
   }
   return arr;
}

/* -------------------------------- division -------------------------------- */
// 1d array /= value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>&> operator/=(std::array<T, N>& arr /*ref*/, const TT& d) noexcept {
   for (auto& a : arr) a /= d;
   return arr;
}
// 1d array /= 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N>& operator/=(std::array<T, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
// 2d array /= 1d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator/=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
// 2d array /= 2d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator/=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
/* -------------------------------------------------------------------------- */
// 1d array / value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator/(std::array<T, N> arr /*copy*/, const TT& d) noexcept {
   for (auto& a : arr) a /= d;
   return arr;
}
// value / 1d array
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator/(const TT& d, std::array<T, N> arr /*copy*/) noexcept {
   for (auto& a : arr) a = d / a;
   return arr;
}
// 1d array / 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N> operator/(std::array<T, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
// 2d array / 1d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator/(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++]; /*std::array<T, M> /= TT*/
   return arr;
}
// 1d array / 2d array
template <size_t N, std::size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N> operator/(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept {
   for (auto i = 0; auto& a : arr) a = ARR[i++] / a;
   return arr;
}
/* -------------------------------------------------------------------------- */
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& stream, const std::array<T, N>& arr) noexcept {
   bool first = true;
   stream << "{";
   std::ranges::for_each(arr, [&](const auto& x) { stream << (first ? (first = false, "") : ",") << x; });
   stream << "}";
   return stream;
}
/* -------------------------------------------------------------------------- */

template <size_t N, typename T, typename TT, typename Func>
constexpr void for_each(std::array<T, N>& arr, std::array<TT, N>& ARR, const Func& func) {
   if constexpr (N > 0) {
      [&]<size_t... Is>(std::index_sequence<Is...>) {
         (func(std::get<Is>(arr), std::get<Is>(ARR)), ...);
      }(std::make_index_sequence<N>());
   }
}

template <size_t N, typename T, typename TT, typename Func>
constexpr void for_each(std::array<T, N>& arr, const std::array<TT, N>& ARR, const Func& func) {
   if constexpr (N > 0) {
      [&]<size_t... Is>(std::index_sequence<Is...>) {
         (func(std::get<Is>(arr), std::get<Is>(ARR)), ...);
      }(std::make_index_sequence<N>());
   }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
std::vector<Tddd> operator-(std::vector<Tddd> V, const Tddd& u) noexcept {
   for (auto& v : V) v -= u;
   return V;
};
std::vector<Tddd> operator+(std::vector<Tddd> V, const Tddd& u) noexcept {
   for (auto& v : V) v += u;
   return V;
};
std::vector<Tddd> operator*(std::vector<Tddd> V, const Tddd& u) noexcept {
   for (auto& v : V) v *= u;
   return V;
};
std::vector<Tddd> operator/(std::vector<Tddd> V, const Tddd& u) noexcept {
   for (auto& v : V) v /= u;
   return V;
};

/* =========================================================================== */
/*                       Mathematical Vector Operations                       */
/* =========================================================================== */
template <size_t N0, std::size_t N1, typename T>
constexpr std::array<std::array<T, N0>, N1> Transpose(const std::array<std::array<T, N1>, N0>& M) noexcept {
   std::array<std::array<T, N0>, N1> ret;
   std::size_t i = 0, j = 0;
   for (i = 0; const auto& Mi : M) {
      for (j = 0; const auto& Mij : Mi)
         ret[j++][i] = Mij;
      i++;
   }
   return ret;
}

template <size_t N, typename T>
constexpr std::array<std::vector<T>, N> Transpose(const std::vector<std::array<T, N>>& M) noexcept {
   std::array<std::vector<T>, N> ret;
   std::size_t j = 0;
   for (j = 0; j < N; ++j)
      ret[j].reserve(M.size());
   for (const auto& Mi : M) {
      j = 0;
      for (const auto& Mij : Mi) {
         ret[j++].push_back(Mij);
      }
   }
   return ret;
}

template <typename T, std::size_t N>
constexpr std::array<T, 2> MinMax(const std::array<T, N>& arr) noexcept {
   auto minmax = std::minmax_element(arr.begin(), arr.end());
   return {*minmax.first, *minmax.second};
}
constexpr Tdd MinMax(const Tdd& A) noexcept { return ((A[0] < A[1]) ? Tdd{A[0], A[1]} : Tdd{A[1], A[0]}); };

template <typename T, std::size_t N0, std::size_t N1>
constexpr auto MinMaxColumns(const std::array<std::array<T, N0>, N1>& arr) {
   std::array<std::array<T, 2>, N0> minmax_array;
   std::size_t i = 0;
   for (i = 0; i < N0; ++i)
      minmax_array[i] = {arr[0][i], arr[0][i]};

   for (const auto& row : arr) {
      for (i = 0; const auto& a : row) {
         if (minmax_array[i][0] >= a)
            minmax_array[i][0] = a;
         else if (minmax_array[i][1] <= a)
            minmax_array[i][1] = a;
         i++;
      }
   }

   return minmax_array;
}

constexpr T3Tdd MinMaxColumns(const T2Tddd& A) noexcept {
   return {MinMax(Tdd{A[0][0], A[1][0]}), MinMax(Tdd{A[0][1], A[1][1]}), MinMax(Tdd{A[0][2], A[1][2]})};
};

template <typename T, std::size_t N>
constexpr auto MinMaxColumns(const std::vector<std::array<T, N>>& vec_arr) noexcept {
   std::array<std::array<T, 2>, N> minmax_array;
   std::size_t i = 0;
   for (i = 0; i < N; ++i)
      minmax_array[i] = {vec_arr[0][i], vec_arr[0][i]};
   for (const auto& arr : vec_arr) {
      for (i = 0; const auto& a : arr) {
         if (minmax_array[i][0] >= a)
            minmax_array[i][0] = a;
         else if (minmax_array[i][1] <= a)
            minmax_array[i][1] = a;
         i++;
      }
   }
   return minmax_array;
}

template <size_t N, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, T>::type
Total(const std::array<T, N>& arr) noexcept {
   T ret = 0;
   for (const auto& a : arr) ret += a;
   return ret;
}

template <typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, T>::type
Total(const std::array<T, 2>& arr) noexcept {
   return std::get<0>(arr) + std::get<1>(arr);
}

template <typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, T>::type
Total(const std::array<T, 3>& arr) noexcept {
   return std::get<0>(arr) + std::get<1>(arr) + std::get<2>(arr);
}

template <size_t N0, std::size_t N1, typename T>
constexpr std::array<T, N1> Total(const std::array<std::array<T, N1>, N0>& arr) noexcept {
   std::array<T, N1> ret{};
   std::size_t i = 0, j = 0;
   for (i = 0; i < N0; ++i)
      ret += arr[i];
   return ret;
}

template <std::size_t N1, typename T>
constexpr std::array<T, N1> Total(const std::array<std::array<T, N1>, 3>& arr) noexcept {
   return std::get<0>(arr) + std::get<1>(arr) + std::get<2>(arr);
}

/* -------------------------------------------------------------------------- */

template <typename T, std::size_t N1, std::size_t N2>
constexpr std::array<std::array<T, N2>, N1> TensorProduct(const T scalar, const std::array<std::array<T, N2>, N1>& matrix) noexcept {
   std::array<std::array<T, N2>, N1> result{};
   std::size_t i = 0, j = 0;
   for (i = 0; i < N1; ++i) {
      for (j = 0; j < N2; ++j) {
         result[i][j] = scalar * matrix[i][j];
      }
   }
   return result;
}

template <typename T, std::size_t N1, std::size_t N2>
constexpr std::array<std::array<T, N2>, N1> TensorProduct(const std::array<T, N1>& vec1, const std::array<T, N2>& vec2) noexcept {
   std::array<std::array<T, N2>, N1> ret{};
   std::size_t m = 0, j = 0;
   for (m = 0; m < N1; ++m) {
      for (j = 0; j < N2; ++j) {
         ret[m][j] = vec1[m] * vec2[j];
      }
   }
   return ret;
}

// special case for 3x3 matrix
template <typename T>
constexpr std::array<std::array<T, 3>, 3> TensorProduct(const std::array<T, 3>& vec1, const std::array<T, 3>& vec2) noexcept {
   return {{{std::get<0>(vec1) * std::get<0>(vec2), std::get<0>(vec1) * std::get<1>(vec2), std::get<0>(vec1) * std::get<2>(vec2)},
            {std::get<1>(vec1) * std::get<0>(vec2), std::get<1>(vec1) * std::get<1>(vec2), std::get<1>(vec1) * std::get<2>(vec2)},
            {std::get<2>(vec1) * std::get<0>(vec2), std::get<2>(vec1) * std::get<1>(vec2), std::get<2>(vec1) * std::get<2>(vec2)}}};
}

// Base case for scalar and vector
template <typename T, std::size_t N>
constexpr auto TensorProduct(const T scalar, const std::array<T, N>& vec) noexcept {
   std::array<T, N> result{};
   for (size_t i = 0; i < N; ++i) {
      result[i] = scalar * vec[i];
   }
   return result;
}

// Recursive case for tensors
template <typename T, std::size_t N, typename... Arrays>
constexpr auto TensorProduct(const std::array<T, N>& tensor1, const Arrays&... tensors) noexcept {
   std::array<decltype(TensorProduct(tensor1[0], tensors...)), N> result{};
   for (size_t i = 0; i < N; ++i)
      result[i] = TensorProduct(tensor1[i], tensors...);
   return result;
}

/* -------------------------------------------------------------------------- */

template <typename T, std::size_t N>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, T>::type
Dot(const std::array<T, N>& arr, const std::array<T, N>& ARR) noexcept {
   T ret = 0;
   std::size_t i = 0;
   for (const auto& a : arr)
      ret = std::fma(a, ARR[i++], ret);
   return ret;
}

double Dot(const std::array<double, 3>& arr, const std::array<double, 3>& ARR) noexcept {
   return {std::fma(std::get<0>(arr), std::get<0>(ARR), std::fma(std::get<1>(arr), std::get<1>(ARR), std::get<2>(arr) * std::get<2>(ARR)))};
}

template <typename T, std::size_t N0, std::size_t N1, std::size_t N2>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<T, N2>, N1>>::type
Dot(const std::array<std::array<T, N0>, N1>& arr, const std::array<std::array<T, N2>, N0>& ARR) noexcept {
   std::array<std::array<T, N2>, N1> ret{};
   T sum = 0;
   std::size_t i = 0, j = 0, k = 0;
   for (i = 0; i < N1; ++i) {
      for (j = 0; j < N2; ++j) {
         sum = 0;
         for (k = 0; k < N0; ++k) {
            // sum += arr[i][k] * ARR[k][j];
            sum = std::fma(arr[i][k], ARR[k][j], sum);
         }
         ret[i][j] = sum;
      }
   }
   return ret;
}

template <size_t N, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, void>::type
Fill(std::array<T, N>& arr, const T d) noexcept {
   for (auto& a : arr)
      a = d;
}

template <size_t N0, std::size_t N1, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, void>::type
Fill(std::array<std::array<T, N1>, N0>& arr, const T d) noexcept {
   for (auto& a : arr)
      Fill(a, d);
}

template <size_t N0, std::size_t N1, std::size_t N2, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, void>::type
Fill(std::array<std::array<std::array<T, N2>, N1>, N0>& arr, const T d) noexcept {
   for (auto& a : arr)
      Fill(a, d);
}

template <size_t N0, std::size_t N1, std::size_t N2, std::size_t N3, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, void>::type
Fill(std::array<std::array<std::array<std::array<T, N3>, N2>, N1>, N0>& arr, const T d) noexcept {
   for (auto& a : arr)
      Fill(a, d);
}

/* -------------------------------------------------------------------------- */
template <size_t N0, std::size_t N1, typename T>
constexpr std::array<T, N1> Dot(const std::array<T, N0>& arr, const std::array<std::array<T, N1>, N0>& ARRARR) noexcept {
   /*
   {a,b,c}.{{A0,A1,A2,A3},{B0,B1,B2,B3},{C0,C1,C2,C3}}
   */
   std::array<T, N1> ret{};
   // for_each(ARRARR, arr, [&](auto &ARR, auto &a) {
   //    for (size_t j = 0; j < N1; ++j) {
   //       ret[j] += a * ARR[j];
   //    }
   // });
   //
   std::size_t i = 0, j = 0;
   for (i = 0; const auto& ARR : ARRARR) {
      for (j = 0; j < N1; ++j)
         ret[j] = std::fma(arr[i], ARR[j], ret[j]);
      // ret[j] += arr[i] * ARR[j];
      i++;
   }

   return ret;
}

template <size_t N0, std::size_t N1, std::size_t N2, typename T>
constexpr std::array<std::array<T, N2>, N1> Dot(const std::array<T, N0>& arr, const std::array<std::array<std::array<T, N2>, N1>, N0>& ARRARR) noexcept {
   std::array<std::array<T, N2>, N1> ret = {};
   std::size_t i = 0, j = 0;
   for (i = 0; const auto& ARR : ARRARR) {
      for (j = 0; j < N1; ++j)
         ret[j] += arr[i] * ARR[j];
      // ret[j] = std::fma(arr[i], ARR[j], ret[j]);
      i++;
   }
   return ret;
}

template <size_t N0, std::size_t N1, std::size_t N2, std::size_t N3, typename T>
constexpr std::array<std::array<std::array<T, N3>, N2>, N1> Dot(const std::array<T, N0>& arr, const std::array<std::array<std::array<std::array<T, N3>, N2>, N1>, N0>& ARRARR) noexcept {
   std::array<std::array<std::array<T, N3>, N2>, N1> ret = {};
   std::size_t i = 0, j = 0;
   for (i = 0; const auto& ARR : ARRARR) {
      for (j = 0; j < N1; ++j)
         ret[j] += arr[i] * ARR[j];
      // ret[j] = std::fma(arr[i], ARR[j], ret[j]);
      i++;
   }
   return ret;
}

template <size_t N0, std::size_t N1, std::size_t N2, std::size_t N3, std::size_t N4, typename T>
constexpr std::array<std::array<std::array<std::array<T, N4>, N3>, N2>, N1> Dot(const std::array<T, N0>& arr, const std::array<std::array<std::array<std::array<std::array<T, N4>, N3>, N2>, N1>, N0>& ARRARR) noexcept {
   std::array<std::array<std::array<std::array<T, N4>, N3>, N2>, N1> ret = {};
   std::size_t i = 0, j = 0;
   for (i = 0; const auto& ARR : ARRARR) {
      for (j = 0; j < N1; ++j)
         ret[j] += arr[i] * ARR[j];
      // ret[j] = std::fma(arr[i], ARR[j], ret[j]);
      i++;
   }
   return ret;
}

/* -------------------------------------------------------------------------- */

template <size_t N0, std::size_t N1, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<T, N0>>::type
Dot(const std::array<std::array<T, N1>, N0>& AA, const std::array<T, N1>& arr) noexcept {
   std::array<T, N0> ret = {};
   std::size_t i = 0;
   for (const auto& A : AA)
      ret[i++] = Dot(A, arr);
   return ret;
}

template <size_t N0, std::size_t N1, std::size_t N2, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<T, N1>, N0>>::type
Dot(const std::array<std::array<std::array<T, N2>, N1>, N0>& AAA, const std::array<T, N2>& arr) noexcept {
   std::array<std::array<T, N1>, N0> ret = {};
   std::size_t i = 0, j = 0;
   for (i = 0; const auto& AA : AAA) {
      for (j = 0; const auto& A : AA) {
         ret[i][j++] = Dot(A, arr);
      }
      i++;
   }
   return ret;
}

template <size_t N0, std::size_t N1, std::size_t N2, std::size_t N3, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<std::array<T, N3>, N1>, N0>>::type
Dot(const std::array<std::array<std::array<T, N2>, N1>, N0>& AAA, const std::array<std::array<T, N3>, N2>& arr) noexcept {
   std::array<std::array<std::array<T, N3>, N1>, N0> ret = {};
   int i = 0, j = 0;
   for (i = 0; const auto& AA : AAA) {
      for (j = 0; const auto& A : AA) {
         ret[i][j++] = Dot(A, arr);
      }
      i++;
   }
   return ret;
}

template <size_t N0, std::size_t N1, std::size_t N2, std::size_t N3, std::size_t N4, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<std::array<std::array<T, N4>, N3>, N1>, N0>>::type
Dot(const std::array<std::array<std::array<T, N2>, N1>, N0>& AAA, const std::array<std::array<std::array<T, N4>, N3>, N2>& arr) noexcept {
   std::array<std::array<std::array<std::array<T, N4>, N3>, N1>, N0> ret = {};
   std::size_t i = 0, j = 0;
   for (i = 0; const auto& AA : AAA) {
      for (j = 0; const auto& A : AA) {
         ret[i][j++] = Dot(A, arr);
      }
      i++;
   }
   return ret;
}

template <size_t N0, std::size_t N1, std::size_t N2, std::size_t N3, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<std::array<T, N2>, N1>, N0>>::type
Dot(const std::array<std::array<std::array<std::array<T, N3>, N2>, N1>, N0>& AAAA, const std::array<T, N3>& arr) noexcept {
   std::array<std::array<std::array<T, N2>, N1>, N0> ret = {};
   std::size_t i = 0, j = 0, k = 0;
   for (i = 0; const auto& AAA : AAAA) {
      for (j = 0; const auto& AA : AAA) {
         for (k = 0; const auto& A : AA) {
            ret[i][j][k++] = Dot(A, arr);
         }
         j++;
      }
      i++;
   }
   return ret;
}

/* -------------------------------------------------------------------------- */
// template <size_t N0, std::size_t N1, typename T>
// constexpr std::array<T, N1> Dot(const std::array<std::array<T, N0>, N1>& arr, const std::array<T, N0>& ARR) noexcept {
//    /*
//    {{A0,A1,A2},{B0,B1,B2},{C0,C1,C2},{D0,D1,D2}}.{a,b,c}={a*A0+b*A1+c*A2,a*B0+b*B1+c*B2,a*C0+b*C1+c*C2,a*D0+b*D1+c*D2}
//    */
//    std::array<T, N1> ret;
//    for_each(ret, arr, [&](auto& r, auto& ar) { r = Dot(ar, ARR); });
//    return ret;
// }

template <typename... Ts, typename... Is, std::size_t... Indices>
constexpr void ToArrayImpl(const std::tuple<Ts...>& tup, std::array<std::common_type_t<Ts...>, sizeof...(Ts)>& arr, std::index_sequence<Indices...>) noexcept {
   ((arr[Indices] = std::get<Indices>(tup)), ...);
}

template <typename... Ts>
constexpr auto ToArray(const std::tuple<Ts...>& tup) {
   using T = std::common_type_t<Ts...>;
   constexpr std::size_t N = sizeof...(Ts);
   std::array<T, N> arr{};
   ToArrayImpl(tup, arr, std::make_index_sequence<N>());
   return arr;
}

template <typename T, std::size_t N>
constexpr std::array<T, N> ToArray(const std::array<T, N>& arr) noexcept { return arr; }

// constexpr std::array<T, N1> Dot(const std::array<T, N0> &arr, const std::array<std::array<T, N1>, N0> &ARR) {
//    std::array<T, N1> ret = {};
//    for_each(arr, [&](const auto &ar) {
//       for_each(ARR, [&](const auto &AR) {
//          for_each(ret, AR, [&](auto &r, auto &A) {
//             r += ar * A;
//          });
//       });
//    });

//    return ret;
// }
//
//

template <size_t N, typename T>
constexpr T NormSquared(const std::array<T, N>& arr) noexcept {
   T sum = {};
   for (const auto& a : arr)
      sum = std::fma(a, a, sum);
   return sum;
}

template <size_t N>
constexpr double NormSquared(const std::array<double, N>& arr) noexcept {
   double sum = 0.;
   for (const auto& a : arr)
      sum = std::fma(a, a, sum);
   return sum;
}

template <size_t N, typename T>
constexpr T Norm(const std::array<T, N>& arr) noexcept {
   return std::sqrt(Dot(arr, arr));
}

template <size_t N, std::size_t M, typename T>
constexpr auto Norm(const std::array<std::array<T, M>, N>& mat) noexcept {
   //$ Mathematica implementation returns the maximum of the singular values of the matrix, which is not the same as the Frobenius norm implemented here.
   static_assert(!std::is_array_v<T>, "Norm is not defined for arrays deeper than 2 dimensions.");
   T sum = {};
   for (const auto& row : mat)
      sum += NormSquared(row);
   return std::sqrt(sum);
}

template <size_t N>
constexpr std::array<double, N> Normalize(std::array<double, N> arr) noexcept {
   double norm = 0.;
   for (const auto& a : arr)
      norm = std::fma(a, a, norm);

   // Check if norm is zero to avoid division by zero
   if (norm == 0.)
      return arr;  // or return {}; for an array of zeros

   norm = std::sqrt(norm);

   for (auto& a : arr)
      a /= norm;
   return arr;
   // static_assert(N > 0, "Array must have at least one element.");
   // T norm = Norm(arr);
   // auto c = (norm == static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(1) / norm;
   // std::ranges::transform(arr, arr.begin(), [&c](T a) { return a * c; });
   // return arr;
}

double Norm(const std::array<double, 2>& arr) noexcept {
   return std::hypot(std::get<0>(arr), std::get<1>(arr));
}

double Norm(const std::array<double, 3>& arr) noexcept {
   return std::hypot(std::get<0>(arr), std::get<1>(arr), std::get<2>(arr));
}

constexpr double Norm(const double value) noexcept {
   return std::abs(value);
}

double Distance(const Tddd& A, const Tddd& B) {
   return Norm(A - B);
};

template <size_t M = 0, std::size_t N, typename T>
constexpr T RootMeanSquare(const std::array<T, N>& arr) noexcept { return std::sqrt(Dot(arr, arr) / N); }

// template <typename T>
// constexpr std::array<T, 3> Cross(const std::array<T, 3>& A, const std::array<T, 3>& B) noexcept {
//    return {{std::fma(-std::get<2>(A), std::get<1>(B), std::get<1>(A) * std::get<2>(B)),
//             std::fma(-std::get<0>(A), std::get<2>(B), std::get<2>(A) * std::get<0>(B)),
//             std::fma(-std::get<1>(A), std::get<0>(B), std::get<0>(A) * std::get<1>(B))}};
// }

#include <type_traits>

// template <typename T>
// inline constexpr std::array<T, 3> Cross(const std::array<T, 3>& A, const std::array<T, 3>& B) noexcept {
//    static_assert(std::is_arithmetic_v<T>, "Arithmetic type required.");
//    return {{std::fma(-std::get<2>(A), std::get<1>(B), std::get<1>(A) * std::get<2>(B)),
//             std::fma(-std::get<0>(A), std::get<2>(B), std::get<2>(A) * std::get<0>(B)),
//             std::fma(-std::get<1>(A), std::get<0>(B), std::get<0>(A) * std::get<1>(B))}};
// }

template <typename T>
inline constexpr std::array<T, 3> Cross(const std::array<T, 3>& A, const std::array<T, 3>& B) noexcept {
   static_assert(std::is_arithmetic_v<T>, "Arithmetic type required.");
   return {{std::fma(std::get<1>(A), std::get<2>(B), -std::get<2>(A) * std::get<1>(B)),
            std::fma(std::get<2>(A), std::get<0>(B), -std::get<0>(A) * std::get<2>(B)),
            std::fma(std::get<0>(A), std::get<1>(B), -std::get<1>(A) * std::get<0>(B))}};
}

/* -------------------------------------------------------------------------- */

// void FusedMultiplyIncrement(std::array<double, 3>& ret, const double arr, const std::array<double, 3>& ARR) noexcept {
//    std::get<0>(ret) = std::fma(arr, std::get<0>(ARR), std::get<0>(ret));
//    std::get<1>(ret) = std::fma(arr, std::get<1>(ARR), std::get<1>(ret));
//    std::get<2>(ret) = std::fma(arr, std::get<2>(ARR), std::get<2>(ret));
// }

template <size_t N>
std::array<double, N> FusedMultiplyAdd(const double W, const std::array<double, N>& ARR, std::array<double, N> ret) noexcept {
   for (size_t i = 0; i < N; ++i)
      ret[i] = std::fma(W, ARR[i], ret[i]);
   return ret;
}

std::array<double, 3> FusedMultiplyAdd(const double W, const std::array<double, 3>& ARR, const std::array<double, 3>& ret) noexcept {
   return {std::fma(W, std::get<0>(ARR), std::get<0>(ret)),
           std::fma(W, std::get<1>(ARR), std::get<1>(ret)),
           std::fma(W, std::get<2>(ARR), std::get<2>(ret))};
}

template <size_t N>
void FusedMultiplyIncrement(const double W, const std::array<double, N>& ARR, std::array<double, N>& ret) noexcept {
   for (size_t i = 0; i < N; ++i)
      ret[i] = std::fma(W, ARR[i], ret[i]);
}

void FusedMultiplyIncrement(const double W, const std::array<double, 3>& ARR, std::array<double, 3>& ret) noexcept {
   std::get<0>(ret) = std::fma(W, std::get<0>(ARR), std::get<0>(ret));
   std::get<1>(ret) = std::fma(W, std::get<1>(ARR), std::get<1>(ret));
   std::get<2>(ret) = std::fma(W, std::get<2>(ARR), std::get<2>(ret));
}

template <size_t N, size_t N1>
void FusedMultiplyIncrement(const double W, const std::array<std::array<double, N1>, N>& ARR, std::array<std::array<double, N1>, N>& ret) noexcept {
   for (size_t i = 0; i < N1; ++i)
      FusedMultiplyIncrement(W, ARR[i], ret[i]);
}

// void FusedMultiplyIncrement(const double W, const std::vector<double>& ARR, std::vector<double>& ret) noexcept {
//    for (std::size_t i = 0; auto& r : ret)
//       r = std::fma(W, ARR[i++], r);
// }

#include <execution>

void FusedMultiplyIncrement(const double W, const std::vector<double>& ARR, std::vector<double>& ret) noexcept {
   std::transform(std::execution::unseq, ARR.begin(), ARR.end(), ret.begin(), ret.begin(),
                  [W](const double& a, double& r) {
                     return std::fma(a, W, r);
                  });
}

void FusedMultiplyIncrement(const double W, const std::vector<float>& ARR, std::vector<float>& ret) noexcept {
   std::transform(std::execution::unseq, ARR.begin(), ARR.end(), ret.begin(), ret.begin(),
                  [W](const float& a, float& r) {
                     return std::fma(a, W, r);
                  });
}

// void FusedMultiplyIncrement(const std::vector<double>& ARR, const double W, std::vector<double>& ret) noexcept {
//    for (std::size_t i = 0; auto& r : ret)
//       r = std::fma(W, ARR[i++], r);
// }

void FusedMultiplyIncrement(const std::vector<double>& ARR, const double W, std::vector<double>& ret) noexcept {
   std::transform(std::execution::unseq, ARR.begin(), ARR.end(), ret.begin(), ret.begin(),
                  [W](const double& a, double& r) {
                     return std::fma(a, W, r);
                  });
}

void FusedMultiplyIncrement(const std::vector<float>& ARR, const float W, std::vector<float>& ret) noexcept {
   std::transform(std::execution::unseq, ARR.begin(), ARR.end(), ret.begin(), ret.begin(),
                  [W](const float& a, float& r) {
                     return std::fma(a, W, r);
                  });
}

void FusedMultiplyIncrement(const std::vector<double>& ARR, const float W, std::vector<float>& ret) noexcept {
   std::transform(std::execution::unseq, ARR.begin(), ARR.end(), ret.begin(), ret.begin(),
                  [W](const float& a, float& r) {
                     return std::fma(a, W, r);
                  });
}

void FusedMultiplyIncrement(const std::vector<double>& ARR, const float W, std::vector<double>& ret) noexcept {
   std::transform(std::execution::unseq, ARR.begin(), ARR.end(), ret.begin(), ret.begin(),
                  [W](const double& a, double& r) {
                     return std::fma(a, W, r);
                  });
}

void FusedMultiplyIncrement(const std::vector<float>& ARR, const float W, std::vector<double>& ret) noexcept {
   std::transform(std::execution::unseq, ARR.begin(), ARR.end(), ret.begin(), ret.begin(),
                  [W](const double& a, double& r) {
                     return std::fma(a, W, r);
                  });
}

template <size_t N>
void FusedMultiplyIncrement(const std::array<double, N>& ARR, const double W, std::array<double, N>& ret) noexcept {
   for (std::size_t i = 0; i < N; ++i)
      ret[i] = std::fma(W, ARR[i], ret[i]);
}

// template <size_t N>
// void FusedMultiplyIncrement(const std::array<double, N>& ARR, const double W, std::array<double, N>& ret) noexcept {
//    if constexpr (N > 0) {
//       ret[N - 1] = std::fma(W, ARR[N - 1], ret[N - 1]);
//       FusedMultiplyIncrement<N - 1>(ARR, W, ret);
//    }
// }

template <>
void FusedMultiplyIncrement<0>(const std::array<double, 0>&, const double, std::array<double, 0>&) noexcept {
   // Do nothing
}

void FusedMultiplyIncrement(const double arr, const double ARR, double& ret) noexcept {
   ret = std::fma(arr, ARR, ret);
}

/* -------------------------------------------------------------------------- */

// double Det(const T3Tddd& M) {
//    return -(std::get<0>(std::get<2>(M)) * std::get<1>(std::get<1>(M)) * std::get<2>(std::get<0>(M))) +
//           std::get<0>(std::get<1>(M)) * std::get<1>(std::get<2>(M)) * std::get<2>(std::get<0>(M)) +
//           std::get<0>(std::get<2>(M)) * std::get<1>(std::get<0>(M)) * std::get<2>(std::get<1>(M)) -
//           std::get<0>(std::get<0>(M)) * std::get<1>(std::get<2>(M)) * std::get<2>(std::get<1>(M)) -
//           std::get<0>(std::get<1>(M)) * std::get<1>(std::get<0>(M)) * std::get<2>(std::get<2>(M)) +
//           std::get<0>(std::get<0>(M)) * std::get<1>(std::get<1>(M)) * std::get<2>(std::get<2>(M));
// };

double Det(const T3Tddd& M) {
   // // Unpacking each row of the matrix
   // auto [a, b, c] = std::get<0>(M);  // First row
   // auto [d, e, f] = std::get<1>(M);  // Second row
   // auto [g, h, i] = std::get<2>(M);  // Third row

   // // Calculating the determinant using std::fma
   // return std::fma(a, std::fma(e, i, -f * h),
   //                 std::fma(b, std::fma(f, g, -d * i),
   //                          std::fma(c, std::fma(d, h, -e * g), 0.)));

   return M[0][0] * (std::fma(M[1][1], M[2][2], -M[1][2] * M[2][1])) -
          M[0][1] * (std::fma(M[1][0], M[2][2], -M[1][2] * M[2][0])) +
          M[0][2] * (std::fma(M[1][0], M[2][1], -M[1][1] * M[2][0]));
}

/* -------------------------------------------------------------------------- */

template <std::size_t N>
constexpr std::array<double, N + 1> Subdivide(const double xmin, const double xmax) noexcept {
   static_assert(N > 0, "The number of divisions N must be a positive integer.");
   std::array<double, N + 1> ret;
   const double dx = (xmax - xmin) / N;
   for (std::size_t i = 0; i < N + 1; i++)
      ret[i] = i * dx + xmin;
   return ret;
};

/* -------------------------------------------------------------------------- */

// struct SubdivideDomainIntoTriangles {
//    const int division;
//    const double max;
//    const double min;
//    std::vector<std::array<int, 3>> parameterOnTriangle_index, parameterOnSquare_index;
//    std::vector<std::array<std::array<double, 2>, 3>> parameterOnTriangle, parameterOnSquare;
//    SubdivideDomainIntoTriangles(const int division, const double min = 0., const double max = 1.) : division(division), max(max), min(min) {
//       if (division < 1)
//          throw std::invalid_argument("The number of divisions must be a positive integer.");

//       if (min >= max)
//          throw std::invalid_argument("The minimum value must be less than the maximum value.");
//       //
//       std::vector<double> vec(division + 1);
//       const double dx = (max - min) / (double)division;
//       for (std::size_t i = 0; i < division + 1; i++)
//          vec[i] = i * dx + min;
//       double dt = 1. / (double)division / 2.;
//       //
//       std::vector<std::array<std::array<int, 2>, 3>> tmp_index;
//       int index = 0;
//       for (auto i = 0; i < vec.size() - 1; ++i)
//          for (auto j = 0; j < vec.size(); ++j) {
//             parameterOnTriangle.push_back({{{vec[i], vec[j]}, {vec[i + 1], vec[j]}, {vec[i], vec[j + 1]}}});
//             tmp_index.push_back({{{i, j}, {i + 1, j}, {i, j + 1}}});
//             if (std::abs(1 - vec[i + 1] - vec[j]) < dt && std::abs(1 - vec[i] - vec[j + 1]) < dt)
//                break;
//             parameterOnTriangle.push_back({{{vec[i], vec[j + 1]}, {vec[i + 1], vec[j]}, {vec[i + 1], vec[j + 1]}}});
//             tmp_index.push_back({{{i, j + 1}, {i + 1, j}, {i + 1, j + 1}}});
//          }

//       auto convert = [&](Tii ij) { return ij[0] + (division + 1) * ij[1]; };

//       for (auto i = 0; i < tmp_index.size(); ++i) {
//          auto [ij0, ij1, ij2] = tmp_index[i];
//          parameterOnTriangle_index.push_back({convert(ij0), convert(ij1), convert(ij2)});
//       }

//       tmp_index.clear();
//       for (auto i = 0; i < vec.size() - 1; ++i)
//          for (auto j = 0; j < vec.size() - 1; ++j) {
//             parameterOnSquare.push_back({{{vec[i], vec[j]}, {vec[i + 1], vec[j]}, {vec[i], vec[j + 1]}}});
//             tmp_index.push_back({{{i, j}, {i + 1, j}, {i, j + 1}}});
//             parameterOnSquare.push_back({{{vec[i], vec[j + 1]}, {vec[i + 1], vec[j]}, {vec[i + 1], vec[j + 1]}}});
//             tmp_index.push_back({{{i, j + 1}, {i + 1, j}, {i + 1, j + 1}}});
//          }

//       for (auto i = 0; i < tmp_index.size(); ++i) {
//          auto [ij0, ij1, ij2] = tmp_index[i];
//          parameterOnSquare_index.push_back({convert(ij0), convert(ij1), convert(ij2)});
//       }

//       std::cout << "parameterOnTriangle.size() = " << parameterOnTriangle.size() << std::endl;
//       std::cout << "parameterOnTriangle_index.size() = " << parameterOnTriangle_index.size() << std::endl;
//       std::cout << "parameterOnSquare.size() = " << parameterOnSquare.size() << std::endl;
//       std::cout << "parameterOnSquare_index.size() = " << parameterOnSquare_index.size() << std::endl;
//    };
// };

// std::vector<std::array<std::array<double, 2>, 3>> SubdivideTriangleIntoTriangles(const int divide) {
std::vector<std::array<std::array<double, 2>, 3>> SymmetricSubdivisionOfTriangle_00_10_01(const int divide) {
   std::vector<double> vec(divide + 1);
   const double max = 1., min = 0.;
   const double dt = (max - min) / (double)divide;
   //
   for (std::size_t i = 0; i < divide + 1; i++)
      vec[i] = i * dt + min;
   //
   double resolution = 0.5 / (double)divide;
   //
   std::vector<std::array<std::array<double, 2>, 3>> t0t1_triangles;
   t0t1_triangles.reserve(divide * divide);
   //
   for (auto i = 0; i < vec.size() - 1; ++i)
      for (auto j = 0; j < vec.size() - 1; ++j) {
         t0t1_triangles.push_back({{{vec[i], vec[j]}, {vec[i + 1], vec[j]}, {vec[i], vec[j + 1]}}});

         //@
         //@ |\
         //@ | \
         //@ ----

         if (std::abs(1. - vec[i + 1] - vec[j]) < 1E-10 && std::abs(1. - vec[i] - vec[j + 1]) < 1E-10)
            break;

         t0t1_triangles.push_back({{{vec[i], vec[j + 1]}, {vec[i + 1], vec[j]}, {vec[i + 1], vec[j + 1]}}});

         //@ ----
         //@ |\ |
         //@ | \|
         //@ ----
      }

   return t0t1_triangles;
};

std::vector<T3Tddd> SymmetricSubdivisionOfTriangle(const T3Tddd& ABC, const int divide) {
   auto parametersOnTriangle = SymmetricSubdivisionOfTriangle_00_10_01(divide);
   std::vector<T3Tddd> triangle_vertices(parametersOnTriangle.size());
   auto convert = [&](const Tdd& t0t1) -> Tddd { return std::get<0>(t0t1) * std::get<0>(ABC) +
                                                        std::get<1>(t0t1) * std::get<1>(ABC) +
                                                        (1. - std::get<0>(t0t1) - std::get<1>(t0t1)) * std::get<2>(ABC); };
   int i = 0;
   for (const auto& t3tdd : parametersOnTriangle)
      triangle_vertices[i++] = {convert(std::get<0>(t3tdd)), convert(std::get<1>(t3tdd)), convert(std::get<2>(t3tdd))};
   return triangle_vertices;
};

std::vector<T3Tdd> SymmetricSubdivisionOfTriangle(const T3Tdd& ABC, const int divide) {
   auto parametersOnTriangle = SymmetricSubdivisionOfTriangle_00_10_01(divide);
   std::vector<T3Tdd> triangle_vertices(parametersOnTriangle.size());
   auto convert = [&](const Tdd& t0t1) -> Tdd { return std::get<0>(t0t1) * std::get<0>(ABC) +
                                                       std::get<1>(t0t1) * std::get<1>(ABC) +
                                                       (1. - std::get<0>(t0t1) - std::get<1>(t0t1)) * std::get<2>(ABC); };
   int i = 0;
   for (const auto& t3tdd : parametersOnTriangle)
      triangle_vertices[i++] = {convert(std::get<0>(t3tdd)), convert(std::get<1>(t3tdd)), convert(std::get<2>(t3tdd))};
   return triangle_vertices;
};

std::vector<std::array<std::array<double, 2>, 3>> SubdivideSquareIntoTriangles(const int x_divide, const int y_divide) {
   std::vector<double> vec_x(x_divide + 1);
   std::vector<double> vec_y(y_divide + 1);
   const double xmax = 1., xmin = 0.;
   const double ymax = 1., ymin = 0.;
   const double dx = (xmax - xmin) / (double)x_divide;
   const double dy = (ymax - ymin) / (double)y_divide;
   //
   for (std::size_t i = 0; i < x_divide + 1; i++)
      vec_x[i] = i * dx + xmin;
   for (std::size_t i = 0; i < y_divide + 1; i++)
      vec_y[i] = i * dy + ymin;
   //
   std::vector<std::array<std::array<double, 2>, 3>> t0t1_triangles;
   t0t1_triangles.reserve(2 * x_divide * y_divide);
   //
   for (auto i = 0; i < vec_x.size() - 1; ++i)
      for (auto j = 0; j < vec_y.size() - 1; ++j) {
         t0t1_triangles.push_back({{{vec_x[i], vec_y[j]}, {vec_x[i + 1], vec_y[j]}, {vec_x[i], vec_y[j + 1]}}});
         t0t1_triangles.push_back({{{vec_x[i], vec_y[j + 1]}, {vec_x[i + 1], vec_y[j]}, {vec_x[i + 1], vec_y[j + 1]}}});
      }
   return t0t1_triangles;
};

std::vector<std::array<std::array<double, 2>, 3>> SubdivideSquareIntoTriangles(const int divide) {
   return SubdivideSquareIntoTriangles(divide, divide);
};

/* -------------------------------------------------------------------------- */
template <size_t N, typename T>
constexpr std::array<T, N> Reverse(const std::array<T, N>& vecs) {
   std::array<T, N> reversed_vecs{};
   std::reverse_copy(vecs.begin(), vecs.end(), reversed_vecs.begin());
   return reversed_vecs;
}
template <size_t N, typename T>
constexpr bool DuplicateFreeQ(const std::array<T, N>& vecs) {
   for (size_t i = 0; const auto& v : vecs) {
      for (size_t j = i + 1; j < N; ++j) {
         if (v == vecs[j])
            return false;
      }
      i++;
   }
   return true;
}

template <size_t N1, std::size_t N2, typename T>
std::unordered_set<T> Intersection(const std::array<T, N1>& arr1, const std::array<T, N2>& arr2) {
   std::unordered_set<T> result;
   for (const auto& elem1 : arr1)
      if (std::find(arr2.begin(), arr2.end(), elem1) != arr2.end())
         result.emplace(elem1);
   return result;
}

template <size_t N, typename T>
std::vector<T> ToVector(const std::array<T, N>& arr) {
   return std::vector<T>(arr.begin(), arr.end());
}

template <size_t N0, std::size_t N1, typename T>
constexpr std::array<T, N0 + N1> Join(const std::array<T, N0>& a, const std::array<T, N1>& b) {
   std::array<T, N0 + N1> result{};
   std::ranges::copy(a, result.begin());
   std::ranges::copy(b, result.begin() + N0);
   return result;
}
/* -------------------------------------------------------------------------- */

template <typename T, std::size_t N>
constexpr std::array<T, N> RotateLeft(const std::array<T, N>& arr, const int n = 1) {
   std::array<T, N> result = arr;
   std::ranges::rotate(result, result.begin() + (n % N));
   return result;
}

template <typename T, std::size_t N>
constexpr std::array<T, N> RotateRight(const std::array<T, N>& arr, const int n = 1) {
   std::array<T, N> result = arr;
   std::ranges::rotate(result, result.begin() + (N - n % N));
   return result;
}

/* =========================================================================== */
/*                                Interpolation                               */
/* =========================================================================== */

/*DOC_EXTRACT interpolation:TriShape

### 三角形形状関数

線形の三角形形状関数は，$`\xi_2 = 1-\xi_0-\xi_1`$として，

```math
{\bf N}_{\rm l}({\boldsymbol \xi})=
\left(
\begin{array}{c}
\xi_0\\
\xi_1\\
\xi_2
\end{array}
\right)
```

2次の三角形形状関数は，$`\xi_2 = 1-\xi_0-\xi_1`$として，

```math
{\bf N}_{\rm q}({\boldsymbol \xi})=
\left(
\begin{array}{c}
\xi_0 (2  \xi_0 - 1)\\
\xi_1 (2 \xi_1 - 1)\\
\xi_2 (2 \xi_2 - 1)\\
4 \xi_0 \xi_1\\
4 \xi_1 \xi_2\\
4 \xi_2 \xi_0
\end{array}
\right)
```

<img src="triangle_vertcies_order_for_shape_functions.png" width="400">

この形状関数の使い方は，節点上の変数が$`(v_0,v_1,v_2,v_4,v_5,v_6)^{\intercal}`$の場合，
次のようにして，三角形上の値を補間できる．

```math
v({\boldsymbol \xi}) =
\left(
\begin{array}{c}
\xi_0 (2  \xi_0 - 1)\\
\xi_1 (2 \xi_1 - 1)\\
\xi_2 (2 \xi_2 - 1)\\
4 \xi_0 \xi_1\\
4 \xi_1 \xi_2\\
4 \xi_2 \xi_0
\end{array}
\right)^{\intercal}
\left(
\begin{array}{c}
v_0\\
v_1\\
v_2\\
v_3\\
v_4\\
v_5
\end{array}
\right)
=
\left(
\begin{array}{c}
v_0\\
v_1\\
v_2\\
v_3\\
v_4\\
v_5
\end{array}
\right)^{\intercal}
\left(
\begin{array}{c}
\xi_0 (2  \xi_0 - 1)\\
\xi_1 (2 \xi_1 - 1)\\
\xi_2 (2 \xi_2 - 1)\\
4 \xi_0 \xi_1\\
4 \xi_1 \xi_2\\
4 \xi_2 \xi_0
\end{array}
\right)
```

$v_0$がベクトルであっても，この計算は可能である．

Mathematicaを使った上の計算を考えてみる．Mathematicaのルールとして，

- Mathematicaでは，ベクトルの転置はない
- Row majorで行列を表現する
- `Dot[A,B]`の計算は，`B`がベクトルの場合，（Mathematicaにはベクトルの転置はないが）`B`が転置されたかのような計算結果を返す

ということを考えると，上の計算は次のようになる．
${\bf N}^{\intercal}V$の転置はいらなず，`Dot[N, V]`で計算できる．
また，$V^{\intercal}{\bf N}$の計算は，`Dot[V, N]`でそのまま計算できる．

```Mathematica
In[1]:= V = {{a, b}, {c, d}, {e, f}};
Dot[{N0, N1, N2}, V]　(*第１引数に，転置が要らないことに注意．既に転置されていると考える．*)
Dot[Transpose[V], {N0, N1, N2}]　(*この順番なら一般的な方法でOK*)

Out[2]= {a N0 + c N1 + e N2, b N0 + d N1 + f N2}

Out[3]= {a N0 + c N1 + e N2, b N0 + d N1 + f N2}
```

ちなみに，節点3と節点5の線上のパラメタは，$`t_0 = 1/2`$である．
これを2次補間の形状関数に代入すると，

```math
(N_0, N_1, N_2, N_3, N_4, N_5) = (0, t_1(2 t_1-1), t_1 (2 t_1-1), 2 t_1, 2 (1 - 2 t_1) t_1, 1 - 2 t_1)
```

となり，この線上では，節点0の影響を受けず，補間値はそれ以外の（内部）の情報からのみ決まる．

*/

template <std::size_t N>
constexpr std::array<double, N> TriShape(double t0, double t1, const auto& p0p1p2) noexcept {
   static_assert(N == 3 || N == 6, "Unsupported shape function size. Only 3 or 6 are supported.");
   auto t2 = 1 - t0 - t1;
   if constexpr (N == 3) {
      return {t0, t1, t2};
   } else if constexpr (N == 6) {
      std::array<double, 6> deflt{t0 * (2. * t0 - 1.),
                                  t1 * (2. * t1 - 1.),
                                  t2 * (2. * t2 - 1.),
                                  4. * t0 * t1,
                                  4. * t1 * t2,
                                  4. * t2 * t0};
      constexpr Tddd TriShape3 = {1., 1., -1.};
      if (!std::get<0>(p0p1p2)) {
         auto [M0_0, M0_1, M0_2] = deflt[0] * TriShape3;
         deflt[0] = 0;
         deflt[3] += M0_1;
         deflt[4] += M0_2;
         deflt[5] += M0_0;
      }
      if (!std::get<1>(p0p1p2)) {
         auto [M1_0, M1_1, M1_2] = deflt[1] * TriShape3;
         deflt[1] = 0;
         deflt[3] += M1_0;
         deflt[4] += M1_1;
         deflt[5] += M1_2;
      }
      if (!std::get<2>(p0p1p2)) {
         auto [M2_0, M2_1, M2_2] = deflt[2] * TriShape3;
         deflt[2] = 0;
         deflt[4] += M2_0;
         deflt[5] += M2_1;
         deflt[3] += M2_2;
      }
      return deflt;
   }
}

template <std::size_t N>
constexpr std::array<double, N> TriShape(std::array<double, 2> t0t1, const auto& p0p1p2) noexcept { return TriShape<N>(std::get<0>(t0t1), std::get<1>(t0t1), p0p1p2); }

template <std::size_t N>
constexpr std::array<double, N> TriShape(double t0, double t1) noexcept {
   constexpr std::array<bool, 3> ooo{true, true, true};
   return TriShape<N>(t0, t1, ooo);
}

template <std::size_t N>
constexpr std::array<double, N> TriShape(std::array<double, 2> t0t1) noexcept { return TriShape<N>(std::get<0>(t0t1), std::get<1>(t0t1)); }

template <std::size_t N, int i, int j>
constexpr std::array<double, N> D_TriShape(double t0, double t1, const auto& p0p1p2) noexcept {
   static_assert(N == 3 || N == 6, "Unsupported shape function size. Only 3 or 6 are supported.");
   auto t2 = 1 - t0 - t1;
   if constexpr (N == 3) {
      if (i == 0 && j == 0)
         return {t0, t1, t2};
      else if (i == 1 && j == 0)
         return {1., 0., -1.};
      else if (i == 0 && j == 1)
         return {0., 1., -1.};
      else
         return {0., 0., 0.};
   } else if constexpr (N == 6) {
      std::array<double, 6> deflt;
      if (i == 0 && j == 0)
         return TriShape<6>(t0, t1);
      else if (i == 1 && j == 0)
         deflt = {-1. + 4. * t0, 0., 1. - 4. * (1. - t0 - t1), 4. * t1, -4. * t1, -4. * t0 + 4. * (1 - t0 - t1)};
      else if (i == 2 && j == 0)
         deflt = {4., 0, 4., 0, 0, -8.};
      else if (i == 0 && j == 1)
         deflt = {0, -1 + 4. * t1, 1. - 4. * (1. - t0 - t1), 4 * t0, 4 * (1 - t0 - t1) - 4 * t1, -4 * t0};
      else if (i == 0 && j == 2)
         deflt = {0., 4., 4., 0., -8., 0.};
      else if (i == 1 && j == 1)
         deflt = {0, 0, 4, 4, -4, -4};
      else
         deflt = {0., 0., 0., 0., 0., 0.};

      constexpr Tddd TriShape3 = {1., 1., -1.};
      if (!std::get<0>(p0p1p2)) {
         auto [M0_0, M0_1, M0_2] = deflt[0] * TriShape3;
         deflt[0] = 0;
         deflt[3] += M0_1;
         deflt[4] += M0_2;
         deflt[5] += M0_0;
      }
      if (!std::get<1>(p0p1p2)) {
         auto [M1_0, M1_1, M1_2] = deflt[1] * TriShape3;
         deflt[1] = 0;
         deflt[3] += M1_0;
         deflt[4] += M1_1;
         deflt[5] += M1_2;
      }
      if (!std::get<2>(p0p1p2)) {
         auto [M2_0, M2_1, M2_2] = deflt[2] * TriShape3;
         deflt[2] = 0;
         deflt[4] += M2_0;
         deflt[5] += M2_1;
         deflt[3] += M2_2;
      }
      return deflt;
   }
}

template <std::size_t N, int i, int j>
constexpr std::array<double, N> D_TriShape(double t0, double t1) noexcept {
   constexpr std::array<bool, 3> ooo{true, true, true};
   return D_TriShape<N, i, j>(t0, t1, ooo);
}

/*DOC_EXTRACT interpolation:ModTriShape

### 範囲 {t_0,t_1} = {[0,1],[0,1]} -> [t0,t1]=[0,1],[0,1-t0]

前で示した通り，$`(\xi_0,\xi_1)=(\xi_0,\eta(1-\xi_0))`$と置き換えることで，矩形領域を三角形領域に変換できる．

これを利用して，形状関数の引数を$`N(\xi_0,\eta(1-\xi_0))`$として与え，
$`(\xi_0,\eta)`$をそれぞれ$`[0,1]`$の範囲で変化させることで，
三角形の形状を補間できる．

このプログムでは，
`N(x,y)=TriShape<3>(x,y)`と定義し，
`N(x,y(1-x))=ModTriShape<3>(x,y)`と定義している．

*/

template <std::size_t N>
constexpr std::array<double, N> ModTriShape(double t0, double t1, const auto& p0p1p2) noexcept {
   static_assert(N == 3 || N == 6, "Unsupported shape function size. Only 3 or 6 are supported.");
   const double t2 = 1. - t0 - t1;
   const double t0m1 = t0 - 1.;
   const double t1m1 = t1 - 1.;
   if constexpr (N == 3) {
      return {t0,
              t1 * (1. - t0),
              t0m1 * t1m1};
   } else if constexpr (N == 6) {
      std::array<double, N> deflt{{t0 * (2. * t0 - 1.), t0m1 * t1 * (1 + 2 * t0m1 * t1), t1m1 * t0m1 * (1 + 2 * t1m1 * t0 - 2 * t1), -4 * t0m1 * t0 * t1, -4 * std::pow(t0m1, 2) * t1m1 * t1, 4 * t1m1 * t0m1 * t0}};
      constexpr std::array<double, 3> TriShape3 = std::array<double, 3>{1., 1., -1.};
      if (!std::get<0>(p0p1p2)) {
         auto [M0_0, M0_1, M0_2] = deflt[0] * TriShape3;
         deflt[0] = 0;
         deflt[3] += M0_1;
         deflt[4] += M0_2;
         deflt[5] += M0_0;
      }
      if (!std::get<1>(p0p1p2)) {
         auto [M1_0, M1_1, M1_2] = deflt[1] * TriShape3;
         deflt[1] = 0;
         deflt[3] += M1_0;
         deflt[4] += M1_1;
         deflt[5] += M1_2;
      }
      if (!std::get<2>(p0p1p2)) {
         auto [M2_0, M2_1, M2_2] = deflt[2] * TriShape3;
         deflt[2] = 0;
         deflt[4] += M2_0;
         deflt[5] += M2_1;
         deflt[3] += M2_2;
      }
      return deflt;
   }
}

template <std::size_t N>
constexpr std::array<double, N> ModTriShape(double t0, double t1) noexcept { return ModTriShape<N>(t0, t1, std::array<bool, 3>{true, true, true}); }

template <std::size_t N>
constexpr std::array<double, N> ModTriShape(std::array<double, 2> t0t1) noexcept {
   return ModTriShape<N>(std::get<0>(t0t1), std::get<1>(t0t1), std::array<bool, 3>{true, true, true});
}

template <int i, int j>
constexpr std::array<double, 6> D_TriShape_Quadratic(double t0, double t1, const auto& p0p1p2) noexcept {
   std::array<double, 6>
       deflt{};
   if constexpr (i == 0 && j == 0)
      deflt = {t0 * (-1 + 2 * t0), t1 * (-1 + 2 * t1), (-1 + t0 + t1) * (-1 + 2 * t0 + 2 * t1), 4 * t0 * t1, -4 * t1 * (-1 + t0 + t1), -4 * t0 * (-1 + t0 + t1)};
   else if constexpr (i == 0 && j == 1)
      deflt = {0, -1 + 4 * t1, -3 + 4 * t0 + 4 * t1, 4 * t0, -4 * (-1 + t0 + 2 * t1), -4 * t0};
   else if constexpr (i == 0 && j == 2)
      deflt = {0, 4, 4, 0, -8, 0};
   else if constexpr (i == 1 && j == 0)
      deflt = {-1 + 4 * t0, 0, -3 + 4 * t0 + 4 * t1, 4 * t1, -4 * t1, -4 * (-1 + 2 * t0 + t1)};
   else if constexpr (i == 1 && j == 1)
      deflt = {0, 0, 4, 4, -4, -4};
   else if constexpr (i == 1 && j == 2)
      deflt = {0, 0, 0, 0, 0, 0};
   else if constexpr (i == 2 && j == 0)
      deflt = {4, 0, 4, 0, 0, -8};
   else if constexpr (i == 2 && j == 1)
      deflt = {0, 0, 0, 0, 0, 0};
   else if constexpr (i == 2 && j == 2)
      deflt = {0, 0, 0, 0, 0, 0};
   else
      deflt = {0., 0., 0., 0., 0., 0.};

   constexpr Tddd TriShape3 = {1., 1., -1.};
   if (!std::get<0>(p0p1p2)) {
      auto [M0_0, M0_1, M0_2] = deflt[0] * TriShape3;
      deflt[0] = 0;
      deflt[3] += M0_1;
      deflt[4] += M0_2;
      deflt[5] += M0_0;
   }
   if (!std::get<1>(p0p1p2)) {
      auto [M1_0, M1_1, M1_2] = deflt[1] * TriShape3;
      deflt[1] = 0;
      deflt[3] += M1_0;
      deflt[4] += M1_1;
      deflt[5] += M1_2;
   }
   if (!std::get<2>(p0p1p2)) {
      auto [M2_0, M2_1, M2_2] = deflt[2] * TriShape3;
      deflt[2] = 0;
      deflt[4] += M2_0;
      deflt[5] += M2_1;
      deflt[3] += M2_2;
   }
   return deflt;
}

template <int i, int j>
constexpr std::array<double, 6> D_TriShape_Quadratic(double t0, double t1) noexcept {
   return D_TriShape_Quadratic<i, j>(t0, t1, std::array<bool, 3>{true, true, true});
}

template <int i, int j>
constexpr std::array<double, 6> D_TriShape_Quadratic(const std::array<double, 2>& t0t1) noexcept {
   return D_TriShape_Quadratic<i, j>(std::get<0>(t0t1), std::get<1>(t0t1), std::array<bool, 3>{true, true, true});
}

/* -------------------------------------------------------------------------- */
constexpr std::array<double, 6> QuadraticTriShapeCenter(double t0, double t1, const auto& p0p1p2) noexcept {
   std::array<double, 6> deflt{((-1 + t0) * t0) / 2.,
                               ((-1 + t0) * t1 * (1 + (-1 + t0) * t1)) / 2.,
                               ((-1 + t0) * (t0 * (-1 + t1) - t1) * (-1 + t1)) / 2.,
                               (1 - t0) * (1 + (-1 + t0) * t1),
                               -((t0 * (-1 + t1) - t1) * (1 + (-1 + t0) * t1)),
                               (-1 + t0) * (t0 * (-1 + t1) - t1)};
   constexpr Tddd TriShape3 = {1., 1., -1.};
   if (!std::get<0>(p0p1p2)) {
      auto [M0_0, M0_1, M0_2] = deflt[0] * TriShape3;
      deflt[0] = 0;
      deflt[3] += M0_1;
      deflt[4] += M0_2;
      deflt[5] += M0_0;
   }
   if (!std::get<1>(p0p1p2)) {
      auto [M1_0, M1_1, M1_2] = deflt[1] * TriShape3;
      deflt[1] = 0;
      deflt[3] += M1_0;
      deflt[4] += M1_1;
      deflt[5] += M1_2;
   }
   if (!std::get<2>(p0p1p2)) {
      auto [M2_0, M2_1, M2_2] = deflt[2] * TriShape3;
      deflt[2] = 0;
      deflt[4] += M2_0;
      deflt[5] += M2_1;
      deflt[3] += M2_2;
   }
   return deflt;
}

template <int i, int j>
constexpr std::array<double, 6> D_QuadraticTriShapeCenter(double t0, double t1, const auto& p0p1p2) noexcept {
   std::array<double, 6> deflt{};
   if constexpr (i == 0 && j == 0)
      deflt = {((-1 + t0) * t0) / 2., ((-1 + t0) * t1 * (1 + (-1 + t0) * t1)) / 2., ((-1 + t0) * (t0 * (-1 + t1) - t1) * (-1 + t1)) / 2., (1 - t0) * (1 + (-1 + t0) * t1), -((t0 * (-1 + t1) - t1) * (1 + (-1 + t0) * t1)), (-1 + t0) * (t0 * (-1 + t1) - t1)};
   else if constexpr (i == 0 && j == 1)
      deflt = {0, ((-1 + t0) * (1 + 2 * (-1 + t0) * t1)) / 2., ((-1 + t0) * (1 + 2 * t0 * (-1 + t1) - 2 * t1)) / 2., -std::pow(-1 + t0, 2), -(std::pow(-1 + t0, 2) * (-1 + 2 * t1)), std::pow(-1 + t0, 2)};
   else if constexpr (i == 0 && j == 2)
      deflt = {0, std::pow(-1 + t0, 2), std::pow(-1 + t0, 2), 0, -2 * std::pow(-1 + t0, 2), 0};
   else if constexpr (i == 1 && j == 0)
      deflt = {-0.5 + t0, (t1 * (1 + 2 * (-1 + t0) * t1)) / 2., ((1 + 2 * t0 * (-1 + t1) - 2 * t1) * (-1 + t1)) / 2., -1 - 2 * (-1 + t0) * t1, 1 + 2 * t1 * (-1 + t0 + t1 - t0 * t1), 1 + 2 * t0 * (-1 + t1) - 2 * t1};
   else if constexpr (i == 1 && j == 1)
      deflt = {0, 0.5 + 2 * (-1 + t0) * t1, 1.5 + 2 * t0 * (-1 + t1) - 2 * t1, 2 - 2 * t0, -2 * (-1 + t0) * (-1 + 2 * t1), 2 * (-1 + t0)};
   else if constexpr (i == 1 && j == 2)
      deflt = {0, 2 * (-1 + t0), 2 * (-1 + t0), 0, 4 - 4 * t0, 0};
   else if constexpr (i == 2 && j == 0)
      deflt = {1, std::pow(t1, 2), std::pow(-1 + t1, 2), -2 * t1, -2 * (-1 + t1) * t1, 2 * (-1 + t1)};
   else if constexpr (i == 2 && j == 1)
      deflt = {0, 2 * t1, 2 * (-1 + t1), -2, 2 - 4 * t1, 2};
   else if constexpr (i == 2 && j == 2)
      deflt = {0, 2, 2, 0, -4, 0};
   else
      deflt = {0., 0., 0., 0., 0., 0.};

   constexpr Tddd TriShape3 = {1., 1., -1.};
   if (!std::get<0>(p0p1p2)) {
      auto [M0_0, M0_1, M0_2] = deflt[0] * TriShape3;
      deflt[0] = 0;
      deflt[3] += M0_1;
      deflt[4] += M0_2;
      deflt[5] += M0_0;
   }
   if (!std::get<1>(p0p1p2)) {
      auto [M1_0, M1_1, M1_2] = deflt[1] * TriShape3;
      deflt[1] = 0;
      deflt[3] += M1_0;
      deflt[4] += M1_1;
      deflt[5] += M1_2;
   }
   if (!std::get<2>(p0p1p2)) {
      auto [M2_0, M2_1, M2_2] = deflt[2] * TriShape3;
      deflt[2] = 0;
      deflt[4] += M2_0;
      deflt[5] += M2_1;
      deflt[3] += M2_2;
   }
   return deflt;
}

template <int i, int j>
constexpr std::array<double, 6> D_QuadraticTriShapeCenter(double t0, double t1) noexcept {
   constexpr std::array<bool, 3> ooo{true, true, true};
   return D_QuadraticTriShapeCenter<i, j>(t0, t1, ooo);
}

constexpr std::array<double, 6> QuadraticTriShapeCenter(double t0, double t1) noexcept {
   return QuadraticTriShapeCenter(t0, t1, std::array<bool, 3>{true, true, true});
}

/* -------------------------------------------------------------------------- */

namespace std {

template <typename T, std::size_t N>
struct hash<std::array<T, N>> {
   std::size_t operator()(std::array<T, N> const& arr) const {
      std::size_t seed = 0;
      for (const auto& elem : arr) {
         seed ^= hash<T>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
   }
};
}  // namespace std
/* -------------------------------------------------------------------------- */

// T2Tdd Inverse(const T2Tdd& M) {
//    const auto [x00, x01] = std::get<0>(M);
//    const auto [x10, x11] = std::get<1>(M);
//    const double det = std::fma(-x01, x10, x00 * x11);
//    return {{{x11 / det, -x01 / det}, {-x10 / det, x00 / det}}};
// };

// T3Tddd Inverse(const T3Tddd& mat) {
//    const auto [x00, x01, x02] = std::get<0>(mat);
//    const auto [x10, x11, x12] = std::get<1>(mat);
//    const auto [x20, x21, x22] = std::get<2>(mat);
//    const double inv_det = 1. / (std::fma(std::fma(-x02, x11, x01 * x12), x20, std::fma(std::fma(x02, x10, -x00 * x12), x21, std::fma(std::fma(-x01, x10, x00 * x11), x22, 0))));
//    return {Tddd{inv_det * std::fma(-x12, x21, x11 * x22),
//                 inv_det * std::fma(x02, x21, -x01 * x22),
//                 inv_det * std::fma(-x02, x11, x01 * x12)},
//            Tddd{inv_det * std::fma(x12, x20, -x10 * x22),
//                 inv_det * std::fma(-x02, x20, x00 * x22),
//                 inv_det * std::fma(x02, x10, -x00 * x12)},
//            Tddd{inv_det * std::fma(-x11, x20, x10 * x21),
//                 inv_det * std::fma(x01, x20, -x00 * x21),
//                 inv_det * std::fma(-x01, x10, x00 * x11)}};
//    // auto bc = Cross(std::get<1>(mat), std::get<2>(mat));
//    // return T3Tddd{bc, Cross(std::get<2>(mat), std::get<0>(mat)), Cross(std::get<0>(mat), std::get<1>(mat))} / (Dot(std::get<0>(mat), bc));
// };

template <typename T>
void IdentityMatrix(T& M) {
   std::size_t i = 0;
   for (auto& m : M) {
      m.fill(0.);
      m[i] = 1.;
      i++;
   }
};

#endif