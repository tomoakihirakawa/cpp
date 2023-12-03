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
// 2d array + 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator+=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) += (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}
// 2d array + 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator+=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) += (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
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
template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator+(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept {
   return arr += ARR;
}
// 1d array + 2d array
template <size_t N, size_t M, typename T, typename TT>
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
// 2d array -= 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator-=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) -= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}
// 2d array -= 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator-=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) -= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
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
template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator-(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr -= ARR; }

// 1d array - 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator-(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept {
   for (auto i = 0; auto& a : arr)
      a = ARR[i++] - a;
   return arr;
}

template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<T>::value && !is_std_array<TT>::value, std::array<std::array<TT, M>, N>> operator-(const std::array<std::array<TT, M>, N>& ARR, const std::array<std::array<T, M>, N>& arr) noexcept {
   std::array<std::array<TT, M>, N> result{};
   for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
         result[i][j] = ARR[i][j] - arr[i][j];
      }
   }
   return result;
}

template <size_t N0, size_t N1, size_t N2, typename T, typename TT>
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
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator*=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) *= (std::get<Is>(ARR))), ...); }(std::make_index_sequence<N>());
   return arr;
}
// 2d array *= 2d array
template <size_t N, size_t M, typename T, typename TT>
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
template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator*(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr *= ARR; }
// 1d array * 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N> operator*(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept {
   for (auto i = 0; auto& a : arr) a *= ARR[i++];
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
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator/=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
// 2d array /= 2d array
template <size_t N, size_t M, typename T, typename TT>
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
template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator/(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++]; /*std::array<T, M> /= TT*/
   return arr;
}
// 1d array / 2d array
template <size_t N, size_t M, typename T, typename TT>
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
template <size_t N0, size_t N1, typename T>
constexpr std::array<std::array<T, N0>, N1> Transpose(const std::array<std::array<T, N1>, N0>& M) noexcept {
   std::array<std::array<T, N0>, N1> ret;
   for (auto i = 0; const auto& Mi : M) {
      for (auto j = 0; const auto& Mij : Mi)
         ret[j++][i] = Mij;
      i++;
   }
   return ret;
}

template <size_t N, typename T>
constexpr std::array<std::vector<T>, N> Transpose(const std::vector<std::array<T, N>>& M) noexcept {
   std::array<std::vector<T>, N> ret;
   for (size_t j = 0; j < N; ++j)
      ret[j].reserve(M.size());
   for (const auto& Mi : M) {
      size_t j = 0;
      for (const auto& Mij : Mi) {
         ret[j++].push_back(Mij);
      }
   }
   return ret;
}

template <typename T, size_t N>
constexpr T Min(const std::array<T, N>& arr) noexcept { return *std::min_element(arr.begin(), arr.end()); }

template <typename T, size_t N>
constexpr T Max(const std::array<T, N>& arr) noexcept { return *std::max_element(arr.begin(), arr.end()); }

template <typename T, size_t N>
constexpr std::array<T, 2> MinMax(const std::array<T, N>& arr) noexcept {
   auto minmax = std::minmax_element(arr.begin(), arr.end());
   return {*minmax.first, *minmax.second};
}
constexpr Tdd MinMax(const Tdd& A) noexcept { return ((A[0] < A[1]) ? Tdd{A[0], A[1]} : Tdd{A[1], A[0]}); };

template <typename T, size_t N0, size_t N1>
constexpr auto MinMaxColumns(const std::array<std::array<T, N0>, N1>& arr) {
   std::array<std::array<T, 2>, N0> minmax_array;
   for (size_t i = 0; i < N0; ++i)
      minmax_array[i] = {arr[0][i], arr[0][i]};

   for (const auto& row : arr) {
      for (size_t i = 0; const auto& a : row) {
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

template <typename T, size_t N>
constexpr auto MinMaxColumns(const std::vector<std::array<T, N>>& vec_arr) noexcept {
   std::array<std::array<T, 2>, N> minmax_array;
   for (size_t i = 0; i < N; ++i)
      minmax_array[i] = {vec_arr[0][i], vec_arr[0][i]};
   for (const auto& arr : vec_arr) {
      for (size_t i = 0; const auto& a : arr) {
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

template <size_t N0, size_t N1, typename T>
constexpr std::array<T, N1> Total(const std::array<std::array<T, N1>, N0>& arr) noexcept {
   std::array<T, N1> ret{};
   for (size_t i = 0; i < N0; ++i) {
      for (size_t j = 0; j < N1; ++j) {
         ret[j] += arr[i][j];
      }
   }
   return ret;
}

/* -------------------------------------------------------------------------- */

template <typename T, size_t N1, size_t N2>
constexpr std::array<std::array<T, N2>, N1> TensorProduct(const T scalar, const std::array<std::array<T, N2>, N1>& matrix) noexcept {
   std::array<std::array<T, N2>, N1> result{};
   for (size_t i = 0; i < N1; ++i) {
      for (size_t j = 0; j < N2; ++j) {
         result[i][j] = scalar * matrix[i][j];
      }
   }
   return result;
}

template <typename T, size_t N1, size_t N2>
constexpr std::array<std::array<T, N2>, N1> TensorProduct(const std::array<T, N1>& vec1, const std::array<T, N2>& vec2) noexcept {
   std::array<std::array<T, N2>, N1> ret{};
   for (size_t m = 0; m < N1; ++m) {
      for (size_t j = 0; j < N2; ++j) {
         ret[m][j] = vec1[m] * vec2[j];
      }
   }
   return ret;
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

template <typename T, size_t N>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, T>::type
Dot(const std::array<T, N>& arr, const std::array<T, N>& ARR) noexcept {
   T ret = 0;
   // int i = 0;
   // std::ranges::for_each(arr, [&](const auto& a) { ret = std::fma(a, ARR[i++], ret); });
   // rewrite this using for
   for (size_t i = 0; i < N; ++i)
      ret += arr[i] * ARR[i];
   return ret;
}

template <typename T, size_t N0, size_t N1, size_t N2>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<T, N2>, N1>>::type
Dot(const std::array<std::array<T, N0>, N1>& arr, const std::array<std::array<T, N2>, N0>& ARR) noexcept {
   std::array<std::array<T, N2>, N1> ret{};
   T sum = 0;
   for (size_t i = 0; i < N1; ++i) {
      for (size_t j = 0; j < N2; ++j) {
         sum = 0;
         for (size_t k = 0; k < N0; ++k) {
            sum += arr[i][k] * ARR[k][j];
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

template <size_t N0, size_t N1, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, void>::type
Fill(std::array<std::array<T, N1>, N0>& arr, const T d) noexcept {
   for (auto& a : arr)
      Fill(a, d);
}

template <size_t N0, size_t N1, size_t N2, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, void>::type
Fill(std::array<std::array<std::array<T, N2>, N1>, N0>& arr, const T d) noexcept {
   for (auto& a : arr)
      Fill(a, d);
}

template <size_t N0, size_t N1, size_t N2, size_t N3, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, void>::type
Fill(std::array<std::array<std::array<std::array<T, N3>, N2>, N1>, N0>& arr, const T d) noexcept {
   for (auto& a : arr)
      Fill(a, d);
}
/* -------------------------------------------------------------------------- */
template <size_t N0, size_t N1, typename T>
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
   for (size_t i = 0; const auto& ARR : ARRARR) {
      for (size_t j = 0; j < N1; ++j)
         ret[j] += arr[i] * ARR[j];
      i++;
   }

   return ret;
}

template <size_t N0, size_t N1, size_t N2, typename T>
constexpr std::array<std::array<T, N2>, N1> Dot(const std::array<T, N0>& arr, const std::array<std::array<std::array<T, N2>, N1>, N0>& ARRARR) noexcept {
   std::array<std::array<T, N2>, N1> ret = {};
   for (size_t i = 0; const auto& ARR : ARRARR) {
      for (size_t j = 0; j < N1; ++j)
         ret[j] += arr[i] * ARR[j];
      i++;
   }
   return ret;
}

template <size_t N0, size_t N1, size_t N2, size_t N3, typename T>
constexpr std::array<std::array<std::array<T, N3>, N2>, N1> Dot(const std::array<T, N0>& arr, const std::array<std::array<std::array<std::array<T, N3>, N2>, N1>, N0>& ARRARR) noexcept {
   std::array<std::array<std::array<T, N3>, N2>, N1> ret = {};
   for (int i = 0; const auto& ARR : ARRARR) {
      for (size_t j = 0; j < N1; ++j)
         ret[j] += arr[i] * ARR[j];
      i++;
   }
   return ret;
}

template <size_t N0, size_t N1, size_t N2, size_t N3, size_t N4, typename T>
constexpr std::array<std::array<std::array<std::array<T, N4>, N3>, N2>, N1> Dot(const std::array<T, N0>& arr, const std::array<std::array<std::array<std::array<std::array<T, N4>, N3>, N2>, N1>, N0>& ARRARR) noexcept {
   std::array<std::array<std::array<std::array<T, N4>, N3>, N2>, N1> ret = {};
   for (int i = 0; const auto& ARR : ARRARR) {
      for (size_t j = 0; j < N1; ++j)
         ret[j] += arr[i] * ARR[j];
      i++;
   }
   return ret;
}

/* -------------------------------------------------------------------------- */

template <size_t N0, size_t N1, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<T, N0>>::type
Dot(const std::array<std::array<T, N1>, N0>& AA, const std::array<T, N1>& arr) noexcept {
   std::array<T, N0> ret = {};
   for (int i = 0; const auto& A : AA) {
      ret[i] = Dot(A, arr);
      i++;
   }
   return ret;
}

template <size_t N0, size_t N1, size_t N2, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<T, N1>, N0>>::type
Dot(const std::array<std::array<std::array<T, N2>, N1>, N0>& AAA, const std::array<T, N2>& arr) noexcept {
   std::array<std::array<T, N1>, N0> ret = {};
   for (int i = 0; const auto& AA : AAA) {
      for (int j = 0; const auto& A : AA) {
         ret[i][j] = Dot(A, arr);
         j++;
      }
      i++;
   }
   return ret;
}

template <size_t N0, size_t N1, size_t N2, size_t N3, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<std::array<T, N3>, N1>, N0>>::type
Dot(const std::array<std::array<std::array<T, N2>, N1>, N0>& AAA, const std::array<std::array<T, N3>, N2>& arr) noexcept {
   std::array<std::array<std::array<T, N3>, N1>, N0> ret = {};
   for (int i = 0; const auto& AA : AAA) {
      for (int j = 0; const auto& A : AA) {
         ret[i][j++] = Dot(A, arr);
      }
      i++;
   }
   return ret;
}

template <size_t N0, size_t N1, size_t N2, size_t N3, size_t N4, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<std::array<std::array<T, N4>, N3>, N1>, N0>>::type
Dot(const std::array<std::array<std::array<T, N2>, N1>, N0>& AAA, const std::array<std::array<std::array<T, N4>, N3>, N2>& arr) noexcept {
   std::array<std::array<std::array<std::array<T, N4>, N3>, N1>, N0> ret = {};
   for (int i = 0; const auto& AA : AAA) {
      for (int j = 0; const auto& A : AA) {
         ret[i][j++] = Dot(A, arr);
      }
      i++;
   }
   return ret;
}

template <size_t N0, size_t N1, size_t N2, size_t N3, typename T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, std::array<std::array<std::array<T, N2>, N1>, N0>>::type
Dot(const std::array<std::array<std::array<std::array<T, N3>, N2>, N1>, N0>& AAAA, const std::array<T, N3>& arr) noexcept {
   std::array<std::array<std::array<T, N2>, N1>, N0> ret = {};
   for (int i = 0; const auto& AAA : AAAA) {
      for (int j = 0; const auto& AA : AAA) {
         for (int k = 0; const auto& A : AA) {
            ret[i][j][k] = Dot(A, arr);
            k++;
         }
         j++;
      }
      i++;
   }
   return ret;
}

/* -------------------------------------------------------------------------- */
// template <size_t N0, size_t N1, typename T>
// constexpr std::array<T, N1> Dot(const std::array<std::array<T, N0>, N1>& arr, const std::array<T, N0>& ARR) noexcept {
//    /*
//    {{A0,A1,A2},{B0,B1,B2},{C0,C1,C2},{D0,D1,D2}}.{a,b,c}={a*A0+b*A1+c*A2,a*B0+b*B1+c*B2,a*C0+b*C1+c*C2,a*D0+b*D1+c*D2}
//    */
//    std::array<T, N1> ret;
//    for_each(ret, arr, [&](auto& r, auto& ar) { r = Dot(ar, ARR); });
//    return ret;
// }

template <typename... Ts, typename... Is, size_t... Indices>
constexpr void ToArrayImpl(const std::tuple<Ts...>& tup, std::array<std::common_type_t<Ts...>, sizeof...(Ts)>& arr, std::index_sequence<Indices...>) noexcept {
   ((arr[Indices] = std::get<Indices>(tup)), ...);
}

template <typename... Ts>
constexpr auto ToArray(const std::tuple<Ts...>& tup) {
   using T = std::common_type_t<Ts...>;
   constexpr size_t N = sizeof...(Ts);
   std::array<T, N> arr{};
   ToArrayImpl(tup, arr, std::make_index_sequence<N>());
   return arr;
}

template <typename T, size_t N>
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
      sum += a * a;
   return sum;
}

template <size_t N, typename T>
constexpr T Norm(const std::array<T, N>& arr) noexcept {
   return std::sqrt(Dot(arr, arr));
}

template <size_t N, size_t M, typename T>
constexpr auto Norm(const std::array<std::array<T, M>, N>& mat) noexcept {
   //$ Mathematica implementation returns the maximum of the singular values of the matrix, which is not the same as the Frobenius norm implemented here.
   static_assert(!std::is_array_v<T>, "Norm is not defined for arrays deeper than 2 dimensions.");
   T sum = {};
   for (const auto& row : mat)
      sum += NormSquared(row);
   return std::sqrt(sum);
}

template <size_t N, typename T>
constexpr std::array<T, N> Normalize(std::array<T, N> arr) noexcept {
   static_assert(N > 0, "Array must have at least one element.");
   T norm = Norm(arr);
   auto c = (norm == static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(1) / norm;
   std::ranges::transform(arr, arr.begin(), [&c](T a) { return a * c; });
   return arr;
}

constexpr double Norm(const double value) noexcept {
   return std::abs(value);
}

template <size_t M = 0, size_t N, typename T>
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

// double Det(const T3Tddd& M) {
//    return -(std::get<0>(std::get<2>(M)) * std::get<1>(std::get<1>(M)) * std::get<2>(std::get<0>(M))) +
//           std::get<0>(std::get<1>(M)) * std::get<1>(std::get<2>(M)) * std::get<2>(std::get<0>(M)) +
//           std::get<0>(std::get<2>(M)) * std::get<1>(std::get<0>(M)) * std::get<2>(std::get<1>(M)) -
//           std::get<0>(std::get<0>(M)) * std::get<1>(std::get<2>(M)) * std::get<2>(std::get<1>(M)) -
//           std::get<0>(std::get<1>(M)) * std::get<1>(std::get<0>(M)) * std::get<2>(std::get<2>(M)) +
//           std::get<0>(std::get<0>(M)) * std::get<1>(std::get<1>(M)) * std::get<2>(std::get<2>(M));
// };

double Det(const T3Tddd& M) {
   // Unpacking each row of the matrix
   auto [a, b, c] = std::get<0>(M);  // First row
   auto [d, e, f] = std::get<1>(M);  // Second row
   auto [g, h, i] = std::get<2>(M);  // Third row

   // Calculating the determinant using std::fma
   return std::fma(a, std::fma(e, i, -f * h),
                   std::fma(b, std::fma(f, g, -d * i),
                            std::fma(c, std::fma(d, h, -e * g), 0)));
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

template <size_t N1, size_t N2, typename T>
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

template <size_t N0, size_t N1, typename T>
constexpr std::array<T, N0 + N1> Join(const std::array<T, N0>& a, const std::array<T, N1>& b) {
   std::array<T, N0 + N1> result{};
   std::ranges::copy(a, result.begin());
   std::ranges::copy(b, result.begin() + N0);
   return result;
}
/* -------------------------------------------------------------------------- */

template <typename T, size_t N>
constexpr std::array<T, N> RotateLeft(const std::array<T, N>& arr, const int n = 1) {
   std::array<T, N> result = arr;
   std::ranges::rotate(result, result.begin() + (n % N));
   return result;
}

template <typename T, size_t N>
constexpr std::array<T, N> RotateRight(const std::array<T, N>& arr, const int n = 1) {
   std::array<T, N> result = arr;
   std::ranges::rotate(result, result.begin() + (N - n % N));
   return result;
}

/* =========================================================================== */
/*                                Interpolation                               */
/* =========================================================================== */
// 沢山ヘッダーを作るとごちゃごちゃしてくるので，まずはarrayのユーティリティーについては，このヘッダーにまとめていく
template <typename T, size_t N>
constexpr std::array<T, N> TriShape(T t0, T t1) noexcept {
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
constexpr std::array<double, N> TriShape(double t0, double t1) noexcept { return TriShape<double, N>(t0, t1); }

/*DOC_EXTRACT interpolation:ModTriShape

## 範囲を修正した三角形形状関数

普通の三角形形状関数は，$`{\mathbf N}=(N_0,N_1,N_2) = (t_0,t_1,1-t_0-t_1)`$．
これを使った，$`{\rm Dot}({\mathbf N},\{{\mathbf X_0},{\mathbf X_1},{\mathbf X_2}\})`$は，$`t_0,t_1=[0,1]`$で平行四辺形を作る．
$`t_0,t_1=[0,1]`$の範囲で，三角形を形成するように変数変換したいことがある．
そのたびに，変数変換をプログラムするのは面倒なので，予め形状関数自体を変更しておく．
変更した形状関数は，`ModTriShape`にあるように，
3点の場合は，

```math
\begin{align}
N_0 &= t_0 \\
N_1 &= -t_1(t_0-1) \\
N_2 &= (t_0-1)(t_1-1)
\end{align}
```

6点の場合は，

```math
\begin{align}
N_0 &= t_0(2t_0-1) \\
N_1 &= t_1(2t_1-1) \\
N_2 &= (1-t_0-t_1)(2(1-t_0-t_1)-1) \\
N_3 &= 4t_0t_1 \\
N_4 &= 4t_1(1-t_0-t_1) \\
N_5 &= 4t_0(1-t_0-t_1)
\end{align}
```

*/

template <typename T, size_t N>
constexpr std::array<T, N> ModTriShape(const T& t0, const T& t1, const auto& p0p1p2) noexcept {
   // b! どのポインターにどれだけの係数を加えるかを得られるようにする
   // b!　これは，線形補間であっても利用できる
   //
   static_assert(N == 3 || N == 6, "Unsupported shape function size. Only 3 or 6 are supported.");
   const double t2 = 1. - t0 - t1;
   const double t0m1 = t0 - 1.;
   const double t1m1 = t1 - 1.;
   if constexpr (N == 3) {
      return {t0,
              -t1 * t0m1,
              t0m1 * t1m1};
   } else if constexpr (N == 6) {
      std::array<T, N> deflt{{t0 * (-1 + 2 * t0),
                              t0m1 * t1 * (1 + 2 * t0m1 * t1),
                              t1m1 * t0m1 * (1 + 2 * t1m1 * t0 - 2 * t1),
                              -4 * t0m1 * t0 * t1,
                              -4 * std::pow(t0m1, 2) * t1m1 * t1,
                              4 * t1m1 * t0m1 * t0}};
      const auto TriShape3 = TriShape<3>(static_cast<T>(1), static_cast<T>(1));
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

template <size_t N>
constexpr std::array<double, N> ModTriShape(double t0, double t1) noexcept { return ModTriShape<double, N>(t0, t1, std::array<bool, 3>{true, true, true}); }

template <typename T, size_t N>
constexpr std::array<T, N> ModTriShape(T t0, T t1) noexcept { return ModTriShape<T, N>(t0, t1, std::array<bool, 3>{true, true, true}); }

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

namespace std {

template <typename T, size_t N>
struct hash<std::array<T, N>> {
   size_t operator()(std::array<T, N> const& arr) const {
      size_t seed = 0;
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
   size_t i = 0;
   for (auto& m : M) {
      m.fill(0.);
      m[i++] = 1.;
   }
};

#endif