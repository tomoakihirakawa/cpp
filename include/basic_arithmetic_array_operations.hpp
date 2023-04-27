#ifndef basic_arithmetic_array_operations_H
#define basic_arithmetic_array_operations_H

#include <algorithm>
#include <array>
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

template <typename>
struct is_std_array : std::false_type {};

template <template <typename, std::size_t> class Array, typename T, std::size_t N>
struct is_std_array<Array<T, N>> : std::true_type {};
/* -------------------------------------------------------------------------- */
template <size_t N, typename T>
constexpr std::array<T, N> operator-(std::array<T, N> arr) noexcept {
   for (auto& a : arr)
      a = -a;
   return arr;
}
/* -------------------------------- += ------------------------------- */
// 1d array - value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>&> operator+=(std::array<T, N>& arr /*ref*/, const TT& d) noexcept {
   for (auto& a : arr) a += d;
   return arr;
}
// 1d array - 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N>& operator+=(std::array<T, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) += (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
   return arr;
}
// 2d array - 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator+=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) += (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
   return arr;
}
// 2d array - 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator+=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) += (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
   return arr;
}
/* -------------------------------------------------------------------------- */
// 1d array - value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator+(std::array<T, N> arr /*copy*/, const TT& d) noexcept { return arr += d; }
// value - 1d array
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator+(const TT& d, std::array<T, N> arr /*copy*/) noexcept { return arr + d; }
// 1d array - 1d array
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator+(std::array<T, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr += ARR; }
// 2d array - 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator+(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept {
   return arr += ARR;
}
// 1d array - 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N> operator+(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept { return arr += ARR; }
/* -------------------------------- -= ------------------------------- */
// 1d array - value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>&> operator-=(std::array<T, N>& arr /*ref*/, const TT& d) noexcept {
   for (auto& a : arr) a -= d;
   return arr;
}
// 1d array - 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N>& operator-=(std::array<T, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) -= (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
   return arr;
}
// 2d array - 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator-=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) -= (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
   return arr;
}
// 2d array - 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator-=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) -= (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
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
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator-(const TT& d, std::array<T, N> arr /*copy*/) noexcept {
   for (auto& a : arr) a = d - a;
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
constexpr std::array<std::array<T, M>, N> operator-(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept {
   for (auto i = 0; auto& a : arr) a = ARR[i++] - a;
   return arr;
}
/* -------------------------------- *= ------------------------------- */
// 1d array - value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>&> operator*=(std::array<T, N>& arr /*ref*/, const TT& d) noexcept {
   for (auto& a : arr) a *= d;
   return arr;
}
// 1d array - 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N>& operator*=(std::array<T, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) *= (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
   return arr;
}
// 2d array - 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator*=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) *= (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
   return arr;
}
// 2d array - 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator*=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   if constexpr (N > 0)
      [&]<size_t... Is>(std::index_sequence<Is...>) { ((std::get<Is>(arr) *= (std::get<Is>(ARR))), ...); }
   (std::make_index_sequence<N>());
   return arr;
}
/* -------------------------------------------------------------------------- */
// 1d array - value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator*(std::array<T, N> arr /*copy*/, const TT& d) noexcept {
   for (auto& a : arr) a *= d;
   return arr;
}
// value - 1d array
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator*(const TT& d, std::array<T, N> arr /*copy*/) noexcept {
   for (auto& a : arr) a *= d;
   return arr;
}
// 1d array - 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N> operator*(std::array<T, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr *= ARR; }
// 2d array - 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator*(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept { return arr *= ARR; }
// 1d array - 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N> operator*(const std::array<TT, N>& ARR, std::array<std::array<T, M>, N> arr /*copy*/) noexcept {
   for (auto i = 0; auto& a : arr) a *= ARR[i++];
   return arr;
}
/* -------------------------------- /= ------------------------------- */
// 1d array - value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>&> operator/=(std::array<T, N>& arr /*ref*/, const TT& d) noexcept {
   for (auto& a : arr) a /= d;
   return arr;
}
// 1d array - 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N>& operator/=(std::array<T, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
// 2d array - 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator/=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
// 2d array - 2d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::array<std::array<T, M>, N>& operator/=(std::array<std::array<T, M>, N>& arr /*ref*/, const std::array<std::array<TT, M>, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
/* -------------------------------------------------------------------------- */
// 1d array - value
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator/(std::array<T, N> arr /*copy*/, const TT& d) noexcept {
   for (auto& a : arr) a /= d;
   return arr;
}
// value - 1d array
template <size_t N, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<T, N>> operator/(const TT& d, std::array<T, N> arr /*copy*/) noexcept {
   for (auto& a : arr) a = d / a;
   return arr;
}
// 1d array - 1d array
template <size_t N, typename T, typename TT>
constexpr std::array<T, N> operator/(std::array<T, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++];
   return arr;
}
// 2d array - 1d array
template <size_t N, size_t M, typename T, typename TT>
constexpr std::enable_if_t<!is_std_array<TT>::value, std::array<std::array<T, M>, N>> operator/(std::array<std::array<T, M>, N> arr /*copy*/, const std::array<TT, N>& ARR) noexcept {
   for (auto i = 0; auto& a : arr) a /= ARR[i++]; /*std::array<T, M> += TT*/
   return arr;
}
// 1d array - 2d array
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
      }
      (std::make_index_sequence<N>());
   }
}

template <size_t N, typename T, typename TT, typename Func>
constexpr void for_each(std::array<T, N>& arr, const std::array<TT, N>& ARR, const Func& func) {
   if constexpr (N > 0) {
      [&]<size_t... Is>(std::index_sequence<Is...>) {
         (func(std::get<Is>(arr), std::get<Is>(ARR)), ...);
      }
      (std::make_index_sequence<N>());
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
constexpr T Total(const std::array<T, N>& arr) noexcept {
   T ret = 0;
   std::ranges::for_each(arr, [&ret](const auto& a) { ret += a; });
   return ret;
}

template <typename T, size_t N>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, T>::type
Dot(const std::array<T, N>& arr, const std::array<T, N>& ARR) noexcept {
   T ret = 0;
   for (size_t i = 0; i < N; ++i) {
      ret += arr[i] * ARR[i];
   }
   return ret;
}

template <typename T, size_t N0, size_t N1, size_t N2>
constexpr std::array<std::array<T, N2>, N1> Dot(const std::array<std::array<T, N0>, N1>& arr, const std::array<std::array<T, N2>, N0>& ARR) noexcept {
   std::array<std::array<T, N2>, N1> ret{};
   for (size_t i = 0; i < N1; ++i) {
      for (size_t j = 0; j < N2; ++j) {
         T sum = 0;
         for (size_t k = 0; k < N0; ++k) {
            sum += arr[i][k] * ARR[k][j];
         }
         ret[i][j] = sum;
      }
   }
   return ret;
}

template <size_t N0, size_t N1, typename T>
constexpr std::array<T, N1> Dot(const std::array<T, N0>& arr, const std::array<std::array<T, N1>, N0>& ARRARR) noexcept {
   /*
   {a,b,c}.{{A0,A1,A2,A3},{B0,B1,B2,B3},{C0,C1,C2,C3}}
   */
   std::array<T, N1> ret = {};
   // for_each(ARRARR, arr, [&](auto &ARR, auto &a) {
   //    for (size_t j = 0; j < N1; ++j) {
   //       ret[j] += a * ARR[j];
   //    }
   // });
   //
   for (int i = 0; const auto& ARR : ARRARR) {
      for (size_t j = 0; j < N1; ++j)
         ret[j] += arr[i] * ARR[j];
      i++;
   }

   return ret;
}

template <size_t N0, size_t N1, typename T>
constexpr std::array<T, N1> Dot(const std::array<std::array<T, N0>, N1>& arr, const std::array<T, N0>& ARR) noexcept {
   /*
   {{A0,A1,A2},{B0,B1,B2},{C0,C1,C2},{D0,D1,D2}}.{a,b,c}={a*A0+b*A1+c*A2,a*B0+b*B1+c*B2,a*C0+b*C1+c*C2,a*D0+b*D1+c*D2}
   */
   std::array<T, N1> ret;
   for_each(ret, arr, [&](auto& r, auto& ar) { r = Dot(ar, ARR); });
   return ret;
}

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

template <size_t N, typename T>
constexpr T Norm(const std::array<T, N>& arr) noexcept {
   T norm = static_cast<T>(0);
   for (const auto& a : arr)
      norm += a * a;
   return std::sqrt(norm);
}

template <size_t N, typename T>
constexpr std::array<T, N> Normalize(std::array<T, N> arr) noexcept {
   static_assert(N > 0, "Array must have at least one element.");
   T norm = Norm(arr);
   if (norm == static_cast<T>(0))
      return arr;
   else {
      for (auto& a : arr)
         a /= norm;
      return arr;
   }
}

template <size_t M = 0, size_t N, typename T>
constexpr T RootMeanSquare(const std::array<T, N>& arr) noexcept { return std::sqrt(Dot(arr, arr) / N); }

template <typename T>
constexpr std::array<T, 3> Cross(const std::array<T, 3>& A, const std::array<T, 3>& B) noexcept {
   return {{A[1] * std::get<2>(B) - std::get<2>(A) * std::get<1>(B),
            std::get<2>(A) * std::get<0>(B) - A[0] * std::get<2>(B),
            A[0] * std::get<1>(B) - A[1] * std::get<0>(B)}};
};
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
   std::copy(a.begin(), a.end(), result.begin());
   std::copy(b.begin(), b.end(), result.begin() + N0);
   return result;
};
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

template <typename T, size_t N>
constexpr std::array<T, N> ModTriShape(const T& t0, const T& t1, const auto& p0p1p2) noexcept {
   // b! どのポインターにどれだけの係数を加えるかを得られるようにする
   // b!　これは，線形補間であっても利用できる
   //
   static_assert(N == 3 || N == 6, "Unsupported shape function size. Only 3 or 6 are supported.");
   auto t2 = 1 - t0 - t1;
   auto t0m1 = t0 - 1;
   auto t1m1 = t1 - 1;
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
      static const auto TriShape3 = TriShape<3>(static_cast<T>(1), static_cast<T>(1));
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

#endif