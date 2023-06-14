#ifndef basic_arithmetic_vector_operations_H
#define basic_arithmetic_vector_operations_H

#include <algorithm>  //transformなど
#include <cmath>
#include <numeric>  //数値のシーケンスの処理に特化したアルゴリズム
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
   // std::transform(v.begin(), v.end(), v.begin(), [&din](const T c) { return c * din; });
   for (auto &u : v)
      u *= din;
   return v;
};
//! matrix * scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator*(std::vector<std::vector<T>> v, const T din) {
   // std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &c) { return c * din; });
   for (auto &u : v)
      u *= din;
   return v;
};
//@ scaler * vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator*(const T din, std::vector<T> v) {
   // std::transform(v.begin(), v.end(), v.begin(), [&din](const T c) { return c * din; });
   for (auto &u : v)
      u *= din;
   return v;
};
//! scaler * matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator*(const T din, std::vector<std::vector<T>> v) {
   // std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &c) { return c * din; });
   for (auto &u : v)
      u *= din;
   return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<std::vector<T>>> operator*(std::vector<std::vector<std::vector<T>>> v, const T &din) {
   // std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<std::vector<T>> &c) { return c * din; });
   for (auto &u : v)
      u *= din;
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
// template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<T> operator-(std::vector<T> v, const T din) {
//    std::transform(v.begin(), v.end(), v.begin(), [&din](const T tmp) { return tmp - din; });
//    return v;
// };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(std::vector<T> v, const T din) {
   for (auto &u : v)
      u -= din;
   return v;
};

//! matrix - scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> v, const T &din) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp) { return tmp - din; });
   return v;
};
//@ scaler - vector
// template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<T> operator-(const T din, std::vector<T> v) {
//    std::transform(v.begin(), v.end(), v.begin(), [&din](const T tmp) { return din - tmp; });
//    return v;
// };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(const T din, std::vector<T> v) {
   for (auto &u : v)
      u = din - u;
   return v;
};

//! scaler - matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(const T din, std::vector<std::vector<T>> v) {
   std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp) { return din - tmp; });
   return v;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator-=(std::vector<T> &v, const std::vector<T> &w) {
   for (size_t i = 0; i < v.size(); ++i)
      v[i] -= w[i];
   return v;
};

//@ vector - vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(std::vector<T> v, const std::vector<T> &w) {
   return v -= w;
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
   std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b) { return b - a; });  // 逆にする
   return v;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator-=(std::vector<T> &v, const T &w) { return (v = v - w); };

// template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<T> &operator-=(std::vector<T> &v, const std::vector<T> &w) {
//    return (v = v - w);
// };

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

template <typename... TT>
struct hash<std::tuple<TT...>> {
   size_t operator()(std::tuple<TT...> const &tt) const {
      size_t seed = 0;
      std::apply([&](auto const &...args) {
         ((seed ^= hash<TT>()(args) + 0x9e3779b9 + (seed << 6) + (seed >> 2)), ...);
      },
                 tt);
      return seed;
   }
};

}  // namespace std

#include <unordered_map>

#endif