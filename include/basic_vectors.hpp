#ifndef basic_vectors_H
#define basic_vectors_H

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "basic_arithmetic_vector_operations.hpp"
#include "basic_exception.hpp"

/* -------------------------------------------------------------------------- */

template <>
void IdentityMatrix<VV_d>(VV_d &mat) {
   int i = 0, j = 0;
   for (auto &m : mat) {
      j = 0;
      for (auto &n : m) {
         if (i == j)
            n = 1.;
         else
            n = 0.;
         j++;
      }
      i++;
   }
};
/* -------------------------------------------------------------------------- */
template <typename T>
std::array<T, 2> Sort(std::array<T, 2> &ab) {
   if (std::get<0>(ab) <= std::get<1>(ab))
      return ab;
   else
      return {std::get<1>(ab), std::get<0>(ab)};
};
/* -------------------------------------------------------------------------- */
template <typename T>
void Swap(std::array<T, 2> &ab) {
   auto a = std::get<0>(ab);
   std::get<0>(ab) = std::get<1>(ab);
   std::get<1>(ab) = a;
};
/* -------------------------------------------------------------------------- */

// template <typename T>
// tuple_of<std::tuple_size<T>::value, tuple_of<std::tuple_size<T>::value, double>> Inverse(const tuple_of<std::tuple_size<T>::value, tuple_of<std::tuple_size<T>::value, double>> &M) {
//    tuple_of<std::tuple_size<T>::value, tuple_of<std::tuple_size<T>::value, double>> tup_mat;
//    std::vector<std::vector<T>> vec_mat(std::tuple_size<T>::value, std::vector<T>(std::tuple_size<T>::value));
//    int i = 0, j = 0;
//    for_each(tup_mat, [&](const auto &tup_mat_i) {
//       j = 0;
//       for_each(tup_mat_i, [&](const auto &tup_mat_ij) { vec_mat[i][j++] = tup_mat_ij; });
//       i++;
//    });
//    vec_mat = Inverse(vec_mat);
//    for_each(tup_mat, [&](const auto &tup_mat_i) {
//       j = 0;
//       for_each(tup_mat_i, [&](const auto &tup_mat_ij) { tup_mat_ij = vec_mat[i][j++]; });
//       i++;
//    });
//    return;
// };

// template <typename T>
// tuple_of<std::tuple_size<T>::value, tuple_of<std::tuple_size<T>::value, double>> Inverse(const tuple_of<std::tuple_size<T>::value, tuple_of<std::tuple_size<T>::value, double>> &M) {
//    tuple_of<std::tuple_size<T>::value, tuple_of<std::tuple_size<T>::value, double>> tup_mat;
//    std::vector<std::vector<T>> vec_mat(std::tuple_size<T>::value, std::vector<T>(std::tuple_size<T>::value));
//    int i = 0, j = 0;
//    for_each(tup_mat, [&](const auto &tup_mat_i) {
//       j = 0;
//       for_each(tup_mat_i, [&](const auto &tup_mat_ij) { vec_mat[i][j++] = tup_mat_ij; });
//       i++;
//    });
//    vec_mat = Inverse(vec_mat);
//    for_each(tup_mat, [&](const auto &tup_mat_i) {
//       j = 0;
//       for_each(tup_mat_i, [&](const auto &tup_mat_ij) { tup_mat_ij = vec_mat[i][j++]; });
//       i++;
//    });
//    return;
// };

/* -------------------------------------------------------------------------- */
T4T4d Inverse(const T4T4d &mat) {
   auto [x00, x01, x02, x03] = std::get<0>(mat);
   auto [x10, x11, x12, x13] = std::get<1>(mat);
   auto [x20, x21, x22, x23] = std::get<2>(mat);
   auto [x30, x31, x32, x33] = std::get<3>(mat);
   double det = (x01 * x13 * x22 * x30 - x01 * x12 * x23 * x30 - x00 * x13 * x22 * x31 + x00 * x12 * x23 * x31 - x01 * x13 * x20 * x32 + x00 * x13 * x21 * x32 + x01 * x10 * x23 * x32 - x00 * x11 * x23 * x32 + x03 * (x12 * x21 * x30 - x11 * x22 * x30 - x12 * x20 * x31 + x10 * x22 * x31 + x11 * x20 * x32 - x10 * x21 * x32) + (x01 * x12 * x20 - x00 * x12 * x21 - x01 * x10 * x22 + x00 * x11 * x22) * x33 + x02 * (-(x13 * x21 * x30) + x11 * x23 * x30 + x13 * x20 * x31 - x10 * x23 * x31 - x11 * x20 * x33 + x10 * x21 * x33));
   //
   T4T4d inv = {{{-(x13 * x22 * x31) + x12 * x23 * x31 + x13 * x21 * x32 - x11 * x23 * x32 - x12 * x21 * x33 + x11 * x22 * x33, x03 * x22 * x31 - x02 * x23 * x31 - x03 * x21 * x32 + x01 * x23 * x32 + x02 * x21 * x33 - x01 * x22 * x33, -(x03 * x12 * x31) + x02 * x13 * x31 + x03 * x11 * x32 - x01 * x13 * x32 - x02 * x11 * x33 + x01 * x12 * x33, x03 * x12 * x21 - x02 * x13 * x21 - x03 * x11 * x22 + x01 * x13 * x22 + x02 * x11 * x23 - x01 * x12 * x23},
                 {x13 * x22 * x30 - x12 * x23 * x30 - x13 * x20 * x32 + x10 * x23 * x32 + x12 * x20 * x33 - x10 * x22 * x33, -(x03 * x22 * x30) + x02 * x23 * x30 + x03 * x20 * x32 - x00 * x23 * x32 - x02 * x20 * x33 + x00 * x22 * x33, x03 * x12 * x30 - x02 * x13 * x30 - x03 * x10 * x32 + x00 * x13 * x32 + x02 * x10 * x33 - x00 * x12 * x33, -(x03 * x12 * x20) + x02 * x13 * x20 + x03 * x10 * x22 - x00 * x13 * x22 - x02 * x10 * x23 + x00 * x12 * x23},
                 {-(x13 * x21 * x30) + x11 * x23 * x30 + x13 * x20 * x31 - x10 * x23 * x31 - x11 * x20 * x33 + x10 * x21 * x33, x03 * x21 * x30 - x01 * x23 * x30 - x03 * x20 * x31 + x00 * x23 * x31 + x01 * x20 * x33 - x00 * x21 * x33, -(x03 * x11 * x30) + x01 * x13 * x30 + x03 * x10 * x31 - x00 * x13 * x31 - x01 * x10 * x33 + x00 * x11 * x33, x03 * x11 * x20 - x01 * x13 * x20 - x03 * x10 * x21 + x00 * x13 * x21 + x01 * x10 * x23 - x00 * x11 * x23},
                 {x12 * x21 * x30 - x11 * x22 * x30 - x12 * x20 * x31 + x10 * x22 * x31 + x11 * x20 * x32 - x10 * x21 * x32, -(x02 * x21 * x30) + x01 * x22 * x30 + x02 * x20 * x31 - x00 * x22 * x31 - x01 * x20 * x32 + x00 * x21 * x32, x02 * x11 * x30 - x01 * x12 * x30 - x02 * x10 * x31 + x00 * x12 * x31 + x01 * x10 * x32 - x00 * x11 * x32, -(x02 * x11 * x20) + x01 * x12 * x20 + x02 * x10 * x21 - x00 * x12 * x21 - x01 * x10 * x22 + x00 * x11 * x22}}};
   return inv / det;
};

double Det(const T3Tddd &M) {
   return -(std::get<0>(std::get<2>(M)) * std::get<1>(std::get<1>(M)) * std::get<2>(std::get<0>(M))) +
          std::get<0>(std::get<1>(M)) * std::get<1>(std::get<2>(M)) * std::get<2>(std::get<0>(M)) +
          std::get<0>(std::get<2>(M)) * std::get<1>(std::get<0>(M)) * std::get<2>(std::get<1>(M)) -
          std::get<0>(std::get<0>(M)) * std::get<1>(std::get<2>(M)) * std::get<2>(std::get<1>(M)) -
          std::get<0>(std::get<1>(M)) * std::get<1>(std::get<0>(M)) * std::get<2>(std::get<2>(M)) +
          std::get<0>(std::get<0>(M)) * std::get<1>(std::get<1>(M)) * std::get<2>(std::get<2>(M));
};
/* ------------------------------------------------------ */
template <typename T>
bool IntersectingQ(const std::vector<T> &A, const std::vector<T> &B) {
   for (const auto &a : A)
      for (const auto &b : B)
         if (b == a)
            return true;
   return false;
};
template <typename T>
bool IntersectingQ(const std::unordered_set<T> &A, const std::unordered_set<T> &B) {
   for (const auto &a : A)
      for (const auto &b : B)
         if (b == a)
            return true;
   return false;
};
template <typename T>
bool IntersectingQ(const std::unordered_set<T> &A, const std::vector<T> &B) {
   for (const auto &a : A)
      for (const auto &b : B)
         if (b == a)
            return true;
   return false;
};
template <typename T>
bool IntersectingQ(const std::vector<T> &A, const std::unordered_set<T> &B) {
   for (const auto &a : A)
      for (const auto &b : B)
         if (b == a)
            return true;
   return false;
};
template <typename T>
bool IntersectingQ(const std::unordered_set<T> &A, const T &b) {
   for (const auto &a : A)
      if (b == a)
         return true;
   return false;
};
template <typename T>
bool IntersectingQ(const std::vector<T> &A, const T &b) {
   for (const auto &a : A)
      if (b == a)
         return true;
   return false;
};
/* ------------------------------------------------------ */

template <typename T>
class Chain {
  public:
   int isComplete;
   int isReversed;
   int addedBack;
   int addedFront;
   std::vector<T *> chained;
   Chain() : chained({}) {
      addedBack = false;
      addedFront = false;
      isComplete = false;
      isReversed = false;
   };
   ////////////////
   void display() {
      std::cout << "--------------------------------------" << std::endl;
      std::cout << "isComplete = " << isComplete << std::endl;
      std::cout << "isReversed = " << isReversed << std::endl;
      std::cout << " addedBack = " << addedBack << std::endl;
      std::cout << " addedFront = " << addedFront << std::endl;
      std::cout << "   chained = " << chained << std::endl;
   };
   ///////
   void checkComplete(bool oneside_match) {
      if (oneside_match) {
         // A--B--Aでは鎖とは言えない
         // A--B--C--AならOK
         if ((*chained.begin()) == (*chained.rbegin()) && chained.size() > 3) {
            isComplete = 2;
         } else
            isComplete = 1;
      } else {
         this->chained = {};
         isComplete = 0;
         isReversed = 0;
      }
   };
   /////
   void join_front(const std::vector<T *> &base, const std::vector<T *> &pc) {
      this->chained = base;
      bool oneside_match = false;
      if ((*chained.rbegin()) == (*pc.begin())) {  // right direction
         chained.insert(chained.end(), /*後尾に挿入*/ pc.begin() + 1, pc.end());
         oneside_match = true;
         isReversed = false;
         addedBack = true;
         addedFront = false;
      } else if (*(chained.rbegin()) == *pc.rbegin()) {  // oppsite direction
         chained.insert(chained.end(), /*後尾に反転して挿入*/ pc.rbegin() + 1, pc.rend());
         oneside_match = true;
         isReversed = true;
         addedBack = true;
         addedFront = false;
      }
      // else if (*(chained.begin()) == *pc.rbegin())
      // { //right direction
      //   chained.insert(chained.begin(), /*前方に挿入*/ pc.begin(), pc.end() - 1);
      //   oneside_match = true;
      //   isReversed = false;
      //   addedBack = false;
      //   addedFront = true;
      // }
      // else if (*(chained.begin()) == *pc.begin())
      // { //oppsite direction
      //   chained.insert(chained.begin(), /*前方に反転して挿入*/ pc.rbegin(), pc.rend() - 1);
      //   oneside_match = true;
      //   isReversed = true;
      //   addedBack = false;
      //   addedFront = true;
      // }
      checkComplete(oneside_match);
   };
   ////////////////
   void join_front_fix_order(const std::vector<T *> &base, const std::vector<T *> &pc) {
      this->chained = base;
      bool oneside_match = false;
      if ((*chained.rbegin()) == (*pc.begin())) {  // right direction
         chained.insert(chained.end(), /*後尾に挿入*/ pc.begin() + 1, pc.end());
         oneside_match = true;
         isReversed = false;
         addedBack = true;
         addedFront = false;
      }
      // else if (*(chained.rbegin()) == *pc.rbegin())
      // { //oppsite direction
      //   chained.insert(chained.end(), /*後尾に反転して挿入*/ pc.rbegin() + 1, pc.rend());
      //   oneside_match = true;
      //   isReversed = true;
      //   addedBack = true;
      //   addedFront = false;
      // }
      // else if (*(chained.begin()) == *pc.rbegin())
      // { //right direction
      //   chained.insert(chained.begin(), /*前方に挿入*/ pc.begin(), pc.end() - 1);
      //   oneside_match = true;
      //   isReversed = false;
      //   addedBack = false;
      //   addedFront = true;
      // }
      // else if (*(chained.begin()) == *pc.begin())
      // { //oppsite direction
      //   chained.insert(chained.begin(), /*前方に反転して挿入*/ pc.rbegin(), pc.rend() - 1);
      //   oneside_match = true;
      //   isReversed = true;
      //   addedBack = false;
      //   addedFront = true;
      // }
      checkComplete(oneside_match);
   };
   // 前後ろに適当にくっつけていくのは，うまくいかないことがわかった，
   // なぜなら，cutlineどうしがくっついたりする場合がでるためだ，
   void join(const std::vector<T *> &base, const std::vector<T *> &pc) {
      this->chained = base;
      bool oneside_match = false;
      if ((*chained.rbegin()) == (*pc.begin())) {  // right direction
         chained.insert(chained.end(), /*後尾に挿入*/ pc.begin() + 1, pc.end());
         oneside_match = true;
         isReversed = false;
         addedBack = true;
         addedFront = false;
      } else if (*(chained.rbegin()) == *pc.rbegin()) {  // oppsite direction
         chained.insert(chained.end(), /*後尾に反転して挿入*/ pc.rbegin() + 1, pc.rend());
         oneside_match = true;
         isReversed = true;
         addedBack = true;
         addedFront = false;
      } else if (*(chained.begin()) == *pc.rbegin()) {  // right direction
         chained.insert(chained.begin(), /*前方に挿入*/ pc.begin(), pc.end() - 1);
         oneside_match = true;
         isReversed = false;
         addedBack = false;
         addedFront = true;
      } else if (*(chained.begin()) == *pc.begin()) {  // oppsite direction
         chained.insert(chained.begin(), /*前方に反転して挿入*/ pc.rbegin(), pc.rend() - 1);
         oneside_match = true;
         isReversed = true;
         addedBack = false;
         addedFront = true;
      }
      checkComplete(oneside_match);
   };
};

///////////////////////////////////////
template <typename T>
std::vector<T *> Flatten(const std::vector<std::vector<T *>> &mat) {
   // std::vector<T *> ret;
   // for (const auto &m : mat)
   // 	for (const auto &n : m)
   // 		ret.emplace_back(n);
   // return ret;
   std::vector<T *> ret;
   for (const auto &part : mat)
      ret.insert(ret.end(), part.cbegin(), part.cend());
   return ret;
};
template <typename T>
std::vector<T> Flatten(const std::vector<std::vector<T>> &mat) {
   std::vector<T> ret;
   for (const auto &part : mat)
      ret.insert(ret.end(), part.cbegin(), part.cend());
   return ret;
};
template <typename T>
std::vector<T> Flatten(const std::vector<std::unordered_set<T>> &mat) {
   std::vector<T> ret(0);
   ret.reserve(mat.size() * mat[0].size());
   for (const auto &m : mat)
      for (const auto &n : m)
         ret.emplace_back(n);
   return ret;
};
std::vector<Tddd> Flatten(const std::vector<std::vector<Tddd>> &mat) {
   std::vector<Tddd> ret(0);
   ret.reserve(mat.size() * mat[0].size());
   for (const auto &m : mat)
      for (const auto &n : m)
         ret.emplace_back(n);
   return ret;
};
// 鎖状につなげていく．
// 1つ目のベクトルの順番を入れ替えることもあり得る．
// 2つ目のベクトルは，同じ方向を向いていなければならない．
template <typename T>
std::vector<T> FlattenAsChain(/*not const*/ std::vector<std::vector<T>> Vps) {
   int count = 0;  // 無限ループをさける
   std::vector<T> ret = Vps[0];
   Vps.erase(Vps.begin());
   while (!Vps.empty()) {
      for (auto it = std::begin(Vps); it < std::end(Vps); it++) {
         // 前か後につぎつぎにくっつけていく
         if (*ret.rbegin() == *std::begin(*it)) {
            ret.emplace_back(*std::rbegin(*it) /*next*/);
            Vps.erase(it);
            break;
         } else if (*ret.begin() == *std::rbegin(*it)) {
            ret.insert(ret.begin(), *std::begin(*it) /*next*/);
            Vps.erase(it);
            break;
         }
         if (it == std::end(Vps) - 1)
            return ret;
      }

      if (count++ > 10000)
         throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "loop infinity　多分，点の周りの面が重複したりして，数珠繋ぎにできない"));
   };
   return ret;
};
//////////////////////////////////////
template <typename T>
std::tuple<T, T, T, T> Flatten(const std::tuple<std::tuple<T, T>, std::tuple<T, T>> &v) {
   return {std::get<0>(std::get<0>(v)), std::get<1>(std::get<0>(v)),
           std::get<0>(std::get<1>(v)), std::get<1>(std::get<1>(v))};
};
template <typename T>
std::tuple<T, T, T, T, T, T> Flatten(const std::tuple<std::tuple<T, T, T>, std::tuple<T, T, T>> &v) {
   return {std::get<0>(std::get<0>(v)), std::get<1>(std::get<0>(v)), std::get<2>(std::get<0>(v)),
           std::get<0>(std::get<1>(v)), std::get<1>(std::get<1>(v)), std::get<2>(std::get<1>(v))};
};
template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T> Flatten(const std::tuple<std::tuple<T, T, T>, std::tuple<T, T, T>, std::tuple<T, T, T>> &v) {
   return {std::get<0>(std::get<0>(v)), std::get<1>(std::get<0>(v)), std::get<2>(std::get<0>(v)),
           std::get<0>(std::get<1>(v)), std::get<1>(std::get<1>(v)), std::get<2>(std::get<1>(v)),
           std::get<0>(std::get<2>(v)), std::get<1>(std::get<2>(v)), std::get<2>(std::get<2>(v))};
};
template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T, T> Flatten(const std::tuple<std::tuple<T, T, T>, std::tuple<T, T, T>, std::tuple<T, T, T>, std::tuple<T, T, T>> &v) {
   return {std::get<0>(std::get<0>(v)), std::get<1>(std::get<0>(v)), std::get<2>(std::get<0>(v)),
           std::get<0>(std::get<1>(v)), std::get<1>(std::get<1>(v)), std::get<2>(std::get<1>(v)),
           std::get<0>(std::get<2>(v)), std::get<1>(std::get<2>(v)), std::get<2>(std::get<2>(v)),
           std::get<0>(std::get<3>(v)), std::get<1>(std::get<3>(v)), std::get<2>(std::get<3>(v))};
};
template <typename T>
std::tuple<T, T, T,
           T, T, T,
           T, T, T,
           T, T, T,
           T, T, T,
           T, T, T,
           T, T, T,
           T, T, T>
Flatten(const std::tuple<std::tuple<T, T, T, T>,
                         std::tuple<T, T, T, T>,
                         std::tuple<T, T, T, T>,
                         std::tuple<T, T, T, T>,
                         std::tuple<T, T, T, T>,
                         std::tuple<T, T, T, T>> &v) {
   return {std::get<0>(std::get<0>(v)), std::get<1>(std::get<0>(v)), std::get<2>(std::get<0>(v)), std::get<3>(std::get<0>(v)),
           std::get<0>(std::get<1>(v)), std::get<1>(std::get<1>(v)), std::get<2>(std::get<1>(v)), std::get<3>(std::get<1>(v)),
           std::get<0>(std::get<2>(v)), std::get<1>(std::get<2>(v)), std::get<2>(std::get<2>(v)), std::get<3>(std::get<2>(v)),
           std::get<0>(std::get<3>(v)), std::get<1>(std::get<3>(v)), std::get<2>(std::get<3>(v)), std::get<3>(std::get<3>(v)),
           std::get<0>(std::get<4>(v)), std::get<1>(std::get<4>(v)), std::get<2>(std::get<4>(v)), std::get<3>(std::get<4>(v)),
           std::get<0>(std::get<5>(v)), std::get<1>(std::get<5>(v)), std::get<2>(std::get<5>(v)), std::get<3>(std::get<5>(v))};
};
template <typename T>
std::tuple<T, T, T,
           T, T, T,
           T, T, T,
           T, T, T,
           T, T, T,
           T, T, T,
           T, T, T,
           T, T, T>
Flatten(const std::tuple<std::tuple<T, T, T>,
                         std::tuple<T, T, T>,
                         std::tuple<T, T, T>,
                         std::tuple<T, T, T>,
                         std::tuple<T, T, T>,
                         std::tuple<T, T, T>,
                         std::tuple<T, T, T>,
                         std::tuple<T, T, T>> &v) {
   return {std::get<0>(std::get<0>(v)), std::get<1>(std::get<0>(v)), std::get<2>(std::get<0>(v)),
           std::get<0>(std::get<1>(v)), std::get<1>(std::get<1>(v)), std::get<2>(std::get<1>(v)),
           std::get<0>(std::get<2>(v)), std::get<1>(std::get<2>(v)), std::get<2>(std::get<2>(v)),
           std::get<0>(std::get<3>(v)), std::get<1>(std::get<3>(v)), std::get<2>(std::get<3>(v)),
           std::get<0>(std::get<4>(v)), std::get<1>(std::get<4>(v)), std::get<2>(std::get<4>(v)),
           std::get<0>(std::get<5>(v)), std::get<1>(std::get<5>(v)), std::get<2>(std::get<5>(v)),
           std::get<0>(std::get<6>(v)), std::get<1>(std::get<6>(v)), std::get<2>(std::get<6>(v)),
           std::get<0>(std::get<7>(v)), std::get<1>(std::get<7>(v)), std::get<2>(std::get<7>(v))};
};
//////////////////////////////////////
template <class T>
std::vector<std::vector<T>> Transpose(const std::vector<std::vector<T>> &mat) {
   if (mat.empty())
      return mat;
   std::vector<std::vector<T>> ans(mat[0].size(), std::vector<T>(mat.size()));
   for (size_t i = 0; i < mat.size(); i++)
      for (size_t j = 0; j < mat[i].size(); j++)
         ans[j][i] = mat[i][j];
   return ans;
};

VVV_d Transpose(const std::vector<VV_d> &mat) {
   VVV_d ans(mat[0][0].size(), VV_d(mat[0].size(), V_d(mat.size())));
   for (size_t i = 0; i < mat.size(); i++)
      for (size_t j = 0; j < mat[i].size(); j++)
         for (size_t k = 0; k < mat[i][j].size(); k++)
            ans[k][j][i] = mat[i][j][k];
   return ans;
};

VV_d TensorProduct(const V_d &vec1, const V_d &vec2) {
   VV_d ret(vec1.size(), V_d(vec2.size()));
   for (auto m = 0; m < vec1.size(); m++)
      for (auto j = 0; j < vec2.size(); j++)
         ret[m][j] = vec1[m] * vec2[j];
   return ret;
};

VVV_d TensorProductSet(const V_d &vec1, const V_d &vec2) {
   VVV_d ret(vec1.size(), VV_d(vec2.size(), V_d(2, 0)));
   for (size_t m = 0; m < vec1.size(); m++)
      for (size_t j = 0; j < vec2.size(); j++)
         ret[m][j] = {vec1[m], vec2[j]};
   return ret;
};

template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Dot(const std::vector<T> &vec1, const std::vector<T> &vec2) {
   T ret = 0;
   for (size_t i = 0; const auto &v1 : vec1) {
      ret = std::fma(v1, vec2[i++], ret);
   }
   return ret;
}

template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> Dot(const std::vector<std::vector<T>> &mat, const std::vector<T> &vec) {
   std::vector<T> ans(mat.size());
   int i = 0;
   for (const auto &m : mat)
      ans[i++] = Dot(m, vec);
   return ans;
};

template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> Dot(const std::vector<T> &vec, const std::vector<std::vector<T>> &mat) {
   std::vector<T> ans(mat[0].size(), 0.);
   int j = 0, i = 0;
   T tmp;
   for (const auto &mj : mat) {
      tmp = vec[j];
      i = 0;
      for (const auto &mi : mj) {
         ans[i] = std::fma(mi, tmp, ans[i]);
         i++;
      }
      j++;
   }
   return ans;
}

template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> Dot(const std::vector<std::vector<T>> &mat1, const std::vector<std::vector<T>> &mat2) {
   std::vector<std::vector<T>> ans(mat1.size(), std::vector<T>(mat2[0].size(), 0.));
   for (size_t x = 0; x < mat1.size(); x++) {
      for (size_t y = 0; y < mat2[0].size(); y++) {
         for (size_t j = 0; j < mat2.size(); j++) {
            ans[x][y] = std::fma(mat1[x][j], mat2[j][y], ans[x][y]);
         }
      }
   }
   return ans;
}

template <typename T>
std::vector<T> ToVector(const std::unordered_set<T> &uo) {
   return std::vector<T>(uo.begin(), uo.end());
};
template <typename T>
std::unordered_set<T> ToUnorderedSet(const std::vector<T> &v) { return std::unordered_set<T>(v.begin(), v.end()); };

template <typename T>
std::vector<T> ToVector(const std::tuple<T, T> &v) {
   return {std::get<0>(v), std::get<1>(v)};
};

template <typename T>
std::vector<T> ToVector(const std::tuple<T, T, T> &v) {
   return {std::get<0>(v), std::get<1>(v), std::get<2>(v)};
};

template <typename T>
std::vector<T> ToVector(const std::tuple<T, T, T, T, T, T, T> &v) { return {std::get<0>(v),
                                                                            std::get<1>(v),
                                                                            std::get<2>(v),
                                                                            std::get<3>(v),
                                                                            std::get<4>(v),
                                                                            std::get<5>(v),
                                                                            std::get<6>(v)}; };
template <typename T>
std::vector<T> ToVector(const std::tuple<T, T, T, T, T, T, T, T> &v) { return {std::get<0>(v),
                                                                               std::get<1>(v),
                                                                               std::get<2>(v),
                                                                               std::get<3>(v),
                                                                               std::get<4>(v),
                                                                               std::get<5>(v),
                                                                               std::get<6>(v),
                                                                               std::get<7>(v)}; };

std::vector<double> ToVector(const Tddd &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v)}; };
std::vector<double> ToVector(const T4d &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v)}; };
std::vector<double> ToVector(const T6d &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v), std::get<4>(v), std::get<5>(v)}; };
std::vector<double> ToVector(const T7d &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v), std::get<4>(v), std::get<5>(v), std::get<6>(v)}; };

VV_d ToVector(const std::vector<Tddd> &v) {
   VV_d ret(v.size(), {0, 0, 0});
   for (auto i = 0; i < v.size(); ++i)
      ret[i] = ToVector(v[i]);
   return ret;
};
VV_d ToVector(const T3Tddd &v) {
   return {ToVector(std::get<0>(v)),
           ToVector(std::get<1>(v)),
           ToVector(std::get<2>(v))};
};
VV_d ToVector(const T4T4d &v) {
   return {ToVector(std::get<0>(v)),
           ToVector(std::get<1>(v)),
           ToVector(std::get<2>(v)),
           ToVector(std::get<3>(v))};
};
VV_d ToVector(const T6T6d &v) {
   return {ToVector(std::get<0>(v)),
           ToVector(std::get<1>(v)),
           ToVector(std::get<2>(v)),
           ToVector(std::get<3>(v)),
           ToVector(std::get<4>(v)),
           ToVector(std::get<5>(v))};
};
VV_d ToVector(const T7T7d &v) {
   return {ToVector(std::get<0>(v)),
           ToVector(std::get<1>(v)),
           ToVector(std::get<2>(v)),
           ToVector(std::get<3>(v)),
           ToVector(std::get<4>(v)),
           ToVector(std::get<5>(v)),
           ToVector(std::get<6>(v))};
};
T6d ToT6d(const Tddd tmp) {
   return {std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), 0., 0., 0.};
};
Tddd ToTddd(const V_d &v) { return {v[0], v[1], v[2]}; };
Tddd ToTddd(const T6d &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v)}; };
Tdd ToTdd(const V_d &v) { return {v[0], v[1]}; };
// std::vector<double> ToVector(const Tdd &v) { return {std::get<0>(v), std::get<1>(v)}; };
// double Norm(const T4d &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) + std::get<1>(t) * std::get<1>(t) + std::get<2>(t) * std::get<2>(t) + std::get<3>(t) * std::get<3>(t)); };
// double Norm(const T6d &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) +
//                                              std::get<1>(t) * std::get<1>(t) +
//                                              std::get<2>(t) * std::get<2>(t) +
//                                              std::get<3>(t) * std::get<3>(t) +
//                                              std::get<4>(t) * std::get<4>(t) +
//                                              std::get<5>(t) * std::get<5>(t)); };
// double Norm(const T7d &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) +
//                                              std::get<1>(t) * std::get<1>(t) +
//                                              std::get<2>(t) * std::get<2>(t) +
//                                              std::get<3>(t) * std::get<3>(t) +
//                                              std::get<4>(t) * std::get<4>(t) +
//                                              std::get<5>(t) * std::get<5>(t) +
//                                              std::get<6>(t) * std::get<6>(t)); };
// double Norm(const Tddd &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) + std::get<1>(t) * std::get<1>(t) + std::get<2>(t) * std::get<2>(t)); };
// double Norm(const Tdd &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) + std::get<1>(t) * std::get<1>(t)); };
// T4d Normalize(const T4d &X) { return X / Norm(X); };
// T7d Normalize(const T7d &X) { return X / Norm(X); };
// Tddd Normalize(const Tddd &X) { return X / Norm(X); };
// Tdd Normalize(const Tdd &X) { return X / Norm(X); };

/* -------------------------------------------------------------------------- */

// template <typename T>
// std::vector<T> Dot(const std::vector<std::vector<double>> &A, const std::vector<T> &B) {
//    // (M x N) . (N x 3) = (M x 3)
//    const int N = std::tuple_size<T>::value;
//    std::vector<T> ret(A.size());
//    for (auto i = 0; i < A.size(); ++i) {
//       for_each(ret[i], [](auto &r) { r = 0; });
//       for (auto j = 0; j < N; ++j)
//          std::get<j>(ret[i]) += A[i][j] * std::get<j>(B[i]);
//    }
//    return ret;
// };

// template <typename T, typename = std::enable_if_t<std::is_same_v<T, std::tuple<double, double>> ||
//                                                   std::is_same_v<T, std::tuple<double, double, double>> ||
//                                                   std::is_same_v<T, std::tuple<double, double, double, double>> ||
//                                                   std::is_same_v<T, std::tuple<double, double, double, double, double>> ||
//                                                   std::is_same_v<T, std::tuple<double, double, double, double, double, double>> ||
//                                                   std::is_same_v<T, std::tuple<double, double, double, double, double, double, double>> ||
//                                                   std::is_same_v<T, std::tuple<double, double, double, double, double, double, double, double>> ||
//                                                   std::is_same_v<T, std::tuple<double, double, double, double, double, double, double, double, double>>>>
// std::vector<T> Dot(const std::vector<std::vector<double>> &A, const std::vector<T> &B) {
//    const int N = std::tuple_size<T>::value;
//    std::vector<T> ret(A.size());
//    int j = 0;
//    for (auto i = 0; i < A.size(); ++i) {
//       for_each(ret[i], [&](auto &r) { r = 0; });
//       j = 0;
//       for_each10(ret[i], B[i], [&](auto &r, auto &b) { r += A[i][j++] * b; });
//    }
//    return ret;
// };

// double Dot(const T6d &v, const T6d &u) {
//    return (std::get<0>(v) * std::get<0>(u) +
//            std::get<1>(v) * std::get<1>(u) +
//            std::get<2>(v) * std::get<2>(u) +
//            std::get<3>(v) * std::get<3>(u) +
//            std::get<4>(v) * std::get<4>(u) +
//            std::get<5>(v) * std::get<5>(u));
// };
// Tddd Dot(const T6d &v, const T6Tddd &A) {
//    return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<0>(std::get<1>(A)) * std::get<1>(v) + std::get<0>(std::get<2>(A)) * std::get<2>(v) + std::get<0>(std::get<3>(A)) * std::get<3>(v) + std::get<0>(std::get<4>(A)) * std::get<4>(v) + std::get<0>(std::get<5>(A)) * std::get<5>(v),
//            std::get<1>(std::get<0>(A)) * std::get<0>(v) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<1>(std::get<2>(A)) * std::get<2>(v) + std::get<1>(std::get<3>(A)) * std::get<3>(v) + std::get<1>(std::get<4>(A)) * std::get<4>(v) + std::get<1>(std::get<5>(A)) * std::get<5>(v),
//            std::get<2>(std::get<0>(A)) * std::get<0>(v) + std::get<2>(std::get<1>(A)) * std::get<1>(v) + std::get<2>(std::get<2>(A)) * std::get<2>(v) + std::get<2>(std::get<3>(A)) * std::get<3>(v) + std::get<2>(std::get<4>(A)) * std::get<4>(v) + std::get<2>(std::get<5>(A)) * std::get<5>(v)};
// };
// double Dot(const T4d &v, const T4d &u) {
//    return (std::get<0>(v) * std::get<0>(u) +
//            std::get<1>(v) * std::get<1>(u) +
//            std::get<2>(v) * std::get<2>(u) +
//            std::get<3>(v) * std::get<3>(u));
// };
// double Dot(const Tddd &v, const Tddd &u) {
//    return (std::get<0>(v) * std::get<0>(u) +
//            std::get<1>(v) * std::get<1>(u) +
//            std::get<2>(v) * std::get<2>(u));
// };

// double Dot(const Tdd &v, const Tdd &u) {
//    return (std::get<0>(v) * std::get<0>(u) + std::get<1>(v) * std::get<1>(u));
// };

// Tdd Dot(const Tdd &t, const T2Tdd &u) {
//    return {std::get<0>(t) * std::get<0>(std::get<0>(u)) + std::get<1>(t) * std::get<0>(std::get<1>(u)),
//            std::get<0>(t) * std::get<1>(std::get<0>(u)) + std::get<1>(t) * std::get<1>(std::get<1>(u))};
// };

// Tdd Dot(const T2Tdd &u, const Tdd &t) {
//    return {std::get<0>(t) * std::get<0>(std::get<0>(u)) + std::get<1>(t) * std::get<1>(std::get<0>(u)),
//            std::get<0>(t) * std::get<0>(std::get<1>(u)) + std::get<1>(t) * std::get<1>(std::get<1>(u))};
// };

// Tddd Dot(const Tdd &t, const T2Tddd &u) {
//    return {std::get<0>(t) * std::get<0>(std::get<0>(u)) + std::get<1>(t) * std::get<0>(std::get<1>(u)),
//            std::get<0>(t) * std::get<1>(std::get<0>(u)) + std::get<1>(t) * std::get<1>(std::get<1>(u)),
//            std::get<0>(t) * std::get<2>(std::get<0>(u)) + std::get<1>(t) * std::get<2>(std::get<1>(u))};
// };

// Tddd Dot(const T3Tddd &A, const Tddd &v) {
//    return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<1>(std::get<0>(A)) * std::get<1>(v) + std::get<2>(std::get<0>(A)) * std::get<2>(v),
//            std::get<0>(std::get<1>(A)) * std::get<0>(v) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<2>(std::get<1>(A)) * std::get<2>(v),
//            std::get<0>(std::get<2>(A)) * std::get<0>(v) + std::get<1>(std::get<2>(A)) * std::get<1>(v) + std::get<2>(std::get<2>(A)) * std::get<2>(v)};
// };

// /*
// In[12]:= Dot[{t0,t1,t2},{{X0x,X0y},{X1x,X1y},{X2x,X2y}}]
// Out[12]= {t0 X0x+t1 X1x+t2 X2x,t0 X0y+t1 X1y+t2 X2y
// */
// Tdd Dot(const Tddd &v, const T3Tdd &A) {
//    return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<0>(std::get<1>(A)) * std::get<1>(v) + std::get<0>(std::get<2>(A)) * std::get<2>(v),
//            std::get<1>(std::get<0>(A)) * std::get<0>(v) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<1>(std::get<2>(A)) * std::get<2>(v)};
// };

// Tddd Dot(const Tddd &v, const T3Tddd &A) {
//    return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<0>(std::get<1>(A)) * std::get<1>(v) + std::get<0>(std::get<2>(A)) * std::get<2>(v),
//            std::get<1>(std::get<0>(A)) * std::get<0>(v) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<1>(std::get<2>(A)) * std::get<2>(v),
//            std::get<2>(std::get<0>(A)) * std::get<0>(v) + std::get<2>(std::get<1>(A)) * std::get<1>(v) + std::get<2>(std::get<2>(A)) * std::get<2>(v)};
// };

// T4d Dot(const T4T4d &A, const T4d &v) {
//    return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<1>(std::get<0>(A)) * std::get<1>(v) + std::get<2>(std::get<0>(A)) * std::get<2>(v) + std::get<3>(std::get<0>(A)) * std::get<3>(v),
//            std::get<0>(std::get<1>(A)) * std::get<0>(v) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<2>(std::get<1>(A)) * std::get<2>(v) + std::get<3>(std::get<1>(A)) * std::get<3>(v),
//            std::get<0>(std::get<2>(A)) * std::get<0>(v) + std::get<1>(std::get<2>(A)) * std::get<1>(v) + std::get<2>(std::get<2>(A)) * std::get<2>(v) + std::get<3>(std::get<2>(A)) * std::get<3>(v),
//            std::get<0>(std::get<3>(A)) * std::get<0>(v) + std::get<1>(std::get<3>(A)) * std::get<1>(v) + std::get<2>(std::get<3>(A)) * std::get<2>(v) + std::get<3>(std::get<3>(A)) * std::get<3>(v)};
// };

// T4d Dot(const T4d &v, const T4T4d &A) {
//    return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<0>(std::get<1>(A)) * std::get<1>(v) + std::get<0>(std::get<2>(A)) * std::get<2>(v) + std::get<0>(std::get<3>(A)) * std::get<3>(v),
//            std::get<0>(v) * std::get<1>(std::get<0>(A)) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<1>(std::get<2>(A)) * std::get<2>(v) + std::get<1>(std::get<3>(A)) * std::get<3>(v),
//            std::get<0>(v) * std::get<2>(std::get<0>(A)) + std::get<1>(v) * std::get<2>(std::get<1>(A)) + std::get<2>(std::get<2>(A)) * std::get<2>(v) + std::get<2>(std::get<3>(A)) * std::get<3>(v),
//            std::get<0>(v) * std::get<3>(std::get<0>(A)) + std::get<1>(v) * std::get<3>(std::get<1>(A)) + std::get<2>(v) * std::get<3>(std::get<2>(A)) + std::get<3>(std::get<3>(A)) * std::get<3>(v)};
// };

// T5d Dot(const T5d &v, const T5T5d &A) {
//    return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<0>(std::get<1>(A)) * std::get<1>(v) + std::get<0>(std::get<2>(A)) * std::get<2>(v) + std::get<0>(std::get<3>(A)) * std::get<3>(v) + std::get<0>(std::get<4>(A)) * std::get<4>(v),
//            std::get<0>(v) * std::get<1>(std::get<0>(A)) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<1>(std::get<2>(A)) * std::get<2>(v) + std::get<1>(std::get<3>(A)) * std::get<3>(v) + std::get<1>(std::get<4>(A)) * std::get<4>(v),
//            std::get<0>(v) * std::get<2>(std::get<0>(A)) + std::get<1>(v) * std::get<2>(std::get<1>(A)) + std::get<2>(std::get<2>(A)) * std::get<2>(v) + std::get<2>(std::get<3>(A)) * std::get<3>(v) + std::get<2>(std::get<4>(A)) * std::get<4>(v),
//            std::get<0>(v) * std::get<3>(std::get<0>(A)) + std::get<1>(v) * std::get<3>(std::get<1>(A)) + std::get<2>(v) * std::get<3>(std::get<2>(A)) + std::get<3>(std::get<3>(A)) * std::get<3>(v) + std::get<3>(std::get<4>(A)) * std::get<4>(v),
//            std::get<0>(v) * std::get<4>(std::get<0>(A)) + std::get<1>(v) * std::get<4>(std::get<1>(A)) + std::get<2>(v) * std::get<4>(std::get<2>(A)) + std::get<3>(v) * std::get<4>(std::get<3>(A)) + std::get<4>(std::get<4>(A)) * std::get<4>(v)};
// };

/* -------------------------------------------------------------------------- */

// Tddd Cross(const Tddd &A, const Tddd &X) {
//    return {{std::get<1>(A) * std::get<2>(X) - std::get<2>(A) * std::get<1>(X),
//             std::get<2>(A) * std::get<0>(X) - std::get<0>(A) * std::get<2>(X),
//             std::get<0>(A) * std::get<1>(X) - std::get<1>(A) * std::get<0>(X)}};
// };

// Tddd Cross(const Tdd &A, const Tdd &X) {
//    return {{0., 0., std::get<0>(A) * std::get<1>(X) - std::get<1>(A) * std::get<0>(X)}};
// };

// T2Tdd Transpose(const T2Tdd &A) {
//    return {{{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A))},
//             {std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A))}}};
// };
// T2Tddd Transpose(const T3Tdd &A) {
//    return {{{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A)), std::get<0>(std::get<2>(A))},
//             {std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A)), std::get<1>(std::get<2>(A))}}};
// };

// T3Tdd Transpose(const T2Tddd &A) {
//    return {{{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A))},
//             {std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A))},
//             {std::get<2>(std::get<0>(A)), std::get<2>(std::get<1>(A))}}};
// };

// T3Tddd Transpose(const T3Tddd &A) {
//    return {{{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A)), std::get<0>(std::get<2>(A))},
//             {std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A)), std::get<1>(std::get<2>(A))},
//             {std::get<2>(std::get<0>(A)), std::get<2>(std::get<1>(A)), std::get<2>(std::get<2>(A))}}};
// };

// T3T4d Transpose(const T4Tddd &A) {
//    return {{{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A)), std::get<0>(std::get<2>(A)), std::get<0>(std::get<3>(A))},
//             {std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A)), std::get<1>(std::get<2>(A)), std::get<1>(std::get<3>(A))},
//             {std::get<2>(std::get<0>(A)), std::get<2>(std::get<1>(A)), std::get<2>(std::get<2>(A)), std::get<2>(std::get<3>(A))}}};
// };

// std::tuple<V_d, V_d> Transpose(const std::vector<Tdd> &V_Tdd) {
//    V_d v0(V_Tdd.size()), v1(V_Tdd.size());
//    for (auto i = 0; i < V_Tdd.size(); ++i) {
//       v0[i] = std::get<0>(V_Tdd[i]);
//       v1[i] = std::get<1>(V_Tdd[i]);
//    }
//    return {v0, v1};
// };
// std::tuple<V_d, V_d, V_d> Transpose(const std::vector<Tddd> &V_Tdd) {
//    V_d v0(V_Tdd.size()), v1(V_Tdd.size()), v2(V_Tdd.size());
//    for (auto i = 0; i < V_Tdd.size(); ++i) {
//       v0[i] = std::get<0>(V_Tdd[i]);
//       v1[i] = std::get<1>(V_Tdd[i]);
//       v2[i] = std::get<2>(V_Tdd[i]);
//    }
//    return {v0, v1, v2};
// };

/* -------------------------------------------------------------------------- */
// bool isfinite(const V_d &v_IN)
// {
//   for (const auto &v : v_IN)
//     if (!std::isfinite(v))
//       return false;
//   return true;
// };

bool isfinite(const V_d &v_IN) {
   for (const auto &v : v_IN)
      if (v < -1E+20 || v > 1E+20)
         return false;
   return true;
};

bool myIsfinite(const double v) {
   if (v < -1E+20 || v > 1E+20)
      return false;
   else
      return true;
};

bool isFinite(const double v, const double eps = 1E+20) {
   if (v < -eps || v > eps || v != v || std::isnan(v))
      return false;
   else
      return true;
};

bool isFinite(const V_d &v_IN) {
   for (const auto &v : v_IN)
      if (!isFinite(v))
         return false;
   return true;
};

bool isFinite(const VV_d &vv_IN) {
   for (const auto &v : vv_IN)
      if (!isFinite(v))
         return false;
   return true;
};

bool isFinite(const Tdd &v) {
   if (isFinite(std::get<0>(v)) && isFinite(std::get<1>(v)))
      return true;
   else
      return false;
};
bool isFinite(const T4d &v) {
   if (isFinite(std::get<0>(v)) && isFinite(std::get<1>(v)) && isFinite(std::get<2>(v)) && isFinite(std::get<3>(v)))
      return true;
   else
      return false;
};
bool isFinite(const Tddd &v, const double eps = 1E+20) {
   if (isFinite(std::get<0>(v), eps) && isFinite(std::get<1>(v), eps) && isFinite(std::get<2>(v), eps))
      return true;
   else
      return false;
};
bool isFinite(const T6d &v, const double eps = 1E+20) {
   if (isFinite(std::get<0>(v), eps) && isFinite(std::get<1>(v), eps) && isFinite(std::get<2>(v), eps) && isFinite(std::get<3>(v), eps) && isFinite(std::get<4>(v), eps) && isFinite(std::get<5>(v), eps))
      return true;
   else
      return false;
};

bool isFinite(const T3Tddd &V) {
   auto [X, Y, Z] = V;
   if (isFinite(X) && isFinite(Y) && isFinite(Z))
      return true;
   else
      return false;
};
/* ------------------------------------------------------ */
// template <typename T, std::size_t i>
// typename std::enable_if<std::tuple_size<T>::value == i, double>::type
// Total(const T &t) {
//    return 0.0;
// }

// template <typename T, std::size_t i>
// typename std::enable_if<std::is_same<typename std::tuple_element<i, T>::type, double>::value, double>::type
// Total(const T &t) {
//    return std::get<i>(t) + Total<T, i + 1>(t);
// }

// template <typename T, std::size_t i>
// typename std::enable_if<std::is_same<typename std::tuple_element<i, T>::type, int>::value, int>::type
// Total(const T &t) {
//    return std::get<i>(t) + Total<T, i + 1>(t);
// }

double Total(const T4d &v) { return std::get<0>(v) + std::get<1>(v) + std::get<2>(v) + std::get<3>(v); };
double Total(const Tddd &v) { return std::get<0>(v) + std::get<1>(v) + std::get<2>(v); };
double Total(const Tdd &v) { return std::get<0>(v) + std::get<1>(v); };

Tdd Total(const T3Tdd &v) {
   return {std::get<0>(std::get<0>(v)) + std::get<0>(std::get<1>(v)) + std::get<0>(std::get<2>(v)),
           std::get<1>(std::get<0>(v)) + std::get<1>(std::get<1>(v)) + std::get<1>(std::get<2>(v))};
};
Tddd Total(const T2Tddd &X) {
   return {(std::get<0>(std::get<0>(X)) + std::get<0>(std::get<1>(X))),
           (std::get<1>(std::get<0>(X)) + std::get<1>(std::get<1>(X))),
           (std::get<2>(std::get<0>(X)) + std::get<2>(std::get<1>(X)))};
};
Tddd Total(const T3Tddd &X) {
   return {(std::get<0>(std::get<0>(X)) + std::get<0>(std::get<1>(X)) + std::get<0>(std::get<2>(X))),
           (std::get<1>(std::get<0>(X)) + std::get<1>(std::get<1>(X)) + std::get<1>(std::get<2>(X))),
           (std::get<2>(std::get<0>(X)) + std::get<2>(std::get<1>(X)) + std::get<2>(std::get<2>(X)))};
};
Tddd Total(const T4Tddd &X) {
   return {(std::get<0>(std::get<0>(X)) + std::get<0>(std::get<1>(X)) + std::get<0>(std::get<2>(X)) + std::get<0>(std::get<3>(X))),
           (std::get<1>(std::get<0>(X)) + std::get<1>(std::get<1>(X)) + std::get<1>(std::get<2>(X)) + std::get<1>(std::get<3>(X))),
           (std::get<2>(std::get<0>(X)) + std::get<2>(std::get<1>(X)) + std::get<2>(std::get<2>(X)) + std::get<2>(std::get<3>(X)))};
};

Tddd Total(const T8Tddd &X) {
   return {
       (std::get<0>(std::get<0>(X)) + std::get<0>(std::get<1>(X)) + std::get<0>(std::get<2>(X)) + std::get<0>(std::get<3>(X)) + std::get<0>(std::get<4>(X)) + std::get<0>(std::get<5>(X)) + std::get<0>(std::get<6>(X)) + std::get<0>(std::get<7>(X))),
       (std::get<1>(std::get<0>(X)) + std::get<1>(std::get<1>(X)) + std::get<1>(std::get<2>(X)) + std::get<1>(std::get<3>(X)) + std::get<1>(std::get<4>(X)) + std::get<1>(std::get<5>(X)) + std::get<1>(std::get<6>(X)) + std::get<1>(std::get<7>(X))),
       (std::get<2>(std::get<0>(X)) + std::get<2>(std::get<1>(X)) + std::get<2>(std::get<2>(X)) + std::get<2>(std::get<3>(X)) + std::get<2>(std::get<4>(X)) + std::get<2>(std::get<5>(X)) + std::get<2>(std::get<6>(X)) + std::get<2>(std::get<7>(X)))};
};

Tddd Total(const std::vector<Tddd> &V) {
   Tddd ret = {0, 0, 0};
   for (const auto &v : V)
      ret += v;
   return ret;
};
T4d Total(const std::vector<T4d> &V) {
   T4d ret = {0, 0, 0, 0};
   for (const auto &v : V)
      ret += v;
   return ret;
};
T6d Total(const std::vector<T6d> &V) {
   T6d ret = {0, 0, 0, 0, 0, 0};
   for (const auto &v : V)
      ret += v;
   return ret;
};
double Total(const std::vector<double> &V) {
   double ret = 0;
   for (const auto &v : V)
      ret += v;
   return ret;
};
/* -------------------------------------------------------------------------- */

double Mean(const Tdd &v) { return Total(v) / 2.; };
double Mean(const T4d &v) { return Total(v) / 4.; };
double Mean(const Tddd &v) { return Total(v) / 3.; };

Tddd Mean(const Tddd &A, const Tddd &B, const Tddd &C) { return (A + B + C) / 3.; };

Tddd Mean(const T2Tddd &X) { return Total(X) / 2.; };
Tddd Mean(const T3Tddd &X) { return Total(X) / 3.; };
Tddd Mean(const T4Tddd &X) { return Total(X) / 4.; };
Tddd Mean(const T8Tddd &X) { return Total(X) / 8.; };

Tddd Mean(const std::vector<Tddd> &X) {
   Tddd ret = {0., 0., 0.};
   for (const auto &x : X) {
      std::get<0>(ret) += std::get<0>(x);
      std::get<1>(ret) += std::get<1>(x);
      std::get<2>(ret) += std::get<2>(x);
   }
   return ret / ((double)X.size());
};

double Mean(const std::vector<double> &X) { return Total(X) / X.size(); };

// Tddd Mean(const std::unordered_set<Tddd> &X)
// {
// 	Tddd ret = {0., 0., 0.};
// 	for (const auto &x : X)
// 	{
// 		std::get<0>(ret) += std::get<0>(x);
// 		std::get<1>(ret) += std::get<1>(x);
// 		std::get<2>(ret) += std::get<2>(x);
// 	}
// 	return ret / ((double)X.size());
// };

template <class T>
T Min(const std::vector<T> &v) {
   return *std::min_element(std::begin(v), std::end(v));
};
template <class T>
T Min(const std::vector<std::vector<T>> &v) {
   T tmp, ret(0.);
   for (size_t i = 0; i < v.size(); i++) {
      tmp = Min(v[i]);
      if (ret > tmp)
         ret = tmp;
   }
   return ret;
};

Tdd MinMax(const std::vector<double> &v) {
   auto it = std::minmax_element(std::begin(v), std::end(v));
   return {*it.first, *it.second};
};

template <class T>
T Max(const std::vector<T> &v) { return *std::max_element(std::begin(v), std::end(v)); };
template <class T>
T Max(const std::vector<std::vector<T>> &v) {
   T tmp, ret(0.);
   for (size_t i = 0; i < v.size(); i++) {
      tmp = Max(v[i]);
      if (ret < tmp)
         ret = tmp;
   }
   return ret;
};

double Min(const Tdd &A) { return ((std::get<0>(A) < std::get<1>(A)) ? std::get<0>(A) : std::get<1>(A)); };
double Max(const Tdd &A) { return ((std::get<0>(A) > std::get<1>(A)) ? std::get<0>(A) : std::get<1>(A)); };
// Tdd MinMax(const Tdd &A) { return ((std::get<0>(A) < std::get<1>(A)) ? Tdd{std::get<0>(A), std::get<1>(A)} : Tdd{std::get<1>(A), std::get<0>(A)}); };
// T3Tdd MinMaxTranspose(const T2Tddd &A) {
//    return {MinMax(Tdd{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A))}),
//            MinMax(Tdd{std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A))}),
//            MinMax(Tdd{std::get<2>(std::get<0>(A)), std::get<2>(std::get<1>(A))})};
// };
double Min(const Tddd &A) {
   // イコールが重要
   const auto [X, Y, Z] = A;
   return X <= Y ? (X <= Z ? X : Z) : (Y <= Z ? Y : Z);
};
double Max(const Tddd &A) {
   auto [X, Y, Z] = A;
   return X >= Y ? (X >= Z ? X : Z) : (Y >= Z ? Y : Z);
};
double Min(const T4d &A) {
   const auto [X, Y, Z, W] = A;
   if (X <= Y && X <= Z && X <= W)
      return X;
   else if (Y <= Z && Y <= W)
      return Y;
   else if (Z <= W)
      return Z;
   else
      return W;
};
double Max(const T4d &A) {
   auto [X, Y, Z, W] = A;
   return X >= Y ? (X >= Z ? (X >= W ? X : W) : (Z >= W ? Z : W)) : (Y >= Z ? (Y >= W ? Y : W) : (Z >= W ? Z : W));
};
Tdd MinMax(const T4d &A) {
   auto [X, Y, Z, W] = A;
   return {X >= Y ? (X >= Z ? (X >= W ? X : W) : (Z >= W ? Z : W)) : (Y >= Z ? (Y >= W ? Y : W) : (Z >= W ? Z : W)),
           X <= Y ? (X <= Z ? (X <= W ? X : W) : (Z <= W ? Z : W)) : (Y <= Z ? (Y <= W ? Y : W) : (Z <= W ? Z : W))};
};
double FiniteMin(const Tddd &A) {
   // イコールが重要
   auto [X, Y, Z] = A;
   if (isFinite(X) && !isFinite(Y) && !isFinite(Z)) {
      Y = X;
      Z = X;
   } else if (!isFinite(X) && !isFinite(Y) && isFinite(Z)) {
      Y = Z;
      X = Z;
   } else if (!isFinite(X) && isFinite(Y) && !isFinite(Z)) {
      Z = Y;
      X = Y;
   } else if (isFinite(X) && isFinite(Y) && !isFinite(Z))
      Z = X;
   else if (isFinite(X) && !isFinite(Y) && isFinite(Z))
      Y = Z;
   else if (!isFinite(X) && isFinite(Y) && isFinite(Z))
      X = Z;
   //
   return Min(Tddd{X, Y, Z});
};

// T3Tdd MinMaxTranspose(const T3Tddd &A) {
//    const auto [Xs, Ys, Zs] = Transpose(A);
//    return {{{Min(Xs), Max(Xs)}, {Min(Ys), Max(Ys)}, {Min(Zs), Max(Zs)}}};
// };
// T3Tdd MinMaxTranspose(const T4Tddd &A) {
//    const auto [Xs, Ys, Zs] = Transpose(A);
//    return {{{Min(Xs), Max(Xs)}, {Min(Ys), Max(Ys)}, {Min(Zs), Max(Zs)}}};
// };

double Max(const std::vector<Tddd> &A) {
   double ret = -1E-10;
   for (const auto &a : A)
      if (ret < Max(a))
         ret = Max(a);
   return ret;
};

T3Tdd MinMax(const std::tuple<V_d, V_d, V_d> &A) {
   return {{MinMax(std::get<0>(A)),
            MinMax(std::get<1>(A)),
            MinMax(std::get<2>(A))}};
};
// T3Tdd MinMaxTranspose(const std::vector<Tddd> &A) { return MinMax(Transpose(A)); };
/* ------------------------------------------------------ */
double Rot(const V_d vec1, const V_d vec2) {
   return vec1[0] * vec2[1] - vec1[1] * vec2[0];
};
std::vector<V_d> Inv(const std::vector<V_d> &mat) {
   std::vector<V_d> ans(mat.size(), V_d(mat[0].size(), 0.));
   double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
   ans[1][1] = mat[0][0] / det;
   ans[0][1] = -mat[0][1] / det;
   ans[1][0] = -mat[1][0] / det;
   ans[0][0] = mat[1][1] / det;
   return ans;
};
//==========================================================
// vector operators
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> Cross(const std::vector<T> &A) {
   return {-A[1], A[0]};
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> Cross(const std::vector<T> &A, const std::vector<T> &X) {
   if (A.size() == 3)
      return {A[1] * X[2] - A[2] * X[1],
              A[2] * X[0] - A[0] * X[2],
              A[0] * X[1] - A[1] * X[0]};
   else if (A.size() == 2)
      return Cross(std::vector<T>{A[0], A[1], 0.}, std::vector<T>{X[0], X[1], 0.});
   else {
      std::stringstream ss;
      ss << "Invalid vectors are passed\n";
      ss << "A " << A << " X " << X;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   }
};
//==========================================================
V_d log10(const V_d &vec) {
   V_d ret(vec.size());
   for (size_t i = 0; i < vec.size(); i++)
      ret[i] = std::log10(vec[i]);
   return ret;
};
V_d log(const V_d &vec) {
   V_d ret(vec.size());
   for (size_t i = 0; i < vec.size(); i++)
      ret[i] = std::log(vec[i]);
   return ret;
};
/* ------------------------------------------------------ */
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Norm(const T &x) { return std::abs(x); };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Norm(const std::vector<T> &vec) {
   // return std::sqrt(std::inner_product(vec.cbegin(), vec.cend(), vec.cbegin(), 0.));
   T ret = 0;
   for (const auto &v : vec)
      ret += v * v;
   return std::sqrt(ret);
};
double Norm(const std::vector<Tddd> &vec) {
   double ret = 0;
   for (const auto &v : vec)
      ret += Dot(v, v);
   return std::sqrt(ret);
};
/* ------------------------------------------------------ */
V_d Normalize(const V_d &X) { return X / Norm(X); };
VV_d Orthogonalize(VV_d VV) {
   // VVは正方行列に限る
   for (auto i = 0; i < VV.size(); ++i) {
      for (auto j = 0; j < i; ++j)
         VV[i] -= Dot(VV[i], VV[j]) * VV[j];  // VV[j]は直行化が終わったベクトル
      VV[i] /= Norm(VV[i]);
   }
   return VV;
};
double Norm3d(const V_d &vec) {
   if (vec.size() != 3) {
      std::stringstream ss;
      ss << vec;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   }
   return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Abs(const V_d &vec) {
   T tmp(0);
   for (size_t i = 0; i < vec.size(); i++)
      tmp += vec[i] * vec[i];
   return sqrt(tmp);
};
/* ------------------------------------------------------ */
// template <typename T>
// std::vector<T> RotateLeft(const std::vector<T> &vecs, const int n = 1) {  // 2020/03/22
//    std::vector<T> ret(vecs);
//    std::rotate(ret.begin(), ret.begin() + n, ret.end());
//    return ret;
// }
// template <typename T>
// std::tuple<T, T, T> Reverse(const std::tuple<T, T, T> &vecs) {  // 2022年3月21日
//    return {std::get<2>(vecs), std::get<1>(vecs), std::get<0>(vecs)};
// };
// template <typename T>
// std::tuple<T, T, T> RotateLeft(const std::tuple<T, T, T> &vecs, int n = 1) {  // 2021/12/05
//    n = n % 3;
//    if (n == 0)
//       return vecs;
//    else if (n == 1)
//       return {std::get<1>(vecs), std::get<2>(vecs), std::get<0>(vecs)};
//    else
//       return {std::get<2>(vecs), std::get<0>(vecs), std::get<1>(vecs)};
// }

// template <typename T, size_t N>
// std::array<T, N> RotateLeft(const std::array<T, N> &arr, int n = 1) {
//    std::array<T, N> result;
//    n = n % N;

//    for (size_t i = 0; i < N; ++i) {
//       result[(i + N - n) % N] = arr[i];
//    }

//    return result;
// }

// T3Tddd RotateLeft(const T3Tddd &vecs, int n = 1) {  // 2021/12/05
//    n = n % 3;
//    if (n == 0)
//       return vecs;
//    else if (n == 1)
//       return {{std::get<1>(vecs), std::get<2>(vecs), std::get<0>(vecs)}};
//    else
//       return {{std::get<2>(vecs), std::get<0>(vecs), std::get<1>(vecs)}};
// }
// T3Tddd RotateRight(const T3Tddd &vecs, int n = 1) {  // 2021/12/05
//    n = n % 3;
//    if (n == 0)
//       return vecs;
//    else if (n == 1)
//       return {{std::get<2>(vecs), std::get<0>(vecs), std::get<1>(vecs)}};
//    else
//       return {{std::get<1>(vecs), std::get<2>(vecs), std::get<0>(vecs)}};
// }
// Tddd RotateRight(const Tddd &vecs, int n = 1) {  // 2021/12/05
//    n = n % 3;
//    if (n == 0)
//       return vecs;
//    else if (n == 1)
//       return {std::get<2>(vecs), std::get<0>(vecs), std::get<1>(vecs)};
//    else
//       return {std::get<1>(vecs), std::get<2>(vecs), std::get<0>(vecs)};
// }
// template <class Type>
// std::vector<Type> RotateRight(const std::vector<Type> &vecs, const int n) {  // 2021/04/07
//    std::vector<Type> ret(vecs);
//    std::rotate(ret.rbegin(), ret.rbegin() + n, ret.rend());
//    return ret;
// }
/* ------------------------------------------------------ */
// 2021/06/14
struct Quaternion {
   Tddd v;
   double a, b, c, d;
   T4d q;

   // cos(q/2) + (ux*i + uy*j+ uz*k) * sin(q/2)
   Quaternion() : a(1), b(0), c(0), d(0), v({0, 0, 0}), q({1, 0, 0, 0}){};
   Quaternion(const double aIN, const double bIN, const double cIN, const double dIN) : a(aIN),
                                                                                        b(bIN),
                                                                                        c(cIN),
                                                                                        d(dIN),
                                                                                        v({bIN, cIN, dIN}),
                                                                                        q({aIN, bIN, cIN, dIN}){};
   Quaternion(const T4d &qIN) : a(std::get<0>(qIN)),
                                b(std::get<1>(qIN)),
                                c(std::get<2>(qIN)),
                                d(std::get<3>(qIN)),
                                v({std::get<1>(qIN), std::get<2>(qIN), std::get<3>(qIN)}),
                                q(qIN){};
   Quaternion(const Tddd &axis, const double angle) : v(Normalize(axis) * sin(angle / 2.)),
                                                      a(cos(angle / 2.)),
                                                      b(std::get<0>(v)),
                                                      c(std::get<1>(v)),
                                                      d(std::get<2>(v)),
                                                      q({a, b, c, d}){
                                                          // 空間回転　q = cos(theta/2) + n*sin(theta/2)
                                                      };
   /* ------------------------------------------------------ */
   T3Tddd Rv() const {
      //%固定した座標系(global座標)における．位置ベクトルの回転をするために使う．
      // ノーマライズされていなくていい
      double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
      return {{{a2 + b2 - c2 - d2, 2. * b * c - 2. * a * d, 2. * a * c + 2. * b * d},
               {2. * b * c + 2. * a * d, a2 - b2 + c2 - d2, -2. * a * b + 2. * c * d},
               {-2. * a * c + 2. * b * d, 2. * a * b + 2. * c * d, a2 - b2 - c2 + d2}}};
      // OK
   }
   T3Tddd Rs() const {
      //%物体座標系を回転することで，global座標が物体座標にとってどのように移動するかを計算するために使う．
      // ノーマライズされていなくていい
      double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
      return {{{a2 + b2 - c2 - d2, 2. * b * c + 2. * a * d, -2. * a * c + 2. * b * d},
               {2. * b * c - 2. * a * d, a2 - b2 + c2 - d2, 2. * a * b + 2. * c * d},
               {2. * a * c + 2. * b * d, -2. * a * b + 2. * c * d, a2 - b2 - c2 + d2}}};
      // OK
   }
   Tddd Rv(const Tddd &uIN) const { return Dot(this->Rv(), uIN); }
   Tddd Rs(const Tddd &uIN) const { return Dot(this->Rs(), uIN); }

   // ノーマライズされていなくていい
   T3Tddd R() const { return this->Rv(); }
   Tddd R(const Tddd &uIN) const { return Dot(this->Rv(), uIN); }
   /* ------------------------------------------------------ */
   // ドローン工学入門(1.74) or https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
   //  double yaw() const { return std::atan2(2. * (a * d + b * c), 1. - 2. * (c * c + d * d)); }
   //  double pitch() const { return std::asin(2. * (a * c - b * d)); }
   //  double roll() const { return std::atan2(2. * (a * b + c * d), 1. - 2. * (a * a + d * d)); }

   // double yaw() const { return std::atan2(2.0 * (c * d + a * b), a * a - b * b - c * c + d * d); };
   // double pitch() const { return std::asin(-2.0 * (b * d - a * c)); };
   // double roll() const { return std::atan2(2.0 * (b * c + a * d), a * a + b * b - c * c - d * d); };

   double roll() const { return std::atan2(2.0 * (d * c + a * b), 1.0 - 2.0 * (b * b + c * c)); }
   double pitch() const { return std::asin(2.0 * (c * a - d * b)); }
   double yaw() const { return std::atan2(2.0 * (d * a + b * c), -1.0 + 2.0 * (a * a + b * b)); }

   // yaw，pitch，rollは，回転の順序でもある
   Tddd YPR() const { return {this->yaw(), this->pitch(), this->roll()}; }

   T3Tddd Ryaw() const {
      double t = this->yaw();
      return {{{cos(t), sin(t), 0.}, {-sin(t), cos(t), 0.}, {0., 0., 1.}}};
   }

   T3Tddd Rpitch() const {
      double t = this->pitch();
      return {{{cos(t), 0., -sin(t)}, {0., 1., 0.}, {sin(t), 0., cos(t)}}};
   }

   T3Tddd Rroll() const {
      double t = this->roll();
      return {{{1., 0., 0.}, {0., cos(t), sin(t)}, {0., -sin(t), cos(t)}}};
   }

   Tddd Ryaw(const Tddd &uIN) const { return Dot(this->Ryaw(), uIN); }
   Tddd Rpitch(const Tddd &uIN) const { return Dot(this->Rpitch(), uIN); }
   Tddd Rroll(const Tddd &uIN) const { return Dot(this->Rroll(), uIN); }
   /* ------------------------------------------------------ */

   const T4d &operator()() const { return q; }

   Quaternion conjugate() const {
      return Quaternion(T4d{a, -b, -c, -d});
   };

   Quaternion approxNextQuaternion(const Tddd &w, const double dt) const {
      return Quaternion(this->q + this->d_dt(w * dt)());
   };

   Quaternion d_dt(const Tddd &w /*angular velocity*/) const {
      // 回転変化率w=dtheta/dtは，クォータニオン変化はQ.d_dt(w)と同じである．
      /*
      W = {{0, -w0, -w1, -w2},
              {w0, 0, w2, -w1},
              {w1, -w2, 0, w0},
              {w2, w1, -w0, 0}}

      inverse[W] = {{0, w0, w1, w2},
                              {-w0, 0, -w2, w1},
                              {-w1, w2, 0, -w0},
                              {-w2, -w1, w0,0}}/(w0^2 + w1^2 + w2^2)
                              = -W/(w0^2 + w1^2 + w2^2)
      */

      return Quaternion((-this->b * std::get<0>(w) - this->c * std::get<1>(w) - this->d * std::get<2>(w)) * 0.5, /*{0,-w0,-w1,-w2}.{a,b,c,d}/2*/
                        (this->a * std::get<0>(w) - this->d * std::get<1>(w) + this->c * std::get<2>(w)) * 0.5,  /*{w0,0,w2,-w1}.{a,b,c,d}/2*/
                        (this->d * std::get<0>(w) + this->a * std::get<1>(w) - this->b * std::get<2>(w)) * 0.5,  /*{w1,-w2,0,w0}.{a,b,c,d}/2*/
                        (-this->c * std::get<0>(w) + this->b * std::get<1>(w) + this->a * std::get<2>(w)) * 0.5 /*{w2,w1,-w0,0}.{a,b,c,d}/2*/);
   };

   void set(const T4d &qIN) {
      double n = Norm(qIN);
      this->a = std::get<0>(qIN) / n;
      std::get<0>(this->v) = this->b = std::get<1>(qIN) / n;
      std::get<1>(this->v) = this->c = std::get<2>(qIN) / n;
      std::get<2>(this->v) = this->d = std::get<3>(qIN) / n;
      this->q = {this->a, this->b, this->c, this->d};
   };

   void set(const Quaternion &qIN) {
      this->set(qIN());
   };
};

Quaternion operator*(const Quaternion &A, const Quaternion &B) {
   // https://en.wikipedia.org/wiki/Quaternion
   // Tddd v = A.a * B.v + B.a * A.v + Cross(A.v, B.v);
   // return Quaternion(T4d{A.a * B.a - Dot(A.v, B.v),
   // 					  std::get<0>(v),
   // 					  std::get<1>(v),
   // 					  std::get<2>(v)}); // ok
   auto [a1, b1, c1, d1] = A.q;
   auto [a2, b2, c2, d2] = B.q;
   return Quaternion(T4d{a1 * a2 + -b1 * b2 + -c1 * c2 + -d1 * d2,
                         a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                         a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                         a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2});
};
Quaternion &operator*=(Quaternion &A, const Quaternion &B) {
   return A = (A * B);  // ok
};
Quaternion operator*(Quaternion A, const double dt) {
   A.set(A.q * dt);  // ok
   return A;
};
/* ------------------------- 足し算 ------------------------ */
Quaternion operator+(Quaternion B, const double A) {
   //()はただのT4d
   B.set(A + B.q);
   return B;
};
Quaternion operator+(const double A, Quaternion B) {
   //()はただのT4d
   B.set(A + B.q);
   return B;
};
Quaternion operator+(Quaternion A, const Quaternion &B) {
   //()はただのT4d
   A.set(A.q + B.q);
   return A;
};
Quaternion operator+(Quaternion A, const T4d &B) {
   //()はただのT4d
   A.set(A.q + B);
   return A;
};
Quaternion operator+(const T4d &B, Quaternion A) {
   //()はただのT4d
   A.set(A.q + B);
   return A;
};
Quaternion &operator+=(Quaternion &A, const Quaternion &B) {
   A.set(A.q + B.q);
   return A;
};
double Norm(const Quaternion &A) {
   return Norm(A.q);
};
/* ------------------------------------------------------ */
// double VectorAngle(const Tddd &V1, const Tddd &V2) { return std::atan2(Norm(Cross(V1, V2)), Dot(V1, V2)); };
double VectorAngle(const Tddd &V1, const Tddd &V2) { return std::acos(Dot(V1, V2) / (Norm(V1) * Norm(V2))); };
double VectorAngle(const Tddd &X1, const Tddd &X2, const Tddd &axis) {
   // r = RotationMatrix[VectorAngle[A, {1, 0, 0}] * If[Dot[cross, axis] >= 0, 1, -1], axis];
   return VectorAngle(X1, X2) * (Dot(Cross(X1, X2), axis) >= 0. ? 1. : -1.);
   // return std::atan2(Norm(Cross(X1 - axis, X2 - axis)), Dot(X1 - axis, X2 - axis));
};

double MyVectorAngle(const Tddd &v0, const Tddd &v1, const Tddd &frontdir /*右手系のzとなる*/) {
   //      /|
   //     / |
   //    /  |y
   //   /   |
   //  /q___|
   //     x
   // q = atan2(double y, double x), q = [-pi, pi]
   //
   //
   // this can distingish ccw(positive) or cw(negative)

   // auto Y = Cross(frontdir, v0); //右手系
   // return atan2(Dot(v1, Y / Norm(Y)), Dot(v1, v0 / Norm(v0)));

   // V_d z = Cross(v0/*x*/, v1);
   // V_d z = frontdir;
   auto z = frontdir / Norm(frontdir);
   auto y = Cross(z, v0 /*x*/);
   y = y / Norm(y);
   auto x = v0 / Norm(v0);
   // {x,y,z}右手系

   //      /|
   //     / |
   //  v0/  |Dot(v1,y)
   //   /   |
   //  /q___|
   //     Dot(v1,v0/Norm(v0)

   return std::atan2(Dot(v1, y), Dot(v1, x));
};
double MyVectorAngle(const V_d &v0, const V_d &v1, const V_d &frontdir /*右手系のzとなる*/) {
   //      /|
   //     / |
   //    /  |y
   //   /   |
   //  /q___|
   //     x
   // q = atan2(double y, double x), q = [-pi, pi]
   //
   //
   // this can distingish ccw(positive) or cw(negative)

   // auto Y = Cross(frontdir, v0); //右手系
   // return atan2(Dot(v1, Y / Norm(Y)), Dot(v1, v0 / Norm(v0)));

   // V_d z = Cross(v0/*x*/, v1);
   // V_d z = frontdir;
   V_d z = frontdir / Norm(frontdir);
   V_d y = Cross(z, v0 /*x*/);
   y = y / Norm(y);
   V_d x = v0 / Norm(v0);
   // {x,y,z}右手系

   //      /|
   //     / |
   //  v0/  |Dot(v1,y)
   //   /   |
   //  /q___|
   //     Dot(v1,v0/Norm(v0)

   return std::atan2(Dot(v1, y), Dot(v1, x));
};

double MyVectorAngle(const V_d &v0, const V_d &v1) {
   // cannot distingish ccw(positive) or cw(negative)
   // return MyVectorAngle(v0, v1, Cross(v0, v1));

   return std::acos(Dot(v0, v1) / (Norm(v0) * Norm(v1)));
   // これはMathematicaの定義と同じ:
   //  a = {a0, a1, a2}
   //  b = {b0, b1, b2}
   //  Refine[VectorAngle[a, b],
   //   Assumptions -> {a0 \[Element] Reals && a1 \[Element] Reals &&
   //      a2 \[Element] Reals && b0 \[Element] Reals &&
   //      b1 \[Element] Reals && b2 \[Element] Reals}]
};
double MyVectorAngle(const Tddd &v0, const Tddd &v1) {
   return std::acos(Dot(v0, v1) / (Norm(v0) * Norm(v1)));
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T MathematicaVectorAngle(const std::vector<T> &V1, const std::vector<T> &V2) {
   if (V1.size() > 1) {
      return std::atan2(Norm(Cross(V1, V2)), Dot(V1, V2));
   } else {
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, " size is invalid"));
      return 0.;
   }
};

// point base
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T VectorAngle(const std::vector<T> &X1, const std::vector<T> &X2, const std::vector<T> &X0) {
   if (X1.size() > 1) {
      return std::atan2(Norm(Cross(X1 - X0, X2 - X0)), Dot(X1 - X0, X2 - X0));
   } else {
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, " size is invalid"));
      return 0.;
   }
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T VectorAngle(const std::vector<T> &X1, const std::vector<T> &X2) {
   return VectorAngle(X1, X2, {0., 0.});
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T VectorAngleDirected(const std::vector<T> &X1, const std::vector<T> &X2) {
   T a = VectorAngle(X1, X2, {0., 0.});
   return (a < 0) ? (2 * M_PI - a) : a;
};
// template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// T DirectedArea(const std::vector<T> &vec1,
//                const std::vector<T> &vec2,
//                const std::vector<T> &n /*the direction*/) {
//    //  vec1
//    // --->*
//    //     | vec2
//    //     v
//    return Dot(Cross(vec1, vec2), n);
// };

// template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// T DirectedArea(const std::vector<T> &a,
//                const std::vector<T> &b,
//                const std::vector<T> &c,
//                const std::vector<T> &n) {
//    // a --> b
//    //       |
//    //       v
//    //	   c
//    return Dot(Cross(b - a, c - b), n);
// };

double TriangleArea(const Tddd &a, const Tddd &b, const Tddd &c) {
   double A = Norm(a - c);
   double B = Norm(b - a);
   double C = Norm(c - b);
   // this->area = std::abs(Norm(c) / Dot(a, b));
   double s = 0.5 * (A + B + C);
   return std::sqrt(s * (s - A) * (s - B) * (s - C));
};

double TriangleArea(const T3Tddd &abc) {
   auto [a, b, c] = abc;
   double A = Norm(a - c);
   double B = Norm(b - a);
   double C = Norm(c - b);
   // this->area = std::abs(Norm(c) / Dot(a, b));
   double s = 0.5 * (A + B + C);
   return std::sqrt(s * (s - A) * (s - B) * (s - C));
};
static const std::tuple<Tiii, Tiii, Tiii, Tiii> TetrahedronPolygons = {{1, 2, 3}, {0, 1, 3}, {0, 2, 1}, {0, 3, 2}};
T4Tddd TetrahedronNormals(const T4Tddd &X) {
   auto [x0, x1, x2, x3] = X;
   Tddd c = (x0 + x1 + x2 + x3) * 0.25;
   Tddd n0 = Normalize(Cross(x1 - x0, x3 - x0)), m0 = (x0 + x1 + x3) / 3. - c;
   Tddd n1 = Normalize(Cross(x3 - x0, x2 - x0)), m1 = (x0 + x2 + x3) / 3. - c;
   Tddd n2 = Normalize(Cross(x2 - x0, x1 - x0)), m2 = (x0 + x1 + x2) / 3. - c;
   Tddd n3 = Normalize(Cross(x2 - x1, x3 - x1)), m3 = (x1 + x2 + x3) / 3. - c;
   return {(Dot(n0, m0) >= 0 ? 1. : -1.) * n0,
           (Dot(n1, m1) >= 0 ? 1. : -1.) * n1,
           (Dot(n2, m2) >= 0 ? 1. : -1.) * n2,
           (Dot(n3, m3) >= 0 ? 1. : -1.) * n3};
};

double TetrahedronVolume(const T4Tddd &X) {
   return 0.16666666666666666667 * std::abs(Det(T3Tddd{std::get<1>(X) - std::get<0>(X),
                                                       std::get<2>(X) - std::get<0>(X),
                                                       std::get<3>(X) - std::get<0>(X)}));
};

Tddd TetrahedronCircumCenter(const Tddd &a, const Tddd &b, const Tddd &c, const Tddd &d) {
   double a2 = Dot(a, a);
   return Dot(Inverse(T3Tddd{b - a, c - a, d - a}),
              0.5 * Tddd{Dot(b, b) - a2, Dot(c, c) - a2, Dot(d, d) - a2});
};

Tddd TetrahedronCircumCenter(const T4Tddd &abcd) {
   return TetrahedronCircumCenter(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd), std::get<3>(abcd));
};

Tddd TriangleAngles(const T3Tddd &abc) {
   auto [a, b, c] = abc;
   return {VectorAngle(b - a, c - a),
           VectorAngle(c - b, a - b),
           VectorAngle(a - c, b - c)};
};

Tddd TriangleNormal(const Tddd &a, const Tddd &b, const Tddd &c) {
   return Normalize(Cross((b - a), (c - a)));
};
Tddd TriangleNormal(const T3Tddd &abc) {
   auto [a, b, c] = abc;
   return Normalize(Cross((b - a), (c - a)));
};

bool isFlat(const Tddd &a, const Tddd &b, const double lim_rad) {
   return Dot(Normalize(a), Normalize(b)) > cos(lim_rad);  // if this satisfied which means falt
};

// // 多角形の頂点における隣り合う辺が作る面積を計算する．法線と逆の場合はマイナス．範囲は[pi,-pi]
// V_d DirectedArea(const VV_d &abcd,
//                  const V_d &n) {
//    V_d ret(abcd.size());
//    ret[0] = DirectedArea(abcd[abcd.size() - 1], abcd[0], abcd[1], n);
//    for (size_t i = 0; i < abcd.size() - 2; i++)
//       ret[i + 1] = DirectedArea(abcd[i], abcd[i + 1], abcd[i + 2], n);
//    ret[abcd.size() - 1] = DirectedArea(abcd[abcd.size() - 2], abcd[abcd.size() - 1], abcd[0], n);
//    return ret;
// };

double InteriorAngle(const V_d &x, const V_d &b, const V_d &z) {
   // this can distingish ccw(positive) or cw(negative)
   auto Y = Cross(z, x);  // 右手系
   return std::atan2(Dot(b, Y / Norm(Y)), Dot(b, x / Norm(x)));
};

// V_d RotateVector(const double angle, const V_d &axis, const V_d &v) {
//    V_d k = axis / Norm(axis);
//    return v * std::cos(angle) + Cross(k, v) * std::sin(angle) + k * (Dot(k, v)) * (1. - std::cos(angle));
// };

// V_d RotateVectorPI2(const V_d &axis, const V_d &v) {
//    V_d k = axis / Norm(axis);
//    return Cross(k, v) + k * (Dot(k, v));
// };

// double SphericalInteriorAngle(const V_d &org, const V_d &p0, const V_d &p1, const V_d &p2) {
//    auto v0 = (p0 - org);
//    auto v1 = (p1 - org);
//    auto v2 = (p2 - org);
//    auto s01 = RotateVectorPI2(Cross(v1, v0), v1);
//    auto s02 = RotateVectorPI2(Cross(v1, v2), v1);
//    return MyVectorAngle(s01, s02);
// };

T3Tddd RotationMatrix(const double theta, const Tddd &V) {
   // // Euler Rodrigues
   // double c = cos(theta), s = sin(theta), e = (1. - c);
   // return {{std::get<0>(V) * std::get<0>(V) * e + c, std::get<0>(V) * std::get<1>(V) * e - std::get<2>(V) * s, std::get<0>(V) * std::get<2>(V) * e + std::get<1>(V) * s},
   // 		{std::get<0>(V) * std::get<1>(V) * e + std::get<2>(V) * s, std::get<1>(V) * std::get<1>(V) * e + c, std::get<1>(V) * std::get<2>(V) * e - std::get<0>(V) * s},
   // 		{std::get<0>(V) * std::get<2>(V) * e - std::get<1>(V) * s, std::get<2>(V) * std::get<1>(V) * e + std::get<0>(V) * s, std::get<2>(V) * std::get<2>(V) * e + c}};
   Quaternion q(V, theta);
   return q.Rv();
};
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/*                        vector modification operators                       */
/* -------------------------------------------------------------------------- */
Tddd Projection(const Tddd &v, Tddd n) {
   /* the component in n direction of v will be returned */
   n = Normalize(n);
   return Dot(v, n) * n;
};

Tddd Chop(const Tddd &v, const Tddd &n) {
   /* the component in n direction of v will be chopped */
   return v - Projection(v, n);
};

Tddd Reflect(const Tddd &v, const Tddd &n) {
   /* n is a normal vector of a surface*/
   return v - 2. * Projection(v, n);
};

Tddd Scaled(const Tddd &v, const double d) {
   /* Returns the scaled vector */
   return d * Normalize(v);
};

Tddd Mirror(const Tddd &position, const Tddd &a_point_on_mirror, const Tddd &normal_vector) {
   return position + 2 * Projection(a_point_on_mirror - position, normal_vector);
};
/* ------------------------------------------------------ */
/*                        線形方程式の解法                   */
/* ------------------------------------------------------ */

V_d diagonal_scaling_vector(VV_d mat) {
   V_d ret(mat.size());
   int s = mat.size();
   for (auto i = 0; i < s; ++i)
      ret[i] = 1. / mat[i][i];
   return ret;
};

/* ------------------------------------------------------ */
template <typename T>
bool Between(const T &x, const std::array<T, 2> &minmax) {
   const auto &[min, max] = minmax;
   return (min <= x && x <= max && isFinite(x));
}

/* -------------------------------------------------------------------------- */
/////////////////////////////////////////////////////
template <typename T>
std::vector<T> Range(const T xmin, const T xmax, const T di) {
   std::vector<T> ret;
   T i = 0.;
   while (di * i + xmin <= xmax) {
      ret.push_back(di * i + xmin);
      i = i + 1.;
   }
   return ret;
};
// template <typename T>
// std::vector<T> Subdivide(const T xmin, const T xmax, const int n)
// {
// 	T dx = (xmax - xmin) / n;
// 	std::vector<T> ret(n + 1);
// 	for (int i = 0; i < n + 1; i++)
// 		ret[i] = i * dx + xmin;
// 	return ret;
// };
std::vector<double> Subdivide(const double xmin, const double xmax, const int n) {
   if (n <= 0) {
      std::stringstream ss;
      ss << "Subdivide(" << xmin << "," << xmax << "," << n << ")の位置3に与えられた分割数は正の整数でなければなりません．";
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   }
   std::vector<double> ret(n + 1);
   const double dx = (xmax - xmin) / n;
   for (int i = 0; i < n + 1; i++)
      ret[i] = i * dx + xmin;
   return ret;
};
std::vector<double> Subdivide(const Tdd &xminxmax, const int n) { return Subdivide(std::get<0>(xminxmax), std::get<1>(xminxmax), n); };
template <class T>
std::vector<std::vector<T>> Subdivide(const std::vector<T> &xmin, const std::vector<T> &xmax, const int n) {
   std::vector<T> dx = (xmax - xmin) / n;
   std::vector<std::vector<T>> ret(xmin.size(), std::vector<T>(n + 1, 0.));
   for (size_t i = 0; i < xmin.size(); i++)
      ret[i] = Subdivide(xmin[i], xmax[i], n);
   return Transpose(ret);
};
template <typename T>
std::vector<T> SubdivideByStep(const T xmin, const T xmax, const T di) {
   return Subdivide(xmin, xmax, (int)((xmax - xmin) / di));
};
template <typename T>
std::vector<std::vector<T>> SubdivideByStep(const std::vector<T> &xmin, const std::vector<T> &xmax, const T di) {
   return Subdivide(xmin, xmax, (int)((Norm(xmax - xmin)) / di));
};
template <typename T>
std::vector<std::vector<T>> SubdivideByStep(const std::vector<std::vector<T>> &x, const T di) {
   std::vector<std::vector<T>> ret;
   for (size_t i = 0; i < x.size() - 1; i++) {
      if (!ret.empty())
         ret.pop_back();
      ret = Join(ret, SubdivideByStep(x[i], x[i + 1], di));
   }
   return ret;
};
template <typename T>
std::vector<std::vector<T>> SubdivideByStepExclude(const std::vector<std::vector<T>> &x_in, const T di) {
   std::vector<std::vector<T>> x = x_in;
   std::vector<std::vector<T>> ret;
   T n;
   for (size_t i = 0; i < x.size() - 1; i++) {
      if (x.size() == 1) {
         std::cout << "invalid length" << std::endl;
         std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorOff << std::endl;
         return x;
      } else if (x.size() == 2) {
         return SubdivideByStep(x[i], x[i + 1], di);
      } else if (i == 0) {
         n = Norm(x[i + 1] - x[i]) / di;
         std::vector<std::vector<T>> tmp = SubdivideByStep(x[i], x[i + 1] - (x[i + 1] - x[i]) / (2. * n), di);
         ret = Join(ret, tmp);
      } else if (i == x.size() - 2) {
         n = Norm(x[i + 1] - x[i]) / di;
         std::vector<std::vector<T>> tmp = SubdivideByStep(x[i] + (x[i + 1] - x[i]) / (2. * n), x[i + 1], di);
         ret = Join(ret, tmp);
      } else if (i > 0 && i < x.size() - 2) {
         n = Norm(x[i + 1] - x[i]) / di;
         std::vector<std::vector<T>> tmp = SubdivideByStep(x[i] + (x[i + 1] - x[i]) / (2. * n), x[i + 1] - (x[i + 1] - x[i]) / (2 * n), di);
         ret = Join(ret, tmp);
      }
   }
   return ret;
};
template <typename T>
std::vector<std::vector<T>> Subdivide(const std::vector<std::vector<T>> &x, const int n) {
   T len(0.);
   for (size_t i = 0; i < x.size() - 1; i++)
      len += Norm(x[i + 1] - x[i]);
   T di = len / n;
   int counter(0);
   while (counter < 100000) {
      std::vector<std::vector<T>> tmp = SubdivideByStep(x, di);
      if (tmp.size() == n)
         return SubdivideByStep(x, di);
      else if (tmp.size() > n) {
         di = di * 1.01;
      } else if (tmp.size() < n) {
         di = di * 0.99;
      };
      counter++;
   }
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorOff << std::endl;
   std::cout << Red << "can not subdivide by " << n << " counter = " << counter << colorOff << std::endl;
   abort();
   return SubdivideByStep(x, di);
};
template <typename T>
std::vector<std::vector<T>> SubdivideExclude(const std::vector<std::vector<T>> &x, const int n) {
   T len(0.);
   for (size_t i = 0; i < x.size() - 1; i++)
      len += Norm(x[i + 1] - x[i]);
   T di = len / n;
   int counter(0);
   while (counter < 100000) {
      std::vector<std::vector<T>> tmp = SubdivideByStepExclude(x, di);
      if (tmp.size() == n)
         return tmp;
      else if (tmp.size() > n) {
         di = di * 1.01;
      } else if (tmp.size() < n) {
         di = di * 0.999;
      };
      counter++;
   }
   std::cout << Red << __func__ << ":" << colorOff << std::endl;
   std::cout << Red << "can not subdivide by " << n << " counter = " << counter << colorOff << std::endl;
   abort();
   return SubdivideByStepExclude(x, di);
};

#endif