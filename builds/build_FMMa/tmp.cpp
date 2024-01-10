//////////////////////////////////////////////////////////
#ifndef INCL_fundamental
#define INCL_fundamental

#include <stdlib.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#ifndef _USE_MATH_DEFINES
   #define _USE_MATH_DEFINES
   #include <cmath>
   #ifndef M_PI
      #define M_PI (3.14159265358979323846)
   #endif
#endif

/////// emscripten //////
#ifdef __EMSCRIPTEN__
   #include <emscripten/bind.h>
#endif
/////////////////////////
// #ifndef INCL_lapacke
// #define INCL_lapacke
// #include "lapacke.h"
// #include <map>// for OPENCL
// #endif
///////////////////////////////////////////////////////////
// manipulater
const std::string black("\033[0;30m");
const std::string red("\033[0;31m");
const std::string green("\033[0;32m");
const std::string yellow("\033[0;33m");
const std::string blue("\033[0;34m");
const std::string magenta("\033[0;35m");
const std::string cyan("\033[0;36m");
const std::string white("\033[0;37m");

const std::string Black("\033[1;30m");
const std::string Red("\033[1;31m");
const std::string Green("\033[1;32m");
const std::string Yellow("\033[1;33m");
const std::string Blue("\033[1;34m");
const std::string Magenta("\033[1;35m");
const std::string Cyan("\033[1;36m");
const std::string White("\033[1;37m");

const std::string _black("\033[4;0;30m");
const std::string _red("\033[4;0;31m");
const std::string _green("\033[4;0;32m");
const std::string _yellow("\033[4;0;33m");
const std::string _blue("\033[4;0;34m");
const std::string _magenta("\033[4;0;35m");
const std::string _cyan("\033[4;0;36m");
const std::string _white("\033[4;0;37m");

const std::string _Black("\033[4;1;30m");
const std::string _Red("\033[4;1;31m");
const std::string _Green("\033[4;1;32m");
const std::string _Yellow("\033[4;1;33m");
const std::string _Blue("\033[4;1;34m");
const std::string _Magenta("\033[4;1;35m");
const std::string _Cyan("\033[4;1;36m");
const std::string _White("\033[4;1;37m");
const std::string reset("\033[0m");

///////////////////////////////////////////////////////
// 20200526
std::vector<int> stoi(const std::vector<std::string> &vec) {
   std::vector<int> ret(vec.size());
   std::transform(vec.begin(), vec.end(), ret.begin(), [](const std::string &str) { return std::stoi(str); });
   return ret;
};
std::vector<std::vector<int>> stoi(const std::vector<std::vector<std::string>> &vec) {
   std::vector<std::vector<int>> ret(0);
   for (const auto &v : vec)
      ret.push_back(stoi(v));
   return ret;
};
std::vector<std::vector<std::vector<int>>> stoi(const std::vector<std::vector<std::vector<std::string>>> &vec) {
   std::vector<std::vector<std::vector<int>>> ret(0);
   for (const auto &v : vec)
      ret.push_back(stoi(v));
   return ret;
};
std::vector<double> stod(const std::vector<std::string> &vec) {
   std::vector<double> ret(vec.size());
   std::transform(vec.begin(), vec.end(), ret.begin(), [](const std::string &str) { return std::stod(str); });
   return ret;
};
std::vector<std::vector<double>> stod(const std::vector<std::vector<std::string>> &vec) {
   std::vector<std::vector<double>> ret(0);
   for (const auto &v : vec)
      ret.push_back(stod(v));
   return ret;
};
std::vector<std::vector<std::vector<double>>> stod(const std::vector<std::vector<std::vector<std::string>>> &vec) {
   std::vector<std::vector<std::vector<double>>> ret(0);
   for (const auto &v : vec)
      ret.push_back(stod(v));
   return ret;
};
///////////////////////////////////////////////////////
int InverseQuotientMod(const std::vector<int> &row_col, const int b) {
   return row_col[0] * b + row_col[1];
};
std::vector<int> QuotientMod(const int a, const int b) {
   return {(int)a / b, a % b};
};
//==========================================================
template <class T>
std::vector<T> Take(std::vector<T> &vec, const std::vector<int> &beg_end_skip) {
   std::vector<T> ret;
   switch (beg_end_skip.size()) {
      case 2:
         if (beg_end_skip[1] < 0) {
            for (auto i = beg_end_skip[0]; i < vec.size() + beg_end_skip[1]; i++)
               ret.push_back(vec[i]);
            return ret;
         } else {
            for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
               ret.push_back(vec[i]);
            return ret;
         }
      case 3:
         if (beg_end_skip[1] < 0) {
            for (auto i = beg_end_skip[0]; i < vec.size() + beg_end_skip[1]; i++)
               if (i % beg_end_skip[2] == 0)
                  ret.push_back(vec[i]);
            return ret;
         } else {
            for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
               if (i % beg_end_skip[2] == 0)
                  ret.push_back(vec[i]);
            return ret;
         }
      default:
         std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
         abort();
         return vec;
   };
};
template <class T>
std::vector<std::vector<T>> Take(std::vector<std::vector<T>> &mat, const std::vector<int> &beg_end_skip) {
   std::vector<std::vector<T>> ret;
   switch (beg_end_skip.size()) {
      case 2:
         if (beg_end_skip[1] < 0) {
            for (auto i = beg_end_skip[0]; i < mat.size() + beg_end_skip[1]; i++) {
               std::vector<T> vec;
               for (auto j = 0; j < mat[i].size(); j++)
                  vec.push_back(mat[i][j]);
               ret.push_back(vec);
            };
            return ret;
         } else {
            for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++) {
               std::vector<T> vec;
               for (auto j = 0; j < mat[i].size(); j++)
                  vec.push_back(mat[i][j]);
               ret.push_back(vec);
            };
            return ret;
         }
      case 3:
         if (beg_end_skip[1] < 0) {
            for (auto i = beg_end_skip[0]; i < mat.size().beg_end_skip[1]; i++) {
               if (i % beg_end_skip[2] == 0) {
                  std::vector<T> vec;
                  for (auto j = 0; j < mat[i].size(); j++)
                     vec.push_back(mat[i][j]);
                  ret.push_back(vec);
               }
            };
            return ret;
         } else {
            for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++) {
               if (i % beg_end_skip[2] == 0) {
                  std::vector<T> vec;
                  for (auto j = 0; j < mat[i].size(); j++)
                     vec.push_back(mat[i][j]);
                  ret.push_back(vec);
               }
            };
            return ret;
         }
      default:
         std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
         abort();
         return mat;
   };
};
template <class T>
std::vector<std::vector<T>> Take(std::vector<std::vector<T>> &mat, const std::vector<int> &beg_end_skip, const std::vector<int> &beg_end_skip2) {
   std::vector<std::vector<T>> ret;
   switch (beg_end_skip.size()) {
      case 2:
         if (beg_end_skip[1] < 0) {
            for (auto i = beg_end_skip[0]; i < mat.size() + beg_end_skip[1]; i++) {
               std::vector<T> vec;
               for (auto j = beg_end_skip[0]; j < beg_end_skip2[1]; j++)
                  vec.push_back(mat[i][j]);
               ret.push_back(vec);
            };
            return ret;
         } else {
            for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++) {
               std::vector<T> vec;
               for (auto j = beg_end_skip[0]; j < beg_end_skip2[1]; j++)
                  vec.push_back(mat[i][j]);
               ret.push_back(vec);
            };
            return ret;
         }
      case 3:
         if (beg_end_skip[1] < 0) {
            for (auto i = beg_end_skip[0]; i < mat.size() + beg_end_skip[1]; i++)
               if (i % beg_end_skip[2] == 0)
                  ret.push_back(Take(mat[i], beg_end_skip2));
            return ret;
         } else {
            for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
               if (i % beg_end_skip[2] == 0)
                  ret.push_back(Take(mat[i], beg_end_skip2));
            return ret;
         }
      default:
         std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
         abort();
         return mat;
   };
};
//==========================================================
std::vector<std::vector<double>> Import(const std::string &fname) {
   std::ifstream fin(fname);
   std::vector<std::vector<double>> mat;
   std::string str;
   double tmp;
   while (getline(fin, str)) {
      std::stringstream ss;
      ss << str;
      std::vector<double> row;
      while (ss >> tmp)
         row.push_back(tmp);
      mat.push_back(row);
   }
   fin.close();
   return mat;
};
// example
// std::vector<std::vector<double>> mat=Import("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");
//=========================================================
template <typename T>
std::vector<T> Flatten(const std::vector<std::vector<T>> &mat) {
   std::vector<T> ret(0);
   for (size_t i = 0; i < mat.size(); i++)
      for (size_t j = 0; j < mat[i].size(); j++) {
         if (!mat[i].empty())
            ret.push_back(mat[i][j]);
      };
   return ret;
};
//==========================================================
// Template <typename T>
// std::ostream &operator<<(std::ostream &stream, const std::map<std::string, std::string> &v){
//   stream << "{";
//   for(size_t i=0; i<v.size()-1; i++)
//     stream << v[i] << "," << v[i]->second;
//   stream << *v.rbegin() << "}";
//   return stream;
// };
// std::ostream &operator<<(std::ostream &stream, const std::vector<std::string> &v){
//   #if defined FULL_DEBUG
//   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
//   #endif
//   std::string str="{";
//   for(auto i=0; i<v.size()-1; i++)
//     str = str + v[i] + ",";
//   str = str + *v.rbegin() + "}";
//   stream << str;
//   return stream;
// };
std::ostream &operator<<(std::ostream &stream, const std::vector<std::string> &v) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   stream << "{";
   for (auto i = 0; i < v.size() - 1; i++)
      stream << v[i] << ",";
   stream << *v.rbegin();
   stream << "}";
   return stream;
};
template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::vector<T> &v) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   stream << "{";
   for (auto i = 0; i < v.size() - 1; i++)
      stream << v[i] << ",";
   stream << *v.rbegin();
   stream << "}";
   return stream;
};
///////////////////////////////////////////////////////////////
template <class T>
T Min(const std::vector<T> &v) {
   return *std::min_element(std::begin(v), std::end(v));
};
template <class T>
T Min(const std::vector<std::vector<T>> &v) {
   T tmp, ret(0.);
   for (auto i = 0; i < v.size(); i++) {
      tmp = Min(v[i]);
      if (ret > tmp)
         ret = tmp;
   }
   return ret;
};

template <class T>
T Max(const std::vector<T> &v) {
   return *std::max_element(std::begin(v), std::end(v));
};
template <class T>
T Max(const std::vector<std::vector<T>> &v) {
   T tmp, ret(0.);
   for (auto i = 0; i < v.size(); i++) {
      tmp = Max(v[i]);
      if (ret < tmp)
         ret = tmp;
   }
   return ret;
};

bool MemberQ(const std::vector<int> &list, const int form) {  // Mathematica like
   bool ret;
   for (auto &v : list)
      if (v == form)
         return true;
   return false;
};

bool MemberQ(const std::vector<std::vector<int>> &list, const std::vector<int> &form) {  // Mathematica like
   bool ret;
   for (auto &v : list) {
      ret = ret || (v == form);
      if (ret)
         return ret;
   }
   return false;
};

double Sum(const std::vector<double>::iterator first, const std::vector<double>::iterator second) {
   return std::accumulate(first, second, 0.);
};
double Sum(const std::vector<double> &v) {
   return std::accumulate(v.cbegin(), v.cend(), 0.);
};
//////////////////////////////////////////////////////////
std::vector<double> operator-(const std::vector<double> &v) {
   std::vector<double> ret(v);
   std::transform(ret.cbegin(), ret.cend(), ret.begin(), std::negate<double>());
   return ret;
};
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>> &v) {
   std::vector<std::vector<double>> ret(v);
   std::transform(ret.cbegin(), ret.cend(), ret.begin(), [](std::vector<double> tmp) { return -tmp; });
   return ret;
};
///////////////////////////////////////////////////
// vector
template <class T>
std::vector<T> operator*(const std::vector<T> &v, const T din) {
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](T tmp) { return tmp * din; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>> &v, const T din) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<T> tmp) { return tmp * din; });
   return ret;
};
// vector
template <class T>
std::vector<T> operator*(const T din, const std::vector<T> &v) {
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](T tmp) { return tmp * din; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<T>> operator*(const T din, const std::vector<std::vector<T>> &v) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<T> tmp) { return tmp * din; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<std::vector<T>>> operator*(const std::vector<std::vector<std::vector<T>>> &v, const T din) {
   std::vector<std::vector<std::vector<T>>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<std::vector<T>> tmp) { return tmp * din; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<std::vector<T>>> operator*(const T din, const std::vector<std::vector<std::vector<T>>> &v) {
   std::vector<std::vector<std::vector<T>>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<std::vector<T>> tmp) { return tmp * din; });
   return ret;
};
// vector x vector
template <class T>
std::vector<T> operator*(const std::vector<T> &v, const std::vector<T> &w) {
#if defined debug_fundamental
   if (v.size() != w.size()) {
      std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
      std::cout << Red << ": vectors have different sizes" << reset << std::endl;
      std::cout << Red << std::setprecision(5) << v << reset << std::endl;
      std::cout << Red << std::setprecision(5) << w << reset << std::endl;
      abort();
   }
#endif
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](T a, T b) { return a * b; });
   return ret;
};
template <class T>
std::vector<T> &operator*=(std::vector<T> &v, const T w) {
   v = v * w;
   return v;
};
template <class T>
std::vector<T> &operator*=(std::vector<T> &v, const std::vector<T> &w) {
   v = v * w;
   return v;
};
// vector x matrix
template <class T>
std::vector<std::vector<T>> operator*(const std::vector<T> &w, const std::vector<std::vector<T>> &v) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, T b) { return b * a; });
   return ret;
};
// matrix x vector
template <class T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, T b) { return a * b; });
   return ret;
};
// matrix x matrix
template <class T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, std::vector<T> b) { return a * b; });
   return ret;
};
template <class T>
void operator*=(std::vector<T> &v, const std::vector<T> &w) {
   v = v * w;
};
template <class T>
std::vector<std::vector<T>> &operator*=(std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   v = v * w;
   return v;
};
template <class T>
std::vector<std::vector<T>> &operator*=(std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   v = v * w;
   return v;
};
/////////////////////////////////////////////////
// vector
template <class T>
std::vector<T> operator/(const std::vector<T> &v, const T din) {
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](T tmp) { return tmp / din; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<T>> operator/(const std::vector<std::vector<T>> &v, const T din) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<T> tmp) { return tmp / din; });
   return ret;
};
// vector
template <class T>
std::vector<T> operator/(const T din, const std::vector<T> &v) {
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](T tmp) { return din / tmp; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<T>> operator/(const T din, const std::vector<std::vector<T>> &v) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<T> tmp) { return din / tmp; });
   return ret;
};
// vector x vector
template <class T>
std::vector<T> operator/(const std::vector<T> &v, const std::vector<T> &w) {
#if defined debug_fundamental
   if (v.size() != w.size()) {
      std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
      std::cout << Red << ": vectors have different sizes" << reset << std::endl;
      std::cout << Red << std::setprecision(5) << v << reset << std::endl;
      std::cout << Red << std::setprecision(5) << w << reset << std::endl;
      abort();
   }
#endif
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](T a, T b) { return a / b; });
   return ret;
};
template <class T>
std::vector<T> &operator/=(std::vector<T> &v, const T w) {
   v = v / w;
   return v;
};
template <class T>
std::vector<T> &operator/=(std::vector<T> &v, const std::vector<T> &w) {
   v = v / w;
   return v;
};
// vector x matrix
template <class T>
std::vector<std::vector<T>> operator/(const std::vector<T> &w, const std::vector<std::vector<T>> &v) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, T b) { return b / a; });
   return ret;
};
// matrix x vector
template <class T>
std::vector<std::vector<T>> operator/(const std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, T b) { return a / b; });
   return ret;
};
// matrix x matrix
template <class T>
std::vector<std::vector<T>> operator/(const std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, std::vector<T> b) { return a / b; });
   return ret;
};
template <class T>
std::vector<std::vector<T>> &operator/=(std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   v = v / w;
   return v;
};
template <class T>
std::vector<std::vector<T>> &operator/=(std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   v = v / w;
   return v;
};
/////////////////////////////////////////////////
// vector
template <class T>
std::vector<T> operator-(const std::vector<T> &v, const T din) {
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](T tmp) { return tmp - din; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>> &v, const T din) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<T> tmp) { return tmp - din; });
   return ret;
};
// vector
template <class T>
std::vector<T> operator-(const T din, const std::vector<T> &v) {
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](T tmp) { return din - tmp; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<T>> operator-(const T din, const std::vector<std::vector<T>> &v) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<T> tmp) { return din - tmp; });
   return ret;
};
// vector x vector
template <class T>
std::vector<T> operator-(const std::vector<T> &v, const std::vector<T> &w) {
#if defined debug_fundamental
   if (v.size() != w.size()) {
      std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
      std::cout << Red << ": vectors have different sizes" << reset << std::endl;
      std::cout << Red << std::setprecision(5) << v << reset << std::endl;
      std::cout << Red << std::setprecision(5) << w << reset << std::endl;
      abort();
   }
#endif

   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](T a, T b) { return a - b; });
   return ret;
};
template <class T>
std::vector<T> &operator-=(std::vector<T> &v, const T w) {
   v = v - w;
   return v;
};
template <class T>
std::vector<T> &operator-=(std::vector<T> &v, const std::vector<T> &w) {
   v = v - w;
   return v;
};
// vector x matrix
template <class T>
std::vector<std::vector<T>> operator-(const std::vector<T> &w, const std::vector<std::vector<T>> &v) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, T b) { return b - a; });
   return ret;
};
// matrix x vector
template <class T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, T b) { return a - b; });
   return ret;
};
// matrix x matrix
template <class T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, std::vector<T> b) { return a - b; });
   return ret;
};
template <class T>
std::vector<std::vector<T>> &operator-=(std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   v = v - w;
   return v;
};
template <class T>
std::vector<std::vector<T>> &operator-=(std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   v = v - w;
   return v;
};
/////////////////////////////////////////////////
// vector
template <class T>
std::vector<T> operator+(const std::vector<T> &v, const T din) {
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](T tmp) { return tmp + din; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>> &v, const T din) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<T> tmp) { return tmp + din; });
   return ret;
};
// vector
template <class T>
std::vector<T> operator+(const T din, const std::vector<T> &v) {
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](T tmp) { return din + tmp; });
   return ret;
};
// matrix
template <class T>
std::vector<std::vector<T>> operator+(const T din, const std::vector<std::vector<T>> &v) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), std::back_inserter(ret), [&din](std::vector<T> tmp) { return din + tmp; });
   return ret;
};
// vector x vector
template <class T>
std::vector<T> operator+(const std::vector<T> &v, const std::vector<T> &w) {
#if defined debug_fundamental
   if (v.size() != w.size()) {
      std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
      std::cout << Red << ": vectors have different sizes" << reset << std::endl;
      std::cout << Red << std::setprecision(5) << v << reset << std::endl;
      std::cout << Red << std::setprecision(5) << w << reset << std::endl;
      abort();
   }
#endif
   std::vector<T> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](T a, T b) { return a + b; });
   return ret;
};
template <class T>
std::vector<T> &operator+=(std::vector<T> &v, const T w) {
   v = v + w;
   return v;
};
template <class T>
std::vector<T> &operator+=(std::vector<T> &v, const std::vector<T> &w) {
   v = v + w;
   return v;
};
// vector x matrix
template <class T>
std::vector<std::vector<T>> operator+(const std::vector<T> &w, const std::vector<std::vector<T>> &v) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, T b) { return b + a; });
   return ret;
};
// matrix x vector
template <class T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, T b) { return a + b; });
   return ret;
};
// matrix x matrix
template <class T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   std::vector<std::vector<T>> ret;
   std::transform(v.cbegin(), v.cend(), w.cbegin(), std::back_inserter(ret), [](std::vector<T> a, std::vector<T> b) { return a + b; });
   return ret;
};
template <class T>
std::vector<std::vector<T>> &operator+=(std::vector<std::vector<T>> &v, const std::vector<T> &w) {
   v = v + w;
   return v;
};
template <class T>
std::vector<std::vector<T>> &operator+=(std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w) {
   v = v + w;
   return v;
};
/////////////////////////////////////////////////////////////
double Mean(const std::vector<double> &v) {
   return accumulate(v.cbegin(), v.cend(), 0.0) / v.size();
};
//==========================================================
struct Timer {
   std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
   Timer() {
      start = std::chrono::high_resolution_clock::now();
   };
   std::chrono::duration<double> elapsed() {
      end = std::chrono::high_resolution_clock::now();
      return end - start;
   };
};
//==========================================================
double sgn(const double x) { return (x > 0.) ? 1. : (x < 0. ? -1. : 0.); };
// template<class T>
// std::vector<T> Transpose(const std::vector< T >& mat){
//   std::vector<T> ans;
//   for(auto j=0; j<mat[i].size(); j++)
//     {
//       for(auto i=0; i<mat.size(); i++)
// 	{
// 	  tmp.push_back(mat[i][j]);
// 	}
//       ans.push_back(tmp);
//     }
//   return ans;
// };
std::vector<std::vector<double>> Transpose(const std::vector<std::vector<double>> &mat) {
   std::vector<std::vector<double>> ans(mat[0].size(), std::vector<double>(mat.size()));
   for (auto i = 0; i < mat.size(); i++)
      for (auto j = 0; j < mat[i].size(); j++)
         ans[j][i] = mat[i][j];
   return ans;
};
std::vector<std::vector<std::vector<double>>> Transpose(const std::vector<std::vector<std::vector<double>>> &mat) {
   std::vector<std::vector<std::vector<double>>> ans(mat[0][0].size(), std::vector<std::vector<double>>(mat[0].size(), std::vector<double>(mat.size())));
   for (auto i = 0; i < mat.size(); i++)
      for (auto j = 0; j < mat[i].size(); j++)
         for (auto k = 0; k < mat[i][j].size(); k++)
            ans[k][j][i] = mat[i][j][k];
   return ans;
};

std::vector<std::vector<double>> TensorProduct(const std::vector<double> &vec1, const std::vector<double> &vec2) {
   std::vector<std::vector<double>> ret(vec1.size(), std::vector<double>(vec2.size()));
   for (auto m = 0; m < vec1.size(); m++)
      for (auto j = 0; j < vec2.size(); j++)
         ret[m][j] = vec1[m] * vec2[j];
   return ret;
};

std::vector<std::vector<std::vector<double>>> TensorProductSet(const std::vector<double> &vec1, const std::vector<double> &vec2) {
   std::vector<std::vector<std::vector<double>>> ret(vec1.size(), std::vector<std::vector<double>>(vec2.size(), std::vector<double>(2, 0)));
   for (auto m = 0; m < vec1.size(); m++)
      for (auto j = 0; j < vec2.size(); j++)
         ret[m][j] = {vec1[m], vec2[j]};
   return ret;
};

template <class T>
T Dot(const std::vector<T> &vec1, const std::vector<T> &vec2) {
#if defined debug_fundamental
   if (vec1.size() == 2 && vec2.size() == 3) {
      /*この場合は，vec1の３要素目を0パディングする*/
   } else if (vec1.size() != vec2.size()) {
      std::cout << Red << "Dot: vectors have different sizes" << reset << std::endl;
      std::cout << Red << std::setprecision(5) << vec1 << reset << std::endl;
      std::cout << Red << std::setprecision(5) << vec2 << reset << std::endl;
      abort();
   }
#endif
   return std::inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0.);
};
//// for pointer
template <class T>
T Dot(const std::vector<T *> &vec1, const std::vector<T *> &vec2) {
#if defined debug_fundamental
   if (vec1.size() == 2 && vec2.size() == 3) {
      /*この場合は，vec1の３要素目を0パディングする*/
   } else if (vec1.size() != vec2.size()) {
      std::cout << Red << "Dot: vectors have different sizes" << reset << std::endl;
      std::cout << Red << std::setprecision(5) << vec1 << reset << std::endl;
      std::cout << Red << std::setprecision(5) << vec2 << reset << std::endl;
      abort();
   }
#endif

   T ans(0.);
   for (auto i = 0; i < vec1.size(); i++)
      ans += *vec1[i] * *vec2[i];
   return ans;
};
// std::vector<double> Dot(const std::vector< std::vector<double> > mat, const std::vector<double> vec){
//   std::vector<double> ans(vec.size(), 0.);
//   for(size_t i=0; i<mat.size(); i++)
//     for(size_t j=0; j<vec.size(); j++)
//       ans[i] += mat[i][j] * vec[j];
//   return ans;
// };
template <class T>
std::vector<T> Dot(const std::vector<std::vector<T>> &mat, const std::vector<T> &vec) {
   std::vector<T> ans(mat.size());
   for (auto i = 0; i < mat.size(); i++)
      ans[i] = Dot(mat[i], vec);
   return ans;
};
template <class T>
std::vector<T> Dot(const std::vector<T> &vec, const std::vector<std::vector<T>> &mat) {
   return Dot(Transpose(mat), vec);
};
template <class T>
std::vector<std::vector<T>> Dot(const std::vector<std::vector<T>> &mat1, const std::vector<std::vector<T>> &mat2) {
   if (mat1[0].size() != mat2.size()) {
      std::cout << __func__ << ": passed variables have dimensions that can not be computed" << std::endl;
      abort();
   }
   std::vector<std::vector<T>> ans(mat1.size(), std::vector<T>(mat2[0].size(), 0.));
   for (auto x = 0; x < mat1.size(); x++)
      for (auto y = 0; y < mat2[0].size(); y++)
         for (auto j = 0; j < mat2.size(); j++)
            ans[x][y] += mat1[x][j] * mat2[j][y];
   return ans;
};
double Rot(const std::vector<double> vec1, const std::vector<double> vec2) {
   return vec1[0] * vec2[1] - vec1[1] * vec2[0];
};
std::vector<std::vector<double>> Inv(const std::vector<std::vector<double>> &mat) {
   std::vector<std::vector<double>> ans(mat.size(), std::vector<double>(mat[0].size(), 0.));
   double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
   ans[1][1] = mat[0][0] / det;
   ans[0][1] = -mat[0][1] / det;
   ans[1][0] = -mat[1][0] / det;
   ans[0][0] = mat[1][1] / det;
   return ans;
};
//==========================================================
// vector operators
template <class T>
std::vector<T> Cross(const std::vector<T> &A) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   return {-A[1], A[0]};
};
template <class T>
std::vector<T> Cross(const std::vector<T> &A, const std::vector<T> &X) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   if (A.size() == 3)
      return {A[1] * X[2] - A[2] * X[1],
              A[2] * X[0] - A[0] * X[2],
              A[0] * X[1] - A[1] * X[0]};
   else if (A.size() == 2)
      return Cross(std::vector<T>{A[0], A[1], 0.}, std::vector<T>{X[0], X[1], 0.});
   else {
      std::cout << Red << "Invalid Vector is passed" << reset << std::endl;
      std::cout << Red << std::setprecision(5) << A << reset << std::endl;
      std::cout << Red << std::setprecision(5) << X << reset << std::endl;
      abort();
   }
};
//==========================================================
// template <class Ttype> Ttype Abs(const std::vector<Ttype>& vec){
//   Ttype tmp(0);
//   for(size_t i=0; i<vec.size(); i++)
//     tmp += vec[i]*vec[i];
//   return (Ttype)std::sqrt(tmp);
// };
/*std::vector<double> abs(const std::vector<double>& vec){
  std::vector<double> ret(vec.size());
  for(size_t i=0; i<vec.size(); i++)
  ret[i] = abs(vec[i]);
  return ret;
  };*/
std::vector<double> log10(const std::vector<double> &vec) {
   std::vector<double> ret(vec.size());
   for (auto i = 0; i < vec.size(); i++)
      ret[i] = std::log10(vec[i]);
   return ret;
};
std::vector<double> log(const std::vector<double> &vec) {
   std::vector<double> ret(vec.size());
   for (size_t i = 0; i < vec.size(); i++)
      ret[i] = std::log(vec[i]);
   return ret;
};
template <class T>
T Norm(const std::vector<T> &vec) {
   return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.));
};
template <class T>
T Abs(const std::vector<double> &vec) {
   T tmp(0);
   for (size_t i = 0; i < vec.size(); i++)
      tmp += vec[i] * vec[i];
   return std::sqrt(tmp);
};
//==========================================================
template <class T>
T VectorAngle(const std::vector<T> &X1, const std::vector<T> &X2, const std::vector<T> &X0) {
   return std::atan2(Rot(std::vector<T>{X1[0] - X0[0], X1[1] - X0[1]}, std::vector<T>{X2[0] - X0[0], X2[1] - X0[1]}),
                     Dot(std::vector<T>{X1[0] - X0[0], X1[1] - X0[1]}, std::vector<T>{X2[0] - X0[0], X2[1] - X0[1]}));
   // std::vector<T> vec1 = {x1,y1}, vec2 = {x2,y2};
   // T tmp = Dot(vec1,vec2) / ( Abs(vec1)*Abs(vec2) );
   // return sgn( Rot(vec1,vec2) )*acos(tmp);};
};
template <class T>
T VectorAngle(const std::vector<T> &X1, const std::vector<T> &X2) {
   return VectorAngle(X1, X2, {0., 0.});
};
template <class T>
T VectorAngleDirected(const std::vector<T> &X1, const std::vector<T> &X2) {
   T a = VectorAngle(X1, X2, {0., 0.});
   return (a < 0) ? (M_2_PI - a) : a;
};
//==========================================================
double SgLog(const double u, const double hsig, const double maxmin) {
   return hsig + sgn(u) * exp(maxmin * (std::abs(u) - 1.));
};
double InvSgLog(const double h, const double hsig, const double maxmin) {
   return sgn(h - hsig) * std::log(std::abs(h - hsig)) / maxmin + 1.;
};
std::vector<double> SgLog(const std::vector<double> &h, const double hsig, const double maxmin) {
   std::vector<double> ret;
   std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [hsig, maxmin](double tmp) { return SgLog(tmp, hsig, maxmin); });
   return ret;
};
double DSgLog(const double u, const double hsig, const double maxmin) {
   return maxmin * sgn(u) * exp(maxmin * (std::abs(u) - 1.));
};  // 単調増加
std::vector<double> DSgLog(const std::vector<double> &h, const double hsig, const double maxmin) {
   std::vector<double> ret;
   std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [hsig, maxmin](double tmp) { return DSgLog(tmp, hsig, maxmin); });
   return ret;
};
//============================================================
double Sg(const double h, const double h_sig, const double beta) {
   return sgn(h) * pow(std::abs(h), beta) + h_sig;
};
double InvSg(const double h, const double h_sig, const double beta) {
   return sgn(h - h_sig) * pow(std::abs(h - h_sig), 1. / beta);
};
std::vector<double> Sg(const std::vector<double> &h, const double h_sig, const double beta) {
   std::vector<double> ret;
   std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [h_sig, beta](double tmp) { return Sg(tmp, h_sig, beta); });
   return ret;
};
double DSg(const double h, const double h_sig, const double beta) {
   return beta * pow(std::abs(h), beta - 1.);
};  // 単調増加
std::vector<double> DSg(const std::vector<double> &h, const double h_sig, const double beta) {
   std::vector<double> ret;
   std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [h_sig, beta](double tmp) { return DSg(tmp, h_sig, beta); });
   return ret;
};
//==========================================================
double linspace(const std::vector<double> &v, int size, int i) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   return v[0] + (v[1] - v[0]) * i / (double)(size - 1.);
};
/* double init_position(int size, int i, std::vector<double> v, double c, double beta){ */
/*   cout << "max = " << v[1] << ", min = " << v[0] << ", d = " << v[1]-v[0] << endl; */
/*   return Sg( v[0] + (v[1] - v[0]) * i/(double)(size - 1.), c, beta);; */
/* }; */
double linspace(const std::vector<double> &v, int size, int i, double c, double beta) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   std::vector<double> w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
   return Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
};
void linspace(const std::vector<double> &v, int size, std::vector<double> &O_X) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   O_X.resize(size);
   for (int i = 0; i < size; i++)
      O_X[i] = v[0] + (v[1] - v[0]) * i / ((double)size - 1.);
};
void linspace(const std::vector<double> &v, int size, std::vector<double> &O_X, double c, double beta) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   std::vector<double> w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
   O_X.resize(size);
   for (int i = 0; i < size; i++)
      O_X[i] = Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
};
std::vector<double> linspace(const std::vector<double> &v, int size) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   std::vector<double> O_X(size);
   for (int i = 0; i < size; i++)
      O_X[i] = v[0] + (v[1] - v[0]) * i / (double)(size - 1.);
   return O_X;
};
std::vector<double> linspace(const std::vector<double> &v, int size, double c, double beta) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   std::vector<double> w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
   std::vector<double> O_X(size);
   for (int i = 0; i < size; i++)
      O_X[i] = Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
   return O_X;
};
void linspace(const std::vector<std::vector<double>> &mat, const std::vector<int> size, std::vector<double> &O_X) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   std::vector<double> v;
   for (size_t n = 0; n < mat.size(); n++) {
      linspace(mat[n], size[n], v);
      O_X.insert(std::end(O_X), std::begin(v), std::end(v));
   }
};
////////// Join has been modified on 20200411
template <typename T>
std::vector<T> Join(std::vector<T> a, const std::vector<T> &b) {
   a.reserve(a.size() + b.size());
   a.insert(a.end(), b.begin(), b.end());
   return a;
};
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
template <typename T>
std::vector<T> Subdivide(const T xmin, const T xmax, const int n) {
   T dx = (xmax - xmin) / n;
   std::vector<T> ret(n + 1);
   for (int i = 0; i < n + 1; i++)
      ret[i] = i * dx + xmin;
   return ret;
};
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
         std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
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
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
   std::cout << Red << "can not subdivide by " << n << " counter = " << counter << reset << std::endl;
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
   std::cout << Red << __func__ << ":" << reset << std::endl;
   std::cout << Red << "can not subdivide by " << n << " counter = " << counter << reset << std::endl;
   abort();
   return SubdivideByStepExclude(x, di);
};
//==========================================================
void gauleg(const double x1, const double x2, std::vector<double> &x, std::vector<double> &w) {
#if defined FULL_DEBUG
   std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
   const double EPS = 1.0e-14;
   double z1, z, xm, xl, pp, p3, p2, p1;
   int n = x.size();
   int m = (n + 1) / 2;
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);
   for (int i = 0; i < m; i++) {
      z = std::cos(3.1415926535897932385 * (i + 0.75) / (n + 0.5));
      do {
         p1 = 1.0;
         p2 = 0.0;
         for (int j = 0; j < n; j++) {
            p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
         }
         pp = n * (z * p1 - p2) / (z * z - 1.0);
         z1 = z;
         z = z1 - p1 / pp;
      } while (std::abs(z - z1) > EPS);
      x[i] = xm - xl * z;
      x[n - 1 - i] = xm + xl * z;
      w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
      w[n - 1 - i] = w[i];
   }
};
//===========================================================
// iはNコのデータがある場合，0<=i<=N-1となる
// i=N-1の場合，q[N-1]=1.なのでh=1.はq[N-1]=1.<= h && h < q[N]

//============= 20191222 Interpolation =============
double Bspline(const double h, const std::vector<double> &q, const int i, const int K) {
   // switch(K)
   //   {
   //   case 1 :
   //     return (q[i] <= h && h < q[i+1]) ? 1. : 0.;
   //   default :
   //     return (h - q[i])/(q[i + K - 1] - q[i]) * Bspline(h,q,i,K-1)
   // 	+(q[i+K] - h)/( q[i+K] - q[i+1]) * Bspline(h,q,i+1,K-1);
   //   };

   switch (K) {
      case 1:
         return (q[i] <= h && h < q[i + 1]) ? 1. : 0.;
      case 2:
         return ((q[i] <= h && h < q[i + 1]) ? (h - q[i]) / (q[i + 1] - q[i]) : 0.) + ((q[i + 1] <= h && h < q[i + 2]) ? (q[i + 2] - h) / (q[i + 2] - q[i + 1]) : 0.);
      case 3:
         return ((q[i] <= h && h < q[i + 1]) ? (h - q[i]) / (q[i + 2] - q[i]) * (h - q[i]) / (q[i + 1] - q[i]) : 0.) + ((q[i + 1] <= h && h < q[i + 2]) ? ((h - q[i]) / (q[i + 2] - q[i]) * (q[i + 2] - h) + (q[i + 3] - h) / (q[i + 3] - q[i + 1]) * (h - q[i + 1])) / (q[i + 2] - q[i + 1]) : 0.) + ((q[i + 2] <= h && h < q[i + 3]) ? (q[i + 3] - h) / (q[i + 3] - q[i + 1]) * (q[i + 3] - h) / (q[i + 3] - q[i + 2]) : 0.);
      default:
         return (h - q[i]) / (q[i + K - 1] - q[i]) * Bspline(h, q, i, K - 1) + (q[i + K] - h) / (q[i + K] - q[i + 1]) * Bspline(h, q, i + 1, K - 1);
   };

   // switch(K)
   //   {
   //   case 1 :
   //     return (q[i] <= h && h < q[i+1]) ? 1. : 0.;
   //   case 2 :
   //     return ((q[i] <= h && h < q[i+1]) ? (h - q[i])/(q[i+1] - q[i]) : 0.)
   // 	+ ((q[i+1] <= h && h < q[i+2]) ? (q[i+2] - h)/( q[i+2] - q[i+1]) : 0.);
   //   case 3 :
   //     return (h - q[i])/(q[i+2] - q[i]) * Bspline(h,q,i,2)
   // 	+(q[i+3] - h)/( q[i+3] - q[i+1]) * Bspline(h,q,i+1,2);
   //   default :
   //     return (h - q[i])/(q[i + K - 1] - q[i]) * Bspline(h,q,i,K-1)
   // 	+(q[i+K] - h)/( q[i+K] - q[i+1]) * Bspline(h,q,i+1,K-1);
   //   };
};

std::vector<double> Bspline(const std::vector<double> &h, const std::vector<double> &q, const int i, const int K) {
   std::vector<double> ret(h.size());
   for (auto k = 0; k < h.size(); k++)
      ret[k] = Bspline(h[k], q, i, K);
   return ret;
};

std::vector<std::vector<double>> Bspline(const std::vector<double> &h, const std::vector<double> &q, const int K) {
   std::vector<std::vector<double>> ret(h.size(), std::vector<double>(h.size() + K, 0));
   for (auto k = 0; k < ret.size(); k++)
      for (auto i = 0; i < ret[k].size(); i++)
         ret[k][i] = Bspline(h[k], q, i, K);
   return ret;
};

double D_Bspline(const double h, const std::vector<double> &q, const int i, const int K) {
   return (K - 1.) * (1. / (q[i + K - 1] - q[i]) * Bspline(h, q, i, K - 1) - 1. / (q[i + K] - q[i + 1]) * Bspline(h, q, i + 1, K - 1));
};
double D_Bspline(const double h, const std::vector<double> &q, const int i, const int K, const int n) {
   if (K == 1 && n > 0) {
      return 0.;
   } else {
      switch (n) {
         case 0:
            return Bspline(h, q, i, K);
         default:
            return (K - 1) * (1. / (q[i + K - 1] - q[i]) * D_Bspline(h, q, i, K - 1, n - 1) - 1. / (q[i + K] - q[i + 1]) * D_Bspline(h, q, i + 1, K - 1, n - 1));
      };
   };
};

struct ludcmp {
   int n;
   std::vector<std::vector<double>> lu;
   std::vector<int> indx;
   double d;
   std::vector<std::vector<double>> aref;

   ludcmp(const std::vector<std::vector<double>> &a) : n(a.size()), lu(a), aref(a), indx(n) {
      const double TINY = 1.0e-40;
      int i, imax, j, k;
      double big, temp;
      std::vector<double> vv(n);
      d = 1.0;
      for (i = 0; i < n; i++) {
         big = 0.0;
         for (j = 0; j < n; j++)
            if ((temp = std::abs(lu[i][j])) > big)
               big = temp;
         if (big == 0.0)
            throw("Singular matrix in LUdcmp");
         vv[i] = 1.0 / big;
      }
      for (k = 0; k < n; k++) {
         big = 0.0;
         for (i = k; i < n; i++) {
            temp = vv[i] * std::abs(lu[i][k]);
            if (temp > big) {
               big = temp;
               imax = i;
            }
         }
         if (k != imax) {
            for (j = 0; j < n; j++) {
               temp = lu[imax][j];
               lu[imax][j] = lu[k][j];
               lu[k][j] = temp;
            }
            d = -d;
            vv[imax] = vv[k];
         }
         indx[k] = imax;
         if (lu[k][k] == 0.0)
            lu[k][k] = TINY;
         for (i = k + 1; i < n; i++) {
            temp = lu[i][k] /= lu[k][k];
            for (j = k + 1; j < n; j++)
               lu[i][j] -= temp * lu[k][j];
         }
      }
   };
   void solve(const std::vector<double> &b, std::vector<double> &x) {
      int i, ii = 0, ip, j;
      double sum;
      if (b.size() != n || x.size() != n)
         throw("solve bad sizes");
      for (i = 0; i < n; i++)
         x[i] = b[i];
      for (i = 0; i < n; i++) {
         ip = indx[i];
         sum = x[ip];
         x[ip] = x[i];
         if (ii != 0)
            for (j = ii - 1; j < i; j++)
               sum -= lu[i][j] * x[j];
         else if (sum != 0.0)
            ii = i + 1;
         x[i] = sum;
      }
      for (i = n - 1; i >= 0; i--) {
         sum = x[i];
         for (j = i + 1; j < n; j++)
            sum -= lu[i][j] * x[j];
         x[i] = sum / lu[i][i];
      }
   };

   void solve(std::vector<std::vector<double>> &b, std::vector<std::vector<double>> &x) {
      int i, j, m = b[0].size();
      if (b.size() != n || x.size() != n || b[0].size() != x.size())
         throw("solve bad sizes");
      std::vector<double> xx(n);
      for (j = 0; j < m; j++) {
         for (i = 0; i < n; i++)
            xx[i] = b[i][j];
         solve(xx, xx);
         for (i = 0; i < n; i++)
            x[i][j] = xx[i];
      }
   };
   void inverse(std::vector<std::vector<double>> &ainv) {
      int i, j;
      ainv.resize(n, std::vector<double>(n, 0));
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++)
            ainv[i][j] = 0.;
         ainv[i][i] = 1.;
      }
      solve(ainv, ainv);
   };

   std::vector<std::vector<double>> Inverse() {
      std::vector<std::vector<double>> ainv(n, std::vector<double>(n, 0));
      int i, j;
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++)
            ainv[i][j] = 0.;
         ainv[i][i] = 1.;
      }
      solve(ainv, ainv);
      return ainv;
   };

   double det() {
      double dd = d;
      for (int i = 0; i < n; i++)
         dd *= lu[i][i];
      return dd;
   };
   void mprove(std::vector<double> &b, std::vector<double> &x) {
      int i, j;
      std::vector<double> r(n);
      for (i = 0; i < n; i++) {
         double sdp = -b[i];
         for (j = 0; j < n; j++)
            sdp += (double)aref[i][j] * (double)x[j];
         r[i] = sdp;
      }
      solve(r, r);
      for (i = 0; i < n; i++)
         x[i] -= r[i];
   };
};

std::vector<std::vector<double>> Inverse(const std::vector<std::vector<double>> &mat) {
   ludcmp LU(mat);
   std::vector<std::vector<double>> ret;
   LU.inverse(ret);
   return ret;
};

// std::vector<double> PeriodicKnot(const std::vector<double>& h, const int K){
//   int s(h.size());
//   std::vector<double> q(s + K);
//   double EPS = 1.E-15;
//   for(auto i=-K; i <= -1; i++)
//     q[i] = h[i+s-1]-(s-1);
//   for(auto i=s; i <= s+K-2; i++)
//     q[i] = h[i-s+1]+(s-1);
//   for(auto i=0; i < s; i++)
//     q[i + K] = h[i];
//   return q;
// };
std::vector<double> UniformKnots(const std::vector<double> &h, const int K) {
   return Subdivide(-.5, .5, h.size() + K - 1);
};
std::vector<double> UniformKnots(const int size, const int K) {
   return UniformKnots(Subdivide(double(-1), double(1), size - 1), K);
};
std::vector<double> OpenUniformKnots(const std::vector<double> &h, const int K) {
   int s(h.size());
   std::vector<double> q(s + K);
   double EPS = 1.E-15;
   for (auto i = 0; i < K; i++) {
      q[i] = h[0] - EPS * (i - K + 1.);
      q[i + s] = h[s - 1] + EPS * i;
   }
   for (auto i = 0; i < s - K; i++)
      q[i + K] = (h[i] + h[(i + K)]) / 2.;

   return q;
};
std::vector<double> OpenUniformKnots(const int size, const int K) {
   return OpenUniformKnots(Subdivide(double(-1), double(1), size - 1), K);
};

std::vector<double> Bspline_knot(const int size, const int K) {
   std::vector<double> h = Subdivide(double(-1), double(1), size - 1), q(size + K);
   double EPS = 1.E-15;
   // for(auto i=0; i < K; i++)
   //   q[i] = h[0]-EPS*(i-K+1.);
   // for(auto i=0; i < size - K; i++)
   //   q[i + K] = (h[i] + h[(i + K)])/2.;
   // for(auto i=size; i < size + K; i++)
   //   q[i] = h[size - 1]+EPS*(i-size);
   for (auto i = 0; i < K; i++) {
      q[i] = h[0] - EPS * (i - K + 1.);
      q[i + size] = h[size - 1] + EPS * i;
   }
   for (auto i = 0; i < size - K; i++)
      q[i + K] = (h[i] + h[i + K]) / 2.;
   return q;
};

std::vector<double> Bspline_knot_periodic(const int size, const int K) {
   std::vector<double> h = Subdivide(double(-1), double(1), size - 1), q(size + 4 * K - 2);
   for (auto i = 0; i < size; i++)
      q[i] = h[i];
   for (auto i = -2 * K; i < 0; i++)
      q[i] = h[i + size - 1] - (size - 1);
   for (auto i = size; i < size + 2 * K - 1; i++)
      q[i] = h[i - size + 1] + (size - 1);
   return q;
};
std::vector<double> Bspline_knot(const std::vector<double> &h, const int K) {
   int s(h.size());
   std::vector<double> q(s + K);
   double EPS = 1.E-15;
   for (auto i = 0; i < K; i++)
      q[i] = h[0] - EPS * (i - K + 1.);
   for (auto i = 0; i < s - K; i++)
      q[i + K] = (h[i] + h[(i + K)]) / 2.;
   for (auto i = s; i < s + K; i++)
      q[i] = h[s - 1] + EPS * (i - s);
   return q;
};
std::vector<std::vector<double>> Bspline_knot(const std::vector<std::vector<double>> &h, const int K) {
   std::vector<std::vector<double>> q(h.size());
   for (auto i = 0; i < h.size(); i++)
      q[i] = Bspline_knot(h[i], K);
   return q;
};
std::vector<double> Bspline_vector(const double h, const int s, const int K) {
   //  --------------  s --------------->
   // {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)}
   std::vector<double> ret(s, 0.), q = OpenUniformKnots(s, K);
   for (auto j = 0; j < s; j++)
      ret[j] = Bspline(h, q, j, K);
   return ret;
};  // parametric -1 to 1
std::vector<double> Bspline_vector(const double h, const int s, const int K, const std::vector<double> &q) {
   std::vector<double> ret(s, 0.);
   for (auto j = 0; j < s; j++)
      ret[j] = Bspline(h, q, j, K);
   return ret;
};
std::vector<double> D_Bspline_vector(const double h, const int s, const int K, int n) {
   //  --------------  s --------------->
   // {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)}
   std::vector<double> ret(s, 0.), q = OpenUniformKnots(s, K);
   for (auto j = 0; j < s; j++)
      ret[j] = D_Bspline(h, q, j, K, n);
   return ret;
};
std::vector<double> D_Bspline_vector(const double h, const int s, const int K, int n, const std::vector<double> &q) {
   std::vector<double> ret(s, 0.);
   for (auto j = 0; j < s; j++)
      ret[j] = D_Bspline(h, q, j, K, n);
   return ret;
};
//           --------------  s  --------------->
// |        {
// |        {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)},
// |        {B0(h1),B1(h1),B2(h1),B3(h1),B4(h1)},
// |        {B0(h2),B1(h2),B2(h2),B3(h2),B4(h2)},
// h.size() {B0(h3),B1(h3),B2(h3),B3(h3),B4(h3)},
// |        {B0(h4),B1(h4),B2(h4),B3(h4),B4(h4)},
// |        {B0(h5),B1(h5),B2(h5),B3(h5),B4(h5)},
// |        {B0(h6),B1(h6),B2(h6),B3(h6),B4(h6)},
// V        }
std::vector<std::vector<double>> Bspline_matrix(const std::vector<double> &h, const int s, const int K) {
   std::vector<std::vector<double>> ret(h.size(), std::vector<double>(s, 0.));
   for (auto i = 0; i < h.size(); i++)
      ret[i] = Bspline_vector(h[i], s, K);
   return ret;
};
std::vector<std::vector<double>> Bspline_matrix(const std::vector<double> &h, const int s, const int K, const std::vector<double> &q) {
   std::vector<std::vector<double>> ret(h.size(), std::vector<double>(s, 0.));
   for (auto i = 0; i < h.size(); i++)
      ret[i] = Bspline_vector(h[i], s, K, q);
   return ret;
};
std::vector<std::vector<double>> Bspline_matrix(const std::vector<std::vector<double>> &h, const int s, const int K, const std::vector<std::vector<double>> &q) {
   std::vector<std::vector<double>> ret(h.size(), std::vector<double>(s, 0.));
   for (auto i = 0; i < h.size(); i++)
      for (auto j = 0; j < s; j++)
         ret[i][j] = Bspline(h[i][j], q[i], s, K);
   return ret;
};
std::vector<std::vector<double>> Bspline_matrix(const std::vector<std::vector<double>> &h, const int s, const int K, const std::vector<double> &q) {
   std::vector<std::vector<double>> ret(h.size(), std::vector<double>(s, 0.));
   for (auto i = 0; i < h.size(); i++)
      for (auto j = 0; j < s; j++)
         ret[i][j] = Bspline(h[i][j], q, i, K);
   return ret;
};
std::vector<std::vector<double>> D_Bspline_matrix(const std::vector<double> &h, const int s, const int K, const int n) {
   std::vector<std::vector<double>> ret(h.size(), std::vector<double>(s, 0.));
   for (auto i = 0; i < h.size(); i++)
      ret[i] = D_Bspline_vector(h[i], s, K, n);
   return ret;
};
std::vector<std::vector<double>> D_Bspline_matrix(const std::vector<double> &h, const int s, const int K, const int n, const std::vector<double> &q) {
   std::vector<std::vector<double>> ret(h.size(), std::vector<double>(s, 0.));
   for (auto i = 0; i < h.size(); i++)
      ret[i] = D_Bspline_vector(h[i], s, K, n, q);
   return ret;
};
//    ---------------  s --------------->
// | {
// | {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)},
// | {B0(h1),B1(h1),B2(h1),B3(h1),B4(h1)},
// s {B0(h2),B1(h2),B2(h2),B3(h2),B4(h2)},
// | {B0(h3),B1(h3),B2(h3),B3(h3),B4(h3)},
// | {B0(h4),B1(h4),B2(h4),B3(h4),B4(h4)}
// V }
std::vector<std::vector<double>> Bspline_matrix(const int s, const int K) {
   return Bspline_matrix(Subdivide(-1., 1., s - 1) /*h*/, s, K);
};
std::vector<std::vector<double>> D_Bspline_matrix(const int s, const int K, const int n) {
   return D_Bspline_matrix(Subdivide(-1., 1., s - 1) /*h*/, s, K, n);
};
///////////////////////////////////////////////////////////
struct ParametricInterpolation {
  public:
   std::vector<std::vector<double>> samp3;
   std::vector<double> samp2;
   std::vector<double> q, xi, eta;
   std::vector<double> q1, q2;
   std::vector<std::vector<double>> InvMatB, InvMatBT;
   int s, K;
   int s1, s2;
   ///////////// 1 parm  ///////////
   ParametricInterpolation(const std::vector<double> &samp2_IN, const int K_IN) : samp2(samp2_IN),
                                                                                  s(samp2_IN.size()),
                                                                                  xi(Subdivide(double(-1), double(1), samp2_IN.size() - 1)),
                                                                                  q(OpenUniformKnots(Subdivide(double(-1), double(1), samp2_IN.size() - 1), K_IN)),
                                                                                  K(K_IN),
                                                                                  InvMatBT(Inverse(Transpose(Bspline_matrix(samp2_IN.size(), K_IN)))),
                                                                                  InvMatB(Inverse(Bspline_matrix(samp2_IN.size(), K_IN))){};
   std::vector<double> N(const double a) { return Dot(InvMatBT, Bspline_vector(a, s, K)); };
   double operator()(const double a) { return Dot(N(a), samp2); };
   //===========================
   std::vector<double> N(const double a, const int n) { return Dot(InvMatBT, D_Bspline_vector(a, s, K, n)); };
   double operator()(const double a, const int n) { return Dot(Dot(InvMatB, D_Bspline_vector(a, s, K, n)), samp2); };
   //========= specific ========
   double Basis(const double h, const int i) { return Bspline(h, q, i, K); };
   std::vector<double> Basis(const std::vector<double> &h, const int i) { return Bspline(h, q, i, K); };
   double Shape(const double h, const int i) { return Dot(InvMatBT[i], Bspline_vector(h, s, K)); };
   double DShape(const double h, const int i, const int n) { return Dot(InvMatBT[i], D_Bspline_vector(h, s, K, n)); };
   /////////// 2 parm  ////////////
   /*
   s1 is the size of a column vector : points in y direction
   s2 is the size of a row vector : points in x direction
 */
   ParametricInterpolation(const std::vector<std::vector<double>> &samp3_IN, const int K_IN) : samp3(samp3_IN),
                                                                                               s1(samp3_IN.size()) /*y*/,
                                                                                               s2(samp3_IN[0].size()) /*x*/,
                                                                                               q1(OpenUniformKnots(samp3_IN.size(), K_IN)),
                                                                                               q2(OpenUniformKnots(samp3_IN[0].size(), K_IN)),
                                                                                               K(K_IN),
                                                                                               InvMatB(Inverse(Bspline_matrix(samp3_IN.size() /*y*/, K_IN))),
                                                                                               InvMatBT(Inverse(Transpose(Bspline_matrix(samp3_IN[0].size() /*x*/, K_IN)))){};
   ParametricInterpolation(){};
   void reset(const std::vector<std::vector<double>> &samp3_IN, const int K_IN) {
      q1.clear();
      q2.clear();
      samp3.clear();
      InvMatB.clear();
      InvMatBT.clear();
      // insert
      samp3 = samp3_IN;
      s1 = samp3_IN.size() /*y*/;
      s2 = samp3_IN[0].size() /*x*/;
      q1 = OpenUniformKnots(samp3_IN.size(), K_IN);
      q2 = OpenUniformKnots(samp3_IN[0].size(), K_IN);
      K = K_IN;
      InvMatB = Inverse(Bspline_matrix(samp3_IN.size(), K_IN));
      InvMatBT = Inverse(Transpose(Bspline_matrix(samp3_IN[0].size(), K_IN)));
   };
   ////////////////////////////////////////////////////
   std::vector<double> Nx(const double b /*y*/) {
      return Dot(InvMatBT, Bspline_vector(b, s2, K, q2));
   };
   std::vector<double> Ny(const double a /*x*/) {
      return Dot(Bspline_vector(a, s1, K, q1), InvMatB);
   };
   std::vector<std::vector<double>> N(const std::vector<double> &ba) {
      return TensorProduct(Ny(ba[1]), Nx(ba[0]));
   };
   double operator()(const double b /*x*/, const double a /*y*/) {
      return Dot(Ny(a), Dot(samp3, Nx(b)));
   };
   double operator()(const std::vector<double> &ba) {
      return Dot(Ny(ba[1]), Dot(samp3, Nx(ba[0])));
   };
   ////////////////// DERIVATIVES /////////////////
   std::vector<double> DNx(const double b /*y*/, const int n) {
      return Dot(InvMatBT, D_Bspline_vector(b, s2, K, n, q2));
   };
   std::vector<double> DNy(const double a /*x*/, const int n) {
      return Dot(D_Bspline_vector(a, s1, K, n, q1), InvMatB);
   };
   std::vector<std::vector<double>> DN(const std::vector<double> &ba, const std::vector<int> &n) {
      return TensorProduct(DNy(ba[1], n[1]), DNx(ba[0], n[0]));
   };
   double operator()(const double b /*x*/, const double a /*y*/, const std::vector<int> &n) {
      return Dot(DNy(a, n[1]), Dot(samp3, DNx(b, n[0])));
   };
   double operator()(const std::vector<double> &ba, const std::vector<int> &n) {
      return Dot(DNy(ba[1], n[1]), Dot(samp3, DNx(ba[0], n[0])));
   };
   ///////////////// specific shape function ////////////
   double Nx(const double b /*y*/, const int j) {
      return Dot(InvMatBT[j], Bspline_vector(b, s2, K, q2));
   };
   double Ny(const double a /*x*/, const int m) {
      std::vector<double> By = Bspline_vector(a, s1, K, q1);
      double ret(0.);
      for (auto l = 0; l < s1; l++)
         ret += InvMatB[l][m] * By[l];
      return ret;
   };
   double N(const std::vector<double> &ba, const std::vector<int> &jm) { /*2D vector access*/
      return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
   };
   double N(const std::vector<double> &ba, const int mj) { /*1D access*/
      return Ny(ba[1], (int)(mj / s2) /*row*/) * Nx(ba[0], mj % s2 /*col*/);
   };
   // double Basis(const std::vector<double>& ba, const std::vector<int> &jm){
   //   return samp3[jm[1]][jm[0]] * Nx(ba[0], jm[0]) * Ny(ba[1], jm[1]);
   // };
   double Basis(const std::vector<double> &ba, const std::vector<int> &jm) {
      return Bspline(ba[0], q2, jm[0], K) * Bspline(ba[1], q1, jm[1], K);
   };
   double Shape(const std::vector<double> &ba, const std::vector<int> &jm) {
      return Nx(ba[0], jm[0]) * Ny(ba[1], jm[1]);
   };
};
///////////////////////////////////////////////////////////
class node_complex {
   //    *
   //    |
   // *--*--*
   //    |
   //    *
  public:
   std::vector<std::vector<double>> v_complex;
   std::vector<node_complex *> nodes;
   std::vector<bool> hits;
   node_complex(const std::vector<std::vector<double>> &v_complex_IN) : v_complex(v_complex_IN){};

   void connect(node_complex &node_IN) {
      node_complex *p = std::addressof(node_IN);
      for (auto &n : nodes)
         if (p == n)
            return;
      nodes.push_back(p);
   };
};

///////////////////////////////////////////////////////////
class glLINES {
  public:
   std::vector<std::vector<double>> v_complex;
   std::vector<std::vector<float>> s_complex;
   std::vector<std::vector<int>> f_v_complex;  // これのstd::vector<int>が１面 -> 3頂点を指定
   std::vector<bool> hits;

   std::vector<std::vector<double>> samp3X, samp3Y, samp3Z;
   // sample点がない場合
   glLINES(int row = 20, int col = 20, double scale = 5.) {
      /* sampleを作成 */
      std::vector<std::vector<double>> samp3X(row, std::vector<double>(col, 0));
      std::vector<std::vector<double>> samp3Y(row, std::vector<double>(col, 0));
      std::vector<std::vector<double>> samp3Z(row, std::vector<double>(col, 0));
      for (auto i = 0; i < row; i++) {
         for (auto j = 0; j < col; j++) {
            double x = (-1. + 2. / (col - 1) * j);
            double y = (-1. + 2. / (row - 1) * i);
            samp3X[i][j] = x * scale;
            samp3Y[i][j] = y * scale;
            samp3Z[i][j] = std::sin(2. * M_PI * x) * std::sin(2. * M_PI * y);
         }
      }
      /* gl.LINESでメッシュがかけるような形式で格納 */
      set(samp3X, samp3Y, samp3Z);
   };
   // sample点がある場合
   glLINES(const std::vector<std::vector<double>> &samp3X_IN,
           const std::vector<std::vector<double>> &samp3Y_IN,
           const std::vector<std::vector<double>> &samp3Z_IN) : samp3X(samp3X_IN), samp3Y(samp3Y_IN), samp3Z(samp3Z_IN) {
      /* gl.LINESでメッシュがかけるような形式で格納 */
      set(samp3X_IN, samp3Y_IN, samp3Z_IN);
   };
   // verticesの格納
   void set(const std::vector<std::vector<double>> &samp3X_IN,
            const std::vector<std::vector<double>> &samp3Y_IN,
            const std::vector<std::vector<double>> &samp3Z_IN) {
      auto col = samp3X_IN[0].size();
      auto row = samp3X_IN.size();
      v_complex.clear();
      s_complex.clear();
      f_v_complex.clear();
      /* gl.LINESでメッシュがかけるような形式で格納 */
      int l = 0;
      //    8 points
      //    0   1
      //   7*---*2
      //    |   |
      //    |   |
      //   6*---*3
      //    5   4
      for (auto i = 0; i < row - 1; i++)
         for (auto j = 0; j < col - 1; j++) {
            v_complex.push_back({samp3X_IN[i][j], samp3Z_IN[i][j], samp3Y_IN[i][j]});              // 0
            v_complex.push_back({samp3X_IN[i][j + 1], samp3Z_IN[i][j + 1], samp3Y_IN[i][j + 1]});  // 1
            f_v_complex.push_back({l++, l++});
            v_complex.push_back({samp3X_IN[i][j + 1], samp3Z_IN[i][j + 1], samp3Y_IN[i][j + 1]});              // 2
            v_complex.push_back({samp3X_IN[i + 1][j + 1], samp3Z_IN[i + 1][j + 1], samp3Y_IN[i + 1][j + 1]});  // 3
            f_v_complex.push_back({l++, l++});
            v_complex.push_back({samp3X_IN[i + 1][j + 1], samp3Z_IN[i + 1][j + 1], samp3Y_IN[i + 1][j + 1]});  // 4
            v_complex.push_back({samp3X_IN[i + 1][j], samp3Z_IN[i + 1][j], samp3Y_IN[i + 1][j]});              // 5
            f_v_complex.push_back({l++, l++});
            v_complex.push_back({samp3X_IN[i + 1][j], samp3Z_IN[i + 1][j], samp3Y_IN[i + 1][j]});  // 6
            v_complex.push_back({samp3X_IN[i][j], samp3Z_IN[i][j], samp3Y_IN[i][j]});              // 7
            f_v_complex.push_back({l++, l++});
         }

      s_complex.resize(v_complex.size(), std::vector<float>(4, 0));  // RGBS
      for (auto &tmp : s_complex)
         tmp = {.9, .0, .0, 1.};

      this->hits.resize(f_v_complex.size());
   };
};
///////////////////////////////////////////////////////////
struct ParametricInterpolation3D {
  public:
   std::vector<std::vector<double>> v_complex;
   std::vector<std::vector<float>> s_complex;
   std::vector<std::vector<int>> f_v_complex;  // これのstd::vector<int>が１面 -> 3頂点を指定

   std::vector<std::vector<double>> samp3X, samp3Y, samp3Z;
   std::vector<double> q1, q2;
   std::vector<std::vector<double>> inv_coefmat, inv_coefmatT;
   int s1, s2, K;
   ////////// 2 parm  //////////
   /*
   s1 is the size of a column vector : points in y direction
   s2 is the size of a row vector : points in x direction
 */
   ParametricInterpolation3D(const std::vector<std::vector<double>> &samp3X_IN,
                             const std::vector<std::vector<double>> &samp3Y_IN,
                             const std::vector<std::vector<double>> &samp3Z_IN,
                             const int K_IN) : K(K_IN),
                                               s1(samp3X_IN.size()) /*y*/,
                                               s2(samp3X_IN[0].size()) /*x*/,
                                               q1(OpenUniformKnots(samp3X_IN.size(), K_IN)),
                                               q2(OpenUniformKnots(samp3X_IN[0].size(), K_IN)),
                                               inv_coefmat(Inverse(Bspline_matrix(samp3X_IN.size() /*y*/, K_IN))),
                                               inv_coefmatT(Inverse(Transpose(Bspline_matrix(samp3X_IN[0].size() /*x*/, K_IN)))),
                                               samp3X(samp3X_IN),
                                               samp3Y(samp3Y_IN),
                                               samp3Z(samp3Z_IN) {
      if (samp3X_IN.size() != samp3Y_IN.size() || samp3X_IN.size() != samp3Z_IN.size()) {
         std::cerr << Red << "sizes of input samples are differen" << reset << std::endl;
         std::cerr << Red << samp3X_IN.size() << reset << std::endl;
         std::cerr << Red << samp3Y_IN.size() << reset << std::endl;
         std::cerr << Red << samp3Z_IN.size() << reset << std::endl;
         abort();
      }
   };
   ParametricInterpolation3D(){};
   void Reset(const std::vector<std::vector<double>> &samp3X_IN,
              const std::vector<std::vector<double>> &samp3Y_IN,
              const std::vector<std::vector<double>> &samp3Z_IN,
              const int K_IN) {
      if (samp3X_IN.size() != samp3Y_IN.size() || samp3X_IN.size() != samp3Z_IN.size()) {
         std::cerr << Red << "sizes of input samples are differen" << reset << std::endl;
         std::cerr << Red << samp3X_IN.size() << reset << std::endl;
         std::cerr << Red << samp3Y_IN.size() << reset << std::endl;
         std::cerr << Red << samp3Z_IN.size() << reset << std::endl;
         abort();
      }
      q1.clear();
      q2.clear();
      inv_coefmat.clear();
      inv_coefmatT.clear();
      samp3X.clear();
      samp3Y.clear();
      samp3Z.clear();
      // insert
      K = K_IN;
      s1 = samp3X_IN.size() /*y*/;
      s2 = samp3X_IN[0].size() /*x*/;
      q1 = OpenUniformKnots(samp3X_IN.size(), K_IN);
      q2 = OpenUniformKnots(samp3X_IN[0].size(), K_IN);
      inv_coefmat = Inverse(Bspline_matrix(samp3X_IN.size(), K_IN));
      inv_coefmatT = Inverse(Transpose(Bspline_matrix(samp3X_IN[0].size(), K_IN)));
      samp3X = samp3X_IN;
      samp3Y = samp3Y_IN;
      samp3Z = samp3Z_IN;
   };
   ////////////////////////////////////////////////////
   std::vector<double> Nx(const double b /*y*/) {
      return Dot(inv_coefmatT, Bspline_vector(b, s2, K, q2));
   };
   std::vector<double> Ny(const double a /*x*/) {
      return Dot(Bspline_vector(a, s1, K, q1), inv_coefmat);
   };
   std::vector<std::vector<double>> N(const std::vector<double> &ba) {
      return TensorProduct(Ny(ba[1]), Nx(ba[0]));
   };
   std::vector<double> operator()(const double b /*x*/, const double a /*y*/) {
      std::vector<double> shapey = Ny(a), shapex = Nx(b), ret(3, 0.);
      for (auto m = 0; m < s1; m++) {
         ret[0] += shapey[m] * std::inner_product(samp3X[m].begin(), samp3X[m].end(), shapex.begin(), 0.);
         ret[1] += shapey[m] * std::inner_product(samp3Y[m].begin(), samp3Y[m].end(), shapex.begin(), 0.);
         ret[2] += shapey[m] * std::inner_product(samp3Z[m].begin(), samp3Z[m].end(), shapex.begin(), 0.);
      }
      return ret;
   };
   std::vector<double> operator()(const std::vector<double> &ba) {
      std::vector<double> shapey = Ny(ba[1]), shapex = Nx(ba[0]), ret(3, 0.);
      for (auto m = 0; m < s1; m++) {
         ret[0] += shapey[m] * std::inner_product(samp3X[m].begin(), samp3X[m].end(), shapex.begin(), 0.);
         ret[1] += shapey[m] * std::inner_product(samp3Y[m].begin(), samp3Y[m].end(), shapex.begin(), 0.);
         ret[2] += shapey[m] * std::inner_product(samp3Z[m].begin(), samp3Z[m].end(), shapex.begin(), 0.);
      }
      return ret;
   };
   ////////////////// DERIVATIVES /////////////////
   std::vector<double> DNx(const double b /*y*/, const int n) {
      return Dot(inv_coefmatT, D_Bspline_vector(b, s2, K, n, q2));
   };
   std::vector<double> DNy(const double a /*x*/, const int n) {
      std::vector<double> By = D_Bspline_vector(a, s1, K, n, q1);
      std::vector<double> ret(s1, 0.);
      for (auto m = 0; m < s1; m++)
         for (auto l = 0; l < s1; l++)
            ret[m] += inv_coefmat[l][m] * By[l];
      return ret;
   };
   std::vector<std::vector<double>> DN(const std::vector<double> &ba, const std::vector<int> &n) {
      std::vector<double> shapey = DNy(ba[1], n[1]), shapex = DNx(ba[0], n[0]);
      std::vector<std::vector<double>> ret(s1, std::vector<double>(s2, 0.));
      for (auto m = 0; m < s1; m++)
         for (auto j = 0; j < s2; j++)
            ret[m][j] = shapey[m] * shapex[j];
      return ret;
   };
   std::vector<double> operator()(const double b /*x*/, const double a /*y*/, const std::vector<int> &n) {
      std::vector<double> shapey = DNy(a, n[1]), shapex = DNx(b, n[0]), ret(3, 0.);
      for (auto m = 0; m < s1; m++)
         for (auto j = 0; j < s2; j++) {
            ret[0] += samp3X[m][j] * shapey[m] * shapex[j];
            ret[1] += samp3Y[m][j] * shapey[m] * shapex[j];
            ret[2] += samp3Z[m][j] * shapey[m] * shapex[j];
         }
      return ret;
   };
   std::vector<double> operator()(const std::vector<double> &ba, const std::vector<int> &n) {
      std::vector<double> shapey = DNy(ba[1], n[1]), shapex = DNx(ba[0], n[0]), ret(3, 0.);
      for (auto m = 0; m < s1; m++)
         for (auto j = 0; j < s2; j++) {
            ret[0] += samp3X[m][j] * shapey[m] * shapex[j];
            ret[1] += samp3Y[m][j] * shapey[m] * shapex[j];
            ret[2] += samp3Z[m][j] * shapey[m] * shapex[j];
         }
      return ret;
   };
   ///////////////// specific shape function ////////////
   double Nx(const double b /*y*/, const int j) {
      return std::inner_product(inv_coefmatT[j].begin(), inv_coefmatT[j].end(), Bspline_vector(b, s2, K, q2).begin(), 0.);
   };
   double Ny(const double a /*x*/, const int m) {
      std::vector<double> By = Bspline_vector(a, s1, K, q1);
      double ret(0.);
      for (auto l = 0; l < s1; l++)
         ret += inv_coefmat[l][m] * By[l];
      return ret;
   };
   double N(const std::vector<double> &ba, const std::vector<int> &jm) { /*2D vector access*/
      return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
   };
   double N(const std::vector<double> &ba, const int mj) { /*1D access*/
      return Ny(ba[1], (int)(mj / s2) /*row*/) * Nx(ba[0], mj % s2 /*col*/);
   };
   std::vector<double> basis(const std::vector<double> &ba, const std::vector<int> &jm) {
      return std::vector<double>{samp3X[jm[1]][jm[0]], samp3Y[jm[1]][jm[0]], samp3Z[jm[1]][jm[0]]} * Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
   };
   std::vector<double> cross(const std::vector<double> &ba) {
      std::vector<double> Dshapey = DNy(ba[1], 1), Dshapex = DNx(ba[0], 1), shapey = DNy(ba[1], 0), shapex = DNx(ba[0], 0), ret_dxi(3, 0.), ret_deta(3, 0.);
      for (auto m = 0; m < s1; m++)
         for (auto j = 0; j < s2; j++) {
            ret_dxi[0] += samp3X[m][j] * shapey[m] * Dshapex[j];
            ret_dxi[1] += samp3Y[m][j] * shapey[m] * Dshapex[j];
            ret_dxi[2] += samp3Z[m][j] * shapey[m] * Dshapex[j];
            ret_deta[0] += samp3X[m][j] * Dshapey[m] * shapex[j];
            ret_deta[1] += samp3Y[m][j] * Dshapey[m] * shapex[j];
            ret_deta[2] += samp3Z[m][j] * Dshapey[m] * shapex[j];
         }
      return Cross(ret_dxi, ret_deta);
   };
};
////////////////////////////////////////////////////////
struct RBF_fn {
   virtual double rbf(double r) = 0;
};
struct RBF_interp {
   int dim, n;
   std::vector<std::vector<double>> pts;
   std::vector<double> vals;
   std::vector<double> w;
   RBF_fn &fn;
   bool norm;

   RBF_interp(const std::vector<std::vector<double>> &ptss /* x for f(x) */,
              const std::vector<double> &valss /*f(x)*/,
              RBF_fn &func,
              bool nrbf = false)
       : dim(ptss[0].size()), n(ptss.size()), pts(ptss), vals(valss), w(n), fn(func), norm(nrbf) {
      double sum;
      std::vector<std::vector<double>> rbf(n, std::vector<double>(n, 0.));
      std::vector<double> rhs(n);
      for (auto i = 0; i < n; i++) {
         sum = 0.;
         for (auto j = 0; j < n; j++)
            sum += (rbf[i][j] = fn.rbf(rad(&pts[i][0], &pts[j][0])));
         rhs[i] = norm ? sum * vals[i] : vals[i];
      }
      ludcmp lu(rbf);
      lu.solve(rhs, w);
   }
   double operator()(const std::vector<double> &pt) {
      double fval, sum = 0., sumw = 0.;
      if (pt.size() != dim)
         throw("RBF_interp bad pt size");
      for (auto i = 0; i < n; i++) {
         fval = fn.rbf(rad(&pt[0], &pts[i][0]));
         sumw += w[i] * fval;
         sum += fval;
      }
      return norm ? sumw / sum : sumw;
   }

   double rad(const double *p1, const double *p2) {
      double sum = 0.;
      for (int i = 0; i < dim; i++)
         sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
      return std::sqrt(sum);
   }
};
struct RBF_multiquadric : RBF_fn {
   double r02;
   RBF_multiquadric(double scale = 1.) : r02(scale * scale) {}
   double rbf(double r) { return std::sqrt(r * r + r02); }
};
struct RBF_thinplate : RBF_fn {
   double r0;
   RBF_thinplate(double scale = 1.) : r0(scale) {}
   double rbf(double r) { return r <= 0. ? 0. : r * r * std::log(r / r0); }
};
struct RBF_gauss : RBF_fn {
   double r0;
   RBF_gauss(double scale = 1.) : r0(scale) {}
   double rbf(double r) { return exp(-0.5 * r * r / (r0 * r0)); }
};
struct RBF_inversemultiquadric : RBF_fn {
   double r02;
   RBF_inversemultiquadric(double scale = 1.) : r02(scale * scale) {}
   double rbf(double r) { return 1. / std::sqrt(r * r + r02); }
};
struct Shep_interp {
   int dim, n;
   const std::vector<std::vector<double>> &pts;
   const std::vector<double> &vals;
   double pneg;

   Shep_interp(std::vector<std::vector<double>> &ptss, std::vector<double> &valss, double p = 2.)
       : dim(ptss[0].size()), n(ptss.size()), pts(ptss), vals(valss), pneg(-p) {}

   double interp(std::vector<double> &pt) {
      double r, w, sum = 0., sumw = 0.;
      if (pt.size() != dim)
         throw("RBF_interp bad pt size");
      for (auto i = 0; i < n; i++) {
         if ((r = rad(&pt[0], &pts[i][0])) == 0.)
            return vals[i];
         sum += (w = pow(r, pneg));
         sumw += w * vals[i];
      }
      return sumw / sum;
   }

   double rad(const double *p1, const double *p2) {
      double sum = 0.;
      for (auto i = 0; i < dim; i++)
         sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
      return std::sqrt(sum);
   }
};

// #ifdef INCL_lapacke
// struct LU_LAPACK{
// public:
//   char TRANS = 'N';
//   int INFO;
//   int LDA;
//   int LDB;
//   int N;
//   int NRHS = 1;
//   std::vector<int> IPIV;
//   std::vector<double> Flatten_A;

//   LU_LAPACK(const std::vector<std::vector<double>>& A/* must be square matrix */):
//     Flatten_A(Flatten(A)),
//     INFO(A.size()),
//     LDA(A.size()),
//     LDB(A.size()),
//     N(A.size()),
//     IPIV(A.size())
//   {
//     std::cout << "compute the LU factorization..." << std::endl;
//     //void LAPACK_dgetrf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda, lapack_int* ipiv, lapack_int *info );
//     LAPACK_dgetrf(&N, &N, &Flatten_A[0], &LDA, &IPIV[0], &INFO);
//     if(INFO)
//       std::cout << "an error occured : "<< INFO << std::endl;
//   };
//   /////////////////////////////////////////////////////////////
//   std::vector<double> solve(const std::vector<double>& B){
//     std::vector<double> ret(B);
//     // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
//     std::cout << "solving the system..."<< std::endl;
//     // void LAPACK_dgetrs( char* trans, lapack_int* n, lapack_int* nrhs, const double* a, lapack_int* lda, const lapack_int* ipiv,double* b, lapack_int* ldb, lapack_int *info );
//     dgetrs_(&TRANS, &N, &NRHS, &Flatten_A[0], &LDA, &IPIV[0], &ret[0], &LDB,&INFO);
//     if(INFO)
//       {
// 	std::cout << "an error occured : "<< INFO << std::endl;
// 	abort();
//       }
//     return ret;
//   };
//   /////////////////////////////////////////////////////////////
//   void solve(const std::vector<double>& B, std::vector<double>& ret){
//     // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
//     std::cout << "solving the system..."<< std::endl;
//     // void LAPACK_dgetrs( char* trans, lapack_int* n, lapack_int* nrhs, const double* a, lapack_int* lda, const lapack_int* ipiv,double* b, lapack_int* ldb, lapack_int *info );
//     dgetrs_(&TRANS, &N, &NRHS, &Flatten_A[0], &LDA,&IPIV[0], &ret[0], &LDB,&INFO);
//     if(INFO)
//       {
// 	std::cout << "an error occured : "<< INFO << std::endl;
// 	abort();
//       }
//   };
// };
// #endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//=========================================================
// void Load(const std::string& filename, const std::vector<std::string>& sep, std::vector<std::vector<std::vector<double>>>& mat){
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     std::cout << Red << "The file can not be opened" << reset << std::endl;
//   std::string str, S, cell;
//   double v;
//   std::size_t found, s(sep.size()), min_len, left, right;
//   bool loop;
//   int row_counter(0), col_counter(0);
//   while(!strm.eof( ) && mat.size() > row_counter)
//     {
//       std::getline(strm,str);
//       if(!str.empty())
// 	{
// 	  if((left = str.find("{"))!=std::string::npos);
// 	  if((right = str.find("}"))!=std::string::npos);
// 	  cell = str.substr(left,right);
// 	  str = str.substr(right+1);

// 	  std::vector<double> row_vector;
// 	  loop = true;
// 	  col_counter = 0;

// 	  while(loop && mat[row_counter].size() > col_counter)
// 	    {
// 	      min_len = 999999;
// 	      for(auto i=0; i<s; i++)
// 		if((found = cell.find(sep[i]))!=std::string::npos)
// 		  if(min_len > found)
// 		    {
// 		      min_len = found;
// 		      S = sep[i];
// 		    }
// 	      col_counter++;
// 	    }
// 	  row[row_counter].push_back(row_vector);
// 	  row_counter++;
// 	}
//     }
//   strm.close();
// };
//=========================================================
template <class Ttype>
void load(std::ifstream &in, std::vector<Ttype> &vec, const int row_size) {
   vec.clear();
   if (!in) {
      std::cout << "入力ファイルが開けません" << std::endl;
   }
   std::string str;
   int i;
   Ttype tmp;
   i = 0;
   getline(in, str);
   std::cout << std::string(40, '=') << std::endl;
   std::cout << "この行は，" << str << std::endl;
   std::cout << std::string(40, '-') << std::endl;
   std::cout << "読み込み結果：" << std::endl;
   for (size_t j = 0; j < row_size; j++) {

      if (str.find(',') != std::string::npos) {
         i = str.find(',');
         //	std::cout << "\',\'は" << i << "番目にあります" << endl;
         //	std::cout << str << "　を" << i << "まで読み込みます -> " << str.substr(0,i) << " -> これを数値に変換すると ->";
         tmp = (Ttype)stod(str.substr(0, i));
         //	std::cout << std::setprecision(10) << tmp << endl;
         str = str.substr(i + 1);
         //	std::cout << i+1 << "番目以降を抜き出し文字列を置き換えます -> " << str << endl;
      } else {
         //	std::cout << "\',\'は" << "見つかりませんでした -> ";
         tmp = (Ttype)stod(str);
         //	std::cout << str << "　を" << "最後まで読み込みます -> " << std::setprecision(10) << tmp << endl;
      }
      vec.push_back(tmp);
      //      std::cout << tmp << "を代入 -> ";
      std::cout << "vec[" << j << "] = " << std::setprecision(5) << vec[j] << std::endl;
   }
   std::cout << std::string(40, '-') << std::endl;
};
//===========================================================
template <class Ttype>
void load(std::ifstream &in, std::vector<std::vector<Ttype>> &vec, const int row_size) {
   vec.clear();
   vec.resize(row_size, std::vector<Ttype>(0));
   if (!in) {
      std::cout << "入力ファイルが開けません" << std::endl;
   }
   std::string str;
   int i;
   Ttype tmp;
   for (int line = 0;; line++) {
      i = 0;
      getline(in, str);
      if (in.eof())
         break;
      std::cout << std::string(40, '=') << std::endl;
      std::cout << "この行は，" << str << std::endl;
      std::cout << std::string(40, '-') << std::endl;
      std::cout << "読み込み結果：" << std::endl;
      for (int j = 0; j < row_size; j++) {

         if (str.find(',') != std::string::npos) {
            i = str.find(',');
            //	std::cout << "\',\'は" << i << "番目にあります" << endl;
            //	std::cout << str << "　を" << i << "まで読み込みます -> " << str.substr(0,i) << " -> これを数値に変換すると ->";
            tmp = (Ttype)stod(str.substr(0, i));
            //	std::cout << std::setprecision(10) << tmp << endl;
            str = str.substr(i + 1);
            //	std::cout << i+1 << "番目以降を抜き出し文字列を置き換えます -> " << str << endl;
         } else {
            //	std::cout << "\',\'は" << "見つかりませんでした -> ";
            tmp = (Ttype)stod(str);
            //	std::cout << str << "　を" << "最後まで読み込みます -> " << std::setprecision(10) << tmp << endl;
         }
         vec[j].push_back(tmp);
         //      std::cout << tmp << "を代入 -> ";
         std::cout << "vec[" << j << "][" << line << "] = " << std::setprecision(5) << vec[j][line] << std::endl;
      }
      std::cout << std::string(40, '-') << std::endl;
   };
};
//=====================================================================
//=================== Mathematica output loader =======================
//=====================================================================
std::vector<std::string> StringSplit(const std::string &strIN, const std::vector<std::string> &SEP) {
   std::string str(strIN), foundFirstSep; /* expexting "1,2,3,4,5" or "1, 2,3, 4,5"*/
   std::vector<std::string> ret(0), tmp(0);
   while (!str.empty()) {
      size_t foundFirst(999999);
      bool bfound(false);
      for (const auto &sep : SEP) {
         size_t found = str.find(sep);
         if (!sep.empty()                                       /* ignore sep, "" */
             && found != std::string::npos                      /* if sep is found */
             && (found < foundFirst)) /* choose closest one */  // do not use "<=" instead
         {
            foundFirst = found;
            foundFirstSep = sep;
            bfound = true;
         }
      }
      if (bfound) {
         tmp.push_back(str.substr(0, foundFirst));
         str = str.substr(foundFirst + foundFirstSep.length());
      } else {
         tmp.push_back(str);

         // remove first and last cells if they are empty strings
         while (tmp.size() > 0)
            if ((*tmp.begin()).empty())
               tmp.erase(tmp.begin());
            else if ((*tmp.rbegin()).empty() ||
                     (*tmp.rbegin()).find_first_not_of(" ") == std::string::npos /* there is no string except spaces */ ||
                     (*tmp.rbegin()).find_first_not_of("\r") == std::string::npos /* there is no string except \r */ ||
                     (*tmp.rbegin()).find_first_not_of("\n") == std::string::npos /* there is no string except \n */ ||
                     (*tmp.rbegin()).find_first_not_of("\t") == std::string::npos /* there is no string except \t */)
               tmp.pop_back();
            else
               break;

         return tmp;
      }
   }

   // remove first and last cells if they are empty strings

   while (tmp.size() > 0)
      if ((*tmp.begin()).empty())
         tmp.erase(tmp.begin());
      else if ((*tmp.rbegin()).empty() ||
               (*tmp.rbegin()).find_first_not_of(" ") == std::string::npos /* there is no string except spaces */ ||
               (*tmp.rbegin()).find_first_not_of("\r") == std::string::npos /* there is no string except \r */ ||
               (*tmp.rbegin()).find_first_not_of("\n") == std::string::npos /* there is no string except \n */ ||
               (*tmp.rbegin()).find_first_not_of("\t") == std::string::npos /* there is no string except \t */)
         tmp.pop_back();
      else
         break;

   return tmp;
};
std::string StringTrim(const std::string &strIN, const std::vector<std::string> &SEP) {
   std::vector<std::string> str = StringSplit(strIN, SEP);
   std::string ret;
   for (const auto &s : str)
      ret += s;
   return ret;
};
std::vector<std::string> StringTrim(const std::vector<std::string> &strIN, const std::vector<std::string> &SEP) {
   std::vector<std::string> strOUT(strIN.size());
   for (auto i = 0; i < strOUT.size(); i++)
      strOUT[i] = StringTrim(strIN[i], SEP);
   return strOUT;
};
std::string StringJoin(const std::vector<std::string> &strIN, const std::string sep = "") {
   std::string ret(*strIN.begin());
   for (auto i = 1; i < strIN.size(); i++)
      ret += sep + strIN[i];
   return ret;
};
std::string Directory(const std::string &f) {
   auto dir = StringSplit(f, {"/"});
   dir.pop_back();
   return "/" + StringJoin(dir, "/") + "/";
};
//============================================================
//========================= Load  ============================
//============================================================
// std::vector<std::vector<double>> Load(const std::string& filename, const std::vector<std::string>& sep){
//   std::vector<std::vector<double>> mat;
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     std::cout << Red << filename << " can not be opened" << reset << std::endl;
//   else
//     std::cout << Blue << filename << " is opened" << reset << std::endl;

//   std::string str, S;
//   double v;
//   std::size_t found, s(sep.size()), min_len;
//   bool loop;
//   while(!strm.eof( ))
//     {
//       std::getline(strm,str);
//       if(!str.empty())
// 	{
// 	  std::vector<double> row_vector;
// 	  loop = true;
// 	  while(loop)
// 	    {
// 	      min_len = 999999;
// 	      for(auto i=0; i<s; i++)
// 		if((found = str.find(sep[i]))!=std::string::npos)
// 		  if(min_len > found)
// 		    {
// 		      min_len = found;
// 		      S = sep[i];
// 		    }

// 	      if(min_len!=999999/*if not 0*/)
// 		{
// 		  row_vector.push_back((double)stod(str.substr(0,min_len)));
// 		  str = str.substr(min_len + S.length());
// 		}
// 	      else
// 		{
// 		  row_vector.push_back((double)stod(str));
// 		  loop = false;
// 		}
// 	    }
// 	  mat.push_back(row_vector);
// 	}
//     }
//   strm.close();
//   return mat;
// };
//////////////////////////////
void Load(const std::string &filename, std::vector<std::vector<std::string>> &ret_mat, const std::vector<std::string> &SEP) {
   ret_mat.clear();
   std::ifstream strm(filename, std::ios::in);
   if (!strm)
      std::cout << Red << filename << " can not be opened" << reset << std::endl;
   else
      std::cout << Blue << filename << " is opened" << reset << std::endl;

   int row_counter(0);
   while (!strm.eof())  // loop until the end of file
   {
      std::string read;
      std::getline(strm, read);

      if (!StringTrim(read, {" "}).empty()) {
         ret_mat.push_back(StringTrim(StringSplit(read, SEP), {" ", "\t", "\n", "\r"}));
      }
   }
   strm.close();
   std::cout << Blue << filename << " is closed" << reset << std::endl;
};

std::vector<std::vector<std::string>> Load(const std::string &filename, const std::vector<std::string> &SEP) {
   std::vector<std::vector<std::string>> ret;
   Load(filename, ret, SEP);
   return ret;
};
// //////////////////////////////
// void Load(const std::string& filename, std::vector<std::vector<std::vector<double>>>& vvv)
// {
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     {
//       std::cout << Red << "The file can not be opened: " << filename << reset << std::endl;
//       std::cin.ignore();
//       abort();
//     }
//   std::string str;
//   while(!strm.eof( ))
//     {
//       std::getline(strm,str);

//       // std::cout <<  str << std::endl;
//       // str = str.substr(str.find_last_of("{")+1,  str.find_last_of("}")-1);
//       // std::cout <<  str << std::endl;

//       std::vector<std::string> str_vec = to_cell(str);/* expecting {"{1,2,3,5,5}","{1,2,3,5,5}","{1,2,3,5,5}"} */
//       std::vector<std::vector<double>> row_vec(0);
//       for(const auto& v: str_vec)
// 	{
// 	  std::string tmp = v.substr(v.find_first_of("{")+1,  v.find_last_of("}")-1);
// 	  tmp.find("");
// 	  std::vector<double> vec_doub = string_to_vector_double(/* expecting {1,2,3,5,5} */tmp,{",",", "});
// 	  row_vec.push_back(vec_doub);
// 	}
//       vvv.push_back(row_vec);
//     }
//   strm.close();
// };
// /////////////////////////////////
// void Load(const std::string& filename, std::vector<std::vector<double>>& vv)
// {
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     {
//       std::cout << Red << "The file can not be opened: " << filename << reset << std::endl;
//       std::cin.ignore();
//       abort();
//     }
//   std::string str;
//   while(!strm.eof( ))
//     {
//       std::getline(strm,str);

//       // std::cout <<  str << std::endl;
//       // str = str.substr(str.find_last_of("{")+1,  str.find_last_of("}")-1);
//       // std::cout <<  str << std::endl;

//       std::vector<double> vec_doub = string_to_vector_double(/* expecting {1,2,3,5,5} */str,{",",", "});
//       vv.push_back(vec_doub);
//     }
//   strm.close();
// };
//============================================================
//============================================================
template <class T>
std::vector<T> PointsToSurface(T t0, T t1, const std::vector<std::vector<T>> &a) {
   T EPS = 1E-15;

   if (a.size() == 3) {
      T t2 = 1. - t1 - t0;
      T x0 = 1. - t0 - 2. * t1, x1 = 1. - t1 - 2. * t2, x2 = 1. - t2 - 2. * t0;

      if (t0 < EPS) {
         return a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
      } else if (t1 < EPS) {
         return a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
      } else if (t2 < EPS) {
         return a[2] * (1. - x2) / 2. + a[0] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
      } else {
         std::vector<T> l0 = a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
         std::vector<T> l1 = a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
         std::vector<T> l2 = a[2] * (1. - x2) / 2. + a[0] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
         std::vector<T> tmp = {1 / t0, 1 / t1, 1 / t2};
         return Dot(tmp / (1 / t0 + 1 / t1 + 1 / t2), {l0, l1, l2});
      }
   } else {
      T t2 = 1. - t0, t3 = 1. - t1;
      T x0 = 2. * t3 - 1., x1 = 2. * t0 - 1., x2 = 2. * t1 - 1., x3 = 2. * t2 - 1.;

      if (t0 < EPS) {
         return a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
      } else if (t1 < EPS) {
         return a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
      } else if (t2 < EPS) {
         return a[2] * (1. - x2) / 2. + a[3] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
      } else if (t3 < EPS) {
         return a[3] * (1. - x3) / 2. + a[0] * (1. + x3) / 2. /*+p2A3*(x3^2-1)*/;
      } else {
         std::vector<T> l0 = a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
         std::vector<T> l1 = a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
         std::vector<T> l2 = a[2] * (1. - x2) / 2. + a[3] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
         std::vector<T> l3 = a[3] * (1. - x3) / 2. + a[0] * (1. + x3) / 2. /*+p2A3*(x3^2-1)*/;
         std::vector<T> tmp = {1 / t0, 1 / t1, 1 / t2, 1 / t3};
         return Dot(tmp / (1 / t0 + 1 / t1 + 1 / t2 + 1 / t3), {l0, l1, l2, l3});
      }
   }
};
//============================================================

std::vector<double> string_to_vector_double(const std::string &strIN, const std::vector<std::string> &SEP) {
   std::vector<std::string> str = StringSplit(strIN, SEP);
   std::vector<double> ret;
   for (auto i = 0; i < str.size(); i++)
      if (!str[i].empty())
         ret.push_back(stod(str[i]));

   return ret;
   // std::string str(strIN);/* expexting "1,2,3,4,5" or "1, 2,3, 4,5"*/
   // std::vector<double> ret(0);
   // size_t found;
   // int i(0);
   // while(!str.empty())
   //   {
   //     for(auto sep: SEP)
   //       if((found = str.find(sep))!=std::string::npos)
   // 	  {
   // 	    // std::cout << str << std::endl;
   // 	    // std::cout << found << std::endl;
   // 	    // std::cout << sep << std::endl;
   // 	    // std::cout << str << std::endl;
   // 	    ret.push_back(stod(
   // 			       str.substr((i==0 ? 0 : found+sep.length()), found)
   // 			       ));
   // 	    str = str.substr(found+sep.length());
   // 	    break;
   // 	  }
   // 	else
   // 	  {
   // 	    ret.push_back(stod(
   // 			       str.substr((i==0 ? 0 : found+sep.length()))
   // 			       ));
   // 	    return ret;
   // 	  }
   //     i++;
   //   }
   // return ret;
};
///////////////
std::vector<std::string> to_cell(const std::string &strIN) {
   std::string str(strIN); /* expexting "{1,2,3,4,5} {1,2,3,4,5} {1,2,3,4,5}"*/
   std::vector<std::string> ret(0);
   size_t found, right, left;
   int i(0);
   while (!str.empty()) {
      if ((left = str.find("{")) != std::string::npos) {
         if ((right = str.find("}")) != std::string::npos) {
            std::string cell = str.substr(left, right - left + 1);
            ret.push_back(cell);
            str = str.substr(right + 1);
            // std::cout << str << std::endl;
         } else
            return ret;
      } else
         return ret;
   }
   return ret;
};
///////////////////////////

template <class Type>
std::vector<Type> RotateLeft(const std::vector<Type> &vecs, const int n) {  // 2020/03/22
   std::vector<Type> ret(vecs);
   std::rotate(ret.begin(), ret.begin() + n, ret.end());
   return ret;
}

template <class Type>
std::vector<Type> Drop(const std::vector<Type> &vecs, std::vector<int> vec_n) {  // 2020/03/22
   std::vector<Type> ret(vecs);
   std::sort(vec_n.begin(), vec_n.end());
   for (auto i = 0; i < vec_n.size(); i++)
      ret.erase(ret.begin() + (vec_n[i] - i));
   return ret;
};

int Position(std::vector<int> vecIN, const int n) {  // 2020/03/22
   std::vector<int>::iterator it = std::find(vecIN.begin(), vecIN.end(), n);
   if (it == vecIN.end())
      std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;

   return std::distance(vecIN.begin(), it);
};

std::vector<int> Union(const std::vector<std::vector<int>> &vecs) {  // Union without sorting!
   std::vector<int> ret(0);

   for (const auto &vec : vecs)
      for (const auto &v : vec)
         if (!MemberQ(ret, v))
            ret.push_back(v);

   return ret;
};

template <class Type>
std::vector<std::vector<Type>> SortVectorChain(const std::vector<std::vector<Type>> &vecs) {  // 2020/03/22
   std::vector<Type> v = *vecs.begin();
   std::vector<std::vector<Type>> rest = Drop(vecs, {0});
   Type knot = *v.rbegin();
   int n = 0;
   std::vector<std::vector<Type>> ret = {v};
   while (rest.size() != 0) {
      std::vector<Type> tmp = rest[n];
      // std::cout << "knot= " << knot << std::endl;
      // std::cout << "tmp= " << tmp << std::endl;
      // std::cout << "MemberQ("<< tmp << "," << knot << ") = " << MemberQ(tmp, knot) << std::endl;
      if (MemberQ(tmp, knot)) {
         // std::cout << tmp << std::endl;
         // std::cout << Green << "RotateLeft(tmp, Position(tmp, knot))= " << RotateLeft(tmp, Position(tmp, knot)) << reset << std::endl;
         //	  v.clear();
         v = RotateLeft(tmp, Position(tmp, knot));
         ret.push_back(v);
         knot = *v.rbegin();

         rest = Drop(rest, {n});

         n = 0;
      } else
         n++;

      // std::cout << Blue << "ret= " << ret << reset << std::endl;
   }
   return ret;
}

template <class T>
std::vector<T> Reverse(std::vector<T> vec) {
   std::reverse(vec.begin(), vec.end());
   return vec;
};
//===========================================================
template <class T>
int isIntersectingSurface(const std::vector<T> &p0,
                          const std::vector<T> &p1,
                          const std::vector<T> &p2,
                          const std::vector<T> &a,
                          const std::vector<T> &b) {
   /* 0:頂点の最大最小の範囲の外で，片方にa,bgがある */
   /* 1:拡大した面には入れているが，多角形の頂点の最大最小範囲にすら入れていない */
   /* 2:a,bは多角形の面と交差していないが，かなり惜しい */
   /* 3:a,bは多角形の面と交差 */
   T maxmin;
   for (auto i = 0; i < 3; i++) {
      maxmin = Max(std::vector<T>{p0[i], p1[i], p2[i]});
      if (maxmin < a[i] && maxmin < b[i])
         return 0;
      maxmin = Min(std::vector<T>{p0[i], p1[i], p2[i]});
      if (maxmin > a[i] && maxmin > b[i])
         return 0;
   }

   std::vector<T> n = Cross(p1 - p0, p2 - p0);
   T d = Dot(Norm(p0 - a) < Norm(p1 - a) ? (p1 - a) : (p0 - a), n) / Dot(b - a, n);
   if (d < 0 || d > 1)
      return 0; /*面に到達できていない*/

   std::vector<T> ps = a + (b - a) * d;

   // ポリゴン頂点の最大最小でチェック
   if (Min(std::vector<T>{p0[0], p1[0], p2[0]}) > ps[0] || ps[0] > Max(std::vector<T>{p0[0], p1[0], p2[0]}))
      return 1; /*面の最大最小範囲にすら入れていない*/
   if (Min(std::vector<T>{p0[1], p1[1], p2[1]}) > ps[1] || ps[1] > Max(std::vector<T>{p0[1], p1[1], p2[1]}))
      return 1; /*面の最大最小範囲にすら入れていない*/
   if (Min(std::vector<T>{p0[2], p1[2], p2[2]}) > ps[2] || ps[2] > Max(std::vector<T>{p0[2], p1[2], p2[2]}))
      return 1; /*面の最大最小範囲にすら入れていない*/

   if (Dot(Cross(p0 - ps, p1 - ps), n) >= 0 && Dot(Cross(p1 - ps, p2 - ps), n) >= 0 && Dot(Cross(p2 - ps, p0 - ps), n) >= 0)
      return 3; /*a,bは面と交差*/
   else
      return 2; /*a,bは面と交差していないが，かなり惜しい*/
   ;
}

template <class T>
std::vector<T> pOnSurface(const std::vector<T> &p0,
                          const std::vector<T> &p1,
                          const std::vector<T> &p2,
                          const std::vector<T> &a,
                          const std::vector<T> &b) {
   std::vector<T> n = Cross(p1 - p0, p2 - p0);
   return a + (b - a) * Dot(/*tangential vector*/ Norm(p0 - a) < Norm(p1 - a) ? (p1 - a) : (p0 - a),
                            /*normal vector*/ n) /
                  Dot(b - a, n);
}
//===========================================================
struct LoadObj {

   std::vector<std::vector<double>> v_complex;
   std::vector<std::vector<float>> s_complex;
   std::vector<std::vector<int>> f_v_complex;  // これのstd::vector<int>が１面 -> 3頂点を指定
   std::vector<bool> hit;

   std::vector<std::vector<float>> s;
   std::vector<std::vector<double>> v, vn;

   std::vector<std::vector<int>> f_v;  // start from 0
   std::vector<std::vector<int>> f_t;
   std::vector<std::vector<int>> f_vn;

   std::vector<int> checkAllIntersection(const std::vector<double> &a,
                                         const std::vector<double> &b) {
      std::vector<int> check0123(f_v_complex.size());
      int i = 0;
      for (const auto &f : f_v_complex) {
         check0123[i] = isIntersectingSurface(v_complex[f[0]],
                                              v_complex[f[1]],
                                              v_complex[f[2]], a, b);

         if (check0123[i] == 3) {
            s_complex[f[0]] = {0., 1., 0., 1.};
            s_complex[f[1]] = {0., 1., 0., 1.};
            s_complex[f[2]] = {0., 1., 0., 1.};
         }
         // if(check0123[i] == 2){
         // 	s_complex[f[0]] += {0.,1.,0.,1.};
         // 	s_complex[f[1]] += {0.,1.,0.,1.};
         // 	s_complex[f[2]] += {0.,1.,0.,1.};
         // }

         i++;
      }
      return check0123;
   }

   std::vector<int> checkAllIntersection(const std::vector<std::vector<double>> &ab) {
      return checkAllIntersection(ab[0], ab[1]);
   }

   std::vector<double> surface(int i, double t0, double t1) {
      std::vector<int> indexes = f_v[i];

      std::vector<std::vector<double>> a(indexes.size());
      for (auto n = 0; n < indexes.size(); n++)
         a[n] = v[indexes[n] - 1];

      return PointsToSurface(t0, t1, a);
   };

   LoadObj(){};
   LoadObj(const std::string &filename) {
      std::vector<std::vector<std::string>> read_line;
      Load(filename, read_line, {"    ", "   ", "  ", " "});
      this->load(read_line);
      generateComplex();
   };

   void generateComplex() {
      int j = 0;
      if (s.empty()) {
         v_complex.clear();
         f_v_complex.clear();
         for (const auto &Ind : f_v) {
            for (const auto &i : Ind) {
               v_complex.push_back(v[i]);
            }
            f_v_complex.push_back(std::vector<int>{j, j + 1, j + 2});
            j = j + 3;
         }
      } else {
         v_complex.clear();
         s_complex.clear();
         f_v_complex.clear();
         for (const auto &Ind : f_v) {
            //	std::vector<float> tmp = {std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),1.};
            std::vector<float> tmp = {.4, .4, .4, 1.};
            for (const auto &i : Ind) {
               v_complex.push_back(v[i]);
               s_complex.push_back(tmp);
            }
            f_v_complex.push_back(std::vector<int>{j, j + 1, j + 2});
            j = j + 3;
         }
      }
      // }

      // s_complex.resize(f_v.size(),std::vector<float>(4,0));
      // int j=0;
      // for(const auto& Ind:f_v)
      //   for(const auto& i:Ind){
      //     s_complex[j] = v[i];
      // 	j++;
      //   }
   };

   // stringを変換しv, f_nに格納
   void load(std::vector<std::vector<std::string>> &read_line) {
      for (auto &line : read_line)
         if (line[0] == "v") {
            line.erase(line.begin());
            v.push_back(stod(line));
         } else if (line[0] == "f") {
            line.erase(line.begin());
            std::vector<int> f_v_tmp, f_t_tmp, f_vn_tmp;
            for (auto &l : line) {
               std::vector<std::string> v_t_vn = StringSplit(l, {"/"});
               if (v_t_vn.size() == 1) {
                  f_v_tmp.push_back(stoi(v_t_vn[0]) - 1);
               } else if (v_t_vn.size() == 2) {
                  f_v_tmp.push_back(stoi(v_t_vn[0]) - 1);
                  if (!v_t_vn[1].empty())
                     f_t_tmp.push_back(stoi(v_t_vn[1]) - 1);
               } else if (v_t_vn.size() == 3) {
                  f_v_tmp.push_back(stoi(v_t_vn[0]) - 1);
                  if (!v_t_vn[1].empty())
                     f_t_tmp.push_back(stoi((v_t_vn[1])) - 1);
                  if (!v_t_vn[2].empty())
                     f_vn_tmp.push_back(stoi(v_t_vn[2]) - 1);
               }
            }
            if (f_v_tmp.size() > 0)
               f_v.push_back(f_v_tmp);
            if (f_t_tmp.size() > 0)
               f_t.push_back(f_t_tmp);
            if (f_vn_tmp.size() > 0)
               f_vn.push_back(f_vn_tmp);
         } else if (line[0] == "vn") {
            line.erase(line.begin());
            vn.push_back(stod(line));
         }
      s.resize(v.size(), std::vector<float>(4, 0));  // RGBS
      for (auto &tmp : s)
         tmp = {.9, .9, .9, 1.};

      // for(auto& tmp:s)
      //   tmp = {std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),1.};
   };

   void ConstructFromString(const std::string &strIN) {
      std::stringstream strm = std::stringstream(strIN);
      std::vector<std::vector<std::string>> read_line;
      std::string read;
      while (std::getline(strm, read))  // loop until the end of file
      {
         if (!StringTrim(read, {" "}).empty()) {
            read_line.push_back(StringTrim(StringSplit(read, {"    ", "   ", "  ", " "}), {"    ", "   ", "  ", " ", "\t", "\n", "\r"}));
         }
      }
      load(read_line);
      generateComplex();
#ifdef _OPENMP
   #pragma omp parallel for
#endif
      glLINES lines;
      for (auto i = 0; i < lines.f_v_complex.size(); i++) {
         checkAllIntersection({lines.v_complex[lines.f_v_complex[i][0]],
                               lines.v_complex[lines.f_v_complex[i][1]]});
      }
   };

   // std::string JSON(){
   //   std::string vertices = "{\n\"vertices\": [\n";
   //   for(auto i = 0; i<v.size(); i++){
   //     for(auto j = 0; j<v[i].size(); j++)
   //   	vertices += std::to_string(float(v[i][j])) + ", ";
   //     if(i != v.size()-1)
   // 	vertices += "\n";
   //     if(i == v.size()-1)
   //       vertices.pop_back();
   //   }
   //   vertices.pop_back();
   //   vertices += "],\n";
   //   std::string indices = "\"indices\": [\n";
   //   for(auto i = 0; i<f_v.size(); i++){
   //     for(auto j = 0; j<f_v[i].size(); j++)
   //       indices += std::to_string(f_v[i][j]-1) + ", ";
   //     if(i != f_v.size()-1)
   //       indices += "\n";
   //     if(i == f_v.size()-1)
   //       indices.pop_back();
   //   }
   //   indices.pop_back();
   //   indices += "]\n}";
   //   return vertices + indices;
   // };

   ////////// 自動的に三角形にする
   // objindex をなおそう
   std::string JSON() {
      ///////////////// vertices
      std::string str = "{\n\"vertices\": [\n";
      for (auto i = 0; i < v.size(); i++) {
         for (auto j = 0; j < v[i].size(); j++)
            str += std::to_string(float(v[i][j])) + ", ";
         if (i != v.size() - 1)
            str += "\n";
         if (i == v.size() - 1)
            str.pop_back();
      }
      str.pop_back();
      str += "],\n";
      /////////////////// scalars
      str += "\"scalars\": [\n";
      for (auto i = 0; i < s.size(); i++) {
         for (auto j = 0; j < s[i].size(); j++)
            str += std::to_string(float(s[i][j])) + ", ";
         if (i != s.size() - 1)
            str += "\n";
         if (i == s.size() - 1)
            str.pop_back();
      }
      str.pop_back();
      str += "],\n";
      ///////////////// vertices_complex
      str += "\n\"vertices_complex\": [\n";
      for (auto i = 0; i < v_complex.size(); i++) {
         for (auto j = 0; j < v_complex[i].size(); j++)
            str += std::to_string(float(v_complex[i][j])) + ", ";
         if (i != v_complex.size() - 1)
            str += "\n";
         if (i == v_complex.size() - 1)
            str.pop_back();
      }
      str.pop_back();
      str += "],\n";
      /////////////////// scalars_complex
      str += "\"scalars_complex\": [\n";
      for (auto i = 0; i < s_complex.size(); i++) {
         for (auto j = 0; j < s_complex[i].size(); j++)
            str += std::to_string(float(s_complex[i][j])) + ", ";
         if (i != s_complex.size() - 1)
            str += "\n";
         if (i == s_complex.size() - 1)
            str.pop_back();
      }
      str.pop_back();
      str += "],\n";
      ////////////////// index_complex
      str += "\"indices_complex\": [\n";
      for (auto i = 0; i < f_v_complex.size(); i++) {
         if (f_v_complex[i].size() == 3) {
            for (auto j = 0; j < f_v_complex[i].size(); j++)
               str += std::to_string(f_v_complex[i][j]) + ", ";
            if (i != f_v_complex.size() - 1)
               str += "\n";
            if (i == f_v_complex.size() - 1)
               str.pop_back();
         } else if (f_v_complex[i].size() == 4) {
            str += std::to_string(f_v_complex[i][0]) + ", " + std::to_string(f_v_complex[i][1]) + ", " + std::to_string(f_v_complex[i][2]) + ", ";
            str += std::to_string(f_v_complex[i][2]) + ", " + std::to_string(f_v_complex[i][3]) + ", " + std::to_string(f_v_complex[i][0]) + ", ";
            ;
            if (i != f_v_complex.size() - 1)
               str += "\n";
            if (i == f_v_complex.size() - 1)
               str.pop_back();
         }
      }
      str.pop_back();
      str += "],\n";
      ////////////////// indices
      str += "\"indices\": [\n";
      for (auto i = 0; i < f_v.size(); i++) {
         if (f_v[i].size() == 3) {
            for (auto j = 0; j < f_v[i].size(); j++)
               str += std::to_string(f_v[i][j]) + ", ";
            if (i != f_v.size() - 1)
               str += "\n";
            if (i == f_v.size() - 1)
               str.pop_back();
         } else if (f_v[i].size() == 4) {
            str += std::to_string(f_v[i][0]) + ", " + std::to_string(f_v[i][1]) + ", " + std::to_string(f_v[i][2]) + ", ";
            str += std::to_string(f_v[i][2]) + ", " + std::to_string(f_v[i][3]) + ", " + std::to_string(f_v[i][0]) + ", ";
            ;
            if (i != f_v.size() - 1)
               str += "\n";
            if (i == f_v.size() - 1)
               str.pop_back();
         }
      }
      str.pop_back();
      str += "]\n}";
      ////////////////////
      return str;
   };
   std::string JSON(const std::string &strIN) {
      ConstructFromString(strIN);
      return JSON();
   };
};
std::string obj2json(const std::string &filename) {
   LoadObj loaded(filename);
   return loaded.JSON();
};
//============================================================
// class MyArray{
// public:
//   std::vector<double> vec;

//   MyArray(){
//   }

//   std::vector<double> getvec(){
//     return vec;
//   }

//   void setDoubleVec(double& address, int len){
//     vec.resize(len);
//     for(auto i=0; i<len; i++)
//       vec[i] = *(address + i);
//   }

//   std::vector<double> show(){
//     return vec;
//   }
// };

//===========================================================
float divide(float x, float y) {
   return x / y;
};

/////// emscripten //////
#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(stl_wrappers) {
   // ベクトルのバインドを登録する必要がある．
   emscripten::register_vector<int>("VectorInt");
   emscripten::register_vector<double>("VectorDouble");
   emscripten::register_vector<float>("VectorFloat");
   emscripten::register_vector<std::vector<double>>("VVDouble");
   emscripten::register_vector<std::vector<int>>("VVInt");
   emscripten::register_vector<std::vector<float>>("VVFloat");
}
// チェイン形式で，emscripten::class_のメンバー関数らによって，関数，属性などをバインドしていく
// 必要なだけバインドし，コードのサイズを少なくするように心がける．（https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html）
EMSCRIPTEN_BINDINGS(my_module) {

   emscripten::class_<glLINES>("glLINES")
       .constructor<>()
       .property("v_complex", &glLINES::v_complex)
       .property("s_complex", &glLINES::s_complex)
       .property("f_v_complex", &glLINES::f_v_complex);

   emscripten::class_<LoadObj>("LoadObj")
       .constructor<>()
       .function("ConstructFromString", &LoadObj::ConstructFromString)
       .function("JSON", emscripten::select_overload<std::string()>(&LoadObj::JSON))
       .function("OBJ2JSON", emscripten::select_overload<std::string(const std::string &)>(&LoadObj::JSON))
       .property("v_complex", &LoadObj::v_complex)
       .property("s_complex", &LoadObj::s_complex)
       .property("f_v_complex", &LoadObj::f_v_complex)
       .property("v", &LoadObj::v)
       .property("s", &LoadObj::s)
       .property("f_v", &LoadObj::f_v);
   /*
   オーバーロードがある場合：
   emscripten::select_overload<>
 */

   emscripten::function("isIntersectingSurfaceVec",
                        emscripten::select_overload<int(const std::vector<float> &,
                                                        const std::vector<float> &,
                                                        const std::vector<float> &,
                                                        const std::vector<float> &,
                                                        const std::vector<float> &)>(&isIntersectingSurface));
   emscripten::function("divide", divide);

   emscripten::function("FlattenVVD",
                        emscripten::select_overload<std::vector<double>(const std::vector<std::vector<double>> &)>(&Flatten));
   emscripten::function("FlattenVVI",
                        emscripten::select_overload<std::vector<int>(const std::vector<std::vector<int>> &)>(&Flatten));
   emscripten::function("FlattenVVF",
                        emscripten::select_overload<std::vector<float>(const std::vector<std::vector<float>> &)>(&Flatten));

   // emscripten::class_<MyArray>("MyArray")
   //   .constructor<>()
   //   .function("getvector", &MyArray::getvec)
   //   .function("setDoubleVec", &MyArray::setDoubleVec)
   //   .function("show", &MyArray::show)
   //   ;

   /*
   vectorを使いたい場合
   https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html?highlight=memory
 */
}
#endif

#endif
