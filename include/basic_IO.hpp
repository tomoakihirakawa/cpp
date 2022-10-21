#ifndef basic_IO_H
#define basic_IO_H

#include <iostream>
#include <sstream>
#include <unordered_set>
#include <vector>
#include "basic_constants.hpp"

std::ostream &operator<<(std::ostream &stream, const std::vector<std::string> &v) {
   stream << "{";
   if (!v.empty()) {
      for (size_t i = 0; i < v.size() - 1; i++)
         stream << v[i] << ",";
      stream << *v.rbegin();
   }
   stream << "}";
   return stream;
};

template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::vector<T *> &v) {
   std::stringstream ss;
   ss << "{";
   if (!v.empty()) {
      for (size_t i = 0; i < v.size() - 1; i++)
         ss << v[i] << ",";
      ss << *v.rbegin();
   }
   ss << "}";
   stream << ss.str();
   return stream;
   // stream << "{";
   // if (!v.empty())
   // {
   //     for (size_t i = 0; i < v.size() - 1; i++)
   //         stream << v[i] << ",";
   //     stream << *v.rbegin();
   // }
   // stream << "}";
   // return stream;
};

template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::unordered_set<T *> &v) {
   int i = 0;
   std::vector<T *> w(v.size());
   for (const auto &u : v)
      w[i++] = u;
   stream << w;
   return stream;
};

std::ostream &operator<<(std::ostream &stream, const std::vector<bool> &v) {
   stream << "{";
   if (!v.empty()) {
      for (size_t i = 0; i < v.size() - 1; i++)
         stream << v[i] << ",";
      stream << *v.rbegin();
   }
   stream << "}";
   return stream;
};

std::ostream &operator<<(std::ostream &stream, const std::vector<int> &v) {
   stream << "{";
   if (!v.empty()) {
      for (size_t i = 0; i < v.size() - 1; i++)
         stream << v[i] << ",";
      stream << *v.rbegin();
   }
   stream << "}";
   return stream;
};

std::ostream &operator<<(std::ostream &stream, const std::vector<double> &v) {
   stream << "{";
   if (!v.empty()) {
      for (size_t i = 0; i < v.size() - 1; i++)
         stream << v[i] << ",";
      stream << *v.rbegin();
   }
   stream << "}";
   return stream;
};

std::ostream &operator<<(std::ostream &stream, const std::vector<std::vector<int>> &v) {
   stream << "{";
   if (!v.empty()) {
      for (size_t i = 0; i < v.size() - 1; i++)
         stream << v[i] << ",";
      stream << *v.rbegin();
   }
   stream << "}";
   return stream;
};
std::ostream &operator<<(std::ostream &stream, const std::vector<std::vector<double>> &v) {
   stream << "{";
   if (!v.empty()) {
      for (size_t i = 0; i < v.size() - 1; i++)
         stream << v[i] << ",";
      stream << *v.rbegin();
   }
   stream << "}";
   return stream;
};
// std::stringstream &operator<<(std::stringstream &stream, const std::tuple<double, double, double> &v)
// {
// 	stream << "{" << std::get<0>(v) << "," << std::get<1>(v) << "," << std::get<2>(v) << "}";
// 	return stream;
// };
// std::stringstream &operator<<(std::stringstream &stream, const std::vector<std::tuple<double, double, double>> &v)
// {
// 	for (const auto &u : v)
// 		stream << "{" << u << "}";
// 	return stream;
// };

template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T> &v) { return (stream << "{" << std::get<0>(v) << "," << std::get<1>(v) << "}"); };

template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T, T> &v) {
   return (stream << "{"
                  << std::get<0>(v) << ","
                  << std::get<1>(v) << ","
                  << std::get<2>(v) << "}");
};

template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T, T, T> &v) {
   return (stream << "{"
                  << std::get<0>(v) << ","
                  << std::get<1>(v) << ","
                  << std::get<2>(v) << ","
                  << std::get<3>(v) << "}");
};

template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T, T, T, T> &v) {
   return (stream << "{"
                  << std::get<0>(v) << ","
                  << std::get<1>(v) << ","
                  << std::get<2>(v) << ","
                  << std::get<3>(v) << ","
                  << std::get<4>(v) << "}");
};

template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T, T, T, T, T> &v) {
   return (stream << "{"
                  << std::get<0>(v) << ","
                  << std::get<1>(v) << ","
                  << std::get<2>(v) << ","
                  << std::get<3>(v) << ","
                  << std::get<4>(v) << ","
                  << std::get<5>(v) << "}");
};

std::ostream &operator<<(std::ostream &stream, const T3Tdd &v) {
   auto [v0, v1, v2] = v;
   stream << "{{" << std::get<0>(v0) << "," << std::get<1>(v0) << "},";
   stream << "{" << std::get<0>(v1) << "," << std::get<1>(v1) << "},";
   stream << "{" << std::get<0>(v2) << "," << std::get<1>(v2) << "}}";
   return stream;
};
std::ostream &operator<<(std::ostream &stream, const T3Tddd &v) {
   auto [v0, v1, v2] = v;
   stream << "{{" << std::get<0>(v0) << "," << std::get<1>(v0) << "," << std::get<2>(v0) << "},";
   stream << "{" << std::get<0>(v1) << "," << std::get<1>(v1) << "," << std::get<2>(v1) << "},";
   stream << "{" << std::get<0>(v2) << "," << std::get<1>(v2) << "," << std::get<2>(v2) << "}}";
   return stream;
};

/* ------------------------------------------------------ */

template <typename T>
void Print(const T *const v_IN) { std::cout << v_IN << std::endl; };
template <typename T>
void Print(const T &v_IN) { std::cout << v_IN << std::endl; };
template <typename T, typename U>
void Print(const T &v_IN, const U &color) { std::cout << color << v_IN << colorOff << std::endl; };
//
template <typename T>
void MatrixForm(const std::vector<std::vector<T>> &mat) {
   std::stringstream ss;
   ss << Red << "{" << colorOff;
   for (auto it = mat.begin(); it != mat.end(); ++it) {
      ss << Green << "{" << colorOff;
      for (auto jt = (*it).begin(); jt != ((*it).end() - 1); ++jt)
         ss << (*jt) << Green << "," << colorOff;
      ss << (*((*it).end() - 1)) << Green << "}" << colorOff;
      if (it != mat.end() - 1)
         ss << ",\n";
   }
   ss << Red << "}" << colorOff;
   std::cout << ss.str() << std::endl;
};

template <typename T, typename U>
void MatrixForm(const std::vector<std::vector<T>> &mat, const U &w) {
   std::stringstream ss;
   ss << Red << "{" << colorOff;
   for (auto it = mat.begin(); it != mat.end(); ++it) {
      ss << Green << "{" << colorOff;
      for (auto jt = (*it).begin(); jt != ((*it).end() - 1); ++jt)
         ss << w << (*jt) << Green << "," << colorOff;
      ss << w << (*((*it).end() - 1)) << Green << "}" << colorOff;
      if (it != mat.end() - 1)
         ss << ",\n";
   }
   ss << Red << "}" << colorOff;
   std::cout << ss.str() << std::endl;
};

#endif
