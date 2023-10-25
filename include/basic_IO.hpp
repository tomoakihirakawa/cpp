#ifndef basic_IO_H
#define basic_IO_H

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <unordered_set>
#include <vector>
#include "basic_constants.hpp"

// template <typename... T>
// std::ostream &operator<<(std::ostream &stream, const std::tuple<T...> &V) {
//    constexpr std::size_t n = std::tuple_size<std::tuple<T...>>::value;  // std::tuple_size<decltype(V)>::value;はエラー　
//    std::size_t i = 0;
//    stream << "{";
//    for_each(V, [&](const auto &v) {
//       if constexpr (i++ == n - 1)
//          stream << v;
//       else
//          stream << v << ",";
//    });
//    stream << "}";
//    return stream;
// };

template <size_t N, typename T, typename STREAM>
STREAM &operator<<(STREAM &stream, const std::array<T, N> &V) {
   stream << "{";
   bool first = true;
   std::ranges::for_each(V, [&](const auto &v) {
      stream << (first ? "" : ",") << v;  // fold expressionを使えば，別にタプルの長さをチェック必要はない
      first = false;
   });
   stream << "}";
   return stream;
}

template <size_t N, typename T, typename STREAM>
STREAM &operator<<(STREAM &stream, const std::vector<std::array<T, N>> &V) {
   stream << "{";
   bool first = true;
   std::ranges::for_each(V, [&](const auto &v) {
      stream << (first ? "" : ",") << v;  // fold expressionを使えば，別にタプルの長さをチェック必要はない
      first = false;
   });

   stream << "}";
   return stream;
}

template <typename... T>
std::ostream &operator<<(std::ostream &stream, const std::tuple<T...> &V) {
   stream << "{";
   std::apply([&stream](const auto &...v) {
      bool first = true;
      ((stream << (first ? "" : ",") << v, first = false), ...);  // fold expressionを使えば，別にタプルの長さをチェック必要はない
   },
              V);
   stream << "}";
   return stream;
}

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
   /* ----------------------------------------------- */
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

// template <typename T, typename U>
// std::ostream &operator<<(std::ostream &stream, const tuple_of<std::tuple_size<T>::value,
//                                                               tuple_of<std::tuple_size<T>::value, U>> &V) {
//    stream << "{";
//    for_each(v, [&stream](const auto &U) { stream << u << ","; });
//    stream << "}";
//    return stream;
// };

// template <typename T>
// std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T> &v) { return (stream << "{" << std::get<0>(v) << "," << std::get<1>(v) << "}"); };

// template <typename T>
// std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T, T> &v) {
//    return (stream << "{"
//                   << std::get<0>(v) << ","
//                   << std::get<1>(v) << ","
//                   << std::get<2>(v) << "}");
// };

// template <typename T>
// std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T, T, T> &v) {
//    return (stream << "{"
//                   << std::get<0>(v) << ","
//                   << std::get<1>(v) << ","
//                   << std::get<2>(v) << ","
//                   << std::get<3>(v) << "}");
// };

// template <typename T>
// std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T, T, T, T> &v) {
//    return (stream << "{"
//                   << std::get<0>(v) << ","
//                   << std::get<1>(v) << ","
//                   << std::get<2>(v) << ","
//                   << std::get<3>(v) << ","
//                   << std::get<4>(v) << "}");
// };

// template <typename T>
// std::ostream &operator<<(std::ostream &stream, const std::tuple<T, T, T, T, T, T> &v) {
//    return (stream << "{"
//                   << std::get<0>(v) << ","
//                   << std::get<1>(v) << ","
//                   << std::get<2>(v) << ","
//                   << std::get<3>(v) << ","
//                   << std::get<4>(v) << ","
//                   << std::get<5>(v) << "}");
// };

// template <typename T>
// std::ostream &operator<<(std::ostream &stream, const std::vector<T> &v) {
//    stream << "{";
//    if (!v.empty()) {
//       for (size_t i = 0; i < v.size() - 1; i++)
//          stream << v[i] << ",";
//       stream << *v.rbegin();
//    }
//    stream << "}";
//    return stream;
// };

// std::ostream &operator<<(std::ostream &stream, const T3Tdd &v) {
//    auto [v0, v1, v2] = v;
//    stream << "{{" << std::get<0>(v0) << "," << std::get<1>(v0) << "},";
//    stream << "{" << std::get<0>(v1) << "," << std::get<1>(v1) << "},";
//    stream << "{" << std::get<0>(v2) << "," << std::get<1>(v2) << "}}";
//    return stream;
// };
// std::ostream &operator<<(std::ostream &stream, const T3Tddd &v) {
//    auto [v0, v1, v2] = v;
//    stream << "{{" << std::get<0>(v0) << "," << std::get<1>(v0) << "," << std::get<2>(v0) << "},";
//    stream << "{" << std::get<0>(v1) << "," << std::get<1>(v1) << "," << std::get<2>(v1) << "},";
//    stream << "{" << std::get<0>(v2) << "," << std::get<1>(v2) << "," << std::get<2>(v2) << "}}";
//    return stream;
// };

/* ------------------------------------------------------ */

// template <typename T>
// void Print(const T *const v_IN) { std::cout << v_IN << std::endl; };
template <typename T>
void Print(const T &v_IN) {
   std::cout << v_IN << std::endl;
};
template <typename T, typename U>
void Print(const T &v_IN, const U &color) {
   std::cout << color << v_IN << colorOff << std::endl;
};
template <typename... Args>
void Print(Args... args) {
   std::stringstream stream;
   (stream << ... << args);
   std::cout << stream.str() << colorOff << std::endl;
}
//
template <typename T>
void DebugPrint(const T &v_IN) {
#if defined(_debugging_)
   std::cout << v_IN << std::endl;
#endif
};
template <typename T, typename U>
void DebugPrint(const T &v_IN, const U &color) {
#if defined(_debugging_)
   std::cout << color << v_IN << colorOff << std::endl;
#endif
};
//
template <typename... Args>
void DebugPrint(Args... args) {
#if defined(_debugging_)
   std::stringstream stream;
   (stream << ... << args);
   std::cout << stream.str() << colorOff << std::endl;
#endif
}

// template <typename T>
// void MatrixForm(const std::vector<std::vector<T>> &mat, const int n = 4, const int m = 7) {
//    std::stringstream ss;
//    ss << std::setprecision(n) << Red << "{" << colorOff;
//    for (auto it = mat.begin(); it != mat.end(); ++it) {
//       ss << Green << (it == mat.begin() ? "{" : " {") << colorOff;
//       for (auto jt = (*it).begin(); jt != ((*it).end() - 1); ++jt)
//          ss << std::setw(m) << std::setfill(' ') << (*jt) << Green << "," << colorOff;
//       ss << std::setw(m) << std::setfill(' ') << (*((*it).end() - 1)) << Green << "}" << colorOff;
//       if (it != mat.end() - 1)
//          ss << ",\n";
//    }
//    ss << Red << "}" << colorOff;
//    std::cout << ss.str() << std::endl;
// };

#include <functional>

template <typename T, typename Func>
std::string MatrixForm(const T &mat,
                       Func color,
                       const int prec,
                       const int width) {

   auto COL = [](const auto &x, const auto &Red, const int w = 0) {
      std::stringstream ss;
      ss << Red << std::setw(w) << std::setfill(' ') << x << colorOff;
      return ss.str();
   };

   std::stringstream ss;
   int i = 0;
   ss << std::setprecision(prec) << COL("{", Magenta);
   for (auto it = mat.begin(); it != mat.end(); ++it) {
      int j = 0;
      ss << Green << (it == mat.begin() ? "{" : " {") << colorOff;
      for (auto jt = (*it).begin(); jt != ((*it).end() - 1); ++jt) {
         if (color(i, j))
            ss << COL(*jt, Red, width) << COL(",", Green);
         else
            ss << std::setw(width) << std::setfill(' ') << (*jt) << COL(",", Green);
         j++;
      }

      if (color(i, j))
         ss << COL((*((*it).end() - 1)), Red, width) << COL("}", Green);
      else
         ss << std::setw(width) << std::setfill(' ') << (*((*it).end() - 1)) << COL("}", Green);

      if (it != mat.end() - 1) {
         ss << ",\n";
         j++;
      }
      i++;
   }
   ss << COL("}", Magenta);
   return ss.str();
}

template <typename T>
std::string MatrixForm(const T &mat,
                       const int prec = 4,
                       const int width = 7) {
   return MatrixForm(
       mat, [](const auto &i, const auto &j) { return false; },
       prec, width);
}

template <typename T, typename U>
std::string MatrixForm(const T &mat, const U &w) {
   return MatrixForm(
       mat, [](const auto &i, const auto &j) { return false; }, 4, w);
};

/* -------------------------------------------------------------------------- */
#include <sys/utsname.h>
#include <unistd.h>

void logMachineInformation(std::ofstream &ofs) {
   char hostname[256];
   if (gethostname(hostname, sizeof(hostname)) != 0) {
      ofs << "Failed to get hostname\n";
      return;
   }

   struct utsname unameData;
   if (uname(&unameData) != 0) {
      ofs << "Failed to get system information\n";
      return;
   }

   ofs << "Machine Information:" << std::endl;
   ofs << "  Hostname: " << hostname << std::endl;
   ofs << "  Operating System: " << unameData.sysname << " " << unameData.release << std::endl;
   ofs << "  Working Directory: " << std::filesystem::current_path() << std::endl;  // Added current directory
   ofs << std::endl;
}

bool initializeLogFile(const std::string &logfilename, int argc, char **argv) {
   std::vector<std::string> args(argv, argv + argc);
   // Read existing content
   std::ifstream ifs(logfilename);
   std::string existingContent;
   if (ifs.is_open()) {
      existingContent.assign(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
   }
   ifs.close();

   // Clear and open log file
   std::ofstream ofs(logfilename);
   if (!ofs.is_open()) {
      std::cerr << "Failed to open " << logfilename << std::endl;
      return false;
   }
   ofs << std::setw(10) << std::setfill('-') << "" << std::endl;

   // Log time of execution
   auto now = std::chrono::system_clock::now();
   std::time_t time = std::chrono::system_clock::to_time_t(now);
   ofs << "Execution time: " << std::ctime(&time) << std::endl;
   // Log machine info
   logMachineInformation(ofs);

   ofs << "Command line arguments:" << std::endl;
   for (size_t i = 0; i < args.size() && i < 5; ++i) {  // Changed '||' to '&&' here
      ofs << "  arg" << i << ": " << args[i] << std::endl;
   }
   // Add existing content back
   ofs << existingContent;
   ofs << std::setw(10) << std::setfill('-') << "" << std::endl;
   ofs.close();
   return true;
}

std::string getUserName() {
   const char *username_cstr = std::getenv("USER");
   if (username_cstr == nullptr) {
      username_cstr = std::getenv("USERNAME");
   }

   std::string username;
   if (username_cstr != nullptr) {
      username = std::string(username_cstr);
   } else {
      throw std::runtime_error("Failed to get username");
   }

   return username;
}

std::string getMachineName() {
   char hostname[1024];
   hostname[1023] = '\0';
   if (gethostname(hostname, 1023) == -1) {
      throw std::runtime_error("Failed to get machine name");
   }

   return std::string(hostname);
}
/* -------------------------------------------------------------------------- */

#endif
