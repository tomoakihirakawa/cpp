#ifndef basic_arithmetic_array_operations_H
#define basic_arithmetic_array_operations_H

#include <array>
#include "basic_alias.hpp"

// b! -------------------------------------------------------------------------- */
// b!                                 + - * /                                    */
// b! -------------------------------------------------------------------------- */
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator-(std::array<T, N> arr /*copy*/) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) *= -1;
      operator-<M + 1, N>(arr);
   }
   return arr;
};
//@ ---------------------------- 1D array vs T--------------------------- */
template <size_t M = 0, size_t N, typename T>
std::array<T, N> &operator+=(std::array<T, N> &arr /*reference*/, const T &d) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) += d;
      operator+=<M + 1, N>(arr, d);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> &operator-=(std::array<T, N> &arr /*reference*/, const T &d) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) -= d;
      operator-=<M + 1, N>(arr, d);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> &operator*=(std::array<T, N> &arr /*reference*/, const T &d) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) *= d;
      operator*=<M + 1, N>(arr, d);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> &operator/=(std::array<T, N> &arr /*reference*/, const T &d) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) /= d;
      operator/=<M + 1, N>(arr, d);
   }
   return arr;
};
/* -------------------------------------------------------------------------- */
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator+(std::array<T, N> arr /*copy*/, const T &d) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) += d;
      operator+=<M + 1, N>(arr, d);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator-(std::array<T, N> arr /*copy*/, const T &d) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) -= d;
      operator-=<M + 1, N>(arr, d);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator*(std::array<T, N> arr /*copy*/, const T &d) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) *= d;
      operator*=<M + 1, N>(arr, d);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator/(std::array<T, N> arr /*copy*/, const T &d) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) /= d;
      operator/=<M + 1, N>(arr, d);
   }
   return arr;
};
/* -------------------------------------------------------------------------- */
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator+(const T &d, std::array<T, N> arr /*copy*/) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) += d;
      operator+=<M + 1, N>(arr, d);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator-(const T &d, std::array<T, N> arr /*copy*/) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) = d - std::get<M>(arr);
      operator-=<M + 1, N>(d, arr);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator*(const T &d, std::array<T, N> arr /*copy*/) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) *= d;
      operator*=<M + 1, N>(arr, d);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> operator/(const T &d, std::array<T, N> arr /*copy*/) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) = d / std::get<M>(arr);
      operator/=<M + 1, N>(d, arr);
   }
   return arr;
};
//@ -------------------------------------------------------------------------- */
//@ -------------------------------------------------------------------------- */
//@ -------------------------------------------------------------------------- */
//! --------------------------- 1d array - 1d array -------------------------- */
template <size_t M = 0, size_t N, typename T>
std::array<T, N> &operator+=(std::array<T, N> &arr /*reference*/, const std::array<T, N> &ARR) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) += std::get<M>(ARR);
      operator+=<M + 1, N>(arr, ARR);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> &operator-=(std::array<T, N> &arr /*reference*/, const std::array<T, N> &ARR) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) -= std::get<M>(ARR);
      operator-=<M + 1, N>(arr, ARR);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> &operator*=(std::array<T, N> &arr /*reference*/, const std::array<T, N> &ARR) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) *= std::get<M>(ARR);
      operator*=<M + 1, N>(arr, ARR);
   }
   return arr;
};
template <size_t M = 0, size_t N, typename T>
std::array<T, N> &operator/=(std::array<T, N> &arr /*reference*/, const std::array<T, N> &ARR) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      std::get<M>(arr) /= std::get<M>(ARR);
      operator/=<M + 1, N>(arr, ARR);
   }
   return arr;
};

template <size_t N, typename T>
std::array<T, N> operator+(std::array<T, N> arr /*copy*/, const std::array<T, N> &ARR) { return arr += ARR; };
template <size_t N, typename T>
std::array<T, N> operator-(std::array<T, N> arr /*copy*/, const std::array<T, N> &ARR) { return arr -= ARR; };
template <size_t N, typename T>
std::array<T, N> operator*(std::array<T, N> arr /*copy*/, const std::array<T, N> &ARR) { return arr *= ARR; };
template <size_t N, typename T>
std::array<T, N> operator/(std::array<T, N> arr /*copy*/, const std::array<T, N> &ARR) { return arr /= ARR; };

//! ---------------------------- 2d array - T --------------------------- */

template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator+=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const T &d) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) += d;  //! which means std::array<T, N1> += d
      operator+=<M + 1>(arrarr, d);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator-=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const T &d) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) -= d;
      operator-=<M + 1>(arrarr, d);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator*=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const T &d) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) *= d;
      operator*=<M + 1>(arrarr, d);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator/=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const T &d) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) /= d;
      operator/=<M + 1>(arrarr, d);
   }
   return arrarr;
};

template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator+(std::array<std::array<T, N1>, N0> arrarr /*copy*/, const T &d) { return arrarr += d; };
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator-(std::array<std::array<T, N1>, N0> arrarr /*copy*/, const T &d) { return arrarr -= d; };
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator*(std::array<std::array<T, N1>, N0> arrarr /*copy*/, const T &d) { return arrarr *= d; };
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator/(std::array<std::array<T, N1>, N0> arrarr /*copy*/, const T &d) { return arrarr /= d; };

template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator+(const T &d, std::array<std::array<T, N1>, N0> arrarr /*copy*/) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) += d;  //! which means std::array<T, N1> += d
      operator+=<M + 1>(arrarr, d);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator-(const T &d, std::array<std::array<T, N1>, N0> arrarr /*copy*/) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) = std::get<M>(arrarr) - d;
      operator-=<M + 1>(arrarr, d);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator*(const T &d, std::array<std::array<T, N1>, N0> arrarr /*copy*/) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) *= d;
      operator*=<M + 1>(arrarr, d);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator/(const T &d, std::array<std::array<T, N1>, N0> arrarr /*copy*/) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) = std::get<M>(arrarr) / d;
      operator/=<M + 1>(arrarr, d);
   }
   return arrarr;
};
//! --------------------------- 2d array - 1d array -------------------------- */
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator+=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const std::array<T, N1> &ARR) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) += ARR;
      operator+=<M + 1>(arrarr, ARR);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator-=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const std::array<T, N1> &ARR) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) -= ARR;
      operator-=<M + 1>(arrarr, ARR);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator*=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const std::array<T, N1> &ARR) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) *= ARR;
      operator*=<M + 1>(arrarr, ARR);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator/=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const std::array<T, N1> &ARR) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) /= ARR;
      operator/=<M + 1>(arrarr, ARR);
   }
   return arrarr;
};
template <size_t M, size_t N, typename T>
std::array<std::array<T, N>, M> operator+(std::array<std::array<T, N>, M> arrarr /*copy*/, const std::array<T, N> &ARR) { return arrarr += ARR; };
template <size_t M, size_t N, typename T>
std::array<std::array<T, N>, M> operator-(std::array<std::array<T, N>, M> arrarr /*copy*/, const std::array<T, N> &ARR) { return arrarr -= ARR; };
template <size_t M, size_t N, typename T>
std::array<std::array<T, N>, M> operator*(std::array<std::array<T, N>, M> arrarr /*copy*/, const std::array<T, N> &ARR) { return arrarr *= ARR; };
template <size_t M, size_t N, typename T>
std::array<std::array<T, N>, M> operator/(std::array<std::array<T, N>, M> arrarr /*copy*/, const std::array<T, N> &ARR) { return arrarr /= ARR; };
//! --------------------------- 2d array - 2d array -------------------------- */
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator+=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const std::array<std::array<T, N1>, N0> &ARRARR) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) += std::get<M>(ARRARR);
      operator+=<M + 1>(arrarr, ARRARR);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator-=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const std::array<std::array<T, N1>, N0> &ARRARR) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) -= std::get<M>(ARRARR);
      operator-=<M + 1>(arrarr, ARRARR);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator*=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const std::array<std::array<T, N1>, N0> &ARRARR) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) *= std::get<M>(ARRARR);
      operator*=<M + 1>(arrarr, ARRARR);
   }
   return arrarr;
};
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> &operator/=(std::array<std::array<T, N1>, N0> &arrarr /*reference*/, const std::array<std::array<T, N1>, N0> &ARRARR) {
   if constexpr (M < N0) {
      std::get<M>(arrarr) /= std::get<M>(ARRARR);
      operator/=<M + 1>(arrarr, ARRARR);
   }
   return arrarr;
};

template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator+(std::array<std::array<T, N1>, N0> arrarr /*copy*/, const std::array<std::array<T, N1>, N0> &ARRARR) { return arrarr += ARRARR; };
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator-(std::array<std::array<T, N1>, N0> arrarr /*copy*/, const std::array<std::array<T, N1>, N0> &ARRARR) { return arrarr -= ARRARR; };
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator*(std::array<std::array<T, N1>, N0> arrarr /*copy*/, const std::array<std::array<T, N1>, N0> &ARRARR) { return arrarr *= ARRARR; };
template <size_t M = 0, size_t N0, size_t N1, typename T>
std::array<std::array<T, N1>, N0> operator/(std::array<std::array<T, N1>, N0> arrarr /*copy*/, const std::array<std::array<T, N1>, N0> &ARRARR) { return arrarr /= ARRARR; };

/* -------------------------------------------------------------------------- */
/*                       Mathematical Vector Operations                       */
/* -------------------------------------------------------------------------- */

template <size_t N, typename T>
T Total(const std::array<T, N> &arr) {
   T ret = 0;
   for_each(arr, [&ret](const auto &a) { ret += a; });
   return ret;
}

template <size_t M = 0, size_t N, typename T>
T Dot(const std::array<T, N> &arr, const std::array<T, N> &ARR) {
   T ret = 0;
   if constexpr (M < N) {
      ret += std::get<M>(arr) * std::get<M>(ARR);
      ret += Dot<M + 1>(arr, ARR);
   }
   return ret;
}

template <size_t N0, size_t N1, typename T>
std::array<T, N1> Dot(const std::array<T, N0> &arr, const std::array<std::array<T, N1>, N0> &ARRARR) {
   /*
   {a,b,c}.{{A0,A1,A2,A3},{B0,B1,B2,B3},{C0,C1,C2,C3}}
   */
   std::array<T, N1> ret = {};

   // for (size_t i = 0; i < N0; ++i) {
   //    for (size_t j = 0; j < N1; ++j) {
   //       ret[j] += arr[i] * ARR[i][j];
   //    }
   // }

   // for_each(arr, ARR, [&](const auto &ar, const auto &AR) {
   //    // for (size_t j = 0; j < N1; ++j) {
   //    //    ret[j] += ar * AR[j];
   //    // }
   //    for_each(ret, AR, [&](auto &r, auto &A) {
   //       r += ar * A;
   //    });
   // });

   for_each(ARRARR, [&](const auto &ARR) {
      for_each(ret, arr, ARR, [&](auto &r, auto &a, auto &A) {
         r += a * A;
      });
   });

   return ret;
}

template <size_t N0, size_t N1, typename T>
std::array<T, N1> Dot(const std::array<std::array<T, N0>, N1> &arr, const std::array<T, N0> &ARR) {
   /*
   {{A0,A1,A2},{B0,B1,B2},{C0,C1,C2},{D0,D1,D2}}.{a,b,c}={a*A0+b*A1+c*A2,a*B0+b*B1+c*B2,a*C0+b*C1+c*C2,a*D0+b*D1+c*D2}
   */
   std::array<T, N1> ret = {};
   for_each(ret, arr, [&](auto &r, auto &ar) { r = Dot(ar, ARR); });

   return ret;
}

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

template <size_t M = 0, size_t N, typename T>
T Norm(const std::array<T, N> &arr) { return std::sqrt(Dot(arr, arr)); }

template <size_t M = 0, size_t N, typename T>
T RootMeanSquare(const std::array<T, N> &arr) { return std::sqrt(Dot(arr, arr) / N); }

template <typename T>
std::array<T, 3> Cross(const std::array<T, 3> &A, const std::array<T, 3> &B) {
   return {{std::get<1>(A) * std::get<2>(B) - std::get<2>(A) * std::get<1>(B),
            std::get<2>(A) * std::get<0>(B) - std::get<0>(A) * std::get<2>(B),
            std::get<0>(A) * std::get<1>(B) - std::get<1>(A) * std::get<0>(B)}};
};

#endif