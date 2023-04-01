#ifndef basic_alias_H
#define basic_alias_H

#include <string>
#include <tuple>
#include <vector>
/* ------------------------------------------------------ */
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
/* ------------------------------------------------------ */
using V_s = std::vector<std::string>;
using VV_s = std::vector<std::vector<std::string>>;
using VVV_s = std::vector<std::vector<std::string>>;
/* ------------------------------------------------------ */
using V_i = std::vector<int>;
using VV_i = std::vector<std::vector<int>>;
using VVV_i = std::vector<std::vector<std::vector<int>>>;
/* ------------------------------------------------------ */
using V_Tdd = std::vector<std::tuple<double, double>>;
using V_Tddd = std::vector<std::tuple<double, double, double>>;
using VV_Tdd = std::vector<V_Tdd>;
using VV_Tddd = std::vector<V_Tddd>;
/* ------------------------------------------------------ */
using Tdd = std::tuple<double, double>;
using Tddd = std::tuple<double, double, double>;
using T4d = std::tuple<double, double, double, double>;
using T5d = std::tuple<double, double, double, double, double>;
using T6d = std::tuple<double, double, double, double, double, double>;
using T7d = std::tuple<double, double, double, double, double, double, double>;
using T8d = std::tuple<double, double, double, double, double, double, double, double>;
using T9d = std::tuple<double, double, double, double, double, double, double, double, double>;
using T10d = std::tuple<double, double, double, double, double, double, double, double, double, double>;
using T11d = std::tuple<double, double, double, double, double, double, double, double, double, double, double>;
using T12d = std::tuple<double, double, double, double, double, double, double, double, double, double, double, double>;
using T13d = std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double>;
using T14d = std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double>;
using T15d = std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double, double>;
using T16d = std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double>;
using T17d = std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double>;
//

using Tii = std::tuple<int, int>;
using Tiii = std::tuple<int, int, int>;
using T4i = std::tuple<int, int, int, int>;
using T5i = std::tuple<int, int, int, int, int>;
using T6i = std::tuple<int, int, int, int, int, int>;
/* ------------------------------------------------------ */
using T2Tdd = std::tuple<Tdd, Tdd>;
using T3Tdd = std::tuple<Tdd, Tdd, Tdd>;
using T4Tdd = std::tuple<Tdd, Tdd, Tdd, Tdd>;
using T5Tdd = std::tuple<Tdd, Tdd, Tdd, Tdd, Tdd>;
using T6Tdd = std::tuple<Tdd, Tdd, Tdd, Tdd, Tdd, Tdd>;
//
using T2Tddd = std::tuple<Tddd, Tddd>;
using T3Tddd = std::tuple<Tddd, Tddd, Tddd>;
using T4Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd>;
using T5Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd>;
using T6Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T7Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T8Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T9Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T10Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T11Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T12Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T13Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T14Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T15Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T16Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T17Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T18Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T19Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T20Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T21Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T22Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T23Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T24Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T25Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T64Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd,
                           Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd,
                           Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd,
                           Tddd, Tddd, Tddd, Tddd>;
//
using T2T6d = std::tuple<T6d, T6d>;
using T3T4d = std::tuple<T4d, T4d, T4d>;
using T4T4d = std::tuple<T4d, T4d, T4d, T4d>;
using T4Tiii = std::tuple<Tiii, Tiii, Tiii, Tiii>;
using T6T6d = std::tuple<T6d, T6d, T6d, T6d, T6d, T6d>;
using T7T7d = std::tuple<T7d, T7d, T7d, T7d, T7d, T7d, T7d>;
//
using T4T3Tddd = std::tuple<T3Tddd, T3Tddd, T3Tddd, T3Tddd>;
//
using T6T2Tddd = std::tuple<T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd>;
using T6T4Tddd = std::tuple<T4Tddd, T4Tddd, T4Tddd, T4Tddd, T4Tddd, T4Tddd>;
using T12T3Tddd = std::tuple<T3Tddd, T3Tddd, T3Tddd, T3Tddd, T3Tddd, T3Tddd, T3Tddd, T3Tddd, T3Tddd, T3Tddd, T3Tddd, T3Tddd>;
using T12T2Tddd = std::tuple<T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd, T2Tddd>;
//
using T2b = std::tuple<bool, bool>;
using T3b = std::tuple<bool, bool, bool>;
using T4b = std::tuple<bool, bool, bool, bool>;
using T5b = std::tuple<bool, bool, bool, bool, bool>;
using T6b = std::tuple<bool, bool, bool, bool, bool, bool>;
using T7b = std::tuple<bool, bool, bool, bool, bool, bool, bool>;
using T8b = std::tuple<bool, bool, bool, bool, bool, bool, bool, bool>;

/* -------------------------------------------------------------------------- */
/*                                    array                                   */
/* -------------------------------------------------------------------------- */

template <size_t M = 0, size_t N, typename T, typename Func>
void for_each(std::array<T, N> &arr, const Func &func) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      func(std::get<M>(arr));
      for_each<M + 1>(arr, func);
   }
}

template <size_t M = 0, size_t N, typename T, typename Func>
void for_each(const std::array<T, N> &arr, const Func &func) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      func(std::get<M>(arr));
      for_each<M + 1>(arr, func);
   }
}

template <size_t M = 0, size_t N, typename T, typename Func>
void for_each(std::array<T, N> &arr, std::array<T, N> &ARR, const Func &func) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      func(std::get<M>(arr), std::get<M>(ARR));
      for_each<M + 1>(arr, ARR, func);
   }
}

template <size_t M = 0, size_t N, typename T, typename Func>
void for_each(const std::array<T, N> &arr, const std::array<T, N> &ARR, const Func &func) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      func(std::get<M>(arr), std::get<M>(ARR));
      for_each<M + 1>(arr, ARR, func);
   }
}

template <size_t M = 0, size_t N, typename T, typename Func>
void for_each(const std::array<T, N> &arr0, const std::array<T, N> &arr1, const std::array<T, N> &arr2, const Func &func) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      func(std::get<M>(arr0), std::get<M>(arr1), std::get<M>(arr2));
      for_each<M + 1>(arr0, arr1, arr2, func);
   }
}

template <size_t M = 0, size_t N, typename T, typename Func>
void for_each(std::array<T, N> &arr0, std::array<T, N> &arr1, std::array<T, N> &arr2, const Func &func) {
   if constexpr (M < std::tuple_size<std::array<T, N>>::value) {
      func(std::get<M>(arr0), std::get<M>(arr1), std::get<M>(arr2));
      for_each<M + 1>(arr0, arr1, arr2, func);
   }
}

/* ------------------------------------------------------ */
/*                     タプルのベクトル演算                  */
/* ------------------------------------------------------ */
template <size_t N = 0, typename T>
bool all_of(const T &t) {
   if constexpr (N < std::tuple_size<T>::value) {
      return std::get<N>(t) && all_of<N + 1>(t);
   } else
      return true;
}
template <size_t N = 0, typename T, typename U>
bool all_of(const T &t, const U &func) {
   if constexpr (N < std::tuple_size<T>::value) {
      return func(std::get<N>(t)) && all_of<N + 1>(t, func);
   } else
      return true;
}
template <size_t N = 0, typename T, typename U>
bool any_of(const T &t, const U &func) {
   if constexpr (N < std::tuple_size<T>::value) {
      return func(std::get<N>(t)) || any_of<N + 1>(t, func);
   } else
      return false;
}
template <size_t N = 0, typename T, typename U>
bool none_of(const T &t, const U &func) {
   if constexpr (N < std::tuple_size<T>::value) {
      return !func(std::get<N>(t)) && none_of<N + 1>(t, func);
   } else
      return true;
}
//
template <size_t N = 0, typename T, typename U>
void for_each(T &t, const U &func) {
   if constexpr (1 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
   } else if constexpr (2 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
   } else if constexpr (3 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
   } else if constexpr (4 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
      func(std::get<3>(t));
   } else if constexpr (5 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
      func(std::get<3>(t));
      func(std::get<4>(t));
   } else if constexpr (N < std::tuple_size<T>::value) {
      func(std::get<N>(t));
      for_each<N + 1>(t, func);
   }
}
template <size_t N = 0, typename T, typename U>
void for_each(const T &t, const U &func) {
   if constexpr (1 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
   } else if constexpr (2 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
   } else if constexpr (3 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
   } else if constexpr (4 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
      func(std::get<3>(t));
   } else if constexpr (5 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
      func(std::get<3>(t));
      func(std::get<4>(t));
   } else if constexpr (N < std::tuple_size<T>::value) {
      func(std::get<N>(t));
      for_each<N + 1>(t, func);
   }
}

template <size_t N = 0, typename T, typename U>
void for_each1(T &t, const U &func) {
   if constexpr (1 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
   } else if constexpr (2 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
   } else if constexpr (3 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
   } else if constexpr (4 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
      func(std::get<3>(t));
   } else if constexpr (5 == std::tuple_size<T>::value) {
      func(std::get<0>(t));
      func(std::get<1>(t));
      func(std::get<2>(t));
      func(std::get<3>(t));
      func(std::get<4>(t));
   } else if constexpr (N < std::tuple_size<T>::value) {
      func(std::get<N>(t));
      for_each1<N + 1>(t, func);
   }
}
/* -------------------------------------------------------- */

template <size_t N = 0, typename T, typename H, typename U>
void for_each(T &t, H &h, const U &func) {
   if constexpr (N < std::tuple_size<T>::value && N < std::tuple_size<H>::value) {
      func(std::get<N>(t), std::get<N>(h));
      for_each<N + 1>(t, h, func);
   }
}

template <size_t N = 0, typename T, typename H, typename U>
void for_each01(T t, H &h, const U &func) {
   if constexpr (N < std::tuple_size<T>::value && N < std::tuple_size<H>::value) {
      func(std::get<N>(t), std::get<N>(h));
      for_each01<N + 1>(t, h, func);
   }
}

template <size_t N = 0, typename T, typename H, typename U>
void for_each10(T &t, const H &h, const U &func) {
   if constexpr (N < std::tuple_size<T>::value && N < std::tuple_size<H>::value) {
      func(std::get<N>(t), std::get<N>(h));
      for_each10<N + 1>(t, h, func);
   }
}

// template <size_t N = 0, typename T, typename H, typename U>
// void for_each(T t, H &h, const U &func) {
//    if constexpr (N < std::tuple_size<T>::value && N < std::tuple_size<H>::value) {
//       func(std::get<N>(t), std::get<N>(h));
//       for_each01<N + 1>(t, h, func);
//    }
// }
/* -------------------------------------------------------------------------- */
template <size_t N = 0, typename T0, typename T1, typename T2, typename U>
void for_each(T0 &t0, T1 &t1, T2 &t2, const U &func) {
   if constexpr (N < std::tuple_size<T0>::value &&
                 N < std::tuple_size<T1>::value &&
                 N < std::tuple_size<T2>::value) {
      func(std::get<N>(t0), std::get<N>(t1), std::get<N>(t2));
      for_each<N + 1>(t0, t1, t2, func);
   }
}
template <size_t N = 0, typename T0, typename T1, typename T2, typename U>
void for_each000(const T0 &t0, const T1 &t1, const T2 &t2, const U &func) {
   if constexpr (N < std::tuple_size<T0>::value &&
                 N < std::tuple_size<T1>::value &&
                 N < std::tuple_size<T2>::value) {
      func(std::get<N>(t0), std::get<N>(t1), std::get<N>(t2));
      for_each000<N + 1>(t0, t1, t2, func);
   }
}
template <size_t N = 0, typename T0, typename T1, typename T2, typename U>
void for_each011(const T0 &t0, T1 &t1, T2 &t2, const U &func) {
   if constexpr (N < std::tuple_size<T0>::value &&
                 N < std::tuple_size<T1>::value &&
                 N < std::tuple_size<T2>::value) {
      func(std::get<N>(t0), std::get<N>(t1), std::get<N>(t2));
      for_each011<N + 1>(t0, t1, t2, func);
   }
}
/* -------------------------------------------------------------------------- */
template <size_t N = 0, typename T0, typename T1, typename T2, typename T3, typename U>
void for_each(T0 &t0, T1 &t1, T2 &t2, T3 &t3, const U &func) {
   if constexpr (N < std::tuple_size<T0>::value &&
                 N < std::tuple_size<T1>::value &&
                 N < std::tuple_size<T2>::value &&
                 N < std::tuple_size<T3>::value) {
      func(std::get<N>(t0), std::get<N>(t1), std::get<N>(t2), std::get<N>(t3));
      for_each<N + 1>(t0, t1, t2, t3, func);
   }
}
template <size_t N = 0, typename T0, typename T1, typename T2, typename T3, typename U>
void for_each0111(const T0 &t0, T1 &t1, T2 &t2, T3 &t3, const U &func) {
   if constexpr (N < std::tuple_size<T0>::value &&
                 N < std::tuple_size<T1>::value &&
                 N < std::tuple_size<T2>::value &&
                 N < std::tuple_size<T3>::value) {
      func(std::get<N>(t0), std::get<N>(t1), std::get<N>(t2), std::get<N>(t3));
      for_each0111<N + 1>(t0, t1, t2, t3, func);
   }
}
/* -------------------------------------------------------------------------- */
template <size_t N = 0, typename T0, typename T1, typename T2, typename T3, typename T4, typename U>
void for_each(T0 &t0, T1 &t1, T2 &t2, T3 &t3, T4 &t4, const U &func) {
   if constexpr (N < std::tuple_size<T0>::value &&
                 N < std::tuple_size<T1>::value &&
                 N < std::tuple_size<T2>::value &&
                 N < std::tuple_size<T3>::value &&
                 N < std::tuple_size<T4>::value) {
      func(std::get<N>(t0), std::get<N>(t1), std::get<N>(t2), std::get<N>(t3), std::get<N>(t4));
      for_each<N + 1>(t0, t1, t2, t3, t4, func);
   }
}
template <size_t N = 0, typename T0, typename T1, typename T2, typename T3, typename T4, typename U>
void for_each01111(const T0 &t0, T1 &t1, T2 &t2, T3 &t3, T4 &t4, const U &func) {
   if constexpr (N < std::tuple_size<T0>::value &&
                 N < std::tuple_size<T1>::value &&
                 N < std::tuple_size<T2>::value &&
                 N < std::tuple_size<T3>::value &&
                 N < std::tuple_size<T4>::value) {
      func(std::get<N>(t0), std::get<N>(t1), std::get<N>(t2), std::get<N>(t3), std::get<N>(t4));
      for_each01111<N + 1>(t0, t1, t2, t3, t4, func);
   }
}
/* -------------------------------------------------------------------------- */
// example from https://stackoverflow.com/questions/38885406/produce-stdtuple-of-same-type-in-compile-time-given-its-length-by-a-template-a
// typenameはコンパイラーに，以下が型であることを伝えるためのもの．
template <size_t I, typename T>
struct tuple_n {
   template <typename... Args>
   using type = typename tuple_n<I - 1, T>::template type<T, Args...>;
   // Iを減らすと同時に，引数を増やしていき，Iがゼロとなった場合，最終的に詰め込まれた引数を持つタプルとなる
};

template <typename T>
struct tuple_n<0, T> {
   template <typename... Args>
   using type = std::tuple<Args...>;
};

template <size_t I, typename T>
using tuple_of = typename tuple_n<I, T>::template type<>;

template <size_t I, typename T>
using tuple_tuple_of = typename tuple_n<I, tuple_of<I, T>>::template type<>;

/*
tuple_of<3>
なら
tuple_n<3>::template type<>
が呼ばれ，さらに
tuple_n<3 - 1>::template type<double>
が呼ばれ，さらに
tuple_n<2 - 1>::template type<double,double>
が呼ばれ，さらに
tuple_n<1-1>::template type<double,double,double>
が呼ばれ，これは，
tuple_n<0>::template std::tuple<double,double,double>
となる
*/

/* -------------------------------------------------------------------------- */
#endif