/*
メタプログラミングを十分に理解するために，
https://riptutorial.com/cplusplus/topic/462/metaprogramming
を全て実行する．
*/

#define GenericMinMax
/* ------------------------------------------------------ */
#if defined(calculating_factorials)

#define No1

#if defined(No1)

// https://riptutorial.com/cplusplus/example/1525/calculating-factorials
#include <iostream>
template <unsigned int n> //コンパイル時にnを使いまわせることに驚き
struct factorial
{
    // enum // enumは必要なのか？
    // {
    //     value = n * factorial<n - 1>::value // nが使い回されている．再起的
    // };

    static const int A = n * factorial<n - 1>::A; // この方法でもOK
    static const int B = n * factorial<n - 1>::B; // この方法でもOK
};

template <>
struct factorial<0>
{
    // enum
    // {
    //     value = 1
    // };

    static const int A = 1; // この方法でもOK
    static const int B = 1; // この方法でもOK
};

int main()
{
    std::cout << factorial<7>::A << std::endl;
    std::cout << factorial<7>::B << std::endl;
}

#elif defined(No2)

// constexprを使った例
#include <iostream>
constexpr long long factorial(long long n)
{
    return (n == 0) ? 1 : n * factorial(n - 1);
}

int main()
{
    char test[factorial(3)];
    std::cout << factorial(7) << std::endl;
}

#elif defined(selfmade)

#include <iostream>
constexpr long long factorial(long long n)
{
    return (n == 0) ? 1 : n * factorial(n - 1);
}

constexpr long long combination(long long n, long long r)
{
    return factorial(n) / (factorial(r) * factorial(n - r));
}

int main()
{
    std::cout << combination(7, 3) << std::endl;
}

#endif
/* ------------------------------------------------------ */
#elif defined(calculating_power)

#include <functional>
#include <type_traits>
#include <utility>

template <class, class = std::void_t<>>
struct has_hash
    : std::false_type
{
};

template <class T>
struct has_hash<T, std::void_t<decltype(std::hash<T>()(std::declval<T>()))>>
    : std::true_type
{
};

int main(){};
/* ------------------------------------------------------ */
#elif defined(GenericMinMax)
#include <iostream>
#include <tuple>
template <typename T1, typename T2>
auto min(const T1 &a, const T2 &b) -> typename std::common_type<const T1 &, const T2 &>::type
{
    return a < b ? a : b;
}

template <typename T1, typename T2, typename... Args>
auto min(const T1 &a, const T2 &b, const Args &...args) -> typename std::common_type<const T1 &, const T2 &, const Args &...>::type
{
    return min(min(a, b), args...);
}

template <typename T>
auto min(const std::tuple<T, T> &a)
{
    static const int N = std::tuple_size(a);
    return std::get<0>(a) < std::get<1>(a) ? std::get<0>(a) : std::get<1>(a);
}

int main()
{
    std::cout << min(1, 2) << std::endl;
}

#endif