//https://theolizer.com/cpp-school2/cpp-school2-2/#typename.h
#include <string>
//
#if defined(__GNUC__)
#include <cxxabi.h>
std::string getNameByTypeInfo(std::type_info const &iTypeInfo)
{
    char *aName;
    int status = 0;
    aName = abi::__cxa_demangle(iTypeInfo.name(), 0, 0, &status);
    std::string ret(aName);
    std::free(aName);
    return ret;
}
#else
std::string getNameByTypeInfo(std::type_info const &iTypeInfo)
{
    return iTypeInfo.name();
}
#endif

#define TYPENAME(dType) getNameByTypeInfo(typeid(dType))

///////////////////////////////////////////////////////////////

//可変個引数テンプレート（variadictemplate）
//https://theolizer.com/cpp-school2/cpp-school2-14/#parameter_pack

#define test1
#ifdef test1
#include <iostream>
#include <iomanip>
#include <ctime>

template <typename tType>
struct Type
{
};

//再起定義の最後を処理する通常のオーバーロード関数
void printLogImpl() { std::cout << std::endl; }

//再起定義の途中を処理する関数テンプレート
template <typename tFirst, typename... tRest>
void printLogImpl(tFirst first, tRest... rest)
{
    std::cout << first << "  [" << TYPENAME(Type<decltype(first)>) << "]";
    auto i = sizeof...(rest);
    std::cout << "     sizeof...(rest) = " << i << std::endl;
    printLogImpl(rest...);
}

//最初の入口となる関数テンプレート
template <typename... T>
void printLog(T... t)
{
    std::time_t now = std::time(nullptr);
    std::cout << std::put_time(std::localtime(&now), "%T") << " : ";
    printLogImpl(t...);
}
#endif
#ifdef test2
#include <iostream>
#include <iomanip>
#include <ctime>

template <typename... tTypes>
void printLogImpl(tTypes const &...){};

template <typename... tTypes>
void printLog(tTypes const &... a)
{
    printLogImpl(std::cout << a << std::endl...);
};

#endif

int main()
{
    printLog(1, 1.2, "three", 10, 3, 2);
}
