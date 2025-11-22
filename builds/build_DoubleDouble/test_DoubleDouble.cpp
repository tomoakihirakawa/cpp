#include <quadmath.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdfloat>
#include <string>
#include <vector>

/*DOC_EXTRACT 0_error

DoubleDoubleを作成する前に，数値がどのタイミングで丸められるかを知っておく必要がある．

```cpp
sh clean
cmake -DCMAKE_BUILD_TYPE=Debug ../ -DSOURCE_FILE=test_DoubleDouble.cpp
make
./test_DoubleDouble
```

この実装は間違いがある．

```cpp
DoubleDouble operator+(const DoubleDouble& other) const {
    DoubleDouble result;
    result.a = a + other.a;//!ここで丸められてしまい誤差が生じる
    result.b = b + other.b;
    return result;
}
```

なので，このdouble同士の演算においても誤差を取りこぼさない工夫が必要となる．

*/

struct DoubleDouble {
   double a, b;

   DoubleDouble() : a(0.0), b(0.0) {}

   DoubleDouble(const std::float128_t& xy) {
      b = static_cast<double>(xy - (a = static_cast<double>(xy)));
   }

   // Example method to add two DoubleDouble numbers
   DoubleDouble operator+(const DoubleDouble& other) const {
      DoubleDouble result;
      result.a = a + other.a;  //! ここで丸められてしまい誤差が生じる
      result.b = b + other.b;
      return result;
   }

   // Example method to add two DoubleDouble numbers
   DoubleDouble operator-(const DoubleDouble& other) const {
      DoubleDouble result;
      result.a = a - other.a;
      result.b = b - other.b;
      return result;
   }

   // std::float128_t conversion operator
   operator std::float128_t() const {
      return std::float128_t(a) + std::float128_t(b);
   }
};

template <typename STREAM>
auto operator<<(STREAM& stream, const std::same_as<std::float128_t> auto& q) -> STREAM& {
   char buf[128];
   quadmath_snprintf(buf, sizeof(buf), "%.50Qg", q);
   stream << buf;
   return stream;
}

int main() {
   // Example of using double and double precision
   double pi = M_PI;
   double PI = 3.1415926535897932384626433832795028841971693993751;
   double PI2 = 6.2831853071795864769252867665590057683943387987502;
   DoubleDouble dd_X(3.14159265358979323846264338328q);

   std::cout.precision(50);  // Set precision for output
   auto w = std::setw(60);
   std::cout << w << "" << "3.1415926535897932384626433832795028841971693993751" << std::endl;
   std::cout << w << "Double precision pi: " << pi << std::endl;
   //! 元々doubleなので，演算はdouble精度で当然行われる
   std::cout << w << "dd_X.a + dd_X.b =  " << dd_X.a + dd_X.b << std::endl;
   //! 以下は，doubleの演算の後に，float128_tへと変換するので，精度は変わらない．double以上はIEEE754の丸め規則に従う
   std::cout << w << "std::float128_t(dd_X.a + dd_X.b) =  " << std::float128_t(dd_X.a + dd_X.b) << std::endl;
   //! float128_tへと変換され，演算はfloat128精度で行われる
   std::cout << w << "std::float128_t(dd_X.a) + std::float128_t(dd_X.b) =  " << std::float128_t(dd_X.a) + std::float128_t(dd_X.b) << std::endl;
   //! doubleの精度より小さい桁は，IEEE754の丸め規則に従って決まり，ランダムではない
   std::cout << w << "std::float128_t(dd_X.a) =  " << std::float128_t(dd_X.a) << std::endl;
   //! なので，その決まった方法で残った桁と実際の値の差をbとして保持し，後に必要なときに足し合わせる
   std::cout << w << "std::float128_t(dd_X.b) =  " << std::float128_t(dd_X.b) << std::endl;

   std::cout << "=========================== 足し算でも精度を維持するか ==============================" << std::endl;

   std::cout << w << "" << "6.2831853071795864769252867665590057683943387987502" << std::endl;
   std::cout << w << "(dd_X + dd_X).a =  " << (dd_X + dd_X).a << std::endl;
   std::cout << w << "(dd_X + dd_X).b =  " << (dd_X + dd_X).b << std::endl;
   std::cout << w << "std::float128_t(dd_X + dd_X) =  " << std::float128_t(dd_X + dd_X) << std::endl;

   std::cout << "=========================== 引き算でも精度を維持するか ===============================" << std::endl;
   DoubleDouble dd_Y = dd_X + dd_X;
   std::cout << w << "" << "3.1415926535897932384626433832795028841971693993751" << std::endl;
   std::cout << w << "(dd_X - dd_X).a =  " << (dd_Y - dd_X).a << std::endl;
   std::cout << w << "(dd_X - dd_X).b =  " << (dd_Y - dd_X).b << std::endl;
   std::cout << w << "std::float128_t(dd_X - dd_X) =  " << std::float128_t(dd_Y - dd_X) << std::endl;

   std::cout << "ただし，この実装は間違いがある．" << std::endl;
}
