#include "lib_quadmath.hpp"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <quadmath.h>
#include <stdfloat>
#include <string>
#include <vector>

/*DOC_EXTRACT 0_error

DoubleDoubleを作成する前に，数値がどのタイミングで丸められるかを知っておく必要がある．

```cpp
sh clean
cmake -DCMAKE_BUILD_TYPE=Debug ../ -DSOURCE_FILE=test_DoubleDouble_correct.cpp
make
./test_DoubleDouble_correct
```

*/

int main() {
  // Example of using double and double precision
  double pi = M_PI;
  std::float128_t PI = 3.14159265358979323846264338327950288419716939937510582097494q;
  std::float128_t PI2 = 6.28318530717958647692528676655900576839433879875021164194989q;
  DoubleDouble PI_dd(PI);
  DoubleDouble PI2_dd(PI2);

  std::cout.precision(50); // Set precision for output
  auto w = std::setw(70);
  std::cout << "=========================== 丸めに関する基礎知識 ==============================" << std::endl;
  std::cout << w << "True pi: " << "3.1415926535897932384626433832795028841971693993751" << std::endl;
  std::cout << w << "Double precision pi: " << double(PI) << std::endl;
  std::cout << "倍精度(double)では，この時点で既に丸められている．" << std::endl;
  std::cout << "ただし，丸めはIEEE754の規則に従って行われるため，丸められた値は決まっている．" << std::endl;
  //! 元々doubleなので，演算はdouble精度で当然行われる
  std::cout << "=========================== DoubleDoubleの内部表現について ==============================" << std::endl;
  std::cout << w << "DoubleDoubleは，highとlowの2つのdoubleで表現される:" << std::endl;
  std::cout << w << "PI_dd.a (high) =  " << PI_dd.a << std::endl;
  std::cout << w << "PI_dd.b (low)  =  " << PI_dd.b << std::endl;
  std::cout << w << "highとlowをfloat128_tに変換して足し合わせれば，元の値と同じになる:" << std::endl;
  std::cout << w << "std::float128_t(PI_dd.a) + std::float128_t(PI_dd.b) =  " << std::float128_t(PI_dd.a) + std::float128_t(PI_dd.b) << std::endl;
  std::cout << "=========================== 足し算で精度を維持するか ==============================" << std::endl;

  std::cout << w << "True PI+PI: " << "6.2831853071795864769252867665590057683943387987502" << std::endl;
  std::cout << w << "倍精度 double(PI) + double(PI) = " << static_cast<std::float128_t>(double(PI) + double(PI)) << std::endl;
  std::cout << w << "倍倍精度 PI_dd + PI_dd = " << static_cast<std::float128_t>(PI_dd + PI_dd) << std::endl;
  std::cout << w << "std::float128_t(PI + PI) = " << static_cast<std::float128_t>(PI + PI) << std::endl;

  std::cout << "=========================== 引き算でも精度を維持するか ===============================" << std::endl;
  std::cout << w << "True 2PI-PI: " << "3.1415926535897932384626433832795028841971693993751" << std::endl;
  std::cout << w << "倍精度 double(PI2) - double(PI) = " << static_cast<std::float128_t>(double(PI2) - double(PI)) << std::endl;
  std::cout << w << "倍倍精度 PI2_dd - PI_dd = " << static_cast<std::float128_t>(PI2_dd - PI_dd) << std::endl;
  std::cout << w << "std::float128_t(PI2 - PI) = " << static_cast<std::float128_t>(PI2 - PI) << std::endl;

  /* -------------------------------------------------------------------------- */
  //
  std::cout << w << "" << "-9.86960440108935861883449099987615113531" << std::endl;
  std::complex<DoubleDouble> dd_complex(DoubleDouble(PI), (DoubleDouble)(PI));
  std::complex<DoubleDouble> dd2_complex((DoubleDouble)(PI), (DoubleDouble)(PI2));
  std::cout << w << "dd_complex: " << (std::float128_t)(dd_complex).real() << std::endl;
  std::cout << w << "dd_complex*dd2_complex: " << (std::float128_t)(dd_complex * dd2_complex).real() << std::endl;
}
