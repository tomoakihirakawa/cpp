#include <quadmath.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <stdfloat>
#include <string>
#include <vector>
#include "lib_quadmath.hpp"

/*

```cpp
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_DoubleDouble_correct2.cpp
make
./test_DoubleDouble_correct2
```

*/

// テスト結果をチェックし、表示するためのヘルパー関数
void check(const std::string& test_name, bool passed) {
   std::cout << "Test: " << std::left << std::setw(45) << test_name << " | "
             << (passed ? "\033[32mPASSED\033[0m" : "\033[31mFAILED\033[0m") << std::endl;
}

// DoubleDoubleの結果を__float128の期待値と比較するヘルパー関数
template <typename TYPE>
void check_near(const std::string& test_name, const TYPE& result, const std::float128_t& expected, const std::float128_t& tolerance = 1e-31q) {
   auto diff = fabsq(static_cast<std::float128_t>(result) - expected);
   bool passed = diff < tolerance;
   check(test_name, passed);
   if (!passed) {
      std::cout.precision(50);
      std::cout << "  Result:   " << static_cast<std::float128_t>(result) << std::endl;
      std::cout << "  Expected: " << expected << std::endl;
      std::cout << "  Difference: " << static_cast<std::float128_t>(diff) << std::endl;
   }
}

// --- Test Cases ---

// 1. 基本的な精度とコンストラクタのテスト
void test_precision() {
   std::cout << "\n--- 1. Testing Basic Precision and Constructors ---\n";
   auto PI_q = 3.14159265358979323846264338327950288q;
   DoubleDouble dd_pi(PI_q);
   check_near("Constructor from __float128", dd_pi, PI_q);

   // doubleでは丸められる値をテスト
   double d_val = 1.0 + 1e-17;
   std::float128_t q_val = 1.0q + 1e-17q;
   DoubleDouble dd_val(q_val);
   check("Constructor preserves small values", static_cast<double>(dd_val.a) == 1.0 && dd_val.b != 0.0);
   check_near("Value correctness for small additions", dd_val, q_val);
}

//! 1.5 accumulation errorのテスト
void test_accumulation_error() {
   std::cout << "\n--- 1.5 Testing Accumulation Error ---\n";
   std::float128_t PI_q = 3.14159265358979323846264338327950288q;
   std::float128_t float128_sum;
   LongDoubleLongDouble long_dd_sum;
   DoubleDouble dd_sum;
   int N = 100000;  // 足し合わせる回数
   {
      DoubleDouble dd_pi(PI_q);
      // 1000回足し合わせる
      for (int i = 0; i < N; ++i) {
         dd_sum += dd_pi * dd_pi;
      }
   }
   {
      LongDoubleLongDouble long_dd_pi(PI_q);
      // 1000回足し合わせる
      for (int i = 0; i < N; ++i) {
         // long_dd_sum += long_dd_pi * long_dd_pi;
         long_dd_sum += (PI_q * PI_q);
      }
   }
   {
      // 1000回足し合わせる
      for (int i = 0; i < N; ++i) {
         float128_sum += PI_q * PI_q;
      }
   }
   // 結果を検証
   check_near("Accumulation of pi", dd_sum, N * PI_q * PI_q);
   check_near("Accumulation of pi (LongDoubleLongDouble)", long_dd_sum, N * PI_q * PI_q);
   check_near("Accumulation of pi (__float128)", float128_sum, N * PI_q * PI_q);
   // 丸め誤差の確認
}

// 2. 四則演算のテスト
void test_arithmetic() {
   std::cout << "\n--- 2. Testing Arithmetic Operations ---\n";
   auto PI_q = 3.14159265358979323846264338327950288q;
   auto E_q = 2.71828182845904523536028747135266250q;
   DoubleDouble dd_pi(PI_q);
   DoubleDouble dd_e(E_q);

   // 足し算
   check_near("Addition (pi + e)", dd_pi + dd_e, PI_q + E_q);

   // 引き算
   check_near("Subtraction (pi - e)", dd_pi - dd_e, PI_q - E_q);

   // 掛け算
   check_near("Multiplication (pi * e)", dd_pi * dd_e, PI_q * E_q);

   // 割り算
   std::float128_t E = 2.7182818284590452353602874713526624977572470937000q;
   check_near("Division (pi / e)", dd_pi / E, PI_q / E_q);

   // 複合代入演算子
   DoubleDouble dd_sum = dd_pi;
   dd_sum += dd_e;
   check_near("Compound Addition (+=)", dd_sum, PI_q + E_q);
}

// 3. 情報落ち（Catastrophic Cancellation）のテスト
void test_cancellation() {
   std::cout << "\n--- 3. Testing Catastrophic Cancellation Avoidance ---\n";
   // ほぼ等しい2つの値を用意
   auto val1_q = 1.0q;
   auto val2_q = 1.0q + 1e-25q;
   DoubleDouble dd1(val1_q);
   DoubleDouble dd2(val2_q);

   // doubleでは情報落ちして0になるか、精度が大幅に悪化する
   double d1 = 1.0;
   double d2 = 1.0 + 1e-16;  // 1e-25はdoubleでは表現不可
   check("Standard double loses precision", (d2 - d1) != 1e-16);

   // DoubleDoubleでは正しく計算できるはず
   DoubleDouble diff = dd2 - dd1;
   check_near("Subtraction of nearly-equal numbers", diff, 1e-25q);
}

// 4. 比較演算子のテスト
void test_comparison() {
   std::cout << "\n--- 4. Testing Comparison Operators ---\n";
   DoubleDouble dd1(1.0);
   DoubleDouble dd2(1.0);
   DoubleDouble dd3(1.0 + 1e-30);

   check("Equality (==)", dd1 == dd2);
   check("Inequality (!=)", dd1 != dd3);
   check("Less than (<)", dd1 < dd3);
   check("Greater than (>)", dd3 > dd1);
   check("Less than or equal (<=)", dd1 <= dd2);
   check("Greater than or equal (>=)", dd2 >= dd1);

   // doubleとの比較
   // check("Comparison with double (==)", dd1 == 1.0);
   check("Comparison with double (<)", dd1 < (1.0 + 1e-20));
}

// 5. constexprのテスト
void test_constexpr() {
   std::cout << "\n--- 5. Testing Constexpr Evaluation ---\n";
   // コンパイル時に計算を実行
   constexpr DoubleDouble a(1.23);
   std::cout << a.a << std::endl;
   std::cout << a.b << std::endl;
   constexpr DoubleDouble b(4.56);
   std::cout << a.a << std::endl;
   std::cout << a.b << std::endl;
   constexpr DoubleDouble sum = a + b;
   constexpr DoubleDouble prod = a * b;

   // 結果を検証
   check_near("Constexpr addition", sum, 1.23q + 4.56q);
   check_near("Constexpr multiplication", prod, 1.23q * 4.56q);

   // static_assertでコンパイル時にチェックすることも可能
   // static_assert((a + b) == DoubleDouble(5.79));
   // check("static_assert for constexpr", true);
}

// 6. std::complexとの連携テスト
void test_complex() {
   std::cout << "\n--- 6. Testing std::complex Integration ---\n";
   auto PI_q = 3.14159265358979323846264338327950288q;
   DoubleDouble dd_pi(PI_q);
   DoubleDouble dd_one(1.0);

   std::complex<DoubleDouble> c1(dd_pi, dd_one);  // pi + 1i
   std::complex<DoubleDouble> c2(dd_one, dd_pi);  // 1 + pi*i

   std::complex<DoubleDouble> c_sum = c1 + c2;
   std::complex<std::float128_t> q_sum = {PI_q + 1.0q, 1.0q + PI_q};

   check_near("Complex addition (real part)", c_sum.real(), q_sum.real());
   check_near("Complex addition (imag part)", c_sum.imag(), q_sum.imag());

   // (pi + i) * (1 + pi*i) = (pi - pi) + (pi^2 + 1)i = (pi^2 + 1)i
   std::complex<DoubleDouble> c_prod = c1 * c2;
   std::complex<std::float128_t> q_prod = {0.0q, PI_q * PI_q + 1.0q};

   check_near("Complex multiplication (real part)", c_prod.real(), q_prod.real());
   check_near("Complex multiplication (imag part)", c_prod.imag(), q_prod.imag());
}

int main() {
   std::cout << std::fixed;  // 固定小数点表記に設定
   std::cout.precision(10);  // デフォルトの表示精度

   std::cout << "=================================================" << std::endl;
   std::cout << "      Running Tests for DoubleDouble Class       " << std::endl;
   std::cout << "=================================================" << std::endl;

   test_precision();
   test_accumulation_error();
   test_arithmetic();
   test_cancellation();
   test_comparison();
   test_constexpr();
   test_complex();

   std::cout << "\n-------------------------------------------------" << std::endl;
   std::cout << "                 All tests finished.             " << std::endl;
   std::cout << "-------------------------------------------------" << std::endl;

   return 0;
}
