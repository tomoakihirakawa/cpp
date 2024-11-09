/*DOC_EXTRACT 0_0_0_CRS

\insert{compressed_row_storage}

### CRSの使用例

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test3_CRS.cpp
make
./test3_CRS
```

ここには，`A`かける`V`をCRSを使って高速に計算する例を示している．

*/

/*DOC_EXTRACT 0_0_1_CRS

#### CRSは，ある行ベクトルを格納するクラスと考える

私のプログラムでは，Row-major orderで行列を格納しており，次のように行列を定義している．

```cpp
std::vector<std::vector<double>> Mat; // <- std::vector<ROW VECTOR> Mat;
```

```math
\begin{pmatrix}
\{a_{11} & a_{12} & a_{13} & \cdots & a_{1n}\}&\leftarrow {\text{a ROW VECTOR}} \\
\{a_{21} & a_{22} & a_{23} & \cdots & a_{2n}\}&\\
\{a_{31} & a_{32} & a_{33} & \cdots & a_{3n}\}&\\
\{\vdots & \vdots & \vdots & \ddots & \vdots \}&\\
\{a_{m1} & a_{m2} & a_{m3} & \cdots & a_{mn}\}&
\end{pmatrix}
```

CRSは，このROW VECTORを格納するクラスであり，CRSのベクトルが行列となる．

```cpp
std::vector<CRS*> Mat_CRS(A.size());
```
*/

/*DOC_EXTRACT 0_0_2_CRS

#### `setIndexCRS`

CRSは，`CRS->setIndexCRS(i)`のようにして，自身の行番号を保持しておく．
このインデックスは，`std::vector<CRS*>`と掛け算をする相手である`V`の行番号に対等させておく必要がある．

**掛け算`Dot(A,V)`において，CRS（これは行ベクトルと考える）は，自分に保存されている{row index,value}のセットを元に，`V[row index]*value`のようにして足し合わせていく．**

*/

/*DOC_EXTRACT 0_0_3_CRS

#### 値を格納：`set`と`increment`

値を格納するには，２つの方法があり，`set`と`increment`がある．
このように，インデックスと値を指定するようにしているのは，
値がゼロの場合は，何もせず，インデックスも保存しないようにしているためである．

値を設定する，`set`と`increment`の第一引数は，CRSのポインタである．

*/

/*DOC_EXTRACT 0_0_3_CRS

#### `selfDot`

`selfDot`は，CRSに保存した`A`と`V`を掛け合わせる関数である．

*/

#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

using V_d = std::vector<double>;

/* -------------------------------------------------------------------------- */

int main() {

   std::vector<std::array<double, 5>> n_dot_pardot_crsdot_crsseldot;

   const int s = 100000;
   VV_d A(s, V_d(s, 0.));
   V_d V(A.size(), 1.);

   std::vector<CRS*> A_CRS(A.size());

   std::array<double, 2> time_dot, time_pardot, time_crsdot, time_crsseldot;

#pragma omp parallel for
   for (size_t i = 0; i < A.size(); ++i) {
      A_CRS[i] = new CRS();
      A_CRS[i]->setIndexCRS(i);
      A_CRS[i]->value = V[i];
   }

   for (auto n = 2000; n < 20000; n += 2000) {

      std::cout << "n = " << n << std::endl;
      std::cout << "make A" << std::endl;
#pragma omp parallel
      for (auto& a : A)
#pragma omp single nowait
      {
         auto i = 0;
         for (auto& ai : a) {
            if (i < n)
               ai = 1.;
            else
               ai = 0;
            i++;
         }
      }

      /* ------------------------------------------------------------ */
      /*                          CRS行列の作成                         */
      /* ------------------------------------------------------------ */

      std::cout << "make CRS" << std::endl;
#pragma omp parallel for
      for (size_t i = 0; i < A.size(); ++i) {
         A_CRS[i]->setIndexCRS(i);
         A_CRS[i]->value = V[i];
      }

      std::cout << "set" << std::endl;
#pragma omp parallel for
      for (size_t i = 0; i < A.size(); ++i) {
         for (size_t j = 0; j < n; ++j)
            A_CRS[i]->set(A_CRS[j], A[i][j]);
      }

      std::cout << "setVectorCRS" << std::endl;
#pragma omp parallel for
      for (size_t i = 0; i < A.size(); ++i)
         A_CRS[i]->setVectorCRS();

      /* ------------------------------------------------------------ */

      std::cout << "calculation" << std::endl;

      TimeWatch watch;

      // std::cout << Total(Dot(A, V)) << std::endl;
      // std::cout << "Total(Dot(A, V))" << std::endl;
      // std::cout << Yellow << " Original" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      time_dot = watch();

      std::cout << Total(ParallelDot(A, V)) << std::endl;
      std::cout << "Total(ParallelDot(A, V))" << std::endl;
      std::cout << Yellow << " Original" << Blue << "\nElapsed time: " << Red << (time_pardot = watch()) << colorReset << " s\n";

      std::cout << Total(PreparedDot(A_CRS, V)) << std::endl;
      std::cout << "Total(Dot(A_CRS, V))" << std::endl;
      std::cout << Yellow << " CRS" << Blue << "\nElapsed time: " << Red << (time_crsdot = watch()) << colorReset << " s\n";

      std::cout << Total(Dot(A_CRS)) << std::endl;
      std::cout << "Total(selfDot(A_CRS))" << std::endl;
      std::cout << Yellow << " CRS" << Blue << "\nElapsed time: " << Red << (time_crsseldot = watch()) << colorReset << " s\n";

      n_dot_pardot_crsdot_crsseldot.push_back({n, time_dot[0], time_pardot[0], time_crsdot[0], time_crsseldot[0]});
   }

   for (size_t i = 0; i < A.size(); ++i)
      delete A_CRS[i];

   std::ofstream ofs("n_dot_pardot_crsdot_crsseldot.dat");
   ofs << "n time_dot time_pardot time_crsdot time_crsseldot\n";
   for (const auto [n, time_dot, time_pardot, time_crsdot, time_crsseldot] : n_dot_pardot_crsdot_crsseldot)
      ofs << n << " " << time_dot << " " << time_pardot << " " << time_crsdot << " " << time_crsseldot << "\n";

   ofs.close();
}
