/*DOC_EXTRACT 0_0_0_CRS

\insert{compressed_row_storage}

### CRSの使用例

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test3_CRS.cpp
make
./test3_CRS
```

ここには，`A`かける`V`をCRSを使って高速に計算する例を示している．

*/

#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

using V_d = std::vector<double>;

const int s = 50000;
VV_d A(s, V_d(s, 0.));

struct MatrixCRS {
   std::vector<CRS*> vec_CRS;
   ~MatrixCRS() {
      for (auto& crs : vec_CRS)
         if (crs != nullptr)
            delete crs;
   }

   MatrixCRS(const VV_d& A, const V_d& V) {
      vec_CRS.resize(A.size(), nullptr);
      for (size_t i = 0; i < V.size(); ++i) {
         auto p = new CRS();
         vec_CRS[i] = p;
         p->setIndexCRS(i);
         p->value = V[i];
      }
      for (size_t i = 0; i < A.size(); ++i)
         for (size_t j = 0; j < A.size(); ++j) {
            vec_CRS[i]->set(vec_CRS[j], A[i][j]);
         }
      for (auto& crs : vec_CRS)
         crs->setVectorCRS();
   }

   V_d dot() {
      return selfDot(vec_CRS);
   };
};

int main() {

   // make A
   int n = 1000;
#pragma omp parallel
   for (size_t i = 0; i < A.size(); ++i)
#pragma omp single nowait
      for (size_t j = 0; j < n && j < A.size(); ++j) {
         A[i][j] = 1.;
      }

   V_d V;
   for (size_t i = 0; i < A.size(); ++i)
      V.emplace_back(i + 1);

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
   このインデックスは，`std::vector<VRS*>`と掛け算をする相手である`V`の行番号に対等させておく必要がある．

   **掛け算`Dot(A,V)`において，CRS（これは行ベクトルと考える）は，自分に保存されている{row index,value}のセットを元に，`V[row index]*value`のようにして足し合わせていく．**

  */

   MatrixCRS Mat_CRS(A, V);

   /*DOC_EXTRACT 0_0_3_CRS

   #### 値を格納：`set`と`increment`

   値を格納するには，２つの方法があり，`set`と`increment`がある．
   このように，インデックスと値を指定するようにしているのは，
   値がゼロの場合は，何もせず，インデックスも保存しないようにしているためである．

   値を設定する，`set`と`increment`の第一引数は，CRSのポインタである．

  */

   TimeWatch watch;
   watch();

   std::cout << Total(ParallelDot(A, V)) << std::endl;
   std::cout << Yellow << " Original" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

   watch();
   std::cout << Total(Dot(Mat_CRS.vec_CRS, V)) << std::endl;
   std::cout << Yellow << " CRS" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

   /*DOC_EXTRACT 0_0_3_CRS

   #### `selfDot`

   `selfDot`は，CRSに保存した`A`と`V`を掛け合わせる関数である．

  */

   std::cout << Total(Mat_CRS.dot()) << std::endl;
   std::cout << Yellow << " CRS" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

   return 0;
}
