/*DOC_EXTRACT

## Compressed Sparse Row (CSR)

CSRは行列を表現する方法の一つである．
このCSRクラスは，std::unordered_mapを用いて，行列の非ゼロ要素を表現する．
std::unordered_mapのkeyはポインタであり，valueはdoubleである．
CSRクラス自身が，行列の行番号を保存しており，keyであるCSRクラスは行列の列番号を保存している．
*/

#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

int main() {

   /* -------------------------------------------------------------------------- */
   const VV_d A = {{4., -1, 0, -1},
                   {-1., 4, -1, 0},
                   {0., -1, 4, -1},
                   {-1., 0, -1, 4}};
   const V_d b = {15., 10., 10, 15};
   // 初期値が大事
   V_d x0(b.size(), 0.);
   V_d ans = {6.875, 5.625, 5.625, 6.875};
   /* -------------------------------------------------------------------------- */

   std::unordered_set<CSR*> V_CSR;

   {
      auto csr = new CSR();
      csr->setIndexCSR(0);
      csr->value = 15.;
      V_CSR.emplace(csr);
   }
   {
      auto csr = new CSR();
      csr->setIndexCSR(1);
      csr->value = 10.;
      V_CSR.emplace(csr);
   }
   {
      auto csr = new CSR();
      csr->setIndexCSR(2);
      csr->value = 10.;
      V_CSR.emplace(csr);
   }
   {
      auto csr = new CSR();
      csr->setIndexCSR(3);
      csr->value = 15.;
      V_CSR.emplace(csr);
   }
   for (const auto& crs0 : V_CSR)
      for (const auto& crs1 : V_CSR) {
         int i = crs0->getIndexCSR();
         int j = crs1->getIndexCSR();
         if (A[i][j] != static_cast<double>(0.))
            crs0->increment(crs1, A[i][j]);
      }

   std::cout << " Dot : " << Dot(V_CSR, V_CSR) << std::endl;
   std::cout << " Dot : " << Dot(V_CSR) << std::endl;
   std::cout << " Dot : " << Dot(A, b) << std::endl;

   Timer timer;
   std::cout << "time:" << timer() << std::endl;
   bool finished = false;
   double error;
   for (auto restart = 0; restart < 1; ++restart) {
      for (auto i = 1; i <= 10; ++i) {
         gmres gm(V_CSR, b, x0, i);
         //   std::cout << "gm.x = " << gm.x << std::endl;
         std::cout << "time:" << timer() << std::endl;
         error = gm.err;  // Norm(Dot(A, gm.x) - b);
         x0 = gm.x;
         std::cout << "       Restart : " << restart << std::endl;
         std::cout << "             i : " << i << std::endl;
         std::cout << "estimate error : " << gm.err << std::endl;
         std::cout << "   actual error: " << Norm(Dot(A, x0) - b) << std::endl;
         std::cout << "      Solution : " << x0 << std::endl;
         std::cout << Red << "--------------------------------" << colorOff << std::endl;
      }
   }
   std::cout << "time:" << timer() << std::endl;
};
