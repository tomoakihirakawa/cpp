#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include "lib_Fourier.hpp"

/*

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_2D_DFT.cpp
make
./example1_2D_DFT
```

*/

int main() {
   //! generate 2d data
   const int N = 20, M = 20;
   std::vector<std::vector<double>> data2D(N, std::vector<double>(M, 0.));

   auto fx = [](double t) {
      const int k = 49;
      return std::cos(k * 2. * M_PI / N * t);
   };
   auto fy = [](double t) {
      const int k = 20;
      return std::sin(k * 2. * M_PI / M * t);
   };

   for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
         data2D[i][j] = fx(i) * fy(j);

   std::cout << "2D data" << std::endl;
   std::ofstream ofs("2D_data.dat");
   for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
         ofs << data2D[i][j] << (M - 1 == j ? "\n" : " ");
   ofs.close();

   /* ----------------- 2D フーリエ変換 --------------- */

   auto cn2D = DFT(data2D);

   std::cout << "2D DFT" << std::endl;
   {
      std::ofstream ofs2("2D_DFT_real.dat");
      for (int i = 0; i < N; ++i)
         for (int j = 0; j < M; ++j)
            ofs2 << cn2D[i][j].real() << (M - 1 == j ? "\n" : " ");
      ofs2.close();
   }
   {
      std::ofstream ofs2("2D_DFT_imag.dat");
      for (int i = 0; i < N; ++i)
         for (int j = 0; j < M; ++j)
            ofs2 << cn2D[i][j].imag() << (M - 1 == j ? "\n" : " ");
      ofs2.close();
   }

   /* -------------- 逆フーリエ変換を行うと元のデータに戻るか確認する． -------------- */

   auto inverse_c = InverseDFT(cn2D);

#pragma omp parallel
   for (auto i = 0; i < 100000; ++i)
#pragma omp single nowait
      InverseDFT(cn2D);  //! 何秒くらいかかるか．チェック

   std::cout << "Inverse DFT" << std::endl;
   {
      std::ofstream ofs2("Inverse_DFT.dat");
      for (int i = 0; i < N; ++i)
         for (int j = 0; j < M; ++j)
            ofs2 << inverse_c[i][j].real() << (M - 1 == j ? "\n" : " ");
      ofs2.close();
   }
}
