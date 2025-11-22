/*DOC_EXTRACT 0FourierTransform

## FMMのM2L変換として作成する係数$L_k^j$の計算

対象とする積分：

$$
\begin{aligned}
L_{j}^{k}&=\mathop{\sum}\limits_{{n}={0}}\limits^{N}{\mathop{\sum}\limits_{{m}=-{n}}\limits^{n}{{f}_{n}^{m}{g}_{{n}+{j}}^{{m}-{k}}}}\\
&=\mathop{\sum}\limits_{{n}={0}}\limits^{N}{\mathop{\sum}\limits_{{m}=-{n}+{N}}\limits^{{n}+{N}}{{f}_{n}^{{m}-{N}}{g}_{{n}+{j}}^{-\left({{k}-{m}}\right)-{N}}}}\\
&=\mathop{\sum}\limits_{{n}={0}}\limits^{N}{\mathop{\sum}\limits_{{m}={0}}\limits^{2N}{{f}_{n}^{{m}-{N}}{g}_{{n}+{j}}^{-\left({{k}-{m}}\right)-{N}}}}\\
&={\left\{{{\mathcal{F}}_{2D}^{-{1}}\left[{{\mathcal{F}}_{2D}\left[{{\mathbf{f}}_{p}^{{q}-{N}}}\right]\odot{\mathcal{F}}_{2D}\left[{{\mathbf{g}}_{p}^{-{q}-{N}}}\right]}\right]}\right\}}_{j}^{k}
\end{aligned}
$$

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Debug ../ -DSOURCE_FILE=example1_2D_DFT_complex_convolution_check.cpp
make
./example1_2D_DFT_complex_convolution_check
```

*/

#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>
#include "lib_Fourier.hpp"

const int N = 10;  // 行列のサイズ

std::complex<double> f(const int n, const int m) {
   if (std::abs(m) > n || n < 0 || n > N)
      return std::complex<double>(0, 0);
   else
      return std::complex<double>(n + m + 1, n - m + 1);
}

std::complex<double> g(const int n, const int m) {
   if (std::abs(m) > n || n < 0 || n > N)
      return std::complex<double>(0, 0);
   else
      return std::complex<double>(n - m + 1, n + m + 1);
}

void write2Dcsv(const std::string& filename, const std::vector<std::vector<double>>& data) {
   std::ofstream ofs(filename);
   for (const auto& row : data) {
      for (size_t j = 0; j < row.size(); ++j) {
         ofs << std::setprecision(15) << row[j];
         if (j < row.size() - 1)
            ofs << ",";
      }
      ofs << std::endl;
   }
   std::cout << "write2Dcsv: " << filename << ", size = {" << data.size() << ", " << data[0].size() << "}" << std::endl;
   ofs.close();
};

void write2Dcsv(const std::string& filename, const std::vector<std::vector<std::complex<double>>>& data) {
   std::ofstream ofs(filename);
   for (const auto& row : data) {
      for (size_t j = 0; j < row.size(); ++j) {
         ofs << std::setprecision(15) << row[j].real();
         if (j < row.size() - 1)
            ofs << ",";
      }
      ofs << std::endl;
   }
   std::cout << "write2Dcsv: " << filename << ", size = {" << data.size() << ", " << data[0].size() << "}" << std::endl;
   ofs.close();
};

std::vector<std::vector<double>> Re(const auto& data) {
   std::vector<std::vector<double>> result(data.size(), std::vector<double>(data[0].size()));
   for (size_t i = 0; i < data.size(); ++i)
      for (size_t j = 0; j < data[0].size(); ++j)
         result[i][j] = data[i][j].real();
   return result;
}
std::vector<std::vector<double>> Im(const auto& data) {
   std::vector<std::vector<double>> result(data.size(), std::vector<double>(data[0].size()));
   for (size_t i = 0; i < data.size(); ++i)
      for (size_t j = 0; j < data[0].size(); ++j)
         result[i][j] = data[i][j].imag();
   return result;
}

int main() {

   //! シンプルな畳み込みの計算 convolutionのサイズは，
   auto conv = [](const int j, const int k) {
      std::complex<double> conv_jk(0, 0);
      for (int n = 0; n <= N; ++n)
         for (int m = -2 * n; m <= 2 * n; ++m)
            conv_jk += f(n, m) * g(j + n, k - m);
      return conv_jk;
   };

   //! f, gで行列を作成し，シンプルな畳み込み相互相関を計算，DFTを使った計算とで比較する．

   std::array<std::array<std::complex<double>, 2 * N + 1>, N + 1> F, G;

   for (int p = 0; p <= N; ++p) {
      F[p].fill(std::complex<double>(0, 0));
      G[p].fill(std::complex<double>(0, 0));
      for (int q = 0; q <= 2 * N; ++q) {
         F[p][q] = f(p, q - N);
         G[p][q] = g(p, q - N);
      }
   }

   write2Dcsv("./Re_F.csv", Re(F));
   write2Dcsv("./Im_F.csv", Im(F));
   write2Dcsv("./Re_G.csv", Re(G));
   write2Dcsv("./Im_G.csv", Im(G));

   /* ------------------------------------------------ */
   Fourier2D<std::complex<double>> fourier2D(G.size(), G[0].size(), F.size(), F[0].size());
   fourier2D.add(G, {true, false}, F);
   auto convFFT = InverseDFT(fourier2D);
   auto reFFT = Re(convFFT);
   auto ImFFT = Im(convFFT);
   write2Dcsv("./Re_DFTconvolution.csv", reFFT);
   write2Dcsv("./Im_DFTconvolution.csv", ImFFT);

   std::cout << "FFT size = {" << convFFT.size() << ", " << convFFT[0].size() << "}" << std::endl;

   std::array<std::array<std::complex<double>, 4 * N + 1>, 2 * N + 1> convolution, convolutionFFT;

   for (int j = 0; j <= N; ++j)
      for (int k = -2 * N; k <= 2 * N; ++k) {
         convolution[j][k + 2 * N] = conv(j, k);
         convolutionFFT[j][k + 2 * N] = convFFT[N - j][k + 2 * N];
      }

   write2Dcsv("./Re_convolution.csv", Re(convolution));
   write2Dcsv("./Im_convolution.csv", Im(convolution));
   write2Dcsv("./Re_convolutionFFT.csv", Re(convolutionFFT));
   write2Dcsv("./Im_convolutionFFT.csv", Im(convolutionFFT));

   for (int j = 0; j <= N; ++j) {
      for (int k = -N; k <= N; ++k)
         std::cout << std::abs(conv(j, k) - convFFT[N - j][k + 2 * N]) << ", ";
      std::cout << std::endl;
   }
   /* -------------------------------------------------------------------------- */
}
