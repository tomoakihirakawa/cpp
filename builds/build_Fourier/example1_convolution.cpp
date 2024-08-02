/*DOC_EXTRACT FourierTransform

## 畳み込み積分

畳み込み積分は，信号解析や画像処理などの分野でよく使われる演算．
直感的には，2つの関数をスライドさせながら積分することで，2つの関数の類似度を評価する．


*/

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

std::complex<double> coeff(const std::vector<double>& sample, const int n) {
   int N = sample.size();
   std::complex<double> sum = 0;
   for (int k = 0; k <= N - 1; ++k) {
      sum += std::polar(sample[k], -n * 2 * M_PI / N * k);
   }
   // sample[k] means f(k * dt) or f(k * T / N). T is the maximum time.
   return sum / static_cast<double>(N);
};

std::vector<std::complex<double>> DFT(const std::vector<double>& sample) {
   int N = sample.size();
   std::vector<std::complex<double>> result(N);
   for (int n = 0; n < N; ++n)
      result[n] = coeff(sample, n);
   return result;
}

int main() {
   const std::vector<double> list = {1, 1, 2, 2, 1, 1, 0, 0};
   for (int n = 0; n < list.size(); ++n)
      std::cout << coeff(list, n) << std::endl;

   for (auto&& c : DFT(list))
      std::cout << c << std::endl;

   std::vector<double> list2(1000);
   auto f = [](double t) {
      double T = 15.;
      return std::sin(2 * M_PI / T * t);
   };

   double T0 = 150.;
   double N = list2.size();
   {
      std::ofstream ofs("original.dat");
      for (int k = 0; k < N; ++k) {
         list2[k] = f(k * T0 / N);
         ofs << list2[k] << std::endl;
      }
      ofs.close();
   }

   {
      std::ofstream ofs("dft.dat");
      int n = 0;
      for (auto&& c : DFT(list2)) {
         ofs << (double)(n / T0) << " " << c.real() << " " << c.imag() << std::endl;
         n++;
         // the index of c is n = 0, 1, ..., N-1
         // c[n] shows the component of frequency n * 2 * pi / T0
      }
      ofs.close();
   }
   return 0;
}
