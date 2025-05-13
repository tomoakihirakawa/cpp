#ifndef lib_Fourier_H
#define lib_Fourier_H

#include <cmath>
#include <complex>
#include <execution>  // C++17以上
#include <numeric>    // for std::transform_reduce
#include <vector>

/* -------------------------------------------------------------------------- */

// FFT and IFFT implementation for complex data
template <typename T>
void FFT(std::vector<std::complex<T>>& data, bool inverse = false) {
   size_t n = data.size();
   if (n <= 1) return;

   // Split data into even and odd parts
   std::vector<std::complex<T>> even(n / 2), odd(n / 2);
   for (size_t i = 0; i < n / 2; ++i) {
      even[i] = data[i * 2];
      odd[i] = data[i * 2 + 1];
   }

   FFT(even, inverse);
   FFT(odd, inverse);

   T angle = 2 * M_PI / n * (inverse ? 1 : -1);
   std::complex<T> w(1), wn(std::cos(angle), std::sin(angle));
   for (size_t i = 0; i < n / 2; ++i) {
      data[i] = even[i] + w * odd[i];
      data[i + n / 2] = even[i] - w * odd[i];
      if (inverse) {
         data[i] /= 2;
         data[i + n / 2] /= 2;
      }
      w *= wn;
   }
}

// FFT (Forward Transform) for real data
template <typename T>
std::vector<std::complex<T>> FFT(const std::vector<T>& sample) {
   if (sample.empty()) return {};

   size_t n = 1;
   while (n < sample.size()) n *= 2;  // Ensure size is a power of 2

   std::vector<std::complex<T>> data(n, 0);
   for (size_t i = 0; i < sample.size(); ++i) {
      data[i] = static_cast<std::complex<T>>(sample[i]);
   }

   FFT(data, false);
   return data;
}

// FFT (Forward Transform) for complex data
template <typename T>
std::vector<std::complex<T>> FFT(const std::vector<std::complex<T>>& sample) {
   if (sample.empty()) return {};

   size_t n = 1;
   while (n < sample.size()) n *= 2;  // Ensure size is a power of 2

   std::vector<std::complex<T>> data(sample.begin(), sample.end());
   data.resize(n, 0);  // Zero-pad to power of 2

   FFT(data, false);
   return data;
}

// IFFT (Inverse Transform) for complex data
template <typename T>
std::vector<std::complex<T>> IFFT(const std::vector<std::complex<T>>& sample) {
   if (sample.empty()) return {};

   size_t n = 1;
   while (n < sample.size()) n *= 2;

   std::vector<std::complex<T>> data(sample.begin(), sample.end());
   data.resize(n, 0);  // Zero-pad to power of 2

   FFT(data, true);
   return data;
}

// Convolution using FFT
template <typename T>
std::vector<std::complex<T>> ConvolveFFT(const std::vector<T>& f, const std::vector<T>& g) {
   if (f.empty() || g.empty()) return {};

   size_t len = 1;
   while (len < f.size() + g.size() - 1) len *= 2;  // Ensure size is a power of 2

   std::vector<std::complex<T>> F(len, 0), G(len, 0);
   for (size_t i = 0; i < f.size(); ++i) F[i] = f[i];
   for (size_t i = 0; i < g.size(); ++i) G[i] = g[i];

   // FFT
   FFT(F);
   FFT(G);

   // Multiply in frequency domain
   for (size_t i = 0; i < len; ++i) {
      F[i] *= G[i];
   }

   // IFFT
   auto result = IFFT(F);

   // Resize result to the expected convolution length
   result.resize(f.size() + g.size() - 1);
   return result;
}

// Overloaded Convolution for complex data
template <typename T>
std::vector<std::complex<T>> ConvolveFFT(const std::vector<std::complex<T>>& f, const std::vector<std::complex<T>>& g) {
   if (f.empty() || g.empty()) return {};

   size_t len = 1;
   while (len < f.size() + g.size() - 1) len *= 2;  // Ensure size is a power of 2

   std::vector<std::complex<T>> F(len, 0), G(len, 0);
   for (size_t i = 0; i < f.size(); ++i) F[i] = f[i];
   for (size_t i = 0; i < g.size(); ++i) G[i] = g[i];

   // FFT
   FFT(F);
   FFT(G);

   // Multiply in frequency domain
   for (size_t i = 0; i < len; ++i) {
      F[i] *= G[i];
   }

   // IFFT
   auto result = IFFT(F);

   // Resize result to the expected convolution length
   result.resize(f.size() + g.size() - 1);
   return result;
}

/* -------------------------------------------------------------------------- */

template <typename T>
std::vector<std::complex<double>> DFT(const std::vector<T>& sample) {
   //! sample can be real or complex
   int s = sample.size();
   std::vector<std::complex<double>> result(s);
   std::complex<double> sum = 0;
   double c;
   for (int n = 0; n < s; ++n) {
      c = -n * 2. * M_PI / s;

      // sum = 0;
      // for (int k = 0; k <= s - 1; ++k)
      //    sum += sample[k] * std::polar(1.0, c * k);  //! sampleが複素数の場合に対応するためpolarの引数としてsampleは渡せない

      std::complex<double> sum = std::transform_reduce(std::execution::par_unseq, sample.begin(), sample.end(), std::complex<double>{0},
                                                       std::plus<>(),
                                                       [&](const auto& val) {
                                                          int k = &val - &sample[0];            //! これはポインタの差分でkを求める方法
                                                          return val * std::polar(1.0, c * k);  //! sampleが複素数の場合に対応するためpolarの引数としてsampleは渡せない
                                                       });
      result[n] = sum / static_cast<double>(s);
   }
   return result;
}

template <typename T>
std::vector<std::vector<std::complex<double>>> DFT(const std::vector<std::vector<T>>& sample2D) {
   //! sample can be real or complex
   int N = sample2D.size();     //! 行
   int M = sample2D[0].size();  //! 列
   std::vector<std::vector<std::complex<double>>> result(N, std::vector<std::complex<double>>(M));
   std::complex<double> sum = 0;
   double cx, cy;
   /*
   for (int n = 0; n < N; ++n)
      for (int m = 0; m < M; ++m) {
         sum = 0;
         cx = -n * 2 * M_PI / N;
         cy = -m * 2 * M_PI / M;
         for (int k = 0; k <= N - 1; ++k)
            for (int j = 0; j <= M - 1; ++j)
               sum += std::polar(sample2D[k][j], cx * k + cy * j);
         result[n][m] = sum / static_cast<double>(N * M);
      }
   return result;
   */
   // 上のようにまとめてDFTを行うと，O(N^2M^2)の計算量になり，計算効率が悪い．
   // このように2Dの場合は，行ごとにDFTを行い，次に列ごとにDFTを行うのが良い．
   std::vector<std::vector<std::complex<double>>> cn2D(N, std::vector<std::complex<double>>(M));
   std::vector<std::vector<std::complex<double>>> cn2Dtmp(N, std::vector<std::complex<double>>(M));
   // ! x方向のフーリエ変換
   for (int i = 0; i < N; ++i)
      cn2D[i] = DFT(sample2D[i]);  //! 長さM
   // ! y方向のフーリエ変換 cn2D[;;][j]
   for (int j = 0; j < M; ++j) {
      std::vector<std::complex<double>> cn_j(N);
      for (int i = 0; i < N; ++i)
         cn_j[i] = cn2D[i][j];  //! 長さN
      cn_j = DFT(cn_j);         //! 長さN
      for (int i = 0; i < N; ++i)
         cn2Dtmp[i][j] = cn_j[i];
   }
   // ! cn2Dtmp[;;][j]をcn2D[;;][j]にコピー
   return cn2Dtmp;
}

template <typename T, size_t N>
std::array<std::complex<double>, N> DFT(const std::array<T, N>& sample) {
   //! sample can be real or complex
   int s = sample.size();
   double S = s;
   std::array<std::complex<double>, N> result;
   std::complex<double> sum = 0;
   double c;
   for (int n = 0; n < s; ++n) {
      sum = 0;
      c = -n * 2 * M_PI / s;
      for (int k = 0; k <= s - 1; ++k)
         sum += sample[k] * std::polar(1.0, c * k);
      result[n] = sum / S;
   }
   return result;
}

/* -------------------------------------------------------------------------- */

template <typename T>
std::vector<std::complex<double>> InverseDFT(const std::vector<T>& sample) {
   //! sample can be real or complex
   int s = sample.size();
   std::vector<std::complex<double>> result(s);
   std::complex<double> sum = 0;
   double c;
   for (int n = 0; n < s; ++n) {
      sum = 0;
      c = n * 2 * M_PI / s;
      for (int k = 0; k <= s - 1; ++k)
         sum += sample[k] * std::polar(1.0, c * k);
      result[n] = sum * (double)s;
   }
   return result;
}

template <typename T>
std::vector<std::vector<std::complex<double>>> InverseDFT(const std::vector<std::vector<T>>& sample2D) {
   //! sample can be real or complex
   int N = sample2D.size();     //! 行
   int M = sample2D[0].size();  //! 列
   std::vector<std::vector<std::complex<double>>> result(N, std::vector<std::complex<double>>(M));
   std::complex<double> sum = 0;
   double cx, cy;
   std::vector<std::vector<std::complex<double>>> cn2D(N, std::vector<std::complex<double>>(M));
   std::vector<std::vector<std::complex<double>>> cn2Dtmp(N, std::vector<std::complex<double>>(M));
   // ! x方向のフーリエ変換
   for (int i = 0; i < N; ++i)
      cn2D[i] = InverseDFT(sample2D[i]);  //! 長さM
   // ! y方向のフーリエ変換 cn2D[;;][j]
   for (int j = 0; j < M; ++j) {
      std::vector<std::complex<double>> cn_j(N);
      for (int i = 0; i < N; ++i)
         cn_j[i] = cn2D[i][j];  //! 長さN
      cn_j = InverseDFT(cn_j);  //! 長さN
      for (int i = 0; i < N; ++i)
         cn2Dtmp[i][j] = cn_j[i];
   }
   // ! cn2Dtmp[;;][j]をcn2D[;;][j]にコピー
   return cn2Dtmp;
}

template <typename T, size_t N>
std::array<std::complex<double>, N> InverseDFT(const std::array<T, N>& sample) {
   //! sample can be real or complex
   int s = sample.size();
   std::array<std::complex<double>, N> result;
   std::complex<double> sum = 0;
   double c;
   for (int n = 0; n < s; ++n) {
      sum = 0;
      c = n * 2 * M_PI / s;
      for (int k = 0; k <= s - 1; ++k)
         sum += sample[k] * std::polar(1.0, c * k);
      result[n] = sum * (double)s;
   }
   return result;
}

/* -------------------------------------------------------------------------- */
template <typename T>
std::vector<std::complex<double>> DiscreteConvolve(const std::vector<T>& f, const std::vector<T>& g) {
   int len = f.size() + g.size() - 1;
   std::vector<std::complex<double>> F(len, 0), G(len, 0);
   for (int i = 0; i < f.size(); ++i) F[i] = f[i];
   for (int i = 0; i < g.size(); ++i) G[i] = g[i];

   std::vector<std::complex<double>> FourierGF(len, 0);
   F = DFT(F);
   G = DFT(G);
   for (int n = 0; n < len; ++n)
      FourierGF[n] += F[n] * G[n];

   return InverseDFT(FourierGF);
}

template <typename T, size_t N>
std::array<std::complex<double>, N> DiscreteConvolve(const std::array<T, N>& f, const std::array<T, N>& g) {
   int len = f.size() + g.size() - 1;
   std::array<std::complex<double>, N> F, G;
   for (int i = 0; i < f.size(); ++i) F[i] = f[i];
   for (int i = 0; i < g.size(); ++i) G[i] = g[i];

   std::array<std::complex<double>, N> FourierGF;
   F = DFT(F);
   G = DFT(G);
   for (int n = 0; n < len; ++n)
      FourierGF[n] += F[n] * G[n];

   return InverseDFT(FourierGF);
}

template <typename T>
struct DiscreteConvolveClass {

   int len;
   std::vector<std::complex<double>> FourierF, FourierG, FourierGF;
   std::vector<std::complex<double>> InverseFourierGF;

   DiscreteConvolveClass(const std::vector<T>& f, const std::vector<T>& g) : len(f.size() + g.size() - 1),
                                                                             FourierF(len, 0),
                                                                             FourierG(len, 0),
                                                                             FourierGF(len, 0) {
      for (int i = 0; i < f.size(); ++i)
         FourierF[i] = f[i];
      for (int i = 0; i < g.size(); ++i)
         FourierG[i] = g[i];

      FourierF = DFT(FourierF);
      FourierG = DFT(FourierG);

      for (int n = 0; n < len; ++n)
         FourierGF[n] += FourierF[n] * FourierG[n];

      InverseFourierGF = InverseDFT(FourierGF);
   };
};

#endif