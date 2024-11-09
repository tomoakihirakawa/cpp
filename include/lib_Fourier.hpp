#ifndef lib_Fourier_H
#define lib_Fourier_H

#include <cmath>
#include <complex>
#include <vector>

template <typename T>
std::vector<std::complex<double>> DFT(const std::vector<T>& sample) {
   //! sample can be real or complex
   int s = sample.size();
   std::vector<std::complex<double>> result(s);
   std::complex<double> sum = 0;
   double c;
   for (int n = 0; n < s; ++n) {
      sum = 0;
      c = -n * 2 * M_PI / s;
      for (int k = 0; k <= s - 1; ++k)
         sum += sample[k] * std::polar(1.0, c * k);
      result[n] = sum / static_cast<double>(s);
   }
   return result;
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
      result[n] = sum;
   }
   return result;
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
      result[n] = sum;
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
      FourierGF[n] += F[n] * G[n] * (double)len;

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
      FourierGF[n] += F[n] * G[n] * (double)len;

   return InverseDFT(FourierGF);
}

#endif