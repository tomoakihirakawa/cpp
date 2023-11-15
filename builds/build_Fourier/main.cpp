/*DOC_EXTRACT FourierTransform

# フーリエ変換

## 複素フーリエ級数展開

```math
f(t) = \sum_{n=-\infty}^{\infty} c_n \exp(i n \omega^\ast t), \quad c_n = \frac{1}{T^\ast} \int_{-\frac{T^\ast}{2}}^{\frac{T^\ast}{2}} f(t) \exp(-i n \omega^\ast t) \, dt, \quad \omega^\ast = \frac{2\pi}{T^\ast}
```

省略

## 離散フーリエ変換

サンプル数が$`N+1`$，$`(k=0,1,...N+1)`$のとき，台形則を使った関数$`g(t)`$の数値積分は，

```math
\int_0^{T^\ast} g(t) dt = \left[\frac{g(0) + g(N \delta t)}{2} + \sum_{k=1}^{N-1} g(k \delta t) \right] \delta t, \quad \delta t = \frac{T^\ast}{N+1}
```

この台形則を使って，
サンプル数が$`N`$，$`(k=0,1,...N-1)`$，$`g(t) = f(t) \exp(-i n \omega^\ast t)`$として．
複素フーリエ係数$`c_n`$を数値積分で計算すると，

```math
\begin{align}
c_n &= \frac{1}{T^\ast} \left[ \frac{g(0) + g((N-1)\delta t)}{2} + \sum_{k=1}^{N-2} g(k \delta t) \right] \delta t, \quad \delta t = \frac{T^\ast}{N}, \quad g(0) = g((N-1)\delta t)\\
&= \frac{1}{N} \sum_{k=0}^{N-2} g(k \delta t)\\
&= \frac{1}{N} \sum_{k=0}^{N-2} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{T^\ast} k \frac{T^\ast}{N} \right) \right]\\
&= \frac{1}{N} \sum_{k=0}^{N-2} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{N} k \right) \right]
\end{align}
```

となる．ここで$`\delta t`$は区間$`T^\ast`$を$`N`$等分したときの各小区間の長さを表す．
そして$`g(0) = g((N-1)\delta t)`$は，関数$`g(t)`$が周期$`T^\ast`$であることを意味する．
この式は，離散フーリエ変換（DFT）における$`c_n`$の近似計算に相当する．

$`c_n`$は，$`\omega _n = \frac{2 \pi n}{T^\ast}`$の角周波数成分を表す．つまり，$`f_n = \frac{n}{T^\ast}`$の周波数成分を表し．$`T_n = \frac{T^\ast}{n}`$の周期成分を表す．
このことから，周波数分解能は$`\Delta f = \frac{1}{T^\ast}`$，周期分解能は$`\Delta T = T^\ast`$とわかる．
また，基本周期$`T^\ast`$と同じ周期でサンプリングしていては，$`T^\ast`$の周期成分は分解できない．
少なくとも，$`T^\ast`$の周期成分を分解するには，$`T^\ast`$の半分の周期でサンプリングする必要がある．

---

Mathematicaの組み込み関数と比較して確かめてみる．
Mathematicaの`Fourier`関数の`FourierParameters`オプションが，`{-1,-1}`の場合に上記の式と一致する．

```Mathematica
list = N@{1, 1, 2, 2, 1, 1, 0, 0}

MyFourier[list_, n_] := With[{len = Length[list]},
   Sum[list[[k + 1]]*Exp[-I*n*2 \[Pi]/len*k], {k, 0, len - 2}]/len
   ];

Column@Fourier[list, FourierParameters -> {-1, -1}]
Column@Table[MyFourier[list, n], {n, 0, Length[list] - 1}]
```

*/

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

std::complex<double> coeff(const std::vector<double>& sample, const int n) {
   int N = sample.size();
   std::complex<double> sum = 0;
   for (int k = 0; k <= N - 2; ++k) {
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
   //    const std::vector<double> list = {1, 1, 2, 2, 1, 1, 0, 0};
   //    for (int n = 0; n < list.size(); ++n)
   //       std::cout << coeff(list, n) << std::endl;

   //    for (auto&& c : DFT(list))
   //       std::cout << c << std::endl;

   // check DFT of sin cos wave and output

   std::vector<double> list2(1000);
   auto f = [](double t) {
      double T = 15.;
      return std::sin(2 * M_PI / T * t);
   };

   double T0 = 15.;
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
      std::ofstream("dft.dat");
      std::ofstream ofs("dft.dat");
      for (auto&& c : DFT(list2)) {
         ofs << c.real() << " " << c.imag() << std::endl;
         // the index of c is n = 0, 1, ..., N-1
         // c[n] shows the component of frequency n * 2 * pi / T0
      }
      ofs.close();
   }
   return 0;
}
