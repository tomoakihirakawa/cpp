/*DOC_EXTRACT FourierTransform

複素フーリエ級数展開は

```math
f(t) = \sum_{n=-\infty}^{\infty} c_n \exp(i n \omega_0 t), \quad c_n = \frac{1}{T} \int_{-\frac{T}{2}}^{\frac{T}{2}} f(t) \exp(-i n \omega_0 t) \, dt, \quad \omega_0 = \frac{2\pi}{T}
```

サンプル数が$N+1$，$(k=0,1,...N+1)$のとき，台形則を使った関数$g(t)$の数値積分は，

```math
\int_0^T g(t) dt = \left[\frac{g(0) + g(N \delta t)}{2} + \sum_{k=1}^{N-1} g(k \delta t) \right] \delta t, \quad \delta t = \frac{T}{N+1}
```

ここではサンプル数が$N$，$(k=0,1,...N-1)$のとき，
$g(t) = f(t) \exp(-i n \omega_0 t)$として，$c_n$を台形則を使って計算する．

```math
\begin{align}
c_n &= \frac{1}{T} \left[ \frac{g(0) + g((N-1)\delta t)}{2} + \sum_{k=1}^{N-2} g(k \delta t) \right] \delta t, \quad \delta t = \frac{T}{N}, \quad g(0) = g((N-1)\delta t)\\
&= \frac{1}{N} \sum_{k=0}^{N-2} g(k \delta t)\\
&= \frac{1}{N} \sum_{k=0}^{N-2} \left[ f(k \delta t) \exp\left( -i n \frac{2 \pi}{T} k \frac{T}{N} \right) \right]\\
&= \frac{1}{N} \sum_{k=0}^{N-2} \left[ f_k \exp\left( -i n \frac{2 \pi}{N} k \right) \right]
\end{align}
```

ここで$\delta t$は区間$T$を$N$等分したときの各小区間の長さを表す．
そして$g(0) = g((N-1)\delta t)$は，関数$g(t)$が周期$T$であることを意味する．
この式は，離散フーリエ変換（DFT）における$c_n$の近似計算に相当する．

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
#include <iostream>
#include <vector>

const std::vector<double> list = {1, 1, 2, 2, 1, 1, 0, 0};

std::complex<double> coeff(int n) {
   int N = list.size();
   std::complex<double> sum = 0;
   for (int k = 0; k < N; ++k) {
      sum += list[k] * std::polar(1.0, -n * 2 * M_PI / N * k);
   }
   return sum / static_cast<double>(N);
};

std::vector<std::complex<double>> DFT(const std::vector<double>& list) {
   int N = list.size();
   std::vector<std::complex<double>> result(N);
   for (int n = 0; n < N; ++n) 
      result[n] = coeff(n);   
   return result;
}

int main() {
   for (int n = 0; n < list.size(); ++n)
      std::cout << coeff(n) << std::endl;

   for (auto&& c : DFT(list))
      std::cout << c << std::endl;

   return 0;
}
