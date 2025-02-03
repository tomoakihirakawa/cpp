/*DOC_EXTRACT FourierTransform

## 畳み込み積分

畳み込み積分は，関数`g`をスライドさせながら`f`と掛け合わせ和を求めることである．

```math
(f \ast g)(t) = \int_{-\infty}^{\infty} dx f(x) g(t-x)
```

の形の積分である．ここで，$`\ast`$は畳み込み積分を表す．離散データの畳み込み積分は，次のように計算できる．

```math
(f \ast g)_j = \sum_{k=0}^{N-1} f_k g_{j-k}
```

畳み込み積分の値は，フーリエ変換された`f`と`g`の積を逆フーリエ変換した結果と等しい．

```math
(f \ast g)_j = \mathcal{F}^{-1}[\mathcal{F}[f] \cdot \mathcal{F}[g]]_j
```

$`\mathcal{F}^{-1}[\mathcal{F}[f] \cdot \mathcal{F}[g]]`$はデータ列であって，`_j`は，そのデータ列の`j`番目をとることを意味する．

全て逆フーリエ変換してから一つだけを抜き出す必要はなく，`j`だけが必要なら`j`番目のデータだけを取り出すよう和を取ればいい．

```math
(f \ast g)_j = \sum_{n=0}^{N-1} (\mathcal{F}[f] \cdot \mathcal{F}[g])_n \exp\left(i n \frac{2\pi}{N} j\right)
```

以下は，離散フーリエ変換，逆フーリエ変換，畳み込み積分を行うMatheamticaのコードである．

データをスライドさせて掛け合わせた結果できるデータ列の和は，スライド`len=Length[f] + Length[g] - 1`まで値を持ちえる．
これ以上のスライド＆掛け算の結果はゼロになる．なので，畳み込み和が返すデータ列は，`Length[f] + Length[g] - 1`となる．

離散フーリエを使った，畳み込み和の内部では，まず`f`と`g`の長さをゼロ埋めして`len`の長さに揃えた後，それぞれのフーリエ変換を計算し，掛け合わせる．
その結果を逆フーリエ変換して，畳み込み和を求める．

```Mathematica
MyFourier[list_, n_] := With[{len = Length[list], c = -I*n*2*\[Pi]/Length[list]},
Sum[list[[k + 1]]*Exp[c*k], {k, 0, len - 1}]/len
];

MyInverseFourier[list_, n_] := With[{len = Length[list], c = I*n*2*\[Pi]/Length[list]},
Sum[list[[k + 1]]*Exp[c*k], {k, 0, len - 1}]
];

MyDiscreteConvolve[f_, g_] := Module[{len, F, G, FourierGF},
len = Length[f] + Length[g] - 1;
F = PadRight[f, len];
G = PadRight[g, len];
FourierGF = Table[MyFourier[F, n]*MyFourier[G, n], {n, 0, len - 1}];
Return[N@Table[len*MyInverseFourier[FourierGF, n], {n, 0, len - 1}]];
]
```

![sample_conv.png](sample_conv.png)

(see `example0.nb`)

*/

/*
```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_convolution.cpp
make
./example1_convolution
```
*/

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include "lib_Fourier.hpp"
#include "lib_measurement.hpp"

std::ostream& operator<<(std::ostream& os, const std::vector<std::complex<double>>& cn) {
   for (auto&& c : cn)
      os << c << std::endl;
   return os;
}

int main() {
   // const std::vector<double> f = {1, 2, 3, 4, 5, 4, 3, 2, 1};
   // const std::vector<double> g = {1, -1, 1, -1, 1, -1, 3, 0, 0};
   // auto len = f.size() + g.size() - 1;
   // std::vector<double> convolution(len, 0);
   // double sum = 0;

   // const std::vector<std::complex<double>> f = {1, 2, 3, 4, 5, 4, 3, 2, 1, 1};
   // const std::vector<std::complex<double>> g = {1, -1, 1, -1, 1, -1, 3, 0, 0, 1};

   std::vector<std::complex<double>> f = {1, 2, 3, 4, 5, 4, 3, 2, 1, 1};
   std::vector<std::complex<double>> g = {1, -1, 1, -1, 1, -1, 3, 0, 0, 1};

   f.resize(16 * 2, 0);
   g.resize(16 * 2, 0);

   auto len = f.size() + g.size() - 1;
   std::vector<std::complex<double>> convolution(len, 0);
   std::complex<double> sum = 0;

   TimeWatch tw;

   for (int x = 0; x < len; ++x) {  // Correct range
      sum = 0;
      for (int i = 0; i < len; ++i) {
         auto a = (0 <= i && i < f.size()) ? f[i] : 0;
         auto b = (0 <= x - i && x - i < g.size()) ? g[x - i] : 0;
         sum += a * b;
      }
      convolution[x] = sum;
   }

   std::cout << "convolution" << std::endl;
   // for (auto i = 0; i < len; ++i)
   //    std::cout << i << " " << convolution[i] << std::endl;

   std::cout << "Elapsed time : " << tw()[0] << std::endl;

   std::cout << "DiscreteConvolve" << std::endl;
   DiscreteConvolveClass DCC(f, g);
   // std::cout << DCC.FourierF << std::endl;
   // std::cout << DCC.FourierG << std::endl;
   // std::cout << DCC.FourierGF << std::endl;
   // std::cout << DCC.InverseFourierGF << std::endl;

   std::cout << "Elapsed time : " << tw()[0] << std::endl;

   ConvolveFFT(f, g);
   // std::cout << ConvolveFFT(f, g) << std::endl;

   std::cout
       << "Elapsed time : " << tw()[0] << std::endl;
   return 0;
   // {1., 1., 2., 2., 3., 1., 4., 5., 8., 9., 13., 10., 8., 5., 3.}
}
