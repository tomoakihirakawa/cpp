/*DOC_EXTRACT FourierTransform

# フーリエ変換

## 複素フーリエ級数展開

```math
f(t) = \sum_{n=-\infty}^{\infty} c_n \exp(i n \omega^\ast t), \quad c_n = \frac{1}{T^\ast} \int_{-\frac{T^\ast}{2}}^{\frac{T^\ast}{2}} f(t) \exp(-i n \omega^\ast t) \, dt, \quad \omega^\ast = \frac{2\pi}{T^\ast}
```

$`\exp({i \theta}) = \cos \theta + i \sin \theta`$なので，
フーリエ係数の実部には，$`\cos \theta`$の係数が，虚部には，$`\sin \theta`$の係数が含まれる．

$`c_n=\frac{a_n - i \mathrm{sgn}(n) b_n}{2}`$

## 離散フーリエ変換（インデックス周期$`N`$のフーリエ変換）

次のようなで$`N`$個の離散データがあるとする．

```cpp
{1, 1, 2, 2, 1, 1, 0, 0}
```

これが，周期的に繰り返すとする．

```cpp
{1, 1, 2, 2, 1, 1, 0, 0},{1, 1, 2, 2, 1, 1, 0, 0},{1, 1, 2, 2, 1, 1, 0, 0},...
```
初めのデータを$`0`$番として数えると，$`N`$番目のデータは$`0`$番目のデータと等しいことになる．
この無限に続く数字をフーリエ級数で表現するなら，$`0`$番目と$`N`$番目のデータは，級数を構成する三角関数の$`0`$と$`2\pi`$に対応させるのが自然だろう．
つまり，dataとindex，angle，periodの対応は次のようになる．

```cpp
data  : {1, 1, 2, ...,                               0, 0}, {1, 1, ...
index : {0, 1, 2, ...,                           N-2, N-1}, {N, N+1, ...
angle : {0, 2pi/N, 2pi*2/N, ..., 2pi*(N-2)/N, 2pi*(N-1)/N}, {2pi, 2pi*(N+1)/N, ...
period: {0,   T/N,    2T/N, ...,   T*(N-2)/N,   T*(N-1)/N}, {T,     T*(N+1)/N, ...
```

複素フーリエ係数$`c_n`$を台形則で数値積分すると，

```math
\begin{align}
c_n &= \frac{1}{T^\ast} \left[ \frac{g_n(0) + g_n(N\delta t)}{2} + \sum_{k=1}^{N-1} g_n(k \delta t) \right] \delta t, \quad \delta t = \frac{T^\ast}{N}, \quad g_n(0) = g_n(N\delta t),\quad g_n(t) = f(t) \exp(-i n \omega^\ast t)\\
&= \frac{1}{N} \sum_{k=0}^{N-1} g_n(k \delta t) {\quad\text{became simple additions}}\\
&= \frac{1}{N} \sum_{k=0}^{N-1} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{T^\ast} k \frac{T^\ast}{N} \right) \right]\\
&= \frac{1}{N} \sum_{k=0}^{N-1} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{N} k \right) \right]
\end{align}
```

これからわかるように，$`c_n`$は周期$`T^\ast`$に依存しておらず（$`f(kT^\ast/N)`$は，$`T^\ast`$によらず常に$`k`$番めデータ値を指しているので，$`T^\ast`$に依存していない），データの数$`N`$に依存している．

この結果は，複素フーリエ係数$`c_n`$の式において，$`T^\ast`$を$`N`$として数値積分したものとも考えられる．つまり，時間軸ではなく，インデックス軸で積分していることと同じになっている．

NOTE: $`c_n`$を変形すると，

$$
c_n = \frac{1}{N} \sum_{k=0}^{N-1} \left[ f\left(k\frac{T^\ast}{N}\right)
\exp\left(k\right)\right]\exp\left( -i n \frac{2 \pi}{N}\right)
$$

となっており，$`n`$に関して周期$`N`$の周期関数となっている．
また$`\cos(\theta)=\cos(-\theta)`$であるため，$`\Re[c_n]=\Re[c_{-n}]`$で
$`\sin(\theta)=-\sin(-\theta)`$であるため，$`\Im[c_n]=-\Im[c_{-n}]`$である．
１周期分つまり$`N`$分の係数ではなく，$`N/2`$分の係数さえわかれば元の関数を復元できる．

---

Mathematicaの組み込み関数と比較して確かめてみる．
Mathematicaの`Fourier`関数の`FourierParameters`オプションが，`{-1,-1}`の場合に上記の式と一致する．
`MyFourier`は全く同じ結果を返す．

```Mathematica
(*example0.nb*)
list = {1., 1., 2., 2., 1., 1., 0., 0.};
MyFourier[list_, n_] := With[{len = Length[list]}, Sum[list[[k + 1]]*Exp[-I*n*2. \[Pi]/len*k], {k, 0, len - 1}]/len];
Column[Fourier[list, FourierParameters -> {-1, -1}], Frame -> All]
Column[cn = Table[MyFourier[list, n], {n, 0, Length[list] - 1}], Frame -> All]
```

c++での`MyFourier`と同じ関数を作って実行してみる．実行方法：

```cpp
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_simple.cpp
./example0_simple
```

結果は，Matheamticaと同じになる．

| list | 元データ| 離散フーリエ変換 |
|:---:|:---:|:---:|
| random list | ![sample_original_mathematica.png](sample_original_mathematica.png) | ![sample_ReIm_cn_mathematica.png](sample_ReIm_cn_mathematica.png) |
| cos wave | ![sample_original_cosWave.png](sample_original_cosWave.png) | ![sample_ReIm_cn_cosWave.png](sample_ReIm_cn_cosWave.png) |
| square wave | ![sample_original_squareWave.png](sample_original_squareWave.png) | ![sample_ReIm_cn_squareWave.png](sample_ReIm_cn_squareWave.png) |
| triangle wave | ![sample_original_triangleWave.png](sample_original_triangleWave.png) | ![sample_ReIm_cn_triangleWave.png](sample_ReIm_cn_triangleWave.png) |


## 逆離散フーリエ変換

フーリエ係数$`c_n`$から元の関数$`f(t)`$を復元することを考える．
三角関数を掛けて積分することで係数を抽出できたので，その方法で関数を抽出する．

```math
f\left(k\frac{T^\ast}{N}=k\delta t\right) = \frac{1}{N} \sum_{k=0}^{N-1} \left[ c_n \exp\left( i n \frac{2 \pi}{N} k \right) \right]
```

```Mathematica
list = {1., 1., 2., 2., 1., 1., 0., 0.};
MyInverseFourier[list_, n_] := With[{len = Length[list]}, Sum[list[[k + 1]]*Exp[I*n*2 \[Pi]/len*k], {k, 0, len - 1}]/len];
Column[InverseFourier[cn, FourierParameters -> {-1, -1}], Frame -> All]
Column[Table[MyInverseFourier[cn, n]*(Length[list]), {n, 0, Length[list] - 1,1}], Frame -> All]
```

1周期分積分するので，$`T^\ast`$で割っている．$`\delta t`$とかけるので，結果として$`N`$で割ることになる．

## 離散フーリエ変換によるデータの補間

| list | 逆フーリエ変換 | フーリエ級数展開 |
|:---:|:---:|:---:|
| random list | ![sample_Re_inv_mathematica.png](sample_Re_inv_mathematica.png) | ![sample_interpolation_mathematica.png](sample_interpolation_mathematica.png) |
| cos wave | ![sample_Re_inv_cosWave.png](sample_Re_inv_cosWave.png) | ![sample_interpolation_cosWave.png](sample_interpolation_cosWave.png) |
| square wave | ![sample_Re_inv_squareWave.png](sample_Re_inv_squareWave.png) | ![sample_interpolation_squareWave.png](sample_interpolation_squareWave.png) |
| triangle wave | ![sample_Re_inv_triangleWave.png](sample_Re_inv_triangleWave.png) | ![sample_interpolation_triangleWave.png](sample_interpolation_triangleWave.png) |

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
