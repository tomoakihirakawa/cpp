## ⛵ FMMのM2L変換として作成する係数$L _k^j$の計算 

対象とする積分：

$$
\begin{aligned}
L _{j}^{k}&=\mathop{\sum}\limits _{{n}={0}}\limits^{N}{\mathop{\sum}\limits _{{m}=-{n}}\limits^{n}{{f} _{n}^{m}{g} _{{n}+{j}}^{{m}-{k}}}}\\
&=\mathop{\sum}\limits _{{n}={0}}\limits^{N}{\mathop{\sum}\limits _{{m}=-{n}+{N}}\limits^{{n}+{N}}{{f} _{n}^{{m}-{N}}{g} _{{n}+{j}}^{-\left({{k}-{m}}\right)-{N}}}}\\
&=\mathop{\sum}\limits _{{n}={0}}\limits^{N}{\mathop{\sum}\limits _{{m}={0}}\limits^{2N}{{f} _{n}^{{m}-{N}}{g} _{{n}+{j}}^{-\left({{k}-{m}}\right)-{N}}}}\\
&={\left\{{{\mathcal{F}} _{2D}^{-{1}}\left[{{\mathcal{F}} _{2D}\left[{{\mathbf{f}} _{p}^{{q}-{N}}}\right]\odot{\mathcal{F}} _{2D}\left[{{\mathbf{g}} _{p}^{-{q}-{N}}}\right]}\right]}\right\}} _{j}^{k}
\end{aligned}
$$

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Debug ../ -DSOURCE_FILE=example1_2D_DFT_complex_convolution_check.cpp
make
./example1_2D_DFT_complex_convolution_check
```

[:material-microsoft-visual-studio-code:example1_2D_DFT_complex_convolution_check.cpp#L1](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Fourier/example1_2D_DFT_complex_convolution_check.cpp:1)

---
# 🐋 DFTを使った畳み込み・相互相関の計算 

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_2D_DFT_external_image.cpp
make
./example1_2D_DFT_external_image
```

２つの画像の畳み込み・相互相関をDFTを使って計算してみる．

<figure>
<img src="./img_kernel_shift.png" alt="img_kernel_shift.png" width="500px">
<figcaption>kernel image, shift image注意：縦軸が反転して表示させているので</figcaption>
</figure>

畳み込み&畳み込みと相互相関&畳み込みを行って結果を比較する．
２つの違いは$j-n$か$j+n$かだけ．
$j=0$として$\sum$の中でのインデックスの取り方をイメージすると，
$f _0g _0+f _1g _{-1}+f _2g _{-2}+\dots$というように，インデックスは互いに反対方向に進んでいくのが，畳み込み．
なので，$j=2N _{max}-1$でも

* $j$方向の畳み込みは，$j\geq 0$で値をもつ
* 一方で相互相関は，$j< 0$でも値をもつ

このことは，次のように図でイメージすると明確で，$j$によらない固定された核となる行列に対して，シフトさせるもう片方の行列
をまず$j$だけ進める．次に，相互相関の場合は，２つの行列が重なり合う要素だけで内積をとる．
畳み込みの場合は，$j$を軸にして，シフト行列を反転させてから内積をとる．
明らかに，$j<0$でも相互相関は値を持ち，畳み込みの場合は$j<0$では値を持たないものの，
行列サイズよりも大きな$j$でも値を持ち得る．

<figure>
<img src="./img_conv_corr.png" alt="img_conv_corr.png" width="500px">
<figcaption>畳み込み&畳み込み（左）と相互相関&畳み込み（右）の結果</figcaption>
</figure>

<figure>
<img src="./img_process_conv_corr.png" alt="img_process_conv_corr.png" width="500px">
<figcaption>畳み込みと相互相関の処理の捉え方</figcaption>
</figure>

<figure>
<img src="./img_process_flipped_conv.png" alt="img_process_flipped_conv.png" width="500px">
<figcaption>図を$j$方向に反転させて畳み込むと結果は相互相関と一致する</figcaption>
</figure>

[:material-microsoft-visual-studio-code:example1_2D_DFT_external_image.cpp#L1](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Fourier/example1_2D_DFT_external_image.cpp:1)

---
# 🐋 複素フーリエ級数展開（周期関数） 

## ⛵ 周期関数なら複素フーリエ級数展開 

周期$T$で周波数$W$の周期関数$f(t)$をフーリエ級数展開すると次のようになる．

$$
f(t) = \sum _{n=-\infty}^{\infty} c _n \exp(i n W t),\quad
c _n = \frac{1}{T} \int _{-\frac{T}{2}}^{\frac{T}{2}} f(t) \exp(-i n W t) \, dt, \quad W = \frac{2\pi}{T}
$$

$c _n=\frac{a _n - i \mathrm{sgn}(n) b _n}{2}$

## ⛵ 非周期関数:周期無限大ならフーリエ変換 

非周期関数は，$T^\ast \to \infty$の関数として考えることができるので，極限をとる．

例えば，$f(t)=\cos(px)$の場合を考えると，以下のようになり，

$$
\begin{align*}
\int _{-\infty}^{\infty} \cos(x) e^{-i\omega x} \, dx
&= \int _{-\infty}^{\infty} \frac{e^{-ipx} + e^{ipx}}{2} \cdot e^{-i\omega x} \, dx \\
&= \frac{1}{2} \int _{-\infty}^{\infty} \left( e^{-i(p + \omega)x} + e^{i(p - \omega)x} \right) dx \\
&= \frac{1}{2} \left( 2\pi \delta(p + \omega) + 2\pi \delta(p - \omega) \right)
\end{align*}
$$

# 🐋 離散複素フーリエ級数展開とその係数（離散フーリエ変換） 

離散フーリエ変換と逆離散フーリエ変換を端的に示すと次のようになる．
結果は，Mathematicaの`Fourier`関数，`InverseFourier`関数の`FourierParameters`オプションが，`{-1,-1}`の場合と一致する．

```Mathematica
MyFourier[list_, n_] := With[{len = Length[list]}, Sum[list[[k + 1]]*Exp[-I*n*2 \[Pi]/len*k], {k, 0, len - 1}]/len];
MyInverseFourier[list_, n_] := With[{len = Length[list]}, Sum[list[[k + 1]]*Exp[I*n*2 \[Pi]/len*k], {k, 0, len - 1}]];

(*離散データ*)
u = N@{1, 2, 3, 4, 5, 4, 3, 2, 1};

Grid[
Transpose@{cn = Table[MyFourier[u, n], {n, 0, Length[u] - 1}],
Fourier[u, FourierParameters -> {-1, -1}]}
, Frame -> All]

Grid[
Transpose@{Table[MyInverseFourier[cn, n], {n, 0, Length[cn] - 1}],
InverseFourier[cn, FourierParameters -> {-1, -1}]}
, Frame -> All]
```

## ⛵ 離散フーリエ変換（インデックス周期$N$のフーリエ変換） 

$N$個の離散データ`{1, 1, 2, 2, 1, 1, 0, 0}`があるとする．
これが，周期的に繰り返すとする．

```cpp
{1, 1, 2, 2, 1, 1, 0, 0},{1, 1, 2, 2, 1, 1, 0, 0},{1, 1, 2, 2, 1, 1, 0, 0},...
```

与えられたデータは離散的だが，
これを連続で周期的な関数$f(t)$として表したい．
また，データインデックス上では，対応する与えられたデータ値と一致させたい．
つまり，$i=2$を関数に与えたら，2を返してほしい．これはいわば補間のようなものである．
複素フーリエ級数展開によって，この願いを叶えるためには，
$c _0,...,c _{N-1}$の$N$コの係数があればよい．

$$
f(t) = \sum _{n=0}^{N-1} c _n \exp(i n \omega^\ast t), \quad c _n =\frac{1}{T^\ast} \int _{0}^{T^\ast} f(t) \exp(-i n \omega^\ast t) \, dt, \quad \omega^\ast = \frac{2\pi}{T^\ast}
$$

ここでの添字の振り方は，初めのデータを$0$番として数えることにする．
そう数えると周期性を仮定したので，$N$番目のデータは$0$番目のデータと等しいことになる．
この無限に続く数字をフーリエ級数で表現すると，
$0$番目と$N$番目のデータは，級数を構成する三角関数の$0$と$2\pi$に対応している．
dataとindex，angle，periodの対応は次のようなイメージである．

```cpp
data  : {1, 1, 2, ...,                               0, 0}, {1, 1, ...
index : {0, 1, 2, ...,                           N-2, N-1}, {N, N+1, ...
angle : {0, 2pi/N, 2pi*2/N, ..., 2pi*(N-2)/N, 2pi*(N-1)/N}, {2pi, 2pi*(N+1)/N, ...
period: {0,   T/N,    2T/N, ...,   T*(N-2)/N,   T*(N-1)/N}, {T,     T*(N+1)/N, ...
```

複素フーリエ係数$c _n$を台形則で数値積分すると，

$$
\begin{align}
c _n &= \frac{1}{T^\ast} \left[ \frac{g _n(0) + g _n(N\delta t)}{2} + \sum _{k=1}^{N-1} g _n(k \delta t) \right] \delta t, \quad \delta t = \frac{T^\ast}{N}, \quad g _n(0) = g _n(N\delta t),\quad g _n(t) = f(t) \exp(-i n \omega^\ast t)\\
&= \frac{1}{N} \sum _{k=0}^{N-1} g _n(k \delta t) {\quad\text{became simple additions}}\\
&= \frac{1}{N} \sum _{k=0}^{N-1} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{T^\ast} k \frac{T^\ast}{N} \right) \right]\\
&= \frac{1}{N} \sum _{k=0}^{N-1} \left[ f _k \exp\left( -i n \frac{2 \pi}{N} k \right) \right], \quad f _k = f\left(k\frac{T^\ast}{N}\right)
\end{align}
$$

これからわかるように，$c _n$は周期$T^\ast$に依存しておらず，データの数$N$に依存している．
（$f(kT^\ast/N)$は，$T^\ast$によらず常に$k$番目データ値`data[k]`を指しているので，$T^\ast$に依存していない）
離散フーリエ係数は，周期とは無関係なのである．

$c _n$が大きさを表す波の周波数は，数式から$n/T^\ast$であるとわかる．

最後の式は，連続した関数のフーリエ係数を抽出するための式と照らし合わせると，
$T^\ast$を$N$と置き換えた形になっている．
時間軸ではなく，インデックス軸で積分しているようなものである．

$c _n$は，$c _n=c _{n+N}$であり$n$に関して周期$N$の周期関数となっている．
ただし，実部は偶関数で，虚部は奇関数である．
なので，$c _n=c _{-n}$ではなく，$c _n=\overline{c _{-n}}$である．
ただし，これは複素フーリエ変換の定義そのものの性質である．

$$
c _n = \frac{a _0}{2}, \quad c _n = \frac{a _n - i b _n}{2}, \quad c _{-n} = \frac{a _n + i b _n}{2}, \quad c _{-n} = \overline{c _n}
$$

なので，$c _n$を求めるさい，$n=N/2$より大きいインデックスの係数の計算は省略できる．

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


## ⛵ 逆離散フーリエ変換 

フーリエ係数$c _n$から元の関数$f _\kappa=f(t=\kappa \delta t)$を復元することを考える．
三角関数を掛けて積分することで係数を抽出できたので，その方法で関数を抽出する．

フーリエ変換は次のように定義している．

$$
c _n = \frac{1}{N} \sum _{k=0}^{N-1} \left[ f _k \exp\left( -i n \frac{2 \pi}{N} k \right) \right]
$$

$\exp\left( -i n \frac{2 \pi}{N} k \right)$ではなく，$N\exp\left( i n \frac{2 \pi}{N} \kappa \right)$を掛けて積分する．ここで，$k$と区別するために$\kappa$を使っている．

$$
\begin{equation}
\begin{aligned}
\sum _{n=0}^{N-1}{c _n} \exp\left( i n \frac{2 \pi}{N} \kappa \right)&=N\sum _{n=0}^{N-1}{\frac{1}{N} \sum _{k=0}^{N-1} \left[ f _k \exp\left( -i n \frac{2 \pi}{N} k \right) \right]
} \exp\left( i n \frac{2 \pi}{N} \kappa \right)\\
&=N\sum _{n=0}^{N-1}{\frac{1}{N} \sum _{k=0}^{N-1} \left[ f _k \exp\left( -i n \frac{2 \pi}{N} k \right) \right]
} \exp\left( i n \frac{2 \pi}{N} \kappa \right)\\
&=N\frac{1}{N} f _\kappa \\
&=f\left(\kappa\frac{T^\ast}{N}=\kappa\delta t\right)
\end{aligned}
\end{equation}
$$

このように，フーリエ係数を使って元の関数を復元できる．
フーリエ変換の定義上，$N$を掛ける必要がある．

```Mathematica
list = {1., 1., 2., 2., 1., 1., 0., 0.};
MyInverseFourier[list_, n_] := With[{len = Length[list]}, Sum[list[[k + 1]]*Exp[I*n*2 \[Pi]/len*k], {k, 0, len - 1}]];
Column[InverseFourier[cn, FourierParameters -> {-1, -1}], Frame -> All]
Column[Table[MyInverseFourier[cn, n], {n, 0, Length[list] - 1}],
Frame -> All]
```

## ⛵ 離散フーリエ変換によるデータの補間 

| list | 逆フーリエ変換 | フーリエ級数展開 |
|:---:|:---:|:---:|
| random list | ![sample_Re_inv_mathematica.png](sample_Re_inv_mathematica.png) | ![sample_interpolation_mathematica.png](sample_interpolation_mathematica.png) |
| cos wave | ![sample_Re_inv_cosWave.png](sample_Re_inv_cosWave.png) | ![sample_interpolation_cosWave.png](sample_interpolation_cosWave.png) |
| square wave | ![sample_Re_inv_squareWave.png](sample_Re_inv_squareWave.png) | ![sample_interpolation_squareWave.png](sample_interpolation_squareWave.png) |
| triangle wave | ![sample_Re_inv_triangleWave.png](sample_Re_inv_triangleWave.png) | ![sample_interpolation_triangleWave.png](sample_interpolation_triangleWave.png) |

[:material-microsoft-visual-studio-code:example0_simple.cpp#L1](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Fourier/example0_simple.cpp:1)

## ⛵ 畳み込み積分 

畳み込み積分は，関数`g`をスライドさせながら`f`と掛け合わせ和を求めることである．

$$
(f \ast g)(t) = \int _{-\infty}^{\infty} dx f(x) g(t-x)
$$

ここで，$`\ast`$は畳み込み積分を表す．
畳み込み積分は，$`f`$と$`g`$のフーリエ変換を掛け合わせて逆フーリエ変換したものと等しい．
畳み込み積分は，shiftに関する関数と捉えることができる．

```Mathematica
ClearAll["Global`*"];

(*離散信号の定義*)
f[x_] := Piecewise[{{1/2, 0 <= x}, {-1/2, True}}]
g[x_] := Piecewise[{{1, -1 <= x <= 1}, {0, True}}]

(*面白い例*)
f[x_] := Sech[x]^2;
g[x_] := HeavisideTheta[x] - 1/2;

convolve[x_] = Integrate[f[T]*g[x - T], {T, -Infinity, Infinity}];

(*convolve[x_]=Convolve[g[t],f[t],t,x];*)

F[w_] = FourierTransform[f[x], x, w, FourierParameters -> {1, -1}];
G[w_] = FourierTransform[g[x], x, w, FourierParameters -> {1, -1}];
(*convolve[x_]=InverseFourierTransform[F[w]*G[w],w,x,\
FourierParameters->{1,-1}];*)

Manipulate[
Plot[{f[x], g[shift - x], f[x]*g[shift - x], convolve[x]}, {x, -5,
5}
, AxesLabel -> {"x", ""}
, PlotLegends -> {"f[x]", "g[shift-x]", "f[x]g[shift-x]",
"convolve[x]"}
, PlotLabel -> "shift=" <> ToString[shift]
, PlotStyle -> {{Orange, Dashed}, {Orange, DotDashed}, Orange,
Blue}
, Filling -> {3 -> Axis}
, FillingStyle -> Directive[Opacity[0.5], Blue]
, PlotRange -> {{-5, 5}, {-1, 1}}
, BaseStyle -> Black
, Epilog -> {Black, Line[{{shift, -6}, {shift, 6}}]}]
, {shift, -5, 5}]

Export[FileNameJoin[{NotebookDirectory[],
"sample_convolve.gif"}], %, "GIF", {"ImageSize" -> 700,
"AnimationRepetitions" -> \[Infinity], "AnimationRate" -> 10}]
```

<img src="sample_convolve.gif" alt="sample_convolve.gif" width="700">

### 🪼 離散データの畳み込み積分 

離散データの畳み込み積分は，次のように計算できる．

$$
(f \ast g) _j = \sum _{k=0}^{N-1} f _k g _{j-k}
$$

畳み込み積分の値は，フーリエ変換された`f`と`g`の積を逆フーリエ変換した結果と等しい．

$$
(f \ast g) _j = \mathcal{F}^{-1}[\mathcal{F}[f] \cdot \mathcal{F}[g]] _j
$$

$\mathcal{F}^{-1}[\mathcal{F}[f] \cdot \mathcal{F}[g]]$はデータ列であって，`_j`は，そのデータ列の`j`番目をとることを意味する．

全て逆フーリエ変換してから一つだけを抜き出す必要はなく，`j`だけが必要なら`j`番目のデータだけを取り出すよう和を取ればいい．

$$
(f \ast g) _j = \sum _{n=0}^{N-1} (\mathcal{F}[f] \cdot \mathcal{F}[g]) _n \exp\left(i n \frac{2\pi}{N} j\right)
$$

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

<img src="sample_descrete_convolve.gif" alt="sample_discrete_convolve.gif" width="700">

(see `example0.nb`)

[:material-microsoft-visual-studio-code:example1_convolution.cpp#L1](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Fourier/example1_convolution.cpp:1)

---
