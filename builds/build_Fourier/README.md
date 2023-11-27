# Contents
- [🐋 フーリエ変換](#🐋-フーリエ変換)
    - [⛵ 複素フーリエ級数展開](#⛵-複素フーリエ級数展開)
    - [⛵ 離散フーリエ変換](#⛵-離散フーリエ変換)
    - [⛵ 逆離散フーリエ変換](#⛵-逆離散フーリエ変換)


---
# 🐋 フーリエ変換 

## ⛵ 複素フーリエ級数展開 

```math
f(t) = \sum _{n=-\infty}^{\infty} c _n \exp(i n \omega^\ast t), \quad c _n = \frac{1}{T^\ast} \int _{-\frac{T^\ast}{2}}^{\frac{T^\ast}{2}} f(t) \exp(-i n \omega^\ast t) \, dt, \quad \omega^\ast = \frac{2\pi}{T^\ast}
```

$`\exp({i \theta}) = \cos \theta + i \sin \theta`$なので，
フーリエ係数の実部には，$`\cos \theta`$の係数が，虚部には，$`\sin \theta`$の係数が含まれる．

$`c _n=\frac{a _n - i \mathrm{sgn}(n) b _n}{2}`$

## ⛵ 離散フーリエ変換 

関数が周期関数であるとする．多くの場合その周期つまり基本周期はわかるだろうから，その周期を$`T^\ast`$とする．
その周期の積分区間で積分することで，フーリエ係数を求めことができる．

元のデータが$`f(0)=f(T^\ast)`$を満たす周期関数であることを考慮して，有効なサンプリングを$`N`$回行うなら
（１周期区間にびっしりとサンプリングするだろう）サンプリング時刻は$`t=k\delta t=k\frac{T^\ast}{N},k=0,1,2,...,N-1`$とすべき（等間隔とする）．

台形則を使って，
サンプル数が$`N`$，$`(k=0,1,...N-1)`$，$`g _n(t) = f(t) \exp(-i n \omega^\ast t)`$として．
複素フーリエ係数$`c _n`$を数値積分で計算すると，

```math
\begin{align}
c _n &= \frac{1}{T^\ast} \left[ \frac{g _n(0) + g _n((N-1)\delta t)}{2} + \sum _{k=1}^{N-2} g _n(k \delta t) \right] \delta t, \quad \delta t = \frac{T^\ast}{N}, \quad g _n(0) = g _n((N-1)\delta t)\\
&= \frac{1}{N} \sum _{k=0}^{N-2} g _n(k \delta t)\\
&= \frac{1}{N} \sum _{k=0}^{N-2} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{T^\ast} k \frac{T^\ast}{N} \right) \right]\\
&= \frac{1}{N} \sum _{k=0}^{N-2} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{N} k \right) \right]
\end{align}
```

となる．

この係数$`c _n`$を使って，元の関数$`f(t)`$を表すと

```math
f(t) = \sum _{n=-\infty}^{\infty} c _n \exp(i n \omega^\ast t)
```

関数$`g _n(t)`$が周期$`T^\ast`$の関数と仮定し$`g _n(0) = g _n((N-1)\delta t)`$とした．
また，$`\delta t`$は区間$`T^\ast`$を$`N`$等分したときの各小区間の長さである．
これが，離散フーリエ変換である．

💡 $`c _n`$は，$`\omega _n = \frac{2 \pi n}{T^\ast}`$の角周波数成分を表す．つまり，$`f _n = \frac{n}{T^\ast}`$の周波数成分を表し．$`T _n = \frac{T^\ast}{n}`$の周期成分を表す．
このことから，周波数分解能は$`\Delta f = \frac{1}{T^\ast}`$，周期分解能は$`\Delta T = T^\ast`$とわかる．
また，基本周期$`T^\ast`$と同じ周期でサンプリングしていては，$`T^\ast`$の周期成分は分解できない．
少なくとも，$`T^\ast`$の周期成分を分解するには，$`T^\ast`$の半分以下の周期でサンプリングする必要がある．

💡 $`c _n`$を変形すると，

$$
c _n = \frac{1}{N} \sum _{k=0}^{N-2} \left[ f\left(k\frac{T^\ast}{N}\right)
\exp\left(k\right)\right]\exp\left( -i n \frac{2 \pi}{N}\right)
$$

となっており，$`n`$に関して周期$`N`$の周期関数となっている．
また$`\cos(\theta)=\cos(-\theta)`$であるため，$`\Re[c _n]=\Re[c _{-n}]`$で
$`\sin(\theta)=-\sin(-\theta)`$であるため，$`\Im[c _n]=-\Im[c _{-n}]`$である．
１周期分つまり$`N`$分の係数ではなく，$`N/2`$分の係数さえわかれば元の関数を復元できる．

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

## ⛵ 逆離散フーリエ変換

[./example0_simple.cpp#L1](./example0_simple.cpp#L1)

---
