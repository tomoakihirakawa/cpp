# Contents
- [🐋 フーリエ変換](#🐋-フーリエ変換)
    - [⛵ 複素フーリエ級数展開](#⛵-複素フーリエ級数展開)
    - [⛵ 離散フーリエ変換](#⛵-離散フーリエ変換)


---
# 🐋 フーリエ変換 

## ⛵ 複素フーリエ級数展開 

```math
f(t) = \sum _{n=-\infty}^{\infty} c _n \exp(i n \omega^\ast t), \quad c _n = \frac{1}{T^\ast} \int _{-\frac{T^\ast}{2}}^{\frac{T^\ast}{2}} f(t) \exp(-i n \omega^\ast t) \, dt, \quad \omega^\ast = \frac{2\pi}{T^\ast}
```

省略

## ⛵ 離散フーリエ変換 

サンプル数が$`N+1`$，$`(k=0,1,...N+1)`$のとき，台形則を使った関数$`g(t)`$の数値積分は，

```math
\int _0^{T^\ast} g(t) dt = \left[\frac{g(0) + g(N \delta t)}{2} + \sum _{k=1}^{N-1} g(k \delta t) \right] \delta t, \quad \delta t = \frac{T^\ast}{N+1}
```

この台形則を使って，
サンプル数が$`N`$，$`(k=0,1,...N-1)`$，$`g(t) = f(t) \exp(-i n \omega^\ast t)`$として．
複素フーリエ係数$`c _n`$を数値積分で計算すると，

```math
\begin{align}
c _n &= \frac{1}{T^\ast} \left[ \frac{g(0) + g((N-1)\delta t)}{2} + \sum _{k=1}^{N-2} g(k \delta t) \right] \delta t, \quad \delta t = \frac{T^\ast}{N}, \quad g(0) = g((N-1)\delta t)\\
&= \frac{1}{N} \sum _{k=0}^{N-2} g(k \delta t)\\
&= \frac{1}{N} \sum _{k=0}^{N-2} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{T^\ast} k \frac{T^\ast}{N} \right) \right]\\
&= \frac{1}{N} \sum _{k=0}^{N-2} \left[ f\left(k\frac{T^\ast}{N}\right) \exp\left( -i n \frac{2 \pi}{N} k \right) \right]
\end{align}
```

となる．
関数$`g(t)`$が周期$`T^\ast`$の関数と仮定し$`g(0) = g((N-1)\delta t)`$とした．
また，$`\delta t`$は区間$`T^\ast`$を$`N`$等分したときの各小区間の長さである．
これが，離散フーリエ変換である．

$`c _n`$は，$`\omega _n = \frac{2 \pi n}{T^\ast}`$の角周波数成分を表す．つまり，$`f _n = \frac{n}{T^\ast}`$の周波数成分を表し．$`T _n = \frac{T^\ast}{n}`$の周期成分を表す．
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

[./main.cpp#L1](./main.cpp#L1)

---
