# Contents
- [🐋 数値積分](#🐋-数値積分)
    - [⛵ 台形則](#⛵-台形則)
    - [⛵ ルジャンドル多項式，ルジャンドル補間，ガウス・ルジャンドル積分](#⛵-ルジャンドル多項式，ルジャンドル補間，ガウス・ルジャンドル積分)


---
# 🐋 数値積分 

## ⛵ 台形則 

台形則は，関数の積分を台形の面積の和で近似する方法である．

```math
\int _a^b f(x) dx \approx \left(\frac{f(a)+f(b)}{2} + \sum _{i=1}^{N-1} f(a+i\Delta x)\right)\Delta x, \quad \Delta x = \frac{b-a}{N}
```

[./TrapezoidalRule.cpp#L4](./TrapezoidalRule.cpp#L4)

---
## ⛵ ルジャンドル多項式，ルジャンドル補間，ガウス・ルジャンドル積分 

🚧 このファイルは，ルジャンドル多項式，ルジャンドル補間，ガウス・ルジャンドル積分を扱う．

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=LegendrePolynomials.cpp
make
./LegendrePolynomials > LegendrePolynomials.dat
```

```sh
gnuplot
file = 'LegendrePolynomials.dat'
plot for [i=2:7] file using 1:i title sprintf('order %d', i-2)
```

[./LegendrePolynomials.cpp#L19](./LegendrePolynomials.cpp#L19)

---
