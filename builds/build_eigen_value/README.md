# Contents

- [🐋 固有値問題](#🐋-固有値問題)
    - [⛵ 固有値の計算](#⛵-固有値の計算)


---
# 🐋 固有値問題 

## ⛵ 固有値の計算 

行列$`A`$をQR分解$`A=QR`$し，
$`A _k = Q _k^{-1} A Q _k`$の計算を繰り返すことで，
$`A _k`$の対角成分が$`A`$の固有値に収束することを確認する．

わかりやすいように$`\cdot`$で行列の積を表す．

```math
Q _2^{-1} \cdot (Q _1^{-1} \cdot (Q _0^{-1} \cdot (A = Q _0R _0) \cdot Q _0=Q _1R _1) \cdot Q _1=Q _2 \cdot R _2) \cdot R _2
```


[./testEigenValues.cpp#L11](./testEigenValues.cpp#L11)


---
