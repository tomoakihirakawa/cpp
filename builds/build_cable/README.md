# Contents
- [🐋 ケーブルの動的解析](#🐋-ケーブルの動的解析)
    - [⛵ 直線要素を用いたシミュレーション](#⛵-直線要素を用いたシミュレーション)
    - [⛵ 実行方法](#⛵-実行方法)
    - [⛵ ⛵ 浮体係留用に`Network`の派生クラスを作成](#⛵-⛵-浮体係留用に`Network`の派生クラスを作成)


---
# 🐋 ケーブルの動的解析 

## ⛵ 直線要素を用いたシミュレーション 

💡 弦の振動を支配する方程式として，波動方程式$`\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}`$よく紹介される．
この方程式は，ある固定した点$`x`$における弦の変位$`u`$の加速度が，弦の曲げ剛性$`c^2`$かける曲率に比例することを表している．

直線で結ばれた節点上にケーブルの自重を集中させ，その節点に働く張力や重力から，節点の運動を追っていく．

剛性は，ヤング率$`E`$と断面積$`A`$から$`EA`$．
張力$`T`$は，$`T = EA \frac{\Delta L}{L}`$となる．

オイラー法，Leap-Frog法，Runge-Kutta法を用いて，弾性体の動きをシミュレーション．

* 剛性$`[N/m]`$:$`1400 \times 10^6`$
* 減衰$`[N/(m/s^2)]`$:$`0.9`$
* 自然長$`[m]`$:$`1`$

![sample.gif](sample.gif)

`const double stiffness = 10000;`の場合
![sample_2.gif](sample_2.gif)

チェーン

## ⛵ 実行方法 

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

## ⛵ ⛵ 浮体係留用に`Network`の派生クラスを作成  

`networkLine`には，`natural_length`，`stiffness`，`damping`，`weight_per_unit_length`の4つのパラメータを持たせる．

`natural_length`は，`moorinLine`の`total_length`と`MooringLine`の`getPoints().size()`から決まる．

それをまとめる`MooringLine`は，`total_length`を持つ．

📝 内部で，自動でタイムステップを細かく取り，与えられた時間までシミュレーションするようにしたい．
[../../include/MooringLine.hpp#L6](../../include/MooringLine.hpp#L6)

[./main.cpp#L16](./main.cpp#L16)

---
