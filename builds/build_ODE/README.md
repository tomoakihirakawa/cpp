# Contents

- [🐋 ODEの初期値問題](#🐋-ODEの初期値問題)
    - [⛵️ 減衰調和振動子/Damped Harmonic Oscillatorの例](#⛵️-減衰調和振動子/Damped-Harmonic-Oscillatorの例)
    - [⛵️ Runge-Kutta Integration of ODE](#⛵️-Runge-Kutta-Integration-of-ODE)


---
# 🐋 ODEの初期値問題

## ⛵️ 減衰調和振動子/Damped Harmonic Oscillatorの例

減衰調和振動子の式から，
次のような加速度$`a(x,v)=\frac{d^2x}{dt^2}`$を
[プログラム中で宣言](../../builds/build_ODE/example_DampedHrmonicOscillator.cpp#L47)し，

$$
\begin{align*}
m \frac{d^2x}{dt^2} + b \frac{dx}{dt} + k x &= 0\\
\rightarrow a(x,v) &= -\gamma v - \omega^2 x, \quad v=\frac{dx}{dt},\quad \gamma=\frac{b}{m}, \quad \omega^2=\frac{k}{m}
\end{align*}
$$

$`\gamma = 1, \omega = 10`$として，初期値問題をといてみる．
加速度の評価回数$`N`$を合わせて比較した例：

| ![](figN25.png) | ![](figN50.png) |  ![](figError.png) |
|:---:|:---:|:---:|
|$`N=25`$ evaluations|$`N=50`$ evaluations|the sum of differences|

* [後退オイラー](../../builds/build_ODE/example_DampedHrmonicOscillator.cpp#L77)の１回の計算で溜まる誤差は$`O(\Delta t^2)`$．次時刻における速度と加速度が正確に計算できなければ使えない．
* [リープフロッグ](../../builds/build_ODE/example_DampedHrmonicOscillator.cpp#L98)の１回の計算で溜まる誤差は$`O({\Delta t}^3)`$となる．[LeapFrogのクラス](../../include/integrationOfODE.hpp#L280)
* [4次のルンゲクッタ](../../builds/build_ODE/example_DampedHrmonicOscillator.cpp#L116)の１回の計算で溜まる誤差は$`O({\Delta t}^5)`$となる．しかし，加速度を4階も計算する必要がある．[RungeKuttaのクラス](../../include/integrationOfODE.hpp#L11)


** LeapFrog **

時間間隔$`\Delta t`$が変化する場合でも使える形でプログラムしている（[LeapFrogのクラス](../../include/integrationOfODE.hpp#L280)）．
$`\Delta t`$が変化する場合，"半分蹴って-移動-半分蹴って"，"半分蹴って-移動-半分蹴って"の手順を繰り返す．


[./example_DampedHrmonicOscillator.cpp#L4](./example_DampedHrmonicOscillator.cpp#L4)


**Runge-Kutta**
このように，ルンゲクッタを使って２階微分方程式を解く場合，
２階微分を２つの1階微分にわけて考え，互いに独立した２つのルンゲクッタを用意し，変数を独立して更新する必要がある．


[./example_DampedHrmonicOscillator.cpp#L117](./example_DampedHrmonicOscillator.cpp#L117)


## ⛵️ Runge-Kutta Integration of ODE

![](RK.png)


[./example_RungeKutta.cpp#L1](./example_RungeKutta.cpp#L1)


---
