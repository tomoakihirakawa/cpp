# Contents

    - [⛵️ロボットの節の位置をLightHillの曲線上に乗せる](#⛵️ロボットの節の位置をLightHillの曲線上に乗せる)
    - [⛵️ニュートン法](#⛵️ニュートン法)
    - [⛵️準ニュートン法](#⛵️準ニュートン法)


---
## ⛵️ロボットの節の位置をLightHillの曲線上に乗せる 

LightHillの式は以下のようになる．

$$
{\bf x}^{\rm LH}(x,t) = (x,y^{\rm LH}(x,t)),\quad
y^{\rm LH}(x,t) = \left( \frac{c _1}{L} x + {c _2} \left(\frac{x}{L}\right)^2 \right) \sin \left( \frac{2 \pi}{L} x - \omega t \right)
$$

ここで，$`c _1, c _2, L, \omega`$は定数である．

ロボットの$`i`$番目の節の位置は，$`{\bf x} _{i}^{\rm rb} = {\bf x} _{i-1}^{\rm rb} + r \left( \cos \theta _i, \sin \theta _i \right)`$である．
次の関数を使って表すことにする．

$$
{\bf x} _{i}^{\rm rb} = X^{\rm rb}({\bf x} _{i-1}^{\rm rb},r,\theta _i),
\quad X^{\rm rb}({\bf a},r,\theta) = {\bf a} + r \left( \cos \theta, \sin \theta \right)
$$

ここで，$`r`$はロボットの節の長さ，$`\theta _i`$はロボットの節の角度である．
頭の位置を$`X _{0}^{\rm rb}=(0,0)`$とする．
次の節の位置は，$`X _{1}^{\rm rb} = r (\cos \theta _1, \sin \theta _1)`$である．
さらに次の節の位置は，$`X _{2}^{\rm rb} = X _{1}^{\rm rb} + r (\cos \theta _2, \sin \theta _2)`$となる．

目的関数$`f`$の微分は，

$$
\frac{df}{d\theta} = -r \sin\theta\frac{d y^{\rm LH} }{dx}-r\cos\theta
$$

💡 この目的関数$`f`$には，前の節の位置を与える必要がある．節の位置は，後ろの節の位置によって変わらないので，この目的関数を先頭から順番に最適化することは問題ない．

![./output_lighthill/sample.gif](output_lighthill/sample.gif)


[./example01_NewtonRaphson.cpp#L73](./example01_NewtonRaphson.cpp#L73)


## ⛵️ニュートン法 

**最適か否かを判断するための関数**（目的関数）のヤコビ行列を使う場合とヘッセ行列を使う場合がある．
目的関数の根を見つける場合は，ヤコビ行列を使う．
最適化の問題の多くは，目的関数の最大最小を求めることなので，ヘッセ行列を利用したニュートン法を用いる．


[./example0_NewtonRaphson.cpp#L1](./example0_NewtonRaphson.cpp#L1)


## ⛵️準ニュートン法 

ニュートン法で使うヤコビ行列などを別のものに置き換えた方法．


[./example1_Broyden.cpp#L1](./example1_Broyden.cpp#L1)


---
