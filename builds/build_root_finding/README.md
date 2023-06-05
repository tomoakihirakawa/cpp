# Contents

- [🐋ニュートン法](#🐋ニュートン法)
    - [⛵️ニュートン法](#⛵️ニュートン法)
    - [⛵️例）ロボットの節をLightHillの曲線上に乗せる](#⛵️例）ロボットの節をLightHillの曲線上に乗せる)
        - [🪸目的関数$`f`$](#🪸目的関数$`f`$)
    - [⛵️準ニュートン法](#⛵️準ニュートン法)


---
# 🐋ニュートン法 

## ⛵️ニュートン法 

**最適か否かを判断するための関数**（目的関数）のヤコビ行列を使う場合とヘッセ行列を使う場合がある．
目的関数の根を見つける場合は，ヤコビ行列を使う．
最適化の問題の多くは，目的関数の最大最小を求めることなので，ヘッセ行列を利用したニュートン法を用いる．


[./example0_NewtonRaphson_0.cpp#L1](./example0_NewtonRaphson_0.cpp#L1)


## ⛵️例）ロボットの節をLightHillの曲線上に乗せる 

LightHillの式：

$$
{\bf x}^{\rm LH}(x,t) = (x,y^{\rm LH}(x,t)),\quad
y^{\rm LH}(x,t) = \left( \frac{c _1}{L} x + {c _2} \left(\frac{x}{L}\right)^2 \right) \sin \left( \frac{2 \pi}{L} x - \omega t \right)
$$

ロボットの$`i`$番目の節の位置：

$$
{\bf x} _{i}^{\rm rb} = {\bf x} _{i-1}^{\rm rb} + r \left( \cos \theta _i, \sin \theta _i \right)
$$

ここで，変数の意味は以下の通り．

| variable | meaning |
|:---:|:---:|
| $`L`$ | 全長 |
| $`\omega`$ | 角周波数 |
| $`k`$ | 波数 |
| $`c _1`$ | 振幅1 |
| $`c _2`$ | 振幅2 |
| $`n`$ | ロボットの関節の数 |
| $`r`$ | ロボットの関節間の長さ |
| $`\theta _i`$ | $`i`$番目の関節が進行方向となす角度 |

### 🪸目的関数$`f`$ 

LightHillの式にこの節を乗せるには，どのような目的関数$`f`$を用いればよいだろうか．
最適化する節の一つ前の節の位置を$`{\bf a}=(a _x,a _y)`$とすると，次の目的関数$`f`$が考えられる．

$$
f(\theta) = y^{\rm LH}(x,t) - a _y - r \sin \theta
$$

ニュートン法には微分が必要．

$$
\frac{df}{d\theta} = -r \sin\theta\frac{d y^{\rm LH} }{dx}-r\cos\theta
$$

💡 この目的関数$`f`$には，前の節の位置が含まれているが，この目的関数を使って，先頭から順番に角度を決めていけば，各最適化において見積もる角度は常に１つだけとなる．

| $`n=5`$ | $`n=10`$ | $`n=50`$ |
|:---:|:---:|:---:|
| ![sample5.gif](sample5.gif)  | ![sample10.gif](sample10.gif) | ![sample50.gif](sample50.gif) |

💡 ただし，$`f`$を目的関数とすると根への収束が良くなかったので，$`f^2/2`$を目的関数として計算した．目的関数の微分は，$`f \frac{df}{d\theta}`$としている．


[./example0_NewtonRaphson_1.cpp#L6](./example0_NewtonRaphson_1.cpp#L6)


## ⛵️準ニュートン法 

ニュートン法で使うヤコビ行列などを別のものに置き換えた方法．


[./example1_Broyden.cpp#L1](./example1_Broyden.cpp#L1)


---
