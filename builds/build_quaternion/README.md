# Contents

- [🐋 クォータニオンによる回転と平行移動](#🐋-クォータニオンによる回転と平行移動)


---
# 🐋 クォータニオンによる回転と平行移動 

クォータニオンは，3D回転を効率的に計算するために便利な表現．

ある回転軸$v$に対して，角度$\theta$だけ回転するクォータニオン$q$は以下のように表される．

```math
\begin{aligned}
q &= a + bi + cj + dk = \cos(\theta/2) +  \sin(\theta/2) \cdot \dfrac{v _x i + v _y k + v _z k}{\|v\|}\\
v &= (v _x, v _y, v _z)
\end{aligned}
```

`Rv()` a rotation transformation in the global coordinate system.

```math
Rv = \begin{bmatrix}
a^2 + b^2 - c^2 - d^2 & 2 \cdot b \cdot c - 2 \cdot a \cdot d & 2 \cdot a \cdot c + 2 \cdot b \cdot d \\
2 \cdot b \cdot c + 2 \cdot a \cdot d & a^2 - b^2 + c^2 - d^2 & -2 \cdot a \cdot b + 2 \cdot c \cdot d \\
-2 \cdot a \cdot c + 2 \cdot b \cdot d & 2 \cdot a \cdot b + 2 \cdot c \cdot d & a^2 - b^2 - c^2 + d^2 \\
\end{bmatrix}
```

The `Rs()` function looks like it's used to calculate how the global coordinates move relative to the object's coordinates, through rotation.
The matrix representation is:

```math
Rs = \begin{bmatrix}
a^2 + b^2 - c^2 - d^2 & 2 \cdot b \cdot c + 2 \cdot a \cdot d & -2 \cdot a \cdot c + 2 \cdot b \cdot d \\
2 \cdot b \cdot c - 2 \cdot a \cdot d & a^2 - b^2 + c^2 - d^2 & 2 \cdot a \cdot b + 2 \cdot c \cdot d \\
2 \cdot a \cdot c + 2 \cdot b \cdot d & -2 \cdot a \cdot b + 2 \cdot c \cdot d & a^2 - b^2 - c^2 + d^2 \\
\end{bmatrix}
```

[../../include/basic_vectors.hpp#L1345](../../include/basic_vectors.hpp#L1345)



* 以下のコードは、3Dオブジェクトの回転と平行移動を実行します．
* translate関数：ネットワークの全点を指定した量だけ平行移動します．
* rotate関数：指定したクォータニオンと中心点を使用して、ネットワークの全点を回転します．
* main関数：bunny、cow、camelオブジェクトをロードして、それぞれを回転させ、結果をファイルに出力します．

![sample.gif](sample.gif)

$x$軸に対して回転 -> $y$軸に対して回転 -> $z$軸に対して回転

$(0.1,0,0)$を中心にして$x$軸に対して回転 -> $y$軸に対して回転 -> $z$軸に対して回転

時計回りが正である．


[./main.cpp#L5](./main.cpp#L5)


---
