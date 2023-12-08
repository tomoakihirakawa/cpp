# Contents
- [🐋 クォータニオンを使った物体の３次元回転](#🐋-クォータニオンを使った物体の３次元回転)
    - [⛵ クォータニオンを使ったシンプルな回転](#⛵-クォータニオンを使ったシンプルな回転)
    - [⛵ クォータニオンの時間微分，角速度](#⛵-クォータニオンの時間微分，角速度)
    - [⛵ ⛵ クォータニオンの微分](#⛵-⛵-クォータニオンの微分)
    - [⛵ ⛵ 角加速度からクォータニオンの微分を計算](#⛵-⛵-角加速度からクォータニオンの微分を計算)
        - [🪼 クォータニオンの正規化](#🪼-クォータニオンの正規化)
    - [⛵ 剛体の回転と平行移動](#⛵-剛体の回転と平行移動)


---
# 🐋 クォータニオンを使った物体の３次元回転 

## ⛵ クォータニオンを使ったシンプルな回転 

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
[../../include/basic_vectors.hpp#L1290](../../include/basic_vectors.hpp#L1290)


* 以下のコードは、3Dオブジェクトの回転と平行移動を実行します．
* translate関数：ネットワークの全点を指定した量だけ平行移動します．
* rotate関数：指定したクォータニオンと中心点を使用して、ネットワークの全点を回転します．
* main関数：bunny、cow、camelオブジェクトをロードして、それぞれを回転させ、結果をファイルに出力します．

```
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=validateRotation.cpp
make
./validateRotation
```

![sample.gif](sample.gif)

$x$軸に対して回転 -> $y$軸に対して回転 -> $z$軸に対して回転

$(0.1,0,0)$を中心にして$x$軸に対して回転 -> $y$軸に対して回転 -> $z$軸に対して回転

時計回りが正である．

[./validateRotation.cpp#L5](./validateRotation.cpp#L5)

---
## ⛵ クォータニオンの時間微分，角速度 

以下を実行して，ルンゲクッタを使いクォータニオンの時間微分を時間積分する．

```
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=validateAngularVelocity.cpp
make
./validateAngularVelocity
```

## ⛵ ⛵ クォータニオンの微分  

クォータニオンは姿勢を表し，クォータニオンをかけることで姿勢を変化させることができる．

上のoperator*に定義されるように
クォータニオン同士の積$`{\boldsymbol q} _1 * {\boldsymbol q} _2`$は次のように計算される．

```math
\begin{aligned}
{\boldsymbol q} _1 * {\boldsymbol q} _2 =
\begin{bmatrix}
a1 * a2 + -b1 * b2 + -c1 * c2 + -d1 * d2\\
a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2\\
a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2\\
a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2
\end{bmatrix}
\end{aligned}
```

## ⛵ ⛵ 角加速度からクォータニオンの微分を計算  

初期値問題を解くための数値計算は，$`x _{\text next} = x + \frac{dx}{dt} dt`$
という形に対して考えられているため，これに合わせるために，クォータニオンの微分を考える必要がある．

次の計算の$`\lim _{dt \to 0}`$を取ると，角加速度からクォータニオンの微分を計算できる．

```cpp
auto nextQ = Quaternion({1., 0., 0.}, w[0] * dt) * Quaternion({0., 1., 0.}, w[1] * dt) * Quaternion({0., 0., 1.}, w[2] * dt);
return ((nextQ * q)() - q()) / dt;
```

結果として以下が得られる．

```math
\begin{aligned}
\begin{bmatrix}
q0\\
q1\\
q2\\
q3
\end{bmatrix}
=
\begin{bmatrix}
-q1 * w0 - q2 * w1 - q3 * w2\\
q0 * w0 + q3 * w1 - q2 * w2\\
-q3 * w0 + q0 * w1 + q1 * w2\\
q2 * w0 - q1 * w1 + q0 * w2
\end{bmatrix}
\end{aligned}
```

これを使えば，$`q _{\text next} = q + \frac{dq}{dt} dt`$という形で初期値問題を解くことができる．
[../../include/basic_vectors.hpp#L1513](../../include/basic_vectors.hpp#L1513)


![sample_dQdt.gif](sample_dQdt.gif)

👀 青がRK1，ピンクがRK2．緑がRK4．
RK1は時間が進むにつれて大きくなっている．
RK4も初期の状態と比べて若干大きくなっているように見える．

### 🪼 クォータニオンの正規化 

以上のような誤差は．回転行列を作成する際に，クォータニオンを正規化していないため，
回転行列の行列式が1になっていないことが原因である．
回転行列の行列式は，スケーリング（拡大縮小）を表しておリ，体積や長さを保つためには，行列式は1でなければならない．

**必ず回転行列を計算する際は，正規化したクォータニオンを使うべきである．**

## ⛵ 剛体の回転と平行移動 

剛体上の点$`X`$の位置を剛体とともに移動させるとき，
剛体の重心に関する回転と平行移動を使って，次のように計算できる．

```math
\begin{aligned}
R _{\rm new}\cdot (X-X _{\rm initial COM}) + X _{\rm new COM}
\end{aligned}
```

ここの回転行列$`R _{\rm new}`$は，「初期姿勢からの更新された姿勢までの回転」を施すものである．
初期姿勢に対する更新された姿勢を表すクォータニオン$`Q _{\rm new}`$から計算する．
[../../include/basic_vectors.hpp#L1629](../../include/basic_vectors.hpp#L1629)

[./validateAngularVelocity.cpp#L5](./validateAngularVelocity.cpp#L5)

---
