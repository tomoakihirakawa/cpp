#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"

auto w = std::setw(20);

/*DOC_EXTRACT 0_1_1_quasiNewton

## ロボットの節をLighthillの曲線上に乗せる

### レビュー

\insert{./REVIEW.md}

### Lighthillの式

Lighthillの式：

```math
{\bf x}^{\rm LH}(x,t) = (x,y^{\rm LH}(x,t)),\quad
y^{\rm LH}(x,t) = \left( \frac{c_1}{L} x + {c_2} \left(\frac{x}{L}\right)^2 \right) \sin \left( \frac{2 \pi}{L} x - \omega t \right)
```

ロボットの$i$番目の節の位置ベクトル：

```math
{\bf x}_{i}^{\rm rb} = {\bf x}_{i-1}^{\rm rb} + r \left( \cos \theta_i, \sin \theta_i \right)
```

ここで，変数の意味は以下の通り．

| variable | meaning |
|:---:|:---:|
| $L$ | 全長 |
| $\omega$ | 角周波数 |
| $k$ | 波数 |
| $c_1$ | 振幅1 |
| $c_2$ | 振幅2 |
| $n$ | ロボットの関節の数 |
| $r$ | ロボットの関節間の長さ |
| $\theta_i$ | $i$番目の関節が進行方向となす角度 |

### 目的関数$f$

Lighthillの式にこの節を乗せるには，どのような目的関数$f$を用いればよいだろうか．
最適化する節の一つ前の節の位置を${\bf a}=(a_x,a_y)$とすると，次の目的関数$f$が考えられる．

```math
f(\theta) = y^{\rm LH}(x,t) - a_y - r \sin \theta
```

ニュートン法には微分が必要．

```math
\frac{df}{d\theta} = -r \sin\theta\frac{d y^{\rm LH} }{dx}-r\cos\theta
```


ただ，$f$を目的関数とすると根への収束が良くなかったので，$f^2/2$を目的関数として計算する．目的関数の微分は，$f \frac{df}{d\theta}$としている．

NOTE: この目的関数$f$には，前の節の位置が含まれているが，この目的関数を使って，先頭から順番に角度を決めていけば，各最適化において見積もる角度は常に１つだけとなる．

| $n=5$ | $n=10$ | $n=50$ |
|:---:|:---:|:---:|
| ![sample5.gif](sample5.gif)  | ![sample10.gif](sample10.gif) | ![sample50.gif](sample50.gif) |

### 工夫点

愚直にニュートン法を適用すると，比較的振幅が大きい場合，正しい角度が得られない．
例えば以下のケース．

```cpp
double L = 0.71;
double w = 2. * M_PI * 1.0;
double k = 2. * M_PI * 2.0;
double c1 = 0.1;
double c2 = 0.1;
int nodes = 10;
int steps = 20;
```

そのような場合，\ref{LighthillRobot:scale}{ここ}のニュートン法のステップ幅を小さくすることで，正しい角度が得られる場合がある．


| `scale` | $n=5$ | $n=10$ | $n=50$ |
|:---:|:---:|:---:|:---:|
| `scale=1.0` | ![sample_5_bad.gif](sample_5_bad.gif)  | ![sample_10_bad.gif](sample_10_bad.gif) | ![sample_50_bad.gif](sample_50_bad.gif) |
| `scale=0.1` | ![sample_5_bad_mod.gif](sample_5_bad_mod.gif) | ![sample_10_bad_mod.gif](sample_10_bad_mod.gif) | ![sample_50_bad_mod.gif](sample_50_bad_mod.gif) |


LighthillRobotのクラスは，\ref{newton:LighthillRobot}{ここ}で宣言している．

### ロボットのエネルギー効率について

話がNewton法から離れるが，ロボットのエネルギー効率について．この内容は後で移動しておく．

ロボットの運動エネルギーは，$`\frac{1}{2}m v^2`$．
ロボットの運動エネルギーがロボットの出力だけから得られるとすると，
ロボットの出力は，このロボットの運動エネルギーの時間変化，$m v\frac{dv}{dt}$となる．
供給電力$P$は，電流$I$と電圧$V$の積$P = I V$なので，ロボットのエネルギー効率は，

```math
\eta = \frac{m v a}{I V}
```

*/

int main() {

   double L = 0.71;
   double w = 2. * M_PI * 1.0;
   double k = 2. * M_PI * 2.0;
   double c1 = 0.1;
   double c2 = 0.1;
   int nodes = 5;
   int steps = 20;

   LighthillRobot lhr(L, w, k, c1, c2, nodes);

   for (auto i = 0; i < steps; i++) {
      std::ofstream outFile("./output_lighthill/lighthill" + std::to_string(i) + ".txt");
      TimeWatch time;
      double t = (double)i / steps;
      auto Q = lhr.getAngles(t);
      auto xy = lhr.anglesToX(Q);
      std::cout << "time : " << time() << std::endl;
      for (auto j = 0; j < xy.size(); j++) {
         auto [x, y] = xy[j];
         outFile << j << " " << x << " " << y << " " << Q[j] << std::endl;
         std::cout << j << " " << x << " " << y << " " << Q[j] << std::endl;
      }
      outFile.close();
   }

   // true value
   for (auto i = 0; i < steps; i++) {
      std::ofstream outFile("./output_lighthill/lighthill" + std::to_string(i) + "_analitical.txt");
      double t = (double)i / steps;
      for (auto j = 0; j < 1000; j++) {
         double x = (double)j / 1000;
         outFile << x << " " << lhr.yLH(x, t) << std::endl;
      }
      outFile.close();
   }
}
