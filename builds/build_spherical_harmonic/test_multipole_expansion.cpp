#include <array>
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT 0_1_spherical_harmonics

\insert{Multipole_Expansion}

## 精度の確認

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_multipole_expansion.cpp
make
./test_multipole_expansion
```

### $`G_{\rm apx}`$の精度

$`{\bf c}=(x,y,0)`$を変化させてプロットした結果：

|      | **n=4** | **n=5** | **n=6** | **n=7** | **n=8** |
|:----:|:---:|:---:|:---:|:---:|:---:|
| **x = (0,0,0), a = (5,5,5)**    | ![n4_A_5_5_5](./output_n4_A_5_5_5.png)       | ![n5_A_5_5_5](./output_n5_A_5_5_5.png)        | ![n6_A_5_5_5](./output_n6_A_5_5_5.png)        | ![n7_A_5_5_5](./output_n7_A_5_5_5.png)       | ![n8_A_5_5_5](./output_n8_A_5_5_5.png)       |
| **x = (0,0,0), a = (10,10,10)** | ![n4_A_10_10_10](./output_n4_A_10_10_10.png) | ![n5_A_10_10_10](./output_n5_A_10_10_10.png)  | ![n6_A_10_10_10](./output_n6_A_10_10_10.png)  | ![n7_A_10_10_10](./output_n7_A_10_10_10.png) | ![n8_A_10_10_10](./output_n8_A_10_10_10.png) |

この結果からわかるように，Green関数の実際の値は，$`{\bf c}`$によって変わらないが，$`G_{\rm apx}`$の値は$`{\bf c}`$によって変化し，
$`{\bf c}`$が$`{\bf x}`$に近いところでは，$`G_{\rm apx}`$の値は$`G`$の値に近づく．

$`a_{near},b_{near}`$は，より小さければ精度が良く，
また，$`a_{far},b_{far}`$は，より大きければ精度が良くなる．

### $`G_{\rm apx}`$の勾配$`\nabla G_{\rm apx}`$の精度

$`\nabla G_{\rm apx}`$は，$`\nabla_{\rm \circ}=(\frac{\partial}{\partial r},\frac{\partial}{\partial a},\frac{\partial}{\partial b})`$とすると，

```math
\nabla G_{\rm apx} =
\nabla_{\rm \circ} G_{\rm apx}
\begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix}
```

具体的には`gradGapx`のように

```math
\begin{align*}
\nabla_{\circ} G_{\rm apx}(n, {\bf x},{\bf a},{\bf c})
& = \sum_{k=0}^{n} \sum_{m=-k}^{k}\nabla_{\circ}\left(r^k Y(k, -m, a, b)\right)_{(r,a,b)=(r_{near},a_{near},b_{near})}
\frac{1}{r_{far}^{k+1}} Y(k, m, a_{far}, b_{far})\\
\nabla_{\circ}\left(r^k Y(k, -m, a, b)\right)
&= \left(k r^{k-1} Y, r^k \frac{\partial Y}{\partial a}, r^k \frac{\partial Y}{\partial b},
\right)\\
\frac{\partial Y}{\partial a} &= \sqrt{\frac{(k - |m|)!}{(k + |m|)!}} \frac{d P_k^{|m|}}{d x}(x)_{x=\cos(a) } e^{i mb}\\
\frac{\partial Y}{\partial b} &= \sqrt{\frac{(k - |m|)!}{(k + |m|)!}} P_k^{|m|}(\cos(a)) i m e^{i mb}\\
\frac{d P_k^{m}}{d x}(x) &= \frac{(-1)^m}{\sqrt{1-x^2}} \left( \frac{m x}{\sqrt{1-x^2}} P_k^{m}(x) + P_k^{m+1}(x) \right)
\end{align*}
```

勾配の座標変換は，$`Y(k,m,a_{far},b_{far})`$には影響しない．

```math
\begin{align*}
\nabla G_{\rm apx}
&= \nabla_{\circ} G_{\rm apx} \begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix}\\
& = \sum_{k=0}^{n} \sum_{m=-k}^{k}\nabla_{\circ}\left(r^k Y(k, -m, a, b)\right)_{(r,a,b)=(r_{near},a_{near},b_{near})}
\begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix}
\frac{1}{r_{far}^{k+1}} Y(k, m, a_{far}, b_{far})
\end{align*}
```

$`{\bf c}=(x,y,0)`$を変化させてプロットした結果：

| | **n=4** | **n=5** | **n=6** | **n=7** | **n=8** |
|:----:|:---:|:---:|:---:|:---:|:---:|
| **x = (0,0,0), a = (5,5,5)** | ![n4_A_5_5_5](./output_n4_A_5_5_5_grad.png) | ![n5_A_5_5_5](./output_n5_A_5_5_5_grad.png) | ![n6_A_5_5_5](./output_n6_A_5_5_5_grad.png) | ![n7_A_5_5_5](./output_n7_A_5_5_5_grad.png) | ![n8_A_5_5_5](./output_n8_A_5_5_5_grad.png) |
| **x = (0,0,0), a = (10,10,10)** | ![n4_A_10_10_10](./output_n4_A_10_10_10_grad.png) | ![n5_A_10_10_10](./output_n5_A_10_10_10_grad.png) | ![n6_A_10_10_10](./output_n6_A_10_10_10_grad.png) | ![n7_A_10_10_10](./output_n7_A_10_10_10_grad.png) | ![n8_A_10_10_10](./output_n8_A_10_10_10_grad.png) |

*/
int main() {

   for (auto a = -10.; a <= 10.; a += 1.) {
      std::array<double, 3> A = {-10, -10, 1};
      // std::array<double, 3> A = {5, 5, 5};
      std::array<double, 3> X = {10, 10, 1};

      const double dx = 0.5;

      for (int n : {4, 5, 6, 7, 8}) {
         std::string name = "./output/output_n" + std::to_string(n);
         name += "_A_" + std::to_string((int)A[0]) + "_" + std::to_string((int)A[1]) + "_" + std::to_string((int)A[2]);
         name += "_X_" + std::to_string((int)X[0]) + "_" + std::to_string((int)X[1]) + "_" + std::to_string((int)X[2]) + ".txt";
         std::cout << name << std::endl;
         std::ofstream ofs(name);
         // std::ofstream ofs("output_n" + std::to_string(n) + "_A_5_5_5.txt");
         for (double x = -20.0; x <= 20.0; x += dx) {
            for (double y = -20.0; y <= 20.0; y += dx) {
               double z = 0;
               std::array<double, 3> center = {x, y, z};
               auto error = std::log10(std::abs(1 - Gapx(n, X, A, center) / G(X, A)));
               ofs << x << " " << y << " " << 0 << " " << error << "\n";
            }
         }
      }

      for (int n : {4, 5, 6, 7, 8}) {
         std::string name = "./output/output_n" + std::to_string(n);
         name += "_A_" + std::to_string((int)A[0]) + "_" + std::to_string((int)A[1]) + "_" + std::to_string((int)A[2]);
         name += "_X_" + std::to_string((int)X[0]) + "_" + std::to_string((int)X[1]) + "_" + std::to_string((int)X[2]) + "_switch.txt";
         std::cout << name << std::endl;
         std::ofstream ofs(name);
         // std::ofstream ofs("output_n" + std::to_string(n) + "_A_5_5_5.txt");
         for (double x = -20.0; x <= 20.0; x += dx) {
            for (double y = -20.0; y <= 20.0; y += dx) {
               double z = 0;
               std::array<double, 3> center = {x, y, z};
               if (Norm(A - center) < Norm(X - center)) {
                  auto error = std::log10(std::abs(1 - Gapx(n, A, X, center) / G(X, A)));
                  ofs << x << " " << y << " " << 0 << " " << error << "\n";
               } else {
                  auto error = std::log10(std::abs(1 - Gapx(n, X, A, center) / G(X, A)));
                  ofs << x << " " << y << " " << 0 << " " << error << "\n";
               }
            }
         }
      }

      for (int n : {3, 4, 5, 6, 7, 8}) {
         std::string name = "./output/output_n" + std::to_string(n);
         name += "_A_" + std::to_string((int)A[0]) + "_" + std::to_string((int)A[1]) + "_" + std::to_string((int)A[2]);
         name += "_X_" + std::to_string((int)X[0]) + "_" + std::to_string((int)X[1]) + "_" + std::to_string((int)X[2]) + "_grad.txt";
         std::cout << name << std::endl;
         std::ofstream ofs(name);
         // std::ofstream ofs("output_n" + std::to_string(n) + "_A_5_5_5_grad.txt");
         for (double x = -20.0; x <= 20.0; x += dx) {
            for (double y = -20.0; y <= 20.0; y += dx) {
               double z = 0;
               std::array<double, 3> center = {x, y, z};
               auto g = gradG(X, A);
               auto v = Norm(g - gradGapx(n, X, A, center)) / Norm(g);
               auto error = std::log10(v);
               ofs << x << " " << y << " " << 0 << " " << error << "\n";
            }
         }
      }
   }
}

/*DOC_EXTRACT spherical_harmonics

## 境界要素法への応用

境界要素法で最も計算時間を要するのは，連立１次方程式の**係数行列の作成**と**それを解く**ことである．

反復法を使えば，方程式を早く解けそうだが，実際そこまで速く解けない．
その理由は，BEMの係数行列が密行列であるために，反復法で最も時間を要する行列-ベクトル積の時間が短縮できないためである．
ナイーブなBEMでは，反復解法の利点を十分に活かせない．

しかし，
多重極展開を使えば，
**BEMの係数行列をあたかも疎行列のように，行列-ベクトル積が実行でき，
反復解法を高速に実行できる．**

<img src="./fmm.png" width="400px">

### 境界積分方程式

ラプラス方程式とグリーンの定理を合わせて，境界積分方程式が得られる．
これのグリーン関数$G$を多重極展開によって$`G_{\rm apx}`$で置き換えると，

```math
\alpha ({\bf{a}})\phi ({\bf{a}}) = \iint _\Gamma {\left( {G_{\rm apx}({\bf{x}},{\bf a},{\bf c})\phi_n ({\bf{x}}) - \phi ({\bf{x}})\nabla G_{\rm apx}({\bf{x}},{\bf a},{\bf c})\cdot {\bf{n}}(\bf x)} \right)dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t)
```

となり，原点$`{\bf a}`$と積分変数$`{\bf x}`$が分離できる．

```math
\alpha ({\bf{a}})\phi ({\bf{a}})={\bf Y}({\bf a},{\bf c})\cdot\iint _\Gamma {\left( {{{\bf Y}^\ast}({\bf x},{\bf c})\phi_n ({\bf{x}}) - \phi ({\bf{x}}){{\bf Y}_n^\ast}({\bf x},{\bf c})} \right) dS}\quad\text{on}\quad{\bf x} \in \Gamma(t).
```

ここで，$`{\bf Y}({\bf a},{\bf c})`$は，
$`{\bf Y}=\{\frac{1}{r_{far}^{-k+1}}Y(0,-k,a,b),\frac{1}{r_{far}^{-k+1+1}}Y(0,-k+1,a,b),\frac{1}{r_{far}^{-k+2+1}}Y(0,-k+2,a,b),...,\frac{1}{r_{far}^{k+1}}Y(n,k,a,b)\}`$
のようなベクトル．

```math
{\bf n}({\bf x})\cdot\nabla G_{\rm apx}({\bf x},{\bf a},{\bf c})=\sum_{k=0}^n \sum_{m=-k}^k
{\bf n}({\bf x}) \cdot \left( \nabla_{\circ}(r^k Y(k, -m, a, b))_{(r,a,b)=(r_{near},a_{near},b_{near})}
\begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix} \right)
\frac{1}{r_{far}^{k+1}} Y(k,m,a_{far}, b_{far})={\bf Y}_n^\ast({\bf x},{\bf c})\cdot{\bf Y}({\bf a},{\bf c})
```

ただ，十分な精度でグリーン関数を近似するためには，
$`\|{\bf x - \bf c}\|`$が$`\|{\bf a - \bf c}\|`$よりも十分に小さい必要がある．

### 空間分割

$`\bf c`$を一つに固定するのではなく，空間を分割して，それぞれのセルの中心において$`{\bf c}`$を固定する．
各セルのインデックスを$`\square i`$として，その中心座標を$`{\bf c}_{\square i}`$のように表す．
そうすると，

```math
\alpha ({\bf a})\phi ({\bf a})=\sum_{\square i} {\bf Y}({\bf a},{\bf c}_{\square i})\cdot\iint_{\Gamma_{\square i}}{( {{{\bf Y}^\ast}({\bf x},{\bf c}_{\square i})\phi_n ({\bf x}) - \phi ({\bf x}){{\bf Y}_n^\ast}({\bf x},{\bf c}_{\square i})} ) dS}
```

さらに，原点の近傍セルの積分は，多重極展開を使わずに，元々のグリーン関数を使って計算することにすると，

```math
\begin{align*}
\alpha ({\bf{a}})\phi ({\bf{a}})=& \iint_{\Gamma_{\rm near-fields}}( {G({\bf x},{\bf a})\phi_n ({\bf x}) - \phi (\bf x) G_n({\bf x},{\bf a})})dS\\
& + \sum_{\square i}\{{\bf Y}({\bf a},{\bf c}_{\square i})\cdot\iint _{\Gamma _{\square i}}{({{{\bf Y}^\ast}({\bf x},{\bf c}_{\square i})\phi_n ({\bf{x}}) - \phi ({\bf{x}}){{\bf Y}_n^\ast}({\bf x},{\bf c}_{\square i})})dS}\}
\end{align*}
```

### 局所展開

Graf's Addition Theoremを使って，$`{\bf Y}^\ast({\bf x},{\bf c}_{\square i})`$を$`{\bf Y}^\ast({\bf x},{\bf c})`$の線形結合で表す．

```math
{\bf Y}^\ast({\bf x},{\bf c}_{\square i}) = \sum_{\square j} {\bf Y}^\ast({\bf x},{\bf c}_{\square j}){\bf Y}({\bf c}_{\square j},{\bf c}_{\square i})
```

*/

// ### モーメントの移動

// ```math
// \iint _{\Gamma _{\square i}} {{{\bf Y}^\ast}({\bf x},{\bf c}_{\square i})\phi_n ({\bf{x}})}dS
// ```

// は，近似の展開中心を$`\square i`$番目のセルの中心$`{\bf c}_{\square i}`$した場合のモーメントである．
// $`{\bf c}_{\square j}`$を中心としたモーメントを，$`{\bf c}_{\square i}`$におけるモーメントを使って表す．

// ```math
// \iint _{\Gamma _{\square j}} {{{\bf Y}^\ast}({\bf x},{\bf c}_{\square j})\phi_n ({\bf{x}})}dS
// = \iint _{\Gamma _{\square j}} {{{\bf Y}^\ast}({\bf x},{\bf c}_{\square j})\phi_n ({\bf{x}})}dS
// ```

// これは，$`\phi_n`$と$`\phi`$の積のモーメントである．
