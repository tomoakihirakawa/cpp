# Contents

- [🐋多重極展開(Multipole Expansion)](#🐋多重極展開(Multipole-Expansion))
    - [⛵️Green関数の多重極展開](#⛵️Green関数の多重極展開)
        - [🪸球面座標系への変換](#🪸球面座標系への変換)
        - [🪸$`G _{\rm apx}`$の精度](#🪸$`G-_{\rm-apx}`$の精度)
        - [🪸$`G _{\rm apx}`$の勾配$`\nabla G _{\rm apx}`$](#🪸$`G-_{\rm-apx}`$の勾配$`\nabla-G-_{\rm-apx}`$)
        - [🪸$`\nabla G _{\rm apx}`$の精度](#🪸$`\nabla-G-_{\rm-apx}`$の精度)
    - [⛵️境界要素法への応用](#⛵️境界要素法への応用)


---
# 🐋多重極展開(Multipole Expansion) 

## ⛵️Green関数の多重極展開 

次のGreen関数を考える．

$$
G({\bf x},{\bf a}) = \frac{1}{\|{\bf x}-{\bf a}\|},
\quad \nabla G({\bf x},{\bf a}) = -\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}
$$

近似解 $`G _{\rm apx}({\bf x- \bf c},{\bf a - \bf c})`$ を以下の式で定義する：

$$
G _{\rm apx}(n, {\bf x- \bf c},{\bf a - \bf c}) \approx \sum _{k=0}^{n} \sum _{m=-k}^{k} \left( \frac{r _{near}}{r _{far}} \right)^k \frac{1}{r _{far}} Y(k, -m, a _{near}, b _{near}) Y(k, m, a _{far}, b _{far})
$$

$$
Y(k, m, a, b) = \sqrt{\frac{(k - |m|)!}{(k + |m|)!}}(-1)^m P _k^{|m|}(\cos(a)) e^{i mb}
$$

ここで，

- $`Y(k, m, a, b)`$ は球面調和関数
- $`r _{near}`$ と $`r _{far}`$ はベクトル $`{\bf x - c}`$ と $`{\bf a - c}`$ のノルム
- $`a _{near}`$, $`b _{near}`$, $`a _{far}`$, $`b _{far}`$ はベクトル $`{\bf x - c}`$ と $`{\bf a - c}`$ の球面座標


[./test_multipole_expansion.cpp#L8](./test_multipole_expansion.cpp#L8)


### 🪸球面座標系への変換 

$`{\bf x}=(x,y,z)`$から球面座標$`(r,a,b)`$への変換は次のように行う．

$$
r = \|{\bf x}\|, \quad a = \arctan \frac{\sqrt{x^2 + y^2}}{z}, \quad b = \arctan \frac{y}{x}
$$

$`r _\parallel=\sqrt{x^2+y^2}`$とする．$`\frac{\partial}{\partial t}(\arctan(f(t))) = \frac{f'(t)}{1 + f(t)^2}`$なので，

$$
\nabla r = \frac{\bf x}{r},\quad
\nabla a = \frac{1}{r^2r _\parallel} \left(xz,yz,-r _\parallel^2\right),\quad
\nabla b = \frac{1}{r _\parallel^2} \left(-y,x,0\right)
$$


[./test_multipole_expansion.cpp#L45](./test_multipole_expansion.cpp#L45)


### 🪸$`G _{\rm apx}`$の精度 

$`{\bf c}=(x,y,0)`$を変化させてプロットした結果：

| | **n=4** | **n=5** | **n=6** | **n=7** | **n=8** |
|:----:|:---:|:---:|:---:|:---:|:---:|
| **$`{\bf x} = (0,0,0),{\bf a} = (5,5,5)`$** | ![n4_A_5_5_5](output_n4_A_5_5_5.png) | ![n5_A_5_5_5](output_n5_A_5_5_5.png) | ![n6_A_5_5_5](output_n6_A_5_5_5.png) | ![n7_A_5_5_5](output_n7_A_5_5_5.png) | ![n8_A_5_5_5](output_n8_A_5_5_5.png) |
| **$`{\bf x} = (0,0,0),{\bf a} = (10,10,10)`$** | ![n4_A_10_10_10](output_n4_A_10_10_10.png) | ![n5_A_10_10_10](output_n5_A_10_10_10.png)  | ![n6_A_10_10_10](output_n6_A_10_10_10.png)  | ![n7_A_10_10_10](output_n7_A_10_10_10.png) | ![n8_A_10_10_10](output_n8_A_10_10_10.png) |

この結果からわかるように，Green関数の実際の値は，$`{\bf c}`$によって変わらないが，$`G _{\rm apx}`$の値は$`{\bf c}`$によって変化し，
$`{\bf c}`$が$`{\bf x}`$に近いところでは，$`G _{\rm apx}`$の値は$`G`$の値に近づく．

$`a _{near},b _{near}`$は，より小さければ精度が良く，
また，$`a _{far},b _{far}`$は，より大きければ精度が良くなる．


[./test_multipole_expansion.cpp#L81](./test_multipole_expansion.cpp#L81)


### 🪸$`G _{\rm apx}`$の勾配$`\nabla G _{\rm apx}`$ 

$`\nabla G _{\rm apx}`$は，$`\nabla _{\rm \circ}=(\frac{\partial}{\partial r},\frac{\partial}{\partial a},\frac{\partial}{\partial b})`$とすると，

$$
\nabla G _{\rm apx} =
\nabla _{\rm \circ} G _{\rm apx}
\begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix}
$$

### 🪸$`\nabla G _{\rm apx}`$の精度 

$`{\bf c}=(x,y,0)`$を変化させてプロットした結果：

| | **n=4** | **n=5** | **n=6** | **n=7** | **n=8** |
|:----:|:---:|:---:|:---:|:---:|:---:|
| **$`{\bf x} = (0,0,0),{\bf a} = (5,5,5)`$** | ![n4_A_5_5_5](output_n4_A_5_5_5_grad.png) | ![n5_A_5_5_5](output_n5_A_5_5_5_grad.png) | ![n6_A_5_5_5](output_n6_A_5_5_5_grad.png) | ![n7_A_5_5_5](output_n7_A_5_5_5_grad.png) | ![n8_A_5_5_5](output_n8_A_5_5_5_grad.png) |
| **$`{\bf x} = (0,0,0),{\bf a} = (10,10,10)`$** | ![n4_A_10_10_10](output_n4_A_10_10_10_grad.png) | ![n5_A_10_10_10](output_n5_A_10_10_10_grad.png) | ![n6_A_10_10_10](output_n6_A_10_10_10_grad.png) | ![n7_A_10_10_10](output_n7_A_10_10_10_grad.png) | ![n8_A_10_10_10](output_n8_A_10_10_10_grad.png) |


[./test_multipole_expansion.cpp#L141](./test_multipole_expansion.cpp#L141)


## ⛵️境界要素法への応用 

境界要素法で最も計算時間を要するのは，連立１次方程式の**係数行列の作成**と**それを解く**ことである．

反復法を使えば，方程式を早く解けそうだが，実際そこまで速く解けない．
その理由は，BEMの係数行列が密行列であるために，反復法で最も時間を要する行列-ベクトル積の時間が短縮できないためである．
ナイーブなBEMでは，反復解法の利点を十分に活かせない．

しかし，
多重極展開を使えば，
**BEMの係数行列をあたかも疎行列のように，行列-ベクトル積が実行でき，
反復解法を高速に実行できる．**


[./test_multipole_expansion.cpp#L254](./test_multipole_expansion.cpp#L254)


---
