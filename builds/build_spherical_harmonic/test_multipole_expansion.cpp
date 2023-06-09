#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include "basic_arithmetic_array_operations.hpp"

/*DOC_EXTRACT BEM

# 多重極展開(Multipole Expansion)

## Green関数の多重極展開

次のGreen関数を考える．

$$
G({\bf x},{\bf a}) = \frac{1}{\|{\bf x}-{\bf a}\|},
\quad \nabla G({\bf x},{\bf a}) = -\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}
$$

近似解 $G_{\rm apx}({\bf x},{\bf a},{\bf c})$ を以下の式で定義する：

$$
G_{\rm apx}(n, {\bf x},{\bf a},{\bf c}) &\approx \sum_{k=0}^n \sum_{m=-k}^k \left( \frac{r_{near}}{r_{far}} \right)^k \frac{1}{r_{far}} Y(k, -m, a_{near}, b_{near}) Y(k, m, a_{far}, b_{far})={\bf Y}^*({\bf x},{\bf c})\cdot{\bf Y}({\bf a},{\bf c})
$$

ここで，$(r_{near},a_{near},b_{near})$は，球面座標系に${\bf x}-{\bf c}$を変換したものであり，
$(r_{far},a_{far},b_{far})$は，球面座標系に${\bf a}-{\bf c}$を変換したもの．$Y(k, m, a, b)$は球面調和関数：

$$
Y(k, m, a, b) = \sqrt{\frac{(k - |m|)!}{(k + |m|)!}} P_k^{|m|}(\cos(a)) e^{i mb}
$$

$P_k^m(x)$はルジャンドル陪関数：

$$
P_k^m(x) = \frac{(-1)^m}{2^k k!} (1-x^2)^{m/2} \frac{d^{k+m}}{dx^{k+m}}(x^2-1)^k
$$

*/

double G(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) { return 1 / Norm(X_near - X_far); }

std::array<double, 3> gradG(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) {
   return -(X_near - X_far) / std::pow(Norm(X_near - X_far), 3.);
}

/*DOC_EXTRACT BEM

### 球面座標系への変換

${\bf x}=(x,y,z)$から球面座標$(r,a,b)$への変換は次のように行う．

$$
r = \|{\bf x}\|, \quad a = \arctan \frac{\sqrt{x^2 + y^2}}{z}, \quad b = \arctan \frac{y}{x}
$$

$r_\parallel=\sqrt{x^2+y^2}$とする．$\frac{\partial}{\partial t}(\arctan(f(t))) = \frac{f'(t)}{1 + f(t)^2}$なので，
$(r,a,b)$の$(x,y,z)$に関する勾配は次のようになる．

$$
\nabla r = \frac{\bf x}{r},\quad
\nabla a = \frac{1}{r^2r_\parallel} \left(xz,yz,-r_\parallel^2\right),\quad
\nabla b = \frac{1}{r_\parallel^2} \left(-y,x,0\right)
$$

*/

// \label{ToSphericalCoordinates}
std::array<double, 3> ToSphericalCoordinates(const std::array<double, 3>& X) {
   return {Norm(X),
           std::atan2(std::sqrt(std::get<0>(X) * std::get<0>(X) + std::get<1>(X) * std::get<1>(X)), std::get<2>(X)),
           std::atan2(std::get<1>(X), std::get<0>(X))};
};

std::array<std::array<double, 3>, 3> gradSphericalCoordinates(const std::array<double, 3>& X) {
   auto [x, y, z] = X;
   auto r = Norm(X);
   auto R = std::sqrt(x * x + y * y);
   return {std::array<double, 3>{x / r, y / r, z / r},
           std::array<double, 3>{x * z / (r * r * R), y * z / (r * r * R), -R / (r * r)},
           std::array<double, 3>{-y / (R * R), x / (R * R), 0.}};
};

/*DOC_EXTRACT BEM

### $G_{\rm apx}$の精度

${\bf c}=(x,y,0)$を変化させてプロットした結果：

| | **n=4** | **n=5** | **n=6** | **n=7** | **n=8** |
|:----:|:---:|:---:|:---:|:---:|:---:|
| **$`{\bf x} = (0,0,0),{\bf a} = (5,5,5)`$** | ![n4_A_5_5_5](./output_n4_A_5_5_5.png) | ![n5_A_5_5_5](./output_n5_A_5_5_5.png) | ![n6_A_5_5_5](./output_n6_A_5_5_5.png) | ![n7_A_5_5_5](./output_n7_A_5_5_5.png) | ![n8_A_5_5_5](./output_n8_A_5_5_5.png) |
| **$`{\bf x} = (0,0,0),{\bf a} = (10,10,10)`$** | ![n4_A_10_10_10](./output_n4_A_10_10_10.png) | ![n5_A_10_10_10](./output_n5_A_10_10_10.png)  | ![n6_A_10_10_10](./output_n6_A_10_10_10.png)  | ![n7_A_10_10_10](./output_n7_A_10_10_10.png) | ![n8_A_10_10_10](./output_n8_A_10_10_10.png) |

この結果からわかるように，Green関数の実際の値は，${\bf c}$によって変わらないが，$G_{\rm apx}$の値は${\bf c}$によって変化し，
${\bf c}$が${\bf x}$に近いところでは，$G_{\rm apx}$の値は$G$の値に近づく．

$a_{near},b_{near}$は，より小さければ精度が良く，
また，$a_{far},b_{far}$は，より大きければ精度が良くなる．

*/

// Compute the factorial of a_near given number
constexpr int factorial(const int n) {
   if (n < 0)
      throw std::runtime_error("factorial of negative number");
   if (n == 0 || n == 1)
      return 1;
   return n * factorial(n - 1);
}

double P(const double k, const double m, const double x) { return std::pow(-1., m) * std::assoc_legendre(k, m, x); };

// Compute the spherical harmonic function sph(k, m, a_far, b_far)
constexpr std::complex<double> Y(const int k, const int m, const double a_far, const double b_far) {
   if (k < 0 || abs(m) > k) {
      return std::complex<double>(0.0, 0.0);
   }
   const double assocLegendre = std::sqrt(static_cast<double>(factorial(k - abs(m))) / static_cast<double>(factorial(k + abs(m)))) * P(k, abs(m), std::cos(a_far));
   return std::polar(assocLegendre, m * b_far);
}

double Gapx(unsigned p,
            std::array<double, 3> X_near_IN,
            std::array<double, 3> X_far_IN,
            const std::array<double, 3>& center) {
   auto X_near = X_near_IN - center;
   auto X_far = X_far_IN - center;
   auto [r_near, a_near, b_near] = ToSphericalCoordinates(X_near);
   auto [r_far, a_far, b_far] = ToSphericalCoordinates(X_far);

   auto inv_r_far = 1. / r_far;
   auto r_near_inv_r = r_near / r_far;
   double c;
   std::complex<double> accum = 0;
   for (int k = 0; k <= p; ++k) {
      c = std::pow(r_near_inv_r, k) * inv_r_far;
      for (int m = -k; m <= k; ++m)
         accum += c * Y(k, -m, a_near, b_near) * Y(k, m, a_far, b_far);
   }
   return accum.real();
}

/*DOC_EXTRACT BEM

### $G_{\rm apx}$の勾配$\nabla G_{\rm apx}$の精度

$\nabla G_{\rm apx}$は，$\nabla_{\rm \circ}=(\frac{\partial}{\partial r},\frac{\partial}{\partial a},\frac{\partial}{\partial b})$とすると，

$$
\nabla G_{\rm apx} =
\nabla_{\rm \circ} G_{\rm apx}
\begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix}
$$

具体的には`gradGapx`のように

$$
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
$$

勾配の座標変換は，$Y(k,m,a_{far},b_{far})$には影響しない．

$$
\begin{align*}
\nabla G_{\rm apx}
&= \nabla_{\circ} G_{\rm apx} \begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix}\\
& = \sum_{k=0}^{n} \sum_{m=-k}^{k}\nabla_{\circ}\left(r^k Y(k, -m, a, b)\right)_{(r,a,b)=(r_{near},a_{near},b_{near})}
\begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix}
\frac{1}{r_{far}^{k+1}} Y(k, m, a_{far}, b_{far})
\end{align*}
$$

${\bf c}=(x,y,0)$を変化させてプロットした結果：

| | **n=4** | **n=5** | **n=6** | **n=7** | **n=8** |
|:----:|:---:|:---:|:---:|:---:|:---:|
| **$`{\bf x} = (0,0,0),{\bf a} = (5,5,5)`$** | ![n4_A_5_5_5](./output_n4_A_5_5_5_grad.png) | ![n5_A_5_5_5](./output_n5_A_5_5_5_grad.png) | ![n6_A_5_5_5](./output_n6_A_5_5_5_grad.png) | ![n7_A_5_5_5](./output_n7_A_5_5_5_grad.png) | ![n8_A_5_5_5](./output_n8_A_5_5_5_grad.png) |
| **$`{\bf x} = (0,0,0),{\bf a} = (10,10,10)`$** | ![n4_A_10_10_10](./output_n4_A_10_10_10_grad.png) | ![n5_A_10_10_10](./output_n5_A_10_10_10_grad.png) | ![n6_A_10_10_10](./output_n6_A_10_10_10_grad.png) | ![n7_A_10_10_10](./output_n7_A_10_10_10_grad.png) | ![n8_A_10_10_10](./output_n8_A_10_10_10_grad.png) |

*/

double dPdx(const double k, const double m, const double x) {
   double tmp = 1. / std::sqrt(1 - x * x);
   return std::pow(-1., m) * tmp * (m * x * tmp * std::assoc_legendre(k, m, x) + std::assoc_legendre(k, m + 1, x));
}

constexpr std::complex<double> dYda(const int k, const int m, const double a_far, const double b_far) {
   if (k < 0 || abs(m) > k) {
      return std::complex<double>(0.0, 0.0);
   }
   const double assocLegendre = std::sqrt(static_cast<double>(factorial(k - abs(m))) / static_cast<double>(factorial(k + abs(m)))) * dPdx(k, abs(m), std::cos(a_far)) * (-std::sin(a_far));
   return std::polar(assocLegendre, m * b_far);
}

constexpr std::complex<double> dYdb(const int k, const int m, const double a_far, const double b_far) {
   if (k < 0 || abs(m) > k) {
      return std::complex<double>(0.0, 0.0);
   }
   const double assocLegendre = std::sqrt(static_cast<double>(factorial(k - abs(m))) / static_cast<double>(factorial(k + abs(m)))) * P(k, abs(m), std::cos(a_far));
   return std::polar(assocLegendre * m, m * b_far) * std::complex<double>(0.0, 1.0);
}

std::array<double, 3> gradGapx(unsigned p,
                               std::array<double, 3> X_near_IN,
                               std::array<double, 3> X_far_IN,
                               const std::array<double, 3>& center) {
   auto X_near = X_near_IN - center;
   auto X_far = X_far_IN - center;

   auto [r_near, a_near, b_near] = ToSphericalCoordinates(X_near);
   auto [r_far, a_far, b_far] = ToSphericalCoordinates(X_far);

   auto inv_r_far = 1. / r_far;
   auto r_near_inv_r = r_near * inv_r_far;
   std::array<double, 3> grad_sphe_Gapx = {0., 0., 0.};
   std::complex<double> Y_far;
   for (int k = 0; k <= p; ++k)
      for (int m = -k; m <= k; ++m) {
         Y_far = Y(k, m, a_far, b_far);
         std::get<0>(grad_sphe_Gapx) += (k * std::pow(r_near_inv_r, k - 1) * std::pow(inv_r_far, 2.) * Y(k, -m, a_near, b_near) * Y_far).real();
         std::get<1>(grad_sphe_Gapx) += (std::pow(r_near_inv_r, k) * inv_r_far * dYda(k, -m, a_near, b_near) * Y_far).real();
         std::get<2>(grad_sphe_Gapx) += (std::pow(r_near_inv_r, k) * inv_r_far * dYdb(k, -m, a_near, b_near) * Y_far).real();
      }
   return Dot(grad_sphe_Gapx, gradSphericalCoordinates(X_near));
}

// int main() {
//    std::array<double, 3> X_near = {1, 1, 0};
//    std::array<double, 3> X_far = {100, 100, 100};
//    std::cout << "  G  ,   Gapx " << std::endl;
//    for (int i = 1; i < 40; i++) {
//       std::cout << G(X_near, X_far) << " , " << Gapx(i, X_near, X_far) << "\n";
//    }
// }

int main() {
   std::array<double, 3> A = {10, 10, 10};
   // std::array<double, 3> A = {5, 5, 5};
   std::array<double, 3> X = {0, 0, 0};

   for (int n : {4, 5, 6, 7, 8}) {
      std::ofstream ofs("output_n" + std::to_string(n) + "_A_10_10_10.txt");
      // std::ofstream ofs("output_n" + std::to_string(n) + "_A_5_5_5.txt");
      for (double x = -20.0; x <= 20.0; x += .1) {
         for (double y = -20.0; y <= 20.0; y += .1) {
            double z = 0;
            std::array<double, 3> center = {x, y, z};
            auto error = std::log10(std::abs(1 - Gapx(n, X, A, center) / G(X, A)));
            ofs << x << " " << y << " " << 0 << " " << error << "\n";
         }
         ofs << "\n";
      }
   }

   for (int n : {4, 5, 6, 7, 8}) {
      std::ofstream ofs("output_n" + std::to_string(n) + "_A_10_10_10_grad.txt");
      // std::ofstream ofs("output_n" + std::to_string(n) + "_A_5_5_5_grad.txt");
      for (double x = -20.0; x <= 20.0; x += .1) {
         for (double y = -20.0; y <= 20.0; y += .1) {
            double z = 0;
            std::array<double, 3> center = {x, y, z};
            auto g = gradG(X, A);
            auto v = Norm(g - gradGapx(n, X, A, center)) / Norm(g);
            auto error = std::log10(v);
            ofs << x << " " << y << " " << 0 << " " << error << "\n";
         }
         ofs << "\n";
      }
   }
}

/*DOC_EXTRACT BEM

## 境界要素法への応用

境界要素法で最も計算時間を要するのは，連立１次方程式の**係数行列の作成**と**それを解く**ことである．

反復法を使えば，方程式を早く解けそうだが，実際そこまで速く解けない．
その理由は，BEMの係数行列が密行列であるために，反復法で最も時間を要する行列-ベクトル積の時間が短縮できないためである．
ナイーブなBEMでは，反復解法の利点を十分に活かせない．

しかし，
多重極展開を使えば，
**BEMの係数行列をあたかも疎行列のように，行列-ベクトル積が実行でき，
反復解法を高速に実行できる．**

### 境界積分方程式

ラプラス方程式とグリーンの定理を合わせて，境界積分方程式が得られる．
これのグリーン関数$G$を多重極展開によって$G_{\rm apx}$で置き換えると，

$$
\alpha ({\bf{a}})\phi ({\bf{a}}) = \iint _\Gamma {\left( {G_{\rm apx}({\bf{x}},{\bf a},{\bf c})\phi_n ({\bf{x}}) - \phi ({\bf{x}})\nabla G_{\rm apx}({\bf{x}},{\bf a},{\bf c})\cdot {\bf{n}}(\bf x)} \right)dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t)
$$

となり，原点${\bf a}$と積分変数${\bf x}$が分離できる．

$$
\alpha ({\bf{a}})\phi ({\bf{a}})
= {\bf Y}({\bf a},{\bf c})\cdot\iint _\Gamma {\left( {{\bf Y^*}({\bf x},{\bf c})\phi_n ({\bf{x}}) - \phi ({\bf{x}}){{\bf Y}_n^*}({\bf x},{\bf c})} \right) dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
$$

ここで，${\bf Y}({\bf a},{\bf c})$は，
${\bf Y}=\{\frac{1}{r_{far}^{-k+1}}Y(0,-k,a,b),\frac{1}{r_{far}^{-k+1+1}}Y(0,-k+1,a,b),\frac{1}{r_{far}^{-k+2+1}}Y(0,-k+2,a,b),...,\frac{1}{r_{far}^{k+1}}Y(n,k,a,b)\}$
のようなベクトル．

$$
\begin{align*}
{\bf n}({\bf x})\cdot\nabla G_{\rm apx}({\bf x},{\bf a},{\bf c}) & = \sum_{k=0}^{n} \sum_{m=-k}^{k} {\bf n}({\bf x}) \cdot \left\{ \nabla_{\circ}\left(r^k Y(k, -m, a, b)\right)_{(r,a,b)=(r_{near},a_{near},b_{near})} \begin{bmatrix} \nabla r \\ \nabla a \\ \nabla b \end{bmatrix} \right\} \frac{1}{r_{far}^{k+1}} Y(k, m, a_{far}, b_{far})\\
&={\bf Y}_n^*({\bf x},{\bf c})\cdot{\bf Y}({\bf a},{\bf c})
\end{align*}
$$

ただ，十分な精度でグリーン関数を近似するためには，
$\|{\bf x - \bf c}\|$が$\|{\bf a - \bf c}\|$よりも十分に小さい必要がある．

### 空間分割

$\bf c$を一つに固定するのではなく，空間を分割して，それぞれのセルの中心において${\bf c}$を固定する．
各セルのインデックスを$\square i$として，その中心座標を${\bf c}_{\square i}$のように表す．
そうすると，

$$
\alpha ({\bf a})\phi ({\bf a})=\sum_{\square i} {\bf Y}({\bf a},{\bf c}_{\square i})\cdot\iint_{\Gamma_{\square i}}{\left( {{\bf Y^*}({\bf x},{\bf c}_{\square i})\phi_n ({\bf x}) - \phi ({\bf x}){{\bf Y}_n^*}({\bf x},{\bf c}_{\square i})} \right) dS}
$$

さらに，原点の近傍セルの積分は，多重極展開を使わずに，元々のグリーン関数を使って計算することにすると，

$$
\begin{align*}
\alpha ({\bf{a}})\phi ({\bf{a}})
=& \iint_{\Gamma_{\rm near-filed}}
\left( {G({\bf x},{\bf a})\phi_n ({\bf x}) - \phi (\bf x) G_n({\bf x},{\bf a})} \right)dS\\
& + \sum_{\square i}
\left\{
{\bf Y}({\bf a},{\bf c}_{\square i})\cdot\iint _{\Gamma _{\square i}} {\left( {{\bf Y^*}({\bf x},{\bf c}_{\square i})\phi_n ({\bf{x}}) - \phi ({\bf{x}}){{\bf Y}_n^*}({\bf x},{\bf c}_{\square i})} \right) dS}
\right\}
\end{align*}
$$

*/