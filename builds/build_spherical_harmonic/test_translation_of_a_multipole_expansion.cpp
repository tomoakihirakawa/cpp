#include <array>
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

// python3.11 ../../extract_comments.py README.md -source ./ -include ../../

/*

## `multipole_expansion`クラスのチェック

## 多重極展開とその移動

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion.cpp
make
./test_translation_of_a_multipole_expansion
```

ここで示す，多重極展開は次式を近似する．

```math
(G,\nabla G\cdot {\bf n})=\left(\frac{1}{\|{\bf x}-{\bf a}\|}, -\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}\cdot{\bf n}\right)
```

ガウス・ルジャンドル積分を使う際には，これに重みをかけて足し合わせる．その重みは別に計算し，保存しておく．元々の関数の近似の係数と重みの係数を混同しないように注意する．
現在のところ，以下のような値を与えて多重極展開を計算している．

* カーネル$G$には，位置と数値積分のための重みを与えている．
* カーネル$\nabla G\cdot {\bf n}$には，位置と数値積分のための重み，そして法線ベクトルを与えている．

```cpp
void increment(const Tddd& XIN, const std::array<double, 2> weights, const Tddd& normal) {
   const Tddd R = XIN - this->X;
   auto set_coeffs = [&](int n, int m) -> std::array<std::complex<double>, 2> {
      return {SolidHarmonicR(n, m, R) * weights[0],
               Dot(normal, Grad_SolidHarmonicR(n, m, R)) * weights[1]};
   };
   this->set(set_coeffs);
};
```

*/

int main() {
   std::array<double, 3> Xc = {0., 0., 0.};
   std::array<double, 3> Xc_ = Xc + Tddd{.5, .5, .5};
   std::array<double, 3> OL = Xc_ + Tddd{7., 7., 7.};
   std::array<double, 3> OL2 = OL + Tddd{-.5, .5, .5};
   std::array<double, 3> O = OL2 + Tddd{-.5, .5, .5};
   //
   ExpCoeffs<15> M0(Xc);
   auto X012 = std::array<Tddd, 3>{Tddd{-1., 0., -1}, Tddd{1., 0., -1}, Tddd{0., 1., -1}};
   auto q012 = std::array<double, 3>{1, 1., 1.};
   auto phi012 = std::array<double, 3>{1., 1., 1.};
   //
   double integral0 = 0.0, integral1 = 0.0;
   std::array<double, 2> weights;
   auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
   std::array<double, 3> normal = Normalize(cross);
   Tddd X, R;
   double ig = 0, ign = 0, nr, value;

   for (const auto& [xi0, xi1, ww] : __array_GW10xGW10__) {
      X = Dot(ModTriShape<3>(xi0, xi1), X012);
      value = Dot(ModTriShape<3>(xi0, xi1), q012);
      R = X - O;
      weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
      nr = Norm(R);
      //$ 直接積分
      ig += weights[0] / nr;
      ign += -weights[1] * Dot(R / (nr * nr * nr), normal);
      //$ 多重極展開

      M0.increment(X, weights, normal);
   }

   auto M2M_0 = M2M(M0, Xc_);

   auto M2M_M2L_0 = M2L(M2M_0, OL);

   auto M2M_M2L_L2L_0 = L2L(M2M_M2L_0, OL2);

   auto accuracy = [](std::complex<double> a, double b) { return (a.real() / b - 1.) * 100; };

   int n = 3, m = -2;
   for (auto x : {0., 0.5, 1.0}) {
      for (auto y : {0., 0.5, 1.0}) {
         for (auto z : {0., 0.5, 1.0}) {
            auto [rho, theta, phi] = ToSphericalCoordinates(x, y, z);
            std::cout << ToSphericalCoordinates(x, y, z) << SolidHarmonicR(n, m, std::array<double, 3>{x, y, z}) << ", " << std::polar((std::pow(rho, n) / factorial(n + m)) * std::assoc_legendre(n, m, std::cos(theta)), m * phi) << std::endl;
         }
      }
   }

   //! 　直接積分
   std::cout << Red << ig << colorReset << std::endl;
   std::cout << Green << ign << colorReset << std::endl;

   //! 　P2M後のモーメントを使った積分結果
   std::cout << "check moment expansion" << std::endl;
   std::cout << Red << M0.IG_using_M(O) << " accuracy: " << accuracy(M0.IG_using_M(O), ig) << "\%" << colorReset << std::endl;
   std::cout << Green << M0.IGn_using_M(O) << " accuracy: " << accuracy(M0.IGn_using_M(O), ign) << "\%" << colorReset << std::endl;

   //! 　M2M後のモーメントを使った積分結果
   std::cout << "check M2M" << std::endl;
   std::cout << Red << M2M_0.IG_using_M(O) << " accuracy: " << accuracy(M2M_0.IG_using_M(O), ig) << "\%" << colorReset << std::endl;
   std::cout << Green << M2M_0.IGn_using_M(O) << " accuracy: " << accuracy(M2M_0.IGn_using_M(O), ign) << "\%" << colorReset << std::endl;

   //! 　M2L後のモーメントを使った積分結果
   std::cout << "check M2L" << std::endl;
   std::cout << Red << M2M_M2L_0.IG_using_L(O) << " accuracy: " << accuracy(M2M_M2L_0.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
   std::cout << Green << M2M_M2L_0.IGn_using_L(O) << " accuracy: " << accuracy(M2M_M2L_0.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;

   //! 　L2L後のモーメントを使った積分結果
   std::cout << "check L2L" << std::endl;
   std::cout << Red << M2M_M2L_L2L_0.IG_using_L(O) << " accuracy: " << accuracy(M2M_M2L_L2L_0.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
   std::cout << Green << M2M_M2L_L2L_0.IGn_using_L(O) << " accuracy: " << accuracy(M2M_M2L_L2L_0.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;
}