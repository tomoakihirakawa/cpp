#include <array>
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT 0_2_1_translation_of_a_multipole_expansion

## `multipole_expansion`クラスのチェック

FMMアルゴリズムでは，展開中心から遠くにある遠方原点の値は，モーメントを計算した後に渡される．
ここでチェックするのは，その計算過程を行うクラス`multipole_expansion`が問題なく動作するかどうかである．

NOTE:境界要素法におけるモーメントは，極そのものではなく，極の面積分（３D）である．

NOTE:多重極モーメントを計算するために，極の値を与えられなければならない．
\cite{Liu_2009}

NOTE:効率化するために要求されるオペレーションは，極の値が変化した際に，できるだけ少ない計算でモーメントを更新することである．

1. モーメントの計算（近傍にある複数の極を変数分離し足し合わせる）
2. 遠方の原点を決めて渡し，計算しておいたモーメントと積和を計算する
3. この計算結果と，展開しない計算結果との差をプロット

一つ前の例では，展開位置を変えることで，多重極展開の精度がどのように変化するかを調べた．
原点位置の移動による展開精度の変化は，展開中心の移動による展開精度の変化と同じである．
展開精度は，（多分）相対距離を規格化した上での，展開中心と極と原点との相対的位置関係で決まっているからである．

## 展開中心の移動（M2M）

多数の極を空間的にグループ分けして，
グループの中心位置を展開中心として多重極展開したとする．

次に，そのグループをさらにまとめて新たな多重極展開を行うことを考える．
この操作は，１ステップ目で得られた各グループの多重極展開係数を利用することで効率的に行うことができる．
各極に対する多重極展開は計算せずに済むからである．

変更されるのは，多重極係数ではなく，球面調和関数自体と，少しの係数のみである．

ここでは，始めに，１ステップ目として座標原点を中心とした多重極展開を行い，
次に，様々な場所での多重極展開を行って，前回同様に精度を検証する．

もし，２ステップ目において，展開中心が１ステップ目同様に原点であれば，
前回と同じ結果が得られるはずである．


```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion.cpp
make
./test_translation_of_a_multipole_expansion
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