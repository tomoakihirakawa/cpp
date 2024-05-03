#include <array>
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT 0_2_1_translation_of_a_multipole_expansion

## ツリー構造を用いた高速多重極展開のテスト

1. まずobjファイルを読み込み`Network`クラスのインスタンスを作成
2. 三角格子をツリーの最下層のバケットに割り振る（空間分割して，バケットに保存）
3. ３次元のポテンシャル問題として，まずは各セル中心において，保存された三角形のポテンシャル面積分を多重極展開しモーメントを保存
4. M2M演算を行い，親のセルにモーメントを伝播．最上層のセルに到達するまで繰り返す
5. M2L演算を行い，最上層のセルでのモーメントを局所展開
6. L2L演算を行い，局所展開されたモーメントを子のセルに伝播
7. 最下層の局所展開係数を使って，ポテンシャルを計算し，直接積分と比較

*/

int main() {

   auto obj = new Network("./bunny.obj");
   obj->makeBucketPoints(obj->getScale() / 4.);
   obj->makeBucketFaces(obj->getScale() / 4.);

   std::array<double, 3> Xc = {0., 0., 0.};
   std::array<double, 3> Xc_ = Xc + Tddd{0., 0., 0.};
   std::array<double, 3> OL = Xc_ + Tddd{10., 10., 10.};
   std::array<double, 3> OL2 = OL + Tddd{-1., 0., 0.};
   std::array<double, 3> O = OL2 + Tddd{-1., 0., 0.};
   //
   ExpCoeffs<3> M0(Xc);
   auto X012 = std::array<Tddd, 3>{Tddd{0., 0., -1.}, Tddd{1., 0., -1.}, Tddd{0., 1., -1.}};
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
      ign += -weights[1] * Dot(R / nr, cross) / (nr * nr);
      //$ 多重極展開
      M0.increment(X, weights, normal);
   }

   auto M2M_0 = M2M(M0, Xc_);

   auto M2M_M2L_0 = M2L(M2M_0, OL);

   auto M2M_M2L_L2L_0 = L2L(M2M_M2L_0, OL2);

   auto accuracy = [](std::complex<double> a, double b) { return (a.real() / b - 1.); };

   //! 　直接積分
   std::cout << Red << ig << colorReset << std::endl;
   std::cout << Green << ign << colorReset << std::endl;

   //! 　P2M後のモーメントを使った積分結果
   std::cout << "check moment expansion" << std::endl;
   std::cout << Red << M0.IG_using_M(O) << " accuracy:" << accuracy(M0.IG_using_M(O), ig) << colorReset << std::endl;
   std::cout << Green << M0.IGn_using_M(O) << " accuracy:" << accuracy(M0.IGn_using_M(O), ign) << colorReset << std::endl;

   //! 　M2M後のモーメントを使った積分結果
   std::cout << "check M2M" << std::endl;
   std::cout << Red << M2M_0.IG_using_M(O) << " accuracy:" << accuracy(M2M_0.IG_using_M(O), ig) << colorReset << std::endl;
   std::cout << Green << M2M_0.IGn_using_M(O) << " accuracy:" << accuracy(M2M_0.IGn_using_M(O), ign) << colorReset << std::endl;

   //! 　M2L後のモーメントを使った積分結果
   std::cout << "check M2L" << std::endl;
   std::cout << Red << M2M_M2L_0.IG_using_L(O) << " accuracy:" << accuracy(M2M_M2L_0.IG_using_L(O), ig) << colorReset << std::endl;
   std::cout << Green << M2M_M2L_0.IGn_using_L(O) << " accuracy:" << accuracy(M2M_M2L_0.IGn_using_L(O), ign) << colorReset << std::endl;

   //! 　L2L後のモーメントを使った積分結果
   std::cout << "check L2L" << std::endl;
   std::cout << Red << M2M_M2L_L2L_0.IG_using_L(O) << " accuracy:" << accuracy(M2M_M2L_L2L_0.IG_using_L(O), ig) << colorReset << std::endl;
   std::cout << Green << M2M_M2L_L2L_0.IGn_using_L(O) << " accuracy:" << accuracy(M2M_M2L_L2L_0.IGn_using_L(O), ign) << colorReset << std::endl;
}