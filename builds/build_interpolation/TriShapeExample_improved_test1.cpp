/*DOC_EXTRACT 1_3_0_interpolation

## 接続関係を利用した補間精度の向上（擬2次補間）

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=TriShapeExample_improved_test1.cpp
make
./TriShapeExample_improved_test1
```

2次補間を利用する，要素は，2次要素と呼ばれ，
一般的には，三角形の頂点に加え，辺上にもサンプル点を配置する．

<img src="pseudo_quad.png" width="700">


ここで紹介する擬2次補間要素は，辺上に存在しない節点を周辺の要素を使って近似し，その値を使って三角形上に2次要素を作る方法である．
擬2次補間は，線形要素と同じメッシュを使ったとしても節点数を増やす必要はないが，
メッシュの接続関係を線形補間よりも多く考慮しており，また高次の補間であるため，精度の向上が期待できる．

<!-- この方法の実装には，要素同士の接続情報を計算中に効率的に取得する必要がある． -->

辺上の節点の補間は，2次補間を使って行う．
この辺に隣接する三角形を中央にもつ2次補間は２通り考えられ，この２通りの補間の平均値を辺上の節点の値とする．
一度，辺上の節点の値が決まれば，擬2次補間要素は一般的な2次補間と全く同じである．

ただし，次の章で示すが，擬2次補間要素を方程式の離散化に適用するためには，
一方的に値を補間する機能だけでなく，補間された値がどの節点の値の線形結合で決まるかが取得できる機能もプログラムに実装する必要がある．

擬2次補間は，よく知られている2次補間の形状関数を基本としている．

```math
\begin{aligned}
v({\boldsymbol \xi}) =
N({\xi_0,\xi_1})^{\intercal}
V\end{aligned}
,\quad
N({\xi_0,\xi_1})=\left(
\begin{array}{c}
\xi _0 (2 \xi _0 - 1)\\
\xi _1 (2 \xi _1 - 1)\\
\xi _2 (2 \xi _2 - 1)\\
4 \xi _0 \xi _1\\
4 \xi _1 \xi _2\\
4 \xi _2 \xi _0
\end{array}
\right),\quad
V=\left(
\begin{array}{c}
v_0\\v_1\\v_2\\v_3\\v_4\\v_5
\end{array}
\right)
```

Fig. \ref{fig:pseudo_quad_schematic}に示すように，
この形状関数の係数を，対応する節点の値に掛けて足し合わせることで，
三角形要素の内部の任意の点における値を補間することができる．

<img src="pseudo_quad_white.png" width="700">

ただし，辺上の節点$3,4,5$は設定していないので，
隣接する三角形の頂点の値を使った2次補間の平均で近似する：

```math
\begin{aligned}
v({\boldsymbol \xi}) =
N({\xi_0,\xi_1})^{\intercal}
\left(
\begin{array}{c}
v_0\\v_1\\v_2\\
\frac{1}{2}\left({N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right) V_{\rm {\ell}01in} + N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right) V_{\rm {\ell}01out}}\right)\\
\frac{1}{2}\left({N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right) V_{\rm {\ell}12in} + N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right) V_{\rm {\ell}12out}}\right)\\
\frac{1}{2}\left({N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right) V_{\rm {\ell}20in} + N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right) V_{\rm {\ell}20out}}\right)
\end{array}
\right)
\end{aligned}
```

この式を\eqref{eq:general_usage_of_shape_function}の形に書き直すために．
次のような関係を使う：

```math
\begin{aligned}
N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right) V_{\rm {\ell}12in} &=N_{\rm q}\left(\frac{1}{4},\frac{1}{2}\right) V_{\rm {\ell}01in},\\
N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right) V_{\rm {\ell}20in} &=N_{\rm q}\left(\frac{1}{2},\frac{1}{4}\right) V_{\rm {\ell}01in},\\
N({\xi_0,\xi_1})^{\intercal}
\left(
\begin{array}{c}
v_0\\v_1\\v_2\\0\\0\\0
\end{array}
\right)
&=
\left(
\begin{array}{c}
0\\0\\0\\\xi_2 (2\xi_2 - 1)\\\xi_0 (2\xi_0 - 1)\\\xi_1 (2\xi_1 - 1)
\end{array}
\right)
V_{\rm {\ell}01in}
\end{aligned}
```

これを使って，$V_{\rm {\ell}12in}$と$V_{\rm {\ell}21in}$の代わりに，$V_{\rm {\ell}01in}$を使った式に置き換える．

```math
\begin{aligned}
{\bf x}({\boldsymbol \xi})&=
\left(
\left(\begin{array}{c}
0\\0\\0\\\xi_2 (2\xi_2 - 1)\\\xi_0 (2\xi_0 - 1)\\\xi_1 (2\xi_1 - 1)\\
\end{array}
\right)
+2 \xi_0 \xi_1 N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right)
+2 \xi_1 \xi_2 N_{\rm q}\left(\frac{1}{2},\frac{1}{4}\right)
+2 \xi_2 \xi_0 N_{\rm q}\left(\frac{1}{4},\frac{1}{2}\right)
\right)V_{\rm {\ell}01in}\\
&+2 \xi_0 \xi_1 N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right)V_{\rm {\ell}01out}\\
&+2 \xi_1 \xi_2 N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right)V_{\rm {\ell}12out}\\
&+2 \xi_2 \xi_0 N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right)V_{\rm {\ell}21out}
\end{aligned}
\label{eq:pseudo_quadratic_interpolation}
```

このように，\eqref{eq:general_usage_of_shape_function}の形の式を４つ足し合わせることで，擬2次補間を実装することができる．

ただし，辺が角を成している場合，この補間では，角にはならず，滑らかに補間されてしまう．
そのため，角を成している辺上の節点は，線形補間を使って近似することにする．つまり，辺が繋ぐ２節点の平均で近似する．
例えば，辺01が角となっている場合，\eqref{eq:pseudo_quadratic_interpolation}は次のように書き換える．

```math
\begin{aligned}
{\bf x}({\boldsymbol \xi})&=
\left(
\left(\begin{array}{c}
0\\0\\0\\\xi_2 (2\xi_2 - 1)\\\xi_0 (2\xi_0 - 1)\\\xi_1 (2\xi_1 - 1)\\
\end{array}
\right)
+
2\xi_0\xi_1
\left(\begin{array}{c}
0\\0\\0\\0\\1\\1
\end{array}
\right)
+2 \xi_1 \xi_2 N_{\rm q}\left(\frac{1}{2},\frac{1}{4}\right)
+2 \xi_2 \xi_0 N_{\rm q}\left(\frac{1}{4},\frac{1}{2}\right)
\right)V_{\rm {\ell}01in}\\
&+2 \xi_1 \xi_2 N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right)V_{\rm {\ell}12out}\\
&+2 \xi_2 \xi_0 N_{\rm q}\left(\frac{1}{4},\frac{1}{4}\right)V_{\rm {\ell}21out}
\end{aligned}
```

0,1節点は，$`V_{\rm {\ell}01in}`$における4,5節点であるため，$`2\xi_0\xi_1(0,0,0,1,1)^{\intercal}`$の項に

<img src="peak_function_interpolations.png">

積分結果を比較すると，線形補間よりも擬2次補間の方が精度が向上していることがわかる．
規則的なメッシュでは，２次補間と擬2次補間の結果は同程度であるが，
不規則なメッシュでは，擬2次補間の精度は２次補間により悪いが，線形補間よりも良い．

<img src="peak_function_interpolation_integration.png">

* [pseudo_quad要素の精度peaks_function積分.nb](pseudo_quad要素の精度peaks_function積分.nb)
* [説明.key](説明.key)

*/

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Network.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "kernelFunctions.hpp"

bool checker(const networkLine *line) { return useOppositeFace(line, 0.7); };

double peaksFunction(double x, double y) {
   return 3 * (1 - x) * (1 - x) * std::exp(-x * x - (y + 1) * (y + 1)) -
          10 * (x / 5 - x * x * x - y * y * y * y * y) * std::exp(-x * x - y * y) -
          1. / 3 * std::exp(-(x + 1) * (x + 1) - y * y);
}

double D_peaksFunction_dx(double x, double y) {
   return -(1.0 / 3.0) * std::exp((-1 - x) * (1 + x) - y * y) * (-2 - 2 * x) - 6 * std::exp(-x * x - (1 + y) * (1 + y)) * (1 - x) - 6 * std::exp(-x * x - (1 + y) * (1 + y)) * (1 - x) * (1 - x) * x - 10 * std::exp(-x * x - y * y) * (1.0 / 5.0 - 3 * x * x) + 20 * std::exp(-x * x - y * y) * x * (x / 5.0 - x * x * x - std::pow(y, 5));
}

double D_peaksFunction_dy(double x, double y) {
   return 2.0 / 3.0 * std::exp((-1 - x) * (1 + x) - y * y) * y + 50 * std::exp(-x * x - y * y) * std::pow(y, 4) - 6 * std::exp(-x * x - (1 + y) * (1 + y)) * (1 - x) * (1 - x) * (1 + y) + 20 * std::exp(-x * x - y * y) * y * (x / 5.0 - x * x * x - std::pow(y, 5));
}

double exact_norm_cross(double x, double y) {
   return Norm(Cross(Tddd{1., 0., D_peaksFunction_dx(x, y)}, Tddd{0., 1., D_peaksFunction_dy(x, y)}));
}
int main() {

   auto get_triangle_verticies = [](const int divide) {
      T3Tdd domain = {{{3., -3.}, {-3., 3.}, {-3., -3.}}};
      // T3Tdd domain = {{{1., -1.}, {-1., 1.}, {-1., -1.}}};
      auto subdivide = SubdivideSquareIntoTriangles(divide);
      std::vector<T3Tddd> triangle_vertices(subdivide.size());
      int k = 0;
      T3Tddd tmp;
      for (auto t01 : subdivide) {
         int i = 0;
         for (auto [t0, t1] : t01) {
            auto [x, y] = Dot(TriShape<3>(t0, t1), domain);
            tmp[i++] = {x, y, 0.};
         }
         triangle_vertices[k++] = tmp;
      }
      return triangle_vertices;
   };

   auto DistortMesh = [](Network *net) {
      auto points = ToVector(net->getPoints());
      auto lines = ToVector(net->getLines());
      for (int k = 0; k < 5; ++k) {
         lines = RandomSample(lines);
         if (k % 3 == 0)
            for (int j = 0; auto &l : lines)
               if (j++ % 3 == 0)
                  if (l->canFlip())
                     l->flip();

         if (k % 3 == 1)
            for (const auto &l : lines)
               l->flipIfTopologicallyBetter(10.0, 10.0, 4);

         for (const auto &l : lines)
            l->flipIfBetter(10.0, 10.0, 4);

         for (int j = 0; j < 10; ++j) {
#pragma omp parallel
            for (const auto &p : points)
#pragma omp nowait
            {
               if (std::ranges::none_of(p->getLines(), [](const auto &l) { return l->getFaces().size() == 1; })) {
                  if (j < 5)
                     p->setX(p->X + 0.5 * NeighborAverageSmoothingVector(p, p->X));
                  else
                     p->setX(p->X + 0.2 * DistorsionMeasureWeightedSmoothingVector(p, p->X));
               }
            }
            net->setGeometricProperties();
         }
      }
   };

   Tddd X;

   for (const auto use_distorted_mesh : {true, false}) {

      std::string id = use_distorted_mesh ? "_distorted" : "_regular";

      //@ -------------------------------------------------------------------------- */
      //@                             peaksFunctionを使った例                          */
      //@ -------------------------------------------------------------------------- */
      for (auto divide : {10, 20, 30}) {
         auto name = _HOME_DIR_ + "/output/peaksFunction_linear" + id + std::to_string(divide) + ".dat";
         std::ofstream file_peaks_linear(name);
         std::cout << "file_peaks_linear : " << name << std::endl;
         name = _HOME_DIR_ + "/output/peaksFunction_pseudo_quad" + id + std::to_string(divide) + ".dat";
         std::ofstream file_peaks_pseudo_quad(name);
         std::cout << "file_peaks_pseudo_quad : " << name << std::endl;
         /* ------------------------------------------------- */
         auto net_peaks = new Network();
         net_peaks->setFaces(get_triangle_verticies(divide));
         if (use_distorted_mesh)
            DistortMesh(net_peaks);
         for (const auto &f : net_peaks->getFaces()) {
            for (const auto &p : f->getPoints()) {
               auto [x, y, z] = p->X;
               p->setX(Tddd{x, y, peaksFunction(x, y)});
            }
         }
         net_peaks->setGeometricProperties();
         /* ------------------------------------------------- */
         for (const auto &f : net_peaks->getFaces()) {
            for (const auto &p : f->getPoints()) {
               auto [x, y, z] = p->X;
               file_peaks_linear << x << " " << y << " " << peaksFunction(x, y) << " ";
               // DodecaPoints dodecapoints(f, p);
               for (auto t0t1 : SymmetricSubdivisionOfTriangle_00_10_01(6)) {
                  for (auto [t0, t1] : t0t1) {
                     auto [x, y, z] = f->dodecaPoints[0]->interpolate(t0, t1, [&](networkPoint *p) { return p->X; });
                     file_peaks_pseudo_quad << x << " " << y << " " << z << " ";
                  }
                  file_peaks_pseudo_quad << std::endl;
               }
            }
            file_peaks_linear << std::endl;
         }

         file_peaks_linear.close();
         file_peaks_pseudo_quad.close();
         delete net_peaks;
      }

      // b% -------------------------------------------------------------------------- */
      // b%                                積分のチェック                                 */
      // b% -------------------------------------------------------------------------- */

      int end_divide = 50, divide = 5;
      std::vector<std::tuple<int, double>> results_quad(end_divide - divide + 1), results_linear(end_divide - divide + 1), result_original(end_divide - divide + 1);
      std::vector<std::tuple<int, double>> error_quad(end_divide - divide + 1), error_linear(end_divide - divide + 1);
#pragma omp parallel
      for (auto divide = 5; divide <= end_divide; ++divide)
#pragma omp nowait
      {
         /* ------------------------------------------------ */
         auto net_peaks = new Network();
         net_peaks->setFaces(get_triangle_verticies(divide));
         if (use_distorted_mesh)
            DistortMesh(net_peaks);
         for (const auto &f : net_peaks->getFaces()) {
            for (const auto &p : f->getPoints()) {
               auto [x, y, z] = p->X;
               p->setX(Tddd{x, y, peaksFunction(x, y)});
            }
         }
         net_peaks->setGeometricProperties();
         /* ------------------------------------------------ */

         double integrate_linear = 0., integrate_pseudo = 0., integrate_original = 0., error_integrate_linear = 0., error_integrate_pseudo = 0.;
         double tmp0, tmp1, tmp2, tmp3;
         Tddd N, X;

         for (const auto &f : net_peaks->getFaces()) {
            auto X012 = ToX(f->getPoints());
            // auto quadpoints = DodecaPoints(f, f->getPoints()[0], [&](const networkLine *l) { return l->getFaces().size() >= 2; });
            auto quadpoints = DodecaPoints(f, f->getPoints()[0], [&](const networkLine *l) { return l->getFaces().size() >= 2; });
            for (const auto &[t0, t1, ww] : __array_GW10xGW10__) {

               //@ 線形補間
               N = ModTriShape<3>(t0, t1);
               tmp0 = ww * Norm(Cross(X012[1] - X012[0], X012[2] - X012[0])) * (1 - t0);

               //@ 擬2次補間
               auto [xi0, xi1, xi2] = N;
               tmp1 = ww * Norm(quadpoints.cross(xi0, xi1)) * (1 - t0);

               //@ 厳密な外積を使って積分した場合
               {
                  auto [x, y, z] = Dot(N, X012);
                  auto DX0 = Dot(D_TriShape<3, 1, 0>(xi0, xi1), X012);
                  auto DX1 = Dot(D_TriShape<3, 0, 1>(xi0, xi1), X012);
                  DX0[2] = 0.;
                  DX1[2] = 0.;
                  tmp2 = ww * exact_norm_cross(x, y) * Norm(Cross(DX0, DX1)) * (1 - t0);
               }
               //
               {
                  auto [x, y, z] = quadpoints.X(xi0, xi1);
                  auto DX0 = quadpoints.D_X<1, 0>(xi0, xi1);
                  auto DX1 = quadpoints.D_X<0, 1>(xi0, xi1);
                  DX0[2] = 0.;
                  DX1[2] = 0.;
                  tmp3 = ww * exact_norm_cross(x, y) * Norm(Cross(DX0, DX1)) * (1 - t0);
               }
               integrate_linear += tmp0;
               integrate_pseudo += tmp1;
               integrate_original += tmp2;
               error_integrate_linear += std::abs(tmp0 - tmp2);
               error_integrate_pseudo += std::abs(tmp1 - tmp3);
            }
         }
         results_quad[divide - 5] = {divide, integrate_pseudo};
         results_linear[divide - 5] = {divide, integrate_linear};
         result_original[divide - 5] = {divide, integrate_original};
         error_quad[divide - 5] = {divide, error_integrate_pseudo};
         error_linear[divide - 5] = {divide, error_integrate_linear};
         delete net_peaks;
         std::cout << "id : " << id << "divide : " << divide << " is done." << std::endl;
      }

      std::ofstream file_peaks_integral_linear(_HOME_DIR_ + "/output/peaks_integral_linear" + id + ".dat");
      for (const auto &[divide, integrate_linear] : results_linear) {
         file_peaks_integral_linear << divide << " " << integrate_linear << std::endl;
         std::cout << "integrate_linear : " << divide << " : " << integrate_linear << std::endl;
      }
      file_peaks_integral_linear.close();
      std::ofstream file_peaks_integral_quad(_HOME_DIR_ + "/output/peaks_integral_pseudo_quad" + id + ".dat");
      for (const auto &[divide, integrate_pseudo] : results_quad) {
         file_peaks_integral_quad << divide << " " << integrate_pseudo << std::endl;
         std::cout << "integrate_pseudo : " << divide << " : " << integrate_pseudo << std::endl;
      }
      file_peaks_integral_quad.close();

      std::ofstream file_peaks_integral_original(_HOME_DIR_ + "/output/peaks_integral_original" + id + ".dat");
      for (const auto &[divide, integrate_original] : result_original) {
         file_peaks_integral_original << divide << " " << integrate_original << std::endl;
         std::cout << "integrate_original : " << divide << " : " << integrate_original << std::endl;
      }
      file_peaks_integral_original.close();

      std::ofstream file_peaks_integral_error_linear(_HOME_DIR_ + "/output/peaks_integral_linear_error" + id + ".dat");
      for (const auto &[divide, error_integrate_linear] : error_linear) {
         file_peaks_integral_error_linear << divide << " " << error_integrate_linear << std::endl;
         std::cout << "error_integrate_linear : " << divide << " : " << error_integrate_linear << std::endl;
      }
      file_peaks_integral_error_linear.close();

      std::ofstream file_peaks_integral_error_quad(_HOME_DIR_ + "/output/peaks_integral_pseudo_quad_error" + id + ".dat");
      for (const auto &[divide, error_integrate_pseudo] : error_quad) {
         file_peaks_integral_error_quad << divide << " " << error_integrate_pseudo << std::endl;
         std::cout << "error_integrate_pseudo : " << divide << " : " << error_integrate_pseudo << std::endl;
      }
      file_peaks_integral_error_quad.close();
   }
}
