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

$$
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
$$

Fig. \ref{fig:pseudo_quad_schematic}に示すように，
この形状関数の係数を，対応する節点の値に掛けて足し合わせることで，
三角形要素の内部の任意の点における値を補間することができる．

<img src="pseudo_quad_white.png" width="700">

ただし，辺上の節点$3,4,5$は設定していないので，
隣接する三角形の頂点の値を使った2次補間の平均で近似する：

$$
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
$$

この式を\eqref{eq:general_usage_of_shape_function}の形に書き直すために．
次のような関係を使う：

$$
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
$$

これを使って，$V_{\rm {\ell}12in}$と$V_{\rm {\ell}21in}$の代わりに，$V_{\rm {\ell}01in}$を使った式に置き換える．

$$
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
$$

このように，\eqref{eq:general_usage_of_shape_function}の形の式を４つ足し合わせることで，擬2次補間を実装することができる．

ただし，辺が角を成している場合，この補間では，角にはならず，滑らかに補間されてしまう．
そのため，角を成している辺上の節点は，線形補間を使って近似することにする．つまり，辺が繋ぐ２節点の平均で近似する．
例えば，辺01が角となっている場合，\eqref{eq:pseudo_quadratic_interpolation}は次のように書き換える．

$$
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
$$

0,1節点は，$`V_{\rm {\ell}01in}`$における4,5節点であるため，$`2\xi_0\xi_1(0,0,0,1,1)^{\intercal}`$の項に

<img src="peak_function_interpolations.png">

積分結果を比較すると，線形補間よりも擬2次補間の方が精度が向上していることがわかる．
規則的なメッシュでは，２次補間と擬2次補間の結果は同程度であるが，
不規則なメッシュでは，擬2次補間の精度は２次補間により悪いが，線形ｆ補間よりも良い．

<img src="peak_function_interpolation_integration.png">

* [pseudo_quad要素の精度peaks_function積分.nb](pseudo_quad要素の精度peaks_function積分.nb)
* [説明.key](説明.key)


擬２次補間は，端部では線形補間を使うため，端部の精度は線形補間と同程度である．端部の精度を上げるには分割数を増やす他ない．
BEMに擬２次補間を適用した際に，造波精度が思ったよりも向上しなかったのは，このためだと考えられる．
ただし，喫水から離れた場所での，波の伝播に関する精度は，1.5倍以上の節点数での計算と同程度の精度が得られると考えられる．


*/

#include <array>
#include <cmath>
#include <concepts>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Network.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "kernelFunctions.hpp"
#include "basic_alias.hpp"          
#define GW4xGW4
#define bad_placement
#define peaks_function

bool checker(const networkLine *line) { return useOppositeFace(line, 0.7); };

double peaksFunction(double x, double y) {
#ifdef peaks_function
   const double x3 = x * x * x;
   const double y5 = y * y * y * y * y;
   const double term_A = 9.0 * std::exp(2.0 * x) * std::pow(x - 1.0, 2);
   const double term_B_inner = -x + 5.0 * x3 + 5.0 * y5;
   const double term_B = std::exp(2.0 * y) * (-1.0 + 6.0 * std::exp(1.0 + 2.0 * x) * term_B_inner);
   const double outer_factor = (1.0 / 3.0) * std::exp(-x * (2.0 + x) - std::pow(1.0 + y, 2));
   return outer_factor * (term_A + term_B);
#else
   return 0.;
#endif
}

double D_peaksFunction_dx(double x, double y) {
#ifdef peaks_function
   return -(1.0 / 3.0) * std::exp((-1 - x) * (1 + x) - y * y) * (-2 - 2 * x) - 6 * std::exp(-x * x - (1 + y) * (1 + y)) * (1 - x) - 6 * std::exp(-x * x - (1 + y) * (1 + y)) * (1 - x) * (1 - x) * x - 10 * std::exp(-x * x - y * y) * (1.0 / 5.0 - 3 * x * x) + 20 * std::exp(-x * x - y * y) * x * (x / 5.0 - x * x * x - std::pow(y, 5));
#else
   return 0.;
#endif
}

double D_peaksFunction_dy(double x, double y) {
#ifdef peaks_function
   return 2.0 / 3.0 * std::exp((-1 - x) * (1 + x) - y * y) * y + 50 * std::exp(-x * x - y * y) * std::pow(y, 4) - 6 * std::exp(-x * x - (1 + y) * (1 + y)) * (1 - x) * (1 - x) * (1 + y) + 20 * std::exp(-x * x - y * y) * y * (x / 5.0 - x * x * x - std::pow(y, 5));
#else
   return 0.;
#endif
}

T3Tdd domain = {{{2., -2.}, {-2., 2.}, {-2., -2.}}};

double exact_norm_cross(double x, double y) {
   return Norm(Cross(Tddd{1., 0., D_peaksFunction_dx(x, y)}, Tddd{0., 1., D_peaksFunction_dy(x, y)}));
}

// T3Tdd domain = {{{1., -1.}, {-1., 1.}, {-1., -1.}}};

int main() {

   auto get_triangle_verticies = [](const int divide) {
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
      std::cout << "get_triangle_verticies : " << triangle_vertices.size() << std::endl;
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
               l->flipIfTopologicallyBetter(10.0, 10.0, 5);

         for (const auto &l : lines)
            l->flipIfBetter(10.0, 10.0, 4);

         for (int j = 0; j < 5; ++j) {
// #pragma omp parallel
            for (const auto &p : points)
// #pragma omp nowait
            {
               if (std::ranges::none_of(p->getLines(), [](const auto &l) { return l->getFaces().size() == 1; })) {
                  // if (j < 5)
                  //    p->setX(p->X + 0.5 * NeighborAverageSmoothingVector(p, p->X));
                  // else
                  auto V = DistorsionMeasureWeightedSmoothingVector(p, p->X);
                  if (isFinite(V))
                     p->setX(p->X + 0.05 * V);
               }
            }
            net->setGeometricProperties();
         }
      }
   };

   Tddd X;

#if defined(GW4xGW4)
   auto GWGW_linear = __array_GW4xGW4__;
   auto GWGW_quad = __array_GW8xGW8__;

   std::string foldername = "GW4GW4";
#else
   auto GWGW_linear = __array_GW7xGW7__;
   auto GWGW_quad = __array_GW14xGW14__;

   std::string foldername = "GW7GW7";
#endif

#if defined(peaks_function)
   foldername += "_peak";
#else
   foldername += "_flat";
#endif
#if defined(bad_placement)
   foldername += "_bad";
#endif

   for (const auto use_distorted_mesh : {false,true}) {

      std::string id = use_distorted_mesh ? "_distorted" : "_regular";

      //@ -------------------------------------------------------------------------- */
      //@                             peaksFunctionを使った例                          */
      //@ -------------------------------------------------------------------------- */
      if (false)
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
                  int N;
                  if (std::abs(x) < 2 && std::abs(y) < 2)
                     N = 10;
                  else
                     N = 10;
                  for (auto t0t1 : SymmetricSubdivisionOfTriangle_00_10_01(N)) {
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

      // continue;

      // b% -------------------------------------------------------------------------- */
      // b%                                積分のチェック                                  */
      // b% -------------------------------------------------------------------------- */

      std::vector<std::tuple<int, int, int, double> > results_pseudo_quad, error_pseudo_quad;
      std::vector<std::tuple<int, int, int, double> > results_pseudo_quad_NoEdge, error_pseudo_quad_NoEdge;
      std::vector<std::tuple<int, int, int, double> > result_quad, error_quad;
      std::vector<std::tuple<int, int, int, double> > result_quad_NoEdge, error_quad_NoEdge;  // 追加
      std::vector<std::tuple<int, int, int, double> > results_linear, error_linear;
      std::vector<std::tuple<int, int, int, double> > results_linear_NoEdge, error_linear_NoEdge;  // 追加
      std::vector<std::tuple<int, int, int, double> > result_original;

      std::vector<std::tuple<int, int, int, double> > results_pseudo_quad_G, error_pseudo_quad_G;
      std::vector<std::tuple<int, int, int, double> > results_pseudo_quad_NoEdge_G, error_pseudo_quad_NoEdge_G;
      std::vector<std::tuple<int, int, int, double> > result_quad_G, error_quad_G;
      std::vector<std::tuple<int, int, int, double> > result_quad_NoEdge_G, error_quad_NoEdge_G;  // 追加
      std::vector<std::tuple<int, int, int, double> > results_linear_G, error_linear_G;
      std::vector<std::tuple<int, int, int, double> > results_linear_NoEdge_G, error_linear_NoEdge_G;  // 追加
      std::vector<std::tuple<int, int, int, double> > result_original_G;

      auto run = [&](const int divide) {
         std::cout << Green << "Generating Network. divide : " << divide << colorReset << std::endl;
         /* ------------------------------------------------------ */
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
         std::cout << Red << "Network is generated. node size : " << net_peaks->getPoints().size() << colorReset << std::endl;
         net_peaks->setGeometricProperties();

         /* ------------------------------------------------ */
         for (auto fundamental_solution_Q : {true, false})  // fundamental_solution_Q = true の場合は、擬2次補間の精度が向上する
         {
            const Tddd A = {2., 2., peaksFunction(2., 2.)};
            auto G = [&](const Tddd &X) {
               if (fundamental_solution_Q) {
                  return 1. / Norm(X - A);
               } else
                  return 1.;
            };
            auto G_exact = [&](double x, double y) {
               if (fundamental_solution_Q) {
                  Tddd X = {x, y, peaksFunction(x, y)};
                  return 1. / Norm(X - A);
               } else
                  return 1.;
            };
            // ...existing code...
            double I_original = 0.;
            double I_linear = 0., error_I_linear = 0.;
            double I_linear_NoEdge = 0., error_I_linear_NoEdge = 0.;
            double I_quad = 0., error_I_quad = 0.;
            double I_quad_NoEdge = 0., error_I_quad_NoEdge = 0.;
            double I_pseudo = 0., error_I_pseudo = 0.;
            double I_pseudo_NoEdge = 0., error_I_pseudo_NoEdge = 0.;
            // ...existing code...
            Tddd N_xi, X_xi, DX1, DX0;
            double dXdXi, dXidZeta, tmp, exact, error, dxydxi0xi1;
            //
            auto weight_exact_norm_cross_J = [&](double quadrature_weight, Tddd xyz, Tddd DX0, Tddd DX1) {
               DX0[2] = DX1[2] = 0.;
               auto x = xyz[0];
               auto y = xyz[1];
               return quadrature_weight * G_exact(x, y) * exact_norm_cross(x, y) /*constant 1.*/ * Norm(Cross(DX0, DX1));
            };
            int gauss_points_linear = 0;
            int gauss_points_linear_NoEdge = 0;
            int gauss_points_standard_quad = 0;
            int gauss_points_standard_quad_NoEdge = 0;
            T3Tddd X012;
            networkPoint *pOrg;
            for (const auto &f : net_peaks->getFaces()) {
               auto [p0, p1, p2] = f->getPoints();
               int index = 0;
               pOrg = p0;

               auto a = A;
#ifdef bad_placement
               a = Tddd{-2., -2., peaksFunction(-2., -2.)};
#endif
               if (Norm(p1->X - a) < Norm(p0->X - a) && Norm(p1->X - a) < Norm(p2->X - a)) {
                  index = 1;
                  pOrg = p1;
               } else if (Norm(p2->X - a) < Norm(p0->X - a) && Norm(p2->X - a) < Norm(p1->X - a)) {
                  index = 2;
                  pOrg = p2;
               }
               X012 = ToX(f->getPoints(pOrg));
               //{2,2,peaksFunction(2.,2.)}に最も近いものを[0]とする
               /* --------------------------------------------------- */
               std::array<Tddd, 6> X012345;
               // P0, P1, P2 は頂点そのもの
               X012345[0] = X012[0];
               X012345[1] = X012[1];
               X012345[2] = X012[2];
               // P3, P4, P5 は辺の中点
               auto P0 = X012[0];
               auto P1 = X012[1];
               auto P2 = X012[2];
               auto P3 = 0.5 * (X012[0] + X012[1]);
               auto P4 = 0.5 * (X012[1] + X012[2]);
               auto P5 = 0.5 * (X012[2] + X012[0]);
               X012345[3] = {P3[0], P3[1], peaksFunction(P3[0], P3[1])};
               X012345[4] = {P4[0], P4[1], peaksFunction(P4[0], P4[1])};
               X012345[5] = {P5[0], P5[1], peaksFunction(P5[0], P5[1])};
               /* ---------------------------------------------------- */
               auto dodecapoints = DodecaPoints(f, f->getPoints(pOrg)[0], [&](const networkLine *l) { return l->getFaces().size() >= 2; });
               /* ----------------------- 数値積分 --------------------- */
               bool is_edge = isEdge(f);
               for (const auto &[zeta0, zeta1, ww] : GWGW_linear) {
                  N_xi = ModTriShape<3>(zeta0, zeta1);
                  auto [xi0, xi1, xi2] = N_xi;
                  dXidZeta = (1 - zeta0);
                  /* ----------------------------------------------------- */
                  {
                     //@ 線形補間
                     DX0 = Dot(D_TriShape<3, 1, 0>(xi0, xi1), X012);
                     DX1 = Dot(D_TriShape<3, 0, 1>(xi0, xi1), X012);
                     dXdXi = Norm(Cross(DX0, DX1));  //! Jacobian (xi0,xi1)->(x,y)
                     X_xi = Dot(N_xi, X012);
                     tmp = ww * G(X_xi) * dXdXi;
                     I_linear += tmp* dXidZeta;
                     //
                     //@ 厳密な外積を使って積分した場合 for 線形補間
                     DX0[2] = DX1[2] = 0.;
                     dxydxi0xi1 = Norm(Cross(DX0, DX1));
                     exact = ww * G_exact(X_xi[0], X_xi[1]) * exact_norm_cross(X_xi[0], X_xi[1]) * dxydxi0xi1;
                     //
                     error = std::pow(tmp - exact, 2) * dXidZeta;
                     error_I_linear += error;
                     if (!is_edge) {
                        I_linear_NoEdge += tmp;
                        error_I_linear_NoEdge += error;
                     }
                     //
                     I_original += exact;
                  }
                  /* ----------------------------------------------------- */
                  {
                     //@ 擬2次補間
                     dXdXi = Norm(dodecapoints.cross(xi0, xi1));  //! Jacobian (xi0,xi1)->(x,y)
                     X_xi = dodecapoints.X(xi0, xi1);
                     tmp = ww * G(X_xi) * dXdXi;
                     I_pseudo += tmp * dXidZeta;
                     //
                     DX0 = dodecapoints.D_X<1, 0>(xi0, xi1);
                     DX1 = dodecapoints.D_X<0, 1>(xi0, xi1);
                     DX0[2] = DX1[2] = 0.;
                     dxydxi0xi1 = Norm(Cross(DX0, DX1));
                     exact = ww * G_exact(X_xi[0], X_xi[1]) * exact_norm_cross(X_xi[0], X_xi[1]) * dxydxi0xi1 ;
                     error = std::pow(tmp - exact, 2)* dXidZeta;
                     error_I_pseudo += error;
                     if (!is_edge) {
                        I_pseudo_NoEdge += tmp;
                        error_I_pseudo_NoEdge += error;
                        gauss_points_linear_NoEdge++;
                     }
                     gauss_points_linear++;
                  }
                  /* ----------------------------------------------------- */
               }
               /* -------------------------------------------------------------------------- */
               if (divide <= 10000)
                  for (const auto &[zeta0, zeta1, ww] : GWGW_quad) {
                     auto [xi0, xi1, xi2] = ModTriShape<3>(zeta0, zeta1);
                     /* ----------------------------------------------------- */
                     {
                        //@ 標準2次補間
                        //! １．辺の中点のx,yを計算
                        //! ２．そのx3=(x01,y01),x4=(x12,y12),x5=(x20,y20)におけz3,z4,z5=peak()を計算し X012345を作る．
                        dXdXi = Norm(dodecapoints.cross(xi0, xi1, X012345));  //! Jacobian (xi0,xi1)->(x,y)
                        dXidZeta = (1 - zeta0);
                        X_xi = dodecapoints.X(xi0, xi1, X012345);
                        tmp = ww * G(X_xi) * dXdXi ;
                        I_quad += tmp* dXidZeta;
                        //
                        DX0 = dodecapoints.D_X<1, 0>(xi0, xi1, X012345);
                        DX1 = dodecapoints.D_X<0, 1>(xi0, xi1, X012345);
                        DX0[2] = DX1[2] = 0.;
                        dxydxi0xi1 = Norm(Cross(DX0, DX1));
                        exact = ww * G_exact(X_xi[0], X_xi[1]) * exact_norm_cross(X_xi[0], X_xi[1]) * dxydxi0xi1 ;
                        error = std::pow(tmp - exact, 2)* dXidZeta;
                        error_I_quad += error;
                        if (!is_edge) {
                           I_quad_NoEdge += tmp;
                           error_I_quad_NoEdge += error;
                           gauss_points_standard_quad_NoEdge++;
                        }
                        gauss_points_standard_quad++;
                     }
                     /* ----------------------------------------------------- */
                  }
            }

            //! quad辺の上の節点の数が重複するので面を２つ持つ辺の数を数えてその分をsから引く

            int Np = net_peaks->getPoints().size();
            int Np_NoEdge = net_peaks->getPoints([](networkPoint *p) { return !isEdge(p); }).size();
            int Np4Quad = Np + net_peaks->getLines().size();
            //! さらに端部と関係するedgeも消す
            int Np4Quad_NoEdge = Np_NoEdge + net_peaks->getLines([](networkLine *l) { auto p01 = l->getPoints(); return !isEdge(p01[0]) & !isEdge(p01[1]); }).size();

            std::cout << "I_original : " << I_original << std::endl;
            std::cout << "Np : " << Np << " Np_NoEdge : " << Np_NoEdge << " Np4Quad : " << Np4Quad << " Np4Quad_NoEdge : " << Np4Quad_NoEdge << std::endl;
            // openmpによる重複に注意

            error_I_linear = std::sqrt(error_I_linear);
            error_I_linear_NoEdge = std::sqrt(error_I_linear_NoEdge);
            error_I_pseudo = std::sqrt(error_I_pseudo);
            error_I_pseudo_NoEdge = std::sqrt(error_I_pseudo_NoEdge);
            error_I_quad = std::sqrt(error_I_quad);
            error_I_quad_NoEdge = std::sqrt(error_I_quad_NoEdge);

#pragma omp critical
            {
               if (fundamental_solution_Q) {
                  results_linear_G.push_back({divide, Np, gauss_points_linear, I_linear});
                  results_linear_NoEdge_G.push_back({divide, Np_NoEdge, gauss_points_linear_NoEdge, I_linear_NoEdge});
                  results_pseudo_quad_G.push_back({divide, Np, gauss_points_linear, I_pseudo});
                  results_pseudo_quad_NoEdge_G.push_back({divide, Np_NoEdge, gauss_points_linear_NoEdge, I_pseudo_NoEdge});
                  error_linear_G.push_back({divide, Np, gauss_points_linear, error_I_linear});
                  error_linear_NoEdge_G.push_back({divide, Np_NoEdge, gauss_points_linear_NoEdge, error_I_linear_NoEdge});
                  error_pseudo_quad_G.push_back({divide, Np, gauss_points_linear, error_I_pseudo});
                  error_pseudo_quad_NoEdge_G.push_back({divide, Np_NoEdge, gauss_points_linear_NoEdge, error_I_pseudo_NoEdge});

                  if (divide <= 10000) {
                     result_quad_G.push_back({divide, Np4Quad, gauss_points_standard_quad, I_quad});
                     result_quad_NoEdge_G.push_back({divide, Np4Quad_NoEdge, gauss_points_standard_quad_NoEdge, I_quad_NoEdge});
                     error_quad_G.push_back({divide, Np4Quad, gauss_points_standard_quad, error_I_quad});
                     error_quad_NoEdge_G.push_back({divide, Np4Quad_NoEdge, gauss_points_standard_quad_NoEdge, error_I_quad_NoEdge});
                  }
               } else {
                  results_linear.push_back({divide, Np, gauss_points_linear, I_linear});
                  results_linear_NoEdge.push_back({divide, Np_NoEdge, gauss_points_linear_NoEdge, I_linear_NoEdge});
                  results_pseudo_quad.push_back({divide, Np, gauss_points_linear, I_pseudo});
                  results_pseudo_quad_NoEdge.push_back({divide, Np_NoEdge, gauss_points_linear_NoEdge, I_pseudo_NoEdge});
                  error_linear.push_back({divide, Np, gauss_points_linear, error_I_linear});
                  error_linear_NoEdge.push_back({divide, Np_NoEdge, gauss_points_linear_NoEdge, error_I_linear_NoEdge});
                  error_pseudo_quad.push_back({divide, Np, gauss_points_linear, error_I_pseudo});
                  error_pseudo_quad_NoEdge.push_back({divide, Np_NoEdge, gauss_points_linear_NoEdge, error_I_pseudo_NoEdge});

                  if (divide <= 10000) {
                     result_quad.push_back({divide, Np4Quad, gauss_points_standard_quad, I_quad});
                     result_quad_NoEdge.push_back({divide, Np4Quad_NoEdge, gauss_points_standard_quad_NoEdge, I_quad_NoEdge});
                     error_quad.push_back({divide, Np4Quad, gauss_points_standard_quad, error_I_quad});
                     error_quad_NoEdge.push_back({divide, Np4Quad_NoEdge, gauss_points_standard_quad_NoEdge, error_I_quad_NoEdge});
                  }
               }
            }
         }
         // result_original.push_back({divide, Np, I_original});
         delete net_peaks;
      };

#pragma omp parallel
      for (auto i : {1., 1.25, 1.75, 2.75, 4.75, 7.75})
#pragma omp single nowait
         for (auto order : {40000., 10000., 1000., 100., 10.}) {
            run((int)(std::sqrt(i * order)));
         }

      /* ----------------------------------------------------------- */
      {
         auto write_and_print_results = [&](const std::string &file_basename, const auto &results_vector) {
            std::ofstream file_stream("./output_fundamental_" + foldername + "/" + file_basename + id + ".dat");
            std::cout << "--- Writing " << file_basename << id << ".dat ---" << std::endl;
            for (const auto &[divide, s, gausspoints, value] : results_vector) {
               file_stream << divide << " " << s << " " << gausspoints << " " << std::setprecision(16) << value << std::endl;
               std::cout << file_basename << " : " << divide << " : " << value << std::endl;
            }
         };

         write_and_print_results("integral_linear", results_linear_G);
         write_and_print_results("integral_linear_NoEdge", results_linear_NoEdge_G);
         write_and_print_results("integral_pseudo_quad", results_pseudo_quad_G);
         write_and_print_results("integral_pseudo_quad_NoEdge", results_pseudo_quad_NoEdge_G);
         write_and_print_results("integral_quad", result_quad_G);  // 2次補間を追加した場合
         write_and_print_results("integral_quad_NoEdge", result_quad_NoEdge_G);
         // write_and_print_results("integral_original_G", "I_original_G", result_original_G);
         write_and_print_results("integral_linear_error", error_linear_G);
         write_and_print_results("integral_linear_NoEdge_error", error_linear_NoEdge_G);
         write_and_print_results("integral_pseudo_quad_error", error_pseudo_quad_G);
         write_and_print_results("integral_pseudo_quad_NoEdge_error", error_pseudo_quad_NoEdge_G);
         write_and_print_results("integral_quad_error", error_quad_G);  // 2次補間を追加した場合
         write_and_print_results("integral_quad_NoEdge_error", error_quad_NoEdge_G);
      }
      /* ----------------------------------------------------------- */
      {
         auto write_and_print_results = [&](const std::string &file_basename, const auto &results_vector) {
            std::ofstream file_stream("./output_" + foldername + "/" + file_basename + id + ".dat");
            std::cout << "--- Writing " << file_basename << id << ".dat ---" << std::endl;
            for (const auto &[divide, s, gausspoints, value] : results_vector) {
               file_stream << divide << " " << s << " " << gausspoints << " " << std::setprecision(16) << value << std::endl;
               std::cout << file_basename << " : " << divide << " : " << value << std::endl;
            }
         };

         write_and_print_results("integral_linear", results_linear);
         write_and_print_results("integral_linear_NoEdge", results_linear_NoEdge);
         write_and_print_results("integral_pseudo_quad", results_pseudo_quad);
         write_and_print_results("integral_pseudo_quad_NoEdge", results_pseudo_quad_NoEdge);
         write_and_print_results("integral_quad", result_quad);  // 2次補間を追加した場合
         write_and_print_results("integral_quad_NoEdge", result_quad_NoEdge);
         // write_and_print_results("integral_original", "I_original", result_original);

         write_and_print_results("integral_linear_error", error_linear);
         write_and_print_results("integral_linear_NoEdge_error", error_linear_NoEdge);
         write_and_print_results("integral_pseudo_quad_error", error_pseudo_quad);
         write_and_print_results("integral_pseudo_quad_NoEdge_error", error_pseudo_quad_NoEdge);
         write_and_print_results("integral_quad_error", error_quad);  // 2次補間を追加した場合
         write_and_print_results("integral_quad_NoEdge_error", error_quad_NoEdge);
      }
      /* ----------------------------------------------------------- */
   }
}
