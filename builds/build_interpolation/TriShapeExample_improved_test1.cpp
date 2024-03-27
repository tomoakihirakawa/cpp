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
         std::ofstream file_peaks_linear(_HOME_DIR_ + "/output/peaksFunction_linear" + id + std::to_string(divide) + ".dat");
         std::ofstream file_peaks_pseudo_quad(_HOME_DIR_ + "/output/peaksFunction_pseudo_quad" + id + std::to_string(divide) + ".dat");
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
               for (auto t0t1 : SubdivideTriangleIntoTriangles(6)) {
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

      int end_divide = 100, divide = 5;
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
