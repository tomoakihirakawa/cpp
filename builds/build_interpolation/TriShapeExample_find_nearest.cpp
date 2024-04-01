/*DOC_EXTRACT 1_3_1_interpolation

### 複雑な3Dオブジェクトの形状補間

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=TriShapeExample_find_nearest.cpp
make
./TriShapeExample_find_nearest
```


*/

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Network.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "kernelFunctions.hpp"
#include "minMaxOfFunctions.hpp"

namespace fs = std::filesystem;

template <std::size_t N>
std::array<double, N> ToColor(const std::array<networkPoint *, N> &points) {
   std::array<double, N> colors;
   for (int i = 0; auto p : points)
      colors[i++] = p->value;
   return colors;
}

bool checker(const networkLine *line) {
   // return false;
   auto faces = line->getFaces();
   if (faces.size() != 2)
      return false;
   return isFlat(faces[0]->normal, faces[1]->normal, M_PI / 4.);
};

constexpr std::array<std::array<double, 2>, 3> M_sub{{{0., 0.5}, {0.5, 0.}, {0.5, 0.5}}};

int main() {
   fs::path outputDir = fs::path(_HOME_DIR_) / "output";
   fs::create_directories(outputDir);
   Tddd X;

   /* -------------------------------- OBJ を使った例 ------------------------------- */

   std::string name = "torrus0d3remesh";

   // for (std::string name : {"plate0d2", "plate0d1", "plate0d075"}) {
   Network net(fs::path("./") / (name + ".obj"));
   std::vector<std::array<std::array<double, 2>, 3>> t0t1_triangles = SymmetricSubdivisionOfTriangle_00_10_01(9);
   int i = 0;

   std::vector<Tddd> nearestXs, nearestXs_using_findMinium;
   std::vector<double> min0, min1;

   double max = 100.;
   for (auto k = 0; k < max; k++) {
      double t = k / max;
      const double r = 0.5;
      Tddd a = {r * std::cos(2. * M_PI * t), r * std::sin(2. * M_PI * t), 0.3};

      Tddd nearestX = {1E+10, 1E+10, 1E+10};
      Tddd nearestX_using_findMinium = {1E+10, 1E+10, 1E+10};

      for (auto f : net.getFaces()) {
         auto baseFileName = outputDir / name;
         std::ofstream file_linear(baseFileName.string() + "_linear" + std::to_string(i) + ".dat");
         std::ofstream file_quad(baseFileName.string() + "_quad" + std::to_string(i) + ".dat");
         std::ofstream file_pseudo_quadratic(baseFileName.string() + "_pseudo_quad" + std::to_string(i) + ".dat");
         std::ofstream file_pseudo_quadratic_cornertreated(baseFileName.string() + "_pseudo_quad_cornertreated" + std::to_string(i) + ".dat");

         DodecaPoints dodecapoints(f, f->getPoints()[0], checker);

         //   NewtonRaphson<Tdd> nr(1. / 3., 1. / 3.);
         //   for (auto i = 0; i < 100; i++) {
         //      double[t0, t1] = nr.X;
         //      Tdd gradF = {dodecapoints.D_interpolate<1, 0>(t0, t1, [&](networkPoint *p) { return p->X - a; }),
         //                   dodecapoints.D_interpolate<0, 1>(t0, t1, [&](networkPoint *p) { return p->X - a; })};
         //      T2Tdd Hessian = {{dodecapoints.D_interpolate<2, 0>(t0, t1, [&](networkPoint *p) { return p->X - a; }),
         //                        dodecapoints.D_interpolate<1, 1>(t0, t1, [&](networkPoint *p) { return p->X - a; })},
         //                       {dodecapoints.D_interpolate<1, 1>(t0, t1, [&](networkPoint *p) { return p->X - a; }),
         //                        dodecapoints.D_interpolate<0, 2>(t0, t1, [&](networkPoint *p) { return p->X - a; })}};
         //      nr.update(gradF, Hessian);
         //      if (Norm(nr.dX) < 1e-10)
         //         break;
         //   }

         /* -------------------------------------------------------------------------- */

         std::function<Tddd(networkPoint *)> conversion([&](networkPoint *p) -> Tddd { return p->X; });
         auto [t0t1, min_value_on_this_surface] = dodecapoints.findMinimum(conversion, [&](const Tddd &X) -> double { return Norm(X - a); });

         if (Norm(nearestX_using_findMinium - a) > Norm(min_value_on_this_surface - a))
            nearestX_using_findMinium = dodecapoints.interpolate(t0t1, [&](networkPoint *p) { return p->X; });

         auto [t0t1_, nearest] = dodecapoints.Nearest(a);

         if (Norm(nearestX - a) > Norm(nearest - a))
            nearestX = nearest;

         /* -------------------------------------------------------------------------- */
         //  for (auto t0t1 : t0t1_triangles) {
         //     for (auto [t0, t1] : t0t1) {
         //        auto X012 = f->vertices;
         //        // auto N3 = ModTriShape<3>(t0, t1);
         //        auto N3 = TriShape<3>(t0, t1);
         //        auto centered_N3 = Dot(N3, T3Tdd{{{0., 0.5}, {0.5, 0.}, {0.5, 0.5}}});
         //        {
         //           auto [x, y, z] = X = Dot(N3, X012);
         //           file_linear << x << " " << y << " " << z << " ";
         //        }
         //        {
         //           auto [p0, p1, p2] = f->getPoints();
         //           auto quadpoint = QuadPoints(f, p0, checker);
         //           auto N6 = TriShape<6>(centered_N3, quadpoint.available);
         //           auto [x, y, z] = X = Dot(N6, ToX(quadpoint.points));
         //           file_quad << x << " " << y << " " << z << " ";
         //        }
         //        {
         //           auto [p0, p1, p2] = f->getPoints();
         //           DodecaPoints dodecapoints(f, p0);
         //           auto [x, y, z] = dodecapoints.interpolate(t0, t1, [&](networkPoint *p) { return p->X; });
         //           file_pseudo_quadratic << x << " " << y << " " << z << " ";
         //        }
         //        {
         //           auto [p0, p1, p2] = f->getPoints();
         //           DodecaPoints dodecapoints(f, p0, checker);
         //           auto [x, y, z] = dodecapoints.interpolate(t0, t1, [&](networkPoint *p) { return p->X; });
         //           file_pseudo_quadratic_cornertreated << x << " " << y << " " << z << " ";
         //        }
         //     }
         //     file_linear << std::endl;
         //     file_quad << std::endl;
         //     file_pseudo_quadratic << std::endl;
         //     file_pseudo_quadratic_cornertreated << std::endl;
         //  }
         //  i++;
         //  file_linear.close();
         //  file_quad.close();
         //  file_pseudo_quadratic.close();
         //  file_pseudo_quadratic_cornertreated.close();
      }

      nearestXs.push_back(nearestX);
      nearestXs_using_findMinium.push_back(nearestX_using_findMinium);

      min0.push_back(Norm(nearestX - a));
      min1.push_back(Norm(nearestX_using_findMinium - a));
   }
   std::ofstream file_nearest(outputDir / (name + "_nearest.dat"));
   for (auto [x, y, z] : nearestXs)
      file_nearest << x << " " << y << " " << z << std::endl;
   file_nearest.close();

   std::ofstream file_nearest_using_findMinium(outputDir / (name + "_nearest_using_findMinium.dat"));
   for (auto [x, y, z] : nearestXs_using_findMinium)
      file_nearest_using_findMinium << x << " " << y << " " << z << std::endl;
   file_nearest_using_findMinium.close();

   std::cout << "total min0 (Nearest): " << Sum(min0) << std::endl;
   std::cout << "total min1 (findMinium): " << Sum(min1) << std::endl;
   std::cout << name << " is done. The file size is " << i << std::endl;
}
