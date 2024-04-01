/*DOC_EXTRACT 1_3_1_interpolation

### 複雑な3Dオブジェクトの形状補間

この例では，`obj`ファイルを読み込んで，その面を補間する．

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=TriShapeExample_improved_test2.cpp
make
./TriShapeExample_improved_test2
```

<img src="TriShapeExample_improved_test2_torrus0d1remesh_linear0.png" width="700">

<img src="TriShapeExample_improved_test2_torrus0d1remesh_linear0_highreso.png" width="700">

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
   //
   Tddd X;
   /* -------------------------------- OBJ を使った例 ------------------------------- */
   for (std::string name : {"torrus0d1remesh", "torrus0d2remesh", "torrus0d3remesh"}) {
      // for (std::string name : {"plate0d2", "plate0d1", "plate0d075"}) {
      Network net(fs::path("./") / (name + ".obj"));
      std::vector<std::array<std::array<double, 2>, 3>> t0t1_triangles = SymmetricSubdivisionOfTriangle_00_10_01(9);
      int i = 0;
      for (auto f : net.getFaces()) {
         auto baseFileName = outputDir / name;
         std::ofstream file_linear(baseFileName.string() + "_linear" + std::to_string(i) + ".dat");
         std::ofstream file_quad(baseFileName.string() + "_quad" + std::to_string(i) + ".dat");
         std::ofstream file_pseudo_quadratic(baseFileName.string() + "_pseudo_quad" + std::to_string(i) + ".dat");
         std::ofstream file_pseudo_quadratic_cornertreated(baseFileName.string() + "_pseudo_quad_cornertreated" + std::to_string(i) + ".dat");
         for (auto t0t1 : t0t1_triangles) {
            for (auto [t0, t1] : t0t1) {
               auto X012 = f->vertices;
               // auto N3 = ModTriShape<3>(t0, t1);
               auto N3 = TriShape<3>(t0, t1);
               auto centered_N3 = Dot(N3, T3Tdd{{{0., 0.5}, {0.5, 0.}, {0.5, 0.5}}});
               {
                  auto [x, y, z] = X = Dot(N3, X012);
                  file_linear << x << " " << y << " " << z << " ";
               }
               {
                  auto [p0, p1, p2] = f->getPoints();
                  auto quadpoint = QuadPoints(f, p0, checker);
                  auto N6 = TriShape<6>(centered_N3, quadpoint.available);
                  auto [x, y, z] = X = Dot(N6, ToX(quadpoint.points));
                  file_quad << x << " " << y << " " << z << " ";
               }
               {
                  auto [p0, p1, p2] = f->getPoints();
                  DodecaPoints dodecapoints(f, p0);
                  auto [x, y, z] = dodecapoints.interpolate(t0, t1, [&](networkPoint *p) { return p->X; });
                  file_pseudo_quadratic << x << " " << y << " " << z << " ";
               }
               {
                  auto [p0, p1, p2] = f->getPoints();
                  DodecaPoints dodecapoints(f, p0, checker);
                  auto [x, y, z] = dodecapoints.interpolate(t0, t1, [&](networkPoint *p) { return p->X; });
                  file_pseudo_quadratic_cornertreated << x << " " << y << " " << z << " ";
               }
            }
            file_linear << std::endl;
            file_quad << std::endl;
            file_pseudo_quadratic << std::endl;
            file_pseudo_quadratic_cornertreated << std::endl;
         }
         i++;
         file_linear.close();
         file_quad.close();
         file_pseudo_quadratic.close();
         file_pseudo_quadratic_cornertreated.close();
      }
      /* -------------------------------------------------------------------------- */
      for (auto p : net.getPoints()) {
         auto [x, y, z] = p->X;
         p->value = std::cos(2 * 2 * M_PI * std::sqrt(x * x + y * y));
      }
      i = 0;
      for (auto f : net.getFaces()) {
         auto baseFileName = outputDir / name;
         std::ofstream file_linear(baseFileName.string() + "_linear_color" + std::to_string(i) + ".dat");
         std::ofstream file_quad(baseFileName.string() + "_quad_color" + std::to_string(i) + ".dat");
         std::ofstream file_pseudo_quadratic(baseFileName.string() + "_pseudo_quad_color" + std::to_string(i) + ".dat");
         for (auto t0t1 : t0t1_triangles) {
            for (auto [t0, t1] : t0t1) {
               auto X012 = f->vertices;
               // auto N3 = ModTriShape<3>(t0, t1);
               auto N3 = TriShape<3>(t0, t1);
               auto [x, y, z] = X = Dot(N3, X012);
               auto centered_N3 = Dot(N3, T3Tdd{{{0., 0.5}, {0.5, 0.}, {0.5, 0.5}}});
               {
                  file_linear << x << " " << y << " " << z << " " << Dot(N3, ToColor(f->getPoints())) << " ";
               }
               {
                  auto [p0, p1, p2] = f->getPoints();
                  auto quadpoint = QuadPoints(f, p0, checker);
                  auto N6 = TriShape<6>(centered_N3, quadpoint.available);
                  file_quad << x << " " << y << " " << z << " " << Dot(N6, ToColor(quadpoint.points)) << " ";
               }
               {
                  auto [p0, p1, p2] = f->getPoints();
                  DodecaPoints dodecapoints(f, p0, checker);
                  auto [x, y, z] = dodecapoints.interpolate(t0, t1, [&](networkPoint *p) { return p->X; });
                  file_pseudo_quadratic << x << " " << y << " " << z << " " << dodecapoints.interpolate(t0, t1, [&](networkPoint *p) { return p->value; }) << " ";
               }
            }
            file_linear << std::endl;
            file_quad << std::endl;
            file_pseudo_quadratic << std::endl;
         }
         i++;
         file_linear.close();
         file_quad.close();
         file_pseudo_quadratic.close();
      }
      std::cout << name << " is done. The file size is " << i << std::endl;
   }
}
