/*DOC_EXTRACT 0_1_quadratic_interpolation

## ２次補間

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_quadratic_interpolation.cpp
make
./example0_quadratic_interpolation
```

*/

#include "Network.hpp"
#include "vtkWriter.hpp"

int main() {

   PVDWriter pvd("./semisub.pvd");
   PVDWriter pvd_interp_points_middle("./interp_points_middle.pvd");
   PVDWriter pvd_interp_points_full("./interp_points_full.pvd");

   auto obj = new Network("./semisub.obj");
   auto interp_points_middle = new Network();
   auto interp_points_full = new Network();
   auto points_for_line = new Network();
   constexpr std::array<std::array<double, 2>, 3> parameters{{{0., 0.5}, {0.5, 0.}, {0.5, 0.5}}};

   // set CORNER
   for (const auto& l : obj->getLines()) {
      auto fs = l->getFaces();
      l->CORNER = !isFlat(fs[0]->normal, fs[1]->normal, M_PI / 180. * 30.);
   }
   int i = 0;
   for (const auto& f : obj->getFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      QuadPoints quadpoint(p0, f);

      {
         auto getPoint = [&](const auto& l, const auto& f) {
            Tddd X = {0., 0., 0.};
            if (l->CORNER)
               X = l->X;
            else {
               auto fs = l->getFaces();
               {
                  QuadPoints quadpoint(l, fs[0]);
                  X += Dot(TriShape<6>(0.25, 0.25, quadpoint.ignore), ToX(quadpoint.points));
               }
               {
                  QuadPoints quadpoint(l, fs[1]);
                  X += Dot(TriShape<6>(0.25, 0.25, quadpoint.ignore), ToX(quadpoint.points));
               }
               X /= 2.;
            }
            new networkPoint(points_for_line, X);
         };
         for (const auto& l : f->getLines())
            getPoint(l, f);
      }

      if (i < 10) {
         for (const auto& [t0, t1, ww] : __array_GW6xGW6__) {
            auto N = TriShape<6>(Dot(ModTriShape<3>(t0, t1), parameters), quadpoint.ignore);
            new networkPoint(interp_points_middle, Dot(N, ToX(quadpoint.points)));
            auto M = TriShape<6>(t0, t1, quadpoint.ignore);
            new networkPoint(interp_points_full, Dot(M, ToX(quadpoint.points)));
         }
         {
            auto filename = "./points_middle" + std::to_string(i) + ".vtp";
            std::ofstream ofs(filename);
            vtkPolygonWrite(ofs, interp_points_middle->getPoints());
            std::cout << filename << std::endl;
            pvd_interp_points_middle.push(filename, i);
         }
         {
            auto filename = "./points_full" + std::to_string(i) + ".vtp";
            std::ofstream ofs(filename);
            vtkPolygonWrite(ofs, interp_points_full->getPoints());
            std::cout << filename << std::endl;
            pvd_interp_points_full.push(filename, i);
         }
      }
      i++;
   }
   auto filename = "./points_for_line.vtp";
   std::ofstream ofs(filename);
   vtkPolygonWrite(ofs, points_for_line->getPoints());
   std::cout << filename << std::endl;
   pvd_interp_points_middle.output();
   pvd_interp_points_full.output();

   //    for (auto i = 0; i < 18; ++i) {
   //       obj->rotate(20. / 180. * M_PI, {0., 0., 1.});
   //       auto filename = "./semisub" + std::to_string(i) + ".vtp";
   //       std::ofstream ofs(filename);
   //       vtkPolygonWrite(ofs, obj->getFaces());
   //       pvd.push(filename, i);
   //    }
   //    pvd.output();
}