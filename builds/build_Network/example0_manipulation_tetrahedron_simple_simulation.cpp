/*DOC_EXTRACT 0_1_manipulation_tetrahedron

## 四面体の操作

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_manipulation_tetrahedron.cpp
make
./example0_manipulation_tetrahedron
```

*/

#include "tetgen1.6.0/tetgen.h"
//
#include "Network.hpp"
#include "vtkWriter.hpp"

using Tddd = std::array<double, 3>;
using DataVariant = std::variant<double, Tddd>;
using DataMap = std::unordered_map<networkPoint*, DataVariant>;

int main() {

   /* -------------------------------------------------------------------- */
   /*                           四面体の生成と出力                            */
   /* -------------------------------------------------------------------- */

   auto coil = new Network("./input/coil.off");
   auto box = new Network("./input/box.off");

   coil->tetrahedralize();
   box->tetrahedralize();

   {
      std::ofstream ofs("./outptut/coil.vtu");
      vtkUnstructuredGridWrite(ofs, coil->getTetras());
      ofs.close();
   }
   {
      std::ofstream ofs("./outptut/box.vtu");
      vtkUnstructuredGridWrite(ofs, box->getTetras());
      ofs.close();
   }

   /* -------------------------------------------------------------------------- */

   coil->makeBucketPoints();
   coil->makeBucketFaces();
   coil->makeBucketTetras();
   box->makeBucketPoints();
   box->makeBucketFaces();
   box->makeBucketTetras();

   /* -------------------------------------------------------------------------- */

   std::vector<std::tuple<std::string, DataMap>> data = {{"fliped", data1}, {"xyz", data2}, {"tetra_size", data3}};

   // time stel
   const double dt = 0.01;
   const double simulation_time = 0;
   auto points = coil->getPoints();
   for (const auto& p : points)
      p->RK_X.initialize(dt, simulation_time, p->X, 4);

   std::ofstream ofs("./output/simu_box" + std::to_string(step) + ".vtu");
   vtkUnstructuredGridWrite(ofs, coil->getTetras(), data);
   ofs.close();

   std::array<double, 3> acceleration = {0., 0., 9.8};

   PVDWriter pvd("./output/simu_coil.pvd");

   for (auto step = 0; step < 100; ++step) {
      do {

         for (auto p : points) {
            p->RK_X.push(p->velocity);
            p->RK_Velocity.push(acceleration);
            p->velocity = p->RK_Velocity.get();
         }

         {
            std::strign name = "./output/simu_coil" + std::to_string(step) + ".vtu";
            std::ofstream ofs(name);
            vtkUnstructuredGridWrite(ofs, coil->getTetras(), data);
            ofs.close();
            pvd.push(name, step);
         }

      } while (std::any_of(points, [](const networkPoint* p) { return p->RK_X.finished; }));
   }
   pvd.output();
}