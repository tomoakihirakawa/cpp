/*DOC_EXTRACT 0_1_manipulation_tetrahedron

## 四面体の操作

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_manipulation_tetrahedron_simple_simulation.cpp
make
./example0_manipulation_tetrahedron_simple_simulation
```

*/

#include "tetgen1.6.0/tetgen.h"
//
#include "Network.hpp"
#include "vtkWriter.hpp"

using Tddd = std::array<double, 3>;
using DataVariant = std::variant<double, Tddd>;
using DataMap = std::unordered_map<networkPoint*, DataVariant>;

int main(const int argc, const char** argv) {

   std::string name0 = argc > 1 ? argv[1] : "./input/coil.obj";
   std::string name1 = argc > 2 ? argv[2] : "./input/box.obj";

   /* -------------------------------------------------------------------- */
   /*                           四面体の生成と出力                            */
   /* -------------------------------------------------------------------- */

   auto coil = new Network(name0);
   auto box = new Network(name1);

   if (argc > 3) {
      coil->tetrahedralize(argv[3]);
      box->tetrahedralize(argv[3]);
   } else {
      coil->tetrahedralize();
      box->tetrahedralize();
   }

   {
      std::ofstream ofs("./output/coil.vtu");
      vtkUnstructuredGridWrite(ofs, coil->getTetras());
      ofs.close();
   }

   {
      std::ofstream ofs("./output/box.vtu");
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

   std::cout << "tetras.size() : " << coil->getTetras().size() << std::endl;
   std::cout << "tetras.size() : " << box->getTetras().size() << std::endl;
   /* -------------------------------------------------------------------------- */

   //% 1. 衝突を検知したいオブジェクトのバケツを取得
   //% 2. 衝突を探査する距離を設定．表面からその距離までの範囲で衝突を探査する
   //% 3. 互いの表面節点をトラバースし，衝突を検知する

   // std::vector<std::tuple<std::string, DataMap>> data = {{"fliped", data1}, {"xyz", data2}, {"tetra_size", data3}};

   // time stel
   const double dt = 0.01;
   double simulation_time = 0;

   std::ofstream ofs("./output/simu_box.vtu");
   vtkUnstructuredGridWrite(ofs, coil->getTetras());
   ofs.close();

   std::array<double, 6> acceleration = {0., 0., 9.8, 0., 0., 0.};

   PVDWriter pvd("./output/simu_coil.pvd");
   PVDWriter pvd_collision("./output/collision.pvd");

   for (auto step = 0; step < 100; ++step) {
      for (const auto& p : coil->getPoints()) {
         p->RK_X.initialize(dt, simulation_time, p->X, 4);
         p->RK_generalized_velocity.initialize(dt, simulation_time, p->velocity, 4);
      }
      std::cout << "step : " << step << std::endl;
      do {
         for (auto p : coil->getPoints()) {
            p->RK_X.push(p->velocityTranslational());
            p->RK_generalized_velocity.push(acceleration);

            auto tmp = p->RK_generalized_velocity.getX();
            p->velocity[0] = tmp[0];
            p->velocity[1] = tmp[1];
            p->velocity[2] = tmp[2];
            p->setXSingle(p->RK_X.getX());
         }

         {
            std::string name = "./output/simu_coil" + std::to_string(step) + ".vtu";
            std::ofstream ofs(name);
            vtkUnstructuredGridWrite(ofs, coil->getTetras());
            ofs.close();
            pvd.push(name, step);
         }

         // collision detection
         double range = 0.1;
         std::unordered_set<networkFace*> faces;
         for (auto p : coil->getSurfacePoints()) {
            auto tmp = box->BucketFaces.getData(p->X, range);
            faces.insert(tmp.begin(), tmp.end());
         }
         for (auto p : box->getSurfacePoints()) {
            auto tmp = coil->BucketFaces.getData(p->X, range);
            faces.insert(tmp.begin(), tmp.end());
         }

         std::string name = "collision" + std::to_string(step) + ".vtu";
         std::ofstream ofs("./output/" + name);
         vtkUnstructuredGridWrite(ofs, faces);
         ofs.close();
         pvd_collision.push(name, step);

      } while (std::ranges::any_of(coil->getPoints(), [](const networkPoint* p) { return p->RK_X.finished; }));
      simulation_time = (*(coil->getPoints().begin()))->RK_X.gett();
   }
   pvd.output();
   pvd_collision.output();
}