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
   double k = argc > 3 ? std::stod(argv[3]) : 1.;

   /* -------------------------------------------------------------------- */
   /*                           四面体の生成と出力                            */
   /* -------------------------------------------------------------------- */

   auto coil = new Network(name0);
   auto box = new Network(name1);

   if (argc > 4) {
      // argv[3] to sting
      auto command = std::string(argv[4]);
      coil->tetrahedralize(command);
      std::cout << "tetrahedralize coil" << std::endl;
      box->tetrahedralize(command);
      std::cout << "tetrahedralize box" << std::endl;
   } else {
      coil->tetrahedralize();
      std::cout << "tetrahedralize coil" << std::endl;
      box->tetrahedralize();
      std::cout << "tetrahedralize box" << std::endl;
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

   PVDWriter pvd("./output/simu_coil.pvd");
   PVDWriter pvd_collision("./output/collision.pvd");

   box->makeBuckets();
   /* -------------------------------------------------------------------------- */
   auto data1 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
   auto data2 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
   auto data3 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
   for (const auto& p : coil->getPoints()) {
      data1[p] = 0.;
      data2[p] = Tddd{0., 0., 0.};
      data3[p] = 0.;
   }
   /* -------------------------------------------------------------------------- */

   for (auto step = 0; step < 3000; ++step) {
      for (const auto& p : coil->getPoints()) {
         p->RK_X.initialize(dt, simulation_time, p->X, 4);
         p->RK_generalized_velocity.initialize(dt, simulation_time, p->velocity, 4);
      }

      std::cout << "step : " << step << std::endl;

      coil->makeBuckets();
      coil->setContactFaces({box});

      do {

#pragma omp parallel
         for (auto p : coil->getPoints())
#pragma omp single nowait
         {
            p->RK_X.push(p->velocityTranslational());

            p->acceleration = {0., 0., -9.8, 0., 0., 0.};
            p->stress = {0., 0., 0.};

            /* --------------------------------------------------- */

            for (const auto& [f, _, __] : p->ContactFaces) {
               auto n = f->normal;
               auto a = Dot(p->velocityTranslational() / dt, n) * n;
               p->acceleration[0] -= a[0];
               p->acceleration[1] -= a[1];
               p->acceleration[2] -= a[2];
            }

            /* --------------------------------------------------- */

            for (auto t : p->Tetras) {
               auto [p0, p1, p2, p3] = t->getPoints(p);

               //! 初期位置
               auto n_init = TriangleNormal(p0->initialX, p1->initialX, p2->initialX);
               auto V_init = p0->initialX - (p1->initialX + p2->initialX + p3->initialX) / 3;
               auto d_N_init = Dot(V_init, n_init);               //! 初期位置のN成分
               auto d_H_init = Norm(V_init - d_N_init * n_init);  //! 初期位置のH成分

               //! 現在の位置
               auto n = TriangleNormal(p0->X, p1->X, p2->X);
               auto V = p0->X - (p1->X + p2->X + p3->X) / 3;
               auto d_N = Dot(V, n);          //! 現在のN成分
               auto d_H = Norm(V - d_N * n);  //! 現在のH成分
               //
               p->stress += k * std::pow(Norm(d_N_init - d_N), 2) * n;
               p->stress += k * std::pow(Norm(d_H_init - d_H), 2) * Normalize(V - d_N * n);
            }
            p->acceleration[0] += p->stress[0];
            p->acceleration[1] += p->stress[1];
            p->acceleration[2] += p->stress[2];
            /* --------------------------------------------------- */
            // damping
            p->acceleration[0] -= p->velocity[0] / dt * 0.05;
            p->acceleration[1] -= p->velocity[1] / dt * 0.05;
            p->acceleration[2] -= p->velocity[2] / dt * 0.05;

            p->RK_generalized_velocity.push(p->acceleration);
         }

#pragma omp parallel
         for (auto p : coil->getPoints())
#pragma omp single nowait
         {
            auto tmp = p->RK_generalized_velocity.getX();
            p->velocity[0] = tmp[0];
            p->velocity[1] = tmp[1];
            p->velocity[2] = tmp[2];
            p->setXSingle(p->RK_X.getX());
         }

         coil->setGeometricProperties();
      } while (std::ranges::any_of(coil->getPoints(), [](const networkPoint* p) { return p->RK_X.finished; }));
      simulation_time = (*(coil->getPoints().begin()))->RK_X.gett();

      /* ---------------------------------------------------------- */
      std::cout << "output" << std::endl;
#pragma omp parallel
      for (const auto& p : coil->getPoints())
#pragma omp single nowait
      {
         // data1[p] = p->X[0];
         data1[p] = (double)(p->ContactFaces.size());
         data2[p] = Tddd{p->acceleration[0], p->acceleration[1], p->acceleration[2]};
         data3[p] = p->stress;
      }
      std::vector<std::tuple<std::string, DataMap>> data = {{"x", data1}, {"accel", data2}, {"stress", data3}};
      std::string name = "simu_coil" + std::to_string(step) + ".vtu";
      std::ofstream ofs("./output/" + name);

      /* ---------------------------------------------------------- */

      vtkUnstructuredGridWrite(ofs, coil->getTetras(), data);
      ofs.close();
      pvd.push(name, step);

      /* ---------------------------------------------------------- */

      pvd.output();
      pvd_collision.output();
   }
}