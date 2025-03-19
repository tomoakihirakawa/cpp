/*

## 四面体の操作

To check if the tetrahedra are correctly generated and can be accessed, we will manipulate the tetrahedra in this example.

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_manipulation_tetrahedron_simple_simulation_messed.cpp
make
./example0_manipulation_tetrahedron_simple_simulation_messed
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
   double dt = argc > 4 ? std::stod(argv[4]) : 0.001;

   /* -------------------------------------------------------------------- */
   /*                           四面体の生成と出力                            */
   /* -------------------------------------------------------------------- */

   auto coil = new Network(name0);
   auto box = new Network(name1);

   if (argc > 5) {
      // argv[3] to sting
      auto command = std::string(argv[5]);
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
         p->RK_defGrad.initialize(dt, simulation_time, p->defGrad, 4);
         p->RK_generalized_velocity.initialize(dt, simulation_time, p->velocity, 4);
      }

      std::cout << "step : " << step << std::endl;
      do {

         /* ------------------------------ velocityGrad ------------------------------ */
         for (auto p : coil->getPoints()) {
            p->velocityGrad.fill({0., 0., 0.});
            double volume = 0.;
            for (auto tet : p->Tetras) {
               p->velocityGrad += tet->volume * tet->grad([](const networkPoint* p) { return p->velocityTranslational(); });
               volume += tet->volume;
            }
            p->velocityGrad /= volume;
         }
         /* ------------------------------ CauchyStress ------------------------------ */

         for (auto p : coil->getPoints()) {

            p->RightCauchyGreen = Dot(Transpose(p->defGrad), p->defGrad);
            std::array<std::array<double, 3>, 3> I = {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}};
            p->GreenLagrangeStrain = 0.5 * (p->RightCauchyGreen - I);
            // ２階Piola-Kirchhoff応力
            double lambda = 100;
            double mu = 0;
            // S = lambda * tr(E) * I + 2 * mu * E
            std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> C;
            p->SecondPiolaKirchhoffStress.fill({0., 0., 0.});
            for (int i = 0; i < 3; ++i)
               for (int j = 0; j < 3; ++j)
                  for (int k = 0; k < 3; ++k)
                     for (int l = 0; l < 3; ++l) {
                        auto c = 0.;
                        if (i == j && k == l)
                           c = lambda;
                        if (i == k && j == l)
                           c += mu;
                        if (i == l && j == k)
                           c += mu;
                        p->SecondPiolaKirchhoffStress[i][j] += c * p->GreenLagrangeStrain[k][l];
                     }
            // P = F . S
            p->FirstPiolaKirchhoffStress = Dot(p->defGrad, p->SecondPiolaKirchhoffStress);
            // Cauchy応力
            auto J = Det(p->defGrad);
            p->CauchyStress = Dot(p->FirstPiolaKirchhoffStress, Transpose(p->defGrad)) / J;
         }

         /* ------------------------------ stressDiv ------------------------------ */

         for (auto p : coil->getPoints()) {
            p->stressDiv.fill(0.);
            double volume = 0.;
            for (auto tet : p->Tetras) {
               p->stressDiv += tet->volume * tet->div([](networkPoint* p) -> std::array<std::array<double, 3>, 3> { return p->CauchyStress; });
               volume += tet->volume;
            }
            p->stressDiv /= volume;
         }

#pragma omp parallel
         for (auto p : coil->getPoints())
#pragma omp single nowait
         {
            // std::cout << "p->RightCauchyGreen = " << p->RightCauchyGreen << std::endl;
            // std::cout << "p->GreenLagrangeStrain = " << p->GreenLagrangeStrain << std::endl;
            // std::cout << "p->CauchyStress = " << p->CauchyStress << std::endl;
            // std::cout << "p->defGrad = " << p->defGrad << std::endl;
            // std::cout << "p->velocityGrad = " << p->velocityGrad << std::endl;
            // std::cout << "p->stressDiv = " << p->stressDiv << std::endl;
            p->RK_X.push(p->velocityTranslational());
            p->RK_defGrad.push(Dot(p->velocityGrad, p->defGrad));

            p->acceleration = {0., 0., -9.8, 0., 0., 0.};

            // 衝突を模擬
            // if (p->X[2] < 0.55 && p->velocity[2] < 0) {
            //    double dist = 0.55 - p->X[2];
            //    p->acceleration[2] += 10000 * std::pow(dist, 2);
            // }

            p->acceleration[0] += p->stressDiv[0];
            p->acceleration[1] += p->stressDiv[1];
            p->acceleration[2] += p->stressDiv[2];

            // damping
            // double total_w = 0;
            // std::array<double, 3> v = {0., 0., 0.};
            // for (auto q : p->getNeighbors()) {
            //    auto w = 1. / Norm(q->X - p->X);
            //    v += w * q->velocityTranslational();
            //    total_w += w;
            // }
            // v /= total_w;
            // v -= p->velocityTranslational();

            // p->acceleration[0] += 0.1 * v[0];
            // p->acceleration[1] += 0.1 * v[1];
            // p->acceleration[2] += 0.1 * v[2];

            // damping
            p->acceleration[0] -= p->velocity[0] / dt * 0.5;
            p->acceleration[1] -= p->velocity[1] / dt * 0.5;
            p->acceleration[2] -= p->velocity[2] / dt * 0.5;

            // auto r = box->siginedDistance(p->X);
            // double k = 1. / std::pow(Norm(r), 2);
            // auto r_normalize = Normalize(r);
            // acceleration[0] -= k * r_normalize[0];
            // acceleration[1] -= k * r_normalize[1];
            // acceleration[2] -= k * r_normalize[2];
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
            p->defGrad = p->RK_defGrad.getX();
         }
         // std::cout << "collision detection" << std::endl;
         // double range = 0.1;
         // std::unordered_set<networkFace*> faces;
         // for (auto p : coil->getSurfacePoints()) {
         //    auto tmp = box->BucketFaces.getData(p->X, range);
         //    faces.insert(tmp.begin(), tmp.end());
         // }
         // for (auto p : box->getSurfacePoints()) {
         //    auto tmp = coil->BucketFaces.getData(p->X, range);
         //    faces.insert(tmp.begin(), tmp.end());
         // }

         coil->setGeometricProperties();
      } while (std::ranges::any_of(coil->getPoints(), [](const networkPoint* p) { return p->RK_X.finished; }));
      simulation_time = (*(coil->getPoints().begin()))->RK_X.gett();

#pragma omp parallel
      for (auto p : coil->getPoints())
#pragma omp single nowait
      {
         auto a = 0.8;
         std::array<std::array<double, 3>, 3> I = {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}};
         p->defGrad = a * p->RK_defGrad.getX() + (1. - a) * I;

         double restitution = 0.5;
         if (p->X[2] < 0.55 && p->velocity[2] < 0) {
            p->velocity[2] *= -restitution;
         }
      }

      /* ---------------------------------------------------------- */
      std::cout << "output" << std::endl;
#pragma omp parallel
      for (const auto& p : coil->getPoints())
#pragma omp single nowait
      {
         data1[p] = p->X[0];
         data2[p] = p->stressDiv;
         double total_vol = 0;
         p->velocityGrad.fill({0., 0., 0.});
         for (auto tet : p->Tetras) {
            auto vol = tet->volume;
            p->velocityGrad += vol * tet->grad([](const networkPoint* p) { return std::array<double, 3>{p->velocity[0], p->velocity[1], p->velocity[2]}; });
            total_vol += vol;
         }
         p->velocityGrad /= total_vol;
         // data3[p] = p->velocityGrad[0];

         data3[p] = Det(p->defGrad);
      }
      std::vector<std::tuple<std::string, DataMap>> data = {{"x", data1}, {"stressDiv", data2}, {"det F", data3}};
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