/*DOC_EXTRACT 0_0_load_output_3d_model

# `Network`

#### 実行方法

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

Tddd DistorsionMeasureWeightedSmoothingVector_modified(const networkPoint* p, std::function<Tddd(const networkPoint*)> position) {
   const int max_sum_depth = 20;
   auto faces = p->getFaces();
   std::vector<double> weights;
   weights.reserve(faces.size());
   std::vector<Tddd> positions;
   positions.reserve(faces.size());
   Tddd X0, X1, X2, Xmid, vertical, X_ideal, To_ideal;
   double W, height, normalized_discrepancy;
   for (const auto& f : faces) {
      auto [p0, p1, p2] = f->getPoints(p);
      X0 = position(p0);
      X1 = position(p1);
      X2 = position(p2);
      W = std::pow(CircumradiusToInradius(X0, X1, X2), 2);  // CircumradiusToInradius(X0, X1, X2)は最小で2
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3.) * 0.5;
      X_ideal = height * vertical + Xmid;
      To_ideal = X_ideal - position(p0);
      normalized_discrepancy = Norm(To_ideal) / height;
      weights.push_back(normalized_discrepancy * W);
      positions.emplace_back(To_ideal);
   }
   weights /= Sum(weights);
   Tddd V = {0., 0., 0.};
   for (int i = 0; i < weights.size(); ++i)
      V += weights[i] * positions[i];
   return V;
};

int main() {

   /* -------------------------------------------------------------------- */
   /*                           四面体の生成と出力                            */
   /* -------------------------------------------------------------------- */
   {
      auto obj = new Network("./input/bunny.off");
      obj->tetrahedralize();
      //@ 点に値を付与\label{add_data_custom_name}
      auto data1 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
      auto data2 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
      auto data3 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
      for (const auto& p : obj->getPoints()) {
         data1[p] = p->X[0];
         data2[p] = p->X;
         data3[p] = (double)p->Tetras.size();
      }
      std::vector<std::tuple<std::string, DataMap>> data = {{"x", data1}, {"xyz", data2}, {"tetra_size", data3}};
      std::ofstream ofs("./outptut/tetras.vtu");
      vtkUnstructuredGridWrite(ofs, obj->getTetras(), data);
      ofs.close();
      std::cout << "Wrote tetras.vtu" << std::endl;

      /* ---------------------------- 四面体を持つ表面のエッジフリップテスト ---------------------------- */

      {
         //@ 一部の辺をflipしたあと，四面体を出力
         int count = 0;
         for (auto i = 0; i < 100; ++i) {
            for (const auto& f : RandomSample(obj->getSurfaces())) {
               bool found = false;
               for (const auto& l : f->getLines()) {
                  auto [p, q] = l->getPoints();
                  data1[p] = -10000.;  // わかりやすいように色をつけるフリップ前-10000，フリップ後10000
                  data1[q] = -10000.;
                  if (l->flip()) {
                     auto [p, q] = l->getPoints();
                     data1[p] = 10000.;
                     data1[q] = 10000.;
                     std::cout << Green << "Flipped!" << colorReset << std::endl;
                     found = true;
                     count++;
                     break;
                  }
               }
               if (found)
                  break;
            }
         }
         std::vector<std::tuple<std::string, DataMap>> data = {{"fliped", data1}, {"xyz", data2}, {"tetra_size", data3}};

         std::ofstream ofs("./output/tetras_after_flip.vtu");
         vtkUnstructuredGridWrite(ofs, obj->getTetras(), data);
         ofs.close();
         std::cout << "Wrote tetras_after_flip.vtu" << std::endl;
      }

      /* ---------------------------- 座標が四面体の内部か外部かの判定テスト ---------------------------- */
      {
         PVDWriter pvd("./output/tetras_InsideQ.pvd");
         auto bounds = obj->scaledBounds(0.7);
         auto [xyz0, xyz1] = Transpose(bounds);
         int N = 1000;
         for (auto i = 0; i < N; ++i) {
            std::vector<networkTetra*> tetras;
            auto a = (double)i / N;
            auto X = a * xyz1 + (1 - a) * xyz0;
            for (auto tet : obj->getTetras())
               if (tet->InsideQ(X))
                  tetras.emplace_back(tet);
            if (tetras.size() > 0) {
               auto filename = "tetras_InsideQ" + std::to_string(i) + ".vtu";
               std::ofstream ofs("./output/" + filename);
               vtkUnstructuredGridWrite(ofs, tetras, data);
               pvd.push(filename, i);
               ofs.close();
               if (i == 0)
                  std::cout << "Wrote " << filename << std::endl;
            }
         }
         pvd.output();
      }

      /* ---------------------------- 四面体の表面の点を移動させるテスト ---------------------------- */
      {
         /*
            表面の点を移動させて，四面体を調整していくことができるかのテスト
         */
         PVDWriter pvd("./output/tetras_surface_move_points.pvd");
         auto normal = [&](const networkPoint* p) {
            Tddd N = {0., 0., 0.};
            for (const auto& f : p->getSurfaces())
               N += f->normal;
            return Normalize(N);
         };

         double a = 0.1;
         for (auto i = 0; i < 100; ++i) {
            //@ 表面の点を移動
            for (auto p : obj->getSurfacePoints())
               p->setXSingle(p->X + normal(p) * 0.001);

            //@ 四面体の節点の調整
            for (auto j = 0; j < 5; ++j)
               for (auto p : obj->getPoints()) {
                  if (!p->SurfaceQ()) {
                     p->setXSingle(p->X + a * DistorsionMeasureWeightedSmoothingVector_modified(p, [&](const networkPoint* p) { return p->X; }));
                  }
               }

            obj->setGeometricProperties();

            auto filename = "tetras_surface_move_points" + std::to_string(i) + ".vtu";
            std::ofstream ofs("./output/" + filename);
            vtkUnstructuredGridWrite(ofs, obj->getTetras(), data);
            pvd.push(filename, i);
            ofs.close();
            std::cout << "Wrote " << "./output/" + filename << std::endl;
         }
      }
   }
}