/*DOC_EXTRACT 0_1_manipulation_tetrahedron

## 四面体の操作

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_manipulation_tetrahedron.cpp
make
./example0_manipulation_tetrahedron
```

* \ref{make_bucket_tetras}{これは}，空間分割のためのバケットを作成し，四面体をバケットに登録する例
* \ref{edge_flip}{これは}，四面体を持つ表面のエッジフリップテストの例
* \ref{test_insideQ}{これは}，座標が四面体の内部か外部かの判定テストの例

*/

#include "tetgen1.6.0/tetgen.h"
//
#include "Network.hpp"
#include "vtkWriter.hpp"

Tddd DistorsionMeasureWeightedSmoothingVector_modified(const networkPoint *p, std::function<Tddd(const networkPoint *)> position) {
  const int max_sum_depth = 20;
  auto faces = p->getFaces();
  std::vector<double> weights;
  weights.reserve(faces.size());
  std::vector<Tddd> positions;
  positions.reserve(faces.size());
  Tddd X0, X1, X2, Xmid, vertical, X_ideal, To_ideal;
  double W, height, normalized_discrepancy;
  for (const auto &f : faces) {
    auto [p0, p1, p2] = f->getPoints(p);
    X0 = position(p0);
    X1 = position(p1);
    X2 = position(p2);
    W = std::pow(CircumradiusToInradius(X0, X1, X2), 2); // CircumradiusToInradius(X0, X1, X2)は最小で2
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

using Tddd = std::array<double, 3>;
using DataVariant = std::variant<double, Tddd>;
using DataMap = std::unordered_map<networkPoint *, DataVariant>;

auto lines = std::string(50, '=');

int main(int argc, char *argv[]) {
  std::string fileName = "./input/tetrahedron_input.stl";
  if (argc >= 2) {
    fileName = argv[1];
  }

  /* -------------------------------------------------------------------- */
  /*                           四面体の生成と出力                            */
  /* -------------------------------------------------------------------- */
  std::cout << Red << lines << colorReset << std::endl;
  std::cout << Red << "1. 四面体の生成と出力" << colorReset << std::endl;

  auto obj = new Network(fileName);
  obj->tetrahedralize();
  //@ 点に値を付与\label{add_data_custom_name}
  std::unordered_map<networkPoint *, double> data1;
  std::unordered_map<networkPoint *, Tddd> data2;
  std::unordered_map<networkPoint *, std::variant<double, Tddd>> data3;
  for (const auto &p : obj->getPoints()) {
    data1[p] = p->X[0];
    data2[p] = p->X;
    data3[p] = (double)p->Tetras.size();
  }
  using NamedData = std::tuple<std::string, DataVariant3>;
  std::vector<NamedData> data;
  data.emplace_back("x", DataVariant3{data1});
  data.emplace_back("xyz", DataVariant3{data2});
  data.emplace_back("tetra_size", DataVariant3{data3});
  std::ofstream ofs("./outptut/tetras.vtu");
  vtkUnstructuredGridWrite(ofs, obj->getTetras(), data);
  ofs.close();
  std::cout << "Wrote tetras.vtu" << std::endl;

  /* --------------------- 空間分割のためのバケットを作成し，四面体をバケットに登録する --------------------- */
  std::cout << Red << lines << colorReset << std::endl;
  std::cout << Red << "2. 空間分割のためのバケットを作成し，四面体をバケットに登録する" << colorReset << std::endl;

  // \label{make_bucket_tetras}
  obj->makeBucketTetras(obj->getScale() / 5);

  int i = 0;
  PVDWriter pvd("./output/tetras_in_bucket.pvd");
  for (const auto &VVVdata : obj->BucketTetras.data)
    for (const auto &VVdata : VVVdata)
      for (const auto &Vdata : VVdata) {
        auto name = "cell" + std::to_string(i) + ".vtu";
        std::ofstream ofs("./output/" + name);
        vtkUnstructuredGridWrite(ofs, Vdata);
        ofs.close();
        pvd.push(name, i);
        std::cout << "Vdata.size() : " << Vdata.size() << ", i : " << i << std::endl;
        i++;
      }

  pvd.output();

  /* ---------------------------- 四面体を持つ表面のエッジフリップテスト ---------------------------- */
  std::cout << Red << lines << colorReset << std::endl;
  std::cout << Red << "3. 四面体を持つ表面のエッジフリップテスト" << colorReset << std::endl;

  // \label{edge_flip}
  {
    {
      std::ofstream ofs("./output/tetras_before_flip.vtu");
      vtkUnstructuredGridWrite(ofs, obj->getTetras());
      ofs.close();
      std::cout << "Wrote tetras_before_flip.vtu" << std::endl;
    }
    int count = 0;
    for (auto i = 0; i < 50; ++i) {
      for (const auto &f : RandomSample(obj->getSurfaces())) {
        bool found = false;
        for (const auto &l : f->getLines()) {
          auto [p, q] = l->getPoints();
          data1[p] = -10000.; // わかりやすいように色をつけるフリップ前-10000，フリップ後10000
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

    obj->tetrahedralize(); // flip後は必ず，四面体の情報を再構築する必要がある．
    std::vector<NamedData> data;
    data.emplace_back("fliped", DataVariant3{data1});
    data.emplace_back("xyz", DataVariant3{data2});
    data.emplace_back("tetra_size", DataVariant3{data3});

    std::ofstream ofs("./output/tetras_after_flip.vtu");
    vtkUnstructuredGridWrite(ofs, obj->getTetras(), data);
    ofs.close();
    std::cout << "Wrote tetras_after_flip.vtu" << std::endl;
  }

  /* --------------------------- エッジのマージテストも同様に実装可能 --------------------------- */

  std::cout << Red << lines << colorReset << std::endl;
  std::cout << Red << "4. エッジのマージテストも同様に実装可能" << colorReset << std::endl;

  {
    // delete all inner faces;
    std::vector<networkLine *> inner_lines;
    for (const auto &l : obj->getLines())
      if (std::ranges::none_of(l->Faces, [&](const auto &f) { return f->SurfaceQ(); }))
        inner_lines.push_back(l);
    for (const auto &l : inner_lines)
      delete l;

    obj->setGeometricProperties();
  }

  // とにかく出力
  {
    obj->setGeometricProperties();
    std::ofstream ofs("./output/tetras_after_merge.vtp");
    vtkPolygonWrite(ofs, obj->getSurfaces());
    ofs.close();
    std::cout << "Wrote tetras_after_merge.vtp" << std::endl;
  }

  {

    int count = 0;
    for (auto i = 0; i < 2000; ++i) {
      for (const auto &f : RandomSample(obj->getSurfaces())) {
        bool found = false;
        for (const auto &l : f->getLines()) {
          auto p = l->merge();
          if (p) {
            std::cout << Green << "Merged!" << colorReset << std::endl;
            obj->checkConnectivity();
            found = true;
            count++;
            break;
          }
        }
        if (found) {
          break;
        }
      }
    }

    // とにかく出力
    {
      obj->setGeometricProperties();
      std::ofstream ofs("./output/tetras_after_merge.vtp");
      vtkPolygonWrite(ofs, obj->getSurfaces());
      ofs.close();
      std::cout << "Wrote tetras_after_merge.vtp" << std::endl;
    }

    obj->tetrahedralize(); // flip後は必ず，四面体の情報を再構築する必要がある．
    std::vector<NamedData> data;
    data.emplace_back("fliped", DataVariant3{data1});
    data.emplace_back("xyz", DataVariant3{data2});
    data.emplace_back("tetra_size", DataVariant3{data3});

    std::ofstream ofs("./output/tetras_after_merge.vtu");
    vtkUnstructuredGridWrite(ofs, obj->getTetras(), data);
    ofs.close();
    std::cout << "Wrote tetras_after_merge.vtu" << std::endl;
  }

  /* ---------------------------- 座標が四面体の内部か外部かの判定テスト ---------------------------- */
  std::cout << Red << lines << colorReset << std::endl;
  std::cout << Red << "5. 座標が四面体の内部か外部かの判定テスト" << colorReset << std::endl;

  // \label{test_insideQ}
  {
    PVDWriter pvd("./output/tetras_InsideQ.pvd");
    auto bounds = obj->scaledBounds(0.7);
    auto [xyz0, xyz1] = Transpose(bounds);
    int N = 100;
    for (auto i = 0; i < N; ++i) {
      std::vector<networkTetra *> tetras;
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
  std::cout << Red << lines << colorReset << std::endl;
  std::cout << Red << "6. 四面体の表面の点を移動させるテスト" << colorReset << std::endl;

  // \label{surface_tetras_move_points}
  {
    /*
       表面の点を移動させて，四面体を調整していくことができるかのテスト
    */
    PVDWriter pvd("./output/tetras_surface_move_points.pvd");
    auto normal = [&](const networkPoint *p) {
      Tddd N = {0., 0., 0.};
      for (const auto &f : p->getSurfaces())
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
            p->setXSingle(p->X + a * DistorsionMeasureWeightedSmoothingVector_modified(p, [&](const networkPoint *p) { return p->X; }));
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
    pvd.output();
  }
}