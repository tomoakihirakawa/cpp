/*DOC_EXTRACT 0_0_load_output_3d_model

# `Network`

節点に隣接する節点や辺や要素の情報を効率的に取得するためには，接続関係を管理し続ける必要がある．`Network`クラスは，接続関係の情報を保持し，その情報をもとに相互にアクセスするための機能を提供する．

* 節点や辺や面の相互アクセス
* メッシュの細分化
* `obj`や`off`オブジェクトデータの読み込みと出力

## 点・線・面の接続関係とその整理

1. `networkFace->Lines`を設定
2. `networkFace->setPoints()`は，`networkFace->Lines`が設定されていることを前提として，`networkFace->Points`と`networkFace->PLPLPL`を設定する．
3. `Network::setGeometricProperties()`は，`f->setGeometricProperties(ToX(f->setPoints()))`を実行している．

## 3Dファイルの読み込みと出力

### 読み込み `Network`

\ref{Network::constructor}{Networkのコンストラクタ}では，引数として，**OFFファイル**または**OBJファイル**をあたえることができる．
`Load3DFile`クラスを使ってデータを読み込み，`Network`クラスを作成している．

```cpp
auto obj = new Network("./bunny.obj");//ファイルからNetworkオブジェクトを作成
```

### 出力 `vtkPolygonWrite`，`vtkUnstructuredGridWrite`

`Network`クラスは，`getFaces`メンバ関数を使って簡単に面の情報を取得できる．

`vtkPolygonWrite`を使うと，`Network`クラスの面の情報を，`vtp`ファイルとして出力できる．
`vtkPolygonWrite`には，`ofstream`と，`std::vector<networkFace*>`や`std::vector<networkLine*>`などを渡し，出力できる．

#### 面の出力

面や線や点の情報を出力する場合は，`vtkPolygonWrite`を使う．

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getFaces());
ofs.close();
```

#### 線の出力

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getEdges());
ofs.close();
```

#### 点の出力

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getPoints());
ofs.close();
```

#### 四面体の出力

四面体のような内部構造を持つデータを出力する場合は，`vtkUnstructuredGridWrite`を使う．

```cpp
auto obj = new Network("./input/bunny.off");
obj->tetrahedralize();
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
```

\ref{add_data_default_name}{このようにして}，点に値を付与し，vtpとして出力することもできる．
また，\ref{add_data_custom_name}{カスタム名}を付けることもできる．

#### 実行方法

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_load_3d_file.cpp
make
./example0_load_3d_file
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

   {
      auto off = new Network("./input/bunny.off");
      std::ofstream ofs("./output/bunny_off.vtp");
      vtkPolygonWrite(ofs, off->getFaces());
   }

   {
      auto obj = new Network("./input/bunny.obj");
      std::ofstream ofs("./output/bunny_obj.vtp");
      vtkPolygonWrite(ofs, obj->getFaces());
   }

   {
      auto line_obj = new Network("./input/line.obj");
      std::ofstream ofs("./output/line_obj.vtp");
      vtkPolygonWrite(ofs, line_obj->getLines());
   }

   /*DOC_EXTRACT 0_1_load_output_3d_model

   ### `PVDWriter`を使ったpvdファイルの作成方法

   pvdファイルは，ファイルと時間をセットにしてまとめ，paraview上で，3Dファイルのアニメーションを再生するためのファイル．

   ```cpp
   PVDWriter pvd("./bunny_obj.pvd");//出力するpvdファイル名を指定しクラスを作成
   pvd.push(filename, time);//`filename`には，`vtp`ファイルなどの3Dファイル名を，`time`には，そのファイルの時間を指定
   pvd.output();//最後にpvdファイルを出力
   ```

   | 面のアニメーション | 線のアニメーション |
   |:---------------:|:---------------:|
   | <img src="sample.gif" width="500px"> | <img src="sample_line.gif" width="500px"> |

   NOTE: QuickTimeで作成したmovファイルをgifに変換するには，次のようにする．

   ```sh
   ffmpeg -i line.mov -filter_complex "[0:v] fps=30, scale=iw*0.5:ih*0.5 [v]" -map "[v]" sample_line.gif
   ```

    */

   /* -------------------------------------------------------------------------- */
   /*                                    面の出力                                 */
   /* -------------------------------------------------------------------------- */
   {
      PVDWriter pvd("./output/bunny_obj.pvd");
      auto obj = new Network("./input/bunny.obj");
      for (auto i = 0; i < 18; ++i) {
         obj->rotate(20. / 180. * M_PI, {0., 0., 1.});
         auto filename = "bunny_obj" + std::to_string(i) + ".vtp";
         std::ofstream ofs("./output/" + filename);
         vtkPolygonWrite(ofs, obj->getFaces());
         pvd.push(filename, i);
      }
      pvd.output();
   }

   /* -------------------------------------------------------------------------- */
   /*                                    線の出力                                  */
   /* -------------------------------------------------------------------------- */
   {
      PVDWriter pvd("./output/line_obj.pvd");
      auto obj = new Network;
      const double L = 10.;
      std::vector<networkPoint*> points;
      for (const auto& x : Subdivide(0, L, 10))
         points.push_back(new networkPoint(obj, {x, 0, 0}));
      for (auto i = 0; i < points.size() - 1; ++i)
         new networkLine(obj, points[i], points[i + 1]);

      for (int i = 0; const auto& t : Subdivide(0, 2 * M_PI, 100)) {
         for (auto& p : points) {
            auto [x, y, z] = p->X;
            auto c1 = 0.5, c2 = 0.1;
            auto c = c1 * x / L + c2 * std::pow(x / L, 2);
            p->setXSingle({x, c * std::cos(t + 4 * M_PI / L * x), c * std::sin(t + 4 * M_PI / L * x)});
         }
         auto filename = "line_obj" + std::to_string(i++) + ".vtp";
         std::ofstream ofs("./output/" + filename);
         vtkPolygonWrite(ofs, obj->getLines());
         pvd.push(filename, t);
         ofs.close();
      }
      pvd.output();
   }

   /* -------------------------------------------------------------------------- */
   /*                                    点の出力                                  */
   /* -------------------------------------------------------------------------- */
   using Tddd = std::array<double, 3>;
   using DataVariant = std::variant<double, Tddd>;
   using DataMap = std::unordered_map<networkPoint*, DataVariant>;
   {
      auto obj = new Network("./input/bunny.obj");
      {
         std::string filename = "point_obj.vtp";
         std::ofstream ofs("./output/" + filename);
         vtkPolygonWrite(ofs, obj->getPoints());
         ofs.close();
      }
      {
         //@ 点に値を付与\label{add_data_default_name}
         std::unordered_map<networkPoint*, double> data;
         for (const auto& p : obj->getPoints())
            data[p] = p->X[0];
         std::ofstream ofs("./output/point_obj_with_data_default.vtp");
         vtkPolygonWrite(ofs, obj->getPoints(), data);
         ofs.close();
      }
      {
         //@ 点に値を付与\label{add_data_custom_name}
         auto data1 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
         auto data2 = std::unordered_map<networkPoint*, std::variant<double, Tddd>>();
         for (const auto& p : obj->getPoints()) {
            data1[p] = p->X[0];
            data2[p] = p->X;
         }

         std::vector<std::tuple<std::string, DataMap>> data = {{"x", data1}, {"xyz", data2}};
         std::ofstream ofs("./output/point_obj_with_data_custom.vtp");
         vtkPolygonWrite(ofs, obj->getPoints(), data);
         ofs.close();
      }
   }

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
   }
}