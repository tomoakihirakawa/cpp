/*DOC_EXTRACT 0_0_load_output_3d_model

# `Network`

数値シミュレーションの多くは，三角形や四面体の３Dメッシュを利用して行われる．

ある節点に隣接する節点や辺や要素を素早く効率的に取得するためには，節点や辺や面の接続関係をわかりやすく管理する必要がある．
`Network`クラスは，節点や辺や面の接続関係を保持し，その接続関係から相互にアクセスできるようにするために作ったクラスである．
また，メッシュを細分化したり，一度作成した`Network`クラスのオブジェクトから`vtk`や`obj`ファイルなどを出力したりもできる．

* 節点や辺や面の相互アクセス
* メッシュの細分化
* データの読み込みと出力

## 点・線・面の接続関係とその整理

1. `networkFace->Lines`を設定
2. `networkFace->setPoints()`は，`networkFace->Lines`が設定されていることを前提として，`networkFace->Points`と`networkFace->PLPLPL`を設定する．
3. `Network::setGeometricProperties()`は，`f->setGeometricProperties(ToX(f->setPoints()))`を実行している．

## 3Dファイルの読み込みと出力

### 読み込み `Network`

\ref{Network::constructor}{Networkのコンストラクタ}では，引数として，**OFFファイル**または**OBJファイル**をあたえることができる．
`Load3DFile`クラスを使ってデータを読み込み，`Network`クラスを作成する．

```cpp
auto obj = new Network("./bunny.obj");//ファイルからNetworkオブジェクトを作成
```

### 出力 `vtkPolygonWrite`

`Network`クラスは，`getFaces`メンバ関数を使って簡単に面の情報を取得できる．

`vtkPolygonWrite`を使うと，`Network`クラスの面の情報を，`vtp`ファイルとして出力できる．
`vtkPolygonWrite`には，`ofstream`と，`std::vector<networkFace*>`や`std::vector<networkLine*>`などを渡し，出力できる．

#### 面の出力

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getFaces());
ofs.close();
```

<img src="sample.png" width="500px">

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

#include "Network.hpp"
#include "vtkWriter.hpp"

int main() {

   {
      auto off = new Network("./bunny.off");
      std::ofstream ofs("./bunny_off.vtp");
      vtkPolygonWrite(ofs, off->getFaces());
   }

   {
      auto obj = new Network("./bunny.obj");
      std::ofstream ofs("./bunny_obj.vtp");
      vtkPolygonWrite(ofs, obj->getFaces());
   }

   {
      auto line_obj = new Network("./line.obj");
      std::ofstream ofs("./line_obj.vtp");
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
      PVDWriter pvd("./bunny_obj.pvd");
      auto obj = new Network("./bunny.obj");
      for (auto i = 0; i < 18; ++i) {
         obj->rotate(20. / 180. * M_PI, {0., 0., 1.});
         auto filename = "./bunny_obj" + std::to_string(i) + ".vtp";
         std::ofstream ofs(filename);
         vtkPolygonWrite(ofs, obj->getFaces());
         pvd.push(filename, i);
      }
      pvd.output();
   }

   /* -------------------------------------------------------------------------- */
   /*                                    線の出力                                  */
   /* -------------------------------------------------------------------------- */
   {
      PVDWriter pvd("./line_obj.pvd");
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
         auto filename = "./line_obj" + std::to_string(i++) + ".vtp";
         std::ofstream ofs(filename);
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
      auto obj = new Network("./bunny.obj");
      {
         std::ofstream ofs("./point_obj.vtp");
         vtkPolygonWrite(ofs, obj->getPoints());
         ofs.close();
      }
      {
         //@ 点に値を付与\label{add_data_default_name}
         std::unordered_map<networkPoint*, double> data;
         for (const auto& p : obj->getPoints())
            data[p] = p->X[0];
         std::ofstream ofs("./point_obj_with_data_default.vtp");
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
         std::ofstream ofs("./point_obj_with_data_custom.vtp");
         vtkPolygonWrite(ofs, obj->getPoints(), data);
         ofs.close();
      }
   }
}
