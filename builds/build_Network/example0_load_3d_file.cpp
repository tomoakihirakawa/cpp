/*DOC_EXTRACT 0_0_load_output_3d_model

# `Network`

数値シミュレーションの多くは，三角形や四面体の３Dメッシュを利用して行われる．
単純に配列にメッシュ情報を格納していては，シミュレーションの質を高めるための工夫を加えることが難しいと思われる．
例えば，ある節点に隣接する節点や辺や要素を取得するのは，配列に格納しているだけでは効率的に行うことができない．
`Network`クラスは，節点や辺や面の接続関係を保持し，接続関係から相互にアクセスできるようにするためのクラスである．
また，メッシュを細分化することもできる．
また，一度作成した`Network`クラスのオブジェクトから，`vtk`ファイルや`obj`ファイルなどを出力することができる．

* データの読み込みと出力
* 節点や辺や面の相互アクセス
* メッシュの細分化

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
```

![sample.png](sample.png)

#### 線の出力

```cpp
auto obj = new Network("./bunny.obj");
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getEdges());
```

#### 実行方法

```shell
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
   | ![sample.gif](sample.gif) | ![sample_line.gif](sample_line.gif) |

   NOTE: QuickTimeで作成したmovファイルをgifに変換するには，次のようにする．

   ```sh
   ffmpeg -i line.mov -filter_complex "[0:v] fps=30, scale=iw*0.5:ih*0.5 [v]" -map "[v]" sample_line.gif
   ```

    */

   // pvd for bunny.obj
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

   // pvd for line.obj
   {
      PVDWriter pvd("./line_obj.pvd");
      auto obj = new Network;
      const double L = 10.;
      auto X = Subdivide({0, L}, 10);
      std::vector<networkPoint*> points;
      for (const auto& x : X)
         points.push_back(new networkPoint(obj, {x, 0, 0}));
      for (auto i = 0; i < points.size() - 1; ++i)
         new networkLine(obj, points[i], points[i + 1]);

      int i = 0;
      for (const auto& t : Subdivide({0, 2 * M_PI}, 100)) {
         for (auto& p : points) {
            auto [x, y, z] = p->X;
            auto c1 = 0.5, c2 = 0.1;
            auto c = c1 * x / L + c2 * std::pow(x / L, 2);
            p->setXSingle({x, c * cos(t + 4 * M_PI / L * x), c * sin(t + 4 * M_PI / L * x)});
         }
         auto filename = "./line_obj" + std::to_string(i++) + ".vtp";
         std::ofstream ofs(filename);
         vtkPolygonWrite(ofs, obj->getLines());
         pvd.push(filename, t);
      }
      pvd.output();
   }
}
