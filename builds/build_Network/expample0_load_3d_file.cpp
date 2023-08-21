#include "Network.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT load_3d_file

# `Network`

## 3Dファイルを読み込み，`vtkPolygonWrite`を使った出力方法

### 読み込み `Network`

\ref{Network::constructor}{Networkのコンストラクタ}では，拡張子から，
与えられたファイルが，

* OFFファイル
* OBJファイル

かチェクして，`Load3DFile`クラスを使ってデータを読み込み`Network`クラスとして読み込む．

### 出力 `vtkPolygonWrite`

`vtkPolygonWrite`には，`ofstream`と，`std::vector<networkFace*>`や`std::vector<networkTetra*>`などを渡し，出力できる．

```cpp
std::ofstream ofs("./bunny_obj.vtp");
vtkPolygonWrite(ofs, obj->getFaces());
```

![sample.png](sample.png)

```shell
$ cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=load_3d_file.cpp
$ make
$ ./load_3d_file
```

*/

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

   /*DOC_EXTRACT pvd

   ### `PVDWriter`を使ったpvdファイルの作成方法

   pvdファイルは，ファイルと時間をセットにしてまとめ，paraview上で，3Dファイルのアニメーションを再生するためのファイルである．

   次のようにして，出力するpvdファイル名を指定しクラスを作成する．

   ```cpp
   PVDWriter pvd("./bunny_obj.pvd");
   ```

   `filename`には，`vtp`ファイルなどの3Dファイル名を，`time`には，そのファイルの時間を指定する：

   ```cpp
   pvd.push(filename, time);
   ```

   最後にpvdファイルを出力する．

   ```cpp
   pvd.output();
   ```

   ![sample.gif](sample.gif)

    */

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
