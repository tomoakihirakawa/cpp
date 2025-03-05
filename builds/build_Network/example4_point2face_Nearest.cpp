#include "Network.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 0_4_0_face_to_face_contact

### 点から面までの最短ベクトル `Nearest`

\ref{Nearest_}{`Nearest`}関数は，点から面までの最短ベクトルを求める関数である．

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example4_point2face_Nearest.cpp
make
./example4_point2face_Nearest
```

<img src="example4.png" style="display: block; margin: 0 auto; height: 300px;">

*/

std::string outdir = "./example4/";
PVDWriter pvd(outdir + "line_face_intersection.pvd");

// use filesystem

int main() {

   PVDWriter pvd(pvdfile);
   std::filesystem::path outdir = "example4";
   std::filesystem::path pvdfile = outdir / "point2face.pvd";
   auto obj = new Network("./input/sphere.obj");
   //! 平均を原点に移動
   obj->translate(-Mean(ToX(obj->getPoints())));
   obj->resetInitialX();
   for (const auto &f : obj->getFaces())
      f->setGeometricProperties(ToX(f->setPoints()));  // @ face->Lines must have been determined
   std::ofstream ofs(outdir / "sphere.vtp");
   vtkPolygonWrite(ofs, obj->getFaces());
   ofs.close();
   std::cout << outdir / "sphere.vtp" << std::endl;

   //! バケツを作成し，点と面をバケツに保存する．
   auto edge2edge = obj->getScale();
   obj->makeBuckets();

   int n = 1000;
   for (auto i = 0; i < n; i++) {
      double q = 2 * M_PI / n * i;
      auto fromX = std::array<double, 3>{edge2edge / 2.5 * std::cos(q), 0., edge2edge / 2.5 * std::sin(q)};
      std::string name = "points" + std::to_string(i) + ".vtp";
      std::ofstream ofs(outdir / name);
      vtkPolygonWrite(ofs, {(T2Tddd){fromX, Nearest(fromX, ToX(obj->getFaces()))}});
      ofs.close();
      pvd.push(name, i);
   }
   pvd.output();
   std::cout << pvdfile << std::endl;
}