#include "Network.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 0_4_face_to_face_contact

### 面と面の接触判定

\ref{basic_geometry:IntersectQ}{`IntersectQ`}関数は，交差判定には使えるが，接触判定には使えない．

**オブジェクト同士の接触**をプログラム上で定義するなら，
２面の最短距離が，ある閾値以下にある，とするのが自然な定義だろう．

#### ２面の最短距離

２つのポリゴン面上において最短距離にある２点の片方はある三角形の頂点である．
ただし，三角形が曲面を成している場合は違う．
これには，$N_{vertex}*M_{triangle} + M_{vertex}*N_{triangle}$の計算量がかかり，
また，この一つひとつの計算において，\ref{Nearest(const Tddd &X, const T3Tddd &abc)}{Nearest}のような計算を行う．
この計算は，空間分割を使って，調べる面の数を減らせば，多くの場合，実用上問題とはならない時間内で終わる．


もう一つの方法は，よりナイーブな方法で，

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example4_point2face.cpp
make
./example4_point2face
```

![./example4/anim.gif](./example4/anim.gif)

*/

// std::string outdir = "./example4/";
// PVDWriter pvd(outdir + "line_face_intersection.pvd");

// use filesystem
std::filesystem::path outdir = "example4";
std::filesystem::path pvdfile = outdir / "point2face.pvd";
PVDWriter pvd(pvdfile);
int main() {

   auto obj = new Network("./bunny.obj");
   //! 平均を原点に移動
   obj->translate(-Mean(ToX(obj->getPoints())));
   obj->resetInitialX();
   for (const auto &f : obj->getFaces())
      f->setGeometricProperties(ToX(f->setPoints()));  // @ face->Lines must have been determined
   std::ofstream ofs(outdir / "sphere.vtp");
   vtkPolygonWrite(ofs, obj->getFaces());
   ofs.close();

   //! バケツを作成し，点と面をバケツに保存する．
   auto edge2edge = obj->getScale();
   auto bucket_spacing = edge2edge / 10.;
   obj->makeBucketPoints(bucket_spacing);
   obj->makeBucketFaces(bucket_spacing);

   int n = 1000;
   for (auto i = 0; i < n; i++) {
      double q = 2 * M_PI / n * i;
      auto fromX = std::array<double, 3>{edge2edge / 2.5 * std::cos(q), 0., edge2edge / 2.5 * std::sin(q)};
      auto closestX = Nearest(fromX, ToX(obj->getFaces()));
      std::string name = "points" + std::to_string(i) + ".vtp";
      std::ofstream ofs(outdir / name);
      vtkPolygonWrite(ofs, {(T2Tddd){fromX, closestX}});
      ofs.close();
      pvd.push(name, i);
   }
   pvd.output();
}