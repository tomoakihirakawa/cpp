#include "Network.hpp"
#include "vtkWriter.hpp"

/*
DOC_EXTRACT 0_4_face_to_face_contact

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

ただし，頂点が面上にある場合など，特異な場合の処理には注意が必要である．

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example4_point2face_interaction.cpp
make
./example4_point2face_interaction
```


![./example4/anim.gif](./example4/anim.gif)

*/

std::string outdir = "./example4/";
int main() {

   auto obj = new Network("./bunny.obj");
   //! 平均を原点に移動
   obj->translate(-Mean(ToX(obj->getPoints())));
   obj->resetInitialX();
   for (const auto &f : obj->getFaces())
      f->setGeometricProperties(ToX(f->setPoints()));  // @ face->Lines must have been determined

   //! バケツを作成し，点と面をバケツに保存する．
   auto edge2edge = obj->getScale();
   auto bucket_spacing = edge2edge / 10.;
   obj->makeBucketPoints(bucket_spacing);
   obj->makeBucketFaces(bucket_spacing);

   int n = 20;
   for (auto i = 0; i < n; i++) {
      double q = 2 * M_PI / n * i;
      auto formX = std::array<double, 3>{std::cos(q), std::sin(q), 0.};
      auto closestX = obj->BucketFaces.findClosestPoint(formX);
      std::string name = "points" + std::to_string(i) + ".vtp";
      std::ofstream ofs(outdir + name);
      vtkPolygonWrite(ofs, {(T2Tddd){fromX, closestX}});
      ofs.close();
      pvd2.push(name, i);
   }

   std::ofstream ofs(outdir + "all_intersected_faces.vtp");
   vtkPolygonWrite(ofs, all_intersected_faces);
   ofs.close();

   pvd.output();
   pvd2.output();
}