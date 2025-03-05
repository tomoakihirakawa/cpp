#include "Network.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 0_3_line_face_interaction

## 空間分割の応用例：オブジェクトの接触や交差の判定

### 線分と面の交差判定

`Network`クラスは，`makeBucketPoints`でバケツ`BucketPoints`を準備し，内部に保存している点をバケツに保存する．
同様に，`makeBucketFaces`でバケツを`BucketFaces`を準備し，内部に保存している面をバケツに保存する．

要素の接触や交差の判定には，\ref{basic_geometry:IntersectQ}{`IntersectQ`}関数を使う．
また，接触判定の高速化のために，空間分割を使う．

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example3_line_face_interaction.cpp
make
./example3_line_face_interaction
```

<gif src="./example3/anim_faster.gif" width="500px">

*/

std::string outdir = "./example3/";
PVDWriter pvd(outdir + "line_face_intersection.pvd");
PVDWriter pvd2(outdir + "line.pvd");

int main() {

   auto obj = new Network("./bunny.obj");
   obj->translate(-Mean(ToX(obj->getPoints())));
   obj->resetInitialX();
   for (const auto &f : obj->getFaces())
      f->setGeometricProperties(ToX(f->setPoints()));  // @ face->Lines must have been determined
   {
      std::ofstream ofs(outdir + "bunny.vtp");
      vtkPolygonWrite(ofs, obj->getFaces());
      ofs.close();
   }
   //! バケツを作成し，点をバケツに保存する．
   auto spacing = obj->getScale() / 10.;
   obj->makeBucketPoints(spacing);
   //! バケツを作成し，面をバケツに保存する．
   obj->makeBucketFaces(spacing);

   const double r = 0.1;
   std::unordered_set<networkFace *> all_intersected_faces;
   const int N = 1000;
   for (auto i = 0; i < N; i++) {
      auto theta = 4 * M_PI * i / N;
      std::unordered_set<networkFace *> intersect_faces;
      Tddd A = {r * std::cos(theta), -0.06 + 0.15 / N * i, r * std::sin(theta)};
      Tddd B = {-r * std::cos(theta), -0.06 + 0.15 / N * i, -r * std::sin(theta)};
      obj->BucketFaces.apply(obj->BucketFaces.line2indices(A, B), [&](networkFace *f) {
         if (IntersectQ(T2Tddd{A, B}, (T3Tddd)(*f))) {
            std::cout << "intersect" << std::endl;
            intersect_faces.emplace(f);
            all_intersected_faces.emplace(f);
         }
      });
      {
         std::string name = "intersected_faces" + std::to_string(i) + ".vtp";
         std::ofstream ofs(outdir + name);
         vtkPolygonWrite(ofs, all_intersected_faces);
         ofs.close();
         pvd.push(name, i);
      }
      {
         std::string name = "line" + std::to_string(i) + ".vtp";
         std::ofstream ofs(outdir + name);
         vtkPolygonWrite(ofs, {(T2Tddd){A, B}});
         ofs.close();
         pvd2.push(name, i);
      }
   }

   std::ofstream ofs(outdir + "all_intersected_faces.vtp");
   vtkPolygonWrite(ofs, all_intersected_faces);
   ofs.close();

   pvd.output();
   pvd2.output();
}