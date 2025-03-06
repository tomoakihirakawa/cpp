#include "Network.hpp"
#include "vtkWriter.hpp"
/*DOC_EXTRACT 0_4_1_face_to_face_contact

### 面と面の接触判定（面と面の最短距離の計算）

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example5_face2face_contact.cpp
make
./example5_face2face_contact
```

\ref{basic_geometry:IntersectQ}{`IntersectQ`}関数は，交差判定には使えるが，接触判定には使えない．
オブジェクト同士の**接触**をプログラム上で定義するなら，互いの面において最も近くにある面同士の最短距離を計算が，ある閾値以下にあるときに接触しているとみなす方法が自然である．

\ref{Nearest(const T3Tddd &XYZ, const T3Tddd &ABC)}{Nearest}関数は，面と面の最短距離を求める関数．
重心座標付近を徐々に細かく探索して，より正確な最短距離を求める．探査範囲を狭めながらこれを繰り返す．
`N=3`，4x4点で最短となる重心座標を求める．次に範囲を0.7倍に狭め，これを数十回繰り返す．

<img src="example5_face2face_contact.png" style="display: block; margin: 0 auto; height: 300px;">

<img src="example5_face2face_contact.gif" style="display: block; margin: 0 auto; height: 300px;">

<img src="example5_face2face_contact2.gif" style="display: block; margin: 0 auto; height: 300px;">

> [!WARNING]
> 捕捉した側の三角形の節点に1を，捕捉された側の三角形の節点には-1の値を与えて出力している．捕捉される側の点が-1ではないからといって，接触していないわけではない．
> 接触判定は，２つ同時に行うことは難しい．オブジェクト毎に接触判定を行う必要がある．

#### 接触判定の条件の制限

場合によっては，面同士が向き合っている場合にのみ接触したとみなしたい場合もある．例えば次の２つが考えられる．

- `isFlat(X1nearest - X0nearest, F->normal, 10. * M_PI / 180.)`
- `isFlat(f->normal, F->normal, 10. * M_PI / 180.)`

*/

int main() {

   std::filesystem::path outdir = "example5";
   PVDWriter pvd(outdir / "face_face_interaction.pvd");
   std::filesystem::path pvdfile = outdir / "face2face.pvd";
   auto obj = new Network("./input/sphere.obj");
   //! 平均を原点に移動
   obj->translate(-Mean(ToX(obj->getPoints())));
   obj->resetInitialX();

   for (const auto &f : obj->getSurfaces())
      f->setGeometricProperties(ToX(f->setPoints()));  // @ face->Lines must have been determined
   std::ofstream ofs(outdir / "sphere.vtp");
   vtkPolygonWrite(ofs, obj->getSurfaces());
   ofs.close();
   std::cout << outdir / "sphere.vtp" << std::endl;

   //! バケツを作成し，点と面をバケツに保存する．
   auto edge2edge = obj->getScale();
   obj->makeBuckets();

   int n = 100;

   for (auto i = 0; i < n; i++) {
      double q = 2 * M_PI / n * i;
      auto fromC = std::array<double, 3>{edge2edge / 2.5 * std::cos(q), 0., edge2edge / 2.5 * std::sin(q)};
      auto fromX = std::array<double, 3>{0., 0.1, 0.} + fromC;
      auto fromY = std::array<double, 3>{0., -0.1, 0.} + fromC;
      auto fromZ = std::array<double, 3>{0.1, 0., 0.1} + fromC;
      T3Tddd XYZ = {fromX, fromY, fromZ};
      Tddd X0, X1, X0nearest, X1nearest;
      double min_distance = 1E+10;
      for (const auto &f : obj->getSurfaces()) {
         auto ABC = ToX(f);
         auto [X0, X1] = Nearest(XYZ, ABC);
         double d = Norm(X0 - X1);
         if (min_distance > d) {
            X0nearest = X0;
            X1nearest = X1;
            min_distance = d;
         }
      }
      std::string name = "points" + std::to_string(i) + ".vtp";
      std::ofstream ofs(outdir / name);
      // auto [X0, X1] = Nearest(T3Tddd{fromX, fromY, fromZ}, ToX(obj->getSurfaces()));
      T6Tddd V = {fromX, fromY, fromZ, fromX, X0nearest, X1nearest};
      vtkPolygonWrite(ofs, V);
      ofs.close();
      pvd.push(name, i);
   }
   pvd.output();
   std::cout << pvdfile << std::endl;

   {
      PVDWriter pvd(outdir / "object_object_interaction.pvd");
      auto OBJ = new Network("./input/sphere.obj");
      std::vector<networkFace *> faces;
      faces.reserve(obj->getSurfaces().size() + OBJ->getSurfaces().size());
      for (const auto &f : obj->getSurfaces())
         faces.emplace_back(f);
      for (const auto &f : OBJ->getSurfaces())
         faces.emplace_back(f);
      std::unordered_map<networkPoint *, double> hit;
      hit.reserve(obj->getPoints().size() + OBJ->getPoints().size());

      Tddd mean = Mean(ToX(OBJ->getPoints()));
      OBJ->translate(-mean);
      OBJ->rotate(-M_PI, Tddd{0., 1., 0.});
      OBJ->translate(Tddd{0.0065, 0., -0.05});
      OBJ->resetInitialX();
      OBJ->makeBuckets();

      for (auto &p : obj->getPoints())
         hit[p] = 0.;
      for (auto &p : OBJ->getPoints())
         hit[p] = 0.;

      Tddd X0, X1, X0nearest, X1nearest;
      const int N = 100;
      double range = 0.002;
      const Tdd ranges = {-1.5 * range, 1.5 * range};
      const T3Tdd ranges3 = {ranges, ranges, ranges};
      for (auto i = 0; i < N; i++) {
         std::cout << "i = " << i << std::endl;
         OBJ->translate(0.0002 * Tddd{0., std::cos(M_PI / N), 0.});
         // OBJ->rotate(M_PI / N, Tddd{1., 1., 1.});
         OBJ->resetInitialX();
         OBJ->makeBuckets();
         std::cout << "object is translated" << std::endl;
         for (auto &p : hit)
            p.second = 0.;

         std::cout << "traverse faces and find nearest faces" << std::endl;
         for (const auto &F : OBJ->getSurfaces()) {
            for (const auto &f : obj->getSurfaces()) {
               if (IntersectQ(F->bounds + ranges3, f->bounds + ranges3)) {
                  auto [X0, X1] = Nearest(ToX(F), ToX(f));
                  if (Norm(X0 - X1) < range && isFlat(Normalize(X1 - X0), F->normal, 10. * M_PI / 180.)) {
                     for (const auto &p : F->getPoints()) hit[p] = 1.;
                     for (const auto &p : f->getPoints()) hit[p] = -1.;
                     std::cout << "hit" << std::endl;
                     break;
                  }
               }
            }
         }
         std::string name = "faces" + std::to_string(i) + ".vtp";
         std::ofstream ofs(outdir / name);
         std::vector<std::tuple<std::string, decltype(hit)>> data{{"hit", hit}};
         vtkPolygonWrite(ofs, faces, data);
         ofs.close();
         std::cout << outdir / name << std::endl;
         pvd.push(name, i);
      }
      pvd.output();
      std::cout << pvdfile << std::endl;
   }
}