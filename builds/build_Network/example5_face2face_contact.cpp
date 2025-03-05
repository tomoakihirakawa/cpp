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

\ref{Nearest(const T3Tddd &XYZ, const T3Tddd &ABC)}{Nearest}関数は，面と面の最短距離を求める関数であり，次の流れで計算される．

1. 準ニュートン法を使って，凡その最短距離となる２面の重心座標（計4変数）を求める．
2. その重心座標付近を細かく探索して，より正確な最短距離を求める．探査範囲を狭めながらこれを繰り返す．

始めの，準ニュートン法で`1/5`の精度で重心座標が求まると仮定している．分割探査だけを使うと見落としが生じる**気がするので**含めている．
これ以降は，`N=2`，３x3点点で最短となる重心座標を求める．8回繰り返すと，焼く1E-10の精度で最短距離が求まる．

<img src="example5_face2face_contact.png" style="display: block; margin: 0 auto; height: 300px;">

*/

int main() {

   std::filesystem::path outdir = "example5";
   PVDWriter pvd(outdir / "face_face_interaction.pvd");
   std::filesystem::path pvdfile = outdir / "face2face.pvd";
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
      for (const auto &f : obj->getFaces()) {
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
      // auto [X0, X1] = Nearest(T3Tddd{fromX, fromY, fromZ}, ToX(obj->getFaces()));
      T6Tddd V = {fromX, fromY, fromZ, fromX, X0nearest, X1nearest};
      vtkPolygonWrite(ofs, V);
      ofs.close();
      pvd.push(name, i);
   }
   pvd.output();
   std::cout << pvdfile << std::endl;
}