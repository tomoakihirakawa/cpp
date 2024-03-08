#include "Network.hpp"
#include "basic.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 0_0_0_quaternion

# クォータニオンを使った物体の３次元回転

## クォータニオンを使ったシンプルな回転

クォータニオンは，3D回転を効率的に計算するために便利な表現．

\insert{0_0_0_Quaternion}

 * 以下のコードは、3Dオブジェクトの回転と平行移動を実行します．
 * translate関数：ネットワークの全点を指定した量だけ平行移動します．
 * rotate関数：指定したクォータニオンと中心点を使用して、ネットワークの全点を回転します．
 * main関数：bunny、cow、camelオブジェクトをロードして、それぞれを回転させ、結果をファイルに出力します．

```
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=validateRotation.cpp
make
./validateRotation
```

 ![sample.gif](sample.gif)

$x$軸に対して回転 -> $y$軸に対して回転 -> $z$軸に対して回転

$(0.1,0,0)$を中心にして$x$軸に対して回転 -> $y$軸に対して回転 -> $z$軸に対して回転

時計回りが正である．

*/

void translate(Network* const net, const Tddd& shift) {
   for (auto& p : net->getPoints())
      p->setXSingle(p->initialX + shift);
   net->setGeometricProperties();
};

void rotate(Network* const net, const Quaternion& Q, const Tddd& c = {0, 0, 0}) {
   for (auto& p : net->getPoints())
      p->setXSingle(Q.Rv(p->initialX - c) + c);
   net->setGeometricProperties();
};

void rotate_axis(Network* const net, const Quaternion& Q, const Tddd& c = {0, 0, 0}) {
   for (auto& p : net->getPoints())
      p->setXSingle(Q.Rs(p->initialX - c) + c);
   net->setGeometricProperties();
};

int main() {
   std::string outdir = "./output/";
   int t = 0;
   if (true) {
      PVDWriter PVD(outdir + "rotating_bunny.pvd");
      auto net = new Network("bunny.obj", "bunny");
      auto mean = Mean(ToX(net->getPoints()));
      translate(net, -mean);
      net->resetInitialX();
      int l = 0;
      t = 0;
      for (const auto& center : std::vector<Tddd>{{0, 0, 0}, {0.1, 0., 0.}}) {
         // 各軸に対して回転させる
         for (const auto& axis : std::vector<Tddd>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}) {
            for (auto i = 0; i < 50; ++i) {
               auto Q = Quaternion(axis, 2 * M_PI / 50. * i);
               vtkPolygonWriter<networkPoint*> vtp;
               rotate_axis(net, Q, center);
               for (const auto& f : net->getFaces()) {
                  vtp.add(f->getPoints());
                  vtp.addPolygon(f->getPoints());
               }
               std::cout << l << ", {yaw,pitch,roll} = " << Q.YPR() << std::endl;
               std::string name = "bunny" + std::to_string(l++) + ".vtp";
               PVD.push(name, t);
               t += 1;
               std::ofstream ofs(outdir + name);
               vtp.write(ofs);
            }
         }
      }
      PVD.output();
   }

   if (true) {
      PVDWriter PVD(outdir + "rotating_cow.pvd");
      auto net = new Network("cow.obj", "cow");
      auto mean = Mean(ToX(net->getPoints()));
      translate(net, -mean);
      net->resetInitialX();
      int l = 0;
      Quaternion Q;
      t = 0;
      for (auto i = 0; i < 50; ++i) {
         vtkPolygonWriter<networkPoint*> vtp;
         // Q *= Quaternion({ 1, 0, 0 }, 2 * M_PI / 50.);
         Q *= Quaternion({0, 0, 1}, 2 * M_PI / 50.);
         rotate(net, Q);
         for (const auto& f : net->getFaces()) {
            vtp.add(f->getPoints());
            vtp.addPolygon(f->getPoints());
         }
         std::cout << l << ", {yaw,pitch,roll} = " << Q.YPR() << std::endl;
         std::string name = "armadillo" + std::to_string(l++) + ".vtp";
         PVD.push(name, t);
         t += 1;
         std::ofstream ofs(outdir + name);
         vtp.write(ofs);
      }
      PVD.output();
   }

   if (true) {
      PVDWriter PVD(outdir + "rotating_camel.pvd");
      auto net = new Network("camel.obj", "camel");
      auto mean = Mean(ToX(net->getPoints()));
      translate(net, -mean);
      net->resetInitialX();
      int l = 0;
      t = 0;
      Quaternion Qyaw, Qpitch, Qroll;
      for (auto i = 0; i < 50; ++i) {
         vtkPolygonWriter<networkPoint*> vtp;
         Qyaw *= Quaternion({0, 0, 1}, 2 * M_PI / 50.);
         Qpitch *= Quaternion({0, 1, 0}, 2 * M_PI / 50.);
         Qroll *= Quaternion({1, 0, 0}, 2 * M_PI / 50.);
         auto Q = Quaternion({0, 0, 1}, Qyaw.yaw());
         Q *= Quaternion({0, 1, 0}, Qpitch.pitch());
         Q *= Quaternion({1, 0, 0}, Qroll.roll());
         rotate(net, Q);
         for (const auto& f : net->getFaces()) {
            vtp.add(f->getPoints());
            vtp.addPolygon(f->getPoints());
         }
         std::cout << l << ", {yaw,pitch,roll} = " << Q.YPR() << std::endl;
         std::string name = "camel" + std::to_string(l++) + ".vtp";
         PVD.push(name, t);
         t += 1;
         std::ofstream ofs(outdir + name);
         vtp.write(ofs);
      }
      PVD.output();
   }
};
