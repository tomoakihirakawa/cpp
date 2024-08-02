#include "Network.hpp"
#include "basic.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 0_0_0_quaternion

# クォータニオンを使った物体の３次元回転

\insert{0_0_0_Quaternion}

## クォータニオンを使った物体の３次元回転の例

3Dの物体の回転と平行移動を行う例．

 * translate：ネットワークの全点を指定した量だけ平行移動します．
 * rotate：指定したクォータニオンと中心点を使用して、ネットワークの全点を回転します．
 * main：bunny、cow、camelオブジェクトをロードして、それぞれを回転させ、結果をファイルに出力します．

```
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=validateRotation.cpp
make
./validateRotation
```

<img src="./sample.gif" alt="sample" width="300" height="200">

クォータニオンから作られた回転行列`Rv`を使って，

1. 物体を$`x`$軸に回転 -> $`y`$軸に回転 -> $`z`$軸の順に回転させる．
2. 次に，$`(0.1,0,0)`$を中心にして$`x`$軸に対して回転 -> $`y`$軸に対して回転 -> $`z`$軸に対して回転させる．

`Rv`は自身の目線は変えないまま（自身の座標系を変えないまま），物体をその座標系において回転させる．

*/

void translate(Network* const net, const Tddd& shift) {
   for (auto& p : net->getPoints())
      p->setXSingle(p->initialX + shift);
   net->setGeometricProperties();
};

void rotate(Network* const net, const Quaternion& Q, const Tddd& c = {0, 0, 0}) {
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
               rotate(net, Q, center);
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
