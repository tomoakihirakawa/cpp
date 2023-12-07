#include "Network.hpp"
#include "basic.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 1_0_0_quaternion

## クォータニオンの時間微分，角速度

以下を実行して，ルンゲクッタを使いクォータニオンの時間微分を時間積分する．

```
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=validateAngularVelocity.cpp
make
./validateAngularVelocity
```

\insert{0_1_0_Quaternion}

![sample_dQdt.gif](sample_dQdt.gif)

SEE: 青がRK1，ピンクがRK2．緑がRK4．
RK1は時間が進むにつれて大きくなっている．
RK4も初期の状態と比べて若干大きくなっているように見える．

### クォータニオンの正規化

以上のような誤差は．回転行列を作成する際に，クォータニオンを正規化していないため，
回転行列の行列式が1になっていないことが原因である．
回転行列の行列式は，スケーリング（拡大縮小）を表しておリ，体積や長さを保つためには，行列式は1でなければならない．

**必ず回転行列を計算する際は，正規化したクォータニオンを使うべきである．**

## 剛体の回転と平行移動

\insert{rigidTransformation}

*/

void translate(Network* const net, const Tddd& shift) {
   for (auto& p : net->getPoints())
      p->setXSingle(p->initialX + shift);
   net->setGeometricProperties();
};

int main() {
   const double dt = 1.;
   const int iteration = 50;
   std::string outdir = "./output/";
   for (const auto normaliaze : {true, false})
      for (const auto RK_order : {1, 2, 4}) {

         std::string id = "RK_" + std::to_string(RK_order);
         if (normaliaze)
            id = id + "_normalized";
         else
            id = id + "_not_normalized";

         auto net = new Network("cow.obj", "cow");
         auto mean = Mean(ToX(net->getPoints()));
         translate(net, -mean);
         net->resetInitialX();
         net->COM = net->initial_center_of_mass = {0., 0., 0.};

         double simulation_time = 1.;
         //
         int l = 0;
         // 各軸に対して回転させる
         PVDWriter PVD(outdir + "cow_time_integrated" + id + ".pvd");
         for (const auto& axis : std::vector<Tddd>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}) {
            for (auto i = 0; i < iteration; ++i) {
               net->RK_Q.initialize(dt, simulation_time, net->Q(), RK_order);
               /* -------------------------------------------------------------------------- */
               // RK法で時間積分
               for (auto j = 0; j < RK_order; ++j) {
                  net->RK_Q.push(AngularVelocityTodQdt(axis * 2 * M_PI / iteration, net->Q));  // クォータニオン->T4dとしてプッシュ
                  net->Q = net->RK_Q.getX();
                  std::cout << "Q = " << net->Q() << ", Det(net->Q.R()) = " << Det(net->Q.R()) << std::endl;
                  //! ここに正規化を入れると，時間積分がうまくいかない．
               }
               if (normaliaze) net->Q.normalize();
               /* -------------------------------------------------------------------------- */
               for (auto& p : net->getPoints())
                  p->setXSingle(rigidTransformation(net->ICOM, net->COM, net->Q.R(), p->initialX));
               // p->setXSingle(net->Q.Rv(p->initialX - net->ICOM) + net->ICOM);
               net->setGeometricProperties();
               //
               vtkPolygonWriter<networkPoint*> vtp;
               for (const auto& f : net->getFaces()) {
                  vtp.add(f->getPoints());
                  vtp.addPolygon(f->getPoints());
               }

               std::string name = "cow_time_integrated" + id + "_" + std::to_string(l++) + ".vtp";
               PVD.push(name, simulation_time);
               std::ofstream ofs(outdir + name);
               vtp.write(ofs);

               simulation_time += dt;
            }
         }
         PVD.output();
         std::cout << "Done." << std::endl;
         std::cout << "paraview " + outdir + "cow_time_integrated" + id + ".pvd" << std::endl;
      }
};
