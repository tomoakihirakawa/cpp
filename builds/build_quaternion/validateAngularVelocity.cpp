#include "Network.hpp"
#include "basic.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 1_0_0_quaternion

\insert{0_1_0_Quaternion}

## クォータニオンの微分の数値的な時間積分の例

以下を実行して，ルンゲクッタを使いクォータニオンの時間微分を時間積分する．

```
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=validateAngularVelocity.cpp
make
./validateAngularVelocity
```

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
   const int iteration = 10;
   std::string outdir = "./output/";
   for (const auto normaliaze : {1, 2, 3})
      for (const auto RK_order : {1, 2, 4}) {

         std::string id = "RK_" + std::to_string(RK_order);
         if (normaliaze == 1)
            id = id + "_normalized_at_end";
         else if (normaliaze == 2)
            id = id + "_normalized_at_middle";
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
                  auto w = axis * 2 * M_PI / iteration;
                  auto dqdt = AngularVelocityTodQdt(w, net->RK_Q.getX());
                  net->RK_Q.push(dqdt);  // クォータニオン->T4dとしてプッシュ
                  net->Q = normaliaze == 2 ? Normalize(net->RK_Q.getX()) : net->RK_Q.getX();
                  std::cout << "Q = " << net->Q() << ", Det(net->Q.R()) = " << Det(net->Q.R()) << std::endl;
               }
               if (normaliaze == 1)
                  net->Q.normalize();
               /* -------------------------------------------------------------------------- */
               for (auto& p : net->getPoints())
                  p->setXSingle(net->rigidTransformation(p->initialX));
               // p->setXSingle(rigidTransformation(net->ICOM, net->COM, net->Q.R(), p->initialX));
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
   std::cout << "paraview ./output/cow_time_integratedRK_4_not_normalized.pvd ./output/cow_time_integratedRK_4_normalized_at_end.pvd ./output/cow_time_integratedRK_4_normalized_at_middle.pvd" << std::endl;
};
