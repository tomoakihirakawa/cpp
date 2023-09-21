#include "Network.hpp"
#include "basic.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 1_0_0_quaternion

## クォータニオンの時間微分，角速度

\insert{0_1_0_Quaternion}

![sample_dQdt.gif](sample_dQdt.gif)

青がRK1，ピンクがRK2．緑がRK4．
RK1は時間が進むにつれて大きくなっている．

*/

void translate(Network* const net, const Tddd& shift) {
   for (auto& p : net->getPoints())
      p->setXSingle(p->initialX + shift);
   net->setGeometricProperties();
};

int main() {

   auto net = new Network("cow.obj", "cow");
   auto mean = Mean(ToX(net->getPoints()));
   translate(net, -mean);
   net->resetInitialX();
   net->COM = net->initial_center_of_mass = {0., 0., 0.};

   const int iteration = 50;
   const double dt = 1.;
   double simulation_time = 1.;
   const int RK_order = 1;

   int l = 0;
   // 各軸に対して回転させる
   for (const auto& axis : std::vector<Tddd>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}) {
      for (auto i = 0; i < iteration; ++i) {
         net->RK_Q.initialize(dt, simulation_time, net->Q(), RK_order);
         /* -------------------------------------------------------------------------- */
         // RK法で時間積分
         for (auto j = 0; j < RK_order; ++j) {
            net->RK_Q.push(AngularVelocityTodQdt(axis * 2 * M_PI / iteration, net->Q));  // クォータニオン->T4dとしてプッシュ
            net->Q = net->RK_Q.getX();
         }
         /* -------------------------------------------------------------------------- */
         for (auto& p : net->getPoints())
            p->setXSingle(net->Q.Rv(p->initialX - net->ICOM) + net->ICOM);
         net->setGeometricProperties();
         //
         vtkPolygonWriter<networkPoint*> vtp;
         for (const auto& f : net->getFaces()) {
            vtp.add(f->getPoints());
            vtp.addPolygon(f->getPoints());
         }

         std::ofstream ofs("./output/cow_time_integrated" + std::to_string(RK_order) + "_" + std::to_string(l++) + ".vtp");

         vtp.write(ofs);
      }
   }
};
