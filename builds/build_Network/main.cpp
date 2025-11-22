#include "Network.hpp"

int main() {
   auto net = new Network("../../obj/bunny.obj");
   mk_vtu("./vtu/obj.vtu", net->getFaces());
   for (const auto &f : net->getFaces()) {
      auto h = std::sqrt(f->area) / 5. /*粒子間隔*/;
      for (const auto &X : particlize(f, 0.001 / 2))
         new networkPoint(net, net, X);
   }
   mk_vtu("./vtu/obj_particlized.vtu", {net->getPoints()});
   delete net;
}
