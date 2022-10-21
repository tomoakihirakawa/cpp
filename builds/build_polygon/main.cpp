#include "Network.hpp"
#include "vtkWriter.hpp"
int main() {

   std::string name = "bunny";
   auto net = new Network("./bunny.obj", name);
   vtkPolygonWrite("./bunny.vtp", ToX(net->getFaces()));

   vtkPolygonWriter<networkPoint *> vtp;
   auto [x01, y01, z01] = net->getBounds();

   for (const auto &x : Subdivide(std::get<0>(x01), std::get<1>(x01), 100))
      for (const auto &y : Subdivide(std::get<0>(y01), std::get<1>(y01), 100))
         for (const auto &z : Subdivide(std::get<0>(z01), std::get<1>(z01), 100)) {
            Tddd X = {x, y, z};
            if (net->isInside(X)) {
               auto p = new networkPoint(net, net, X);
               vtp.add(p);
            }
         }

   std::ofstream ofs("./out.vtp");
   vtp.write(ofs);
   ofs.close();
};
