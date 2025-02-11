#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <vector>
#include "tetgen1.6.0/tetgen.h"
//
#include "Network.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT

点が，楕円の内部にあるかどうかを判定する．

<img src="example5_is_in_an_ellipsoid.png" width="500px">

*/

int main() {

   std::string outdir = "./example5/";
   PVDWriter pvd(outdir + "ellipsoidPoints.pvd");

   int gridCount = 100;               // Number of points along each axis
   double gridSize = 2. / gridCount;  // Grid size

   Ellipsoid ellipsoid_front(1.0, 0.5, 0.3);
   Ellipsoid ellipsoid_back(0.5, 0.25, 0.15);

   int total_step = 50;
   const Tddd initial_vector = {1., 0., 0.};
   for (int step = 0; step < total_step; ++step) {
      std::cout << Green << "step = " << step << colorReset << std::endl;
      Network* net = new Network();
      double theta = 2. * step * M_PI / total_step;
      Tddd center = {0.1 * cos(theta), 0.1 * sin(theta), 0.1 * sin(theta)};

      double r = 1.;
      double eps = 0.1;
      Tddd n_front = r * Tddd{cos(theta), std::sin(theta), 0.};
      Quaternion Q_front(initial_vector, n_front);
      ellipsoid_front.setProperties(r + eps, r / 1.5, r / 1.5, center + n_front, Q_front);

      r = 0.5;
      Tddd n_back = r * Tddd{-cos(theta), -std::sin(theta), 0.};
      Quaternion Q_back(initial_vector, {-cos(theta), -std::sin(theta), 0.});
      ellipsoid_back.setProperties(r + eps, r / 1.5, r / 1.5, center + n_back, Q_back);

      int count = 0;
      for (int i = -gridCount; i < gridCount; ++i)
         for (int j = -gridCount; j < gridCount; ++j)
            for (int k = -gridCount; k < gridCount; ++k) {
               Tddd point = {i * gridSize, j * gridSize, k * gridSize};
               if (ellipsoid_front.InsideQ(point) || ellipsoid_back.InsideQ(point)) {
                  new networkPoint(net, point);
                  count++;
               }
            }

      std::cout << "count = " << count++ << std::endl;
      std::string filename = "ellipsoidStep" + std::to_string(step) + ".vtp";
      std::ofstream ofs(outdir + filename);
      vtkPolygonWrite(ofs, net->getPoints());
      ofs.close();
      pvd.push(filename, step);
      delete net;
   }

   pvd.output();
   return 0;
}
