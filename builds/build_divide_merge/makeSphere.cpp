#include "InterpolationRBF.hpp"
#include "Network.hpp"

std::string home_dir = std::getenv("HOME");

bool refine(netLp l, double len) {
   if (l->length() > len) {
      l->divide();
      return true;
   }
   return false;
};

#define sphare

#if defined(sphare)
int main() {
   Timer time;
   // https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/
   const double X = .525731112119133606f;
   const double Z = .850650808352039932f;
   const double N = 0.;

   const VV_d vertices =
       {{-X, N, Z}, {X, N, Z}, {-X, N, -Z}, {X, N, -Z}, {N, Z, X}, {N, Z, -X}, {N, -Z, X}, {N, -Z, -X}, {Z, X, N}, {-Z, X, N}, {Z, -X, N}, {-Z, -X, N}};

   const VV_i connection =
       {{0, 4, 1}, {0, 9, 4}, {9, 5, 4}, {4, 5, 8}, {4, 8, 1}, {8, 10, 1}, {8, 3, 10}, {5, 3, 8}, {5, 2, 3}, {2, 7, 3}, {7, 10, 3}, {7, 6, 10}, {7, 11, 6}, {11, 0, 6}, {0, 1, 6}, {6, 1, 10}, {9, 0, 11}, {9, 11, 2}, {9, 2, 5}, {7, 2, 11}};

   Network net;
   net.setFaces(connection, net.setPoints(vertices));
   net.displayStates();

   const double radius = 0.15;
   for (const auto &p : net.getPoints())
      p->setX(radius * Normalize(p->X));

   //* ------------------------------------------------------ */
   //*                          線の分割                        */
   //* ------------------------------------------------------ */
   Histogram Histo;
   for (auto count = 0; count < 100; ++count) {
      //! ------------------------------------------------------ */
      if (count % 1 == 0) {
         mk_vtu("./output/sphere_divide" + std::to_string(count) + ".vtu", net.getFaces());
         std::ofstream ofs("./output/sphere_divide" + std::to_string(count) + ".obj");
         creteOBJ(ofs, net);
         ofs.close();
      }
      //! ------------------------------------------------------ */
      Histogram h(extLength(net.getLines()));
      //   std::cout << Grid({"count", h.count}, 50) << std::endl;
      //   std::cout << Grid({"cumulative_count", h.cumulative_count}, 50) << std::endl;
      //   std::cout << Grid({"diff", h.diff}, 50) << std::endl;
      //   std::cout << Grid({"interval", h.interval}, 50) << std::endl;
      int num = 0;
      /* ------------------------*/
      //* | 0 | 1 | 2 |[3]| 4 |     diff, mid_interval
      //!   0   1   2  [3]  4   5   interval, cumulative_count
      /* ------------------------*/
      Histo.set(extLength(net.getLines()));
      std::cout << "time:" << time() << std::endl;
      for (auto j = 0; j < Histo.cumulative.size(); ++j)
         if (Histo.cumulative[j] >= 0.8) {
            num = j;
            break;
         }

      std::cout << "Histo.interval[" << num << "]" << Histo.interval[num] << std::endl;
      auto tmp = net.Lines;
      if (!num == 0) {
         bool found = false;
         do {
            found = false;
            for (const auto &l : net.getLines()) {
               //    std::cout << l->length() << ", " << h.interval[num + 1] << std::endl;
               if (l->length() >= h.interval[num]) {
                  auto p = l->divide();
                  p->setX(radius * Normalize(p->X));
                  found = true;
                  break;
               }
            }
         } while (found);

         const double small = M_PI / 180 * 0.1;
         //  for (auto i = 0; i < 1; i++) {
         //     AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         //     for (const auto &l : net.Lines)
         //        l->flipIfTopologicalyBetter(0.5 * M_PI / 180., 0.5 * M_PI / 180.);
         //     AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         //     for (const auto &l : net.Lines)
         //        l->flipIfBetter(M_PI / 180.);
         //  }

         AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         for (const auto &l : net.Lines)
            l->flipIfBetter(M_PI / 180.);
         AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         LaplacianSmoothingPreserveShape(net.getPoints(), small);
         //
         for (const auto &p : net.getPoints())
            p->setX(radius * Normalize(p->X));
      }
   }
}
#endif