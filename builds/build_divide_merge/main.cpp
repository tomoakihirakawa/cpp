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

#define bem
#if defined(bem)
int main(int arg, char **argv) {
   Timer time;
   // JSON json(std::ifstream("./input.json"));
   // Network net(json["objfile"][0], json["name"][0]);
   // auto name = json["name"][0];
   // int remesh = stoi(json["remesh"])[0];
   /* ------------------------------------------------------ */
   // JSON json(std::ifstream("./input.json"));
   std::string input_file{argv[1]};        // input
   std::string output_directory{argv[2]};  // output dir
   std::string output_name{argv[3]};       // output name
   int remesh = std::atoi(argv[4]);
   Network net(input_file, output_name);
   mk_vtu(output_directory + "/" + output_name + ".vtu", {net.getFaces()});
   /* ------------------------------------------------------ */
   net.displayStates();

   auto extractLength = [](const V_netLp &lines) {
      V_d ret;
      for (auto l : lines)
         ret.emplace_back(l->length());
      return ret;
   };
   //* ------------------------------------------------------ */
   //*                          線の分割                       */
   //* ------------------------------------------------------ */
   V_netLp lines;
   Histogram Histo;
   for (auto count = 0; count <= remesh; ++count) {
      //! ------------------------------------------------------ */
      mk_vtu(output_directory + "/" + output_name + std::to_string(count) + ".vtu", net.getFaces());
      std::ofstream ofs(output_directory + "/" + output_name + std::to_string(count) + ".obj");
      creteOBJ(ofs, net);
      ofs.close();
      std::cout << "output " << time() << std::endl;
      //! ------------------------------------------------------ */
      std::cout << Grid({"count", Histo.count}, 50) << std::endl;
      std::cout << Grid({"cumulative_count", Histo.cumulative_count}, 50) << std::endl;
      std::cout << Grid({"diff", Histo.diff}, 50) << std::endl;
      std::cout << Grid({"interval", Histo.interval}, 50) << std::endl;
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
      Divide(net.Lines, [&](auto l) { return l->length() > Histo.interval[num]; });
      //
      for (auto i = 0; i < 10; i++)
         Merge(net.Lines, [&](auto l) { return l->length() < Histo.interval[num]; });
      // (*net.Lines.begin())->merge();
      //
      std::cout << net << std::endl;
      std::cout << net.getPoints().size() << std::endl;
      std::cout << net.Lines.size() << std::endl;
      /*
      ! 注意
      @ 強い制限のため，湾曲した境界面上を点は移動することができない．
      @ しかし，円筒の内面などは滑らかにしたいができない．
      @ 湾曲した面の辺をフリップしてしまうと，そのまま動かなくなる．これは面の法線方向を変化させるような平滑化を禁止しているためで，フリップにはその制限がかかっていないためである．
      */
      const double small = M_PI / 180 * 0.01;
      for (auto i = 0; i < 1; i++) {
         AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         for (const auto &l : net.Lines)
            l->flipIfTopologicalyBetter(0.5 * M_PI / 180., 0.5 * M_PI / 180.);
         AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         for (const auto &l : net.Lines)
            l->flipIfBetter(M_PI / 180.);
      }

      std::cout << "1. time:" << time() << std::endl;
      AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
      for (const auto &l : net.Lines)
         l->flipIfBetter(M_PI / 180.);
      std::cout << "2. time:" << time() << std::endl;
      AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
      std::cout << "3. time:" << time() << std::endl;
      LaplacianSmoothingPreserveShape(net.getPoints(), small);
   }
}
#elif defined(sphare)
int main() {
   // https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/
   const double X = .525731112119133606f;
   const double Z = .850650808352039932f;
   const double N = 0.;

   const VV_d vertices =
       {{-X, N, Z},
        {X, N, Z},
        {-X, N, -Z},
        {X, N, -Z},
        {N, Z, X},
        {N, Z, -X},
        {N, -Z, X},
        {N, -Z, -X},
        {Z, X, N},
        {-Z, X, N},
        {Z, -X, N},
        {-Z, -X, N}};

   const VV_i connection =
       {{0, 4, 1},
        {0, 9, 4},
        {9, 5, 4},
        {4, 5, 8},
        {4, 8, 1},
        {8, 10, 1},
        {8, 3, 10},
        {5, 3, 8},
        {5, 2, 3},
        {2, 7, 3},
        {7, 10, 3},
        {7, 6, 10},
        {7, 11, 6},
        {11, 0, 6},
        {0, 1, 6},
        {6, 1, 10},
        {9, 0, 11},
        {9, 11, 2},
        {9, 2, 5},
        {7, 2, 11}};

   Network net;
   net.generateNetwork(vertices, connection);
   net.setBounds();
   net.displayStates();

   double lim_len = .2;
   auto extractLength = [](const V_netLp &lines) {
      V_d ret;
      for (auto l : lines)
         ret.emplace_back(l->length());
      return ret;
   };
   //* ------------------------------------------------------ */
   //*                          線の分割                        */
   //* ------------------------------------------------------ */
   for (auto count = 0; count < 10; ++count) {
      //! ------------------------------------------------------ */
      mk_vtu("./vtu/sphere_divide" + std::to_string(count) + ".vtu", net.getFaces());
      std::ofstream ofs(home_dir + "/vtu/sphere_divide" + std::to_string(count) + ".obj");
      creteOBJ(ofs, net);
      ofs.close();

      //! ------------------------------------------------------ */
      Histogram h(extLength(net.getLines()));
      std::cout << Grid({"count", h.count}, 50) << std::endl;
      std::cout << Grid({"cumulative_count", h.cumulative_count}, 50) << std::endl;
      std::cout << Grid({"diff", h.diff}, 50) << std::endl;
      std::cout << Grid({"interval", h.interval}, 50) << std::endl;

      int num = 0;
      /* ------------------------*/
      //* | 0 | 1 | 2 |[3]| 4 |   diff, mid_interval
      //! 0   1   2  [3]  4   5   interval, cumulative_count
      /* ------------------------*/
      for (auto j = 0; j < h.diff.size(); ++j)
         if (h.cumulative_count[j + 1] / h.data_size >= 0.8 && h.diff[j] >= h.diff[num == 0 ? j : num])
            num = j;
      std::cout << "num " << num << std::endl;
      std::cout << h.interval[num] << std::endl;
      if (!num == 0) {

         bool found = false;
         do {
            found = false;
            for (const auto &l : net.getLines()) {
               std::cout << l->length() << ", " << h.interval[num + 1] << std::endl;
               if (l->length() >= h.interval[num]) {
                  std::cout << 1 << std::endl;
                  auto p = l->divide();
                  std::cout << 2 << std::endl;
                  p->setX(p->getX() / Norm(p->getX()));
                  std::cout << 3 << std::endl;
                  found = true;
                  std::cout << 4 << std::endl;
                  for (auto &p : Flatten(BFS(p, 3)))
                     for (auto &L : p->getLines())
                        L->flipIfIllegal();
                  std::cout << 5 << std::endl;
               }
            }
         } while (found);
      }
   }
   //* ------------------------------------------------------ */
   //*                          面のマージ                      */
   //* ------------------------------------------------------ */
   int cc = 0;
   for (auto count = 0; count < 10; ++count) {
      //! ------------------------------------------------------ */
      Histogram h(extLength(net.getLines()));
      std::cout << Grid({"count", h.count}, 50) << std::endl;
      std::cout << Grid({"cumulative_count", h.cumulative_count}, 50) << std::endl;
      std::cout << Grid({"diff", h.diff}, 50) << std::endl;
      std::cout << Grid({"interval", h.interval}, 50) << std::endl;

      int num = 0;
      /* ------------------------*/
      //* | 0 | 1 | 2 |[3]| 4 |   diff, mid_interval
      //! 0   1   2  [3]  4   5   interval, cumulative_count
      /* ------------------------*/
      for (auto j = 0; j < h.diff.size(); ++j)
         if (h.cumulative_count[j] / h.data_size <= 0.8 && h.diff[j] >= std::abs(h.diff[num]))
            num = j;
      if ((h.diff[num] && !num == h.interval.size()) > 0 || num == 0)
         num++;
      std::cout << "num " << num << std::endl;
      std::cout << h.interval[num] << std::endl;
      try {
         if (num != h.interval.size()) {
            bool found = false;
            do {
               found = false;
               std::cout << "size = " << net.getLines().size() << std::endl;
               net.displayStates();
               /* ------------------------------------------------------ */
               // デバッグのために，適当に削除していく　
#define good_merge
#if defined(bad_merge)
               for (const auto &l : net.getLines()) {
                  auto p = l->merge();
                  p->setX(p->getX() / Norm(p->getX()));
                  for (auto &p : Flatten(BFS(p, 3)))
                     for (auto &L : p->getLines())
                        L->flipIfIllegal();  // フリップがないと歪な三角形ができる．
                  found = true;
                  mk_vtu("./vtu/sphere_merge_0" + std::to_string(cc++) + ".vtu", net.getFaces());
                  break;
               }
#elif defined(good_merge)
               /* ------------------------------------------------------ */
               // 比較的大きい方から
               for (const auto &l : net.getLines()) {
                  if (!l) {
                     std::cout << "l = " << l << std::endl;
                     std::cout << l->length() << ", " << h.interval[num + 1] << std::endl;
                  }
                  if (l->length() <= h.interval[num]) {
                     auto p = l->merge();
                     p->setX(p->getX() / Norm(p->getX()));
                     found = true;
                     for (auto &p : Flatten(BFS(p, 3)))
                        for (auto &L : p->getLines())
                           L->flipIfIllegal();
                     //! ------------------------------------------------------ */
                     mk_vtu(home_dir + "/vtu/sphere_merge_0" + std::to_string(cc++) + ".vtu", net.getFaces());
                     // std::ofstream ofs("./vtu/sphere_merge" + std::to_string(count) + ".obj");
                     // creteOBJ(ofs, net);
                     // ofs.close();
                     break;
                  }
                  if (found)
                     break;
               }
#endif
            } while (found);
         }
      } catch (std::exception e) {
         std::cout << e.what() << "\n";
      }
   }
}
#endif