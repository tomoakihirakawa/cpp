#include "Network.hpp"

std::string home_dir = std::getenv("HOME");

bool refine(netLp l, double len) {
   if (l->length() > len) {
      l->divide();
      return true;
   }
   return false;
};

int main(int arg, char **argv) {
   Timer time;
   std::string input_file{argv[1]};        // input
   std::string output_directory{argv[2]};  // output dir
   std::string output_name{argv[3]};       // output name
   int remesh = std::atoi(argv[4]);
   //
   Network net(input_file, output_name);
   net.setGeometricProperties();
   mk_vtu(output_directory + "/" + output_name + ".vtu", {net.getFaces()});
   net.displayStates();
   //
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
      Histo.set(extLength(net.getLines()));
      std::cout << "time:" << time() << std::endl;
      for (auto j = 0; j < Histo.cumulative.size(); ++j)
         if (Histo.cumulative[j] >= 0.8) {
            num = j;
            break;
         }
      std::cout << "Histo.interval[" << num << "]" << Histo.interval[num] << std::endl;
      Divide(net.Lines, [&](auto l) { return l->length() > Histo.interval[num]; });
      // for (auto i = 0; i < 10; i++)
      //    Merge(net.Lines, [&](auto l) { return l->length() < Histo.interval[num]; });

      const double small = M_PI / 180 * 0.01;
      for (auto i = 0; i < 1; i++) {
         AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         for (const auto &l : net.Lines)
            l->flipIfTopologicalyBetter(0.5 * M_PI / 180., 0.5 * M_PI / 180.);
         AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         for (const auto &l : net.Lines)
            l->flipIfBetter(M_PI / 180.);
      }
   }
}