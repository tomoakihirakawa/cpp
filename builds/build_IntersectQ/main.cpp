#include "Network.hpp"
#include "vtkWriter.hpp"

int main() {
   std::vector<Tddd> vec;
   for (const auto &x : Subdivide(0, 1, 3))
      for (const auto &y : Subdivide(0, 1, 3))
         for (const auto &z : Subdivide(0.5, 0.51, 3))
            vec.emplace_back(Tddd{x + RandomReal({-1, 1}), y + RandomReal({-1, 1}), z + RandomReal({-.1, .1})});

   Network cube("./cube.obj", "cube");

   for (auto i = 0; i < 100; ++i) {
      cube.translate(Tddd{RandomReal({-.2, .2}), RandomReal({-.2, .2}), 0.});
      T2Tddd AB;
      std::vector<T2Tddd> lines, not_lines;
      for (const auto &A : vec)
         for (const auto &B : vec) {
            AB = {A, B};
            if (IntersectQ(cube.bounds, AB))
               lines.emplace_back(AB);
            else
               not_lines.emplace_back(AB);
         }
      {
         std::ofstream ofs("./output/intersecting_lines" + std::to_string(i) + ".vtp");
         vtkPolygonWrite(ofs, lines);
      }
      {
         std::ofstream ofs("./output/not_intersecting_lines" + std::to_string(i) + ".vtp");
         vtkPolygonWrite(ofs, not_lines);
      }
      {
         std::ofstream ofs("./output/cube" + std::to_string(i) + ".vtp");
         vtkPolygonWrite(ofs, cube.getFaces());
         ofs.close();
      }
   }
};
