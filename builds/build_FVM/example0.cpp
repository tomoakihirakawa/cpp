/*DOC_EXTRACT exampo0

# 有限体積法

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0.cpp
make
./example0 pq
```

*/

#include "tetgen1.6.0/tetgen.h"
//
#include "Network.hpp"
#include "vtkWriter.hpp"

using Tddd = std::array<double, 3>;
using DataVariant = std::variant<double, Tddd>;
using DataMap = std::unordered_map<networkPoint*, DataVariant>;

int main(const int argc, const char** argv) {

   std::string obj_dir = "../../obj/FVM_sample/";
   auto container = new Network(obj_dir + "container.obj", "container");
   auto inlet = new Network(obj_dir + "inlet.obj", "inlet");
   auto outlet = new Network(obj_dir + "outlet.obj", "outlet");
   auto water = new Network(obj_dir + "water.obj", "water");

   /* -------------------------------------------------------------------- */
   /*                           四面体の生成と出力                            */
   /* -------------------------------------------------------------------- */

   if (argc > 1) {
      auto command = std::string(argv[1]);
      for (auto net : {container, inlet, outlet, water}) {
         //! 四面体の生成
         net->tetrahedralize(command);
         net->makeBuckets();
      }
   } else
      throw std::runtime_error("Please input the command for tetrahedralization.");

   water->setContactFaces({container, inlet, outlet});

   /* --------------------------------------------------------------------- */

   //! 初期条件の設定
   for (const auto& p : water->getPoints())
      for (const auto f : p->getContactFaces()) {
         if (f->getNetwork() == inlet)
            p->U_FVM = 10. * f->normal;
         if (f->getNetwork() == outlet)
            p->U_FVM = -10. * f->normal;
         if (f->getNetwork() == container)
            p->U_FVM = Chop(p->U_FVM, f->normal);
      }

   //! 移流項の計算
   for (const auto& tet : water->getTetras()) {
      //$ 方程式は，各節点に対して立てる
   }

   /* --------------------------------------------------------------------- */
   /* ------------------------- 接触の確認の為の出力 ------------------------- */

   auto data1 = std::unordered_map<networkPoint*, double>();
   auto data2 = std::unordered_map<networkPoint*, Tddd>();
   for (const auto& p : water->getPoints())
      for (const auto f : p->getContactFaces()) {
         if (data1.find(p) == data1.end())
            data1[p] = 0.;
         if (data2.find(p) == data2.end())
            data2[p] = {0., 0., 0.};
         auto v = data1[p];
         if (f->getNetwork() == inlet)
            v -= 3.;
         if (f->getNetwork() == outlet)
            v += 3.;
         if (f->getNetwork() == container)
            v += 1.;
         data1[p] = v;
         data2[p] = p->U_FVM;
      }

   for (auto net : {container, inlet, outlet, water}) {
      auto name = "./output/" + net->getName() + ".vtu";
      std::ofstream ofs(name);
      vtkUnstructuredGridWrite(ofs, net->getTetras(), {{"data1", data1}, {"data2", data2}});
      ofs.close();
      std::cout << "paraview " << name << std::endl;
   }
}