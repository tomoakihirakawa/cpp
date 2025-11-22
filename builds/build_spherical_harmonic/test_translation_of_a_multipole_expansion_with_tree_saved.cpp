#include <array>
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

// python3.11 ../../extract_comments.py README.md -source ./ -include ../../

/*

## ツリー構造を使った多重極展開の移動

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion_with_tree.cpp
make
./test_translation_of_a_multipole_expansion_with_tree
```

*/

#include "Network.hpp"
#include "vtkWriter.hpp"

auto obj = std::make_unique<Network>("./bunny.obj");

// Function to write polygons to a file
void writePolygonToFile(const std::string& name, auto polygon) {
   std::ofstream ofs(name);
   vtkPolygonWrite(ofs, polygon);
   ofs.close();
}
test_translation_of_a_multipole_expansion_with_tree_saved.cpp
int main() {

   std::array<double, 3> center = {0., 0., 0.};
   double num = 0.;
   for (const auto& p : obj->getPoints()) {
      center += p->X;
      num += 1.;
   }
   center /= num;
   obj->translate(-center);

   std::ofstream ofs("./bunny_obj.vtp");
   vtkPolygonWrite(ofs, obj->getFaces());

   //    obj->makeBucketPoints(obj->getScale() / 3.);
   //    obj->makeBucketFaces(obj->getScale() / 3.);

   Buckets<networkFace*> bucket_faces(obj->bounds, obj->getScale() / 15.);

   bucket_faces.add(obj->getFaces(), [](auto F) { return F->X; });
   bucket_faces.setLevel(0, 2);
   bucket_faces.generateTree();
   int l = 0, global_l = 0;

   //    int i = 0;
   //    for (const auto& B : bucket_faces.getAllBucket()) {
   //       std::ofstream ofs("./bunny_obj" + std::to_string(i++) + ".vtp");
   //       std::unordered_set<networkFace*> tmp;
   //       for (const auto& F : B->_data_)
   //          tmp.insert(F);
   //       vtkPolygonWrite(ofs, tmp);
   //       ofs.close();
   //    }

   /* -------------------------------------------------------------------------- */

   std::array<double, 3> O = Tddd{0., 0., 0.};
   //
   std::array<double, 2> weights;
   std::array<double, 3> normal, cross;
   Tddd X, R;
   double ig = 0, ign = 0, nr, value;

   // * ある特定のレベルの全バケツにアクセス
   // * 自分の下，または上にある全てのバケツにアクセス
   //  * バケツは，自分にExpCoeffsを保存する
   auto accuracy = [](std::complex<double> a, double b) { return (a.real() / b - 1.) * 100; };

   /* -------------------------------------------------------------------------- */

   int level = 1;

   /* --------------------------- Multipole Expansion -------------------------- */

   std::cout << "bucket_faces.getAll() = " << bucket_faces.getAll().size() << std::endl;
   int count_faces = 0;
   //
   auto ME = [&](std::shared_ptr<Buckets<networkFace*>> bucket_of_face) {
      std::array<double, 3> X, R;
      std::array<double, 2> weights;
      double value, nr;
      auto q012 = std::array<double, 3>{1, 1., 1.};
      for (const auto& F : bucket_of_face->_data_) {
         count_faces += 1;
         auto X012 = ToX(F->getPoints());
         auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
         auto normal = F->normal;
         for (const auto& [xi0, xi1, ww] : __array_GW10xGW10__) {
            X = Dot(ModTriShape<3>(xi0, xi1), X012);
            value = Dot(ModTriShape<3>(xi0, xi1), q012);
            R = X - O;
            weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
            nr = Norm(R);
            //$ 直接積分
            ig += weights[0] / nr;
            ign += -weights[1] * Dot(R / (nr * nr * nr), normal);
            //$ 多重極展開
            bucket_of_face->multipole_expansion.increment(X, weights, normal);
         }
      }
   };
   bucket_faces.forEachAtLevelParallel(2, ME);
   std::cout << "count_faces = " << count_faces << std::endl;

   /* ------------------------- Multipole to Multipole ------------------------- */
   bucket_faces.forEachAtLevel(1, [&](std::shared_ptr<Buckets<networkFace*>> bucket_of_face) {
      bucket_of_face->multipole_expansion.initialize(bucket_of_face->X);
      for (auto B : bucket_of_face->getAllBucket())
         M2M(B->multipole_expansion, bucket_of_face->multipole_expansion);
   });

   /* ------------------------- Multipole to Local Expansion ------------------------- */

   //! 　直接積分
   std::cout << Red << ig << colorReset << std::endl;
   std::cout << Green << ign << colorReset << std::endl;

   {
      auto expansion_M2L_from_level1 = bucket_faces.getBucketAtLevel(1, O)->local_expansion;
      expansion_M2L_from_level1.initialize(O);
      //
      bucket_faces.forEachAtLevel(1, [&](std::shared_ptr<Buckets<networkFace*>> far_bucket) {
         M2L(far_bucket->multipole_expansion, expansion_M2L_from_level1);
      });

      std::cout << "check L2L" << std::endl;
      std::cout << Red << expansion_M2L_from_level1.IG_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level1.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
      std::cout << Green << expansion_M2L_from_level1.IGn_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level1.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;

      auto expansion_M2L_from_level2 = bucket_faces.getBucketAtLevel(2, O)->local_expansion;
      expansion_M2L_from_level2.initialize(O);
      //
      bucket_faces.forEachAtLevel(2, [&](std::shared_ptr<Buckets<networkFace*>> far_bucket) {
         M2L(far_bucket->multipole_expansion, expansion_M2L_from_level2);
      });

      std::cout << "check L2L" << std::endl;
      std::cout << Red << expansion_M2L_from_level2.IG_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level2.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
      std::cout << Green << expansion_M2L_from_level2.IGn_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level2.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;
   }

   {
      /* ------------------------- Multipole to Local Expansion ------------------------- */
      ExpCoeffs<10> expansion_M2L_from_level1(O);
      ExpCoeffs<10> expansion_M2L_from_level2(O);

      bucket_faces.forEachAtLevel(1, [&](std::shared_ptr<Buckets<networkFace*>> bucket) {
         M2L(bucket->multipole_expansion, expansion_M2L_from_level1);
      });

      bucket_faces.forEachAtLevel(2, [&](std::shared_ptr<Buckets<networkFace*>> bucket) {
         M2L(bucket->multipole_expansion, expansion_M2L_from_level2);
      });

      std::cout << "check L2L" << std::endl;
      std::cout << Red << expansion_M2L_from_level1.IG_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level1.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
      std::cout << Green << expansion_M2L_from_level1.IGn_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level1.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;

      std::cout << "check L2L" << std::endl;
      std::cout << Red << expansion_M2L_from_level2.IG_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level2.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
      std::cout << Green << expansion_M2L_from_level2.IGn_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level2.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;
   }
}