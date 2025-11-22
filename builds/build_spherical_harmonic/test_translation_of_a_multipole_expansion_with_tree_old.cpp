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

   Buckets<networkFace*> bucket_faces(obj->bounds, obj->getScale() / 7.);
   bucket_faces.add(obj->getFaces(), [](auto F) { return F->X; });
   bucket_faces.setLevel(0, 3);
   bucket_faces.generateTree();
   int l = 0, global_l = 0;

   bucket_faces.forEachAll([&](Buckets<networkFace*>* bucket_of_face) {
      bucket_of_face->multipole_expansion.initialize(bucket_of_face->X);
      bucket_of_face->local_expansion.initialize(bucket_of_face->X);
   });

   /* -------------------------------------------------------------------------- */

   auto& expansion_M2L_from_level1 = bucket_faces.getBucketAtLevel(1, Tddd{0., 0., 0.})->local_expansion;
   std::array<double, 3> O = expansion_M2L_from_level1.X, normal, cross;
   std::array<double, 2> weights;
   Tddd X, R;
   double ig = 0, ign = 0, nr, value;

   // * ある特定のレベルの全バケツにアクセス
   // * 自分の下，または上にある全てのバケツにアクセス
   //  * バケツは，自分にExpCoeffsを保存する

   auto accuracy = [](std::complex<double> a, double b) { return (a.real() / b - 1.) * 100; };
   int level = 1;

   /* ------------------------------- バケツの可視化 ----------------------------- */

   std::vector<T8Tddd> cube_level, cube_level_deepest;

   {
      bucket_faces.forEachAtDeepest([&](Buckets<networkFace*>* bucket_of_face) {
         T8Tddd t8tddd = (CoordinateBounds)(bucket_of_face->bounds);
         cube_level_deepest.push_back(t8tddd);
      });
      std::ofstream ofs("./output/cube_level_deepest.vtp");
      vtkPolygonWrite(ofs, cube_level_deepest);
      ofs.close();
   }

   for (int i = 0; i < 4; i++) {
      cube_level.clear();
      bucket_faces.forEachAtLevel(i, [&](Buckets<networkFace*>* bucket_of_face) {
         T8Tddd t8tddd = (CoordinateBounds)(bucket_of_face->bounds);
         cube_level.push_back(t8tddd);
      });
      std::ofstream ofs("./output/cube_level" + std::to_string(i) + ".vtp");
      vtkPolygonWrite(ofs, cube_level);
      ofs.close();
   }

   /* --------------------------- Multipole Expansion -------------------------- */

   std::vector<T2Tddd> lines_ME;
   std::vector<std::vector<T2Tddd>> lines_M2M, lines_M2L, lines_L2L, lines_shifting_in_order;

   std::cout << "bucket_faces.getAll() = " << bucket_faces.getAll().size() << std::endl;
   int count_faces = 0;
   //
   auto ME = [&](Buckets<networkFace*>* bucket_of_face) {
      if (bucket_of_face == nullptr)
         return;
      std::array<double, 3> X, R;
      std::array<double, 2> weights;
      double value, nr;
      auto q012 = std::array<double, 3>{1, 1., 1.};
      for (const auto& F : bucket_of_face->_data_) {
         count_faces += 1;
         auto X012 = ToX(F->getPoints());
         auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
         auto normal = F->normal;
         for (const auto& [xi0, xi1, ww] : __array_GW5xGW5__) {
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
            //
            lines_ME.push_back({bucket_of_face->X, X});
            lines_shifting_in_order[0].push_back({bucket_of_face->X, X});
         }
      }
   };

   lines_shifting_in_order.push_back({});
   bucket_faces.forEachAtDeepest(ME);
   std::cout << "count_faces = " << count_faces << std::endl;

   // 一番下から一番上までのM2Mできるようにする
   //@ ------------------------- Multipole to Multipole ------------------------- */

   lines_M2M.resize(3);
   bucket_faces.forEachAtLevel({2, 1, 0}, [&](Buckets<networkFace*>* bucket_of_face) {
      for (auto& B : bucket_of_face->getAllBucket()) {
         M2M(B->multipole_expansion, bucket_of_face->multipole_expansion);
         lines_M2M[bucket_of_face->level].push_back({bucket_of_face->X, B->X});
      }
   });

   int j = 1;
   for (auto i : {2, 1, 0}) {
      lines_shifting_in_order.push_back({});
      bucket_faces.forEachAtLevel(i, [&](Buckets<networkFace*>* bucket_of_face) {
         for (auto& B : bucket_of_face->getAllBucket()) {
            lines_shifting_in_order[j].push_back({bucket_of_face->X, B->X});
         }
      });
      j += 1;
   }

   //! 　直接積分
   std::cout << Red << ig << colorReset << std::endl;
   std::cout << Green << ign << colorReset << std::endl;

   // expansion_M2L_from_level1.initialize(expansion_M2L_from_level1.X);

   //% ------------------------------ Dual Tree Traversal ----------------------------- */

   std::cout << "M2L" << std::endl;
   lines_M2L.resize(2);

   std::cout << expansion_M2L_from_level1.coeffs << std::endl;

   bucket_faces.forEachAtLevel({0, 1}, [&](Buckets<networkFace*>* the_bucket) {
      bucket_faces.forEachAtLevel({the_bucket->level}, [&](Buckets<networkFace*>* far_bucket) {
         if (the_bucket == far_bucket)
            return;
         M2L(far_bucket->multipole_expansion, the_bucket->local_expansion);
         lines_M2L[the_bucket->level].push_back({the_bucket->X, far_bucket->X});
      });
   });

   std::cout << expansion_M2L_from_level1.coeffs << std::endl;

   for (auto i : {0, 1}) {
      lines_shifting_in_order.push_back({});
      bucket_faces.forEachAtLevel(j, [&](Buckets<networkFace*>* the_bucket) {
         bucket_faces.forEachAtLevel({the_bucket->level}, [&](Buckets<networkFace*>* far_bucket) {
            M2L(far_bucket->multipole_expansion, the_bucket->local_expansion);
            lines_shifting_in_order[j].push_back({the_bucket->X, far_bucket->X});
         });
      });
      j += 1;
   }

   {
      std::cout << Red << expansion_M2L_from_level1.IG_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level1.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
      std::cout << Green << expansion_M2L_from_level1.IGn_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level1.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;
   }

   /* --------------------------------------------------------------------------------- */

   std::cout << "L2L" << std::endl;
   lines_L2L.resize(2);
   bucket_faces.forEachAtLevel({1, 2}, [&](Buckets<networkFace*>* bucket_of_face) {
      for (auto& B : bucket_of_face->getAllBucket()) {
         L2L(bucket_of_face->local_expansion, B->local_expansion);
         lines_L2L[bucket_of_face->level - 1].push_back({bucket_of_face->X, B->X});
      }
   });

   for (auto i : {1, 2}) {
      lines_shifting_in_order.push_back({});
      bucket_faces.forEachAtLevel(i, [&](Buckets<networkFace*>* bucket_of_face) {
         for (auto& B : bucket_of_face->getAllBucket()) {
            L2L(bucket_of_face->local_expansion, B->local_expansion);
            lines_shifting_in_order[j].push_back({bucket_of_face->X, B->X});
         }
      });
      j += 1;
   }

   /* -------------------------------------------------------------------------- */

   {
      std::ofstream ofs("output/lines_ME.vtp");
      vtkPolygonWrite(ofs, lines_ME);
      ofs.close();
   }
   for (auto i = 0; i < 3; i++) {
      std::ofstream ofs("output/lines_M2M" + std::to_string(i) + ".vtp");
      vtkPolygonWrite(ofs, lines_M2M[i]);
      ofs.close();
   }
   for (auto i = 0; i < 2; i++) {
      std::ofstream ofs("output/lines_M2L" + std::to_string(i) + ".vtp");
      vtkPolygonWrite(ofs, lines_M2L[i]);
      ofs.close();
   }
   for (auto i = 0; i < 2; i++) {
      std::ofstream ofs("output/lines_L2L" + std::to_string(i) + ".vtp");
      vtkPolygonWrite(ofs, lines_L2L[i]);
      ofs.close();
   }

   for (auto i = 0; i < lines_shifting_in_order.size(); i++) {
      std::ofstream ofs("output/lines_shifting_in_order" + std::to_string(i) + ".vtp");
      vtkPolygonWrite(ofs, lines_shifting_in_order[i]);
      ofs.close();
   }

   // 表示させて計算をチェック

   {
      std::cout << Red << expansion_M2L_from_level1.IG_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level1.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
      std::cout << Green << expansion_M2L_from_level1.IGn_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level1.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;
   }

   // auto expansion_M2L_from_level2 = bucket_faces.getBucketAtLevel(2, O)->local_expansion;
   // expansion_M2L_from_level2.initialize(O);
   // //
   // bucket_faces.forEachAtLevel(2, [&](std::shared_ptr<Buckets<networkFace*>> far_bucket) {
   //    M2L(far_bucket->multipole_expansion, expansion_M2L_from_level2);
   // });

   // std::cout << "check L2L" << std::endl;
   // std::cout << Red << expansion_M2L_from_level2.IG_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level2.IG_using_L(O), ig) << "\%" << colorReset << std::endl;
   // std::cout << Green << expansion_M2L_from_level2.IGn_using_L(O) << " accuracy: " << accuracy(expansion_M2L_from_level2.IGn_using_L(O), ign) << "\%" << colorReset << std::endl;
}