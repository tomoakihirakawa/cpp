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

/* -------------------------------------------------------------------------- */

const double scale = 3.;
bool isFar(Buckets<networkFace*>* A, Buckets<networkFace*>* B) { return !isInside(B->X, A->scaledBounds(scale)); };
bool isNear(Buckets<networkFace*>* A, Buckets<networkFace*>* B) { return isInside(B->X, A->scaledBounds(scale)); };
bool isFar(const std::shared_ptr<Buckets<networkFace*>>& A, const std::shared_ptr<Buckets<networkFace*>>& B) { return !isInside(B->X, A->scaledBounds(scale)); };
bool isNear(const std::shared_ptr<Buckets<networkFace*>>& A, const std::shared_ptr<Buckets<networkFace*>>& B) { return isInside(B->X, A->scaledBounds(scale)); };

// まだ問題がある．
void checkAndAddBuckets(std::shared_ptr<Buckets<networkFace*>> A, std::shared_ptr<Buckets<networkFace*>> B) {
   if (isNear(A, B)) {
      for (auto& A_c : A->getAllBucket()) {
         for (auto& B_c : B->getAllBucket()) {

            if (A_c != B_c && isNear(A_c, B_c))
               A_c->buckets_near.emplace(B_c);

            // if (B_c->data1D.empty())
            //    continue;
            if (isFar(A_c, B_c) && A_c != B_c)
               A_c->buckets_for_M2L.emplace(B_c);
            else
               checkAndAddBuckets(A_c, B_c);
         }
      }
   }
};

/* -------------------------------------------------------------------------- */

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
   ofs.close();

   //    obj->makeBucketPoints(obj->getScale() / 3.);
   //    obj->makeBucketFaces(obj->getScale() / 3.);

   Buckets<networkFace*> bucket_faces(obj->bounds, obj->getScale() / 5.);
   bucket_faces.add(obj->getFaces(), [](auto F) { return F->X; });
   int max_level = 3;
   bucket_faces.setLevel(0, max_level);
   bucket_faces.generateTree();
   int l = 0, global_l = 0;

   bucket_faces.forEachAll([&](Buckets<networkFace*>* bucket_of_face) {
      bucket_of_face->multipole_expansion.initialize(bucket_of_face->X);
      bucket_of_face->local_expansion.initialize(bucket_of_face->X);
   });

   /* -------------------------------------------------------------------------- */

   //! store M2L buckets for the bucket at level 1
   bucket_faces.forEachAtLevel({1}, [&](Buckets<networkFace*>* A) {
      A->buckets_for_M2L.clear();
      bucket_faces.forEachAtLevel({1}, [&](Buckets<networkFace*>* B) {
         //$ check if another bucket is inside the bucket. If it is not, add it to the list of buckets for M2L
         if (isFar(A, B))
            A->buckets_for_M2L.emplace(B);
         else
            for (auto& A_c : A->getAllBucket())
               for (auto& B_c : B->getAllBucket()) {
                  if (A_c != B_c && isNear(A_c, B_c))
                     A_c->buckets_near.emplace(B_c);

                  if (A_c != B_c && isFar(A_c, B_c))
                     A_c->buckets_for_M2L.emplace(B_c);
                  else
                     checkAndAddBuckets(A_c, B_c);
               }
      });
   });

   /* -------------------------------------------------------------------------- */

   // expansion_M2L_from_level1場所が悪い？
   // Tddd target_X = {-0.0012, -0.00376, 0.03518};
   Tddd target_X = {0., 0., 0.};
   const auto& target_bucket = bucket_faces.getBucketAtDeepest(target_X);
   auto& expansion_M2L_from_level1 = target_bucket->local_expansion;
   std::array<double, 3> O = expansion_M2L_from_level1.X, normal, cross;
   std::array<double, 2> weights;
   Tddd X, R;
   double ig = 0, ign = 0, nr, value;

   // * ある特定のレベルの全バケツにアクセス
   // * 自分の下，または上にある全てのバケツにアクセス
   // * バケツは，自分にExpCoeffsを保存する

   auto accuracy = [](std::complex<double> a, double b) { return (a.real() / b - 1.) * 100; };

   /* ------------------------------- バケツの可視化 ----------------------------- */

   std::vector<T8Tddd> cube_level, cube_level_deepest, cube_near;
   std::vector<std::vector<T8Tddd>> cube_M2L;
   std::vector<std::vector<T3Tddd>> face_in_bucket;
   std::vector<std::vector<T2Tddd>> line_M2L;

   {
      bucket_faces.forEachAtDeepest([&](Buckets<networkFace*>* bucket_of_face) {
         T8Tddd t8tddd = (CoordinateBounds)(bucket_of_face->bounds);
         cube_level_deepest.push_back(t8tddd);
      });
      std::ofstream ofs("./output/cube_level_deepest.vtp");
      vtkPolygonWrite(ofs, cube_level_deepest);
      ofs.close();
   }

   for (int i = 0; i < max_level; i++) {
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

   std::cout << "bucket_faces.getAll() = " << bucket_faces.getAll().size() << std::endl;
   int count_faces = 0;
   //
   auto ME = [&](Buckets<networkFace*>* bucket_of_face) -> std::array<double, 2> {
      // if (bucket_of_face == nullptr)
      //    return std::array<double, 2>{0., 0.};
      auto q012 = std::array<double, 3>{1, 1., 1.};
      for (const auto& F : bucket_of_face->data1D) {
         count_faces += 1;
         auto X012 = ToX(F->getPoints());
         auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
         auto normal = F->normal;
         for (const auto& [xi0, xi1, ww] : __array_GW5xGW5__) {
            auto X = Dot(ModTriShape<3>(xi0, xi1), X012);
            auto value = Dot(ModTriShape<3>(xi0, xi1), q012);
            auto R = X - O;
            auto weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
            auto nr = Norm(R);
            //$ 直接積分
            ig += weights[0] / nr;
            ign += -weights[1] * Dot(R / (nr * nr * nr), normal);
            //$ 多重極展開
            bucket_of_face->multipole_expansion.increment(X, weights, normal);
         }
      }
      return std::array<double, 2>{ig, ign};
   };

   TimeWatch tw;

   std::cout << Magenta << "ME ..." << colorReset << std::endl;
   bucket_faces.forEachAtDeepest(ME);
   std::cout << "count_faces = " << count_faces << std::endl;
   std::cout << Magenta << "ME" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
   // 一番下から一番上までのM2Mできるようにする
   //@ ------------------------- Multipole to Multipole ------------------------- */

   std::cout << Magenta << "M2M ..." << colorReset << std::endl;
   for (int i = max_level - 1; i >= 0; i--)
      bucket_faces.forEachAtLevel({i}, [&](Buckets<networkFace*>* bucket_of_face) {
         for (auto& B : bucket_of_face->getAllBucket())
            M2M(B->multipole_expansion, bucket_of_face->multipole_expansion);
      });
   std::cout << Magenta << "M2M" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   //! 　直接積分
   std::cout << Red << ig << colorReset << std::endl;
   std::cout << Green << ign << colorReset << std::endl;

   // expansion_M2L_from_level1.initialize(expansion_M2L_from_level1.X);

   /* -------------------------------------------------------------------------- */

   std::cout << Magenta << "M2L ..." << colorReset << std::endl;

   // std::cout << expansion_M2L_from_level1.coeffs << std::endl;
   // for (int i = 0; i < max_level; i++)
   //    bucket_faces.forEachAtLevel({i}, [&](Buckets<networkFace*>* A) {
   //       bucket_faces.forEachAtLevel({A->level}, [&](Buckets<networkFace*>* far_bucket) {
   //          if (A == far_bucket)
   //             return;
   //          M2L(far_bucket->multipole_expansion, A->local_expansion);
   //       });
   //    });

   for (int i = 0; i <= max_level; i++)
      bucket_faces.forEachAtLevel({i}, [&](Buckets<networkFace*>* A) {
#pragma omp parallel
         for (auto& B : A->buckets_for_M2L)
#pragma omp single nowait
         {
            M2L(A->multipole_expansion, B->local_expansion);
         }
      });

   {
      cube_M2L.resize(max_level + 1);
      line_M2L.resize(max_level + 1);
      face_in_bucket.resize(max_level + 1);
      int l = bucket_faces.getBucketAtDeepest(O)->level;
      std::cout << "l = " << l << std::endl;

      for (int i = 0; i <= l; i++) {
         auto A = bucket_faces.getBucketAtLevel(i, O);
         // std::cout << "A->level = " << A->level << std::endl;
         for (auto& B : A->buckets_for_M2L) {

            line_M2L[B->level].push_back(T2Tddd{B->X, A->X});

            for (auto& F : B->data1D)
               face_in_bucket[B->level].push_back((T3Tddd)(*F));

            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            // std::cout << "B->level = " << B->level << std::endl;
            cube_M2L[B->level].push_back(t8tddd);
         }
      }

      //! 出力
      for (int i = 0; i <= max_level; i++) {
         {
            std::ofstream ofs("./output/cube_M2L" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, cube_M2L[i]);
            ofs.close();
         }
         {
            std::ofstream ofs("./output/line_M2L" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, line_M2L[i]);
            ofs.close();
         }
         {
            std::ofstream ofs("./output/face_in_bucket" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, face_in_bucket[i]);
            ofs.close();
         }
      }
   }

   {
      cube_near.clear();
      for (auto& B : target_bucket->buckets_near)
         cube_near.push_back((CoordinateBounds)(B->bounds));
      std::ofstream ofs("./output/cube_near.vtp");
      vtkPolygonWrite(ofs, cube_near);
      ofs.close();
   }

   std::cout << Magenta << "M2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
   std::cout << expansion_M2L_from_level1.coeffs << std::endl;

   /* -------------------------------------------------------------------------- */
   std::cout << "L2Lで，どの程度の精度が得られるか" << std::endl;
   for (int i = 0; i < max_level; i++) {
      bucket_faces.forEachAtLevel({i}, [&](Buckets<networkFace*>* bucket_of_face) {
         for (auto& B : bucket_of_face->getAllBucket())
            L2L(bucket_of_face->local_expansion, B->local_expansion);
      });
      // auto ig_at_target = expansion_M2L_from_level1.IG_using_L(O);
      // auto ign_at_target = expansion_M2L_from_level1.IGn_using_L(O);
      // std::cout << Red << ig_at_target << " accuracy: " << accuracy(ig_at_target, ig) << "\%" << colorReset << std::endl;
      // std::cout << Green << ign_at_target << " accuracy: " << accuracy(ign_at_target, ign) << "\%" << colorReset << std::endl;
      // //
      // std::cout << Red << ig_at_target << " accuracy: " << accuracy(ig_at_target, ig) << "\%" << colorReset << std::endl;
      // std::cout << Green << ign_at_target << " accuracy: " << accuracy(ign_at_target, ign) << "\%" << colorReset << std::endl;
   }
   std::cout << Magenta << "L2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */

   std::cout << Magenta << "L2P ..." << colorReset << std::endl;

   double IG = 0;
   double IGn = 0;

   auto accum = [&](const auto& b) {
      for (const auto& F : b->all_stored_objects) {
         auto q012 = std::array<double, 3>{1, 1., 1.};
         count_faces += 1;
         auto X012 = ToX(F->getPoints());
         auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
         auto normal = F->normal;
         for (const auto& [xi0, xi1, ww] : __array_GW5xGW5__) {
            auto X = Dot(ModTriShape<3>(xi0, xi1), X012);
            auto value = Dot(ModTriShape<3>(xi0, xi1), q012);
            auto R = X - O;
            auto weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
            auto nr = Norm(R);
            //$ 直接積分
            IG += weights[0] / nr;
            IGn += -weights[1] * Dot(R / (nr * nr * nr), normal);
         }
      }
      auto ig_at_target = expansion_M2L_from_level1.IG_using_L(O) + IG;
      auto ign_at_target = expansion_M2L_from_level1.IGn_using_L(O) + IGn;
      std::cout << Red << ig_at_target << " accuracy: " << accuracy(ig_at_target, ig) << "\%" << colorReset << std::endl;
      std::cout << Green << ign_at_target << " accuracy: " << accuracy(ign_at_target, ign) << "\%" << colorReset << std::endl;
      //
      std::cout << Red << ig_at_target << " accuracy: " << accuracy(ig_at_target, ig) << "\%" << colorReset << std::endl;
      std::cout << Green << ign_at_target << " accuracy: " << accuracy(ign_at_target, ign) << "\%" << colorReset << std::endl;
      //
      std::cout << red << IG << " accuracy: " << accuracy(IG, ig) << "\%" << colorReset << std::endl;
      std::cout << green << IGn << " accuracy: " << accuracy(IGn, ign) << "\%" << colorReset << std::endl;
   };

   accum(target_bucket);
   for (auto& B : target_bucket->buckets_near)
      accum(B);

   std::cout << Magenta << "L2P" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
}