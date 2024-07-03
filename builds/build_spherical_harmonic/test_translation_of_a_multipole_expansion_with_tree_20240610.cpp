#include <array>
#include <memory>
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

// python3.11 ../../extract_comments.py README.md -source ./ -include ../../

/*DOC_EXTRACT 0_2_1_translation_of_a_multipole_expansion

## ツリー構造を使った多重極展開の移動

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion_with_tree_20240610.cpp
make
./test_translation_of_a_multipole_expansion_with_tree_20240610
 paraview check_M2L.pvsm
```

*/

int num_DO_NOT_EXPAND = 0;

#include "Network.hpp"
#include "vtkWriter.hpp"

auto obj = std::make_unique<Network>("./bunny.obj");
// auto obj = std::make_unique<Network>("./Armadillo.obj");

// Function to write polygons to a file
void writePolygonToFile(const std::string& name, auto polygon) {
   std::ofstream ofs(name);
   vtkPolygonWrite(ofs, polygon);
   ofs.close();
}

/* -------------------------------------------------------------------------- */

using sp_pole4FMM = std::shared_ptr<pole4FMM>;

const double scale = 3;
bool isFar(Buckets<sp_pole4FMM>* A, Buckets<sp_pole4FMM>* B) { return !isInside(B->X, A->scaledBounds(scale)); };
bool isNear(Buckets<sp_pole4FMM>* A, Buckets<sp_pole4FMM>* B) { return isInside(B->X, A->scaledBounds(scale)); };
bool isFar(const std::shared_ptr<Buckets<sp_pole4FMM>>& A, const std::shared_ptr<Buckets<sp_pole4FMM>>& B) { return !isInside(B->X, A->scaledBounds(scale)); };
bool isNear(const std::shared_ptr<Buckets<sp_pole4FMM>>& A, const std::shared_ptr<Buckets<sp_pole4FMM>>& B) { return isInside(B->X, A->scaledBounds(scale)); };

// まだ問題がある．
bool isFar(const Buckets<sp_pole4FMM>* A, const Buckets<sp_pole4FMM>* B);

template <typename T>
void checkAndAddBucketsImpl(T A, T B) {
   if (isNear(A, B)) {
      for (auto& A_c : A->getAllBucket()) {
         if (!A_c->all_stored_objects.empty()) {
            for (auto& B_c : B->getAllBucket()) {
               if (!B_c->all_stored_objects.empty()) {
                  bool isnear = false;
                  if (A_c != B_c && isNear(A_c, B_c)) {
                     A_c->buckets_near.emplace_back(B_c);
                     isnear = true;
                  }

                  if (isFar(A_c, B_c) && A_c != B_c) {
                     if (A_c->all_stored_objects.size() < num_DO_NOT_EXPAND && !isnear)
                        A_c->buckets_for_DI.emplace_back(B_c);
                     else
                        A_c->buckets_for_M2L.emplace_back(B_c);
                  } else {
                     checkAndAddBucketsImpl(A_c, B_c);
                  }
               }
            }
         }
      }
   }
}

void checkAndAddBuckets(std::shared_ptr<Buckets<sp_pole4FMM>> A, std::shared_ptr<Buckets<sp_pole4FMM>> B) {
   checkAndAddBucketsImpl(A, B);
}

void checkAndAddBuckets(Buckets<sp_pole4FMM>* A, Buckets<sp_pole4FMM>* B) {
   checkAndAddBucketsImpl(A, B);
}

/* -------------------------------------------------------------------------- */

int main() {

   /* -------------------------------------------------------------------------- */
   /*                  境界面の読み込み．境界面を原点に移動する．                       */
   /* -------------------------------------------------------------------------- */

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

   /* -------------------------------------------------------------------------- */
   /*                           バケットの作成．極の追加                             */
   /* -------------------------------------------------------------------------- */
   std::cout << "バケットの作成．極の追加" << std::endl;
   int count_faces = 0;
   TimeWatch tw;

   Buckets<sp_pole4FMM> B_poles(obj->bounds, obj->getScale() / 6);

   double ig = 0, ign = 0;

   for (auto& F : obj->getFaces()) {
      auto q012 = std::array<double, 3>{1, 1., 1.};
      count_faces += 1;
      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      auto normal = F->normal;
      for (const auto& [xi0, xi1, ww] : __array_GW10xGW10__) {
         auto X = Dot(ModTriShape<3>(xi0, xi1), X012);
         auto value = Dot(ModTriShape<3>(xi0, xi1), q012);
         auto weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
         B_poles.add(X, std::make_shared<pole4FMM>(X, weights, normal));  //$ 極の追加
      }
   };

   点に誤差を持たせ，プロットする．

   /* -------------------------------------------------------------------------- */
   /*                                ツリー構造を生成                               */
   /* -------------------------------------------------------------------------- */
   std::cout << "ツリー構造を生成" << std::endl;
   int max_level = 4;
   B_poles.setLevel(0, max_level);
   auto condition = [](auto bucket) {
      // return (bucket->level < 1 || (bucket->all_stored_objects.size() > 3000 && bucket->level + 1 <= bucket->max_level));
      return (bucket->level < 1 || (bucket->all_stored_objects.size() > 0 && bucket->level + 1 <= bucket->max_level));
   };

   // ここを修正する．
   //! どのように修正するかは要検討．
   //! 途中で展開が不要になるくらい，ソースが少ない場合，その内部のソースは直接積分で計算する．

   B_poles.generateTree(condition);

   std::cout << Magenta << "Tree" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
   int l = 0, global_l = 0;

   B_poles.forEachAll([&](Buckets<sp_pole4FMM>* B) {
      B->multipole_expansion.initialize(B->X);
      B->local_expansion.initialize(B->X);
   });

   /* ---------------------------------- 極の展開 ---------------------------------- */

   std::cout << "極の展開" << std::endl;

   B_poles.forEachAtDeepestParallel([&](Buckets<sp_pole4FMM>* B) {
      B->multipole_expansion.increment(B->all_stored_objects);
   });

   std::cout << Magenta << "Multipole Expansion" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */
   /*                      各レベルの各セルのM2Lの相手を保存する                       */
   /* -------------------------------------------------------------------------- */

   //! store M2L buckets for the bucket at level 1
   //! ここは，Aの展開係数を，M2Lすべきバケツに保存する．Aが空なら，M2LすべきバケツはAにとってないことになる．
   B_poles.forEachAtLevel({1}, [&](Buckets<sp_pole4FMM>* A) {
      A->buckets_for_M2L.clear();
      //! (1) この設定では，Mする場所が，ソース点があるバケツに限られる．これは，常に妥当な設定である．
      if (!A->all_stored_objects.empty())
         B_poles.forEachAtLevel({1}, [&](Buckets<sp_pole4FMM>* B) {  //$ check if another bucket is inside the bucket. If it is not, add it to the list of buckets for M2L
            if (!B->all_stored_objects.empty()) {
               //! (2) この設定は，Lできる場所が，節点などに限られる．しかも，nearにすら含めない．．．．

               if (A != B && isNear(A, B))
                  A->buckets_near.emplace_back(B);

               if (isFar(A, B))
                  A->buckets_for_M2L.emplace_back(B);
               else
                  checkAndAddBuckets(A, B);  // % 次のレベルに移動する．
            }
         });
   });

   //! -------------------------------------------------------------------------- */
   //!                                   直接積分                                   */
   //! -------------------------------------------------------------------------- */

   // expansion_M2L_from_level1場所が悪い？
   // Tddd O = {-0.0012, -0.00376, 0.03518};
   // Tddd O = {-0.05, -0., 0.};

   Tddd O = ToX(RandomSample(ToVector(obj->getPoints()))[0]);
   const auto& target_bucket = B_poles.getBucketAtDeepest(O);

   tw();

   for (auto& F : obj->getFaces()) {
      auto q012 = std::array<double, 3>{1, 1., 1.};
      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      auto normal = F->normal;
      for (const auto& [xi0, xi1, ww] : __array_GW10xGW10__) {
         auto X = Dot(ModTriShape<3>(xi0, xi1), X012);
         auto value = Dot(ModTriShape<3>(xi0, xi1), q012);
         auto R = X - O;
         auto weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
         auto nr = Norm(R);
         //$ 直接積分
         ig += weights[0] / nr;
         if (std::ranges::none_of(X012, [&](const auto& x) { return Norm(x - O) < 1e-10; }))
            ign += -weights[1] * Dot(R / (nr * nr * nr), normal);
      }
   };

   auto elapsed_time = tw();
   std::cout << "Elasped time : " << elapsed_time << std::endl;
   std::cout << "Approximate elapse time for the whole calculation : " << elapsed_time[0] * obj->getPoints().size() << ", obj->getPoints().size():" << obj->getPoints().size() << std::endl;

   // * ある特定のレベルの全バケツにアクセス
   // * 自分の下，または上にある全てのバケツにアクセス
   // * バケツは，自分にExpCoeffsを保存する

   auto accuracy = [](std::complex<double> a, double b) { return (a.real() / b - 1.) * 100; };

   // $ -------------------------------------------------------------------------- */
   // $                          バケツの可視化のための出力                            */
   // $ -------------------------------------------------------------------------- */

   std::vector<T8Tddd> cube_level, cube_level_deepest, cube_near;
   std::vector<std::vector<T8Tddd>> cube_M2L;
   std::vector<std::vector<Tddd>> poles_in_bucket;
   std::vector<std::vector<T2Tddd>> line_M2L;

   {
      B_poles.forEachAtDeepest([&](Buckets<sp_pole4FMM>* B) {
         T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
         cube_level_deepest.push_back(t8tddd);
      });
      std::ofstream ofs("./output/cube_level_deepest.vtp");
      vtkPolygonWrite(ofs, cube_level_deepest);
      ofs.close();
   }

   for (int i = 0; i < max_level; i++) {
      cube_level.clear();
      B_poles.forEachAtLevel(i, [&](Buckets<sp_pole4FMM>* B) {
         T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
         cube_level.push_back(t8tddd);
      });
      std::ofstream ofs("./output/cube_level" + std::to_string(i) + ".vtp");
      vtkPolygonWrite(ofs, cube_level);
      ofs.close();
   }

   // 一番下から一番上までのM2Mできるようにする
   //@ ------------------------- Multipole to Multipole ------------------------- */

   std::cout << Magenta << "M2M ..." << colorReset << std::endl;
   for (int i = max_level - 1; i >= 0; i--)
      B_poles.forEachAtLevelParallel({i}, [&](Buckets<sp_pole4FMM>* B) {
         // for (auto& b : B->getAllBucket())
         //    M2M(b->multipole_expansion, B->multipole_expansion);
         B->forEachBuckets([&](Buckets<sp_pole4FMM>* b) {
            M2M(b->multipole_expansion, B->multipole_expansion);
         });
      });
   std::cout << Magenta << "M2M" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   //! 直接積分
   std::cout << Red << ig << colorReset << std::endl;
   std::cout << Green << ign << colorReset << std::endl;

   // expansion_M2L_from_level1.initialize(expansion_M2L_from_level1.X);

   /* -------------------------------------------------------------------------- */
   /*                         Multipole 2 Local expansion                        */
   /* -------------------------------------------------------------------------- */

   std::cout << Magenta << "M2L ..." << colorReset << std::endl;

   int count_M2L = 0;

   std::vector<int> levels_zero2max;
   for (int i = 0; i <= max_level; i++)
      levels_zero2max.push_back(i);

   B_poles.forEachAtLevel(levels_zero2max, [&](Buckets<sp_pole4FMM>* A) {
#pragma omp parallel
      for (auto& B : A->buckets_for_M2L)
#pragma omp single nowait
         M2L(A->multipole_expansion, B->local_expansion);
   });

   std::cout << Magenta << "M2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   // 根本的にM2Lを高速化する方法があるはず

   /* ----------------------------------- 出力 ----------------------------------- */

   {
      cube_M2L.resize(max_level + 1);
      line_M2L.resize(max_level + 1);
      poles_in_bucket.resize(max_level + 1);
      int l = B_poles.getBucketAtDeepest(O)->level;
      std::cout << "l = " << l << std::endl;

      for (int i = 0; i <= l; i++) {
         auto A = B_poles.getBucketAtLevel(i, O);
         // std::cout << "A->level = " << A->level << std::endl;
         for (auto& B : A->buckets_for_M2L) {

            line_M2L[B->level].push_back(T2Tddd{B->X, A->X});

            for (auto& pole : B->all_stored_objects)
               poles_in_bucket[B->level].push_back(pole->X);

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
            std::ofstream ofs("./output/poles_in_bucket" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, poles_in_bucket[i]);
            ofs.close();
         }
      }

      //

      {
         poles_in_bucket[0].clear();
         for (auto& B : target_bucket->buckets_near)
            for (auto& pole : B->all_stored_objects)
               poles_in_bucket[0].push_back(pole->X);

         std::ofstream ofs("./output/poles_in_bucket_near.vtp");
         vtkPolygonWrite(ofs, poles_in_bucket[0]);
         ofs.close();
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
   std::cout << Magenta << "出力" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
   /* -------------------------------------------------------------------------- */
   /*                           Local 2 Local expansion                          */
   /* -------------------------------------------------------------------------- */
   std::cout << "L2Lで，どの程度の精度が得られるか" << std::endl;
   B_poles.forEachAtLevelParallel(levels_zero2max, [&](Buckets<sp_pole4FMM>* B) {
      // for (auto& b : B->getAllBucket())
      //    L2L(B->local_expansion, b->local_expansion);
      B->forEachBuckets([&](Buckets<sp_pole4FMM>* b) {
         L2L(B->local_expansion, b->local_expansion);
      });
   });

   std::cout << Magenta << "L2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */

   std::cout << Magenta << "L2P ..." << colorReset << std::endl;

   double IG = 0;
   double IGn = 0;

   auto direct_integration = [&](const auto& b) {
      for (const auto& pole : b->all_stored_objects) {
         auto R = pole->X - O;
         auto nr = Norm(R);
         //$ 直接積分
         IG += pole->weights[0] / nr;
         IGn += -pole->weights[1] * Dot(R / (nr * nr * nr), pole->normal);
      }
   };

   direct_integration(target_bucket);
   for (auto& B : target_bucket->buckets_near)
      direct_integration(B);

   // このバケツだけでなく，このバケツに至るまでの，各層のbuckets_for_DIについて，直接積分を行う必要がある．M2Lを行っていないため．
   for (auto& B : target_bucket->buckets_for_DI)
      direct_integration(B);

   double ig_at_target = target_bucket->local_expansion.IG_using_L(O).real() + IG;
   double ign_at_target = target_bucket->local_expansion.IGn_using_L(O).real() + IGn;
   std::cout << Red << "直接積分:" << ig << ", FMM:" << ig_at_target << ", accuracy: " << accuracy(ig_at_target, ig) << "\%" << colorReset << std::endl;
   std::cout << Green << "直接積分:" << ign << ", FMM:" << ign_at_target << ", accuracy: " << accuracy(ign_at_target, ign) << "\%" << colorReset << std::endl;
   std::cout << Red << IG << " accuracy: " << accuracy(IG, ig) << "\%" << colorReset << std::endl;
   std::cout << Green << IGn << " accuracy: " << accuracy(IGn, ign) << "\%" << colorReset << std::endl;
   std::cout << Magenta << "L2P" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */
   //@ 全体の総積分の速度

   tw();
   int count_origins = 0;
   {
      for (auto& p : obj->getPoints()) {
         IG = 0;
         IGn = 0;
         auto b = B_poles.getBucketAtDeepest(p->X);
         auto O = p->X;

         auto direct_integration = [&](const auto& b) {
            for (const auto& pole : b->all_stored_objects) {
               auto R = pole->X - O;
               auto nr = Norm(R);
               //$ 直接積分
               IG += pole->weights[0] / nr;
               IGn += -pole->weights[1] * Dot(R / (nr * nr * nr), pole->normal);
            }
         };

         direct_integration(b);
         for (auto& B : b->buckets_near)
            direct_integration(B);

         IG += b->local_expansion.IG_using_L(O).real();
         IGn += b->local_expansion.IGn_using_L(O).real();

         count_origins++;
         if (count_origins % 1000 == 0)
            std::cout << "count_origins = " << count_origins << std::endl;
      }
   }
   std::cout << "count_origins = " << count_origins << std::endl;
   std::cout << Magenta << "Total" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
}
