#include <array>
#include <memory>
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

// ignを正しく計算する．どうやら符号逆．
// python3.11 ../../extract_comments.py README.md -source ./ -include ../../

/*

\insert{Multipole_Expansion}

## ツリー構造を使った多重極展開の移動

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion_with_tree_20240812.cpp
make
./test_translation_of_a_multipole_expansion_with_tree_20240812
paraview check_M2L.pvsm
```

*/

int num_DO_NOT_EXPAND = 0;

#include "Network.hpp"
#include "vtkWriter.hpp"

// auto obj = std::make_unique<Network>("./bunny.obj");
auto obj = std::make_unique<Network>("./pumpkin.obj");
// auto obj = std::make_unique<Network>("./Armadillo.obj");

// Function to write polygons to a file
void writePolygonToFile(const std::string& name, auto polygon) {
   std::ofstream ofs(name);
   vtkPolygonWrite(ofs, polygon);
   ofs.close();
}

/* -------------------------------------------------------------------------- */

using sp_pole4FMM = std::shared_ptr<pole4FMM>;

// const double scale = 3;
// bool isFar(Buckets<sp_pole4FMM>* A, Buckets<sp_pole4FMM>* B) { return !isInside(B->X, A->scaledBounds(scale)); };
// bool isNear(Buckets<sp_pole4FMM>* A, Buckets<sp_pole4FMM>* B) { return isInside(B->X, A->scaledBounds(scale)); };
// bool isFar(const std::shared_ptr<Buckets<sp_pole4FMM>>& A, const std::shared_ptr<Buckets<sp_pole4FMM>>& B) { return !isInside(B->X, A->scaledBounds(scale)); };
// bool isNear(const std::shared_ptr<Buckets<sp_pole4FMM>>& A, const std::shared_ptr<Buckets<sp_pole4FMM>>& B) { return isInside(B->X, A->scaledBounds(scale)); };

// まだ問題がある．
// bool isFar(const Buckets<sp_pole4FMM>* A, const Buckets<sp_pole4FMM>* B);

// template <typename T>
// void checkAndAddBucketsImpl(T A, T B) {
//    if (isNear(A, B)) {
//       for (auto& A_c : A->getAllBucket()) {
//          if (!A_c->all_stored_objects.empty()) {
//             for (auto& B_c : B->getAllBucket()) {
//                if (!B_c->all_stored_objects.empty()) {
//                   bool isnear = false;
//                   if (A_c != B_c && isNear(A_c, B_c)) {
//                      A_c->buckets_near.emplace_back(B_c);
//                      isnear = true;
//                   }

//                   if (isFar(A_c, B_c) && A_c != B_c) {
//                      if (A_c->all_stored_objects.size() < num_DO_NOT_EXPAND && !isnear)
//                         A_c->buckets_for_DI.emplace_back(B_c);
//                      else
//                         A_c->buckets_for_M2L.emplace_back(B_c);
//                   } else {
//                      checkAndAddBucketsImpl(A_c, B_c);
//                   }
//                }
//             }
//          }
//       }
//    }
// }

// void checkAndAddBuckets(std::shared_ptr<Buckets<sp_pole4FMM>> A, std::shared_ptr<Buckets<sp_pole4FMM>> B) {
//    checkAndAddBucketsImpl(A, B);
// }

// void checkAndAddBuckets(Buckets<sp_pole4FMM>* A, Buckets<sp_pole4FMM>* B) {
//    checkAndAddBucketsImpl(A, B);
// }

/* -------------------------------------------------------------------------- */

// std::array<double, 2> direct_integration(const Buckets<sp_pole4FMM>& b, const Tddd& O) {
//    std::array<double, 2> igign = {0, 0};
//    for (const auto& pole : b.all_stored_objects) {
//       auto R = pole->X - O;
//       auto nr = Norm(R);
//       //$ 直接積分
//       igign[0] += pole->weights[0] / nr;
//       if (nr > 1)
//          igign[1] += -(pole->weights[1] / nr) * Dot((R / nr) / nr, pole->normal);
//    }
//    return igign;
// };

// std::array<double, 2> direct_integration(const std::shared_ptr<Buckets<sp_pole4FMM>>& b, const Tddd& O) {
//    std::array<double, 2> igign = {0, 0};
//    for (const auto& pole : b->all_stored_objects) {
//       auto R = pole->X - O;
//       auto nr = Norm(R);
//       //$ 直接積分
//       igign[0] += pole->weights[0] / nr;
//       if (nr > 1)
//          igign[1] += -(pole->weights[1] / nr) * Dot((R / nr) / nr, pole->normal);
//    }
//    return igign;
// };

// std::array<double, 2> direct_integration(const std::vector<std::shared_ptr<Buckets<sp_pole4FMM>>>& buckets, const Tddd& O) {
//    std::array<double, 2> igign = {0, 0};
//    for (const auto& b : buckets)
//       igign += direct_integration(b, O);
//    return igign;
// };

/* -------------------------------------------------------------------------- */

int main() {

   /* -------------------------------------------------------------------------- */
   /*                  境界面の読み込み．境界面を原点に移動する．                       */
   /* -------------------------------------------------------------------------- */

   std::array<double, 3> center = {0., 0., 0.};  //{5., 5., 5.};
   double num = 0.;
   for (const auto& p : obj->getPoints()) {
      center += p->X;
      num += 1.;
   }
   center /= num;
   obj->translate(-center);
   // obj->scale(10.);

   std::ofstream ofs("./bunny_obj.vtp");
   vtkPolygonWrite(ofs, obj->getFaces());
   ofs.close();

   /* -------------------------------------------------------------------------- */
   /*                           バケットの作成．極の追加                             */
   /* -------------------------------------------------------------------------- */
   std::cout << "バケットの作成．極の追加" << std::endl;
   int count_faces = 0;
   TimeWatch tw;

   Buckets<sp_pole4FMM> B_poles(obj->scaledBounds(1.1), obj->getScale() / 5);
   for (auto& F : obj->getFaces()) {
      auto q012 = std::array<double, 3>{1, 1., 1.};
      count_faces += 1;
      auto [p0, p1, p2] = F->getPoints();
      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      auto normal = F->normal;
      for (const auto& [xi0, xi1, ww] : __array_GW6xGW6__) {
         auto N = ModTriShape<3>(xi0, xi1);
         auto X = Dot(N, X012);
         auto value = Dot(N, q012);
         auto weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
         std::function<Tdd()> get_values = [&]() {
            return Tdd{Dot(N, std::array<double, 3>{1., 1., 1.}),
                       Dot(N, std::array<double, 3>{1., 1., 1.})};
         };
         B_poles.add(X, std::make_shared<pole4FMM>(X, weights, normal, get_values));  //$ 極の追加
      }
   };

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
      B->multipole_expansion.increment_moments(B->all_stored_objects);
   });

   std::cout << Magenta << "Multipole Expansion" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */
   /*                      各レベルの各セルのM2Lの相手を保存する                       */
   /* -------------------------------------------------------------------------- */
   std::cout << "各レベルの各セルのM2Lの相手を保存する" << std::endl;
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

   std::cout << Magenta << "M2L buckets for the bucket at level 1" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   //! -------------------------------------------------------------------------- */
   //!                                   直接積分                                   */
   //! -------------------------------------------------------------------------- */

   // expansion_M2L_from_level1場所が悪い？
   // Tddd O = {-0.0012, -0.00376, 0.03518};
   // Tddd O = {-0.05, -0., 0.};

   std::cout << "Direct Integration..." << std::endl;

   for (auto& p : obj->getPoints())
      p->igign = direct_integration(B_poles, p->X);

   std::cout << Magenta << "Direct Integration" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

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
            B->multipole_expansion.M2M(b->multipole_expansion);
         });
      });
   std::cout << Magenta << "M2M" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

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
         B->local_expansion.M2L(A->multipole_expansion);
   });

   std::cout << Magenta << "M2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */
   /*                           Local 2 Local expansion                          */
   /* -------------------------------------------------------------------------- */
   std::cout << "L2Lで，どの程度の精度が得られるか" << std::endl;
   B_poles.forEachAtLevelParallel(levels_zero2max, [&](Buckets<sp_pole4FMM>* B) {
      // for (auto& b : B->getAllBucket())
      //    L2L(B->local_expansion, b->local_expansion);
      B->forEachBuckets([&](Buckets<sp_pole4FMM>* b) {
         b->local_expansion.L2L(B->local_expansion);
      });
   });

   std::cout << Magenta << "L2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   std::cout << Magenta << "L2P ..." << colorReset << std::endl;
   /* -------------------------------------------------------------------------- */
   //@ 全体の総積分の速度

   tw();
   int count_origins = 0;

   for (auto& p : obj->getPoints()) {
      p->IgPhi_IgnPhin_far.fill(0);
      p->IgPhi_IgnPhin_near.fill(0);
      auto b = B_poles.getBucketAtDeepest(p->X);

      p->IgPhi_IgnPhin_near += direct_integration(b, p->X);
      for (auto& B : b->buckets_near)
         p->IgPhi_IgnPhin_near += direct_integration(B, p->X);

      auto IGIGn = b->local_expansion.L2P(p->X);

      p->IgPhi_IgnPhin_far[0] = IGIGn[0];
      p->IgPhi_IgnPhin_far[1] = IGIGn[1];

      count_origins++;
      if (count_origins % 1000 == 0)
         std::cout << "count_origins = " << count_origins << std::endl;

      p->IgPhi_IgnPhin_FMM = p->IgPhi_IgnPhin_far + p->IgPhi_IgnPhin_near;
   }

   std::cout << "count_origins = " << count_origins << std::endl;
   std::cout << Magenta << "Total" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   //@ -------------------------------------------------------------------------- */
   //@                                     出力                                    */
   //@ -------------------------------------------------------------------------- */

   {
      Tddd O = ToX(RandomSample(ToVector(obj->getPoints()))[0]);
      const auto& target_bucket = B_poles.getBucketAtDeepest(O);
      //
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
         std::unordered_map<networkPoint*, double> data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13;
         double max_abs_ign = 0;
         for (const auto& p : obj->getPoints()) {
            auto d_igign = p->igign - p->IgPhi_IgnPhin_FMM;
            data1[p] = d_igign[0];
            data2[p] = d_igign[1];
            data3[p] = p->igign[0];
            data4[p] = p->igign[1];
            data5[p] = std::abs(d_igign[0] / p->igign[0]);
            data6[p] = std::abs(d_igign[1] / p->igign[1]);
            data7[p] = p->IgPhi_IgnPhin_FMM[0];
            data8[p] = p->IgPhi_IgnPhin_FMM[1];
            data9[p] = p->IgPhi_IgnPhin_near[0];
            data10[p] = p->IgPhi_IgnPhin_near[1];
            data11[p] = p->IgPhi_IgnPhin_far[0];
            data12[p] = p->IgPhi_IgnPhin_far[1];
            max_abs_ign = std::max(max_abs_ign, std::abs(p->igign[1]));
         }

         for (const auto& p : obj->getPoints()) {
            auto d_igign = p->igign - p->IgPhi_IgnPhin_FMM;
            data13[p] = std::abs(d_igign[1] / max_abs_ign);
         }

         std::vector<std::tuple<std::string, std::unordered_map<networkPoint*, double>>> data = {{"diff_ig", data1},
                                                                                                 {"diff_ign", data2},
                                                                                                 {"ig", data3},
                                                                                                 {"ign", data4},
                                                                                                 {"rel_err_ig", data5},
                                                                                                 {"rel_err_ign", data6},
                                                                                                 {"ig_FMM", data7},
                                                                                                 {"ign_FMM", data8},
                                                                                                 {"ig_near", data9},
                                                                                                 {"ign_near", data10},
                                                                                                 {"ig_far", data11},
                                                                                                 {"ign_far", data12},
                                                                                                 {"rel_err_ign_max", data13}};
         std::ofstream ofs("./output/nodes_igign.vtp");
         vtkPolygonWrite(ofs, obj->getPoints(), data);
         ofs.close();
      }

      {
         poles_in_bucket[0].clear();
         for (auto& B : target_bucket->buckets_near)
            for (auto& pole : B->all_stored_objects)
               poles_in_bucket[0].push_back(pole->X);

         std::ofstream ofs("./output/poles_in_bucket_near.vtp");
         vtkPolygonWrite(ofs, poles_in_bucket[0]);
         ofs.close();
      }

      {
         cube_near.clear();
         for (auto& B : target_bucket->buckets_near)
            cube_near.push_back((CoordinateBounds)(B->bounds));
         std::ofstream ofs("./output/cube_near.vtp");
         vtkPolygonWrite(ofs, cube_near);
         ofs.close();
      }
   }
   std::cout << Magenta << "出力" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
}
