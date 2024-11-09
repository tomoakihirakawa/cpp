/*

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion_with_tree_ExaFmm.cpp
make
./test_translation_of_a_multipole_expansion_with_tree_ExaFmm
 paraview check_M2L.pvsm
```

*/

#include <array>
#include <memory>
#include "Network.hpp"
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"
#include "vtkWriter.hpp"

auto obj = std::make_unique<Network>("./bunny.obj");

using sp_pole4FMM = std::shared_ptr<pole4FMM>;

const double scale = 3;

// Determine if two buckets are near or far based on their scaled bounds
bool isFar(Buckets<sp_pole4FMM>* A, Buckets<sp_pole4FMM>* B) { return !isInside(B->X, A->scaledBounds(scale)); }
bool isNear(Buckets<sp_pole4FMM>* A, Buckets<sp_pole4FMM>* B) { return isInside(B->X, A->scaledBounds(scale)); }
bool isFar(const std::shared_ptr<Buckets<sp_pole4FMM>>& A, const std::shared_ptr<Buckets<sp_pole4FMM>>& B) { return !isInside(B->X, A->scaledBounds(scale)); }
bool isNear(const std::shared_ptr<Buckets<sp_pole4FMM>>& A, const std::shared_ptr<Buckets<sp_pole4FMM>>& B) { return isInside(B->X, A->scaledBounds(scale)); }

void checkAndAddBuckets(std::shared_ptr<Buckets<sp_pole4FMM>> A, std::shared_ptr<Buckets<sp_pole4FMM>> B) {
   if (isNear(A, B)) {
      for (auto& A_c : A->getAllBucket()) {
         for (auto& B_c : B->getAllBucket()) {
            if (A_c != B_c) {
               auto size = B_c->all_stored_objects.size();
               if (size == 0)
                  continue;
               else if (isFar(A_c, B_c)) {
                  A_c->buckets_for_M2L.emplace_back(B_c);
               } else {
                  A_c->buckets_near.emplace_back(B_c);
                  checkAndAddBuckets(A_c, B_c);
               }
            }
         }
      }
   }
}

int main() {
   /* -------------------------------------------------------------------------- */
   /*                  Load and center the object at the origin                   */
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
   /*                           Create buckets and add poles                      */
   /* -------------------------------------------------------------------------- */
   std::cout << "Creating buckets and adding poles" << std::endl;
   int count_faces = 0;
   TimeWatch tw;

   Buckets<sp_pole4FMM> B_poles(obj->bounds, obj->getScale() / 5);

   double ig = 0, ign = 0;

   for (auto& F : obj->getFaces()) {
      auto q012 = std::array<double, 3>{1, 1., 1.};
      count_faces += 1;
      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      auto normal = F->normal;
      for (const auto& [xi0, xi1, ww] : __array_GW5xGW5__) {
         auto X = Dot(ModTriShape<3>(xi0, xi1), X012);
         auto value = Dot(ModTriShape<3>(xi0, xi1), q012);
         auto weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
         B_poles.add(X, std::make_shared<pole4FMM>(X, weights, normal));  // Add poles
      }
   }

   /* -------------------------------------------------------------------------- */
   /*                                Generate the tree                            */
   /* -------------------------------------------------------------------------- */
   std::cout << "Generating tree structure" << std::endl;
   int max_level = 5;
   B_poles.setLevel(0, max_level);
   auto subdivide_condition = [](auto bucket) {
      return (bucket->level < 1 || (bucket->all_stored_objects.size() > 0 && bucket->level + 1 <= bucket->max_level));
   };

   B_poles.generateTree(subdivide_condition);

   std::cout << Magenta << "Tree" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   B_poles.forEachAll([&](Buckets<sp_pole4FMM>* B) {
      B->multipole_expansion.initialize(B->X);
      B->local_expansion.initialize(B->X);
   });

   /* ---------------------------------- Perform multipole expansion ---------------------------------- */

   B_poles.forEachAtDeepest([&](Buckets<sp_pole4FMM>* B) {
      B->multipole_expansion.increment(B->all_stored_objects);
   });

   std::cout << Magenta << "Multipole Expansion" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */
   /*                      Store M2L partners for each bucket at level 1         */
   /* -------------------------------------------------------------------------- */

   B_poles.forEachAtLevel({1}, [&](Buckets<sp_pole4FMM>* A) {
      A->buckets_for_M2L.clear();
      B_poles.forEachAtLevel({1}, [&](Buckets<sp_pole4FMM>* B) {
         if (isFar(A, B)) {
            if (!B->all_stored_objects.empty()) {
               A->buckets_for_M2L.emplace_back(B);
            }
         } else {
            for (auto& A_c : A->getAllBucket())
               for (auto& B_c : B->getAllBucket())
                  if (A_c != B_c) {
                     if (!B_c->all_stored_objects.empty())
                        continue;
                     else if (isFar(A_c, B_c)) {
                        A_c->buckets_for_M2L.emplace_back(B_c);
                     } else {
                        A_c->buckets_near.emplace_back(B_c);
                        checkAndAddBuckets(A_c, B_c);
                     }
                  }
         }
      });
   });

   //! -------------------------------------------------------------------------- */
   //!                                   Direct Integration                       */
   //! -------------------------------------------------------------------------- */

   Tddd target_X = ToX(RandomSample(ToVector(obj->getPoints()))[0]);
   const auto& target_bucket = B_poles.getBucketAtDeepest(target_X);
   auto& expansion_M2L_from_level1 = target_bucket->local_expansion;
   std::array<double, 3> O = expansion_M2L_from_level1.X;

   tw();

   for (auto& F : obj->getFaces()) {
      auto q012 = std::array<double, 3>{1, 1., 1.};
      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      auto normal = F->normal;
      for (const auto& [xi0, xi1, ww] : __array_GW5xGW5__) {
         auto X = Dot(ModTriShape<3>(xi0, xi1), X012);
         auto value = Dot(ModTriShape<3>(xi0, xi1), q012);
         auto R = X - O;
         auto weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);
         auto nr = Norm(R);
         ig += weights[0] / nr;
         if (std::ranges::none_of(X012, [&](const auto& x) { return Norm(x - O) < 1e-10; }))
            ign += -weights[1] * Dot(R / (nr * nr * nr), normal);
      }
   }

   auto elapsed_time = tw();
   std::cout << "Elapsed time : " << elapsed_time << std::endl;
   std::cout << "Approximate elapsed time for the whole calculation : " << elapsed_time[0] * obj->getPoints().size() << ", obj->getPoints().size():" << obj->getPoints().size() << std::endl;

   auto accuracy = [](std::complex<double> a, double b) { return (a.real() / b - 1.) * 100; };

   // $ -------------------------------------------------------------------------- */
   // $                          Visualize the buckets                            */
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

   // Perform multipole to multipole (M2M)
   std::cout << Magenta << "M2M ..." << colorReset << std::endl;
   for (int i = max_level - 1; i >= 0; i--) {
      B_poles.forEachAtLevelParallel({i}, [&](Buckets<sp_pole4FMM>* B) {
         B->forEachBuckets([&](Buckets<sp_pole4FMM>* b) {
            M2M(b->multipole_expansion, B->multipole_expansion);
         });
      });
   }
   std::cout << Magenta << "M2M" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   //! Direct integration
   std::cout << Red << ig << colorReset << std::endl;
   std::cout << Green << ign << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */
   /*                         Multipole 2 Local expansion                        */
   /* -------------------------------------------------------------------------- */

   std::cout << Magenta << "M2L ..." << colorReset << std::endl;

   int count_M2L = 0;

   std::vector<int> levels_zero2max;
   for (int i = 0; i <= max_level; i++)
      levels_zero2max.push_back(i);

   B_poles.forEachAtLevel(levels_zero2max, [&](Buckets<sp_pole4FMM>* A) {
      int count_M2L_of_A = 0;
#pragma omp parallel
      for (auto& B : A->buckets_for_M2L)
#pragma omp single nowait
      {
         M2L(A->multipole_expansion, B->local_expansion);
      }
   });

   std::cout << Magenta << "M2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
   std::cout << expansion_M2L_from_level1.coeffs << std::endl;

   /* ----------------------------------- Output ----------------------------------- */

   {
      cube_M2L.resize(max_level + 1);
      line_M2L.resize(max_level + 1);
      poles_in_bucket.resize(max_level + 1);
      int l = B_poles.getBucketAtDeepest(O)->level;
      std::cout << "l = " << l << std::endl;

      for (int i = 0; i <= l; i++) {
         auto A = B_poles.getBucketAtLevel(i, O);
         for (auto& B : A->buckets_for_M2L) {
            line_M2L[B->level].push_back(T2Tddd{B->X, A->X});

            for (auto& pole : B->all_stored_objects)
               poles_in_bucket[B->level].push_back(pole->X);

            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            cube_M2L[B->level].push_back(t8tddd);
         }
      }

      //! Output
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

   /* -------------------------------------------------------------------------- */
   /*                           Local 2 Local expansion                          */
   /* -------------------------------------------------------------------------- */

   std::cout << "L2Lで，どの程度の精度が得られるか" << std::endl;

   B_poles.forEachAtLevelParallel(levels_zero2max, [&](Buckets<sp_pole4FMM>* B) {
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
         IG += pole->weights[0] / nr;
         IGn += -pole->weights[1] * Dot(R / (nr * nr * nr), pole->normal);
      }
   };

   direct_integration(target_bucket);
   for (auto& B : target_bucket->buckets_near)
      direct_integration(B);
   for (auto& B : target_bucket->buckets_for_DI)
      direct_integration(B);

   double ig_at_target = target_bucket->local_expansion.IG_using_L(O).real() + IG;
   double ign_at_target = target_bucket->local_expansion.IGn_using_L(O).real() + IGn;

   std::cout << Red << "Direct Integration:" << ig << ", FMM:" << ig_at_target << ", accuracy: " << accuracy(ig_at_target, ig) << "\%" << colorReset << std::endl;
   std::cout << Green << "Direct Integration:" << ign << ", FMM:" << ign_at_target << ", accuracy: " << accuracy(ign_at_target, ign) << "\%" << colorReset << std::endl;
   std::cout << Red << IG << " accuracy: " << accuracy(IG, ig) << "\%" << colorReset << std::endl;
   std::cout << Green << IGn << " accuracy: " << accuracy(IGn, ign) << "\%" << colorReset << std::endl;

   std::cout << Magenta << "L2P" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   /* -------------------------------------------------------------------------- */
   //@ Overall integration speed

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
               IG += pole->weights[0] / nr;
               IGn += -pole->weights[1] * Dot(R / (nr * nr * nr), pole->normal);
            }
         };

         direct_integration(b);
         for (auto& B : b->buckets_near)
            direct_integration(B);
         for (auto& B : b->buckets_for_DI)
            direct_integration(B);

         IG += b->local_expansion.IG_using_L(O).real();
         IGn += b->local_expansion.IGn_using_L(O).real();
      }

      count_origins++;
   }
   std::cout << "count_origins = " << count_origins << std::endl;
   std::cout << Magenta << "Total" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
}
