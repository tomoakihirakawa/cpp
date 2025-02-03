#include <array>
#include <memory>

bool _PSEUDO_QUADRATIC_ELEMENT_ = false;

#include "/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_bem/BEM_setBoundaryTypes.hpp"
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT 0_2_1_translation_of_a_multipole_expansion

\insert{Multipole_Expansion}

## ツリー構造を使った多重極展開の移動

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion_with_tree_20240818.cpp
make
./test_translation_of_a_multipole_expansion_with_tree_20240818 ./pumpkin.obj
paraview check_M2L.pvsm
```

*/

/*DOC_EXTRACT 0_2_1_translation_of_a_multipole_expansion

# Fast Multipole Method

## pole class

pole class has the following attributes:

- position
- weights
- normal vector
- updater function (to update the intensity, that is the potential, of the pole)

## Buckets class

Buckets class stores specified objects as `Buckets<T>`, and generates tree structure until the number of objects in a bucket is less than or equal to the specified number of objects per bucket.

The step to generate the tree structure should be as follows:

1. add objects to the bucket
2. set the maximum level of the tree using `setLevel`
3. generate the tree structure using `generateTree` while specifying the condition to stop the generation of the tree structure


# Fast Multipole Method

The Fast Multipole Method (FMM) is an algorithm for the efficient calculation of the integration of the pole/potential using the tree structure, the multipole expansion, shifting expansion, and the local expansion. Since FMM calculates integration/summation, such as BIE and does not make the coefficient matrix, solver for the simultaneous linear equations should be iterative methods. GMRES is commonly used for the solver with FMM.

| First steps | GRMES iterative step | description | | |
| --- | --- | --- | --- | --- |
| 1 | | add poles to the root bucket | | |
| 2 | | generate the tree structure from the root bucket | | |
| 3 (before M2M) | | expansion of the poles | | |
| 4 | 1 | **update the intensity of the poles** | | |
| 5 | 2 | Multipole to Multipole (M2M): shift the multipole expansion at each center, from the deeper level to the upper level | about 8 🪣 -> 1 parent 🪣 | use pre-computed SPH |
| 6 | 3 |  Multipole to Local (M2L)| every 🪣 -> (only same level) -> many local 🪣 | use pre-computed SPH |
| 7 | 4 | Local to Local (L2L) | 1 🪣 -> about 8 children 🪣 | use pre-computed SPH |
| 8 | 5 | Add direct integration for the near field and the integration using the local expansion for the far field | | |

Many part of process are dependent on relative position of the poles and the buckets. Therefore, many part of the first steps are saved and reused in the following iterative steps. Remaining part for iterative steps are the update of the intensity of the poles, and simple incrementatation in four-fold for-loops. However, the number of incrementation is not negligible, and the direct integration for the near field also takes time. FMM is surely faster than the direct summation when the number of poles is more than about 10000, but the calculation time is already long when the number of poles is about 10000.
*/

#define _BEM_

#include "Network.hpp"
#include "vtkWriter.hpp"

// auto obj = std::make_unique<Network>("./bunny.obj");
// auto obj = std::make_unique<Network>("./pumpkin.obj");
// auto obj = std::make_unique<Network>("./HumanFace.obj");
// auto obj = std::make_unique<Network>("./HumanBrain.obj");
// auto obj = std::make_unique<Network>("./Armadillo.obj");

// Function to write polygons to a file
void writePolygonToFile(const std::string& name, auto polygon) {
   std::ofstream ofs(name);
   vtkPolygonWrite(ofs, polygon);
   ofs.close();
}

/* -------------------------------------------------------------------------- */

bool bool_output = true;

int main(int argc, char* argv[]) {

   if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
      return 1;
   }
   if (argc > 2)
      bool_output = argv[2];

   auto obj = std::make_unique<Network>(std::string(argv[1]));

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

   TimeWatch tw;

   /*

   1. add poles to the bucket
   2. generate tree
   3. multipole expansion

   */

   /* -------------- BEM ----------------- */
#ifdef _BEM_
   obj->setGeometricProperties();
   for (const auto& f : obj->getFaces())
      f->setIntegrationInfo();
   setBoundaryTypes(obj.get());
   setPhiPhinOnFace(obj.get());
   auto size = setNodeFaceIndices(obj.get());
   std::cout << "size = " << size << std::endl;
#endif
   /* ------------------------------------- */

   //@ -------------------------------------------------------------------------- */
   //@                       バケットの作成．極の追加 add                             */
   //@ -------------------------------------------------------------------------- */

   Buckets<sp_pole4FMM> B_poles(obj->scaledBounds(1.1), obj->getScale() / 6.);

   std::cout << "バケットの作成．極の追加" << std::endl;

   /*
      phiOnFaceやphinOnFaceを更新することで，get_valuesの結果は自動的に更新される．
   */

   for (auto& F : obj->getFaces()) {
      auto [p0, p1, p2] = F->getPoints();
      auto closest_p_to_origin = p0;
      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);

      for (const auto& [t0t1, ww, shape3, X, cross, norm_cross] : F->map_Point_LinearIntegrationInfo_vector[0].at(closest_p_to_origin)) {
         auto [xi0, xi1] = t0t1;
         auto weights = Tdd{norm_cross * ww, norm_cross * ww};
         //$ 極の追加
         auto pole = std::make_shared<pole4FMM>(X,
                                                weights,
                                                F->normal,
                                                [p0, p1, p2, F,
                                                 f0 = std::get<1>(pf2ID(p0, F)),
                                                 f1 = std::get<1>(pf2ID(p1, F)),
                                                 f2 = std::get<1>(pf2ID(p2, F)),
                                                 shape3](pole4FMM* self) -> void {
                                                   //! 再度計算する updater
                                                   std::get<0>(self->values) = Dot(shape3, Tddd{p0->meanPhiOnFace(), p1->meanPhiOnFace(), p2->meanPhiOnFace()});
                                                   std::get<1>(self->values) = Dot(shape3, Tddd{p0->phinOnFace.at(f0), p1->phinOnFace.at(f1), p2->phinOnFace.at(f2)});
                                                });
         B_poles.add(X, pole);
         pole->update();
      }
   }

   std::cout << Magenta << "Add poles" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   //@ -------------------------------------------------------------------------- */
   //@                                ツリー構造を生成                               */
   //@ -------------------------------------------------------------------------- */

   std::cout << "ツリー構造を生成" << std::endl;
   int max_level = 10;
   B_poles.setLevel(0, max_level);
   B_poles.generateTree([](auto bucket) {
      if (bucket->all_stored_objects_vector.empty())
         return false;
      else
         return bucket->all_stored_objects_vector.size() > 1000 && bucket->level < bucket->max_level;
   });
   std::cout << Magenta << "Tree" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   // show info of tree
   for (auto i = 0; i < B_poles.level_buckets.size(); ++i) {
      int mean_M2L_size = 0;
      for (auto m2l : B_poles.level_buckets[i])
         mean_M2L_size += m2l->buckets_for_M2L.size();
      mean_M2L_size /= B_poles.level_buckets[i].size();
      std::cout << "level = " << i << ", size = " << B_poles.level_buckets[i].size() << ", mean M2L size = " << mean_M2L_size << std::endl;
   }

   /* -------------------------------------------------------------------------- */

   double area = 0.;
   for (auto& F : obj->getFaces())
      area += F->area;
   const double too_close_distance = 0.1 * std::sqrt(area / obj->getPoints().size());

   //@ -------------------------------------------------------------------------- */
   //@                                  FMM                                      */
   //@ -------------------------------------------------------------------------- */

   //! 実際に問題を解く場合は，計算の流れにそって，phiOnFaceやphinOnFaceを更新する．
   //! ここでは，適当な値を入れている．
   for (const auto& p : obj->getPoints()) {
      for (const auto& [f, i] : p->f2Index) {
         p->phiOnFace.at(f) = 1.;   //! this is known value to calculate b
         p->phinOnFace.at(f) = 1.;  //! this is known value to calculate b
      }
   }

   std::cout << "極の展開" << std::endl;
   MultipoleExpansion(B_poles);
   std::cout << Magenta << "Multipole Expansion" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   TimeWatch twFMM;

   /* -------------------------------------------------------------------------- */

   setM2M(B_poles);
   setM2L(B_poles);
   setL2L(B_poles);

   /* -------------------------------------------------------------------------- */

   std::cout << Red << "Total Elapsed time : " << twFMM() << colorReset << std::endl;

   /* ----------------------------------- 比較 ----------------------------------- */

   for (const auto& p : obj->getPoints()) {
      p->phiOnFace_copy = p->phiOnFace;
      p->phinOnFace_copy = p->phinOnFace;
      for (const auto& [f, i] : p->f2Index) {
         if (isDirichletID_BEM(p, f)) {
            p->phinOnFace.at(f) = 1;
            p->phiOnFace.at(f) = 1;
         } else if (isNeumannID_BEM(p, f)) {
            p->phinOnFace.at(f) = 1;
            p->phiOnFace.at(f) = 1;
         } else
            throw std::runtime_error("Error: Boundary type is not defined.");
      }
   }

   updatePole_ME_M2M_M2L_L2L(B_poles);

   TimeWatch twFMM2;
#pragma omp parallel
   for (auto& p : obj->getPoints())
#pragma omp single nowait
   {
      double A = 0, n = 0, eps = 0;
      for (auto& f : p->getFaces())
         A += f->area;
      eps = std::sqrt(A / M_PI) * 0.01;
      for (const auto& [f, i] : p->f2Index) {
         auto [IgPhin_IgnPhi_near, IgPhin_IgnPhi_far] = integrate(B_poles, p->X, eps);
         std::get<1>(IgPhin_IgnPhi_near) += p->solid_angle * p->phiOnFace.at(f);
         p->IgPhi_IgnPhin_near = IgPhin_IgnPhi_near;
         p->IgPhi_IgnPhin_far = IgPhin_IgnPhi_far;
         p->IgPhi_IgnPhin_FMM = p->IgPhi_IgnPhin_near + p->IgPhi_IgnPhin_far;
      }
   }

   std::cout << Magenta << "integration including L2P" << Green << ", Elapsed time : " << twFMM2() << colorReset << std::endl;

   for (const auto& p : obj->getPoints()) {
      p->phiOnFace = p->phiOnFace_copy;
      p->phinOnFace = p->phinOnFace_copy;
   }

   //! -------------------------------------------------------------------------- */
   //!                                   直接積分                                   */
   //! -------------------------------------------------------------------------- */

   if (true) {
      TimeWatch twDirect;
      std::cout << "Direct Integration..." << std::endl;

      for (const auto& p : obj->getPoints())
         p->phiphin = {std::get<2>(p->X), std::get<2>(p->X)};

      auto v = ToVector(obj->getPoints());
#pragma omp parallel for
      for (int i = 0; i < v.size(); ++i) {
         auto& p = v[i];  // Safe as 'p' is now private to each thread
         double A = 0, n = 0, eps = 0;
         for (auto& f : p->getFaces())
            A += f->area;
         eps = std::sqrt(A / M_PI) * 0.01;
         p->igign = direct_integration_rigid_mode_technique(&B_poles, p->X, eps);
         // p->igign = direct_integration(&B_poles, p->X);
         // std::get<1>(p->igign) -= p->almost_solid_angle * p->phiphin[0];
         std::get<1>(p->igign) -= p->solid_angle * p->phiphin[0];
      }

      std::cout << Magenta << "Direct Integration" << Green << ", Elapsed time : " << twDirect() << colorReset << std::endl;
   }

   tw.reset();

   // % -------------------------------------------------------------------------- */
   // %                                     出力                                    */
   // % -------------------------------------------------------------------------- */
   if (bool_output) {

      //! バケツの可視化のための出力

      std::vector<T8Tddd> cube_level, cube_level_deepest, cube_near;
      std::vector<std::vector<T8Tddd>> cube_M2L;
      std::vector<std::vector<Tddd>> poles_in_bucket;
      std::vector<std::vector<T2Tddd>> line_M2L;

      std::cout << "output" << std::endl;

      {
         B_poles.forEachAtDeepest([&](Buckets<sp_pole4FMM>* B) {
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            cube_level_deepest.push_back(t8tddd);
         });
         std::ofstream ofs("./output/cube_level_deepest.vtp");
         vtkPolygonWrite(ofs, cube_level_deepest);
         ofs.close();
      }

      std::cout << "paraview ./output/cube_level_deepest.vtp" << std::endl;

      for (int i = 0; i < B_poles.max_level; i++) {
         cube_level.clear();
         B_poles.forEachAtLevel(i, [&](Buckets<sp_pole4FMM>* B) {
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            cube_level.push_back(t8tddd);
         });
         std::ofstream ofs("./output/cube_level" + std::to_string(i) + ".vtp");

         vtkPolygonWrite(ofs, cube_level);
         ofs.close();
      }

      std::cout << "paraview ./output/cube_level*.vtp" << std::endl;

      Tddd O = ToX(RandomSample(ToVector(obj->getPoints()))[0]);
      const auto& target_bucket = B_poles.getBucketAtDeepest(O);
      //
      cube_M2L.resize(B_poles.max_level + 1);
      line_M2L.resize(B_poles.max_level + 1);
      poles_in_bucket.resize(B_poles.max_level + 1);
      int l = B_poles.getBucketAtDeepest(O)->level;
      std::cout << "l = " << l << std::endl;

      for (int i = 0; i <= l; i++) {
         auto A = B_poles.getBucketAtLevel(i, O);
         // std::cout << "A->level = " << A->level << std::endl;
         for (auto& B : A->buckets_for_M2L) {

            line_M2L[B->level].push_back(T2Tddd{B->X, A->X});

            for (auto& pole : B->all_stored_objects_vector)
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

      std::cout << "paraview ./output/cube_M2L*.vtp" << std::endl;
      std::cout << "paraview ./output/line_M2L*.vtp" << std::endl;

      {
         std::unordered_map<networkPoint*, double> data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13, data14;
#ifdef _BEM_
         std::unordered_map<networkPoint*, double> data1_BEM, data2_BEM, data3_BEM;
#endif
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
            data14[p] = p->almost_solid_angle;
            data1_BEM[p] = p->CORNER ? 0 : (p->Neumann ? 1 : (p->Dirichlet ? 2 : 3));
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
                                                                                                 {"rel_err_ign_max", data13},
                                                                                                 {"almost_solid_angle", data14},
                                                                                                 {"boundary_type", data1_BEM}};

         std::ofstream ofs("./output/nodes_igign.vtp");
         vtkPolygonWrite(ofs, obj->getPoints(), data);
         ofs.close();
      }

      {
         poles_in_bucket[0].clear();
         for (auto& B : target_bucket->buckets_near)
            for (auto& pole : B->all_stored_objects_vector)
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
