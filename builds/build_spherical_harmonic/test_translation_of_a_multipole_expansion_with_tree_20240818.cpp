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

## 要素法特有の話

そもそものガウス点がすくなければツリーを伸ばしていけて，直接積分の量を減らせるため，早くなる
直接積分の部分のガウス点は減らしたくない．これを両立することが大事だ．

このような議論はされていない．

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

int main(int argc, char* argv[]) {

   if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
      return 1;
   }

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
   setBoundaryTypes(obj.get(), {});
   setPhiPhinOnFace(obj.get());
   auto size = setNodeFaceIndices(obj.get());
   std::cout << "size = " << size << std::endl;
#endif
   /* ------------------------------------- */

   //@ -------------------------------------------------------------------------- */
   //@                       バケットの作成．極の追加 add                             */
   //@ -------------------------------------------------------------------------- */

   //! 実際に問題を解く場合は，計算の流れにそって，phiOnFaceやphinOnFaceを更新する．
   //! ここでは，適当な値を入れている．
   auto targets_points = obj->getPoints();
   for (const auto& p : targets_points) {
      for (const auto& [f, i] : p->f2Index) {
         p->phiphin[0] = p->phiOnFace[f] = RandomReal({-1, 1});
         p->phiphin[1] = p->phinOnFace[f] = RandomReal({-1, 1});
      }
   }

   /* -------------------------------------------------------------------------- */

   auto bounds = obj->getUniformBounds(0.6);
   Buckets<std::shared_ptr<source4FMM<target4FMM>>> B_poles(bounds, (bounds[0][1] - bounds[0][0]) / 3.);
   std::cout << "バケットの作成．極の追加" << std::endl;
   /* phiOnFaceやphinOnFaceを更新することで，get_valuesの結果は自動的に更新される */
   for (auto& F : obj->getSurfaces()) {
      auto [p0, p1, p2] = F->getPoints();
      auto closest_p_to_origin = p0;
      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      auto key0 = std::get<1>(pf2ID(p0, F));
      auto key1 = std::get<1>(pf2ID(p1, F));
      auto key2 = std::get<1>(pf2ID(p2, F));
      std::array<double*, 2> pair_pointer_phiphin0 = {&p0->phiOnFace.at(key0), &p0->phinOnFace.at(key0)};
      std::array<double*, 2> pair_pointer_phiphin1 = {&p1->phiOnFace.at(key1), &p1->phinOnFace.at(key1)};
      std::array<double*, 2> pair_pointer_phiphin2 = {&p2->phiOnFace.at(key2), &p2->phinOnFace.at(key2)};
      for (const auto& [t0t1, quadrature_weight, shape3, X, cross, J_det] : F->map_Point_LinearIntegrationInfo_vector[0].at(closest_p_to_origin)) {
         // b% [2] ソース点を直接積分に置き換えたあとの，積分

         auto get_weighted_source_densities = [shape3,
                                               W = quadrature_weight * J_det,
                                               pair_pointer_phiphin0,
                                               pair_pointer_phiphin1,
                                               pair_pointer_phiphin2]() -> std::array<double, 2> { return {W * Dot(shape3, Tddd{*pair_pointer_phiphin0[0], *pair_pointer_phiphin1[0], *pair_pointer_phiphin2[0]}),
                                                                                                           W * Dot(shape3, Tddd{*pair_pointer_phiphin0[1], *pair_pointer_phiphin1[1], *pair_pointer_phiphin2[1]})}; };

         auto use_this_soruce_when_set_direct_integration = [p0, p1, p2, shape3, cross, X012,
                                                             W = quadrature_weight * J_det,
                                                             pair_pointer_phiphin0,
                                                             pair_pointer_phiphin1,
                                                             pair_pointer_phiphin2](const target4FMM* origin) -> std::vector<std::tuple<std::array<double*, 2>, std::array<double, 2>>> {
            std::array<double, 3> cross = Cross(X012[1] - X012[0], X012[2] - X012[0]), WGN = {0., 0., 0.}, WGnN = {0., 0., 0.}, N012_geometry, R;
            double norm_cross = Norm(cross);
            double tmp, nr;
            for (const auto& [t0, t1, ww] : __array_GW1xGW1__) {
               N012_geometry = ModTriShape<3>(t0, t1);
               nr = Norm(R = (Dot(N012_geometry, X012) - origin->Xtarget));
               tmp = (ww * (1. - t0) / nr);
               WGN += (norm_cross * tmp) * N012_geometry;
               WGnN += (-Dot(R / nr, cross) * tmp / nr) * N012_geometry;
            }
            if (p0 == origin)
               std::get<0>(WGnN) = 0.;  //! リジッドモードテクニック
            if (p1 == origin)
               std::get<1>(WGnN) = 0.;  //! リジッドモードテクニック
            if (p2 == origin)
               std::get<2>(WGnN) = 0.;  //! リジッドモードテクニック
            return {{pair_pointer_phiphin0, std::array<double, 2>{WGN[0], WGnN[0]}},
                    {pair_pointer_phiphin1, std::array<double, 2>{WGN[1], WGnN[1]}},
                    {pair_pointer_phiphin2, std::array<double, 2>{WGN[2], WGnN[2]}}};
         };

         auto pole = std::make_shared<source4FMM<target4FMM>>(X, F->normal,
                                                              get_weighted_source_densities,
                                                              use_this_soruce_when_set_direct_integration);
         B_poles.add(X, pole);
      }
   }
   std::cout << Magenta << "Add poles" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
   //@ -------------------------------------------------------------------------- */
   //@                                ツリー構造を生成                               */
   //@ -------------------------------------------------------------------------- */
   //
   std::cout << "ツリー構造を生成" << std::endl;
   int max_level = 5;
   B_poles.setLevel(0, max_level);
   B_poles.generateTree([](auto bucket) {
      if (bucket->all_stored_objects_vector.empty())
         return false;
      else
         return bucket->all_stored_objects_vector.size() > 500 && bucket->level < bucket->max_level;
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

   //@ -------------------------------------------------------------------------- */
   //@                                   FMM                                      */
   //@ -------------------------------------------------------------------------- */

   //! 繰り返し計算によって問題が生じないかもチェック
   initializeFMM(B_poles, targets_points);

   // for (auto i = 0; i < 3; ++i)
   {
      updateFMM(B_poles);
      TimeWatch twFMM2;
#pragma omp parallel
      for (auto& p : targets_points)
#pragma omp single nowait
      {
         for (const auto& [f, i] : p->f2Index) {
            // auto [IgPhin_IgnPhi_near, IgPhin_IgnPhi_far] = integrate(B_poles, p->X, eps);
            // std::get<1>(IgPhin_IgnPhi_near) -= p->getSolidAngle() * p->phiOnFace.at(f);  //! ここにphiOnFaceは使わない独立させる
            p->wGPhin_wGnPhi_near = p->integrate();
            p->wGPhin_wGnPhi_far = p->integrateFMM();
            p->wGPhin_wGnPhi_FMM = p->wGPhin_wGnPhi_near + p->wGPhin_wGnPhi_far;
            p->diagonal_coefficient = p->wGPhin_wGnPhi_FMM[0] - p->wGPhin_wGnPhi_FMM[1];
         }
      }

      std::cout << Magenta << "integration including L2P" << Green << ", Elapsed time : " << twFMM2() << colorReset << std::endl;
   }

   //! -------------------------------------------------------------------------- */
   //!                                   直接積分                                   */
   //! -------------------------------------------------------------------------- */

   if (true) {
      TimeWatch twDirect;
      std::cout << "Direct Integration..." << std::endl;
#pragma omp parallel
      for (auto& p : obj->getPoints())
#pragma omp single nowait
      {
         p->igign = p->integrate(B_poles);
      }

      std::cout << Magenta << "Direct Integration" << Green << ", Elapsed time : " << twDirect() << colorReset << std::endl;
   }

   tw.reset();

   // % -------------------------------------------------------------------------- */
   // %                                     出力                                    */
   // % -------------------------------------------------------------------------- */
   {

      //! バケツの可視化のための出力

      std::vector<T8Tddd> cube_level, cube_level_deepest, cube_near;
      std::vector<std::vector<T8Tddd>> cube_M2L;
      std::vector<std::vector<Tddd>> poles_in_bucket;
      std::vector<std::vector<T2Tddd>> line_M2L;
      PVDWriter cubePVD("./output/cubes.pvd");
      PVDWriter lineM2LPVD("./output/line_M2L.pvd");

      std::cout << "output" << std::endl;

      {
         B_poles.forEachAtDeepest([&](Buckets<std::shared_ptr<source4FMM<target4FMM>>>* B) {
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            cube_level_deepest.push_back(t8tddd);
         });
         std::ofstream ofs("./output/cube_level_deepest.vtp");
         vtkPolygonWrite(ofs, cube_level_deepest);
         ofs.close();
      }

      std::cout << "paraview ./output/cube_level_deepest.vtp" << std::endl;

      //! Oはランダムに選ぶのではなく，例としてふさわしい，最も深いレベルの中にある極を選ぶ
      Tddd sample_Origin;
      int sample_max_level = 0;

      for (int i = 0; i < B_poles.max_level; ++i) {
         cube_level.clear();
         B_poles.forEachAtLevel(i, [&](Buckets<std::shared_ptr<source4FMM<target4FMM>>>* B) {
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            cube_level.push_back(t8tddd);
            if (sample_max_level < i && !B->all_stored_objects_vector.empty())
               sample_Origin = B->all_stored_objects_vector[0]->X;
         });
         auto name = "./output/cube_level" + std::to_string(i) + ".vtp";
         std::ofstream ofs(name);
         std::cout << "name = " << name << std::endl;
         vtkPolygonWrite(ofs, cube_level);
         ofs.close();
      }

      std::cout << "paraview ./output/cube_level*.vtp" << std::endl;

      const auto& target_bucket = B_poles.getBucketAtDeepest(sample_Origin);
      //
      cube_M2L.resize(B_poles.max_level + 1);
      line_M2L.resize(B_poles.max_level + 1);
      poles_in_bucket.resize(B_poles.max_level + 1);
      int l = B_poles.getBucketAtDeepest(sample_Origin)->level;
      std::cout << "l = " << l << std::endl;

      for (int i = 0; i <= l; i++) {
         auto A = B_poles.getBucketAtLevel(i, sample_Origin);
         // std::cout << "A->level = " << A->level << std::endl;
         for (auto& B : A->buckets_for_M2L) {
            line_M2L[B->level].push_back(T2Tddd{B->X, A->X});
            for (auto& pole : B->all_stored_objects_vector)
               poles_in_bucket[B->level].push_back(pole->X);
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            std::cout << "B->level = " << B->level << std::endl;
            cube_M2L[B->level].push_back(t8tddd);
         }
      }

      //! 出力
      for (int i = 0; i <= max_level; i++) {
         {
            std::ofstream ofs("./output/cube_M2L" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, cube_M2L[i]);
            cubePVD.push("cube_level" + std::to_string(i) + ".vtp", i);
            ofs.close();
         }
         {
            std::ofstream ofs("./output/line_M2L" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, line_M2L[i]);
            lineM2LPVD.push("line_M2L" + std::to_string(i) + ".vtp", i);
            ofs.close();
         }
         {
            std::ofstream ofs("./output/poles_in_bucket" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, poles_in_bucket[i]);
            ofs.close();
         }
      }

      cubePVD.output();
      lineM2LPVD.output();

      std::cout << "paraview ./output/cube_M2L*.vtp" << std::endl;
      std::cout << "paraview ./output/line_M2L*.vtp" << std::endl;

      {
         std::unordered_map<networkPoint*, double> data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13, data14, data15;
#ifdef _BEM_
         std::unordered_map<networkPoint*, double> data1_BEM, data2_BEM, data3_BEM;
#endif
         double max_abs_ign = 0;
         for (const auto& p : obj->getPoints()) {
            auto d_igign = p->igign - p->wGPhin_wGnPhi_FMM;
            auto d_igign_far = p->igign - p->wGPhin_wGnPhi_far;
            auto d_igign_near = p->igign - p->wGPhin_wGnPhi_near;
            data1[p] = d_igign[0];
            data2[p] = d_igign[1];
            data3[p] = p->igign[0];
            data4[p] = p->igign[1];
            data5[p] = std::abs(d_igign[0] / p->igign[0]);
            data6[p] = std::abs(d_igign[1] / p->igign[1]);
            data7[p] = p->wGPhin_wGnPhi_FMM[0];
            data8[p] = p->wGPhin_wGnPhi_FMM[1];
            data9[p] = p->wGPhin_wGnPhi_near[0];
            data10[p] = p->wGPhin_wGnPhi_near[1];
            data11[p] = p->wGPhin_wGnPhi_far[0];
            data12[p] = p->wGPhin_wGnPhi_far[1];
            data14[p] = d_igign_far[1];
            data15[p] = d_igign_near[1];
            data1_BEM[p] = p->CORNER ? 0 : (p->Neumann ? 1 : (p->Dirichlet ? 2 : 3));
            max_abs_ign = std::max(max_abs_ign, std::abs(p->igign[1]));
         }

         for (const auto& p : obj->getPoints()) {
            auto d_igign = p->igign - p->wGPhin_wGnPhi_FMM;
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
                                                                                                 {"ign_org-far", data14},
                                                                                                 {"ign_org-near", data15},
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
