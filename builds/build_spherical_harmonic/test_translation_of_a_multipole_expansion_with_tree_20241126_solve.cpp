#include <array>
#include <memory>
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

// ignを正しく計算する．どうやら符号逆．
// python3.11 ../../extract_comments.py README.md -source ./ -include ../../

/*DOC_EXTRACT 0_2_1_translation_of_a_multipole_expansion

\insert{Multipole_Expansion}

## ツリー構造を使った多重極展開の移動

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion_with_tree_20240818.cpp
make
./test_translation_of_a_multipole_expansion_with_tree_20240818
paraview check_M2L.pvsm
```

*/

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

   auto obj = std::make_unique<Network>(std::string(argv[1]));

   /* -------------------------------------------------------------------------- */

   auto setAlmostSolidAngle = [&](auto& B_poles, auto& points) {
      TimeWatch tw;

      // initialize coeffs
      B_poles.forEachAll([&](Buckets<sp_pole4FMM>* B) {
         B->multipole_expansion.initialize();
         B->local_expansion.initialize();
      });

#pragma omp parallel
      for (auto& p : points)
#pragma omp single nowait
         p->phiphin = {1., 0.};

      std::cout << "極の展開 reuse" << std::endl;
      MultipoleExpansionReuse(B_poles);
      std::cout << Magenta << "Multipole Expansion" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "M2M ..." << colorReset << std::endl;
      M2M(B_poles);
      std::cout << Magenta << "M2M" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "M2L ..." << colorReset << std::endl;
      M2L(B_poles);
      std::cout << Magenta << "M2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "L2L ..." << colorReset << std::endl;
      L2L(B_poles);
      std::cout << Magenta << "L2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "L2P ..." << colorReset << std::endl;

#pragma omp parallel
      for (auto& p : points)
#pragma omp single nowait
      {
         auto [IgPhi_IgnPhin_near, IgPhi_IgnPhin_far] = integrate(B_poles, p->X, 1e-4);
         p->almost_solid_angle = IgPhi_IgnPhin_near.first;
      }
      std::cout << Magenta << "L2P" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
   };

   auto updatePole_M2M_M2L_L2L_integrate = [&](auto& B_poles, auto& points) {
      TimeWatch tw;

#pragma omp parallel
      for (auto& p : points)
#pragma omp single nowait
         p->phiphin = {2., 2.};  //!  これは適当に与えている

      std::cout << "極の展開 reuse" << std::endl;
      MultipoleExpansionReuse(B_poles);
      std::cout << Magenta << "Multipole Expansion" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "M2M ..." << colorReset << std::endl;
      M2M(B_poles);
      std::cout << Magenta << "M2M" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "M2L ..." << colorReset << std::endl;
      M2L(B_poles);
      std::cout << Magenta << "M2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "L2L ..." << colorReset << std::endl;
      L2L(B_poles);
      std::cout << Magenta << "L2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "L2P ..." << colorReset << std::endl;

#pragma omp parallel
      for (auto& p : points)
#pragma omp single nowait
      {
         auto [IgPhi_IgnPhin_near, IgPhi_IgnPhin_far] = integrate(B_poles, p->X, 1e-4);
         p->IgPhi_IgnPhin_near = IgPhi_IgnPhin_near;
         p->IgPhi_IgnPhin_far = IgPhi_IgnPhin_far;
         p->IgPhi_IgnPhin_FMM = IgPhi_IgnPhin_near + IgPhi_IgnPhin_far;
      }
      std::cout << Magenta << "L2P" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
   };

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

   Buckets<sp_pole4FMM> B_poles(obj->scaledBounds(1.1), obj->getScale() / 5.);

   //@ -------------------------------------------------------------------------- */
   //@                           バケットの作成．極の追加                             */
   //@ -------------------------------------------------------------------------- */

   std::cout << "バケットの作成．極の追加" << std::endl;
   int count_faces = 0;
   for (auto& F : obj->getFaces()) {

      count_faces += 1;
      auto [p0, p1, p2] = F->getPoints();

      p0->phiphin = {1., 1.};
      p1->phiphin = {1., 1.};
      p2->phiphin = {1., 1.};

      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      for (const auto& [xi0, xi1, ww] : __array_GW6xGW6__) {
         auto N = ModTriShape<3>(xi0, xi1);
         auto X = Dot(N, X012);
         auto weights = Tdd{1., 1.} * Norm(cross) * ww * (1. - xi0);

         //% このようにポインターを使って極を計算するようにしておけば，ポインタ先の値が変わると，極の値は自動的に変わる．
         //% GMRESでは，Krylov空間の基底ベクトルを変えて，行列ベクトル積を計算する．これを高速化するのが，高速多重極展開である．
         //% 高速多重極展開は，行列ベクトル積を直接計算しないし，行列や基底ベクトルを配列として保存することもしない．
         //% 基底の各値をポインタ先に保存し，M2M, M2L, L2L, L2Pを計算することで，結果的に行列ベクトル積を計算する．
         std::function<Tdd()> get_values = [&]() {
            return Tdd{Dot(N, std::array<double, 3>{p0->phiphin[0], p1->phiphin[0], p2->phiphin[0]}),
                       Dot(N, std::array<double, 3>{p0->phiphin[1], p1->phiphin[1], p2->phiphin[1]})};
         };

         B_poles.add(X, std::make_shared<pole4FMM>(X, weights, F->normal, get_values));  //$ 極の追加
      }
   };

   std::cout << Magenta << "Add poles" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   //@ -------------------------------------------------------------------------- */
   //@                                ツリー構造を生成                               */
   //@ -------------------------------------------------------------------------- */
   //! どのように修正するかは要検討．
   //! 途中で展開が不要になるくらい，ソースが少ない場合，その内部のソースは直接積分で計算する．
   std::cout << "ツリー構造を生成" << std::endl;
   int max_level = 7;
   B_poles.setLevel(0, max_level);
   B_poles.generateTree([](auto bucket) {
      if (bucket->all_stored_objects_vector.empty())
         return false;
      else
         return bucket->all_stored_objects_vector.size() > 1250 && bucket->level < bucket->max_level;
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
   //@                                 極の展開                                    */
   //@ -------------------------------------------------------------------------- */

   std::cout << "極の展開" << std::endl;
   MultipoleExpansion(B_poles);
   std::cout << Magenta << "Multipole Expansion" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   TimeWatch twFMM;

   setAlmostSolidAngle(B_poles, obj->getPoints());

   updatePole_M2M_M2L_L2L_integrate(B_poles, obj->getPoints());

   std::cout << Red << "Total Elapsed time : " << twFMM() << colorReset << std::endl;

   //! -------------------------------------------------------------------------- */
   //!                                   直接積分                                   */
   //! -------------------------------------------------------------------------- */

   if (true) {
      TimeWatch twDirect;
      std::cout << "Direct Integration..." << std::endl;

      // for (auto& p : obj->getPoints())
      //    p->igign = direct_integration(B_poles, p->X);

      auto v = ToVector(obj->getPoints());
#pragma omp parallel for
      for (int i = 0; i < v.size(); ++i) {
         auto& p = v[i];  // Safe as 'p' is now private to each thread
         p->igign = direct_integration(B_poles, p->X);
      }

      std::cout << Magenta << "Direct Integration" << Green << ", Elapsed time : " << twDirect() << colorReset << std::endl;
   }

   tw.reset();

   // % -------------------------------------------------------------------------- */
   // %                                     出力                                    */
   // % -------------------------------------------------------------------------- */
   if (true) {

      //! バケツの可視化のための出力

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
