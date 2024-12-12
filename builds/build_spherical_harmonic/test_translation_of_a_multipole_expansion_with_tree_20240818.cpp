#include <array>
#include <memory>
#include "/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_bem/BEM_setBoundaryTypes.hpp"
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
./test_translation_of_a_multipole_expansion_with_tree_20240818 ./pumpkin.obj
paraview check_M2L.pvsm
```

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
   setBoundaryTypes(obj.get());
   initializePhiPhinOnFace(obj.get());
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

      for (const auto& [t0t1, ww, shape3, X, cross, norm_cross] : F->map_Point_LinearIntegrationInfo.at(closest_p_to_origin)) {
         auto [xi0, xi1] = t0t1;
         auto weights = Tdd{1., 1.} * norm_cross * ww;

         auto id0 = pf2ID(p0, F);
         auto id1 = pf2ID(p1, F);
         auto id2 = pf2ID(p2, F);

         auto phi0 = &p0->phiOnFace[std::get<1>(id0)];
         auto phi1 = &p1->phiOnFace[std::get<1>(id1)];
         auto phi2 = &p2->phiOnFace[std::get<1>(id2)];
         auto phin0 = &p0->phinOnFace[std::get<1>(id0)];
         auto phin1 = &p1->phinOnFace[std::get<1>(id1)];
         auto phin2 = &p2->phinOnFace[std::get<1>(id2)];

         std::function<Tdd()> get_values = [phi0, phi1, phi2, phin0, phin1, phin2, shape3]() {
            return Tdd{Dot(shape3, std::array<double, 3>{*phi0, *phi1, *phi2}),
                       Dot(shape3, std::array<double, 3>{*phin0, *phin1, *phin2})};
         };
         B_poles.add(X, std::make_shared<pole4FMM>(X, weights, F->normal, get_values));  //$ 極の追加
      }
   };

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

   /*DOC_EXTRACT 0_2_1_translation_of_a_multipole_expansion

   1. 立体角と特異的な計算を含む係数を，積分を使って計算する（リジッドモードテクニック）　ただ，直接解法とは違って，phiの係数行列を完全に抜き出す必要はない．
   2. 極の追加：各面に対して極を追加し，バケットに格納する．
   3. ツリー構造の生成：バケットに格納された極を基にツリー構造を生成する．
   4. 多重極展開：ツリー構造を用いて多重極展開を行う．
   5. 特異的な積分計算を省くために，BIEを使って特異でない部分を使って計算する（FMMを利用）．
   6. 線形連立方程式の右辺bを計算する（FMMを利用）．
   7. GMRESに与える，行列ベクトル積を返す関数を作成する．
   8. GMRESクラスに，AdotV関数，b，first guessを与えて解く．

   */

   std::cout << "極の展開" << std::endl;
   MultipoleExpansion(B_poles);
   std::cout << Magenta << "Multipole Expansion" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   TimeWatch twFMM;

   /* ----------------------------- almost solid angleの計算 ----------------------------- */

   TimeWatch tw_solid_angle;
   std::size_t count = 0;
   for (const auto& p : obj->getPoints()) {
      count += p->f2Index.size();
      //! copy
      p->phiOnFace_copy = p->phiOnFace;
      p->phinOnFace_copy = p->phinOnFace;
      for (const auto& [f, i] : p->f2Index) {
         p->phiOnFace.at(f) = 1.;   //! this is known value to calculate b
         p->phinOnFace.at(f) = 0.;  //! this is known value to calculate b
      }
   }

   MEreuse_M2M_M2L_L2L(B_poles);

#pragma omp parallel
   for (auto& p : obj->getPoints())
#pragma omp single nowait
   {
      auto [IgPhin_IgnPhi_near, IgPhin_IgnPhi_far] = integrate(B_poles, p->X, too_close_distance);
      p->almost_solid_angle = -(IgPhin_IgnPhi_near[1] + IgPhin_IgnPhi_far[1]);
   }

   std::cout << Magenta << "L2P" << Green << ", Elapsed time : " << tw_solid_angle() << colorReset << std::endl;

   for (const auto& p : obj->getPoints()) {
      count += p->f2Index.size();
      p->phiOnFace = p->phiOnFace_copy;
      p->phinOnFace = p->phinOnFace_copy;
   }

   /* ---------------------------------- bの計算 ---------------------------------- */

   auto MatrixVectorProduct = [too_close_distance, &obj, &B_poles](const bool LHS = true) -> V_d {
      /*
         基本とする形　(G,G,G,G).(phin,phin,phin,phin) = (Gn,Gn,Gn,Gn).(phi,phi,phi,phi)
         境界条件に応じて変形　(G,-Gn,G,G).(phin,phi,phin,phin) = (Gn,-G,Gn,Gn).(phi,phin,phi,phi)
      */
      std::size_t count = 0;
      for (const auto& p : obj->getPoints())
         count += p->f2Index.size();
      std::vector<double> V(count, 0.);

      MEreuse_M2M_M2L_L2L(B_poles);
#pragma omp parallel
      for (const auto& p : obj->getPoints())
#pragma omp single nowait
      {
         for (const auto& [f, i] : p->f2Index) {
            /*

            元々は，IgnPhiのIgnの一部に特異的な計算が含まれているが，それを除いている．

            */
            auto [IgPhin_IgnPhi_near, IgPhin_IgnPhi_far] = integrate(B_poles, p->X, too_close_distance);
            std::get<1>(IgPhin_IgnPhi_near) -= p->almost_solid_angle * p->phiOnFace.at(f);
            auto [IgPhin, IgnPhi] = IgPhin_IgnPhi_near + IgPhin_IgnPhi_far;
            if (p->CORNER && isNeumannID_BEM(p, f) /*行の変更*/) {
               if (LHS)
                  V[i] = p->phiOnFace.at(f);  // unknown
               else
                  V[pf2Index(p, nullptr)] = p->phiOnFace.at(nullptr);  // known
            } else {
               if (LHS) {
                  if (isDirichletID_BEM(p, f))
                     V[i] = IgPhin;  // unknown
                  else if (isNeumannID_BEM(p, f))
                     V[i] = -IgnPhi;  // unknown
                  else
                     throw std::runtime_error("Error: Boundary type is not defined.");
               } else {
                  if (isDirichletID_BEM(p, f))
                     V[i] = IgnPhi;  // known
                  else if (isNeumannID_BEM(p, f))
                     V[i] = -IgPhin;  // known
                  else
                     throw std::runtime_error("Error: Boundary type is not defined.");
               }
            }
         }
      }
      return V;
   };

   auto return_A_dot_v = [&](const V_d& V) -> V_d {
      //! 値を更新
      for (const auto& p : obj->getPoints()) {
         for (const auto& [f, i] : p->f2Index) {
            if (isDirichletID_BEM(p, f))
               p->phinOnFace.at(f) = V[i];  //! this is unknown value that will be calculated
            else if (isNeumannID_BEM(p, f))
               p->phiOnFace.at(f) = V[i];  //! this is unknown value that will be calculated
            else
               throw std::runtime_error("Error: Boundary type is not defined.");
         }
      }
      return MatrixVectorProduct(true);
   };

   std::vector<double> b = MatrixVectorProduct(false);

   std::cout << Red << "Total Elapsed time : " << twFMM() << colorReset << std::endl;

   /* ------------------------------ GMRES ------------------------------------- */

   // std::cout << "use gmres" << std::endl;
   std::vector<int> list = {3, 3};
   std::vector<double> error;
   std::unordered_map<networkPoint*, double> data_gmres_ans, data_b;
   std::vector<double> x0(count, 1.);
   for (auto gmres_size : list) {
      gmres* GMRES = new gmres(return_A_dot_v, b, x0, gmres_size);
      for (auto& p : obj->getPoints())
         for (const auto& [f, i] : p->f2Index) {
            data_gmres_ans[p] = GMRES->x[i];
            data_b[p] = b[i];
         }
      std::cout << "gmres size = " << gmres_size << std::endl;
      std::cout << "gmres error = " << GMRES->err << std::endl;
      error.push_back(GMRES->err);
      x0 = GMRES->x;
   }

   std::cout << "gmres size list = " << list << std::endl;
   std::cout << "gmres error = " << error << std::endl;

   /* ----------------------------------- 比較 ----------------------------------- */

   for (auto& p : obj->getPoints())
      p->phiphin = {1., 1.};

   MEreuse_M2M_M2L_L2L(B_poles);

   TimeWatch twFMM2;
#pragma omp parallel
   for (auto& p : obj->getPoints())
#pragma omp single nowait
   {
      auto [IgPhin_IgnPhi_near, IgPhin_IgnPhi_far] = integrate(B_poles, p->X, too_close_distance);
      p->IgPhi_IgnPhin_near = IgPhin_IgnPhi_near;
      std::get<1>(p->IgPhi_IgnPhin_near) -= p->almost_solid_angle * p->phiphin[0];
      p->IgPhi_IgnPhin_far = IgPhin_IgnPhi_far;
      p->IgPhi_IgnPhin_FMM = p->IgPhi_IgnPhin_near + p->IgPhi_IgnPhin_far;
   }
   std::cout << Magenta << "L2P" << Green << ", Elapsed time : " << twFMM2() << colorReset << std::endl;

   //! -------------------------------------------------------------------------- */
   //!                                   直接積分                                   */
   //! -------------------------------------------------------------------------- */

   if (true) {
      TimeWatch twDirect;
      std::cout << "Direct Integration..." << std::endl;

      auto v = ToVector(obj->getPoints());
#pragma omp parallel for
      for (int i = 0; i < v.size(); ++i) {
         auto& p = v[i];  // Safe as 'p' is now private to each thread
         p->igign = direct_integration_rigid_mode_technique(&B_poles, p->X, too_close_distance);
         // p->igign = direct_integration(&B_poles, p->X);
         std::get<1>(p->igign) -= p->almost_solid_angle * p->phiphin[0];
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
                                                                                                 {"boundary_type", data1_BEM},
                                                                                                 {"gmres_ans", data_gmres_ans},
                                                                                                 {"b", data_b}};

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
