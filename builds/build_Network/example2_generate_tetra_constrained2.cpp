/*DOC_EXTRACT 1_1_tetra


## 四面体の生成（制約付き四面分割 constrained tetrahedralization）

* PLC: piecewise linear complex
* CDT: constrained Delaunay triangulation

CDTの生成法には，主に２つの方法がある\ref{Schewchuk2002}：

* naive gift wrapping algorithm (これはadvancing front algorithmとも呼ばれるものと同じだろう)
* sweep algorithm

[杉原厚吉,計算幾何学](杉原厚吉,計算幾何学)によれば，ドロネー四面体分割以外に，綺麗な四面体分割を作成する方法はほとんど知られていないらしい．
四面体分割は，三角分割の場合のように，最小内角最大性が成り立たたず，スリーバー（sliver）と呼ばれる，外接円が大きくないものの潰れた悪い四面体が作られる可能性がある．
このスリーバーをうまく削除することが重要となる．

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example2_generate_tetra_constrained2.cpp
make
./example2_generate_tetra_constrained2
```

`bunny.obj`のような複雑なポリゴンには，この方法ではうまくいかない．

*/

#include <utility>
#include "Network.hpp"
#include "kernelFunctions.hpp"
#include "minMaxOfFunctions.hpp"
#include "vtkWriter.hpp"

int main(int arg, char **argv) {
   //! read obj file from command line argument
   if (arg < 2)
      throw std::invalid_argument("no input file");
   std::string name = argv[1];
   std::cout << "filename : " << name << std::endl;
   /* -------------------------------------------------------------------------- */
   // 四面体の許容サイズは適切か確認（max_radius_of_sphere）
   Timer timer;
   // std::string name = "./cube2.obj";
   // std::string name = "./water100.obj";
   // std::string name = "./bunny.obj";
   /* -------------------------------------------------------------------------- */
   int n = 10;
   auto object = new Network(name, "object");
   auto candidate = new Network(name, "candidate");
   auto initialFaces = object->getFaces();
   object->makeBucketFaces(0.5 * Mean(extLength(object->getLines())));
   object->BucketFaces.setVector();
   // object->genOctreeOfFaces({1, 5}, 1);
   auto ofs = std::ofstream(_HOME_DIR_ + "/output/polygon.vtu");
   vtkUnstructuredGridWrite(ofs, object->getFaces());
   ofs.close();
   /* -------------------------------------------------------------------------- */
   auto [rangeX, rangeY, rangeZ] = object->getBounds();
   auto vecX = Subdivide(rangeX, n);
   auto vecY = Subdivide(rangeY, n);
   auto vecZ = Subdivide(rangeZ, n);
   auto [x0, x1] = rangeX;
   auto [y0, y1] = rangeY;
   auto [z0, z1] = rangeZ;
   double particle_spacing = Mean(Tddd{vecX[1] - vecX[0], vecY[1] - vecY[0], vecZ[1] - vecZ[0]});
   double depth = particle_spacing;
   /* -------------------------------------------------------------------------- */
   Tdd r = {-0.001, 0.001};
   std::vector<Tddd> XYZ;
   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   int count = 0;
   for (const auto &xyz : XYZ)
      new networkPoint(object, xyz);

   object->setGeometricProperties();
   // object->genOctreeOfPoints({minDepth, maxDepth}, 1);
   //
   {
      auto ofs = std::ofstream(_HOME_DIR_ + "/output/points.vtp");
      vtkPolygonWrite(ofs, object->getPoints());
      ofs.close();
   }
   /* --------------------------------- バケツの準備 --------------------------------- */
   Buckets<networkPoint *> buckets(object->getBounds(), particle_spacing);
   Buckets<networkFace *> face_buckets(object->getBounds(), particle_spacing);
   Buckets<networkTetra *> tetra_buckets(object->getBounds(), particle_spacing);

   Buckets<networkPoint *> candidate_buckets(object->getBounds(), particle_spacing);
   std::unordered_map<networkPoint *, networkFace *> candidate_p2f;
   std::unordered_map<networkFace *, networkPoint *> candidate_f2p;

   std::cout << "timer : " << timer() << std::endl;
   buckets.add(object->getPoints());
   std::cout << "all_stored_objects.size() = " << buckets.all_stored_objects.size() << std::endl;
   std::cout << "timer : " << timer() << std::endl;

   // b* ------------------ ボロノイ図，テトラが問題なくできているか確認 ------------------- */

   auto check = [&]() {
      for (const auto &tet : object->getTetras()) {
         bool error_found = buckets.any_of(tet->circumcenter, 5 * tet->circumradius, /*少し半径を小さくし，ほぼ半径上にあるものはスルーする*/
                                           [&](const auto &p) { return !MemberQ(tet->Points, p) && tet->circumradius - 1E-10 >= Norm(ToX(p) - tet->circumcenter); });
         if (error_found)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error found");
      }
   };

   // b% -------------------------------------------------------------------------- */

   Print("形成された面の内，テトラを持たない面の法線方向において，最も近い点を候補とし，新たにテトラを作成できないか調べる");
   PVDWriter pvd(_HOME_DIR_ + "/output/tetras.pvd");
   std::vector<networkPoint *> points_candidates = ToVector(object->getPoints());

   int last_tet_size = 0;
   int i = 0;
   auto output = [&]() {
      last_tet_size = tetra_buckets.all_stored_objects.size();
      std::cout << "points.size() = " << object->getPoints().size() << std::endl;
      std::cout << " faces.size() = " << object->getFaces().size() << std::endl;
      std::cout << "tetras.size() = " << object->getTetras().size() << std::endl;
      std::cout << Magenta << i << Green << ", Elapsed time : " << timer() << colorReset << std::endl;
      auto filename = _HOME_DIR_ + "/output/tetras" + std::to_string(i) + ".vtu";
      std::ofstream ofs(filename);
      // vtkPolygonWrite(ofs, object->getTetras());
      vtkUnstructuredGridWrite(ofs, object->getTetras());
      pvd.push(filename, count);
      pvd.output();
      ofs.close();
      i++;
   };

   int n_loop = 30;
   auto initial_face = object->getFaces();
   for (auto k = 0; k <= n_loop; k++) {
      auto faces = object->getFaces();
      std::cout << "faces.size() = " << faces.size() << std::endl;
      for (const auto &f : k < 20 ? initial_face : faces) {
         const auto [p0, p1, p2] = f->Points;
         const auto f_center = f->centroid;
         Tddd first_X, empty_direction = {0., 0., 0.};

         /*
         探査方法を決定する．
         */

         Tddd direction_to_search = f->normal;
         auto [t0, t1] = f->Tetras;
         if (t0 != nullptr && t1 != nullptr)
            continue;
         else if (t0 != nullptr)
            direction_to_search = Dot(direction_to_search, t0->incenter - f_center) > 0. ? -direction_to_search : direction_to_search;
         else if (t1 != nullptr)
            direction_to_search = Dot(direction_to_search, t1->incenter - f_center) > 0. ? -direction_to_search : direction_to_search;
         else
            direction_to_search.fill(0.);

         if (!object->isInside_MethodBucket(f->center + 1E-5 * direction_to_search))
            continue;

         /*DOC_EXTRACT 1_1_tetra

         ## スコアリングと選択

         四面体の外接球の中心に点が近いほどスコアは低くなる．

         外接球の半径が小さすぎる場合は四面体の候補から外す．

         */

         std::tuple<networkPoint *, Tddd, std::array<double, 3>> best_p = {nullptr, {0., 0., 0.}, {-1., 0., 0.}};
         bool first = true;
         //! 既にある点を候補とする
         const double too_small_height = 1E-8;
         const double too_far_from_face = 1.;
         const double too_long_edge = 1.;

         for (const auto &p : points_candidates) {

            // if (k >= n_loop - 5) {
            //    int num = (int)(isFace(p, p0, p1) != nullptr) + (int)(isFace(p, p1, p2) != nullptr) + (int)(isFace(p, p2, p0) != nullptr);
            //    if (num == 2) {
            //       auto [is_generated, t] = genTetra(object, Append(f->getPoints(), p));
            //       if (is_generated) {
            //          tetra_buckets.add(t->incenter, t);
            //          break;
            //       }
            //    }
            // }

            //! 明らかに四面体を作成できない点を除外

            if (f->MemberQ(p))
               continue;  //! この点を使って四面体は作成できない

            if (Norm(direction_to_search) != 0.)
               if (Dot(direction_to_search, p->X - f_center) < 0.)
                  continue;  //! この点を使って四面体は作成できない

            if (Norm(p->X - f_center) > too_far_from_face)
               continue;  //! この点を使って四面体は作成できない

            if (std::abs(Dot(p->X - f_center, f->normal)) < too_small_height)
               continue;  //! この点を使って四面体は作成できない

            if (Norm(p->X - p0->X) > too_long_edge || Norm(p->X - p1->X) > too_long_edge || Norm(p->X - p2->X) > too_long_edge)
               continue;  //! この点を使って四面体は作成できない

            /* ----------------------------------- */

            //! 四面体を作成してみて，形状が悪いものは除外

            Tetrahedron testing_Tetra(Append(ToX(f), p->X));

            auto [p0, p1, p2] = f->getPoints();

            const auto tet_circum_c = testing_Tetra.circumcenter;
            const auto tet_circum_r = testing_Tetra.circumradius;
            const auto tet_in_r = testing_Tetra.inradius;

            const double shape_criteria = tet_circum_r / tet_in_r;

            // if (testing_Tetra.volume < 1E-10)
            //    continue;
            // if (0. > tet_in_r || shape_criteria > 1E+5 || tet_in_r < 1E-5 || tet_circum_r < 1E-5 ||
            if (!isFinite(tet_in_r) || !isFinite(tet_circum_r) || !isFinite(1. / tet_in_r) || !isFinite(1. / tet_circum_r))
               continue;  //! この点を使って四面体は作成できない

            // if (!object->isInside_MethodBucket(p->X) || !object->isInside_MethodBucket(testing_Tetra.centroid) || !object->isInside_MethodBucket(testing_Tetra.incenter))
            //    continue;  //! この点を使って四面体は作成できない

            bool intersection_found = false;
            const double e = 1E-10;
            for (const auto &pq : std::vector<std::array<networkPoint *, 2>>{{p, p0}, {p, p1}, {p, p2}, {p0, p1}, {p1, p2}, {p2, p0}}) {
               auto [A, B] = pq;
               for (auto &l : object->getLines()) {
                  auto [C, D] = l->getPoints();
                  if (A != C && A != D && B != C && B != D) {
                     auto s = 0.5, t = 0.5, ds = 0., dt = 0.;
                     auto justToCheckIfIntersects = [&](const double t, const double s) -> Tddd {
                        return C->X + s * (D->X - C->X) - (A->X + t * (B->X - A->X));
                     };
                     auto gradFF = [&](const double t, const double s) -> Tdd {
                        auto F = C->X + s * (D->X - C->X) - (A->X + t * (B->X - A->X));
                        auto dFdt = -(B->X - A->X);
                        auto dFds = (D->X - C->X);
                        return {Dot(F, dFdt), Dot(F, dFds)};
                     };
                     BroydenMethod<std::array<double, 2>> BM(gradFF(t, s), gradFF(t, s + 1E-10));
                     // std::cout << "{A,B,C,D} = " << A->X << ", " << B->X << ", " << C->X << ", " << D->X << std::endl;
                     for (auto i = 0; i < 20; i++) {
                        BM.updateGoodBroyden(gradFF(t, s), gradFF(std::clamp(t - dt, 0., 1.), std::clamp(s - ds, 0., 1.)));
                        t = std::clamp(BM.X[0], 0., 1.);
                        s = std::clamp(BM.X[1], 0., 1.);
                        dt = BM.dX[0];
                        ds = BM.dX[1];
                        // std::cout << "i" << i << ": dt = " << dt << ", ds = " << ds << ", t = " << t << ", s = " << s << std::endl;
                        // std::cin.ignore();
                        // if (Norm(BM.dX) < 1E-12)
                        {
                           // std::cout << "break" << std::endl;
                           if (Norm(justToCheckIfIntersects(t, s)) < 1E-10 && 0. - e <= t && t <= 1. + e && 0. - e <= s && s <= 1. + e) {
                              // std::cout << "is intersected" << justToCheckIfIntersects(t, s) << ", " << Norm(justToCheckIfIntersects(t, s)) << std::endl;
                              intersection_found = true;
                              break;
                           }
                        }
                     }
                     if (intersection_found)
                        break;
                  }
                  if (intersection_found)
                     break;
               }
               if (intersection_found)
                  break;
            }
            if (intersection_found)
               continue;

            std::array<double, 3> scores = {0., 0., shape_criteria};
            double diff;

            buckets.apply(tet_circum_c, tet_circum_r, [&](const auto &other_p) {
               // if ((std::get<0>(scores) < std::get<2>(best_p)[0] && std::get<1>(scores) < std::get<2>(best_p)[1]) || f->MemberQ(other_p))
               //    return;  //! この点を使って四面体は作成できない
               std::get<0>(scores) -= std::clamp(tet_circum_r - Norm(other_p->X - tet_circum_c), 0., tet_circum_r);
            });

            if (isFace(p, p0, p1))
               std::get<1>(scores) += 1.;
            if (isFace(p, p1, p2))
               std::get<1>(scores) += 1.;
            if (isFace(p, p2, p0))
               std::get<1>(scores) += 1.;

            if ((first && std::get<0>(best_p)) ||
                std::get<2>(best_p)[0] < std::get<0>(scores) /*higher score*/
                                                             //  || (std::get<2>(best_p)[0] == std::get<0>(scores) && std::get<2>(best_p)[1] <= std::get<1>(scores))
            ) {
               //  (std::get<2>(best_p)[0] == std::get<0>(scores) && std::get<2>(best_p)[1] == std::get<1>(scores) && std::get<2>(best_p)[2] < shape_criteria)
               /*if equal score compare second score*/
               best_p = {p, p->X, scores};
               first = false;
            }
         }

         // //! 新たに作成する点候補
         // for (auto i = 10; i > 1; --i) {
         //    auto Xtmp = good_X(0.2 * i);
         //    auto sc = score(Xtmp) - 1.;
         //    if (sc > s) {
         //       s = sc;
         //       X = Xtmp;
         //       p = nullptr;
         //    }
         // }

         // print info about p, X, and s
         if (std::get<0>(best_p) == nullptr)
            continue;
         else {
            auto p = std::get<0>(best_p);
            const auto X = std::get<0>(best_p)->X;
            const Tddd center = Centroid(X, p0->X, p1->X, p2->X);
            const T4Tddd tet_t4tddd = {X, p0->X, p1->X, p2->X};
            Tetrahedron tet_org(tet_t4tddd);
            const Tddd mean = tet_org.incenter;
            const double scale = 1E-3;
            const T4Tddd tet = tet_t4tddd - scale * (T4Tddd{X - mean, p0->X - mean, p1->X - mean, p2->X - mean});
            int num = (int)(isFace(p, p0, p1) != nullptr) + (int)(isFace(p, p1, p2) != nullptr) + (int)(isFace(p, p2, p0) != nullptr);
            if (tetra_buckets.none_of(center, Norm(center - p0->X), [&](const auto &other_tet) {
                   return IntersectQ((T4Tddd)(*other_tet), tet);
                })) {
               auto [is_generated, t] = genTetra(object, Append(f->getPoints(), std::get<0>(best_p)));
               if (is_generated)
                  tetra_buckets.add(t->incenter, t);
            } else if (num >= 2) {
               auto [is_generated, t] = genTetra(object, Append(f->getPoints(), std::get<0>(best_p)));
               if (is_generated)
                  tetra_buckets.add(t->incenter, t);
            }
         }

         if (count++ % 5 == 0 && last_tet_size != tetra_buckets.all_stored_objects.size())
            output();
      }
   };

   output();
};

/*DOC_EXTRACT 10_1_vtk

# vtk, vtp, vtu

* VTK (Visualization Toolkit)
   VTKは，3次元データを可視化するためのライブラリでフォーマットという意味ではない．
* VTU (VTK Unstructured Grid Format)
   VTUは，内部構造や体積データの解析の場合に適している．体積のある非構造格子データを扱う際はこれを使う．
* VTP (VTK PolyData Format)
   VTPは，表面のみの表示や表面の特性に焦点を当てる場合に適している．


以下は，どちらも四面体を表現している．

<details>
<summary>VTUフォーマット</summary>

```xml
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
    <UnstructuredGrid>
        <Piece NumberOfPoints="4" NumberOfCells="1">
            <Points>
                <DataArray type="Float32" NumberOfComponents="3" format="ascii">
                    0.157726 -0.00244936 -0.15 0.140393 -0.05 -0.15 0.123855 -0.0239571 -0.15
                    0.162817 -0.05 -0.2
                </DataArray>
            </Points>
            <Cells>
                <DataArray type="Int32" Name="connectivity" format="ascii">
                    0 1 2 3
                </DataArray>
                <DataArray type="Int32" Name="offsets" format="ascii">
                    4
                </DataArray>
                <DataArray type="UInt8" Name="types" format="ascii">
                    10
                </DataArray>
            </Cells>
        </Piece>
    </UnstructuredGrid>
</VTKFile>

```

</details>


<details>
<summary>VTPフォーマット</summary>

```xml
<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
    <PolyData>
        <Piece NumberOfLines="0" NumberOfPoints="4" NumberOfPolys="4" NumberOfStrips="0"
            NumberOfVerts="0">
            <Points>
                <DataArray NumberOfComponents="3" format="ascii" type="Float32">
                    0.157726 -0.00244936 -0.15 0.140393 -0.05 -0.15 0.123855 -0.0239571 -0.15
                    0.162817 -0.05 -0.2
                </DataArray>
            </Points>
            <PointData>
            </PointData>
            <CellData Normals="cell_normals" Scalars="cell_scalars">
            </CellData>
            <Polys>
                <DataArray Name="connectivity" format="ascii" type="Int32">
                    2 3 1 0 3 2 0 1 3 0 2 1
                </DataArray>
                <DataArray Name="offsets" format="ascii" type="Int32">
                    3 6 9 12
                </DataArray>
            </Polys>
            <Lines>
                <DataArray Name="connectivity" format="ascii" type="Int32">

                </DataArray>
                <DataArray Name="offsets" format="ascii" type="Int32">

                </DataArray>
            </Lines>
        </Piece>
    </PolyData>
</VTKFile>
```

</details>

*/