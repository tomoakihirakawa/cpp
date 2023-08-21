
/*DOC_EXTRACT

## 四面体を生成（制約付き四面分割 constrained tetrahedralization）

* PLC: piecewise linear complex
* CDT: constrained Delaunay triangulation

CDTの生成法には，主に２つの方法がある[Schewchuk 2002](Schewchuk 2002)：

* naive gift wrapping algorithm (これはadvancing front algorithmとも呼ばれるものと同じだろう)
* sweep algorithm

### Advancing Front Algorithm


*/

#include <utility>
#include "Network.hpp"
#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

int main(int arg, char **argv) {
   /* -------------------------------------------------------------------------- */
   // 四面体の許容サイズは適切か確認（max_radius_of_sphere）
   Timer timer;
   std::string name = "./bunny.obj";
   /* -------------------------------------------------------------------------- */
   int n = 10;
   auto object = new Network(name, "object");
   auto initialFaces = object->getFaces();
   object->makeBucketFaces(0.5 * Mean(extLength(object->getLines())));
   object->BucketFaces.setVector();
   // object->genOctreeOfFaces({1, 5}, 1);
   auto ofs = std::ofstream(_HOME_DIR_ + "/output/polygon.vtp");
   vtkPolygonWrite(ofs, object->getFaces());
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
   for (const auto &x : vecX)
      for (const auto &y : vecY)
         for (const auto &z : vecZ) {
            Tddd xyz = {x + RandomReal(r), y + RandomReal(r), z + RandomReal(r)};
            if (object->isInside_MethodBucket(xyz))
               XYZ.push_back(xyz);
         }
   /* -------------------------------------------------------------------------- */
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
   //
   std::cout << "timer : " << timer() << std::endl;
   buckets.add(object->getPoints());
   std::cout << "all_stored_objects.size() = " << buckets.all_stored_objects.size() << std::endl;
   std::cout << "timer : " << timer() << std::endl;
   /* -------------------------------------------------------------------------- */
   // b* ------------------ ボロノイ図，テトラが問題なくできているか確認 ------------------- */
   auto check = [&]() {
      for (const auto &tet : object->getTetras()) {
         bool error_found = buckets.any_of(tet->circumcenter, 5 * tet->circumradius, /*少し半径を小さくし，ほぼ半径上にあるものはスルーする*/
                                           [&](const auto &p) { return !MemberQ(tet->Points, p) && tet->circumradius - 1E-10 >= Norm(ToX(p) - tet->circumcenter); });
         if (error_found)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error found");
      }
   };
   /* -------------------------------------------------------------------------- */
   /*                           最短距離にある２点を繋ぐ                             */
   /* -------------------------------------------------------------------------- */
   Print("最短距離にある２点を繋ぐ", Red);
   // randomに分散させると繋ぐ２点がとても少なくなる
   std::unordered_map<networkPoint *, std::vector<networkPoint *>> p_neighbors;
   for (const auto &a : object->getPoints()) {
      networkPoint *p_closest = nullptr;
      double distance = 1E+10;
      std::vector<networkPoint *> alredy_stored_points;
      buckets.apply(a->X, 2 * depth, [&](const auto &b) {
         if (a != b) {
            bool found_one_in_similar_direction = false;
            for (auto &p : alredy_stored_points) {
               if (Dot(p->X - a->X, b->X - a->X) > 0 /*この方向で確認*/) {
                  found_one_in_similar_direction = true;
                  if (Distance(a, p) >= Distance(a, b)) {
                     distance = Distance(a, b);
                     p = b;
                  }
               }
            }
            if (!found_one_in_similar_direction)
               alredy_stored_points.emplace_back(b);
         }
      });

      auto &nei = p_neighbors[a];
      for (auto &p : alredy_stored_points) {
         if (Distance(a, p) < 3 * particle_spacing) {
            link(a, p);
            nei.emplace_back(p);
         }
      }
   }
   std::cout << Green << "Elapsed time : " << timer() << colorOff << std::endl;
   {
      std::ofstream ofs(_HOME_DIR_ + "/output/first_lines.vtp");
      vtkPolygonWrite(ofs, object->getLines());
   }
   check();

   // b% -------------------------------------------------------------------------- */

   auto accumulate = [&](const T_4P &ABCD, const T_2P &pq) {
      double ACCUM = 0.;
      if (!(std::ranges::any_of(ABCD, [&](const auto &P) { return P == std::get<0>(pq); }) &&
            std::ranges::any_of(ABCD, [&](const auto &P) { return P == std::get<1>(pq); })))
         std::ranges::for_each(ToT4T3P(ABCD), [&](const auto &tri) {
            auto t = Triangle(ToX(tri));
            if (IntersectQ(t.incenter, t.inradius, ToX(pq)))
               ACCUM += std::abs(t.inradius - Norm(t.incenter - Nearest(t.incenter, ToX(pq))));
         });
      return ACCUM;
   };

   /* -------------------------------------------------------------------------- */

   auto criterion = [&](const T_4P &abcd) {
      if (DuplicateFreeQ(abcd)) {
         auto [a, b, c, d] = abcd;

         Tetrahedron T(T4Tddd{a->X, b->X, c->X, d->X});

         double count = 0, cr;
         const double small = 1E-6;
         // 四面体の形に関する条件
         if (T.inradius > particle_spacing / 20       /*内接円半径に関する条件*/
             && T.circumradius < 2 * particle_spacing /*外接円半径に関する条件*/
             && object->isInside_MethodBucket(T.incenter)
             //  && std::get<0>(object->isInside_MethodOctree(T.circumcenter)) /**/
             //  && tetra_buckets.none_of([&](const auto &t) { return IntersectQ(T, t->scaled(0.9999)); })
             && tetra_buckets.none_of(T.incenter, 3 * particle_spacing, [&](const auto &t) { return IntersectQ(T, t->scaled(1 - 1E-10)); })) {
            /* ----------------------------------------------- */
            //% 外接級の内部に点がある場合はペナルティー
            const double range = 3 * particle_spacing;
            double ACCUM = 0., excess = 0;
            buckets.apply(T.incenter, range, [&](const auto &p) {
               auto tmp = (T.circumradius - Distance(p->X, T.circumcenter));  // / T.circumradius;
               if (tmp > 0)
                  ACCUM += tmp;
            });
            // if (ACCUM > 0.5)
            //    return 1E+50;

            //% 歪な形状の場合ぺネルティー
            ACCUM += Norm(T.circumcenter - T.incenter) / T.inradius;

            // //% 歪な形状の場合ぺネルティー
            ACCUM = ACCUM * log10(1 / Min(T.solidangles));

            if (!isFinite(ACCUM))
               return 1E+50;
            else
               return ACCUM;
         } else if (
             tetra_buckets.none_of(T.incenter, 3 * particle_spacing, [&](const auto &t) { return IntersectQ(T, t->scaled(1 - 1E-10)); }) &&
             std::ranges::all_of(ConnectedQ(a, b, c), [](auto v) { return v; }) &&
             std::ranges::all_of(ConnectedQ(a, b, d), [](auto v) { return v; }) &&
             std::ranges::all_of(ConnectedQ(a, c, d), [](auto v) { return v; }) &&
             std::ranges::all_of(ConnectedQ(b, c, d), [](auto v) { return v; }))
            return 0.;
      }
      return 1E+50;
   };
   // b% -------------------------------------------------------------------------- */
   auto gen1 = [&](const networkFace *f,
                   const auto &candidates_points,
                   const Tddd &n = {0., 0., 0.} /*serach direction*/) {
      double SCORE = 1E+20;
      int count = 0;
      auto [a, b, c] = f->getPoints();
      T_4P abcd = {nullptr, nullptr, nullptr, nullptr}, ABCD = {nullptr, nullptr, nullptr, nullptr};
      Tddd f_center = Mean(ToX(f->getPoints()));
      auto normal = TriangleNormal(ToX(f->getPoints()));
      auto condition = [&](const Tddd &X, const Tddd &search_direction) {
         return Norm(search_direction) < 0.1 || Dot(X - f_center, search_direction) > 0;
      };
      auto update_score = [&](const T_4P &abcd) {
         auto accum = criterion(abcd);
         if (SCORE >= accum) {
            SCORE = accum;
            ABCD = abcd;
            count++;
         }
      };
      if (f) {
         auto [t0, t1] = f->Tetras;
         if (t0 && !t1) {
            auto empty_direction = Projection(Normalize(f_center - t0->center), normal);
            for (const auto &d : candidates_points)
               if (condition(d->X, n) && condition(d->X, empty_direction))
                  update_score({a, b, c, d});
         } else if (!t0 && t1) {
            auto empty_direction = Projection(Normalize(f_center - t1->center), normal);
            for (const auto &d : candidates_points)
               if (condition(d->X, n) && condition(d->X, empty_direction))
                  update_score({a, b, c, d});
         } else if (!t0 && !t1) {
            for (const auto &d : candidates_points)
               if (condition(d->X, n))
                  update_score({a, b, c, d});
         }
      }
      if (count) {
         auto ret = genTetra(object, ABCD);
         if (std::get<0>(ret)) {
            tetra_buckets.add(std::get<1>(ret)->incenter, std::get<1>(ret));
            return ret;
         }
      }
      return std::tuple<bool, networkTetra *>{false, nullptr};
   };
   // b% -------------------------------------------------------------------------- */
   Print("形成された面の内，テトラを持たない面の法線方向において，最も近い点を候補とし，新たにテトラを作成できないか調べる");

   int i = 0;
   PVDWriter pvd(_HOME_DIR_ + "/output/tetras.pvd");
   std::vector<networkPoint *> points_inside;
   std::unordered_map<networkPoint *, std::tuple<Tddd, double>> map_points_inside;
   for (auto k = 0; k < 10; k++) {
      auto faces = object->getFaces();
      for (const auto &f : faces) {

         auto get_points_inside = [&](const Tddd &X, const networkPoint *except_p = nullptr) {
            points_inside.clear();
            auto r = Circumradius(Append(ToX(f), X));
            auto Xc = Circumcenter(Append(ToX(f), X));
            buckets.apply(Xc, r, [&](const auto &p) {
               if (Norm(p->X - Xc) < r * 0.999 && except_p != p)
                  points_inside.emplace_back(p);
            });
         };

         auto gen_point = [&](const Tddd &X) -> networkPoint * {
            auto p = new networkPoint(object, X);
            buckets.add(p->X, p);
            return p;
         };

         auto gen_tetra = [&](const auto &p) {
            auto [is_generated, t] = genTetra(object, Append(f->getPoints(), p));
            if (is_generated)
               tetra_buckets.add(t->incenter, t);
         };

         Tddd first_X, empty_direction = {0., 0., 0.};
         double mean_len = Mean(extLength(f->getLines()));  // 適当な長さ

         auto [t0, t1] = f->Tetras;
         if (!t0 && !t1) {
            // Print("どちら側でもいいのでとりあえず");
            first_X = Circumcenter(ToX(f)) - mean_len * f->normal;
            empty_direction = -f->normal;
            if (!object->isInside_MethodBucket(first_X)) {
               first_X = Circumcenter(ToX(f)) + mean_len * f->normal;
               empty_direction = f->normal;
            }
         } else if (t0 && !t1) {
            // Print("t1 direction");
            empty_direction = Projection(Normalize(f->center - t0->center), f->normal);
            first_X = Circumcenter(ToX(f)) + mean_len * Normalize(empty_direction);
         } else if (!t0 && t1) {
            // Print("t0 direction");
            empty_direction = Projection(Normalize(f->center - t1->center), f->normal);
            first_X = Circumcenter(ToX(f)) + mean_len * Normalize(empty_direction);
         } else if (t0 && t1)
            continue;

         // 対象となる点の候補を取得
         {
            points_inside.clear();
            auto X = first_X;
            auto r = Circumradius(Append(ToX(f), X));
            auto Xc = Circumcenter(Append(ToX(f), X));
            // 異常に大きい半径のものは，排除する
            if (r < 10 * Mean(extLength(f->getLines())))
               buckets.apply(Xc, 1.5 * r, [&](const auto &p) {
                  if (Dot(p->X - f->center, empty_direction) > 0.)
                     if (!MemberQ(f->getPoints(), p))
                        points_inside.emplace_back(p);
               });
         }

         auto condition = [&](const Tddd &X) {
            Tetrahedron T(Append(ToX(f), X));
            auto [p0, p1, p2, p3] = T.vertices;
            auto S = InSphere(p0, p1, p2, p3);
            auto s0 = InSphere(p0, (p0 + p1) / 2, (p0 + p2) / 2, (p0 + p3) / 2);
            auto s1 = InSphere(p1, (p1 + p0) / 2, (p1 + p2) / 2, (p1 + p3) / 2);
            auto s2 = InSphere(p2, (p2 + p0) / 2, (p2 + p1) / 2, (p2 + p3) / 2);
            auto s3 = InSphere(p3, (p3 + p0) / 2, (p3 + p1) / 2, (p3 + p2) / 2);
            if (object->isInside_MethodBucket(X) &&
                object->isInside_MethodBucket(T.centroid) &&
                // 周辺のテトラの辺が，新たに作成するテトラの内接球に入っていないかチェック
                tetra_buckets.none_of(X, T.circumradius, [&](const auto &t) {
                   return IntersectQ(S, T6T2Tddd(*t)) ||
                          IntersectQ(s0, T6T2Tddd(*t)) ||
                          IntersectQ(s1, T6T2Tddd(*t)) ||
                          IntersectQ(s2, T6T2Tddd(*t)) ||
                          IntersectQ(s3, T6T2Tddd(*t));
                }) &&
                // 新たに作成するテトラの辺が，周辺のテトラの内接球に入っていないかチェック
                tetra_buckets.none_of(X, T.circumradius, [&](const auto &t) {
                   auto [p0, p1, p2, p3] = t->vertices;
                   auto S = InSphere(p0, p1, p2, p3);
                   auto s0 = InSphere(p0, (p0 + p1) / 2, (p0 + p2) / 2, (p0 + p3) / 2);
                   auto s1 = InSphere(p1, (p1 + p0) / 2, (p1 + p2) / 2, (p1 + p3) / 2);
                   auto s2 = InSphere(p2, (p2 + p0) / 2, (p2 + p1) / 2, (p2 + p3) / 2);
                   auto s3 = InSphere(p3, (p3 + p0) / 2, (p3 + p1) / 2, (p3 + p2) / 2);
                   return IntersectQ(S, T6T2Tddd(T)) ||
                          IntersectQ(s0, T6T2Tddd(T)) ||
                          IntersectQ(s1, T6T2Tddd(T)) ||
                          IntersectQ(s2, T6T2Tddd(T)) ||
                          IntersectQ(s3, T6T2Tddd(T));
                })
                //  &&
                //  buckets.none_of(X, T.circumradius, [&](const auto &p) { return Norm(X - T.circumcenter) < T.circumradius * 0.9 && !MemberQ(f->getPoints(), p); })
            ) {
               if (Norm(empty_direction) < 0.1)
                  return true;
               else if (t0 && Dot(Projection(t0->center - f->center, f->normal), Projection(X - f->center, f->normal)) > 0.)
                  return false;
               else if (t1 && Dot(Projection(t1->center - f->center, f->normal), Projection(X - f->center, f->normal)) > 0.)
                  return false;
               else if (t0 && t1)
                  return false;
               else
                  return Dot(X - f->center, Dot(empty_direction, f->normal) * f->normal) > 0;
            }
            return false;
         };

         /*

         ## score

         作成するテトラの外接球に他の点がどれほど食い込むかを評価し，最もスコアが良いものを選択する．

         1. `Tetrahedron`を作成する．

         作成するテトラが，ポリゴンの外部にある場合は，排除する．
         作成するテトラが，既存のテトラと干渉する場合は，排除する．

         2. 四面体の中心座標と外接球を使って，Bspline関数定義し，周囲の点上の値だけ引いていく．このスコアが最も良いものを採用する．

         */

         auto score = [&](const Tddd &X) -> double {
            Tetrahedron T(Append(ToX(f), X));

            // 排除する条件1
            // 外接円が大きすぎる場合は排除する
            if (T.circumradius > 0.05)
               return -1E+10;

            // 排除する条件2
            // ポリゴンの外部にある場合は排除する
            if (!object->isInside_MethodBucket(X) || !object->isInside_MethodBucket(T.centroid))
               return -1E+10;

            // 排除する条件3
            // 周辺のテトラと干渉する場合は排除する
            auto [p0, p1, p2, p3] = T.vertices;
            auto S = InSphere(p0, p1, p2, p3);
            auto s0 = InSphere(p0, (p0 + p1) / 2, (p0 + p2) / 2, (p0 + p3) / 2);
            auto s1 = InSphere(p1, (p1 + p0) / 2, (p1 + p2) / 2, (p1 + p3) / 2);
            auto s2 = InSphere(p2, (p2 + p0) / 2, (p2 + p1) / 2, (p2 + p3) / 2);
            auto s3 = InSphere(p3, (p3 + p0) / 2, (p3 + p1) / 2, (p3 + p2) / 2);
            if (tetra_buckets.any_of(X, T.circumradius, [&](const auto &t) {
                   auto [p0, p1, p2, p3] = t->vertices;
                   auto S_ = InSphere(p0, p1, p2, p3);
                   auto s0_ = InSphere(p0, (p0 + p1) / 2, (p0 + p2) / 2, (p0 + p3) / 2);
                   auto s1_ = InSphere(p1, (p1 + p0) / 2, (p1 + p2) / 2, (p1 + p3) / 2);
                   auto s2_ = InSphere(p2, (p2 + p0) / 2, (p2 + p1) / 2, (p2 + p3) / 2);
                   auto s3_ = InSphere(p3, (p3 + p0) / 2, (p3 + p1) / 2, (p3 + p2) / 2);
                   return IntersectQ(S_, T6T2Tddd(T)) || IntersectQ(s0_, T6T2Tddd(T)) || IntersectQ(s1_, T6T2Tddd(T)) || IntersectQ(s2_, T6T2Tddd(T)) || IntersectQ(s3_, T6T2Tddd(T)) ||
                          IntersectQ(S, T6T2Tddd(*t)) || IntersectQ(s0, T6T2Tddd(*t)) || IntersectQ(s1, T6T2Tddd(*t)) || IntersectQ(s2, T6T2Tddd(*t)) || IntersectQ(s3, T6T2Tddd(*t));
                   //   IntersectQ(T, *t);
                }))
               return -1E+10;

            // 新たに作成するテトラの外接球に他の点がどれほど食い込むかを評価し，最もスコアが良いものを選択する．
            double score = 0.;
            buckets.apply(T.circumcenter, T.circumradius, [&](const auto &p) {
               score -= w_Bspline(Norm(p->X - T.circumcenter), T.circumradius) / T.circumradius;
            });

            score -= T.circumradius / std::sqrt(f->area / 2.);

            return score;
         };

         /*

         まずは理想的な位置に点を置いてみる．
         ただし，近くに既に点がある場合は，それを使う．

         * 四面体の干渉チェック
         * 四面体が外部に飛び出していないかチェック

         */

         map_points_inside.clear();
         map_points_inside[nullptr] = {first_X, score(first_X) - 1};

         for (const auto &p : points_inside)
            map_points_inside[p] = {p->X, score(p->X)};

         auto find_best_score = [&]() {
            std::tuple<networkPoint *, Tddd, double> best_p = {nullptr, {0., 0., 0.}, -1E+10};
            networkPoint *min_p = nullptr;
            for (const auto &[p, X_s] : map_points_inside) {
               auto [X, s] = X_s;
               // closer to zero (maximum) is better
               if (std::get<2>(best_p) < s) {
                  std::get<0>(best_p) = p;
                  std::get<1>(best_p) = X;
                  std::get<2>(best_p) = s;
               }
            }
            return best_p;
         };

         auto [p, X, s] = find_best_score();
         // print info about p, X, and s
         if (s < -1E+5)
            continue;
         else {
            if (p == nullptr)
               p = gen_point(X);

            std::cout << "p = " << p << ", X = " << X << ", s = " << s << std::endl;
            gen_tetra(p);
         }

         if (count++ % 50 == 0) {
            std::cout << Magenta << i << Green << ", Elapsed time : " << timer() << colorOff << std::endl;
            auto filename = _HOME_DIR_ + "/output/faces" + std::to_string(i) + ".vtp";
            std::ofstream ofs(filename);
            vtkPolygonWrite(ofs, object->getTetras());
            pvd.push(filename, count);
            pvd.output();
            ofs.close();

            i++;
         }
      }
   };
};
