#include <utility>
#define DEM
#include "Network.hpp"
#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

std::vector<T4Tddd> toCubeFaces(const auto &accum) {
   std::vector<T4Tddd> cube(6 * accum.size());
   int i = 0;
   for (const auto &atree : accum) {
      auto [X0, X1] = std::get<0>(atree->bounds);
      auto [Y0, Y1] = std::get<1>(atree->bounds);
      auto [Z0, Z1] = std::get<2>(atree->bounds);
      cube[i++] = T4Tddd{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}};
      cube[i++] = T4Tddd{{X0, Y0, Z1}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}};
      cube[i++] = T4Tddd{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}};
      cube[i++] = T4Tddd{{X0, Y1, Z0}, {X1, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}};
      cube[i++] = T4Tddd{{X0, Y0, Z0}, {X0, Y0, Z1}, {X0, Y1, Z1}, {X0, Y1, Z0}};
      cube[i++] = T4Tddd{{X1, Y0, Z0}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X1, Y1, Z0}};
   }
   return cube;
};

void Add(std::unordered_set<T_PP> &lines, const T_PP &l) {
   if (lines.find({std::get<0>(l), std::get<1>(l)}) == lines.end() && lines.find({std::get<1>(l), std::get<0>(l)}) == lines.end())
      lines.emplace(l);
};
/* -------------------------------------------------------------------------- */
double Distance(const Tddd &X, const T3Tddd &X012) {
   auto itx = IntersectionSphereTriangle(X, 1E+20, X012);
   return Norm(itx.X - X);
};
double Distance(const T3Tddd &X012, const Tddd &X) {
   auto itx = IntersectionSphereTriangle(X, 1E+20, X012);
   return Norm(itx.X - X);
};
/* -------------------------------------------------------------------------- */
/**
 * ./main bunny.obj 7 1 3
 */
/* -------------------------------------------------------------------------- */
#define check_voronoi0
#if defined(check_voronoi0)
int main(int arg, char **argv) {
   /* -------------------------------------------------------------------------- */
   // 四面体の許容サイズは適切か確認（max_radius_of_sphere）
   Timer timer;
   std::string name = "./bunny.obj";
   /* -------------------------------------------------------------------------- */
   int n = 10;
   auto object = new Network(name, "object");
   auto initialFaces = object->getFaces();
   object->genOctreeOfFaces({1, 5}, 1);
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
            // if (std::get<0>(object->isInside_MethodOctree(xyz)))
            //    XYZ.push_back(xyz);
            /* ---------------------------------------------------------- */
            auto [isinside, cell, f] = object->isInside_MethodOctree(xyz);
            if (isinside) {
               if (!cell->faces_.empty()) {
                  if (Norm(Nearest(xyz, cell->faces_) - xyz) > particle_spacing)
                     XYZ.push_back(xyz);
               } else
                  XYZ.push_back(xyz);
            }
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
   std::cout << Green << "Elapsed time : " << timer() << colorReset << std::endl;
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
             && std::get<0>(object->isInside_MethodOctree(T.incenter))
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
            auto empty_direction = Projection(Normalize(f_center - t0->centroid), normal);
            for (const auto &d : candidates_points)
               if (condition(d->X, n) && condition(d->X, empty_direction))
                  update_score({a, b, c, d});
         } else if (!t0 && t1) {
            auto empty_direction = Projection(Normalize(f_center - t1->centroid), normal);
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
   std::vector<networkTetra *> tetras;
   for (const auto &t : object->getTetras())
      tetras.emplace_back(t);

   /*

   外接球に他の点が入らないか確認しながら，面に対してある方向に徐々に点を移動させていく．
   比較的綺麗な四面体を作れる範囲に，既存の点が存在しない場合，その点を新たに追加し，四面体を作成する．

   */
   std::vector<networkPoint *> points_inside;
   for (auto k = 0; k < 10; k++) {
      auto faces = object->getFaces();
      for (const auto &f : faces) {

         /*
         四面体を生成できる条件

         IntersectQ()がfalseである場合
         -> 生成可能．

         IntersectQ()がtrueである場合
         -> 面が共通しており，接触している面の反対側に四面体が作られる場合は生成可能．



         */

         auto [t0, t1] = f->Tetras;
         if (t0 && t1)
            continue;

         double mean_len = Mean(extLength(f->getLines()));  // 適当な長さ

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

         Tddd first_X = Circumcenter(ToX(f)) - mean_len * f->normal, empty_direction = {0., 0., 0.};

         if (t0 && !t1) {
            empty_direction = Projection(Normalize(f->center - t0->centroid), f->normal);
            first_X = Circumcenter(ToX(f)) + mean_len * Normalize(empty_direction);
         } else if (!t0 && t1) {
            empty_direction = Projection(Normalize(f->center - t1->centroid), f->normal);
            first_X = Circumcenter(ToX(f)) + mean_len * Normalize(empty_direction);
         }

         auto condition = [&](const Tddd &X) {
            if (Norm(empty_direction) < 0.1)
               return true;
            Tddd n = f->normal;
            if (Dot(empty_direction, f->normal) < 0.)
               n *= -1;
            return Dot(X - f->center, n) > 0;
         };

         if (points_inside.empty() && std::get<0>(object->isInside_MethodOctree(first_X))) {
            auto p = gen_point(first_X);
            gen_tetra(p);
         } else {
            auto tmp = points_inside;
            for (const auto &p : tmp) {
               get_points_inside(p->X, p);
               if (condition(p->X))
                  if (points_inside.empty() && !MemberQ(f->getPoints(), p)) {
                     auto [isinside, _, __] = object->isInside_MethodOctree(Incenter(Append(ToX(f), p->X)));
                     if (isinside)
                        gen_tetra(p);
                  }
            }
         }

         if (count++ % 50 == 0) {
            std::cout << Magenta << i << Green << ", Elapsed time : " << timer() << colorReset << std::endl;
            std::ofstream ofs(_HOME_DIR_ + "/output/tetras" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, object->getTetras());
            ofs.close();
            i++;
         }
      }
   };
   if (false)
      for (auto I = 0; I < 5000; I++) {
         bool found = false;
         int count = 0;
         /* -------------------------------------------------------------------------- */
         if (I == 0) {
            for (const auto &f : initialFaces) {
               auto [ismade, tet] = gen1(f, buckets.getBucket(Mean(ToX(f->getPoints())), 1.5 * particle_spacing), -f->normal);
               if (ismade) {
                  tetras.emplace_back(tet);
                  std::cout << "add tetra " << count << ", object->getTetras().size() = " << object->getTetras().size() << std::endl;
                  // if (count++ > 100)
                  //    break;
               }
            }
            std::cout << Magenta << i << Green << ", Elapsed time : " << timer() << colorReset << std::endl;
            std::ofstream ofs(_HOME_DIR_ + "/output/tetras" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, object->getTetras());
            ofs.close();
            i++;
         }
         count = 0;
         /* -------------------------------------------------------------------------- */
         for (const auto &f : object->getFaces()) {
            auto [ismade, tet] = gen1(f, buckets.getBucket(Mean(ToX(f->getPoints())), 2. * particle_spacing));
            if (ismade) {
               tetras.emplace_back(tet);
               std::cout << "add tetra " << count << ", object->getTetras().size() = " << object->getTetras().size() << std::endl;
               if (count++ > 100)
                  break;
            }
         }
         std::cout << Magenta << i << Green << ", Elapsed time : " << timer() << colorReset << std::endl;
         std::ofstream ofs(_HOME_DIR_ + "/output/tetras" + std::to_string(i) + ".vtp");
         vtkPolygonWrite(ofs, object->getTetras());
         ofs.close();
         i++;
      }
};
#elif defined(check_voronoi)
int main(int arg, char **argv) {
   Timer timer;
   /* ----------------------- 引数の読み込み ---------------------- */
   std::string name{argv[1]};  //"/Users/tomoaki/Dropbox/markdown/cpp/obj/cow.obj";
   int maxDepth = std::atoi(argv[2]);
   int minNumber = std::atoi(argv[3]);
   double spacing = std::atof(argv[4]);
   std::cout << "maxDepth = " << maxDepth << std::endl;
   std::cout << "minNumber = " << minNumber << std::endl;
   /* ------------------------------------------------------ */
   auto net = new Network;
   double particle_spacing = spacing;
   auto [xb0, xb1] = Tdd{-10., 10.};
   auto [yb0, yb1] = Tdd{-10., 10.};
   auto [zb0, zb1] = Tdd{-10., 10.};
   auto X = Subdivide(xb0, xb1, (int)std::round((xb1 - xb0) / particle_spacing));  // 粒子のX座標
   auto Y = Subdivide(yb0, yb1, (int)std::round((yb1 - yb0) / particle_spacing));  // 粒子のY座標
   auto Z = Subdivide(zb0, zb1, (int)std::round((zb1 - zb0) / particle_spacing));  // 粒子のZ座標
   Tdd r = {-0.1, 0.1};
   for (const auto &x : X)
      for (const auto &y : Y)
         for (const auto &z : Z)
            new networkPoint(
                net, net,
                {x + RandomReal(r), y + RandomReal(r), z + RandomReal(r)});
   net->setGeometricProperties();
   std::cout << "net points size:" << net->getPoints().size() << std::endl;
   vtkPolygonWrite("./output/points.vtp", extX(net->getPoints()));
   /* ------------------------------------------------------ */
   octree<networkPoint *> tree(net->getBounds(), ToVector(net->getPoints()), {maxDepth, minNumber});  // 八分木構造を生成．ポリゴン外部のキューブは自動で内部で削除せれる
   std::cout << "timer : " << timer() << std::endl;
   auto tmp = tree.getDepth(maxDepth);
   std::cout << "deepest cells :" << tmp.size() << std::endl;
   vtkPolygonWrite("./output/cubeCenter.vtp", toX(tmp));
   /* ------------------------------------------------------ */
   std::unordered_map<networkPoint *, std::unordered_set<networkPoint *>> p_neighbors;  // テンプレートはさまざまな制約を生み，エラーの原因となる
   for (const auto &c : tmp) {
      if (c->faces_.size() == 2) {
         auto p0 = c->faces_[0];
         auto p1 = c->faces_[1];
         p_neighbors[p0].emplace(p1);
         p_neighbors[p1].emplace(p0);
      } else if (c->faces_.size() == 3) {
         auto p0 = c->faces_[0];
         auto p1 = c->faces_[1];
         auto p2 = c->faces_[2];
         p_neighbors[p0].emplace(p1);
         p_neighbors[p1].emplace(p2);
         p_neighbors[p2].emplace(p0);
         p_neighbors[p0].emplace(p2);
         p_neighbors[p1].emplace(p0);
         p_neighbors[p2].emplace(p1);
      }
   }
   std::vector<T2Tddd> lines;
   /* ----------------------- 引数の読み込み ---------------------- */
   std::string name{
       argv[1]};  //"/Users/tomoaki/Dropbox/markdown/cpp/obj/cow.obj";
   int maxDepth = std::atoi(argv[2]);
   int minNumber = std::atoi(argv[3]);
   std::cout << "maxDepth = " << maxDepth << std::endl;
   std::cout << "minNumber = " << minNumber << std::endl;
   auto net = new Network;
   /* ------------------------------------------------------ */
   /*              八分木構造と球体の干渉のチェックと出力           */
   /* ------------------------------------------------------ */
   double particle_spacing = 0.4;
   auto [xb0, xb1] = Tdd{-1., 1.};
   auto [yb0, yb1] = Tdd{-1., 1.};
   auto [zb0, zb1] = Tdd{-1., 1.};
   auto X = Subdivide(xb0, xb1, (int)std::round((xb1 - xb0) / particle_spacing * 2));  // 粒子のX座標
   auto Y = Subdivide(yb0, yb1, (int)std::round((yb1 - yb0) / particle_spacing));      // 粒子のY座標
   auto Z = Subdivide(zb0, zb1, (int)std::round((zb1 - zb0) / particle_spacing));      // 粒子のZ座標
   Tdd r = {-0.1, 0.1};
   for (const auto &x : X)
      for (const auto &y : Y)
         for (const auto &z : Z)
            new networkPoint(net, net, {x + RandomReal(r), y + RandomReal(r), z + RandomReal(r)});
   net->setGeometricProperties();
   octree tree(net->getBounds(), {maxDepth, minNumber}, ToX(net->getPoints()));  // 八分木構造を生成．ポリゴン外部のキューブは自動で内部で削除せれる
   vtkPolygonWrite("./output/points.vtp", ToX(net->getPoints()));
   vtkPolygonWrite("./output/cubeAll.vtp", toCubeFaces(tree.getAllDeepest()));
};
#else
int main(int arg, char **argv) {
   Timer timer;
   int minNumber = std::atoi(argv[3]);
   /* ---------------------- ポリゴンの読み込み --------------------- */
   auto net = new Network(name, "object");
   std::cout << "{maxDepth, minNumber} = " << Tii{maxDepth, minNumber} << std::endl;
   octree tree(net->getBounds(), {maxDepth, minNumber}, extVertices(net->getFaces()));  // 八分木構造を生成．ポリゴン外部のキューブは自動で内部で削除せれる
   tree.deleteOuside();
   std::cout << "timer : " << timer() << std::endl;
   {
      // vtkPolygonWriter<key here>;
      vtkPolygonWriter</*specidy key*/ networkPoint *> vtp;
      std::unordered_map<networkPoint *, double> values, solidangles;
      std::unordered_map<networkPoint *, Tddd> normals;
      double v = 0;
      for (const auto &f : net->getFaces()) {
         auto [p0, p1, p2] = f->getPointsTuple();
         vtp.add({p0, p1, p2});
         vtp.addPolygon({p0, p1, p2});
         for_each(f->getPointsTuple(), [&](auto &p) { values[p] = Norm(p->getXtuple()); });
         for_each(f->getPointsTuple(), [&](auto &p) { normals[p] = p->getNormalTuple(); });
         for_each(f->getPointsTuple(), [&](auto &p) { solidangles[p] = p->getSolidAngle(); });
      }
      vtp.addPointData("values", values);
      vtp.addPointData("normals", normals);
      vtp.addPointData("solidangles", solidangles);
      std::ofstream ofs("./output/bunny.vtp");
      vtp.write(ofs);
   }

   /* ------------------------- ボクセルの出力 ------------------------- */
   vtkPolygonWrite("./output/cubeAll.vtp", toCubeFaces(tree.getAllDeepestInside()));
   auto trees = tree.getAllDeepestInside();
   /* ------------------------------------------------------ */
   /*              八分木構造と球体の干渉のチェックと出力           */
   /* ------------------------------------------------------ */
   double particle_spacing = 0.04;
   auto [xb0, xb1] = Tdd{-0.05, 0.05};
   auto [yb0, yb1] = Tdd{0.05, 0.15};
   auto [zb0, zb1] = Tdd{0.06, 0.06 + 0.1};
   auto X = Subdivide(xb0, xb1, (int)std::round((xb1 - xb0) / particle_spacing));  // 粒子のX座標
   auto Y = Subdivide(yb0, yb1, (int)std::round((yb1 - yb0) / particle_spacing));  // 粒子のY座標
   auto Z = Subdivide(zb0, zb1, (int)std::round((zb1 - zb0) / particle_spacing));  // 粒子のZ座標
   int i = 0;
   for (const auto &x : X)
      for (const auto &y : Y)
         for (const auto &z : Z) {
            Tddd xyz = {x, y, z};
            std::cout << "timer : " << timer() << std::endl;
            geometry::Sphere sp(xyz, 0.05);  // 球体オブジェク
            vtkPolygonWrite("./output/cube" + std::to_string(i) + ".vtp", toCubeFaces(tree.getIntersectInside(sp)));
            vtkPolygonWrite("./output/particle" + std::to_string(i) + ".vtp", xyz);
            i++;
         }
   //! ------------------------------------------------------ */
   //!                         vtpの出力                       */
   //! ------------------------------------------------------ */
   std::cout << "timer : " << timer() << std::endl;
   {
      vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
      double v = 0;
      Tddd u = {0., 0., 0.};
      for (const auto &t : trees) {
         auto [X0, X1, X2, X3, X4, X5, X6, X7] = t->getVertices();
         std::shared_ptr<Tddd> x0(new Tddd(X0));
         std::shared_ptr<Tddd> x1(new Tddd(X1));
         std::shared_ptr<Tddd> x2(new Tddd(X2));
         std::shared_ptr<Tddd> x3(new Tddd(X3));
         std::shared_ptr<Tddd> x4(new Tddd(X4));
         std::shared_ptr<Tddd> x5(new Tddd(X5));
         std::shared_ptr<Tddd> x6(new Tddd(X6));
         std::shared_ptr<Tddd> x7(new Tddd(X7));
         vtp.add({x0, x1, x2, x3, x4, x5, x6, x7});
         vtp.addPolygon({{x0, x1, x2, x3},
                         {x4, x5, x6, x7},
                         {x0, x1, x5, x4},
                         {x2, x3, x7, x6},
                         {x0, x4, x7, x3},
                         {x1, x2, x6, x5}});
         vtp.addPointData("scaler0", {{x0, v},
                                      {x1, v},
                                      {x2, v},
                                      {x3, v},
                                      {x4, v},
                                      {x5, v},
                                      {x6, v},
                                      {x7, v}});
         v += 1.;
         vtp.addPointData("scaler1", {{x0, v},
                                      {x1, v},
                                      {x2, v},
                                      {x3, v},
                                      {x4, v},
                                      {x5, v},
                                      {x6, v},
                                      {x7, v}});
         v += 1.;
         vtp.addPointData("vector", {{x0, X0 - u},
                                     {x1, X1 - u},
                                     {x2, X2 - u},
                                     {x3, X3 - u},
                                     {x4, X4 - u},
                                     {x5, X5 - u},
                                     {x6, X6 - u},
                                     {x7, X7 - u}});
      }
      std::cout << vtp.verticies.size() << std::endl;
      std::cout << "outputting " << std::endl;
      std::ofstream ofs("./output/output2.vtp");
      vtp.write(ofs);
   }
};

#endif