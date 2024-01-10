#include <utility>
#define DEM
#include "Network.hpp"
#include "SPH_weightingFunctions.hpp"
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

void Add(std::unordered_set<std::tuple<networkPoint *, networkPoint *>> &lines,
         const std::tuple<networkPoint *, networkPoint *> &l) {
   std::tuple<networkPoint *, networkPoint *> l0 = {std::get<0>(l), std::get<1>(l)};
   std::tuple<networkPoint *, networkPoint *> l1 = {std::get<1>(l), std::get<0>(l)};
   if (lines.find(l0) == lines.end() && lines.find(l1) == lines.end())
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
/* -------------------------------------------------------------------------- */
#define check_SPH3
#ifdef check_SPH3
int main(int argc, char **argv) {
   if (argc < 4) {
      std::stringstream ss;
      ss << "number of arguments: " << argc << Green << " -> $./main bunny.obj 5 0";
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   }
   /* ----------------------------------------------------------- */
   Timer timer;
   /* ----------------------- 引数の読み込み ---------------------- */
   std::string name{argv[1]};  //"/Users/tomoaki/Dropbox/markdown/cpp/obj/cow.obj";
   int minDepth = std::atoi(argv[2]);
   int maxDepth = std::atoi(argv[3]);
   std::cout << "minDepth = " << minDepth << std::endl;
   std::cout << "maxDepth = " << maxDepth << std::endl;
   std::cout << "timer : " << timer() << std::endl;
   /* ---------------------- ポリゴンの読み込み --------------------- */
   auto net = new Network(name, "object");
   net->translate({0, -0.1, 0});
   net->rotate(M_PI / 2., {1, 0, 0});
   octree tree(net->getBounds(), {minDepth, maxDepth}, 1, ToX(net->getFaces()));  // 八分木構造を生成．ポリゴン外部のキューブは自動で内部で削除せれる
   tree.deleteOuside();
   std::cout << "timer : " << timer() << ", 八分木構造を生成" << std::endl;
   if (false) {
      vtkPolygonWriter<networkPoint *> vtp;
      for (const auto &f : net->getFaces()) {
         auto [p0, p1, p2] = f->getPointsTuple();
         vtp.add(p0, p1, p2);
         vtp.addPolygon(p0, p1, p2);
      }
      std::ofstream ofs("./output/bunny.vtp");
      vtp.write(ofs);
      ofs.close();
      std::cout << "timer : " << timer() << std::endl;
      mk_vtu("./output/cubeAll.vtu", toCubeFaces(tree.getAllDeepestInside()));
      //* ----------------------------------------------------------- */
      //*                         隣接するセルを保存                     */
      //* ----------------------------------------------------------- */
      std::cout << "timer : " << timer() << std::endl;
      tree.setNeighbors();
      std::cout << "timer : " << timer() << ", tree.setNeighbors()" << std::endl;
      tree.setVectorsToTriangle();
      std::cout << "timer : " << timer() << ", tree.setVectorsToTriangle()" << std::endl;
   }
   // b! ------------------------------------------------------ */
   // b!           signed distance の結果を出力                   */
   // b! ------------------------------------------------------ */
   if (false) {
      vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
      for (const auto &cell : tree.getAllDeepestInside()) {
         auto [X0, X1, X2, X3, X4, X5, X6, X7] = cell->getVertices();
         std::shared_ptr<Tddd> p0(new Tddd(X0)), p1(new Tddd(X1)), p2(new Tddd(X2)), p3(new Tddd(X3)), p4(new Tddd(X4)), p5(new Tddd(X5)), p6(new Tddd(X6)), p7(new Tddd(X7));
         vtp.add({p0, p1, p2, p3, p4, p5, p6, p7});
         vtp.addPolygon({{p4, p6, p7, p5}, {p0, p2, p6, p4}, {p2, p3, p7, p6}, {p3, p1, p5, p7}, {p0, p4, p5, p1}, {p0, p1, p3, p2}});
         auto [d0, d1, d2, d3, d4, d5, d6, d7] = cell->scalers;
         vtp.addPointData("distance", {{p0, d0}, {p1, d1}, {p2, d2}, {p3, d3}, {p4, d4}, {p5, d5}, {p6, d6}, {p7, d7}});
         auto [v0, v1, v2, v3, v4, v5, v6, v7] = cell->vectors;
         vtp.addPointData("vectors", {{p0, v0}, {p1, v1}, {p2, v2}, {p3, v3}, {p4, v4}, {p5, v5}, {p6, v6}, {p7, v7}});
         //
         auto X = Mean(cell->getVertices());
         std::shared_ptr<Tddd> p(new Tddd(X));
         vtp.add(p);
         auto [x0x1, y0y1, z0z1] = cell->bounds;
         vtp.addPointData("interpolated vector", {{p, cell->Interpolate(X, cell->vectors) /*可視化の際は注意*/}});
      }
      std::cout << vtp.verticies.size() << std::endl;
      std::ofstream ofs("./output/output2.vtp");
      vtp.write(ofs);
      ofs.close();
      std::cout << "timer : " << timer() << ", output2.vtp" << std::endl;
   }
   //  b* ------------------------------------------------------ */
   //  b*              八分木構造と球体の干渉のチェックと出力          */
   //  b* ------------------------------------------------------ */
   if (false) {
      //! 干渉のチェックを行い，球に入るキューブを抜き出す
      int i = 0;
      auto [x01, y01, z01] = net->getBounds();
      for (const auto &x : Subdivide(x01, 3))
         for (const auto &y : Subdivide(y01, 3))
            for (const auto &z : Subdivide(z01, 3)) {
               Tddd xyz = {x, y, z};  // 平滑化距離の範囲
               std::cout << "timer : " << timer() << std::endl;
               {
                  double h = net->getScale() / 20.;
                  Sphere sp(xyz, h);  // 球体オブジェクト
                  vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
                  for (const auto &cell : tree.getIntersectInside(sp)) {
                     for_each((T6T4Tddd)(*cell), [&](const auto &aT4Tddd) {
                        auto [X0, X1, X2, X3] = aT4Tddd;
                        std::shared_ptr<Tddd> x0(new Tddd(X0)), x1(new Tddd(X1)), x2(new Tddd(X2)), x3(new Tddd(X3));
                        vtp.add(x0, x1, x2, x3);
                        vtp.addPolygon(x0, x1, x2, x3);
                        //
                        auto X = Mean(cell->getVertices());
                        std::shared_ptr<Tddd> p(new Tddd(X));
                        vtp.add(p);
                        auto [x0x1, y0y1, z0z1] = cell->bounds;
                        auto v = cell->Interpolate(X, cell->vectors);
                        vtp.addPointData("interpolated vector", {{p, v}});
                     });
                  }
                  std::ofstream ofs("./output/cube" + std::to_string(i) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
               {
                  vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
                  std::shared_ptr<Tddd> x(new Tddd(xyz));
                  vtp.add(x);
                  std::ofstream ofs("./output/smoothing_length" + std::to_string(i) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
               i++;
            }
   }
   if (false) {
      // b$ ------------------------------------------------------ */
      // b$       SPHで，Rigid Bodyとして使うような点群を生成           */
      // b$ ------------------------------------------------------ */
      auto RigidBody = new Network;
      auto [x01, y01, z01] = net->getBounds();
      // auto [x01, y01, z01] = T3Tdd{{-0.1, 0.1}, {-0.1, 0.1}, {-0.1, 0.1}};
      for (const auto &x : Subdivide(x01, 10))
         for (const auto &y : Subdivide(y01, 10))
            for (const auto &z : Subdivide(z01, 10)) {
               Tddd xyz = {x, y, z};
               if (tree.isIntersectInside(xyz))
                  new networkPoint(RigidBody, xyz);
            }
      std::cout << "timer : " << timer() << ", networkPointの生成" << std::endl;
      std::cout << "RigidBody->getPoints().size()=" << RigidBody->getPoints().size() << std::endl;
      vtkPolygonWriter<networkPoint *> vtp;
      for (const auto &p : RigidBody->getPoints())
         vtp.add(p);
      std::ofstream ofs("./output/particle_inside.vtp");
      vtp.write(ofs);
      ofs.close();
      /* --------------------------------------------------------------
         1) 点を内包するセルを，ツリーから抽出
         2) 点上のベクトルを，セルのベクトルで補間
      ----------------------------------------------------------------- */
      {
         // signed_distance_vectorを設定
         // SPHを実行
         std::unordered_map<networkPoint *, Tddd> data;
         for (const auto &p : RigidBody->getPoints()) {
            auto cell = tree.getIntersect(ToX(p));
            p->signed_distance_vector = cell->Interpolate(ToX(p), cell->vectors);
            // if (Between(Norm(p->signed_distance_vector), {-1E-5, 1E-5}))
            //    data[p] = {0., 0., 0.};
            // else
            data[p] = p->signed_distance_vector;  //{0., 0., 0.};
         }
         vtkPolygonWriter<networkPoint *> vtp;
         for (const auto &p : ToVector(RigidBody->getPoints()))
            vtp.add(p);
         vtp.addPointData("vector", data);
         std::ofstream ofs("./output/particle_inside_with_signed_vector.vtp");
         vtp.write(ofs);
         ofs.close();
      }
      std::cout << "timer : " << timer() << ", networkPointの生成" << std::endl;
      /* --------------------------------------------------------------
         1) RigidBodyでmakeBucketPoints
         2) バケツを利用して効率的に指定したバケツ内の粒子を取得
      ----------------------------------------------------------------- */
      {
         Timer timer;
         std::cout << "timer : " << timer() << std::endl;
         // RigidBody->makeBucketPoints(0.005);
         RigidBody->makeBucketPoints(net->getScale() / 10.);
         std::cout << "timer : " << timer() << "ある間隔でバケツを作成し，点を保存" << std::endl;
         auto [x01, y01, z01] = RigidBody->getBounds();
         int i = 0;
         for (const auto &x : Subdivide(x01, 2))
            for (const auto &y : Subdivide(y01, 2))
               for (const auto &z : Subdivide(z01, 2)) {
                  vtkPolygonWriter<networkPoint *> vtp;
                  vtp.add(RigidBody->BucketPoints.get(Tddd{x, y, z}, Tdd{0., 0.0005 * i}));
                  std::ofstream ofs("./output/pointsInBuckets" + std::to_string(i++) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
                  std::cout << "timer : " << timer() << std::endl;
               }
      }
   }
   /* -----------------------------------------------------------------------
      1) Fluid作成とRigidBodyを作成
      2)
   -------------------------------------------------------------------------- */
   {
      auto polyFluid = new Network("./OC5water.obj", "water");
      auto ofs0 = std::ofstream("./output/polyFluid.vtp");
      vtkPolygonWrite(ofs0, ToX(polyFluid->getFaces()));

      auto polyWaveTank = new Network("./OC5tank2.obj", "tank");
      auto ofs1 = std::ofstream("./output/polyWaveTank.vtp");
      vtkPolygonWrite(ofs1, ToX(polyWaveTank->getFaces()));

      auto polyWall = new Network("./OC5tank2.obj", "wall");
      auto ofs2 = std::ofstream("./output/polyWall.vtp");
      vtkPolygonWrite(ofs2, ToX(polyWall->getFaces()));
      //
      octree treeFluid(polyFluid->getBounds(), {minDepth, maxDepth}, 1, ToX(polyFluid->getFaces()));
      treeFluid.deleteOuside();
      octree treeWaveTank(polyWaveTank->getBounds(), {minDepth, maxDepth}, 1, ToX(polyWaveTank->getFaces()));
      treeWaveTank.deleteOuside();
      treeWaveTank.setNeighbors();
      treeWaveTank.setVectorsToTriangle();
      //
      auto ofs3 = std::ofstream("./output/treeWaveTank.vtp");
      vtkPolygonWrite(ofs3, toCubeFaces(treeWaveTank.getAllDeepestInside()));
      //
      octree treeWall(polyWall->getBounds(), {minDepth, maxDepth}, 1, ToX(polyWall->getFaces()));
      treeWall.deleteOuside();
      //
      auto Fluid = new Network;
      auto WaveTank = new Network;
      auto Wall = new Network;
      //
      int n = 15;
      
      // Tdd rangeX = {-100, 300}, rangeY = {-10, 10}, rangeZ = {-10, 50};
      Tdd rangeX = {-0.5, 19.5}, rangeY = {-1.5, 1.5}, rangeZ = {-0.5, 2};
      auto vecX = Subdivide(rangeX, std::round(std::get<1>(rangeX) - std::get<0>(rangeX)) * (n - 1));
      auto vecY = Subdivide(rangeY, std::round(std::get<1>(rangeY) - std::get<0>(rangeY)) * (n - 1));
      auto vecZ = Subdivide(rangeZ, std::round(std::get<1>(rangeZ) - std::get<0>(rangeZ)) * (n - 1));
      //
      double particle_spacing = 1. / n;
      double volume = std::pow(particle_spacing, 3);
      double C_SML = 3.1;
      double max_dt = 0.001;
      double initial_surface_z_position = 0.5;

      auto isInside = [](const auto &tree, const Network *net, const Tddd &X) {
         // セルに含まれるかどうかではなく，ポリゴンの内部にあるかどうかで判定するための関する
         auto cells = tree.getIntersectInside(X);
         if (cells.empty())
            return false;
         else {
            if ((cells[0]->faces_.empty()))
               return true;
            else {
               // 最も近い面とdotが同じかどうか
               double distance = 1E+20, dot;
               for (const auto &f : cells[0]->faces_) {
                  auto itx = IntersectionSphereTriangle(X, 1E+20, f);
                  if (distance >= itx.distance) {
                     distance = itx.distance;
                     dot = Dot(itx.X - X, TriangleNormal(f));
                  }
               }
               return dot >= 0.;
            }
         }
      };

      for (const auto &x : vecX)
         for (const auto &y : vecY)
            for (const auto &z : vecZ) {
               Tddd xyz = {x, y, z};
               if (isInside(treeWaveTank, polyWaveTank, xyz) || isInside(treeWall, polyWall, xyz)) {
                  if (isInside(treeWaveTank, polyWaveTank, xyz)) {
                     auto p = new networkPoint(WaveTank, xyz);
                     p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
                     p->radius_SPH = C_SML * particle_spacing;
                     p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
                  }
               } else if (isInside(treeFluid, polyFluid, xyz)) {
                  auto p = new networkPoint(Fluid, xyz);
                  p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
                  p->radius_SPH = C_SML * particle_spacing;
                  p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
               }
            }

      for (const auto &cell : treeWaveTank.getAllDeepestInside()) {
         auto p = new networkPoint(Wall, cell->center);
         p->setDensityVolume(_WATER_DENSITY_, cell->getVolume());  // 質量(mass)は関数内部で自動で決められる
         p->radius_SPH = C_SML * particle_spacing;
         p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
      }

      for (const auto &p : WaveTank->getPoints()) {
         auto cell = treeWaveTank.getIntersect(p->X);
         p->normal_SPH = cell->Interpolate(p->X, cell->vectors);
      }
      //
      Fluid->setGeometricProperties();
      WaveTank->setGeometricProperties();
      Wall->setGeometricProperties();
      std::cout << "timer : " << timer() << ", networkPointの生成" << std::endl;
      std::cout << "Fluid->getPoints().size()=" << Fluid->getPoints().size() << std::endl;
      std::cout << "WaveTank->getPoints().size()=" << WaveTank->getPoints().size() << std::endl;
      std::cout << "Wall->getPoints().size()=" << Wall->getPoints().size() << std::endl;
      {
         auto ofs = std::ofstream("./output/Fluid.vtp");
         vtkPolygonWrite(ofs, ToX(Fluid->getPoints()));
      }
      {
         auto ofs = std::ofstream("./output/Wall.vtp");
         vtkPolygonWrite(ofs, ToX(WaveTank->getPoints()));
      }
      {
         auto ofs = std::ofstream("./output/Wall.vtp");
         vtkPolygonWrite(ofs, ToX(Wall->getPoints()));
      }
      /* -------------------------------------------------------------------------- */
      double real_time = 0.;
      PVDWriter pvdFluidSPH("./output/FluidSPH.pvd");
      PVDWriter pvdWaveTankSPH("./output/WaveTankSPH.pvd");
      PVDWriter pvdWallSPH("./output/Wall.pvd");
      int count = 0, i = 0;
      //
      std::vector<Network *> RigidBodies = {WaveTank};
      for (auto count = 0; count < 10000; ++count) {
         auto bound = CoordinateBounds(WaveTank->bounds) + CoordinateBounds(Fluid->bounds);
         auto points = Fluid->getPoints();
         for (const auto &p : points) {
            if (!bound.isInside(p->X)) {
               std::cout << Blue << "delete" << std::endl;
               delete p;
            }
         }
         std::cout << Blue << "count = " << count << std::endl;
         // if (count < 100) {
         //    for (const auto &p : Fluid->getPoints())
         //       p->mu_SPH = 0.1;
         //    for (const auto &p : Wall->getPoints())
         //       p->mu_SPH = 0.1;
         // } else {
         //    for (const auto &p : Fluid->getPoints())
         //       p->mu_SPH = 0.001138;
         //    for (const auto &p : Wall->getPoints())
         //       p->mu_SPH = 0.001138;
         // }

         setContactPoints(Fluid, {Fluid}, RigidBodies, real_time, C_SML, particle_spacing);
         // 出力
         if (count % 3 == 0) {
            {
               vtkPolygonWriter<networkPoint *> vtp;
               vtp.add(ToVector(Fluid->getPoints()));
               setData(vtp, Fluid);
               std::ofstream ofs("./output/FluidSPH" + std::to_string(i) + ".vtp");
               vtp.write(ofs);
               ofs.close();
            }
            if(count==0)
            {
               vtkPolygonWriter<networkPoint *> vtp;
               vtp.add(ToVector(WaveTank->getPoints()));
               setData(vtp, WaveTank);
               std::ofstream ofs("./output/WaveTankSPH" + std::to_string(i) + ".vtp");
               vtp.write(ofs);
               ofs.close();
            }
            // {
            //    vtkPolygonWriter<networkPoint *> vtp;
            //    vtp.add(ToVector(Wall->getPoints()));
            //    setData(vtp, Wall);
            //    std::ofstream ofs("./output/WallSPH" + std::to_string(i) + ".vtp");
            //    vtp.write(ofs);
            //    ofs.close();
            // }
            pvdFluidSPH.push("./FluidSPH" + std::to_string(i) + ".vtp", real_time);
            pvdFluidSPH.output();
            pvdWaveTankSPH.push("./WaveTankSPH" + std::to_string(i) + ".vtp", real_time);
            pvdWaveTankSPH.output();
            pvdWallSPH.push("./WallSPH" + std::to_string(i) + ".vtp", real_time);
            pvdWallSPH.output();

            i++;
         }
         double alpha;
         // if (count < 200)
         //    alpha = 0.;
         // else
         alpha = 1.;
         developByEISPH(Fluid, {Fluid}, RigidBodies, real_time, C_SML, particle_spacing, alpha);
         std::cout << "real_time = " << real_time << std::endl;
      }
   }
};
#elif defined(check_voronoi0)
int main(int arg, char **argv) {
   Timer timer;
   std::string name(argv[1]);
   int maxDepth(std::atoi(argv[2]));
   int maxNumber(std::atoi(argv[3]));
   double spacing(std::atof(argv[4]));
   std::cout << Grid({"maxDepth", "maxNumber", "spacing"}, 10) << std::endl;
   std::cout << Grid({maxDepth, maxNumber, spacing}, 10) << std::endl;
   /* ------------------------------make------------------------ */
   // auto net = new Network(name, "object");
   //
   double particle_spacing = spacing;
   auto [xb0, xb1] = Tdd{-5., 5.};
   auto [yb0, yb1] = Tdd{-5., 5.};
   auto [zb0, zb1] = Tdd{-5., 5.};
   auto X = Subdivide(xb0, xb1, (int)std::round((xb1 - xb0) / particle_spacing));  // 粒子のX座標
   auto Y = Subdivide(yb0, yb1, (int)std::round((yb1 - yb0) / particle_spacing));  // 粒子のY座標
   auto Z = Subdivide(zb0, zb1, (int)std::round((zb1 - zb0) / particle_spacing));  // 粒子のZ座標
   Tdd r = {-.2, .2};
   std::vector<Tddd> XYZ;
   for (const auto &x : X)
      for (const auto &y : Y)
         for (const auto &z : Z)
            XYZ.push_back({x + RandomReal(r), y + RandomReal(r), z + RandomReal(r)});
   /* ----------------------------------------------------------- */
   int I = 0;
   for (auto J = 0; J < XYZ.size(); J++) {
      auto net = new Network;
      int count = 0;
      for (const auto &xyz : XYZ) {
         new networkPoint(net, net, xyz);
         if (count++ > J)
            break;
      }
      net->setGeometricProperties();

      Buckets<networkPoint *> buckets(net->getBounds(), net->getScale() / 10);
      std::cout << "timer : " << timer() << std::endl;
      buckets.add(net->getPoints());
      std::cout << "all_stored_objects.size() = " << buckets.all_stored_objects.size() << std::endl;
      std::cout << "timer : " << timer() << std::endl;
      int ii = 0;
      for (auto i = 0; i < buckets.buckets.size(); ++i)
         for (auto j = 0; j < buckets.buckets[i].size(); ++j)
            for (auto k = 0; k < buckets.buckets[i][j].size(); ++k) {
               {
                  auto list = buckets.buckets[i][j][k];
                  // vtkPolygonWrite("./output/bucket" + std::to_string(ii++) +
                  // ".vtp", ToX(list));
               }
               {
                  auto list = buckets.get(buckets.itox({i, j, k}), 1);
                  // vtkPolygonWrite("./output/bucket" + std::to_string(ii++) +
                  // ".vtp", ToX(list));
               }
               {
                  auto list = buckets.get(buckets.itox({i, j, k}), 2);
                  // vtkPolygonWrite("./output/bucket" + std::to_string(ii++) +
                  // ".vtp", ToX(list));
               }
               {
                  auto list = buckets.get(buckets.itox({i, j, k}), 3);
                  // vtkPolygonWrite("./output/bucket" + std::to_string(ii++) +
                  // ".vtp", ToX(list));
               }
            }
      //
      std::cout << "timer : " << timer() << std::endl;
      //
      std::unordered_set<T_PP> lines;
      networkPoint *p_nearest = nullptr;
      for (const auto &p : buckets.getAll()) {
         p_nearest = nullptr;
         buckets.getNearest(p, {0, 10}, p_nearest);
         if (p_nearest) Add(lines, {p, p_nearest});
      }

      std::cout << "timer : " << timer() << std::endl;
      double search_radius = 3.;
      double max_radius_of_sphere = 1;
      /* -------------------------------------------------------------------------- */
      auto apply = [&](const T_3P &abc, const Tddd &X, const Tii &depth_range) {
         bool found = false;
         auto [a, b, c] = abc;
         if (a != b && a != c && b != c) {
            auto s = net->getTetras().size();
            Tddd center;
            double r;
            buckets.apply([&](const auto &d) {
               if (found)
                  return false /*means stop loops*/;
               else if (d != a && d != b && d != c) {
                  center = circumcenter(a, b, c, d);
                  r = Norm(ToX(a) - center);
                  if (isFinite(r) && r < max_radius_of_sphere && !buckets.anyInside(center, r, T_4P{a, b, c, d})) {
                     genTetra(net, a, b, c, d);
                     if (s != net->getTetras().size()) {
                        found = true;
                        return false;
                     }
                  }
               }
               return true;
            },
                          X, depth_range);
         }
         return found;
      };

      auto apply2 = [&](const T_2P &ab, const Tddd &X, const Tii &depth_range) {
         bool found = false;
         auto [a, b] = ab;
         if (a != b) {
            buckets.apply([&](const auto &c) {
               if (found)
                  return false;
               else if (a != c && b != c)
                  if (apply({a, b, c}, (ToX(a) + ToX(b) + ToX(c)) / 3., depth_range)) {
                     found = true;
                     return false;
                  }
               return true;
            },
                          X, depth_range);
         }
         return found;
      };
      /* -------------------------------------------------------------------------- */
      //  ２点から四面体を生成
      /**
       * cを使って四面体を生成するかどうか，高速にチェックするアルゴリズムが必要．
       */
      std::cout << "timer : " << timer() << std::endl;
      for (const auto &[a, b] : lines) {
         apply2({a, b}, 0.5 * (ToX(a) + ToX(b)), {0, 5});
         //  auto tmp = net->getPoints();  // buckets.getObjectsFlattened(0.5 * (ToX(a) + ToX(b)), Norm(ToX(a) - ToX(b)));
         //  for (const auto &c : buckets.getObjectsFlattened(0.5 * (ToX(a) + ToX(b)), Norm(ToX(a) - ToX(b))))
         //     if (apply({a, b, c}, (ToX(a) + ToX(b) + ToX(c)) / 3., {0, 1}))
         //        break;
      }
      std::cout << "２点から四面体を生成 : " << timer() << std::endl;
      // 既存の面の３点から四面体を生成
      for (auto k = 0; k < 10; k++) {
         for (const auto &f : net->getFaces()) {
            auto [t0, t1] = f->Tetras;
            if (!t1) {
               apply(f->getPoints(), f->X, {0, 8});
            }
         }
      }
      std::cout << "既存の面の３点から四面体を生成 : " << timer() << std::endl;
      /* -------------------------------------------------------------------------- */
      vtkPolygonWriter<networkPoint *> vtp;
      for (const auto &f : net->getFaces()) {
         auto [p0, p1, p2] = f->getPoints();
         vtp.add({p0, p1, p2});
         vtp.addPolygon({p0, p1, p2});
      }

      std::unordered_map<netPp, double> data;
      for (const auto &p : net->getPoints()) {
         auto [x, y, z] = ToX(p);
         auto v = y * x;
         data[p] = v / std::abs(v) * Norm(ToX(p));
      }

      if (J % 10 == 0) {
         I++;
         std::cout << "net points size:" << net->getPoints().size() << std::endl;
         std::cout << "net faces size:" << net->getFaces().size() << std::endl;
         std::cout << "net tetras size:" << net->getTetras().size() << std::endl;
         vtkPolygonWrite("./output/points" + std::to_string(I) + ".vtp", net->getPoints(), data);
         /* -------------------------------------------------------------------------- */
         net->setGeometricProperties();
         std::unordered_map<netPp, double> data_voroni;
         auto voronoi = new Network;
         std::vector<std::tuple<netPp, netPp>> voronoi_lines;
         for (const auto &t : net->getTetras()) {
            for_each(t->Faces, [&](const auto &f) {
               auto [t0, t1] = f->Tetras;
               if (t0 && t1) {
                  auto p0 = new networkPoint(voronoi, voronoi, t0->circumcenter);
                  auto p1 = new networkPoint(voronoi, voronoi, t1->circumcenter);
                  voronoi_lines.push_back({p0, p1});
                  data_voroni[p0] = std::get<1>(ToX(p0));
                  data_voroni[p1] = std::get<1>(ToX(p1));
               }
            });
         }
         std::cout << "timer : " << timer() << std::endl;
         vtkPolygonWrite("./output/voronoiPoints" + std::to_string(I) + ".vtp", voronoi->getPoints(), data_voroni);
         vtkPolygonWrite("./output/voronoiLines" + std::to_string(I) + ".vtp", voronoi_lines, data_voroni);
         std::ofstream ofs("./output/tetras" + std::to_string(I) + ".vtp");
         vtp.write(ofs);
         ofs.close();
         std::cout << Red << "I:" << I << colorReset << std::endl;
      }
   }
};
#elif defined(check_voronoi)
int main(int argc, char **argv) {
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
   auto X = Subdivide(
       xb0, xb1,
       (int)std::round((xb1 - xb0) / particle_spacing));  // 粒子のX座標
   auto Y = Subdivide(
       yb0, yb1,
       (int)std::round((yb1 - yb0) / particle_spacing));  // 粒子のY座標
   auto Z = Subdivide(
       zb0, zb1,
       (int)std::round((zb1 - zb0) / particle_spacing));  // 粒子のZ座標
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
   for (const auto &[p, n] : p_neighbors)
      for (const auto &q : n)
         lines.push_back(T2Tddd{p->getXtuple(), q->getXtuple()});
   vtkPolygonWrite("./output/cubeLines.vtp", lines);
};
#elif defined(check_points)
int main(int argc, char **argv) {
   Timer timer;
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
#elif defined(check_triangles)
int main(int argc, char **argv) {
   Timer timer;
   /* ----------------------- 引数の読み込み ---------------------- */
   std::string name{argv[1]};  //"/Users/tomoaki/Dropbox/markdown/cpp/obj/cow.obj";
   int maxDepth = std::atoi(argv[2]);
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