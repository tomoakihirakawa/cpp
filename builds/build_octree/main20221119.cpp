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
#if defined(check_SPH3)
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
      setVectorsToTriangle(tree);
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
   /* ------------------------------------------------------------------------
      1) Fluid作成とRigidBodyを作成
      2)
   -------------------------------------------------------------------------- */
   {
      // auto polyFluid = new Network("./OC5water.obj", "water");
      auto polyFluid = new Network("../../obj/2022Arai/simple_case/water.obj", "water");
      {
         auto ofs0 = std::ofstream("./output/polyFluid.vtp");
         vtkPolygonWrite(ofs0, ToX(polyFluid->getFaces()));
         ofs0.close();
         polyFluid->genOctreeOfFaces({minDepth, maxDepth}, 1);
      }
      // auto polyWaveTank = new Network("./OC5tank2.obj", "tank");
      auto polyWaveTank = new Network("../../obj/2022Arai/simple_case/tank.obj", "tank");
      {
         auto ofs1 = std::ofstream("./output/polyWaveTank.vtp");
         vtkPolygonWrite(ofs1, ToX(polyWaveTank->getFaces()));
         ofs1.close();
         polyWaveTank->genOctreeOfFaces({minDepth, maxDepth}, 1);
         polyWaveTank->octreeOfFaces->setNeighbors();
         setVectorsToTriangle(polyWaveTank->octreeOfFaces);
         auto ofs3 = std::ofstream("./output/treeWaveTank.vtp");
         vtkPolygonWrite(ofs3, toCubeFaces(polyWaveTank->octreeOfFaces->getAllDeepestInside()));
         ofs3.close();
      }
      // auto polyWall = new Network("./OC5wavemaker.obj", "wall");
      auto polyWall = new Network("../../obj/2022Arai/simple_case/tank.obj", "wall");
      {
         auto ofs2 = std::ofstream("./output/polyWall.vtp");
         vtkPolygonWrite(ofs2, ToX(polyWall->getFaces()));
         ofs2.close();
         polyWall->genOctreeOfFaces({minDepth, maxDepth}, 1);
         polyWall->octreeOfFaces->setNeighbors();
         setVectorsToTriangle(polyWall->octreeOfFaces);
         auto ofs3 = std::ofstream("./output/treeWall.vtp");
         vtkPolygonWrite(ofs3, toCubeFaces(polyWall->octreeOfFaces->getAllDeepestInside()));
         ofs3.close();
      }
      //
      auto Fluid = new Network;
      auto WaveTank = new Network;
      auto Wall = new Network;
      //
      std::vector<std::tuple<Network *, Network *>> RigidBodies = {
          std::tuple<Network *, Network *>{WaveTank, polyWaveTank}
          // , std::tuple<Network *, Network *>{Wall, polyWall}
      };
      //
      int n = 50;
      //
      // Tdd rangeX = {-0.5, 19.5}, rangeY = {-1.5, 1.5}, rangeZ = {-0.5, 2};
      // Tdd rangeX = {-0.5, 2.5}, rangeY = {-0.5, 1.5}, rangeZ = {-0.5, 1};
      CoordinateBounds range(polyFluid->getBounds());
      range += CoordinateBounds(polyWaveTank->getBounds());
      range += CoordinateBounds(polyWall->getBounds());
      auto [rangeX, rangeY, rangeZ] = range.bounds;
      auto vecX = Subdivide(rangeX, std::round(std::get<1>(rangeX) - std::get<0>(rangeX)) * (n - 1));
      auto vecY = Subdivide(rangeY, std::round(std::get<1>(rangeY) - std::get<0>(rangeY)) * (n - 1));
      auto vecZ = Subdivide(rangeZ, std::round(std::get<1>(rangeZ) - std::get<0>(rangeZ)) * (n - 1));
      //
      double particle_spacing = 1. / n;
      double volume = std::pow(particle_spacing, 3);
      double C_SML = 3.;
      double max_dt = 0.001;
      double initial_surface_z_position = 0.5;
      //
      for (const auto &x : vecX)
         for (const auto &y : vecY)
            for (const auto &z : vecZ) {
               Tddd xyz = {x, y, z};
               auto [isInsideOfTank, cell0] = polyWaveTank->isInside_MethodOctree(xyz);
               auto [isInsideOfWall, cell1] = polyWall->isInside_MethodOctree(xyz);
               auto [isInsideOfFluid, cell2] = polyFluid->isInside_MethodOctree(xyz);
               if (isInsideOfTank) {
                  auto p = new networkPoint(WaveTank, xyz);
                  p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
                  p->radius_SPH = C_SML * particle_spacing;
                  p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
               } else if (isInsideOfWall) {
                  auto p = new networkPoint(Wall, xyz);
                  p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
                  p->radius_SPH = C_SML * particle_spacing;
                  p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
               } else if (isInsideOfFluid) {
                  auto p = new networkPoint(Fluid, xyz);
                  p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
                  p->radius_SPH = C_SML * particle_spacing;
                  p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
               }
            }
      // for (const auto &cell : polyWaveTank->octreeOfFaces->getAllDeepestInside()) {
      //    auto p = new networkPoint(Wall, cell->center);
      //    p->setDensityVolume(_WATER_DENSITY_, cell->getVolume());  // 質量(mass)は関数内部で自動で決められる
      //    p->radius_SPH = C_SML * particle_spacing;
      //    p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
      // }
      std::cout << "timer : " << timer() << std::endl;
      for (const auto &[obj, poly] : RigidBodies)
         for (const auto &p : obj->getPoints())
            p->normal_SPH = poly->interpolateVector(p->X);
      //
      std::cout << "timer : " << timer() << ", networkPointの生成" << std::endl;
      std::cout << "Fluid->getPoints().size()=" << Fluid->getPoints().size() << std::endl;
      std::cout << "WaveTank->getPoints().size()=" << WaveTank->getPoints().size() << std::endl;
      std::cout << "Wall->getPoints().size()=" << Wall->getPoints().size() << std::endl;
      {
         auto ofs = std::ofstream("./output/Fluid.vtp");
         vtkPolygonWrite(ofs, ToX(Fluid->getPoints()));
         ofs.close();
      }
      {
         auto ofs = std::ofstream("./output/WaveTank.vtp");
         vtkPolygonWrite(ofs, ToX(WaveTank->getPoints()));
         ofs.close();
      }
      {
         auto ofs = std::ofstream("./output/Wall.vtp");
         vtkPolygonWrite(ofs, ToX(Wall->getPoints()));
         ofs.close();
      }
      std::cout << "timer : " << timer() << std::endl;
      Fluid->setGeometricProperties();
      WaveTank->setGeometricProperties();
      Wall->setGeometricProperties();
      /* -------------------------------------------------------------------------- */
      double real_time = 0.;
      PVDWriter pvdFluidSPH("./output/FluidSPH.pvd");
      PVDWriter pvdWaveTankSPH("./output/WaveTankSPH.pvd");
      PVDWriter pvdWallSPH("./output/Wall.pvd");
      int count = 0, i = 0;
      //
      for (auto count = 0; count < 5000; ++count) {
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

         // setContactPoints(Fluid, RigidBodies, real_time, C_SML, particle_spacing);
         // 出力
         if (count % 5 == 0) {
            {
               vtkPolygonWriter<networkPoint *> vtp;
               vtp.add(ToVector(Fluid->getPoints()));
               setData(vtp, Fluid);
               std::ofstream ofs("./output/FluidSPH" + std::to_string(i) + ".vtp");
               vtp.write(ofs);
               ofs.close();
            }
            // {
            //    vtkPolygonWriter<networkPoint *> vtp;
            //    vtp.add(ToVector(WaveTank->getPoints()));
            //    setData(vtp, WaveTank);
            //    std::ofstream ofs("./output/WaveTankSPH" + std::to_string(i) + ".vtp");
            //    vtp.write(ofs);
            //    ofs.close();
            // }
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
         developByEISPH(Fluid, RigidBodies, real_time, C_SML, particle_spacing);
         std::cout << "real_time = " << real_time << std::endl;
      }
   }
};
#elif defined(check_voronoi0)
int main(int arg, char **argv) {
   //! NOTE
   // 四面体の許容サイズは適切か確認（max_radius_of_sphere）
   Timer timer;
   std::string name(argv[1]);
   int minDepth(std::atoi(argv[2]));
   int maxDepth(std::atoi(argv[3]));
   double spacing(std::atof(argv[4]));
   std::cout << Grid({"minDepth", "maxDepth", "spacing"}, 10) << std::endl;
   std::cout << Grid({minDepth, maxDepth, spacing}, 10) << std::endl;
   /* ------------------------------make------------------------ */
   // auto [xb0, xb1] = Tdd{-5., 5.};
   // auto [yb0, yb1] = Tdd{-5., 5.};
   // auto [zb0, zb1] = Tdd{-5., 5.};
   // auto X = Subdivide(xb0, xb1, (int)std::round((xb1 - xb0) / particle_spacing));  // 粒子のX座標
   // auto Y = Subdivide(yb0, yb1, (int)std::round((yb1 - yb0) / particle_spacing));  // 粒子のY座標
   // auto Z = Subdivide(zb0, zb1, (int)std::round((zb1 - zb0) / particle_spacing));  // 粒子のZ座標
   /* -------------------------------------------------------------------------- */
   int n = 20;
   auto object = new Network(name, "object");
   object->genOctreeOfFaces({0, 5}, 1);
   auto ofs = std::ofstream("./output/polygon.vtp");
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
   Tdd r = {-0.0001, 0.0001};
   std::vector<Tddd> XYZ;
   for (const auto &x : vecX)
      for (const auto &y : vecY)
         for (const auto &z : vecZ) {
            Tddd xyz = {x + RandomReal(r), y + RandomReal(r), z + RandomReal(r)};
            if (std::get<0>(object->isInside_MethodOctree(xyz)))
               XYZ.push_back(xyz);
         }
   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   int I = 0;
   for (auto J = XYZ.size() - 1; J < XYZ.size(); J++) {
      auto net = new Network;
      int count = 0;
      for (const auto &xyz : XYZ) {
         new networkPoint(net, xyz);
         if (count++ > J)
            break;
      }
      net->setGeometricProperties();
      //

      // net->genOctreeOfPoints({minDepth, maxDepth}, 1);
      // auto ofs = std::ofstream("./output/octreeOfPoints" + std::to_string(I) + ".vtp");
      // vtkPolygonWrite(ofs, toCubeFaces(net->octreeOfPoints->getAllDeepestInside()));
      // ofs.close();
      /* --------------------------------- バケツの準備 --------------------------------- */
      Buckets<networkPoint *> buckets(object->getBounds(), particle_spacing);
      Buckets<networkFace *> face_buckets(object->getBounds(), particle_spacing);
      Buckets<networkTetra *> tetra_buckets(object->getBounds(), particle_spacing);
      //
      std::cout << "timer : " << timer() << std::endl;
      buckets.add(net->getPoints());
      std::cout << "all_stored_objects.size() = " << buckets.all_stored_objects.size() << std::endl;
      std::cout << "timer : " << timer() << std::endl;
      /* -------------------------------------------------------------------------- */
      // b* ------------------ ボロノイ図，テトラが問題なくできているか確認 ------------------- */
      auto check = [&]() {
         for (const auto &tet : net->getTetras()) {
            bool error_found = buckets.any_of(tet->circumcenter, 5 * tet->circumradius,
                                              /*少し半径を小さくし，ほぼ半径上にあるものはスルーする*/
                                              [&](const auto &p) { return !MemberQ(tet->Points, p) && tet->circumradius - 1E-10 >= Norm(ToX(p) - tet->circumcenter); });
            if (error_found)
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error found");
         }
      };
      /* -------------------------------------------------------------------------- */
      auto conditions = [&](const T_4P &abcd) {
         if (!DuplicateFreeQ(abcd)) {
            return false;
         }
         auto [a, b, c, d] = abcd;
         /* ---------------------------------------------------------- */
         // 3点でなく2点だけを共有し，法線方向が同じ面がある場合，四面体は作れない
         // そのような面が一つでもあれば．
         auto isFaceOverlap = [&](const auto &f) {
            std::vector<networkPoint *> sharing, not_sharing;
            for_each(abcd, [&](const auto &a) {
               if (MemberQ(f->getPoints(), a))
                  sharing.emplace_back(a);
               else
                  not_sharing.emplace_back(a);
            });
            if (sharing.size() == 2) {
               for (const auto &P : not_sharing) {
                  auto n = TriangleNormal(ToX(sharing[0]), ToX(sharing[1]), ToX(P));
                  if (Between(std::abs(Dot(f->normal, n)), {1. - 1E-10, 1. + 1E-10})) {
                     // 同じ平面
                     auto mean = (ToX(sharing[0]) + ToX(sharing[1])) / 2.;
                     if (Dot(ToX(P) - mean, f->center - mean) > 0.)
                        return true;  // 重なっている
                  }
               }
               // auto [t0, t1] = f->Tetras;
               // if (t0) {
               //    if (IntersectQ(ToX(abcd), ToX(t0->Points)))
               //       return true;
               // }
               // if (t1) {
               //    if (IntersectQ(ToX(abcd), ToX(t1->Points)))
               //       return true;
               // }
            }
            return false;  // 問題なし
         };
         /* ---------------------------------------------------------- */
         bool nogoood = any_of(abcd, [&](const auto &a) {
            return std::any_of(a->Faces.begin(), a->Faces.end(), isFaceOverlap);
         });
         if (nogoood)
            return false;
         /* ---------------------------------------------------------- */
         Tetrahedron tet(ToX(abcd));
         // if (Max(tet.solidangles) >= M_PI / 2.)
         //    return false;
         // else
         {
            auto [isinside, f] = object->isInside_MethodOctree(tet.circumcenter);
            if (!isinside) {
               return false;
            } else
               return buckets.none_of(tet.circumcenter,
                                      4 * tet.circumradius,
                                      [&](const auto &p) {
                                         /*少し半径を小さくし，ほぼ半径上にあるものはスルーする*/
                                         if (!MemberQ(abcd, p) && tet.circumradius >= Norm(ToX(p) - tet.circumcenter))
                                            return true;
                                         //   for (const auto &T : p->getTetras())
                                         //      if (T.circumradius >= Norm(ToX(p) - T.circumcenter)) {
                                         //         return true;
                                         //      }
                                         return false;
                                      });
         }
      };
      // b% -------------------------------------------------------------------------- */
      Print("最短距離にある２点を繋ぐ", Red);
      // randomに分散させると繋ぐ２点がとても少なくなる
      std::unordered_map<networkPoint *, std::vector<networkPoint *>> p_neighbors;
      for (const auto &a : net->getPoints()) {
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
         std::ofstream ofs("./output/first_lines.vtp");
         vtkPolygonWrite(ofs, net->getLines());
      }
      check();
      // b% -------------------------------------------------------------------------- */
      Print("繋がった点は，四面体の候補である", Red);

      auto subtetras = [](const Tetrahedron &tet) {
         auto [a, b, c, d] = tet.verticies;
         Tetrahedron tet0({a, a / 2. + b / 2., a / 2. + c / 2., a / 2. + d / 2.});
         Tetrahedron tet1({b, b / 2. + a / 2., b / 2. + c / 2., b / 2. + d / 2.});
         Tetrahedron tet2({c, c / 2. + a / 2., c / 2. + b / 2., c / 2. + d / 2.});
         Tetrahedron tet3({d, d / 2. + a / 2., d / 2. + b / 2., d / 2. + c / 2.});
         return std::vector<Tetrahedron>{tet, tet0, tet1, tet2, tet3};
      };

      auto edges = [](const Tetrahedron &tet) {
         auto [a, b, c, d] = tet.verticies;
         return std::vector<T2Tddd>{{a, b}, {a, c}, {a, d}, {b, c}, {b, d}, {c, d}};
      };
      /* -------------------------------------------------------------------------- */
      auto gen = [&](const T_4P &abcd) {
         if (!DuplicateFreeQ(abcd))
            return std::tuple<bool, networkTetra *>{false, nullptr};
         Tetrahedron tet(ToX(abcd));
         auto [isinside, f] = object->isInside_MethodOctree(tet.circumcenter);
         if (isinside) {
            int count = 0;
            int on_surface_count = 0;
            auto [a, b, c, d] = abcd;

            std::vector<Tetrahedron> creating_tets = subtetras(tet);
            std::vector<T2Tddd> creating_edges = edges(tet);

            if (tet.inradius > 1E-4 &&
                buckets.none_of(tet.circumcenter, 3 * depth,
                                [&](const auto &p) {
                                   if (!MemberQ(abcd, p) && Norm(ToX(p) - tet.circumcenter) - (tet.circumradius - 1E-10) <= 0.)
                                      return true;
                                   else
                                      return false;
                                }) &&
                tetra_buckets.none_of(  // tet.incenter, 3 * depth,
                    [&](const auto &loop_TET) {
                       auto loop_edges = edges(*loop_TET);
                       auto loop_subtetras = subtetras(*loop_TET);
                       if (std::any_of(std::begin(creating_tets), std::end(creating_tets), [&](const auto &t) {
                              if (std::any_of(std::begin(loop_edges), std::end(loop_edges),
                                              [&](const auto &ab) { return IntersectQ(t.incenter, t.inradius + 2E-5, ab); }))
                                 return true;
                              if (std::any_of(std::begin(loop_subtetras), std::end(loop_subtetras),
                                              [&](const auto &T) { return IntersectQ(t.incenter, t.inradius, T.incenter, T.inradius + 2E-5); }))
                                 return true;
                              return false;
                           }))
                          return true;
                       if (std::any_of(std::begin(loop_subtetras), std::end(loop_subtetras), [&](const auto &t) {
                              if (std::any_of(std::begin(creating_edges), std::end(creating_edges),
                                              [&](const auto &ab) { return IntersectQ(t.incenter, t.inradius + 2E-5, ab); }))
                                 return true;
                              if (std::any_of(std::begin(creating_tets), std::end(creating_tets),
                                              [&](const auto &T) { return IntersectQ(t.incenter, t.inradius, T.incenter, T.inradius + 2E-5); }))
                                 return true;
                              return false;
                           }))
                          return true;
                       return false;
                    })) {
               //
               auto ret = genTetra(net, a, b, c, d);
               tetra_buckets.add(std::get<1>(ret)->incenter, std::get<1>(ret));
               return ret;
            }
         }
         return std::tuple<bool, networkTetra *>{false, nullptr};
      };
      /* -------------------------------------------------------------------------- */
      int counter = 0;
      for (const auto &a : net->getPoints()) {
         auto vec = p_neighbors[a];
         for (auto it = vec.begin(); it != vec.end(); ++it) {
            for (auto jt = std::next(it, 1); jt != vec.end(); ++jt) {
               for (auto kt = std::next(jt, 1); kt != vec.end(); ++kt)
                  gen(T_4P{a, *it, *jt, *kt});
            }
         }
      }
      // b% -------------------------------------------------------------------------- */
      Print("形成された面の内，テトラを持たない面の法線方向において，最も近い点を候補とし，新たにテトラを作成できないか調べる");
      for (auto i = 0; i < 30; i++) {
         int count = 0, generated = 0, good_condition = 0;
         for (const auto &f : net->getFaces()) {
            bool found = false;
            auto [t0, t1] = f->Tetras;
            if (t0 && !t1) {
               // 候補点を探す
               std::vector<networkPoint *> candidates;
               auto [a, b, c] = f->getPoints();
               Tddd dir = Normalize(f->center - t0->centroid);
               buckets.apply(f->center, 5. * depth, [&](const auto &d) {
                  if (DuplicateFreeQ(T_4P{a, b, c, d}) && Dot(d->X - f->center, dir) > 0.) {
                     candidates.emplace_back(d);
                     if (candidates.size() > 5) {
                        sortByDistance(candidates, f->center);
                        candidates.pop_back();
                     }
                  }
               });
               for (const auto &d : candidates) {
                  gen(T_4P{a, b, c, d});
               }
            }
         }
         std::cout << Magenta << i << Green << ", Elapsed time : " << timer() << colorReset << std::endl;
         {
            std::unordered_map<std::shared_ptr<Tddd>, double> data;
            vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
            for (const auto &tet : net->getTetras()) {
               std::shared_ptr<Tddd> x(new Tddd(tet->incenter));
               //
               for (const auto &t : subtetras(*tet)) {
                  std::shared_ptr<Tddd> x(new Tddd(t.incenter));
                  data[x] = t.inradius;
                  vtp.add(x);
               }
               vtp.addPointData("inradius", data);
               {
                  std::ofstream ofs("./output/tetra_insphere" + std::to_string(i) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
            }
            /* ------------------------------------------- */
            std::ofstream ofs("./output/first_tetra" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, net->getTetras());
            ofs.close();
         }
      }
      break;
      /* -------------------------------------------------------------------------- */
      // for (const auto &a : net->getPoints()) {
      //    std::unordered_set<networkPoint *> vec;
      //    buckets.apply(a->X, 3 * depth, [&](const auto &b) {
      //       if (a != b)
      //          vec.emplace(b);
      //    });
      //    for (auto it = vec.begin(); it != vec.end(); ++it) {
      //       for (auto jt = std::next(it, 1); jt != vec.end(); ++jt) {
      //          for (auto kt = std::next(jt, 1); kt != vec.end(); ++kt)
      //             gen(T_4P{a, *it, *jt, *kt});
      //       }
      //    }
      // }
      /* -------------------------------------------------------------------------- */
      // int counter = 0;
      // for (const auto &a : net->getPoints()) {
      //    auto vec = p_neighbors[a];
      //    for (auto it = vec.begin(); it != vec.end(); ++it) {
      //       for (auto jt = std::next(it, 1); jt != vec.end(); ++jt) {
      //          for (auto kt = std::next(jt, 1); kt != vec.end(); ++kt)
      //             gen(T_4P{a, *it, *jt, *kt});
      //       }
      //    }
      // }
      /* -------------------------------------------------------------------------- */
      {
         int counter = 0;
         std::unordered_map<std::shared_ptr<Tddd>, double> data;
         vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
         std::unordered_set<networkTetra *> tetras;
         for (const auto &tet : net->getTetras()) {
            //
            std::shared_ptr<Tddd> x(new Tddd(tet->incenter));
            // std::cout << "incenter = " << tet->incenter << ", inradius = " << tet->inradius << std::endl;
            //
            for (const auto &t : subtetras(*tet)) {
               std::shared_ptr<Tddd> x(new Tddd(t.incenter));
               data[x] = t.inradius;
               vtp.add(x);
            }
            //
            vtp.addPointData("inradius", data);
            {
               std::ofstream ofs("./output/tetra_insphere" + std::to_string(counter) + ".vtp");
               vtp.write(ofs);
               ofs.close();
            }
            //
            {
               std::ofstream ofs("./output/first_tetra" + std::to_string(counter) + ".vtp");
               tetras.emplace(tet);
               vtkPolygonWrite(ofs, tetras);
               ofs.close();
            }
            counter++;
         }
      }

      std::cout << Green << "Elapsed time : " << timer() << colorReset << std::endl;
      {
         std::ofstream ofs("./output/first_tetras.vtp");
         vtkPolygonWrite(ofs, net->getTetras());
      }
      std::cin.ignore();
      check();
      // b% -------------------------------------------------------------------------- */
      // Print("面の探査");
      // std::vector<T_4P> abcd;
      // for (const auto &a : net->getPoints()) {
      //    for (const auto &l0 : a->getLines()) {
      //       auto b = (*l0)(a);
      //       for (const auto &l1 : b->getLines())
      //          if (l0 != l1) {
      //             auto c = (*l1)(b);
      //             for (const auto &l2 : c->getLines())
      //                if (l1 != l2) {
      //                   {
      //                      auto d = (*l2)(c);
      //                      if (DuplicateFreeQ(T_4P{a, b, c, d})) {  // triangle found
      //                         abcd.push_back({a, b, c, d});
      //                      }
      //                   }
      //                }
      //          }
      //    }
      // }
      // for (const auto &[a, b, c, d] : abcd)
      //    genTetra(net, a, b, c, d);
      // std::cout << Green << "Elapsed time : " << timer() << colorReset << std::endl;
      // b% -------------------------------------------------------------------------- */
      // Print("予め繋いだline(２点)から四面体を生成");
      // //! 適当な３点目を選択し，さらに適当な４点目を選ぶ．そして，その四面体の内部に点があるかどうかをチェックする（すでに造られた四面体の外心が内部にあるかチェックすべきだろう）．
      // for (const auto &l : net->getLines()) {
      //    auto [a, b] = l->getPoints();
      //    buckets.any_of((a->X + b->X) / 2., 2 * depth, [&](const auto &c) {
      //       if (c != a && c != b)
      //          return generateTetra({a, b, c}, (a->X + b->X + c->X) / 3., 2 * depth);
      //       return false;
      //    });
      // }
      // std::cout << Green << "Elapsed time : " << timer() << colorReset << std::endl;
      // b% -------------------------------------------------------------------------- */
      Print("形成された面の内，テトラを持たない面の法線方向において，最も近い点を候補とし，新たにテトラを作成できないか調べる");
      for (auto i = 0; i < 30; i++) {
         int count = 0, generated = 0, good_condition = 0;
         for (const auto &f : net->getFaces()) {
            bool found = false;
            auto [t0, t1] = f->Tetras;
            if (t0 && !t1 || !t0 && t1) {
               // 候補点を探す
               std::vector<networkPoint *> candidates;
               auto [a, b, c] = f->getPoints();
               Tddd dir = Normalize(f->center - t0->centroid);
               buckets.apply(f->center, 5. * depth, [&](const auto &d) {
                  if (DuplicateFreeQ(T_4P{a, b, c, d}) && Dot(d->X - f->center, dir) > 0.) {
                     candidates.emplace_back(d);
                     if (candidates.size() > 5) {
                        sortByDistance(candidates, f->center);
                        candidates.pop_back();
                     }
                  }
               });
               for (const auto &d : candidates) {
                  gen(T_4P{a, b, c, d});
               }
            }
         }
         //
         // for (auto k = 0; k < 10; k++) {
         //    for (const auto &f : net->getFaces()) {
         //       auto [t0, t1] = f->Tetras;
         //       if (!t1) {
         //          auto abc = f->getPoints();
         //          auto X = f->X;
         //          generateTetra(abc, X, depth);
         //       }
         //    }
         // }
         std::cout << Green << "Elapsed time : " << timer() << colorReset << std::endl;
         {
            std::ofstream ofs("./output/tetra" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, net->getTetras());
         }
         // check();
      }
      // b! ---------------------------------------------------------------------------- */
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

      if (J % 10 == 0 || J == XYZ.size() - 1) {
         I++;
         std::cout << "net points size:" << net->getPoints().size() << std::endl;
         std::cout << "net lines size:" << net->getLines().size() << std::endl;
         std::cout << "net faces size:" << net->getFaces().size() << std::endl;
         std::cout << "net tetras size:" << net->getTetras().size() << std::endl;
         {
            std::ofstream ofs("./output/points" + std::to_string(I) + ".vtp");
            vtkPolygonWrite(ofs, net->getPoints(), data);
         }
         {
            std::ofstream ofs("./output/faces" + std::to_string(I) + ".vtp");
            vtkPolygonWrite(ofs, net->getFaces(), data);
         }
         /* -------------------------------------------------------------------------- */
         net->setGeometricProperties();
         std::unordered_map<netPp, double> data_voroni;
         auto voronoi = new Network;
         std::vector<std::tuple<netPp, netPp>> voronoi_lines;
         for (const auto &t : net->getTetras()) {
            for_each(t->Faces, [&](const auto &f) {
               auto [t0, t1] = f->Tetras;
               if (t0 && t1) {
                  auto p0 = new networkPoint(voronoi, t0->circumcenter);
                  auto p1 = new networkPoint(voronoi, t1->circumcenter);
                  voronoi_lines.push_back({p0, p1});
                  data_voroni[p0] = std::get<1>(ToX(p0));
                  data_voroni[p1] = std::get<1>(ToX(p1));
               }
            });
         }

         // int i = 0;
         // for (const auto &t : net->getTetras()) {
         //    std::ofstream ofs("./output/tetras" + std::to_string(i++) + ".vtp");
         //    vtkPolygonWrite(ofs, t->Faces);
         //    ofs.close();
         // }
         {
            std::ofstream ofs("./output/tetrasAll.vtp");
            vtkPolygonWrite(ofs, net->getTetras());
            ofs.close();
         }
         std::cout << "timer : " << timer() << std::endl;
         {
            std::ofstream ofs("./output/voronoiPoints" + std::to_string(I) + ".vtp");
            vtkPolygonWrite(ofs, voronoi->getPoints(), data_voroni);
         }
         {
            std::ofstream ofs("./output/voronoiLines" + std::to_string(I) + ".vtp");
            vtkPolygonWrite(ofs, voronoi_lines, data_voroni);
         }
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