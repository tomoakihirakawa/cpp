// #define _debugging_

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
int main() {
   Timer timer;
   int minDepth = 4, maxDepth = 5;
   std::cout << "minDepth = " << minDepth << std::endl;
   std::cout << "maxDepth = " << maxDepth << std::endl;
   std::cout << "timer : " << timer() << std::endl;
   /* ------------------------------------------------------------------------
      1) Fluid作成とRigidBodyを作成
      2)
   -------------------------------------------------------------------------- */
   auto polyFluid = new Network("./OC5water.obj", "water");
   // auto polyFluid = new Network("../../obj/2022Arai/simple_case/water.obj", "water");
   // auto polyFluid = new Network("../../obj/2022Arai/static_pressure/water.obj", "water");
   {
      auto ofs0 = std::ofstream("./output/polyFluid.vtp");
      vtkPolygonWrite(ofs0, ToX(polyFluid->getFaces()));
      ofs0.close();
      polyFluid->genOctreeOfFaces({minDepth, maxDepth}, 1);
   }
   auto polyWaveTank = new Network("./OC5tank.obj", "tank");
   // auto polyWaveTank = new Network("../../obj/2022Arai/simple_case/wall.obj", "tank");
   // auto polyWaveTank = new Network("../../obj/2022Arai/static_pressure/tank.obj", "tank");
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
   auto polyWall = new Network("./OC5tank.obj", "wall");
   // auto polyWall = new Network("../../obj/2022Arai/simple_case/wall.obj", "wall");
   // auto polyWall = new Network("../../obj/2022Arai/static_pressure/tank.obj", "wall");
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
   polyFluid->translate(Tddd{1E-5, 1E-5, 1E-5} * M_PI);
   polyWaveTank->translate(Tddd{1E-5, 1E-5, 1E-5} * M_PI);
   polyWall->translate(Tddd{1E-5, 1E-5, 1E-5} * M_PI);
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
   // Tdd rangeX = {-0.5, 19.5}, rangeY = {-1.5, 1.5}, rangeZ = {-0.5, 2};
   // Tdd rangeX = {-0.5, 2.5}, rangeY = {-0.5, 1.5}, rangeZ = {-0.5, 1};
   //
   double particle_spacing = 0.04;
   double &ps = particle_spacing;
   double volume = std::pow(particle_spacing, 3);
   //
   CoordinateBounds range(polyFluid->getBounds());
   range += CoordinateBounds(polyWaveTank->getBounds());
   range += CoordinateBounds(polyWall->getBounds());
   auto [rangeX, rangeY, rangeZ] = range.bounds;
   //
   auto [rX0, rX1] = rangeX;
   auto [rY0, rY1] = rangeY;
   auto [rZ0, rZ1] = rangeZ;
   auto vecX = Subdivide({ps * (int)(rX0 / ps), ps * (int)(rX1 / ps)}, (int)(rX1 / ps) - (int)(rX0 / ps));
   auto vecY = Subdivide({ps * (int)(rY0 / ps), ps * (int)(rY1 / ps)}, (int)(rY1 / ps) - (int)(rY0 / ps));
   auto vecZ = Subdivide({ps * (int)(rZ0 / ps), ps * (int)(rZ1 / ps)}, (int)(rZ1 / ps) - (int)(rZ0 / ps));
   //
   double C_SML = 3.000001;
   double max_dt = 0.001;
   double initial_surface_z_position = 0.1;
   //
   double rand = RandomReal({1E-10, 1E-10});
   for (const auto &x : vecX)
      for (const auto &y : vecY)
         for (const auto &z : vecZ) {
            Tddd xyz = {x, y, z};
            // xyz += 1E-4;  // M_PI / 1000.;
            //
            auto [isInsideOfTank, cell0, f0] = polyWaveTank->isInside_MethodOctree(xyz);
            auto [isInsideOfWall, cell1, f1] = polyWall->isInside_MethodOctree(xyz);
            auto [isInsideOfFluid, cell2, f2] = polyFluid->isInside_MethodOctree(xyz);
            if (isInsideOfTank) {
               auto p = new networkPoint(WaveTank, xyz);
               p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
               p->radius_SPH = C_SML * ps;
               p->C_SML = C_SML;
               p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
            } else if (isInsideOfWall) {
               auto p = new networkPoint(Wall, xyz);
               p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
               p->radius_SPH = C_SML * ps;
               p->C_SML = C_SML;
               p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
            } else if (isInsideOfFluid) {
               auto p = new networkPoint(Fluid, xyz);
               p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
               p->radius_SPH = C_SML * ps;
               p->C_SML = C_SML;
               p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
            }
         }
   // b# -------------------------------------------------------------------------- */
   // b#                             外向きベクトルの設定                               */
   // b# -------------------------------------------------------------------------- */
   WaveTank->makeBucketPoints(particle_spacing);
   //
   std::cout << "timer : " << timer() << std::endl;
   for (const auto &[obj, poly] : RigidBodies)
      for (const auto &p : obj->getPoints()) {
         // p->normal_SPH = poly->interpolateVector(p->X);
         auto [X, f] = Nearest_(p->X, poly->getFaces());
         p->normal_SPH = X - p->X;
         p->mirroring_face = f;
      }
   //
   {
      vtkPolygonWriter<networkPoint *> vtp;
      vtp.add(WaveTank->getPoints());
      std::unordered_map<networkPoint *, Tddd> normal_SPH;
      for (const auto &p : WaveTank->getPoints())
         normal_SPH[p] = p->normal_SPH;
      vtp.addPointData("normal_SPH", normal_SPH);
      std::ofstream ofs("./output/WaveTank_check_normal.vtp");
      vtp.write(ofs);
      ofs.close();
   }
   //
   for (const auto &[obj, poly] : RigidBodies)
      for (const auto &p : obj->getPoints()) {
         auto [isinside, cell, f] = poly->isInside_MethodOctree(p->X);
         if (isinside && f) {
            // auto intsp = IntersectionSphereTriangle(p->X, 1E+10, ToX(f));
            // p->normal_SPH = intsp.X - p->X;
            auto nearest = Nearest(p->X, ToX(f));
            p->normal_SPH = nearest - p->X;
            p->mirroring_face = f;
         }
      };
   // b# -------------------------------------------------------------------------- */
   // b# -------------------------------------------------------------------------- */
   // b# -------------------------------------------------------------------------- */
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
      double mu = 0.001138;
      for (const auto &p : Fluid->getPoints())
         p->mu_SPH = mu;
      for (const auto &p : WaveTank->getPoints())
         p->mu_SPH = mu;
      for (const auto &p : Wall->getPoints())
         p->mu_SPH = mu;

      std::cout << Blue << "count = " << count << ", mu = " << mu << std::endl;

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
         if (count % 100 == 0) {
            vtkPolygonWriter<networkPoint *> vtp;
            vtp.add(ToVector(WaveTank->getPoints()));
            setDataOmitted(vtp, WaveTank);
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
         // pvdWaveTankSPH.push("./WaveTankSPH" + std::to_string(i) + ".vtp", real_time);
         // pvdWaveTankSPH.output();
         // pvdWallSPH.push("./WallSPH" + std::to_string(i) + ".vtp", real_time);
         // pvdWallSPH.output();

         i++;
      }
      developByEISPH(Fluid, RigidBodies, real_time, C_SML, particle_spacing);
      std::cout << "real_time = " << real_time << std::endl;
   }
};