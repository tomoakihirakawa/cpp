/*DOC_EXTRACT 1_0_0_SPH

# 実行方法

ファイルをダウンロードして，`build_sph`ディレクトリに移動．
⚠️上書きされるので注意．

```sh
git clone https://github.com/tomoakihirakawa/cpp.git
cd ./cpp/builds/build_sph
```

`clean`でCMake関連のファイルを削除して（ゴミがあるかもしれないので），
`cmake`で`Makefile`を生成して，`make`でコンパイルする．

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

次に，入力ファイルを生成．

```sh
python3 input_generator.py
```

例えば，`./input_files/static_pressure_PS0d0125_CSML2d4_RK1`が生成される．
入力ファイルを指定して実行．

```sh
./main ./input_files/static_pressure_PS0d0125_CSML2d4_RK1
```

*/

double delta_t;
const double too_small_total_w = 1E-3;

#include <filesystem>
#include <utility>
#define DEM
#include "Network.hpp"

std::unordered_set<Network *> _ALL_NET_;

#include "SPH.hpp"
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

JSONoutput jsonout;

int main(int argc, char **argv) {
   /* -------------------------------------------------------------------------- */
   /*                           Set up logging to file                           */
   /* -------------------------------------------------------------------------- */
   if (!initializeLogFile("./log_" + getUserName() + ".txt", argc, argv))
      return 1;
   if (!initializeLogFile("./log.txt", argc, argv))
      return 1;
   /* -------------------------------------------------------------------------- */

   if (argc <= 1)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "argv <= 1. write input json file directory!\\ex.\\$ ./main ./input");
   std::string input_directory{argv[1]};  // input directory
   input_directory += "/";
   //
   std::string id = "";
   // if (argc >= 3)
   //    id = argv[2];  // input directory
   // b! -------------------------------------------------------------------------- */
   std::cout << "input_directory : " << input_directory << std::endl;
   std::string input_main_file = "setting.json";
   JSON settingJSON(std::ifstream(input_directory + input_main_file));
   for (const auto &line : settingJSON())
      std::cout << Red << std::setw(30) << line.first << " : " << line.second << colorReset << std::endl;

   // Check for the existence of required keys in the JSON file
   const std::vector<std::string> required_keys = {"output_directory", "CSML", "end_time", "end_time_step", "max_dt", "RK_order", "initial_surface_z_position"};
   for (const auto &key : required_keys)
      if (!settingJSON.find(key))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "setting.json does not have " + key);

   const auto output_directory = settingJSON["output_directory"][0] + id + "/";
   const double CSML = stod(settingJSON["CSML"])[0];
   const double end_time_step = stod(settingJSON["end_time_step"])[0];
   const double end_time = stod(settingJSON["end_time"])[0];
   const double max_dt = stod(settingJSON["max_dt"])[0];
   const int RK_order = stoi(settingJSON["RK_order"])[0];
   const double initial_surface_z_position = stod(settingJSON["initial_surface_z_position"])[0];
   //
   std::filesystem::path output_path = output_directory;
   if (!std::filesystem::exists(output_path)) {
      try {
         std::filesystem::create_directory(output_path);
      } catch (const std::filesystem::filesystem_error &e) {
         std::cerr << e.what() << '\n';
      }
   }
   /* --------------------- copy input directory and files --------------------- */
   {
      std::filesystem::path output_directory_path(output_directory);
      std::cout << "copy input directory and files" << std::endl;
      for (const auto &file : Append(settingJSON["input_files"], input_main_file)) {
         std::filesystem::path source(input_directory + file);
         auto new_path = output_directory_path / std::filesystem::relative(source, ".");
         if (source != new_path) {
            std::filesystem::create_directories(new_path.parent_path());
            if (!std::filesystem::exists(new_path) || std::filesystem::is_regular_file(new_path)) {
               std::cout << Green << "Copying from " << source << Magenta << " to " << new_path << colorReset << std::endl;
               std::filesystem::copy(source, new_path, std::filesystem::copy_options::overwrite_existing);
            }
         }
      }
      //
      std::filesystem::path source_directory = std::filesystem::current_path();
      for (const auto &entry : std::filesystem::directory_iterator(source_directory)) {
         auto filename = entry.path().filename().string();
         auto new_path = output_directory_path / filename;
         if (entry.path() != new_path && (entry.path().extension() == ".hpp" || filename == "main.cpp")) {
            if (!std::filesystem::exists(new_path) || std::filesystem::is_regular_file(new_path)) {
               std::cout << Blue << "Copying from " << entry.path() << Magenta << " to " << new_path << colorReset << std::endl;
               std::filesystem::copy(entry.path(), new_path, std::filesystem::copy_options::overwrite_existing);
            }
         }
      }
   }

   // b! -------------------------------------------------------------------------- */
   TimeWatch watch;

   int minDepth = 1, maxDepth = 5;
   std::cout << "minDepth = " << minDepth << std::endl;
   std::cout << "maxDepth = " << maxDepth << std::endl;
   /* ------------------------------------------------------------------------ */
   std::vector<std::tuple<Network *, Network *, JSON>> all_objects;
   std::vector<std::tuple<Network *, Network *, JSON>> RigidBodies;
   //
   std::vector<std::tuple<Network *, JSON>> probes;
   //
   std::unordered_map<Network *, PVDWriter *> net2PVD;
   auto PVD_SPP = new PVDWriter(output_directory + "SPP.pvd");
   //
   for (const auto &file : settingJSON["input_files"]) {
      JSON J(std::ifstream(input_directory + file));
      if (J["type"][0] == "probe") {
         auto tmp = new Network();
         tmp->setName(J["name"][0]);
         new networkPoint(tmp, {stod(J["location"][0]), stod(J["location"][1]), stod(J["location"][2])});
         probes.push_back({tmp, J});
         net2PVD[tmp] = new PVDWriter(output_directory + J["name"][0] + ".pvd");
      } else {
         auto polyNet = new Network(J["objfile"][0], J["name"][0]);
         auto particlesNet = new Network;
         //
         auto ofs = std::ofstream(output_directory + J["name"][0] + "_init.vtp");
         vtkPolygonWrite(ofs, ToX(polyNet->getFaces()));
         ofs.close();
         //
         polyNet->genOctreeOfFaces({minDepth, maxDepth}, 1);
         all_objects.push_back({particlesNet, polyNet, J});
         if (J["type"][0] == "RigidBody") {
            RigidBodies.push_back({particlesNet, polyNet, J});
            particlesNet->isFluid = polyNet->isFluid = false;
            particlesNet->isRigidBody = polyNet->isRigidBody = true;
         } else if (J["type"][0] == "Fluid") {
            particlesNet->isFluid = polyNet->isFluid = true;
            particlesNet->isRigidBody = polyNet->isRigidBody = false;
         }
         net2PVD[particlesNet] = new PVDWriter(output_directory + J["name"][0] + "_particle.pvd");
         net2PVD[polyNet] = new PVDWriter(output_directory + J["name"][0] + "_polygon.pvd");
      }
   }

   double particle_spacing = stod(settingJSON["particle_spacing"])[0];
   double &ps = particle_spacing;
   double volume = std::pow(particle_spacing, 3);
   //
   CoordinateBounds range;
   for (const auto &[object, polygon, _] : all_objects) {
      polygon->setGeometricProperties();
      range += CoordinateBounds(polygon->getBounds());
   }
   auto [rangeX, rangeY, rangeZ] = range.bounds;
   //
   auto [rX0, rX1] = rangeX;
   auto [rY0, rY1] = rangeY;
   auto [rZ0, rZ1] = rangeZ;
   auto vecX = Subdivide(ps * (int)(rX0 / ps), ps * (int)(rX1 / ps), (int)(rX1 / ps) - (int)(rX0 / ps));
   auto vecY = Subdivide(ps * (int)(rY0 / ps), ps * (int)(rY1 / ps), (int)(rY1 / ps) - (int)(rY0 / ps));
   auto vecZ = Subdivide(ps * (int)(rZ0 / ps), ps * (int)(rZ1 / ps), (int)(rZ1 / ps) - (int)(rZ0 / ps));

   Tddd closest = {0., 0., 0.};
   double min = 1E+20;
   for (const auto &x : vecX)
      for (const auto &y : vecY)
         for (const auto &z : vecZ)
            if (min >= Norm(Tddd{x, y, z})) {
               min = Norm(Tddd{x, y, z});
               closest = {x, y, z};
            }

   // double rand = RandomReal({1E-10, 1E-10});
   for (const auto &x : vecX)
      for (const auto &y : vecY)
         for (const auto &z : vecZ) {
            Tddd xyz = {x, y, z};
            xyz -= closest;
            xyz += particle_spacing / 2.;
            xyz += particle_spacing * 1E-12;  // std::numbers::pi / 1000.;
            for (const auto &[object, polygon, _] : all_objects) {
               // 優先順位で粒子を配置
               auto [isInside, cell, f] = polygon->isInside_MethodOctree(xyz, particle_spacing * 1E-10);
               if (isInside) {
                  auto p = new networkPoint(object, xyz);
                  p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
                  p->particle_spacing = particle_spacing;
                  p->C_SML_next = p->C_SML = 3.5;  // CSML;
                  p->lap_U.fill(0.);
                  p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
                  p->isFluid = object->isFluid;
                  p->isNeumannSurface = false;
                  break;
               }
            }
         }

   auto WaveTank = std::get<0>(all_objects[0]);
   // auto Wall = std::get<0>(all_objects[1]);
   auto Fluid = std::get<0>(all_objects[1]);
   for (const auto &p : Fluid->getPoints()) {
      // auto rand = RandomRealTddd(-particle_spacing * 1E-14, particle_spacing * 1E-14);
      // p->setX(p->X + rand);
      p->setDensityVolume(_WATER_DENSITY_, volume);
      p->setX(p->X);
   }

   /*DOC_EXTRACT 0_1_0_preparation_wall_and_freesurface

   ## N.S.方程式を解く前の準備

   壁粒子の法線ベクトル`p->v_to_surface_SPH`を計算する．

   */

   // b# -------------------------------------------------------------------------- */
   // b#                             外向きベクトルの設定                               */
   // b# -------------------------------------------------------------------------- */
   for (const auto &[particlesNet, poly, J] : all_objects) {
      particlesNet->makeBucketPoints(3. * particle_spacing);
      poly->makeBucketFaces(3. * particle_spacing);
      //
      std::cout << Yellow << poly->getName() << " makeBucketFaces" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

      double ps = particle_spacing;

      if (particlesNet->isRigidBody)
         for (const auto &p : particlesNet->getPoints()) {
            // p->v_to_surface_SPH = poly->interpolateVector(p->X);
            auto [X, f] = Nearest_(p->X, poly->getFaces());
            // p->v_to_surface_SPH = ((int)(Distance(X, p) / ps) + 1 / 2.) * ps * Normalize(X - p->X);
            // p->v_to_surface_SPH = p->v_to_surface_SPH = ((int)(Distance(X, p) / ps)) * ps * Normalize(X - p->X);

            p->v_to_surface_SPH = p->normal_SPH = ((int)(Distance(X, p) / ps) + .5) * ps * Normalize(X - p->X);

            // p->v_to_surface_SPH = p->normal_SPH = X - p->X;

            // p->v_to_surface_SPH = p->normal_SPH = ((int)(Distance(X, p) / ps) + 1E-12) * ps * Normalize(X - p->X);
            p->mirroring_face = f;
         }
      // vtkPolygonWriter<networkPoint *> vtp;
      // vtp.add(particlesNet->getPoints());
      // std::unordered_map<networkPoint *, Tddd> normal_SPH;
      // for (const auto &p : particlesNet->getPoints())
      //    v_to_surface_SPH[p] = p->v_to_surface_SPH;
      // vtp.addPointData("v_to_surface_SPH", normal_SPH);
      // std::ofstream ofs(output_directory + J.at("name")[0] + "_check_normal.vtp");
      // vtp.write(ofs);
      // ofs.close();
   }

   // b# -------------------------------------------------------------------------- */
   // b# -------------------------------------------------------------------------- */
   // b# -------------------------------------------------------------------------- */
   for (const auto &[WaveTank, _, J] : all_objects) {
      std::cout << "WaveTank->getPoints().size()=" << WaveTank->getPoints().size() << std::endl;
      {
         auto ofs = std::ofstream(output_directory + J.at("name")[0] + ".vtp");
         vtkPolygonWrite(ofs, ToX(WaveTank->getFaces()));
         ofs.close();
      }
      WaveTank->setGeometricProperties();
   }
   /* -------------------------------------------------------------------------- */
   double simulation_time = 0.;
   // PVDWriter pvdWallSPH(output_directory + "Wall.pvd");q
   int time_step = 0, k = 0, i = 0, j = 0, l = 0;

   auto active_RigidBodies = RigidBodies;
   auto active_all_objects = all_objects;

   auto pack_actives = [&](auto &ACTIVES, const auto &ALL) {
      ACTIVES.clear();
      for (const auto &PART_POLY_JSON : ALL) {
         auto [particlesNet, poly, J] = PART_POLY_JSON;
         if (J.find("inactivate") || J.find("inactive")) {
            if (J.find("inactivate")) {
               if (J.at("inactivate").size() != 2)
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "inactivate size is not 2");
               double start = stod(J.at("inactivate")[0]);
               double end = stod(J.at("inactivate")[1]);
               if (!Between(simulation_time, {start, end}))
                  ACTIVES.emplace_back(PART_POLY_JSON);
            }
            if (J.find("inactive")) {
               if (J.find("inactive") && J.at("inactive").size() != 2)
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "inactive size is not 2");
               double start = stod(J.at("inactive")[0]);
               double end = stod(J.at("inactive")[1]);
               if (!Between(simulation_time, {start, end}))
                  ACTIVES.emplace_back(PART_POLY_JSON);
            }
         } else
            ACTIVES.emplace_back(PART_POLY_JSON);
      }
   };

   for (auto time_step = 0; time_step < end_time_step; ++time_step) {

      Print(Magenta, "time_step = ", time_step, ", simulation_time = ", simulation_time, ", end_time = ", end_time, ", max_dt = ", max_dt, ", RK_order = ", RK_order);

      {
         auto points = Fluid->getPoints();
         for (auto p : points)
            if (p->isAuxiliary) {
               p->surfacePoint->auxPoint = nullptr;
               delete p;
            }
      }

      {
         auto points = Fluid->getPoints();
         for (auto p : points) {
            p->isAuxiliary = false;
         }
      }

      pack_actives(active_all_objects, all_objects);
      pack_actives(active_RigidBodies, RigidBodies);

      _ALL_NET_.clear();
      for (const auto &[particlesNet, _, __] : all_objects)
         _ALL_NET_.emplace(particlesNet);

      if (end_time < simulation_time)
         break;

      // int N = 100;
      // if (time_step == N) {
      //    for (const auto &[object, _, __] : all_objects)
      //       for (const auto &p : object->getPoints())
      //          p->mu_SPH = _WATER_MU_10deg_;
      // } else if (time_step < N) {
      //    for (const auto &[object, _, __] : all_objects)
      //       for (const auto &p : object->getPoints())
      //          p->mu_SPH = _WATER_MU_10deg_ * 10.;
      // }

      // developByEISPH(Fluid, RigidBodies, simulation_time, CSML, particle_spacing, time_step < 50 ? 1E-12 : max_dt);

      developByEISPH(Fluid,
                     active_RigidBodies,
                     simulation_time,
                     CSML,
                     particle_spacing,
                     time_step < 50 ? max_dt / 100 : max_dt,
                     // max_dt,
                     RK_order);

      // freeze particle a while
      for (const auto &p : Fluid->getPoints()) {
         p->vec_time_SPH.emplace_back(simulation_time);
         p->vec_U_SPH.emplace_back(p->U_SPH);
      }

      /*DOC_EXTRACT 0_3_0_SPH

      ## 出力

      */

      if (time_step % 5 == 0) {

         int num_aux = 0, num_fluid = 0, num_neuman_surface = 0, num_surface = 0;
         for (const auto &p : Fluid->getPoints()) {
            if (p->isAuxiliary)
               num_aux++;
            else if (p->isFluid)
               num_fluid++;
            else if (p->isNeumannSurface)
               num_neuman_surface++;
            else if (p->isSurface)
               num_surface++;
         }

         std::cout << std::setw(10) << "auxiliary : " << num_aux << std::endl;
         std::cout << std::setw(10) << "fluid : " << num_fluid << std::endl;
         std::cout << std::setw(10) << "wall : " << num_neuman_surface << std::endl;
         std::cout << std::setw(10) << "surface : " << num_surface << std::endl;

         auto Output = [&](Network *Fluid, const std::string &name, const int i) {
            if (Fluid != nullptr) {
               std::vector<networkPoint *> points;
               for (const auto &p : Fluid->getPoints())
                  if (!p->isAuxiliary && p->isCaptured)
                     points.emplace_back(p);

               vtkPolygonWriter<networkPoint *> vtp;
               vtp.add(points);
               setDataOmitted(vtp, Fluid);
               std::ofstream ofs(output_directory + name + std::to_string(i) + ".vtp");
               vtp.write(ofs);
               ofs.close();

               if (net2PVD.find(Fluid) != net2PVD.end()) {
                  auto pvd = net2PVD[Fluid];
                  pvd->push("./" + name + std::to_string(i) + ".vtp", simulation_time);
                  pvd->output();
               } else {
                  if (net2PVD.find(nullptr) == net2PVD.end())
                     net2PVD.insert({nullptr, new PVDWriter(output_directory + "other.pvd")});

                  net2PVD[nullptr]->push("./" + name + std::to_string(i) + ".vtp", simulation_time);
                  net2PVD[nullptr]->output();
               }
            }
         };

         // pallarelize Output using OpenMP

#pragma omp parallel sections
         {
#pragma omp section
            {
               Output(Fluid, "FluidSPH", i);
            }
#pragma omp section
            {
               Output(WaveTank, "WaveTankSPH", i);
            }
         }
         i++;

         // int count = 0;
         // std::vector<networkPoint *> a_wall_p = {};
         // if (time_step % 5 == 0) {
         //    vtkPolygonWriter<networkPoint *> vtp;
         //    for (const auto &p : WaveTank->getPoints())
         //       if (p->isCaptured)
         //       // if (p->isChecked)
         //       {
         //          a_wall_p.emplace_back(p);
         //          vtp.add(p);
         //       }
         //    setDataOmitted(vtp, WaveTank);
         //    std::ofstream ofs(output_directory + "WaveTankSPH" + std::to_string(j) + ".vtp");
         //    vtp.write(ofs);
         //    ofs.close();
         //    //
         //    net2PVD[WaveTank]->push("./WaveTankSPH" + std::to_string(j) + ".vtp", simulation_time);
         //    net2PVD[WaveTank]->output();
         //    j++;
         // }

         /* -------------------------------------------------------------------------- */

         // count = 0;
         // for (const auto &p : a_wall_p) {
         //    vtkPolygonWriter<networkPoint *> vtp;
         //    auto X = p->X + 2 * p->v_to_surface_SPH;
         //    for (const auto &[particlesNet, _, __] : all_objects) {
         //       particlesNet->BucketPoints.apply(X, p->SML(), [&](const auto &q) {
         //          if (Distance(q, X) < p->SML())
         //             vtp.add(q);
         //       });
         //    }
         //    setDataOmitted(vtp, WaveTank);
         //    setDataOmitted(vtp, Fluid);
         //    std::ofstream ofs(output_directory + "referenced_points_by_a_wall_p_" + std::to_string(count++) + "_" + std::to_string(i) + ".vtp");
         //    vtp.write(ofs);
         //    ofs.close();
         //    break;
         // }

         //
         // pvdWaveTankSPH.push("./WaveTankSPH" + std::to_string(i) + ".vtp", simulation_time);
         // pvdWaveTankSPH.output();
         // pvdWallSPH.push("./WallSPH" + std::to_string(i) + ".vtp", simulation_time);
         // pvdWallSPH.output();

         // b# -------------------------------------------------------------------------- */
         // b#                              output JSON files                             */
         // b# -------------------------------------------------------------------------- */

         // ポリゴンの節点上で圧力を計算する．

         //! -------------------------------------------------------------------------- */
         //!                                    プローブ                                  */
         //! -------------------------------------------------------------------------- */
         Print("プローブ出力");
         jsonout.push("time", simulation_time);
         for (const auto &[probe, J] : probes) {
            vtkPolygonWriter<networkPoint *> vtp;
            for (const auto &p : probe->getPoints())
               vtp.add(p);

            std::unordered_map<networkPoint *, double> PRESSURE, RHO, PRESSURE_N, RHO_N;
            for (const auto &p : probe->getPoints()) {
               double pressure = 0, rho = 0, w, total = 0.;
               double pressure_normalized = 0, rho_normalized = 0;
               int c = 0;
               for (const auto &[particlesNet, _, __] : all_objects) {
                  particlesNet->BucketPoints.apply(p->X, particle_spacing * CSML, [&](const auto &q) {
                     if (isFinite(q->p_SPH) && isFinite(q->volume)) {
                        w = q->volume * w_Bspline(Norm(q->X - p->X), particle_spacing * CSML);
                        pressure += q->p_SPH * w;
                        rho += q->rho * w;
                        if (particlesNet->isFluid) {
                           total += w;
                           pressure_normalized += q->p_SPH * w;
                           rho_normalized += q->rho * w;
                        }
                        c++;
                     }
                  });
               }
               if (total > 1E-15) {
                  PRESSURE[p] = pressure;
                  RHO[p] = rho;
                  PRESSURE_N[p] = pressure_normalized / total;
                  RHO_N[p] = rho_normalized / total;
               } else {
                  PRESSURE[p] = 0.;
                  RHO[p] = 0.;
                  PRESSURE_N[p] = 0.;
                  RHO_N[p] = 0.;
               }
            }
            vtp.addPointData("pressure", PRESSURE);
            vtp.addPointData("density", RHO);
            auto name = output_directory + J.at("name")[0] + "_" + std::to_string(k) + ".vtp";
            std::ofstream ofs(name);
            vtp.write(ofs);
            ofs.close();
            //
            net2PVD[probe]->push(name, simulation_time);
            net2PVD[probe]->output();
            // b# -------------------------------------------------------------------------- */
            // b#                              output JSON files                             */
            // b# -------------------------------------------------------------------------- */
            jsonout.push(probe->getName() + "_pressure", std::get<1>(*PRESSURE.begin()));
            jsonout.push(probe->getName() + "_rho", std::get<1>(*RHO.begin()));
            jsonout.push(probe->getName() + "_pressure_normalized", std::get<1>(*PRESSURE_N.begin()));
            jsonout.push(probe->getName() + "_rho_normalized", std::get<1>(*RHO_N.begin()));
            std::ofstream os(output_directory + "/result.json");
            jsonout.output(os);
            os.close();
         }

         //* -------------------------------------------------------------------------- */
         //*                                   オブジェクト                               */
         //* -------------------------------------------------------------------------- */

         /*DOC_EXTRACT 0_3_0_SPH

         ## 出力（ポリゴン）

         */

         Print("オブジェクト出力");
         for (const auto &[_, poly, J] : all_objects) {
            vtkPolygonWriter<networkPoint *> vtp;
            for (const auto &f : poly->getFaces()) {
               auto [p0, p1, p2] = f->getPoints();
               vtp.add(p0);
               vtp.add(p1);
               vtp.add(p2);
               vtp.addPolygon(std::array<networkPoint *, 4>{p0, p1, p2});
            }
            std::unordered_map<networkPoint *, double> PRESSURE, RHO;
            for (const auto &p : poly->getPoints()) {
               double pressure = 0, rho = 0;
               int c = 0;
               for (const auto &[particlesNet, _, __] : all_objects) {
                  particlesNet->BucketPoints.apply(p->X, particle_spacing * CSML, [&](const auto &q) {
                     pressure += q->p_SPH * q->volume * w_Bspline(Norm(q->X - p->X), particle_spacing * CSML);
                     rho += q->mass * w_Bspline(Norm(q->X - p->X), particle_spacing * CSML);
                     c++;
                  });
               }
               if (c > 10) {
                  PRESSURE[p] = pressure;
                  RHO[p] = rho;
               } else {
                  PRESSURE[p] = 0.;
                  RHO[p] = rho;
               }
            }
            vtp.addPointData("pressure", PRESSURE);
            vtp.addPointData("density", RHO);
            auto name = output_directory + J.at("name")[0] + "_polygon" + std::to_string(k) + ".vtp";
            std::ofstream ofs(name);
            vtp.write(ofs);
            ofs.close();
            //
            net2PVD[poly]->push(name, simulation_time);
            net2PVD[poly]->output();
         }
         /* -------------------------------------------------------------------------- */
         auto bound = CoordinateBounds(WaveTank->bounds) + CoordinateBounds(Fluid->bounds);
         auto points = Fluid->getPoints();
         for (const auto &p : points) {
            if (!bound.isInside(p->X)) {
               std::cout << Blue << "delete" << std::endl;
               delete p;
            }
         }

         k++;
      }
   }
};