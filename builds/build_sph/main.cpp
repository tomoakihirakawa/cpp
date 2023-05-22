#define _debugging_
#include <filesystem>
#include <utility>
#define DEM
#include "Network.hpp"
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
/**
 * ./main bunny.obj 7 1 3
 */
/* -------------------------------------------------------------------------- */
JSONoutput jsonout;
int main(int arg, char **argv) {
   if (arg <= 1)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "argv <= 1. write input json file directory!\\ex.\\$ ./main ./input");
   std::string input_directory{argv[1]};  // input directory
   input_directory += "/";
   //
   std::string id = "";
   if (arg >= 3)
      id = argv[2];  // input directory
   // b! -------------------------------------------------------------------------- */
   std::cout << "input_directory : " << input_directory << std::endl;
   std::string input_main_file = "setting.json";
   JSON settingJSON(std::ifstream(input_directory + input_main_file));
   for (const auto &line : settingJSON())
      std::cout << Red << std::setw(30) << line.first << " : " << line.second << colorOff << std::endl;

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
         std::filesystem::create_directories(new_path.parent_path());
         std::filesystem::copy(source, new_path, std::filesystem::copy_options::overwrite_existing);
      }
   }
   // b! -------------------------------------------------------------------------- */
   Timer timer;
   int minDepth = 1, maxDepth = 5;
   std::cout << "minDepth = " << minDepth << std::endl;
   std::cout << "maxDepth = " << maxDepth << std::endl;
   std::cout << "timer : " << timer() << std::endl;
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
         new networkPoint(tmp, {stod(J["location"][0]),
                                stod(J["location"][1]),
                                stod(J["location"][2])});
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
   //
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
   auto vecX = Subdivide({ps * (int)(rX0 / ps), ps * (int)(rX1 / ps)}, (int)(rX1 / ps) - (int)(rX0 / ps));
   auto vecY = Subdivide({ps * (int)(rY0 / ps), ps * (int)(rY1 / ps)}, (int)(rY1 / ps) - (int)(rY0 / ps));
   auto vecZ = Subdivide({ps * (int)(rZ0 / ps), ps * (int)(rZ1 / ps)}, (int)(rZ1 / ps) - (int)(rZ0 / ps));

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
            // xyz += 1E-4;  // std::numbers::pi / 1000.;
            for (const auto &[object, polygon, _] : all_objects) {
               // 優先順位で粒子を配置
               auto [isInside, cell, f] = polygon->isInside_MethodOctree(xyz, particle_spacing * 1E-10);
               if (isInside) {
                  auto p = new networkPoint(object, xyz);
                  p->setDensityVolume(_WATER_DENSITY_, volume);  // 質量(mass)は関数内部で自動で決められる
                  p->radius_SPH = CSML * ps;
                  p->C_SML = CSML;
                  p->lap_U.fill(0.);
                  p->pressure_SPH = _WATER_DENSITY_ * _GRAVITY_ * (initial_surface_z_position - std::get<2>(p->X));
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
      p->setDensityVolume(_WATER_DENSITY_, std::pow(particle_spacing, 3));
      p->setX(p->X);
   }
   // b# -------------------------------------------------------------------------- */
   // b#                             外向きベクトルの設定                               */
   // b# -------------------------------------------------------------------------- */
   for (const auto &[particlesNet, poly, J] : all_objects) {
      particlesNet->makeBucketPoints(particle_spacing);
      //
      std::cout << "timer : " << timer() << std::endl;

      if (particlesNet->isRigidBody)
         for (const auto &p : particlesNet->getPoints()) {
            // p->normal_SPH = poly->interpolateVector(p->X);
            auto [X, f] = Nearest_(p->X, poly->getFaces());
            // p->normal_SPH = ((int)(Distance(X, p) / particle_spacing) + 1 / 2.) * particle_spacing * Normalize(X - p->X);
            p->vector_to_wall_SPH = p->normal_SPH = ((int)(Distance(X, p) / particle_spacing) + 1E-20) * particle_spacing * Normalize(X - p->X);
            // p->normal_SPH = ((int)(Distance(X, p) / particle_spacing) + 1 / 4.) * particle_spacing * Normalize(X - p->X);
            // p->normal_SPH = X - p->X;
            p->mirroring_face = f;
         }
      vtkPolygonWriter<networkPoint *> vtp;
      vtp.add(particlesNet->getPoints());
      std::unordered_map<networkPoint *, Tddd> normal_SPH;
      for (const auto &p : particlesNet->getPoints())
         normal_SPH[p] = p->normal_SPH;
      vtp.addPointData("normal_SPH", normal_SPH);
      std::ofstream ofs(output_directory + J.at("name")[0] + "_check_normal.vtp");
      vtp.write(ofs);
      ofs.close();
   }
   // b# -------------------------------------------------------------------------- */
   // b# -------------------------------------------------------------------------- */
   // b# -------------------------------------------------------------------------- */
   //
   std::cout << "timer : " << timer() << ", networkPointの生成" << std::endl;
   for (const auto &[WaveTank, _, J] : all_objects) {
      std::cout << "WaveTank->getPoints().size()=" << WaveTank->getPoints().size() << std::endl;
      {
         auto ofs = std::ofstream(output_directory + J.at("name")[0] + ".vtp");
         vtkPolygonWrite(ofs, ToX(WaveTank->getPoints()));
         ofs.close();
      }
      WaveTank->setGeometricProperties();
   }
   /* -------------------------------------------------------------------------- */
   double real_time = 0.;
   // PVDWriter pvdWallSPH(output_directory + "Wall.pvd");q
   int time_step = 0, k = 0, i = 0, j = 0, l = 0;
   for (auto time_step = 0; time_step < end_time_step; ++time_step) {
      if (end_time < real_time)
         break;

      // int N = 1000;

      // if (time_step == N) {
      //    for (const auto &[object, _, __] : all_objects)
      //       for (const auto &p : object->getPoints())
      //          p->mu_SPH = _WATER_MU_10deg_;
      // } else if (time_step < N) {
      //    for (const auto &[object, _, __] : all_objects)
      //       for (const auto &p : object->getPoints())
      //          p->mu_SPH = _WATER_MU_10deg_ * 10;
      // }

      // developByEISPH(Fluid, RigidBodies, real_time, CSML, particle_spacing, time_step < 50 ? 1E-12 : max_dt);
      developByEISPH(Fluid,
                     RigidBodies,
                     real_time,
                     CSML,
                     particle_spacing,
                     time_step < 5 ? max_dt / 100 : max_dt,
                     RK_order);

      std::cout << "real_time = " << real_time << std::endl;

      // freeze particle a while
      // for (const auto &p : Fluid->getPoints()) {
      //    if (time_step < 10)
      //       p->U_SPH.fill(0.);
      // }

      // 出力
      if (time_step % 5 == 0) {
         auto Output = [&](Network *Fluid, const std::string &name, const int i) {
            if (Fluid != nullptr) {
               vtkPolygonWriter<networkPoint *> vtp;
               vtp.add(ToVector(Fluid->getPoints()));
               setDataOmitted(vtp, Fluid);
               std::ofstream ofs(output_directory + name + std::to_string(i) + ".vtp");
               vtp.write(ofs);
               ofs.close();

               if (net2PVD.find(Fluid) != net2PVD.end()) {
                  auto pvd = net2PVD[Fluid];
                  pvd->push("./" + name + std::to_string(i) + ".vtp", real_time);
                  pvd->output();
               } else {
                  if (net2PVD.find(nullptr) == net2PVD.end())
                     net2PVD.insert({nullptr, new PVDWriter(output_directory + "other.pvd")});

                  net2PVD[nullptr]->push("./" + name + std::to_string(i) + ".vtp", real_time);
                  net2PVD[nullptr]->output();
               }
            }
         };

         Output(Fluid, "FluidSPH", i);
         Output(Fluid->surfaceNet, "surdaceNet", i);

         i++;

         int count = 0;
         std::vector<networkPoint *> a_wall_p = {};
         if (time_step % 5 == 0) {
            vtkPolygonWriter<networkPoint *> vtp;
            for (const auto &p : WaveTank->getPoints())
               if (p->isCaptured) {
                  if ((count++) % 500 == 0)
                     a_wall_p.emplace_back(p);
                  vtp.add(p);
               }
            setDataOmitted(vtp, WaveTank);
            std::ofstream ofs(output_directory + "WaveTankSPH" + std::to_string(j) + ".vtp");
            vtp.write(ofs);
            ofs.close();
            //
            net2PVD[WaveTank]->push("./WaveTankSPH" + std::to_string(j) + ".vtp", real_time);
            net2PVD[WaveTank]->output();
            j++;
         }

         count = 0;
         for (const auto &p : a_wall_p) {
            vtkPolygonWriter<networkPoint *> vtp;
            auto X = p->X + 2 * p->normal_SPH;
            for (const auto &[particlesNet, _, __] : all_objects) {
               particlesNet->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
                  if (Distance(q, X) < p->radius_SPH)
                     vtp.add(q);
               });
            }
            setDataOmitted(vtp, WaveTank);
            setData(vtp, Fluid, X);
            std::ofstream ofs(output_directory + "referenced_points_by_a_wall_p_" + std::to_string(count++) + "_" + std::to_string(i) + ".vtp");
            vtp.write(ofs);
            ofs.close();
            break;
         }

         //
         // pvdWaveTankSPH.push("./WaveTankSPH" + std::to_string(i) + ".vtp", real_time);
         // pvdWaveTankSPH.output();
         // pvdWallSPH.push("./WallSPH" + std::to_string(i) + ".vtp", real_time);
         // pvdWallSPH.output();
         // b# -------------------------------------------------------------------------- */
         // b#                              output JSON files                             */
         // b# -------------------------------------------------------------------------- */
         // ポリゴンの節点上で圧力を計算する．
         //! -------------------------------------------------------------------------- */
         //!                                    プローブ                                  */
         //! -------------------------------------------------------------------------- */
         Print("プローブ出力");
         jsonout.push("time", real_time);
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
            net2PVD[probe]->push(name, real_time);
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
            net2PVD[poly]->push(name, real_time);
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