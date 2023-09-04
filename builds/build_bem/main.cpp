/*DOC_EXTRACT 0_0_BEM
\insert{REVIEW.md}
*/

// #define _debugging_

#define BEM
#define use_lapack
int time_step;
double real_time = 0;

#define simulation
#include <sys/utsname.h>
#include <unistd.h>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <regex>
#include "Network.hpp"
#include "integrationOfODE.hpp"
#include "kernelFunctions.hpp"
#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"

// pvd cpg_pvd("./vtu/bem.pvd");

#include "BEM.hpp"
#include "svd.hpp"

void logMachineInformation(std::ofstream &ofs) {
   char hostname[256];
   if (gethostname(hostname, sizeof(hostname)) != 0) {
      ofs << "Failed to get hostname\n";
      return;
   }

   struct utsname unameData;
   if (uname(&unameData) != 0) {
      ofs << "Failed to get system information\n";
      return;
   }

   ofs << "Machine Information:" << std::endl;
   ofs << "  Hostname: " << hostname << std::endl;
   ofs << "  Operating System: " << unameData.sysname << " " << unameData.release << std::endl;
   ofs << std::endl;  // Add a blank line for separation
}

int main(int argc, char **argv) {

   /* -------------------------------------------------------------------------- */
   /*                           Set up logging to file                           */
   /* -------------------------------------------------------------------------- */

   // Read the contents of the existing log file
   std::ifstream ifs("log.txt");
   std::string existingContent;
   if (ifs.is_open()) {
      existingContent.assign(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
   }

   // Open log file in trunc mode to clear its contents
   std::ofstream ofs("log.txt");
   if (!ofs.is_open()) {
      std::cerr << "Failed to open log file\n";
      return 1;
   }
   ofs << std::setw(10) << std::setfill('-') << "" << std::endl;

   // Get current time
   auto now = std::chrono::system_clock::now();
   std::time_t time = std::chrono::system_clock::to_time_t(now);
   ofs << "Execution time: " << std::ctime(&time) << std::endl;

   // Log machine information
   logMachineInformation(ofs);

   // Check command line arguments
   if (argc <= 1) {
      ofs << "No arguments provided. Write input JSON file directory.\n";
      return 1;
   }

   // Write command line arguments to log file
   ofs << "Command line arguments:" << std::endl;
   for (int i = 1; i < argc; ++i)
      ofs << "  arg" << i << ": " << argv[i] << std::endl;

   // Write the existing log content back to the log file after the new log entry
   ofs << existingContent;
   ofs << std::setw(10) << std::setfill('-') << "" << std::endl;

   /* -------------------------------------------------------------------------- */

   if (argc <= 1)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "argv <= 1. write input json file directory!\nex.\n$ ./main ./input");

   std::string input_directory(argv[1]);  // input directory
   input_directory += "/";
   std::cout << input_directory << std::endl;

   /* -------------------------------------------------------------------------- */
   /*                           read setting.json                                */
   /* -------------------------------------------------------------------------- */

   // Read settings from a JSON file and store values
   JSON settingJSON(input_directory + "/setting.json");
   for (auto &[key, value] : settingJSON()) std::cout << key << ": " << value << std::endl;

   // Check for the existence of required keys in the JSON file
   if (!std::ranges::all_of(settingJSON.find({"end_time_step", "end_time", "output_directory", "max_dt", "input_files"}), [](auto v) { return v; }))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "there is missing key in setting.json");

   // Store values from the JSON file into variables
   const std::string output_directory = settingJSON.at("output_directory")[0];
   const double max_dt = stod(settingJSON.at("max_dt"))[0];
   const int end_time_step = stoi(settingJSON.at("end_time_step"))[0];
   const double end_time = stod(settingJSON.at("end_time"))[0];

   std::cout << "" << std::endl;

   if (settingJSON.find("WATER_DENSITY", [](auto STR_VEC) { _WATER_DENSITY_ = stod(STR_VEC[0]); }))
      std::cout << "WATER_DENSITY: " << _WATER_DENSITY_ << std::endl;

   if (settingJSON.find("GRAVITY", [](auto STR_VEC) { _GRAVITY_ = stod(STR_VEC[0]); }))
      std::cout << "GRAVITY: " << _GRAVITY_ << std::endl;

   // Set default values for optional keys and update if found in the JSON file
   double stop_remesh_time = 1E+10;
   double force_remesh_time = 0;
   int grid_refinement = 0;

   settingJSON.find("stop_remesh_time", [&](auto STR_VEC) { stop_remesh_time = stod(STR_VEC[0]); });
   settingJSON.find("force_remesh_time", [&](auto STR_VEC) { force_remesh_time = stod(STR_VEC[0]); });
   settingJSON.find("grid_refinement", [&](auto STR_VEC) { grid_refinement = stoi(STR_VEC[0]); });
   /* -------------------------------------------------------------------------- */
   /*                read JSON files, input_files in setting.json                */
   /* -------------------------------------------------------------------------- */

   // Define containers and helper functions for network objects and output information
   std::map<Network *, outputInfo> NetOutputInfo;
   auto setup_network_output_info = [&NetOutputInfo](Network *net, const JSON &injson, const std::string &output_directory) {
      std::cout << "setup_network_output_info" << std::endl;
      auto output_name = injson.at("name")[0];
      NetOutputInfo[net].pvd_file_name = output_name;
      NetOutputInfo[net].vtu_file_name = output_name + "_";
      NetOutputInfo[net].PVD = new PVDWriter(output_directory + "/" + output_name + ".pvd");
   };

   std::vector<Network *> FluidObject, RigidBodyObject, SoftBodyObject;
   auto initialize_network_objects = [&FluidObject, &RigidBodyObject, &SoftBodyObject](Network *net, const JSON &injson) {
      std::cout << "initialize_network_objects" << std::endl;
      auto type = injson.at("type")[0];
      if (type == "RigidBody") {
         RigidBodyObject.emplace_back(net);
         net->isRigidBody = true;
         net->isSoftBody = false;
      } else if (type == "SoftBody" || type == "FixedBody") {
         SoftBodyObject.emplace_back(net);
         net->isRigidBody = false;
         net->isSoftBody = true;
      } else if (type == "Fluid") {
         FluidObject.emplace_back(net);
         net->isRigidBody = net->isSoftBody = false;
      }
      net->isFixed = (injson.find("isFixed") && stob(injson.at("isFixed"))[0]);
      for (auto i = 0; i < 10; ++i)
         AreaWeightedSmoothingPreserveShape(net->getPoints(), 0.1);

      net->isFloatingBody = (injson.find("velocity") && injson.at("velocity")[0] == "floating");
      net->isFloatingBody = net->isFloatingBody || (injson.find("acceleration") && injson.at("acceleration")[0] == "floating");

      net->resetInitialX();
      net->setGeometricProperties();
   };

   for (auto input_file_name : settingJSON["input_files"]) {
      std::cout << Green << input_directory + input_file_name << colorOff << std::endl;
      JSON injson(input_directory + input_file_name);

      // check required keys
      if (!std::ranges::all_of(injson.find({"name", "objfile", "type"}), [](auto v) { return v; }))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, input_directory + input_file_name + " has missing key");

      // display contents of the JSON file
      for (auto &[key, value] : injson()) std::cout << Green << key << colorOff << ": " << value << std::endl;

      auto object_name = injson.at("name")[0];
      if (!injson.find("ignore") || !stob(injson["ignore"])[0]) {
         auto net = new Network(injson.at("objfile")[0], object_name);
         net->inputJSON = injson;
         net->applyTransformations(injson);
         setup_network_output_info(net, injson, output_directory);
         initialize_network_objects(net, injson);
         std::filesystem::copy_file(input_directory + input_file_name, output_directory + "/" + input_file_name, std::filesystem::copy_options::overwrite_existing);
         mk_vtu(output_directory + "/" + object_name + "_init.vtu", {net->getFaces()});
      } else {
         Print("skipped");
      }
   }

   /* ----------------- Create output directory and copy files ----------------- */

   std::filesystem::create_directories(output_directory);
   std::filesystem::copy_file(input_directory + "/setting.json", output_directory + "/setting.json", std::filesystem::copy_options::overwrite_existing);
   std::filesystem::copy_file("./main.cpp", output_directory + "/main.cpp", std::filesystem::copy_options::overwrite_existing);
   //
   std::regex pattern("^BEM.*\\.hpp$");
   for (auto &entry : std::filesystem::directory_iterator("."))
      if (std::regex_match(entry.path().filename().string(), pattern))
         std::filesystem::copy_file(entry.path(), output_directory / entry.path().filename(), std::filesystem::copy_options::overwrite_existing);

   /* -------------------------------------------------------------------------- */

   auto water = FluidObject[0];
   PVDWriter cornerPointsPVD(output_directory + "/cornerPointsPVD.pvd");
   PVDWriter cornerPVD(output_directory + "/corner.pvd");
   Print("setting done");

   /* -------------------------------------------------------------------------- */
   try {
      //  b* ------------------------------------------------------ */
      //  b*                         メインループ                      */
      //  b* ------------------------------------------------------ */
      TimeWatch watch;
      for (time_step = 0; time_step < end_time_step; time_step++) {
         if (end_time < real_time) break;
         show_info(*water);
         //! 体積を保存するようにリメッシュする必要があるだろう．
         setBoundaryTypes(*water, Join(RigidBodyObject, SoftBodyObject));
         double rad = M_PI / 180;

         flipIf(*water, {10 * rad, rad}, {10 * rad, rad}, false);

         // b# ------------------------------------------------------ */
         // b#                       刻み時間の決定                     */
         // b# ------------------------------------------------------ */

         const auto Points = water->getPoints();
         const auto Faces = water->getFaces();
         double dt = dt_CFL(*water, max_dt, .5);
         Print("===========================================================================");
         Print("       dt :" + Red + std::to_string(dt) + colorOff);
         Print("time_step :" + Red + std::to_string(time_step) + colorOff);
         Print("real time :" + Red + std::to_string(real_time) + colorOff);
         Print("---------------------------------------------------------------------------");

         double spacing = Mean(extLength(water->getLines())) * 3;
         Buckets<networkFace *> FMM_BucketsFaces(water->bounds, spacing);
         Buckets<networkPoint *> FMM_BucketsPoints(water->bounds, spacing);
         for (const auto &f : water->getFaces())
            FMM_BucketsFaces.add(f->getXtuple(), f);

         for (const auto &p : water->getPoints())
            FMM_BucketsPoints.add(ToX(p), p);

         // b@ ------------------------------------------------------ */
         // b@        初期値問題を解く（時間微分方程式を数値積分する）           */
         // b@ ------------------------------------------------------ */
         int RK_order = 4;
         for (const auto &p : Points) {
            p->RK_phi.initialize(dt, real_time, std::get<0>(p->phiphin), RK_order);
            p->RK_X.initialize(dt, real_time, ToX(p), RK_order);
         }
         for (const auto &net : RigidBodyObject) {
            net->RK_COM.initialize(dt, real_time, net->COM, RK_order);
            net->RK_Q.initialize(dt, real_time, net->Q(), RK_order);
            net->RK_Velocity.initialize(dt, real_time, net->velocity, RK_order);
            //
            if (net->interp_accel.size() > 10)
               net->interp_accel.pop();
            net->interp_accel.push(real_time, net->acceleration);
         }
         for (const auto &net : SoftBodyObject) {
            // !いらないはずのもの
            for (const auto &p : net->getPoints())
               p->RK_X.initialize(dt, real_time, ToX(p), RK_order);
         }
         // b@ ----------------------------------------------------- */
         /*DOC_EXTRACT 0_1_BEM

         ## 計算プログラムの概要

         | 項目 | 詳細|
         |---:|:---|
         | 要素 | 線形三角要素 |
         | 時間発展方法 | 4次のルンゲクッタ |
         | 解析領域 | 時間領域 |
         | 境界条件 | 水面の境界条件は非線形であるが，非線形のまま解く |

         ### 計算の流れ

         1. 境界条件の設定
         2. 境界値問題（BIE）を解き，$`\phi`$と$`\phi_n`$を求める
         3. 三角形の線形補間を使って節点の流速を計算する
         4. 次時刻の$`\Omega(t+\Delta t)`$がわかるので，修正流速を計算する
         5. 浮体の加速度を計算する．境界値問題（BIE）を解き，$`\phi_t`$と$`\phi_{nt}`$を求め，浮体面上の圧力$`p`$を計算する必要がある
         6. 全境界面の節点の位置を更新．ディリクレ境界では$`\phi`$を次時刻の値へ更新

         */

         int RK_step = 0;
         BEM_BVP BVP;
         do {
            auto RK_time = (*Points.begin())->RK_X.gett();  //%各ルンゲクッタの時刻を使う
            std::cout << "RK_step = " << ++RK_step << "/" << RK_order << ", RK_time = " << RK_time << ", real_time = " << real_time << std::endl;

            setBoundaryTypes(*water, Join(RigidBodyObject, SoftBodyObject));
            std::cout << Green << "setBoundaryTypes" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

            setNeumannVelocity(Join(RigidBodyObject, SoftBodyObject));

            BVP.solve(*water, FMM_BucketsPoints, FMM_BucketsFaces);
            std::cout << Green << "BVP.solve -> {Φ,Φn}が決まる" << Blue << "\nElapsed time: " << Red << watch() << " s\n";

            calculateCurrentVelocities(*water);
            std::cout << Green << "U_BEMを計算" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

            // calculateCurrentUpdateVelocities(*water, 100, (RK_step == 4));
            // if (RK_step == 4)
            //    std::cout << Green << "do shift" << colorOff << " s\n";
            //
            calculateCurrentUpdateVelocities(*water, 50);

            std::cout << Green << "U_update_BEMを計算" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

            BVP.solveForPhiPhin_t(*water, RigidBodyObject);
            std::cout << Green << "BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる" << Blue << "\nElapsed time: " << Red << watch() << colorOff << " s\n";

            // b$ --------------------------------------------------- */

            /*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

            ### 浮体の重心位置・姿勢・速度の更新

            浮体の重心位置は，重心に関する運動方程式を解くことで求める．
            姿勢は，角運動量に関する運動方程式などを使って，各加速度を求める．姿勢はクオータニオンを使って表現する．

            */

            for (const auto &net : RigidBodyObject) {

               // \label{BEM:impose_velocity}
               //  重心位置と姿勢の時間発展
               if (net->inputJSON.find("velocity")) {
                  if (net->inputJSON.at("velocity")[0] != "fixed") {
                     std::cout << "updating " << net->getName() << "'s (RigidBodyObject) velocity" << std::endl;
                     auto COM_old = net->COM;
                     net->RK_COM.push(net->velocityTranslational());
                     net->COM = net->RK_COM.getX();
                     Quaternion q;
                     q = q.d_dt(net->velocityRotational());  // w->クォータニオン
                     net->RK_Q.push(q());                    // クォータニオン->T4dとしてプッシュ
                     net->Q = net->RK_Q.getX();

                     std::cout << "name = " << net->getName() << std::endl;
                     std::cout << "net->COM - COM_old= " << net->COM - COM_old << std::endl;
                  }
                  //
                  // if (net->inputJSON.find("angle") && net->inputJSON.at("angle")[0] == "Hadzic2005") {
                  //    auto axis = stod(net->inputJSON.at("angle")[1]);
                  //    double start = stod(net->inputJSON.at("angle")[2]);
                  //    Hadzic2005 hadzic(start);
                  //    Quaternion q({0., 1., 0.}, hadzic.getAngle(net->RK_Q.get_t()));
                  //    net->Q = q;
                  // }
               }

               bool use_given_velocity = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] != "update" && net->inputJSON.at("velocity")[0] != "floating");
               bool update_velocity_using_predetermined_accel = (net->inputJSON.find("velocity") && net->inputJSON.find("acceleration") && net->inputJSON.at("velocity")[0] == "update");
               bool update_velocity_using_solved_accel = (net->inputJSON.find("acceleration") && net->inputJSON.at("acceleration")[0] == "floating");
               bool update_velocity_using_solved_accel2 = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating");
               if (use_given_velocity) {
                  std::cout << "use " << net->getName() << "'s (RigidBodyObject) predetermiend velocity" << std::endl;
               } else if (update_velocity_using_solved_accel || update_velocity_using_predetermined_accel || update_velocity_using_solved_accel2) {
                  std::cout << "updating " << net->getName() << "'s (RigidBodyObject) velocity from acceleration" << std::endl;
                  net->RK_Velocity.push(net->acceleration);
                  net->velocity = net->RK_Velocity.getX();
                  std::cout << "acceleration = " << net->acceleration << std::endl;
                  std::cout << "velocity = " << net->velocity << std::endl;
               } else {
                  std::cout << net->getName() << "'s (RigidBodyObject) velocity is not updated" << std::endl;
               }

               net->RigidBodyMovePoints();
            }

            // b$ --------------------------------------------------- */

            for (const auto &net : SoftBodyObject) {
               std::cout << "updating " << net->getName() << "'s (SoftBodyObject) position" << std::endl;
               for (const auto &p : net->getPoints()) {
                  p->RK_X.push(p->velocityTranslational());  //@ 位置xの時間発展
                  // p->setXSingle(p->RK_X.getX());
               }
               net->setGeometricProperties();
            }

            // b$ --------------------------------------------------- */

            //@ Φの時間発展，Φnの時間発展はない
            for (const auto &p : Points) {
               if (!p->Neumann /*Neumannを変更しても，あとでBIEによって上書きされるので，からわない．*/) {
                  p->RK_phi.push(p->DphiDt(p->U_update_BEM, 0.));
                  std::get<0>(p->phiphin) = p->phi_Dirichlet = p->RK_phi.getX();  // 角点の法線方向はわからないので，ノイマンの境界条件phinを与えることができない．
               }
               //@ 位置xの時間発展
               p->RK_X.push(p->U_update_BEM);
               p->setXSingle(p->RK_X.getX());
            }
            // b$ --------------------------------------------------- */
            std::cout << Green << "name:" << water->getName() << ": setBounds" << colorOff << std::endl;

            for (auto &nets : {FluidObject, RigidBodyObject, SoftBodyObject})
               for (auto &net : nets)
                  net->setGeometricProperties();

            std::ofstream ofs(output_directory + "/water" + std::to_string(RK_step) + ".obj");
            createOBJ(ofs, *water);
            ofs.close();

            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
         } while (!((*Points.begin())->RK_X.finished));

         /* -------------------------------------------------------------------------- */
         for (const auto &net : SoftBodyObject)
            AreaWeightedSmoothingPreserveShape(net->getPoints(), 1);

         /* ------------------------------------------------------ */

         {
            double mean_phi = 0.;
            for (const auto &p : water->getPoints())
               mean_phi += std::get<0>(p->phiphin);
            mean_phi /= water->getPoints().size();
            for (const auto &p : water->getPoints()) {
               p->phi_Dirichlet -= mean_phi;
               // p->phi_Neumann -= mean_phi;
               std::get<0>(p->phiphin) -= mean_phi;
            }
         }

         /* ------------------------------------------------------ */

         std::cout << Green << "real_timeを取得" << colorOff << std::endl;
         real_time = (*Points.begin())->RK_X.gett();

         for (const auto &p : Points) {
            p->U_BEM_last = p->U_BEM;
            p->U_tangential_BEM_last = p->U_tangential_BEM;
         }

         /* ------------------------------------------------------ */

         std::ofstream ofs(output_directory + "/water_current.obj");
         createOBJ(ofs, *water);
         ofs.close();
         // b# -------------------------------------------------------------------------- */
         // b#                              output JSON files                             */
         // b# -------------------------------------------------------------------------- */
         std::cout << "output JSON files" << std::endl;
         jsonout.push("time", real_time);
         // water
         jsonout.push(water->getName() + "_volume", water->getVolume());
         jsonout.push(water->getName() + "_EK", KinematicEnergy(water->getFaces()));
         jsonout.push(water->getName() + "_EP", PotentialEnergy(water->getFaces()));
         jsonout.push(water->getName() + "_E", TotalEnergy(water->getFaces()));
         // bodies
         for (const auto &net : Join(RigidBodyObject, SoftBodyObject)) {
            if ((net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating") ||
                (net->inputJSON.find("output") && std::ranges::any_of(net->inputJSON.at("output"), [](const auto &s) { return s == "json"; }))) {
               auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
               jsonout.push(net->getName() + "_pitch", net->quaternion.pitch());
               jsonout.push(net->getName() + "_yaw", net->quaternion.yaw());
               jsonout.push(net->getName() + "_roll", net->quaternion.roll());
               jsonout.push(net->getName() + "_force", tmp.surfaceIntegralOfPressure());
               jsonout.push(net->getName() + "_torque", tmp.getFroudeKrylovTorque(net->COM));
               jsonout.push(net->getName() + "_accel", net->acceleration);
               jsonout.push(net->getName() + "_velocity", net->velocity);
               jsonout.push(net->getName() + "_COM", net->COM);
               jsonout.push(net->getName() + "_area", tmp.area);
               jsonout.push(net->getName() + "_EK", Dot(net->velocity, net->velocity) * net->mass / 2.);
               jsonout.push(net->getName() + "_EP", net->mass * _GRAVITY_ * (net->COM[2] - net->ICOM[2]));
            }
         }
         std::ofstream os(output_directory + "/result.json");
         jsonout.output(os);
         os.close();

         // b# -------------------------------------------------------------------------- */
         // b#                            output Paraview files                           */
         // b# -------------------------------------------------------------------------- */

         // 流体
         for (const auto &net : FluidObject) {
            auto filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, net->getFaces(), dataForOutput(net, dt));
            NetOutputInfo[net].PVD->push(filename, real_time);
            NetOutputInfo[net].PVD->output();
         }

         for (const auto &net : Join(RigidBodyObject, SoftBodyObject)) {
            VV_VarForOutput data;
            if (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating") {
               auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
               uomap_P_Tddd P_COM, P_COM_p, P_accel, P_velocity, P_rotational_velocity, P_rotational_accel, P_FroudeKrylovTorque;
               uomap_P_d P_Pitch, P_Yaw, P_Roll, P_pressure, P_radius;
               //    for (const auto &p : water->getPoints()) {
               //       P_accel[p] = net->accelRigidBody(ToX(p));
               //       P_velocity[p] = net->velocityRigidBody(ToX(p));
               //       P_rotational_velocity[p] = net->velocityRotational();
               //       // P_FroudeKrylovTorque[p] = tmp.getFroudeKrylovTorque(net->COM);
               //       P_pressure[p] = p->pressure;
               //    }

               for (const auto &p : net->getPoints()) {
                  P_COM[p] = net->COM;
                  P_COM_p[p] = net->COM - ToX(p);
                  P_Pitch[p] = net->quaternion.pitch();
                  P_Yaw[p] = net->quaternion.yaw();
                  P_Roll[p] = net->quaternion.roll();
                  P_accel[p] = net->accelRigidBody(ToX(p));
                  P_velocity[p] = net->velocityRigidBody(ToX(p));
                  P_rotational_velocity[p] = net->velocityRotational();
               }

               data = {
                   {"vector to COM", P_COM_p},
                   {"COM", P_COM},
                   {"pitch", P_Pitch},
                   {"yaw", P_Yaw},
                   {"roll", P_Roll},
                   {"velocity", P_accel},
                   {"acceleration", P_velocity},
                   {"rotational valocity", P_rotational_velocity},
                   {"rotational acceleration", P_rotational_accel},
                   {"pressure", P_pressure},
                   {"FroudeKrylovTorque", P_FroudeKrylovTorque}};
               mk_vtu(output_directory + "/actingFacesOn" + NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu", tmp.actingFaces, data);
            }
            auto filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, net->getFaces(), data);
            NetOutputInfo[net].PVD->push(filename, real_time);
            NetOutputInfo[net].PVD->output();
         }
         // 流体
         // {
         //    std::vector<Tddd> points;
         //    for (const auto &p : water->getPoints())
         //       if (p->CORNER)
         //          points.emplace_back(ToX(p));

         //    auto filename = "cornerPoints" + std::to_string(time_step) + ".vtu";
         //    mk_vtu(output_directory + "/" + filename, points);
         //    cornerPointsPVD.push(filename, real_time);
         //    cornerPointsPVD.output();
         // }
         // {
         //    std::unordered_set<networkFace *> faces;
         //    for (const auto &f : water->getFaces()) {
         //       // for (const auto p : f->getPoints())
         //       // 	if (p->CORNER)
         //       // 	{
         //       // 		faces.emplace(f);
         //       // 		break;
         //       // 	}
         //       for_each(f->getPoints(), [&](const auto &p) {
         //          if (p->CORNER)
         //             faces.emplace(f);
         //       });
         //    }
         //    auto filename = "corner" + std::to_string(time_step) + ".vtu";
         //    mk_vtu(output_directory + "/" + filename, faces, data);
         //    cornerPVD.push(filename, real_time);
         //    cornerPVD.output();
         // }
         // b# -------------------------------------------------------------------------- */
         //
         // mk_vtu(output_directory + "/" + obj.getName() + std::to_string(time_step) + ".vtu", obj.getFaces(), datacpg);
         // mk_vtu(output_directory + "/" + name + std::to_string(time_step) + ".vtu", cpg.well_Faces, datacpg);
         // cpg_pvd.push(name, name + std::to_string(time_step) + ".vtu", time_step * t_rep * dt);
         // cpg_pvd.output_();
         // b* ------------------------------------------------------ */
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
   return 0;
};

/*DOC_EXTRACT 2_0_HOW_TO_RUN

# 実行方法

ファイルをダウンロードして，`build_bem`ディレクトリに移動．

```shell
$ git clone https://github.com/tomoakihirakawa/cpp.git
$ cd ./cpp/builds/build_bem
```

`clean`でCMake関連のファイルを削除して（ゴミがあるかもしれないので），
`cmake`で`Makefile`を生成して，`make`でコンパイルする．

```shell
$ sh clean
$ cmake -DCMAKE_BUILD_TYPE=Release ../
$ make
```

次に，入力ファイルを生成．

```shell
$ python3 input_generator.py
```

例えば，`./input_files/Hadzic2005`が生成される．入力ファイルを指定して実行．

```shell
$ ./main ./input_files/Hadzic2005
```

*/

/*DOC_EXTRACT 3_0_EXAMPLES

# Examples

**[See the Examples here!](EXAMPLES.md)**

*/
