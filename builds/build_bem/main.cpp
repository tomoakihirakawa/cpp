// #define _debugging_

#define BEM
#define use_lapack
int time_step;
double real_time = 0;

#define simulation
#include <filesystem>
#include <regex>
#include "Network.hpp"
#include "integrationOfODE.hpp"
#include "kernelFunctions.hpp"
#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"

pvd cpg_pvd("./vtu/bem.pvd");

#include "BEM.hpp"
#include "svd.hpp"

int main(int arg, char **argv) {
   if (arg <= 1)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "argv <= 1. write input json file directory!\nex.\n$ ./main ./input");
   std::string input_directory(argv[1]);  // input directory
   input_directory += "/";
   std::cout << input_directory << std::endl;
   try {
      /* -------------------------------------------------------------------------- */
      /*                           read setting.json                                */
      /* -------------------------------------------------------------------------- */
      // Read settings from a JSON file and store values
      JSON settingJSON(input_directory + "/setting.json");
      for (auto &[key, value] : settingJSON())
         std::cout << key << ": " << value << std::endl;

      // Check for the existence of required keys in the JSON file
      const std::vector<std::string> required_keys = {"end_time_step", "end_time", "output_directory", "max_dt", "input_files"};
      for (const auto &key : required_keys)
         if (!settingJSON.find(key))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "setting.json does not have " + key);

      // Store values from the JSON file into variables
      std::string output_directory = settingJSON.at("output_directory")[0];
      double max_dt = stod(settingJSON.at("max_dt"))[0];
      const int end_time_step = stoi(settingJSON.at("end_time_step"))[0];
      const double end_time = stod(settingJSON.at("end_time"))[0];

      // Set default values for optional keys and update if found in the JSON file
      double stop_remesh_time = 1E+10;
      double force_remesh_time = 0;
      int grid_refinement = 0;
      if (settingJSON.find("stop_remesh_time"))
         stop_remesh_time = stod(settingJSON.at("stop_remesh_time")[0]);
      if (settingJSON.find("force_remesh_time"))
         force_remesh_time = stod(settingJSON.at("force_remesh_time")[0]);
      if (settingJSON.find("grid_refinement"))
         grid_refinement = stoi(settingJSON.at("grid_refinement")[0]);

      /* -------------------------------------------------------------------------- */
      /*                read JSON files, input_files in setting.json                */
      /* -------------------------------------------------------------------------- */
      // Define containers and helper functions for network objects and output information
      std::map<Network *, outputInfo> NetOutputInfo;
      auto setup_network_output_info = [&NetOutputInfo](Network *net, const JSON &injson, const std::string &output_directory) {
         auto output_name = injson.at("name")[0];
         NetOutputInfo[net].pvd_file_name = output_name;
         NetOutputInfo[net].vtu_file_name = output_name + "_";
         NetOutputInfo[net].PVD = new PVDWriter(output_directory + "/" + output_name + ".pvd");
      };

      std::vector<Network *> FluidObject, RigidBodyObject, SoftBodyObject;
      auto initialize_network_objects = [&FluidObject, &RigidBodyObject, &SoftBodyObject](Network *net, const JSON &injson) {
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
         net->resetInitialX();
         net->setGeometricProperties();
      };

      for (auto input_file_name : settingJSON["input_files"]) {
         std::cout << input_directory + input_file_name << std::endl;
         JSON injson(input_directory + input_file_name);

         const std::vector<std::string> required_keys = {"name", "objfile", "type"};
         for (const auto &key : required_keys)
            if (!injson.find(key))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, input_directory + input_file_name + " does not have " + key);

         auto object_name = injson["name"][0];

         for (auto &[key, value] : injson())
            std::cout << key << ": " << value << std::endl;

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

      /* -------------------------------------------------------------------------- */

      // Create output directory and copy files
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
      //  b* ------------------------------------------------------ */
      //  b*                         メインループ                      */
      //  b* ------------------------------------------------------ */
      TimeWatch watch;
      for (time_step = 0; time_step < end_time_step; time_step++) {
         if (end_time < real_time)
            break;
         show_info(*water);
         //! 体積を保存するようにリメッシュする必要があるだろう．
         // auto radius = Mean(extLength(water->getLines()));
         setBoundaryConditions(*water, Join(RigidBodyObject, SoftBodyObject));
         double rad = M_PI / 180;
         // flipIf(*water, {10 * rad, rad}, {5 * rad /*結構小さく*/, rad}, false);
         flipIf(*water, {10 * rad, rad}, {10 * rad /*結構小さく*/, rad}, false);
         // b# ------------------------------------------------------ */
         // b#                       刻み時間の決定                     */
         // b# ------------------------------------------------------ */
         const auto Points = water->getPoints();
         const auto Faces = water->getFaces();
         double dt = max_dt;
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

         mk_vtu(output_directory + "/FMM_BucketsFaces.vtu", FMM_BucketsFaces.getT4Tddd());
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
         }
         // b$ --------------------------------------------------- */
         for (const auto &net : SoftBodyObject) {
            // !いらないはずのもの
            for (const auto &p : net->getPoints())
               p->RK_X.initialize(dt, real_time, ToX(p), RK_order);
         }

         int RK_step = 0;
         BEM_BVP BVP;
         do {
            //! 壁面の動きは，マイステップ更新することにした．この結果はphin()で参照される
            auto RK_time = (*Points.begin())->RK_X.gett();  //%各ルンゲクッタの時刻を使う
            std::cout << "RK_step = " << ++RK_step << "/" << RK_order << ", RK_time = " << RK_time << ", real_time = " << real_time << std::endl;
            // b# ------------------------------------------------------ */
            // b#      物体のノイマン境界の速度 u(t) at Neumann を設定         */
            // b# ------------------------------------------------------ */
            for (const auto &net : RigidBodyObject) {
               std::cout << "----------------" << std::endl;
               std::cout << net->getName() << "　の流速の計算方法" << std::endl;
               if (net->isFixed) {
                  net->mass = 1E+20;
                  net->inertia.fill(1E+20);
                  net->COM.fill(0.);
                  net->initial_center_of_mass.fill(0.);
               }

               if (net->inputJSON.find("velocity")) {
                  std::string move_name = net->inputJSON["velocity"][0];
                  std::cout << "move_name = " << move_name << std::endl;
                  if (move_name == "fixed") {
                     net->velocity.fill(0.);
                     net->acceleration.fill(0.);
                  } else if (move_name != "floating") {
                     net->velocity = velocity(move_name, net->inputJSON["velocity"], RK_time);  // T6d //@ Φnを計算するために，物体表面の速度forced_velocityは，保存しておく必要がある
                                                                                                // net->acceleration = forced_motion::acceleration(RK_time); // T6d //@ 圧力を計算するために，物体表面の加速度は，保存しておく必要がある
                  } else if (move_name == "floating") {
                     std::cout << "floatingの場合は，加速度の時間積分によってシミュレートされる" << std::endl;
                  }
               } else {
                  std::cout << "指定がないので速度はゼロ" << std::endl;
                  net->velocity.fill(0.);
                  net->acceleration.fill(0.);
               }
               std::cout << "----------------" << std::endl;
            }
            // b$ --------------------------------------------------- */
            for (const auto &net : SoftBodyObject) {
               std::cout << "----------------" << std::endl;
               std::cout << net->getName() << "　の流速の計算方法．soft bodyの場合，各節点に速度を与える．" << std::endl;
               net->velocity.fill(0.);
               net->acceleration.fill(0.);
               if (net->inputJSON.find("velocity")) {
                  std::string move_name = net->inputJSON["velocity"][0];
                  std::cout << "move_name = " << move_name << std::endl;
                  if (move_name == "fixed") {
                     for (const auto &p : net->getPoints()) {
                        p->velocity.fill(0.);
                        p->acceleration.fill(0.);
                     }
                  } else {
                     for (const auto &p : net->getPoints())
                        p->velocity = velocity(move_name, net->inputJSON["velocity"], p, RK_time);  // T6d //@ Φnを計算するために，物体表面の速度forced_velocityは，保存しておく必要がある
                  }
               } else {
                  std::cout << "指定がないので速度はゼロ" << std::endl;
                  for (const auto &p : net->getPoints()) {
                     p->velocity.fill(0.);
                     p->acceleration.fill(0.);
                  }
               }
               std::cout << "----------------" << std::endl;
            }
            // やはり，ルンゲクッタは3回目に同じ場所での計算を行うわけなので，角点に関しては，毎回の時間発展＆修正はまずいだろう．
            // 予測するΩ(t+δt)は，ルンゲクッタに合わせたものでないとおかしいだろう．つまり，RK4なら　Ω(t0),Ω(t1),Ω(t2=t1),Ω(t3)でないといけないだろう．
            // b% -------------------------------------------------------- */
            // b%            境界条件（角点・ディリクレ・ノイマン）の決定               */
            // b% -------------------------------------------------------- */
            setBoundaryConditions(*water, Join(RigidBodyObject, SoftBodyObject));
            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
            // b! ------------------------------------------------------ */
            // b!           　境界値問題を解く-> {Φ,Φn}が決まる                */
            // b! ------------------------------------------------------ */
            std::cout << Green << "境界値問題を解く-> {Φ,Φn}が決まる" << colorOff << std::endl;
            BVP.solve(*water, FMM_BucketsPoints, FMM_BucketsFaces);
            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
            // b* ------------------------------------------------------ */
            // b*                    微分∇ΦやDUDtを計算                    */
            // b* ------------------------------------------------------ */
            std::cout << Green << "微分∇ΦやDUDtを計算" << colorOff << std::endl;
            calculateVelocities(*water);
            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
            // b* ------------------------------------------------------ */
            // b*           　境界値問題を解く-> {Φt,Φtn}が決まる               */
            // b* ------------------------------------------------------ */
            std::cout << Green << "境界値問題を解く-> {Φt,Φtn}が決まる" << colorOff << std::endl;
            BVP.solveForPhiPhin_t(water, RigidBodyObject);
            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
            // @ -------------------------------------------------------- */
            // @                 ディリクレ境界ではΦを時間積分                 */
            // @ -------------------------------------------------------- */
            std::cout << Green << "ディリクレ境界ではΦを時間積分，ノイマン境界ではΦnを陽に与える" << colorOff << std::endl;
            // b$ --------------------------------------------------- */
            for (const auto &net : RigidBodyObject) {
               std::cout << "name:" << net->getName() << std::endl;
               if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] != "fixed") {
                  net->RK_COM.push(net->velocityTranslational());
                  net->COM = net->RK_COM.getX();
                  Quaternion q;
                  q = q.d_dt(net->velocityRotational());  // w->クォータニオン
                  net->RK_Q.push(q());                    // クォータニオン->T4dとしてプッシュ
                  net->Q = net->RK_Q.getX();
               }
               if (!net->inputJSON.find("acceleration") || (net->inputJSON.find("acceleration") && net->inputJSON["acceleration"][0] != "fixed")) {
                  net->RK_Velocity.push(net->acceleration);
                  net->velocity = net->RK_Velocity.getX();
               }
               net->RigidBodyMovePoints();
            }
            // b$ --------------------------------------------------- */
            for (const auto &net : SoftBodyObject) {
               std::cout << "name:" << net->getName() << std::endl;
               for (const auto &p : net->getPoints()) {
                  p->RK_X.push(p->velocityTranslational());  //@ 位置xの時間発展
                  // p->setXSingle(p->RK_X.getX());
               }
               net->setGeometricProperties();
            }
            // b$ --------------------------------------------------- */
            for (const auto &p : Points) {
               //@ Φの時間発展，Φnの時間発展はない
               if (!p->Neumann /*Neumannを変更しても，あとでBIEによって上書きされるので，何からわない．*/) {
                  p->RK_phi.push(p->DphiDt(p->U_update_BEM, 0.));
                  std::get<0>(p->phiphin) = p->phi_Dirichlet = p->RK_phi.getX();  // 角点の法線方向はわからないので，ノイマンの境界条件phinを与えることができない．
               }
               //@ 位置xの時間発展
               p->RK_X.push(p->U_update_BEM);
               p->setXSingle(p->RK_X.getX());
            }
            // b$ --------------------------------------------------- */
            std::cout << Green << "name:" << water->getName() << ": setBounds" << colorOff << std::endl;
            water->setGeometricProperties();

            std::ofstream ofs(output_directory + "/water" + std::to_string(RK_step) + ".obj");
            creteOBJ(ofs, *water);
            ofs.close();

            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
         } while (!((*Points.begin())->RK_X.finished));

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
         creteOBJ(ofs, *water);
         ofs.close();

         // b# -------------------------------------------------------------------------- */
         // b#                              output JSON files                             */
         // b# -------------------------------------------------------------------------- */
         jsonout.push("time", real_time);
         jsonout.push(water->getName() + "_volume", water->getVolume());
         for (const auto &net : Join(RigidBodyObject, SoftBodyObject)) {
            if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] == "floating") {
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
            }
         }
         std::ofstream os(output_directory + "/result.json");
         jsonout.output(os);
         os.close();

         // b# -------------------------------------------------------------------------- */
         // b#                            output Paraview files                           */
         // b# -------------------------------------------------------------------------- */

         auto data = dataForOutput(*water, dt);
         // 流体
         for (const auto &net : FluidObject) {
            auto filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, net->getFaces(), data);
            NetOutputInfo[net].PVD->push(filename, real_time);
            NetOutputInfo[net].PVD->output();
         }

         for (const auto &net : Join(RigidBodyObject, SoftBodyObject)) {
            VV_VarForOutput data;
            if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] == "floating") {
               auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
               //    uomap_P_Tddd P_COM, P_COM_p, P_accel, P_velocity, P_rotational_velocity, P_rotational_accel, P_FroudeKrylovTorque;
               //    uomap_P_d P_Pitch, P_Yaw, P_Roll, P_pressure;

               //    for (const auto &p : water->getPoints()) {
               //       P_accel[p] = net->accelRigidBody(ToX(p));
               //       P_velocity[p] = net->velocityRigidBody(ToX(p));
               //       P_rotational_velocity[p] = net->velocityRotational();
               //       // P_FroudeKrylovTorque[p] = tmp.getFroudeKrylovTorque(net->COM);
               //       P_pressure[p] = p->pressure;
               //    }

               //    for (const auto &p : net->getPoints()) {
               //       P_COM[p] = net->COM;
               //       P_COM_p[p] = net->COM - ToX(p);
               //       P_Pitch[p] = net->quaternion.pitch();
               //       P_Yaw[p] = net->quaternion.yaw();
               //       P_Roll[p] = net->quaternion.roll();
               //       P_accel[p] = net->accelRigidBody(ToX(p));
               //       P_velocity[p] = net->velocityRigidBody(ToX(p));
               //       P_rotational_velocity[p] = net->velocityRotational();
               //    }

               //    data = {
               //        {"vector to COM", P_COM_p},
               //        {"COM", P_COM},
               //        {"pitch", P_Pitch},
               //        {"yaw", P_Yaw},
               //        {"roll", P_Roll},
               //        {"velocity", P_accel},
               //        {"acceleration", P_velocity},
               //        {"rotational valocity", P_rotational_velocity},
               //        {"rotational acceleration", P_rotational_accel},
               //        {"pressure", P_pressure},
               //        {"FroudeKrylovTorque", P_FroudeKrylovTorque}};
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