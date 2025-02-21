#pragma once

#include "BEM.hpp"
#include "basic.hpp"

std::string toLowerCase(const std::string &str) {
   std::string strLower = str;
   std::transform(strLower.begin(), strLower.end(), strLower.begin(),
                  [](unsigned char c) { return std::tolower(c); });
   return strLower;
}

struct SimulationSettings {
   const std::string setting_filename = "setting.json";
   std::filesystem::path input_directory;
   JSON settingJSON;
   //
   std::filesystem::path output_directory;
   int end_time_step;
   double end_time;
   double max_dt;
   int ALEPERIOD;
   double stop_remesh_time;
   double force_remesh_time;
   int grid_refinement;
   bool _LINEAR_ELEMENT_ = true;
   bool _PSEUDO_QUADRATIC_ELEMENT_ = false;
   bool _ALE_ON_LINEAR_ELEMENT_ = false;
   bool _ALE_ON_PSEUDO_QUADRATIC_ELEMENT_ = false;
   //    double _WATER_DENSITY_ = 1000;
   //    double _GRAVITY_ = 9.81;

   std::map<Network *, outputInfo> NetOutputInfo;

   std::vector<Network *> FluidObject, RigidBodyObject, SoftBodyObject, AbsorberObject;
   std::vector<JSON> MeasurementJSONs;

   /* -------------------------------------------------------------------------- */

   SimulationSettings(std::filesystem::path input_directory)
       : input_directory(input_directory),
         settingJSON(input_directory / setting_filename),
         output_directory(settingJSON.at("output_directory")[0]),
         end_time_step(stoi(settingJSON.at("end_time_step"))[0]),
         end_time(stod(settingJSON.at("end_time"))[0]),
         max_dt(stod(settingJSON.at("max_dt"))[0]),
         ALEPERIOD(stoi(settingJSON.at("ALEPERIOD"))[0]) {
      {

         std::cout << "input_directory : " << input_directory << std::endl;

         /* -------------------------------------------------------------------------- */
         /*                           read setting.json                                */
         /* -------------------------------------------------------------------------- */

         for (auto &[key, value] : settingJSON()) std::cout << key << ": " << value << std::endl;

         for (const auto &key : {"end_time_step", "end_time", "output_directory", "max_dt", "input_files"}) {
            if (!settingJSON.find(key)) {
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Missing key in setting.json: " + std::string(key));
            }
         }

         if (settingJSON.find("element") && toLowerCase(settingJSON.at("element")[0]) == "linear")
            _LINEAR_ELEMENT_ = true;
         else if (settingJSON.find("element") && (toLowerCase(settingJSON.at("element")[0]).contains("quad") && toLowerCase(settingJSON.at("element")[0]).contains("pseudo")))
            _PSEUDO_QUADRATIC_ELEMENT_ = true;
         else
            _LINEAR_ELEMENT_ = true;

         if (_LINEAR_ELEMENT_)
            std::cout << "LINEAR_ELEMENT" << std::endl;
         if (_PSEUDO_QUADRATIC_ELEMENT_)
            std::cout << "PSEUDO_QUADRATIC_ELEMENT" << std::endl;

         if (settingJSON.find("ALE") && toLowerCase(settingJSON.at("ALE")[0]) == "linear")
            _ALE_ON_LINEAR_ELEMENT_ = true;
         else if (settingJSON.find("ALE") && (toLowerCase(settingJSON.at("ALE")[0]).contains("quad") || toLowerCase(settingJSON.at("ALE")[0]).contains("pseudo")))
            _ALE_ON_PSEUDO_QUADRATIC_ELEMENT_ = true;
         else
            _ALE_ON_PSEUDO_QUADRATIC_ELEMENT_ = true;

         if (_ALE_ON_LINEAR_ELEMENT_)
            std::cout << "ALE_ON_LINEAR_ELEMENT" << std::endl;
         if (_ALE_ON_PSEUDO_QUADRATIC_ELEMENT_)
            std::cout << "ALE_ON_PSEUDO_QUADRATIC_ELEMENT" << std::endl;

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

         for (auto input_file_name : settingJSON["input_files"]) {
            std::cout << Green << input_directory / input_file_name << colorReset << std::endl;
            JSON injson(input_directory / input_file_name);

            /* -------------------- display contents of the JSON file ------------------- */
            for (auto &[key, value] : injson()) std::cout << Green << key << colorReset << ": " << value << std::endl;
            auto object_name = injson.at("name")[0];

            /* --------------------------- check required keys -------------------------- */
            if (!injson.find("name"))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"name\" key");
            if (!injson.find("type"))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"type\" key");
            auto type = injson.at("type")[0];
            if ((type.contains("Fluid") || type.contains("Body")) && !injson.find("objfile"))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"objfile\" key");
            else if (type.contains("Measurement") || type.contains("gauge")) {
               MeasurementJSONs.emplace_back(injson);
               std::cout << "type = " << type << std::endl;
               std::cout << "skipped" << std::endl;
               continue;
            }

            /* ------------------------ create Network object --------------------------- */

            if (!injson.find("ignore") || !stob(injson["ignore"])[0]) {

               auto net = new Network(injson.at("objfile")[0], object_name);
               net->inputJSON = injson;
               net->applyTransformations(injson);
               setOutputInfo(net, injson.at("name")[0], output_directory);

               /* ----------------------------- 境界面データの回転と平行移動 ----------------------------- */
               if (injson.find("rotation")) {
                  const auto angle_theta = stod(injson.at("rotation"));
                  std::array<double, 3> axis = {angle_theta[0], angle_theta[1], angle_theta[2]};
                  const auto angle = angle_theta[3];
                  net->rotate(angle, axis);
               }
               if (injson.find("translation")) {
                  const auto translation = stod(injson.at("translation"));
                  net->translate(std::array<double, 3>{translation[0], translation[1], translation[2]});
               }

               setTypes(net);
               std::filesystem::copy_file(input_directory / input_file_name, output_directory / input_file_name, std::filesystem::copy_options::overwrite_existing);
               mk_vtu(output_directory / (object_name + "_init.vtu"), {net->getFaces()});
               /* -------------------------------------------------------------------------- */
               for (auto &[key, value] : injson()) {
                  if (key.contains("mooring")) {
                     //$ ------------------------ MOORING ---------------------- */
                     //$ extract key the contains "mooring" as key name. extract them as vector
                     std::cout << "/* ------------------------ MOORING ---------------------- */" << std::endl;
                     std::cout << "initialize mooring" << std::endl;

                     //! check if the value contains 12 elements
                     if (value.size() != 13)  //! show contents of the value nicely
                        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mooring line must have 12 elements");

                     const auto name = value[0];
                     const int n_points = std::stoi(value[8]);  //$ number of points
                     const std::array<double, 3> X_begin = {std::stod(value[1]), std::stod(value[2]), std::stod(value[3])};
                     const std::array<double, 3> X_end = {std::stod(value[4]), std::stod(value[5]), std::stod(value[6])};
                     const double total_length = std::stod(value[7]);  //! [m]
                     const double w = std::stod(value[9]);             //! [kg/m]
                     const double stiffness = std::stod(value[10]);    //! [N/m]
                     const double damp = std::stod(value[11]);         //! [N/(m/s^2)]
                     const double diam = std::stod(value[12]);         //! [m]

                     std::cout << std::right << std::setw(15) << "name : " << name << std::endl;
                     std::cout << std::right << std::setw(15) << "X_begin : " << X_begin << std::endl;
                     std::cout << std::right << std::setw(15) << "X_end : " << X_end << std::endl;
                     std::cout << std::right << std::setw(15) << "total_length : " << total_length << std::endl;
                     std::cout << std::right << std::setw(15) << "n_points : " << n_points << std::endl;
                     std::cout << std::right << std::setw(15) << "mass per unit length : " << w << std::endl;
                     std::cout << std::right << std::setw(15) << "stiffness : " << stiffness << std::endl;
                     std::cout << std::right << std::setw(15) << "damp : " << damp << std::endl;
                     std::cout << std::right << std::setw(15) << "diam : " << diam << std::endl;
                     std::cout << std::right << std::setw(15) << "total mass" << w * total_length << std::endl;

                     std::cout << std::right << std::setw(15) << "initialize MooringLine" << std::endl;
                     auto mooring_net = new MooringLine(X_begin, X_end, total_length, n_points);
                     mooring_net->setName(name);

                     std::cout << std::right << std::setw(15) << "MooringLine->getPoints().size() = " << mooring_net->getPoints().size() << std::endl;
                     std::cout << std::right << std::setw(15) << "MooringLine->getName() = " << mooring_net->getName() << std::endl;

                     std::cout << std::right << std::setw(15) << "MooringLine->setDensityStiffnessDampingDiameter" << std::endl;
                     mooring_net->setDensityStiffnessDampingDiameter(w, stiffness, damp, diam);

                     net->mooringLines.emplace_back(mooring_net);
                     setOutputInfo(mooring_net, name, output_directory);

                     std::cout << std::right << std::setw(15) << "MooringLine->setEquilibriumState" << std::endl;
                     auto boundary_condition = [&](networkPoint *p) {
                        if (p == mooring_net->lastPoint || p == mooring_net->firstPoint) {
                           p->acceleration.fill(0.);
                           p->velocity.fill(0.);
                        }
                     };

                     mooring_net->setEquilibriumState(boundary_condition);

                     std::cout << "/* ------------------------------------------------------- */" << std::endl;
                  }
               }
               //$ ------------------------------------------------------- */
            } else
               Print("skipped");
         }
      }
   }

   /* -------------------------------------------------------------------------- */

   void setOutputInfo(Network *net, auto output_name, const std::filesystem::path &output_directory) {
      std::cout << "setOutputInfo" << std::endl;
      this->NetOutputInfo[net].pvd_file_name = output_name;
      this->NetOutputInfo[net].vtu_file_name = output_name + "_";
      this->NetOutputInfo[net].PVD = new PVDWriter(output_directory / (output_name + ".pvd"));
   };

   /* -------------------------------------------------------------------------- */

   void setTypes(Network *net) {
      std::cout << "setTypes" << std::endl;
      //$ set type
      auto type = net->inputJSON.at("type")[0];
      if (type.contains("RigidBody")) {
         std::cout << "RigidBody" << std::endl;
         RigidBodyObject.emplace_back(net);
         net->isRigidBody = true;
         net->isSoftBody = net->isFluid = false;
      } else if (type.contains("SoftBody") || type.contains("FixedBody")) {
         std::cout << "SoftBody" << std::endl;
         SoftBodyObject.emplace_back(net);
         net->isSoftBody = true;
         net->isRigidBody = net->isFluid = false;
      } else if (type.contains("Fluid")) {
         std::cout << "Fluid" << std::endl;
         FluidObject.emplace_back(net);
         net->isFluid = true;
         net->isRigidBody = net->isSoftBody = false;
      } else if (type.contains("Absorber") || type.contains("absorb") || type.contains("damping")) {
         std::cout << "Absorber" << std::endl;
         AbsorberObject.emplace_back(net);
         net->isAbsorber = true;
         net->isRigidBody = net->isSoftBody = net->isFluid = false;
         //
         net->inputJSON.find("wave_theory", [&](auto STR_VEC) {
            net->water_wave_theory = WaterWaveTheory();
            auto vec = stod(STR_VEC);
            net->water_wave_theory.A = vec[0];
            net->water_wave_theory.set_T_h(vec[1], vec[2]);
            net->water_wave_theory.bottom_z = vec[3];
            if (vec.size() > 4)
               net->water_wave_theory.theta = vec[4] / 180. * M_PI;
            net->absorb_velocity = [net](const networkPoint *p) { return net->water_wave_theory.gradPhi(p->X, p->RK_X.gett()); };
            net->absorb_gradPhi_t = [net](const networkPoint *p) { return net->water_wave_theory.gradPhi_t(p->X, p->RK_X.gett()); };
            net->absorb_phi = [net](const Tddd &X, const double t) { return net->water_wave_theory.phi(X, t); };
            net->absorb_eta = [net](const Tddd &X, const double t) { return net->water_wave_theory.eta(X, t); };
            net->absorb_gamma = [net](const networkPoint *p) { return std::clamp(p->signed_distance / (2 * net->water_wave_theory.L), 0., 1.); };
         });
         net->inputJSON.find("random_wave_theory", [&](auto STR_VEC) {
            if (STR_VEC.size() > 4)
               throw std::runtime_error("random_wave_theory must have 4 elements");
            auto vec = stod(STR_VEC);
            double H13 = vec[0];
            double T13 = vec[1];
            double h = vec[2];
            double bottom_z = vec[3];
            net->random_water_wave_theory = RandomWaterWaveTheory(H13, T13, h, bottom_z);
            net->absorb_velocity = [net](const networkPoint *p) { return net->random_water_wave_theory.gradPhi(p->X, p->RK_X.gett()); };
            net->absorb_gradPhi_t = [net](const networkPoint *p) { return net->random_water_wave_theory.gradPhi_t(p->X, p->RK_X.gett()); };
            net->absorb_phi = [net](const Tddd &X, const double t) { return net->random_water_wave_theory.phi(X, t); };
            net->absorb_eta = [net](const Tddd &X, const double t) { return net->random_water_wave_theory.eta(X, t); };
            net->absorb_gamma = [net](const networkPoint *p) { return std::clamp(p->signed_distance / (2 * net->random_water_wave_theory.L13), 0., 1.); };
         });
         net->inputJSON.find("wave_theory_L", [&](auto STR_VEC) {
            net->water_wave_theory = WaterWaveTheory();
            auto vec = stod(STR_VEC);
            net->water_wave_theory.A = vec[0];
            net->water_wave_theory.set_L_h(vec[1], vec[2]);
            net->water_wave_theory.bottom_z = vec[3];
            if (vec.size() > 4)
               net->water_wave_theory.theta = vec[4] / 180. * M_PI;
            net->absorb_velocity = [net](const networkPoint *p) { return net->water_wave_theory.gradPhi(p->X, p->RK_X.gett()); };
            net->absorb_gradPhi_t = [net](const networkPoint *p) { return net->water_wave_theory.gradPhi_t(p->X, p->RK_X.gett()); };
            net->absorb_phi = [net](const Tddd &X, const double t) { return net->water_wave_theory.phi(X, t); };
            net->absorb_eta = [net](const Tddd &X, const double t) { return net->water_wave_theory.eta(X, t); };
            net->absorb_gamma = [net](const networkPoint *p) { return std::clamp(p->signed_distance / (2 * net->water_wave_theory.L), 0., 1.); };
         });
      }

      //! isFixedはdefaultでfalse．指定された分だけ順に置き換わる
      //! ただし，指定が１つだけなら，それを全てのに適用する．
      if (net->inputJSON.find("isFixed")) {
         auto isFixed = stob(net->inputJSON.at("isFixed"));
         if (isFixed.size() == 1)
            net->isFixed.fill(isFixed[0]);
         else
            for (auto i = 0; i < isFixed.size(); ++i)
               net->isFixed[i] = isFixed[i];
      }
      // for (auto i = 0; i < 10; ++i)
      //    AreaWeightedSmoothingPreserveShape(net->getPoints(), 0.1);
      //$ set velocity
      std::cout << "set velocity" << std::endl;
      net->isFloatingBody = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating");
      // velocityにfileが指定されている場合は，そのファイルを読み込み，
      // interpolationBsplineであるnet->intpMotionRigidBodyをsetする．
      net->inputJSON.find("velocity", [&](auto STR_VEC) {
         if (STR_VEC[0].contains("file")) {
            if (STR_VEC.size() == 1)
               throw std::runtime_error("Failed to open the input file.");
            std::ifstream file(STR_VEC[1]);
            if (!file.is_open())
               throw std::runtime_error("Failed to open the input file.");
            std::vector<double> T;
            std::vector<std::array<double, 6>> XYZ_Angles;
            std::string line;
            double t, x, y, z, q0, q1, q2;
            while (std::getline(file, line)) {
               if (line[0] == '#') continue;
               replace(line.begin(), line.end(), ',', ' ');
               std::istringstream iss(line);
               iss >> t >> x >> y >> z >> q0 >> q1 >> q2;
               T.push_back(t);
               XYZ_Angles.push_back({x, y, z, q0, q1, q2});
            }
            file.close();
            net->intpMotionRigidBody.set(3, T, XYZ_Angles);
         }
      });
      net->isFloatingBody = net->isFloatingBody || (net->inputJSON.find("acceleration") && net->inputJSON.at("acceleration")[0] == "floating");

      net->resetInitialX();
      net->setGeometricProperties();
   };
};