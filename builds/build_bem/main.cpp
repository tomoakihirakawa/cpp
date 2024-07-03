/*DOC_EXTRACT 0_0_BEM

# BEM-MEL

[README_FOR_STUDENTS.md](README_FOR_STUDENTS.md)

[REVIEW_NOTE0.md](REVIEW_NOTE0.md)

[REVIEW_NOTE1.md](REVIEW_NOTE1.md)

*/

// #define _debugging_

bool _LINEAR_ELEMENT_ = false;
bool _PSEUDO_QUADRATIC_ELEMENT_ = false;
bool _ALE_ON_LINEAR_ELEMENT_ = false;
bool _ALE_ON_PSEUDO_QUADRATIC_ELEMENT_ = false;

#define BEM
#define use_lapack
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
#include "vtkWriter.hpp"

int time_step;
double simulation_time = 0;
JSONoutput jsonout;

// pvd cpg_pvd("./vtu/bem.pvd");

#include "BEM.hpp"
#include "svd.hpp"

std::string toLowerCase(const std::string &str) {
   std::string strLower = str;
   std::transform(strLower.begin(), strLower.end(), strLower.begin(),
                  [](unsigned char c) { return std::tolower(c); });
   return strLower;
}

int main(int argc, char **argv) {

   std::clock_t cpu_clock_start = std::clock();
   auto wall_clock_start = std::chrono::high_resolution_clock::now();

   /*DOC_EXTRACT 0_1_BEM

   ## 入力ファイルの読み込み

   1. 境界条件の設定
   2. 境界値問題（BIE）を解き，$`\phi`$と$`\phi_n`$を求める
   3. 三角形の線形補間を使って節点の流速を計算する

   */

   /* -------------------------------------------------------------------------- */
   /*                           Set up logging to file                           */
   /* -------------------------------------------------------------------------- */
   if (!initializeLogFile("log.txt", argc, argv))
      return 1;

   /* -------------------------------------------------------------------------- */

   if (argc <= 1)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "argv <= 1. write input json file directory!\nex.\n$ ./main ./input");

   std::filesystem::path input_directory(argv[1]);  // input directory
   // input_directory += "/";
   std::cout << input_directory << std::endl;

   /* -------------------------------------------------------------------------- */
   /*                           read setting.json                                */
   /* -------------------------------------------------------------------------- */

   // Read settings from a JSON file and store values
   JSON settingJSON(input_directory / "setting.json");
   for (auto &[key, value] : settingJSON()) std::cout << key << ": " << value << std::endl;

   // Check for the existence of required keys in the JSON file
   if (!std::ranges::all_of(settingJSON.find({"end_time_step", "end_time", "output_directory", "max_dt", "input_files"}), [](auto v) { return v; }))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "there is missing key in setting.json");

   // Store values from the JSON file into variables
   const std::filesystem::path output_directory = settingJSON.at("output_directory")[0];
   const double max_dt = stod(settingJSON.at("max_dt"))[0];

   const int ALEPERIOD = stoi(settingJSON.at("ALEPERIOD"))[0];
   // Then, in your conditional checks, use this function to ensure case-insensitivity

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

   const int end_time_step = stoi(settingJSON.at("end_time_step"))[0];
   const double end_time = stod(settingJSON.at("end_time"))[0];

   // int ALE_period_in_step = 1;
   // if (settingJSON.find("ALE_period_in_step"))
   //    ALE_period_in_step = stod(settingJSON.at("ALE_period_in_step"))[0];

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
   auto setOutputInfo = [&NetOutputInfo](Network *net, auto output_name, const std::filesystem::path &output_directory) {
      std::cout << "setOutputInfo" << std::endl;
      NetOutputInfo[net].pvd_file_name = output_name;
      NetOutputInfo[net].vtu_file_name = output_name + "_";
      NetOutputInfo[net].PVD = new PVDWriter(output_directory / (output_name + ".pvd"));
   };

   std::vector<Network *> FluidObject, RigidBodyObject, SoftBodyObject, AbsorberObject;
   std::vector<JSON> MeasurementJSONs;
   auto setTypes = [&FluidObject, &RigidBodyObject, &SoftBodyObject, &AbsorberObject](Network *net, const JSON &injson) {
      std::cout << "setTypes" << std::endl;

      //$ set type
      auto type = injson.at("type")[0];
      if (type.contains("RigidBody")) {
         RigidBodyObject.emplace_back(net);
         net->isRigidBody = true;
         net->isSoftBody = net->isFluid = false;
      } else if (type.contains("SoftBody") || type.contains("FixedBody")) {
         SoftBodyObject.emplace_back(net);
         net->isSoftBody = true;
         net->isRigidBody = net->isFluid = false;
      } else if (type.contains("Fluid")) {
         FluidObject.emplace_back(net);
         net->isFluid = true;
         net->isRigidBody = net->isSoftBody = false;
      } else if (type.contains("Absorber") || type.contains("absorb") || type.contains("damping")) {
         AbsorberObject.emplace_back(net);
         net->isAbsorber = true;
         net->isRigidBody = net->isSoftBody = net->isFluid = false;
      }

      //! isFixedはdefaultでfalse．指定された分だけ順に置き換わる
      //! ただし，指定が１つだけなら，それを全てのに適用する．
      if (injson.find("isFixed")) {
         auto isFixed = stob(injson.at("isFixed"));
         if (isFixed.size() == 1)
            net->isFixed.fill(isFixed[0]);
         else
            for (auto i = 0; i < isFixed.size(); ++i)
               net->isFixed[i] = isFixed[i];
      }
      // for (auto i = 0; i < 10; ++i)
      //    AreaWeightedSmoothingPreserveShape(net->getPoints(), 0.1);
      //$ set velocity
      net->isFloatingBody = (injson.find("velocity") && injson.at("velocity")[0] == "floating");
      // velocityにfileが指定されている場合は，そのファイルを読み込み，
      // interpolationBsplineであるnet->intpMotionRigidBodyをsetする．
      injson.find("velocity", [&](auto STR_VEC) {
         if (STR_VEC[0].contains("file")) {
            if (STR_VEC.size() == 1) throw std::runtime_error("Failed to open the input file.");
            std::ifstream file(STR_VEC[1]);
            if (!file.is_open()) throw std::runtime_error("Failed to open the input file.");
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
      net->isFloatingBody = net->isFloatingBody || (injson.find("acceleration") && injson.at("acceleration")[0] == "floating");

      net->resetInitialX();
      net->setGeometricProperties();
   };

   /* generate network objects from input JSON files */

   for (auto input_file_name : settingJSON["input_files"]) {
      std::cout << Green << input_directory / input_file_name << colorReset << std::endl;
      JSON injson(input_directory / input_file_name);

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
      /* -------------------------------------------------------------------------- */
      /* -------------------- display contents of the JSON file ------------------- */
      for (auto &[key, value] : injson()) std::cout << Green << key << colorReset << ": " << value << std::endl;
      auto object_name = injson.at("name")[0];

      if (!injson.find("ignore") || !stob(injson["ignore"])[0]) {

         auto net = new Network(injson.at("objfile")[0], object_name);
         net->inputJSON = injson;
         net->applyTransformations(injson);
         setOutputInfo(net, injson.at("name")[0], output_directory);
         setTypes(net, injson);
         // output initial state
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

   /* ----------------- Create output directory and copy files ----------------- */

   std::filesystem::create_directories(output_directory);
   std::filesystem::copy_file(input_directory / "setting.json", output_directory / "setting.json", std::filesystem::copy_options::overwrite_existing);
   std::filesystem::copy_file("./main.cpp", output_directory / "main.cpp", std::filesystem::copy_options::overwrite_existing);
   //
   std::regex pattern("^BEM.*\\.hpp$");
   for (auto &entry : std::filesystem::directory_iterator("."))
      if (std::regex_match(entry.path().filename().string(), pattern))
         std::filesystem::copy_file(entry.path(), output_directory / entry.path().filename(), std::filesystem::copy_options::overwrite_existing);

   /* -------------------------------------------------------------------------- */

   // auto water = FluidObject[0];
   PVDWriter cornerPointsPVD(output_directory / "cornerPointsPVD.pvd");
   PVDWriter DirichletSurfacePVD(output_directory / "DirichletSurface.pvd");
   Print("setting done");

   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */
   /* -------------------------------------------------------------------------- */

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

   try {
      //  b* ------------------------------------------------------ */
      //  b*                         メインループ                      */
      //  b* ------------------------------------------------------ */
      TimeWatch watch;
      Buckets<networkFace *> FMM_BucketsFaces;
      Buckets<networkPoint *> FMM_BucketsPoints;
      for (time_step = 0; time_step < end_time_step; time_step++) {
         if (end_time < simulation_time)
            break;

         double dt = 1E+20;
         int RK_order = 4;

         // double spacing = 0.;
         // for (auto &water : FluidObject)
         //    spacing += Mean(extLength(water->getLines())) * 10;
         // spacing /= FluidObject.size();

         CoordinateBounds bounds(FluidObject[0]->bounds);
         for (auto &water : FluidObject)
            bounds += water->bounds;

         FMM_BucketsFaces.initialize(bounds, bounds.getScale() / 10.);
         FMM_BucketsPoints.initialize(bounds, bounds.getScale() / 10.);
         std::vector<Network *> AllObjects = Join(FluidObject, RigidBodyObject, SoftBodyObject);

         for (auto &water : FluidObject) {
            show_info(*water);
            // b# ------------------------------------------------------ */
            // b#                       刻み時間の決定                     */
            // b# ------------------------------------------------------ */
            auto dt_cfl = dt_CFL(*water, max_dt, .4);

            if (dt > dt_cfl)
               dt = dt_cfl;
            if (time_step == 0 && dt > 0.00001)
               dt = 0.00001;
         }
         if (dt < 1E-13)
            dt = 1E-13;
         Print("===========================================================================");
         Print("       dt :", Red, std::setprecision(10), dt, colorReset);
         Print("time_step :", Red, time_step, colorReset);
         Print("real time :", Red, simulation_time, colorReset);
         Print("---------------------------------------------------------------------------");

         Print("makeBucketFaces", Green);

#pragma omp parallel
         for (const auto &net : AllObjects)
#pragma omp single nowait
         {
            // auto spacing = 5 * Mean(extLength(net->getLines()));
            // net->makeBucketFaces(spacing);
            // auto spacing = 5 * Mean(extLength(net->getLines()));
            net->makeBucketFaces(net->getScale() / 7.);
         }

         for (auto &water : FluidObject) {
            //! 体積を保存するようにリメッシュする必要があるだろう．
            Print("setBoundaryTypes", Green);
            setBoundaryTypes(*water, Join(RigidBodyObject, SoftBodyObject));
            Print("setBoundaryTypes", Blue, "\nElapsed time: ", Red, watch(), colorReset, " s\n");
            double rad = M_PI / 180;

            {
               flipIf(*water, {10 * rad /*target n diff*/, 10 * rad /*change n diff*/}, {10 * rad, 10 * rad}, false);
               flipIf(*water, {10 * rad /*target n diff*/, 10 * rad /*change n diff*/}, {10 * rad, 10 * rad}, false);
               flipIf(*water, {10 * rad /*target n diff*/, 10 * rad /*change n diff*/}, {10 * rad, 10 * rad}, false);
            }
            if (time_step < 10 && time_step % 1 == 1) {
               flipIf(*water, {10 * rad, 10 * rad}, {10 * rad, 10 * rad}, true);
               flipIf(*water, {10 * rad, 10 * rad}, {10 * rad, 10 * rad}, true);
               flipIf(*water, {10 * rad, 10 * rad}, {10 * rad, 10 * rad}, true);
            }

            // b@ ------------------------------------------------------ */
            // b@           初期値問題を解く（時間微分方程式を数値積分する）        */
            // b@ ------------------------------------------------------ */
            for (const auto &p : water->getPoints()) {
               p->RK_phi.initialize(dt, simulation_time, std::get<0>(p->phiphin), RK_order);
               p->RK_X.initialize(dt, simulation_time, ToX(p), RK_order);
            }

            for (const auto &f : water->getFaces())
               FMM_BucketsFaces.add(f->getXtuple(), f);
            for (const auto &p : water->getPoints())
               FMM_BucketsPoints.add(ToX(p), p);
         }
         for (const auto &net : RigidBodyObject) {
            net->RK_COM.initialize(dt, simulation_time, net->COM, RK_order);
            net->RK_Q.initialize(dt, simulation_time, net->Q(), RK_order);
            net->RK_Velocity.initialize(dt, simulation_time, net->velocity, RK_order);
            //
            if (net->interp_accel.size() > 10)
               net->interp_accel.pop();
            net->interp_accel.push(simulation_time, net->acceleration);
         }
         for (const auto &net : SoftBodyObject) {
            // !いらないはずのもの
            for (const auto &p : net->getPoints())
               p->RK_X.initialize(dt, simulation_time, ToX(p), RK_order);
         }
         // b@ ----------------------------------------------------- */

         int RK_step = 0;
         BEM_BVP BVP(FluidObject);

         std::vector<double> convergence;
         do {

            auto RK_time = (*(FMM_BucketsPoints.all_stored_objects).begin())->RK_X.gett();  //%各ルンゲクッタの時刻を使う
            std::cout << "RK_step = " << ++RK_step << "/" << RK_order << ", RK_time = " << RK_time << ", simulation_time = " << simulation_time << std::endl;

#pragma omp parallel
            for (const auto &net : AllObjects)
#pragma omp single nowait
               net->makeBucketFaces(net->getScale() / 10.);
            std::cout << Green << "makeBucketFaces" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

            for (auto water : FluidObject) {
               setBoundaryTypes(*water, Join(RigidBodyObject, SoftBodyObject));
               water->setMinDepthFromCORNER();
            }
            std::cout << Green << "setBoundaryTypes" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

            setNeumannVelocity(Join(RigidBodyObject, SoftBodyObject));
            std::cout << Green << "setNeumannVelocity" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

            BVP.solve(FMM_BucketsPoints, FMM_BucketsFaces);
            std::cout << Green << "BVP.solve -> {Φ,Φn}が決まる" << Blue << "\nElapsed time: " << Red << watch() << " s\n";

            for (auto water : FluidObject)
               calculateCurrentVelocities(*water);
            std::cout << Green << "U_BEMを計算" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

            for (auto water : FluidObject) {
               // calculateCurrentUpdateVelocities(*water, 20, (RK_step == 4) ? 0.5 : 0.);
               // calculateCurrentUpdateVelocities(*water, 20, (RK_step == 4) ? 0.5 : 0.);
               bool do_ALE = (RK_step == 4 && (time_step % ALEPERIOD == 0));
               calculateCurrentUpdateVelocities(*water, 20, do_ALE ? 0.5 : 0.);
               if (!do_ALE)
                  std::cout << "ALE is not applied" << std::endl;
               // calculateCurrentUpdateVelocities(*water, 20);
            }
            std::cout << Green << "U_update_BEMを計算" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

            convergence = BVP.solveForPhiPhin_t(Join(RigidBodyObject, SoftBodyObject));
            std::cout << Green << "BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

            // b$ --------------------------------------------------- */

            /*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

            ### 浮体の重心位置・姿勢・速度の更新

            浮体の重心位置は，重心に関する運動方程式を解くことで求める．
            姿勢は，角運動量に関する運動方程式などを使って，各加速度を求める．姿勢はクオータニオンを使って表現する．

            */
#pragma omp parallel
            for (const auto &net : RigidBodyObject)
#pragma omp single nowait
            {  // \label{BEM:impose_velocity}
               //  重心位置と姿勢の時間発展
               if (net->inputJSON.find("velocity")) {
                  //! 時間まで静止させる
                  if ((net->inputJSON.at("velocity")[0].contains("floating") && net->inputJSON.at("velocity").size() >= 2 && std::stod(net->inputJSON.at("velocity")[1]) > simulation_time)) {
                     net->velocity.fill(0);
                     net->acceleration.fill(0);
                  }

                  //! 静止した物体は速度ゼロが与えられているので，下を実行しても動かない
                  if (!net->inputJSON.at("velocity")[0].contains("fixed")) {
                     std::cout << "updating " << net->getName() << "'s (RigidBodyObject) velocity" << std::endl;
                     //
                     net->RK_COM.push(net->velocityTranslational());
                     net->COM = net->RK_COM.getX();
                     //
                     net->RK_Q.push(net->Q.AngularVelocityTodQdt(net->velocityRotational()));
                     net->Q = net->RK_Q.getX();
                     //
                  }
                  std::cout << "name = " << net->getName() << std::endl;
                  std::cout << "net->velocityTranslational() = " << net->velocityTranslational() << std::endl;

                  for (auto &mooring : net->mooringLines) {
                     //! mooring->lastPoint は，浮体ともに動く．
                     auto Xcurrent = mooring->lastPoint->X;
                     mooring->lastPoint->X_last = Xcurrent;
                     Tddd V = (nextPositionOnBody(net, mooring->lastPoint) - Xcurrent) / (net->RK_COM.getTimeAtNextStep() - simulation_time);
                     mooring->simulate(simulation_time, net->RK_COM.getTimeAtNextStep() - simulation_time, [&](networkPoint *p) {
                        if (p == mooring->firstPoint) {
                           p->acceleration.fill(0);
                           p->velocity.fill(0);
                        } else if (p == mooring->lastPoint) {
                           p->acceleration.fill(0);
                           p->velocity[0] = V[0];
                           p->velocity[1] = V[1];
                           p->velocity[2] = V[2];
                        }
                     });

                     if (net->RK_Q.finished)
                        mooring->applyMooringSimulationResult();

                     for (const auto &p : mooring->getPoints())
                        mooring->lastPoint->setX(nextPositionOnBody(net, mooring->lastPoint));
                  }
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

               for (const auto &p : net->getPoints())
                  p->setXSingle(net->Q.R(p->initialX - net->ICOM) + net->COM);
               net->setGeometricProperties();
            }

            // b$ --------------------------------------------------- */

            for (const auto &net : SoftBodyObject) {
               std::cout << "updating " << net->getName() << "'s (SoftBodyObject) position" << std::endl;
               for (const auto &p : net->getPoints()) {
                  auto V = p->velocityTranslational();
                  if (net->isFixed[0])
                     V[0] = 0;
                  if (net->isFixed[1])
                     V[1] = 0;
                  if (net->isFixed[2])
                     V[2] = 0;
                  p->RK_X.push(V);  //@ 位置xの時間発展
                  // p->setXSingle(p->RK_X.getX());
               }

               net->setGeometricProperties();
            }

            // b$ -------------------------------------------------------------------------- */

            /*DOC_EXTRACT 0_4_1_UPDATE_POSITION

            ### 流体の$`\phi`$時間発展，$`\phi_n`$の時間発展はない

            ### 波の吸収（ダンピング領域）

            ```math
            \begin{aligned}
            \gamma &= 1 - 2 \frac{\text{horizontal distance from the center of the absorber}}{\text{width of the absorber}} \\
            \phi_{\rm ref} &= \frac{\sum \phi \cdot \text{area}}{\sum \text{area}}
            \end{aligned}
            ```

            */

            auto isAbsorbed = [&](auto p) -> Network * {
               for (const auto &net : AbsorberObject)
                  if (net->isInside(p->X))
                     return net;
               return nullptr;
            };

            // initialize reference_phi and set absorbedBy
            std::unordered_set<Network *> absorbers;
            for (const auto &p : FMM_BucketsPoints.all_stored_objects) {
               p->absorbedBy = isAbsorbed(p);
               if (p->absorbedBy != nullptr)
                  absorbers.insert(p->absorbedBy);
            }

            std::unordered_map<Network *, std::array<double, 2>> reference_phi;
            for (const auto &absorber : absorbers)
               reference_phi[absorber] = {0, 0};

            // calculate reference_phi and weight
            for (const auto &f : FMM_BucketsFaces.all_stored_objects) {
               auto [p0, p1, p2] = f->getPoints();
               auto absorber = p0->absorbedBy;
               if (absorber != nullptr && absorber == p1->absorbedBy && absorber == p2->absorbedBy) {
                  auto mean_phi = (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3;
                  reference_phi[absorber] += Tdd{mean_phi * f->area, f->area};
               }
            }

            for (const auto &p : FMM_BucketsPoints.all_stored_objects) {
               p->U_absorbed.fill(0.);
               //! damping coeeficient is determined by the horizontal distance from the center of the absorber

               double gamma = 0, ref_phi = 0;

               if (p->absorbedBy != nullptr) {
                  auto [x, y, z] = p->absorbedBy->bounds;
                  auto absorber_center = Tddd{Mean(x), Mean(y), Mean(z)};
                  auto horizontal_dist_from_center = Norm(Chop(p->X - absorber_center, Tddd{0, 0, 1}));
                  double x_scale = std::abs(x[0] - x[1]);
                  gamma = 1. - std::clamp(2 * Norm(horizontal_dist_from_center) / x_scale, 0.1, 1.);
                  ref_phi = reference_phi.at(p->absorbedBy)[0] / reference_phi.at(p->absorbedBy)[1];
               }
               //
               if (!p->Neumann /*Neumannを変更しても，あとでBIEによって上書きされるので，からわない．*/) {
                  p->RK_phi.push(p->DphiDt_damped({gamma, ref_phi}, p->U_update_BEM, 0.));
                  std::get<0>(p->phiphin) = p->phi_Dirichlet = p->RK_phi.getX();  // 角点の法線方向はわからないので，ノイマンの境界条件phinを与えることができない．
               }
               //@ 位置xの時間発展
               p->RK_X.push(p->U_update_BEM - gamma * p->U_update_BEM);
               p->setXSingle(p->RK_X.getX());
               p->phi_tmp = 0;
            }

            // b$ -------------------------------------------------------------------------- */

            for (auto water : FluidObject)
               std::cout << Green << "name:" << water->getName() << ": setBounds" << colorReset << std::endl;

            for (auto net : AllObjects)
               net->setGeometricProperties();

            for (auto water : FluidObject) {
               auto name = water->getName();
               std::ofstream ofs(output_directory / (name + std::to_string(RK_step) + ".obj"));
               createOBJ(ofs, *water);
               ofs.close();
            }

            std::cout << Blue << "Elapsed time: " << Red << watch() << colorReset << " s\n";
         } while (!((*(FMM_BucketsPoints.all_stored_objects).begin())->RK_X.finished));

         /* -------------------------------------------------------------------------- */
         // for (const auto &net : SoftBodyObject)
         //    AreaWeightedSmoothingPreserveShape(net->getPoints(), 1);

         /* ------------------------------------------------------ */

         //! 速度ポテンシャルの平均を0にする
         {
            double mean_phi = 0.;
            for (const auto &p : FMM_BucketsPoints.all_stored_objects)
               mean_phi += std::get<0>(p->phiphin);
            mean_phi /= FMM_BucketsPoints.all_stored_objects.size();
            for (const auto &p : FMM_BucketsPoints.all_stored_objects) {
               p->phi_Dirichlet -= mean_phi;
               // p->phi_Neumann -= mean_phi;
               std::get<0>(p->phiphin) -= mean_phi;
            }
         }

         //! クォータニオンの正規化
         for (const auto &net : RigidBodyObject)
            net->Q.normalize();

         /* ------------------------------------------------------ */

         std::cout << Green << "simulation_timeを取得" << colorReset << std::endl;
         simulation_time = (*(FMM_BucketsPoints.all_stored_objects).begin())->RK_X.gett();

         for (const auto &p : FMM_BucketsPoints.all_stored_objects) {
            p->U_BEM_last = p->U_BEM;
            p->U_tangential_BEM_last = p->U_tangential_BEM;
         }

         /* ------------------------------------------------------ */
         for (auto water : FluidObject) {
            auto name = water->getName();
            std::ofstream ofs(output_directory / (name + ".obj"));
            createOBJ(ofs, *water);
            ofs.close();
         }
         // b# -------------------------------------------------------------------------- */
         // b#                              output JSON files                             */
         // b# -------------------------------------------------------------------------- */

         /*DOC_EXTRACT 1_0_0_JSON_OUTPUT

         ### JSONファイルの出力

         JSONファイルには，計算結果を出力する．

         流体の場合

         | 項目 | 詳細|
         |---:|:---|
         | `simulation_time` | シミュレーション上の時間 |
         | `cpu_time` | CPU時間(CPUがプログラムを実行していた時間の合計) |
         | `wall_clock_time` | 実時間 |
         | `***_volume` | 流体の体積 |
         | `***_EK` | 流体の運動エネルギー |
         | `***_EP` | 流体の位置エネルギー |
         | `***_E` | 流体の全エネルギー |

         剛体などで，浮体か`output`に`json`が指定されている場合

         | 項目 | 詳細|
         |---:|:---|
         | `simulation_time` | シミュレーション上の時間 |
         | `cpu_time` | CPU時間(CPUがプログラムを実行していた時間の合計) |
         | `wall_clock_time` | 実時間 |
         | `***_pitch` | 浮体のピッチ角 |
         | `***_yaw` | 浮体のヨー角 |
         | `***_roll` | 浮体のロール角 |
         | `***_force` | 浮体に働く力 |
         | `***_torque` | 浮体に働くトルク |
         | `***_accel` | 浮体の加速度 |
         | `***_velocity` | 浮体の速度 |
         | `***_COM` | 浮体の重心位置 |
         | `***_area` | 浮体の面積 |
         | `***_EK` | 浮体の運動エネルギー |
         | `***_EP` | 浮体の位置エネルギー |

         */

         std::cout << "output JSON files" << std::endl;

         double worst_iteration = 0;
         double worst_value = 0;
         double worst_grad = 0;
         std::array<double, 3> convergence_info = {0, 0, 0};
         for (auto water : FluidObject) {
            for (const auto &p : water->getPoints()) {
               gradPhi(p, convergence_info);
               if (convergence_info[0] > worst_iteration)
                  worst_iteration = convergence_info[0];
               if (convergence_info[1] > worst_grad)
                  worst_grad = convergence_info[1];
               if (convergence_info[2] > worst_value)
                  worst_value = convergence_info[2];
            }
         }

         jsonout.push("gradPhi_worst_iteration", worst_iteration);
         jsonout.push("gradPhi_worst_value", worst_value);
         jsonout.push("gradPhi_worst_grad", worst_grad);

         jsonout.push("simulation_time", simulation_time);
         // Push CPU time
         double cpu_time_in_seconds = static_cast<double>(std::clock() - cpu_clock_start) / CLOCKS_PER_SEC;
         jsonout.push("cpu_time", cpu_time_in_seconds);

         jsonout.push("eq_of_motion", convergence);

         // Calculate and push wall-clock time
         auto duration = std::chrono::high_resolution_clock::now() - wall_clock_start;
         double wall_clock_time_in_seconds = std::chrono::duration<double>(duration).count();
         jsonout.push("wall_clock_time", wall_clock_time_in_seconds);

         /* 流体 */
         for (const auto &water : FluidObject) {
            jsonout.push(water->getName() + "_volume", water->getVolume());
            jsonout.push(water->getName() + "_EK", KinematicEnergy(water->getFaces()));
            jsonout.push(water->getName() + "_EP", PotentialEnergy(water->getFaces()));
            jsonout.push(water->getName() + "_E", TotalEnergy(water->getFaces()));
            // fluid mesh size faces and point
            jsonout.push(water->getName() + "_face_size", (double)water->getFaces().size());
            jsonout.push(water->getName() + "_point_size", (double)water->getPoints().size());
         }
         /* 物体 */
         for (const auto &net : Join(RigidBodyObject, SoftBodyObject)) {
            if ((net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating") ||
                (net->inputJSON.find("output") && std::ranges::any_of(net->inputJSON.at("output"), [](const auto &s) { return s == "json"; }))) {
               auto tmp = calculateFluidInteraction(FMM_BucketsFaces.all_stored_objects, net);
               jsonout.push(net->getName() + "_pitch", net->quaternion.pitch());
               jsonout.push(net->getName() + "_yaw", net->quaternion.yaw());
               jsonout.push(net->getName() + "_roll", net->quaternion.roll());
               auto [force, torque] = tmp.surfaceIntegralOfPressure();
               jsonout.push(net->getName() + "_force", force);
               jsonout.push(net->getName() + "_torque", torque);
               auto [drag_force, drag_torque] = tmp.surfaceIntegralOfVerySimplifiedDrag();
               jsonout.push(net->getName() + "_drag_force", drag_force);
               jsonout.push(net->getName() + "_drag_torque", drag_torque);
               jsonout.push(net->getName() + "_accel", net->acceleration);
               jsonout.push(net->getName() + "_velocity", net->velocity);
               jsonout.push(net->getName() + "_COM", net->COM);
               jsonout.push(net->getName() + "_area", tmp.area);
               jsonout.push(net->getName() + "_EK", Dot(net->velocity, net->velocity) * net->mass / 2.);
               jsonout.push(net->getName() + "_EP", net->getMass3D() * _GRAVITY3_ * (net->COM[2] - net->ICOM[2]));
            }
         }

         /* 計測 */
         for (const auto &json : MeasurementJSONs) {
            if (json.find("position")) {
               const auto position = json.at("position");
               if (position.size() != 6)
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "position must have 6 elements");
               const double x0 = std::stod(position[0]);
               const double y0 = std::stod(position[1]);
               const double z0 = std::stod(position[2]);
               const double x1 = std::stod(position[3]);
               const double y1 = std::stod(position[4]);
               const double z1 = std::stod(position[5]);

               Tddd nearest_intersection = {1E+20, 1E+20, 1E+20};
               for (const auto &water : FluidObject) {
                  water->BucketFaces.apply(water->BucketFaces.line2indices({x0, y0, z0}, {x1, y1, z1}), [&](networkFace *f) {
                     Tddd X0 = {x0, y0, z0}, X1 = {x1, y1, z1};
                     auto [isintersect, X, t0t1T] = IntersectQ_(T2Tddd{X0, X1}, (T3Tddd)(*f));
                     if (isintersect && (Norm(X - X0) < Norm(nearest_intersection - X0)))
                        nearest_intersection = X;
                  });
               }
               jsonout.push(json.at("name")[0] + "_intersection", nearest_intersection);
            }
         }

         std::ofstream os(output_directory / "result.json");
         jsonout.output(os);
         os.close();

         // b# -------------------------------------------------------------------------- */
         // b#                            output Paraview files                           */
         // b# -------------------------------------------------------------------------- */

         // check dirichlet surface to check interpolated surface
         {
            auto dirichletsurface = new Network();
            std::filesystem::path filename_dir_surface = "DirichletSurface" + std::to_string(time_step) + ".vtu";
            std::vector<std::array<std::array<double, 2>, 3>> t0t1_triangles = SymmetricSubdivisionOfTriangle_00_10_01(4);
            for (const auto &water : FluidObject) {
               for (const auto &f : water->getFaces()) {
                  if (f->Dirichlet) {
                     DodecaPoints dodecapoint(f, f->getPoints()[0], [](const networkLine *line) -> bool { return !line->CORNER; });
                     for (auto t0t1 : t0t1_triangles) {
                        new networkFace(dirichletsurface, {new networkPoint(dirichletsurface, dodecapoint.X(t0t1[0])),
                                                           new networkPoint(dirichletsurface, dodecapoint.X(t0t1[1])),
                                                           new networkPoint(dirichletsurface, dodecapoint.X(t0t1[2]))});
                     }
                  }
               }
            }

            std::unordered_map<networkPoint *, Tddd> data;
            for (auto p : dirichletsurface->getPoints())
               data[p] = ToX(p);
            mk_vtu(output_directory / filename_dir_surface, dirichletsurface->getFaces(), {{"position", data}});
            DirichletSurfacePVD.push(filename_dir_surface, simulation_time);
            DirichletSurfacePVD.output();
            delete dirichletsurface;
         }

         // 流体
         for (const auto &net : FluidObject) {
            std::filesystem::path filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory / filename, net->getFaces(), dataForOutput(net, dt));
            NetOutputInfo[net].PVD->push(filename, simulation_time);
            NetOutputInfo[net].PVD->output();
         }

         for (const auto &net : Join(RigidBodyObject, SoftBodyObject)) {
            VV_VarForOutput data;
            if (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating") {
               auto tmp = calculateFluidInteraction(FMM_BucketsFaces.all_stored_objects, net);
               uomap_P_Tddd P_COM, P_COM_p, P_accel, P_velocity, P_rotational_velocity, P_rotational_accel;
               uomap_P_d P_Pitch, P_Yaw, P_Roll, P_radius;
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
                   {"velocity", P_velocity},
                   {"acceleration", P_accel},
                   {"rotational velocity", P_rotational_velocity},
                   {"rotational acceleration", P_rotational_accel}};
               std::filesystem::path name = "actingFacesOn" + NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
               mk_vtu(output_directory / name, tmp.actingFaces, data);
            }
            std::filesystem::path filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory / filename, net->getFaces(), data);
            NetOutputInfo[net].PVD->push(filename, simulation_time);
            NetOutputInfo[net].PVD->output();

            //$ --------------------------------------------------- */
            //$                  MOORING LINE                       */
            //$ --------------------------------------------------- */
            /*DOC_EXTRACT 0_6_MOORING_LINE

            ### 係留索の出力


            */
            //! rigid body includes mooring lines output them
            for (const auto &mooring : net->mooringLines) {

               std::unordered_map<networkPoint *, Tddd> P_total_force, P_dragforce, P_tension, P_gforce, P_velocity, P_accel;
               P_total_force.reserve(mooring->getPoints().size());
               P_dragforce.reserve(mooring->getPoints().size());
               P_gforce.reserve(mooring->getPoints().size());
               P_tension.reserve(mooring->getPoints().size());
               P_velocity.reserve(mooring->getPoints().size());
               P_accel.reserve(mooring->getPoints().size());
               for (auto p : mooring->getPoints()) {
                  P_total_force[p] = p->getForce();
                  P_dragforce[p] = p->getDragForce();
                  P_tension[p] = p->getTension();
                  P_gforce[p] = p->getGravitationalForce();
                  P_velocity[p] = p->velocityTranslational();
                  P_accel[p] = p->accelTranslational();
               }
               std::vector<std::tuple<std::string, decltype(P_dragforce)>> data = {std::make_tuple("drag force", P_dragforce),
                                                                                   std::make_tuple("tension", P_tension),
                                                                                   std::make_tuple("gravitational force", P_gforce),
                                                                                   std::make_tuple("total force", P_total_force),
                                                                                   std::make_tuple("velocity", P_velocity),
                                                                                   std::make_tuple("acceleration", P_accel)};
               std::filesystem::path filename = NetOutputInfo[mooring].vtu_file_name + std::to_string(time_step) + ".vtp";
               std::ofstream ofs(output_directory / filename);
               vtkPolygonWrite(ofs, mooring->getLines(), data);
               NetOutputInfo[mooring].PVD->push(filename, simulation_time);
               NetOutputInfo[mooring].PVD->output();
            }
         }

         // b* ------------------------------------------------------ */
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
   return 0;
};

/*DOC_EXTRACT 2_0_0_HOW_TO_RUN

# 実行方法

## ファイルのダウンロード

上書きされるので注意．ダウンロードしたら，`build_bem`ディレクトリに移動．

```sh
git clone https://github.com/tomoakihirakawa/cpp.git
cd ./cpp/builds/build_bem
```

## 入力ファイルの生成．

```sh
python3 input_generator.py
```

例えば，`./input_files/Hadzic2005`が生成される．

## プログラムのコンパイルと実行

`clean`でCMake関連のファイルを削除して（ゴミがあるかもしれないので），
`cmake`で`Makefile`を生成して，`make`でコンパイルする．

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

実行

```sh
./main ./input_files/Hadzic2005
```

*/

/*DOC_EXTRACT 3_0_EXAMPLES

# Examples

**[See the Examples here!](EXAMPLES.md)**

*/
