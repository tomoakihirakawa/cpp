#include "pch.hpp"

/*DOC_EXTRACT 0_0_BEM

# BEM-MEL

<img src="./sample_Goring1979.gif">

*/

// #define _debugging_

bool _LINEAR_ELEMENT_ = false;
bool _PSEUDO_QUADRATIC_ELEMENT_ = false;
bool _ALE_ON_LINEAR_ELEMENT_ = false;
bool _ALE_ON_PSEUDO_QUADRATIC_ELEMENT_ = false;
int time_step;
double simulation_time = 0;

#define BEM

#define simulation
// #include <sys/utsname.h>
// #include <unistd.h>
// #include <chrono>
// #include <ctime>
// #include <filesystem>
// #include <regex>
// #include "tetgen1.6.0/tetgen.h"
//
#include "Network.hpp"
// #include "integrationOfODE.hpp"
// #include "kernelFunctions.hpp"
// #include "minMaxOfFunctions.hpp"
// #include "rootFinding.hpp"
// #include "vtkWriter.hpp"

JSONoutput jsonout;

// pvd cpg_pvd("./vtu/bem.pvd");

#include "BEM.hpp"
#include "BEM_inputfile_reader.hpp"
// #include "svd.hpp"
// 追加
#include "OutputCommon.hpp"
#include "OutputJSON.hpp"
#include "OutputParaview.hpp"

int main(int argc, char **argv) {
  std::clock_t cpu_clock_start = std::clock();
  auto wall_clock_start = std::chrono::high_resolution_clock::now();

  /*DOC_EXTRACT 0_1_BEM

  ## 入力ファイルの読み込み

  1. 境界条件の設定
  2. 境界値問題（BIE）を解き，$\phi$と$\phi_n$を求める
  3. 三角形の線形補間を使って節点の流速を計算する

  */

  /* --------------------------------------------------------------------------
   */
  /*                           Set up logging to file */
  /* --------------------------------------------------------------------------
   */
  if (!initializeLogFile("log.txt", argc, argv))
    return 1;
  if (argc <= 1)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "argv <= 1. write input json file directory!\nex.\n$ ./main ./input");
  SimulationSettings setting(argv[1]);
  std::filesystem::path input_directory = setting.input_directory;
  JSON settingJSON = setting.settingJSON;
  const double max_dt = setting.max_dt;
  const int ALEPERIOD = setting.ALEPERIOD;
  const int end_time_step = setting.end_time_step;
  const double end_time = setting.end_time;
  const double stop_remesh_time = setting.stop_remesh_time;
  const double force_remesh_time = setting.force_remesh_time;
  const bool tetrahedralize = setting.tetrahedralize;
  const bool surface_flip = setting.surface_flip;
  const int grid_refinement = setting.grid_refinement;
  const std::filesystem::path output_directory = setting.output_directory;
  _LINEAR_ELEMENT_ = setting._LINEAR_ELEMENT_;
  _PSEUDO_QUADRATIC_ELEMENT_ = setting._PSEUDO_QUADRATIC_ELEMENT_;
  _ALE_ON_LINEAR_ELEMENT_ = setting._ALE_ON_LINEAR_ELEMENT_;
  _ALE_ON_PSEUDO_QUADRATIC_ELEMENT_ = setting._ALE_ON_PSEUDO_QUADRATIC_ELEMENT_;
  std::map<std::string, outputInfo> NetOutputInfo = setting.NetOutputInfo;
  std::vector<Network *> FluidObject = setting.FluidObject;
  std::vector<Network *> RigidBodyObject = setting.RigidBodyObject;
  std::vector<Network *> SoftBodyObject = setting.SoftBodyObject;
  std::vector<Network *> AbsorberObject = setting.AbsorberObject;
  std::vector<JSON> MeasurementJSONs = setting.MeasurementJSONs;
  /* --------------------------------------------------------------------------
   */

  if (tetrahedralize)
    for (auto &network : FluidObject)
      network->tetrahedralize();

  /* --------------------------------------------------------------------------
   */
  /* ----------------- Create output directory and copy files -----------------
   */

  std::filesystem::create_directories(output_directory);
  std::filesystem::copy_file(input_directory / "setting.json", output_directory / "setting.json", std::filesystem::copy_options::overwrite_existing);
  std::filesystem::copy_file("./main.cpp", output_directory / "main.cpp", std::filesystem::copy_options::overwrite_existing);
  //
  std::regex pattern("^BEM.*\\.hpp$");
  for (auto &entry : std::filesystem::directory_iterator("."))
    if (std::regex_match(entry.path().filename().string(), pattern))
      std::filesystem::copy_file(entry.path(), output_directory / entry.path().filename(), std::filesystem::copy_options::overwrite_existing);

  /* --------------------------------------------------------------------------*/

  // auto water = FluidObject[0];
  PVDWriter cornerPointsPVD(output_directory / "cornerPointsPVD.pvd");
  PVDWriter DirichletSurfacePVD(output_directory / "DirichletSurface.pvd");
  Print("setting done");

  /* --------------------------------------------------------------------------*/
  /* --------------------------------------------------------------------------*/
  /* --------------------------------------------------------------------------*/

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
  2. 境界値問題（BIE）を解き，$\phi$と$\phi_n$を求める
  3. 三角形の線形補間を使って節点の流速を計算する
  4. 次時刻の$\Omega(t+\Delta t)$がわかるので，修正流速を計算する
  5.
  浮体の加速度を計算する．境界値問題（BIE）を解き，$\phi_t$と$\phi_{nt}$を求め，浮体面上の圧力$p$を計算する必要がある
  6. 全境界面の節点の位置を更新．ディリクレ境界では$\phi$を次時刻の値へ更新

  */

  try {
    //  * ------------------------------------------------------ */
    //  *                         メインループ                   */
    //  * ------------------------------------------------------ */
    Buckets<networkFace *> Buckets_Surfaces;
    Buckets<networkPoint *> Buckets_AllPoints;

    double time_setIGIGn, time_solve, unknownsizse;

    for (time_step = 0; time_step < end_time_step; time_step++) {
      if (end_time < simulation_time)
        break;

      double dt = 1E+20;
      int RK_order = 4;

      double spacing = 0.;
      for (auto &water : FluidObject)
        spacing += Mean(extLength(water->getLines())) * 10;
      spacing /= FluidObject.size();
      double rad = M_PI / 180;

      for (auto &water : FluidObject) {
        int count = 0;

        /* --------------------------------- 四面体の削除 --------------------------------- */

        double mean_tetra_vol = 0;

        {
          auto tetras = water->getTetras();
          for (const auto &t : tetras)
            mean_tetra_vol += t->getVolume();
          mean_tetra_vol /= tetras.size();
        }

        // delete inner lines and faces
        {
          std::vector<networkFace *> inner_faces;
          for (const auto &f : water->getFaces())
            if (!f->SurfaceQ())
              inner_faces.push_back(f);
          std::vector<networkLine *> inner_lines;
          for (const auto &l : water->getLines())
            if (std::ranges::none_of(l->getFaces(), [&](const auto &f) { return f->SurfaceQ(); }))
              inner_lines.push_back(l);
          std::vector<networkPoint *> inner_points;
          for (const auto &p : water->getPoints()) {
            bool isInner = true;
            if (std::ranges::none_of(p->getFaces(), [&](const auto &f) { return f->SurfaceQ(); }))
              inner_points.push_back(p);
          }

          for (const auto &f : inner_faces)
            delete f;

          for (const auto &l : inner_lines)
            delete l;

          for (const auto &p : inner_points)
            delete p;
        }

        auto tetras = water->getTetras();
        if (!tetras.empty()) {
          std::cout << "tetrahedra are not empty after deleting inner lines and faces!" << std::endl;
          for (const auto &t : tetras)
            delete t;
        }

        water->setGeometricProperties();
        water->checkConnectivity();

        for (const auto &l : water->getLines())
          if (!l->checkTopology())
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error after deleting inner elements");

        for (auto i = 0; i < 1; i++) {
          /* --------------------------------- divide --------------------------------- */

          for (auto iter = 0; iter < 10; iter++) {
            bool divided_any = false;
            auto SurfaceFaces = water->getSurfaces();
            for (const auto &f : SurfaceFaces) {
              auto lines = f->getLines();
              for (const auto &l : lines) {
                auto len = l->length();

                auto mean_line_len = 0.;
                {
                  count = 0;
                  for (const auto &f : SurfaceFaces)
                    for (const auto &l : f->getLines()) {
                      auto len = l->length();
                      mean_line_len += len;
                      count++;
                    }
                  mean_line_len /= count;
                }

                if (len > 2.5 * mean_line_len) {
                  l->divide();
                  divided_any = true;
                  std::cout << Red << "time_step " << time_step << ": line divided due to large length. length = " << len << ", mean length = " << mean_line_len << colorReset << std::endl;
                  break;
                }
              }
              if (divided_any)
                break;
            }
            if (!divided_any)
              break;
          }

          water->setGeometricProperties();
          water->checkConnectivity();
          for (const auto &l : water->getLines())
            if (!l->checkTopology())
              throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error detected in line after division");

          /* ---------------------------------- merge --------------------------------- */

          bool any_found_small_line = false;
          for (auto iter = 0; iter < 10; iter++) {
            bool found_small_line = false;
            auto Surfaces = water->getSurfaces();
            std::unordered_set<networkLine *> SurfaceLines;
            for (const auto &f : Surfaces)
              if (f->SurfaceQ()) {
                auto lines = f->getLines();
                for (const auto &l : lines) {
                  auto surfaces = l->getSurfaces();
                  if (surfaces.size() != 2) {
                    std::cout << "surfaces.size() != 2 but " << surfaces.size() << std::endl;
                    std::cout << "f: " << f << ", l: " << l << std::endl;
                    std::cout << "surfaces: " << surfaces << std::endl;
                    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "surfaces.size() != 2");
                  }
                  SurfaceLines.emplace(l);
                }
              }

            for (const auto &l : SurfaceLines) {
              auto surfaces = l->getSurfaces();
              if (surfaces.size() != 2) {
                std::cout << "surfaces: " << surfaces << std::endl;
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "surfaces.size() != 2");
              }
              auto f0 = surfaces[0];
              auto f1 = surfaces[1];
              auto len = l->length();
              auto mean_line_len = 0.;
              {
                count = 0;
                for (const auto &f : surfaces)
                  for (const auto &l : f->getLines()) {
                    auto len = l->length();
                    mean_line_len += len;
                    count++;
                  }
                mean_line_len /= count;
              }

              // 線の長さが平均の1/5以下ならマージ
              if (len < mean_line_len / 5.) {
                l->merge();
                found_small_line = true;
                any_found_small_line = true;
                std::cout << "time_step " << time_step << ": line merged due to small length. length = " << len << ", mean length = " << mean_line_len << std::endl;
                break;
              }

              auto [p0, p1, p2] = f0->getPoints(l);
              auto ci0 = CircumradiusToInradius(p0->X, p1->X, p2->X);
              auto [q0, q1, q2] = f1->getPoints(l);
              auto ci1 = CircumradiusToInradius(q0->X, q1->X, q2->X);

              // 接する面の法線が反対方向でかつ近ければマージ
              if (Dot(f0->normal, -f1->normal) > std::cos(20. * rad)) {
                l->flip(true); // できればでいい
                l->merge();
                found_small_line = true;
                any_found_small_line = true;
                std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
                break;
              } else if (Norm(p2->X - q2->X) < std::sqrt(TriangleArea(p0->X, p1->X, p2->X) + TriangleArea(q0->X, q1->X, q2->X)) / 5.) {
                l->flip(true); // できればでいい
                l->merge();
                found_small_line = true;
                any_found_small_line = true;
                std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
                break;
              } else if (ci0 > 50. || ci1 > 50.) {
                l->merge();
                found_small_line = true;
                any_found_small_line = true;
                std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
                break;
              }
            }
            if (!found_small_line)
              break;
          }

          water->setGeometricProperties();
          water->checkConnectivity();

          for (const auto &l : water->getLines())
            if (!l->checkTopology())
              throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error detected in line after merge");

          /* --------------------------- エッジフリップによる表面メッシュの改善-------------------------- */

          if (surface_flip) {
            flipIf(*water, {20 * rad /*target n diff*/, 20 * rad /*change n diff*/}, {20 * rad, 20 * rad}, false);
            flipIf(*water, {20 * rad /*target n diff*/, 20 * rad /*change n diff*/}, {20 * rad, 20 * rad}, false);
            flipIf(*water, {20 * rad /*target n diff*/, 20 * rad /*change n diff*/}, {20 * rad, 20 * rad}, false);
            flipIf(*water, {10 * rad /*target n diff*/, 10 * rad /*change n diff*/}, {10 * rad, 10 * rad}, false);
            flipIf(*water, {10 * rad /*target n diff*/, 10 * rad /*change n diff*/}, {10 * rad, 10 * rad}, false);
            flipIf(*water, {10 * rad /*target n diff*/, 10 * rad /*change n diff*/}, {10 * rad, 10 * rad}, false);
          }

          water->setGeometricProperties();
          water->checkConnectivity();

          for (const auto &l : water->getLines())
            if (!l->checkTopology())
              throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error after flip");
        }

        water->tetrahedralize();
        water->setGeometricProperties();

        if (tetrahedralize)
          water->improveTetrahedraDelaunay();
      }
      //
      CoordinateBounds bounds_org(FluidObject[0]->bounds);
      for (auto &water : FluidObject)
        bounds_org += water->bounds;

      CoordinateBounds bounds = bounds_org.scaledBounds(1.5);

      Buckets_Surfaces.initialize(bounds, bounds.getScale() / 10.);
      Buckets_AllPoints.initialize(bounds, bounds.getScale() / 10.);
      std::vector<Network *> AllObjects = Join(FluidObject, RigidBodyObject, SoftBodyObject);

      // # ------------------------------------------------------ */
      // #                       刻み時間の決定                   */
      // # ------------------------------------------------------ */
      for (auto &water : FluidObject) {
        show_info(*water);

        auto dt_cfl = dt_CFL(*water, max_dt, .3);
        if (dt > dt_cfl)
          dt = dt_cfl;
        if (time_step <= 2)
          dt = dt / 10.;
      }
      if (dt < 1E-13)
        dt = 1E-13;
      Print("===========================================================================");
      Print("       dt :", Red, std::setprecision(10), dt, colorReset);
      Print("time_step :", Red, time_step, colorReset);
      Print("real time :", Red, simulation_time, colorReset);
      Print("---------------------------------------------------------------------------");

      for (auto &water : FluidObject) {
        // b@ ------------------------------------------------------ */
        // b@           初期値問題を解く（時間微分方程式を数値積分する） */ b@
        // ------------------------------------------------------ */
        for (const auto &p : water->getPoints()) {
          p->RK_phi.initialize(dt, simulation_time, std::get<0>(p->phiphin), RK_order);
          p->RK_X.initialize(dt, simulation_time, ToX(p), RK_order);
        }

        for (const auto &f : water->getSurfaces())
          Buckets_Surfaces.add(f->getXtuple(), f);
        for (const auto &p : water->getPoints())
          Buckets_AllPoints.add(ToX(p), p);
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

      //   正面を向き合い，近い位置にある水面を探査
      //   std::array<networkFace *, 2> finding_pairs;
      //   for (size_t i = 0; i < FluidObject.size(); i++) {
      //     for (size_t j = i + 1; j < FluidObject.size(); j++) {
      //       auto [face1, face2, distance] =
      //       findClosestFacePair(FluidObject[i], FluidObject[j]); if (distance
      //       < 10. * spacing) {
      //         finding_pairs = {face1, face2};
      //         std::cout << "Found facing surfaces between " <<
      //         FluidObject[i]->getName() << " and " <<
      //         FluidObject[j]->getName() << ", distance = " << distance <<
      //         std::endl;
      //       }
      //     }
      //   }

      // b@ ----------------------------------------------------- */

      int RK_step = 0;
      BEM_BVP BVP(FluidObject);
      std::vector<double> convergence;
      TimeWatch watch;

      std::vector<double> ElapsedTimeBIEDiscretization;
      std::vector<double> ElapsedTimeSolve;
      double ElapsedTimeALE;
      double ElapsedTimeTotal;
      do {
        auto RK_time = (*(Buckets_AllPoints.data1D).begin())->RK_X.gett(); //%各ルンゲクッタの時刻を使う
        std::cout << "RK_step = " << ++RK_step << "/" << RK_order << ", RK_time = " << RK_time << ", simulation_time = " << simulation_time << std::endl;

        /* --------------------------------------------------- */

        if (RK_step == 1) {
#pragma omp parallel for
          for (const auto &net : AllObjects)
            net->makeBuckets(net->getScale() / 10.);

          for (auto &water : FluidObject)
            water->setContactFaces(Join(RigidBodyObject, SoftBodyObject));

          /*
        RKで節点が大きく移動する可能性がるので，setBoundaryTypesをここに移動した．
        RKの間，flipはされない．
        setBoundaryTypesは，接触面を更新するためにここに移動した．
        ただ，今は大きく動かないとしてはじめだけにしている．
        */
          std::cout << Green << "makeBuckets" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
          //! 体積を保存するようにリメッシュする必要があるだろう．
          for (auto &water : FluidObject) {
            setBoundaryTypes(water, Join(RigidBodyObject, SoftBodyObject));
            water->setMinDepthFromCORNER();
          }
          std::cout << Green << "setBoundaryTypes" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
        }

        // BIEで解く際も利用するので，吸収体に吸収される点を記録しておく
        for (const auto &p : Buckets_AllPoints.data1D) {
          p->absorbedBy = nullptr;
          for (const auto &net : AbsorberObject)
            if (net->InsideQ(p->X))
              p->absorbedBy = net;
        }

        /* --------------------------------------------------------------------------
         */
        /*                       境界値問題を解き，速度を計算 */
        /* --------------------------------------------------------------------------
         */

        setNeumannVelocity(Join(RigidBodyObject, SoftBodyObject));
        std::cout << Green << "setNeumannVelocity" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
        auto [time_setIGIGn_, time_solve_, unknownsize_] = BVP.solve();
        time_setIGIGn = time_setIGIGn_;
        time_solve = time_solve_;
        ElapsedTimeBIEDiscretization.push_back(time_setIGIGn);
        ElapsedTimeSolve.push_back(time_solve);
        unknownsizse = unknownsize_;
        std::cout << Green << "BVP.solve -> {Φ,Φn}が決まる" << Blue << "\nElapsed time: " << Red << watch() << " s\n";

        for (auto water : FluidObject)
          calculateCurrentVelocities(*water);
        std::cout << Green << "U_BEMとU_update_BEMを計算" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

        bool do_ALE = (RK_step == RK_order && (time_step % ALEPERIOD == 0));
        if (do_ALE) {
          for (auto water : FluidObject)
            calculateCurrentUpdateVelocities(*water, 30, 0.01);
          ElapsedTimeALE = watch()[0];
          std::cout << Green << "ALEのU_update_BEMを計算" << Blue << "\nElapsed time: " << Red << ElapsedTimeALE << colorReset << " s\n";
        }

        /* --------------------------------------------------------------------------
         */
        /*                                加速度の計算 */
        /* --------------------------------------------------------------------------
         */

#pragma omp parallel for
        for (const auto &net : RigidBodyObject) { // \label{BEM:impose_velocity}
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
              net->COM = net->RK_COM.getX();
              net->Q = Normalize(net->RK_Q.getX());
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

          // 人工的な粘性で結果が一致するようになるかどうかチェックする．

          bool use_given_velocity = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] != "update" && net->inputJSON.at("velocity")[0] != "floating");
          bool update_velocity_using_predetermined_accel = (net->inputJSON.find("velocity") && net->inputJSON.find("acceleration") && net->inputJSON.at("velocity")[0] == "update");
          bool update_velocity_using_solved_accel = (net->inputJSON.find("acceleration") && net->inputJSON.at("acceleration")[0] == "floating");
          bool update_velocity_using_solved_accel2 = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating");
          if (use_given_velocity) {
            std::cout << "use " << net->getName() << "'s (RigidBodyObject) predetermiend velocity" << std::endl;
          } else if (update_velocity_using_solved_accel || update_velocity_using_predetermined_accel || update_velocity_using_solved_accel2) {
            std::cout << "updating " << net->getName() << "'s (RigidBodyObject) velocity from acceleration" << std::endl;
            net->velocity = net->RK_Velocity.getX();
            std::cout << "acceleration = " << net->acceleration << std::endl;
            std::cout << "velocity = " << net->velocity << std::endl;
          } else {
            std::cout << net->getName() << "'s (RigidBodyObject) velocity is not updated" << std::endl;
          }
        }

        /* --------------------------------------------------------------------------
         */

        std::vector<Network *> movableObjects;
        for (auto &net : Join(RigidBodyObject, SoftBodyObject))
          if (std::ranges::none_of(net->isFixed, [](const auto &v) { return v == true; })) {
            movableObjects.push_back(net);
          }

        convergence = BVP.solveForPhiPhin_t(movableObjects);
        std::cout << Green << "BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

        // b$
        // -----------------------------------------------------------------------
        // */

        /*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

        ### 浮体の重心位置・姿勢・速度の更新

        浮体の重心位置は，重心に関する運動方程式を解くことで求める．
        姿勢は，角運動量に関する運動方程式などを使って，各加速度を求める．姿勢はクオータニオンを使って表現する．

        */

#pragma omp parallel for
        for (const auto &net : RigidBodyObject) { // \label{BEM:impose_velocity}
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
              net->RK_COM.push(net->velocityTranslational());
              net->RK_Q.push(net->Q.AngularVelocityTodQdt(net->velocityRotational()));
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

          // 人工的な粘性で結果が一致するようになるかどうかチェックする．

          bool use_given_velocity = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] != "update" && net->inputJSON.at("velocity")[0] != "floating");
          bool update_velocity_using_predetermined_accel = (net->inputJSON.find("velocity") && net->inputJSON.find("acceleration") && net->inputJSON.at("velocity")[0] == "update");
          bool update_velocity_using_solved_accel = (net->inputJSON.find("acceleration") && net->inputJSON.at("acceleration")[0] == "floating");
          bool update_velocity_using_solved_accel2 = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating");
          if (use_given_velocity) {
            std::cout << "use " << net->getName() << "'s (RigidBodyObject) predetermiend velocity" << std::endl;
          } else if (update_velocity_using_solved_accel || update_velocity_using_predetermined_accel || update_velocity_using_solved_accel2) {
            std::cout << "updating " << net->getName() << "'s (RigidBodyObject) velocity from acceleration" << std::endl;
            net->RK_Velocity.push(net->acceleration);
            std::cout << "acceleration = " << net->acceleration << std::endl;
            std::cout << "velocity = " << net->velocity << std::endl;
          } else {
            std::cout << net->getName() << "'s (RigidBodyObject) velocity is not updated" << std::endl;
          }

          //! pushしていなので，上にもっていっても変化しないはず
          for (const auto &p : net->getPoints())
            p->setXSingle(net->rigidTransformation(p->initialX));
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
            p->RK_X.push(V); //@ 位置xの時間発展
                             // p->setXSingle(p->RK_X.getX());
          }

          net->setGeometricProperties();
        }

        // b$
        // --------------------------------------------------------------------------
        // */ b$                                波の吸収（ダンピング領域） b$
        // --------------------------------------------------------------------------
        // */

        /*DOC_EXTRACT 0_4_1_UPDATE_POSITION

        ### 流体の$\phi$時間発展，$\phi_n$の時間発展はない

        ### 波の吸収（ダンピング領域）

        $$
        \begin{aligned}
        \gamma &= 1 - 2 \frac{\text{horizontal distance from the center of the
        absorber}}{\text{width of the absorber}} \\
        \phi_{\rm ref} &= \frac{\sum \phi \cdot \text{area}}{\sum \text{area}}
        \end{aligned}
        $$
        */

        // signed distanceの計算
#pragma omp parallel for
        for (const auto &p : ToVector(Buckets_AllPoints.data1D)) {
          if (p->absorbedBy != nullptr) {
            double min_distance = 1E+20;
            p->signed_distance = Norm(p->X - p->absorbedBy->NearestSurfacePoint(p->X));
          } else
            p->signed_distance = 0;
        }

        // phiの平均
        double mean_phi = 0., total_area = 0;
        for (const auto &f : Buckets_Surfaces.data1D) {
          auto [p0, p1, p2] = f->getPoints();
          mean_phi += (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3 * f->area;
          total_area += f->area;
        }
        mean_phi /= total_area;

        for (const auto &p : Buckets_AllPoints.data1D) {
          p->U_absorbed.fill(0.);
          double gamma = 0, ref_phi = 0;
          Tddd ref_U = {0, 0, 0};
          // U_update_BEMの修正
          if (p->absorbedBy != nullptr) {
            gamma = p->absorbedBy->absorb_gamma(p);
            if (std::ranges::any_of(p->getSurfaces(), [](auto f) { return f->Dirichlet; })) {
              auto nextX = p->RK_X.getX(p->U_update_BEM);
              auto to_eta_in_z = p->absorbedBy->absorb_eta(nextX, p->RK_X.getTimeAtNextStep()) - nextX[2];
              p->U_absorbed[2] = gamma * to_eta_in_z / p->RK_X.getdt();
              p->U_update_BEM[2] += p->U_absorbed[2];
              //
              ref_phi = p->absorbedBy->absorb_phi(p->RK_X.getX(p->U_update_BEM), p->RK_X.getTimeAtNextStep()) + mean_phi;
            }
          }
          if (!p->Neumann /*Neumannを変更しても，あとでBIEによって上書きされるので，からわない．*/) {
            p->RK_phi.push(p->DphiDt_damped({gamma, ref_phi}, p->U_update_BEM, 0.));
            std::get<0>(p->phiphin) = p->phi_Dirichlet = p->RK_phi.getX(); // 角点の法線方向はわからないので，ノイマンの境界条件phinを与えることができない．
          }
          //@ 位置xの時間発展
          p->RK_X.push(p->U_update_BEM);
          p->setXSingle(p->RK_X.getX());
          p->phi_tmp = 0;
        }

        // b$
        // --------------------------------------------------------------------------
        // */

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
      } while (!((*(Buckets_AllPoints.data1D).begin())->RK_X.finished));

      std::cout << Red << "Total elapsed time: " << Red << watch()[1] << colorReset << " s\n";
      ElapsedTimeTotal = watch()[1];

      /* ---------------------------------- 所要時間
       * ---------------------------------- */
      std::cout << "==========================================================="
                   "================="
                << std::endl;
      std::cout << "ElapsedTimeBIEDiscretization=" << ElapsedTimeBIEDiscretization << std::endl;
      std::cout << "ElapsedTimeSolve=" << ElapsedTimeSolve << std::endl;
      std::cout << "ElapsedTimeALE=" << ElapsedTimeALE << std::endl;
      std::cout << "ElapsedTimeTotal=" << ElapsedTimeTotal << std::endl;
      std::cout << "==========================================================="
                   "================="
                << std::endl;
      /* --------------------------------------------------------------------------
       */

      //! 速度ポテンシャルの平均を0にする
      {
        double mean_phi = 0., total_area = 0;
        for (const auto &f : Buckets_Surfaces.data1D) {
          auto [p0, p1, p2] = f->getPoints();
          mean_phi += (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3 * f->area;
          total_area += f->area;
        }
        mean_phi /= total_area;
        // mean_phi /= Buckets_AllPoints.data1D.size();
        for (const auto &p : Buckets_AllPoints.data1D) {
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
      simulation_time = (*(Buckets_AllPoints.data1D).begin())->RK_X.gett();

      /* ------------------------------------------------------ */
      for (auto water : FluidObject) {
        auto name = water->getName();
        std::ofstream ofs(output_directory / (name + ".obj"));
        createOBJ(ofs, *water);
        ofs.close();
      }

      /* --------------------------------------------------------------------------
       */
      /*                                 OUTPUT */
      /* --------------------------------------------------------------------------
       */

      // ここから出力を集約
      OutputContext ctx{.dt = dt, .time_step = time_step, .simulation_time = simulation_time, .output_directory = output_directory, .cpu_clock_start = cpu_clock_start, .wall_clock_start = wall_clock_start};

      // JSON 出力
      OutputJSON::write_step(ctx, jsonout, FluidObject, RigidBodyObject, SoftBodyObject, MeasurementJSONs, Buckets_Surfaces.data1D, convergence, unknownsizse, time_setIGIGn, time_solve);

      // ParaView 出力 (VTU/VTP + PVD)
      OutputParaView::write_step(ctx, NetOutputInfo, FluidObject, RigidBodyObject, SoftBodyObject, Buckets_Surfaces.data1D);

      /* --------------------------------------------------------------------------
       */
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
