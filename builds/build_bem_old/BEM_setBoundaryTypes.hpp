#ifndef BEM_setBoundaryConditions_H
#define BEM_setBoundaryConditions_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

/*DOC_EXTRACT 0_1_BOUNDARY_CONDITIONS

## 境界のタイプを決定する

<img src="./img/schematic_boundary_types_without_float.png" width="600px">

0. 流体と物体の衝突を判定し，流体節点が接触する物体面を保存しておく．

   * \ref{contact_angle}{`networkPoint::contact_angle`}
   * \ref{isInContact}{`networkPoint::isInContact`}
   * \ref{addContactFaces}{`networkPoint::addContactFaces`}

を使って接触判定を行っている．

 \ref{BEM:detection_range}{流体が構造物との接触を感知する半径}の設置も重要．

つぎに，その情報を使って，境界のタイプを次の順で決める．（物理量を与えるわけではない）

1. 面の境界条件：３節点全てが接触している流体面はNeumann面，それ以外はDirichlet面とする．CORNER面は設定しない．
   - Neumann面$`\Gamma^{({\rm N})}`$ : 3点接触流体面
   - Dirichlet面$`\Gamma^{({\rm D})}`$ : それ以外の面

2. 辺の境界条件 : 辺を含む２面がNeumann面ならNeumann辺，２面がDirichlet面ならDirichlet辺，それ以外はCORNERとする．
   - Neumann辺 : 隣接面2面がNeumann面の辺
   - Dirichlet辺 : 隣接面2面がDirichlet面の辺
   - CORNER辺 : それ以外の辺（Neumann面とDirichlet面の間にある辺）

3. 点の境界条件：点を含む面全てがNeumann面ならNeumann点，面全てがDirichlet面ならDirichlet点，それ以外はCORNERとする．
   - Neumann点 : 隣接面全てがNeumann面である点
   - Dirichlet点 : 隣接面全てがDirichlet面である点
   - CORNER点 : それ以外の点（Neumann面とDirichlet面の間にある点）

### 多重節点

NOTE: 面の向き$`\bf n`$がカクッと不連続に変わる節点には，$`\phi`$は同じでも，隣接面にそれぞれ対して異なる$`\phi_n`$を計算できるようにする

NOTE: $`\bf n`$が不連続に変化する節点まわりの要素は，自分のために用意された$`\phi_n`$を選択し補間に用いなければならない

これを多重節点という．

*/

void setRigidBodyVelocityAndAccel_IfPredetermined(Network *net, const double &RK_time) {

   if (net->isFixed) {
      net->mass = 1E+20;
      net->inertia.fill(1E+20);
      net->COM.fill(0.);
      net->initial_center_of_mass.fill(0.);
   }

   std::string move_name_velocity;
   T6d default_acceleration = {0., 0., 0., 0., 0., 0.};

   if (net->inputJSON.find("velocity")) {
      move_name_velocity = net->inputJSON["velocity"][0];
      std::cout << "move_name_velocity = " << move_name_velocity << std::endl;
      if (move_name_velocity == "update")
         std::cout << " velocity is already updated using acceleration" << std::endl;
      else if (move_name_velocity == "fixed")
         net->velocity.fill(0.);
      else if (move_name_velocity == "floating") {
         std::cout << "floatingの場合は，加速度の時間積分によってシミュレートされる" << std::endl;
         net->velocity = net->RK_Velocity.getX();
      } else {
         std::cout << "(RigidBodyObject) velocity is explicityly given as " << move_name_velocity << std::endl;
         double delta_t = 1E-5;
         if (move_name_velocity == "file") {
            net->velocity = net->intpMotionRigidBody.D(RK_time);
            default_acceleration = net->intpMotionRigidBody(RK_time + delta_t / 2.) - net->intpMotionRigidBody(RK_time - delta_t / 2.);
         } else {
            net->velocity = velocity(move_name_velocity, net->inputJSON["velocity"], RK_time);
            default_acceleration = velocity(move_name_velocity, net->inputJSON["velocity"], RK_time + delta_t / 2.) - velocity(move_name_velocity, net->inputJSON["velocity"], RK_time - delta_t / 2.);
         }
         default_acceleration /= delta_t;
      }
   } else {
      std::cout << "指定がないので速度はゼロ" << std::endl;
      net->velocity.fill(0.);
   }

   std::cout << "setting acceleration" << std::endl;
   std::string move_name_accel;
   if (move_name_velocity == "fixed")
      net->acceleration.fill(0.);
   else if (move_name_velocity == "floating")
      std::cout << "floatingの場合は，加速度は計算する" << std::endl;
   else if (net->inputJSON.find("acceleration")) {
      move_name_accel = net->inputJSON["acceleration"][0];
      if (move_name_accel == "fixed")
         net->acceleration.fill(0.);
      else if (move_name_accel == "floating") {
         std::cout << "floatingの場合は，加速度の時間積分によってシミュレートされる" << std::endl;
         // この時点ではわからない
      } else
         net->acceleration = acceleration(move_name_accel, net->inputJSON["acceleration"], RK_time);
   } else {
      std::cout << "指定がないので加速度はdefault_acceleration" << std::endl;
      net->acceleration = default_acceleration;
   }

   for (const auto &p : net->getPoints()) {
      auto tmp = net->velocityRigidBody(p->X);
      p->velocity[0] = tmp[0];
      p->velocity[1] = tmp[1];
      p->velocity[2] = tmp[2];
   }
};

// b# ------------------------------------------------------ */
// b#      物体のノイマン境界の速度 u(t) at Neumann を設定         */
// b# ------------------------------------------------------ */

//\label{BEM:setNeumannVelocity}
void setNeumannVelocity(const std::vector<Network *> &objects) {
   for (auto net : objects) {
      std::cout << Green << "setNeumannVelocity: " << colorReset << net->getName() << std::endl;
      //! 壁面の動きは，マイステップ更新することにした．この結果はphin()で参照される
      net->velocity.fill(0.);
      // net->acceleration.fill(0.);
      for (const auto &p : net->getPoints()) {
         p->velocity.fill(0.);
         // p->acceleration.fill(0.);
      }
      if (net->isRigidBody) {
         auto RK_time = net->RK_COM.gett();  //%各ルンゲクッタの時刻を使う
         setRigidBodyVelocityAndAccel_IfPredetermined(net, RK_time);
      } else if (net->isSoftBody) {
         std::cout << net->getName() << "　の流速の計算方法．soft bodyの場合，各節点に速度を与える．" << std::endl;
         if (net->inputJSON.find("velocity")) {
            std::string move_name = net->inputJSON["velocity"][0];
            std::cout << "move_name = " << move_name << std::endl;
            for (const auto &p : net->getPoints()) {
               auto RK_time = p->RK_X.gett();                                              //%各ルンゲクッタの時刻を使う
               p->velocity = velocity(move_name, net->inputJSON["velocity"], p, RK_time);  // T6d //@ Φnを計算するために，物体表面の速度forced_velocityは，保存しておく必要がある
            }
         }
      }
   }
}

// b# ------------------------------------------------------ */

void setIsMultipleNode(const auto &p) {
   if (p->CORNER)
      p->isMultipleNode = true;
   else if (p->Neumann) {
      auto neumann_N = p->getNormalNeumann_BEM();
      p->isMultipleNode = std::ranges::any_of(p->getFacesNeumann(), [&](const auto &f) { return !isFlat(neumann_N, f->normal, M_PI / 180. * 20); });
   } else
      p->isMultipleNode = false;
};

void setBoundaryTypes(Network &water, const std::vector<Network *> &objects) {
   std::cout << water.getName() << "の境界条件を決定 setBoundaryTypes" << std::endl;
   /* -------------------------------------------------------------------------- */
   /*                             f,l,pの境界条件を決定                             */
   /* -------------------------------------------------------------------------- */

   water.setGeometricProperties();

   std::cout << "step1 点の衝突の判定" << std::endl;
   for (const auto &p : water.getPoints())
      p->clearContactFaces();
   //!!! 衝突の判定がよくエラーが出る箇所
   for (const auto &net : objects) {
#pragma omp parallel
      for (const auto &p : water.getPoints())
#pragma omp single nowait
      {
         // // \label{BEM:detection_range}
         // double mean = 0., w = 0.;

         // for (const auto &f : p->getFaces()) {
         //    mean += f->area;
         //    w += 1.;
         // }

         // p->detection_range = std::sqrt(mean) / 3.;         //! ここも重要：点と面の衝突をどのようにすれば矛盾なく判定できるか．
         // \label{BEM:detection_range}
         p->detection_range = Mean(extLength(p->getLines())) / 4.;
         p->addContactFaces(net->getBucketFaces(), false);
      }
   }

   std::cout << "step2 面の境界条件を判定" << std::endl;
#pragma omp parallel
   for (const auto &f : water.getFaces())
#pragma omp single nowait
   {
      f->Neumann = std::ranges::all_of(f->getPoints(), [&f](const auto &p) { return isInContact(p, f, bfs(p->getContactFaces(), 2)); });
      f->Dirichlet = !f->Neumann;
   }

   std::cout << "step3 線の境界条件を決定" << std::endl;
   // #pragma omp parallel
   for (const auto &l : water.getLines())
   // #pragma omp single nowait
   {
      l->Neumann = std::ranges::all_of(l->getFaces(), [](const auto &f) { return f->Neumann; });
      l->Dirichlet = std::ranges::all_of(l->getFaces(), [](const auto &f) { return f->Dirichlet; });
      l->CORNER = (!l->Neumann && !l->Dirichlet);
   }
   std::cout << "step4 点の境界条件を決定" << std::endl;
#pragma omp parallel
   for (const auto &p : water.getPoints())
#pragma omp single nowait
   {
      p->Neumann = std::ranges::all_of(p->getFaces(), [](const auto &f) { return f->Neumann; });
      p->Dirichlet = std::ranges::all_of(p->getFaces(), [](const auto &f) { return f->Dirichlet; });
      p->CORNER = (!p->Neumann && !p->Dirichlet);
   }

   for (const auto &p : water.getPoints())
      setIsMultipleNode(p);
};

#endif