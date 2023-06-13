#ifndef BEM_setBoundaryConditions_H
#define BEM_setBoundaryConditions_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

/*DOC_EXTRACT BEM

## 境界のタイプを決定する

0. 流体と物体の衝突を判定し，流体節点が接触する物体面を保存しておく．
   \ref{contact_angle}{`networkPoint::contact_angle`}，
   \ref{isInContact}{`networkPoint::isInContact`}，
   \ref{addContactFaces}{`networkPoint::addContactFaces`}
   を使って接触判定を行っている．

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

void setNeumannVelocity(const std::vector<Network *> &objects) {
   // b# ------------------------------------------------------ */
   // b#      物体のノイマン境界の速度 u(t) at Neumann を設定         */
   // b# ------------------------------------------------------ */
   for (const auto &net : objects) {
      //! 壁面の動きは，マイステップ更新することにした．この結果はphin()で参照される
      if (net->isRigidBody) {
         auto RK_time = net->RK_COM.gett();  //%各ルンゲクッタの時刻を使う
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
      if (net->isSoftBody) {
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
               for (const auto &p : net->getPoints()) {
                  auto RK_time = p->RK_X.gett();                                              //%各ルンゲクッタの時刻を使う
                  p->velocity = velocity(move_name, net->inputJSON["velocity"], p, RK_time);  // T6d //@ Φnを計算するために，物体表面の速度forced_velocityは，保存しておく必要がある
               }
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
   }
}

void setIsMultipleNode(const auto &p) {
   if (p->CORNER)
      p->isMultipleNode = true;
   else if (p->Neumann) {
      auto n = p->getNormalNeumann_BEM();
      p->isMultipleNode = std::ranges::any_of(p->getFacesNeumann(), [&n](const auto &f) { return !isFlat(n, f->normal, M_PI / 180. * 20); });
   } else
      p->isMultipleNode = false;
};

void setBoundaryTypes(Network &water, const std::vector<Network *> &objects) {

   /* -------------------------------------------------------------------------- */
   /*                             f,l,pの境界条件を決定                             */
   /* -------------------------------------------------------------------------- */

   water.setGeometricProperties();

   auto radius = Mean(extLength(water.getLines()));
   Print("makeBucketFaces", Green);
   for (const auto &net : objects) {
      radius = Mean(extLength(net->getLines()));
      net->makeBucketFaces(radius);
   }

   std::cout << "step1 点の衝突の判定" << std::endl;
   for (const auto &p : water.getPoints())
      p->clearContactFaces();
   //!!! 衝突の判定がよくエラーが出る箇所
   for (const auto &net : objects) {
#pragma omp parallel
      for (const auto &p : water.getPoints())
#pragma omp single nowait
      {
         //! ここも重要：点と面の衝突をどのようにすれば矛盾なく判定できるか．
         p->radius = (Mean(extLength(p->getLines())) + radius) / 3.;
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
#pragma omp parallel
   for (const auto &l : water.getLines())
#pragma omp single nowait
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