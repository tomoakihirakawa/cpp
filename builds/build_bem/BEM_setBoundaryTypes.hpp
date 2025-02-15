#pragma once

#include "Network.hpp"

//! 境界条件と関連が深いのでここで定義しておく

/*DOC_EXTRACT 0_1_1_2_BOUNDARY_CONDITIONS

## 多重節点

NOTE: 面の向き$`\bf n`$がカクッと不連続に変わる節点には，$`\phi`$は同じでも，隣接面にそれぞれ対して異なる$`\phi_n`$を計算できるようにする

NOTE: $`\bf n`$が不連続に変化する節点まわりの要素は，自分のために用意された$`\phi_n`$を選択し補間に用いなければならない

これを多重節点という．

多重節点を導入すると，未知変数idは，節点idだけではなく，節点と面の組みのidとなる．

### 境界値問題の未知変数ID 多重節点との区別

* `isNeumannID_BEM`と`isDirichletID_BEM`は，節点と面の組みが，境界値問題の未知変数かどうかを判定する．

* `pf2ID`は，節点と面の組みを未知変数IDに変換する．多重節点でない場合は，`{p,nullptr}`が変数のキーとなり，多重節点の場合は，与えられた`{p,f}`が変数のidとなる．

*/

bool isNeumannID_BEM(const netP *p, const netF *f) {
   if (p->Neumann || p->CORNER) {
      if (p->isMultipleNode)
         return p->MemberQ(f) && f->Neumann;
      else
         return (f == nullptr);
   } else
      return false;
};
bool isNeumannID_BEM(const std::tuple<netP *, netF *> &PF) { return isNeumannID_BEM(std::get<0>(PF), std::get<1>(PF)); };
bool isDirichletID_BEM(const auto p, const auto f) { return (p->Dirichlet || p->CORNER) && (f == nullptr); };
bool isDirichletID_BEM(const std::tuple<netP *, netF *> &PF) { return isDirichletID_BEM(std::get<0>(PF), std::get<1>(PF)); };

/*DOC_EXTRACT 0_3_BEM_utilities

## 多重節点を考慮したIDの設定方法

*/

//@ pf2IDは，setNodeFaceIndicesを実行せずとも使える．pf2Indexは，setNodeFaceIndicesを実行してから使う．
std::tuple<networkPoint *, networkFace *> pf2ID(const networkPoint *p, const networkFace *f) {
   if (isNeumannID_BEM(p, f) || isDirichletID_BEM(p, f))
      return {const_cast<networkPoint *>(p), const_cast<networkFace *>(f)};
   else
      return {const_cast<networkPoint *>(p), nullptr};
}

std::vector<std::tuple<networkPoint *, networkFace *>> p2AllIDs(const networkPoint *p) {
   std::vector<std::tuple<networkPoint *, networkFace *>> ret;
   bool nullptr_found = false;
   for (const auto &f : p->getSurfaces()) {
      auto PF = pf2ID(p, f);
      if (nullptr_found && std::get<1>(PF) != nullptr)
         ret.emplace_back(PF);
      else if (!nullptr_found)
         ret.emplace_back(PF);
      if (std::get<1>(PF) == nullptr)
         nullptr_found = true;
   }
   return ret;
};

/*

   係数行列を作成する場合（LU分解など）：
   pf2Index(p0, integ_f)は，積分の重みとに掛かる節点上のある量を指定のに使われる．
   これは，係数行列を作成する際に使うことになる．
   多重節点の場合でも適切にIDを返す．

   係数行列を作成する必要がない場合（GMRESなど）：
   もし，p->f2Index.at(f)が存在するなら，それは多重節点として扱われる．その値を使う．
   つまり，

*/

int pf2Index(const networkPoint *p, networkFace *f) {
   try {
      if (f == nullptr || !p->isMultipleNode || f->Dirichlet)
         return p->f2Index.at(nullptr);
      else
         return p->f2Index.at(f);
   } catch (const std::out_of_range &e) {
      for (const auto &[f, i] : p->f2Index)
         std::cout << f << " " << i << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error");
   }
};

/*

   pf2Index(p0, integ_f)を使えるように刷るためには，setNodeFaceIndicesを実行する必要がある．

*/

std::size_t setNodeFaceIndices(const std::vector<Network *> &objects) {
   std::size_t i = 0;
   for (const auto water : objects)
      for (const auto &q : water->getSurfacePoints()) {
         q->f2Index.clear();
         if (isNeumannID_BEM(q, nullptr) || isDirichletID_BEM(q, nullptr))
            q->f2Index[nullptr] = i++;
         for (const auto &f : q->getSurfaces())
            if (isNeumannID_BEM(q, f) || isDirichletID_BEM(q, f))
               q->f2Index[f] = i++;
      }
   return i;
};

std::size_t setNodeFaceIndices(const Network *objects) {
   return setNodeFaceIndices(std::vector<Network *>{const_cast<Network *>(objects)});
};

/* -------------------------------------------------------------------------- */

#include "BEM_utilities.hpp"

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_1_1_1_BOUNDARY_CONDITIONS

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

*/

void setRigidBodyVelocityAndAccel_IfPredetermined(Network *net, const double &RK_time) {

   if (std::ranges::all_of(net->isFixed, [](const auto &b) { return b; })) {
      net->mass = 1E+20;
      net->inertia.fill(1E+20);
      net->COM.fill(0.);
      net->initial_center_of_mass.fill(0.);
   } else {
      for (auto i = 0; i < net->isFixed.size(); i++) {
         if (net->isFixed[i])
            net->inertia[i] = 1E+20;
      }
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
// b#      物体のノイマン境界の速度 u(t) at Neumann を設定        */
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
         std::cout << net->getName() << "の流速の計算方法．soft bodyの場合，各節点に速度を与える．" << std::endl;
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

/* -------------------------------------------------------------------------- */
/*                             f,l,pの境界条件を決定                             */
/* -------------------------------------------------------------------------- */

void setBoundaryTypes(Network *water, const std::vector<Network *> &objects = {}) {
   std::cout << water->getName() << "の境界条件を決定 setBoundaryTypes" << std::endl;

   water->setGeometricProperties();

   std::cout << "step1 点の衝突の判定" << std::endl;
   for (const auto &p : water->getPoints())
      p->clearContactFaces();
   //!! 衝突の判定がよくエラーが出る箇所
   for (const auto &net : objects) {
#pragma omp parallel
      for (const auto &p : water->getPoints())
#pragma omp single nowait
      {
         //! ここも重要：点と面の衝突をどのようにすれば矛盾なく判定できるか．
         // \label{BEM:detection_range}
         double r = 0., s = 0.;
         for (auto &Q : p->getNeighbors())
            for (auto &q : Q->getNeighbors()) {
               r += Norm(q->X - p->X);
               s += 1.;
            }
         // p->detection_range = r / s / 2.;
         p->detection_range = r / s;
         p->addContactFaces(net->getBucketFaces(), false);
      }
   }

   // 接触判定の問題がか
   // nextvectorの問題？
   // tankを変えてやってみる

   std::cout << "step2 面の境界条件を判定" << std::endl;
   auto faces = water->getSurfaces();

   for (const auto &f : faces) {
      // f->Neumann = std::ranges::all_of(f->getPoints(), [&f](const auto &p) { return isInContact(p, f, bfs(p->getContactFaces(), 2)); });
      //$ faceがNeumannであるための条件は，faceのもつpointがすべて外部の面と接触していることである．
      //$ 以下の設定を使うことで，f->Neumannは自分の頂点pとで，p->getNearestContactFace(f)を使うことができる．
      f->Neumann = std::ranges::all_of(f->getPoints(), [&f](const auto &p) { return p->getNearestContactFace(f) != nullptr; });
      f->Dirichlet = !f->Neumann;
   }

   std::cout << "step3 線の境界条件を決定" << std::endl;
   for (const auto &l : water->getLines()) {
      l->Neumann = std::ranges::all_of(l->getSurfaces(), [](const auto &f) { return f->Neumann; });
      l->Dirichlet = std::ranges::all_of(l->getSurfaces(), [](const auto &f) { return f->Dirichlet; });
      l->CORNER = (!l->Neumann && !l->Dirichlet);
   }
   std::cout << "step4 点の境界条件を決定" << std::endl;

   for (const auto &p : water->getPoints()) {
      p->Neumann = std::ranges::all_of(p->getSurfaces(), [](const auto &f) { return f->Neumann; });
      p->Dirichlet = std::ranges::all_of(p->getSurfaces(), [](const auto &f) { return f->Dirichlet; });
      p->CORNER = (!p->Neumann && !p->Dirichlet);
      // step5 多重節点の判定
      setIsMultipleNode(p);
   }

   //@ ------------------------------------------ */

   for (const auto &f : faces) {
      // f->isPseudoQuadraticElement = _PSEUDO_QUADRATIC_ELEMENT_ && f->Dirichlet;
      f->isPseudoQuadraticElement = _PSEUDO_QUADRATIC_ELEMENT_;
      f->isLinearElement = !f->isPseudoQuadraticElement;
   }

   //@ ------------------------------------------ */

   std::cout << "setBoundaryTypes終了" << std::endl;
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

`phiOnFace`は，各節点`p`における各面`f`に対するポテンシャル`phi`を設定するために使用される．
`phitOnFace`は，各節点`p`における各面`f`に対するポテンシャルの時間微分`dphi/dt`を設定するために使用される．
他も同様である．

*/

/* -------------------------------------------------------------------------- */
/*                         phinOnFace, phintOnFaceの設定                       */
/* -------------------------------------------------------------------------- */
void setPhiPhinOnFace(Network *water) {
   // b! 点
   std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << colorReset << std::endl;

#pragma omp parallel
   for (const auto &p : water->getPoints())
#pragma omp single nowait
   {
      p->phiOnFace.clear();
      p->phitOnFace.clear();
      p->phinOnFace.clear();
      p->phintOnFace.clear();

      auto run = [p](networkFace *const &f) {
         /*Neumann境界のphinは，現在，接触しているオブジェクトの速度contactNormalVelocity(p)またはcontactNormalVelocity(p,f)を使って計算しており，前の時刻のphinは使っていない．*/
         if (p->absorbedBy != nullptr) {
            if (isNeumannID_BEM(p, f)) {
               if (f == nullptr) {
                  p->phiOnFace.insert({nullptr, std::get<0>(p->phiphin)});
                  p->phitOnFace.insert({nullptr, 1E+30});
                  p->phinOnFace.insert({nullptr, std::get<1>(p->phiphin) = Dot(p->absorbedBy->absorb_velocity(p), p->getNormalNeumann_BEM())});
                  p->phintOnFace.insert({nullptr, 1E+30});
               } else {
                  p->phiOnFace.insert({f, std::get<0>(p->phiphin)});
                  p->phitOnFace.insert({f, 1E+30});
                  p->phinOnFace.insert({f, Dot(p->absorbedBy->absorb_velocity(p), f->normal)});
                  p->phintOnFace.insert({f, 1E+30});
               }
            }
            if (isDirichletID_BEM(p, f)) {
               p->phiOnFace.insert({nullptr, std::get<0>(p->phiphin)});
               p->phitOnFace.insert({nullptr, 1E+30});
               p->phinOnFace.insert({nullptr, std::get<1>(p->phiphin)});
               p->phintOnFace.insert({nullptr, 1E+30});
            }
         } else {
            if (isNeumannID_BEM(p, f)) {
               if (f == nullptr) {
                  p->phiOnFace.insert({nullptr, std::get<0>(p->phiphin)});
                  p->phitOnFace.insert({nullptr, 1E+30});
                  p->phinOnFace.insert({nullptr, std::get<1>(p->phiphin) = Dot(contactNormalVelocity(p), p->getNormalNeumann_BEM())});
                  p->phintOnFace.insert({nullptr, 1E+30});
               } else {
                  p->phiOnFace.insert({f, std::get<0>(p->phiphin)});
                  p->phitOnFace.insert({f, 1E+30});
                  p->phinOnFace.insert({f, Dot(contactNormalVelocity(p, f), f->normal)});
                  p->phintOnFace.insert({f, 1E+30});
               }
            }
            if (isDirichletID_BEM(p, f)) {
               p->phiOnFace.insert({nullptr, std::get<0>(p->phiphin)});
               p->phitOnFace.insert({nullptr, 1E+30});
               p->phinOnFace.insert({nullptr, std::get<1>(p->phiphin)});
               p->phintOnFace.insert({nullptr, 1E+30});
            }
         }
      };

      run(nullptr);
      for (const auto &f : p->getSurfaces())
         run(f);
   }

   // b! 面
   std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << colorReset << std::endl;
   for (const auto &f : water->getSurfaces()) {
      auto [p0, p1, p2] = f->getPoints();
      std::get<0>(f->phiphin) = (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3.;
   }
};

void setPhiPhinOnFace(const std::vector<Network *> &objects) {
   for (const auto &water : objects)
      setPhiPhinOnFace(water);
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

// \label{BEM:setPhiPhin_t}
void setPhiPhin_t(std::vector<Network *> WATERS) {
   //    for (const auto water : WATERS)
   // #pragma omp parallel
   //       for (const auto &p : water->getPoints())
   // #pragma omp single nowait
   //          for (const auto &[F, i] : p->f2Index) {
   //             // auto [p, F] = PBF;
   //             //!!ノイマンの場合はこれでDphiDtは計算できませんよ
   //             if (isDirichletID_BEM(p, F))
   //                p->phitOnFace.at(F) = std::get<0>(p->phiphin_t) = p->aphiat(0.);
   //             else if (isNeumannID_BEM(p, F)) {
   //                for (auto &[f, phin_t] : p->phintOnFace) {
   //                   // phin_t = std::get<1>(p->phiphin_t) = (f != nullptr) ? phint_Neumann(f) : phint_Neumann(p);
   //                   phin_t = (f != nullptr) ? phint_Neumann(p, f) : phint_Neumann(p);  // \label{BEM:setphint}
   //                }
   //                std::get<1>(p->phiphin_t) = phint_Neumann(p);
   //             }
   //          }

   /* -------------------------------------------------------------------------- */
   // b! 点
   std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える" << colorReset << std::endl;

   for (const auto water : WATERS)
#pragma omp parallel
      for (const auto &p : water->getPoints())
#pragma omp single nowait
      {
         p->phitOnFace.clear();
         p->phintOnFace.clear();

         auto run = [p](networkFace *const &f) {
            /*Neumann境界のphinは，現在，接触しているオブジェクトの速度contactNormalVelocity(p)またはcontactNormalVelocity(p,f)を使って計算しており，前の時刻のphinは使っていない．*/
            if (p->absorbedBy != nullptr) {
               if (isNeumannID_BEM(p, f)) {
                  if (f == nullptr) {
                     p->phitOnFace.insert({nullptr, 1E+30});
                     p->phintOnFace.insert({nullptr, std::get<1>(p->phiphin_t) = Dot(p->absorbedBy->absorb_gradPhi_t(p), p->getNormalNeumann_BEM())});
                  } else {
                     p->phitOnFace.insert({f, 1E+30});
                     p->phintOnFace.insert({f, Dot(p->absorbedBy->absorb_gradPhi_t(p), f->normal)});
                  }
               }
               if (isDirichletID_BEM(p, f)) {
                  p->phitOnFace.insert({nullptr, p->aphiat(0.)});
                  p->phintOnFace.insert({nullptr, 1E+30});
               }
            } else {
               if (isNeumannID_BEM(p, f)) {
                  if (f == nullptr) {
                     p->phitOnFace.insert({nullptr, 1E+30});
                     p->phintOnFace.insert({nullptr, std::get<1>(p->phiphin_t) = phint_Neumann(p)});
                  } else {
                     p->phitOnFace.insert({f, 1E+30});
                     p->phintOnFace.insert({f, phint_Neumann(p, f)});
                  }
               }
               if (isDirichletID_BEM(p, f)) {
                  p->phitOnFace.insert({nullptr, p->aphiat(0.)});
                  p->phintOnFace.insert({nullptr, 1E+30});
               }
            }
         };

         run(nullptr);
         for (const auto &f : p->getSurfaces())
            run(f);

         // b! 面
         // std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦtを与える．Neumann面にはΦntを与える．" << colorReset << std::endl;
         // for (const auto &f : water->getSurfaces()) {
         //    auto [p0, p1, p2] = f->getPoints();
         //    std::get<0>(f->phiphin_t) = (std::get<0>(p0->phiphin_t) + std::get<0>(p1->phiphin_t) + std::get<0>(p2->phiphin_t)) / 3.;
         // }
      }
   std::cout << "setPhiPhin_t終了" << std::endl;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

template <typename T1, typename T2, typename T3>
void storePhiPhinCommon(const std::vector<Network *> &WATERS, const V_d &ans, T1 phiphinProperty, T2 phiOnFaceProperty, T3 phinOnFaceProperty) {

   for (const auto water : WATERS)
      for (const auto &p : water->getPoints())
         for (const auto &[f, i] : p->f2Index) {
            if (isDirichletID_BEM(p, f)) {
               (p->*phinOnFaceProperty).at(f) = std::get<1>(p->*phiphinProperty) = ans[i];
               (p->*phiOnFaceProperty).at(f) = std::get<0>(p->*phiphinProperty);
            }
            if (isNeumannID_BEM(p, f))
               (p->*phiOnFaceProperty).at(f) = std::get<0>(p->*phiphinProperty) = ans[i];
         }

   //^ 隣接フェイスの面積で重み付けした phi 値の寄与を計算し，結果を phiphinProperty に格納
   for (const auto water : WATERS)
      for (const auto &p : water->getPoints())
         if (p->Neumann) {
            double total = 0;
            std::get<0>(p->*phiphinProperty) = 0;
            for (const auto &f : p->getSurfaces()) {
               total += f->area;
               if ((p->*phiOnFaceProperty).count(f))
                  std::get<0>(p->*phiphinProperty) += (p->*phiOnFaceProperty).at(f) * f->area;
               else
                  std::get<0>(p->*phiphinProperty) += (p->*phiOnFaceProperty).at(nullptr) * f->area;
            }
            std::get<0>(p->*phiphinProperty) /= total;
         }
}

void storePhiPhin(const std::vector<Network *> &WATERS, const V_d &ans) {
   storePhiPhinCommon(WATERS, ans, &networkPoint::phiphin, &networkPoint::phiOnFace, &networkPoint::phinOnFace);
}

void storePhiPhin_t(const std::vector<Network *> &WATERS, const V_d &ans) {
   storePhiPhinCommon(WATERS, ans, &networkPoint::phiphin_t, &networkPoint::phitOnFace, &networkPoint::phintOnFace);
}
