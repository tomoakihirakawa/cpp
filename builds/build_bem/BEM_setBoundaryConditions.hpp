#ifndef BEM_setBoundaryConditions_H
#define BEM_setBoundaryConditions_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

/*DOC_EXTRACT BEM

## 境界条件の設定

1. 流体節点が接触する構造物面を保存する
2. 面の境界条件：３節点全てが接触している流体面はNeumann面，それ以外はDirichlet面とする
3. 辺の境界条件：辺を含む２面がNeumann面ならNeumann辺，２面がDirichlet面ならDirichlet面，それ以外はCORNERとする．
4. 点の境界条件：点を含む面全てがNeumann面ならNeumann点，面全てがDirichlet面ならDirichlet点，それ以外はCORNERとする．

### 多重節点

NOTE: 面の向き$\bf n$がカクッと不連続に変わる節点には，$\phi$は同じでも，隣接面にそれぞれ対して異なる$\phi_n$を計算できるようにする

NOTE: $\bf n$が不連続に変化する節点まわりの要素は，自分のために用意された$\phi_n$を選択し補間に用いなければならない

これを多重節点という．

このループでは，ある面`integ_f`に隣接する節点{p0,p1,p2}の列,IGIGn[origin(fixed),p0],...に値が追加されていく．
（p0が多重接点の場合，適切にp0と同じ位置に別の変数が設定されており，別の面の積分の際にq0が参照される．）
p0は，{面,補間添字}で決定することもできる．
{面,補間添字0}->p0,{面,補間添字1}->p1,{面,補間添字2}->p2というように．

{面A,補間添字},{面B,補間添字},{面C,補間添字}が全て同じ節点p0を指していたとする．
普通の節点なら，IGIGn[origin,{p0,nullptr}]を指す．
多重節点なら，IGIGn[origin,{p0,面A}],IGIGn[origin,{p0,面B}]を指すようにする．
この操作を言葉で言い換えると，

面を区別するかどうかが先にわからないので，face*のまsまかnullptrとすべきかわからないということ．．．．

PBF_index[{p, Dirichlet, ある要素}]
は存在しないだろう．Dirichlet節点は，{p, ある要素}からの寄与を，ある面に

*/

void setIsMultipleNode(const auto &p) {
   if (p->CORNER)
      p->isMultipleNode = true;
   else if (p->Neumann) {
      auto n = p->getNormalNeumann_BEM();
      p->isMultipleNode = std::ranges::any_of(p->getFacesNeumann(), [&n](const auto &f) { return !isFlat(n, f->normal, M_PI / 180. * 20); });
   } else
      p->isMultipleNode = false;
};

void setBoundaryConditions(Network &water, const std::vector<Network *> &objects) {
   water.setGeometricProperties();
   /* -------------------------------------------------------------------------- */
   /*                             f,l,pの境界条件を決定                             */
   /* -------------------------------------------------------------------------- */
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
         // auto toF = extXtuple(ToVector(p->getFaces())) - ToX(p);
         // auto toP = extXtuple(p->getNeighbors()) - ToX(p);
         // double a = Norm(*std::min_element(toP.begin(), toP.end(), [](const auto &a, const auto &b) { return Norm(a) < Norm(b); }));
         // double b = Norm(*std::min_element(toF.begin(), toF.end(), [](const auto &a, const auto &b) { return Norm(a) < Norm(b); }));
         // p->radius = Mean(extLength(extractLines(Flatten(BFS(p, 2))))) / 5.;
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