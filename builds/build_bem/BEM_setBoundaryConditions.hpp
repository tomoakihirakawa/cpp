#ifndef BEM_setBoundaryConditions_H
#define BEM_setBoundaryConditions_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

void setBoundaryConditions(Network &water, const std::vector<Network *> &objects) {
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
         p->radius = (Mean(extLength(p->getLines())) + radius) / 2.;
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
      f->Neumann = all_of(f->getPoints(), [&f](const auto &p) { return isInContact(p, f, bfs(p->getContactFaces(), 2)); });
      f->Dirichlet = !f->Neumann;
   }

   std::cout << "step3 線の境界条件を決定" << std::endl;
#pragma omp parallel
   for (const auto &l : water.getLines())
#pragma omp single nowait
   {
      l->Neumann = all_of(l->getFaces(), [](const auto &f) { return f->Neumann; });
      l->Dirichlet = all_of(l->getFaces(), [](const auto &f) { return f->Dirichlet; });
      l->CORNER = (!l->Neumann && !l->Dirichlet);
   }
   std::cout << "step4 点の境界条件を決定" << std::endl;
#pragma omp parallel
   for (const auto &p : water.getPoints())
#pragma omp single nowait
   {
      p->Neumann = all_of(p->getFaces(), [](const auto &f) { return f->Neumann; });
      p->Dirichlet = all_of(p->getFaces(), [](const auto &f) { return f->Dirichlet; });
      p->CORNER = (!p->Neumann && !p->Dirichlet);
   }
   /* -------------------------------------------------------------------------- */
   /*                         phinOnFace, phintOnFaceの設定                         */
   /* -------------------------------------------------------------------------- */
   // b! 点
   /**
    ## 多重節点
    多重節点という名前は具体性に欠ける．
    普通φnは(節点)にのみ依存する変数だが，nの変化が急なため，不連続性が著しい節点においては，(節点に加え面)にも依存する変数を複数設定する．それらは離散化などで使い分けることになる．
    BIEの離散化における，多重節点扱いについて．
    BIEを数値的に解くために，十分な数の１次方程式を作成する．これは，節点と同じ位置にBIEの原点を取ることで実現できる．
    同じ位置であるにもかかわらず，(節点に加え面)にも依存する変数φnを設定した場合，
    同じ位置であるにもかかわらず，それらを一つ一つを原点として，１次方程式を作成する．
    これらは完全に同じ方程式である．変数の数を節点の数よりも増やしたことによって，方程式の数が増えている．
   */
   std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << colorOff << std::endl;
   auto multiple_node_if = [&](const auto &p, const auto &facesNeuman) {
      return (p->CORNER || any_of(facesNeuman, [&](const auto &f) { return !isFlat(p->getNormalNeumann_BEM(), f->normal, M_PI / 180. * 20); }));
   };
#pragma omp parallel
   for (const auto &p : water.getPoints())
#pragma omp single nowait
   {
      p->phinOnFace.clear();
      p->phintOnFace.clear();

      if (p->Neumann || p->CORNER) {
         auto facesNeuman = p->getFacesNeumann();
         if (multiple_node_if(p, facesNeuman)) {
            for (const auto &f : facesNeuman) {
               p->phinOnFace[f] = Dot(uNeumann(p, f), f->normal);
               p->phintOnFace[f] = 1E+30;
            }
         } else {
            p->phinOnFace[nullptr] = std::get<1>(p->phiphin) = Dot(uNeumann(p), p->getNormalNeumann_BEM());
            p->phintOnFace[nullptr] = 1E+30;
         }
      }
   }

   // b! 面
   std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << colorOff << std::endl;
#pragma omp parallel
   for (const auto &f : water.getFaces())
#pragma omp single nowait
   {
      if (f->Neumann) {
         std::get<1>(f->phiphin) = Dot(uNeumann(f), f->normal);
      } else {
         auto [p0, p1, p2] = f->getPoints();
         std::get<0>(f->phiphin) = (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3.;
      }
   }
};

#endif