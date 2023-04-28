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
};

#endif