#ifndef BEM_derivatives_H
#define BEM_derivatives_H

#include "Network.hpp"

Tddd fitToNeumannVelocity(Tddd VECTOR, const networkPoint *const p) {
   if (p->Neumann || p->CORNER) {
      Tddd nf, vn0, vn1;
      double max_w = kernel_Bspline3(0., p->radius);
      for (const auto &[f, hit_X] : Reverse(p->getContactFacesXCloser()) /*遠い方から*/) {
         nf = f->normal;
         vn0 = Dot(VECTOR, nf) * nf;
         vn1 = Dot(f->getNetwork()->velocityRigidBody(hit_X), nf) * nf;
         VECTOR += kernel_Bspline3(Norm(hit_X - p->getXtuple()), p->radius) / max_w * (vn1 - vn0);
      }
   }
   return VECTOR;
};

Tddd condition_Ua(Tddd VECTOR, const networkPoint *const p) {
   if (p->CORNER) {
      auto cross = Normalize(Cross(p->getNormalNeumann_BEM(), p->getNormalDirichlet_BEM()));
      VECTOR = Dot(VECTOR, cross) * cross;
      // for (const auto &l : p->getLines())
      // 	if (l->CORNER)
      // 	{
      // 		Tddd dir = Normalize((*l)(p)->getXBuffer() - p->getXBuffer());
      // 		VECTOR = Dot(VECTOR, dir) * dir;
      // 	}
      for (const auto &f : p->getFaces())
         if (f->Neumann)
            VECTOR -= Dot(VECTOR, f->normal) * f->normal;
      return VECTOR;
   } else if (p->Dirichlet) {
      return VECTOR - Dot(VECTOR, p->getNormal_BEM()) * p->getNormal_BEM();
   } else {
      for (const auto &f : p->getFaces())
         VECTOR -= Dot(VECTOR, f->normal) * f->normal;
      return VECTOR;
   }
};

Tddd vectorsToSurfaceFromBufferX(const networkPoint *p, const std::vector<T3Tddd> &next_Vrtx) {
   auto closestXFacing = [](const Tddd &p_next_X, const double radius, const std::vector<T3Tddd> &vertices, const Tddd &n) {
      Tddd r = {1E+100, 1E+100, 1E+100};
      for (const auto &vertex : vertices) {
         if (isFacing(TriangleNormal(vertex), n, M_PI / 180 * 20)) {
            auto intxn = IntersectionSphereTriangle(p_next_X, radius, vertex);
            if (intxn.isIntersecting)
               if (Norm(r) >= Norm(intxn.X - p_next_X))
                  r = intxn.X - p_next_X;
         }
      }
      return r;
   };

   if (next_Vrtx.empty())
      return {0., 0., 0.};
   else if (p->Neumann || p->CORNER) {
      std::vector<Tddd> F_clings;
      for (const auto &f : p->getFacesNeumann()) {
         // auto [p0, p1, p2] = f->getPoints();
         // auto n = TriangleNormal(p0->getXBuffer() + p0->U_BUFFER, p1->getXBuffer() + p1->U_BUFFER, p2->getXBuffer() + p2->U_BUFFER);
         auto n = f->normal;
         auto to_closest_X = closestXFacing(p->getXBuffer() + p->U_BUFFER, p->radius, next_Vrtx, n);
         if (isFinite(to_closest_X))
            F_clings.push_back(Dot(to_closest_X, n) * n);
      }
      // Tddd r = Mean(F_clings);
      Tddd r = optimumVector_(F_clings, {0., 0., 0.});
      if (isFinite(r))
         return r;
      else
         return {0., 0., 0.};
   } else
      return {0., 0., 0.};
};

double minViewRatio(const networkPoint *const p) {
   double a = p->getSolidAngle();
   return (2 * M_PI - Min(Tdd{std::abs(a - 2 * M_PI), std::abs(2 * M_PI - a)})) / (2 * M_PI);
};

double normalVariance(const networkPoint *const p) {
   auto n = p->getNormalDirichlet_BEM();
   double m = 0, s = 0;
   for (const auto &f : p->getFacesDirichlet()) {
      m += (M_PI / 2. - VectorAngle(n, f->normal)) / (M_PI / 2.);
      s += 1;
   }
   return m / s;
};

Tddd vectorTangentialShift(const networkPoint *p) {
   auto nextX_U_Ua = [](const networkPoint *p) {
      return p->getXBuffer() + p->U_BUFFER;
   };

   auto next_length = [nextX_U_Ua](const networkLine *const l) {
      auto [p0, p1] = l->getPoints();
      return Norm(nextX_U_Ua(p0) - nextX_U_Ua(p1));
   };

   auto getBaseLength = [next_length](const networkLine *line) {
      // auto [p0, p1] = line->getPoints();
      // std::unordered_set<networkLine *> lc = Join(extLinesCORNER_(p0->getFaces()), extLinesCORNER_(p1->getFaces()));
      // if (!lc.empty()) {
      //    V_d ll;
      //    for (const auto &l : lc)
      //       ll.emplace_back(next_length(l));
      //    return Mean(ll);
      // } else {
      //    // pを引っ張る力は，ノイマン面とディリクレ面で干渉しない
      //    auto [p0, p1] = line->getPoints();
      //    double m = 1, s = 0;
      //    for (const auto &l : Join(p0->getLinesAround(), p1->getLinesAround()))
      //       if (line != l)
      //          if ((line->Dirichlet && (l->Dirichlet || l->CORNER)) || (line->Neumann && (l->Neumann || l->CORNER)) || (line->CORNER && l->CORNER)) {
      //             m *= next_length(l);
      //             s += 1;
      //          }
      //    return std::pow(m, 1. / s);
      // }
      /* ---------------------------------------------- */
      auto [p0, p1] = line->getPoints();
      double m = 1, s = 0;
      std::unordered_set<networkLine *> lines;
      for (const auto &l : Join(p0->getLinesAround(), p1->getLinesAround()))
         lines.emplace(l);
      for (const auto &l : lines) {
         // if ((line->Dirichlet && (l->Dirichlet || l->CORNER)) || (line->Neumann && (l->Neumann || l->CORNER)) || (line->CORNER && l->CORNER)) {
         //    m *= next_length(l);
         //    s += 1;
         // }
         auto [p0, p1] = l->getPoints();
         // auto wiehgt = RootMeanSquare(Tdd{p0->getSolidAngle(), p1->getSolidAngle()} - 2 * M_PI) / (2. * M_PI) + 1.;
         auto wiehgt = 1;
         m *= next_length(l) * wiehgt;
         s += wiehgt;
      }
      return m / s;
   };

   auto vectorToNextNeighborsCenter = [nextX_U_Ua](const networkPoint *const p) {
      double s = 0;
      Tddd ret = {0., 0., 0.};
      Tddd pX = nextX_U_Ua(p);
      /* ------------------------------------------------------ */
      if (p->CORNER) {
         for (const auto &l : p->getLines())
            if (l->CORNER) {
               ret += (nextX_U_Ua(((*l)(p))) - pX);
               s += 1;
            }
      } else {
         for (const auto &l : p->getLines()) {
            ret += (nextX_U_Ua(((*l)(p))) - pX);
            s += 1;
         }
      }
      return ret / s;
   };

   auto vectorToCenter = [nextX_U_Ua](const networkPoint *const p) {
      double s = 0;
      Tddd vec = {0., 0., 0.}, pX = nextX_U_Ua(p);
      if (p->CORNER) {
         for (const auto &l : p->getLinesCORNER()) {
            vec += nextX_U_Ua((*l)(p)) - pX;
            s += 1.;
         }
      } else {
         for (const auto &f : p->getFaces()) {
            auto [p0, p1, p2] = f->getPoints(p);
            // vec += Mean(T3Tddd{nextX_U_Ua(p0), nextX_U_Ua(p1), nextX_U_Ua(p2)}) - pX;
            // s += 1;
            vec += (Incenter(nextX_U_Ua(p0), nextX_U_Ua(p1), nextX_U_Ua(p2)) - pX) * Inradius(nextX_U_Ua(p0), nextX_U_Ua(p1), nextX_U_Ua(p2));
            s += Inradius(nextX_U_Ua(p0), nextX_U_Ua(p1), nextX_U_Ua(p2));
         }
      }
      return vec / s;
   };

   /* ------------------------------------------------------ */
   Tddd V = {0., 0., 0.};
   Tddd v;
   for (const auto &f : p->getFaces()) {
      auto [p0, p1, p2] = f->getPoints(p);
      double d0, d1, d2;
      Tddd v0, v1, v2;
      {
         auto intersect = IntersectionSphereLine(nextX_U_Ua(p0), 1E+20, T2Tddd{nextX_U_Ua(p1), nextX_U_Ua(p2)});
         d0 = intersect.distance;
         v0 = intersect.X - nextX_U_Ua(p0);
      }
      {
         auto intersect = IntersectionSphereLine(nextX_U_Ua(p1), 1E+20, T2Tddd{nextX_U_Ua(p2), nextX_U_Ua(p0)});
         d1 = intersect.distance;
         v1 = intersect.X - nextX_U_Ua(p1);
      }
      {
         auto intersect = IntersectionSphereLine(nextX_U_Ua(p2), 1E+20, T2Tddd{nextX_U_Ua(p0), nextX_U_Ua(p1)});
         d2 = intersect.distance;
         v2 = intersect.X - nextX_U_Ua(p2);
      }
      double mean = (d0 + d1 + d2) / 3.;
      // double var = Norm(Tddd{d0 - mean, d1 - mean, d2 - mean}) / mean;
      V += v0;
      V -= v1;
      V -= v2;
   }
   /* ------------------------------------------------------ */
   Tddd pX = nextX_U_Ua(p);
   double c_LS = 0.5 /*0.1~0.5*/, c_EMT = 0., c_AL = 0.;
   Tddd V_EMT = {0., 0., 0.};
   auto a = minViewRatio(p);
   double tmp = 1, s = 0;
   auto V_LS = vectorToCenter(p);

   return condition_Ua(V_LS * c_LS, p);

   /* --------------------------------------------------------- */
   // if (!p->CORNER) {
   //    return condition_Ua(V_LS * c_LS, p);
   // } else {
   //    for (const auto &l : p->getLines())
   //       V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
   //    return condition_Ua(V_EMT * c_EMT + c_AL * V + V_LS * c_LS, p);
   // }
   /* --------------------------------------------------------- */
   // if (p->CORNER) {
   //    if (a > 1. / 6.)
   //       return condition_Ua(V_LS * c_LS, p);
   //    else
   //       return {0., 0., 0.};
   // } else if (p->Dirichlet) {
   //    for (const auto &l : p->getLines())
   //       V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
   //    return condition_Ua(V_EMT * c_EMT + c_AL * V + V_LS * c_LS, p);
   // } else {
   //    // if (a > 1. / 3.) //比較的滑らかなノイマン面

   //    {
   //       for (const auto &l : p->getLines())
   //          V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
   //       return condition_Ua(V_EMT * c_EMT + c_AL * V + V_LS * c_LS, p);
   //    }
   //    // else //比較的滑らかでないノイマン面
   //    // {
   //    //     V_netFp f;
   //    //     for (const auto &l : p->getLines())
   //    //     {
   //    //         f = l->getFaces();
   //    //         auto b = VectorAngle(f[0]->normal, f[1]->normal) / (2 * M_PI);
   //    //         // std::cout << a << std::endl;
   //    //         V_EMT += (nextX_U_Ua((*l)(p)) - pX) * b;
   //    //     }
   //    //     return condition_Ua(V_EMT, p);
   //    // }
   // }
};

std::tuple<Tddd, double, double, networkLine *, networkFace *> vectorToNearestAjacentSurface(const networkPoint *p) {
   /*
   UartificialClingは，完璧にOmega(t+\delta t)に張り付くようにしなければ，
   面からはなれることで計算の破綻を招く可能性がある．
   */
   std::tuple<Tddd, double, double, networkLine *, networkFace *> ret = {{0, 0, 0}, 0., 0., nullptr, nullptr};
   Tddd pX = p->getXBuffer() + p->U_BUFFER;
   if (p->CORNER) {
      // auto a = minViewRatio(p);
      // if (a > 1 / 3)
      for (const auto &l : p->getLines())
         if (l->CORNER) {
            auto intxn = IntersectionSphereLine(pX, 1E+20, T2Tddd{p->getXBuffer(), (*l)(p)->getXBuffer()});
            if (Norm(std::get<0>(ret) - pX) >= Norm(intxn.X - pX)) {
               std::get<0>(ret) = intxn.X;
               std::get<1>(ret) = intxn.t;
               std::get<2>(ret) = 1 - intxn.t;
               std::get<3>(ret) = l;
               std::get<4>(ret) = nullptr;
            }
         }
   } else {
      for (const auto &f : p->getFaces()) {
         if ((p->Dirichlet && f->Dirichlet) || (p->Neumann && f->Neumann)) {
            auto [p0, p1, p2] = f->getPoints(p);
            auto intxn = IntersectionSphereTriangle(pX, 1E+20, T3Tddd{p0->getXBuffer(), p1->getXBuffer(), p2->getXBuffer()});
            if (Norm(std::get<0>(ret) - pX) >= Norm(intxn.X - pX)) {
               std::get<0>(ret) = intxn.X;
               std::get<1>(ret) = intxn.t0;
               std::get<2>(ret) = intxn.t1;
               std::get<3>(ret) = nullptr;
               std::get<4>(ret) = f;
            }
         }
      }
   }
   return ret;
};

void calculateVectorToSurfaceInBuffer(const Network &net, const bool adjust_dirichlet = true) {
   /*
   @ この方法なら，次の時刻における任意の場所でのポテンシャルを見積もることができる．
   @ このことは，任意のノイマン面上に節点を維持する上で便利である．
   @ Ω(t+δt)をまず見積もり，その面上で最適な格子配置となるように流速を修正する．
   */
   // auto Points = ToVector(net.getPoints());
   auto &Points = (net.getPoints());
   for (const auto &p : Points)
      p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};
   //@ ------------------------------------------------------ */
   //@           次の時刻で最適な格子を目指す修正流速を計算          */
   //@ ------------------------------------------------------ */
   for (auto kk = 0; kk < 20; ++kk) {
      //% ------------------------------------------------------ */
      //%           　　　　 vectorTangentialShift   　 　         */
      //%          ラプラス平滑化と引っ張り合わせた接線方向にシフト      */
      //% ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
         if (!p->Dirichlet || (p->Dirichlet && adjust_dirichlet)) {
            p->U_BUFFER_BUFFER = vectorTangentialShift(p);
         }

      for (const auto &p : Points) {
         if (isFinite(p->U_BUFFER_BUFFER))
            p->U_BUFFER += p->U_BUFFER_BUFFER;
         p->U_BUFFER_BUFFER = {0., 0., 0.};
      }

      //% ------------------------------------------------------ */
      //%              vectorToNearestAjacentSurface             */
      //%           　　　   周辺ディリクレ面に移動        　　　　    */
      //% ------------------------------------------------------ */

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
         if (p->CORNER || (p->Dirichlet && adjust_dirichlet)) {
            if (isFinite(p->U_BUFFER_BUFFER)) {
               p->clungSurface = vectorToNearestAjacentSurface(p);
               p->U_BUFFER_BUFFER = (std::get<0>(p->clungSurface) - (p->getXBuffer() + p->U_BUFFER));
               //! 角のノイマン面から離れるのを防ぐ
               if (p->CORNER)
                  for (const auto &f : p->getFacesNeumann())
                     p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, f->normal) * f->normal;
            }
         }

      for (const auto &p : Points)
         if (p->CORNER || (p->Dirichlet && adjust_dirichlet)) {
            if (isFinite(p->U_BUFFER_BUFFER, 1E+10))
               p->U_BUFFER += p->U_BUFFER_BUFFER;
            p->U_BUFFER_BUFFER = {0., 0., 0.};
         }

      //% ------------------------------------------------------ */
      //%             vectorsToSurfaceFromBufferX                */
      //%              　　　近傍のノイマン面へ移動           　　　   */
      //% ------------------------------------------------------ */
      for (auto ii = 0; ii < 2; ++ii) {
#ifdef _OPENMP
#pragma omp parallel
#endif
         for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
            if (p->CORNER || p->Neumann) {
               {
                  // 接触面候補の次の時刻の位置を予測
                  std::unordered_set<networkFace *> Fs = p->getContactFaces();
                  for (auto &f : p->getContactFaces())
                     for_each(f->getPoints(),
                              [&](const auto &q) {for (auto &F : q->getFaces()){Fs.emplace(F);} });
                  std::vector<T3Tddd> nextBodyVertices;
                  for (auto &f : Fs) {
                     auto [p0, p1, p2] = f->getPoints();
                     // auto X0 = p0->getXtuple() + f->getNetwork()->velocityRigidBody(p0->getXtuple()) * p->RK_X.getdt();
                     // auto X1 = p1->getXtuple() + f->getNetwork()->velocityRigidBody(p1->getXtuple()) * p->RK_X.getdt();
                     // auto X2 = p2->getXtuple() + f->getNetwork()->velocityRigidBody(p2->getXtuple()) * p->RK_X.getdt();
                     /* ------------------------------------------------------ */
                     auto net = f->getNetwork();
                     Quaternion q;
                     q = q.d_dt(net->velocityRotational());
                     //
                     auto COM = net->RK_COM.getX(net->velocityTranslational());
                     auto Q = net->RK_Q.getX(q());
                     auto X0 = velocityRigidBody(p0, COM, Q);
                     auto X1 = velocityRigidBody(p1, COM, Q);
                     auto X2 = velocityRigidBody(p2, COM, Q);
                     /* ------------------------------------------------------ */
                     nextBodyVertices.emplace_back(T3Tddd{X0, X1, X2});
                  }
                  p->U_BUFFER_BUFFER = vectorsToSurfaceFromBufferX(p, nextBodyVertices);
                  //! 角のディリクレ面へのめり込みを防止
                  if (p->CORNER)
                     p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, p->getNormalDirichlet_BEM()) * p->getNormalDirichlet_BEM();
               }
            }
      }
      for (const auto &p : Points) {
         if (isFinite(p->U_BUFFER_BUFFER, 1E+10))
            p->U_BUFFER += p->U_BUFFER_BUFFER;
         p->U_BUFFER_BUFFER = {0., 0., 0.};
      }
   }
};

void calculateVectorFromBufferToContactFaces(const Network &net) {
   /*
   @ ノイマン面に貼り付けるための必要な調整
   */
   // auto Points = ToVector(net.getPoints());
   auto &Points = net.getPoints();

   // for (const auto &p : Points)
   // 	p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};
   for (auto kk = 0; kk < 20; ++kk) {
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
      {
         if (p->CORNER || p->Neumann) {
            // 接触面候補の次の時刻の位置を予測
            std::unordered_set<networkFace *> Fs = p->getContactFaces();
            for (auto &f : p->getContactFaces())
               for_each(f->getPoints(),
                        [&](const auto &q) {for (auto &F : q->getFaces()){Fs.emplace(F);} });
            // for (auto &q : f->getPoints())
            //     for (auto &F : q->getFaces())
            //         Fs.emplace(F);

            std::vector<T3Tddd> nextBodyVertices;
            for (auto &f : Fs) {
               auto [p0, p1, p2] = f->getPoints();

               // auto X0 = p0->getXtuple() + f->getNetwork()->velocityRigidBody(p0->getXtuple()) * p->RK_X.getdt();
               // auto X1 = p1->getXtuple() + f->getNetwork()->velocityRigidBody(p1->getXtuple()) * p->RK_X.getdt();
               // auto X2 = p2->getXtuple() + f->getNetwork()->velocityRigidBody(p2->getXtuple()) * p->RK_X.getdt();
               /* ------------------------------------------------------ */
               auto net = f->getNetwork();
               Quaternion q;
               q = q.d_dt(net->velocityRotational());
               //
               auto COM = net->RK_COM.getX(net->velocityTranslational());
               auto Q = net->RK_Q.getX(q());
               auto X0 = velocityRigidBody(p0, COM, Q);
               auto X1 = velocityRigidBody(p1, COM, Q);
               auto X2 = velocityRigidBody(p2, COM, Q);

               nextBodyVertices.emplace_back(T3Tddd{X0, X1, X2});
            }
            //% ------------------------------------------------------ */
            //%              vectorsToSurfaceFromBufferX              */
            //%              　　　近傍のノイマン面へ移動           　　　   */
            //% ------------------------------------------------------ */
            p->U_BUFFER_BUFFER = vectorsToSurfaceFromBufferX(p, nextBodyVertices);
            //! 角のディリクレ面へのめり込みを防止
            if (p->CORNER)
               p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, p->getNormalDirichlet_BEM()) * p->getNormalDirichlet_BEM();
         }
      }

      for (const auto &p : Points) {
         if (isFinite(p->U_BUFFER_BUFFER))
            p->U_BUFFER += p->U_BUFFER_BUFFER;
         p->U_BUFFER_BUFFER = {0., 0., 0.};
      }
   }

   for (const auto &p : Points) {
      if (!isFinite(p->U_BUFFER)) {
         std::cout << "p->RK_X.getdt() = " << p->RK_X.getdt() << std::endl;
         std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
         std::cout << "p->U_BUFFER = " << p->U_BUFFER << std::endl;
         std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
         std::cout << "p->Neumann = " << p->Neumann << std::endl;
         std::cout << "p->CORNER = " << p->CORNER << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      }
   }
};

#define derivatives_debug
struct derivatives {
   // uomap_P_d P_kappa;
   // std::unordered_set<networkPoint *> Points;
   // std::unordered_set<networkFace *> Faces;
   ~derivatives(){};
   derivatives(const Network &net, bool adjust_dirichlet = false) {
      auto &Points = net.getPoints();
      auto &Faces = net.getFaces();
      // P_laplacian.reserve(Points.size());
#ifdef derivatives_debug
      std::cout << Red << "initialize for parallelization" << colorOff << std::endl;
#endif
      //! initialize for parallelization
      /* ------------------------------------------------------ */
      // #ifdef _OPENMP
      //         std::cout << "曲率などの計算" << std::endl;
      // #pragma omp parallel
      // #endif
      //         for (const auto &p : Points)
      // #ifdef _OPENMP
      // #pragma omp single nowait
      // #endif
      //         {
      //             V_netPp ps = Flatten(BFS(p, 3));
      //             auto interpNormals = InterpolationVectorRBF(ToVector(extX(ps)), ToVector(extNormals(ps)), p->getX());
      //             p->kappa_BEM = interpNormals.div(p->getX()) / 2.; //中心方向法線ベクトルの場合，マイナスをつける．
      //             //! https://en.wikipedia.org/wiki/Mean_curvature
      //         }
      // #endif
      /* ------------------------------------------------------ */

      // b* ------------------------------------------------------ */
      // b*                     U_update_BEMの計算                  */
      // b* ------------------------------------------------------ */
      std::cout << "U_updateBEMを計算👇" << std::endl;
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
      {
         if (!isFinite(p->phiphin)) {
            std::cout << "p->phiphinはfiniteではない！！" << std::endl;
            std::cout << "p->phiphin = " << p->phiphin << std::endl;
            if (p->Neumann)
               std::cout << "p->Neumann" << std::endl;
            if (p->CORNER)
               std::cout << "p->CORNER" << std::endl;
            if (p->Dirichlet)
               std::cout << "p->CORNER" << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         }
         /* ------------------------------------------------------ */
         p->U_BEM = gradPhi(p);
         /* -------------------- おおよそのアップデート流速 ------------------- */
         //@ U_update_BEM は first guess
         // 2022/06/17
         // if (p->CORNER)
         // {
         //     Tddd tang = Normalize(Cross(p->getNormalDirichlet_BEM(), p->getNormalNeumann_BEM()));
         //     p->U_update_BEM = p->U_BEM - Dot(p->U_BEM, tang) * tang;
         // }
         // else

         if (p->Neumann) {
            p->U_update_BEM = uNeumann(p);
            if (!isFinite(p->U_update_BEM))
               p->U_update_BEM = {0., 0., 0.};
         } else
            p->U_update_BEM = p->U_BEM;
      }
      std::cout << "U_updateBEMを計算✅" << std::endl;
      /* ------------------------------------------------------ */

      //@ この後U_update_BEMをclingなどを使って修正する
      // ここの．BUFFER：Ω(t+δt)はルンゲクッタが見積もる時刻の表面と一致しているか？
      // for (const auto &p : Points)
      // 	p->X_BUFFER = p->getXtuple() + p->U_update_BEM * dt;

      // for (const auto &p : Points)
      // 	p->X_BUFFER = p->RK_X.getXinit() + p->U_update_BEM * p->RK_X.getdt();

      // 初期から考えた流速を与える必要があるのでは？
      // ルンゲクッタに従って次のΩ(t+δt)を予測する
      for (const auto &p : Points) {
         p->X_BUFFER = p->RK_X.getX(p->U_update_BEM);
         if (!isFinite(p->X_BUFFER)) {
            std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
            std::cout << "p->X_BUFFER = " << p->X_BUFFER << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         }
      }
      // #ifdef _OPENMP
      //         std::cout << "ラプラシアンを計算" << std::endl;
      // #pragma omp parallel
      // #endif
      //         for (const auto &p : Points)
      // #ifdef _OPENMP
      // #pragma omp single nowait
      // #endif
      //         {
      //             auto ps = Flatten(BFS(p, 2, {p->getNetwork()}));
      //             auto intp = InterpolationVectorRBF(obj3D::extractX(ps), extVelocities(ps), p->getX());
      //             auto tmp = intp.laplacian(p->getX());
      //             p->laplacian_U_BEM = {tmp[0], tmp[1], tmp[2]};
      //         }

      /* ------------------------------------------------------ */

      // double gamma = 72.75 * 1E-3;	  //[N/m] 水20度
      // double gravity = _GRAVITY_;		  //[m/s2]
      // double density = _WATER_DENSITY_; //[kg/m3]
      // double nu = 0.01005 / density;

      // for (auto &[p, v] : this->P_kappa)
      //     v = p->kappa_BEM;

      // for (auto &[p, v] : this->P_laplacian)
      //     v = p->laplacian_U_BEM;
      /*
          この方法は，ノイマン境界条件のclingにおいて，とても自然に無理なく応用できる．
          ディリクレ境界条件に関しては，計算後に補間によってリグリッドしてもいいかもしれない．
      */
      /* ------------------------------------------------------ */
      /*
          X_BUFFERには，U_update_BEMで単純に予測した節点が保存されている．
          予測したディリクレ面は正しいと考える．
          予測したノイマン面よりも，実際に物体を移動させて作った面の方が正しい．
          そこで，ノイマン面と角点に関しては，物体を移動させて作った面に貼り付ける．
      */
      for (const auto &p : Points)
         p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};

      calculateVectorToSurfaceInBuffer(net, adjust_dirichlet);
      std::cout << "calculateVectorToSurfaceInBuffer✅" << std::endl;
      calculateVectorFromBufferToContactFaces(net);
      std::cout << "calculateVectorFromBufferToContactFaces✅" << std::endl;
      // /* ------------------------------------------------------ */

      for (auto &p : Points) {
         p->U_update_BEM += p->U_BUFFER / p->RK_X.getdt();
         // dxdt_correct = p->U_BUFFER / p->RK_X.getdt();
         if (!isFinite(p->U_update_BEM, 1E+10) || !isFinite(p->U_BUFFER, 1E+10)) {
            std::cout << "p->RK_X.getdt() = " << p->RK_X.getdt() << std::endl;
            std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
            std::cout << "p->U_BUFFER = " << p->U_BUFFER << std::endl;
            std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
            std::cout << "p->Neumann = " << p->Neumann << std::endl;
            std::cout << "p->CORNER = " << p->CORNER << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         }
      }

#ifdef derivatives_debug
      std::cout << "φtとφntを一部計算👇" << std::endl;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
      {
         //% ------------------------------------------------------ */
         //%    ノイマン境界面上の加速度から,ノイマン境界面上のφntを計算     */
         //% ------------------------------------------------------ */
         if (p->Neumann || p->CORNER) {
            auto n = p->getNormalNeumann_BEM();
            T6d VW = velocity_from_Neumann_surface(p);
            Tddd angular_velocity = {std::get<3>(VW), std::get<4>(VW), std::get<5>(VW)};
            auto Q = Quaternion();
            auto dQdt = Q.d_dt(angular_velocity);
            for (const auto &[f, _] : p->phinOnFace) {
               if (f == nullptr) {
                  std::get<1>(p->phiphin_t) = accel_normal_from_Neumann_surface(p) - Dot(n, Dot(p->U_BEM, grad_U_LinearElement(p)));
                  std::get<1>(p->phiphin_t) += Dot(n, Dot(velocity_normal_from_Neumann_surface(p) - p->U_BEM, dQdt.Rv()));
                  p->phintOnFace[nullptr] = std::get<1>(p->phiphin_t);
                  /*
                  ∇U=∇∇f={{fxx, fyx, fzx},
                          {fxy, fyy, fzy},
                          {fxz, fyz, fzz}}
                  なので，∇∇f=∇∇f^T
                  */
               } else {
                  double tmp = accel_normal_from_Neumann_surface(p, f) - Dot(n, Dot(p->U_BEM, grad_U_LinearElement(p)));
                  tmp += Dot(n, Dot(Dot(uNeumann(p), f->normal) - p->U_BEM, dQdt.Rv()));
                  p->phintOnFace[f] = tmp;
               }
            }
         }
         //% ------------------------------------------------------ */
         //%                 ディリクレ境界面上のφtを計算                */
         //% ------------------------------------------------------ */
         if (p->Dirichlet || p->CORNER)
            std::get<0>(p->phiphin_t) = p->DphiDt(0.) - Dot(p->U_BEM, p->U_BEM);  //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！

         // this->P_dxdt[p] = p->U_update_BEM; //流速
         // //!この場合マイナスでないと，上の部分が半たんする
         // auto DphiDt = p->DphiDt(p->U_update_BEM, 0.);
         //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！
         // this->P_DphiDt[p] = DphiDt; // update用

         // this->P_aphiat[p] = std::get<0>(p->phiphin_t);  //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！
         // this->P_aphiant[p] = std::get<1>(p->phiphin_t); //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！

         // this->P_U_dot_gradgrad_U[p] = Dot(p->U_BEM, grad_U_LinearElement(p));
         // this->P_pressure[p] = p->pressure_BEM;

         // 10000. * (-1 / 2. * Vphi_Vphi - gravity * (std::get<2>(p->getXtuple())) - aphiat);
         // ここで圧力の項が抜けているが，これは全く流速に関係ないことに気づく．
         // なぜなら，表面上のどこでも同じだけ増加に寄与する大気圧は，
         // grad phiの計算によって，定数のため相殺されるからだ．
      }

#ifdef derivatives_debug
      std::cout << "φtとφntを一部計算✅" << std::endl;
#endif
   }
};

#endif