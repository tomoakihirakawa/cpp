#ifndef BEM_calculateVelocities_H
#define BEM_calculateVelocities_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

//$ -------------------------------------------------------------------------- */
//$                         calculateVecToSurface                              */
//$ -------------------------------------------------------------------------- */

Tddd RK_with_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->U_BEM + p->vecToSurface / p->RK_X.getdt()); };
Tddd RK_with_Ubuff(const networkPoint *p, const Tddd &vecToSurface) { return p->RK_X.getX(p->U_BEM + vecToSurface / p->RK_X.getdt()); };
T3Tddd RK_with_Ubuff(const networkPoint *p0, const networkPoint *p1, const networkPoint *p2) { return {RK_with_Ubuff(p0), RK_with_Ubuff(p1), RK_with_Ubuff(p2)}; };
T3Tddd RK_with_Ubuff(const T_PPP &p012) { return {RK_with_Ubuff(p012[0]), RK_with_Ubuff(p012[1]), RK_with_Ubuff(p012[2])}; };
T3Tddd RK_with_Ubuff(const networkFace *f) { return RK_with_Ubuff(f->getPoints()); };
Tddd RK_with_Ubuff_Normal(const networkFace *f) { return TriangleNormal(RK_with_Ubuff(f->getPoints())); };
double RK_with_Ubuff_Area(const networkFace *f) { return TriangleArea(RK_with_Ubuff(f->getPoints())); };
Tddd RK_with_Ubuff_Normal(const networkPoint *p) {
   Tddd normal = {0., 0., 0.};
   double a = 0, total = 0;
   for (const auto &f : p->getFaces()) {
      a = RK_with_Ubuff_Area(f);
      normal += a * RK_with_Ubuff_Normal(f);
      total += a;
   }
   return Normalize(normal / total);
};

Tddd RK_without_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->U_BEM); };
T3Tddd RK_without_Ubuff(const networkPoint *p0, const networkPoint *p1, const networkPoint *p2) { return {RK_without_Ubuff(p0), RK_without_Ubuff(p1), RK_without_Ubuff(p2)}; };
T3Tddd RK_without_Ubuff(const T_PPP &p012) { return RK_without_Ubuff(std::get<0>(p012), std::get<1>(p012), std::get<2>(p012)); };
T3Tddd RK_without_Ubuff(const networkFace *f) { return RK_without_Ubuff(f->getPoints()); };
T2Tddd RK_without_Ubuff(const networkPoint *p0, const networkPoint *p1) { return {RK_without_Ubuff(p0), RK_without_Ubuff(p1)}; };
Tddd RK_without_Ubuff_Normal(const networkFace *f) { return TriangleNormal(RK_without_Ubuff(f->getPoints())); };
Tddd RK_without_Ubuff_Normal(const networkPoint *p) {
   Tddd normal = {0., 0., 0.};
   double a = 0, total = 0;
   for (const auto &f : p->getFaces()) {
      a = TriangleArea(RK_without_Ubuff(f));
      normal += a * RK_without_Ubuff_Normal(f);
      total += a;
   }
   return Normalize(normal / total);
};
Tddd getNextNormalDirichlet_BEM(const networkPoint *p) {
   Tddd normal = {0., 0., 0.};
   double a = 0, total = 0;
   for (const auto &f : p->getFacesDirichlet()) {
      a = TriangleArea(RK_without_Ubuff(f));
      normal += a * RK_without_Ubuff_Normal(f);
      total += a;
   }
   return Normalize(normal / total);
};
Tddd getNextNormalNeumann_BEM(const networkPoint *p) {
   Tddd normal = {0., 0., 0.};
   double a = 0, total = 0;
   for (const auto &f : p->getFacesNeumann()) {
      a = TriangleArea(RK_without_Ubuff(f));
      normal += a * RK_without_Ubuff_Normal(f);
      total += a;
   }
   return Normalize(normal / total);
};

Tddd condition_Ua(Tddd VECTOR, const networkPoint *const p) {
   /*
   考え方：修正流速は，次の時刻における修正量なので，
   chopする法線方向なども次の時刻における法線方向でないといけない：RK_with_Ubuff_Normal
   */
   if (p->Dirichlet) {
      return Chop(VECTOR, RK_without_Ubuff_Normal(p));
   } else {
      if (p->CORNER)
         VECTOR = Projection(VECTOR, Cross(getNextNormalNeumann_BEM(p), getNextNormalDirichlet_BEM(p)));
      for (const auto &f : p->getFacesNeumann()) {
         VECTOR = Chop(VECTOR, RK_without_Ubuff_Normal(f));
         // VECTOR = Chop(VECTOR, f->normal);  // f->normalでないといけないのか？ 関係なかった
      }
      return VECTOR;
   }
};

void add_vecToSurface_BUFFER_to_vecToSurface(const auto &p) {
   p->vecToSurface += p->vecToSurface_BUFFER;
   if (!isFinite(p->vecToSurface))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not finite");
};

T3Tddd nextBodyVertex(const networkFace *f) {
   auto net = f->getNetwork();
   auto [p0, p1, p2] = f->getPoints();
   if (net->isRigidBody) {
      Quaternion q;
      q = q.d_dt(net->velocityRotational());
      auto COM = net->RK_COM.getX(net->velocityTranslational());
      Quaternion Q(net->RK_Q.getX(q()));
      auto X0 = Q.Rv(p0->initialX - p0->getNetwork()->ICOM) + COM;
      auto X1 = Q.Rv(p1->initialX - p1->getNetwork()->ICOM) + COM;
      auto X2 = Q.Rv(p2->initialX - p2->getNetwork()->ICOM) + COM;
      return T3Tddd{X0, X1, X2};
   } else if (net->isSoftBody) {
      auto X0 = p0->RK_X.getX(p0->velocityTranslational());
      auto X1 = p1->RK_X.getX(p1->velocityTranslational());
      auto X2 = p2->RK_X.getX(p2->velocityTranslational());
      return T3Tddd{X0, X1, X2};
   } else
      return ToX(f);
};

std::vector<T3Tddd> nextBodyVertices(const std::unordered_set<networkFace *> &Fs) {
   std::vector<T3Tddd> ret(Fs.size());
   int i = 0;
   for (auto &f : Fs) {
      auto net = f->getNetwork();
      auto [p0, p1, p2] = f->getPoints();
      if (net->isFixed)
         ret[i] = ToX(f);
      else if (net->isRigidBody) {
         Quaternion q;
         q = q.d_dt(net->velocityRotational());
         auto COM = net->RK_COM.getX(net->velocityTranslational());
         Quaternion Q(net->RK_Q.getX(q()));
         auto X0 = Q.Rv(p0->initialX - p0->getNetwork()->ICOM) + COM;
         auto X1 = Q.Rv(p1->initialX - p1->getNetwork()->ICOM) + COM;
         auto X2 = Q.Rv(p2->initialX - p2->getNetwork()->ICOM) + COM;
         ret[i] = {X0, X1, X2};
      } else if (net->isSoftBody) {
         auto X0 = p0->RK_X.getX(p0->velocityTranslational());
         auto X1 = p1->RK_X.getX(p1->velocityTranslational());
         auto X2 = p2->RK_X.getX(p2->velocityTranslational());
         ret[i] = {X0, X1, X2};
      }
      i++;
   }
   return ret;
};

bool factor_angle(const networkPoint *p, const Tddd &ubuff) {
   for (const auto &f : p->getFaces()) {
      auto [p0, p1, p2] = f->getPoints(p);
      if (0 > Dot(RK_with_Ubuff_Normal(f),
                  TriangleNormal(T3Tddd{RK_with_Ubuff(p0, ubuff), RK_with_Ubuff(p1), RK_with_Ubuff(p2)})))
         return false;
      if (0 > Dot(RK_without_Ubuff_Normal(f),
                  TriangleNormal(T3Tddd{RK_with_Ubuff(p0, ubuff), RK_with_Ubuff(p1), RK_with_Ubuff(p2)})))
         return false;
   }
   return true;
};

double minNextLineLength(const networkPoint *p) {
   double len = 1E+20;
   for (const auto &l : p->getLines()) {
      auto [p0, p1] = l->getPoints();
      double norm = Norm(RK_with_Ubuff(p0) - RK_with_Ubuff(p1));
      if (len > norm)
         len = norm;
   }
   return len;
};

double meanNextLineLength(const networkPoint *p) {
   double len = 0;
   int count = 0;
   for (const auto &l : p->getLines()) {
      auto [p0, p1] = l->getPoints();
      double norm = Norm(RK_with_Ubuff(p0) - RK_with_Ubuff(p1));
      len += norm;
      count++;
   }
   return len / count;
};

Tddd factor(const networkPoint *p, const Tddd &ubuff, const double max_ratio = 0.01) {
   double adjustment = Norm(RK_with_Ubuff(p) - RK_with_Ubuff(p, ubuff));
   double minLenNext = meanNextLineLength(p);
   double ratio = adjustment / minLenNext;
   auto tmp = (ratio < max_ratio) ? ubuff : max_ratio * ubuff / ratio;
   if (factor_angle(p, tmp))
      return tmp;
   else
      return {0, 0, 0};
};

double magicalValue(const networkPoint *p, const networkFace *f) {
   auto [p0, p1, p2] = f->getPoints(p);
   auto nP012 = T3Tddd{RK_with_Ubuff(p0, p1, p2)};
   auto tmp = std::log10(Circumradius(nP012) / Inradius(nP012));
   double max = 10;
   if (tmp > max)
      return max;
   else
      return tmp;
};

double variance2(const networkPoint *p, const networkFace *f) {
   auto [p0, p1, p2] = f->getPoints(p);
   double m = 0, l0, l1, l2;
   m += (l0 = Norm(RK_with_Ubuff(p0) - RK_with_Ubuff(p1)));
   m += (l1 = Norm(RK_with_Ubuff(p1) - RK_with_Ubuff(p2)));
   m += (l2 = Norm(RK_with_Ubuff(p2) - RK_with_Ubuff(p0)));
   m /= 3.;
   auto ret = std::pow(l0 - m, 2) + std::pow(l1 - m, 2) + std::pow(l2 - m, 2);
   return 10 * ret / (m * m);
};

Tddd vectorTangentialShift(const networkPoint *p, const double scale = 1.) {
   Tddd vectorToCenter = {0., 0., 0.};
   double s = 0;
   Tddd pX = RK_with_Ubuff(p);
   if (p->CORNER) {
      for (const auto &l : p->getLinesCORNER()) {
         vectorToCenter += RK_with_Ubuff((*l)(p)) - pX;
         s += 1.;
      }
      vectorToCenter /= s;
      return condition_Ua(scale * vectorToCenter, p);
   } else {
      double min_distance = 1E+10;

      for (const auto &f : p->getFaces()) {
         auto [p0, p1, p2] = f->getPoints(p);
         auto nP012 = T3Tddd{RK_with_Ubuff(p0, p1, p2)};
         auto a = magicalValue(p, f) + variance2(p, f);
         if (min_distance > Norm(Incenter(nP012) - RK_with_Ubuff(p)))
            min_distance = Norm(Incenter(nP012) - RK_with_Ubuff(p));

         auto toIC = Incenter(nP012) - pX;
         auto V = (Norm(toIC) - Inradius(nP012)) * Normalize(toIC);
         auto W = a;  // * TriangleArea(nP012) * Norm(V);
         vectorToCenter += V * W;
         s += W;
      }
      vectorToCenter /= s;  // dimension : [m]
      return condition_Ua(scale * vectorToCenter, p);
   }
};

Tddd vectorTangentialShift2(const networkPoint *p, const double scale = 1.) {
   Tddd vector_to_optimum_X = {0., 0., 0.}, pX = RK_with_Ubuff(p);
   double s = 0;
   if (p->CORNER) {
      for (const auto &l : p->getLinesCORNER()) {
         vector_to_optimum_X += RK_with_Ubuff((*l)(p)) - pX;
         s += 1.;
      }
   } else {
      for (const auto &f : p->getFaces()) {
         auto nP012 = RK_with_Ubuff(f->getPoints(p));
         auto &[np0x, np1x, np2x] = nP012;
         double a = magicalValue(p, f) + variance2(p, f);

         if (std::any_of(f->getLines().begin(), f->getLines().end(), [&](const auto &l) { return l->CORNER; })) {
            a *= std::pow(2., 3);
         } else if (std::any_of(f->getPoints().begin(), f->getPoints().end(), [&](const auto &p) { return p->CORNER; })) {
            a *= std::pow(2., 2);
         } else if (std::any_of(f->getPoints().begin(), f->getPoints().end(), [&](const auto &p) {
                       return std::any_of(p->getFaces().begin(), p->getFaces().end(), [&](const auto &F) {
                          return std::any_of(F->getPoints().begin(), F->getPoints().end(), [&](const auto &q) { return q->CORNER; });
                       });
                    })) {
            a *= 2.;
         }

         auto optimum_position = Norm(np2x - np1x) * Normalize(Chop(np0x - np1x, np2x - np1x)) * sin(M_PI / 3.) + (np2x + np1x) / 2.;
         vector_to_optimum_X += a * (optimum_position - pX);
         s += a;
      }
   }

   vector_to_optimum_X /= s;
   return condition_Ua(scale * vector_to_optimum_X, p);
};

Tddd vectorToNextSurface(const networkPoint *p) {
   /*
   UartificialClingは，完璧にOmega(t+\delta t)に張り付くようにしなければ，
   面からはなれることで計算の破綻を招く可能性がある．
   */
   Tddd pX = RK_with_Ubuff(p), X, ret = {1E+20, 1E+20, 1E+20};
   if (p->Dirichlet) {
      for (const auto &f : p->getFaces()) {
         if ((p->Dirichlet && f->Dirichlet) || (p->Neumann && f->Neumann)) {
            auto actual_corner_face = RK_without_Ubuff(f);
            if (Norm(ret) >= Norm((X = Nearest(pX, actual_corner_face)) - pX))
               ret = X - pX;
         }
      }
      return ret;
   } else if (p->Neumann || p->CORNER) {
      // to_cornerを計算
      Tddd to_corner = {0., 0., 0.};
      if (p->CORNER) {
         std::vector<Tddd> F_clings;
         std::vector<double> distance;
         Tddd to_corner_tmp;
         to_corner_tmp.fill(1E+20);
         for (const auto &l : p->getLines())
            if (l->CORNER) {
               auto actual_corner_line = RK_without_Ubuff(p, (*l)(p));
               auto X = Nearest(pX, actual_corner_line);
               if (Norm(to_corner_tmp) >= Norm(X - pX))
                  to_corner_tmp = X - pX;
            }
         to_corner = 0.8 * to_corner_tmp;
      }

      // to_structure_faceを計算
      auto next_Vrtx = nextBodyVertices(bfs(p->getContactFaces(), 2));
      auto p_next_X = to_corner + RK_with_Ubuff(p);
      Tddd to_structure_face;
      if (!next_Vrtx.empty()) {
         std::vector<Tddd> F_clings;
         F_clings.reserve(5 * next_Vrtx.size());
         std::vector<double> distance, weights;
         Tddd to_closest_X, X;
         for (const auto &f : p->getFacesNeumann() /*pの隣接面*/) {
            // pの周りの面が干渉している構造物面のみを対象として，点に最も近い構造物面上の点を抽出する
            bool found = false;
            to_closest_X.fill(1E+20);
            for (const auto &vertex : next_Vrtx)
               if (isInContact(p_next_X, f->normal, vertex, p->radius)) {
                  X = Nearest(p_next_X, vertex);
                  if (Norm(to_closest_X) >= Norm(X - p_next_X)) {
                     to_closest_X = X - p_next_X;
                     found = true;
                  }
               }
            if (found) {
               F_clings.push_back(Projection(to_closest_X, f->normal));  // この面fにとって，最も近くにある干渉点
               distance.push_back(Norm(to_closest_X));
            }
         }

         auto min = Min(distance);
         for (const auto &d : distance)
            weights.push_back(w_Bspline5(d, p->radius));
         to_structure_face = !F_clings.empty() ? optimumVector(F_clings, {0., 0., 0.}, 1E-12) : Tddd{0., 0., 0.};
      }
      return to_structure_face + to_corner;
   }
   return {0., 0., 0.};
};

void calculateVecToSurface(const Network &net, const int loop = 10) {
   /*
   @ この方法なら，次の時刻における任意の場所でのポテンシャルを見積もることができる．
   @ このことは，任意のノイマン面上に節点を維持する上で便利である．
   @ Ω(t+δt)をまず見積もり，その面上で最適な格子配置となるように流速を修正する．
   */
   for (const auto &p : net.getPoints()) {
      p->vecToSurface_BUFFER.fill(0.);
      p->vecToSurface.fill(0.);
   }
   TimeWatch watch;
   const double scale = 0.1;
   for (auto kk = 0; kk < loop; ++kk) {
      //$ ------------------------------------------------------ */
      //$           　　　　 vectorTangentialShift   　 　         */
      //$          ラプラス平滑化と引っ張り合わせた接線方向にシフト      */
      //$ ------------------------------------------------------ */
      for (auto ii = 0; ii < 30; ++ii) {
         // この計算コストは，比較的やすいので，何度も繰り返しても問題ない．
#pragma omp parallel
         for (const auto &p : net.getPoints())
#pragma omp single nowait
            p->vecToSurface_BUFFER = (p->Neumann ? 0.1 : 1.) * vectorTangentialShift2(p, scale);

         for (const auto &p : net.getPoints()) {
            add_vecToSurface_BUFFER_to_vecToSurface(p);
            p->vecToSurface_BUFFER.fill(0.);
         }
      }
      std::cout << "Elapsed time for 1.vectorTangentialShift : " << watch() << " [s]" << std::endl;
      //$ ------------------------------------------------------ */
      //$              　　　vectorToNextSurface                  */
      //$           　　　   周辺ディリクレ面に移動        　　　　    */
      //$ ------------------------------------------------------ */
#pragma omp parallel
      for (const auto &p : net.getPoints())
#pragma omp single nowait
         p->vecToSurface_BUFFER = p->clungSurface = vectorToNextSurface(p);

      for (const auto &p : net.getPoints()) {
         add_vecToSurface_BUFFER_to_vecToSurface(p);
         p->vecToSurface_BUFFER.fill(0.);
      }
      std::cout << "Elapsed time for 2.vectorToNextSurface: " << watch() << " [s]" << std::endl;
   }
};

//$ -------------------------------------------------------------------------- */
// b! -------------------------------------------------------------------------- */
// b!                               calculateVelocities                          */
// b! -------------------------------------------------------------------------- */
Tddd gradPhi(const networkPoint *const p) {
   try {
      Tddd u, ut, un;
      V_Tddd UW, V, Un, Ut, U;
      V_d weights;
      double w;
      for (const auto &f : p->getFaces()) {
         auto [p0, p1, p2] = f->getPoints(p);
         ut = gradTangential_LinearElement(Tddd{{std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)}}, T3Tddd{{ToX(p0), ToX(p1), ToX(p2)}});
         un.fill(0.);
         if (f->Neumann) {
            if (p->phinOnFace.find(f) != p->phinOnFace.end()) {
               un += f->normal * p->phinOnFace.at(f);
               w = f->area;
            } else if (p->phinOnFace.find(nullptr) != p->phinOnFace.end()) {
               un += f->normal * p->phinOnFace.at(nullptr);
               w = f->area;
            } else
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, std::to_string(p->phinOnFace.size()));
         } else if (f->Dirichlet) {
            un += f->normal * std::get<1>(p->phiphin);
            w = 100 * f->area;  // よりDirichletに合わせるように重みを大きくした
         } else
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         u = ut + un;
         Un.emplace_back(un);
         Ut.emplace_back(ut);
         U.emplace_back(u);
         UW.emplace_back(u * w);
         weights.emplace_back(w);

         if (!isFinite(u, 1E+10)) {
            std::cout << "p->X = " << p->X << std::endl;
            std::cout << "p->U_BEM = " << p->U_BEM << std::endl;
            for (const auto &q : f->getPoints(p))
               std::cout << "q->phiphin = " << q->phiphin << std::endl;
            for (const auto &f : p->getFaces())
               std::cout << "f->area = " << f->area << std::endl;
            std::cout << "un = " << un << std::endl;
            std::cout << "ut = " << ut << std::endl;
            std::cout << "u = " << u << std::endl;
            std::cout << "w = " << w << std::endl;
            std::cout << V_i{f->Dirichlet, f->Neumann} << std::endl;
            std::cout << V_i{p->Dirichlet, p->CORNER, p->Neumann} << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         }
      }

      Tddd grad = Total(UW) / Total(weights);

      //%! --------------------------------------------------------- */
      if (p->CORNER || p->Neumann)
         if (!V.empty()) {
            auto ret = optimumVector(V, grad, weights);
            if (isFinite(ret))
               grad = ret;
         }
      //%! --------------------------------------------------------- */

      return grad;

   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

void calculateVelocities(const Network &net, const int loop = 10) {
#pragma omp parallel
   for (const auto &p : net.getPoints())
#pragma omp single nowait
   {
      p->U_BEM = gradPhi(p);
   }

   std::cout << "p->U_BEM = gradPhi(p)を計算✅" << std::endl;
   //@ この後, U_update_BEMをclingなどを使って修正する
   //@ ルンゲクッタに従って次のΩ(t+δt)を予測する
   calculateVecToSurface(net, loop);
   std::cout << "calculateVecToSurface✅" << std::endl;

   for (const auto &p : net.getPoints()) {
      // dxdt_correct = p->vecToSurface / p->RK_X.getdt();
      p->U_update_BEM = p->U_BEM + p->vecToSurface / p->RK_X.getdt();
      if (!isFinite(p->U_update_BEM, 1E+10) || !isFinite(p->vecToSurface, 1E+10)) {
         std::cout << "p->X = " << p->X << std::endl;
         std::cout << "p->RK_X.getdt() = " << p->RK_X.getdt() << std::endl;
         std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
         std::cout << "p->vecToSurface = " << p->vecToSurface << std::endl;
         std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
         std::cout << "p->Neumann = " << p->Neumann << std::endl;
         std::cout << "p->CORNER = " << p->CORNER << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      }
   }
}

#endif