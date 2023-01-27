#ifndef BEM_derivatives_H
#define BEM_derivatives_H

#include "Network.hpp"

Tddd RK_with_Ubuff(const networkPoint *p) {
   return p->RK_X.getX(p->U_update_BEM + p->U_BUFFER / p->RK_X.getdt());
};

Tddd RK_with_Ubuff(const networkPoint *p, const Tddd &U_BUFFER) {
   return p->RK_X.getX(p->U_update_BEM + U_BUFFER / p->RK_X.getdt());
};

T3Tddd RK_with_Ubuff(const networkPoint *p0, const networkPoint *p1, const networkPoint *p2) {
   return {RK_with_Ubuff(p0), RK_with_Ubuff(p1), RK_with_Ubuff(p2)};
};

T3Tddd RK_with_Ubuff(const T_PPP &p012) {
   return RK_with_Ubuff(std::get<0>(p012), std::get<1>(p012), std::get<2>(p012));
};

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

Tddd RK_without_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->U_update_BEM); };

T3Tddd RK_without_Ubuff(const networkPoint *p0, const networkPoint *p1, const networkPoint *p2) {
   return {RK_without_Ubuff(p0), RK_without_Ubuff(p1), RK_without_Ubuff(p2)};
};

T3Tddd RK_without_Ubuff(const T_PPP &p012) {
   return RK_without_Ubuff(std::get<0>(p012), std::get<1>(p012), std::get<2>(p012));
};

T3Tddd RK_without_Ubuff(const networkFace *f) { return RK_without_Ubuff(f->getPoints()); };

T2Tddd RK_without_Ubuff(const networkPoint *p0, const networkPoint *p1) {
   return {RK_without_Ubuff(p0), RK_without_Ubuff(p1)};
};

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

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

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

void add_U_BUFFER_BUFFER_to_U_BUFFER(const auto &p) { p->U_BUFFER += p->U_BUFFER_BUFFER; };

std::vector<T3Tddd> nextBodyVertices(const std::unordered_set<networkFace *> &Fs) {
   std::vector<T3Tddd> ret;
   for (auto &f : Fs) {
      auto net = f->getNetwork();
      auto [p0, p1, p2] = f->getPoints();
      if (net->isFixed)
         ret.emplace_back(ToX(f));
      else if (net->isRigidBody) {
         Quaternion q;
         q = q.d_dt(net->velocityRotational());
         auto COM = net->RK_COM.getX(net->velocityTranslational());
         Quaternion Q(net->RK_Q.getX(q()));
         auto X0 = Q.Rv(p0->initialX - p0->getNetwork()->ICOM) + COM;
         auto X1 = Q.Rv(p1->initialX - p1->getNetwork()->ICOM) + COM;
         auto X2 = Q.Rv(p2->initialX - p2->getNetwork()->ICOM) + COM;
         ret.emplace_back(T3Tddd{X0, X1, X2});
      } else if (net->isSoftBody) {
         auto X0 = p0->RK_X.getX(p0->velocityTranslational());
         auto X1 = p1->RK_X.getX(p1->velocityTranslational());
         auto X2 = p2->RK_X.getX(p2->velocityTranslational());
         ret.emplace_back(T3Tddd{X0, X1, X2});
      }
   }
   return ret;
};

bool factor_angle(const networkPoint *p, const Tddd &ubuff) {
   for (const auto &f : p->getFaces()) {
      auto [p0, p1, p2] = f->getPoints(p);
      // if (!isFlat(RK_with_Ubuff_Normal(f),TriangleNormal(T3Tddd{RK_with_Ubuff(p0, ubuff), RK_with_Ubuff(p1), RK_with_Ubuff(p2)}),
      //             45. * M_PI / 180.))
      //    return false;
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

Tddd factor(const networkPoint *p, const Tddd &ubuff, const double max_ratio = 0.01) {
   double adjustment = Norm(RK_with_Ubuff(p) - RK_with_Ubuff(p, ubuff));
   double minLenNext = minNextLineLength(p);
   double ratio = adjustment / minLenNext;
   auto tmp = (ratio < max_ratio) ? ubuff : max_ratio * ubuff / ratio;
   if (factor_angle(p, tmp))
      return tmp;
   else
      return {0, 0, 0};
};

//% ------------------------------------------------------ */
//%                   vectorsToStructure                   */
//%              　　　近傍のノイマン面へ移動           　　　   */
//% ------------------------------------------------------ */

Tddd vectorsToStructure(const networkPoint *p, const std::vector<T3Tddd> &next_Vrtx) {
   auto toClosestXFacing = [&](const Tddd &n, const Tddd &p_next_X) {
      Tddd r = {1E+100, 1E+100, 1E+100}, X;
      for (const auto &vertex : next_Vrtx) {
         if (isInContact(p_next_X, n, vertex, p->radius)) {
            X = Nearest(p_next_X, vertex);
            if (Norm(r) >= Norm(X - p_next_X)) {
               r = X - p_next_X;
            }
         }
      }
      return r;
   };

   if (next_Vrtx.empty())
      return {0., 0., 0.};
   else if (p->Neumann || p->CORNER) {
      std::vector<Tddd> F_clings;
      std::vector<double> weights;
      for (const auto &f : p->getFacesNeumann()) {
         auto to_closest_X = toClosestXFacing(f->normal, RK_with_Ubuff(p));
         F_clings.push_back(Projection(to_closest_X, f->normal));
         weights.push_back(w_Bspline5(Norm(to_closest_X), p->radius));
      }
      // Tddd r = 0.1 * optimumVector_(F_clings, {0., 0., 0.}, weights);
      Tddd r = 0.5 * optimumVector_(F_clings, {0., 0., 0.});
      if (isFinite(r))
         return r;
      else
         return {0., 0., 0.};
   } else
      return {0., 0., 0.};
};

void calculateVectorFromBufferToContactFaces(const Network &net) {
   /*
   @ ノイマン面に貼り付けるための必要な調整
   */
   for (auto kk = 0; kk < 10; ++kk) {
#pragma omp parallel
      for (const auto &p : net.getPoints())
#pragma omp single nowait
      {
         p->U_BUFFER_BUFFER = {0., 0., 0.};
         if (p->CORNER || p->Neumann) {
            // 接触面候補の次の時刻の位置を予測
            std::unordered_set<networkFace *> Fs = p->getContactFaces();
            for (auto &f : BFS_Flattened(p->getContactFaces(), 2))
               Fs.emplace(f);
            //! 角のディリクレ面へのめり込みを防止
            p->U_BUFFER_BUFFER = vectorsToStructure(p, nextBodyVertices(Fs));
            if (p->CORNER)
               p->U_BUFFER_BUFFER = Chop(p->U_BUFFER_BUFFER, getNextNormalDirichlet_BEM(p));
         }
      }

      for (const auto &p : net.getPoints()) {
         add_U_BUFFER_BUFFER_to_U_BUFFER(p);
         p->U_BUFFER_BUFFER = {0., 0., 0.};
      }
   }
};

//% ------------------------------------------------------ */
//% ------------------------------------------------------ */
//% ------------------------------------------------------ */

/* -------------------------------------------------------------------------- */
/*                                 接線方向シフト                               */
/* -------------------------------------------------------------------------- */

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
         auto a = magicalValue(p, f);
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
   Tddd vectorToCenter = {0., 0., 0.};
   double s = 0;
   Tddd pX = RK_with_Ubuff(p);
   double min_distance = 1E+10;
   if (p->CORNER) {
      for (const auto &l : p->getLinesCORNER()) {
         vectorToCenter += RK_with_Ubuff((*l)(p)) - pX;
         s += 1.;
      }
      vectorToCenter /= s;
      return condition_Ua(scale * vectorToCenter, p);
   } else {
      for (const auto &f : p->getFaces()) {
         auto nP012 = RK_with_Ubuff(f->getPoints(p));
         auto [np0x, np1x, np2x] = nP012;
         auto mean = (Norm(np1x - np0x) + Norm(np2x - np1x) + Norm(np0x - np2x)) / 3.;
         auto l12 = Norm(np2x - np1x);
         double a = magicalValue(p, f);
         if (any_of(f->getLines(), [&](const auto &l) { return l->CORNER; }))
            a *= 27;
         else if (any_of(f->getPoints(), [&](const auto &p) { return p->CORNER; }))
            a *= 9;
         else if (any_of(f->getPoints(),
                         [&](const auto &p) { return std::any_of(p->getFaces().begin(), p->getFaces().end(),
                                                                 [&](const auto &F) { return any_of(F->getLines(), [&](const auto &l) { return l->CORNER; }); }); }))
            a *= 3;
         auto n = Normalize(Chop(np0x - np1x, np2x - np1x));
         auto X = Norm(np2x - np1x) * n * sin(M_PI / 3.) + (np2x + np1x) / 2.;
         vectorToCenter += a * (X - np0x);  //[m]
         s += a;
      }
      vectorToCenter /= s;
      return condition_Ua(scale * vectorToCenter, p);
   }
   return vectorToCenter;
};

/* -------------------------------------------------------------------------- */

Tddd vectorToNextSurface(const networkPoint *p) {
   /*do not scale this*/
   /*
   UartificialClingは，完璧にOmega(t+\delta t)に張り付くようにしなければ，
   面からはなれることで計算の破綻を招く可能性がある．
   */
   Tddd pX = RK_with_Ubuff(p);
   Tddd X, ret = {1E+20, 1E+20, 1E+20};
   if (p->CORNER) {
      for (const auto &l : p->getLines())
         if (l->CORNER)
            if (Norm(ret) >= Norm((X = Nearest(pX, RK_without_Ubuff(p, (*l)(p)))) - pX))
               ret = X - pX;
      //! 水槽の端のCORNERがノイマン面から離れるのを防ぐ
      for (const auto &f : p->getFacesNeumann())
         ret = Chop(ret, RK_without_Ubuff_Normal(f));
      // このベクトルは次の時刻の点の位置における修正流速なので，RK_without_Ubuff_Normal(f)を使うべき．
      return ret;
   } else {
      for (const auto &f : p->getFaces()) {
         if ((p->Dirichlet && f->Dirichlet) || (p->Neumann && f->Neumann))
            if (Norm(ret) >= Norm((X = Nearest(pX, RK_without_Ubuff(f))) - pX))
               ret = X - pX;
      }
      return ret;
   }
};

//$ -------------------------------------------------------------------------- */
void calculateVectorToSurfaceInBuffer(const Network &net, const int loop = 20, const bool adjust_dirichlet = true) {
   /*
   @ この方法なら，次の時刻における任意の場所でのポテンシャルを見積もることができる．
   @ このことは，任意のノイマン面上に節点を維持する上で便利である．
   @ Ω(t+δt)をまず見積もり，その面上で最適な格子配置となるように流速を修正する．
   */
   const double scale = 0.5;
   for (auto kk = 0; kk < loop; ++kk) {
      //$ ------------------------------------------------------ */
      //$           　　　　 vectorTangentialShift   　 　         */
      //$          ラプラス平滑化と引っ張り合わせた接線方向にシフト      */
      //$ ------------------------------------------------------ */

#pragma omp parallel
      for (const auto &p : net.getPoints()) {
#pragma omp single nowait
         if (p->CORNER)
            p->U_BUFFER_BUFFER = vectorTangentialShift2(p, scale);
         else
            p->U_BUFFER_BUFFER = factor(p, /*0.01 * vectorTangentialShift(p, scale) +*/ vectorTangentialShift2(p, scale));
      }

      for (const auto &p : net.getPoints()) {
         add_U_BUFFER_BUFFER_to_U_BUFFER(p);
         p->U_BUFFER_BUFFER = {0., 0., 0.};
      }

      //$ ------------------------------------------------------ */
      //$                   vectorToNextSurface                  */
      //$           　　　   周辺ディリクレ面に移動        　　　　    */
      //$ ------------------------------------------------------ */

#pragma omp parallel
      for (const auto &p : net.getPoints())
#pragma omp single nowait
         if (p->CORNER || p->Dirichlet)
            p->U_BUFFER_BUFFER = p->clungSurface = vectorToNextSurface(p);

      for (const auto &p : net.getPoints()) {
         add_U_BUFFER_BUFFER_to_U_BUFFER(p);
         p->U_BUFFER_BUFFER = {0., 0., 0.};
      }
      calculateVectorFromBufferToContactFaces(net);
   }
   calculateVectorFromBufferToContactFaces(net);
};

//$ -------------------------------------------------------------------------- */

Tddd gradPhi(const networkPoint *const p) {
   try {
      Tddd u;
      V_Tddd UW, V;
      V_d weights;
      double w;
      for (const auto &f : p->getFaces()) {
         auto [p0, p1, p2] = f->getPoints(p);
         u = gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)}, {ToX(p0), ToX(p1), ToX(p2)});
         if (f->Neumann) {
            if (p->phinOnFace.find(f) != p->phinOnFace.end()) {
               u += f->normal * p->phinOnFace.at(f);
               w = f->area;
            } else if (p->phinOnFace.find(nullptr) != p->phinOnFace.end()) {
               u += f->normal * p->phinOnFace.at(nullptr);
               w = f->area;
            } else
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, std::to_string(p->phinOnFace.size()));
         } else if (f->Dirichlet) {
            u += f->normal * p->phin_Dirichlet;
            w = 100 * f->area;  // よりDirichletに合わせるように重みを大きくした
         } else
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         V.emplace_back(u);
         UW.emplace_back(u * w);
         weights.emplace_back(w);
      }

      Tddd grad = Total(UW) / Total(weights);

      // //%! --------------------------------------------------------- */
      if (p->CORNER || p->Neumann)
         if (!V.empty()) {
            auto ret = optimumVector_(V, {0., 0., 0.}, weights);
            if (isFinite(ret))
               grad = ret;
         }
      // //%! --------------------------------------------------------- */
      return grad;

   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

#define derivatives_debug
struct derivatives {
   ~derivatives(){};
   derivatives(const Network &net, bool adjust_dirichlet = true) {
      auto &Points = net.getPoints();
      auto &Faces = net.getFaces();
#ifdef derivatives_debug
      std::cout << Red << "initialize for parallelization" << colorOff << std::endl;
#endif
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
         p->U_update_BEM = p->U_BEM = gradPhi(p);
         /* -------------------- おおよそのアップデート流速 ------------------- */
         //@ U_update_BEM は first guess
         // 2022/06/17
         // if (p->CORNER) {
         //    p->U_update_BEM = p->U_BEM;  // Chop(p->U_BEM, Normalize(Cross(p->getNormalDirichlet_BEM(), p->getNormalNeumann_BEM())));
         // } else if (p->Neumann)
         //    p->U_update_BEM = Chop(uNeumann(p), p->getNormal_BEM());
         // else
         //    p->U_update_BEM = p->U_BEM;
      }
      std::cout << "U_updateBEMを計算✅" << std::endl;
      /* ------------------------------------------------------ */
      //@ この後U_update_BEMをclingなどを使って修正する
      //@ ルンゲクッタに従って次のΩ(t+δt)を予測する
      // for (const auto &p : Points) {
      //    p->X_BUFFER = p->RK_X.getX(p->U_update_BEM);
      //    if (!isFinite(p->X_BUFFER)) {
      //       std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
      //       std::cout << "p->X_BUFFER = " << p->X_BUFFER << std::endl;
      //       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      //    }
      // }
      /* ------------------------------------------------------ */
      /*
          予測したディリクレ面は正しいと考える．
          予測したノイマン面よりも，実際に物体を移動させて作った面の方が正しい．
          そこで，ノイマン面と角点に関しては，物体を移動させて作った面に貼り付ける．
      */
      for (const auto &p : Points)
         p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};

      calculateVectorToSurfaceInBuffer(net, 20, adjust_dirichlet);
      std::cout << "calculateVectorToSurfaceInBuffer✅" << std::endl;
      // calculateVectorFromBufferToContactFaces(net);
      // std::cout << "calculateVectorFromBufferToContactFaces✅" << std::endl;
      // /* ------------------------------------------------------ */

      for (auto &p : Points) {
         p->U_update_BEM += p->U_BUFFER / p->RK_X.getdt();
         // dxdt_correct = p->U_BUFFER / p->RK_X.getdt();
         if (!isFinite(p->U_update_BEM, 1E+10) || !isFinite(p->U_BUFFER, 1E+10)) {
            std::cout << "p->X = " << p->X << std::endl;
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
            /*
            ∇U=∇∇f={{fxx, fyx, fzx},
                     {fxy, fyy, fzy},
                     {fxz, fyz, fzz}}
            なので，∇∇f=∇∇f^T
            */
            // b* p->phintOnFaceは，std::unordered_map<networkFace *, double>
            // b* 節点のphinを保存する．また，多重節点かどうかも，面がnullptrかどうかで判別できる．
            // b* setBoundaryConditionsで決めている．
            auto n = p->getNormalNeumann_BEM();
            auto Q = Quaternion();
            for (auto &[f, phin_t] : p->phintOnFace) {
               if (f) {
                  auto dQdt = Q.d_dt(NearestContactFace(f)->getNetwork()->velocityRotational());
                  auto n = f->normal;
                  auto [p0, p1, p2] = f->getPoints(p);
                  Tddd phi012 = {std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)};
                  Tddd phin012 = {std::get<1>(p0->phiphin), std::get<1>(p1->phiphin), std::get<1>(p2->phiphin)};
                  Tddd grad_phi = Mean(phin012) * n + gradTangential_LinearElement(phi012, ToX(f));
                  phin_t = std::get<1>(p->phiphin_t) = Dot(n, Dot(uNeumann(p, f) - grad_phi, dQdt.Rv()) + accelNeumann(p, f) - Dot(grad_phi, grad_U_LinearElement(f)));
               } else {
                  auto dQdt = Q.d_dt(NearestContactFace(p)->getNetwork()->velocityRotational());
                  phin_t = std::get<1>(p->phiphin_t) = Dot(n, Dot(uNeumann(p) - p->U_BEM, dQdt.Rv()) + accelNeumann(p) - Dot(p->U_BEM, grad_U_LinearElement(p)));
               }
            }
         }
         //% ------------------------------------------------------ */
         //%                 ディリクレ境界面上のφtを計算                */
         //% ------------------------------------------------------ */
         if (p->Dirichlet || p->CORNER)
            std::get<0>(p->phiphin_t) = p->DphiDt(0.) - Dot(p->U_BEM, p->U_BEM);  //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！
      }

#ifdef derivatives_debug
      std::cout << "φtとφntを一部計算✅" << std::endl;
#endif
   }
};

#endif