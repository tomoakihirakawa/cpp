#ifndef BEM_calculateVelocities_H
#define BEM_calculateVelocities_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

//$ -------------------------------------------------------------------------- */
//$                         calculateVecToSurface                              */
//$ -------------------------------------------------------------------------- */

Tddd RK_without_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->U_update_BEM); };
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
      FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
      total += a;
   }
   return Normalize(normal / total);
};
Tddd getNextNormalDirichlet_BEM(const networkPoint *p) {
   Tddd normal = {0., 0., 0.};
   double a = 0, total = 0;
   for (const auto &f : p->getFacesDirichlet()) {
      a = TriangleArea(RK_without_Ubuff(f));
      // normal += a * RK_without_Ubuff_Normal(f);
      FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
      total += a;
   }
   return Normalize(normal / total);
};
Tddd getNextNormalNeumann_BEM(const networkPoint *p) {
   Tddd normal = {0., 0., 0.};
   double a = 0, total = 0;
   for (const auto &f : p->getFacesNeumann()) {
      a = TriangleArea(RK_without_Ubuff(f));
      // normal += a * RK_without_Ubuff_Normal(f);
      FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
      total += a;
   }
   return Normalize(normal / total);
};

// Tddd RK_with_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->U_update_BEM + p->vecToSurface / p->RK_X.getdt()); };
Tddd RK_with_Ubuff(const networkPoint *p) {
   // return p->RK_X.toReachAtNextTimeQ(p->RK_X.getX(p->U_update_BEM) + p->vecToSurface);
   const auto U = p->RK_X.getVectorToReachAtNextTimeQ(RK_without_Ubuff(p) + p->vecToSurface);
   return p->RK_X.getX(U);
};
Tddd RK_with_Ubuff(const networkPoint *p, const Tddd &vecToSurface) {
   // return p->RK_X.getX(p->U_update_BEM + vecToSurface / p->RK_X.getdt());
   const auto U = p->RK_X.getVectorToReachAtNextTimeQ(RK_without_Ubuff(p) + vecToSurface);
   return p->RK_X.getX(U);
};
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
      // normal += a * RK_with_Ubuff_Normal(f);
      FusedMultiplyIncrement(a, RK_with_Ubuff_Normal(f), normal);
      total += a;
   }
   return Normalize(normal / total);
};

Tddd condition_Ua(Tddd VECTOR, const networkPoint *const p) {
   /*
   考え方：修正流速は，次の時刻における修正量なので，
   chopする法線方向なども次の時刻における法線方向でないといけない：RK_with_Ubuff_Normal
   */

   // if (p->Dirichlet) {
   //    return Chop(VECTOR, RK_without_Ubuff_Normal(p));
   // } else {
   // if (p->CORNER) {
   //    auto next_normal = getNextNormalDirichlet_BEM(p);
   //    for (const auto &f : p->getFacesNeumann())
   //       VECTOR = Projection(VECTOR, Cross(next_normal, RK_without_Ubuff_Normal(f)));
   // }
   // for (const auto &f : p->getFacesNeumann()) {
   //    VECTOR = Chop(VECTOR, RK_without_Ubuff_Normal(f));
   //    // VECTOR = Chop(VECTOR, f->normal);  // f->normalでないといけないのか？ 関係なかった
   // }
   for (const auto &[_, FX] : p->getNearestContactFaces()) {
      auto [F, X] = FX;
      if (F)
         VECTOR = Chop(VECTOR, F->normal);
      // VECTOR = Chop(VECTOR, f->normal);  // f->normalでないといけないのか？ 関係なかった
   }

   return VECTOR;
   // }
};

void add_vecToSurface_BUFFER_to_vecToSurface(const auto &p, const double a = 1.) {
   if (isFinite(p->vecToSurface_BUFFER))
      p->vecToSurface += p->vecToSurface_BUFFER * a;
};

std::array<double, 3> nextPositionOnBody(Network *net, networkPoint *p) {
   auto velocity = net->velocityTranslational();
   auto rotation = net->velocityRotational();
   if (net->isRigidBody) {
      for (auto i = 0; i < 3; ++i)
         if (net->isFixed[i])
            velocity[i] = 0;
      for (auto i = 3; i < 6; ++i)
         if (net->isFixed[i])
            rotation[i - 3] = 0;
   }
   auto next_COM = net->RK_COM.getX(velocity);
   auto next_Q = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(rotation, net->Q)));
   return rigidTransformation(net->ICOM, next_COM, next_Q.Rv(), p->initialX);
   // else return p->X;
}

std::vector<T3Tddd> nextBodyVertices(const std::unordered_set<networkFace *> &Fs) {
   std::vector<T3Tddd> ret(Fs.size());
   int i = 0;
   for (auto &f : Fs) {
      auto net = f->getNetwork();
      auto [p0, p1, p2] = f->getPoints();
      // if (net->isFixed)
      //    ret[i] = ToX(f);
      // else
      if (net->isRigidBody) {
         // 現在のv^nを使って問題ない．加速度はいらない．
         // Quaternion q;
         // q = q.d_dt(net->velocityRotational());
         // Quaternion next_Q(net->RK_Q.getX(net->Q.AngularVelocityTodQdt(net->velocityRotational())));
         // auto next_COM = net->RK_COM.getX(net->velocityTranslational());

         auto velocity = net->velocityTranslational();
         auto rotation = net->velocityRotational();

         for (auto i = 0; i < 3; ++i)
            if (net->isFixed[i])
               velocity[i] = 0;
         for (auto i = 3; i < 6; ++i)
            if (net->isFixed[i])
               rotation[i - 3] = 0;

         auto next_COM = net->RK_COM.getX(velocity);
         auto next_Q = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(rotation, net->Q)));
         auto X_next = [&](const auto &p) { return rigidTransformation(net->ICOM, next_COM, next_Q.Rv(), p->initialX); };

         ret[i] = {X_next(p0), X_next(p1), X_next(p2)};
      } else if (net->isSoftBody) {
         auto v0 = p0->velocityTranslational();
         auto v1 = p1->velocityTranslational();
         auto v2 = p2->velocityTranslational();
         if (std::ranges::all_of(net->isFixed, [](bool isfixed) { return isfixed; })) {
            v0.fill(0.);
            v1.fill(0.);
            v2.fill(0.);
         }
         auto X0 = p0->RK_X.getX(v0);
         auto X1 = p1->RK_X.getX(v1);
         auto X2 = p2->RK_X.getX(v2);
         ret[i] = {X0, X1, X2};
      }
      i++;
   }
   return ret;
};

// \label{BEM:vectorTangentialShift}
// Tddd vectorTangentialShift(const networkPoint *p, double scale = 1.) {
//    Tddd V = {0., 0., 0.};
//    // V += scale * ArithmeticWeightedSmoothingVector(p, [](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
//    // V += scale * AreaWeightedSmoothingVector(p, [](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
//    V += scale * DistorsionMeasureWeightedSmoothingVector2(p, [](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
//    return condition_Ua(V, p);
// };

// \label{BEM:vectorToNextSurface}
Tddd vectorToNextSurface(const networkPoint *p) {
   Tddd pX = RK_with_Ubuff(p);
   Tddd X, ret = {1E+20, 1E+20, 1E+20};
   double dxi = 0.3;
   T3Tdd t0t1_vertices = {{{1., 0.}, {1. - dxi, dxi}, {1. - dxi, 0.}}};
   if (p->Dirichlet) {
      auto t0t1 = SymmetricSubdivisionOfTriangle(t0t1_vertices, 10);
      int index = 0;
      networkFace *closest_face = nullptr;

      auto reset = [&]() {
         ret = {1E+20, 1E+20, 1E+20};
      };

      auto find_vector_using_pseudo_quad = [&](networkFace *f) {
         if (f->Points[0] == p)
            index = 0;
         else if (f->Points[1] == p)
            index = 1;
         else if (f->Points[2] == p)
            index = 2;
         else
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error");
         //// DodecaPoints dodecapoint(f, p, [](const networkLine *line) -> bool { return !line->CORNER; });
         for (const auto &vertices : f->dodecaPoints[index]->interpolate(t0t1, [&](networkPoint *p) -> Tddd { return RK_without_Ubuff(p); })) {
            X = Nearest(pX, vertices);
            if (Norm(ret) >= Norm(X - pX)) {
               ret = X - pX;
               closest_face = f;
            }
         }
      };

      for (const auto &f : p->getFaces()) {
         if ((p->Dirichlet && f->Dirichlet) || (p->Neumann && f->Neumann)) {
            if (_ALE_ON_LINEAR_ELEMENT_) {
               X = Nearest(pX, RK_without_Ubuff(f));
               if (Norm(ret) >= Norm(X - pX))
                  ret = X - pX;
            } else {
               find_vector_using_pseudo_quad(f);
               /* -------------------------------------------------------------------------- */
               // Tdd t0_range = {0., 1.};
               // std::function<Tddd(networkPoint *)> conversion = [&](networkPoint *p) -> Tddd { return RK_without_Ubuff(p); };
               // auto [T01, min] = f->dodecaPoints[index]->findMinimum(conversion, [&](const Tddd &X) { return Norm(X - pX); }, t0_range);
               // ret = f->dodecaPoints[index]->interpolate(T01, [&](networkPoint *p) { return RK_without_Ubuff(p); }) - pX;
            }
         }
      }
      // if (closest_face != nullptr) {
      //    reset();
      //    t0t1 = SymmetricSubdivisionOfTriangle(t0t1_vertices, 20);
      //    find_vector_using_pseudo_quad(closest_face);
      // }

      return ret;
   } else if (p->Neumann || p->CORNER) {
      Tddd to_corner = {0., 0., 0.};
      to_corner.fill(1E+20);
      Tddd to_structure_face = {0., 0., 0.};
      Tddd p_X_on_CORNER = pX;  //! 角にある場合は，正しく修正されるという意味．

      if (p->CORNER) {
         for (const auto &l : p->getLines()) {
            if (l->CORNER) {
               T2Tddd actual_corner_line = {RK_without_Ubuff(p), (1. - dxi) * RK_without_Ubuff(p) + dxi * RK_without_Ubuff((*l)(p))};
               X = Nearest(pX, actual_corner_line);
               if (Norm(to_corner) >= Norm(X - pX)) {
                  to_corner = X - pX;
                  p_X_on_CORNER = X;
               }
            }
         }
      }

      auto next_Vrtx = nextBodyVertices(bfs(p->getContactFaces(), 3));
      if (!next_Vrtx.empty()) {

         //! majorFlatFaces 法線方向が被らない面を抽出する
         std::vector<networkFace *> majorFlatFaces;
         double total_area = 0;
         for (const auto &pf0 : p->getFacesNeumann()) {
            total_area += pf0->area;
            bool largest_among_close_normal_face = true;
            for (const auto &pf1 : p->getFacesNeumann()) {
               if (isFlat(pf0->normal, pf1->normal, M_PI / 180.) && pf0->area < pf1->area)
                  largest_among_close_normal_face = false;
            }
            if (largest_among_close_normal_face)
               majorFlatFaces.push_back(pf0);
         }

         std::vector<double> distances;
         std::vector<Tddd> directions;
         //! 面p_fが干渉する，最も近い構造物面を抽出
         for (const auto &major : majorFlatFaces) {
            Tddd vec_to_closest_struct_face = {1E+20, 1E+20, 1E+20};
            bool found = false;
            for (const auto &struct_vertex : next_Vrtx) {
               if (isInContact(p_X_on_CORNER, major->normal, struct_vertex, p)) {
                  X = Nearest(p_X_on_CORNER, struct_vertex);
                  if (Norm(vec_to_closest_struct_face) >= Norm(X - p_X_on_CORNER)) {
                     //
                     const auto n = TriangleNormal(struct_vertex);
                     const double neglible_range = 1E-3 * std::sqrt(total_area);
                     const double angle = 60 * M_PI / 180.;
                     if (isFlat(n, X - p_X_on_CORNER, angle) || isFlat(n, -(X - p_X_on_CORNER), angle) || Norm(X - p_X_on_CORNER) < neglible_range) {
                        vec_to_closest_struct_face = X - p_X_on_CORNER;
                        found = true;
                     }
                  }
               }
            }
            if (found) {
               distances.push_back(Dot(vec_to_closest_struct_face, major->normal));
               directions.push_back(major->normal);
            }
         }

         to_structure_face = optimalVector(distances, directions, {0., 0., 0.});
         // to_structure_face = optimalVectorSVD(distances, directions);
      }

      return p_X_on_CORNER + to_structure_face - pX;
   }
   return {0., 0., 0.};
};

/*DOC_EXTRACT 0_3_1_INITIAL_VALUE_PROBLEM

### Arbitrary Lagrangian–Eulerian Methods (ALE)

水面を移動する物体が存在する場合，格子が潰れてしまうため，何らかの方法で格子を綺麗に保たなければ，長時間の計算は不可能である．
その方法の一つが，節点の位置を流速$`\nabla \phi`$ではなく任意のベクトル$`{\bf v}`$で移動させる，ALEである．
$`{\bf v}`$で移動する節点位置を$`\boldsymbol{\chi}(t)`$と置くと，
$`\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\nabla\phi\cdot\nabla\phi`$の代わりに，
$`\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\frac{d\boldsymbol\chi}{dt} \cdot\nabla\phi,\frac{d\boldsymbol\chi}{dt} = \bf v`$
を使って$`\phi`$を時間発展させると次時刻の節点位置$`\boldsymbol{\chi}(t+\delta t)`$での$`\phi`$が得られる．

ディリクレ節点（水面）：

求めた流速から，次の時刻の境界面$`\Omega(t+\Delta t)`$を見積もり，その面上で節点を移動させ歪さを解消する．
修正ベクトルは，$`\Delta t`$で割り，求めた流速$`\nabla \phi`$に足し合わせて，節点を時間発展させる．

ノイマン節点：

ノイマン節点も修正流速を加え時間発展させる．
ただし，ノイマン節点の修正流速に対しては，節点が水槽の角から離れないように，工夫を施している．

\ref{BEM:calculateVecToSurface}{`calculateVecToSurface`}で$`\Omega(t+\Delta t)`$上へのベクトルを計算する．

1. まず，\ref{BEM:vectorTangentialShift}{`vectorTangentialShift`}で接線方向にシフトし，
2. \ref{BEM:vectorToNextSurface}{`vectorToNextSurface`}で近くの$`\Omega(t+\Delta t)`$上へのベクトルを計算する．

#### 使っているALEの手法

`CircumradiusToInradius`をもとに，節点に隣接する面の歪みを評価し，その歪みに応じて節点を修正する．
面の歪みを修正するための節点の移動方向は，節点の向かい側にある辺の中点から辺が作る正三角形の高さの方向としている．
歪みを小さくするための方向を歪みの`CircumradiusToInradius`の微分から求めることもできるだろうが，プログラムが複雑になるためこの方法を採用している．

* 重みの最大値と最小値を設定している
* よりディリクレ面の歪みを緩和するように重みを大きくしている
* より喫水線の歪みを緩和するように重みを大きくしている

*/

// \label{BEM:calculateVecToSurface}

Tddd DistorsionMeasureWeightedSmoothingVector_modified(const networkPoint *p,
                                                       const std::array<double, 3> &current_pX,
                                                       std::function<Tddd(const networkPoint *)> position) {
   Tddd V = {0., 0., 0.}, X = {0., 0., 0.}, Xmid, vertical;
   double Wtot = 0, W, height, cuurent_height;
   const int max_sum_depth = 20;

   bool find_neumann_face = std::ranges::any_of(p->getFaces(), [](const auto &f) { return f->Neumann; });
   bool find_dirichlet_face = std::ranges::any_of(p->getFaces(), [](const auto &f) { return f->Dirichlet; });

   // std::vector<double> distances;
   // std::vector<Tddd> directions;

   std::vector<double> weights;
   std::vector<Tddd> positions;

   for (const auto &f : p->getFaces()) {
      if (find_dirichlet_face && f->Neumann)
         continue;
      auto [p0, p1, p2] = f->getPoints(p);
      auto X0 = current_pX;
      auto X1 = position(p1);
      auto X2 = position(p2);
      //! CircumradiusToInradiusの２乗を重みとしている
      // W = std::pow(CircumradiusToInradius(X0, X1, X2), 2);  // CircumradiusToInradius(X0, X1, X2)は最小で2
      W = CircumradiusToInradius(X0, X1, X2);  // CircumradiusToInradius(X0, X1, X2)は最小で2
      // W = CircumradiusToInradius(X0, X1, X2);  // CircumradiusToInradius(X0, X1, X2)は最小で2
      //! 重みの最大値と最小値を設定している
      // W = std::clamp(W, 4., 40.);
      // if ((find_neumann_face && find_dirichlet_face) && f->Neumann)
      //    W *= 0;
      //! より喫水線の歪みを緩和するように重みを大きくする
      // const int min_distance_from_CORNER = p0->minDepthFromCORNER + p1->minDepthFromCORNER + p2->minDepthFromCORNER;
      // W = std::clamp(W, 4., 1000.);

      for (int i = 0; i <= max_sum_depth; ++i)
         if (i == p0->minDepthFromMultipleNode + p1->minDepthFromMultipleNode + p2->minDepthFromMultipleNode) {
            // W *= std::pow(1.25, max_sum_depth - i);
            W *= std::pow(1.1, max_sum_depth - i);
            break;
         }

      //! よりディリクレ面の歪みを緩和するように重みを大きくする
      // if (f->Dirichlet)
      // for (const auto &p : f->getPoints())
      //    if (p->CORNER)
      //       W *= 3.;
      //! 　面それぞれに対する，修正方向は正三角形の高さの方向としている
      Wtot += W;
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3.) * 0.5;
      // cuurent_height = Norm(X0 - Xmid);

      auto error = Norm(height * vertical + Xmid - current_pX) / height;

      weights.push_back(error * W);
      positions.push_back(height * vertical + Xmid - current_pX);
      // FusedMultiplyIncrement(W, height * vertical + Xmid - current_pX, V);
      // directions.push_back(Normalize(height * vertical + Xmid - current_pX));
      // distances.push_back(Norm(height * vertical + Xmid - current_pX));
   }

   weights /= Sum(weights);
   V.fill(0.);
   for (int i = 0; i < weights.size(); ++i)
      FusedMultiplyIncrement(weights[i], positions[i], V);

   return V;
   // return optimalVector(distances, directions, V / Wtot);
};

Tddd ArithmeticWeightedSmoothingVector_modified(const networkPoint *p, const std::array<double, 3> &current_pX,
                                                std::function<Tddd(const networkPoint *)> position,
                                                const double power = 1.) {
   Tddd V = {0., 0, 0.}, qX;
   double Wtot = 0, norm;
   if (!isEdgePoint(p)) {
      if (p->CORNER) {
         for (const auto &l : p->getLines())
            if (l->CORNER) {
               auto q = (*l)(p);
               qX = position(q);
               V += std::pow(Norm(qX - current_pX), power - 1.) * (qX - current_pX);
               Wtot += 1.;
            }
      } else {
         for (const auto &q : p->getNeighbors()) {
            qX = position(q);
            V += std::pow(Norm(qX - current_pX), power - 1.) * (qX - current_pX);
            Wtot += 1.;
         }
      }
      return V / Wtot;
   } else
      return V;
};

void calculateVecToSurface(const Network &net, const int loop, const double coef) {
   auto points = ToVector(net.getPoints());
   for (const auto &p : points) {
      p->vecToSurface_BUFFER.fill(0.);
      p->vecToSurface.fill(0.);
   }

   //! 要素を整えるためのベクトル
   auto addVectorTangentialShift = [&](const int steps = 0) {
      if (coef == 0.) return;
      double scale = coef * ((steps + 1.) / (double)(loop));
#pragma omp parallel
      for (const auto &p : points)
#pragma omp single nowait
      {
         auto X = RK_with_Ubuff(p);
         Tddd V = {0., 0., 0.};
         const double a = (p->CORNER || steps <= 2) ? 0.5 : 1.;
         V = a * DistorsionMeasureWeightedSmoothingVector_modified(p, X, [&](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
         V += (1. - a) * ArithmeticWeightedSmoothingVector_modified(p, X, [&](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
         // V = 0.9 * V + 0.1 * VV;
         // double C = 0, size = 0.;
         // for (const auto &f : p->getFaces()) {
         //    auto [p0, p1, p2] = f->getPoints();
         //    C += Inradius(RK_without_Ubuff(p0), RK_without_Ubuff(p1), RK_without_Ubuff(p2));
         //    size += 1.;
         // }
         // C /= size;
         // C /= 5;
         // if (Norm(V) > C)
         //    V = C * Normalize(V);

         double flatness = p->getMinimalSolidAngle() / (2. * M_PI);
         if (flatness > 0.1)
            flatness = 1.;

         if (p->CORNER)
            V = scale * V * flatness / 2.;
         else
            V = scale * V * flatness;

         if (!isFinite(V))
            V.fill(0.);

         p->vecToSurface_BUFFER = condition_Ua(V, p);

         if (steps == 0)
            p->vecToSurface_BUFFER_BUFFER = p->vecToSurface_BUFFER;
         else
            p->vecToSurface_BUFFER_BUFFER = 0.8 * p->vecToSurface_BUFFER + 0.2 * p->vecToSurface_BUFFER_BUFFER;

         p->vecToSurface_BUFFER = p->vecToSurface_BUFFER_BUFFER;
      }

      for (const auto &p : points) {
         add_vecToSurface_BUFFER_to_vecToSurface(p);
         p->vecToSurface_BUFFER.fill(0.);
      }
   };

   //! 構造物に貼り付けるためのベクトル
   // double scale = 0.999;  //! 完全に壁に貼り付けると，初期状態に戻り調整ができない可能性があるので
   double scale = 1.;  //! 完全に壁に貼り付けると，初期状態に戻り調整ができない可能性があるので
   auto addVectorToNextSurface = [&]() {
#pragma omp parallel
      for (const auto &p : points)
#pragma omp single nowait
         p->vecToSurface_BUFFER = p->clungSurface = scale * vectorToNextSurface(p);

      for (const auto &p : points) {
         add_vecToSurface_BUFFER_to_vecToSurface(p);
         p->vecToSurface_BUFFER.fill(0.);
      }
   };

   TimeWatch watch;
   int steps;
   for (steps = 0; steps < loop; ++steps) {
      addVectorTangentialShift(steps);  // repeating this may led surface detaching
      std::cout << "Elapsed time for 1.vectorTangentialShift : " << watch() << " [s]" << std::endl;
      addVectorToNextSurface();
      std::cout << "Elapsed time for 2.vectorToNextSurface: " << watch() << " [s]" << std::endl;
   }
   addVectorToNextSurface();
};

//$ -------------------------------------------------------------------------- */
// b! -------------------------------------------------------------------------- */
// b!                               calculateVelocities                          */
// b! -------------------------------------------------------------------------- */

void calculateCurrentVelocities(const Network &net) {
#pragma omp parallel
   for (const auto &p : net.getPoints())
#pragma omp single nowait
   {
      p->U_update_BEM = p->U_BEM = gradPhi(p);
      if (p->Neumann)
         p->U_update_BEM = uNeumann(p);
   }
}

void calculateCurrentUpdateVelocities(const Network &net, const int loop, const double coef) {

   for (const auto &p : net.getPoints()) {
      if (!isFinite(p->U_update_BEM)) {
         std::cout << "p->X = " << p->X << std::endl;
         std::cout << "p->U_BEM = " << p->U_BEM << ", gradPhi(p) = " << gradPhi(p) << std::endl;
         std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
         std::cout << "p->vecToSurface = " << p->vecToSurface << std::endl;
         std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
         std::cout << "p->Neumann = " << p->Neumann << std::endl;
         std::cout << "p->CORNER = " << p->CORNER << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      }
   }

   calculateVecToSurface(net, loop, coef);
   const double too_fast = 1E+13;

   //! もしvecToSurfaceが大きすぎる場合は，係数を小さくして再計算
   for (auto i = 0; i < 10; ++i)
      if (std::ranges::any_of(net.getPoints(), [&](const auto &p) { return !isFinite(p->vecToSurface, too_fast); }))
         calculateVecToSurface(net, loop, coef * std::pow(0.1, i));

   //! もしvecToSurfaceが大きすぎる場合は，0にしてあきらめる
   if (std::ranges::any_of(net.getPoints(), [&](const auto &p) { return !isFinite(p->vecToSurface, too_fast); }))
      for (const auto &p : net.getPoints())
         p->vecToSurface.fill(0.);

   for (const auto &p : net.getPoints())
      p->U_update_BEM = p->RK_X.getVectorToReachAtNextTimeQ(p->vecToSurface + RK_without_Ubuff(p));
}

/*DOC_EXTRACT 0_7_OTHERS

### エネルギー保存則（計算精度のチェックに利用できる）

流体全体の運動エネルギーは，ラプラス方程式と発散定理を使うと，次のように境界面に沿った積分で表される．

```math
E_K =\frac{\rho}{2} \iint_\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
```

また，流体の位置エネルギーは，次のように表される．

```math
E_P = \frac{\rho}{2} \iint_\Gamma (0,0,g(z - z_0)^2) \cdot {\bf n} d\Gamma
```

<details>

---

<summary>
NOTE: なぜか？
</summary>

テンソルを使って考えてみると

```math
\begin{align*}
\nabla \cdot (\phi\nabla\phi) &= \frac{\partial\phi}{\partial x_i} \frac{\partial\phi}{\partial x_i} + \phi \frac{\partial^2\phi}{\partial x_i \partial x_i}\\
&= \nabla \phi \cdot \nabla \phi + \phi \nabla^2 \phi\\
&= \nabla \phi \cdot \nabla \phi
\end{align*}
```

よって，

```math
\iiint_\Omega \nabla\phi\cdot\nabla\phi d\Omega = \iiint_\Omega \nabla \cdot (\phi\nabla\phi) d\Omega = \iint_\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
```

---

```math
E_P = \rho g \iiint_\Omega (z - z_0) d\Omega
= \rho g \iiint_\Omega \frac{1}{2} \nabla \cdot (0,0,(z - z_0)^2) d\Omega
= \rho g \iint_\Gamma \frac{1}{2} (0,0,(z - z_0)^2) \cdot {\bf n} d\Gamma
= \frac{1}{2}\rho g \iint_\Gamma (z - z_0)^2 n_z d\Gamma
```

---

</details>

*/

double KinematicEnergy(const std::unordered_set<networkFace *> &faces) {
   double EK = 0;
   for (const auto &f : faces) {
      auto [p0, p1, p2] = f->getPoints();
      EK += Dot((p0->phiphin[0] + p1->phiphin[0] + p2->phiphin[0]) / 3. * gradPhi(f), f->normal) * f->area;
   }
   return _WATER_DENSITY_ * EK / 2.;
};

double PotentialEnergy(const std::unordered_set<networkFace *> &faces) {
   double EP = 0;
   for (const auto &f : faces) {
      auto [p0, p1, p2] = f->getPoints();
      auto intpX = interpolationTriangleLinear0101(T3Tddd{p0->X, p1->X, p2->X});
      for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
         EP += std::pow(intpX(x0, x1)[2], 2) * f->normal[2] * w0w1 * intpX.J(x0, x1);
   }
   return _WATER_DENSITY_ * EP / 2.;
};

double TotalEnergy(const std::unordered_set<networkFace *> &faces) {
   double EK = 0, EP = 0;
   for (const auto &f : faces) {
      auto [p0, p1, p2] = f->getPoints();
      EK += Dot((p0->phiphin[0] + p1->phiphin[0] + p2->phiphin[0]) / 3. * gradPhi(f), f->normal);
      auto intpX = interpolationTriangleLinear0101(T3Tddd{p0->X, p1->X, p2->X});
      for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
         EP += std::pow(intpX(x0, x1)[2], 2) * f->normal[2] * w0w1 * intpX.J(x0, x1);
   }
   return (EK + EP) * _WATER_DENSITY_ / 2.;
};

/*DOC_EXTRACT 0_7_OTHERS

### 内部流速の計算方法（使わなくてもいい）

[Fochesato2005](https://onlinelibrary.wiley.com/doi/10.1002/fld.838)にあるように，
流体内部の流速$`\nabla \phi`$は，BIEを微分して求めることができる．

```math
u({\bf a}) = \nabla\phi({\bf a}) = \int_{\partial \Omega} \frac{\partial Q}{\partial n} ({\bf x})Q({\bf x}, {\bf a}) - \phi({\bf x}) \frac{\partial Q}{\partial n} ({\bf x}, {\bf a}) d\Gamma
```

```math
Q({\bf x},{\bf a}) = \frac{{\bf r}}{4\pi r^3}, \quad \frac{\partial Q}{\partial n} ({\bf x},{\bf a}) = \frac{1}{4\pi r^3} (3 \mathbf{n} - (\mathbf{r} \cdot \mathbf{n}) \frac{\mathbf{r}}{r^2})
```

*/

#endif