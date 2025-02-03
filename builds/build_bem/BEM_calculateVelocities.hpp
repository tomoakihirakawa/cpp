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
   for (const auto &f : p->getSurfaces()) {
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
   for (const auto &f : p->getSurfaces()) {
      a = RK_with_Ubuff_Area(f);
      // normal += a * RK_with_Ubuff_Normal(f);
      FusedMultiplyIncrement(a, RK_with_Ubuff_Normal(f), normal);
      total += a;
   }
   return Normalize(normal / total);
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
   auto COM_next = net->RK_COM.getX(velocity);
   auto Q_next = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(rotation, net->RK_Q.getX())));
   // return rigidTransformation(net->ICOM, COM_next, Q_next.Rv(), p->initialX);
   return rigidTransformation(net->ICOM, COM_next, Q_next.Rv(), p->initialX);
   // else return p->X;
}

std::vector<T3Tddd> nextBodyVertices(const std::unordered_set<networkFace *> &Fs) {
   std::vector<T3Tddd> ret(Fs.size());
   int i = -1;
   for (auto &f : Fs) {
      i++;
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

         if (std::ranges::all_of(velocity, [](double v) { return v == 0; }) && std::ranges::all_of(rotation, [](double v) { return v == 0; })) {
            ret[i] = {p0->initialX, p1->initialX, p2->initialX};
            continue;
         }

         for (auto i = 0; i < 3; ++i) {
            if (net->isFixed[i])
               velocity[i] = 0;
            if (net->isFixed[i + 3])
               rotation[i] = 0;
         }
         auto next_COM = net->RK_COM.getX(velocity);
         auto next_Q = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(rotation, net->RK_Q.getX())));
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

// やはりこの辺りなのだろう

// \label{BEM:vectorToNextSurface}
const double dxi = 0.4;
const T3Tdd t0t1_vertices = {{{1., 0.}, {1. - dxi, dxi}, {1. - dxi, 0.}}};
const auto t0t1 = SymmetricSubdivisionOfTriangle(t0t1_vertices, 3);
const auto t0t1_high = SymmetricSubdivisionOfTriangle(t0t1_vertices, 10);

Tddd vectorToNextSurface(const networkPoint *p) {
   Tddd X_shifted = RK_with_Ubuff(p);
   Tddd X_not_shifted = RK_without_Ubuff(p);
   Tddd X, ret = {1E+20, 1E+20, 1E+20};
   if (p->Dirichlet) {

      auto Nearest_pseudo = [&](const auto X_shifted, const networkFace *f, const auto &t0t1) {
         int index = -1;
         for (auto i = 0; i < 3; ++i)
            if (f->Points[i] == p) {
               index = i;
               break;
            }
         if (index == -1)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error");
         return Nearest(X_shifted, f->dodecaPoints[index]->interpolate(t0t1, [&](networkPoint *p) -> Tddd { return RK_without_Ubuff(p); }));
      };

      networkFace *closest_face = nullptr;

      for (const auto &f : p->getSurfaces()) {
         if ((p->Dirichlet && f->Dirichlet) || (p->Neumann && f->Neumann)) {

            if (_ALE_ON_LINEAR_ELEMENT_)
               X = Nearest(X_shifted, RK_without_Ubuff(f));
            else
               X = Nearest_pseudo(X_shifted, f, t0t1);

            if (Norm(ret) >= Norm(X - X_shifted)) {
               ret = X - X_shifted;
               closest_face = f;
            }
         }
      }
      if (closest_face != nullptr)
         ret = Nearest_pseudo(X_shifted, closest_face, t0t1_high) - X_shifted;

      return ret;
   } else if (p->Neumann || p->CORNER) {
      Tddd to_structure_face = {0., 0., 0.};

      auto next_Vrtx = nextBodyVertices(bfs(p->getContactFaces(), 1));
      if (!next_Vrtx.empty()) {
         auto faces = p->getSurfaces();
         //! facesはDirichet面も含んでいる．これはCORNERの張り付きにおいて重要だと思われる．これがないと，コーナーは無限に上下方向に移動しても条件を満たすことになるだろう．
         //! Dirichet面の法線方向距離も小さくする様にしている．
         double radius = std::sqrt(std::accumulate(faces.begin(), faces.end(), 0., [](double a, const auto &f) { return a + f->area / M_PI; }));
         std::vector<double> distances;
         std::vector<Tddd> directions;

         auto add_vector = [&](const Tddd &V, const Tddd &n) {
            // const double small_angle = 5. * M_PI / 180.;
            // bool found = false;
            // for (auto i = 0; i < directions.size(); ++i)
            //    if (Dot(V, n) < distances[i] && (isFlat(n, directions[i], small_angle) || isFlat(n, -directions[i], small_angle))) {
            //       distances[i] = Dot(V, n);
            //       directions[i] = n;
            //       found = true;
            //       break;
            //    }
            // if (!found)
            {
               distances.push_back(Dot(V, n));
               directions.emplace_back(n);
            }
         };

         //! 面p_fが干渉する，最も近い構造物面を抽出
         const double angle = 50 * M_PI / 180.;
         for (const auto &f : faces /*p->getFacesNeumann()*/) {
            for (const auto &struct_vertex : next_Vrtx) {
               if (isInContact(X_not_shifted /*シフトなし*/, f->normal, struct_vertex, p)) {
                  /*以下はシフトした場所からのベクトルを計算*/
                  Tddd X_nearest = Nearest(X_shifted, struct_vertex), n = TriangleNormal(struct_vertex);
                  Tddd To_nearest = X_nearest - X_shifted;
                  double distance = Norm(To_nearest);
                  if ((isFlat(n, To_nearest, angle) || isFlat(n, -To_nearest, angle)) && (0.3 * radius >= distance))
                     add_vector(To_nearest, n);
                  else if ((0.03 * radius >= distance))
                     add_vector(To_nearest, n);
               }
            }
         }

         if (p->CORNER) {
            /*以下はシフトした場所からのベクトルを計算*/
            Tddd n = getNextNormalDirichlet_BEM(p);
            Tddd To_nearest = RK_without_Ubuff(p) - X_shifted;
            add_vector(To_nearest, n);
         }

         to_structure_face = optimalVector(distances, directions, {0., 0., 0.});
      }

      return to_structure_face;
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
   const int max_sum_depth = 20;
   auto faces = p->CORNER ? p->getFacesDirichlet() : p->getSurfaces();
   std::vector<double> weights;
   weights.reserve(faces.size());
   std::vector<Tddd> positions;
   positions.reserve(faces.size());
   Tddd X0, X1, X2, Xmid, vertical, X_ideal, To_ideal;
   double W, height, normalized_discrepancy;
   for (const auto &f : faces) {
      auto [p0, p1, p2] = f->getPoints(p);
      X0 = current_pX;
      X1 = position(p1);
      X2 = position(p2);
      //! 重みの計算
      W = CircumradiusToInradius(X0, X1, X2);  // CircumradiusToInradius(X0, X1, X2)は最小で2
      for (int i = 0; i <= max_sum_depth; ++i)
         if (i == p0->minDepthFromMultipleNode + p1->minDepthFromMultipleNode + p2->minDepthFromMultipleNode) {
            W *= std::pow(1.1, max_sum_depth - i);
            break;
         }
      //! 　面それぞれに対する，修正方向は正三角形の高さの方向としている
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3.) * 0.5;
      X_ideal = height * vertical + Xmid;
      To_ideal = X_ideal - current_pX;
      normalized_discrepancy = Norm(To_ideal) / height;
      //! 重みと位置を保存
      weights.push_back(normalized_discrepancy * W);
      positions.emplace_back(To_ideal);
   }

   weights /= Sum(weights);
   Tddd V = {0., 0., 0.};
   for (int i = 0; i < weights.size(); ++i)
      V += weights[i] * positions[i];
   return V;
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

   //! 初期化
   for (const auto &p : points) {
      p->vecToSurface.fill(0.);
      p->vecToSurface_BUFFER.fill(0.);
      p->vecToSurface_BUFFER_BUFFER.fill(0.);
   }

   /* ------------------------------ 要素を整えるためのベクトル ----------------------------- */
   auto addVectorTangentialShift = [&](const int steps = 0) {
      auto approxRadius = [&](const networkPoint *p) {
         double area = 0;
         for (const auto &f : p->getSurfaces())
            area += f->area;
         return std::sqrt(area / M_PI);
      };

      if (coef == 0.)
         return;

      double scale = coef * ((steps + 0.1) / (double)(loop));

      Tddd V;
      for (const auto &p : points) {
         auto current_pX = RK_with_Ubuff(p);
         V = scale * DistorsionMeasureWeightedSmoothingVector_modified(p, current_pX, [&](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
         if (!isFinite(V))
            V.fill(0.);
         p->vecToSurface_BUFFER = V;
      }

      for (const auto &p : points) {
         p->vecToSurface_BUFFER = 0.3 * p->vecToSurface_BUFFER_BUFFER + 0.7 * std::min(Norm(p->vecToSurface_BUFFER), approxRadius(p) * 0.3) * Normalize(p->vecToSurface_BUFFER);
         p->vecToSurface_BUFFER_BUFFER = p->vecToSurface_BUFFER;
         add_vecToSurface_BUFFER_to_vecToSurface(p);
         p->vecToSurface_BUFFER.fill(0.);
      }
   };

   /* ---------------------------- 構造物に貼り付けるためのベクトル ---------------------------- */

   auto addVectorToNextSurface = [&]() {
#pragma omp parallel
      for (const auto &p : points)
#pragma omp single nowait
         p->vecToSurface_BUFFER = p->clungSurface = vectorToNextSurface(p);

      for (const auto &p : points) {
         add_vecToSurface_BUFFER_to_vecToSurface(p);
         p->vecToSurface_BUFFER.fill(0.);
      }
   };

   TimeWatch watch;
   for (auto steps = 0; steps < loop; ++steps) {
      addVectorToNextSurface();
      addVectorToNextSurface();
      std::cout << "Elapsed time for 1. addVectorToNextSurface: " << watch() << " [s]" << std::endl;
      addVectorTangentialShift(steps);  // repeating this may led surface detaching
      std::cout << "Elapsed time for 2. addVectorTangentialShift : " << watch() << " [s]" << std::endl;
   }
   addVectorToNextSurface();
   addVectorToNextSurface();
};

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
         p->U_update_BEM = contactNormalVelocity(p);
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
   // const double too_fast = 1E+13;

   // auto nogood = [](auto p) {
   //    auto r = Norm(p->RK_X.getX(p->vecToSurface) - p->X);
   //    return std::ranges::any_of(p->getSurfaces(), [&](const auto &f) { return r > std::sqrt(f->area/M_PI); }); };

   // //! もしvecToSurfaceが大きすぎる場合は，係数を小さくして再計算
   // for (auto i = 0; i < 10; ++i)
   //    // if (std::ranges::any_of(net.getPoints(), [&](const auto &p) { return !isFinite(p->vecToSurface, too_fast); }))
   //    if (std::ranges::any_of(net.getPoints(), [&](const auto &p) { return nogood(p); }))
   //       calculateVecToSurface(net, loop, coef * std::pow(0.1, i));

   //! もしvecToSurfaceが大きすぎる場合は，0にしてあきらめる
   // if (std::ranges::any_of(net.getPoints(), [&](const auto &p) { return !isFinite(p->vecToSurface, too_fast); }))
   // if (std::ranges::any_of(net.getPoints(), [&](const auto &p) { return nogood(p); }))
   //    for (const auto &p : net.getPoints())
   //       p->vecToSurface.fill(0.);

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

double KinematicEnergy(const auto &faces) {
   double EK = 0;
   for (const auto &f : faces) {
      auto [p0, p1, p2] = f->getPoints();
      EK += Dot((p0->phiphin[0] + p1->phiphin[0] + p2->phiphin[0]) / 3. * gradPhi(f), f->normal) * f->area;
   }
   return _WATER_DENSITY_ * EK / 2.;
};

double PotentialEnergy(const auto &faces) {
   double EP = 0;
   for (const auto &f : faces) {
      auto [p0, p1, p2] = f->getPoints();
      auto intpX = interpolationTriangleLinear0101(T3Tddd{p0->X, p1->X, p2->X});
      for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
         EP += std::pow(intpX(x0, x1)[2], 2) * f->normal[2] * w0w1 * intpX.J(x0, x1);
   }
   return _WATER_DENSITY_ * EP / 2.;
};

double TotalEnergy(const auto &faces) {
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