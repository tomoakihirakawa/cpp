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
   if (p->CORNER) {
      auto next_normal = getNextNormalDirichlet_BEM(p);
      for (const auto &f : p->getFacesNeumann())
         VECTOR = Projection(VECTOR, Cross(next_normal, RK_without_Ubuff_Normal(f)));
   }
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

void add_vecToSurface_BUFFER_to_vecToSurface(const auto &p) {
   if (isFinite(p->vecToSurface_BUFFER))
      p->vecToSurface += p->vecToSurface_BUFFER;
};

std::array<double, 3> nextPositionOnBody(Network *net, networkPoint *p) {
   if (!net->isFixed && net->isRigidBody) {
      auto next_COM = net->RK_COM.getX(net->velocityTranslational());
      auto next_Q = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(net->velocityRotational(), net->Q)));
      return rigidTransformation(net->ICOM, next_COM, next_Q.Rv(), p->initialX);
   } else
      return p->X;
}

std::vector<T3Tddd> nextBodyVertices(const std::unordered_set<networkFace *> &Fs) {
   std::vector<T3Tddd> ret(Fs.size());
   int i = 0;
   for (auto &f : Fs) {
      auto net = f->getNetwork();
      auto [p0, p1, p2] = f->getPoints();
      if (net->isFixed)
         ret[i] = ToX(f);
      else if (net->isRigidBody) {
         // 現在のv^nを使って問題ない．加速度はいらない．
         // Quaternion q;
         // q = q.d_dt(net->velocityRotational());
         // Quaternion next_Q(net->RK_Q.getX(net->Q.AngularVelocityTodQdt(net->velocityRotational())));
         // auto next_COM = net->RK_COM.getX(net->velocityTranslational());

         auto next_COM = net->RK_COM.getX(net->velocityTranslational());
         auto next_Q = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(net->velocityRotational(), net->Q)));
         auto X_next = [&](const auto &p) { return rigidTransformation(net->ICOM, next_COM, next_Q.Rv(), p->initialX); };

         ret[i] = {X_next(p0), X_next(p1), X_next(p2)};
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

   if (p->Dirichlet) {
      for (const auto &f : p->getFaces()) {
         if ((p->Dirichlet && f->Dirichlet) || (p->Neumann && f->Neumann)) {
            auto actual_corner_face = RK_without_Ubuff(f);
            X = Nearest(pX, actual_corner_face);
            if (Norm(ret) >= Norm(X - pX))
               ret = X - pX;
         }
      }
      return ret;
   } else if (p->Neumann || p->CORNER) {
      Tddd to_corner = {0., 0., 0.};
      to_corner.fill(1E+20);
      Tddd to_structure_face = {0., 0., 0.};
      Tddd p_X_on_CORNER = pX;  //! 角にある場合は，正しく修正されるという意味．

      if (p->CORNER) {
         for (const auto &l : p->getLines()) {
            if (l->CORNER) {
               auto actual_corner_line = RK_without_Ubuff(p, (*l)(p));
               X = Nearest(pX, actual_corner_line);
               if (Norm(to_corner) >= Norm(X - pX)) {
                  to_corner = X - pX;
                  p_X_on_CORNER = X;
               }
            }
         }
      }

      auto next_Vrtx = nextBodyVertices(bfs(p->getContactFaces(), 4));
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

         std::vector<Tddd> projected_vec_for_each_p_f;
         // 面p_fが干渉する，最も近い構造物面を抽出
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
            if (found)
               projected_vec_for_each_p_f.push_back(Projection(vec_to_closest_struct_face, major->normal));
         }
         to_structure_face = optimumVector(projected_vec_for_each_p_f, {0., 0., 0.}, 1E-12);
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

*/

// \label{BEM:calculateVecToSurface}
void calculateVecToSurface(const Network &net, const int loop, const double coef = 0.1) {
   auto points = ToVector(net.getPoints());
   for (const auto &p : points) {
      p->vecToSurface_BUFFER.fill(0.);
      p->vecToSurface.fill(0.);
   }

   //! 要素を整えるためのベクトル
   auto addVectorTangentialShift = [&](const int k = 0) {
      // この計算コストは，比較的やすいので，何度も繰り返しても問題ない．
      // gradually approching to given a
      if (coef > 0.) {
         double aIN = coef, a;
         double scale = aIN * ((k + 1) / (double)(loop));
         //! ここを0.5とすると角が壊れる
         // const double scale = aIN;
         // if (scale < 0.0001)
         //    scale = 0.0001;
#pragma omp parallel
         for (const auto &p : points)
#pragma omp single nowait
         {
            auto X = RK_with_Ubuff(p);
            // auto V = (NeighborAverageSmoothingVector(p, X) + DistorsionMeasureWeightedSmoothingVector(p, X)) / 10.;
            // auto V = DistorsionMeasureWeightedSmoothingVector(p, X) / 10.;

            auto V = DistorsionMeasureWeightedSmoothingVector(p, X, [&](const networkPoint *p) { return RK_with_Ubuff(p); });

            double C = 0;
            double size = 0.;
            for (const auto &f : p->getFaces()) {
               auto [p0, p1, p2] = f->getPoints();
               C += Inradius(RK_with_Ubuff(p0), RK_with_Ubuff(p1), RK_with_Ubuff(p2));
               size += 1.;
            }
            V /= 10.;
            C /= size;
            C /= 10;

            if (Norm(V) > C)
               V = C * Normalize(V);
            if (!isFinite(V))
               V.fill(0.);

            double flatness = p->getMinimalSolidAngle() / (2. * M_PI);
            if (flatness > 0.1)
               flatness = 1.;
            V = scale * V * flatness;
            /* -------------------------------------------------------------------------- */
            p->vecToSurface_BUFFER = condition_Ua(V, p);
            // p->vecToSurface_BUFFER = vectorTangentialShift(p, scale);
         }
         for (const auto &p : points) {
            add_vecToSurface_BUFFER_to_vecToSurface(p);
            p->vecToSurface_BUFFER.fill(0.);
         }
      }
   };

   //! 構造物に貼り付けるためのベクトル
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
   for (auto kk = 0; kk < loop; ++kk) {
      addVectorTangentialShift(kk);  // repeating this may led surface detaching
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

void calculateCurrentUpdateVelocities(const Network &net, const int loop, const double coef = 0.1) {

   for (const auto &p : net.getPoints()) {
      if (!isFinite(p->U_update_BEM)) {
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

   calculateVecToSurface(net, loop, coef);
   const double too_fast = 1E+13;
   for (auto i = 0; i < 10; ++i) {
      if (std::ranges::any_of(net.getPoints(), [&](const auto &p) { return !isFinite(p->vecToSurface, too_fast); }))
         calculateVecToSurface(net, loop, coef * std::pow(0.1, i));
   }

   if (std::ranges::any_of(net.getPoints(), [&](const auto &p) { return !isFinite(p->vecToSurface, too_fast); }))
      for (const auto &p : net.getPoints())
         p->vecToSurface.fill(0.);

   for (const auto &p : net.getPoints()) {
      // dxdt_correct = p->vecToSurface / p->RK_X.getdt();
      // auto U_update_BEM = p->U_update_BEM + p->vecToSurface / p->RK_X.getdt();
      // if (!isFinite(U_update_BEM, too_fast) || !isFinite(p->vecToSurface, too_fast)) {
      //    std::cout << "p->X = " << p->X << std::endl;
      //    std::cout << "p->RK_X.getdt() = " << p->RK_X.getdt() << std::endl;
      //    std::cout << "p->U_update_BEM = " << U_update_BEM << std::endl;
      //    std::cout << "p->vecToSurface = " << p->vecToSurface << std::endl;
      //    std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
      //    std::cout << "p->Neumann = " << p->Neumann << std::endl;
      //    std::cout << "p->CORNER = " << p->CORNER << std::endl;
      //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // } else
      //    p->U_update_BEM += p->vecToSurface / p->RK_X.getdt();
      p->U_update_BEM = p->RK_X.getVectorToReachAtNextTimeQ(p->vecToSurface + RK_without_Ubuff(p));
   }
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