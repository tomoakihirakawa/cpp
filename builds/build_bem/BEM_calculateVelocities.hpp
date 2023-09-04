#ifndef BEM_calculateVelocities_H
#define BEM_calculateVelocities_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

//$ -------------------------------------------------------------------------- */
//$                         calculateVecToSurface                              */
//$ -------------------------------------------------------------------------- */

Tddd RK_with_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->U_update_BEM + p->vecToSurface / p->RK_X.getdt()); };
Tddd RK_with_Ubuff(const networkPoint *p, const Tddd &vecToSurface) { return p->RK_X.getX(p->U_update_BEM + vecToSurface / p->RK_X.getdt()); };
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
      if (p->CORNER) {
         for (const auto &f : p->getFacesNeumann())
            VECTOR = Projection(VECTOR, Cross(getNextNormalDirichlet_BEM(p), RK_without_Ubuff_Normal(f)));
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
   }
};

void add_vecToSurface_BUFFER_to_vecToSurface(const auto &p) {
   p->vecToSurface += p->vecToSurface_BUFFER;
   if (!isFinite(p->vecToSurface))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not finite");
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
         // 現在のv^nを使って問題ない．加速度はいらない．
         Quaternion q;
         q = q.d_dt(net->velocityRotational());
         auto COM = net->RK_COM.getX(net->velocityTranslational());
         Quaternion Q(net->RK_Q.getX(q()));
         auto X0 = Q.Rv(p0->initialX - net->ICOM) + COM;
         auto X1 = Q.Rv(p1->initialX - net->ICOM) + COM;
         auto X2 = Q.Rv(p2->initialX - net->ICOM) + COM;
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

// \label{BEM:vectorTangentialShift}
Tddd vectorTangentialShift(const networkPoint *p, const double scale = 1.) {
   Tddd vector_to_optimum_X = {0., 0., 0.}, pX = RK_with_Ubuff(p);
   double s = 0;
   // if (p->CORNER) {
   //    for (const auto &l : p->getLinesCORNER()) {
   //       vector_to_optimum_X += RK_with_Ubuff((*l)(p)) - pX;
   //       s += 1.;
   //    }
   // } else

   // auto decrese = [&](const Tddd &V) { return Norm(condition_Ua(V, p)) - Norm(V); };

   auto max_a = 0;
   for (const auto &f : p->getFaces()) {
      auto nP012 = RK_with_Ubuff(f->getPoints(p));
      auto &[np0x, np1x, np2x] = nP012;
      double a = std::log10(CircumradiusToInradius(nP012));
      if (std::ranges::any_of(f->getLines(), [](const auto &l) { return l->CORNER; })) {
         a *= std::pow(2., 3);
      } else if (std::ranges::any_of(f->getPoints(), [](const auto &p) { return p->isMultipleNode; })) {
         a *= std::pow(1.5, 2);
      } else if (std::ranges::any_of(f->getPoints(), [](const auto &p) {
                    return std::ranges::any_of(p->getFaces(), [](const auto &F) {
                       return std::ranges::any_of(F->getPoints(), [](const auto &q) { return q->isMultipleNode; });
                    });
                 })) {
         a *= 1.5;
      }

      auto optimum_position = Norm(np2x - np1x) * Normalize(Chop(np0x - np1x, np2x - np1x)) * sin(M_PI / 3.) + (np2x + np1x) / 2.;
      vector_to_optimum_X += a * condition_Ua(optimum_position - pX, p);
      s += a;
      if (max_a < a)
         max_a = a;
   }

   if (p->isMultipleNode)
      for (const auto &l : p->getLines()) {
         auto q = (*l)(p);
         if (q->isMultipleNode) {
            vector_to_optimum_X += max_a * condition_Ua(RK_with_Ubuff(q) - pX, p);
            s += max_a;
         }
      }
   vector_to_optimum_X /= s;
   // return condition_Ua(scale * vector_to_optimum_X, p);
   return scale * vector_to_optimum_X;
};

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
      Tddd to_structure_face = {0., 0., 0.};

      if (p->CORNER) {
         for (const auto &l : p->getLines()) {
            if (l->CORNER) {
               auto actual_corner_line = RK_without_Ubuff(p, (*l)(p));
               X = Nearest(pX, actual_corner_line);
               if (Norm(to_corner) >= Norm(X - pX))
                  to_corner = X - pX;
            }
         }
      }

      Tddd p_X_on_CORNER = pX + to_corner;
      auto next_Vrtx = nextBodyVertices(bfs(p->getContactFaces(), 3));
      if (!next_Vrtx.empty()) {
         std::vector<networkFace *> pf_to_check;
         for (const auto &pf : p->getFacesNeumann()) {
            if (std::ranges::none_of(pf_to_check, [pf](const auto f) { return VectorAngle(pf->normal, f->normal) < M_PI / 180.; }))
               pf_to_check.push_back(pf);
         }

         std::vector<Tddd> projected_vec_for_each_p_f;
         // 面p_fが干渉する，最も近い構造物面を抽出
         for (const auto &pf : pf_to_check) {
            Tddd vec_to_closest_struct_face = {1E+20, 1E+20, 1E+20};
            bool found = false;
            for (const auto &struct_vertex : next_Vrtx) {
               if (isInContact(p_X_on_CORNER, pf->normal, struct_vertex, p->radius)) {
                  X = Nearest(p_X_on_CORNER, struct_vertex);
                  if (Norm(vec_to_closest_struct_face) >= Norm(X - p_X_on_CORNER)) {
                     vec_to_closest_struct_face = X - p_X_on_CORNER;
                     found = true;
                  }
               }
            }
            if (found)
               projected_vec_for_each_p_f.push_back(Projection(vec_to_closest_struct_face, pf->normal));
         }
         to_structure_face = optimumVector(projected_vec_for_each_p_f, {0., 0., 0.}, 1E-12);
      }

      return p_X_on_CORNER + to_structure_face - pX;
   }
   return {0., 0., 0.};
};

/*DOC_EXTRACT 0_3_1_INITIAL_VALUE_PROBLEM

### 修正流速（激しい波の計算では格子が歪になりやすく，これがないと計算が難しい）

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
void calculateVecToSurface(const Network &net, const int loop, const bool do_shift = true) {
   for (const auto &p : net.getPoints()) {
      p->vecToSurface_BUFFER.fill(0.);
      p->vecToSurface.fill(0.);
   }

   auto addVectorTangentialShift = [&]() {
   // この計算コストは，比較的やすいので，何度も繰り返しても問題ない．
#pragma omp parallel
      for (const auto &p : net.getPoints())
#pragma omp single nowait
      {
         double a = 0.01;
         double scale = 5. * a;
         if (p->isMultipleNode && !p->CORNER)
            scale = 5. * a;
         else if (p->Neumann)
            scale = 5. * a;
         else if (p->CORNER)
            scale = 5. * a;

         p->vecToSurface_BUFFER = vectorTangentialShift(p, scale);
      }
      for (const auto &p : net.getPoints()) {
         add_vecToSurface_BUFFER_to_vecToSurface(p);
         p->vecToSurface_BUFFER.fill(0.);
      }
   };

   auto addVectorToNextSurface = [&]() {
#pragma omp parallel
      for (const auto &p : net.getPoints())
#pragma omp single nowait
         p->vecToSurface_BUFFER = p->clungSurface = vectorToNextSurface(p);

      for (const auto &p : net.getPoints()) {
         add_vecToSurface_BUFFER_to_vecToSurface(p);
         p->vecToSurface_BUFFER.fill(0.);
      }
   };

   TimeWatch watch;
   for (auto kk = 0; kk < loop; ++kk) {
      if (do_shift)
         // for (auto i = 0; i < 10; ++i)
         addVectorTangentialShift();  // repeating this may led surface detaching
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

// Tddd gradPhi(const networkPoint *const p) {
//    try {
//       Tddd u, accum;
//       double w = 0;
//       accum.fill(0.);
//       if (p->CORNER) {
//          for (const auto &f : p->getFacesDirichlet()) {
//             auto [p0, p1, p2] = f->getPoints(p);
//             u = gradTangential_LinearElement(Tddd{{std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)}}, T3Tddd{{ToX(p0), ToX(p1), ToX(p2)}});
//             u += f->normal * p->phinOnFace.at(nullptr);
//             accum += f->area * u;
//             w += f->area;
//          }
//       } else if (p->Dirichlet) {
//          for (const auto &f : p->getFaces()) {
//             auto [p0, p1, p2] = f->getPoints(p);
//             u = gradTangential_LinearElement(Tddd{{std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)}}, T3Tddd{{ToX(p0), ToX(p1), ToX(p2)}});
//             u += f->normal * p->phinOnFace.at(nullptr);
//             accum += f->area * u;
//             w += f->area;
//          }
//       } else {
//          for (const auto &f : p->getFaces()) {
//             auto [p0, p1, p2] = f->getPoints(p);
//             u = gradTangential_LinearElement(Tddd{{std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)}}, T3Tddd{{ToX(p0), ToX(p1), ToX(p2)}});
//             u += f->normal * p->phinOnFace.at(p->phinOnFace.count(f) ? f : nullptr);
//             accum += f->area * u;
//             w += f->area;
//          }
//       }
//       return accum / w;

//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorOff << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };

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

void calculateCurrentUpdateVelocities(const Network &net, const int loop, const bool do_shift = true) {

   calculateVecToSurface(net, loop, do_shift);

   for (const auto &p : net.getPoints()) {
      // dxdt_correct = p->vecToSurface / p->RK_X.getdt();
      p->U_update_BEM += p->vecToSurface / p->RK_X.getdt();

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