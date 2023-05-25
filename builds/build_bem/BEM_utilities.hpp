#ifndef BEM_utilities_H
#define BEM_utilities_H

#include "Network.hpp"

using V_i = std::vector<int>;
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
using V_Netp = std::vector<Network *>;
using V_netFp = std::vector<networkFace *>;
using VV_netFp = std::vector<V_netFp>;

/* -------------------------------------------------------------------------- */
T6d velocity(const std::string &name, const std::vector<std::string> strings, networkPoint *p, double t) {
   if (name.contains("linear") && (name.contains("traveling") || name.contains("wave"))) {
      if (strings.size() == 6) {
         double start = stod(strings[1] /*start*/);
         if (t >= start) {
            double a = std::abs(stod(strings[2] /*a*/));
            double w = std::abs(2 * M_PI / stod(strings[3] /*T*/));
            double h = std::abs(stod(strings[4] /*h*/));
            double z_surface = std::abs(stod(strings[5] /*z_surface*/));
            auto [x, y, z] = p->X - Tddd{0., 0., z_surface};
            DispersionRelation DS(w, h);
            double k = std::abs(DS.k);
            // std::cout << "a = " << a
            //           << ", {w,k} = {" << w << "," << k << "}"
            //           << ", h = " << h
            //           << ", z_surface = " << z_surface
            //           << ", {T, L} = {" << DS.T << ", " << DS.L << "}" << std::endl;
            // t += M_PI / 2. / w;
            return {a * w * cosh(k * (z + h)) / sinh(k * h) * cos(w * (t - start) - k * x) +
                        w * k * a * a / 2 * (cosh(2 * k * (z + h)) - cos(2 * (w * (t - start) - k * x))) / std::pow(sinh(k * h), 2),
                    0.,
                    -a * w * sinh(k * (z + h)) / sinh(k * h) * sin(w * (t - start) - k * x),
                    0.,
                    0.,
                    0.};
         } else
            return {0., 0., 0., 0., 0., 0.};
      } else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be == 6");
   } else if (name.contains("weird")) {
      if (strings.size() == 6) {
         double start = stod(strings[1] /*start*/);
         double a = std::abs(stod(strings[2] /*a*/));
         double w = std::abs(2 * M_PI / stod(strings[3] /*T*/));
         double h = std::abs(stod(strings[4] /*h*/));
         double z_surface = std::abs(stod(strings[5] /*z_surface*/));
         auto [x, y, z] = p->X;
         DispersionRelation DS(w, h);
         double k = std::abs(DS.k);
         std::cout << "a = " << a
                   << ", {w,k} = {" << w << "," << k << "}"
                   << ", h = " << h
                   << ", z_surface = " << z_surface
                   << ", {T, L} = {" << DS.T << ", " << DS.L << "}" << std::endl;
         return {a * w * cos(w * t - k * z),
                 0.,
                 0.,
                 0.,
                 0.,
                 0.};
      } else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be == 6");
   }
   return {0., 0., 0., 0., 0., 0.};
};
T6d velocity(const std::string &name, const std::vector<std::string> strings, const double t) {
   auto g = _GRAVITY_;
   if (name == "Goring1979") {
      double h = 0.25;
      double H = 0.1 * h;  // 造波する初期入射波の波高
      double x = 0;
      double c = std::sqrt(g * (H + h));
      double kappa = std::sqrt(3. * H / (4. * h * h * h));
      double start = stod(strings[1] /*start*/);
      double eta = H * std::pow(1. / cosh(kappa * (-c * (t - start))), 2.);
      return {c * eta / (h + eta), 0., 0., 0, 0, 0};
   } else if (name == "Retzler2000") {
      const std::vector<Tdd> sample = {
          {-0.15000000000000002, 0.},
          {-0.11, 0.},
          {-0.1, 0.},
          {-0.09, 0.},
          {-0.08, 0.},
          {-0.07, 0.},
          {-0.06, 0.},
          {-0.050409836065573754, 0.007920792079208039},
          {-0.02397540983606558, 0.06930693069306948},
          {-0.0006147540983606481, 0.13762376237623775},
          {0.0282786885245902, 0.24851485148514862},
          {0.05040983606557381, 0.3594059405940595},
          {0.05901639344262294, 0.40297029702970305},
          {0.06946721311475412, 0.4425742574257427},
          {0.08299180327868855, 0.4792079207920793},
          {0.10020491803278692, 0.516831683168317},
          {0.11065573770491804, 0.5415841584158416},
          {0.12110655737704923, 0.5663366336633664},
          {0.132172131147541, 0.5910891089108912},
          {0.15000000000000008, 0.6198019801980199},
          {0.1616803278688525, 0.6306930693069308},
          {0.17643442622950822, 0.6445544554455447},
          {0.18872950819672135, 0.6603960396039605},
          {0.20102459016393448, 0.6792079207920793},
          {0.21823770491803285, 0.6801980198019802},
          {0.23053278688524592, 0.6445544554455447},
          {0.25081967213114753, 0.5742574257425743},
          {0.27848360655737703, 0.40297029702970305},
          {0.30000000000000004, 0.23564356435643574},
          {0.319672131147541, 0.10396039603960416},
          {0.3319672131147542, 0.02772277227722786},
          {0.35040983606557385, -0.03960396039603942},
          {0.37069672131147546, -0.085148514851485},
          {0.38913934426229513, -0.08316831683168302},
          {0.4002049180327869, -0.06831683168316827},
          {0.4254098360655738, -0.03069306930693061},
          {0.45000000000000007, -0.009900990099009799},
          {0.5, 0.},
          {0.55, 0.},
          {0.6, 0.},
          {0.65, 0.},
          {0.7, 0.}};
      double start = stod(strings[1] /*start*/);
      const auto intp = InterpolationBspline(3, sample);
      return {intp(t - start), 0., 0., 0., 0., 0.};
   } else if (name == "Chaplin2000") {
      double start = stod(strings[1] /*start*/);
      if (t < start)
         return {0., 0., 0., 0., 0., 0.};
      double h = 0.5;
      // double w = 1.257 / std::sqrt(h / g); /*5.57065*/
      // double A = 0.046 * h;
      // double A = 0.02 * h;
      double w, A;
      if (strings.size() > 3) {
         A = stod(strings[2] /*start*/);
         w = stod(strings[3] /*start*/);
         std::cout << "A = " << A << ", w = " << w << std::endl;
      } else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be > 3. amplitude and frequency");

      auto v = A * w * sin(w * (t - start));
      return {0., v, 0., 0., 0., 0.};
   } else if (name == "flap") {
      // start,A, T, h, l
      // Schaffer,H.A. : Second-order wavemaker theory for irregular waves, Ocean Engineering, 23(1), 47-88, (1996)
      double start = stod(strings[1] /*start*/);
      double A, w, h, l, d, k;
      if (strings.size() > 7) {
         A = std::abs(stod(strings[2] /*A*/));
         w = std::abs(2 * M_PI / stod(strings[3] /*T*/));
         h = std::abs(stod(strings[4] /*h*/));
         l = std::abs(stod(strings[5] /*l*/));
         DispersionRelation DS(w, h);
         k = std::abs(DS.k);
         double d = (l >= 0 ? d : -l);
         std::cout << "A = " << A << ", w = " << w << ", k = " << k << ", h = " << h
                   << ", d = " << d << ", {T, L} = {" << DS.T << ", " << DS.L << "}" << std::endl;
         Tddd axis = {stod(strings[6]), stod(strings[7]), stod(strings[8])};
         auto [wx, wy, wz] = Normalize(axis) * ArcTan((A * g * k * (1 + 2 * h * k * Csch(2 * h * k)) * Sin(t * w)),
                                                      (2. * (-g + (h + l) * Power(w, 2) + g * Cosh(d * k) * Sech(h * k))));
         return {0., 0., 0., wx, wy, wz};
      } else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be > 3. amplitude and frequency");
   } else if (name.contains("sinusoidal") || name.contains("sin")) {
      if (strings.size() == 4) {
         double start = stod(strings[1] /*start*/);
         if (t >= start) {
            double a = std::abs(stod(strings[2] /*a*/));
            double w = std::abs(2 * M_PI / stod(strings[3] /*T*/));
            return {a * w * cos(w * t), 0., 0., 0., 0., 0.};
         }
      } else {
         std::stringstream ss;
         int i = 0;
         for (const auto &s : strings)
            ss << i++ << ":" << s << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }
   }
   return {0., 0., 0., 0., 0., 0.};
};
/* -------------------------------------------------------------------------- */

netFp NearestContactFace(const networkPoint *const p) { return std::get<1>(Nearest_(p->X, p->getContactFaces())); };
netFp NearestContactFace(const networkFace *const f_IN) {
   std::unordered_set<networkFace *> faces;
   std::ranges::for_each(f_IN->getPoints(), [&](const auto &q) { faces.insert(q->getContactFaces().begin(), q->getContactFaces().end()); });
   return std::get<1>(Nearest_(f_IN->center, faces));
};
netFp NearestContactFace(const networkPoint *const p, const networkFace *const f_normal) {
   Tddd r = {1E+100, 1E+100, 1E+100}, X;
   networkFace *ret = nullptr;
   for (const auto &f_target : bfs(p->getContactFaces(), 2))
      if (isInContact(p, f_normal, f_target)) {
         X = Nearest(p->X, ToX(f_target));
         if (Norm(r) >= Norm(X - p->X)) {
            r = X - p->X;
            ret = f_target;
         }
      }
   return ret;
};
std::tuple<netFp, Tddd> NearestContactFace_(const networkPoint *const p, const networkFace *const f_normal) {
   Tddd r = {1E+100, 1E+100, 1E+100}, X;
   networkFace *ret = nullptr;
   for (const auto &f_target : bfs(p->getContactFaces(), 2))
      if (isInContact(p, f_normal, f_target)) {
         X = Nearest(p->X, ToX(f_target));
         if (Norm(r) >= Norm(X - p->X)) {
            r = X - p->X;
            ret = f_target;
         }
      }
   return {ret, r};
};

//$ --------------------------------------------------------------- */

std::tuple<Tddd, double> uNeumann_(const networkPoint *const p, const networkFace *const f_normal) {
   auto [f, vToContact] = NearestContactFace_(p, f_normal);
   if (f) {
      double weight = w_Bspline5(Norm(vToContact), p->radius);
      if (f->getNetwork()->isRigidBody) {
         return {f->getNetwork()->velocityRigidBody(p->X), weight};
      } else if (f->getNetwork()->isSoftBody) {
         auto [p0, p1, p2] = f->getPoints();
         auto v0 = p0->velocityTranslational();
         auto v1 = p1->velocityTranslational();
         auto v2 = p2->velocityTranslational();
         auto [t0, t1, X] = Nearest_(p->X, ToX(f));
         return {v0 * t0 + v1 * t1 + v2 * (1 - t0 - t1), weight};
      }
   }
   return {{0., 0., 0.}, 0.};
};

std::tuple<Tddd, double> accelNeumann_(const networkPoint *const p, const networkFace *const f_normal) {
   auto [f, vToContact] = NearestContactFace_(p, f_normal);
   if (f) {
      double weight = w_Bspline5(Norm(vToContact), p->radius);
      if (f->getNetwork()->isRigidBody) {
         return {f->getNetwork()->accelRigidBody(p->X), weight};
      } else if (f->getNetwork()->isSoftBody) {
         auto [p0, p1, p2] = f->getPoints();
         auto v0 = p0->accelTranslational();
         auto v1 = p1->accelTranslational();
         auto v2 = p2->accelTranslational();
         auto [t0, t1, X] = Nearest_(p->X, ToX(f));
         return {v0 * t0 + v1 * t1 + v2 * (1 - t0 - t1), weight};
      }
   }
   return {{0., 0., 0.}, 0.};
};

Tddd uNeumann(const networkPoint *const p, const networkFace *const f_normal) {
   return std::get<0>(uNeumann_(p, f_normal));
};

Tddd accelNeumann(const networkPoint *const p, const networkFace *const f_normal) {
   return std::get<0>(accelNeumann_(p, f_normal));
};

Tddd uNeumann(const networkPoint *const p) {
   std::vector<Tddd> V;
   std::vector<double> W;
   const Tddd init = {0., 0., 0.};
   for (const auto &f_normal : p->getFacesNeumann()) {
      auto [v, w] = uNeumann_(p, f_normal);
      if (w > 1E-20) {
         V.emplace_back(v);
         W.emplace_back(w);
      }
   }
   if (!V.empty()) {
      auto ret = optimumVector(V, init, W);
      if (isFinite(ret))
         return ret;
   }
   return init;
};

Tddd accelNeumann(const networkPoint *const p) {
   std::vector<Tddd> V;
   std::vector<double> W;
   const Tddd init = {0., 0., 0.};
   for (const auto &f_normal : p->getFacesNeumann()) {
      auto [v, w] = accelNeumann_(p, f_normal);
      if (w > 1E-20) {
         V.emplace_back(v);
         W.emplace_back(w);
      }
   }
   if (!V.empty()) {
      auto ret = optimumVector(V, init, W);
      if (isFinite(ret))
         return ret;
   }
   return init;
};

Tddd uNeumann(const networkFace *const f_normal) {
   auto [p0, p1, p2] = f_normal->getPoints();
   return (uNeumann(p0, f_normal) + uNeumann(p1, f_normal) + uNeumann(p2, f_normal)) / 3.;
};
Tddd accelNeumann(const networkFace *const f_normal) {
   auto [p0, p1, p2] = f_normal->getPoints();
   return (accelNeumann(p0, f_normal) + accelNeumann(p1, f_normal) + accelNeumann(p2, f_normal)) / 3.;
};

//$ --------------------------------------------------------------- */

using map_P_d = std::map<netP *, double>;
using map_P_Vd = std::map<netP *, V_d>;
using map_P_VVd = std::map<netP *, VV_d>;
using map_F_P_Vd = std::map<netF *, map_P_Vd>;
using map_P_P_Vd = std::map<netP *, map_P_Vd>;
using pair_PB = std::pair<netP *, bool>;
using map_pairPB_Tdd = std::unordered_map<std::tuple<netP *, bool, netF *>, Tdd>;
using map_pairPB_pairPB_Tdd = std::unordered_map<std::tuple<netP *, bool, netF *>, map_pairPB_Tdd>;
using map_P_P_Tdd = std::map<netP *, std::map<netP *, Tdd>>;
using map_P_F_P_Vd = std::map<netP *, map_F_P_Vd>;
using VV_SorIorMap = std::vector<std::vector<std::variant<std::string, int, map_P_Vd>>>;

V_netFp takeFaces(const V_Netp &nets) {
   V_netFp ret({});
   for (const auto &n : nets)
      ret.insert(ret.end(), n->getFaces().begin(), n->getFaces().end());
   return DeleteDuplicates(ret);
};

//@ ------------------------------------------------------ */

Tddd gradTangential_LinearElement(const Tddd &phi012, const T3Tddd &X012) {
   // これは{x,y,z}座標系での結果
   auto [X0, X1, X2] = X012;
   auto [phi0, phi1, phi2] = phi012;
   auto n = TriangleNormal(X012);
   return Cross(n, phi0 * (X2 - X1) + phi1 * (X0 - X2) + phi2 * (X1 - X0)) / (2 * TriangleArea(X012));
};

T3Tddd gradTangential_LinearElement(const T3Tddd &V012, const T3Tddd &X012) {
   auto [Vx012, Vy012, Vz012] = Transpose(V012);
   return {gradTangential_LinearElement(Vx012, X012),
           gradTangential_LinearElement(Vy012, X012),
           gradTangential_LinearElement(Vz012, X012)};
};

T3Tddd grad_U_tangential_LinearElement(const networkFace *const f) {
   auto [p0, p1, p2] = f->getPoints();
   return gradTangential_LinearElement({p0->U_BEM, p1->U_BEM, p2->U_BEM}, ToX(f));
};

Tddd grad_LinearElement(const Tddd &F012, const T3Tddd &X012, const Tddd &F_n) {
   //! 三角要素の節点の情報変数F0,F1,F2から，三角要素上でのgrad(F)を計算する．
   return gradTangential_LinearElement(F012, X012) + F_n;
};

T3Tddd grad_LinearElement(const T3Tddd &F012, const T3Tddd &X012, const T3Tddd &F_n) {
   //! 三角要素の節点の情報変数F0,F1,F2から，三角要素上でのgrad(F)を計算する．
   return {grad_LinearElement(std::get<0>(F012), X012, std::get<0>(F_n)),
           grad_LinearElement(std::get<1>(F012), X012, std::get<1>(F_n)),
           grad_LinearElement(std::get<2>(F012), X012, std::get<2>(F_n))};
};

T3Tddd OrthogonalBasis(const Tddd &n) {
   Tddd s0 = Chop(Tddd{1, 0, 0}, n);
   if (Norm(s0) < 1E-1)
      s0 = Chop(Tddd{0, 1, 0}, n);
   s0 = Normalize(s0);
   Tddd s1 = Normalize(Cross(n, s0));
   return {n, s0, s1};
};

//! 注意 成分がx,y,z成分ではないので注意
T3Tddd grad_U_LinearElement(const networkFace *const F, const T3Tddd &orthogonal_basis) {
   auto [s0, s1, s2] = orthogonal_basis;
   auto U_ = [&](const auto &p) { return Tddd{Dot(p->U_BEM, s0), Dot(p->U_BEM, s1), Dot(p->U_BEM, s2)}; };
   auto [P0, P1, P2] = F->getPoints();
   T3Tddd X012 = {P0->X, P1->X, P2->X};
   // 接線方向微分
   auto f_s0_tag = gradTangential_LinearElement(U_(P0) /*s座標の流速*/, X012);
   auto f_s1_tag = gradTangential_LinearElement(U_(P1) /*s座標の流速*/, X012);
   auto f_s2_tag = gradTangential_LinearElement(U_(P2) /*s座標の流速*/, X012);
   // s座標に置き換える
   // auto f_s0s0 = Dot(f_s0_tag, s0);  // n方向なので成分はないはず. n is s0 direction
   auto f_s0s1 = Dot(f_s0_tag, s1);
   auto f_s0s2 = Dot(f_s0_tag, s2);
   //
   // auto f_s1s0 = Dot(f_s1_tag, s0);  // n方向なので成分はないはず
   auto f_s1s1 = Dot(f_s1_tag, s1);
   auto f_s1s2 = Dot(f_s1_tag, s2);
   //
   // auto f_s2s0 = Dot(f_s2_tag, s0);  // n方向なので成分はないはず
   auto f_s2s1 = Dot(f_s2_tag, s1);
   auto f_s2s2 = Dot(f_s2_tag, s2);
   //
   // n方向成分の微分は計算できていないが，接線方向微分でわかる
   auto f_s0s0 = -f_s1s1 - f_s2s2;
   return T3Tddd{{{f_s0s0, f_s0s1, f_s0s2},
                  {f_s0s1, f_s1s1, f_s1s2},
                  {f_s0s2, f_s2s1, f_s2s2}}};
};

T3Tddd grad_U_LinearElement(const networkFace *const F) { return grad_U_LinearElement(F, OrthogonalBasis(F->normal)); };

T3Tddd grad_U_LinearElementNeuamnn(const networkPoint *const p, const T3Tddd &orthogonal_basis) {
   /*
   スカラー量の接線方向勾配を計算することはできるが，法線方向はわからない．
   しかし，連続の式を使えば，phiの法線方向の勾配は，接線方向の勾配から計算することができる．
   ∇U=∇∇f={{fxx, fyx, fzx},{fxy, fyy, fzy},{fxz, fyz, fzz}}, ∇∇f=∇∇f^T
   */
   T3Tddd H{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
   double Atot = 0;
   for (const auto &f : p->getFacesNeumann()) {
      H += f->area * grad_U_LinearElement(f, orthogonal_basis);
      Atot += f->area;
   }
   return H / Atot;
};

#endif