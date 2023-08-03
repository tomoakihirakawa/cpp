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

/*DOC_EXTRACT BEM

## 造波装置など

造波板となるobjectに速度を与えることで，造波装置などを模擬することができる．
\ref{BEM:impose_velocity}{強制運動を課す}

\ref{BEM:Hadzic2005}{ここ}では，Hadzic et al. 2005の造波板の動きを模擬している．
角速度の原点は，板の`COM`としている．

*/

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
   } else if (name.contains("Hadzic2005")) {
      // \label{BEM:Hadzic2005}
      const std::vector<Tdd> sample = {
          {-0.2, 0.023474178403755985},
          {-0.1, 0.023474178403755985},
          {0, 0.023474178403755985},
          {0.19200000000000006, 0.023474178403755985},
          {0.2909090909090909, 0.023474178403755985},
          {0.384, -0.023474178403755985},
          {0.4712727272727274, -0.023474178403755985},
          {0.5760000000000002, 0.07042253521126796},
          {0.674909090909091, 0.16431924882629012},
          {0.7679999999999999, 0.11737089201877993},
          {0.861090909090909, -0.07042253521126796},
          {0.930909090909091, -0.11737089201877993},
          {1.024, 0.023474178403755985},
          {1.1636363636363636, 0.352112676056338},
          {1.2625454545454549, 0.16431924882629012},
          {1.3556363636363638, -0.2112676056338021},
          {1.4545454545454546, -0.25821596244131584},
          {1.5650909090909093, 0.2112676056338021},
          {1.6640000000000001, 0.4929577464788739},
          {1.7687272727272725, 0.2112676056338021},
          {1.8385454545454545, -0.25821596244131584},
          {1.9374545454545453, -0.5399061032863859},
          {2.0305454545454547, -0.07042253521126796},
          {2.082909090909091, 0.39906103286385},
          {2.1759999999999997, 0.8215962441314559},
          {2.274909090909091, 0.44600938967136194},
          {2.356363636363636, -0.44600938967136194},
          {2.432, -0.9624413145539918},
          {2.4727272727272727, -1.056338028169014},
          {2.5658181818181816, -0.5399061032863859},
          {2.6356363636363636, 0.30516431924882603},
          {2.7112727272727275, 1.197183098591549},
          {2.775272727272727, 1.338028169014085},
          {2.8625454545454545, 0.633802816901408},
          {2.9265454545454546, -0.352112676056338},
          {2.9963636363636366, -1.2910798122065739},
          {3.060363636363636, -1.666666666666666},
          {3.159272727272727, -1.009389671361502},
          {3.2290909090909086, 0.25821596244131406},
          {3.298909090909091, 1.619718309859155},
          {3.397818181818182, 2.464788732394366},
          {3.4909090909090903, 1.5727699530516421},
          {3.578181818181818, -0.25821596244131584},
          {3.6363636363636367, -1.431924882629108},
          {3.741090909090909, -2.370892018779342},
          {3.8516363636363637, -1.244131455399062},
          {3.9389090909090907, 0.633802816901408},
          {4.002909090909091, 1.948356807511737},
          {4.061090909090909, 2.793427230046949},
          {4.142545454545454, 2.84037558685446},
          {4.194909090909091, 2.183098591549295},
          {4.2647272727272725, 0.727699530516432},
          {4.340363636363636, -1.1502347417840362},
          {4.4043636363636365, -2.652582159624412},
          {4.491636363636363, -4.061032863849764},
          {4.549818181818182, -4.2018779342723},
          {4.637090909090909, -3.215962441314554},
          {4.724363636363636, -1.056338028169014},
          {4.8, 1.150234741784038},
          {4.893090909090909, 3.356807511737089},
          {5.015272727272727, 4.671361502347417},
          {5.143272727272727, 3.826291079812206},
          {5.224727272727272, 2.323943661971831},
          {5.306181818181818, 0.39906103286385},
          {5.41090909090909, -1.807511737089202},
          {5.498181818181818, -3.122065727699532},
          {5.602909090909091, -3.638497652582158},
          {5.748363636363636, -2.699530516431926},
          {5.876363636363637, -1.1971830985915481},
          {6.010181818181818, 0.5399061032863859},
          {6.120727272727272, 1.619718309859155},
          {6.277818181818182, 2.27699530516432},
          {6.429090909090909, 1.854460093896714},
          {6.597818181818181, 0.8685446009389661},
          {6.737454545454545, -0.023474178403755985},
          {6.836363636363636, -0.5868544600938961},
          {6.999272727272727, -0.9624413145539918},
          {7.1738181818181825, -0.8215962441314559},
          {7.325090909090909, -0.5399061032863859},
          {7.546181818181818, -0.023474178403755985},
          {7.761454545454544, 0.352112676056338},
          {7.982545454545454, 0.44600938967136194},
      };
      double start = stod(strings[1] /*start*/);
      const auto intp = InterpolationBspline(3, sample);
      return {0, 0., 0., 0., intp.D(t - start) / 180. * 2 * M_PI, 0.};
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
      // auto ret = optimumVector(V, init, W);
      auto ret = optimumVector(V, init);
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
      // auto ret = optimumVector(V, init, W);
      auto ret = optimumVector(V, init);
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

Tddd grad_phi_tangential(const networkFace *const f) {
   auto [p0, p1, p2] = f->getPoints();
   return gradTangential_LinearElement(Tddd{{std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)}}, T3Tddd{{ToX(p0), ToX(p1), ToX(p2)}});
};

T3Tddd grad_LinearElement(const T3Tddd &F012, const T3Tddd &X012, const T3Tddd &F_n) {
   //! 三角要素の節点の情報変数F0,F1,F2から，三角要素上でのgrad(F)を計算する．
   return {grad_LinearElement(std::get<0>(F012), X012, std::get<0>(F_n)),
           grad_LinearElement(std::get<1>(F012), X012, std::get<1>(F_n)),
           grad_LinearElement(std::get<2>(F012), X012, std::get<2>(F_n))};
};

T3Tddd OrthogonalBasis(const Tddd &n_IN) {
   auto n = Normalize(n_IN);
   Tddd s0 = Chop(Tddd{1, 0, 0}, n);
   if (Norm(s0) < 1E-3)
      s0 = Chop(Tddd{0, 1, 0}, n);
   s0 = Normalize(s0);
   Tddd s1 = Normalize(Cross(n, s0));
   return {n, s0, s1};
};

double getPhin(const networkPoint *p, const networkFace *f) {
   auto iter = p->phinOnFace.find(const_cast<networkFace *>(f));
   if (iter != p->phinOnFace.end())
      return iter->second;
   else
      return p->phinOnFace.at(nullptr);
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT BEM

## 初期値問題

節点の位置と速度ポテンシャル$`\phi`$に関する初期値問題を解いて行くことが，シミュレーションである．
言い換えると，節点位置$`\frac{d\bf x}{dt}`$と速度ポテンシャル$`\frac{d\phi}{dt}`$を少しずつ$`\Delta t`$ずつ時間積分することが，シミュレーションである．
ちなみに，$`\frac{d\bf x}{dt}`$や$`\frac{d\phi}{dt}`$を計算するには，境界値問題を解く必要がある．

ある時刻において，境界値問題が解けたら，$`\frac{d\bf x}{dt}`$と$`\frac{d\phi}{dt}`$はどのように計算できるだろうか．

### 流速$`\frac{d\bf x}{dt}`$の計算

ある三角要素上の接線流速$`\nabla \phi_{\parallel}`$は，線形三角要素補間を使って次のように計算する．

```math
\nabla \phi_{\parallel} = \frac{\bf n}{2A} \times (({\bf x}_2 - {\bf x}_1) \phi_0 +({\bf x}_0 - {\bf x}_2) \phi_1 + ({\bf x}_1 - {\bf x}_0) \phi_2)
```

三角要素上の流速$`\nabla \phi`$は，次のように計算する．

```math
\nabla \phi = \frac{(\phi_n)_0+(\phi_n)_1+(\phi_n)_2}{3} {\bf n} + \nabla \phi_{\parallel}
```

### $`\frac{d\phi}{dt}`$の計算

ある流体粒子に乗ってみたときの，速度ポテンシャルの時間変化$`\frac{D \phi}{D t}`$は，次のように計算できる．

```math
\frac{D \phi}{D t} = \frac{\partial \phi}{\partial t} + \nabla \phi \cdot \nabla \phi
```

<details>
<summary>
NOTE: オイラー的記述
</summary>

$`\phi=\phi(t,{\bf x})`$のように書き表し，位置と空間を独立させ分けて考える方法を，オイラー的記述という．こう書くと，$`\frac{d \phi}{d t}`$は，$`\frac{\partial \phi}{\partial t}`$であり，これは，速度ポテンシャルの純粋な時間変化ではない．純粋な，ある流体粒子の速度ポテンシャルの時間変化を表すためには，位置が時間によって変わると考え，つまり$`\phi=\phi(t,{\bf x}(t))`$と一時的に考えなおし，そして，時間微分する．そうすると$`\frac{d\phi}{dt} = \frac{\partial \phi}{\partial t} + \frac{d\bf x}{dt}\cdot \nabla \phi`$となる．

</details>

ここの$`\frac{\partial \phi}{\partial t}`$の計算は簡単ではない．そこで，ベルヌーイの式（大気圧と接する水面におけるベルヌーイの式は圧力を含まず簡単）を使って，$`\frac{\partial \phi}{\partial t}`$を消去する．

*/

Tddd gradPhi(const networkFace *const f) {
   double phi_n = 0;
   for (const auto &p : f->getPoints())
      phi_n += getPhin(p, f);
   return grad_phi_tangential(f) + phi_n / 3. * f->normal;
};

Tddd gradPhi(const networkPoint *const p) {
   Tddd u;
   V_Tddd V;
   V_d W;
   for (const auto &f : p->getFaces()) {

      u = grad_phi_tangential(f) + getPhin(p, f) * f->normal;

      V.emplace_back(u);
#if defined(use_angle_weigted_normal)
      W.push_back(f->getAngle(p) * (f->Dirichlet ? 10 : 1.));
#elif defined(use_area_weigted_normal)
      W.push_back(f->area * (f->Dirichlet ? 10 : 1.));
#else
      W.push_back(f->Dirichlet ? 10 : 1.);
#endif
   }
   return optimumVector(V, {0., 0., 0.}, W);
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT BEM

## その他

### 境界値問題の未知変数

`isNeumannID_BEM`と`isDirichletID_BEM`は，節点と面の組みが，境界値問題の未知変数かどうかを判定する．
多重節点でない場合は，`{p,nullptr}`が変数のキーとなり，多重節点の場合は，`{p,f}`が変数のキーとなる．

*/

bool isNeumannID_BEM(const auto p, const auto f) {
   if (p->Neumann || p->CORNER) {
      if (p->isMultipleNode) {
         if (p->MemberQ(f))
            return f->Neumann;
         else
            return false;
      } else
         return (f == nullptr);
   } else
      return false;
};

bool isNeumannID_BEM(const std::tuple<netP *, netF *> &PF) { return isNeumannID_BEM(std::get<0>(PF), std::get<1>(PF)); };

bool isDirichletID_BEM(const auto p, const auto f) {
   if (p->Dirichlet || p->CORNER)
      return (f == nullptr);
   else
      return false;
};

bool isDirichletID_BEM(const std::tuple<netP *, netF *> &PF) { return isDirichletID_BEM(std::get<0>(PF), std::get<1>(PF)); };

std::tuple<networkPoint *, networkFace *> pf2ID(const networkPoint *p, const networkFace *f) {
   /**
   NOTE: non-multiple node ID is {p,nullptr}
   NOTE: Iterating over p->getFaces() and p may not get all IDs since p->getFaces() doesn't contain nullptr which is often used for an ID of a non-multiple node.
    */
   if (f == nullptr || !p->isMultipleNode || f->Dirichlet)
      return {const_cast<networkPoint *>(p), nullptr};
   else
      return {const_cast<networkPoint *>(p), const_cast<networkFace *>(f)};
}

std::unordered_set<std::tuple<networkPoint *, networkFace *>> variableIDs(const networkPoint *p) {
   //{p,f}を変換
   // f cannot be nullptr
   //  {p,f} --o--> {p,nullptr}
   //  {p,f} <--x-- {p,nullptr}

   std::unordered_set<std::tuple<networkPoint *, networkFace *>> ret;
   for (const auto &f : p->getFaces())
      ret.emplace(pf2ID(p, f));
   return ret;
};

/*DOC_EXTRACT BEM

### $`\phi_{nt}`$の計算で必要となる$`{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right)`$について．

$`\nabla`$を，$`(x,y,z)`$の座標系ではなく，
面の法線方向$`{\bf n}`$を$`x`$の代わりにとり，
面に水平な方向を$`t_0,t_1`$とする座標系で考えることにして，$`\nabla^*`$と書くことにする．
$`{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right)`$では，$`{\bf n}`$方向成分だけをとる操作をしているので，
新しい座標系でも同じようにすれば，結果は変わらない．

```math
{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right) =  {(1,0,0)}\cdot\left({\nabla^* \phi \cdot \nabla^* \nabla^* \phi}\right).
\quad
\nabla^* \phi = \left(\phi_n, \phi_{t_0}, \phi_{t_1}\right),
\quad \nabla^* \nabla^* \phi =
\begin{bmatrix}
\phi_{nn} & \phi_{nt_0} & \phi_{nt_1} \\
\phi_{t_0n} & \phi_{t_0t_0} & \phi_{t_0t_1} \\
\phi_{t_1n} & \phi_{t_1t_0} & \phi_{t_1t_1}
\end{bmatrix}
```

最後に第１成分だけが残るので，

```math
{(1,0,0)}\cdot\left({\nabla^* \phi \cdot \nabla^* \nabla^* \phi}\right) = \nabla^* \phi \cdot (\phi_{nn}, \phi_{t_0n}, \phi_{t_1n})
```

$`\phi_{nn}`$は，直接計算できないが，ラプラス方程式から$`\phi_{nn}=- \phi_{t_0t_0}- \phi_{t_1t_1}`$となるので，水平方向の勾配の計算から求められる．

*/

/* -------------------------------------------------------------------------- */
// \label{BEM:grad_U_LinearElement}
T3Tddd HessianOfPhi(auto F, const T3Tddd &basis) {
   auto [P0, P1, P2] = F->getPoints();
   auto X012 = T3Tddd{Dot(basis, P0->X), Dot(basis, P1->X), Dot(basis, P2->X)};

   auto [g0_s0, g0_s1, g0_s2] = Dot(basis, P0->U_BEM);
   auto [g1_s0, g1_s1, g1_s2] = Dot(basis, P1->U_BEM);
   auto [g2_s0, g2_s1, g2_s2] = Dot(basis, P2->U_BEM);

   // auto [g_s0s0, g_s0s1, g_s0s2] = gradTangential_LinearElement(Tddd{getPhin(P0, F), getPhin(P1, F), getPhin(P2, F)}, X012);

   auto [g_s0s0, g_s0s1, g_s0s2] = gradTangential_LinearElement(Tddd{g0_s0, g1_s0, g2_s0}, X012);
   auto [g_s1s0, g_s1s1, g_s1s2] = gradTangential_LinearElement(Tddd{g0_s1, g1_s1, g2_s1}, X012);
   auto [g_s2s0, g_s2s1, g_s2s2] = gradTangential_LinearElement(Tddd{g0_s2, g1_s2, g2_s2}, X012);

   // return T3Tddd{{{-g_s1s1 - g_s2s2, g_s0s1 /*wont be used*/, g_s0s2 /*wont be used*/},
   //                {g_s1s0, g_s1s1 /*wont be used*/, g_s1s2 /*wont be used*/},
   //                {g_s2s0, g_s2s1 /*wont be used*/, g_s2s2 /*wont be used*/}}};

   return T3Tddd{{{-g_s1s1 - g_s2s2, g_s0s1 /*wont be used*/, g_s0s2 /*wont be used*/},
                  {(g_s0s1 + g_s1s0) * 0.5, g_s1s1 /*wont be used*/, g_s1s2 /*wont be used*/},
                  {(g_s0s2 + g_s2s0) * 0.5, g_s2s1 /*wont be used*/, g_s2s2 /*wont be used*/}}};
};

double phint_Neumann(networkFace *F) {
   auto Omega = (NearestContactFace(F)->getNetwork())->velocityRotational();
   auto grad_phi = gradPhi(F);
   auto U_body = uNeumann(F);
   auto ret = Dot(Cross(Omega, F->normal), U_body - grad_phi);
   ret += Dot(F->normal, accelNeumann(F));
   auto basis = OrthogonalBasis(F->normal);
   ret -= Dot(Dot(basis, F->normal) /*=(1,0,0)*/, Dot(Dot(basis, U_body), HessianOfPhi(F, basis)));
   // ret -= Dot(Dot(basis, F->normal) /*=(1,0,0)*/, Dot(Dot(basis, grad_phi), HessianOfPhi(F, basis)));
   return ret;
};

double phint_Neumann(const networkPoint *const p) {
   double phint_acum = 0, W_acum = 0, w, phint;
   V_d Phin, W;
   for (const auto &f : p->getFacesNeumann()) {
#if defined(use_angle_weigted_normal)
      w = f->getAngle(p);
#elif defined(use_area_weigted_normal)
      w = f->area;
#else
      w = 1.;
#endif

      phint = phint_Neumann(f);
      phint_acum += w * phint;
      W_acum += w;

      Phin.push_back(phint);
      W.push_back(w);
   }
   // return phint_acum / W_acum;
   return optimumValue(Phin, 0., W);
};

#endif