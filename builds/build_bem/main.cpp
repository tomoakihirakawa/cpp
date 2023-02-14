// #define _debugging_

#define BEM
#define use_lapack
int time_step;
double real_time = 0;

#define simulation
#include <filesystem>
#include "Network.hpp"
#include "integrationOfODE.hpp"
#include "kernelFunctions.hpp"
#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"
pvd cpg_pvd("./vtu/bem.pvd");

#include "bem.hpp"
#include "svd.hpp"

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
         }
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

/* -------------------------------------------------------------------------- */

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

// b$ -------------------------------------------------------------------------- */
// b$                    流体面と構造物面間の関係．面積重み平均．                         */
// b$ -------------------------------------------------------------------------- */
/*
隣接 (Adjacent)
接触 (Contact)
*/
netFp NearestContactFace(const networkPoint *const p) { return std::get<1>(Nearest_(p->X, p->getContactFaces())); };
netFp NearestContactFace(const networkFace *const f_IN) {
   std::unordered_set<networkFace *> faces;
   for_each(f_IN->getPoints(), [&](const auto &q) { faces.insert(q->getContactFaces().begin(), q->getContactFaces().end()); });
   return std::get<1>(Nearest_(f_IN->center, faces));
};
netFp NearestContactFace(const networkPoint *const p,
                         const networkFace *const f_normal) {
   Tddd r = {1E+100, 1E+100, 1E+100}, X;
   networkFace *ret = nullptr;
   for (const auto &f_target : BFS_Flattened(p->getContactFaces(), 2))
      if (isInContact(p, f_normal, f_target)) {
         X = Nearest(p->X, ToX(f_target));
         if (Norm(r) >= Norm(X - p->X)) {
            r = X - p->X;
            ret = f_target;
         }
      }
   return ret;
};
std::tuple<netFp, Tddd> NearestContactFace_(const networkPoint *const p,
                                            const networkFace *const f_normal) {
   Tddd r = {1E+100, 1E+100, 1E+100}, X;
   networkFace *ret = nullptr;
   for (const auto &f_target : BFS_Flattened(p->getContactFaces(), 2))
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

std::tuple<Tddd, double> uNeumann_(const networkPoint *const p,
                                   const networkFace *const f_normal) {
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

std::tuple<Tddd, double> accelNeumann_(const networkPoint *const p,
                                       const networkFace *const f_normal) {
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
   for (const auto &f_normal : p->getFaces()) {
      auto [v, w] = uNeumann_(p, f_normal);
      if (w > 1E-20) {
         V.emplace_back(v);
         W.emplace_back(w);
      }
   }
   if (!V.empty()) {
      auto ret = optimumVector_(V, init, W);
      if (isFinite(ret))
         return ret;
   }
   return init;
};

Tddd accelNeumann(const networkPoint *const p) {
   std::vector<Tddd> V;
   std::vector<double> W;
   const Tddd init = {0., 0., 0.};
   for (const auto &f_normal : p->getFaces()) {
      auto [v, w] = accelNeumann_(p, f_normal);
      if (w > 1E-20) {
         V.emplace_back(v);
         W.emplace_back(w);
      }
   }
   if (!V.empty()) {
      auto ret = optimumVector_(V, init, W);
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
//$ --------------------------------------------------------------- */
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
/* ------------------------------------------------------ */
V_d volume({});

auto modify = [](Network &water, const JSON &json1) {
   auto isExsist = [&json1](std::string str) { return (json1().find(str) != json1().end() && !json1()[str].empty()); };
   if (isExsist("center_of_mass")) {
      std::get<0>(water.center_of_mass) = stob(json1()["center_of_mass"])[0];
      std::get<1>(water.center_of_mass) = stob(json1()["center_of_mass"])[1];
      std::get<2>(water.center_of_mass) = stob(json1()["center_of_mass"])[2];
   }
   if (isExsist("ignore")) {
      water.IGNORE = stob(json1()["ignore"])[0];
   }
   if (isExsist("rotate")) {
      auto rotate = stod(json1()["rotate"]);
      if (rotate.size() > 1)
         water.rotate(rotate[0], Tddd{rotate[1], rotate[2], rotate[3]});
   }
   if (isExsist("scale")) {
      auto scale = stod(json1()["scale"]);
      if (scale.size() > 1)
         water.scale({scale[0], scale[1], scale[2]});
      else
         water.scale(scale[0]);
   }
   if (isExsist("translate")) {
      auto translate = stod(json1()["translate"]);
      if (translate.size() > 1)
         water.translate({translate[0], translate[1], translate[2]});
   }
   // if (isExsist("remesh"))
   // {
   // 	auto minlen = stod(json1()["remesh"]);
   // 	if (minlen.size() > 0)
   // 		remesh(&water, minlen[0]);
   // }
   // if (isExsist("coarsen"))
   // {
   // 	auto minlen = stod(json1()["coarsen"]);
   // 	if (minlen.size() > 0)
   // 		coarsen(&water, minlen[0]);
   // }
   if (isExsist("reverseNormal")) {
      std::string TorF = json1()["reverseNormal"][0];
      if (TorF.compare("True") == 0 || TorF.compare("true") == 0 || TorF.compare("1") == 0) {
         water.reverseNormal();
         std::cout << "reverse done" << std::endl;
      }
   }
};

//@ ------------------------------------------------------ */

Tddd gradTangential_LinearElement(const Tddd &phi012, const T3Tddd &X012) {
   auto [X0, X1, X2] = X012;
   auto [phi0, phi1, phi2] = phi012;
   auto n = TriangleNormal(X012);
   return Cross(n, phi0 * (X2 - X1) + phi1 * (X0 - X2) + phi2 * (X1 - X0)) / (2 * TriangleArea(X012));
};

T3Tddd gradTangential_LinearElement(const T3Tddd &V012, const T3Tddd &X012) {
   auto [X0, X1, X2] = X012;
   auto [Vx012, Vy012, Vz012] = Transpose(V012);
   return {gradTangential_LinearElement(Vx012, X012),
           gradTangential_LinearElement(Vy012, X012),
           gradTangential_LinearElement(Vz012, X012)};
};

T3Tddd grad_U_tangential_LinearElement(const networkFace *const f) {
   auto [p0, p1, p2] = f->getPoints();
   return gradTangential_LinearElement({p0->U_BEM, p1->U_BEM, p2->U_BEM}, ToX(f));
};

T3Tddd grad_U_LinearElement(const networkPoint *const p) {
   T3Tddd gradUtang = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
   /*
   スカラー量の接線方向勾配を計算することはできるが，法線方向はわからない．
   しかし，連続の式を使えば，phiの法線方向の勾配は，接線方向の勾配から計算することができる．
   */
   double Atot = 0;
   for (const auto &f : p->getFaces()) {
      gradUtang += f->area * grad_U_tangential_LinearElement(f);
      Atot += f->area;
   }
   gradUtang /= Atot;  // 接線方向勾配
   auto q = (p->getNeighbors())[0];
   auto V = q->getXtuple() - ToX(p);
   auto n = p->getNormal_BEM();
   Tddd s0 = Normalize(V - n * Dot(V, n));
   Tddd s1 = Normalize(Cross(n, s0));
   Tddd V0 = Dot(gradUtang, s0);
   Tddd V1 = Dot(gradUtang, s1);
   Tddd V2 = {std::get<2>(V0), std::get<2>(V1), -std::get<0>(V0) - std::get<1>(V1)};
   return T3Tddd{V0, V1, V2};
};

T3Tddd grad_U_LinearElement(const networkFace *const f) {
   T3Tddd gradUtang = grad_U_tangential_LinearElement(f);
   auto [p0, p1, p2] = f->getPoints();
   auto V = ToX(p1) - ToX(p0);
   auto n = f->normal;
   Tddd s0 = Normalize(V - n * Dot(V, n));
   Tddd s1 = Normalize(Cross(n, s0));
   Tddd V0 = Dot(gradUtang, s0);
   Tddd V1 = Dot(gradUtang, s1);
   Tddd V2 = {std::get<2>(V0), std::get<2>(V1), -std::get<0>(V0) - std::get<1>(V1)};
   return T3Tddd{V0, V1, V2};
};

   // b@ ------------------------------------------------------ */

   // bool isFacing(const Tddd &F, const Tddd &f, double rad = 1E-10) {
   //    return (Dot(F, -f) >= cos(rad));
   // };

#include "BEM_derivatives.hpp"

//* ------------------------------------------------------ */
//*                        境界値問題を解く                   */
//* ------------------------------------------------------ */
// #define solve_equations_on_all_points
#define solve_equations_on_all_points_rigid_mode
#define solveBVP_debug
#include "BEM_BVP.hpp"

// b! ------------------------------------------------------ */
// b!           格子のdivide, merge．それに伴うΦ，Φnの付与       */
// b! ------------------------------------------------------ */

Tdd estimate_phiphin(const networkLine *const l) {
   // auto fs = l->getFaces();
   // interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_0(fs[0], l);
   // interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_1(fs[1], l);
   // auto phi0 = Dot(intp_l0_0.N(.5, .5), ToPhi(intp_l0_0.Points));
   // auto phi1 = Dot(intp_l0_1.N(.5, .5), ToPhi(intp_l0_1.Points));
   // auto phin0 = Dot(intp_l0_0.N(.5, .5), ToPhin(intp_l0_0.Points));
   // auto phin1 = Dot(intp_l0_1.N(.5, .5), ToPhin(intp_l0_1.Points));
   // return {(phi0 + phi1) / 2., (phin0 + phin1) / 2.};
   //
   auto [a, b] = l->getPoints();
   return (a->phiphin + b->phiphin) / 2.;
};

/* ------------------------------------------------------ */

void remesh(Network &water,
            const Tdd &limit_angle_D,
            const Tdd &limit_angle_N,
            bool force = false,
            int max_count = 100) {
   std::cout << "remeshing" << std::endl;
   water.setGeometricProperties();
   double mean_length = Mean(extLength(water.getLines()));
   bool isfound = false, ismerged = false;
   int count = 0;
   // double lim_degree_Neumann = limit_angle;
   // double lim_degree = limit_angle;

   networkLine *l;
   Tddd X, V;
   networkPoint *q;
   double meanArea;
   V_netFp Fs;
   Tdd phiphin;
   V_netLp lines;
   double local_mean_length;
   do {
      // なくなるまでやるか？
      isfound = false;
      ismerged = false;
      /* ------------------------------------------------------ */
      for (const auto &p : RandomSample(ToVector(water.getPoints()))) {
         /* ------------------------------------------------------ */
         meanArea = Mean(p->getFaceAreas());
         if (false)
            for (const auto &f : p->getFaces()) {
               if (f->area / meanArea < 1E-3) {
                  p->sortLinesByLength();
                  l = *(p->getLines().rbegin());
                  phiphin = estimate_phiphin(l);
                  auto [a, b] = l->getPoints();
                  X = (a->getXtuple() + b->getXtuple()) / 2.;
                  q = l->divide();
                  q->phiphin = phiphin;
                  q->setX(X);
                  isfound = true;
                  break;
               }
            }
         if (isfound)
            break;
         /* ------------------------------------------------------ */
         //* ------------------------------------------------------ */
         //*                立体角が小さすぎる場合merge                 */
         //* ------------------------------------------------------ */
         if (false)
            if (p->getSolidAngle() < 4 * M_PI / 100. || (4 * M_PI - p->getSolidAngle()) < 4 * M_PI / 100.) {
               p->sortLinesByLength();
               l = p->getLines()[0];
               //@ case2 lの面全体を考慮
               // Tdd phiphin = phiphin_from_faces(l);
               //@ case3 lの点だけを考慮
               // Tdd phiphin = phiphin_from_points(l);
               //@ case4 lの面全体を考慮
               phiphin = estimate_phiphin(l);
               /* ------------------------------------------------------ */
               // auto [a, b] = l->getPoints();
               // if (a->CORNER)
               // 	b = a;
               // else if (b->CORNER)
               // 	a = b;
               // Tddd V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
               // for (const auto &f : a->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               // for (const auto &f : b->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               // Tddd X = V + a->getXtuple();
               // auto q = l->merge();
               // q->phiphin = phiphin;
               // q->setX(X);
               // ismerged = true;
               // break;
               /* ------------------------------------------------------ */
               auto [a, b] = l->getPoints(p);
               X = b->getXtuple();
               q = l->merge();
               q->phiphin = phiphin;
               q->setX(X);
               ismerged = true;
               break;
            }
         //! ------------------------------------------------------ */
         //!             辺の長さが長すぎるまたは短すぎる場合             */
         //! ------------------------------------------------------ */
         local_mean_length = Mean(extLength(extractLines(Flatten(BFS(p, 2)))));
         lines = p->getLines();
         sortByLength(lines);
         for (const auto &l : Reverse(lines)) {
            auto [p0, p1] = l->getPoints();
            Fs = l->getFaces();
            //@ ------------------------------------------------------ */
            // if (l->length() > mean_length * 1.75 /*長すぎる*/)
            if (l->length() > local_mean_length * 1.5 /*長すぎる*/) {
               //@ case2 lの面全体を考慮
               // Tdd phiphin = phiphin_from_faces(l);
               //@ case3 lの点だけを考慮
               // Tdd phiphin = phiphin_from_points(l);
               //@ case4 lの面全体を考慮
               phiphin = estimate_phiphin(l);
               /* ------------------------------------------------------ */
               q = l->divide();
               q->phiphin = phiphin;
               isfound = true;
               break;
            }
            //@ ------------------------------------------------------ */
            if (l->length() < local_mean_length / 20.) {
               /*
               b!マージの原則：マージによってノイマン面は変形してはいけない．
               */
               //@ case5 ２点の平均，移動位置は，ノイマンを崩さない方向：ノイマン面の法線方向成分には移動しない．
               auto [a, b] = l->getPoints();
               if (!(a->CORNER && b->CORNER)) {
                  if (a->CORNER)
                     b = a;
                  else if (b->CORNER)
                     a = b;
               }
               phiphin = (a->phiphin + b->phiphin) / 2.;
               V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
               V -= a->getNormalTuple() * Dot(V, a->getNormalTuple());
               V -= b->getNormalTuple() * Dot(V, b->getNormalTuple());
               // for (const auto &f : a->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               // for (const auto &f : b->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               X = V + a->getXtuple();
               q = l->merge();
               q->phiphin = phiphin;
               q->setX(X);
               /* ------------------------------------------------------ */
               ismerged = true;
               break;
            }
            //@ ------------------------------------------------------ */
            auto [min, max] = MinMax(extAreas(Fs));
            if (false && min / max < 1 / 100.) {
               //@ case5 ２点の平均，移動位置は，ノイマンを崩さない方向：ノイマン面の法線方向成分には移動しない．
               auto [a, b] = l->getPoints();
               if (!(a->CORNER && b->CORNER)) {
                  if (a->CORNER)
                     b = a;
                  else if (b->CORNER)
                     a = b;
               }
               phiphin = (a->phiphin + b->phiphin) / 2.;
               V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
               V -= a->getNormalTuple() * Dot(V, a->getNormalTuple());
               V -= b->getNormalTuple() * Dot(V, b->getNormalTuple());
               // for (const auto &f : a->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               // for (const auto &f : b->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               X = V + a->getXtuple();
               q = l->merge();
               q->phiphin = phiphin;
               q->setX(X);
               /* ------------------------------------------------------ */
               ismerged = true;
               break;
            }
            if (ismerged || isfound)
               break;
         }

         // for (const auto &l : water.getLines())
         // {
         // 	auto [p0, p1] = l->getPoints();
         // 	Fs = l->getFaces();
         // 	if (!l->CORNER)
         // 		if (force)
         // 		{
         // 			isfound = l->flipIfTopologicalyBetter((l->Neumann ? lim_degree_Neumann : lim_degree));
         // 			break;
         // 		}
         // 		else
         // 		{
         // 			isfound = l->flipIfBetter((l->Neumann ? lim_degree_Neumann : lim_degree));
         // 		}
         // }

         if (ismerged || isfound)
            break;
      }
   } while ((ismerged || isfound) && count++ < max_count);

   flipIf(water, limit_angle_D, limit_angle_N, force);
};
// b! ------------------------------------------------------ */
// b! ------------------------------------------------------ */
// b! ------------------------------------------------------ */

void setBoundaryConditions(Network &water, const std::vector<Network *> &objects) {
   auto radius = Mean(extLength(water.getLines()));
   auto Points = ToVector(water.getPoints());
   Print("makeBucketFaces", Green);
   for (const auto &net : objects) {
      radius = Mean(extLength(net->getLines()));
      net->makeBucketFaces(radius);
   }

   // b% -------------------------------------------------------- */
   // b%              境界条件（角点・ディリクレ・ノイマン）の決定              */
   // b% -------------------------------------------------------- */
   // b% step1 点の衝突の判定
   std::cout << "step1 点の衝突の判定" << std::endl;
   for (const auto &p : Points)
      p->clearContactFaces();
   //!!! 衝突の判定がよくエラーが出る箇所
   for (const auto &net : objects) {
#pragma omp parallel
      for (const auto &p : Points)
#pragma omp single nowait
      {
         //! ここも重要：点と面の衝突をどのようにすれば矛盾なく判定できるか．
         p->radius = Mean(extLength(p->getLines())) / 3.;
         // auto toF = extXtuple(ToVector(p->getFaces())) - ToX(p);
         // auto toP = extXtuple(p->getNeighbors()) - ToX(p);
         // double a = Norm(*std::min_element(toP.begin(), toP.end(), [](const auto &a, const auto &b) { return Norm(a) < Norm(b); }));
         // double b = Norm(*std::min_element(toF.begin(), toF.end(), [](const auto &a, const auto &b) { return Norm(a) < Norm(b); }));
         // p->radius = Mean(extLength(extractLines(Flatten(BFS(p, 2))))) / 5.;
         p->addContactFaces(net->getBucketFaces(), false); /**shadowあり*/
      }
   }
   // b% step2 面の境界条件を判定
   std::cout << "step2 面の境界条件を判定" << std::endl;

   // auto isNeumann = [&](const networkFace *const f) {
   //    auto [p0, p1, p2] = f->getPoints();
   //    auto faces_p0 = p0->getContactFaces();
   //    auto faces_p1 = p1->getContactFaces();
   //    auto faces_p2 = p2->getContactFaces();
   //    bool isNeumann = f->isThereAnyFacingFace(faces_p0, M_PI / 9.) &&
   //                     f->isThereAnyFacingFace(faces_p1, M_PI / 9.) &&
   //                     f->isThereAnyFacingFace(faces_p2, M_PI / 9.);
   //    return isNeumann;
   // };

   /*面Aの点が接触している面Bを取得．A,B面が向き合っていればノイマン*/
   for (const auto &f : water.getFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      f->Neumann = isInContact(p0, f, BFS_Flattened(p0->getContactFaces(), 2)) &&
                   isInContact(p1, f, BFS_Flattened(p1->getContactFaces(), 2)) &&
                   isInContact(p2, f, BFS_Flattened(p2->getContactFaces(), 2));
      f->Dirichlet = !f->Neumann;
   }
   // b% step3 線の境界条件を決定
   std::cout << "step3 線の境界条件を決定" << std::endl;
   for (const auto &l : water.getLines()) {
      auto faces = l->getFaces();
      l->Neumann = std::all_of(faces.begin(), faces.end(), [](const auto &f) { return f->Neumann; });
      l->Dirichlet = std::all_of(faces.begin(), faces.end(), [](const auto &f) { return f->Dirichlet; });
      if (!l->Neumann && !l->Dirichlet)
         l->CORNER = true;
      else
         l->CORNER = false;
   }
   // b% step4 点の境界条件を決定
   std::cout << "step4 点の境界条件を決定" << std::endl;
   /*
     周りの面が全てノイマンなら点はノイマン．
     周りの面がディリクレとノイマン両方を持っていれば角点
   　それ以外は，ディリクレ．
   */
   for (const auto &p : Points) {
      auto faces = p->getFacesUO();
      bool isNeumann = std::all_of(faces.begin(), faces.end(), [](const auto &f) { return f->Neumann; }) &&
                       !p->getContactFacesX().empty();
      // getContactFacesXがないと，周りの面がノイマンでも，phinを計算できないことがある．
      bool isDirichlet = std::all_of(faces.begin(), faces.end(), [](const auto &f) { return f->Dirichlet; });

      if (isNeumann)
         p->setN();
      else if (isDirichlet)
         p->setD();
      else
         p->setC();
   }
   // b! ------------------------------------------------------ */
   // b!       　    ノイマン境界の点や面にはΦnを与える                  */
   // b! ------------------------------------------------------ */
   // b! 点
   std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << colorOff << std::endl;
   //@ -------------------------------------------------------------------------------- */
   //@ -------------------------------- 多重節点となる条件 ------------------------------- */
   //@ -------------------------------------------------------------------------------- */
   // b* p->phinOnFaceは，std::unordered_map<networkFace *, double>
   // b* 節点のphinを保存する．また，多重節点かどうかも，面がnullptrかどうかで判別できる．
   auto multiple_node_if = [&](const auto &p, const auto &facesNeuman) {
      // return p->Neumann;
      return (p->CORNER || std::any_of(facesNeuman.begin(), facesNeuman.end(), [&](const auto &f) { return !isFlat(p->getNormalNeumann_BEM(), f->normal, M_PI / 180. * 20); }));
      // return (p->CORNER && std::any_of(facesNeuman.begin(), facesNeuman.end(), [&](const auto &f) { return !isFlat(p->getNormalNeumann_BEM(), f->normal, M_PI / 180. * 20); }));
   };
   //@ -------------------------------------------------------------------------------- */
   for (const auto &p : Points) {
      p->phinOnFace.clear();
      p->phintOnFace.clear();
      if (p->Neumann || p->CORNER) {
         // 角度に応じて変更する
         auto facesNeuman = p->getFacesNeumann();
         if (multiple_node_if(p, facesNeuman)) {
            for (const auto &f : facesNeuman) {
               p->phinOnFace[f] = Dot(uNeumann(p, f), f->normal);
               p->phintOnFace[f] = 1E+30;  // この値は，derivativesクラス内で計算する
            }
         } else {
            std::get<1>(p->phiphin) = Dot(uNeumann(p), p->getNormalNeumann_BEM());
            p->phinOnFace[nullptr] = std::get<1>(p->phiphin);
            p->phintOnFace[nullptr] = 1E+30;  // この値は，derivativesクラス内で計算
         }

         if (!isFinite(p->phiphin)) {
            std::cout << "p->phiphinはfiniteではない！！" << std::endl;
            if (p->Neumann)
               std::cout << "Neumann" << std::endl;
            if (p->Dirichlet)
               std::cout << "Dirichlet" << std::endl;
            if (p->CORNER)
               std::cout << "CORNER" << std::endl;
            std::cout << "p->phiphin = " << p->phiphin << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         }
      }
   }
   // b! 面
   std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << colorOff << std::endl;
   for (const auto &f : water.getFaces()) {
      if (f->Neumann) {
         std::get<1>(f->phiphin) = Dot(uNeumann(f), f->normal);
      } else {
         auto [p0, p1, p2] = f->getPoints();
         std::get<0>(f->phiphin) = (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3.;
      }
   }
};

// b* ------------------------- 出力 ------------------------- */
VV_VarForOutput dataForOutput(const Network &water, const double dt) {
   auto hist = Histogram(extLength(water.getLines()));
   std::stringstream ss;
   ss << "\"cumulative_count\":" << hist.cumulative_count << ","
      << "\"diff\":" << hist.diff << ","
      << "\"count\":" << hist.count << ","
      << "\"interval\":" << hist.interval << ","
      << "\"bin_width\":" << hist.bin_width << ","
      << "\"mid_interval\":" << hist.mid_interval << ",";
   std::string s = ss.str();
   std::replace(s.begin(), s.end(), '{', '[');
   std::replace(s.begin(), s.end(), '}', ']');
   std::cout << "{" << s << "}" << std::endl;
   // if (*hist.data.rbegin() > 5)
   // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "length > 5");

   int ii = 0;
   derivatives ders(water, false);
   //-------------------------------------------
   map_P_Vd P_phiphin;
   // V_d lim_len = Subdivide(0.3, 0.2, 5 - 1);
   /* ------------------------------------------------------ */
   std::cout << "-------------------- 次の時刻の変数の値を得る -------------------- " << std::endl;
   std::cout << "------------------ getImprovedの後に微分を評価 ----------------------- " << std::endl;
   uomap_P_Tddd P_accel_body, P_NearestContactFacesX, P_position,
       P_phin_Neumann, P_phin_Dirichlet, P_velocity_body,
       P_uNeumann, P_normal, P_normal_BEM, P_mirrorPosition, P_U_normal_BEM, P_U_dot_gradgrad_U,
       P_U_tangential_BEM, P_U_BEM, P_adustment_vector, P_U_update_BEM, P_U_cling_to_Neumann, P_phin_vector, P_gradPhi_tangential, P_gradPhi;
   uomap_P_d P_IG, P_IGn, P_isGoodForQuad, P_smin_min, P_s_m, P_minViewRatio, P_DphiDt,
       P_volume, P_phi_Neumann, P_phi_Dirichlet, P_state, P_solidangleBIE, P_height, P_phi, P_phin,
       P_face_size, P_radius, P_lines_size, P_ishit, P_BC, P_Intxn_size, P_ContactFaces,
       P_is_multiple_phiphin, P_min_depth, P_aphiat, P_aphiant, P_pressure, P_update_vs_cling, P_normalVariance, P_isMultipleNode;

   uomap_P_Tddd initial_uomap_P_Tddd;
   uomap_P_d initial_uomap_P_d;
   for (const auto &p : water.getPoints()) {
      initial_uomap_P_Tddd[p] = {1E+30, 1E+30, 1E+30};
      initial_uomap_P_d[p] = 1E+30;
   }
   P_U_dot_gradgrad_U = P_gradPhi = P_gradPhi_tangential = P_phin_vector = P_accel_body = P_position = P_NearestContactFacesX = P_phin_Neumann = P_phin_Dirichlet = P_velocity_body = P_uNeumann = P_normal = P_normal_BEM = P_mirrorPosition = P_U_normal_BEM = P_U_tangential_BEM = P_adustment_vector = P_U_BEM = P_U_update_BEM = P_U_cling_to_Neumann = initial_uomap_P_Tddd;
   P_isMultipleNode = P_DphiDt = P_IG = P_IGn = P_isGoodForQuad = P_smin_min = P_s_m = P_minViewRatio = P_volume = P_phi_Neumann = P_phi_Dirichlet = P_state = P_solidangleBIE = P_height = P_phi = P_phin = P_face_size = P_radius = P_lines_size = P_ishit = P_BC = P_Intxn_size = P_ContactFaces = P_is_multiple_phiphin = P_min_depth = P_aphiat = P_aphiant = P_pressure = P_update_vs_cling = P_normalVariance = initial_uomap_P_d;

   Print("ders.P_phiphin_InnerOuterCornerPを出力");
   try {
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (const auto &p : water.getPoints())
#ifdef _OPENMP
#pragma omp single nowait
#endif
      {
         if (p->Neumann || p->CORNER) {
            P_accel_body[p] = accelNeumann(p);
            P_velocity_body[p] = uNeumann(p);
         }
         // auto [m0, s0, min0, smin0, max0] = distorsion(p, dt);
         // P_s_m[p] = s0 / m0;
         // P_smin_min[p] = smin0 / min0;
         P_volume[p] = water.getVolume();
         // P_phin_Neumann[p] = p->getNormalNeumann_BEM() * p->phin_Neumann;
         P_phin_Dirichlet[p] = p->getNormalDirichlet_BEM() * p->phin_Dirichlet;
         // P_phi_Neumann[p] = p->phi_Neumann;
         P_phi_Dirichlet[p] = p->phi_Dirichlet;
         if (p->Neumann || p->CORNER) {
            auto f = NearestContactFace(p);
            if (f) {
               P_NearestContactFacesX[p] = Nearest(p->X, NearestContactFace(p)) - ToX(p);
            }
         }
         P_adustment_vector[p] = p->U_BUFFER;
         P_U_cling_to_Neumann[p] = p->U_cling_to_Neumann;
         P_update_vs_cling[p] = Norm(p->U_cling_to_Neumann) / Norm(p->U_update_BEM);
         P_U_BEM[p] = p->U_BEM;
         P_U_update_BEM[p] = p->U_update_BEM;
         // P_solidangleBIE[p] = p->getSolidAngle();
         // P_minViewRatio[p] = minViewRatio(p);
         if (p->phinOnFace.empty())
            P_isMultipleNode[p] = 0;
         else if (p->phinOnFace.find(nullptr) != p->phinOnFace.end())
            P_isMultipleNode[p] = 1;
         else
            P_isMultipleNode[p] = 2;
         // P_normalVariance[p] = normalVariance(p);
         P_U_normal_BEM[p] = p->U_normal_BEM;
         // P_U_tangential_BEM[p] = p->U_tangential_BEM;
         // P_state[p] = p->getStatus();
         // P_height[p] = std::get<2>(p->X);
         P_phi[p] = std::get<0>(p->phiphin);
         P_phin[p] = std::get<1>(p->phiphin);
         // P_normal[p] = p->normal;
         P_normal_BEM[p] = p->getNormal_BEM();
         // P_face_size[p] = (double)p->getFaces().size();
         // P_lines_size[p] = (double)p->getLines().size();
         // P_lines_length[p] = extLength(p->getLines());
         // P_Intxn_size[p] = (double)takeIntxn(p->getLines()).size();
         P_ContactFaces[p] = (double)p->getContactFaces().size();
         // P_Intxn_length[p] = extLength(takeIntxn(p->getLines()));
         P_BC[p] = p->Dirichlet ? 0. : (p->Neumann ? 1. : (p->CORNER ? 2. : 1 / 0.));
         // if (!p->getContactFaces().empty())
         // 	P_mirrorPosition[p] = 2. * ((*p->getContactFaces().begin()).second) - ToX(p);
         P_radius[p] = p->radius;
         P_position[p] = ToX(p);
         // P_ishit[p] = (double)(!p->getContactFaces().empty());
         // P_ishit[p] = (double)(p->getStatus());
         // P_is_multiple_phiphin[p] = (double)(p->multiple_phiphin.size());
         // P_min_depth[p] = p->minDepthFromCORNER; //;getMinDepth(p);
         bool isgood = true;
         for (const auto &l : p->getLines())
            isgood = isgood && l->isGoodForQuadInterp();
         // P_isGoodForQuad[p] = isgood;
         // P_aphiat[p] = std::get<0>(p->phiphin_t);
         // P_aphiant[p] = std::get<1>(p->phiphin_t);
         P_pressure[p] = p->pressure_BEM;
         P_uNeumann[p] = uNeumann(p);
         P_DphiDt[p] = p->DphiDt(p->U_update_BEM, 0.);
         P_phin_vector[p] = p->U_normal_BEM;
         // P_gradPhi_tangential[p] = p->U_tangential_BEM;
         P_gradPhi[p] = p->U_BEM;
         P_U_dot_gradgrad_U[p] = Dot(p->U_BEM, grad_U_LinearElement(p));
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
   try {
      // phiのアップデートがされていない
      VV_VarForOutput data = {
          // {"s/m", P_s_m},
          // {"smin/min", P_smin_min},
          {"accel Neumann", P_accel_body},
          {"adustment vector", P_adustment_vector},
          // {"volume", P_volume},
          {"U_cling_to_Neumann", P_U_cling_to_Neumann},
          //  {"absU_cling_to_Neumann/absU_update_BEM", P_update_vs_cling},
          {"velocity Neumann", P_velocity_body},
          {"Nearest face", P_NearestContactFacesX},
          //  {"φn_Neumann", P_phin_Neumann},
          //  {"φn_Dirichlet", P_phin_Dirichlet},
          //  {"φ_Neumann", P_phi_Neumann},
          //  {"φ_Dirichlet", P_phi_Dirichlet},
          //  {"min_depth", P_min_depth},
          {"isMultipleNode", P_isMultipleNode},
          // P_NearestContactFacesXで確かに最寄の構造物までをさすことができている．！
          // これで改善できない？
          //  {"U_update_BEM", P_U_update_BEM},
          {"U_BEM", P_U_BEM},
          //  {"U_normal_BEM", P_U_normal_BEM},
          {"U_tangential_BEM", P_U_tangential_BEM},
          // {"kappa", ders.P_kappa},
          {"ContactFaces", P_ContactFaces},
          // {"dxdt_mod", P_dxdt_mod},
          {"grad_phi", P_gradPhi},
          //  {"phin_vector", P_phin_vector},
          {"gradPhiTangential", P_gradPhi_tangential},
          //  {"z", P_height},
          {"position", P_position},
          {"φ", P_phi},
          {"φn", P_phin},
          //  {"solidangle", P_solidangleBIE},
          //  {"minViewRatio", P_minViewRatio},
          //  {"normalVariance", P_normalVariance},
          //  {"normal", P_normal},
          //  {"uNeumann", P_uNeumann},
          //  {"normal_BEM", P_normal_BEM},
          //  {"all_lines_are_GoodForQuad", P_isGoodForQuad},
          // {"laplacian_phi", P_laplacian},
          // {"face_size", P_face_size},
          // {"line_length", 10, P_lines_length},
          // {"line_size", P_lines_size},
          // {"Intxn_size", P_Intxn_size},
          // {"Intxn_length", 10, P_Intxn_length},
          {"boundary condition", P_BC},
          //  {"radius", P_radius},
          // {"state", P_state},
          // {"correction vector", ders.P_dxdt_correct},
          //  {"DφDt", P_DphiDt},
          //  {"φt", P_aphiat},
          //  {"φnt", P_aphiant},
          {"pressure", P_pressure}
          //  {"U.∇U", P_U_dot_gradgrad_U},
          //  {"IG", P_IG},
          //  {"IGn", P_IGn}
          // {"vector to mirrorPosition", P_mirrorPosition},
          // {"is hit", P_ishit}
      };
      return data;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
   std::cout << __PRETTY_FUNCTION__ << " done" << std::endl;
   // mk_vtu(output_directory + "/" + net.getName() + std::to_string(time_step) + ".vtu", net.getFaces(), datacpg);
   // mk_vtu(output_directory + "/" + name + std::to_string(time_step) + ".vtu", cpg.well_Faces, datacpg);
   // cpg_pvd.push(name, name + std::to_string(time_step) + ".vtu", time_step * t_rep * dt);
   // cpg_pvd.output_();
   return {};
};
/* ------------------------------------------------------ */

double dt_CFL(const Network &water, double min_dt, const Tdd coeff = {0.4, 1}) {
   auto [c0, c1] = coeff;
   for (const auto &p : water.getPoints()) {
      for (const auto &q : p->getNeighbors()) {
         if (min_dt > c0 * Norm(ToX(p) - q->getXtuple()) / Norm(p->U_update_BEM - q->U_update_BEM)) {
            min_dt = c0 * Norm(ToX(p) - q->getXtuple()) / Norm(p->U_update_BEM - q->U_update_BEM);
         }
         if (min_dt > c1 * Norm(ToX(p) - q->getXtuple()) / Norm(p->U_update_BEM)) {
            min_dt = c1 * Norm(ToX(p) - q->getXtuple()) / Norm(p->U_update_BEM);
         }
      }
   }
   return min_dt;
};

void show_info(const Network &net) {
   int total = 0, total_c_face = 0, c = 0, n = 0, d = 0;
   for (const auto &p : net.getPoints()) {
      total++;
      if (p->CORNER) {
         c++;
         total_c_face += p->getFaces().size();
      } else if (p->Neumann)
         n++;
      else if (p->Dirichlet)
         d++;
   }
   std::cout << "net.getPoints() = " << net.getPoints().size() << std::endl;
   std::cout << "Total : " << total << std::endl;
   std::cout << "Total variables: " << d + n + total_c_face << std::endl;
   int doublenode = total - c + total_c_face;
   std::cout << "Total case double-node : " << doublenode << std::endl;
   std::cout << "node reduction : " << (double)(doublenode - total) / (double)doublenode << std::endl;
   std::cout << "CORNER : " << c << std::endl;
   std::cout << "Total CORNER faces : " << total_c_face << std::endl;
   std::cout << "Neumann : " << n << std::endl;
   std::cout << "Dirichlet : " << d << std::endl;
};
JSONoutput jsonout;

// b! ------------------------------------------------------ */

struct calculateFroudeKrylovForce {
   std::vector<networkFace *> actingFaces;
   Tddd force, torque;
   double area;
   T6d acceleration;
   std::vector<std::tuple<Tddd, T3Tddd>> PressureVeticies;
   calculateFroudeKrylovForce(const std::unordered_set<networkFace *> faces /*waterfaces*/,
                              const Network *PasObj)
       : force({0., 0., 0.}),
         torque({0., 0., 0.}),
         area(0.),
         PressureVeticies({}),
         acceleration({0., 0., 0., 0., 0., 0.}) {
      // PasObjと接したfaceの頂点にpressureが設定されている前提
      int count = 0;
      for (const auto &f : faces)
         if (f->Neumann) {
            auto [p0, p1, p2] = f->getPoints();
            if (std::any_of(p0->getContactFaces().begin(), p0->getContactFaces().end(), [&](const auto &F) { return F->getNetwork() == PasObj; }) &&
                std::any_of(p1->getContactFaces().begin(), p1->getContactFaces().end(), [&](const auto &F) { return F->getNetwork() == PasObj; }) &&
                std::any_of(p2->getContactFaces().begin(), p2->getContactFaces().end(), [&](const auto &F) { return F->getNetwork() == PasObj; })) {
               this->PressureVeticies.push_back({{p0->pressure, p1->pressure, p2->pressure}, ToX(f)});
               this->actingFaces.emplace_back(f);
               count++;
            }
         }
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpX = interpolationTriangleLinear0101(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            area += intpX.J(x0, x1) * w0w1;
      }
      std::cout << "接触している面の数:" << count << std::endl;
      std::cout << "表面積:" << area << std::endl;
   };

   Tddd getFroudeKrylovTorque(const Tddd &COM) {
      /*
      crossの引数の順番に注意
      モーメントの計算が，N=rxP
      */
      this->torque = {0., 0., 0.};
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpP = interpolationTriangleLinear0101(P012);
         auto intpX = interpolationTriangleLinear0101(X012);
         auto n = TriangleNormal(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            this->torque += Cross(intpX(x0, x1) - COM, n * intpP(x0, x1)) * intpX.J(x0, x1) * w0w1;
      }
      return this->torque;
   };

   Tddd surfaceIntegralOfPressure() {
      this->force = {0., 0., 0.};
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpP = interpolationTriangleLinear0101(P012);
         auto intpX = interpolationTriangleLinear0101(X012);
         auto n = TriangleNormal(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            this->force += n * intpP(x0, x1) * intpX.J(x0, x1) * w0w1;
      }
      return this->force;
   };
};

// b! ------------------------------------------------------ */

struct outputInfo {
   std::string pvd_file_name;
   std::string vtu_file_name;
   PVDWriter *PVD;
   outputInfo(){};
};

// b! ------------------------------------------------------ */

int main(int arg, char **argv) {
   if (arg <= 1)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "argv <= 1. write input json file directory!\nex.\n$ ./main ./input");
   std::string input_directory(argv[1]);  // input directory
   input_directory += "/";
   std::cout << input_directory << std::endl;
   try {
      //* ------------------------------------------------------ */
      //*                         setting                        */
      //* ------------------------------------------------------ */
      JSON settingJSON(std::ifstream(input_directory + "/setting.json"));
      for (auto &[key, value] : settingJSON())
         std::cout << key << ": " << value << std::endl;
      //
      if (!settingJSON.find("end_time_step"))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "setting.json does not have end_time_step");
      if (!settingJSON.find("end_time"))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "setting.json does not have  end_time");
      if (!settingJSON.find("output_directory"))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "setting.json does not have output_directory");
      if (!settingJSON.find("max_dt"))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "setting.json does not have max_dt");
      if (!settingJSON.find("input_files"))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "setting.json does not have input_files");
      //
      std::map<Network *, outputInfo> NetOutputInfo;
      std::vector<Network *> FluidObject, RigidBodyObject, SoftBodyObject;
      std::string output_directory = settingJSON["output_directory"][0];
      double max_dt = stod(settingJSON["max_dt"])[0];
      const int end_time_step = stoi(settingJSON["end_time_step"])[0];
      const double end_time = stod(settingJSON["end_time"])[0];
      std::filesystem::create_directories(output_directory);
      std::filesystem::copy_file(input_directory + "/setting.json", output_directory + "/setting.json", std::filesystem::copy_options::overwrite_existing);
      /* ------------------------------------------------------ */
      double stop_remesh_time = 1E+10;
      double force_remesh_time = 0;
      if (settingJSON.find("stop_remesh_time"))
         stop_remesh_time = stod(settingJSON["stop_remesh_time"][0]);
      if (settingJSON.find("force_remesh_time"))
         force_remesh_time = stod(settingJSON["force_remesh_time"][0]);
      int grid_refinement = 0;
      if (settingJSON.find("grid_refinement"))
         grid_refinement = stoi(settingJSON["grid_refinement"][0]);
      /* ------------------------------------------------------ */

      for (auto FileName : settingJSON["input_files"]) {
         Print("---------------------------------------------------");
         Print("-----------------" + input_directory + FileName + "------------------");
         Print("---------------------------------------------------");
         JSON J(std::ifstream(input_directory + FileName));
         for (auto &[key, value] : J())
            std::cout << key << ": " << value << std::endl;
         if (!J.find("ignore") || !stob(J["ignore"])[0]) {
            auto net = new Network(J["objfile"][0], J["name"][0]);
            net->inputJSON = J;
            // for (auto &p : net->getPoints())
            // 	p->radius = stod(J["radius"])[0];
            // if (J.find("velocity"))
            // {
            // 	std::get<0>(net->velocity_name_start) = stob(J["velocity"][0] /*name*/)[0];
            // 	std::get<1>(net->velocity_name_start) = stod(J["velocity"][1] /*start*/)[0];
            // }
            if (J.find("mass")) {
               net->mass = stod(J["mass"])[0];
               std::get<0>(net->inertia) = stod(J["mass"])[0];
               std::get<1>(net->inertia) = stod(J["mass"])[0];
               std::get<2>(net->inertia) = stod(J["mass"])[0];
            }
            if (J.find("MOI")) {
               std::get<3>(net->inertia) = stod(J["MOI"])[0];
               std::get<4>(net->inertia) = stod(J["MOI"])[1];
               std::get<5>(net->inertia) = stod(J["MOI"])[2];
            }
            if (J.find("COM"))
               net->COM = net->initial_center_of_mass = Tddd{stod(J["COM"])[0], stod(J["COM"])[1], stod(J["COM"])[2]};

            modify(*net, J);
            NetOutputInfo[net].pvd_file_name = J["name"][0];
            NetOutputInfo[net].vtu_file_name = J["name"][0] + "_";
            NetOutputInfo[net].PVD = new PVDWriter(output_directory + "/" + J["name"][0] + ".pvd");
            std::filesystem::copy_file(input_directory + FileName, output_directory + "/" + FileName, std::filesystem::copy_options::overwrite_existing);
            if (J["type"][0] == "RigidBody") {
               RigidBodyObject.emplace_back(net);
               net->isRigidBody = true;
               net->isSoftBody = false;
            } else if (J["type"][0] == "SoftBody" || J["type"][0] == "FixedBody") {
               SoftBodyObject.emplace_back(net);
               net->isRigidBody = false;
               net->isSoftBody = true;
            } else if (J["type"][0] == "Fluid") {
               FluidObject.emplace_back(net);
               net->isRigidBody = net->isSoftBody = false;
            }
            //
            net->isFixed = (J.find("isFixed") && stob(J["isFixed"])[0]);
            //
            // 2022/12/14
            for (auto i = 0; i < 10; ++i)
               AreaWeightedSmoothingPreserveShape(net->getPoints(), 0.1);
            net->resetInitialX();
            net->setGeometricProperties();
            //
            mk_vtu(output_directory + "/" + J["name"][0] + "_init.vtu", {net->getFaces()});

            //

         } else {
            Print("skipped");
         }
         Print("---------------------------------------------------");
      }
      std::filesystem::copy_file("./main.cpp", output_directory + "/main.cpp", std::filesystem::copy_options::overwrite_existing);
      std::filesystem::copy_file("./BEM_derivatives.hpp", output_directory + "/BEM_derivatives.hpp", std::filesystem::copy_options::overwrite_existing);
      auto water = FluidObject[0];
      /* ------------------------------------------------------ */
      PVDWriter cornerPointsPVD(output_directory + "/cornerPointsPVD.pvd");
      PVDWriter cornerPVD(output_directory + "/corner.pvd");

      Print("setting done");
      //  b* ------------------------------------------------------ */
      //  b*                         メインループ                     */
      //  b* ------------------------------------------------------ */
      TimeWatch watch;
      for (time_step = 0; time_step < end_time_step; time_step++) {
         if (end_time < real_time)
            break;
         show_info(*water);
         //! 体積を保存するようにリメッシュする必要があるだろう．
         // auto radius = Mean(extLength(water->getLines()));
         setBoundaryConditions(*water, Join(RigidBodyObject, SoftBodyObject));
         double rad = M_PI / 180;
         // flipIf(*water, {10 * rad, rad}, {5 * rad /*結構小さく*/, rad}, false);
         flipIf(*water, {10 * rad, rad}, {10 * rad /*結構小さく*/, rad}, false);
         // b# ------------------------------------------------------ */
         // b#                       刻み時間の決定                     */
         // b# ------------------------------------------------------ */
         const auto Points = water->getPoints();
         const auto Faces = water->getFaces();
         //
         double dt = max_dt;
         Print("===========================================================================");
         Print("       dt :" + Red + std::to_string(dt) + colorOff);
         Print("time_step :" + Red + std::to_string(time_step) + colorOff);
         Print("real time :" + Red + std::to_string(real_time) + colorOff);
         Print("---------------------------------------------------------------------------");

         double spacing = Mean(extLength(water->getLines())) * 3;
         Buckets<networkFace *> FMM_BucketsFaces(water->bounds, spacing);
         Buckets<networkPoint *> FMM_BucketsPoints(water->bounds, spacing);
         for (const auto &f : water->getFaces())
            FMM_BucketsFaces.add(f->getXtuple(), f);
         for (const auto &p : water->getPoints())
            FMM_BucketsPoints.add(ToX(p), p);

         mk_vtu(output_directory + "/FMM_BucketsFaces.vtu", FMM_BucketsFaces.getT4Tddd());
         // b@ ------------------------------------------------------ */
         // b@        初期値問題を解く（時間微分方程式を数値積分する）           */
         // b@ ------------------------------------------------------ */
         int RK_order = 4;
         for (const auto &p : Points) {
            p->RK_phi.initialize(dt, real_time, std::get<0>(p->phiphin), RK_order);
            p->RK_X.initialize(dt, real_time, ToX(p), RK_order);
         }
         for (const auto &net : RigidBodyObject) {
            net->RK_COM.initialize(dt, real_time, net->COM, RK_order);
            net->RK_Q.initialize(dt, real_time, net->Q(), RK_order);
            net->RK_Velocity.initialize(dt, real_time, net->velocity, RK_order);
         }
         // b$ --------------------------------------------------- */
         for (const auto &net : SoftBodyObject) {
            // !いらないはずのもの
            for (const auto &p : net->getPoints())
               p->RK_X.initialize(dt, real_time, ToX(p), RK_order);
         }

         int RK_step = 0;

         do {
            //! 壁面の動きは，マイステップ更新することにした．この結果はphin()で参照される
            auto RK_time = (*Points.begin())->RK_X.gett();  //%各ルンゲクッタの時刻を使う
            std::cout << "RK_step = " << ++RK_step << "/" << RK_order << ", RK_time = " << RK_time << ", real_time = " << real_time << std::endl;
            // b# ------------------------------------------------------ */
            // b#      物体のノイマン境界の速度 u(t) at Neumann を設定         */
            // b# ------------------------------------------------------ */
            for (const auto &net : RigidBodyObject) {
               std::cout << "----------------" << std::endl;
               std::cout << net->getName() << "　の流速の計算方法" << std::endl;
               if (net->inputJSON.find("velocity")) {
                  std::string move_name = net->inputJSON["velocity"][0];
                  std::cout << "move_name = " << move_name << std::endl;
                  if (move_name == "fixed") {
                     net->velocity = {0., 0., 0., 0., 0., 0.};
                     net->acceleration = {0., 0., 0., 0., 0., 0.};
                  } else if (move_name != "floating") {
                     net->velocity = velocity(move_name, net->inputJSON["velocity"], RK_time);  // T6d //@ Φnを計算するために，物体表面の速度forced_velocityは，保存しておく必要がある
                                                                                                // net->acceleration = forced_motion::acceleration(RK_time); // T6d //@ 圧力を計算するために，物体表面の加速度は，保存しておく必要がある
                  } else if (move_name == "floating") {
                     std::cout << "floatingの場合は，加速度の時間積分によってシミュレートされる" << std::endl;
                  }
               } else {
                  std::cout << "指定がないので速度はゼロ" << std::endl;
                  net->velocity = {0., 0., 0., 0., 0., 0.};
                  net->acceleration = {0., 0., 0., 0., 0., 0.};
               }
               std::cout << "----------------" << std::endl;
            }
            // b$ --------------------------------------------------- */
            for (const auto &net : SoftBodyObject) {
               std::cout << "----------------" << std::endl;
               std::cout << net->getName() << "　の流速の計算方法．soft bodyの場合，各節点に速度を与える．" << std::endl;
               net->velocity = {0., 0., 0., 0., 0., 0.};
               net->acceleration = {0., 0., 0., 0., 0., 0.};
               if (net->inputJSON.find("velocity")) {
                  std::string move_name = net->inputJSON["velocity"][0];
                  std::cout << "move_name = " << move_name << std::endl;
                  if (move_name == "fixed") {
                     for (const auto &p : net->getPoints()) {
                        p->velocity = {0., 0., 0., 0., 0., 0.};
                        p->acceleration = {0., 0., 0., 0., 0., 0.};
                     }
                  } else {
                     for (const auto &p : net->getPoints())
                        p->velocity = velocity(move_name, net->inputJSON["velocity"], p, RK_time);  // T6d //@ Φnを計算するために，物体表面の速度forced_velocityは，保存しておく必要がある
                  }
               } else {
                  std::cout << "指定がないので速度はゼロ" << std::endl;
                  for (const auto &p : net->getPoints()) {
                     p->velocity = {0., 0., 0., 0., 0., 0.};
                     p->acceleration = {0., 0., 0., 0., 0., 0.};
                  }
               }
               std::cout << "----------------" << std::endl;
            }
            // やはり，ルンゲクッタは3回目に同じ場所での計算を行うわけなので，角点に関しては，毎回の時間発展＆修正はまずいだろう．
            // 予測するΩ(t+δt)は，ルンゲクッタに合わせたものでないとおかしいだろう．つまり，RK4なら　Ω(t0),Ω(t1),Ω(t2=t1),Ω(t3)でないといけないだろう．
            // b% -------------------------------------------------------- */
            // b%            境界条件（角点・ディリクレ・ノイマン）の決定               */
            // b% -------------------------------------------------------- */
            setBoundaryConditions(*water, Join(RigidBodyObject, SoftBodyObject));
            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
            // b! ------------------------------------------------------ */
            // b!           　境界値問題を解く-> {Φ,Φn}が決まる                */
            // b! ------------------------------------------------------ */
            std::cout << Green << "境界値問題を解く-> {Φ,Φn}が決まる" << colorOff << std::endl;
            BEM_BVP BVP;
            BVP.solve(*water, FMM_BucketsPoints, FMM_BucketsFaces);
            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
            // b* ------------------------------------------------------ */
            // b*                    微分∇ΦやDUDtを計算                    */
            // b* ------------------------------------------------------ */
            std::cout << Green << "微分∇ΦやDUDtを計算" << colorOff << std::endl;
            derivatives ders(*water);
            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
            // b* ------------------------------------------------------ */
            // b*           　境界値問題を解く-> {Φt,Φtn}が決まる              */
            // b* ------------------------------------------------------ */
            std::cout << Green << "境界値問題を解く-> {Φt,Φtn}が決まる" << colorOff << std::endl;
            BVP.solveForPhiPhin_t();
            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
            // b# ------------------------------------------------------ */
            // b#         物体のノイマン境界の加速度 accel(t) を計算            */
            // b# ------------------------------------------------------ */
            for (const auto &net : RigidBodyObject) {
               if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] == "floating") {
                  std::cout << net->inputJSON.find("velocity") << std::endl;
                  std::cout << net->inputJSON["velocity"][0] << std::endl;
                  auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
                  auto [mx, my, mz, Ix, Iy, Iz] = net->getInertiaGC();
                  auto force = tmp.surfaceIntegralOfPressure() + _GRAVITY3_ * net->mass;
                  force = force * Tddd{1., 1., 0.5};
                  auto torque = tmp.getFroudeKrylovTorque(net->COM);
                  auto [a0, a1, a2] = force / Tddd{mx, my, mz};
                  auto [a3, a4, a5] = torque / Tddd{Ix, Iy, Iz};
                  net->acceleration = T6d{a0, a1, a2, a3, a4, a5};
                  std::cout << red << net->getName() << "\n"
                            << ", net->getInertiaGC() = " << net->getInertiaGC() << "\n"
                            << ", net->mass = " << net->mass << "\n"
                            << ", tmp.surfaceIntegralOfPressure() = " << tmp.surfaceIntegralOfPressure() << "\n"
                            << ", tmp.getFroudeKrylovTorque(net->COM) = " << tmp.getFroudeKrylovTorque(net->COM) << "\n"
                            << ", tmp.area = " << tmp.area << "\n"
                            << ", net->velocity = " << net->velocity << "\n"
                            << ", net->acceleration = " << net->acceleration << "\n"
                            << ", net->COM = " << net->COM << "\n"
                            << colorOff << std::endl;
               }
            }
            // b* ------------------------------------------------------ */
            // b*                 ディリクレ境界ではΦを時間積分                 */
            // b* ------------------------------------------------------ */
            std::cout << Green << "ディリクレ境界ではΦを時間積分，ノイマン境界ではΦnを陽に与える" << colorOff << std::endl;
            for (const auto &net : RigidBodyObject) {
               std::cout << "name:" << net->getName() << std::endl;
               if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] != "fixed") {
                  net->RK_COM.push(net->velocityTranslational());
                  net->COM = net->RK_COM.getX();
                  Quaternion q;
                  q = q.d_dt(net->velocityRotational());  // w->クォータニオン
                  net->RK_Q.push(q());                    // クォータニオン->T4dとしてプッシュ
                  net->Q = net->RK_Q.getX();
               }
               if (!net->inputJSON.find("acceleration") || (net->inputJSON.find("acceleration") && net->inputJSON["acceleration"][0] != "fixed")) {
                  net->RK_Velocity.push(net->acceleration);
                  net->velocity = net->RK_Velocity.getX();
               }
               net->RigidBodyMovePoints();
            }
            // b$ --------------------------------------------------- */
            for (const auto &net : SoftBodyObject) {
               std::cout << "name:" << net->getName() << std::endl;
               for (const auto &p : net->getPoints()) {
                  p->RK_X.push(p->velocityTranslational());  //@ 位置xの時間発展

                  // p->setXSingle(p->RK_X.getX());
               }
               net->setGeometricProperties();
            }

            for (const auto &p : Points) {
               /* ------------------------ */
               if (!p->Neumann) {
                  // b! ここでノイマン面が変な場所に移動させられてはたまらない
                  //@ Φの時間発展，Φnの時間発展はない
                  p->RK_phi.push(p->DphiDt(p->U_update_BEM, 0.));
                  std::get<0>(p->phiphin) = p->phi_Dirichlet = p->RK_phi.getX();  // 角点の法線方向はわからないので，ノイマンの境界条件phinを与えることができない．
               }
               p->RK_X.push(p->U_update_BEM);  //@ 位置xの時間発展
               p->setXSingle(p->RK_X.getX());
            }

            /* ------------------------------------------------------ */
            std::cout << Green << "name:" << water->getName() << ": setBounds" << colorOff << std::endl;
            water->setGeometricProperties();

            std::ofstream ofs(output_directory + "/water" + std::to_string(RK_step) + ".obj");
            creteOBJ(ofs, *water);
            ofs.close();

            std::cout << Blue << "Elapsed time: " << Red << watch() << colorOff << " s\n";
         } while (!((*Points.begin())->RK_X.finished));

         // b$ ------------------------------------------------------ */
         for (const auto &net : SoftBodyObject)
            AreaWeightedSmoothingPreserveShape(net->getPoints(), 1);

         /* ------------------------------------------------------ */
         {
            double mean_phi = 0.;
            for (const auto &p : water->getPoints())
               mean_phi += std::get<0>(p->phiphin);
            mean_phi /= water->getPoints().size();
            for (const auto &p : water->getPoints()) {
               p->phi_Dirichlet -= mean_phi;
               // p->phi_Neumann -= mean_phi;
               std::get<0>(p->phiphin) -= mean_phi;
            }
         }

         /* ------------------------------------------------------ */
         std::cout << Green << "real_timeを取得" << colorOff << std::endl;
         real_time = (*Points.begin())->RK_X.gett();

         for (const auto &p : Points) {
            p->U_BEM_last = p->U_BEM;
            p->U_tangential_BEM_last = p->U_tangential_BEM;
         }

         /* ------------------------------------------------------ */

         std::ofstream ofs(output_directory + "/water_current.obj");
         creteOBJ(ofs, *water);
         ofs.close();

         // b# -------------------------------------------------------------------------- */
         // b#                              output JSON files                             */
         // b# -------------------------------------------------------------------------- */
         jsonout.push("time", real_time);
         jsonout.push(water->getName() + "_volume", water->getVolume());
         for (const auto &net : Join(RigidBodyObject, SoftBodyObject)) {
            if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] == "floating") {
               auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
               jsonout.push(net->getName() + "_pitch", net->quaternion.pitch());
               jsonout.push(net->getName() + "_yaw", net->quaternion.yaw());
               jsonout.push(net->getName() + "_roll", net->quaternion.roll());
               jsonout.push(net->getName() + "_force", tmp.surfaceIntegralOfPressure());
               jsonout.push(net->getName() + "_torque", tmp.getFroudeKrylovTorque(net->COM));
               jsonout.push(net->getName() + "_accel", net->acceleration);
               jsonout.push(net->getName() + "_velocity", net->velocity);
               jsonout.push(net->getName() + "_COM", net->COM);
               jsonout.push(net->getName() + "_area", tmp.area);
            }
         }
         std::ofstream os(output_directory + "/result.json");
         jsonout.output(os);
         os.close();

         // b# -------------------------------------------------------------------------- */
         // b#                            output Paraview files                           */
         // b# -------------------------------------------------------------------------- */

         auto data = dataForOutput(*water, dt);
         // 流体
         for (const auto &net : FluidObject) {
            auto filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, net->getFaces(), data);
            NetOutputInfo[net].PVD->push(filename, real_time);
            NetOutputInfo[net].PVD->output();
         }

         for (const auto &net : Join(RigidBodyObject, SoftBodyObject)) {
            VV_VarForOutput data;
            if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] == "floating") {
               auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
               //    uomap_P_Tddd P_COM, P_COM_p, P_accel, P_velocity, P_rotational_velocity, P_rotational_accel, P_FroudeKrylovTorque;
               //    uomap_P_d P_Pitch, P_Yaw, P_Roll, P_pressure;

               //    for (const auto &p : water->getPoints()) {
               //       P_accel[p] = net->accelRigidBody(ToX(p));
               //       P_velocity[p] = net->velocityRigidBody(ToX(p));
               //       P_rotational_velocity[p] = net->velocityRotational();
               //       // P_FroudeKrylovTorque[p] = tmp.getFroudeKrylovTorque(net->COM);
               //       P_pressure[p] = p->pressure;
               //    }

               //    for (const auto &p : net->getPoints()) {
               //       P_COM[p] = net->COM;
               //       P_COM_p[p] = net->COM - ToX(p);
               //       P_Pitch[p] = net->quaternion.pitch();
               //       P_Yaw[p] = net->quaternion.yaw();
               //       P_Roll[p] = net->quaternion.roll();
               //       P_accel[p] = net->accelRigidBody(ToX(p));
               //       P_velocity[p] = net->velocityRigidBody(ToX(p));
               //       P_rotational_velocity[p] = net->velocityRotational();
               //    }

               //    data = {
               //        {"vector to COM", P_COM_p},
               //        {"COM", P_COM},
               //        {"pitch", P_Pitch},
               //        {"yaw", P_Yaw},
               //        {"roll", P_Roll},
               //        {"velocity", P_accel},
               //        {"acceleration", P_velocity},
               //        {"rotational valocity", P_rotational_velocity},
               //        {"rotational acceleration", P_rotational_accel},
               //        {"pressure", P_pressure},
               //        {"FroudeKrylovTorque", P_FroudeKrylovTorque}};
               mk_vtu(output_directory + "/actingFacesOn" + NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu", tmp.actingFaces, data);
            }
            auto filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, net->getFaces(), data);
            NetOutputInfo[net].PVD->push(filename, real_time);
            NetOutputInfo[net].PVD->output();
         }
         // 流体
         // {
         //    std::vector<Tddd> points;
         //    for (const auto &p : water->getPoints())
         //       if (p->CORNER)
         //          points.emplace_back(ToX(p));

         //    auto filename = "cornerPoints" + std::to_string(time_step) + ".vtu";
         //    mk_vtu(output_directory + "/" + filename, points);
         //    cornerPointsPVD.push(filename, real_time);
         //    cornerPointsPVD.output();
         // }
         // {
         //    std::unordered_set<networkFace *> faces;
         //    for (const auto &f : water->getFaces()) {
         //       // for (const auto p : f->getPoints())
         //       // 	if (p->CORNER)
         //       // 	{
         //       // 		faces.emplace(f);
         //       // 		break;
         //       // 	}
         //       for_each(f->getPoints(), [&](const auto &p) {
         //          if (p->CORNER)
         //             faces.emplace(f);
         //       });
         //    }
         //    auto filename = "corner" + std::to_string(time_step) + ".vtu";
         //    mk_vtu(output_directory + "/" + filename, faces, data);
         //    cornerPVD.push(filename, real_time);
         //    cornerPVD.output();
         // }
         // b# -------------------------------------------------------------------------- */
         //
         // mk_vtu(output_directory + "/" + obj.getName() + std::to_string(time_step) + ".vtu", obj.getFaces(), datacpg);
         // mk_vtu(output_directory + "/" + name + std::to_string(time_step) + ".vtu", cpg.well_Faces, datacpg);
         // cpg_pvd.push(name, name + std::to_string(time_step) + ".vtu", time_step * t_rep * dt);
         // cpg_pvd.output_();
         // b* ------------------------------------------------------ */
      }
   } catch (error_message &e) {
      e.print();
   };
   return 0;
};
