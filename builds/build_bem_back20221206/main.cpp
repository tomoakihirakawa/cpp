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

// #define debug_bem
#include "bem.hpp"
#include "svd.hpp"

using V_i = std::vector<int>;
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
using V_Netp = std::vector<Network *>;
using V_netFp = std::vector<networkFace *>;
using VV_netFp = std::vector<V_netFp>;

// #define rotation_test
// #define up_down_test
// #define testCase1
// #define experiment_Chaplin1999_Retzler2000
// #define experiment_Retzler2000
// #define experiment_Li2002
#define experiment_Goring1979
// #define Hu2002
// #define WenWenLi2002_MachReflection
// #define experiment_sawai
// #define XueAndLin2011_large_amplitude

/* ------------------------------------------------------ */
T6d velocity(const std::string &name, const std::vector<std::string> strings, networkPoint *p, const double t) {
   /* ------------------------------------------------------ */
   std::cout << "net->inputJSON[\" velocity \"] = " << strings << std::endl;
   auto g = _GRAVITY_;
   return {0., 0., 0., 0., 0., 0.};
};
T6d velocity(const std::string &name, const std::vector<std::string> strings, const double t) {
   /* ------------------------------------------------------ */
   std::cout << "net->inputJSON[\" velocity \"] = " << strings << std::endl;
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

   } else
      return {0., 0., 0., 0., 0., 0.};
};

// T6d rotation(const std::string &name, const std::vector<std::string> strings, const double t) {
//    double start = stod(strings[1] /*start*/);
//    double A, w, h, l, d, k;
//    if (strings.size() > 3) {
//       A = stod(strings[2] /*A*/);
//       w = 2 * M_PI / stod(strings[3] /*T*/);
//       h = stod(strings[4] /*h*/);
//       l = stod(strings[5] /*l*/);
//       DispersionRelation DS(w, h);
//       k = DS.k;
//       std::cout << "A = " << A << ", w = " << w << ", k = " << k << std::endl;
//       std::cout << "{T, L} = {" << DS.T << ", " << DS.L << "}" << std::endl;
//       double d = (l >= 0 ? d : -l);
//       Tddd axis = {stod(strings[6]), stod(strings[7]), stod(strings[8])};
//       double rot = Normalize(axis) * ArcTan((A * g * k * (1 + 2 * h * k * Csch(2 * h * k)) * Sin(t * w)), (2. * (-g + (h + l) * Power(w, 2) + g * Cosh(d * k) * Sech(h * k))));
//       return {0., 0., 0., std::get<0>(rot), std::get<1>(rot), std::get<2>(rot)};
//    } else
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be > 3. amplitude and frequency");
// }

/* ------------------------------------------------------ */

// std::unordered_set<networkFace *> facingFace(const networkFace *const f) {
//    //! getContactFacesは，点の周辺の点の影となった面を含まない．
//    //! しかし，面の接触面を判断するには，そのような除外だけでは十分ではなく，
//    //! さらに，この面と向き合っているかどうかの判定し除外する必要がある．それが以下．
//    std::unordered_set<networkFace *> ret;
//    // for (const auto &q : f->getPoints())
//    // 	for (const auto &F : q->getContactFaces())
//    // 	{
//    // 		// auto angle = MyVectorAngle(f->normal, -F->normal);
//    // 		// if ((angle / M_PI * 180. < 60.) || ((M_PI - angle) / M_PI * 180. < 60.))
//    // 		ret.emplace(F);
//    // 	}

//    // for_each(f->getPoints(),
//    //          [&](const auto &q) {
//    //             for (const auto &F : q->getContactFaces()) {
//    //                ret.emplace(F);
//    //             }
//    //          });
//    for_each(f->getPoints(), [&](const auto &q) { ret.insert(q->getContactFaces().begin(), q->getContactFaces().end()); });
//    return ret;
// };

/* ------------------------------------------------------ */

netFp closestFacingFace(const networkFace *const f_IN) {
   //! もしない場合はnullptrを返すので注意
   // networkFace *closest_face = nullptr;
   // double min_distance = 1E+100;
   // for (const auto &F : facingFace(f_IN)) {
   //    // auto dist = Norm(vectorToTriangle(F, f_IN->getXtuple()));
   //    auto dist = Distance(f_IN->center, F);
   //    if (dist < min_distance) {
   //       closest_face = F;
   //       min_distance = dist;
   //    }
   // }
   // auto [X, closest_face] = Nearest_(f_IN->center, facingFace(f_IN));
   // return closest_face;
   //
   std::unordered_set<networkFace *> faces;
   for_each(f_IN->getPoints(), [&](const auto &q) { faces.insert(q->getContactFaces().begin(), q->getContactFaces().end()); });
   return std::get<1>(Nearest_(f_IN->X, faces));
};
/* ------------------------------------------------------ */

std::unordered_map<networkFace *, networkFace *> getNearestContactFacesOfSurroundedNeumannFace(const networkPoint *const p) {
   std::unordered_map<networkFace *, networkFace *> structure_face;
   networkFace *tmp;
   //@ 点の隣接面のうちNeumann面を抜き出し，それぞれに対してclosestFacingFaceを取得する．
   for (const auto &F : p->getFaces())
      if (F->Neumann && (tmp = closestFacingFace(F) /*nulltrなら入らない*/))
         structure_face[F] = tmp;
   return structure_face;
};

// Tddd getNearestContactFacesX(const Tddd &X,
//                              const double radius,
//                              const std::vector<T3Tddd> &vertices) {
//    Tddd ret = {1E+100, 1E+100, 1E+100};
//    for (const auto &vertex : vertices) {
//       auto intxn = IntersectionSphereTriangle(X, radius, vertex);
//       if (intxn.isIntersecting) {
//          auto r = intxn.X - X;
//          auto nr = Norm(r);
//          if (Norm(ret - X) > nr)
//             ret = intxn.X;
//       }
//    }
//    return ret;
// };
// std::tuple<networkFace *, Tddd> getNearestContactFacesX(const networkPoint *const p) {
//    std::tuple<networkFace *, Tddd> ret = {nullptr, {1E+100, 1E+100, 1E+100}};
//    auto pX = ToX(p);
//    if (!p->getContactFaces().empty())
//       for (const auto &[f, X] : p->getContactFacesX()) {
//          if (isFinite(X - pX))
//             if (Norm(std::get<1>(ret) - pX) > Norm(X - pX))
//                ret = {f, X};
//       }
//    return ret;
// };

Tddd uNeumann(const networkPoint *const p) {
   auto tmpContactFaces = p->getContactFacesXCloser();
   if (tmpContactFaces.empty())
      return {0., 0., 0.};
   std::vector<Tddd> V;
   std::vector<double> W;
   Tddd v;
   double w;
   for (const auto &[f, X] : tmpContactFaces) {
      v = f->getNetwork()->velocityRigidBody(X);
      w = kernel_Bspline3(Norm(X - ToX(p)), p->radius);
      V.emplace_back(v);
      W.emplace_back(w);
   }
   auto [f, X] = tmpContactFaces[0];
   auto ret = optimumVector_(V, {0., 0., 0.}, W);
   if (isFinite(ret))
      return ret;
   else
      return {0., 0., 0.};
};

// Tddd uNeumann(const networkPoint *const p, const Tddd &n) {
//    auto tmpContactFaces = p->getContactFacesXCloser();
//    if (tmpContactFaces.empty())
//       return {0., 0., 0.};
//    std::vector<Tddd> V, N;
//    std::vector<double> W;
//    Tddd v;
//    double w;
//    for (const auto &[f, X] : tmpContactFaces) {
//       v = f->getNetwork()->velocityRigidBody(X);
//       w = kernel_Bspline3(Norm(X - ToX(p)), p->radius);
//       V.emplace_back(v);
//       W.emplace_back(w);
//       N.emplace_back(n);
//    }
//    auto [f, X] = tmpContactFaces[0];
//    return optimumVector_(V, {0., 0., 0.}, N, W);
// };

Tddd accelNeumann(const networkPoint *const p) {
   auto tmpContactFaces = p->getContactFacesXCloser();
   if (tmpContactFaces.empty())
      return {0., 0., 0.};
   std::vector<Tddd> V;
   std::vector<double> W;
   Tddd v;
   double w;
   for (const auto &[f, X] : tmpContactFaces) {
      v = f->getNetwork()->accelRigidBody(X);
      w = kernel_Bspline3(Norm(X - ToX(p)), p->radius);
      V.emplace_back(v);
      W.emplace_back(w);
   }
   auto [f, X] = tmpContactFaces[0];
   auto ret = optimumVector_(V, {0., 0., 0.}, W);
   if (isFinite(ret))
      return ret;
   else
      return {0., 0., 0.};
};

/* ------------------------------------------------------ */

// double phi_stokes2nd(double x, double z, const double t)
// {
// 	const double g = 9.8;
// 	auto tmp = a * w / k / sinh(k * h);
// 	tmp *= (cosh(k * (z + h)) * sin(q) + k * a * 3 * cosh(2 * k * (z + h)) / (8 * pow(sinh(k * h), 3)) * sin(2 * q));
// 	tmp -= pow(k * a, 2) / (2 * sinh(2 * k * h) * g * t / k);
// 	return tmp;
// };
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
double accel_normal_from_Neumann_surface(const networkPoint *const p) {
   /*
    * 面積平均に変更した
    */
   //@ 点の隣接面のうちNeumann面を抜き出し，それぞれに対してclosestFacingFaceを取得する．
   //@ 次に，面それぞれの速度を計算し，点の法線方向速度を計算する．
   double ret = 0., A = 0, Atot = 0;
   int i = 0;
   for (auto &[sF, cF] : getNearestContactFacesOfSurroundedNeumannFace(p) /*この点に隣接する流体面それぞれが最も近くで接している構造物の面*/) {
      A = sF->area;
      Atot += A;
      ret += A * Dot(cF->getNetwork()->accelRigidBody(ToX(p)), sF->normal);
      i++;
   }
   return i == 0 ? ret : ret / Atot;
}
/* ------------------------------------------------------ */
double accel_normal_from_Neumann_surface(const networkPoint *const p, const networkFace *const f) {
   return Dot(closestFacingFace(f)->getNetwork()->accelRigidBody(ToX(p)), f->normal);
}
/* ------------------------------------------------------ */
double velocity_normal_from_Neumann_surface(const networkPoint *const p) {
   /*
    * 面積平均に変更した
    */
   //@ 点の隣接面のうちNeumann面を抜き出し，それぞれに対してclosestFacingFaceを取得する．
   //@ 次に，面それぞれの速度を計算し，点の法線方向速度を計算する．
   double ret = 0., A = 0, Atot = 0;
   int i = 0;
   for (auto &[sF, cF] : getNearestContactFacesOfSurroundedNeumannFace(p) /*この点に隣接する流体面それぞれが最も近くで接している構造物の面*/) {
      A = sF->area;
      Atot += A;
      ret += A * Dot(cF->getNetwork()->velocityRigidBody(ToX(p)), sF->normal);
      i++;
   }
   return i == 0 ? ret : ret / Atot;
}
double velocity_normal_from_Neumann_surface(const networkPoint *const p, const networkFace *const f) {
   return Dot(closestFacingFace(f)->getNetwork()->velocityRigidBody(ToX(p)), f->normal);
}
/* ------------------------------------------------------ */
T6d velocity_from_Neumann_surface(const networkPoint *const p) {
   T6d ret = {0., 0., 0., 0., 0., 0.};
   double A = 0, Atot = 0;
   int i = 0;
   for (auto &[sF, cF] : getNearestContactFacesOfSurroundedNeumannFace(p) /*この点に隣接する流体面それぞれが最も近くで接している構造物の面*/) {
      auto n_cF = cF->normal;
      auto n_sF = sF->normal;
      A = sF->area;
      Atot += A;
      Tddd V = cF->getNetwork()->velocityRigidBody(ToX(p));
      ret += A * Dot(V, p->getNormalTuple());
      i++;
   }
   return i == 0 ? ret : ret / Atot;
}
// Tddd local_velocity_from_Neumann_surface(const networkPoint *const p) {
//    Tddd ret = {0., 0., 0.};
//    double A = 0, Atot = 0;
//    int i = 0;
//    for (auto &[sF, cF] : getNearestContactFacesOfSurroundedNeumannFace(p) /*この点に隣接する流体面それぞれが最も近くで接している構造物の面*/) {
//       auto n_cF = cF->normal;
//       auto n_sF = sF->normal;
//       A = sF->area;
//       Atot += A;
//       Tddd V = cF->getNetwork()->velocityRigidBody(ToX(p));
//       ret += A * V;
//       i++;
//    }
//    return i == 0 ? ret : ret / Atot;
// }
/* ------------------------------------------------------ */
double phin_contact(const netFp f_IN) {
   auto closest_face = closestFacingFace(f_IN);
   if (closest_face) {
      // auto n_point = p->getNormalTuple();
      // auto U_normal_of_face = Dot(f->getNetwork()->velocity, ToT6d(n_face)) * n_face;
      return Dot(closest_face->getNetwork()->velocity,
                 ToT6d(f_IN->getXtuple()));
   } else {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      return 0;
   }
}
/* ------------------------------------------------------ */
V_netLp getLinesAround(netPp p) {
   V_netLp ret = {};
   for (const auto &f : p->getFaces())
      for_each(f->getLines(), [&](const auto &l) { ret.emplace_back(l); });
   return DeleteDuplicates(ret);
};
V_netLp getLinesAround(netLp l) {
   V_netLp ret = {l};
   auto [p, q] = l->getPoints();
   for (const auto &f : p->getFaces())
      for_each(f->getLines(), [&](const auto &l) { ret.emplace_back(l); });
   for (const auto &f : q->getFaces())
      for_each(f->getLines(), [&](const auto &l) { ret.emplace_back(l); });
   return DeleteDuplicates(ret);
};
/* ------------------------------------------------------ */
// using V_RKRK = std::vector<std::shared_ptr<derivativeImprover>>;
// using V_RKRK = std::vector<derivativeImprover *>;
using map_P_d = std::map<netP *, double>;
using map_P_Vd = std::map<netP *, V_d>;
// using map_P_TiiiTdd = std::map<netP *, std::map<Tiii, Tdd>>;
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

VV_d extVelocities(const V_netPp &ps) {
   VV_d ret(0);
   for (const auto &p : ps)
      ret.emplace_back(ToVector(p->U_BEM));
   return ret;
};

//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
Tddd gradTangential_LinearElement(const Tddd &V, const T3Tddd &X012) {
   auto [X0, X1, X2] = X012;
   return Cross(TriangleNormal(X012),
                std::get<0>(V) * (X2 - X1) +
                    std::get<1>(V) * (X0 - X2) +
                    std::get<2>(V) * (X1 - X0)) /
          (2 * TriangleArea(X012));
};

// Tddd gradTangential_LinearElement(const Tddd &V, const networkFace *const f)
// {
// 	return gradTangential_LinearElement(V, f->getXVertices());
// 	// auto [X0, X1, X2] = f->getXVertices();
// 	// return Cross(f->getNormalTuple(),
// 	// 			 std::get<0>(V) * (X2 - X1) +
// 	// 				 std::get<1>(V) * (X0 - X2) +
// 	// 				 std::get<2>(V) * (X1 - X0)) /
// 	// 	   (2 * f->area);
// };

// using T3T2ddd = std::tuple<T2Tddd, T2Tddd, T2Tddd>;

// Tddd gradTangential_LinearElement(const T3T2ddd &V)
// {
// 	auto [XV0, XV1, XV2] = V;
// 	return gradTangential_LinearElement({std::get<0>(XV0), std::get<0>(XV1), std::get<0>(XV2)},
// 										{std::get<1>(XV0), std::get<1>(XV1), std::get<1>(XV2)});
// 	// auto [X0, X1, X2] = f->getXVertices();
// 	// return Cross(f->getNormalTuple(),
// 	// 			 std::get<0>(V) * (X2 - X1) +
// 	// 				 std::get<1>(V) * (X0 - X2) +
// 	// 				 std::get<2>(V) * (X1 - X0)) /
// 	// 	   (2 * f->area);
// };

Tddd gradPhiTangential(const networkFace *const f) {
   auto [p0, p1, p2] = f->getPoints();
   return gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
                                       {ToX(p0), ToX(p1), ToX(p2)});
}

Tddd gradPhiTangential(const networkFace *const f, const networkPoint *p) {
   auto [p0, p1, p2] = f->getPoints(p);
   return gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
                                       {ToX(p0), ToX(p1), ToX(p2)});
};

// Tddd gradPhi(const networkFace *const f)
// {
// 	auto [p0, p1, p2] = f->getPoints();
// 	Tddd ret = gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
// 											{ToX(p0), ToX(p1), ToX(p2)});
// 	if (f->Neumann)
// 		ret += f->getNormalTuple() * (p0->phin_Neumann + p1->phin_Neumann + p2->phin_Neumann) / 3.;
// 	else if (f->Dirichlet)
// 		ret += f->getNormalTuple() * (p0->phin_Dirichlet + p1->phin_Dirichlet + p2->phin_Dirichlet) / 3.;
// 	else
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	return ret;
// }

// Tddd gradPhi(const networkFace *const f,
// 			 const networkPoint *const p,
// 			 const double t0 = 1. /*phinにのみ関するもの*/)
// {
// 	Tddd N = {t0, (1 - t0) / 2., (1 - t0) / 2.};
// 	auto [p0, p1, p2] = f->getPoints(p);
// 	Tddd ret = gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
// 											{ToX(p0), ToX(p1), ToX(p2)});
// 	if (f->Neumann)
// 	{
// 		ret += f->getNormalTuple() * Dot(Tddd{p0->phin_Neumann, p1->phin_Neumann, p2->phin_Neumann}, N);
// 	}
// 	else if (f->Dirichlet)
// 		ret += f->getNormalTuple() * Dot(Tddd{p0->phin_Dirichlet, p1->phin_Dirichlet, p2->phin_Dirichlet}, N);
// 	else
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	return ret;
// }

Tddd gradPhiTangential(const networkFace *const f,
                       const networkPoint *const p,
                       const double t0 = 1. /*phinにのみ関するもの*/) {
   Tddd N = {t0, (1 - t0) / 2., (1 - t0) / 2.};
   auto [p0, p1, p2] = f->getPoints(p);
   Tddd ret = gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
                                           {ToX(p0), ToX(p1), ToX(p2)});
   return ret;
}

double getPhi(const networkLine *const l) {
   auto fs = l->getFaces();
   interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_0(fs[0], l);
   interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_1(fs[1], l);
   return Dot(intp_l0_0.N(.5, .5), ToPhi(intp_l0_0.Points)) * 0.5 +
          Dot(intp_l0_1.N(.5, .5), ToPhi(intp_l0_1.Points)) * 0.5;
}

T3Tddd grad_U_tangential_LinearElement(const networkFace *const f) {
   auto [p0, p1, p2] = f->getPoints();
   /*
           {grad(Ux),
           grad(Uy),
           grad(Uz)}
   */
   return {gradTangential_LinearElement({std::get<0>(p0->U_BEM), std::get<0>(p1->U_BEM), std::get<0>(p2->U_BEM)}, {ToX(p0), ToX(p1), ToX(p2)}),
           gradTangential_LinearElement({std::get<1>(p0->U_BEM), std::get<1>(p1->U_BEM), std::get<1>(p2->U_BEM)}, {ToX(p0), ToX(p1), ToX(p2)}),
           gradTangential_LinearElement({std::get<2>(p0->U_BEM), std::get<2>(p1->U_BEM), std::get<2>(p2->U_BEM)}, {ToX(p0), ToX(p1), ToX(p2)})};
};

// T3Tddd grad_U_LinearElement2(const networkPoint *const p)
// {
// 	T3Tddd grad_U_tangential = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
// 	double Atot = 0;
// 	for (const auto &f : p->getFaces())
// 	{
// 		grad_U_tangential += f->area * grad_U_tangential_LinearElement(f);
// 		Atot += f->area;
// 	}
// 	/*
// 		grad_U_tangential = {grad(Ux),grad(Uy),grad(Uz)}
// 	*/
// 	grad_U_tangential /= Atot;
// 	Tddd x = {1, 0, 0}, y = {0, 1, 0};
// 	auto V00 = Dot(std::get<0>(tmp), x);
// 	auto V10 = Dot(std::get<1>(tmp), x);
// 	auto V20 = Dot(std::get<2>(tmp), x);
// 	auto V01 = Dot(std::get<0>(tmp), y);
// 	auto V11 = Dot(std::get<1>(tmp), y);
// 	auto V21 = Dot(std::get<2>(tmp), y);
// 	return T3Tddd{Tddd{V00, V01, V20},
// 				  Tddd{V10, V11, V21},
// 				  Tddd{V20, V21, -V00 - V11}};
// };

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
/* ------------------------------------------------------ */

Tddd gradPhi(const networkPoint *const p) {
   const double t0 = 1;
   Tddd grad, n = p->getNormal_BEM();
   V_Tddd U, Us, V, N;
   V_d weights;
   double w;

   try {
      for (const auto &f : p->getFaces()) {
         auto [p0, p1, p2] = f->getPoints(p);
         Tddd ret = gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
                                                 {ToX(p0), ToX(p1), ToX(p2)});
         if (f->Neumann) {
            if (p->phinOnFace.find(f) != p->phinOnFace.end())
               ret += f->normal * p->phinOnFace.at(f);
            else if (p->phinOnFace.find(nullptr) != p->phinOnFace.end())
               ret += f->normal * p->phinOnFace.at(nullptr);
            else
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, std::to_string(p->phinOnFace.size()));
         } else if (f->Dirichlet)
            ret += f->normal * p->phin_Dirichlet;
         else
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

         w = f->area;
         U.emplace_back(ret * w);
         weights.emplace_back(w);
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };

   grad = Total(U) / Total(weights);
   return grad;
};

// b@ ------------------------------------------------------ */

Tddd nextX_U(const networkPoint *const p, const double dt) {
   return p->getXBuffer();
};
std::vector<Tddd> nextX_U(const std::vector<networkPoint *> &ps, const double dt) {
   std::vector<Tddd> ret;
   for (const auto &p : ps)
      ret.emplace_back(nextX_U(p, dt));
   return ret;
};
Tddd nextX_U_Ua(const networkPoint *const p, const double dt) {
   return p->getXBuffer() + p->U_BUFFER * dt;
};
std::vector<Tddd> nextX_U_Ua(const std::vector<networkPoint *> &ps, const double dt) {
   std::vector<Tddd> ret;
   for (const auto &p : ps)
      ret.emplace_back(nextX_U_Ua(p, dt));
   return ret;
};
std::tuple<Tddd, Tddd> nextX_U_Ua(const std::tuple<networkPoint *, networkPoint *> ps, const double dt) {
   return {std::get<0>(ps)->getXBuffer() + std::get<0>(ps)->U_BUFFER * dt,
           std::get<1>(ps)->getXBuffer() + std::get<1>(ps)->U_BUFFER * dt};
};
Tddd next_normal_U(const networkFace *const f, const double dt) {
   auto [p0, p1, p2] = f->getPoints();
   return TriangleNormal(f->getXVertices() + dt * T3Tddd{p0->getXBuffer(), p1->getXBuffer(), p2->getXBuffer()});
};

Tddd next_area_weighted_normal(const networkPoint *p, const double dt) {
   T3Tddd VVV;
   Tddd ret = {0., 0., 0.};
   for (const auto &f : p->getFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      VVV = f->getXVertices() + dt * T3Tddd{p0->getXBuffer(), p1->getXBuffer(), p2->getXBuffer()};
      ret += TriangleNormal(VVV) * TriangleArea(VVV);
   }
   return ret;
};

Tddd next_normal_U_Ua(const networkFace *const f, const double dt) {
   auto [p0, p1, p2] = f->getPoints();
   return TriangleNormal(f->getXVertices() + dt * T3Tddd{nextX_U_Ua(p0, dt), nextX_U_Ua(p1, dt), nextX_U_Ua(p2, dt)});
};
/* ------------------------------------------------------ */
bool isFacing(const Tddd &F, const Tddd &f, double rad = 1E-10) {
   return (Dot(F, -f) >= cos(rad));
};
Tddd closestXFacing(const Tddd &p_next_X, const double radius, const std::vector<T3Tddd> &vertices, const Tddd &n) {
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
Tddd closestX(const Tddd &p_next_X, const double radius, const std::vector<T3Tddd> &vertices) {
   Tddd r = {1E+100, 1E+100, 1E+100};
   for (const auto &vertex : vertices) {
      auto intxn = IntersectionSphereTriangle(p_next_X, radius, vertex);
      if (intxn.isIntersecting)
         if (Norm(r) >= Norm(intxn.X - p_next_X))
            r = intxn.X - p_next_X;
   }
   return r;
};

Tddd nextFaceNormal_U_Ua(const networkFace *const f, const double dt) {
   auto [p0, p1, p2] = f->getPoints();
   return TriangleNormal(nextX_U_Ua(p0, dt),
                         nextX_U_Ua(p1, dt),
                         nextX_U_Ua(p2, dt));
};
Tddd nextX_U_Ua(const networkFace *const f, const double dt) {
   // auto [p0, p1, p2] = f->getPoints();
   return (nextX_U_Ua(std::get<0>(f->getPoints()), dt) +
           nextX_U_Ua(std::get<1>(f->getPoints()), dt) +
           nextX_U_Ua(std::get<2>(f->getPoints()), dt)) /
          3.;
};
/* ------------------------------------------------------ */
Tddd nextX_U_Ua_Uc(const networkPoint *const p, const double dt) {
   return p->getXBuffer() + (p->U_BUFFER + p->U_cling_to_Neumann) * dt;
};

Tddd nextFaceNormal_U_Ua_Uc(const networkFace *const f, const double dt) {
   auto [p0, p1, p2] = f->getPoints();
   return TriangleNormal(nextX_U_Ua_Uc(p0, dt),
                         nextX_U_Ua_Uc(p1, dt),
                         nextX_U_Ua_Uc(p2, dt));
};
Tddd nextX_U_Ua_Uc(const networkFace *const f, const double dt) {
   auto [p0, p1, p2] = f->getPoints();
   return (nextX_U_Ua_Uc(p0, dt) + nextX_U_Ua_Uc(p1, dt) + nextX_U_Ua_Uc(p2, dt)) / 3.;
};
/* ------------------------------------------------------ */
std::vector<Tddd> vectorsClingToNeumann_old(const networkPoint *p, const double dt) {
   if (p->Neumann || p->CORNER) {
      // 接触面候補を多く取得しておく
      std::unordered_set<networkFace *> Fs;
      for (auto &[f, X] : p->getContactFacesX()) {
         auto [p0, p1, p2] = f->getPoints();
         for (auto &F : p0->getFaces())
            Fs.emplace(F);
         for (auto &F : p1->getFaces())
            Fs.emplace(F);
         for (auto &F : p2->getFaces())
            Fs.emplace(F);
      }

      // 接触面候補の次の時刻の位置を予測
      std::vector<T3Tddd> next_Vrtx(Fs.size());
      int i = 0;
      for (auto &f : Fs) {
         auto [p0, p1, p2] = f->getPoints();
         auto X0 = ToX(p0) + f->getNetwork()->velocityRigidBody(ToX(p0)) * dt;
         auto X1 = ToX(p1) + f->getNetwork()->velocityRigidBody(ToX(p1)) * dt;
         auto X2 = ToX(p2) + f->getNetwork()->velocityRigidBody(ToX(p2)) * dt;
         next_Vrtx[i++] = {X0, X1, X2};
      }

      if (next_Vrtx.empty())
         return {{0., 0., 0.}};

      // 各面の次の時刻の配置を予測
      // 最短距離でこの面と向き合う面をnext_Vrtxから探しF_clingsに保存
      auto Xp = nextX_U_Ua_Uc(p, dt);
      std::vector<Tddd> F_clings;
      for (const auto &f : p->getFacesNeumann()) {
         auto n = nextFaceNormal_U_Ua_Uc(f, dt);
         auto Xf = nextX_U_Ua_Uc(f, dt);
         auto V = Xf - Xp;
         auto X = Xp;
         auto to_closest_X = closestXFacing(X, p->radius, next_Vrtx, n);
         if (isFinite(to_closest_X)) {
            auto r = Dot(to_closest_X, n) * n;
            if (isFinite(r))
               F_clings.push_back(r);
         }
      }

      Tddd r = Mean(F_clings);
      if (isFinite(r))
         return F_clings;
      else
         return {{0., 0., 0.}};
   } else
      return {{0., 0., 0.}};
};
//$ ------------------------------------------------------ */
//$ ------------------------------------------------------ */
//$ ------------------------------------------------------ */
double next_length(const networkLine *const l, const double dt) {
   auto [p0, p1] = l->getPoints();
   return Norm(nextX_U_Ua(p0, dt) - nextX_U_Ua(p1, dt));
};
V_d next_length(const std::unordered_set<networkLine *> &ls,
                const double dt) {
   V_d ret;
   for (const auto &l : ls)
      ret.emplace_back(next_length(l, dt));
   return ret;
};

double next_max_length(const networkPoint *const p, const double dt) {
   auto p_next_X = nextX_U_Ua(p, dt);
   int size = 0;
   double max_len = -1;
   for (const auto &q : p->getNeighbors()) {
      auto tmp = Norm(nextX_U_Ua(q, dt) - p_next_X);
      if (max_len < 0 || max_len < tmp)
         max_len = Norm(nextX_U_Ua(q, dt) - p_next_X);
   }
   return max_len;
};

double next_min_length(const networkPoint *const p, const double dt) {
   auto p_next_X = nextX_U_Ua(p, dt);
   int size = 0;
   double min_len = -1, tmp;
   for (const auto &q : p->getNeighbors()) {
      tmp = Norm(nextX_U_Ua(q, dt) - p_next_X);
      if (min_len < 0 || min_len > tmp)
         min_len = Norm(nextX_U_Ua(q, dt) - p_next_X);
   }
   return min_len;
};

double next_mean_length(const networkPoint *const p, const double dt) {
   auto p_next_X = nextX_U_Ua(p, dt);
   int size = 0;
   double mean_len = 0.;
   for (const auto &q : p->getNeighbors()) {
      mean_len += Norm(nextX_U_Ua(q, dt) - p_next_X);
      size++;
   }
   return mean_len / (double)size;
};

V_d next_length_shorter(const networkPoint *const p, const double dt) {
   V_d ret(p->getLines().size());
   int i = 0;
   for (const auto &l : p->getLines())
      ret[i++] = next_length(l, dt);
   std::sort(ret.begin(), ret.end(), [](const auto &a, const auto &b) { return a < b; });
   return ret;
};

std::vector<networkLine *> next_lines_shorter(const networkPoint *const p, const double dt) {
   std::vector<networkLine *> ret = p->getLines();
   std::sort(ret.begin(), ret.end(), [dt](const auto &a, const auto &b) { return next_length(a, dt) < next_length(b, dt); });
   return ret;
};

double next_mean_length_except(const networkPoint *const p, const double dt, networkLine *const L) {
   auto p_next_X = nextX_U_Ua(p, dt);
   int size = 0;
   double mean_len = 0.;
   for (const auto &l : p->getLines()) {
      if (L != l) {
         mean_len += Norm(nextX_U_Ua((*l)(p), dt) - p_next_X);
         size++;
      }
   }
   return mean_len / (double)size;
};

T3Tddd next_vertices(const networkFace *const f, const double dt) {
   auto [p0, p1, p2] = f->getPoints();
   return f->getXVertices() + dt * T3Tddd{nextX_U_Ua(p0, dt), nextX_U_Ua(p1, dt), nextX_U_Ua(p2, dt)};
};

double next_area(const networkFace *const f, const double dt) {
   auto [p0, p1, p2] = f->getPoints();
   return TriangleArea(f->getXVertices() + dt * T3Tddd{nextX_U_Ua(p0, dt), nextX_U_Ua(p1, dt), nextX_U_Ua(p2, dt)});
};

Tddd next_face_center(const networkFace *const f, const double dt) {
   return Mean(next_vertices(f, dt));
};

double next_min_length_to_face(const networkPoint *const p, const double dt) {
   auto p_next_X = nextX_U_Ua(p, dt);
   int size = 0;
   double min_len = -1, tmp;
   for (const auto &f : p->getFaces()) {
      tmp = Norm(next_face_center(f, dt) - p_next_X);
      if (min_len < 0 || min_len > tmp)
         min_len = tmp;
   }
   return min_len;
};

double next_neighbors_mean_Area(const networkPoint *const p, const double dt) {
   double ret = 0.;
   int i = 0;
   for (const auto &F : p->getFaces()) {
      ret += next_area(F, dt);
      i++;
   }
   return ret / i;
};

double next_neighbors_mean_Area(const networkFace *const f, const double dt) {
   double ret = 0.;
   int i = 0;
   for (const auto &F : f->getNeighbors()) {
      ret += next_area(F, dt);
      i++;
   }
   return ret / i;
};

Tddd next_neighbors_center_AreaAveraged(const networkPoint *const p, const double dt) {
   int i = 0;
   Tddd ret = {0., 0., 0.};
   if (p->CORNER) {
      for (const auto &l : p->getLines()) {
         if (l->CORNER) {
            ret += nextX_U_Ua(((*l)(p)), dt);
            i++;
         }
      }
      return ret / (double)i;
   } else {
      double A = 0;
      for (const auto &f : p->getFaces()) {
         double a = next_area(f, dt);
         A += a;
         ret += next_face_center(f, dt) * a;
         i++;
      }
      return ret / A;
   }
};

Tddd next_neighbors_center(const networkPoint *const p, const double dt) {
   double s = 0;
   Tddd ret = {0., 0., 0.};
   for (const auto &l : p->getLines()) {
      if (p->CORNER && l->CORNER) {
         ret += nextX_U_Ua(((*l)(p)), dt);
         s += 1.;
      }
      if (!p->CORNER) {
         for (const auto &l : p->getLines()) {
            ret += nextX_U_Ua(((*l)(p)), dt);
            s += 1.;
         }
      }
   }
   return ret / s;
};

int getMinDepth(const networkPoint *const p) {
   int min_depth = 999999;
   for (const auto &[net, depth] : p->net_depth)
      if (min_depth >= depth)
         min_depth = depth;
   return min_depth;
};

std::tuple<Network *, int> getMinDepthAndNetwork(const networkPoint *const p) {
   std::tuple<Network *, int> ret = {nullptr, 999999};
   for (const auto &[net, depth] : p->net_depth)
      if (std::get<1>(ret) >= depth) {
         std::get<0>(ret) = net;
         std::get<1>(ret) = depth;
      }
   return ret;
};

#include "BEM_derivatives.hpp"

//* ------------------------------------------------------ */
//*                        境界値問題を解く                   */
//* ------------------------------------------------------ */
// #define solve_equations_on_all_points
#define solve_equations_on_all_points_rigid_mode
#define solveBVP_debug

void addIG(std::unordered_map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign) {
   std::unordered_map<netP *, Tdd>::iterator it;
   if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
      std::get<0>(it->second) += std::get<0>(igign);  // phiは忘れずに計算
   else
      P_phiphin[P_IN] = {std::get<0>(igign), 0.};  // phiは忘れずに計算
}

void addIGn(std::unordered_map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign) {
   std::unordered_map<netP *, Tdd>::iterator it;
   if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
      std::get<1>(it->second) += std::get<1>(igign);  // phiは忘れずに計算
   else
      P_phiphin[P_IN] = {0., std::get<1>(igign)};  // phiは忘れずに計算
}

void addIG(std::map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign) {
   std::map<netP *, Tdd>::iterator it;
   if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
      std::get<0>(it->second) += std::get<0>(igign);  // phiは忘れずに計算
   else
      P_phiphin[P_IN] = {std::get<0>(igign), 0.};  // phiは忘れずに計算
}

void addIGn(std::map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign) {
   std::map<netP *, Tdd>::iterator it;
   if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
      std::get<1>(it->second) += std::get<1>(igign);  // phiは忘れずに計算
   else
      P_phiphin[P_IN] = {0., std::get<1>(igign)};  // phiは忘れずに計算
}

#include "BEM_BVP.hpp"

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

auto Pd2PVd = [](map_P_d &P_kappa) {
   map_P_Vd ret;
   for (auto &[p, kappa] : P_kappa)
      ret[p] = {kappa};
   return ret;
};

// V_d laplacianV(const netPp p)
// {
// 	auto ps = Flatten(BFS(p, 2, {p->getNetwork()}));
// 	auto intp = InterpolationVectorRBF(extractX(ps), extVelocities(ps), p->getX());
// 	return intp.laplacian(p->getX());
// };

// V_d meanX(const std::unordered_set<networkPoint *> &ps)
// {
// 	V_d ret(3, 0.);
// 	for (const auto &p : ps)
// 		ret += p->getX();
// 	return ret / ((double)(ps.size()));
// };
void normalizePhi(const std::unordered_set<networkPoint *> &ps) {
   double ret = 0.;
   for (const auto &p : ps)
      ret += std::get<0>(p->phiphin);
   ret /= ((double)(ps.size()));
   for (const auto &p : ps)
      std::get<0>(p->phiphin) -= ret;
};
/* ------------------------------------------------------ */

// b! ------------------------------------------------------ */
// b!           格子のdivide, merge．それに伴うΦ，Φnの付与       */
// b! ------------------------------------------------------ */
Tdd phiphin_from_faces(const networkLine *const l) {
   Tdd phiphin = {0., 0.};
   for (const auto &f : l->getFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      phiphin += (p0->phiphin + p1->phiphin + p2->phiphin) / 6.;
   }
   return phiphin;
};
Tdd phiphin_from_points(const networkLine *const l) {
   Tdd phiphin = {0., 0.};
   auto [p, q] = l->getPoints();
   phiphin += p->phiphin / 2.;
   phiphin += q->phiphin / 2.;
   return phiphin;
};

Tdd phiphin_from_points_faces(const networkLine *const l) {
   return phiphin_from_points(l) / 2. + phiphin_from_faces(l) / 2.;
};

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

bool isNeumann(const networkFace *const f) {
   auto [p0, p1, p2] = f->getPoints();
   auto faces_p0 = p0->getContactFaces();
   auto faces_p1 = p1->getContactFaces();
   auto faces_p2 = p2->getContactFaces();
   bool isNeumann = f->isThereAnyFacingFace(faces_p0, M_PI / 9.) &&
                    f->isThereAnyFacingFace(faces_p1, M_PI / 9.) &&
                    f->isThereAnyFacingFace(faces_p2, M_PI / 9.);
   return isNeumann;
};

void setBoundaryConditions(Network &water, const std::vector<Network *> &objects) {
   auto radius = Mean(extLength(water.getLines()));
   auto Points = ToVector(water.getPoints());
   Print("makeBucketFaces", Green);
   for (const auto &net : objects)
      net->makeBucketFaces(radius);

   // b% -------------------------------------------------------- */
   // b%            境界条件（角点・ディリクレ・ノイマン）の決定         */
   // b% -------------------------------------------------------- */
   // b% step1 点の衝突の判定
   std::cout << "step1 点の衝突の判定" << std::endl;
   for (const auto &p : Points)
      p->clearContactFaces();
   //!!! 衝突の判定がよくエラーが出る箇所
   for (const auto &net : objects) {
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
      {
         //! ここも重要：点と面の衝突をどのようにすれば矛盾なく判定できるか．
         // p->radius = radius / 2.5; // Mean(extLength(p->getLines()));
         auto toF = extXtuple(ToVector(p->getFaces())) - ToX(p);
         auto toP = extXtuple(p->getNeighbors()) - ToX(p);
         double a = Norm(*std::min_element(toP.begin(), toP.end(), [](const auto &a, const auto &b) { return Norm(a) < Norm(b); }));
         double b = Norm(*std::min_element(toF.begin(), toF.end(), [](const auto &a, const auto &b) { return Norm(a) < Norm(b); }));
         p->radius = Mean(extLength(extractLines(Flatten(BFS(p, 2))))) / 5.;
         p->addContactFaces(net->getBucketFaces(), false); /**shadowあり*/
      }
   }
   // b% step2 面の境界条件を判定
   std::cout << "step2 面の境界条件を判定" << std::endl;
   /*面Aの点が接触している面Bを取得．A,B面が向き合っていればノイマン*/
   for (const auto &f : water.getFaces()) {
      f->Neumann = isNeumann(f);  //(!faces_p0.empty()) && (!faces_p1.empty()) && (!faces_p2.empty());
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
   for (const auto &p : Points) {
      p->phinOnFace.clear();
      auto n = p->getNormalTuple();
      if (p->Neumann || p->CORNER) {
         // std::get<1>(p->phiphin) = p->phin_Neumann = Dot(uNeumann(p, p->getNormalNeumann_BEM()), p->getNormalNeumann_BEM()); // Dot(uNeumann(p), p->getNormalNeumann_BEM());
         // std::get<1>(p->phiphin) = Dot(uNeumann(p), p->getNormalNeumann_BEM());
         // p->phinOnFace[nullptr] = std::get<1>(p->phiphin);

         bool use_multiple_nodes = false;
         auto faces = p->getFacesNeumann();
         for (auto i = 0; i < faces.size(); ++i) {
            if (std::abs(Dot(n, faces[i]->normal) - 1) > 1E-5) {
               use_multiple_nodes = true;
               break;
            }
         }
         if (use_multiple_nodes) {
            for (const auto &f : p->getFacesNeumann())
               p->phinOnFace[f] = velocity_normal_from_Neumann_surface(p, f);
         } else {
            std::get<1>(p->phiphin) = Dot(uNeumann(p), p->getNormalNeumann_BEM());
            p->phinOnFace[nullptr] = std::get<1>(p->phiphin);
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
      // else if (p->CORNER)
      // {
      // 	bool use_multiple_nodes = false;
      // 	auto faces = p->getFacesNeumann();
      // 	for (auto i = 0; i < faces.size(); ++i)
      // 	{
      // 		if (std::abs(Dot(n, faces[i]->normal) - 1) > 1E-5)
      // 		{
      // 			use_multiple_nodes = true;
      // 			break;
      // 		}
      // 	}
      // 	if (use_multiple_nodes)
      // 	{
      // 		for (const auto &f : p->getFacesNeumann())
      // 			p->phinOnFace[f] = velocity_normal_from_Neumann_surface(f);
      // 	}
      // 	else
      // 	{
      // 		std::get<1>(p->phiphin) = Dot(uNeumann(p), p->getNormalNeumann_BEM());
      // 		p->phinOnFace[nullptr] = std::get<1>(p->phiphin);
      // 	}
      // }
   }
   // b! 面
   std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << colorOff << std::endl;
   for (const auto &f : water.getFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      auto phiphin = (p0->phiphin + p1->phiphin + p2->phiphin) / 3.;

      // if (f->Neumann)
      // {
      // 	std::get<1>(f->phiphin) = phin_contact(f);
      // 	std::get<0>(f->phinTuple) = Dot(uNeumann(p0, f->normal), f->normal);
      // 	std::get<1>(f->phinTuple) = Dot(uNeumann(p1, f->normal), f->normal);
      // 	std::get<2>(f->phinTuple) = Dot(uNeumann(p2, f->normal), f->normal);
      // }
      // else
      // 	std::get<0>(f->phiphin) = std::get<0>(phiphin);

      if (f->Neumann) {
         std::get<1>(f->phiphin) = Dot(uNeumann(p0), f->normal);
         std::get<1>(f->phiphin) = Dot(uNeumann(p1), f->normal);
         std::get<1>(f->phiphin) = Dot(uNeumann(p2), f->normal);

         // if (p0->CORNER)
         // 	p0->phinOnFace[f] = Dot(uNeumann(p0), f->normal);
         // if (p1->CORNER)
         // 	p1->phinOnFace[f] = Dot(uNeumann(p1), f->normal);
         // if (p2->CORNER)
         // 	p2->phinOnFace[f] = Dot(uNeumann(p2), f->normal);
      } else
         std::get<0>(f->phiphin) = std::get<0>(phiphin);
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
   uomap_P_Tddd P_accel_body, P_NearestContactFacesX,
       P_phin_Neumann, P_phin_Dirichlet, P_velocity_body,
       P_uNeumann, P_normal, P_normal_BEM, P_mirrorPosition, P_U_normal_BEM, P_U_dot_gradgrad_U,
       P_U_tangential_BEM, P_U_BEM, P_U_update_BEM, P_U_cling_to_Neumann, P_phin_vector, P_gradPhi_tangential, P_gradPhi;
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
   P_U_dot_gradgrad_U = P_gradPhi = P_gradPhi_tangential = P_phin_vector = P_accel_body = P_NearestContactFacesX = P_phin_Neumann = P_phin_Dirichlet = P_velocity_body = P_uNeumann = P_normal = P_normal_BEM = P_mirrorPosition = P_U_normal_BEM = P_U_tangential_BEM = P_U_BEM = P_U_update_BEM = P_U_cling_to_Neumann = initial_uomap_P_Tddd;
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
            P_accel_body[p] = p->getNormal_BEM() * accel_normal_from_Neumann_surface(p);  //!
            // P_velocity_body[p] = local_velocity_from_Neumann_surface(p);                  //!
         }
         // auto [m0, s0, min0, smin0, max0] = distorsion(p, dt);
         // P_s_m[p] = s0 / m0;
         // P_smin_min[p] = smin0 / min0;
         P_volume[p] = water.getVolume();
         // P_phin_Neumann[p] = p->getNormalNeumann_BEM() * p->phin_Neumann;
         P_phin_Dirichlet[p] = p->getNormalDirichlet_BEM() * p->phin_Dirichlet;
         // P_phi_Neumann[p] = p->phi_Neumann;
         P_phi_Dirichlet[p] = p->phi_Dirichlet;
         // if (p->Neumann || p->CORNER)
         //    if (std::get<0>(getNearestContactFacesX(p))) {
         //       auto tmp = std::get<1>(getNearestContactFacesX(p)) - ToX(p);
         //       if (Norm(tmp) < p->radius)
         //          P_NearestContactFacesX[p] = tmp;
         //    }
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
         P_height[p] = std::get<2>(p->X);
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
          //  {"accel_body", P_accel_body},
          // {"volume", P_volume},
          {"U_cling_to_Neumann", P_U_cling_to_Neumann},
          //  {"absU_cling_to_Neumann/absU_update_BEM", P_update_vs_cling},
          //  {"velocity_body", P_velocity_body},
          //  {"NearestContactFacesX", P_NearestContactFacesX},
          //  {"φn_Neumann", P_phin_Neumann},
          //  {"φn_Dirichlet", P_phin_Dirichlet},
          //  {"φ_Neumann", P_phi_Neumann},
          //  {"φ_Dirichlet", P_phi_Dirichlet},
          //  {"min_depth", P_min_depth},
          //  {"isMultipleNode", P_isMultipleNode},
          // P_NearestContactFacesXで確かに最寄の構造物までをさすことができている．！
          // これで改善できない？
          //  {"U_update_BEM", P_U_update_BEM},
          {"U_BEM", P_U_BEM},
          //  {"U_normal_BEM", P_U_normal_BEM},
          {"U_tangential_BEM", P_U_tangential_BEM},
          // {"kappa", ders.P_kappa},
          //  {"ContactFaces", P_ContactFaces},
          // {"dxdt_mod", P_dxdt_mod},
          {"grad_phi", P_gradPhi},
          //  {"phin_vector", P_phin_vector},
          {"gradPhiTangential", P_gradPhi_tangential},
          {"z", P_height},
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
   calculateFroudeKrylovForce(const std::unordered_set<networkFace *> faces /*waterfaces*/, const Network *PasObj)
       : force({0., 0., 0.}), torque({0., 0., 0.}), area(0.), PressureVeticies({}), acceleration({0., 0., 0., 0., 0., 0.}) {
      // PasObjと接したfaceの頂点にpressureが設定されている前提
      Tddd pressure012 = {1E+60, 1E+60, 1E+60};
      for (const auto &f : faces)
         if (f->Neumann) {
            pressure012 = {1E+60, 1E+60, 1E+60};
            auto [p0, p1, p2] = f->getPoints();
            for (const auto &f0 : p0->getContactFaces())
               if (f0->getNetwork() == PasObj)
                  std::get<0>(pressure012) = p0->pressure;
            for (const auto &f1 : p1->getContactFaces())
               if (f1->getNetwork() == PasObj)
                  std::get<1>(pressure012) = p1->pressure;
            for (const auto &f2 : p2->getContactFaces())
               if (f2->getNetwork() == PasObj)
                  std::get<2>(pressure012) = p2->pressure;
            if (isFinite(pressure012)) {
               this->PressureVeticies.emplace_back(std::tuple<Tddd, T3Tddd>{pressure012, T3Tddd{ToX(p0), ToX(p1), ToX(p2)}});
               this->actingFaces.emplace_back(f);
            }
         }
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpX = interpolationTriangleLinear0101(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
            area += intpX.J(x0, x1) * w0w1;
      }
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
         for (const auto &[x0, x1, w0w1] : __GWGW14__Tuple)
            this->torque += Cross(intpX(x0, x1) - COM, TriangleNormal(X012) * intpP(x0, x1)) * intpX.J(x0, x1) * w0w1;
      }
      return this->torque;
   };

   Tddd getFroudeKrylovForce() {
      this->force = {0., 0., 0.};
      for (const auto &[P012, X012] : this->PressureVeticies) {
         auto intpP = interpolationTriangleLinear0101(P012);
         auto intpX = interpolationTriangleLinear0101(X012);
         for (const auto &[x0, x1, w0w1] : __GWGW14__Tuple)
            this->force += TriangleNormal(X012) * intpP(x0, x1) * intpX.J(x0, x1) * w0w1;
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

int main() {
   try {
      //* ------------------------------------------------------ */
      //*                         setting                        */
      //* ------------------------------------------------------ */
      JSON settingJSON(std::ifstream("./inputs/setting.json"));
      std::map<Network *, outputInfo> NetOutputInfo;
      std::vector<Network *> FluidObject, RigidBodyObject, SoftBodyObject;
      std::string output_directory = settingJSON["output_directory"][0];
      double max_dt = stod(settingJSON["max_dt"])[0];
      std::filesystem::create_directories(output_directory);
      /* ------------------------------------------------------ */
      double stop_remesh_time = 1E+10;
      double force_remesh_time = 0;
      double preparation_time = stod(settingJSON["preparation_time"])[0];
      double preparation_max_dt = stod(settingJSON["preparation_max_dt"])[0];
      if (settingJSON.find("stop_remesh_time"))
         stop_remesh_time = stod(settingJSON["stop_remesh_time"][0]);

      if (settingJSON.find("force_remesh_time"))
         force_remesh_time = stod(settingJSON["force_remesh_time"][0]);

      int grid_refinement = 0;
      if (settingJSON.find("grid_refinement"))
         grid_refinement = stoi(settingJSON["grid_refinement"][0]);
      /* ------------------------------------------------------ */

      for (auto FileName : settingJSON["inputfiles"]) {
         Print("---------------------------------------------------");
         Print("-----------------" + FileName + "------------------");
         Print("---------------------------------------------------");
         JSON J(std::ifstream("./inputs/" + FileName));
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
            NetOutputInfo[net].pvd_file_name = J["output_pvd_file_name"][0];
            NetOutputInfo[net].vtu_file_name = J["output_vtu_file_name"][0];
            NetOutputInfo[net].PVD = new PVDWriter(output_directory + "/" + J["output_pvd_file_name"][0] + ".pvd");
            std::filesystem::copy_file("./inputs/" + FileName, output_directory + "/" + FileName, std::filesystem::copy_options::overwrite_existing);
            if (J["type"][0] == "RigidBody") {
               RigidBodyObject.emplace_back(net);
               net->isRigidBody = true;
               net->isSoftBody = false;
            } else if (J["type"][0] == "SoftBody") {
               SoftBodyObject.emplace_back(net);
               net->isRigidBody = false;
               net->isSoftBody = true;
            } else if (J["type"][0] == "Fluid") {
               FluidObject.emplace_back(net);
               net->isRigidBody = net->isSoftBody = false;
            }
            //
            mk_vtu(output_directory + "/" + J["name"][0] + "_init.vtu", {net->getFaces()});
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
      //  b@ ------------------------------------------------------ */
      //  b@                     グリッドの調整                         */
      //  b@ ------------------------------------------------------ */
      // if (true)
      // {
      // 	setBoundaryConditions(water, RigidBodyObject);
      // 	for (const auto &p : water.getPoints())
      // 	{
      // 		double h = 0.5;
      // 		if (!p->Neumann)
      // 		{
      // 			auto tmp = ToX(p);
      // 			p->setXSingle(Tddd{std::get<0>(tmp), std::get<1>(tmp), h});
      // 		}
      // 	}
      // 	water.setGeometricProperties();
      // }
      //
      {
         PVDWriter preRefinementPVD(output_directory + "/preRefinement.pvd");
         double dt = 1E-6;
         double min_dt = 1E+10;
         for (auto i = 0; i < grid_refinement; i++) {
            std::cout << "refinement i = " << i << std::endl;
            // b% -------------------------------------------------------- */
            // b%            境界条件（角点・ディリクレ・ノイマン）の決定               */
            // b% -------------------------------------------------------- */
            setBoundaryConditions(*water, RigidBodyObject);
            auto data = dataForOutput(*water, dt);
            auto filename = output_directory + "/preRefinement" + std::to_string(i) + ".vtu";
            mk_vtu(filename, {water->getFaces()}, data);
            preRefinementPVD.push(filename, (double)i);
            preRefinementPVD.output();
            std::string name = output_directory + "/water100_refined_" + std::to_string(i) + ".obj";
            // std::string name = output_directory + "/water30_refined_laplacian_0d05_" + std::to_string(i) + ".obj";
            std::ofstream ofs(name);
            std::cout << "name = " << name << std::endl;
            creteOBJ(ofs, *water);
            ofs.close();
            // b* ------------------------------------------------------ */
            // b*                    微分∇ΦやDUDtを計算                    */
            // b* ------------------------------------------------------ */
            std::cout << Green << "微分∇ΦやDUDtを計算" << colorOff << std::endl;
            int loop = 2;
            double rad = M_PI / 180;
            for (auto ii = 0; ii <= loop; ++ii) {
               min_dt = dt_CFL(*water, min_dt);
               if (min_dt < preparation_max_dt)
                  dt = min_dt;
               else
                  dt = preparation_max_dt;

               int RK_order = 3;
               std::map<netPp, RungeKutta<Tddd> *> P_RK_X;
               auto Points = water->getPoints();
               for (const auto &p : Points)
                  P_RK_X[p] = (new RungeKutta(dt, real_time, ToX(p), RK_order));

               do {
                  derivatives ders(*water, true);
                  for (const auto &p : Points) {
                     auto rk = P_RK_X[p];
                     rk->push(p->U_update_BEM * (i < 30 ? 0.1 : 1.));
                     // rk->push(ders.P_dxdt[p]);
                     p->setXSingle(rk->getX());
                  }
                  water->setGeometricProperties();
               } while (!(P_RK_X[*Points.begin()]->finished));

               for (const auto &p : water->getPoints())
                  delete P_RK_X[p];

               // flipIf(*water, {15 * rad, rad}, {10 * rad /*結構小さく*/, rad}, false);
               flipIf(*water, {5 * rad, rad}, {5 * rad /*結構小さく*/, rad}, false);
               // if (time_step > 1 && time_step < 10)
               // {
               // 	flipIf(water, {8 * rad, 30 * rad}, {8 * rad /*結構小さく*/, 10 * rad}, true, 5 /*やはり制限はひつようか？*/);
               // 	flipIf(water, {8 * rad, 10 * rad}, {8 * rad /*結構小さく*/, 30 * rad}, true, 5 /*やはり制限はひつようか？*/);
               // }

               // flipIf(water, {rad, rad}, {10 * rad /*結構小さく*/, rad}, false);
               // if (i > 5)
               // 	if (ii == 0)
               // 	{
               // 		flipIf(water, {10 * rad, 30 * rad}, {10 * rad /*結構小さく*/, 10 * rad}, true, 10 /*やはり制限はひつようか？*/);
               // 		flipIf(water, {10 * rad, rad}, {10 * rad /*結構小さく*/, 30 * rad}, true, 10 /*やはり制限はひつようか？*/);
               // 		// arrangeCORNER(water, M_PI / 180, 10);
               // 		// 	arrangeCORNER(water, rad / 100., 10); //このアレンジで崩れることもあった．
               // 		// }else if (ii == 2)
               // 		// {
               // 		// 	flipIf(water, {rad, 30 * rad}, {5 * rad /*結構小さく*/, 5 * rad}, true, 2 /*やはり制限はひつようか？*/);
               // 		// 	flipIf(water, {rad, 5 * rad}, {5 * rad /*結構小さく*/, 30 * rad}, true, 2 /*やはり制限はひつようか？*/);
               // 	}
               // flipIf(water, {10 * rad, rad}, {10 * rad /*結構小さく*/, rad}, false);
            }
            // if (i % 2 == 0)
            // 	remesh(water, {10 * rad, 10 * rad}, {10 * rad /*結構小さく*/, 10 * rad}, false);
         };
         // 内部の角点に常にいるようにするには？？？？？？？？？？？？？・
         // 面が吸い付く必要がある
         std::ofstream ofs(output_directory + "/water_refined.obj");
         creteOBJ(ofs, *water);
         ofs.close();
      }
      //  b* ------------------------------------------------------ */
      //  b*                         メインループ                     */
      //  b* ------------------------------------------------------ */
      TimeWatch watch;
      for (time_step = 0; time_step < 1200; time_step++) {
         show_info(*water);
         //! 体積を保存するようにリメッシュする必要があるだろう．
         // auto radius = Mean(extLength(water->getLines()));
         setBoundaryConditions(*water, RigidBodyObject);
         // b* ------------------------------------------------------ */
         // if (real_time < stop_remesh_time && time_step > 0)
         // {
         // 	double rad = M_PI / 180;
         // 	flipIf(*water, {15 * rad, rad}, {15 * rad /*結構小さく*/, rad}, false);
         // 	if (time_step > 1 && time_step < 10)
         // 	{
         // 		flipIf(*water, {10 * rad, 30 * rad}, {10 * rad /*結構小さく*/, 20 * rad}, true, 20 /*やはり制限はひつようか？*/);
         // 		flipIf(*water, {10 * rad, 10 * rad}, {10 * rad /*結構小さく*/, 30 * rad}, true, 20 /*やはり制限はひつようか？*/);
         // 	}
         // }
         double rad = M_PI / 180;
         flipIf(*water, {5 * rad, rad}, {5 * rad /*結構小さく*/, rad}, false);
         // flipIf(*water, {10 * rad, 30 * rad}, {10 * rad /*結構小さく*/, 20 * rad}, true, 20 /*やはり制限はひつようか？*/);
         // flipIf(*water, {10 * rad, 10 * rad}, {10 * rad /*結構小さく*/, 30 * rad}, true, 20 /*やはり制限はひつようか？*/);
         // b# ------------------------------------------------------ */
         // b#                       刻み時間の決定                     */
         // b# ------------------------------------------------------ */
         const auto Points = water->getPoints();
         const auto Faces = water->getFaces();

         // double min_dt = 1E+10;
         // double dt = dt_CFL(*water, max_dt, {0.2, 1});
         double dt = max_dt;
         // if (2. < real_time && real_time < 4.)
         // 	dt = 0.005;

         // if (real_time <= preparation_time)
         // 	dt = preparation_max_dt;
         // else if (min_dt < max_dt)
         // 	dt = min_dt;
         // else
         // 	dt = max_dt;
         Print("===========================================================================");
         Print("       dt :" + Red + std::to_string(dt) + colorOff);
         Print("time_step :" + Red + std::to_string(time_step) + colorOff);
         Print("real time :" + Red + std::to_string(real_time) + colorOff);
         Print("---------------------------------------------------------------------------");

         jsonout.push("time", real_time);
         jsonout.push(water->getName() + "_volume", water->getVolume());
         for (const auto net : RigidBodyObject) {
            auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
            auto [mx, my, mz, Ix, Iy, Iz] = net->getInertiaGC();
            auto force = tmp.getFroudeKrylovForce() + _GRAVITY3_ * net->mass;
            auto torque = tmp.getFroudeKrylovTorque(net->COM);
            auto [a0, a1, a2] = force / Tddd{mx, my, mz};
            auto [a3, a4, a5] = torque / Tddd{Ix, Iy, Iz};
            auto acceleration = T6d{a0, a1, a2, a3, a4, a5};
            std::cout << red << net->getName() << "\n"
                      << ", net->getInertiaGC() = " << net->getInertiaGC() << "\n"
                      << ", net->mass = " << net->mass << "\n"
                      << ", tmp.getFroudeKrylovForce() = " << tmp.getFroudeKrylovForce() << "\n"
                      << ", tmp.getFroudeKrylovTorque(net->COM) = " << tmp.getFroudeKrylovTorque(net->COM) << "\n"
                      << ", tmp.area = " << tmp.area << "\n"
                      << ", net->velocity = " << net->velocity << "\n"
                      << ", net->acceleration = " << acceleration << "\n"
                      << ", net->COM = " << net->COM << "\n"
                      << colorOff << std::endl;
            /* ------------------------------------------------------ */
            jsonout.push(net->getName() + "_FK_force", tmp.getFroudeKrylovForce());
            jsonout.push(net->getName() + "_FK_torque", tmp.getFroudeKrylovTorque(net->COM));
            jsonout.push(net->getName() + "_accel", acceleration);
            jsonout.push(net->getName() + "_area", tmp.area);
         }

         double maxUcling = 0;
         for (const auto &p : water->getPoints())
            if (p->Neumann || p->CORNER) {
               auto tmp = Norm(p->U_cling_to_Neumann) / Norm(p->U_update_BEM);
               if (maxUcling < tmp && isFinite(tmp))
                  maxUcling = tmp;
            }
         jsonout.push("maxUcling", maxUcling);
         double maxUcling_minlength = 0;
         for (const auto &p : water->getPoints())
            if (p->Neumann || p->CORNER) {
               double min = Min(extLength(p->getLines()));
               auto tmp = Norm(p->U_cling_to_Neumann) * dt / min;
               if (maxUcling_minlength < tmp && isFinite(tmp))
                  maxUcling_minlength = tmp;
            }
         jsonout.push("maxUcling_minlength", maxUcling_minlength);

         double maxUcling_meanlength = 0;
         for (const auto &p : water->getPoints())
            if (p->Neumann || p->CORNER) {
               double mean = Mean(extLength(p->getLines()));
               auto tmp = Norm(p->U_cling_to_Neumann) * dt / mean;
               if (maxUcling_meanlength < tmp && isFinite(tmp))
                  maxUcling_meanlength = tmp;
            }
         jsonout.push("maxUcling_meanlength", maxUcling_meanlength);

         std::ofstream os(output_directory + "/result.json");
         jsonout.output(os);
         os.close();
         /* ------------------------------------------------------ */

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
            for (const auto &p : net->getPoints())
               p->RK_X.initialize(dt, real_time, ToX(p), RK_order);
         }
         /* ------------------------------------ */
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
            setBoundaryConditions(*water, RigidBodyObject);
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
            derivatives ders(*water, time_step > 1);
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
                  auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
                  auto [mx, my, mz, Ix, Iy, Iz] = net->getInertiaGC();
                  auto force = tmp.getFroudeKrylovForce() + _GRAVITY3_ * net->mass;
                  auto torque = tmp.getFroudeKrylovTorque(net->COM);
                  auto [a0, a1, a2] = force / Tddd{mx, my, mz};
                  auto [a3, a4, a5] = torque / Tddd{Ix, Iy, Iz};
                  net->acceleration = T6d{a0, a1, a2, a3, a4, a5};
                  std::cout << red << net->getName() << "\n"
                            << ", net->getInertiaGC() = " << net->getInertiaGC() << "\n"
                            << ", net->mass = " << net->mass << "\n"
                            << ", tmp.getFroudeKrylovForce() = " << tmp.getFroudeKrylovForce() << "\n"
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
                  /* ------------------------- */
                  Quaternion q;
                  q = q.d_dt(net->velocityRotational());  // w->クォータニオン
                  net->RK_Q.push(q());                    // クォータニオン->T4dとしてプッシュ
                  net->Q = net->RK_Q.getX();
                  /* ------------------------- */
                  net->RK_Velocity.push(net->acceleration);
                  net->velocity = net->RK_Velocity.getX();
                  /* ------------------------- */
                  net->RigidBodyMovePoints();
               }
            }
            // b$ --------------------------------------------------- */
            for (const auto &net : SoftBodyObject) {
               std::cout << "name:" << net->getName() << std::endl;
               if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] != "fixed") {
                  for (const auto &p : net->getPoints()) {
                     p->RK_X.push(p->velocityTranslational());  //@ 位置xの時間発展
                     p->setXSingle(p->RK_X.getX());
                  }
               }
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

         /* ------------------------------------------------------ */
         std::cout << Green << "real_timeを取得" << colorOff << std::endl;
         real_time = (*Points.begin())->RK_X.gett();

         for (const auto &p : Points) {
            p->U_BEM_last = p->U_BEM;
            p->U_tangential_BEM_last = p->U_tangential_BEM;
         }

         /* ------------------------------------------------------ */
         /* ------------------------------------------------------ */
         // for (const auto &p : Points)
         // 	p->X_BUFFER = ToX(p);

         // calculateVectorToSurfaceInBuffer(*water, 0.);

         // for (const auto &p : Points)
         // 	p->calculateBufferPotentialsOnClungSurface();

         // for (const auto &p : Points)
         // 	p->setX(ToX(p) + p->U_BUFFER);

         // for (const auto &p : Points)
         // 	p->copyPotentialsBuffer();

         // water->setGeometricProperties();
         /* ------------------------------------------------------ */
         /* ------------------------------------------------------ */

         // normalizePhi(wavemaker->getPoints());
         //
         std::ofstream ofs(output_directory + "/water_current.obj");
         creteOBJ(ofs, *water);
         ofs.close();

         // b* ------------------------- 出力 ------------------------- */
         auto data = dataForOutput(*water, dt);
         // 流体
         for (const auto &net : FluidObject) {
            auto filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, net->getFaces(), data);
            NetOutputInfo[net].PVD->push(filename, real_time);
            NetOutputInfo[net].PVD->output();
         }
         for (const auto &net : RigidBodyObject) {
            VV_VarForOutput data;
            if (net->inputJSON.find("velocity") && net->inputJSON["velocity"][0] == "floating") {
               auto tmp = calculateFroudeKrylovForce(water->getFaces(), net);
               uomap_P_Tddd P_COM, P_COM_p, P_accel, P_velocity, P_rotational_velocity, P_rotational_accel, P_FroudeKrylovTorque;
               uomap_P_d P_Pitch, P_Yaw, P_Roll;
               for (const auto &p : net->getPoints()) {
                  P_COM[p] = net->COM;
                  P_COM_p[p] = net->COM - ToX(p);
                  P_Pitch[p] = net->quaternion.pitch();
                  P_Yaw[p] = net->quaternion.yaw();
                  P_Roll[p] = net->quaternion.roll();
                  P_accel[p] = net->velocityRigidBody(ToX(p));
                  P_velocity[p] = net->accelRigidBody(ToX(p));
                  P_rotational_velocity[p] = net->velocityRotational();
                  P_FroudeKrylovTorque[p] = tmp.getFroudeKrylovTorque(net->COM);
               }
               data = {{"vector to COM", P_COM_p},
                       {"COM", P_COM},
                       {"pitch", P_Pitch},
                       {"yaw", P_Yaw},
                       {"roll", P_Roll},
                       {"velocity", P_accel},
                       {"acceleration", P_velocity},
                       {"rotational valocity", P_rotational_velocity},
                       {"rotational acceleration", P_rotational_accel},
                       {"FroudeKrylovTorque", P_FroudeKrylovTorque}};
               mk_vtu(output_directory + "/actingFacesOn" + NetOutputInfo[net].vtu_file_name + ".vtu", tmp.actingFaces);
            }
            auto filename = NetOutputInfo[net].vtu_file_name + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, net->getFaces(), data);
            NetOutputInfo[net].PVD->push(filename, real_time);
            NetOutputInfo[net].PVD->output();
         }
         // 流体
         {
            std::vector<Tddd> points;
            for (const auto &p : water->getPoints())
               if (p->CORNER)
                  points.emplace_back(ToX(p));

            auto filename = "cornerPoints" + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, points);
            cornerPointsPVD.push(filename, real_time);
            cornerPointsPVD.output();
         }
         {
            std::unordered_set<networkFace *> faces;
            for (const auto &f : water->getFaces()) {
               // for (const auto p : f->getPoints())
               // 	if (p->CORNER)
               // 	{
               // 		faces.emplace(f);
               // 		break;
               // 	}
               for_each(f->getPoints(), [&](const auto &p) {
                  if (p->CORNER)
                     faces.emplace(f);
               });
            }
            auto filename = "corner" + std::to_string(time_step) + ".vtu";
            mk_vtu(output_directory + "/" + filename, faces, data);
            cornerPVD.push(filename, real_time);
            cornerPVD.output();
         }
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
