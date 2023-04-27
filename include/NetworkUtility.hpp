#ifndef NetworkUtility_H
#define NetworkUtility_H

#include "InterpolationRBF.hpp"
#include "Network.hpp"

/* -------------------------------------------------------------------------- */

void creteOBJ(std::ofstream &ofs, Network &net) {
   std::map<netPp, int> P_i;
   int i = 0;
   for (const auto &p : net.getPoints()) {
      P_i[p] = ++i;
      auto [X0, X1, X2] = ToX(p);
      ofs << "v "
          << X0 << " "
          << X1 << " "
          << X2 << std::endl;
   }
   ofs << std::endl;

   for (const auto &f : net.getFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      ofs << "f "
          << P_i[p0] << " "
          << P_i[p1] << " "
          << P_i[p2] << std::endl;
   }
};
// //% ------------------------------------------------------ */
// //%                         反射点の計算                     */
// //% ------------------------------------------------------ */
bool isInConvexPolygon(const networkPoint *const p) {
   try {
      auto ps = p->getNeighborsSort();  // ソートできない場合エラーとなる
      if (ps.size() < 2)
         return false;
      if (ps.size() == 3)
         return true;
      if (ps.size() != p->getNeighbors().size())
         return false;  // 面が重複する場合を覗くために設置

      return isConvexPolygon(extX(ps), p->getNormalTuple());
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
///////////////////////////////////////////////////////////
bool isStraight(const V_d &v0, const V_d &v1, const double angle = 1E-2 /*閾値*/) {
   //        ^
   //      v1|angle (positive)
   // ---v0-->----
   //      v1|angle (negative)
   //        v
   if (std::abs(MyVectorAngle(v0, v1)) < angle)
      return true;
   return false;
};
///////////////////////////////////////////////////////////
bool isStraight(const netPp p0, const netPp p1, const netPp p2, const double angle = 1E-2 /*閾値*/) {
   return (isStraight(ToVector(ToX(p0) - ToX(p1)),
                      ToVector(ToX(p1) - ToX(p2))) < angle);
};
///////////////////////////////////////////////////////////
bool isFlat(const netPp p, double minangle = M_PI / 180.) {
   if (p->getLines().empty()) {
      // mk_vtu("./vtu/p->getLines().empty().vtu", {{p}});
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty()");
   }

   std::vector<Tddd> normals;
   for (const auto &f : p->getFaces())
      normals.emplace_back(f->normal);
   // 面の法線方向が作る内角が全て1E-2より小さい場合，flat
   for (auto i = 0; i < normals.size(); i++)
      for (auto j = i + 1; j < normals.size(); j++)
         if (Dot(normals[i], normals[j]) < cos(minangle))
            return false;  // not flat
   return true;
};
///////////
bool isFlat(const netLp line, double minangle = M_PI / 180.) {
   auto fs = line->getFaces();
   if (fs.size() != 2)
      return false;
   else if (Dot(fs[0]->normal, fs[1]->normal) < cos(minangle))  // dotはflatなら大きくなる
      return false;
   else
      return true;
};
///////////////////////////////////////////////////////////
Tddd AreaWeightedSmoothingVector(const netPp p) {
   if (!isEdgePoint(p)) {
      Tddd ret = {0., 0., 0.}, p_next = ToX(p), V = {0., 0, 0.};
      T3Tddd t3tdd;
      double Wtot, W;
      // auto lines = p->getLines();
      /* ------------------------------------------- */
      // auto update = [&]() {
      //    Wtot = 0.;
      //    V = {0., 0., 0.};
      //    for (const auto &f : p->getFaces()) {
      //       auto [p0, p1, p2] = f->getPoints(p);
      //       F = {p_next, ToX(p1), ToX(p2)};
      //       W = TriangleArea(F);
      //       Wtot += W;
      //       V += W * (Mean(F) - p_next);
      //    }
      //    V /= Wtot;
      //    p_next += V / 2.;  // このp_nextの位置で引っ張られるかチェック
      // };
      //
      auto update = [&]() {
         Wtot = 0.;
         V = {0., 0., 0.};
         for (const auto &f : p->getFaces()) {
            auto [p0, p1, p2] = f->getPoints(p);
            t3tdd = {p_next, ToX(p1), ToX(p2)};
            // W = TriangleArea(t3tdd) * std::pow(std::log(CircumArea(t3tdd) / InArea(t3tdd)), 2);
            auto inarea = InArea(t3tdd);
            auto carea = CircumArea(t3tdd);
            auto W = TriangleArea(t3tdd);
            //     W *= std::log10(carea / inarea);
            // if (!isFinite(carea) || std::log10(carea / inarea) >= 10.)
            //    W = TriangleArea(t3tdd) * 10.;
            //
            // auto r_i = Inradius(t3tdd);
            // auto r_c = Circumradius(t3tdd);
            // auto a = std::log(r_c / r_i);
            // if (!isFinite(a, 1000.))
            //    a = 1000.;
            // auto W = TriangleArea(t3tdd) * a;
            // auto V = Incenter(t3tdd) - p_next;
            // auto b = r_i / Norm(V);
            // if (!Between(b, {0., 1.}))
            //    b = 1.;
            // V += (1. - b) * V * W;
            //
            // W = TriangleArea(t3tdd) * std::pow(TriangleArea(t3tdd) / InArea(t3tdd), 2);
            Wtot += W;
            // V += W * (Incenter(t3tdd) - p_next);
            // V += W * (Incenter(t3tdd) - p_next) - Inradius(t3tdd) * Normalize((Incenter(t3tdd) - p_next));
            // V += W * (Circumcenter(t3tdd) - p_next) / 3.;
            V += W * (Centroid(t3tdd) - p_next);
         }
         V /= Wtot;
         p_next += V / 4.;  // このp_nextの位置で引っ張られるかチェック
      };
      /* -------------------------------------------- */
      for (auto k = 0; k < 10; ++k) {
         update();
         /* ------------------------------------------------------ */
         // Wtot = 0.;
         // V = {0., 0., 0.};
         // for (const auto &l : lines)
         // {
         // 	auto fs = l->getFaces();
         // 	if (fs.size() > 1)
         // 	{
         // 		W = 1 + Norm(Cross(fs[0]->normal, fs[1]->normal));
         // 		Wtot += W;
         // 		V += W * ((*l)(p)->getXtuple() - p_next);
         // 	}
         // }
         // V /= Wtot;
         // p_next += V / 4.; //このp_nextの位置で引っ張られるかチェック
      }
      return p_next - ToX(p);
   } else
      return {0., 0., 0.};
};
//! -------------------------------------------------------------------------- */

double magicalValue_(const networkPoint *p, const networkFace *f) {
   // auto [p0, p1, p2] = f->getPoints(p);
   auto nP012 = ToX(f->getPoints(p));
   auto tmp = std::log10(Circumradius(nP012) / Inradius(nP012));
   double max = 10;
   if (tmp > max)
      return max;
   else
      return tmp;
};

double variance2_(const networkPoint *p, const networkFace *f) {
   auto [p0, p1, p2] = f->getPoints(p);
   double m = 0, l0, l1, l2;
   m += (l0 = Norm(ToX(p0) - ToX(p1)));
   m += (l1 = Norm(ToX(p1) - ToX(p2)));
   m += (l2 = Norm(ToX(p2) - ToX(p0)));
   m /= 3.;
   auto ret = std::pow(l0 - m, 2) + std::pow(l1 - m, 2) + std::pow(l2 - m, 2);
   return 10 * ret / (m * m);
};

Tddd vectorTangentialShift_(const networkPoint *p, const double scale = 1.) {
   Tddd vectorToCenter = {0., 0., 0.};
   double s = 0;
   Tddd pX = p->X;
   double min_distance = 1E+10;
   if (p->CORNER) {
      for (const auto &l : p->getLinesCORNER()) {
         vectorToCenter += ((*l)(p))->X - pX;
         s += 1.;
      }
      vectorToCenter /= s;
      return scale * vectorToCenter;
   } else {
      for (const auto &f : p->getFaces()) {
         auto nP012 = ToX(f->getPoints(p));
         auto [np0x, np1x, np2x] = nP012;
         auto mean = (Norm(np1x - np0x) + Norm(np2x - np1x) + Norm(np0x - np2x)) / 3.;
         auto l12 = Norm(np2x - np1x);
         double a = magicalValue_(p, f) + variance2_(p, f);
         if (std::ranges::any_of(f->getLines(), [&](const auto &l) { return l->CORNER; }))
            a *= 4;
         else if (std::ranges::any_of(f->getPoints(), [&](const auto &p) { return p->CORNER; }))
            a *= 3;
         else if (std::ranges::any_of(f->getPoints(),
                                      [&](const auto &p) { return std::ranges::any_of(p->getFaces(),
                                                                                      [&](const auto &F) { return std::ranges::any_of(F->getLines(), [&](const auto &l) { return l->CORNER; }); }); }))
            a *= 2;
         auto n = Normalize(Chop(np0x - np1x, np2x - np1x));
         auto X = Norm(np2x - np1x) * n * sin(M_PI / 3.) + (np2x + np1x) / 2.;
         vectorToCenter += a * (X - np0x);  //[m]
         s += a;
      }
      vectorToCenter /= s;
      return scale * vectorToCenter;
   }
   return vectorToCenter;
};

/* ------------------------------------------------------ */
// void LaplacianSmoothing(netPp p) {
//    try {
//       if (!isEdgePoint(p) && isInConvexPolygon(p)) { /*端の点はsmoothingしない*/
//          auto ps = p->getNeighbors();
//          if (ps.size() > 2)  // 2点の場合は2点の中点に動いてしまうので，実行しない
//             p->setX(Mean(extractX(ps)));
//       }
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorOff << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };
// void LaplacianSmoothing(const V_netPp &ps) {
//    try {
//       for (const auto &p : ps)
//          LaplacianSmoothing(p);
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorOff << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };
///////////////////////////////////////////////////////////
Tddd exact_along_surface(const networkPoint *const p, Tddd VECTOR) {
   auto faces = ToVector(p->getFaces());
   std::vector<std::tuple<Tddd, Tddd>> normals;
   normals.reserve(2 * faces.size());
   for (auto i = 0; i < faces.size(); i++)
      for (auto j = i; j < faces.size(); j++)
         if (!isFlat(faces[i]->normal, faces[j]->normal, 1E-2 * M_PI / 180.))
            normals.push_back({faces[i]->normal, faces[j]->normal});
   //
   // std::vector<std::tuple<Tddd, Tddd>> normals;
   // for (const auto &f : p->getFaces())
   // 	for (const auto &F : p->getFaces())
   // 		if (!isFlat(f->normal, F->normal, 1E-2 * M_PI / 180.))
   // 			normals.push_back({f->normal, F->normal});
   Tddd c;
   for (const auto &[N0, N1] : normals) {
      c = Normalize(Cross(N0, N1));
      VECTOR = Dot(VECTOR, c) * c;
   }
   return VECTOR;
};
/* ------------------------------------------------------ */
void AreaWeightedSmoothingPreserveShape(netPp p, const double lim_rad = 1E-10) {
   /*
   AreaWeightedSmoothingVectorの設定に伴って移動が変わる
   */
   try {
      if (!isEdgePoint(p)) {
         auto ps = p->getNeighbors();
         bool isflat = true;
         auto V = AreaWeightedSmoothingVector(p);
         Tddd N, N0, N1;
         bool allflat = true;
         bool isfiniteangles, isfiniteareas, isfinitenormal;
         auto faces = ToVector(p->getFaces());
         for (size_t i = 0; i < faces.size(); i++) {
            auto [p0, p1, p2] = faces[i]->getPoints(p);
            auto triangle = T3Tddd{ToX(p0) + V, ToX(p1), ToX(p2)};
            N0 = TriangleNormal(triangle);
            for (size_t j = i; j < faces.size(); j++) {
               auto [p0, p1, p2] = faces[j]->getPoints(p);
               auto triangle = T3Tddd{ToX(p0) + V, ToX(p1), ToX(p2)};
               N1 = TriangleNormal(triangle);
               if (!isFlat(N0, N1, 1E-2 * M_PI / 180.)) {
                  allflat = false;
                  break;
               }
            }
         }
         // もし面が全てフラットなら調べる必要はない
         for (auto kk = 0; kk < 5; ++kk) {
            if (ps.size() > 2) {
               V = exact_along_surface(p, V);
               for (const auto &f : faces) {
                  auto [p0, p1, p2] = f->getPoints(p);
                  auto triangle = T3Tddd{ToX(p0) + V, ToX(p1), ToX(p2)};
                  isfiniteangles = isFinite(TriangleAngles(triangle));
                  isfiniteareas = isFinite(TriangleArea(triangle));
                  N = TriangleNormal(triangle);
                  isfinitenormal = isFinite(N);
                  if (allflat) {
                     // 全てフラットなので，反転しなければいい．
                     if (!isfiniteangles || !isfiniteareas || !isfinitenormal || !isFlat(N, f->normal, M_PI / 180.)) {
                        isflat = false;
                        break;
                     }
                  } else {
                     if (!isfiniteangles || !isfiniteareas || !isfinitenormal || !isFlat(N, f->normal, lim_rad)) {
                        isflat = false;
                        break;
                     }
                  }
               }
               if (isflat) {
                  p->setX(ToX(p) + V);
                  return;
               }
            }
            V /= 20.;
         }
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
void AreaWeightedSmoothingPreserveShape(const V_netPp &ps /*copy*/, const double lim_rad = 1E-10) {
   for (const auto &p : ps)
      AreaWeightedSmoothingPreserveShape(p, lim_rad);
};
void AreaWeightedSmoothingPreserveShape(const std::unordered_set<networkPoint *> &ps /*copy*/, const double lim_rad = 1E-10) {
   for (const auto &p : ps)
      AreaWeightedSmoothingPreserveShape(p, lim_rad);
};
/* ------------------------------------------------------ */
void LaplacianSmoothingPreserveShape(netPp p, const double lim_rad = 1E-10) {
   try {
      if (!isEdgePoint(p)) { /*端の点はsmoothingしない*/
         auto ps = p->getNeighbors();
         bool isflat = true;
         bool allflat = true;
         Tddd N;
         for (const auto &f : p->getFaces())
            for (const auto &F : p->getFaces())
               if (!isFlat(f->normal, F->normal, 1E-2 * M_PI / 180.)) {
                  allflat = false;
                  break;
               }
         if (ps.size() > 2)  // 2点の場合は2点の中点に動いてしまうので，実行しない
         {
            Tddd V = Mean(extractXtuple(ps)) - ToX(p);
            V = exact_along_surface(p, V);

            for (const auto &f : p->getFaces()) {
               auto [p0, p1, p2] = f->getPoints(p);
               bool isfiniteangles = isFinite(TriangleAngles(T3Tddd{ToX(p0) + V, ToX(p1), ToX(p2)}));
               bool isfiniteareas = isFinite(TriangleArea(T3Tddd{ToX(p0) + V, ToX(p1), ToX(p2)}));
               bool isfinitenormal = isFinite(TriangleNormal(T3Tddd{ToX(p0) + V, ToX(p1), ToX(p2)}));
               // if (!isfiniteangles || !isfiniteareas || !isfinitenormal)
               // {
               // 	isflat = false;
               // 	break;
               // }
               N = TriangleNormal(T3Tddd{ToX(p0) + V, ToX(p1), ToX(p2)});
               if (allflat) {
                  // 全てフラットなので，反転しなければいい．
                  if (!isfiniteangles || !isfiniteareas || !isfinitenormal || !isFlat(N, f->normal, M_PI / 180.)) {
                     isflat = false;
                     break;
                  }
               } else {
                  if (!isfiniteangles || !isfiniteareas || !isfinitenormal || !isFlat(N, f->normal, lim_rad)) {
                     isflat = false;
                     break;
                  }
               }
            }
            if (isflat)
               p->setX(ToX(p) + V);
         }
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
void LaplacianSmoothingPreserveShape(const V_netPp &ps /*copy*/, const double lim_rad = 1E-10) {
   for (const auto &p : ps)
      LaplacianSmoothingPreserveShape(p, lim_rad);
};
void LaplacianSmoothingPreserveShape(const std::unordered_set<networkPoint *> &ps /*copy*/, const double lim_rad = 1E-10) {
   for (const auto &p : ps)
      LaplacianSmoothingPreserveShape(p, lim_rad);
};
/* ------------------------------------------------------ */

void flipIf(Network &water, double limit_angle = M_PI / 180., bool force = false, int times = 0) {
   try {
      // 2022/04/13こっちにBEMのmainから持ってきた
      std::cout << "flipIf" << std::endl;
      water.setGeometricProperties();
      double mean_length = Mean(extLength(water.getLines()));
      bool isfound = false, ismerged = false;
      int count = 0;
      for (const auto &l : water.getLines()) {
         auto [p0, p1] = l->getPoints();
         if (!l->CORNER)
         // if (!p0->CORNER && !p1->CORNER)
         {
            if (force && (times == 0 || count < times)) {
               // p0とp1が角の場合，6という数にこだわる必要がない．
               // つまり6がトポロジカルにベターではない．
               isfound = l->flipIfTopologicalyBetter(limit_angle, M_PI / 180. * 10.);
               if (isfound)
                  count++;
            } else {
               isfound = l->flipIfBetter(limit_angle);
            }
         }
      }
      water.setGeometricProperties();
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
void flipIf(Network &water, const Tdd &limit, bool force = false, int times = 0) {
   try {
      /*
      * フリップによって表面の形状が大きく変わってしまうのはよくない．
      * フリップによって三角形の内角がとても小さくなるのもよくない．
      flipIfでは，このようなフリップをどこまで許容するかを与えることができる．
      * limit_angle：フリップを許容する，辺に隣接する面の法線方向の内角
      * limit_inner_angle：フリップを許容する，フリップによってできる三角形の最小内角の最小
      */
      auto [limit_angle, limit_inner_angle] = limit;
      // 2022/04/13こっちにBEMのmainから持ってきた
      std::cout << "flipIf" << std::endl;
      water.setGeometricProperties();
      double mean_length = Mean(extLength(water.getLines()));
      bool isfound = false, ismerged = false;
      int count = 0;
      for (const auto &l : water.getLines()) {
         auto [p0, p1] = l->getPoints();
         if (!l->CORNER)
         // if (!p0->CORNER && !p1->CORNER)
         {
            if (force && (times == 0 || count < times)) {
               isfound = l->flipIfTopologicalyBetter(limit_angle, limit_inner_angle);
               if (isfound)
                  count++;
            } else {
               isfound = l->flipIfBetter(limit_angle, limit_inner_angle);
               if (isfound)
                  count++;
            }
         }
      }
      water.setGeometricProperties();
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
void flipIf(Network &water, const Tdd &limit_Dirichlet, const Tdd &limit_Neumann, bool force = false, int times = 0) {
   try {
      /*
      * フリップによって表面の形状が大きく変わってしまうのはよくない．
      * フリップによって三角形の内角がとても小さくなるのもよくない．
      flipIfでは，このようなフリップをどこまで許容するかを与えることができる．
      * limit_angle：フリップを許容する，辺に隣接する面の法線方向の内角
      * limit_inner_angle：フリップを許容する，フリップによってできる三角形の最小内角の最小
      */
      auto [limit_angle_D, limit_inner_angle_D] = limit_Dirichlet;
      auto [limit_angle_N, limit_inner_angle_N] = limit_Neumann;
      // 2022/04/13こっちにBEMのmainから持ってきた
      std::cout << "flipIf" << std::endl;
      water.setGeometricProperties();
      double mean_length = Mean(extLength(water.getLines()));
      bool isfound = false, ismerged = false;
      int count = 0;
      for (const auto &l : water.getLines()) {
         auto [p0, p1] = l->getPoints();
         if (!l->CORNER)
         // if (!p0->CORNER && !p1->CORNER)
         {
            if (force && (times == 0 || count < times)) {
               // p0とp1が角の場合，6という数にこだわる必要がない．
               // つまり6がトポロジカルにベターではない．
               if (l->Dirichlet) {
                  isfound = l->flipIfTopologicalyBetter(limit_angle_D, limit_inner_angle_D);
                  if (isfound)
                     count++;
               } else {
                  isfound = l->flipIfTopologicalyBetter(limit_angle_N, limit_inner_angle_N);
                  if (isfound)
                     count++;
               }
            } else {
               if (l->Dirichlet) {
                  isfound = l->flipIfBetter(limit_angle_D, limit_inner_angle_D);
                  if (isfound)
                     count++;
               } else {
                  isfound = l->flipIfBetter(limit_angle_N, limit_inner_angle_N);
                  if (isfound)
                     count++;
               }
            }
         }
      }
      water.setGeometricProperties();
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

/* ------------------------------------------------------ */
// void LaplacianSmoothingIfFlat(netPp p) {
//    try {
//       if (!isEdgePoint(p) && isFlat(p) && isInConvexPolygon(p)) { /*端の点はsmoothingしない*/
//          auto ps = p->getNeighbors();
//          if (ps.size() > 2)  // 2点の場合は2点の中点に動いてしまうので，実行しない
//             p->setX(Mean(extractX(ps)));
//       }
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorOff << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };
//////////////////////////////////////////////////////////
// void LaplacianSmoothingIfFlat(V_netPp ps /*copy*/) {
//    for (const auto &p : ps)
//       LaplacianSmoothingIfFlat(p);
// };
///////////////////////////////////////////////////////////
// void LaplacianSmoothingIfFlat_ExceptIntX(V_netPp ps /*copy*/) {
//    try {
//       std::shuffle(std::begin(ps), std::end(ps), std::default_random_engine());
//       for (const auto &p : ps)
//          if (!isEdgePoint(p) && isFlat(p) && isInConvexPolygon(p)) {  // 端の点はsmoothingしない
//             auto countintx = 0;
//             for (const auto &l : p->getLines())
//                if (l->isIntxn())
//                   countintx++;
//             if (!(countintx > 1)) {
//                auto tmp = Mean(extractX(p->getNeighbors()));
//                p->setX(tmp);
//             }
//          }
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorOff << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };
///////////////////////////////////////////////////////////
void LaplacianSmoothingIfOnStraightLine(V_netPp ps) {
   try {
      std::shuffle(std::begin(ps), std::end(ps), std::default_random_engine());
      for (auto p : ps)
         if (!isEdgePoint(p) && !isFlat(p))  // 端の点はsmoothingしない
         {
            std::map<netPp, V_d> mapP_vec;
            for (auto q : p->getNeighbors())
               mapP_vec[q] = ToVector(q->getXtuple() - ToX(p)) / Norm(q->getXtuple() - ToX(p));

            netPp p0, p1;
            bool isStraight = false;
            auto Ps = TakeFirst(mapP_vec);
            double minarea = 1E-2;
            for (auto i = 0; i < Ps.size() - 1; i++)
               for (auto j = i + 1; j < Ps.size(); j++) {
                  // (Ps[i]) <--------- (p) <-------- (Ps[j]) : angle = 0
                  double abs_angle = std::abs(VectorAngle(Ps[i]->getXtuple() - ToX(p), ToX(p) - Ps[j]->getXtuple()));
                  if (abs_angle < minarea) {
                     p0 = Ps[i];
                     p1 = Ps[j];
                     isStraight = true;
                     minarea = abs_angle;
                  }
               }
            if (isStraight)
               p->setX(ToVector(ToX(p0) + ToX(p1)) / 2.);
         }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
double Distance(const Tddd &A, const Tddd &B) { return Norm(A - B); };

double Distance(const Tddd &X, const networkPoint *const p) {
   return Norm(ToX(p) - X);
};
double Distance(const networkPoint *const p, const Tddd &X) {
   return Norm(ToX(p) - X);
};
double Distance(const networkPoint *const p, const networkPoint *const q) {
   return Norm(p->X - q->X);
};
double Distance(const networkLine *const l, const networkPoint *const p) {
   return Norm(ToX(p) - l->getXtuple());
};
double Distance(const networkPoint *const p, const networkLine *const l) {
   return Norm(ToX(p) - l->getXtuple());
};
V_d Distance(const networkPoint *const p, const std::vector<networkPoint *> &ps) {
   V_d ret(ps.size());
   int i = 0;
   for (const auto &q : ps)
      ret[i++] = Distance(q, p);
   return ret;
};
V_d Distance(const std::vector<networkPoint *> &ps, const networkPoint *p) {
   return Distance(p, ps);
};
/*
sortは，
comp(a,b)->trueなら
comp(b,a)->false
でないといけない
*/
void sortByLength(V_netLp &lines) {
   std::sort(lines.begin(), lines.end(),
             [](const auto &v, const auto &w) {
                if (v->length() - w->length() <= 1E-10)
                   return v < w;
                else
                   return v->length() < w->length();
             });
};

void reverseSortByLength(V_netLp &lines) {
   std::sort(lines.begin(), lines.end(),
             [](const auto &w, const auto &v) {
                if (v->length() - w->length() <= 1E-10)
                   return v < w;
                else
                   return v->length() < w->length();
             });
};

// void reverseSortByLength(const std::unordered_set<networkLine *> &uo_lines) {
//    V_netLp lines(uo_lines.begin(), uo_lines.end());
//    std::sort(lines.begin(), lines.end(),
//              [](const auto &w, const auto &v) {
//                 if (v->length() - w->length() <= 1E-10)
//                    return v < w;
//                 else
//                    return v->length() < w->length();
//              });
// };

void sortByDistance(V_netPp &points, const netPp origin) {
   std::sort(points.begin(), points.end(),
             [&origin](const netPp &v, const netPp &w) {
                if (Norm(ToX(v) - ToX(w)) <= 1E-13)
                   return v < w;
                else
                   return Norm(ToX(v) - ToX(origin)) < Norm(ToX(w) - ToX(origin));
             });
};

void sortByDistance(V_netPp &points, const Tddd &origin) {
   std::sort(points.begin(), points.end(),
             [&](const netPp &v, const netPp &w) {
                if (Norm(ToX(v) - ToX(w)) <= 1E-13)
                   return v < w;
                else
                   return Norm(ToX(v) - origin) < Norm(ToX(w) - origin);
             });
};

// void Merge(const Network *net, const std::function<bool(const networkLine *)> &condition) {
//    bool found = true;
//    int count = 0;
//    while (found && count++ < 5) {
//       found = false;
//       for (const auto &l : net->Lines)
//          if (!found && condition(l) && l->merge()) {
//             found = true;
//             break;
//          }
//    }
// };
void Merge(const std::unordered_set<networkLine *> &uo_lines, const std::function<bool(const networkLine *)> &func) {
   for (const auto &l : uo_lines)
      if (func(l)) {
         l->merge();
         return;
      }
   // std::unordered_set<std::shared_ptr<networkLine>> shared_pointers;
   // std::transform(uo_lines.begin(), uo_lines.end(), std::inserter(shared_pointers, shared_pointers.begin()),
   //                [](networkLine *raw_ptr) { return std::make_shared<networkLine>(*raw_ptr); });
   // std::vector<std::shared_ptr<networkLine>> shared_pointers_copy(shared_pointers.begin(), shared_pointers.end());
   // for (const auto &l : shared_pointers_copy) {
   //    if (func(l.get())) {
   //       l->merge();
   //       shared_pointers.erase(l);
   //    }
   // }
};

void Divide(const std::unordered_set<networkLine *> &uo_lines, const std::function<bool(const networkLine *)> &func) {
   for (const auto &l : uo_lines)
      if (func(l)) {
         l->divide();
         // std::cout << "divide" << std::endl;
      }
};

void Divide(const std::unordered_set<networkLine *> &uo_lines, const double lim_len) {
   std::cout << Red << "|";
   for (const auto &l : uo_lines)
      if (l->length() >= lim_len) {
         std::cout << Green << "|";
         l->divide();
      }
   std::cout << Blue << "|" << colorOff;
   //
   // std::vector<networkLine *> lines(uo_lines.begin(), uo_lines.end());
   // std::cout << Red << "|";
   // reverseSortByLength(lines);
   // for (const auto &l : lines)
   //    if (l->length() >= lim_len) {
   //       std::cout << Green << "|";
   //       l->divide();
   //    }
   // // else
   // //    break;
   // std::cout << Blue << "|" << colorOff;
};

// std::vector<V_netPp> triangulate(const V_netPp &objects, const V_d &normal, const double smallangle = 0.) {
//    // std::cout << __PRETTY_FUNCTION__ << std::endl;
//    try {
//       if (!DuplicateFreeQ(objects)) {
//          std::cout << "objects = " << objects << std::endl;
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Duplication is found !");
//       } else if (objects.size() < 3) {
//          Print(objects);
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "objects.size()<3");
//       } else {
//          std::vector<V_netPp> ret({});
//          geometry::polygon poly(extractX(objects));
//          for (const auto &ind : poly.triangulate(normal, smallangle)) {
//             V_netPp points = {objects[ind[0]], objects[ind[1]], objects[ind[2]]};
//             TriangleNormal(points[0]->getXtuple(), points[1]->getXtuple(), points[2]->getXtuple());
//             ret.emplace_back(points);
//          }
//          return ret;
//       }
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorOff << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };
//////////////////////////////////////////////////////////////
//////////////////////////// 外部関数 /////////////////////////
//////////////////////////////////////////////////////////////
std::vector<netPp> getPointsWithoutLines(const std::vector<netPp> &obj) {
   std::vector<netPp> ret;
   for (const auto &o : obj)
      if (o->getLines().empty())
         ret.emplace_back(o);
   return ret;
};
///////
std::vector<netPp> getPointsWithoutFaces(const std::vector<netPp> &obj) {
   std::vector<netPp> ret;
   V_netLp ls;
   bool facefound;
   for (const auto &o : obj) {
      ls = o->getLines();
      if (ls.empty())
         ret.emplace_back(o);
      else {
         facefound = false;
         for (const auto &l : ls)
            if (!l->getFaces().empty()) {
               facefound = true;
               break;
            }
         if (!facefound)
            ret.emplace_back(o);
      }
   }
   return ret;
};

V_d areas(const V_netFp &fs) {
   V_d ret(fs.size());
   for (auto i = 0; i < fs.size(); i++)
      ret[i] = fs[i]->area;
   return ret;
};
/* ------------------------------------------------------ */
std::vector<std::string> getNames(const V_netFp &fs) {
   std::vector<std::string> ret;
   for (const auto &f : fs)
      ret.emplace_back(f->getNetwork()->getName());
   return ret;
};

std::vector<std::string> getNames(const V_netPp &ps) {
   std::vector<std::string> ret;
   for (const auto &p : ps)
      ret.emplace_back(p->getNetwork()->getName());
   return ret;
};

std::vector<std::string> getNames(const V_Netp &nets) {
   std::vector<std::string> ret;
   for (const auto &n : nets)
      ret.emplace_back(n->getName());
   return ret;
};

V_netLp takeIntxn(const V_netLp &ls) {
   V_netLp ret({});
   for (const auto &l : ls)
      if (l->isIntxn())
         ret.emplace_back(l);
   return ret;
};

template <typename T>
std::vector<T *> takeNetwork(const std::vector<T *> &ps, const std::vector<Network *> net) {
   std::vector<T *> ret = {};
   for (const auto &p : ps)
      if (MemberQ(net, p->getNetwork()))
         ret.emplace_back(p);
   return ret;
};

////////////////////////////////////////////////////////////////
////////////////////////// dislay //////////////////////////////
////////////////////////////////////////////////////////////////

template <typename T>
void displayNames(const std::vector<T *> &ps) {
   std::cout << __PRETTY_FUNCTION__ << std::endl;
   std::map<Network *, int> storages;
   for (const auto &p : ps)
      if (storages.find(p->getNetwork()) != storages.end())
         storages[p->getNetwork()] += 1;
      else
         storages[p->getNetwork()] = 1;
   for (auto [n, i] : storages)
      std::cout << std::setw(8) << n->getName() << "  " << std::setw(8) << i << std::endl;
   if (storages.empty())
      std::cout << Red << "not stored !?" << colorOff << std::endl;
};

template <typename T>
void displayStorages(const std::vector<T *> &ps) {
   std::cout << __PRETTY_FUNCTION__ << std::endl;
   std::map<Network *, int> storages;
   for (const auto &p : ps)
      if (storages.find(p->getStorage()) != storages.end())
         storages[p->getStorage()] += 1;
      else
         storages[p->getStorage()] = 1;
   for (auto [n, i] : storages)
      std::cout << std::setw(8) << n->getName() << "  " << std::setw(8) << i << std::endl;
   if (storages.empty())
      std::cout << Red << "not stored !?" << colorOff << std::endl;
};

void display(const std::unordered_set<networkLine *> &ls) {
   std::cout << __PRETTY_FUNCTION__ << std::endl;
   std::map<size_t, int> data;
   int intx = 0;
   for (const auto &l : ls) {
      auto s = l->getFaces().size();
      if (l->isIntxn())
         intx++;
      if (data.find(s) != data.end())
         data[s] += 1;
      else
         data[s] = 1;
   }
   for (auto [n, i] : data)
      std::cout << std::setw(8) << "Faces size = " << n << "  " << std::setw(8) << i << std::endl;

   std::cout << std::setw(8) << "intx = " << intx << "  " << std::setw(8) << intx << std::endl;
};

void display(Network *net) {
   std::cout << __PRETTY_FUNCTION__ << std::endl;
   V_i num(15, 0), pointsLines(15, 0), linesFaces(15, 0);
   ;
   for (const auto &p : net->getPoints())
      for (auto i = 0; i < 15; i++)
         if (p->getLines().size() == (size_t)i)
            pointsLines[i]++;
   // for (const auto &f : net->getFaces())
   // {
   // 	auto f_Points = f->getPoints();
   // 	for (auto i = 0; i < 15; i++)
   // 	{
   // 		if (f_Points.size() == (size_t)i)
   // 			facesPoints[i]++;
   // 		// if (f->getLines().size() == (size_t)i)
   // 		// 	facesLines[i]++;
   // 	}
   // }
   for (const auto &l : net->getLines()) {
      auto fs = l->getFaces();
      for (auto i = 0; i < 15; i++)
         if (fs.size() == (size_t)i)
            linesFaces[i]++;
   }
   std::cout << "Network name : " << net->getName() << std::endl;
   std::cout << "--------------------Names---------------------------" << std::endl;
   std::cout << Blue << "Points : " << colorOff;
   // displayNames(net->getPoints());
   std::cout << Magenta << " Faces : " << colorOff;
   // displayNames(net->getFaces());
   std::cout << "------------------Storages--------------------------" << std::endl;
   std::cout << Blue << "Points : " << colorOff;
   // displayStorages(net->getPoints());
   std::cout << Magenta << " Faces : " << colorOff;
   // displayStorages(net->getFaces());
   std::cout << "--------------------Lines--------------------------" << std::endl;
   display(net->getLines());
   std::cout << "-----------------------------------------------" << std::endl;
   std::cout << Magenta << "Points.size() : " << net->getPoints().size() << std::endl;
   std::cout << Magenta << " Faces.size() : " << net->getFaces().size() << std::endl;
   std::cout << "-----------------------------------------------" << std::endl;
   int i = 0;
   for (auto &n : num) {
      n = i++;
   };
   std::cout << Red << "                  ";
   for (const auto &n : num) {
      std::cout << std::setw(6) << n;
   };
   std::cout << colorOff << std::endl;

   std::cout << magenta << "Lines of points : ";
   for (const auto &n : pointsLines) {
      std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
   }
   std::cout << colorOff << std::endl;

   // std::cout << magenta << "Points of faces : ";
   // for (const auto &n : facesPoints)
   // {
   // 	std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
   // };
   // std::cout << colorOff << std::endl;

   // std::cout << magenta << " Lines of faces : ";
   // for (const auto &n : facesLines)
   // {
   // 	std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
   // };
   std::cout << colorOff << std::endl;

   std::cout << magenta << " Faces of Lines : ";
   for (const auto &n : linesFaces) {
      std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
   };
   std::cout << colorOff << std::endl;

   V_i connection(3);
   net->setLinesStatus(true);
   for (const auto &p : net->getPoints())
      for (const auto &line : p->getLines_toggle(true)) {
         if ((line->getFaces()).size() == 0)
            connection[0]++;
         else if ((line->getFaces()).size() == 1)
            connection[1]++;
         else if ((line->getFaces()).size() == 2)
            connection[2]++;
      }
   std::cout << magenta << " Connection of faces : " << connection;
   if (connection[0] == 0 && connection[1] == 0 && connection[2] != 0)
      std::cout << Blue << " <- Network is a closed surface";
   std::cout << colorOff << std::endl;
   std::cout << "-----------------------------------------------" << std::endl;
};

//* ------------------------------------------------------ */
//*                     粒子法のためのユーティリティ            */
//* ------------------------------------------------------ */
std::vector<Tddd> fullparticlize(const Tddd &X0, const Tddd &X1, const Tddd &X2, const double dx) {
   /*
   |-o--o--o--o-| n = 4
   x,yはパラメタ
   */
   interpolationTriangleLinearTuple intp({X0, X1, X2});
   auto y_list = [&intp, &dx](const double x) {
      if (1. - x < 1E-10)
         return V_d{0.};
      int n = std::round(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
      if (n > 0)
         return Subdivide(0., 1. - x, n);
      else
         return V_d{};
   };
   /* ------------------------------------------------------ */
   Tddd v = intp(0., 1.) - intp(0., 0.);
   Tddd u = intp(1., 0.) - intp(0., 0.);
   double height = Norm(u - Dot(Normalize(v), u) * Normalize(v));
   int n = std::round(height / dx);
   /* ------------------------------------------------------ */
   std::vector<Tddd> ret = {};
   if (n > 0) {
      for (const auto &x : Subdivide(0., 1., n))
         for (const auto &y : y_list(x))
            ret.emplace_back(intp(x, y));
   }
   return ret;
};
std::vector<Tddd> fullparticlize(networkFace const *f, const double dx) {
   // auto ps = f->getPoints();
   // return fullparticlize(ps[0]->getXtuple(), ps[1]->getXtuple(), ps[2]->getXtuple(), dx);

   auto [p0, p1, p2] = f->getPoints();
   return fullparticlize(ToX(p0), ToX(p1), ToX(p2), dx);
};

/*
        % particlizeは座標とパラメタを返す
        % 面をdx間隔で粒子化し，各点の情報をベクトルで返す
*/

// b% -------------- networkFaceを使ったparticlize ------------- */

std::vector<std::tuple<Tddd, Tdd>> particlize(const networkFace *const f, const double dx) { return particlize(f->getXVertices(), dx); };
std::vector<std::tuple<Tddd, Tdd>> particlize(const V_netFp &faces, const double dx) {
   std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
   for (const auto &f : faces) {
      tmp = particlize(f, dx);
      ret.insert(ret.end(), tmp.begin(), tmp.end());
   }
   return ret;
};
std::vector<std::tuple<Tddd, Tdd>> particlize(const std::unordered_set<networkFace *> &faces, const double dx) {
   std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
   for (const auto &f : faces) {
      tmp = particlize(f, dx);
      ret.insert(ret.end(), tmp.begin(), tmp.end());
   }
   return ret;
};

/* ------------------------ 深さあり ------------------------ */

std::vector<std::tuple<Tddd, Tdd>> particlize(const networkFace *const f, const double dx, const V_d &depth_list /*法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など*/) {
   std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
   double alpha;
   Tddd X0, X1, X2;
   for (const auto &d : depth_list /*double実際の長さ*/) {
      auto [p0, p1, p2] = f->getPoints();
      X0 = ToX(p0) + d / Dot(p0->getNormalTuple(), f->normal) * p0->getNormalTuple();
      X1 = ToX(p1) + d / Dot(p1->getNormalTuple(), f->normal) * p1->getNormalTuple();
      X2 = ToX(p2) + d / Dot(p2->getNormalTuple(), f->normal) * p2->getNormalTuple();
      tmp = particlize({X0, X1, X2}, dx);
      ret.insert(ret.end(), tmp.begin(), tmp.end());
   }
   return ret;
};

Tddd particlize(const networkFace *f, const Tdd &t0t1, const double d) {
   auto [p0, p1, p2] = f->getPoints();
   return Dot(Tddd{std::get<0>(t0t1), std::get<1>(t0t1), 1. - std::get<0>(t0t1) - std::get<1>(t0t1)},
              T3Tddd{ToX(p0) + d / Dot(p0->getNormalTuple(), f->normal) * p0->getNormalTuple(),
                     ToX(p1) + d / Dot(p1->getNormalTuple(), f->normal) * p1->getNormalTuple(),
                     ToX(p2) + d / Dot(p2->getNormalTuple(), f->normal) * p2->getNormalTuple()});
};

// b% ------------------------------------------------------ */
Tddd oppositeX(const std::tuple<networkFace * /*補間に使った三角形の頂点*/,
                                T_PPP /*補間に使った三角形の頂点*/,
                                Tdd /*パラメタt0,t1*/,
                                double /*深さ方向距離*/,
                                double /*粒子間隔*/> &particlize_info) {
   // 粒子化された点の面にたいする反対側の位置を返す．
   // 単純にp+2*Dot(f->p0 - p,n)*nでは，角の面において重なりが生じてしまう．
   // パラメタを使って計算すれば重なりは生じない．
   auto [f, p0p1p2, t0t1, d, dx] = particlize_info;
   if (f) {
      auto [p0, p1, p2] = p0p1p2;  // f->getPoints();
      return Dot(Tddd{std::get<0>(t0t1), std::get<1>(t0t1), 1. - std::get<0>(t0t1) - std::get<1>(t0t1)},
                 T3Tddd{ToX(p0) + d / Dot(p0->getNormalTuple(), f->normal) * p0->getNormalTuple(),
                        ToX(p1) + d / Dot(p1->getNormalTuple(), f->normal) * p1->getNormalTuple(),
                        ToX(p2) + d / Dot(p2->getNormalTuple(), f->normal) * p2->getNormalTuple()});
   } else {
      std::stringstream ss;
      // ss << "particlize_info = " << f << ", " << t0 << ", " << t1 << ", " << d << ", " << dx << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   };
};

Tddd oppositeX(const networkPoint *p) {
   return oppositeX(p->particlize_info);
};
/* ------------------------------------------------------ */

std::vector<networkPoint *> operator+(std::vector<networkPoint *> A /*copy and return*/, const std::vector<networkPoint *> &B) {
   A.insert(A.end(), B.begin(), B.end());
   return A;
};
std::vector<networkPoint *> &operator+=(std::vector<networkPoint *> &v, const std::vector<networkPoint *> &w) {
   return (v = v + w);
};
#endif
