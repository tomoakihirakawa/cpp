#ifndef NetworkUtility_H
#define NetworkUtility_H

#include "Network.hpp"

/* -------------------------------------------------------------------------- */

void writeVertices(std::ofstream &ofs, std::map<networkPoint *, int> &P_i, const auto &points) {
   int i = 0;
   for (const auto &p : points) {
      P_i[p] = ++i;
      auto [X0, X1, X2] = ToX(p);
      ofs << "v " << X0 << " " << X1 << " " << X2 << std::endl;
   }
   ofs << std::endl;
}

void writeFaces(std::ofstream &ofs, const std::map<networkPoint *, int> &P_i, const auto &faces) {
   for (const auto &f : faces) {
      auto [p0, p1, p2] = f->getPoints();
      ofs << "f " << P_i.at(p0) << " " << P_i.at(p1) << " " << P_i.at(p2) << std::endl;
   }
}

void createOBJ(std::ofstream &ofs, Network &net) {
   std::map<networkPoint *, int> pointIndices;
   writeVertices(ofs, pointIndices, net.getPoints());
   writeFaces(ofs, pointIndices, net.getFaces());
}

//% ------------------------------------------------------ */
//%                         反射点の計算                     */
//% ------------------------------------------------------ */

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

bool isStraight(const netPp p0, const netPp p1, const netPp p2, const double angle = 1E-2 /*閾値*/) {
   return (isStraight(ToVector(ToX(p0) - ToX(p1)), ToVector(ToX(p1) - ToX(p2))) < angle);
};

bool isFlat(const netPp p, double minangle = M_PI / 180.) {
   if (p->getLines().empty())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty()");
   auto faces = ToVector(p->getFaces());
   for (auto i = 0; i < faces.size(); i++)
      for (auto j = i; j < faces.size(); j++)
         if (!isFlat(ToX(faces[i]), ToX(faces[j]), minangle))
            return false;
   return true;
};

bool isFlat(const netLp line, double minangle = M_PI / 180.) {
   auto fs = line->getFaces();
   if (fs.size() != 2)
      return false;
   else
      return isFlat(ToX(fs[0]), ToX(fs[1]), minangle);
};

/* -------------------------------------------------------------------------- */

Tddd exact_along_surface(const networkPoint *const p, Tddd VECTOR) {
   for (const auto &f : p->getFaces())
      VECTOR = Chop(VECTOR, f->normal);
   return VECTOR;
   // auto faces = ToVector(p->getFaces());
   // std::vector<std::tuple<Tddd, Tddd>> normals;
   // normals.reserve(2 * faces.size());
   // for (auto i = 0; i < faces.size(); i++)
   //    for (auto j = i; j < faces.size(); j++)
   //       if (!isFlat(faces[i]->normal, faces[j]->normal, 1E-2 * M_PI / 180.))
   //          normals.push_back({faces[i]->normal, faces[j]->normal});
   // for (const auto &[N0, N1] : normals)
   //    VECTOR = Projection(VECTOR, Cross(N0, N1));
   // return VECTOR;
};

bool isFlat(const auto p, const double lim_rad = 1E-2 * M_PI / 180.) {
   auto faces = p->getFaces();
   for (auto i = 0; i < faces.size(); i++)
      for (auto j = i; j < faces.size(); j++)
         if (!isFlat(faces[i]->normal, faces[j]->normal, lim_rad))
            return false;
   return true;
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT smoothing_vector

### 格子の平滑化

境界面を格子を綺麗な状態に整えるために，境界面の点を移動させる．
移動方法は，様々考えられる．

| 移動方法 | 説明 |
|:-------|:---------|
| Laplacian smoothing                     | 境界面の点をその点の近傍の点の重心に移動させる．          |
| Area weighted smoothing                 | 隣接面の面積の重みを掛けて，面の中心に移動させる．        |
| Distorsion measure weighted smoothing   | 隣接面の歪みに関する係数を重みとして掛けて，面の中心に移動させる．|

共通点は，許されない移動を防止する移動前に，移動後の形状をチェックすることである．
`canFlip`と同様，\ref{isValidTriangle}{`isValidTriangle`}を用いる．

徐々に移動させる場合，誤差の蓄積，条件の変化を把握するのが難しい．
大きな変化は防げても，小さな変化には対応できない場合が考えられる．

*/

void SmoothingPreserveShape(netPp p, const std::function<Tddd(const netPp)> &SmoothingVector) {
   try {
      if (!isEdgePoint(p)) { /*端の点はsmoothingしない*/
         auto ps = p->getNeighbors();
         // 無視できる角度
         const double negligible_angle = 1E-2 * M_PI / 180.;
         const bool negligibly_flat = isFlat(p, negligible_angle);
         const double acceptable_change_angle = 1E-5 * M_PI / 180.;
         const auto faces = ToVector(p->getFaces());
         Tddd X0, X1, X2;
         auto check_if_next_faces_valid = [&](Tddd &V) {
            for (const auto &f : faces) {
               auto [p0, p1, p2] = f->getPoints();
               if (p0 != p && p1 != p && p2 != p)
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p0 != p && p1 != p && p2 != p");
               X0 = p0->X;
               X1 = p1->X;
               X2 = p2->X;
               if (p0 == p)
                  X0 += V;
               else if (p1 == p)
                  X1 += V;
               else if (p2 == p)
                  X2 += V;

               T3Tddd vertices = {X0, X1, X2};
               if (!isValidTriangle(vertices, 1E-2 * M_PI / 180.))
                  return false;
               if (!isFlat(Cross(p1->X - p0->X, p2->X - p0->X), Cross(X1 - X0, X2 - X0), negligibly_flat ? negligible_angle / 2. : acceptable_change_angle))
                  return false;
               if (!isFlat(Cross(p2->X - p1->X, p0->X - p1->X), Cross(X2 - X1, X0 - X1), negligibly_flat ? negligible_angle / 2. : acceptable_change_angle))
                  return false;
               if (!isFlat(Cross(Normalize(p1->X - p0->X), Normalize(p2->X - p0->X)), Cross(Normalize(X1 - X0), Normalize(X2 - X0)), negligibly_flat ? negligible_angle / 2. : acceptable_change_angle))
                  return false;
               if (!isFlat(Cross(Normalize(p2->X - p1->X), Normalize(p0->X - p1->X)), Cross(Normalize(X2 - X1), Normalize(X0 - X1)), negligibly_flat ? negligible_angle / 2. : acceptable_change_angle))
                  return false;
            }
            return true;
         };

         if (ps.size() > 2)  // 2点の場合は2点の中点に動いてしまうので，実行しない
         {
            Tddd V = exact_along_surface(p, SmoothingVector(p));
            for (auto i = 0; i < 10; ++i) {
               if (check_if_next_faces_valid(V)) {
                  p->setX(p->X + V);
                  return;
               } else
                  V /= 2.;
            }
         }
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

/* ------------------------------------------------------ */

Tddd DistorsionMeasureWeightedSmoothingVector2(const networkPoint *p, std::function<Tddd(const networkPoint *)> position) {
   if (!isEdgePoint(p)) {
      Tddd ret = {0., 0., 0.}, V = {0., 0, 0.}, X, Xmid, vertical;
      T3Tddd t3tdd;
      double Wtot = 0, W, height;
      for (const auto &f : p->getFaces()) {
         auto [p0, p1, p2] = f->getPoints(p);
         auto [X0, X1, X2] = t3tdd = {position(p0), position(p1), position(p2)};
         W = std::log2(CircumradiusToInradius(t3tdd) - 1.) + 1;
         if (p0->CORNER)
            W *= 1.5;
         if (p1->CORNER)
            W *= 1.5;
         if (p2->CORNER)
            W *= 1.5;
         W *= W;
         Wtot += W;
         Xmid = (X2 + X1) / 2.;
         vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
         height = Norm(X2 - X1) * sin(M_PI / 3.);
         X = height * vertical + Xmid;
         V += W * X;
      }
      double Wtot_first = Wtot;
      for (const auto &l : p->getLines())
         if (l->CORNER) {
            W = Wtot_first / 2.;
            Wtot += W;
            V += W * position((*l)(p));
         }
      return V / Wtot - position(p);
   } else
      return {0., 0., 0.};
};

Tddd DistorsionMeasureWeightedSmoothingVector(const networkPoint *p, std::function<Tddd(const networkPoint *)> position) {
   if (!isEdgePoint(p)) {
      Tddd ret = {0., 0., 0.}, V = {0., 0, 0.}, X, Xmid, vertical;
      T3Tddd t3tdd;
      double Wtot = 0, W, height;
      for (const auto &f : p->getFaces()) {
         auto [p0, p1, p2] = f->getPoints(p);
         auto [X0, X1, X2] = t3tdd = {position(p0), position(p1), position(p2)};
         W = std::log2(CircumradiusToInradius(t3tdd) - 1.);
         Wtot += W;
         Xmid = (X2 + X1) / 2.;
         vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
         height = Norm(X2 - X1) * sin(M_PI / 3.);
         X = height * vertical + Xmid;
         V += W * X;
      }
      return V / Wtot - position(p);
   } else
      return {0., 0., 0.};
};

Tddd DistorsionMeasureWeightedSmoothingVector(const networkPoint *p) {
   if (!isEdgePoint(p)) {
      Tddd ret = {0., 0., 0.}, V = {0., 0, 0.}, X, Xmid, vertical;
      T3Tddd t3tdd;
      double Wtot = 0, W, height;
      for (const auto &f : p->getFaces()) {
         t3tdd = ToX(f->getPoints(p));
         auto [X0, X1, X2] = t3tdd;
         W = std::log2(CircumradiusToInradius(t3tdd) - 1.);
         Wtot += W;
         Xmid = (X2 + X1) / 2.;
         vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
         height = Norm(X2 - X1) * sin(M_PI / 3.);
         X = height * vertical + Xmid;
         V += W * X;
      }
      return V / Wtot - p->X;
   } else
      return {0., 0., 0.};
};

void DistorsionMeasureWeightedSmoothingPreserveShape(const auto &ps, const int iteration = 1) {
   // gradually approching to given a
   const double aIN = 0.05;
   for (auto i = 0; i < iteration; ++i) {
      double a = aIN * ((i + 1) / (double)(iteration));
      for (const auto &p : ps)
         SmoothingPreserveShape(p, [&](const networkPoint *q) -> Tddd { return a * DistorsionMeasureWeightedSmoothingVector(q); });
   }
};

/* ------------------------------------------------------ */

Tddd AreaWeightedSmoothingVector(const networkPoint *p, std::function<Tddd(const networkPoint *)> position) {
   if (!isEdgePoint(p)) {
      Tddd X, V = {0., 0, 0.};
      T3Tddd t3tdd;
      double Wtot = 0, W;
      for (const auto &f : p->getFaces()) {
         auto [p0, p1, p2] = f->getPoints();
         t3tdd = {position(p0), position(p1), position(p2)};
         Wtot += (W = TriangleArea(t3tdd));
         V += W * (X = Centroid(t3tdd));
      }
      return V / Wtot - position(p);
   } else
      return {0., 0., 0.};
}

Tddd AreaWeightedSmoothingVector(const networkPoint *p) {
   if (!isEdgePoint(p)) {
      Tddd X, V = {0., 0, 0.};
      T3Tddd t3tdd;
      double Wtot = 0, W;
      for (const auto &f : p->getFaces()) {
         t3tdd = ToX(f->getPoints());
         Wtot += (W = TriangleArea(t3tdd));
         V += W * (X = Centroid(t3tdd));
      }
      return V / Wtot - p->X;
   } else
      return {0., 0., 0.};
};

void AreaWeightedSmoothingPreserveShape(const auto &ps, const int iteration = 1) {
   const double aIN = 0.1;
   for (auto i = 0; i < iteration; ++i) {
      double a = aIN * ((i + 1) / (double)(iteration));
      for (const auto &p : ps)
         SmoothingPreserveShape(p, [&](const networkPoint *q) -> Tddd { return a * AreaWeightedSmoothingVector(q); });
   }
};
/* ------------------------------------------------------ */

Tddd ArithmeticWeightedSmoothingVector(const networkPoint *p, std::function<Tddd(const networkPoint *)> position) {
   if (!isEdgePoint(p)) {
      Tddd V = {0.0, 0.0, 0.0};
      for (const auto &q : p->getNeighbors())
         V += position(q);
      return V / static_cast<double>(p->getNeighbors().size()) - position(p);
   } else {
      return {0.0, 0.0, 0.0};
   }
}

Tddd ArithmeticWeightedSmoothingVector(const networkPoint *p) {
   if (!isEdgePoint(p)) {
      Tddd V = {0., 0, 0.};
      for (const auto &q : p->getNeighbors()) V += q->X;
      return V / p->getNeighbors().size() - p->X;
   } else
      return {0., 0., 0.};
};

void LaplacianSmoothingPreserveShape(const auto &ps, const int iteration = 1) {
   const double aIN = 0.1;
   for (auto i = 0; i < iteration; ++i) {
      double a = aIN * ((i + 1) / (double)(iteration));
      for (const auto &p : ps)
         SmoothingPreserveShape(p, [&](const networkPoint *q) -> Tddd { return a * ArithmeticWeightedSmoothingVector(q); });
   }
};

/* ------------------------------------------------------ */

void flipIf(Network &water, const Tdd &limit_Dirichlet, const Tdd &limit_Neumann, bool force = false, int iteration = 0) {
   try {
      auto [target_of_max_normal_diffD, acceptable_normal_change_by_flipD] = limit_Dirichlet;
      auto [target_of_max_normal_diffN, acceptable_normal_change_by_flipN] = limit_Neumann;
      // 2022/04/13こっちにBEMのmainから持ってきた
      std::cout << "flipIf" << std::endl;
      water.setGeometricProperties();
      double mean_length = Mean(extLength(water.getLines()));
      bool isfound = false, ismerged = false;
      int count = 0;
      auto V = ToVector(water.getLines());
      for (const auto &l : RandomSample(V)) {
         auto [p0, p1] = l->getPoints();
         if (!l->CORNER) {
            if (force && (iteration == 0 || count < iteration)) {
               // p0とp1が角の場合，6という数にこだわる必要がない．
               // つまり6がトポロジカルにベターではない．
               if (l->Dirichlet) {
                  isfound = l->flipIfTopologicallyBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD);
                  if (isfound) count++;
               } else {
                  isfound = l->flipIfTopologicallyBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN);
                  if (isfound) count++;
               }
            } else {
               if (l->Dirichlet) {
                  isfound = l->flipIfBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD, 5);
                  if (isfound) count++;
               } else {
                  isfound = l->flipIfBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN, 5);
                  if (isfound) count++;
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

void flipIf(Network &water, const Tdd &limit, bool force = false, int iteration = 0) {
   flipIf(water, limit, limit, force, iteration);
};

void flipIf(Network &water, double limit_angle = M_PI / 180., bool force = false, int iteration = 0) {
   flipIf(water, limit_angle, force, iteration);
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
double Distance(const Tddd &A, const Tddd &B) {
   return Norm(A - B);
};

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

void Divide(std::unordered_set<networkLine *> uo_lines, const std::function<bool(const networkLine *)> &func) {
   for (const auto &l : uo_lines)
      if (func(l))
         l->divide();
};

void Divide(const std::unordered_set<networkLine *> &uo_lines, const double lim_len) {
   std::cout << Red << "|";
   for (const auto &l : uo_lines)
      if (l->length() >= lim_len) {
         std::cout << Green << "|";
         l->divide();
      }
   std::cout << Blue << "|" << colorOff;
};

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

std::vector<std::tuple<Tddd, Tdd>> particlize(const networkFace *const f, const double dx) {
   return particlize(f->getXVertices(), dx);
};
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
