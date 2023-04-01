#ifndef networkPoint_H
#define networkPoint_H
#pragma once

#include "InterpolationRBF.hpp"
#include "Network.hpp"

// void networkPoint::calculateBufferPotentialsOnClungSurface()
// {
// 	auto [X, t0, t1, line, face] = this->clungSurface;
// 	if (line)
// 	{
// 		auto p1 = (*line)(this);
// 		Tdd N = {t0, 1. - t0};
// 		this->phiphin_BUFFER = Dot(N, T2Tdd{this->phiphin, p1->phiphin});
// 		this->phiphin_t_BUFFER = Dot(N, T2Tdd{this->phiphin_t, p1->phiphin_t});
// 		this->phi_Neumann_BUFFER = Dot(N, Tdd{this->phi_Neumann, p1->phi_Neumann});
// 		this->phi_Dirichlet_BUFFER = Dot(N, Tdd{this->phi_Dirichlet, p1->phi_Dirichlet});
// 		this->phin_Neumann_BUFFER = Dot(N, Tdd{this->phin_Neumann, p1->phin_Neumann});
// 		this->phin_Dirichlet_BUFFER = Dot(N, Tdd{this->phin_Dirichlet, p1->phin_Dirichlet});
// 	}
// 	else if (face)
// 	{
// 		auto [p0, p1, p2] = face->getPoints();
// 		Tddd N = {t0, t1, 1. - t0 - t1};
// 		this->phiphin_BUFFER = Dot(N, T3Tdd{p0->phiphin, p1->phiphin, p2->phiphin});
// 		this->phiphin_t_BUFFER = Dot(N, T3Tdd{p0->phiphin_t, p1->phiphin_t, p2->phiphin_t});
// 		this->phi_Neumann_BUFFER = Dot(N, Tddd{p0->phi_Neumann, p1->phi_Neumann, p2->phi_Neumann});
// 		this->phi_Dirichlet_BUFFER = Dot(N, Tddd{p0->phi_Dirichlet, p1->phi_Dirichlet, p2->phi_Dirichlet});
// 		this->phin_Neumann_BUFFER = Dot(N, Tddd{p0->phin_Neumann, p1->phin_Neumann, p2->phin_Neumann});
// 		this->phin_Dirichlet_BUFFER = Dot(N, Tddd{p0->phin_Dirichlet, p1->phin_Dirichlet, p2->phin_Dirichlet});
// 	}
// };

/* ------------------------------------------------------ */

inline std::unordered_set<networkLine *> networkPoint::getLinesAround() const {
   std::unordered_set<networkLine *> ret;
   for (const auto &f : this->Faces)
      for_each(f->getLines(), [&](const auto &l) { ret.emplace(l); });
   // for (const auto &l : f->getLines())
   // 	ret.emplace(l);

   return ret;
};

inline std::unordered_set<networkLine *> networkPoint::getLinesOppsoite() const {
   return TakeExcept(getLinesAround(), this->Lines);
};

inline std::vector<std::tuple<networkFace *, Tddd>> networkPoint::getContactFacesX() const {
   // radiusで衝突をチェック
   std::vector<std::tuple<networkFace *, Tddd>> ret;
   Tddd x;
   for (const auto &f : this->getContactFaces()) {
      x = Nearest(this->X, ToX(f));
      if (this->radius > Norm(this->X - x))
         ret.push_back({f, x});
   }
   return ret;
};
inline std::vector<std::tuple<networkFace *, Tddd>> networkPoint::getContactFacesXCloser() const {
   auto c_faces_sorted = this->getContactFacesX();
   if (!c_faces_sorted.empty())
      std::sort(c_faces_sorted.begin(), c_faces_sorted.end(),
                [&](auto a, auto b) { return isFinite(Norm(std::get<1>(a) - this->X)) && isFinite(Norm(std::get<1>(b) - this->X)) &&
                                             (Norm(std::get<1>(a) - this->X) - Norm(std::get<1>(b) - this->X)) < 1E-20; });
   return c_faces_sorted;
};

/* ------------------------------------------------------ */

inline V_netFp networkPoint::getFacesNeumann() const {
   std::unordered_set<networkFace *> tmp;
   for (const auto &l : this->Lines)
      for (const auto &f : l->getFaces())
         if (f->Neumann)
            tmp.insert(f);
   return V_netFp(tmp.begin(), tmp.end());
};
inline V_netFp networkPoint::getFacesDirichlet() const {
   std::unordered_set<networkFace *> tmp;
   for (const auto &l : this->Lines)
      for (const auto &f : l->getFaces())
         if (f->Dirichlet)
            tmp.insert(f);
   return V_netFp(tmp.begin(), tmp.end());
};

/* ------------------------------------------------------ */
/*
この点に隣接する面の中で，f_INと向かいあう面があるかどうかを確認する
*/
// inline bool networkPoint::isThereAnyFacingFace(const networkFace *const f_IN, const double rad) const {
//    for (const auto &f : this->Faces)
//       if (f->isFacing(f_IN, rad))
//          return true;
//    return false;
// };
/* ------------------------------------------------------ */
inline void networkPoint::makeMirroredPoints(const Buckets<networkFace *> &B_face, const double mirroring_distance) {
   /*
           ! 衝突判定のためのmirroring_distanceは，影響半径の半分程度でいい．
           ! それより大きくてもそのミラー粒子は流体粒子の影響半径に入らない．
           ! また，ミラー粒子が流体粒子と同じネットワークにあると，Networkから全粒子を抽出する時に混ざった状態になり困ることがありそうなので，
           ! ミラー粒子は，衝突した面のネットワークに保存することにする．

           b! SPHにおいてミラーリングした粒子を参照する際の，注意：
           ! まず流体の方向を向いている面だけを抜き出し，それらの面でミラーリングされた影響半径内の粒子を取り出す．


           x	|   (x)
                   |  (o)
                   |/      (*)   OUT    |
                   /------------------- V
             o          *    IN

   点
      (x)  |   x
      (o)  |<--o
(*)       | / |    *   IN
                   ----V--------------  A
      (o)     (o)  (*)   OUT    |
                      (x)
(x)
   */
   this->map_Face_MirrorPoint.clear();
   auto x = this->X;
   for (const auto &f : DeleteDuplicates(Flatten(B_face.getObjects(this->X, mirroring_distance)))) {
      if (f->getNetwork() != this->getNetwork())
         if (Dot(f->incenter - x, f->normal) < 0 /*面と点が向き合っているかどうか*/) {
            if (IntersectQ(x, mirroring_distance, ToX(f)))
               this->map_Face_MirrorPoint[f] = new networkPoint(f->getNetwork(), 2. * (Nearest(x, ToX(f)) - x) + x);
         }
   }
};

/* ------------------------------------------------------ */
inline Tddd networkPoint::normalDirichletFace() const {
   // 接触面としての法線方向を計算する良い方法が必要．
   Tddd normal = {0., 0., 0.};
   int count = 0;
   for (auto &f : this->Faces) {
      auto [p0, p1, p2] = f->getPoints();
      if ((p0->Dirichlet || p0->CORNER) && (p1->Dirichlet || p1->CORNER) && (p2->Dirichlet || p2->CORNER))
         normal += f->getAngle(this) * f->normal;
   }
   if (count == 0)
      error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Dirichlet surface not found");
   return Normalize(normal);
};

inline Tddd networkPoint::normalNeumannFace() const {
   Tddd normal = {0., 0., 0.};
   int count = 0;
   for (auto &f : this->Faces) {
      auto [p0, p1, p2] = f->getPoints();
      if ((p0->Neumann || p0->CORNER) && (p1->Neumann || p1->CORNER) && (p2->Neumann || p2->CORNER)) {
         normal += f->getAngle(this) * f->normal;
      }
   }
   if (count == 0)
      error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Dirichlet surface not found");
   return Normalize(normal);
};

inline Tddd networkPoint::normalContanctSurface(const double pw0 = 1., const double pw1 = 1.) const {
   // 接触面としての法線方向を計算する良い方法が必要．
   double len = this->radius;
   Tddd n = {0, 0, 0};
   double totlen = 0;
   n = (1. / pow(len, pw0)) * this->getNormalTuple();
   totlen = 1. / pow(len, pw0);
   if (!this->getContactFaces().empty()) {
      std::unordered_set<networkFace *> contactfaces;
      for (auto &f : this->getContactFaces()) {
         bool duplicate = false;
         for (auto &cface : contactfaces)
            if (M_PI / 180. > MyVectorAngle(cface->normal, f->normal)) {
               duplicate = true;
               break;
            }
         if (!duplicate)
            contactfaces.emplace(f);
      }
      //
      for (const auto &f : contactfaces) {
         len = Norm(f->mirrorPosition(this) - this->X);
         n += (1. / pow(len, pw1)) * sgn(Dot(f->normal, this->getNormalTuple())) * f->normal;
         totlen += 1. / pow(len, pw1);
      }
   }
   return n / totlen;
};

// inline Tddd networkPoint::reflect(const Tddd &v) const {
//    // particlizeした点だけが使える．
//    auto f = std::get<0>(particlize_info);
//    if (f) {
//       auto n = f->normal;
//       return v - 2. * Dot(v, n) * n;
//    } else
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// };

// inline Tddd networkPoint::X_little_inside() const {
//    return this->X - (this->radius / 4.) * this->getNormalTuple();
// };
/*
@ addContactPointsは，引数として半径が与えられない場合，各点のradiusが成す球として，互いの重なりを調べ，重なった点はContactPointsに保存する．
*/

const double contact_angle = 40. * M_PI / 180.;
//

bool isInContact(const networkPoint *p, const T3Tddd &f_target) {
   bool isinradius = p->radius > Norm(p->X - Nearest(p->X, f_target));
   bool isCloseNormal = isFlat(p->getNormal_BEM(), -TriangleNormal(f_target), contact_angle) || isFlat(p->getNormal_BEM(), TriangleNormal(f_target), contact_angle);
   bool anyCloseNormal = std::any_of(p->Faces.begin(), p->Faces.end(), [&](const auto &F) { return isFlat(F->normal, -TriangleNormal(f_target), contact_angle) ||
                                                                                                   isFlat(F->normal, TriangleNormal(f_target), contact_angle); });
   return isinradius && (isCloseNormal || anyCloseNormal);
};

bool isInContact(const networkPoint *p, const Tddd &pX, const T3Tddd &f_target) {
   bool isinradius = p->radius > Norm(pX - Nearest(pX, f_target));
   bool isCloseNormal = isFlat(p->getNormal_BEM(), -TriangleNormal(f_target), contact_angle) || isFlat(p->getNormal_BEM(), TriangleNormal(f_target), contact_angle);
   bool anyCloseNormal = std::any_of(p->Faces.begin(), p->Faces.end(), [&](const auto &F) { return isFlat(F->normal, -TriangleNormal(f_target), contact_angle) ||
                                                                                                   isFlat(F->normal, TriangleNormal(f_target), contact_angle); });
   return isinradius && (isCloseNormal || anyCloseNormal);
};

bool isInContact(const Tddd &X /*base*/, const Tddd &n /*base*/,
                 const T3Tddd &f_target,
                 const double &range) {
   // nとf_targetの法線が近いかどうかを判定する．次に，最も近い点がrange以内にあるかどうかを判定する．
   auto N = TriangleNormal(f_target);
   if (!(isFlat(n, -N, contact_angle) || isFlat(n, N, contact_angle)))
      return false;  // not close normal!
   if (!(range > Norm(X - Nearest(X, f_target))))
      return false;  // not in range!
   return true;
};

// is f contact with f ?
bool isInContact(const networkPoint *p, const networkFace *f_normal, const networkFace *f_target) {
   return isInContact(p->X, f_normal->normal, ToX(f_target), p->radius);
};

// is p contact with f ?
bool isInContact(const networkPoint *p, const networkFace *f_target) {
   return std::any_of(p->Faces.begin(), p->Faces.end(), [&](const auto &f_normal) { return isInContact(p, f_normal, f_target); });
};

bool isInContact(const networkPoint *p, const networkFace *f_normal, const std::unordered_set<networkFace *> &faces_target) {
   return std::any_of(faces_target.begin(), faces_target.end(), [&](const auto &f_target) { return isInContact(p, f_normal, f_target); });
};

std::vector<networkFace *> selectionOfFaces(const networkPoint *const p,
                                            const std::unordered_set<networkFace *> &contactFacesIN,
                                            const int num,
                                            const bool delete_same_direction = true) {
   /*
   節点は，複数の三角形と接触していることがある．
   接触している三角形のうち境界条件として利用する必要がない三角形もある（法線方向が重複する三角形のうち利用すべき三角形は，最も近いものだけでよい）．
   */
   // auto X = p->X;
   // if (contactFacesIN.empty())
   //    return {};
   // else if (contactFacesIN.size() == 1)
   //    return {*contactFacesIN.begin()};
   // std::vector<networkFace *> contactFaces;
   // contactFaces.reserve(contactFacesIN.size());

   // for (const auto &f : contactFacesIN)
   //    if (isInContact(p, f))
   //       contactFaces.emplace_back(f);

   // std::sort(contactFaces.begin(), contactFaces.end(),
   //           [&](const auto &A, const auto &B) { return Norm(X - Nearest(X, ToX(A))) < Norm(X - Nearest(X, ToX(B))); });

   // std::vector<networkFace *> ret = {*contactFaces.begin()};
   // for (const auto &F : contactFaces)
   //    if (std::none_of(ret.begin(), ret.end(),
   //                     [&](const auto &f) { return isFlat(F->normal, -f->normal, M_PI / 180) ||
   //                                                 isFlat(F->normal, f->normal, M_PI / 180); }))
   //       ret.emplace_back(F);
   // return ret;
   /* -------------------------------------------------------------------------- */
   std::vector<std::tuple<networkFace *, double>> f_dist_sort;
   double dist;
   for (const auto &f : contactFacesIN)
      if (f->getNetwork() != p->getNetwork()) {
         // 別のネットワークに属するという条件
         dist = Norm(p->X - Nearest(p->X, ToX(f)));
         if (isInContact(p, f)) {
            std::tuple<networkFace *, double> T = {f, dist};
            if (f_dist_sort.empty())
               f_dist_sort.emplace_back(T);
            else {
               auto inserted = false;
               for (auto it = f_dist_sort.begin(); it != f_dist_sort.end(); ++it) {
                  if (std::get<1>(*it) >= dist) {
                     f_dist_sort.insert(it, T);
                     inserted = true;
                     break;
                  }
               }
               if (!inserted)
                  f_dist_sort.emplace_back(T);
            }
         }
      }

   std::vector<networkFace *> ret;

   if (!f_dist_sort.empty()) {
      for (const auto &[F, D] : f_dist_sort) {
         if (!delete_same_direction)
            ret.emplace_back(F);
         else if (std::none_of(ret.begin(), ret.end(), [&](const auto &f) { return isFlat(F->normal, -f->normal, M_PI / 180) ||
                                                                                   isFlat(F->normal, f->normal, M_PI / 180); }))
            ret.emplace_back(F);
         if (ret.size() >= num)
            return ret;
      };
   }
   return ret;
};

inline void networkPoint::addContactFaces(const Buckets<networkFace *> &B, bool include_self_network = true) {
   /*
    b% この衝突面の追加方法は，境界要素法用のもので，格子の情報を使って点の法線方向を計算し衝突を判断している．
   */
   Tddd x;
   if (this->Lines.empty()) {
      // b!まずは，衝突があり得そうな面を多めに保存する．
      std::unordered_map<networkFace *, Tddd> tmpContactFaces;
      // for (const auto &f : B.getObjects_unorderedset(this->X, 2. * this->radius))
      for (const auto &f : B.getAll()) {
         // auto [p0, p1, p2] = f->getXVertices();
         // auto intxn = IntersectionSphereTriangle(this->X, 3. * this->radius, f->getXVertices());
         // if (Norm(intxn.X - this->X) < 2. * this->radius)
         //    tmpContactFaces[f] = intxn.X;
         //
         x = Nearest(this->X, ToX(f));
         if (3 * this->radius > Norm(this->X - x))
            tmpContactFaces[f] = x;
      }
      // b!各面について衝突があり得るか詳しく調べて判断する．
      for (const auto &[f, intxnX] : tmpContactFaces) {
         //@ 周りの点に隠れて面が見えない状況の場合，その面とは接し得ないので，考慮しない．
         // auto [p0, p1, p2] = f->getPointsTuple();
         // BEMのために導入したもの．周辺の点の方向にある面には衝突し得ないため，無視する．
         // ただし修正面の方向をむくせんを抜き出す必要がある．
         // auto n = f->normal;
         auto thisPointToFace = intxnX - this->X;
         // b$ 点が面を有していない場合（SPH用）
         if (Dot(f->normal, thisPointToFace /* <--この点からfまでのベクトルを向き合いの判定に利用*/) < 0. /*少なくとも向き合っていなければいけない*/)
            this->ContactFaces.emplace(f);
      }
   } else {

      // this->ContactFaces = ToUnorderedSet(selectionOfFaces(this, B.getObjects_unorderedset(this->X, this->radius), 20, false));

      DebugPrint("! まずは，衝突があり得そうな面を多めに保存する．");
      double dist;
      std::unordered_set<networkFace *> faces = B.getObjects_unorderedset(this->X, this->radius);
      faces.insert(this->ContactFaces.begin(), this->ContactFaces.end());
      std::vector<std::tuple<networkFace *, double>> f_dist_sort;
      for (const auto &f : faces)
         // 別のネットワークに属するという条件
         if (include_self_network || f->getNetwork() != this->getNetwork()) {
            dist = Norm(this->X - Nearest(this->X, ToX(f)));
            // 面同士が向き合っているという条件
            if (isInContact(this, f)) {
               std::tuple<networkFace *, double> T = {f, dist};
               if (f_dist_sort.empty())
                  f_dist_sort.push_back(T);
               else {
                  auto inserted = false;
                  for (auto it = f_dist_sort.begin(); it != f_dist_sort.end(); ++it) {
                     if (std::get<1>(*it) >= dist) {
                        f_dist_sort.insert(it, T);
                        inserted = true;
                        break;
                     }
                  }
                  if (!inserted)
                     f_dist_sort.push_back(T);
               }
            }
         }

      if (!f_dist_sort.empty()) {
         DebugPrint("! 異なる方向を向く面の情報だけが欲しいので，同方向の面は無視する");
         for (const auto &[F, D] : f_dist_sort) {
            // if (std::none_of(this->ContactFaces.begin(),
            //                  this->ContactFaces.end(),
            //                  [&](const auto &f) { return isFlat(F->normal, -f->normal, M_PI / 180) ||
            //                                              isFlat(F->normal, f->normal, M_PI / 180); }))
            this->ContactFaces.emplace(F);
            if (this->ContactFaces.size() > 10) {
               return;
            }
         };
      } else
         DebugPrint("faces is empty");
   }
};

// 点なのに，sphereでインターセクションをチェックするのは効率的ではない．
inline void networkPoint::addContactPoints(const Buckets<networkPoint *> &B,
                                           const bool include_self_network = true) {
   for (const auto &q : B.getObjects_unorderedset(this->X, 2. * this->radius /*depth*/))
      if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
         if (this != q) {
            // if (intersection(geometry::Sphere(this->X, this->radius), geometry::Sphere(this->X, q->radius)).isIntersecting)
            if (IntersectQ(this->X, this->radius, q->X, q->radius))
               this->ContactPoints.emplace(q);
         }
};
// radiusのしていがある場合は，その半径内の点を取得する．
inline void networkPoint::addContactPoints(const Buckets<networkPoint *> &B,
                                           const double radius,
                                           const bool include_self_network = true) {
   for (const auto &q : B.getObjects_unorderedset(this->X, 2. * this->radius /*depth*/))
      if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
         if (this != q)
            if (Norm(this->X - q->X) <= radius)
               this->ContactPoints.emplace(q);
};
inline void networkPoint::addContactPoints(const Buckets<networkPoint *> &B,
                                           const int limit_depth,
                                           const int limit_num,
                                           const bool include_self_network = true) {
   for (const auto &q : B.getObjects_unorderedset(this->X, limit_depth, limit_num))
      if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
         if (this != q)
            if (Norm(this->X - q->X) <= radius)
               this->ContactPoints.emplace(q);
};
// radiusのしていがある場合は，その半径内の点を取得する．
inline void networkPoint::addContactPoints(const std::vector<Buckets<networkPoint *>> &Bs,
                                           const double radius,
                                           const bool include_self_network = true) {
   for (const auto &B : Bs)
      addContactPoints(B, radius, include_self_network);
};
inline void networkPoint::addContactPoints(const std::vector<Buckets<networkPoint *>> &Bs,
                                           const int limit_depth,
                                           const int limit_num,
                                           const bool include_self_network = true) {
   for (const auto &B : Bs)
      addContactPoints(B, limit_depth, limit_num, include_self_network);
};

// 追加2021/06/18
V_netFp networkPoint::getFacesSort() const {
   V_netFp ret = {*this->Faces.begin()};
   netFp f;
   netLp l;
   int count = 0;
   /*
    *        *-------*
    *       / \    /  \
    *      / f \  /    \
    *     *-- l -0------*
    *      \    / \    /
    *       \  /   \  /
    *        *-------*
    */
   do {
      f = *ret.rbegin();
      l = std::get<2>(f->getLinesTupleFrom(this));
      if ((*l)(f) != ret[0])
         ret.emplace_back((*l)(f));
      else
         break;
   } while (count++ < 1000);
   return ret;
};
/* ------------------------------------------------------ */
// 追加2021/09/02
inline V_netFp networkPoint::getFacesSort(networkLine *const line) const {
   try {
      auto p1 = (*line)(this);
      if (!p1)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      networkFace *begin_face;
      // 左回りの面を選ぶ．その面からスタート
      for (const auto &f : line->getFaces())
         if (f->getPointFront(line) == p1) {
            begin_face = f;
            break;
         }

      V_netFp ret = {begin_face};
      ret.reserve(10);
      int c = 0;
      networkLine *next_line = line;
      do {
         next_line = begin_face->getLineBack(next_line);
         if (next_line != line)
            ret.emplace_back(begin_face = (*next_line)(begin_face));
         else
            break;
      } while (c++ < 1000);
      if (c >= 1000)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      return ret;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
/* ------------------------------------------------------ */
inline V_netFp networkPoint::getFaces(networkLine *line) const {
   return getFacesSort(line);
};
V_netPp networkPoint::getNeighborsSort2() const {
   V_netPp ret;
   netPp p;
   netLp l;
   for (const auto &f : getFacesSort()) {
      p = std::get<1>(f->getPoints(this));
      ret.emplace_back(p);
      l = f->getLineOpposite(this);
      p = (*l)(f)->getPointOpposite(l);
      ret.emplace_back(p);
   }
   return ret;
};

inline V_d networkPoint::getAngles() const {
   V_d ret(this->Faces.size());
   int i = 0;
   for (const auto &f : this->Faces)
      ret[i++] = f->getAngle(this);
   return ret;
};

inline V_d networkPoint::getAngles(networkLine *const base_line) const {
   if (!base_line)
      return this->getAngles();
   V_d ret(this->Faces.size());
   int i = 0;
   for (const auto &f : this->getFacesSort(base_line))
      ret[i++] = f->getAngle(this);
   return ret;
};

using V_Tup_ddVdPp = std::vector<std::tuple<double, double, V_d, networkPoint *>>;

//! 起点を指定できるようにする．
using V_Tup_ddVd = std::vector<std::tuple<double, double, V_d>>;

using V_Var_VVdVPp = std::vector<std::variant<VV_d, V_netPp>>;

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
inline void networkPoint::resetXinfo() {
   if (this->xline != nullptr)
      this->xline->clearXPoints(this);
   this->xline = nullptr;
   if (this->xface != nullptr)
      this->xface->clearXPoints(this);
   this->xface = nullptr;
};

inline void networkPoint::setX(const Tddd &xyz_IN) {
   try {
      CoordinateBounds::setBounds(xyz_IN);
      for (const auto &l : this->getLines()) {
         // std::cout << "l = " << l << std::endl;
         l->setBoundsSingle();
      }
      for (const auto &f : this->Faces) {
         // std::cout << "f = " << f << std::endl;
         // std::cout << f->Points << std::endl;
         // std::cout << f->Lines << std::endl;
         f->setGeometricProperties(ToX(f->setPoints(f->Lines)));
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
// setXは，接続するlineやfaceのsetBoundsを行う．
inline void networkPoint::setX(const V_d &xyz_IN) {
   this->setX(Tddd{xyz_IN[0], xyz_IN[1], xyz_IN[2]});
};

inline Tddd networkPoint::getNormalTuple() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      normal += f->getAngle(this) * f->normal;
   return Normalize(normal);

   // std::vector<Tddd> normals({});
   // for (const auto &f : this->Faces)
   // 	normals.emplace_back(f->normal);

   // if (normals.empty())
   // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "normals is empty");
   // else
   // 	return Mean(normals);
};
inline std::vector<double> networkPoint::getFaceAreas() const {
   V_d ret(this->Faces.size());
   int i = 0;
   for (const auto &f : this->Faces)
      ret[i++] = f->area;
   return ret;
};

#define linear_area_averaged
inline Tddd networkPoint::getNormalAreaAveraged() const {
#ifdef linear_area_averaged
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      normal += f->area * f->normal;
#elif defined(quad_area_averaged)

   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces) {
      auto [p0, p1, p2] = f->getPointsTuple(this);
      auto l0 = p0->getLineBetween(p1);
      auto l1 = p1->getLineBetween(p2);
      auto l2 = p2->getLineBetween(p0);
      //
      auto fs0 = l0->getFaces();
      auto ps6_l0_f00 = fs0[0]->get6PointsTuple(l0);
      auto ps6_l0_f01 = fs0[1]->get6PointsTuple(l0);
      interpolationTriangleQuadByFixedRange3D intp_l0_0(ToX(ps6_l0_f00));
      interpolationTriangleQuadByFixedRange3D intp_l0_1(ToX(ps6_l0_f01));
      //
      auto fs1 = l1->getFaces();
      auto ps6_l1_f10 = fs1[0]->get6PointsTuple(l1);
      auto ps6_l1_f11 = fs1[1]->get6PointsTuple(l1);
      interpolationTriangleQuadByFixedRange3D intp_l1_0(ToX(ps6_l1_f10));
      interpolationTriangleQuadByFixedRange3D intp_l1_1(ToX(ps6_l1_f11));
      //
      auto fs2 = l2->getFaces();
      auto ps6_l2_f20 = fs2[0]->get6PointsTuple(l2);
      auto ps6_l2_f21 = fs2[1]->get6PointsTuple(l2);
      interpolationTriangleQuadByFixedRange3D intp_l2_0(ToX(ps6_l2_f20));
      interpolationTriangleQuadByFixedRange3D intp_l2_1(ToX(ps6_l2_f21));
      //
      auto [f0p0, f0p1, f0p2] = (*l0)(f)->getPoints();
      auto [f1p0, f1p1, f1p2] = (*l1)(f)->getPoints();
      auto [f2p0, f2p1, f2p2] = (*l2)(f)->getPoints();
      normal += f->area * f->normal;
   }
#endif
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                       Arithmetic                       */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalArithmeticAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      normal += f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletArithmeticAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannArithmeticAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*              Kernel Weighted Averaged                  */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalSplineKernelAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   double radius = 1.5 * Norm(extX(this->getNeighbors()) - this->X);
   for (const auto &f : this->Faces)
      normal += w_Bspline5(Norm(f->getXtuple() - this->X), radius) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletSplineKernelAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   double radius = 1.5 * Norm(extX(this->getNeighbors()) - this->X);
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += w_Bspline5(Norm(f->getXtuple() - this->X), radius) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannSplineKernelAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   double radius = 1.5 * Norm(extX(this->getNeighbors()) - this->X);
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += w_Bspline5(Norm(f->getXtuple() - this->X), radius) * f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                   SubArea Averaged                     */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalSubAreaAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.}, a, b, c;
   for (const auto &f : this->Faces)
      normal += f->getSubArea(this) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletSubAreaAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getSubArea(this) * f->normal;

   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannSubAreaAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->getSubArea(this) * f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                      Area Averaged                     */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalDirichletAreaAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->area * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannAreaAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->area * f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                  Area Averaged Buffer                  */
/* ------------------------------------------------------ */
// inline Tddd networkPoint::getNormalAreaAveraged_Buffer() const {
//    // 角度の重みを掛けた法線ベクトルを足し合わせる
//    Tddd normal = {0., 0., 0.};
//    for (const auto &f : this->Faces)
//       if (f->Dirichlet)
//          normal += f->getAreaBuffer() * f->getNormalBuffer();
//    return Normalize(normal);
// };
// inline Tddd networkPoint::getNormalDirichletAreaAveraged_Buffer() const {
//    // 角度の重みを掛けた法線ベクトルを足し合わせる
//    Tddd normal = {0., 0., 0.};
//    for (const auto &f : this->Faces)
//       if (f->Dirichlet)
//          normal += f->getAreaBuffer() * f->getNormalBuffer();
//    return Normalize(normal);
// };
// inline Tddd networkPoint::getNormalNeumannAreaAveraged_Buffer() const {
//    // 角度の重みを掛けた法線ベクトルを足し合わせる
//    Tddd normal = {0., 0., 0.};
//    for (const auto &f : this->Faces)
//       if (f->Neumann)
//          normal += f->getAreaBuffer() * f->getNormalBuffer();
//    return Normalize(normal);
// };
/* ------------------------------------------------------ */
/*             Inscribed Circle Area Averaged             */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalDirichletInscribedCircleAreaAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getInscribedCircleArea() * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannInscribedCircleAreaAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->getInscribedCircleArea() * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalInscribedCircleAreaAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      normal += f->getInscribedCircleArea() * f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                     Angle Averaged                     */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalDirichlet() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getAngle(this) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumann() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->getAngle(this) * f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                     Angle Averaged                     */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalAngleAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      normal += f->getAngle(this) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletAngleAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getAngle(this) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannAngleAveraged() const {
   // 角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->getAngle(this) * f->normal;
   return Normalize(normal);
}
/* ------------------------------------------------------ */
/*                     Newton method                      */
/* ------------------------------------------------------ */
#include "rootFinding.hpp"
//
Tddd optimumUnitVector(const std::vector<Tddd> &unit_vectors, const Tddd &init_n) {
   //@ 最小にする関数
   //@ F0 = ((x^2+y^2+z^2)-1)^2
   //@ dF0dx = 2*((x^2+y^2+z^2)-1)*2x = 4*x*(x^2+y^2+z^2-1)
   //@ F1 = (nx-x)^2+(ny-y)^2+(nz-z)^2
   //@ dF1dx = -2(nx-x)
   NewtonRaphson NR0(std::get<0>(init_n)), NR1(std::get<1>(init_n)), NR2(std::get<2>(init_n));
   double r, F0, F1, F2, dF0dx, dF1dx, dF2dx, drdx, w;
   for (auto i = 0; i < 100; ++i) {
      // r = NR0.X * NR0.X + NR1.X * NR1.X + NR2.X * NR2.X - 1;
      // F0 = F1 = F2 = r * r;
      // dF0dx = 4 * NR0.X * r;
      // dF1dx = 4 * NR1.X * r;
      // dF2dx = 4 * NR2.X * r;
      F0 = F1 = F2 = dF0dx = dF1dx = dF2dx = 0;
      for (const auto &u_vec : unit_vectors) {
         w = 1;  // kernel_Bspline3(Norm(X - p->getXtuple()), p->radius);
         r = w * (std::get<0>(u_vec) - NR0.X);
         drdx = -w;
         F0 += 2. * r * drdx;  //<- d/dx (d*d)
         dF0dx += 2. * drdx * drdx;
         r = w * (std::get<1>(u_vec) - NR1.X);
         drdx = -w;
         F1 += 2. * r * drdx;  //<- d/dx (d*d)
         dF1dx += 2. * drdx * drdx;
         r = w * (std::get<2>(u_vec) - NR2.X);
         drdx = -w;
         F2 += 2. * r * drdx;  //<- d/dx (d*d)
         dF2dx += 2. * drdx * drdx;
      }
      NR0.update(F0, dF0dx);
      NR1.update(F1, dF1dx);
      NR2.update(F2, dF2dx);
      if (std::abs(NR0.dX) < 1E-8 && std::abs(NR1.dX) < 1E-8 && std::abs(NR2.dX) < 1E-8)
         break;

      auto [a, b, c] = Normalize(Tddd{NR0.X, NR1.X, NR2.X});
      NR0.X = a;
      NR1.X = b;
      NR2.X = c;
   }
   return Normalize(Tddd{NR0.X, NR1.X, NR2.X});
};

inline Tddd networkPoint::getNormalOptimum() const {
   std::vector<Tddd> unit_normals;
   for (const auto &f : this->Faces)
      unit_normals.emplace_back(f->normal);
   return optimumUnitVector(unit_normals, this->getNormalTuple());
};
inline Tddd networkPoint::getNormalDirichletOptimum() const {
   std::vector<Tddd> unit_normals;
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         unit_normals.emplace_back(f->normal);
   return optimumUnitVector(unit_normals, this->getNormalDirichlet());
};
inline Tddd networkPoint::getNormalNeumannOptimum() const {
   std::vector<Tddd> unit_normals;
   for (const auto &f : this->Faces)
      if (f->Neumann)
         unit_normals.emplace_back(f->normal);
   return optimumUnitVector(unit_normals, this->getNormalNeumann());
}
/* ------------------------------------------------------ */
// this->normal_BEM = {0., 0., 0.};
// this->normal_Dir_BEM = {0., 0., 0.};
// this->normal_Neu_BEM = {0., 0., 0.};
// double Qtot = 0, Qtot_Dir = 0, Qtot_Neu = 0;
// for (const auto &f : this->Faces)
// {
// 	auto Q = f->getAngle(this);
// 	Qtot += Q;
// 	auto [l0, l, l1] = f->getLinesTupleFrom(this);
// 	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_f_l(f, l);
// 	auto n = Normalize(intp_f_l.cross(0, .5)) * Q;
// 	this->normal_BEM += n;
// 	if (f->Dirichlet)
// 	{
// 		Qtot_Dir += Q;
// 		this->normal_Dir_BEM += n;
// 	}
// 	else if (f->Neumann)
// 	{
// 		Qtot_Neu += Q;
// 		this->normal_Neu_BEM += n;
// 	}
// }
// this->normal_BEM /= Qtot;
// this->normal_Dir_BEM /= Qtot_Dir;
// this->normal_Neu_BEM /= Qtot_Neu;
inline Tddd networkPoint::getNormalQuadInterpAngleAveraged() const {
   Tddd ret = {0., 0., 0.};
   double Qtot = 0, Qtot_Dir = 0, Qtot_Neu = 0;
   for (const auto &f : this->Faces) {
      auto Q = f->getAngle(this);
      Qtot += Q;
      auto [l0, l, l1] = f->getLinesTupleFrom(this);
      interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_f_l(f, l);
      ret += Normalize(intp_f_l.cross(0, .5)) * Q;
   }
   return ret / Qtot;
};
inline Tddd networkPoint::getNormalDirichletQuadInterpAngleAveraged() const {
   Tddd ret = {0., 0., 0.};
   double Qtot = 0, Qtot_Dir = 0, Qtot_Neu = 0;
   for (const auto &f : this->Faces)
      if (f->Dirichlet) {
         auto Q = f->getAngle(this);
         Qtot += Q;
         auto [l0, l, l1] = f->getLinesTupleFrom(this);
         interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_f_l(f, l);
         ret += Normalize(intp_f_l.cross(0, .5)) * Q;
      }
   return ret / Qtot;
};
inline Tddd networkPoint::getNormalNeumannQuadInterpAngleAveraged() const {
   Tddd ret = {0., 0., 0.};
   double Qtot = 0, Qtot_Dir = 0, Qtot_Neu = 0;
   for (const auto &f : this->Faces)
      if (f->Neumann) {
         auto Q = f->getAngle(this);
         Qtot += Q;
         auto [l0, l, l1] = f->getLinesTupleFrom(this);
         interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_f_l(f, l);
         ret += Normalize(intp_f_l.cross(0, .5)) * Q;
      }
   return ret / Qtot;
};
//@ ------------------------------------------------------ */
//@              BEMでどちらを法線ベクトルとして使うか           */
//@ ------------------------------------------------------ */

// #define use_kernel_weigted_normal
// #define use_angle_weigted_normal
#define use_area_weigted_normal
// #define use_subarea_weigted_normal
// #define use_inscribedcircle_area_weigted_normal
// #define use_arithmetic_averaged_normal
// #define use_optimum_normal

inline Tddd networkPoint::getNormalDirichlet_BEM() const {
   //@ これはEMTに使うので，やはりAngle Averagedがいい
#if defined(use_angle_weigted_normal)
   return getNormalDirichletAngleAveraged();
#elif defined(use_area_weigted_normal)
   return getNormalDirichletAreaAveraged();
#elif defined(use_subarea_weigted_normal)
   return getNormalDirichletSubAreaAveraged();
#elif defined(use_inscribedcircle_area_weigted_normal)
   return getNormalDirichletInscribedCircleAreaAveraged();
#elif defined(use_arithmetic_averaged_normal)
   return getNormalDirichletArithmeticAveraged();
#elif defined(use_kernel_weigted_normal)
   return getNormalDirichletSplineKernelAveraged();
#elif defined(use_optimum_normal)
   return getNormalDirichletOptimum();
#endif
};

inline Tddd networkPoint::getNormalNeumann_BEM() const {
   //@ これはEMTに使うので，やはりAngle Averagedがいい
#if defined(use_angle_weigted_normal)
   return getNormalNeumannAngleAveraged();
#elif defined(use_subarea_weigted_normal)
   return getNormalNeumannSubAreaAveraged();
#elif defined(use_area_weigted_normal)
   return getNormalNeumannAreaAveraged();
#elif defined(use_inscribedcircle_area_weigted_normal)
   return getNormalNeumannInscribedCircleAreaAveraged();
#elif defined(use_arithmetic_averaged_normal)
   return getNormalNeumannArithmeticAveraged();
#elif defined(use_kernel_weigted_normal)
   return getNormalNeumannSplineKernelAveraged();
#elif defined(use_optimum_normal)
   return getNormalNeumannOptimum();
#endif
};
inline Tddd networkPoint::getNormal_BEM() const {
   if (this->CORNER)
      return Normalize((getNormalDirichlet_BEM() + getNormalNeumann_BEM()) / 2.);
   else {
#if defined(use_angle_weigted_normal)
      return getNormalAngleAveraged();
#elif defined(use_subarea_weigted_normal)
      return getNormalSubAreaAveraged();
#elif defined(use_area_weigted_normal)
      return getNormalAreaAveraged();
#elif defined(use_inscribedcircle_area_weigted_normal)
      return getNormalInscribedCircleAreaAveraged();
#elif defined(use_arithmetic_averaged_normal)
      return getNormalArithmeticAveraged();
#elif defined(use_kernel_weigted_normal)
      return getNormalSplineKernelAveraged();
#elif defined(use_optimum_normal)
      return getNormalOptimum();
#endif
   }
};

inline Tddd networkPoint::getNormalDirichlet_BEM_Buffer() const {
   return getNormalDirichletAreaAveraged_Buffer();
};

inline Tddd networkPoint::getNormalNeumann_BEM_Buffer() const {
   return getNormalNeumannAreaAveraged_Buffer();
};
inline Tddd networkPoint::getNormal_BEM_Buffer() const {
   if (this->CORNER)
      return Normalize((getNormalDirichlet_BEM_Buffer() + getNormalNeumann_BEM_Buffer()) / 2.);
   else {
      return getNormalAreaAveraged_Buffer();
   }
};

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
// inline V_d networkPoint::getNormal() const
// {
// 	//角度の重みを掛けた法線ベクトルを足し合わせるls
// 	auto normal = this->getNormalTuple();
// 	return ToVector(normal);
// 	// VV_d normals({});
// 	// for (const auto &f : this->Faces)
// 	// 	normals.emplace_back(f->getNormal());

// 	// if (normals.empty())
// 	// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "normals is empty");
// 	// else
// 	// 	return Mean(normals);
// };

/*SolidAngle_detail

  SolidAngle_detail*/
/*SolidAngle_detail_code*/
/*角点は周囲に点があるので，取得する点を制限しなければならない*/
// inline double networkPoint::getSolidAngle(bool TorF)
// {
// 	auto fs = this->getFaces();
// 	VV_d ns(fs.size(), V_d(3, 0.));
// 	int i = 0;
// 	for (const auto &f : fs)
// 		ns[i++] = f->getNormal();
// 	V_d mean = ns[0], var = {0., 0., 0.};
// 	for (auto i = 0; i < ns.size(); i++)
// 		var += (ns[i] - mean) * (ns[i] - mean);
// 	double s = std::sqrt(Norm(var));
// 	if (s < 1E-5)
// 		return 2. * M_PI;

// 	return SolidAngle(this->getX(), obj3D::extractX(this->getNeighborsSort(TorF)));
// 	;
// };
// inline double networkPoint::getSolidAngle(const V_netFp &faces /*available faces*/)
// {
// 	V_netFp fs({});
// 	for (const auto &f : this->Faces)
// 		if (MemberQ(faces, f))
// 			fs.emplace_back(f);
// 	if (fs.size() < 3)
// 		return 0;

// 	VV_d ns(fs.size(), V_d(3, 0.));
// 	int i = 0;
// 	for (const auto &f : fs)
// 		ns[i++] = f->getNormal();
// 	V_d mean = ns[0], var = {0., 0., 0.};
// 	for (auto i = 0; i < ns.size(); i++)
// 		var += (ns[i] - mean) * (ns[i] - mean);
// 	double s = std::sqrt(Norm(var));
// 	if (s < 1E-5)
// 		return 2. * M_PI;

// 	return SolidAngle(this->getX(), obj3D::extractX(fs));
// 	;
// };
/*SolidAngle_detail_code*/
/*getNeighborsSort_detail
networkFaceのPointsが反時計周りに並んで保存されていれば，
その並びに沿って，`getNeighborsSort`も必ず反時計周りの多角形をなす点ベクトルを返す．
getNeighborsSort_detail*/
/*getNeighborsSort_sort*/
inline V_netPp networkPoint::getNeighborsSort() const {
   try {
      VV_netPp Vps({});
      // V_netPp qs(0);
      for (const auto &f : this->Faces) {
         // qs = f->getPoints();
         // if (qs[0] == this)
         // 	Vps.push_back(V_netPp{qs[1], qs[2]});
         // else if (qs[1] == this)
         // 	Vps.push_back(V_netPp{qs[2], qs[0]});
         // else if (qs[2] == this)
         // 	Vps.push_back(V_netPp{qs[0], qs[1]});
         // else
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "ps is empty");
         // これは以下で実現できる．変更2021/06/18
         auto [q0, q1, q2] = f->getPoints(this);
         Vps.push_back(V_netPp{q1, q2});
      }
      auto ret = FlattenAsChain(Vps);
      if (*ret.begin() == *ret.rbegin())
         ret.pop_back();
      else {
         return ret;
         // この場合，
         //  * thisが端点であるか
         //  * 周りが数珠つなぎとなっていない（エラー）
         //  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not chain !");
      }
      return ret;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
//
inline V_netPp networkPoint::getNeighborsSort(bool TorF) {
   try {
      VV_netPp Vps({});
      // V_netPp qs(0);
      for (const auto &f : this->getFaces(TorF)) {
         auto [q0, q1, q2] = f->getPoints();
         if (q0 == this)
            Vps.push_back(V_netPp{q1, q2});
         else if (q1 == this)
            Vps.push_back(V_netPp{q2, q0});
         else if (q2 == this)
            Vps.push_back(V_netPp{q0, q1});
         else
            throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "ps is empty"));
      }
      auto ret = FlattenAsChain(Vps);
      if (*ret.begin() == *ret.rbegin())
         ret.pop_back();
      else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not chain !");
      return ret;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
/*getNeighborsSort_sort*/
inline V_netPp networkPoint::getXNeighbors() const {
   V_netPp ret({});
   for (const auto &l : this->getLines()) {
      netPp closestXp = nullptr;
      double minlen = 1E+100;
      for (const auto &p : l->getXPoints())
         if (minlen > distance(this, p)) {
            closestXp = p;
            minlen = distance(this, p);
         }
      if (closestXp)
         ret.emplace_back(closestXp);
   }
   return ret;
};

///////////// networkPoint //////////////
// inline networkPoint::networkPoint(Network *network_IN, const V_d &xyz_IN, networkLine *line_IN, networkFace *face_IN)
//     : networkObject(network_IN),
//       CoordinateBounds(xyz_IN),
//       xline(line_IN),
//       xface(face_IN)
// {
//     // std::cout << "new networkPoint ... ";
//     networkObject::storage = network_IN; //なぜか初期化リストに入れれない
//     // std::cout << "done" << std::endl;
// };
// コンストラクタ
inline networkPoint::networkPoint(Network *network_IN, const Tddd &xyz_IN, networkLine *xline_IN, networkFace *xface_IN)
    : network(network_IN),
      CoordinateBounds(xyz_IN /*CoordinateBoundsに暗黙に置換される*/),
      initialX(xyz_IN),
      xline(xline_IN),
      xface(xface_IN),
      force({0., 0., 0., 0., 0., 0.}),
      inertia({0., 0., 0., 0., 0., 0.}),
      acceleration({0., 0., 0., 0., 0., 0.}),
      velocity({0., 0., 0., 0., 0., 0.}),
      mass(0.),
      density(0.),
      volume(0.),
      radius(1.),
      signed_distance_vector({0., 0., 0.}),
      map_Net_ContactPoints({{nullptr, {}}}),
/* ------------------------------------------------------ */
#ifdef BEM
      phiphin({0., 0.}),
      phiphin_t({0., 0.}),
      //   phi_Neumann(0.),
      phi_Dirichlet(0.),
      //   phin_Neumann(0.),
      phin_Dirichlet(0.),
      grad_phi_BEM({0., 0., 0.}),
      U_BEM({0., 0., 0.}),
      U_BEM_last({0., 0., 0.}),
      U_update_BEM({0., 0., 0.}),
      U_tangential_BEM({0., 0., 0.}),
      U_tangential_BEM_last({0., 0., 0.}),
      U_normal_BEM({0., 0., 0.}),
      laplacian_U_BEM({0., 0., 0.}),
      normal_BEM({0., 0., 0.}),
      pressure_BEM(0.),
      U_cling_to_Neumann({0., 0., 0.}),
      U_cling_to_Neumann_({0., 0., 0.}),
      grid_tension({0., 0., 0.}),
      grid_tension_({0., 0., 0.}),
      RK_phi(),
      RK_X(),
#endif
      /* ------------------------------------------------------ */
      particlize_info({nullptr, {nullptr, nullptr, nullptr}, {0., 0.}, 0., 0.}) {

   // std::cout << "new networkPoint ... ";

   if (xline != nullptr)
      xline->addXPoint(this);
   if (xface != nullptr)
      xface->addXPoint(this);

   // this->storage = storage_IN;  //なぜか初期化リストに入れれない
   this->network->add(this);
   // std::cout << "done" << std::endl;
   /* ------------------------------------------------------ */
#ifdef DEM
   this->contactP = {};
   this->mass = 1.;
   this->radius = 1.;
   this->W = 0.;
   //! SPH
   this->a_viscosity = 0.;
   this->mu_SPH = _WATER_MU_10deg_;  // 10 deg
   // 0.001140;  // water 15 deg
   this->DUDt_SPH = {0., 0., 0.};
   this->repulsive_force_SPH = {0., 0., 0.};
   this->DPDt_SPH = 0.;
   this->DrhoDt_SPH = 0.;
   this->density = 0.;
   this->lap_U = {0., 0., 0.};
   this->div_U = 0;
   this->gradP_SPH = {0., 0., 0.};
   this->U_SPH = {0., 0., 0.};
   this->mu_lap_rho_g_SPH = {0., 0., 0.};
   this->tmp_U_SPH = {0., 0., 0.};
   this->interpolated_normal_SPH = {0., 0., 0.};
   this->pn_is_set = false;
   this->radius_SPH = 0.;
   this->pressure_SPH = 0.;
   this->pressure_SPH_ = 0.;
   this->isSurface = false;
   this->isCaptured = false;
#endif

   this->Dirichlet = false;
   this->Neumann = false;
};
// 逆方向の立体角
inline double networkPoint::getSolidAngle() const {
   return SolidAngle(this->X, extX(this->getNeighborsSort()));
};
// inline double networkPoint::getSolidAngleBuffer() const {
//    return SolidAngle(this->getXBuffer(), extXBuffer(this->getNeighborsSort()));
// };

/////////////////////////////////////////

std::vector<bool> isIntxns(const V_netLp &ls) {
   std::vector<bool> ret(ls.size(), false);
   for (auto i = 0; i < ls.size(); i++)
      if (ls[i]->isIntxn())
         ret[i] = true;
   return ret;
};

std::vector<bool> isXPoints(const V_netPp &ps) {
   std::vector<bool> ret(ps.size(), false);
   for (auto i = 0; i < ps.size(); i++)
      if (ps[i]->isXPoint())
         ret[i] = true;
   return ret;
};

std::vector<std::vector<bool>> isXPoints(const VV_netPp &ps) {
   std::vector<std::vector<bool>> ret(ps.size());
   for (auto i = 0; i < ps.size(); i++)
      ret[i] = isXPoints(ps[i]);
   return ret;
};

V_netLp getXLines(const V_netPp &ps) {
   V_netLp ret(ps.size());
   for (auto i = 0; i < ps.size(); i++)
      ret[i] = ps[i]->getXLine();
   return ret;
};

VV_netLp getXLines(const VV_netPp &ps) {
   VV_netLp ret(ps.size());
   for (auto i = 0; i < ps.size(); i++)
      ret[i] = getXLines(ps[i]);
   return ret;
};

#endif