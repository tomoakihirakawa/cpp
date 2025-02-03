#pragma once

/* --------------------------------- SUFACE --------------------------------- */

inline std::vector<networkFace *> networkPoint::getSurfaces() const {
   std::vector<networkFace *> surfaces;
   surfaces.reserve(this->Faces.size());
   for (const auto &f : this->Faces)
      if (f->SurfaceQ())
         surfaces.emplace_back(f);
   return surfaces;
};

inline std::vector<networkFace *> networkPoint::setFaces(const V_netLp &this_lines) {
   std::unordered_set<networkFace *> ret;
   for (const auto &l : this_lines)
      if (l != nullptr)
         for (const auto &f : l->getFaces())
            ret.emplace(f);
   this->Faces.assign(std::make_move_iterator(ret.begin()), std::make_move_iterator(ret.end()));
   return this->Faces;
};

inline std::vector<networkFace *> networkPoint::setFaces() {
   setFaces(this->Lines);
   return this->Faces;
};

/* ------------------------------------------------------ */

inline std::unordered_set<networkLine *> networkPoint::getLinesAround() const {
   std::unordered_set<networkLine *> ret;
   for (const auto &f : this->Faces)
      std::ranges::for_each(f->getLines(), [&](const auto &l) { ret.emplace(l); });
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
      if (this->detection_range > Norm(this->X - x))
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
   for (const auto &f : this->Faces)
      if (f->Neumann)
         tmp.emplace(f);
   return V_netFp(tmp.begin(), tmp.end());
};
inline V_netFp networkPoint::getFacesDirichlet() const {
   std::unordered_set<networkFace *> tmp;
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         tmp.emplace(f);
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
      (*)  | / |    *   IN
           ----V--------------  A
      (o)     (o)  (*)   OUT    |
                      (x)
      (x)
   */
   this->map_Face_MirrorPoint.clear();
   const auto x = this->X;
   std::unordered_set<networkFace *> faces;
   B_face.apply(x, mirroring_distance, [&faces](const auto &f) { faces.emplace(f); });
   for (const auto &f : faces) {
      if (f->getNetwork() != this->getNetwork())
         if (Dot(f->incenter - x, f->normal) < 0 /*面と点が向き合っているかどうか*/) {
            if (IntersectQ(x, mirroring_distance, ToX(f)))
               this->map_Face_MirrorPoint[f] = new networkPoint(f->getNetwork(), 2. * (Nearest(x, ToX(f)) - x) + x);
         }
   }
};

/*DOC_EXTRACT networkPoint::contact_angle

| `networkPoint`のメンバー関数/変数      | 説明                                                                |
|-------------------------|--------------------------------------------------------------------------------|
| \ref{contact_angle}{`contact_angle`}         | ２面の法線ベクトルがこの`contact_angle`大きい場合，接触判定から除外される |
| \ref{isFacing}{`isFacing()`}       | ２面の法線ベクトルが`contact_angle`よりも小さいか判定する．ただし，角度は，向かい合う面がなす最小の角度と考える |
| \ref{isInContact}{`isInContact()`}         | 点の隣接面のいずれかが，与えられた面と接触しているか判定する．範囲内で接触しており，かつ`isFacing`が真である場合`true`を返す． |
| \ref{addContactFaces}{`addContactFaces()`}     | バケツに保存された面を基に，節点が接触した面を`networkPoint::ContactFaces`に登録する．   |

現在の実装方法では，接触判定は`networkPoint::addContactFaces`が起点となる．

`networkPoint::addContactFaces`は，節点と隣接する面の組み合わせに対して，接触判定を行い，
`networkPoint::ContactFaces`，`networkPoint::nearestContactFace`，`networkPoint::f_nearestContactFaces`を追加する．

面がノイマン境界条件であるとは，面の全３節点が，`f_nearestContactFaces`に登錄されていることを意味する．
１つでも，p->f_nearestContactFaces[f]が存在しない場合，fはノイマン境界条件でない（また，同時に，pはノイマン節点でないことになる）．

また，節点がノイマン節点であるためには，隣接する全面がノイマン境界条件である必要がある．
そのため，`p->f_nearestContactFaces[隣接面]`が存在しない場合，pはノイマン境界条件でない．

#### 接触の概念図

![接触の概念図](./contact.png)



*/

// \label{contact_angle}
const double contact_angle = 30. * M_PI / 180.;

// \label{isFacing}
bool isFacing(const Tddd &n1, const Tddd &n2) { return isFacing(n1, n2, contact_angle); };
bool isFacing(const Tddd &n2, const T3Tddd &n1) { return isFacing(n1, n2, contact_angle); };
bool isFacing(const T3Tddd &n1, const Tddd &n2) { return isFacing(n1, n2, contact_angle); };
bool isFacing(const T3Tddd &n1, const T3Tddd &n2) { return isFacing(n1, n2, contact_angle); };

// double detection_range(const networkPoint *p) {
//    double mean_d = 0;
//    int count = 0;
//    std::ranges::for_each(p->getNeighbors(), [&](const auto &q) { mean_d += Norm(p->X - q->X); count++; });
//    return mean_d / count / 3;
// };

// \label{isInContact}
bool isInContact(const networkPoint *p, const T3Tddd &f_target) {
   const auto toNearstX = Nearest(p->X, f_target) - p->X;
   // const double angle = 60 * M_PI / 180.;
   //! toNearstX < 1E-5 * std::sqrt(TriangleArea(f_target))この場合，接触点へのベクトルの向きに関わらず，接触していると判定する．
   //! そうでない場合，接触点へのベクトルが，接触面の法線方向を向いているかどうかを判定する．
   // const double neglible_range = 1E-2 * std::sqrt(TriangleArea(f_target));
   // const bool rough_check = (Norm(toNearstX) < neglible_range) || isFlat(TriangleNormal(f_target), toNearstX, angle) || isFlat(TriangleNormal(f_target), -toNearstX, angle);
   // if (!rough_check)
   //    return false;
   //
   if (p->detection_range < Norm(toNearstX))
      return false;  // not in range!

   bool is_close_normal = isFacing(p->getNormal_BEM(), ToX(f_target));
   bool any_close_normal = std::ranges::any_of(p->Faces, [&](const auto &F) { return isFacing(ToX(F), ToX(f_target)); });
   return is_close_normal || any_close_normal;
};

bool isInContact(const networkPoint *p, const Tddd &pX, const T3Tddd &f_target) {
   const auto toNearstX = Nearest(pX, f_target) - pX;
   // const double angle = 60 * M_PI / 180.;
   //! toNearstX < 1E-5 * std::sqrt(TriangleArea(f_target))この場合，接触点へのベクトルの向きに関わらず，接触していると判定する．
   //! そうでない場合，接触点へのベクトルが，接触面の法線方向を向いているかどうかを判定する．
   // const double neglible_range = 1E-2 * std::sqrt(TriangleArea(f_target));
   // const bool rough_check = (Norm(toNearstX) < neglible_range) || isFlat(TriangleNormal(f_target), toNearstX, angle) || isFlat(TriangleNormal(f_target), -toNearstX, angle);
   // if (!rough_check)
   //    return false;
   //
   if (p->detection_range < Norm(toNearstX))
      return false;  // not in range!

   bool is_close_normal = isFacing(p->getNormal_BEM(), ToX(f_target));
   bool any_close_normal = std::ranges::any_of(p->Faces, [&](const auto &F) { return isFacing(F->normal, ToX(f_target)); });
   return is_close_normal || any_close_normal;
};

bool isInContact(const Tddd &X /*base*/, const Tddd &n /*base*/, const T3Tddd &f_target, const networkPoint *p) {
   // nとf_targetの法線が近いかどうかを判定する．次に，最も近い点がrange以内にあるかどうかを判定する．
   if (!isFacing(n, f_target))
      return false;  // not close normal!
   const auto toNearstX = Nearest(X, f_target) - X;
   //! 新たに追加
   // const double neglible_range = 1E-2 * std::sqrt(TriangleArea(f_target));
   // const double angle = 60 * M_PI / 180.;
   // const bool rough_check = !(isFlat(n, toNearstX, angle) || isFlat(n, -toNearstX, angle) || Norm(toNearstX) < neglible_range);
   // if (!rough_check)
   //    return false;
   //
   //!
   if (p->detection_range < Norm(toNearstX))
      return false;  // not in range!
   return true;
};

// is p contact with f in f_normal direction ?
bool isInContact(const networkPoint *p, const networkFace *f_normal, const networkFace *f_target) {
   return isInContact(p->X, f_normal->normal, ToX(f_target), p);
};

// is p contact with f ?
bool isInContact(const networkPoint *p, const networkFace *f_target) {
   return std::ranges::any_of(p->Faces, [&](const auto &f_normal) { return isInContact(p, f_normal, f_target); });
};

bool isInContact(const networkPoint *p, const networkFace *f_normal, const std::unordered_set<networkFace *> &faces_target) {
   return std::ranges::any_of(faces_target, [&](const auto &f_target) { return isInContact(p, f_normal, f_target); });
};

// \label{addContactFaces}
inline void networkPoint::addContactFaces(const Buckets<networkFace *> &B, bool include_self_network = true) {
   /* ----------------------------- for mesh method ---------------------------- */

   /*DOC_EXTRACT networkPoint::addContactFaces()

   | `networkPoint`のメンバー関数/変数      | 説明                                                                |
   |-------------------------|--------------------------------------------------------------------------------|
   | `addContactFaces()`     | バケツに保存された面を基に，節点が接触した面を`networkPoint::ContactFaces`に登録する．   |
   | `std::unordered_set<networkFace *> ContactFaces`          | 節点が接触した面が登録されている．   |
   | `std::tuple<networkFace *, Tddd> nearestContactFace`    | 節点にとって最も近い面とその座標を登録されている．       |
   | `std::unordered_map<networkFace *, std::tuple<networkFace *, Tddd>> f_nearestContactFaces` | この節点に隣接する各面にとって，最も近い面とその座標をこの変数に登録する．           |

   */

   DebugPrint("! まずは，衝突があり得そうな面を多めに保存する．");

   std::vector<std::tuple<networkFace *, double>> f_dist_sort;
   auto add = [&](const auto &f) {
      if (std::ranges::none_of(f_dist_sort, [&](const auto &F) { return f == std::get<0>(F); })) {
         auto d = Norm(this->X - Nearest(this->X, ToX(f)));
         if (this->detection_range >= d)
            f_dist_sort.emplace_back(f, d);
      }
   };

   std::ranges::for_each(this->ContactFaces, [&](const auto &f) { add(f); });

   B.apply(this->X, 2. * this->detection_range /*広め*/, [&](const auto &f) {
      if (include_self_network || f->getNetwork() != this->getNetwork())
         if (isInContact(this, f))
            add(f);
   });

   std::sort(f_dist_sort.begin(), f_dist_sort.end(), [](const auto &a, const auto &b) { return std::get<1>(a) < (std::get<1>(b) - 1E-15); });

   if (!f_dist_sort.empty()) {
      DebugPrint("! ContactFacesに面を登録する．異なる方向を向く面の情報だけが欲しいので，同方向の面は無視する");
      for (const auto &[F, D] : f_dist_sort) {
         if (std::none_of(this->ContactFaces.begin(), this->ContactFaces.end(), [&](const auto &f) { return isFlat(F->normal, -f->normal, M_PI / 180) || isFlat(F->normal, f->normal, M_PI / 180); }))
            this->ContactFaces.emplace(F);
         if (this->ContactFaces.size() > 7)
            break;
      };
   } else
      DebugPrint("faces is empty");

   double distance = 1E+20, tmp;
   Tddd X_near;
   std::get<0>(this->nearestContactFace) = nullptr;
   std::get<1>(this->nearestContactFace) = {1E+20, 1E+20, 1E+20};
   for (const auto &f : this->getContactFaces())
      if (distance > (tmp = Norm(this->X - (X_near = Nearest(this->X, f))))) {
         distance = tmp;
         std::get<1>(this->nearestContactFace) = X_near;
         std::get<0>(this->nearestContactFace) = f;
      }

   auto NearestContactFace_of_f = [&](const networkFace *const f_normal) -> std::tuple<networkFace *, Tddd> {
      Tddd r = {1E+100, 1E+100, 1E+100};
      networkFace *ret = nullptr;
      for (const auto &f_target : bfs(this->getContactFaces(), 2))
         if (isInContact(this, f_normal, f_target)) {
            X_near = Nearest(this->X, ToX(f_target));
            if (Norm(r - this->X) >= Norm(X_near - this->X)) {
               r = X_near;
               ret = f_target;
            }
         }
      return {ret, r};
   };

   for (const auto &f : this->getFaces())
      this->f_nearestContactFaces[f] = NearestContactFace_of_f(f);
};

/* -------------------------------------------------------------------------- */

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
      std::cerr << e.what() << colorReset << std::endl;
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

/* ------------------------------------------------------ */

inline void networkPoint::setX(const Tddd &xyz_IN) {
   try {
      // this->pre_X = xyz_IN;
      CoordinateBounds::setBounds(xyz_IN);
      for (const auto &l : this->getLines()) {
         // std::cout << "l = " << l << std::endl;
         l->setBoundsSingle();
      }
      for (const auto &f : this->Faces) {
         // std::cout << "f = " << f << std::endl;
         // std::cout << f->Points << std::endl;
         // std::cout << f->Lines << std::endl;
         f->setGeometricProperties(ToX(f->getPoints()));
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
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

//@ ------------------------------------------------------ */
//@              BEMでどちらを法線ベクトルとして使うか           */
//@ ------------------------------------------------------ */
inline Tddd networkPoint::getNormalDirichlet_BEM() const { return getNormalDirichletAreaAveraged(); };
inline Tddd networkPoint::getNormalNeumann_BEM() const { return getNormalNeumannAreaAveraged(); };
inline Tddd networkPoint::getNormal_BEM() const {
   if (this->CORNER)
      return Normalize((getNormalDirichlet_BEM() + getNormalNeumann_BEM()) / 2.);
   else {
      return getNormalAreaAveraged();
   }
};

inline V_netPp networkPoint::getNeighborsSort() const {
   try {
      VV_netPp Vps({});
      for (const auto &f : this->Faces) {
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
      std::cerr << e.what() << colorReset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

// コンストラクタ
inline networkPoint::networkPoint(Network *network_IN, const Tddd &xyz_IN, networkLine *xline_IN, networkFace *xface_IN)
    : network(network_IN),
      CoordinateBounds(xyz_IN /*CoordinateBoundsに暗黙に置換される*/),
      initialX(xyz_IN),
      signed_distance_vector({0., 0., 0.}),
      map_Net_ContactPoints({{nullptr, {}}}),
      particlize_info({nullptr, {nullptr, nullptr, nullptr}, {0., 0.}, 0., 0.}) {
   this->network->add(this);
};
// 逆方向の立体角
inline double networkPoint::getSolidAngle() const {
   std::array<double, 3> normal;
   double length = 0., count = 0.;
   for (const auto &f : this->Faces) {
      normal += f->normal;
      length += std::sqrt(2. * f->area);
      count += 1;
   }

   length /= count;
   normal = length * Normalize(normal);
   double ret = 0.;
   for (const auto &f : this->Faces) {
      auto [p0, p1, p2] = f->getPoints(this);
      ret += SolidAngle(p0->X, p1->X, p2->X, normal + p0->X);
   }
   // std::cout << "ret = " << ret << std::endl;
   return 4. * M_PI - ret;
};

inline double networkPoint::getMinimalSolidAngle() const {
   const double TWO_PI = 2. * M_PI;
   auto ratio = getSolidAngle() / TWO_PI;
   return TWO_PI * (ratio >= 1. ? (2. - ratio) : ratio);
};
