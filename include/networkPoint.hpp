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
      for_each(f->getLinesTuple(), [&](const auto &l) { ret.emplace(l); });
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
   for (const auto &f : this->getContactFaces()) {
      auto intxn = IntersectionSphereTriangle(this->getXtuple(), this->radius, f->getXVertices());
      if (intxn.isIntersecting)
         ret.emplace_back(std::tuple<networkFace *, Tddd>{f, intxn.X});
   }
   return ret;
};
inline std::vector<std::tuple<networkFace *, Tddd>> networkPoint::getContactFacesXCloser() const {
   auto c_faces_sorted = this->getContactFacesX();
   if (c_faces_sorted.empty())
      return c_faces_sorted;
   else {
      auto X = this->getXtuple();
      std::sort(c_faces_sorted.begin(), c_faces_sorted.end(), [X](auto a, auto b) { return isFinite(Norm(std::get<1>(a) - X)) &&
                                                                                           isFinite(Norm(std::get<1>(b) - X)) &&
                                                                                           (Norm(std::get<1>(a) - X) - Norm(std::get<1>(b) - X)) < 1E-10; });
      return c_faces_sorted;
   }
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
inline bool networkPoint::isThereAnyFacingFace(const networkFace *const f_IN, const double rad) const {
   for (const auto &f : this->Faces)
      if (f->isFacing(f_IN, rad))
         return true;
   return false;
};
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
   for (const auto &f : DeleteDuplicates(Flatten(B_face.getObjects(this->getXtuple(), mirroring_distance)))) {
      if (f->getNetwork() != this->getNetwork())
         if (Dot(f->getXtuple() - this->getXtuple(), f->normal) < 0 /*面と点が向き合っているかどうか*/) {
            auto [p0, p1, p2] = f->getPointsTuple();
            auto ITX = intersection(geometry::Sphere(this->getXtuple(), mirroring_distance), geometry::Triangle(p0->getXtuple(), p1->getXtuple(), p2->getXtuple()));
            if (ITX.isIntersecting)
               this->map_Face_MirrorPoint[f] = new networkPoint(f->getNetwork(), 2. * (ITX.X - this->getXtuple()) + this->getXtuple());
         }
   }
};

/* ------------------------------------------------------ */
inline Tddd networkPoint::normalDirichletFace() const {
   //接触面としての法線方向を計算する良い方法が必要．
   Tddd normal = {0., 0., 0.};
   int count = 0;
   for (auto &f : this->Faces) {
      auto [p0, p1, p2] = f->getPointsTuple();
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
      auto [p0, p1, p2] = f->getPointsTuple();
      if ((p0->Neumann || p0->CORNER) && (p1->Neumann || p1->CORNER) && (p2->Neumann || p2->CORNER)) {
         normal += f->getAngle(this) * f->normal;
      }
   }
   if (count == 0)
      error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Dirichlet surface not found");
   return Normalize(normal);
};

inline Tddd networkPoint::normalContanctSurface(const double pw0 = 1., const double pw1 = 1.) const {
   //接触面としての法線方向を計算する良い方法が必要．
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
         len = Norm(f->mirrorPosition(this) - this->getXtuple());
         n += (1. / pow(len, pw1)) * sgn(Dot(f->normal, this->getNormalTuple())) * f->normal;
         totlen += 1. / pow(len, pw1);
      }
   }
   return n / totlen;
};

inline Tddd networkPoint::reflect(const Tddd &v) const {
   // particlizeした点だけが使える．
   auto f = std::get<0>(particlize_info);
   if (f) {
      auto n = f->normal;
      return v - 2. * Dot(v, n) * n;
   } else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
};

inline Tddd networkPoint::X_little_inside() const {
   return this->getXtuple() - (this->radius / 4.) * this->getNormalTuple();
};
/*
@ addContactPointsは，引数として半径が与えられない場合，各点のradiusが成す球として，互いの重なりを調べ，重なった点はContactPointsに保存する．
*/
inline void networkPoint::addContactFaces(const Buckets<networkFace *> &B, bool include_self_network = true) {
   /*
    b% この衝突面の追加方法は，境界要素法用のもので，格子の情報を使って点の法線方向を計算し衝突を判断している．
   */
   if (this->Lines.empty()) {
      // b!まずは，衝突があり得そうな面を多めに保存する．
      std::unordered_map<networkFace *, Tddd> tmpContactFaces;
      // for (const auto &f : B.getObjects_unorderedset(this->getXtuple(), 2. * this->radius))
      for (const auto &f : B.getAll()) {
         auto [p0, p1, p2] = f->getXVertices();
         auto intxn = IntersectionSphereTriangle(this->getXtuple(), 3. * this->radius, f->getXVertices());
         if (Norm(intxn.X - this->getXtuple()) < 2. * this->radius)
            tmpContactFaces[f] = intxn.X;
      }
      // b!各面について衝突があり得るか詳しく調べて判断する．
      for (const auto &[f, intxnX] : tmpContactFaces) {
         //@ 周りの点に隠れて面が見えない状況の場合，その面とは接し得ないので，考慮しない．
         // auto [p0, p1, p2] = f->getPointsTuple();
         // BEMのために導入したもの．周辺の点の方向にある面には衝突し得ないため，無視する．
         // ただし修正面の方向をむくせんを抜き出す必要がある．
         // auto n = f->normal;
         auto thisPointToFace = intxnX - this->getXtuple();
         // b$ 点が面を有していない場合（SPH用）
         if (Dot(f->normal, thisPointToFace /* <--この点からfまでのベクトルを向き合いの判定に利用*/) < 0. /*少なくとも向き合っていなければいけない*/)
            this->ContactFaces.emplace(f);
      }
   } else {
      // b!まずは，衝突があり得そうな面を多めに保存する．
      std::unordered_map<networkFace *, std::tuple<Tddd, Tddd>> tmpContactFaces;
      Tddd r;
      std::unordered_set<networkFace *> faces;
      for (const auto &FF : B.getObjects_unorderedset(this->getXtuple(), 5. * this->radius))
         faces.emplace(FF);

      for (const auto &f : faces) {
         if (include_self_network || !(f->getNetwork() == this->getNetwork())) {
            auto intxn = IntersectionSphereTriangle(this->getXtuple(), this->radius, f->getXVertices());
            if (intxn.isIntersecting) {
               r = intxn.X - this->getXtuple();
               if (Norm(r) < this->radius) {
                  // if (Dot(f->normal, this->getNormalTuple()) <= 0)
                  // if (MyVectorAngle(f->normal, -this->getNormalTuple()) / M_PI * 180. <= 135. ||
                  // 	Norm(r) < this->radius / 50.)
                  {
                     std::get<0>(tmpContactFaces[f]) = intxn.X;
                     std::get<1>(tmpContactFaces[f]) = r;
                  }
               }
            }
         }
      }
      std::vector<std::tuple<networkFace *, double>> faces_sort;
      // b! 各面について衝突があり得るか詳しく調べて判断する．
      for (const auto &[f, intxnX_r] : tmpContactFaces) {
         /*
                  |
         Face |<---thisPointToFace---- thisPoint
                  |
         */
         /*
         _________________________
         |       \     A
         | shadow   \60|
         |      0    \|
         |  q<---------*--------q
         |  |
         |  |
         |  *
         */
         //$ 点が面を有している場合（境界要素法用）
         // auto shadowed = false;
         // auto hitX = std::get<0>(intxnX_r);
         // auto thisPointToFace = std::get<1>(intxnX_r);
         // for (const auto &q : this->getNeighbors())
         // {
         // 	auto thisPointToNeighbor = q->getXtuple() - this->getXtuple();
         // 	auto angle = MyVectorAngle(thisPointToNeighbor, thisPointToFace) * 180. / M_PI;
         // 	if (angle < 10. || Norm(thisPointToFace) > this->radius)
         // 	{
         // 		if (!(Norm(thisPointToFace) <= this->radius / 10.))
         // 		{
         // 			shadowed = true;
         // 			break;
         // 		}
         // 	}
         // }
         // if (!shadowed)
         // if (Dot(f->normal, this->getNormalTuple() /* <--メッシュの場合は法線方向を向き合いの判定に利用*/) <= 0.)
         if (this->isThereAnyFacingFace(f, M_PI / 180. * 20.))
            faces_sort.push_back({f, Norm(std::get<1>(intxnX_r))});
      }
      std::sort(faces_sort.begin(), faces_sort.end(), [](const auto &a, const auto &b) { return Norm(std::get<1>(a)) - Norm(std::get<1>(b)) < 1E-20; });

      std::vector<Tddd> face_normal_vectors;
      for (const auto &f : faces_sort) {
         bool duplication = false;
         for (const auto &n : face_normal_vectors)
            if (isFlat(std::get<0>(f)->normal, n, M_PI / 180))
               duplication = true;
         if (!duplication) {
            this->ContactFaces.emplace(std::get<0>(f));
            face_normal_vectors.emplace_back(std::get<0>(f)->normal);
            if (this->ContactFaces.size() >= 5)
               break;
         }
      };
      /* ------------------------------------------------------ */
      // // b!まずは，衝突があり得そうな面を多めに保存する．
      // std::unordered_map<networkFace *, std::tuple<Tddd, Tddd>> tmpContactFaces;
      // auto X_little_inside = this->X_little_inside(); // this->getXtuple() - (this->radius / 4.) * this->getNormalTuple();
      // for (const auto &f : B.getObjects_unorderedset(X_little_inside, 3. * this->radius))
      // 	if (!(!include_self_network && (f->getNetwork() == this->getNetwork())))
      // 	{
      // 		// auto intxn = intersection(geometry::Sphere(this->getXtuple(), this->radius), geometry::Triangle(f->getXVertices()));

      // 		// 完全に正対している面のみを取得.斜めに見る面は取得しない．
      // 		auto intxn = IntersectionSphereTriangle(X_little_inside, this->radius, f->getXVertices());
      // 		Tddd r = intxn.X - X_little_inside;
      // 		if (intxn.isIntersectingInsideTriangle)
      // 			if (intxn.scale <= 0. /*面の法線方向とは逆向き*/ || Dot(r, f->normal) <= 0.)
      // 			{
      // 				std::get<0>(tmpContactFaces[f]) = intxn.X;
      // 				std::get<1>(tmpContactFaces[f]) = r;
      // 			}
      // 	}
      // //全てを入れる場合
      // for (const auto &[f, intxnX_r] : tmpContactFaces)
      // 	this->ContactFaces.emplace(f);
      /* ------------------------------------------------------ */
      // b!各面について衝突があり得るか詳しく調べて判断する．
      // for (const auto &[f, intxnX_r] : tmpContactFaces)
      // {
      // 	/*
      // 		 |
      // 	Face |<---thisPointToFace---- thisPoint
      // 		 |
      // 	*/
      // 	/*
      // 	_________________________
      // 	|       \     A
      // 	| shadow   \60|
      // 	|      0    \|
      // 	|  q<---------*--------q
      // 	|  |
      // 	|  |
      // 	|  *
      // 	*/
      // 	//$ 点が面を有している場合（境界要素法用）
      // 	auto shadowed = false;
      // 	auto hitX = std::get<0>(intxnX_r);
      // 	auto thisPointToFace = std::get<1>(intxnX_r);
      // 	for (const auto &q : this->getNeighbors())
      // 	{
      // 		auto thisPointToNeighbor = q->getXtuple() - this->getXtuple();
      // 		auto angle = MyVectorAngle(thisPointToNeighbor, thisPointToFace) * 180. / M_PI;
      // 		if (angle < 10. || Norm(thisPointToFace) > this->radius)
      // 		{
      // 			if (!(Norm(thisPointToFace) <= this->radius / 10.))
      // 			{
      // 				shadowed = true;
      // 				break;
      // 			}
      // 		}
      // 	}
      // 	if (!shadowed)
      // 		if (Dot(f->normal, this->getNormalTuple() /* <--メッシュの場合は法線方向を向き合いの判定に利用*/) <= 0. /*少なくとも向き合っていなければいけない*/
      // 			|| Norm(hitX - this->getXtuple()) < this->radius)
      // 			this->ContactFaces.emplace(f);
      // }
   }
   // /*
   //  b% この衝突面の追加方法は，境界要素法用のもので，格子の情報を使って点の法線方向を計算し衝突を判断している．
   // */
   // // b!まずは，衝突があり得そうな面を多めに保存する．
   // std::unordered_map<networkFace *, Tddd> tmpContactFaces;
   // for (const auto &f : B.getObjects_unorderedset(this->getXtuple(), 2. * this->radius))
   // 	if (!(!include_self_network && (f->getNetwork() == this->getNetwork())))
   // 	{
   // 		// auto intxn = intersection(geometry::Sphere(this->getXtuple(), this->radius), geometry::Triangle(f->getXVertices()));
   // 		auto intxn = IntersectionSphereTriangle(this->getXtuple(), 2. * this->radius, f->getXVertices());
   // 		if (intxn.isIntersecting && intxn.scale < 0 /*面の法線方向とは逆向き*/)
   // 			tmpContactFaces[f] = intxn.X;
   // 	}
   // // b!各面について衝突があり得るか詳しく調べて判断する．
   // for (const auto &[f, intxnX] : tmpContactFaces)
   // {
   // 	//@ 周りの点に隠れて面が見えない状況の場合，その面とは接し得ないので，考慮しない．
   // 	// auto [p0, p1, p2] = f->getPointsTuple();
   // 	// BEMのために導入したもの．周辺の点の方向にある面には衝突し得ないため，無視する．
   // 	// ただし修正面の方向をむくせんを抜き出す必要がある．
   // 	// auto n = f->getNormalTuple();
   // 	auto thisPointToFace = intxnX - this->getXtuple();
   // 	if (this->Lines.empty())
   // 	{
   // 		// b$ 点が面を有していない場合（SPH用）
   // 		if (Dot(f->normal, thisPointToFace /* <--この点からfまでのベクトルを向き合いの判定に利用*/) < 0. /*少なくとも向き合っていなければいけない*/)
   // 			this->ContactFaces[f] = intxnX;
   // 	}
   // 	else
   // 	{
   // 		/*
   // 			 |
   // 		Face |<---thisPointToFace---- thisPoint
   // 			 |
   // 		*/
   // 		/*
   // 		_____________________
   // 		|       \  A
   // 		| shadow \ |
   // 		|      60 \|
   // 		|  q<------*--------q
   // 		|  |
   // 		|  |
   // 		|  *
   // 		*/
   // 		//$ 点が面を有している場合（境界要素法用）
   // 		auto shadowed = false;
   // 		for (const auto &q : this->getNeighbors())
   // 		{
   // 			auto thisPointToNeighbor = q->getXtuple() - this->getXtuple();
   // 			auto angle = MyVectorAngle(thisPointToNeighbor, thisPointToFace) * 180. / M_PI;
   // 			if (angle < 60. || Norm(thisPointToFace) > this->radius)
   // 			{
   // 				shadowed = true;
   // 				break;
   // 			}
   // 		}
   // 		if (!shadowed)
   // 			if (Dot(f->normal, this->getNormalTuple() /* <--メッシュの場合は法線方向を向き合いの判定に利用*/) < 0. /*少なくとも向き合っていなければいけない*/)
   // 				this->ContactFaces[f] = intxnX;
   // 	}
   // }
};

//点なのに，sphereでインターセクションをチェックするのは効率的ではない．
inline void networkPoint::addContactPoints(const Buckets<networkPoint *> &B,
                                           const bool include_self_network = true) {
   for (const auto &q : B.getObjects_unorderedset(this->getXtuple(), 2. * this->radius /*depth*/))
      if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
         if (this != q)
            if (intersection(geometry::Sphere(this->getXtuple(), this->radius), geometry::Sphere(this->getXtuple(), q->radius)).isIntersecting)
               this->ContactPoints.emplace(q);
};
// radiusのしていがある場合は，その半径内の点を取得する．
inline void networkPoint::addContactPoints(const Buckets<networkPoint *> &B,
                                           const double radius,
                                           const bool include_self_network = true) {
   for (const auto &q : B.getObjects_unorderedset(this->getXtuple(), 2. * this->radius /*depth*/))
      if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
         if (this != q)
            if (Norm(this->getXtuple() - q->getXtuple()) <= radius)
               this->ContactPoints.emplace(q);
};
inline void networkPoint::addContactPoints(const Buckets<networkPoint *> &B,
                                           const int limit_depth,
                                           const int limit_num,
                                           const bool include_self_network = true) {
   for (const auto &q : B.getObjects_unorderedset(this->getXtuple(), limit_depth, limit_num))
      if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
         if (this != q)
            if (Norm(this->getXtuple() - q->getXtuple()) <= radius)
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

//追加2021/06/18
V_netFp networkPoint::getFacesSort() const {
   V_netFp ret = {(this->getFaces())[0]};
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
//追加2021/09/02
inline V_netFp networkPoint::getFacesSort(networkLine *const line) const {
   try {
      auto p1 = (*line)(this);
      if (!p1)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      networkFace *begin_face;
      //左回りの面を選ぶ．その面からスタート
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
      p = std::get<1>(f->getPointsTuple(this));
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

//!起点を指定できるようにする．
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
      for (const auto &l : this->getLines())
         l->setBoundsSingle();
      for (const auto &f : this->Faces)
         f->setGeometricProperties(ToX(f->setPoints(f->Lines)));
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
   //角度の重みを掛けた法線ベクトルを足し合わせる
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
   //角度の重みを掛けた法線ベクトルを足し合わせる
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
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      normal += f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletArithmeticAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannArithmeticAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
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
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   double radius = 1.5 * Norm(extX(this->getNeighbors()) - this->getXtuple());
   for (const auto &f : this->Faces)
      normal += kernel_Bspline5(Norm(f->getXtuple() - this->getXtuple()), radius) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletSplineKernelAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   double radius = 1.5 * Norm(extX(this->getNeighbors()) - this->getXtuple());
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += kernel_Bspline5(Norm(f->getXtuple() - this->getXtuple()), radius) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannSplineKernelAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   double radius = 1.5 * Norm(extX(this->getNeighbors()) - this->getXtuple());
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += kernel_Bspline5(Norm(f->getXtuple() - this->getXtuple()), radius) * f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                   SubArea Averaged                     */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalSubAreaAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.}, a, b, c;
   for (const auto &f : this->Faces)
      normal += f->getSubArea(this) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletSubAreaAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getSubArea(this) * f->normal;

   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannSubAreaAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
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
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->area * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannAreaAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->area * f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                  Area Averaged Buffer                  */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalAreaAveraged_Buffer() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getAreaBuffer() * f->getNormalBuffer();
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletAreaAveraged_Buffer() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getAreaBuffer() * f->getNormalBuffer();
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannAreaAveraged_Buffer() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->getAreaBuffer() * f->getNormalBuffer();
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*             Inscribed Circle Area Averaged             */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalDirichletInscribedCircleAreaAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getInscribedCircleArea() * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannInscribedCircleAreaAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Neumann)
         normal += f->getInscribedCircleArea() * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalInscribedCircleAreaAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      normal += f->getInscribedCircleArea() * f->normal;
   return Normalize(normal);
};
/* ------------------------------------------------------ */
/*                     Angle Averaged                     */
/* ------------------------------------------------------ */
inline Tddd networkPoint::getNormalDirichlet() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getAngle(this) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumann() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
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
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      normal += f->getAngle(this) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalDirichletAngleAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
   Tddd normal = {0., 0., 0.};
   for (const auto &f : this->Faces)
      if (f->Dirichlet)
         normal += f->getAngle(this) * f->normal;
   return Normalize(normal);
};
inline Tddd networkPoint::getNormalNeumannAngleAveraged() const {
   //角度の重みを掛けた法線ベクトルを足し合わせる
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

// 	return geometry::SolidAngle(this->getX(), obj3D::extractX(this->getNeighborsSort(TorF)));
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

// 	return geometry::SolidAngle(this->getX(), obj3D::extractX(fs));
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
         //これは以下で実現できる．変更2021/06/18
         auto [q0, q1, q2] = f->getPointsTuple(this);
         Vps.push_back(V_netPp{q1, q2});
      }
      auto ret = FlattenAsChain(Vps);
      if (*ret.begin() == *ret.rbegin())
         ret.pop_back();
      else {
         return ret;
         //この場合，
         // * thisが端点であるか
         // * 周りが数珠つなぎとなっていない（エラー）
         // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not chain !");
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
   this->mu_SPH = 0.001005;
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
#endif

   this->Dirichlet = false;
   this->Neumann = false;
};
//逆方向の立体角
inline double networkPoint::getSolidAngle() const {
   return geometry::SolidAngle(this->getXtuple(), extX(this->getNeighborsSort()));
};
inline double networkPoint::getSolidAngleBuffer() const {
   return geometry::SolidAngle(this->getXBuffer(), extXBuffer(this->getNeighborsSort()));
};

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