#ifndef networkFace_H
#define networkFace_H

#include "Network.hpp"
#include "NetworkUtility.hpp"

inline Tddd networkFace::normalVelocityRigidBody(const Tddd &X) const {
   return this->normal * Dot(this->normal, this->network->velocityRigidBody(X));
};

// コンストラクタ
inline networkFace::networkFace(Network *network_IN, const T_LLL &Lines_IN, T_3P points)
    : Triangle(ToX(points)),
      network(network_IN),
      status(false),
      Lines(Lines_IN),
      Tetras({nullptr, nullptr}),
      XPoints(0) {
#ifdef DEM
   this->contactP = {};
#endif
   try {
      // this->network = storage_IN;  //なぜか初期化リストに入れれない
      // this->network->add(this);
      // this->network = network_IN;  //なぜか初期化リストに入れれない
      this->network->add(this);
      std::get<0>(this->Lines)->Add(this);
      std::get<1>(this->Lines)->Add(this);
      std::get<2>(this->Lines)->Add(this);
      // setBounds();
      // this->Points = this->getPointsFromLines();
      setPoints();
      std::ranges::for_each(this->Points, [](auto &p) { p->setFaces(p->Lines); });
      T3Tddd p0p1p2_X = ToX(this->Points);
      CoordinateBounds::setBounds(p0p1p2_X);
      this->area = TriangleArea(p0p1p2_X);
      this->normal = TriangleNormal(p0p1p2_X);
      this->angles = TriangleAngles(p0p1p2_X);
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
// コンストラクタ
inline networkFace::networkFace(const netFp f)
    : Triangle(extractXtuple(f)),
      network(f->network),
      status(f->status),
      //   Lines(f->Lines),
      Lines(f->Lines),
      XPoints(0),
      /* ------------------------------------------------------ */
      // interp_from_p0({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}),
      // interp_from_p1({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}),
      // interp_from_p2({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}),
      /* ------------------------------------------------------ */
      //   force({0., 0., 0., 0., 0., 0.}),
      //   inertia(6, 0.),
      //   acceleration(3, 0.),
      //   velocity(6, 0.),
      //   mass(0.),
      //   center_of_mass(3, 0.),
      map_Net_ContactPoints({{nullptr, {}}}) {
   // this->grid_pull_factor = f->grid_pull_factor;

   this->Dirichlet = f->Dirichlet;
   this->Neumann = f->Neumann;
   /* ------------------------------------------------------ */
   this->xyzInverse = f->xyzInverse;
   this->normal = f->normal;
   this->angles = f->angles;
   this->area = f->area;
   // this->Lines = f->Lines;
   this->network = f->getNetwork();
   // this->storage = f->getStorage();
   this->network->add(this);
   // this->Points = f->getPoints();
   this->Points = f->getPoints();
   auto [l0, l1, l2] = f->getLines();
   // this->Lines = {l0, l1, l2};
#ifdef DEM
   this->contactP = f->contactP;
#endif
};
// b% -------------------------------------------------------------------------- */
// b% particlizeDetailsは普通のparticlizeに詳しい情報を加えて返す．
// b% 深さ毎に，面の頂点をシフトしてm線形補間に利用する．2021/11/17
// b% -------------------------------------------------------------------------- */
using TPPP = std::tuple<networkPoint *, networkPoint *, networkPoint *>;
inline std::unordered_set<networkPoint *> networkFace::particlize(const double dx, const V_d &depth_list) {
   // depth_list: 法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など

   std::unordered_set<networkPoint *> ret;
   // double alpha;
   // T3Tddd X0X1X2;
   int count = 0;
   networkPoint *p0, *p1, *p2;
   for (const auto &d : depth_list /*double実際の長さ*/) {
      auto [p0_, p1_, p2_] = this->getPoints();
      if (count % 3 == 1) {
         p0 = p2_;
         p1 = p0_;
         p2 = p1_;
      } else if (count % 3 == 2) {
         p0 = p1_;
         p1 = p2_;
         p2 = p0_;
      } else {
         p0 = p0_;
         p1 = p1_;
         p2 = p2_;
      }
      T3Tddd X0X1X2 = {p0->getXtuple() + d / Dot(p0->getNormalTuple(), this->normal) * p0->getNormalTuple(),
                       p1->getXtuple() + d / Dot(p1->getNormalTuple(), this->normal) * p1->getNormalTuple(),
                       p2->getXtuple() + d / Dot(p2->getNormalTuple(), this->normal) * p2->getNormalTuple()};
      for (const auto &[xyz, t0t1] : triangleIntoPoints(X0X1X2, dx)) {
         auto p = new networkPoint(this->getNetwork(), xyz);
         p->particlize_info = {this, {p0, p1, p2}, t0t1, d, dx};
         ret.emplace(p);
         this->addParametricPoints(p);
      }
      count++;
   }
   return ret;
};

bool isOnSameLine(const netPp p, const netPp q) {
   if (p && q)
      if (p->getXLine() == q->getXLine())
         return true;
   return false;
};
int chainAndDelete(V_netPp &base, VV_netPp &PC) {
   V_netPp ret;
   for (size_t i = 0; i < PC.size(); i++) {
      Chain<networkPoint> chn;
      chn.join_front(base, PC[i]);
      if (chn.isComplete /*1 or 2*/) {
         PC.erase(PC.begin() + i);  // 使ったPCを消去して，使えなくする
         base = chn.chained;
         return chn.isComplete;
      }
   }
   return 0;
};
// chainAndDeleteFixingOrder ２番目の引数の反転を許さない
int chainAndDeleteFixingOrder(V_netPp &base, VV_netPp &PC) {
   V_netPp ret;
   for (size_t i = 0; i < PC.size(); i++) {
      Chain<networkPoint> chn;
      chn.join_front_fix_order(base, PC[i]);
      if (chn.isComplete /*1 or 2*/) {
         PC.erase(PC.begin() + i);
         base = chn.chained;
         return chn.isComplete;
      }
   }
   return 0;
};

#endif