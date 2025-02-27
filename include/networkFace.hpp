#pragma once

inline void networkFace::setDodecaPoints() {
   int i = 0;
   for (const auto &p : this->Points) {
      if (i < dodecaPoints.size()) {  // Ensure we do not exceed the array bounds
         dodecaPoints[i++] = std::make_shared<DodecaPoints>(this, p, [](const networkLine *line) -> bool { return !line->CORNER && line->getSurfaces().size() >= 2; });
      }
   }
}

/*
   原点を変えて積分する際に，積分のどこが，原点に依存し，どこが依存しないかを把握しておく．
   原点に依存しない部分は，事前に計算しておくことで，計算量を削減できる．

   また，積分で面をトラバースする際に，係数行列のどの行に（どの点に）重みをかけるかを把握しておく．

   積分では，原点（行）・面・点（列）が重要である．

*/

inline void networkFace::setIntegrationInfo() {
   // set linear integration info
   this->map_Point_BEM_IGIGn_info_init.clear();
   this->map_Point_LinearIntegrationInfo_vector.clear();
   this->map_Point_LinearIntegrationInfo_vector.resize(2);
   this->map_Point_PseudoQuadraticIntegrationInfo_vector.clear();
   this->map_Point_PseudoQuadraticIntegrationInfo_vector.resize(2);
   //@ -------------------------------------------------------------------------- */
   auto addIntegrationInfo = [&](networkPoint *const p) {
      DodecaPoints dodecapoint(this, p, [](const networkLine *line) -> bool { return useOppositeFace(line, M_PI / 3); });
      //% -------------------------------------------------------------------------- */
      std::vector<BEM_IGIGn_info_type> temp;
      for (const auto &[p, f] : dodecapoint.quadpoint.points_faces) temp.push_back({p, f, 0., 0.});
      for (const auto &[p, f] : dodecapoint.quadpoint_l0.points_faces) temp.push_back({p, f, 0., 0.});
      for (const auto &[p, f] : dodecapoint.quadpoint_l1.points_faces) temp.push_back({p, f, 0., 0.});
      for (const auto &[p, f] : dodecapoint.quadpoint_l2.points_faces) temp.push_back({p, f, 0., 0.});
      this->map_Point_BEM_IGIGn_info_init[p] = temp;
      //% -------------------------------------------------------------------------- */
      auto X012 = ToX(this->getPoints(p));
      auto add = [&](const int i, const auto &GWGW) {
         std::vector<linear_triangle_integration_info> info_linears;
         std::vector<pseudo_quadratic_triangle_integration_info> info_quadratics;
         for (const auto &[t0, t1, ww] : GWGW) {
            auto N012_geometry = ModTriShape<3>(t0, t1);
            auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);  // constant value for the linear integration
            auto X = Dot(N012_geometry, X012);
            info_linears.emplace_back(linear_triangle_integration_info{Tdd{t0, t1}, ww * (1. - t0), N012_geometry, X, cross, Norm(cross)});
            //$ ------------------------------------ */
            auto [xi0, xi1, xi2] = N012_geometry;
            X = dodecapoint.X(xi0, xi1);
            cross = dodecapoint.cross(xi0, xi1);
            auto Nc_N0_N1_N2 = dodecapoint.N6_new(xi0, xi1);
            info_quadratics.emplace_back(pseudo_quadratic_triangle_integration_info{Tdd{t0, t1}, ww * (1. - t0), N012_geometry, Nc_N0_N1_N2, X, cross, Norm(cross)});
         }
         this->map_Point_LinearIntegrationInfo_vector[i][p] = info_linears;
         this->map_Point_PseudoQuadraticIntegrationInfo_vector[i][p] = info_quadratics;
      };
      add(0, __array_GW1xGW1__);
      add(1, __array_GW7xGW7__);
   };
   //@ -------------------------------------------------------------------------- */

   auto [p0, p1, p2] = this->getPoints();
   addIntegrationInfo(p0);
   addIntegrationInfo(p1);
   addIntegrationInfo(p2);
};

inline Tddd networkFace::normalVelocityRigidBody(const Tddd &X) const {
   return this->normal * Dot(this->normal, this->network->velocityRigidBody(X));
};

// コンストラクタ
inline networkFace::networkFace(Network *network_IN, networkPoint *p0, networkLine *l0, networkPoint *p1, networkLine *l1, networkPoint *p2, networkLine *l2)
    : Triangle(p0->X, p1->X, p2->X), network(network_IN) {
   try {
      this->network->add(this);
      this->setPoints(p0, l0, p1, l1, p2, l2);
      l0->add(this);
      l1->add(this);
      l2->add(this);
      p0->add(this);
      p1->add(this);
      p2->add(this);
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

// コンストラクタ
inline networkFace::networkFace(Network *network_IN, networkPoint *p0, networkPoint *p1, networkPoint *p2)
    : Triangle(p0->X, p1->X, p2->X), network(network_IN) {
   try {
      this->network->add(this);
      auto l0 = link(p0, p1, network_IN);
      auto l1 = link(p1, p2, network_IN);
      auto l2 = link(p2, p0, network_IN);
      this->setPoints(p0, l0, p1, l1, p2, l2);
      l0->add(this);
      l1->add(this);
      l2->add(this);
      p0->add(this);
      p1->add(this);
      p2->add(this);
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
// コンストラクタ
inline networkFace::networkFace(const netFp f)
    : Triangle(extractXtuple(f)),
      network(f->network),
      status(f->status),
      Lines(f->Lines),
      Points(f->Points),
      PLPLPL(f->PLPLPL) {
   this->network->add(this);
   this->Dirichlet = f->Dirichlet;
   this->Neumann = f->Neumann;
   this->normal = f->normal;
   this->angles = f->angles;
   this->area = f->area;
};
// b% -------------------------------------------------------------------------- */
// b% particlizeDetailsは普通のparticlizeに詳しい情報を加えて返す．
// b% 深さ毎に，面の頂点をシフトしてm線形補間に利用する．2021/11/17
// b% -------------------------------------------------------------------------- */

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