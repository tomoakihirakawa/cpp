#ifndef Network_H
#define Network_H

#define BEM  // BEMのメンバー変数関数を有効化する

#include <functional>
#include <typeinfo>
#include <unordered_set>
#include "NetworkCommon.hpp"
#include "basic.hpp"
#include "interpolations.hpp"
#include "rootFinding.hpp"

class networkLine;
class networkTetra;
class networkFace;
class networkPoint;
class Network;
class networkPointX;

using netP = networkPoint;
using netF = networkFace;
using netPp = networkPoint *;
using netFp = networkFace *;
using netTp = networkTetra *;
using V_netLp = std::vector<networkLine *>;
using V_netPp = std::vector<networkPoint *>;
using VV_netPp = std::vector<std::vector<networkPoint *>>;
using V_netFp = std::vector<networkFace *>;
using VV_netFp = std::vector<std::vector<networkFace *>>;
using V_Netp = std::vector<Network *>;
//
using T_LL = std::array<networkLine *, 2>;
using T_LLL = std::array<networkLine *, 3>;
using T_3L = std::array<networkLine *, 3>;
using T_LLLL = std::array<networkLine *, 4>;
using T_4L = std::array<networkLine *, 4>;
using T_5L = std::array<networkLine *, 5>;
using T_6L = std::array<networkLine *, 6>;
//
using T_PP = std::array<networkPoint *, 2>;
using T_2P = std::array<networkPoint *, 2>;
using T_PPP = std::array<networkPoint *, 3>;
using T_3P = std::array<networkPoint *, 3>;
using T_4P = std::array<networkPoint *, 4>;
using T_5P = std::array<networkPoint *, 5>;
using T_6P = std::array<networkPoint *, 6>;
//
//
using T4TPPP = std::array<T_PPP, 4>;
using T4T3P = std::array<T_PPP, 4>;
//
using T_FF = std::array<networkFace *, 2>;
using T_FFF = std::array<networkFace *, 3>;
using T_3F = std::array<networkFace *, 3>;
using T_4F = std::array<networkFace *, 4>;
//
using T_TT = std::array<networkTetra *, 2>;
/* ------------------------------------------------------ */

/*extractFaces_detail
`getFaces()`を持つクラスのポインターベクトルを引数にとる．
extractFaces_detail*/
template <class T>
V_netFp extractFaces(const std::vector<T *> &Ls) {
   // V_netFp ret;
   // for (const auto &l : Ls)
   // 	for (const auto &f : l->getFaces())
   // 		ret.emplace_back(f);
   // return DeleteDuplicates(ret);
   std::unordered_set<networkFace *> ret;
   for (const auto &l : Ls)
      for (const auto &f : l->getPoints())
         ret.emplace(f);
   return V_netFp(ret.begin(), ret.end());
};
/*extractPoints_detail
`getPoint()`を持つクラスのポインターベクトルを引数にとる．
extractPoints_detail*/
// template <class T>
// V_netPp extractPoints(const std::vector<T *> &Ls)
// {
// 	V_netPp ret;
// 	for (const auto &l : Ls)
// 		for (const auto &f : l->getPoints())
// 			ret.emplace_back(f);
// 	return DeleteDuplicates(ret);
// };
template <class T>
V_netPp extractPoints(const std::vector<T *> &Ls) {
   std::unordered_set<networkPoint *> ret;
   for (const auto &l : Ls)
      for (const auto &p : l->getPoints())
         ret.emplace(p);
   return V_netPp(ret.begin(), ret.end());
};

template <class T>
V_netLp extractLines(const std::vector<T *> &points) {
   std::unordered_set<networkLine *> ret;
   for (const auto &p : points)
      for (const auto &l : p->getLines())
         ret.emplace(l);
   return V_netLp(ret.begin(), ret.end());
};
/* ------------------------------------------------------ */
/*     *     */
/*    / \    */
/*   *===*   */
/*    \ /    */
/*     *     */
class networkLine : public CoordinateBounds {

  public:
   std::array<double, 3> X_on_curve;
   /* ----------------------------------------------- */
   /*                 MOORING LINE                    */
   /* ----------------------------------------------- */

  public:
   double diameter = 0.;                //! [m]
   double natural_length = 0.;          //! [m]
   double stiffness = 0.;               //! [N/m]
   double damping = 0.;                 //! [N/(m/s^2)]
   double weight_per_unit_length = 0.;  //! [kg/m]

   /* ------------------------------------------------ */
  public:
   bool public_status = false;
#ifdef DEM
  public:
   double tension;
#endif

  public:
   bool CORNER;
   bool Dirichlet;
   bool Neumann;
   bool isIntxn();
   //
   Tddd X_surface;
   double tension_EMT;
   double tension_EMT_;
#if defined(BEM)
   V_d interpoltedX;
#endif
  public:
   //---------------------------------
   Network *network;
   Network *getNetwork() const { return this->network; };

   bool status;

   bool getStatus() const { return this->status; };
   void setStatus(const bool TorF) { this->status = TorF; };
   bool isStatus(const bool TorF) const { return this->status == TorF; };
   //
   V_netPp XPoints;
   V_netPp getXPoints() const { return this->XPoints; };
   void clearXPoints() { return this->XPoints.clear(); };
   bool clearXPoints(const netP *p) { return network::erase(this->XPoints, p); };
   bool addXPoint(netP *p) { return network::add(this->XPoints, p); };
   bool penetrateQ() const { return !XPoints.empty(); };

   V_netFp getFacesPenetrating() const;
   //---------------------------------
  protected:
   networkPoint *Point_A;
   networkPoint *Point_B;
   // なぜこの型なのか？　いつかうunordered_setに変更する．
   V_netFp Faces;
   //---------------------------------
  public:
   networkLine(Network *network_IN, netP *sPoint_IN, netP *ePoint_IN);
   // コピーコンストラクタ
   networkLine(const networkLine *l)
       : CoordinateBounds(l->getBounds()),
         status(false),
         Point_A(nullptr),
         Point_B(nullptr),
         tension_EMT(0.) {
      std::cout << "copying networkLine ...";
      this->network = l->getNetwork();
      this->Faces = l->getFaces();

      this->Dirichlet = l->Dirichlet;
      this->Neumann = l->Neumann;
      this->CORNER = l->CORNER;

      auto tmpL = l->getPoints();
      set(std::get<0>(tmpL), std::get<1>(tmpL));

      setBoundsSingle();
      std::cout << " done" << std::endl;
   };
   ~networkLine();
   //---------------------------------
   // bool setBounds();
   void setBoundsSingle();
   //---------------------------------
   void set(networkPoint *const sPoint_IN, networkPoint *const ePoint_IN) {
      this->Point_A = sPoint_IN;
      this->Point_B = ePoint_IN;
   };
   void set(networkFace *const Face0_IN, networkFace *const Face1_IN) {
      this->Faces = {Face0_IN, Face1_IN};
   };
   void set(networkFace *const Face0_IN) {
      this->Faces = {Face0_IN};
   };
   //---------------------------------
   netPp divide(const Tddd &midX = {1E+80, 1E+80, 1E+80});
   bool canFlip(const double) const;
   bool flip();
   bool flipIfIllegal();
   bool flipIfBetter(const double min_degree_to_flat = M_PI / 180.,
                     const double min_inner_angle = M_PI / 180.,
                     const int min_n = 5);
   bool flipIfTopologicallyBetter(const double min_degree_of_line = M_PI / 180.,
                                  const double min_degree_of_face = M_PI / 180,
                                  const int s_meanIN = 6);
   void divideIfIllegal();
   bool isAdjacentFacesFlat(const double) const;
   bool islegal() const;
   bool isGoodForQuadInterp() const {
      if (this->CORNER)
         return false;
      else
         return true;
   };
   bool isGoodForQuadInterp_Geo() const {
      // 線の中心位置を決めるために，線が２次補間で近似できるか，
      // 周辺の三角形の状況から判断する
      if (this->Neumann && !this->isAdjacentFacesFlat(M_PI / 3.))
         return false;
      if (this->CORNER)
         return false;
      if (!this->islegal())
         return false;
      return true;
   };
   bool isMergeable() const;
   netPp merge();             // deleteしていない方のpointを返す
   netPp mergeIfMergeable();  // deleteしていない方のpointを返す
   //---------------------------------
   //  netP* operator()(netP* a){return this->Point2Point[a];};
   template <class T>
   T *getTheOther(const std::vector<T *> &PorF, const T *const a) const {
      if (PorF.size() < 2)
         return nullptr;
      else if (a == PorF[0])
         return PorF[1];
      else if (a == PorF[1])
         return PorF[0];
      else
         return nullptr;
   };

   netP *operator()(const netP *const a) const {
      if (this->Point_A == a)
         return this->Point_B;
      else if (this->Point_B == a)
         return this->Point_A;
      else
         return nullptr;
   };

   netF *operator()(const netF *const a) const {
      if (this->Faces.size() < 2)
         return nullptr;
      else if (a == this->Faces[0])
         return this->Faces[1];
      else if (a == this->Faces[1])
         return this->Faces[0];
      else
         return nullptr;
   };

   // VV_d getLocations() const { return obj3D::extractX(this->getPoints()); };
   T2Tddd getLocationsTuple() const;

   // V_d getNormal() const;
   Tddd getNormal() const;

   //---------------------------------
   double length() const;
   //------------------
   // V_netPp getPoints() const { return {this->Point_A, this->Point_B}; };
   T_PP getPoints(const networkPoint *const p) const {
      if (p == this->Point_A)
         return {this->Point_A, this->Point_B};
      else
         return {this->Point_B, this->Point_A};
   };
   T_PP getPoints() const { return {this->Point_A, this->Point_B}; };
   T_PP getPointsTuple() const { return {this->Point_A, this->Point_B}; };
   const V_netFp &getFaces() const { return this->Faces; };
   V_netFp getFacesExcept(const netFp f) const { return TakeExcept(this->Faces, f); };
   V_netFp getFacesExcept(const V_netFp &fs) const { return TakeExcept(this->Faces, fs); };
   V_netFp getFaces(const bool TorF) const;
   //------------------
   // 派生クラスのクラス名で，選ばれる関数
   // int Find(netP *p_IN) const { return network::find(this->getPoints(), p_IN); };
   bool replace(const netP *oldP, netP *newP) {
      if (this->Point_A == oldP) {
         this->Point_A = newP;
         return true;
      } else if (this->Point_B == oldP) {
         this->Point_B = newP;
         return true;
      } else
         return false;
   };
   bool Erase(const netP *p_IN) { return replace(p_IN, nullptr); };
   bool Replace(netP *oldP, netP *newP);
   //------------------
   int Find(netF *f_IN) const { return network::find(this->Faces, f_IN); };
   bool Erase(netF *const f_IN) { return network::erase(this->Faces, f_IN); };
   bool Add(netF *f_IN) { return network::add(this->Faces, f_IN); };
   bool replace(netF *oldF, netF *newF) { return network::myswitch(this->Faces, oldF, newF); };
   // bool Replace(netF *oldF, netF *newF, networkLine *newL = nullptr);
   //------------------
   void Delete();
};

using netL = networkLine;
using netLp = networkLine *;
using V_netLp = std::vector<networkLine *>;
using VV_netLp = std::vector<std::vector<networkLine *>>;
/* ------------------------------------------------------ */
// networkPoint限定のもの
// template <>
// struct Buckets<networkPoint *> : public BaseBuckets<networkPoint *> {
//    Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<networkPoint *>(c_bounds, dL_IN){};
//    Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<networkPoint *>(boundingboxIN, dL_IN){};
// };
/* ------------------------------------------------------ */
// template <>
// struct Buckets<networkFace *> : BaseBuckets<networkFace *> {
//    Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<networkFace *>(c_bounds, dL_IN){};
//    Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<networkFace *>(boundingboxIN, dL_IN){};
// };
/* ------------------------------------------------------ */
/*networkPoint_detail
networkPoint_detail*/
/*networkPoint_code*/
/*    \ /    */
/*   --@--   */
/*    / \    */
#include "integrationOfODE.hpp"
class networkPoint : public CoordinateBounds, public CRS {

   /* ------------------------------ MOORING LINE ------------------------------ */
   /*
      例えば，端部には速度が与えられたとしても，内部の節点の速度は，運動方程式に従って，隣の節点との間の張力や重力，抗力そして相対速度に依存する減衰力によって決まる．
      * RK_X_sub
      * RK_velocity_sub
      の二つを使う
   */

  public:
   std::array<double, 3> getGravitationalForce() { return this->mass * _GRAVITY3_; }

   std::array<double, 3> getTension() {
      std::array<double, 3> force, v, relative_velocity;
      force.fill(0);
      double disp, strain;
      for (const auto &l : this->getLines()) {
         v = (*l)(this)->RK_X_sub.getX() - this->RK_X_sub.getX();
         disp = Norm(v) - l->natural_length;
         strain = disp / l->natural_length;
         if (disp > 0.)
            force += l->stiffness * strain * Normalize(v);
         relative_velocity = this->RK_velocity_sub.getX() - (*l)(this)->RK_velocity_sub.getX();
         force -= l->damping * relative_velocity;  // / this->RK_velocity_sub.dt_fixed;
      }
      return force;
   }

   std::array<double, 3> getDragForce(const double Cd = 0.3) {
      std::array<double, 3> drag_force, relative_velocity, mean_v, fluid_velocity, normalized_relative_velocity, A2B;
      drag_force.fill(0);
      fluid_velocity.fill(0);
      double A;
      for (const auto &l : this->getLines()) {
         // mean_v = 0.5 * ((*l)(this)->RK_velocity_sub.getX() + this->RK_velocity_sub.getX());
         mean_v = this->RK_velocity_sub.getX();
         relative_velocity = fluid_velocity - mean_v;
         // A = M_PI * std::pow(l->diameter, 2);
         A2B = (*l)(this)->RK_X_sub.getX() - this->RK_X_sub.getX();
         normalized_relative_velocity = Normalize(relative_velocity);
         A = l->diameter * Norm(A2B) * 0.5;                                                            //! half area
         A *= Norm(normalized_relative_velocity - Dot(normalized_relative_velocity, Normalize(A2B)));  //! projected area
         drag_force += 0.5 * _WATER_DENSITY_ * Dot(relative_velocity, relative_velocity) * Cd * A * normalized_relative_velocity;

         if (!isFinite(A2B) || !isFinite(Normalize(A2B)) || !isFinite(normalized_relative_velocity) || !isFinite(relative_velocity) || !isFinite(drag_force)) {
            std::cout << "A2B = " << A2B << std::endl;
            std::cout << "Normalize(A2B) = " << Normalize(A2B) << std::endl;
            std::cout << "normalized_relative_velocity = " << normalized_relative_velocity << std::endl;
            std::cout << "relative_velocity = " << relative_velocity << std::endl;
            std::cout << "drag_force = " << drag_force << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         }
      }
      return drag_force;
   }

   std::array<double, 3> getForce() {
      return getTension() + getDragForce() + getGravitationalForce();
   }

   /* -------------------------------------------------------------------------- */
  protected:
   bool status;

  public:
   std::array<double, 3> X_last;
   netPp tmpPoint;
   V_netLp Lines;
   std::unordered_set<networkFace *> Faces;
   bool MemberQ(networkFace *const &f_IN) const { return this->Faces.contains(f_IN); };
   Network *network;
   Network *getNetwork() const { return this->network; };
   // Network *getStorage() const { return this->storage; };
   //
   const V_netLp &getLines() const { return this->Lines; };

   V_netLp getLines_toggle(const bool TorF) const {
      V_netLp ret(0);
      for (const auto &line : this->Lines)
         if (line->getStatus() == TorF) {
            line->setStatus(!TorF);
            ret.emplace_back(line);
         }
      return ret;
   };
   int Find(netL *l_IN) { return network::find(this->Lines, l_IN); };
   bool Erase(netL *l_IN) { return network::erase(this->Lines, l_IN); };
   bool add(netL *l_IN) {
      if (!MemberQ_(this->Lines, l_IN)) {
         return this->Lines.emplace_back(l_IN);
         return true;
      } else
         return false;
   };
   bool Add(netL *l_IN) {
      if (!MemberQ_(this->Lines, l_IN)) {
         return this->Lines.emplace_back(l_IN);
         return true;
      } else
         return false;
   };
   bool getStatus() const { return this->status; };
   void setStatus(bool TorF) { this->status = TorF; };
   bool isStatus(bool TorF) const { return this->status == TorF; };

  public:
   // 今のところBEMのためのもの
   RungeKutta<double> RK_phi;
   RungeKutta<Tddd> RK_X, RK_X_sub;
   //
   // 今のところSoft bodyのためのもの
   RungeKutta<Tddd> RK_velocity, RK_velocity_sub;
   RungeKutta<Tddd> RK_force;
   //
   // 今のところSPHのためのもの
   RungeKutta<Tddd> RK_U;
   RungeKutta<double> RK_rho;
   RungeKutta<double> RK_P;
   //
   // なぜか，SPHでルンゲクッタがうまくいかないので（多分dtがDUDtに入ってくるから整合性が取れない）リープフロッグを使ってみる
   LeapFrog<Tddd> LPFG_X;
   LeapFrog<double> LPFG_rho;
   /* -------------------------------------------------------------------------- */
  public:
   V_netLp getLinesCORNER() const {
      V_netLp ret;
      ret.reserve(this->Lines.size());
      for (const auto &l : this->Lines)
         if (l->CORNER)
            ret.emplace_back(l);
      return ret;
   };
   V_netLp getLinesNeumann() const {
      V_netLp ret;
      ret.reserve(this->Lines.size());
      for (const auto &l : this->Lines)
         if (l->Neumann)
            ret.emplace_back(l);
      return ret;
   };
   V_netLp getLinesDirichlet() const {
      V_netLp ret;
      ret.reserve(this->Lines.size());
      for (const auto &l : this->Lines)
         if (l->Dirichlet)
            ret.emplace_back(l);
      return ret;
   };

   bool replace(const netL *oldL, netL *newL) {
      for (auto &l : this->Lines)
         if (l == oldL) {
            l = newL;
            return true;
         }
      return false;
   };

   bool replace(netF *oldF, netF *newF) {
      for (auto &f : this->Faces)
         if (f == oldF) {
            this->Faces.erase(oldF);
            this->Faces.emplace(newF);
            return true;
         }
      return false;
   };

   void sortLinesByLength() {
      std::sort(this->Lines.begin(), this->Lines.end(),
                [](const netLp l0, const netLp l1) {
                   return (l0->length() < l1->length());
                });
   };
   /* ------------------------------------------------------ */
   /*                       For Physics                      */
   /* ------------------------------------------------------ */
   Tddd initialX;  // 必ず設定される
   //
   T6d force = {0., 0., 0., 0., 0., 0.};
   T6d inertia = {0., 0., 0., 0., 0., 0.};
   T6d velocity = {0., 0., 0., 0., 0., 0.};
   T6d acceleration = {0., 0., 0., 0., 0., 0.};
   double mass = 0.;
   Tddd forceTranslational() const { return {force[0], force[1], force[2]}; };
   Tddd forceRotational() const { return {force[3], force[4], force[5]}; };
   Tddd velocityTranslational() const { return {velocity[0], velocity[1], velocity[2]}; };
   Tddd velocityRotational() const { return {velocity[3], velocity[4], velocity[5]}; };
   Tddd accelTranslational() const { return {acceleration[0], acceleration[1], acceleration[2]}; };
   Tddd accelRotational() const { return {acceleration[3], acceleration[4], acceleration[5]}; };
   /* ------------------------------------------------------ */
   T6d &F = this->force;
   T6d &I = this->inertia;
   T6d &A = this->acceleration;
   T6d &V = this->velocity;
   /* ------------------------------------------------------ */
   double density, density_;
   double &rho = this->density;
   double &rho_ = this->density_;
   double volume = 0., volume_ = 0.;
   double radius = 0.;
   double detection_range = 0.;
   double pressure = 0.;
   /* ------------------------------------------------------ */
   Tddd Fxyz() const { return {force[0], force[1], force[2]}; };
   Tddd Txyz() const { return {force[3], force[4], force[5]}; };
   Tddd Vxyz() const { return {velocity[0], velocity[1], velocity[2]}; };
   Tddd Wxyz() const { return {velocity[3], velocity[4], velocity[5]}; };
   /* ------------------------------------------------------ */

   Tddd signed_distance_vector;

   //! SPH
   std::vector<networkPoint *> points_in_SML;

#ifdef DEM
  public:
   V_netPp contactP;
   V_netPp neighborP;  //! for sph
   void setParticle(double volume_IN, double densityIN) {
      this->density = densityIN;
      this->volume = volume_IN;
      this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
      this->mass = this->volume * this->density;
   };

   // 2023/05/16
   //\label{SPH:auxiliaryPoints}
   std::array<networkPoint *, 0> auxiliaryPoints;
   networkPoint *surfacePoint = nullptr;
   networkPoint *auxPoint = nullptr;
   double W;
   /////////////////////////
   // 物性
   double mu_SPH; /*粘性係数：水は約0.001016 Pa.s*/
   //
   double dt_CFL;
   double total_weight, total_weight_, total_weight__, total_N;
   Tddd normal_SPH, v_to_surface_SPH;
   networkFace *mirroring_face;
   double d_empty_center;
   double pn_SPH;
   bool pn_is_set = false;
   bool isSurface = false, isSurface_next = false;
   bool isSurface_tmp = false, isSurface_next_tmp = false;
   bool isSurface_last_tmp = true;
   const int isSurface_count_init = 3;
   int isSurface_count = isSurface_count_init;
   bool isNearSurface = false;
   bool isNeumannSurface = false;
   bool isInsideOfBody = false;
   bool isCaptured = false, isCaptured_ = false;
   bool isChecked = false;
   bool isFluid = false, isFirstWallLayer = false;
   bool isAuxiliary = false;
   bool hasAuxiliary() { return this->isFluid && this->isCaptured && (this->isSurface || this->isSurface_next); };
   bool isFreeFalling = false;
   // double radius_SPH;
   double particle_spacing;

   double SML() const {
      return C_SML * particle_spacing;
   };
   double SML_next() const {
      return C_SML_next * particle_spacing;
   };

   double SML_lap() const {
      return C_SML * particle_spacing;
   };
   double SML_lap_next() const {
      return C_SML_next * particle_spacing;
   };

   double SML_grad() const {
      return 1.2 * C_SML * particle_spacing;
   };
   double SML_grad_next() const {
      return 1.2 * C_SML_next * particle_spacing;
   };

   networkPoint *nearest_wall_p_next = nullptr;
   networkPoint *nearest_wall_p = nullptr;
   double C_SML, C_SML_next;
   double number_density_SPH;
   double density_interpolated_SPH;
   double pressure_SPH, pressure_SPH_, pressure_SPH__;
   int pressure_equation_index = 99;
   double &p_SPH = this->pressure_SPH;
   double &p_SPH_ = this->pressure_SPH_;
   double &p_SPH__ = this->pressure_SPH__;
   double p_SPH_smoothed = 0.;
   double p_EISPH = 0, p_SPH_last = 0;
   double dp_SPH, dp_SPH_;
   double p_SPH_SPP;
   double DPDt_SPH;
   double div_U_error;
   int checked_points_in_radius_of_fluid_SPH, checked_points_in_radius_SPH, checked_points_SPH;
   int checked_points_in_radius_of_fluid_SPH_next, checked_points_in_radius_SPH_next, checked_points_SPH_next;
   double pressure_Tait(const double rho, double C0 = 1466.) const {
      // double C0 = 1466.; //[m/s]
      // C0 /= 5;
      // double C0 = Norm(this->U_SPH); //[m/s]
      // double r = 7.15;
      double r = 7.;
      double rho_w = 1000.;
      double B = rho_w * C0 * C0 / r;
      return B * (std::pow(rho / rho_w, r) - 1.);
   };
   double pressure_Tait() { return pressure_Tait(this->density); };
   double getStaticPressure() {
      double rho_w = 1000.;
      double g = 9.81;
      return -rho_w * g * std::get<2>(getXtuple());
   };
   /////////////////////////
   T3Tddd tensorproduct_grad_Uij, tensorproduct_grad_Uij_next;
   Tddd lap_U;
   T3Tddd grad_U;
   T3Tddd Mat1, Mat2, Mat3, Mat_B;
   std::unordered_map<networkPoint *, double> map_p_grad;
   std::vector<std::tuple<networkPoint *, double>> vector_p_grad;
   // std::unordered_map<networkPoint *, double> extra_storage;
   // std::vector<std::tuple<networkPoint *, double>> extra_storage_vector;
   // void setLapU(const V_d &lap_U_IN) { this->lap_U = lap_U_IN; };
   void setLapU(const Tddd &lap_U_IN) { this->lap_U = lap_U_IN; };
   void setDensityVolume(const double &den, const double &v) {
      // 質量は保存
      this->density = den;
      this->volume = v;
      this->mass = v * den;
      this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
   };
   void setDensity(const double &den) {
      // 質量は保存
      this->density = den;
      this->volume = this->mass / this->density;
      this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
   };
   void setDensity_ConstantVolume(const double &den) {
      // 質量は保存
      this->density = den;
      this->mass = this->volume * this->density;
      this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
   };
   void setVolume(const double &v) {
      // 質量は保存
      this->volume = v;
      this->density = this->mass / this->volume;
      this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
   };
   double div_U, div_tmpU, div_tmpU_, PoissonRHS, div_U_next, DrhoDt_SPH_next;
   Tddd grad_div_U, grad_div_U_;
   Tddd gradP_SPH, gradP_SPH_;
   std::unordered_map<networkPoint *, Tddd> grad_coeff;
   std::unordered_map<networkPoint *, Tddd> grad_coeff_next;
   //////////////////////////
   netFp face_org;
   double a_viscosity;
   Tddd viscosity_term;  // nu*laplacian(U)
   Tddd U_SPH, U_SPH_, marker_X, marker_U;
   T3Tddd nabla_otimes_U;
   Tddd U_XSPH, U_XSPH_next, U_next;
   InterpolationLagrange<std::array<double, 3>> *interp_U_lag = nullptr;
   InterpolationBspline<std::array<double, 3>> *interp_U_Bspline = nullptr;
   std::vector<double> vec_time_SPH;
   std::vector<Tddd> vec_U_SPH;
   double tmp_density;
   Tddd pre_U_SPH;
   Tddd mu_lap_rho_g_SPH;
   Tddd interp_normal, interp_normal_original, interp_normal_original_modified, interp_normal_original_choped, interp_normal_original_next_choped;
   Tddd interp_normal_water, interp_normal_water_next;
   Tddd intp_normal_Eigen;
   Tddd interp_normal_rigid, interp_normal_rigid_next;
   Tddd X_next;
   double volume_next, mass_next;
   Tddd COM_SPH, vec2COM, vec2COM_next;
   double intp_density, intp_density_next, ddr_intp_density;
   double totalMass_SPH;
   Tddd interp_normal_next, interp_normal_original_next;
   std::vector<Tddd> vector_to_polygon;
   std::vector<Tddd> vector_to_polygon_next;
   Tddd cg_neighboring_particles_SPH;
   Tddd b_vector;
   //
   bool is_grad_corr_M_singular = false;
   bool is_grad_corr_M_next_singular = false;
   std::array<std::array<double, 3>, 3> grad_corr_M = _I3_, grad_corr_M_next = _I3_;
   std::array<std::array<double, 3>, 3> laplacian_corr_M = _I3_, laplacian_corr_M_next = _I3_;
   std::array<std::array<double, 3>, 3> inv_grad_corr_M = _I3_, inv_grad_corr_M_next = _I3_;
   std::array<std::array<double, 3>, 3> Eigenvectors_of_M = _I3_, Eigenvectors_of_M1 = _I3_, Eigenvectors_of_M_next = _I3_, Eigenvectors_of_M1_next = _I3_;
   std::array<double, 3> Eigenvalues_of_M, Eigenvalues_of_M1, Eigenvalues_of_M_next, Eigenvalues_of_M1_next;
   //
   double var_Eigenvalues_of_M = 0.;
   double min_Eigenvalues_of_M = 0.;
   double max_Eigenvalues_of_M = 0.;
   double var_Eigenvalues_of_M1 = 0.;
   double min_Eigenvalues_of_M1 = 0.;
   double max_Eigenvalues_of_M1 = 0.;
   double var_Eigenvalues_of_M_next = 0.;
   double var_Eigenvalues_of_M1_next = 0.;
   // ダミー粒子としての情報

   std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> v_reeDW, v_reeDW_next;
   std::array<std::array<std::array<double, 3>, 3>, 3> v_eeDW, v_eeDW_next;
   std::array<std::array<std::array<double, 3>, 3>, 3> v_rrDW, v_rrDW_next;

   /* ------------------- 多段の時間発展スキームのため ------------------- */
   Tddd DUDt_SPH = {0., 0., 0.};
   Tddd DUDt_update_SPH = {0., 0., 0.};
   Tddd DUDt_modify_SPH = {0., 0., 0.};
   Tddd DUDt_modify_SPH_2 = {0., 0., 0.};
   Tddd ViscousAndGravityForce_, tmp_ViscousAndGravityForce_;
   Tddd ViscousAndGravityForce, tmp_ViscousAndGravityForce;
   Tddd repulsive_force_SPH;
   double DrhoDt_SPH;
#endif
   // std::tuple<networkFace * /*補間に使った三角形の頂点*/,
   // 		   double /*パラメタt0*/,
   // 		   double /*パラメタt1*/,
   // 		   double /*深さ方向距離*/,
   // 		   double /*粒子間隔*/>
   // 	particlize_info;
   std::tuple<networkFace * /*補間に使った三角形の頂点*/,
              T_PPP, /*補間に使った三角形の頂点*/
              Tdd /*パラメタt0,t1*/,
              double /*深さ方向距離*/,
              double /*粒子間隔*/>
       particlize_info;
   // Tddd reflect(const Tddd &v) const;
   /* ------------------------ 境界条件 ------------------------ */
  public:
   bool isMultipleNode = false;
   bool CORNER = false;
   bool Dirichlet = false;
   bool Neumann = false;

   /* -------------------------------------------------------------------------- */
   // std::unordered_map<std::tuple<networkFace *, bool>, int> key2Index;
   // int Index(const networkFace *f) { return key2Index.at(variableID_BEM(f)); };
   // std::tuple<networkFace *, bool> variableID_BEM(const networkFace *f) const;
   // std::unordered_set<std::tuple<networkFace *, bool>> variableIDs_BEM() const;

   /* -------------------------------------------------------------------------- */

   std::map<Network *, int> net_depth;
   int minDepthFromCORNER, minDepthFromCORNER_;              // remeshのために導入
   int minDepthFromMultipleNode, minDepthFromMultipleNode_;  // remeshのために導入
   //
   bool isCorner() const { return this->CORNER; };
   bool isDirichlet() const { return this->Dirichlet; };
   bool isNeumann() const { return this->Neumann; };
   void setC() {
      this->Dirichlet = false;
      this->Neumann = false;
      this->CORNER = true;
   };
   void unsetC() {
      this->CORNER = false;
   };
   //
   void setDirichlet() {
      this->Dirichlet = true;
      this->Neumann = false;
      this->CORNER = false;
   };
   void unsetDirichlet() { this->Dirichlet = false; };
   //
   void setNeumann() {
      this->Dirichlet = false;
      this->Neumann = true;
      this->CORNER = false;
   };
   void unsetNeumann() { this->Neumann = false; };
   void setD() { this->setDirichlet(); };
   void unsetD() { this->unsetDirichlet(); };
   void setN() { this->setNeumann(); };
   void unsetN() { this->unsetNeumann(); };
   //
   // bool isCompleteNeumann()
   // {
   // 	//近傍が全てノイマンであればノイマンです．
   // 	//導入理由は？
   // 	for (const auto &p : getNeighbors())
   // 		if (!p->Neumann)
   // 			return false;
   // 	return true;
   // };

#ifdef BEM
  public:
   std::array<double, 2> igign, igign_FMM;  // spherical harmonics のチェックに利用2024/07/03
   std::array<double, 2> igign_near, igign_far;
   // double phi_n();
   using T_PBF = std::tuple<netP *, bool, netF *>;
   std::unordered_map<T_PBF, std::unordered_map<T_PBF, Tdd>> IGIGn;
   Tdd phiphin;
   Tdd phiphin_t;
   Tddd U_absorbed = {0., 0., 0.};
   Network *absorbedBy = nullptr;
   double phi_tmp = 0;
   double phin_tmp = 0;
   // T2T6d phiphin_t_a6;  // 加速度による微分
   T2T6d phiphin_t_a6;  // 加速度による微分
   double phi_Dirichlet;
   double phin_Dirichlet;
   double dpda;

   std::map<networkFace *, int> face2id;

   std::map<networkFace *, double> phiOnFace, phinOnFace;
   std::map<networkFace *, double> phitOnFace, phintOnFace;
   // std::unordered_map<networkFace *, std::tuple<double, T6d>> phintOnFace;
   // std::unordered_map<networkFace *, double> phintOnFace_a;  // 加速度による微分
   //* ------------------------------------------------- */
   // Tddd X_BUFFER;
   // const Tddd &getXBuffer() const { return X_BUFFER; };
   Tddd vecToSurface;
   double shift_coeff = 0;
   Tddd vecToSurface_BUFFER, vecToSurface_BUFFER_BUFFER;
   Tddd whereToGo = {0., 0., 0.};
   Tddd whereToGo_to_check = {0., 0., 0.};
   // const Tddd &getUBuffer() const { return vecToSurface; };
   //
   // Tdd phiphin_BUFFER;
   // Tdd phiphin_t_BUFFER;
   // double phi_Neumann_BUFFER;
   // double phi_Dirichlet_BUFFER;
   // double phin_Neumann_BUFFER;
   // double phin_Dirichlet_BUFFER;
   //* ------------------------------------------------- */
   std::unordered_map<networkFace *, Tdd> multiple_phiphin;
   V_d nabla_phi() const;

   Tddd grad_phi_BEM;
   Tddd U_BEM, U_BEM_last;
   Tddd grid_tension, grid_tension_, grid_tension__;
   /* ------------------------------------------------------ */
   // clungSurfaceを計算してbufferPotentialsOnClungSurfaceを実行
   //
   //<clingすべき座標，補間パラメタt0，パラメタt1，補間する線，または面>
   // std::tuple<Tddd, double, double, networkLine *, networkFace *> clungSurface;
   Tddd clungSurface;
   // void calculateBufferPotentialsOnClungSurface();
   // void copyPotentialsBuffer()
   // {
   // 	this->phiphin = this->phiphin_BUFFER;
   // 	this->phiphin_t = this->phiphin_t_BUFFER;
   // 	this->phi_Neumann = this->phi_Neumann_BUFFER;
   // 	this->phi_Dirichlet = this->phi_Dirichlet_BUFFER;
   // 	this->phin_Neumann = this->phin_Neumann_BUFFER;
   // 	this->phin_Dirichlet = this->phin_Dirichlet_BUFFER;
   // };
   /* ------------------------------------------------------ */
   Tddd U_update_BEM;  // EMTによってシフトさせた流速，法線方向成分はU_BEMと一致させる必要がある．
   Tddd U_cling_to_Neumann, U_cling_to_Neumann_;
   // Tddd U_mod_BEM;
   Tddd U_tangential_BEM, U_tangential_BEM_last;
   Tddd U_normal_BEM;
   Tddd laplacian_U_BEM;
   Tddd normal_BEM;
   Tddd normal_Dir_BEM;
   Tddd normal_Neu_BEM;
   double pressure_BEM;
   double kappa_BEM;
   // double DphiDt;
   //
   // Tddd normalDirichletFace() const;
   // Tddd normalNeumannFace() const;
   // Tddd normalContanctSurface(const double pw0, const double pw1) const;

   double height(const Tddd &offset = {0., 0., 0.}, const Tddd &g_center = {0., 0., -1E+20}) const {
      // 重力中心から，この点までの方向を高さ方向として，offsetから測った高さを返す
      // return Dot(this->X - offset, Normalize(this->X - g_center));

      if (offset[0] == 0. && offset[1] == 0. && offset[2] == 0. && g_center[0] == 0. && g_center[1] == 0. && g_center[2] == -1E+20)
         return std::get<2>(this->X);
      else
         return Dot(this->X - offset, Normalize(this->X - g_center));
   };

   double aphiat(const double pressure /*zero if on atmosfere*/,
                 const Tddd &offset = {0., 0., 0.},
                 const Tddd &g_center = {0., 0., -1E+20}) const {
      // double g = 9.81;	// [m/s2]
      // double rho = 1000.; // [kg/m3]
      if (offset[0] == 0. && offset[1] == 0. && offset[2] == 0. && g_center[0] == 0. && g_center[1] == 0. && g_center[2] == -1E+20)
         return -0.5 * Dot(this->U_BEM, this->U_BEM) - _GRAVITY_ * std::get<2>(this->X) - pressure / _WATER_DENSITY_;
      else
         return -0.5 * Dot(this->U_BEM, this->U_BEM) - _GRAVITY_ * height(offset, g_center) - pressure / _WATER_DENSITY_;
   };

   double DphiDt(const Tddd &U_modified,
                 const double pressure /*zero if on atmosfere*/,
                 const Tddd &offset = {0., 0., 0.},
                 const Tddd &g_center = {0., 0., -1E+20}) const {
      // [1] C. Wang, B. C. Khoo, and K. S. Yeo, “Elastic mesh technique for 3D BIM simulation with an application to underwater explosion bubble dynamics,” Comput. Fluids, vol. 32, no. 9, pp. 1195–1212, Nov. 2003. equaiton (11)
      // return this->DphiDt(pressure, offset, g_center) + Dot(U_modified - this->U_BEM, this->U_BEM);
      return aphiat(pressure, offset, g_center) + Dot(U_modified, this->U_BEM);  // 以下はと同じ：
      // return Dot(U_modified - 0.5 * this->U_BEM, this->U_BEM) - _GRAVITY_ * height(offset, g_center) - pressure / _WATER_DENSITY_;
      //
   };

   double DphiDt_damped(const std::array<double, 2> gamma_phi,
                        const Tddd &U_modified,
                        const double pressure /*zero if on atmosfere*/,
                        const Tddd &offset = {0., 0., 0.},
                        const Tddd &g_center = {0., 0., -1E+20}) const {
      // [1] C. Wang, B. C. Khoo, and K. S. Yeo, “Elastic mesh technique for 3D BIM simulation with an application to underwater explosion bubble dynamics,” Comput. Fluids, vol. 32, no. 9, pp. 1195–1212, Nov. 2003. equaiton (11)
      // return this->DphiDt(pressure, offset, g_center) + Dot(U_modified - this->U_BEM, this->U_BEM);
      auto [gamma, reference_phi] = gamma_phi;
      return -gamma * (std::get<0>(this->phiphin) - reference_phi) + DphiDt(U_modified, pressure, offset, g_center);  // 以下はと同じ：
      // return Dot(U_modified - 0.5 * this->U_BEM, this->U_BEM) - _GRAVITY_ * height(offset, g_center) - pressure / _WATER_DENSITY_;
      //
   };

   double DphiDt(const double pressure /*zero if on atmosfere*/,
                 const Tddd &offset = {0., 0., 0.},
                 const Tddd &g_center = {0., 0., -1E+20}) const {
      return DphiDt(this->U_BEM, pressure, offset, g_center);
      // // offsetは平均の水面位置
      // double gamma = 72.75 * 1E-3; //[N/m] 水20度
      // double g = 9.81;			 // [m/s2]
      // double rho = 1000.;			 // [kg/m3]
      // // double dphidt = Dot(this->U_BEM, this->U_BEM) / 2. - g * height(offset, g_center) - pressure / rho;
      // // dphidt -= gamma / rho * this->kappa_BEM;
      // return aphiat(pressure, offset, g_center) + Dot(this->U_BEM, this->U_BEM);
   };
      // double phi_t(const double pressure = 0. /*zero if on atmosfere*/,
      // 			 const Tddd &offset = {0., 0., 0.},
      // 			 const Tddd &g_center = {0., 0., -1E+20}) const
      // {
      // 	return this->DphiDt(pressure, offset, g_center) - Dot(this->U_BEM, this->U_BEM);
      // };

      // double pressure(const double DphiDt_IN,
      // 				const Tddd &offset = {0., 0., 0.},
      // 				const Tddd &g_center = {0., 0., -1E+20}) const
      // {
      // 	// offsetは平均の水面位置
      // 	// double gravity = 9.80665; // [m/s2]
      // 	// double density = 997.;	  // [kg/m3]
      // 	double tmp = Dot(this->U_BEM, this->U_BEM) / 2. - _GRAVITY_ * height(offset, g_center) - DphiDt_IN;
      // 	return _WATER_DENSITY_ * tmp;
      // };

      // double pressure_EMT(const Tddd &U_modified,
      //                     const double DphiDt_IN,
      //                     const Tddd &offset = {0., 0., 0.},
      //                     const Tddd &g_center = {0., 0., -1E+20}) const {
      //    // offsetは平均の水面位置
      //    // double gravity = 9.80665; // [m/s2]
      //    // double density = 997.;	  // [kg/m3]
      //    double tmp = Dot(this->U_BEM, this->U_BEM) / 2. - _GRAVITY_ * height(offset, g_center) - DphiDt_IN + Dot(this->U_BEM, U_modified);
      //    return _WATER_DENSITY_ * tmp;
      // };
#endif

   // メッシュの衝突を対し噛める際に使用2022/03/10
   // bool isThereAnyFacingFace(const networkFace *const f, const double rad = 1E-10) const;

  private:
   networkLine *xline;
   networkFace *xface;
   // b% ------------------------------------------------------ */
   // b%       面に対する鏡像位置の粒子．一つの面に対して一点決まる　        */
   // b% ------------------------------------------------------ */
   std::unordered_map<networkFace *, networkPoint *> map_Face_MirrorPoint;

  public:
   void makeMirroredPoints(const Buckets<networkFace *> &B_face, const double mirroring_distance);
   void clearMirroredPoints() {
      for (const auto &[f, p] : this->map_Face_MirrorPoint)
         delete p;
      this->map_Face_MirrorPoint.clear();
   };
   //% ------------------------------------------------------ */
   //%                      接触の判別用　                      */
   //% ------------------------------------------------------ */

  private:
   std::tuple<networkFace *, Tddd> nearestContactFace;
   std::unordered_map<networkFace *, std::tuple<networkFace *, Tddd>> f_nearestContactFaces;
   std::unordered_set<networkFace *> ContactFaces;

   std::unordered_set<networkPoint *> ContactPoints;
   std::unordered_map<Network *, std::unordered_set<networkPoint *>> map_Net_ContactPoints;

  public:
   //% ------------------------------------------------------ */
   //%                          接触の判別                      */
   //% ------------------------------------------------------ */
   /*DOC_EXTRACT networkPoint::getContactFaces()

   * `getContactFaces()`で`ContactFaces`呼び出せる．
   * `getNearestContactFace()`で`nearestContactFace`呼び出せる．
   * `getNearestContactFace(face)`で`f_nearestContactFaces`呼び出せる．

   */

   const std::unordered_set<networkFace *> &getContactFaces() const { return this->ContactFaces; };
   const std::tuple<networkFace *, Tddd> &getNearestContactFace() const { return this->nearestContactFace; };

   const std::unordered_map<networkFace *, std::tuple<networkFace *, Tddd>> &getNearestContactFaces() const { return this->f_nearestContactFaces; };

   const std::tuple<networkFace *, Tddd> getNearestContactFace_(const networkFace *const f) const {
      auto it = this->f_nearestContactFaces.find(const_cast<networkFace *>(f));
      return (it != this->f_nearestContactFaces.end()) ? it->second : std::tuple<networkFace *, Tddd>{nullptr, {0., 0., 0.}};
   };
   networkFace *getNearestContactFace(const networkFace *const f) const { return std::get<0>(getNearestContactFace_(f)); };

   void clearContactFaces() {
      this->nearestContactFace = {nullptr, {0., 0., 0.}};
      this->f_nearestContactFaces.clear();
      this->ContactFaces.clear();
   };

   /* ----------------------------------------------------------- */

   std::vector<std::tuple<networkFace *, Tddd>> getContactFacesX() const;
   std::vector<std::tuple<networkFace *, Tddd>> getContactFacesXCloser() const;
   void addContactFaces(const Buckets<networkFace *> &B, bool);  // 自身と同じfaceを含まない
   // void addContactPoints(const Buckets<networkPoint *> &B, const bool);                                    // 自身と同じnetのpointを含み得る
   // void addContactPoints(const Buckets<networkPoint *> &B, const double radius, const bool);               // 自身と同じnetのpointを含み得る
   // void addContactPoints(const std::vector<Buckets<networkPoint *>> &B, const double radius, const bool);  // 自身と同じnetのpointを含み得る
   // Tddd X_little_inside() const;
   //
   const std::unordered_set<networkPoint *> &getContactPoints() const { return this->ContactPoints; };
   const std::unordered_set<networkPoint *> &getContactPoints(Network *const net) const {
      auto it = this->map_Net_ContactPoints.find(net);
      if (it != this->map_Net_ContactPoints.end())
         return it->second;
      else
         return this->map_Net_ContactPoints.at(nullptr);
   };
   std::unordered_set<networkPoint *> getContactPoints(const std::vector<Network *> &net) const {
      std::unordered_set<networkPoint *> ret;
      for (const auto &n : net) {
         auto tmp = this->getContactPoints(n);
         ret.insert(tmp.begin(), tmp.end());
      }
      return ret;
   };
   void clearContactPoints() {
      this->ContactPoints.clear();
      this->map_Net_ContactPoints.clear();
      this->map_Net_ContactPoints[nullptr] = {};
   };
   // void addContactPoints(networkPoint *const p) {
   //    this->ContactPoints.emplace(p);
   //    this->map_Net_ContactPoints[p->getNetwork()].emplace(p);
   // };
   // void addContactPoints(const std::vector<networkPoint *> &P) {
   //    for (const auto &p : P)
   //       this->addContactPoints(p);
   // };
   // void addContactPoints(const std::vector<std::vector<networkPoint *>> &VVP) {
   //    for (const auto &VP : VVP)
   //       this->addContactPoints(VP);
   // };
   // void addContactPoints(const std::unordered_set<networkPoint *> &P) {
   //    for (const auto &p : P) {
   //       this->ContactPoints.emplace(p);
   //       this->map_Net_ContactPoints[p->getNetwork()].emplace(p);
   //    }
   // };
   // インライン化する．
   // コンタクトポイントをpointとface上で作成，完成させる．
   // 次にRBFを使った力の計算をその点を使って行えるようにする．
   //--------------------
   // コンストラクタ
   // networkPoint(Network *network_IN, const V_d &xyz_IN, networkLine *line_IN = nullptr, networkFace *face_IN = nullptr);
   networkPoint(Network *network_IN,
                const Tddd &xyz_IN,
                networkLine *xline_IN = nullptr,
                networkFace *xface_IN = nullptr);
   ~networkPoint();
   /* -------------------------------------------------------------------------- */
   NewtonRaphson<double> NR_double;
   /* -------------------------------------------------------------------------- */
   V_d getFaceAreas() const;
   std::unordered_set<networkLine *> getLinesAround() const;
   std::unordered_set<networkLine *> getLinesOppsoite() const;
   //---------------------
   networkLine *getXLine() const { return this->xline; };
   networkFace *getXFace() const { return this->xface; };
   bool isXPoint() const { return (xline != nullptr && xface != nullptr); };
   bool intersectQ() const { return this->isXPoint(); };
   /*for a cross point*/
   //--------------
   // void set(const V_d& xyz_IN){object3D::setBounds(xyz_IN);};
   void setX(const V_d &xyz_IN);
   void setX(const Tddd &xyz_IN);
   void setXSingle(const Tddd &xyz_IN) {
      // this->pre_X = this->X;
      CoordinateBounds::setBounds(xyz_IN);
   };
   // bool setXcarefully(const V_d &xyz_IN);
   void resetXinfo();
   void setLinesStatus(const bool TorF) {
      for (const auto &l : Lines)
         l->setStatus(TorF);
   };
   // void setBounds();
   void setBoundsSingle() { CoordinateBounds::setBounds(this->X); };
   V_netPp getXNeighbors() const;
   V_netPp getNeighbors() const {
      V_netPp ret(this->Lines.size());
      int i = 0;
      for (const auto &l : this->Lines)
         ret[i++] = (*l)(this);
      return ret;
   };

   int Find(const netP *lookforobj) const {
      for (int i = 0; i < this->Lines.size(); i++)
         if ((*(this->Lines[i]))(this) == lookforobj)
            return i;
      return -1;
   };

   netL *getLineBetween(const netP *const p) const {
      for (const auto &l : this->Lines)
         if (p == (*l)(this))
            return l;
      return nullptr;
   };

   bool isOpen() const {
      for (const auto &l : this->Lines)
         if (l->getFaces().size() < 2)
            return true;
      return false;
   };
   bool isClosed() const {
      for (const auto &l : this->Lines)
         if (l->getFaces().size() < 2)
            return false;
      return true;
   };

   V_netPp getNeighborsSort() const;
   V_netPp getNeighborsSort2() const;
   V_netPp getNeighborsSort(bool TorF);
   //
   // using V_Tup_ddPp = std::vector<std::tuple<double, double, networkPoint *>>;
   // using V_Tup_ddVdPp = std::vector<std::tuple<double, double, V_d, networkPoint *>>;
   // V_Tup_ddVdPp getNeighbors_Depth1_OnPolarAsTuple(networkLine *base_line) const;
   // V_Tup_ddVdPp getNeighbors_Depth2_OnPolarAsTuple(networkLine *base_line) const;
   // V_Tup_ddVdPp getNeighbors_Depth2_OnPolarAsTuple_parametric(networkLine *base_line) const;
   // using V_Tup_ddVd = std::vector<std::tuple<double, double, V_d>>;
   // V_Tup_ddVd getNeighbors_Depth2_OnPolarAsTuple2() const;
   // // V_Tup_dddPp getNeighborsOnPolar() const;

   // using V_Var_VVdVPp = std::vector<std::variant<VV_d, V_netPp>>;
   // V_Var_VVdVPp getNeighbors_Depth1_OnPolarAsVariant() const;
   // V_Var_VVdVPp getNeighbors_Depth2_OnPolarAsVariant() const;
   /* ------------------------------------------------------ */
   // V_netFp getFaces_intersectQ(const bool TorF) const;
   /* ------------------------------------------------------ */
   V_netFp getFaces(networkLine *line) const;
   /* ------------------------------------------------------ */
   V_netFp getFaces(const bool TorF) const {
      V_netFp ret;
      ret.reserve(this->Lines.size());
      for (const auto &l : this->Lines)
         for (const auto &f : l->getFaces(TorF))
            if (std::find(ret.begin(), ret.end(), f) == ret.end())
               ret.emplace_back(f);
      return ret;
   };

   V_d getAngles() const;
   V_d getAngles(networkLine *const base_line) const;

   Tddd normal_to_be_preserved = {0., 0., 0.};

#define use_binary_search
   /**
   2021/09/02unordered_setを使うよう修正した
   将来的にはunordered setを返す関数に修正すべき
    **/

   const std::unordered_set<networkFace *> &getFaces() const { return this->Faces; };

   std::unordered_set<networkFace *> setFaces(const V_netLp &this_lines) {
      std::unordered_set<networkFace *> ret;
      for (const auto &l : this_lines)
         if (l)
            for (const auto &f : l->getFaces())
               ret.emplace(f);
      return (this->Faces = ret);
   };

   std::unordered_set<networkFace *> setFaces() {
      this->Faces.clear();
      for (const auto &l : this->Lines)
         for (const auto &f : l->getFaces())
            this->Faces.emplace(f);
      return this->Faces;
   };

   std::unordered_set<networkFace *> getFacesFromLines() const {
      std::unordered_set<networkFace *> ret;
      for (const auto &l : this->Lines)
         for (const auto &f : l->getFaces())
            ret.emplace(f);
      return ret;
   };

   bool erase(networkFace *f) {
      auto it = this->Faces.find(f);
      if (it != this->Faces.end()) {
         this->Faces.erase(it);
         return true;
      } else
         return false;
   };

   std::unordered_set<networkFace *> getFaces0() const {
      std::unordered_set<networkFace *> ret;
      for (const auto &l : this->Lines)
         for (const auto &f : l->getFaces())
            ret.emplace(f);
      return ret;
   };
   std::unordered_set<networkFace *> getFaces1() const {
      std::unordered_set<networkFace *> ret;
      for (const auto &p : this->getNeighbors())
         ret.merge(p->getFaces0());
      return ret;
   };
   V_netFp getFacesNeumann() const;
   V_netFp getFacesDirichlet() const;
   const std::unordered_set<networkFace *> &getFacesUO() const {
      // std::unordered_set<networkFace *> ret;
      // for (const auto &l : this->Lines)
      // 	for (const auto &f : l->getFaces())
      // 		ret.emplace(f);
      // return ret;
      return this->Faces;
   };

   V_netFp getFacesSort(networkLine *const l) const;
   V_netFp getFacesSort() const;
   // V_netFp getFacesSort2() const;
   // V_netFp getFacesSort() const;
   //--------------
   // bool Replace(netL *oldL, netL *newL)
   // {
   //   this->replace(oldL, newL); //1
   //   oldL->Erase(this);      //2
   //   newL->Add(this);        //3
   //   return true;
   // };
   //--------------
   // V_d getNormal() const override;
   Tddd getNormalTuple() const;
   Tddd getNormalDirichlet() const;
   Tddd getNormalNeumann() const;
   /* ------------------------------------------------------ */
   Tddd getNormal_BEM() const;
   Tddd getNormalDirichlet_BEM() const;
   Tddd getNormalNeumann_BEM() const;
   /* ------------------------------------------------------ */
   Tddd getNormal_BEM_Buffer() const;
   Tddd getNormalDirichlet_BEM_Buffer() const;
   Tddd getNormalNeumann_BEM_Buffer() const;
   /* ------------------------------------------------------ */
   Tddd getNormalSplineKernelAveraged() const;
   Tddd getNormalDirichletSplineKernelAveraged() const;
   Tddd getNormalNeumannSplineKernelAveraged() const;
   /* ------------------------------------------------------ */
   Tddd getNormalSubAreaAveraged() const;
   Tddd getNormalDirichletSubAreaAveraged() const;
   Tddd getNormalNeumannSubAreaAveraged() const;
   /* ------------------------------------------------------ */
   Tddd getNormalAreaAveraged() const;
   Tddd getNormalDirichletAreaAveraged() const;
   Tddd getNormalNeumannAreaAveraged() const;
   /* ------------------------------------------------------ */
   Tddd getNormalAreaAveraged_Buffer() const;
   Tddd getNormalDirichletAreaAveraged_Buffer() const;
   Tddd getNormalNeumannAreaAveraged_Buffer() const;
   /* ------------------------------------------------------ */
   Tddd getNormalInscribedCircleAreaAveraged() const;
   Tddd getNormalDirichletInscribedCircleAreaAveraged() const;
   Tddd getNormalNeumannInscribedCircleAreaAveraged() const;
   /* ------------------------------------------------------ */
   Tddd getNormalAngleAveraged() const;
   Tddd getNormalDirichletAngleAveraged() const;
   Tddd getNormalNeumannAngleAveraged() const;
   /* ------------------------------------------------------ */
   Tddd getNormalOptimum() const;
   Tddd getNormalDirichletOptimum() const;
   Tddd getNormalNeumannOptimum() const;
   /* ------------------------------------------------------ */
   Tddd getNormalQuadInterpAngleAveraged() const;
   Tddd getNormalNeumannQuadInterpAngleAveraged() const;
   Tddd getNormalDirichletQuadInterpAngleAveraged() const;
   /* ------------------------------------------------------ */
   Tddd getNormalArithmeticAveraged() const;
   Tddd getNormalDirichletArithmeticAveraged() const;
   Tddd getNormalNeumannArithmeticAveraged() const;

   void Delete();
   //--------------

   /*SolidAngle_detail

   SolidAngle_detail*/
   /*SolidAngle_detail_code*/
   double getSolidAngle() const;
   double getSolidAngleBuffer() const;
   double getMinimalSolidAngle() const;
   // double getSolidAngle(bool TorF);
   // double getSolidAngle(const V_netFp &faces);
   /*SolidAngle_detail_code*/
};

template <>
double windingNumber(const Tddd &X, const std::vector<networkPoint *> &V_Points) {
   return 1.;
};

template <std::size_t N>
std::array<Tddd, N> ToX(const std::array<networkPoint *, N> &ps) {
   std::array<Tddd, N> ret;
   for (int i = 0; i < N; i++)
      ret[i] = ps[i]->X;
   return ret;
};

T2Tddd ToX(const networkLine *const l) { return l->getLocationsTuple(); };
const Tddd &ToX(const networkPoint *const p) { return p->X; };

template <typename T>
std::vector<Tddd> ToX(const std::unordered_set<T> &ps) {
   std::vector<Tddd> ret(ps.size());
   int i = 0;
   for (const auto &p : ps)
      ret[i++] = ToX(p);
   return ret;
};

template <typename T, std::size_t N>
std::vector<std::array<Tddd, N>> ToX(const std::unordered_set<std::array<T, N>> &ps) {
   std::vector<std::array<Tddd, N>> ret(ps.size());
   int i = 0;
   for (const auto &p : ps)
      ret[i++] = ToX(p);
   return ret;
};

template <typename T, std::size_t N>
std::vector<std::array<Tddd, N>> ToX(const std::vector<std::array<T, N>> &ps) {
   std::vector<std::array<Tddd, N>> ret(ps.size());
   int i = 0;
   for (const auto &p : ps)
      ret[i++] = ToX(p);
   return ret;
};

template <typename T>
std::vector<Tddd> ToX(const std::vector<T> &ps) {
   std::vector<Tddd> ret(ps.size());
   int i = 0;
   for (const auto &p : ps)
      ret[i++] = ToX(p);
   return ret;
};

template <std::size_t N>
std::array<double, N> ToPhi(const std::array<networkPoint *, N> &ps) {
   std::array<double, N> ret;
   for (int i = 0; i < N; i++)
      ret[i] = std::get<0>(ps[i]->phiphin);
   return ret;
};

template <std::size_t N>
std::array<double, N> ToPhin(const std::array<networkPoint *, N> &ps) {
   std::array<double, N> ret;
   for (int i = 0; i < N; i++)
      ret[i] = std::get<0>(ps[i]->phiphin);
   return ret;
};

//@ ------------------------ 抽出用関数 ----------------------- */
//@ --------------------------------------------------------- */
std::unordered_set<networkFace *> extFaces_(const V_netPp &ps) {
   std::unordered_set<networkFace *> ret;
   for (const auto &p : ps)
      for (const auto &f : p->getFaces())
         ret.emplace(f);
   return ret;
};

// 注意：unordered_setからデータを取る場合．順番は保障されない
std::vector<Tddd> extX(const V_netPp &ps) {
   std::vector<Tddd> ret;
   ret.reserve(ps.size());
   for (const auto &p : ps)
      ret.emplace_back(p->X);
   return ret;
};
// std::vector<Tddd> extX(const std::unordered_set<networkPoint *> &ps)
// {
// 	std::vector<Tddd> ret;
// 	for (const auto &p : ps)
// 		ret.emplace_back(p->X);
// 	return ret;
// };
// std::vector<Tddd> extXBuffer(const std::unordered_set<networkPoint *> &ps) {
//    std::vector<Tddd> ret;
//    for (const auto &p : ps)
//       ret.emplace_back(p->getXBuffer());
//    return ret;
// };
// std::vector<Tddd> extXBuffer(const std::vector<networkPoint *> &ps) {
//    std::vector<Tddd> ret;
//    for (const auto &p : ps)
//       ret.emplace_back(p->getXBuffer());
//    return ret;
// };
std::vector<Tddd> extNormals(const V_netPp &ps) {
   std::vector<Tddd> ret;
   ret.reserve(ps.size());
   for (const auto &p : ps)
      ret.emplace_back(p->getNormalTuple());
   return ret;
};

#ifdef BEM
double extPhi(const networkPoint *const p) {
   return std::get<0>(p->phiphin);
};

double extPhin(const networkPoint *const p) {
   return std::get<1>(p->phiphin);
};

V_d extPhi(const V_netPp &ps) {
   V_d ret(ps.size());
   int i = 0;
   for (const auto &p : ps)
      ret[i++] = std::get<0>(p->phiphin);
   return ret;
};

V_d extPhin(const V_netPp &ps) {
   V_d ret(ps.size());
   int i = 0;
   for (const auto &p : ps)
      ret[i++] = std::get<1>(p->phiphin);
   return ret;
};
std::vector<Tdd> extPhiphin(const V_netPp &ps) {
   std::vector<Tdd> ret(ps.size());
   int i = 0;
   for (const auto &p : ps)
      ret[i++] = p->phiphin;
   return ret;
};

#endif

//@ ------------------------------------------------------ */

/*networkPoint_code*/
double distance(const networkPoint *p0, const networkPoint *p1) { return Norm(p0->X - p1->X); };

/* -------------------------------------------------------------------------- */

netP *Point(const netL *const l0,
            const netL *const l1) {
   auto [a, b] = l0->getPoints();
   auto [p, q] = l1->getPoints();
   if (a && (a == p || a == q))
      return a;
   else if (b && (b == p || b == q))
      return b;
   else
      return nullptr;
};
T_3P Points(const netL *const l0,
            const netL *const l1,
            const netL *const l2) {
   return {Point(l2, l0), Point(l0, l1), Point(l1, l2)};
};

netL *Line(const netP *const p0,
           const netP *const p1) {
   for (const auto &l : p0->Lines)
      if (p1 == (*l)(p0))
         return l;
   return nullptr;
};
T_3L Lines(const netP *const p0,
           const netP *const p1,
           const netP *const p2) {
   return {Line(p0, p1), Line(p1, p2), Line(p2, p0)};
};
/* -------------------------------------------------------------------------- */
netL *link(netP *const p0, netP *const p1, Network *const net) {
   // 2022/09/05
   // もし，点が既に繋がっていた場合，それをつなげる線を返す
   // もし，繋がっていない場合，新たな線を生成し，そのポインターを返す．
   // もし，矛盾した繋ぎ方をしていた場合や，つなげることができない点が与えられた場合は例外処理をする
   try {
      if (!p0 || !p1 || p1 == p0) {
         std::stringstream ss;
         ss << p0 << "<--cannot link-->" << p1;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }

      auto l0 = Line(p0, p1);
      auto l1 = Line(p1, p0);

      if (l0 && l1) {
         if (l0 == l1)
            return l0;
         else
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Points are linked but by differrent lines");
      } else if (!l0 && l1) {
         p0->Add(l1);
         return l1;  // only the other side was linked
      } else if (l0 && !l1) {
         p1->Add(l0);
         return l0;  // only the other side was linked
      } else {
         // std::cout << "lineはコンストラクタで引数のオブジェクトを自身のオブジェクトリスト(Points,Faces)に保存し，さらにオブジェクトのLineリストに自身thisを保存する" << std::endl;
         return new networkLine(net, p0, p1);
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
std::array<netL *, 3> link(netP *const p0,
                           netP *const p1,
                           netP *const p2,
                           Network *const net) {
   return {link(p0, p1, net), link(p1, p2, net), link(p2, p0, net)};
};
netL *link(netP *const p0, netP *const p1) {
   return link(p0, p1, p0->getNetwork());
};

/* -------------------------------------------------------------------------- */
class pathInfo {
  public:
   V_netFp face;
   T3Tddd xyz; /* {start_xyz, cwp_xyz, crosspoint_xyz} */
   double incidentAngle;
   double r;
   pathInfo(V_netFp face_IN,
            const T3Tddd &xyz_IN,
            const double incidentAngle_IN,
            const double r_IN) : face(face_IN), xyz(xyz_IN), incidentAngle(incidentAngle_IN), r(r_IN) {};
   ~pathInfo() {
      //    std::cout << Red << "destructed" << colorReset << std::endl;
   }
   void info() {
      std::cout << Blue << "         face :" << face << std::endl;
      std::cout << Blue << "          xyz :" << xyz << std::endl;
      std::cout << Blue << "incidentAngle :" << incidentAngle << std::endl;
      std::cout << Blue << "            r :" << r << colorReset << std::endl;
   };
};
//^ =====================================================
//^     @
//^    / \
//^   @-->@
//^ =====================================================
struct DodecaPoints;
class networkFace : public Triangle {
  public:
   bool public_status;

  protected:
   bool status;

  public:
   T_LLL Lines;  // 2月28日(月)導入
   T_TT Tetras;

   Network *network;
   Network *getNetwork() const { return this->network; };

   void setLines(const T_LLL &ls) { this->Lines = ls; };
   bool getStatus() const { return this->status; };
   void setStatus(bool TorF) { this->status = TorF; };
   bool isStatus(bool TorF) const { return this->status == TorF; };

  public:
   bool MemberQ(const networkLine *const l) const { return MemberQ_(this->Lines, l) && MemberQ_(this->Lines, l) && MemberQ_(this->Lines, l); };
   bool MemberQ(const networkPoint *const p) const { return MemberQ_(this->Points, p); };
   bool AllMemberQ(const T_3P &points) const { return MemberQ_(this->Points, std::get<0>(points)) && MemberQ_(this->Points, std::get<1>(points)) && MemberQ_(this->Points, std::get<2>(points)); };
   bool AnyMemberQ(const T_3P &points) const { return MemberQ_(this->Points, std::get<0>(points)) || MemberQ_(this->Points, std::get<1>(points)) || MemberQ_(this->Points, std::get<2>(points)); };
   // Tdd grid_pull_factor;
   int grid_pull_depth;

   Tddd normalVelocityRigidBody(const Tddd &X) const;
   //% ==================================================================
   //%                        接続情報の事前計算
   //% ==================================================================

   std::array<std::shared_ptr<DodecaPoints>, 3> dodecaPoints;
   void setDodecaPoints();

   //@ ==================================================================
   //@                               BEM
   //@ ==================================================================

   bool Dirichlet;
   bool Neumann;
   bool isDirichlet() const { return this->Dirichlet; };
   bool isNeumann() const { return this->Neumann; };
   int minDepthToCORNER;
   double tension_EMT;

   //@ ------------------------------------------------------------------
   //@ BEMの積分を効率的に行うためにあらかじめ積分情報を作成しておく
   //@ ------------------------------------------------------------------

   using linear_triangle_integration_info = std::tuple<Tdd /*0: 2D parameter {[0,1], [0,1]} (integration variables)*/,
                                                       double /*1: gaussian weight (integration weight)*/,
                                                       Tddd /*2: 3D parameter {xi0=[0,1], xi1=[0,1-xi0]} to move on a triangle and assscociated with the 2D parameter*/,
                                                       Tddd /*3: 3d position vector using {xi0,xi1,xi2}*/,
                                                       Tddd /*4: cross product dX/dxi0 x dX/dxi1 at the point*/,
                                                       double /*5: the norm of the cross product*/>;

   using pseudo_quadratic_triangle_integration_info = std::tuple<Tdd /*0: 2D parameter {[0,1], [0,1]} (integration variables)*/,
                                                                 double /*1: gaussian weight (integration weight)*/,
                                                                 Tddd /*2: 3D parameter {xi0=[0,1], xi1=[0,1-xi0]} to move on a triangle and assscociated with the 2D parameter*/,
                                                                 std::array<T6d, 4> /*3: shape function for the quadratic element*/,
                                                                 Tddd /*4: 3d position vector using {xi0,xi1,xi2}*/,
                                                                 Tddd /*5: cross product dX/dxi0 x dX/dxi1 at the point*/,
                                                                 double /*6: the norm of the cross product*/>;

   std::map<networkPoint *, std::vector<linear_triangle_integration_info>> map_Point_LinearIntegrationInfo;
   std::map<networkPoint *, std::vector<pseudo_quadratic_triangle_integration_info>> map_Point_PseudoQuadraticIntegrationInfo;

   std::map<networkPoint *, std::vector<linear_triangle_integration_info>> map_Point_LinearIntegrationInfo_HigherResolution;
   std::map<networkPoint *, std::vector<pseudo_quadratic_triangle_integration_info>> map_Point_PseudoQuadraticIntegrationInfo_HigherResolution;

   using BEM_IGIGn_info_type = std::tuple<networkPoint *, networkFace *, double, double>;
   std::map<networkPoint *, std::vector<BEM_IGIGn_info_type>> map_Point_BEM_IGIGn_info_init;

   void setIntegrationInfo();

   bool isPseudoQuadraticElement = false;
   bool isLinearElement = false;

   /* -------------------------------------------------------------------------- */

   Tdd phiphin;
   Tddd phinTuple;

   //@========================================================================

#ifdef DEM
  public:
   V_netPp contactP;
#endif

  private:
   VV_d xyzInverse;
   V_netPp XPoints;
   // V_netPp Points;
  public:
   T_PPP Points;
   std::tuple<networkPoint *, networkLine *, networkPoint *, networkLine *, networkPoint *, networkLine *> PLPLPL;
   //% ------------------------------------------------------ */
   std::unordered_set<networkPoint *> ParametricPoints;
   //! ------------------------------------------------------ */
   std::unordered_set<networkPoint *> ContactPoints;
   std::unordered_map<Network *, std::unordered_set<networkPoint *>> map_Net_ContactPoints;

  public:
   //% ------------------------------------------------------ */
   //%           パーティクライズで作られるParametricPoints        */
   //% ------------------------------------------------------ */
   const std::unordered_set<networkPoint *> &getParametricPoints() const { return this->ParametricPoints; };
   void clearParametricPoints() {
      for (const auto &p : this->ParametricPoints)
         delete p;
      this->ParametricPoints.clear();
   };
   void addParametricPoints(networkPoint *const p) { this->ParametricPoints.emplace(p); };
   void addParametricPoints(const V_netPp &P) { this->ParametricPoints.insert(P.begin(), P.end()); };
   void addParametricPoints(const std::unordered_set<networkPoint *> &P) { this->ParametricPoints.insert(P.begin(), P.end()); };
   //! ------------------------------------------------------ */
   //!                          接触の判別                      */
   //! ------------------------------------------------------ */
   const std::unordered_set<networkPoint *> &getContactPoints() const { return this->ContactPoints; };
   const std::unordered_set<networkPoint *> &getContactPoints(Network *const net) const {
      auto it = this->map_Net_ContactPoints.find(net);
      if (it != this->map_Net_ContactPoints.end())
         return it->second;
      else
         return this->map_Net_ContactPoints.at(nullptr);
   };
   const std::unordered_set<networkPoint *> getContactPoints(const std::vector<Network *> &nets) const {
      std::unordered_set<networkPoint *> ret;
      for (const auto &n : nets) {
         auto points = getContactPoints(n);
         ret.insert(points.begin(), points.end());
      }
      return ret;
   };
   void clearContactPoints() {
      this->ContactPoints.clear();
      this->map_Net_ContactPoints.clear();
      this->map_Net_ContactPoints[nullptr] = {};
   };
   //! ------------------------------------------------------ */
   void reverseNormal() {
      // std::reverse(this->Lines.begin(), this->Lines.end());
      this->Lines = Reverse(this->Lines);
      this->setGeometricProperties(ToX(this->setPoints()));  // setBoundsは，setPointsFromCurrentLines()を実行する．
   };
   /* ------------------------------------------------------ */
   Tddd mirror(const Tddd &v) const {
      // 与えられたvの鏡像を返す
      return v - 2. * Dot(v, this->normal) * this->normal;
   };
   Tddd mirrorPosition(const Tddd &v) const {
      // 与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
      auto u = std::get<0>(this->Points)->X;
      return v + 2. * Dot(u - v, this->normal) * this->normal;
   };
   Tddd mirrorPosition(const Tddd &v, const double scale) const {
      // 与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
      return v + scale * Dot(this->X - v, this->normal) * this->normal;
   };
   Tddd mirrorPosition(const networkPoint *p, const double scale = 2.) const {
      // 与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
      return p->X + scale * Dot(this->X - p->X, this->normal) * this->normal;
   };

   networkFace(Network *network_IN, const T_LLL &Lines_IN, T_3P);
   networkFace(Network *network_IN, std::array<networkPoint *, 3> points);
   /*コピーコンストラクタ*/
   networkFace(const netFp f);
   ~networkFace();
   /* ------------------------------------------------------ */
   std::unordered_set<networkPoint *> particlize(const double dx, const V_d &depth_list);
   /* ------------------------------------------------------ */
   V_netPp getXPoints() const;
   bool addXPoint(netP *p) { return network::add(this->XPoints, p); };
   void clearXPoints() { this->XPoints.clear(); };
   bool clearXPoints(netP *p) { return network::erase(this->XPoints, p); };
   V_netPp getPointsPenetrated() const { return this->XPoints; };
   V_netPp getPointsPenetrate() const {
      V_netPp XPoints;
      std::ranges::for_each(this->getLines(), [&](const auto &l) {for (const auto &p : l->XPoints){XPoints.emplace_back(p);} });
      return XPoints;
   };
   V_netPp getPointsPenetrated(const Network *net) const {
      V_netPp ret({});
      for (const auto &p : this->getPointsPenetrated())
         if (p->getNetwork() == net)
            ret.emplace_back(p);
      return ret;
   };
   V_netPp getPointsPenetrate(const Network *net) const {
      V_netPp ret({});
      for (const auto &p : this->getPointsPenetrate())
         if (p->getNetwork() == net)
            ret.emplace_back(p);
      return ret;
   };

   V_netPp getPointsIntersect() const {
      V_netPp ret = this->getPointsPenetrate();
      ret.insert(ret.end(), this->XPoints.begin(), this->XPoints.end());
      return ret;
   };

   V_netLp getLinesPenetrated() const {
      V_netLp ret;
      for (const auto &p : (this->XPoints))
         if (p->getXLine()->penetrateQ())
            ret.emplace_back(p->getXLine());
      return ret;
   };

   V_netFp getNeighbors() const {
      V_netFp ret;
      networkFace *f;
      std::ranges::for_each(this->Lines, [&](const auto &l) {if (f = (*l)(this)){ret.emplace_back(f);}; });
      return ret;
   };

   netL *getLineBetween(const networkFace *const f) const {
      auto [l0, l1, l2] = this->Lines;
      if (f == (*l0)(this))
         return l0;
      else if (f == (*l1)(this))
         return l1;
      else if (f == (*l2)(this))
         return l2;
      else
         return nullptr;
   };

  private:
   /*
   Pointsから，vectorToXPointとなるPointを探し，線の貫いている面の法線方向と比較する．比較した結果が，正なら与えられたPointsは面の後方に位置するとし，trueを返す．
   注意：`this`ネットワークと異なるネットワークに属する点は，干渉点XPointと判断する．
   */
   bool isXPointsBehind(const V_netPp &p) const {
      int s = p.size();
      V_d vectorToXPoint;
      for (auto i = 0; i < s; i++) {
         if (p[i]->getNetwork() != this->network /* means p[i] is interaction netwrok*/) {
            if (p[(i + 1) % s]->getNetwork() == this->network) {
               vectorToXPoint = ToVector(p[i]->X - p[(i + 1) % s]->X);
               if (Dot(ToVector(p[i]->getXFace()->normal), vectorToXPoint) > 1E-12)
                  return true;
            }
            if (p[(i - 1 + s) % s]->getNetwork() == this->network) {
               vectorToXPoint = ToVector(p[i]->X - p[(i - 1 + s) % s]->X);
               if (Dot(p[i]->getXFace()->normal, p[i]->X - p[(i - 1 + s) % s]->X) > 1E-12)
                  return true;
            }
         }
      }

      return false;
   };

  public:
   VV_d getXInverse() const {
      // this->setBounds();
      return this->xyzInverse;
   };
   V_d parameterize(const V_d &xyz_IN) { return Dot(getXInverse(), xyz_IN); };

   T_PPP getPointsFromLines(const T_LLL &lines) const {
      //! p0 <-> line[0] <-> p1 <-> line[1] <-> p2 <-> line[2] を 作成
      auto A = Intersection(std::get<0>(lines)->getPoints(), std::get<2>(lines)->getPoints());
      auto B = Intersection(std::get<1>(lines)->getPoints(), std::get<0>(lines)->getPoints());
      auto C = Intersection(std::get<2>(lines)->getPoints(), std::get<1>(lines)->getPoints());
      if ((A.size() != 1 || B.size() != 1 || C.size() != 1))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Intersection size is not 1");
      else
         return {*A.begin(), *B.begin(), *C.begin()};
   };

   T_PPP getPointsFromLines() const { return getPointsFromLines(this->Lines); };

   T_PPP setPoints(const T_LLL &lines) {
      this->Points = getPointsFromLines(lines);
      std::get<0>(this->PLPLPL) = std::get<0>(this->Points);
      std::get<1>(this->PLPLPL) = std::get<0>(lines);
      std::get<2>(this->PLPLPL) = std::get<1>(this->Points);
      std::get<3>(this->PLPLPL) = std::get<1>(lines);
      std::get<4>(this->PLPLPL) = std::get<2>(this->Points);
      std::get<5>(this->PLPLPL) = std::get<2>(lines);
      return this->Points;
   };

   T_PPP setPoints() { return setPoints(this->Lines); };

   bool setGeometricProperties(const T3Tddd &toX_this_points) {
      /**
       *! 依存関係の明示方法
       * 下のように使う．これによって，setPointsFromLinesした後に，これを実行する必要があることを印象付けることができる．
       * f->setGeometricProperties(f->setPointsFromLines())
       **/
      //@ Pointsに従って，値を設定する．
      //@ 線の持つ点のインターセクションをチェックすることで，面の持つ点を間接的に取得する．
      try {
         /**
         @ networkFacesの持つ
         @ this->Points = {p0,p1,p2}
         @ this->Lines = {l0,l1,l2}
         @ の関係:
         @ 		    p2
         @         /\
         @ 		   /a2\
         @       /    \
         @  l2  /      \ l1
         @ 	   /a0    a1\
         @    -------------
         @  p0      l0      p1
         */
         /*
         b# setBoundsの伝播．おそらくトポロジカルな変更はなく，値の再計算のみでいい
         @ 点の座標変更：PointのsetBounds()を実行
         @   V
         @  PointのsetBounds()   !!!!!! Facesをどうにかして決めれないか？-->全てのsetBounds()が終わった後に，取り込めばいい．
         @   |
         $   |      線のつなぎ変え：LineのsetBounds()を実行
         $   |       V
         $   +-----> LineのsetBounds()
         !           |
         !           +--->FaceのsetBounds() -(LinesとPointsは確定済み)-> normal,area,angles
         *
         @  PointのsetBoundsSingle() ----> CoordinateBounds::setBounds(this->X);　伝播なし
         $  LineのsetBoundsSingle() ----> CoordinateBounds::setBounds(this->X);　伝播なし
         */
         // this->Points = this->getPointsFromLines();
         // T3Tddd p0p1p2_X = ToX(this->Points);
         // T3Tddd toX_this_points = ToX(this_points);
         CoordinateBounds::setBounds(toX_this_points);
         Triangle::setProperties(toX_this_points);
         if (!isFinite(toX_this_points) || !isFinite(this->area) || !isFinite(this->normal) || !isFinite(this->angles)) {
            std::stringstream ss;
            ss << "線が更新されておらずエラーになる可能性がある" << std::endl;
            ss << "全ての線が更新されている必要があるため" << std::endl;
            ss << "this->Points = " << this->Points << std::endl;
            ss << "this->Lines = " << this->Lines << std::endl;
            ss << "toX_this_points = " << toX_this_points << std::endl;
            ss << "this->area = " << this->area << std::endl;
            ss << "this->normal = " << this->normal << std::endl;
            ss << "this->angles = " << this->angles << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
         } else
            return true;
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
   };

   bool setGeometricProperties() {
      return setGeometricProperties(ToX(this->Points));
   };

   //------------------------
   void Delete();
   //------------------------
   // bool Replace(netL *oldL, netL *newL, netF *newF = nullptr) { return network::replace(this, this->Lines, oldL, newL, newF); };
   bool Replace(netL *oldL, netL *newL) {
      // switchでないと，順番に意味のあるFaceではおかしくなるので注意
      if (this->replace(oldL, newL) && oldL->Erase(this) && newL->Add(this))
         return true;
      else
         return false;
   };

   double getSubArea(const networkPoint *const p) const {
      auto [p0, p1, p2] = this->getPoints(p);
      auto a = (p0->X + p1->X) / 2.;
      auto b = (p0->X + p2->X) / 2.;
      auto c = (p0->X + p1->X + p2->X) / 3.;
      return (TriangleArea(p0->X, a, b) + TriangleArea(a, c, b));
   };
   double getInscribedCircleArea() const {
      double a = std::get<0>(this->Lines)->length();
      double b = std::get<1>(this->Lines)->length();
      double c = std::get<2>(this->Lines)->length();
      // Hellon
      return std::sqrt((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) / 16.);
   };
   //////////////////////
   // このfaceの保存状況に従って，lに対する前後のpointへのインデックスを取得できる
   Tiii point_indicies(const netL *l) const {
      try {
         // for (auto i = 0; i < 3; i++)
         // 	if ((l->getPoints()[0] == this->Points[i] /*back*/ && std::get<1>(l->getPoints()) == this->Points[(i + 1) % 3] /*front*/) ||
         // 		(std::get<1>(l->getPoints()) == this->Points[i] /*back*/ && l->getPoints()[0] == this->Points[(i + 1) % 3] /*front*/))
         // 		return {i, (i + 1) % 3, (i + 2) % 3};

         auto p0 = std::get<0>(l->getPoints());
         auto p1 = std::get<0>(l->getPoints());
         if ((p0 == std::get<0>(this->Points) /*back*/ && p1 == std::get<1>(this->Points) /*front*/) ||
             (p1 == std::get<0>(this->Points) /*back*/ && p0 == std::get<1>(this->Points) /*front*/))
            return {0, 1, 2};
         else if ((p0 == std::get<1>(this->Points) /*back*/ && p1 == std::get<2>(this->Points) /*front*/) ||
                  (p1 == std::get<1>(this->Points) /*back*/ && p0 == std::get<2>(this->Points) /*front*/))
            return {1, 2, 0};
         else if ((p0 == std::get<2>(this->Points) /*back*/ && p1 == std::get<0>(this->Points) /*front*/) ||
                  (p1 == std::get<2>(this->Points) /*back*/ && p0 == std::get<0>(this->Points) /*front*/))
            return {2, 0, 1};

         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         std::stringstream ss;
         ss << "このlineを基準としてインデックスをつくれない：この線はこの面のいっぺんではない" << std::endl;
         ss << "setBoundsを忘れていませんか？ 面の線を変更した際などは，setBoundsを忘れないように" << std::endl;
         ss << "FaceのgetPoints()は，lineから毎回間接的に取得することはやめて，setBoundsの際に保存するように変更しました" << std::endl;
         ss << "input line l = " << l << std::endl;
         // ss << "l->Points = " << l->getPoints() << std::endl;
         // ss << "this face->Points = " << this->Points << std::endl;
         // ss << "this face->Lines = " << this->Lines << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };
   // 与えられたpをindex[0]として，this->Pointsのインデックスを返す
   Tiii point_indicies(const networkPoint *const p) const {
      try {
         if (p == std::get<0>(this->Points))
            return {0, 1, 2};
         else if (p == std::get<1>(this->Points))
            return {1, 2, 0};
         else if (p == std::get<2>(this->Points))
            return {2, 0, 1};
         // for (auto i = 0; i < 3; i++)
         // 	if (p == this->Points[i])
         // 		return {i, (i + 1) % 3, (i + 2) % 3};
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         std::stringstream ss;
         ss << "このpointを基準としてインデックスをつくれない：この線はこの面のいっぺんではない" << std::endl;
         ss << "setBoundsを忘れていませんか？ 面の線を変更した際などは，setBoundsを忘れないように" << std::endl;
         ss << "FaceのgetPoints()は，lineから毎回間接的に取得することはやめて，setBoundsの際に保存するように変更しました" << std::endl;
         ss << "input point p = " << p << std::endl;
         // ss << "this face->Points = " << this->Points << std::endl;
         // ss << "this face->Lines = " << this->Lines << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };

   const Tddd &getAngles() const { return this->angles; };
   Tddd getAngles(const netL *l) const {
      try {
         if (l == std::get<0>(this->Lines))
            return this->angles;
         else if (l == std::get<1>(this->Lines))
            return RotateLeft(this->angles, 1);
         else
            return RotateLeft(this->angles, 2);
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
   };
   Tddd getAngles(const networkPoint *const p) const {
      try {
         auto [i, j, k] = point_indicies(p);
         auto [a0, a1, a2] = this->angles;
         if (i == 0)
            return this->angles;
         else if (i == 1)
            return {a1, a2, a0};
         else
            return {a2, a0, a1};
         // return {this->angles[i], /*ここにこの線lが位置する*/ this->angles[j], this->angles[k]};
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
   };
   double getAngle(const networkPoint *p) const {
      if (p == std::get<0>(this->Points))
         return std::get<0>(this->angles);
      else if (p == std::get<1>(this->Points))
         return std::get<1>(this->angles);
      else if (p == std::get<2>(this->Points))
         return std::get<2>(this->angles);
      else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point is not found");
   };

   /* ------------------------------------------------------ */

   /*
   PointsはsetBoundsと同時にかならず，更新される．
   そのため，問題なくりようできる．
   ただし，位置関係を決めておかないと問題が生じる．
   */
   // const V_netPp &getPoints() const { return this->Points; };

   /* ------------------------------------------------------ */
   //% タプル
   /* memo
   @ this->Points = {p0,p1,p2}
   @ this->Lines = {l0,l1,l2}
   @ [p0,p1,p2] = getPoints()
   @ [l0,l1,l2] = getLines()
   @ の関係:
   @ 		   p2
   @         /\
   @        /a2\
   @   l2  /    \ l1
   @      /a0  a1\
   @      ----------
   @   p0     l0    p1
   */
   /* -------------------------------------------------------------------------- */
   const T_PPP &getPoints() const { return this->Points; };

   T_PPP getPoints(const networkPoint *const p) const {
      if (p == std::get<0>(this->Points))
         return this->Points;  // T_PPP{std::get<1>(this->Points), std::get<2>(this->Points), std::get<0>(this->Points)};
      else if (p == std::get<1>(this->Points))
         return RotateLeft(this->Points, 1);  // T_PPP{std::get<2>(this->Points), std::get<0>(this->Points), std::get<1>(this->Points)};
      else
         return RotateLeft(this->Points, 2);
   };
   T_PPP getPoints(const networkLine *const l) const {
      if (l == std::get<0>(this->Lines))
         return this->Points;
      else if (l == std::get<1>(this->Lines))
         return RotateLeft(this->Points);
      else
         return RotateLeft(this->Points, 2);
   };
   /* -------------------------------------------------------------------------- */
   const T_LLL &getLines() const { return this->Lines; };
   T_LLL getLines(const networkLine *const l) const {
      if (l == std::get<0>(this->Lines))
         return this->Lines;
      else if (l == std::get<1>(this->Lines))
         return RotateLeft(this->Lines, 1);  //{std::get<1>(this->Lines), std::get<2>(this->Lines), std::get<0>(this->Lines)};
      else
         return RotateLeft(this->Lines, 2);
   };
   /* -------------------------------------------------------------------------- */
   bool replace(netL *const oldL, netL *const newL) {
      if (std::get<0>(this->Lines) == oldL) {
         std::get<0>(this->Lines) = newL;
         // this->Lines[0] = newL;
         return true;
      } else if (std::get<1>(this->Lines) == oldL) {
         std::get<1>(this->Lines) = newL;
         // this->Lines[1] = newL;
         return true;
      } else if (std::get<2>(this->Lines) == oldL) {
         std::get<2>(this->Lines) = newL;
         // this->Lines[2] = newL;
         return true;
      } else
         return false;
   };
   bool replace(netP *const oldP, netP *const newP) {
      if (std::get<0>(this->Points) == oldP) {
         std::get<0>(this->Points) = newP;
         return true;
      } else if (std::get<1>(this->Points) == oldP) {
         std::get<1>(this->Points) = newP;
         return true;
      } else if (std::get<2>(this->Points) == oldP) {
         std::get<2>(this->Points) = newP;
         return true;
      } else
         return false;
   };
   std::vector<networkFace *> getNeighbors(const networkFace *const obj) const {
      return {(*std::get<0>(this->Lines))(obj), (*std::get<1>(this->Lines))(obj), (*std::get<2>(this->Lines))(obj)};
   };
   /* ------------------------------------------------------ */
   std::tuple<T_PPP, T_LLL> getPointsLinesTuple(const networkLine *const l) const {
      return {this->getPoints(l), getLines(l)};
   };
   /* ------------------------------------------------------ */
   netP *getPointFront(const netL *l) const { return std::get<1>(getPoints(l)); };
   netP *getPointBack(const netL *l) const { return std::get<0>(getPoints(l)); };
   netP *getPointOpposite(const netL *l) const { return std::get<2>(getPoints(l)); };
   /* ------------------------------------------------------ */
   netL *getLineFront(const networkPoint *p) const { return std::get<0>(getLinesTupleFrom(p)); };
   netL *getLineBack(const networkPoint *p) const { return std::get<2>(getLinesTupleFrom(p)); };
   netL *getLineOpposite(const networkPoint *p) const { return std::get<1>(getLinesTupleFrom(p)); };
   /* ------------------------------------------------------ */
   // netL *getLine(int i0, int i1) const
   // {
   // 	if (std::get<0>(Points) && std::get<1>(Points))
   // 		for (const auto &line : std::get<0>(Points)->Lines)
   // 			if (std::get<1>(Points) == (*line)(std::get<0>(Points)))
   // 				return line;
   // 	return nullptr;
   // };

   // LinesTupleだけにしよう！！えらーがある　switch

   netL *getLine(const netL *l, int j = 0) const {
      // int s = this->Lines.size();
      // if (s != 3)
      // {
      // 	std::stringstream ss;
      // 	ss << "辺の数が３ではない：this->Lines" << this->Lines;
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      // }
      // int s = 3;
      if (std::get<0>(this->Lines) == l) {
         if (j % 3 == 0)
            return std::get<0>(this->Lines);
         else if (j % 3 == 1)
            return std::get<1>(this->Lines);
         else
            return std::get<2>(this->Lines);
      } else if (std::get<1>(this->Lines) == l) {
         if ((1 + j) % 3 == 0)
            return std::get<0>(this->Lines);
         else if ((1 + j) % 3 == 1)
            return std::get<1>(this->Lines);
         else
            return std::get<2>(this->Lines);
      } else if (std::get<2>(this->Lines) == l) {
         if ((2 + j) % 3 == 0)
            return std::get<0>(this->Lines);
         else if ((2 + j) % 3 == 1)
            return std::get<1>(this->Lines);
         else
            return std::get<2>(this->Lines);
      } else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
   netL *getLineFront(const netL *l) const { return getLine(l, 1); };
   netL *getLineBack(const netL *l) const { return getLine(l, -1); };
   netF *getFaceFront(const netL *l) { return (*getLine(l, 1))(this); };
   netF *getFaceBack(const netL *l) { return (*getLine(l, -1))(this); };
   // 与えられたpを0番とするpointsベクトル
   // これがないと，objファイルを生成するときにエラーとなる．
   // というか，引数がないgetPointsを作り，this->Poitnsを出力するようにするとエラーになる．
   // つまり，getPoints()は，getPoints(nullptr)としてデフォルト引数でもう一つのgetPointsが実行されないといけない．
   // 初めの段階では，this->Pointsが設定されていないことがあるからだ．
   // 線を通して点を取得するかどうか選べるように，getPointsThroughLines関数を作るべきだろう
   //  V_netPp getPoints(const networkPoint *p) const
   //  {
   //  	if (p == std::get<0>(this->Points))
   //  		return this->Points;
   //  	if (p == std::get<1>(this->Points))
   //  		return {std::get<1>(this->Points), std::get<2>(this->Points), std::get<0>(this->Points)};
   //  	if (p == std::get<2>(this->Points))
   //  		return {std::get<2>(this->Points), std::get<0>(this->Points), std::get<1>(this->Points)};
   //  	else
   //  		return this->Points;
   //  };

   // 与えられたpから出発するlinesベクトル
   //  V_netLp getLinesFrom(const networkPoint *const p) const
   //  {
   //  	try
   //  	{
   //  		if (p == std::get<0>(this->Points))
   //  			return this->Lines;
   //  		else if (p == std::get<1>(this->Points))
   //  			return RotateLeft(this->Lines, 1);
   //  		else if (p == std::get<2>(this->Points))
   //  			return RotateLeft(this->Lines, 2);
   //  		else
   //  			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   //  	}
   //  	catch (std::exception &e)
   //  	{
   //  		std::cerr << e.what() << colorReset << std::endl;
   //  		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   //  	};
   //  };
   T_LLL getLinesTupleFrom(const networkPoint *const p) const {
      try {
         if (p == std::get<0>(this->Points))
            return this->Lines;
         else if (p == std::get<1>(this->Points))
            return RotateLeft(this->Lines, 1);
         else if (p == std::get<2>(this->Points))
            return RotateLeft(this->Lines, 2);
         else
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
   };
   //------------------------
   // netF *getFaceFront(const networkPoint *p) const { return (*(getLinesFrom(p)[0]))(this); };
   // netF *getFaceBack(const networkPoint *p) const { return (*(getLinesFrom(p)[0]))(this); };
   // netF *getFaceOpposite(const networkPoint *p) const { return (*(getLinesFrom(p)[0]))(this); };
   netF *getFaceFront(const networkPoint *p) const { return (*(std::get<0>(getLinesTupleFrom(p))))(this); };
   netF *getFaceBack(const networkPoint *p) const { return (*(std::get<2>(getLinesTupleFrom(p))))(this); };
   netF *getFaceOpposite(const networkPoint *p) const { return (*(std::get<1>(getLinesTupleFrom(p))))(this); };

   /* -------------------------------------------------------------------------- */

   template <size_t N>
   std::array<networkPoint *, N> getPoints() const {
      static_assert(N == 3 || N == 6, "Unsupported getPoints. Only 3 or 6 are supported.");
      if constexpr (N == 3) {
         return ToArray(this->Points);
      } else if constexpr (N == 6) {
         auto [l0, l1, l2] = this->getLines();
         auto [p0_f0, p1_f0, p2_f0] = (*l1)(this)->getPoints(l1);  // f0
         auto [p0_f1, p1_f1, p2_f1] = (*l2)(this)->getPoints(l2);  // f1
         auto [p0_f2, p1_f2, p2_f2] = (*l0)(this)->getPoints(l0);  // f2
         return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
      }
   }
   //
   template <size_t N, bool networkLine::*BoolMember = nullptr>
   std::array<networkPoint *, N> getPoints(const networkPoint *P = nullptr) const {
      static_assert(N == 3 || N == 6, "Unsupported getPoints. Only 3 or 6 are supported.");
      if constexpr (N == 3) {
         if (P == std::get<0>(this->Points) || P == nullptr)
            return ToArray(this->Points);
         else if (P == std::get<1>(this->Points))
            return ToArray(RotateLeft(this->Points, 1));
         else
            return ToArray(RotateLeft(this->Points, 2));
      } else if constexpr (N == 6) {
         /*              0
          *              0  p2_f0
          *            /  \
          *   p0_f0  3/-l1-\5  p1_f0
          *         /  \  /  \
          *  1    1/---4(p)---\2   2　p2_f2
          */
         auto [l0, l1, l2] = this->getLinesTupleFrom(P);
         auto [p0_f0, p1_f0, p2_f0] = (*l1)(this)->getPoints(l1);  // f0
         auto [p0_f1, p1_f1, p2_f1] = (*l2)(this)->getPoints(l2);  // f1
         auto [p0_f2, p1_f2, p2_f2] = (*l0)(this)->getPoints(l0);  // f2
         if (BoolMember) {
            if (l0->*BoolMember) p2_f2 = nullptr;
            if (l1->*BoolMember) p2_f0 = nullptr;
            if (l2->*BoolMember) p2_f1 = nullptr;
         }
         return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
      }
   }

   template <size_t N, bool networkLine::*BoolMember = nullptr>
   std::array<networkPoint *, N> getPoints(const networkLine *L = nullptr) const {
      static_assert(N == 3 || N == 6, "Unsupported getPoints. Only 3 or 6 are supported.");
      if constexpr (N == 3) {
         return ToArray(this->getPoints(L));
      } else if constexpr (N == 6) {
         auto [l1, l2, l0] = this->getLines(L);
         auto [p0_f0, p1_f0, p2_f0] = (*l1)(this)->getPoints(l1);  // f0
         auto [p0_f1, p1_f1, p2_f1] = (*l2)(this)->getPoints(l2);  // f1
         auto [p0_f2, p1_f2, p2_f2] = (*l0)(this)->getPoints(l0);  // f2
         if (BoolMember) {
            if (l0->*BoolMember) p2_f2 = nullptr;
            if (l1->*BoolMember) p2_f0 = nullptr;
            if (l2->*BoolMember) p2_f1 = nullptr;
         }
         return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
      }
   }

   // template <size_t N, bool networkLine::*BoolMember = nullptr>
   // std::array<std::tuple<networkPoint *, networkFace *>, N> getPointsFaces(const networkLine *L = nullptr) const {
   //    static_assert(N == 3 || N == 6, "Unsupported getPoints. Only 3 or 6 are supported.");
   //    if constexpr (N == 3) {
   //       auto [p0, p1, p2] = this->getPoints(L);
   //       return {{p0, this}, {p1, this}, {p2, this}};
   //    } else if constexpr (N == 6) {
   //       auto [l1, l2, l0] = this->getLines(L);
   //       std::tuple<networkPoint *, networkFace *> p0_f0, p1_f0, p2_f0, p0_f1, p1_f1, p2_f1, p0_f2, p1_f2, p2_f2;
   //       {
   //          auto f0 = (*l1)(this);
   //          auto [p0, p1, p2] = f0->getPoints(l1);  // f0
   //          p0_f0 = {p0, f0};
   //          p1_f0 = {p1, f0};
   //          p2_f0 = {p2, f0};
   //       }
   //       {
   //          auto f1 = (*l2)(this);
   //          auto [p0, p1, p2] = f1->getPoints(l2);  // f1
   //          p0_f1 = {p0, f1};
   //          p1_f1 = {p1, f1};
   //          p2_f1 = {p2, f1};
   //       }
   //       {
   //          auto f2 = (*l0)(this);
   //          auto [p0, p1, p2] = f2->getPoints(l0);  // f2
   //          p0_f2 = {p0, f2};
   //          p1_f2 = {p1, f2};
   //          p2_f2 = {p2, f2};
   //       }
   //       if (BoolMember) {
   //          if (l0->*BoolMember) p2_f2 = nullptr;
   //          if (l1->*BoolMember) p2_f0 = nullptr;
   //          if (l2->*BoolMember) p2_f1 = nullptr;
   //       }
   //       return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
   //    }
   // }
   /* -------------------------------------------------------------------------- */

   std::tuple<networkPoint *, networkLine *, networkPoint *, networkLine *, networkPoint *, networkLine *>
   getPointsAndLines(const networkPoint *origin) const {
      auto [P0, L0, P1, L1, P2, L2] = this->PLPLPL;
      if (P0 == origin || nullptr == origin) {
         return this->PLPLPL;
      } else if (P1 == origin) {
         return {P1, L1, P2, L2, P0, L0};
      } else if (P2 == origin) {
         return {P2, L2, P0, L0, P1, L1};
      } else
         return this->PLPLPL;
   };

   std::tuple<networkPoint *, networkLine *, networkPoint *, networkLine *, networkPoint *, networkLine *>
   getPointsAndLinesFromNearest(const networkPoint *origin) const {
      auto [P0, L0, P1, L1, P2, L2] = this->PLPLPL;
      if (P0 == origin || nullptr == origin) {
         return this->PLPLPL;
      } else if (P1 == origin || (Norm(origin->X - P0->X) > Norm(origin->X - P1->X) && Norm(origin->X - P2->X) > Norm(origin->X - P1->X))) {
         return {P1, L1, P2, L2, P0, L0};
      } else if (P2 == origin || (Norm(origin->X - P0->X) > Norm(origin->X - P2->X) && Norm(origin->X - P1->X) > Norm(origin->X - P2->X))) {
         return {P2, L2, P0, L0, P1, L1};
      } else
         return this->PLPLPL;
   };

   std::tuple<networkPoint *, networkLine *, networkPoint *, networkLine *, networkPoint *, networkLine *>
   getPointsAndLines(const networkLine *origin) const {
      auto [P0, L0, P1, L1, P2, L2] = this->PLPLPL;
      if (L0 == origin || nullptr == origin) {
         return this->PLPLPL;
      } else if (L1 == origin) {
         return {P1, L1, P2, L2, P0, L0};
      } else if (L2 == origin) {
         return {P2, L2, P0, L0, P1, L1};
      } else
         return this->PLPLPL;
   };

   /* -------------------------------------------------------------------------- */

   std::tuple<netPp, netPp, netPp, netPp, netPp, netPp> get6PointsTuple(const networkLine *l /*基準*/) const {
      try {
         auto [l0, l1, l2] = this->getLines(l);  // 修正した2022/03/21
         // if (l1 == l)
         // {
         // 	l0 = std::get<1>(this->Lines);
         // 	l1 = std::get<2>(this->Lines);
         // 	l2 = std::get<0>(this->Lines);
         // }
         // else if (l2 == l)
         // {
         // 	l0 = std::get<2>(this->Lines);
         // 	l1 = std::get<0>(this->Lines);
         // 	l2 = std::get<1>(this->Lines);
         // }

         /*
          *                  0  (p2_f0)
          *                 / \
          *                /f0 \
          *               /     \
          *    (p0_f0)  3/--l-l0-\5 (p0_f2)
          *             / \     / \
          *            /   l1  l2  \
          *           / f1  \ /  f2 \
          * (p2_f1) 1/-------4-------\2  (p2_f2)
          *               (p0_f1)
          */
         auto [p0_f0, p1_f0, p2_f0] = (*l0)(this)->getPoints(l0);  // f0
         auto [p0_f1, p1_f1, p2_f1] = (*l1)(this)->getPoints(l1);  // f1
         auto [p0_f2, p1_f2, p2_f2] = (*l2)(this)->getPoints(l2);  // f2
         return {p2_f0,
                 p2_f1,
                 p2_f2,
                 p0_f0,
                 p0_f1,
                 p0_f2};
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
   };

   std::tuple<netPp, netPp, netPp, netPp, netPp, netPp> get6PointsTuple(const networkPoint *p) const {
      try {
         if (p == nullptr || !MemberQ_(this->Points, p)) {
            p = std::get<0>(this->Points);
         }
         auto [l0, l1, l2] = this->getLinesTupleFrom(p);
         /*              0
          *              0
          *            /  \
          *          3/-l1-\5
          *         /  \  /  \
          *  1    1/---4(p)---\2   2
          */
         return (*l1)(this)->get6PointsTuple(l1);
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
   };

   //@ 2024/02/02本格的に使い始める
   std::array<networkPoint *, 6> get6Points(const networkPoint *p) const {
      if (p == nullptr || !MemberQ_(this->Points, p))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p is nullptr or not found in this->Points");
      /*
       *                  0  (p2_f0)
       *                 / \
       *                /f0 \
       *               /     \
       *    (p0_f0)  3/--l-l0-\5 (p0_f2)
       *             / \     / \
       *            /   l1  l2  \
       *           / f1  \ /  f2 \
       * (p2_f1) 1/-------4-------\2  (p2_f2)
       *               (p0_f1)
       */
      auto [l0, l1, l2] = this->getLinesTupleFrom(p);
      auto [p0_f0, p1_f0, p2_f0] = (*l0)(this)->getPoints(l0);  // f0
      auto [p0_f1, p1_f1, p2_f1] = (*l1)(this)->getPoints(l1);  // f1
      auto [p0_f2, p1_f2, p2_f2] = (*l2)(this)->getPoints(l2);  // f2

      if (l0->CORNER)
         p2_f0 = nullptr;
      if (l1->CORNER)
         p2_f1 = nullptr;
      if (l2->CORNER)
         p2_f2 = nullptr;
      return {p2_f0,
              p2_f1,
              p2_f2,
              p0_f0,
              p0_f1,
              p0_f2};
   };

   //-------------------------
   Tddd getMeanX() const { return Mean(getLocationsTuple()); };

   T3Tddd getLocationsTuple() const {
      return {std::get<0>(this->Points)->X,
              std::get<1>(this->Points)->X,
              std::get<2>(this->Points)->X};
   };
   T3Tddd getXVertices() const {
      return {std::get<0>(this->Points)->X,
              std::get<1>(this->Points)->X,
              std::get<2>(this->Points)->X};
   };

   netFp divide(netLp DivL, netLp newDivL, netLp newMidL, int type);

   /* ------------------------------------------------------ */
   bool isFacing(const networkFace *const F, const double rad = 1E-10) const {
      return isFlat(this->normal, -F->normal, rad) || isFlat(this->normal, F->normal, rad);
   };
};
/* -------------------------------------------------------------------------- */

T3Tddd ToX(const networkFace *const f) { return ToX(f->getPoints()); };
std::vector<T3Tddd> ToX(const std::vector<networkFace *> &fs) {
   std::vector<T3Tddd> ret(fs.size());
   int i = 0;
   for (const auto &f : fs)
      ret[i++] = ToX(f->getPoints());
   return ret;
};
std::vector<T3Tddd> ToX(const std::unordered_set<networkFace *> &fs) {
   std::vector<T3Tddd> ret(fs.size());
   int i = 0;
   for (const auto &f : fs)
      ret[i++] = ToX(f->getPoints());
   return ret;
};
/* -------------------------------------------------------------------------- */
Tddd TriangleNormal(const networkFace *const f) { return TriangleNormal(ToX(f)); };
//@ ------------------------ 抽出用関数など ----------------------- */
std::vector<Tddd> extX(const std::unordered_set<networkFace *> &fs) {
   std::vector<Tddd> ret;
   for (const auto &f : fs)
      ret.emplace_back(f->X);
   return ret;
};
std::vector<T3Tddd> extVertices(const std::unordered_set<networkFace *> &fs) {
   std::vector<T3Tddd> ret(fs.size());
   int i = 0;
   for (const auto &f : fs)
      ret[i++] = f->getXVertices();
   return ret;
};
std::vector<T2Tddd> extX(const std::unordered_set<networkLine *> &ls) {
   std::vector<T2Tddd> ret(ls.size());
   int i = 0;
   for (const auto &l : ls)
      ret[i++] = l->getLocationsTuple();
   return ret;
};
std::vector<Tddd> extNormals(const V_netFp &fs) {
   std::vector<Tddd> ret(fs.size());
   int i = 0;
   for (const auto &f : fs)
      ret[i++] = f->normal;
   return ret;
};
// Tddd vectorToTriangle(const networkFace *f, const Tddd &a) {
//    auto n = f->normal;
//    return n * Dot(n, f->X - a);
// };
//@ ------------------------------------------------------ */

// struct interpolationTriangleQuadByFixedRange3D_use_only_good_lines : public interpolationTriangleQuadByFixedRange3D {
//    bool l0_isGoodForQuadInterp;
//    bool l1_isGoodForQuadInterp;
//    bool l2_isGoodForQuadInterp;
//    std::tuple<networkLine *, networkLine *, networkLine *> Lines;
//    std::tuple<networkPoint *, networkPoint *, networkPoint *, networkPoint *, networkPoint *, networkPoint *> Points;
//    interpolationTriangleQuadByFixedRange3D_use_only_good_lines(const networkFace *const f, const networkLine *const l)
//        : interpolationTriangleQuadByFixedRange3D(),
//          l0_isGoodForQuadInterp(true),
//          l1_isGoodForQuadInterp(true),
//          l2_isGoodForQuadInterp(true),
//          Lines(f->getLines(l)),
//          Points(f->get6PointsTuple(l)) {
//       this->s = ToX(Points);
//       if (!(l0_isGoodForQuadInterp = std::get<0>(Lines)->isGoodForQuadInterp()))
//          this->approxP0();
//       if (!(l1_isGoodForQuadInterp = std::get<1>(Lines)->isGoodForQuadInterp()))
//          this->approxP1();
//       if (!(l2_isGoodForQuadInterp = std::get<2>(Lines)->isGoodForQuadInterp()))
//          this->approxP2();
//       //
//       /*
//        approx 0
//          Q1 for l2
//          Q2 for l1
//       *--Q0 for l0--*
//        \   /   \   /
//         \/      \ /
//       p1*---l0---*p0
//         / \  f  / \
//        /   \   /   \
// Q1*------*------*Q2
//                  p2
//       この修正によって，角ではphi，phinが不連続となる.
//       */
//    };
// };
// struct interpolationTriangleQuadByFixedRange3D_use_only_good_lines_Geo : public interpolationTriangleQuadByFixedRange3D {
//    bool l0_isGoodForQuadInterp;
//    bool l1_isGoodForQuadInterp;
//    bool l2_isGoodForQuadInterp;
//    std::tuple<networkLine *, networkLine *, networkLine *> Lines;
//    std::tuple<networkPoint *, networkPoint *, networkPoint *, networkPoint *, networkPoint *, networkPoint *> Points;
//    interpolationTriangleQuadByFixedRange3D_use_only_good_lines_Geo(const networkFace *const f, const networkLine *const l)
//        : interpolationTriangleQuadByFixedRange3D(),
//          l0_isGoodForQuadInterp(true),
//          l1_isGoodForQuadInterp(true),
//          l2_isGoodForQuadInterp(true),
//          Lines(f->getLines(l)),
//          Points(f->get6PointsTuple(l)) {
//       this->s = ToX(Points);
//       if (!(l0_isGoodForQuadInterp = std::get<0>(Lines)->isGoodForQuadInterp_Geo()))
//          this->approxP0();
//       if (!(l1_isGoodForQuadInterp = std::get<1>(Lines)->isGoodForQuadInterp_Geo()))
//          this->approxP1();
//       if (!(l2_isGoodForQuadInterp = std::get<2>(Lines)->isGoodForQuadInterp_Geo()))
//          this->approxP2();
//       //
//       /*
//        approx 0
//          Q1 for l2
//          Q2 for l1
//       *--Q0 for l0--*
//        \   /   \   /
//         \/      \ /
//       p1*---l0---*p0
//         / \  f  / \
//        /   \   /   \
// Q1*------*------*Q2
//                  p2
//       この修正によって，角ではphi，phinが不連続となる.
//       */
//    };
// };
/* ------------------------------------------------------ */
// inline void Buckets<networkPoint *>::add(const V_netPp &ps)
// {
// 	for (const auto p : ps)
// 		add(p->X, p);
// };
// inline void Buckets<networkPoint *>::add(const std::unordered_set<networkPoint *> &ps)
// {
// 	for (const auto p : ps)
// 		add(p->X, p);
// };
/* ------------------------------------------------------ */
std::vector<netL *> link(const V_netPp &obj, Network *net) {
   try {
      std::vector<netL *> ret;
      int s = obj.size();
      for (int i = 0; i < s; i++) {
         auto l = link(obj[i], obj[(i + 1) % s], net);
         ret.emplace_back(l);
      }
      return ret;
   } catch (const error_message &e) {
      Print(obj);
      if (DuplicateFreeQ(obj))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "no duplication.....???");
      else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "duplication found!");
   }
};
netL *unlink(netP *obj, netP *obj_) {
   if (obj_ == obj)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "a point is trying to unlink itself!");

   if (!obj || !obj_ /*if NULL*/)
      return nullptr;

   auto line = obj->getLineBetween(obj_);
   auto line_ = obj_->getLineBetween(obj);

   if ((line == line_) && line) {
      line->Erase(obj);
      line->Erase(obj_);
      obj->Erase(line);
      obj_->Erase(line);
      return line;
   } else {
      std::cout << Red << obj << " and " << obj_ << " are not linked" << colorReset << std::endl;
      return nullptr;
   }
};
//@ ------------------------------------------------------ */
//@                         extract                        */
//@ ------------------------------------------------------ */
/*
 * under scorer _ means the function returns FLATTEND list
 */
std::unordered_set<networkPoint *> extPointsCORNER_(const std::vector<networkPoint *> &ps) {
   std::unordered_set<networkPoint *> ret;
   for (const auto &p : ps)
      if (p->CORNER)
         ret.emplace(p);
   return ret;
};
std::unordered_set<networkPoint *> extPointsCORNER_(const std::unordered_set<networkPoint *> &ps) {
   std::unordered_set<networkPoint *> ret;
   for (const auto &p : ps)
      if (p->CORNER)
         ret.emplace(p);
   return ret;
};
std::unordered_set<networkPoint *> extPointsCornerOrNeumann_(const std::vector<networkPoint *> &ps) {
   std::unordered_set<networkPoint *> ret;
   for (const auto &p : ps)
      if (p->CORNER || p->Neumann)
         ret.emplace(p);
   return ret;
};
std::unordered_set<networkPoint *> extPointsCornerOrNeumann_(const std::unordered_set<networkPoint *> &ps) {
   std::unordered_set<networkPoint *> ret;
   for (const auto &p : ps)
      if (p->CORNER || p->Neumann)
         ret.emplace(p);
   return ret;
};
// b! ------------------------------------------------------ */
// b! ---------------------- extLines ---------------------- */
// b! ------------------------------------------------------ */
/* ------------------ for unordered_set ----------------- */
std::unordered_set<networkLine *> extLinesCORNER_(const std::unordered_set<networkFace *> &fs) {
   std::unordered_set<networkLine *> ret;
   for (const auto &f : fs)
      std::ranges::for_each(f->getLines(), [&](const auto &l) {if (l->CORNER){ret.emplace(l);}; });
   return ret;
};
std::unordered_set<networkLine *> extLines_(const std::unordered_set<networkFace *> &fs) {
   std::unordered_set<networkLine *> ret;
   for (const auto &f : fs)
      std::ranges::for_each(f->getLines(), [&](const auto &l) { ret.emplace(l); });
   return ret;
};
std::unordered_set<networkLine *> extLinesCORNER_(const std::unordered_set<networkPoint *> &ps) {
   std::unordered_set<networkLine *> ret;
   for (const auto &p : ps)
      for (const auto &l : p->getLines())
         if (l->CORNER)
            ret.emplace(l);
   return ret;
};
std::unordered_set<networkLine *> extLinesCORNER_(const networkPoint *p) {
   std::unordered_set<networkLine *> ret;
   for (const auto &l : p->getLines())
      if (l->CORNER)
         ret.emplace(l);
   return ret;
};
std::unordered_set<networkLine *> extLines_(const std::unordered_set<networkPoint *> &ps) {
   std::unordered_set<networkLine *> ret;
   for (const auto &p : ps)
      for (const auto &l : p->getLines())
         ret.emplace(l);
   return ret;
};
/* --------------------- for vector --------------------- */
std::unordered_set<networkLine *> extLinesCORNER_(const std::vector<networkFace *> &fs) {
   std::unordered_set<networkLine *> ret;
   for (const auto &f : fs)
      std::ranges::for_each(f->getLines(), [&](const auto &l) {if (l->CORNER){ret.emplace(l);}; });
   return ret;
};
std::unordered_set<networkLine *> extLines_(const std::vector<networkFace *> &fs) {
   std::unordered_set<networkLine *> ret;
   for (const auto &f : fs)
      std::ranges::for_each(f->getLines(), [&](const auto &l) { ret.emplace(l); });
   return ret;
};
std::unordered_set<networkLine *> extLinesCORNER_(const std::vector<networkPoint *> &ps) {
   std::unordered_set<networkLine *> ret;
   for (const auto &p : ps)
      for (const auto &l : p->getLines())
         if (l->CORNER)
            ret.emplace(l);
   return ret;
};
std::unordered_set<networkLine *> extLines_(const std::vector<networkPoint *> &ps) {
   std::unordered_set<networkLine *> ret;
   for (const auto &p : ps)
      for (const auto &l : p->getLines())
         ret.emplace(l);
   return ret;
};
// b! ------------------------------------------------------ */
// b! ------------------------------------------------------ */
// b! ------------------------------------------------------ */
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
V_d extLength(const std::unordered_set<networkLine *> &ls) {
   V_d ret(ls.size());
   int i = 0;
   for (const auto &l : ls)
      ret[i++] = l->length();
   return ret;
};
//
V_d extLength(const V_netLp &ls) {
   V_d ret(ls.size());
   int i = 0;
   for (const auto &l : ls)
      ret[i++] = l->length();
   return ret;
};
Tddd extLength(const std::array<networkLine *, 3> &ls) {
   return {std::get<0>(ls)->length(),
           std::get<1>(ls)->length(),
           std::get<2>(ls)->length()};
};
V_d extTensionEMT(const V_netLp &ls) {
   V_d ret;
   for (const auto &l : ls)
      ret.emplace_back(l->tension_EMT);
   return ret;
};
V_d extAreas(const V_netFp &fs) {
   V_d ret(fs.size());
   int i = 0;
   for (const auto &f : fs)
      ret[i++] = f->area;
   return ret;
};

V_d extractAreas(const V_netFp &fs) {
   V_d ret(fs.size());
   int i = 0;
   for (const auto &f : fs)
      ret[i++] = f->area;
   return ret;
};
template <class T>
std::vector<Tddd> extractXtuple(const std::unordered_set<T *> &object) {
   std::vector<Tddd> ret(object.size());
   int i = 0;
   for (auto it = object.begin(); it != object.end(); ++it)
      ret[i++] = (*it)->X;
   return ret;
};
template <class T>
std::vector<Tddd> extractXtuple(const std::vector<T *> &object) {
   std::vector<Tddd> ret(object.size());
   for (auto i = 0; i < object.size(); i++)
      ret[i] = object[i]->X;
   return ret;
};

template <class T>
VV_d extractX(const std::unordered_set<T *> &object) {
   VV_d ret(object.size(), V_d(3));
   int i = 0;
   for (auto it = object.begin(); it != object.end(); ++it)
      ret[i++] = ToVector((*it)->X);
   return ret;
};
template <class T>
VV_d extractX(const std::vector<T *> &object) {
   VV_d ret(object.size(), V_d(3));
   for (auto i = 0; i < object.size(); i++)
      ret[i] = ToVector(object[i]->X);
   return ret;
};
// 2021/09/06追加
std::vector<Tddd> extXtuple(const V_netPp &points) {
   std::vector<Tddd> ret(points.size());
   int i = 0;
   for (const auto &p : points)
      ret[i++] = p->X;
   return ret;
};
std::vector<Tddd> extXtuple(const V_netFp &points) {
   std::vector<Tddd> ret(points.size());
   int i = 0;
   for (const auto &p : points)
      ret[i++] = p->X;
   return ret;
};

// template <class T>
// VVV_d extractX(const std::vector<std::vector<T *>> &object)
// {
// 	VVV_d ret;
// 	ret.reserve(object.size());
// 	for (const auto &obj : object)
// 		ret.emplace_back(obj3D::extractX(obj));
// 	return ret;
// };
// T3Tddd extractX(networkFace const *f)
// {
// 	auto ps = f->getPoints();
// 	return std::make_tuple(ps[0]->X, ps[1]->X, ps[2]->X);
// };
T3Tddd extractXtuple(networkFace const *f) {
   auto [p0, p1, p2] = f->getPoints();
   // return std::make_tuple(p0->X, p1->X, p2->X);
   return {p0->X, p1->X, p2->X};
};
//
std::vector<std::vector<Tddd>> extractXtuple(const std::vector<networkLine *> &lines) {
   std::vector<std::vector<Tddd>> ret;
   for (const auto &l : lines) {
      ret.push_back({});
      auto [p, q] = l->getPoints();
      ret.rbegin()->emplace_back(p->X);
      ret.rbegin()->emplace_back(q->X);
   }
   return ret;
};
/* -------------------------------------------------------------------------- */
#include "searcher.hpp"
/* -------------------------------------------------------------------------- */
/*      @     */
/*     /|\    */
/*    / | \   */
/*   /  |  \  */
/*  @---@---@ */

class networkTetra : public Tetrahedron {
  public:
   Network *network;
   // Network connectivities
   T_4P Points;
   T_6L Lines;
   T_4F Faces;
   networkTetra(Network *const network_IN,
                const T_4P &points,
                const T_6L &lines,
                const T_4F &faces);
};
/* -------------------------------------------------------------------------- */
netF *genFace(Network *const net,
              netL *const l0,
              netL *const l1,
              netL *const l2) {
   try {
      if (l0 && l1 && l2 && l0 != l1 && l0 != l2 && l1 != l2) {
         auto intx_faces = Intersection(l0->getFaces(), l1->getFaces(), l2->getFaces());
         auto s = intx_faces.size();
         if (s == 0)
            return new networkFace(net, T_3L{l0, l1, l2}, Points(l0, l1, l2));
         else if (s == 1)
            return intx_faces[0];
         else {
            std::stringstream ss;
            ss << "too many inteserctions of faces : " << intx_faces;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
         }
      } else {
         std::stringstream ss;
         ss << "{l0,l1,l2} = {" << l0 << "," << l1 << "," << l2 << "}";
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
std::tuple<bool, networkTetra *> genTetra(Network *const net,
                                          netP *const p0,
                                          netP *const p1,
                                          netP *p2,
                                          netP *p3) {
   // 0          3          0
   // *--- 2 ----*---- 2 ---*
   //  \       /   \       /
   //   1  2  4     5  3  0
   //    \   /   0   \   /
   //     \ /         \ /
   //     2*---- 3 ----*1
   //       \         /
   //        1   1   0
   //         \     /
   //          \   /
   //            * 0

   /*
    *このような四面体は作ってはならない
    *      *   *
    *     / \ /  \
    *    /  /\    \
    *   / /   \ ---*
    *  *-------*
    */
   // if (Dot(Cross(ToX(p2) - ToX(p1), ToX(p3) - ToX(p1)), (ToX(p1) + ToX(p2) + ToX(p3)) / 3. - ToX(p0)) < 0)
   //    std::swap(p2, p3);

   if (!DuplicateFreeQ(T_4P{p0, p1, p2, p3}))
      return {false, nullptr};

   if (Dot(TriangleNormal(p1->X, p2->X, p3->X), Mean(p1->X, p2->X, p3->X) - p0->X) < 0)
      std::swap(p2, p3);

   auto [l3, l4, l5] = link(p1, p2, p3, net);
   auto f0 = genFace(net, l3, l4, l5);
   auto [l1, l3_, l0] = link(p0, p2, p1, net);
   auto f1 = genFace(net, l1, l3_, l0);
   auto [l2, l4_, l1_] = link(p0, p3, p2, net);
   auto f2 = genFace(net, l2, l4_, l1_);
   auto [l0_, l5_, l2_] = link(p0, p1, p3, net);
   auto f3 = genFace(net, l0_, l5_, l2_);
   if (l0 == l0_ || l1 == l1_ || l2 == l2_ || l3 == l3_ || l4 == l4_ || l5 == l5_) {
      auto [t0, t1] = f0->Tetras;
      if (t0 && Intersection(t0->Points, T_4P{p0, p1, p2, p3}).size() == 4)
         return {false, t0};
      else if (t1 && Intersection(t1->Points, T_4P{p0, p1, p2, p3}).size() == 4)
         return {false, t1};
      else {
         auto tet = new networkTetra(net, T_4P{p0, p1, p2, p3}, T_6L{l0, l1, l2, l3, l4, l5}, T_4F{f0, f1, f2, f3});
         // std::cout << "Total(tet->solidangles/M_PI) = " << Total(tet->solidangles) / (4. * M_PI) << std::endl;
         // std::cout << "Mean(tet->solidangles/M_PI) = " << Mean(tet->solidangles) / (4. * M_PI) << std::endl;
         // std::cout << "tet->solidangles = " << tet->solidangles / (4. * M_PI) << std::endl;
         return {true, tet};
      }
   } else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "contradictions");
};

std::tuple<bool, networkTetra *> genTetra(Network *const net, const T_4P &abcd) {
   auto [a, b, c, d] = abcd;
   return genTetra(net, a, b, c, d);
};

Tddd Nearest(const Tddd &X, const networkFace *f) {
   return Nearest(X, ToX(f));
   // auto [a, b, c] = ToX(f);
   // Tddd m = (a + b + c) / 3.;
   // Tddd ab = (a + b) / 2., bc = (b + c) / 2., ca = (c + a) / 2.;
   // auto ret = Nearest(X, ToX(f));
   // auto X1 = Nearest(X, T3Tddd{a, ab, ca});
   // auto X2 = Nearest(X, T3Tddd{b, bc, ab});
   // auto X3 = Nearest(X, T3Tddd{c, ca, bc});
   // auto X4 = Nearest(X, T3Tddd{ab, bc, ca});
   // if (Norm(ret - X) > Norm(X1 - X))
   //    ret = X1;
   // if (Norm(ret - X) > Norm(X2 - X))
   //    ret = X2;
   // if (Norm(ret - X) > Norm(X3 - X))
   //    ret = X3;
   // if (Norm(ret - X) > Norm(X4 - X))
   //    ret = X4;
   // return ret;
};
double Distance(const Tddd &X, const networkFace *f) { return Norm(Nearest(X, ToX(f)) - X); };
// Tddd Nearest(const networkPoint *p, const networkFace *f) { return Nearest(ToX(p), ToX(f)); };
std::tuple<Tddd, networkFace *> Nearest_(const Tddd &X, const std::unordered_set<networkFace *> &faces) {
   double nearest_d = 1E+20, d;
   Tddd nearest_x, x;
   networkFace *nearest_f = nullptr;
   std::ranges::for_each(faces, [&](const auto &f) {
      x = Nearest(X, f);
      if (nearest_d >= (d = Norm(X - x))) {
         nearest_d = d;
         nearest_x = x;
         nearest_f = f;
      }
   });
   return {nearest_x, nearest_f};
};
std::tuple<Tddd, networkFace *> Nearest_(const Tddd &X, const std::vector<networkFace *> &faces) {
   double distance = 1E+20, tmp;
   Tddd ret, near;
   networkFace *F = nullptr;
   for (const auto &f : faces)
      if (distance > (tmp = Norm(X - (near = Nearest(X, f))))) {
         distance = tmp;
         ret = near;
         F = f;
      }
   return {ret, F};
};
std::tuple<Tddd, networkFace *> Nearest_(const Tddd &X, const double &r, const std::unordered_set<networkFace *> &faces) {
   double distance = 1E+20, tmp;
   Tddd ret, near;
   networkFace *F = nullptr;
   for (const auto &f : faces)
      if (distance > (tmp = Norm(X - (near = Nearest(X, f))))) {
         distance = tmp;
         ret = near;
         F = f;
      }
   return {ret, F};
};
Tddd Nearest(const Tddd &X, const std::vector<networkFace *> &faces) { return std::get<0>(Nearest_(X, faces)); };
Tddd Nearest(const Tddd &X, const std::unordered_set<networkFace *> &faces) { return std::get<0>(Nearest_(X, faces)); };
Tddd Nearest(const networkPoint *p, const std::unordered_set<networkFace *> &faces) { return Nearest(ToX(p), faces); };

//% ========================================================================== */
//%  Networkは持っているPointsから情報を取り出すメソッドを提供する*/
//%  @-@-@
//%  |\|/|
//%  @-@-@
//%  |/|\|
//%  @-@-@
//% ========================================================================== */
#include "RigidBodyDynamics.hpp"
#include "vtkWriter.hpp"
class MooringLine;
class Network : public CoordinateBounds, public RigidBodyDynamics {
  public:
   virtual ~Network();  // = default;  // Virtual destructor

   vtkPolygonWriter<networkPoint *> vtpWriter;

  public:
   Network *surfaceNet;

   // b# -------------------------------------------------------------------------- */
   // b#                                 octreeに関する                               */
   // b# -------------------------------------------------------------------------- */
  public:
   octree<networkFace *> *octreeOfFaces;
   octree<networkPoint *> *octreeOfPoints;
   void genOctreeOfFaces(const Tii &depthlimit, const int objnum) {
      std::cout << this->getName() << "genOctreeOfFaces(" << depthlimit << ", " << objnum << ")";
      this->setGeometricProperties();
      if (this->octreeOfFaces)
         delete this->octreeOfFaces;
      this->octreeOfFaces = new octree(this->scaledBounds(expand_bounds), depthlimit, objnum, this->Faces);
      Print(this->getName() + "->octreeOfFaces->deleteOuside()");
      this->octreeOfFaces->deleteOuside();
      Print(this->getName() + "->octreeOfFaces->setNeighbors()");
      this->octreeOfFaces->setNeighbors();
   };

   void genOctreeOfPoints(const Tii &depthlimit, const int objnum) {
      Print("octreeOfPointsを生成");
      this->setGeometricProperties();
      if (this->octreeOfPoints)
         delete this->octreeOfPoints;
      this->octreeOfPoints = new octree(this->scaledBounds(expand_bounds), depthlimit, objnum, this->Points);
      this->octreeOfPoints->deleteOuside();
   };
   /* -------------------------------------------------------------------------- */
   std::tuple<bool, octree<networkFace *> *, networkFace *> isInside_MethodOctree(const Tddd &X, const double small = 1E-13) const {
      //! セルに含まれるかどうかではなく，ポリゴンの内部にあるかどうかで判定するための関数
      if (this->octreeOfFaces) {
         auto cells = this->octreeOfFaces->getIntersectInside(X);
         if (cells.empty())
            return {false, nullptr, nullptr};  //! どのセルにもXは入らない
         else {
            auto cell = cells[0];
            if (cell->faces_.empty()) {
               //% cellがfaceを持たない場合，
               //% まずはセルの頂点における回転数を調べることで，
               //% ポリゴンの内部外部のどちらに位置するか調べる
               if (std::ranges::all_of(cell->WNs, [](const auto &w_num) { return w_num > 0.6; }) /*should be inside*/)
                  return {true, cell, nullptr};
               else if (std::ranges::all_of(cell->WNs, [](const auto &w_num) { return w_num < 0.1; }) /*should be outside*/)
                  return {false, cell, nullptr};
               else {
                  //@ 構造物の内部外部に位置するか回転数で判断するのが困難な場合
                  std::unordered_set<networkFace *> faces;
                  for (const auto &N0 : cell->neighbors)
                     for (const auto &nei : N0->neighbors)
                        faces.insert(nei->faces_.begin(), nei->faces_.end());
                  if (faces.empty())
                     return {false, cell, nullptr};
                  else {
                     auto [nearest, f] = Nearest_(X, faces);
                     return {(Dot(nearest - X, f->normal) >= 0. || Norm(nearest - X) < small), cell, f};
                  }
               }
            } else {
               //% faceを持つ場合，最も近い面とdotが同じかどうか
               auto [nearest, f] = Nearest_(X, cell->faces_);
               return {(Dot(nearest - X, f->normal) >= 0. || Norm(nearest - X) < small), cell, f};
            }
         }
      } else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "octreeOfFaces has not been set");
   };

   Tddd interpolateVector(const Tddd &X) const {
      auto [isinside, cell, _] = isInside_MethodOctree(X);
      if (cell)
         return cell->Interpolate(X, cell->vectors);
      else
         return {0., 0., 0.};
   };
   // b# -------------------------------------------------------------------------- */
  public:
   RungeKutta<Tddd> RK_COM;
   RungeKutta<T4d> RK_Q;
   RungeKutta<T6d> RK_Velocity;

  public:
   bool IGNORE;
   // Tdd grid_pull_factor;
   int grid_pull_depth;
   /* ------------------------------------------------------ */
   /*                          体積の計算                      */
   /* ------------------------------------------------------ */
   double getVolume() const {
      // ガウスの定理において，F=(x,y,z)とおいてdivF=3とすると，体積積分の結果は体積の3倍となる．
      // 面積分側は(x*nx+y*ny+z*nz)を核にした面積分と体積の3倍が等しいことになる．
      double ret = 0;
      for (const auto &f : getFaces()) {
         auto intp = interpolationTriangleLinear0101(f->getXVertices());
         for (const auto &[x0, x1, w0w1] : __GWGW5__Tuple)
            ret += w0w1 * Dot(intp(x0, x1), intp.cross(x0, x1));
      }
      return ret / 3.;
   };

   double GaussIntegral2(const Tddd &X) const {
      /*
      integrate(r.n/|r^3|)

       ^^^^
       ||||
      +----+
      | 4pi|-->   outside:0
      +----+
      */
      double ret = 0;
      Tddd r;
      // for (const auto &f : this->getFaces())
      // {
      // 	interpolationTriangleLinear0101 intp(f->getXVertices());
      // 	for (const auto &[x0, x1, w0w1] : __GWGW14__Tuple)
      // 	{
      // 		r = intp(x0, x1) - X;
      // 		ret += Dot(r / std::pow(Norm(r), 3), intp.cross(x0, x1)) * w0w1;
      // 	}
      // }

      for (const auto &f : this->getFaces()) {

         auto [X0, X1, X2] = f->getXVertices();
         X0 -= X2;
         X1 -= X2;
         auto A = X - X2;
         double tmp = 0;
         for (const auto &[t0, t1, w0w1] : __GWGW14__Tuple)
            tmp += (1. - t0) / std::pow(Norm(X0 * t0 + X1 * t1 * (1 - t0) - A), 3) * w0w1;
         ret += Dot(A, Cross(X0, X1)) * tmp;
      }
      return ret;

      // for (const auto &f : this->getFaces())
      // {

      // 	auto [X0, X1, X2] = f->getXVertices();
      // 	X0 -= X2;
      // 	X1 -= X2;
      // 	auto A = X - X2;
      // 	double tmp = 0;
      // 	for (const auto &[t0, t1, w0w1] : __GWGW8__Tuple)
      // 		tmp += 1. / std::pow(Norm(X0 * t0 * std::pow(1 - t0, -1 / 3.) + X1 * t1 * std::pow(1 - t0, 2 / 3.) - A * std::pow(1 - t0, -1 / 3.)), 3) * w0w1;
      // 	ret += Dot(A, Cross(X0, X1)) * tmp;
      // }
      // return ret;
   };

   double GaussIntegral(const Tddd &X) const {
      /*
      integrate(r.n/|r^3|)

       ^^^^
       ||||
      +----+
      | 4pi|-->   outside:0
      +----+
      */
      double ret = 0;
      Tddd r;
      for (const auto &f : this->getFaces()) {
         interpolationTriangleLinear0101 intp(f->getXVertices());
         for (const auto &[x0, x1, w0w1] : __GWGW14__Tuple) {
            r = intp(x0, x1) - X;
            ret += Dot(r / std::pow(Norm(r), 3), intp.cross(x0, x1)) * w0w1;
         }
      }
      return ret;
   };

   double windingNumber(const Tddd &X) const {
      double ret = 0;
      for (const auto &f : this->getFaces()) {
         ret += SolidAngle_VanOosteromAandStrackeeJ1983(X, f->getXVertices());
         // ret += SolidAngle(X, f->getXVertices());
      }
      return ret / (4. * M_PI);
   };

   bool isInside(const Tddd &X) const {
      if (!CoordinateBounds::isInside(X))
         return false;
      else if (this->windingNumber(X) < 0.75)
         return false;
      else
         return true;
   };

   std::vector<double> windingNumber(const std::vector<Tddd> &Xs) const {
      std::vector<double> ret(Xs.size(), 0.);
      T3Tddd V;
      for (const auto &f : this->getFaces()) {
         V = f->getXVertices();
         for (auto i = 0; i < Xs.size(); ++i)
            ret[i] += SolidAngle_VanOosteromAandStrackeeJ1983(Xs[i], V);
      }
      for (auto &r : ret)
         r /= (4. * M_PI);
      return ret;
   };

   bool all_of_isInside(const std::vector<Tddd> &Xs) const {
      if (!CoordinateBounds::isInside(X))
         return false;
      else {
         for (auto &wn : this->windingNumber(Xs)) {
            // ひとつでも小さい値があればfalse
            if (wn < 0.75)
               return false;
         }
         return true;
      }
   };
   //! ------------------------------------------------------ */
   //!                          接触の判別                      */
   //! ------------------------------------------------------ */

   std::unordered_set<networkFace *> getContactFacesOfPoints() const {
      std::unordered_set<networkFace *> ret;
      for (const auto &p : this->getPoints())
         ret.insert(std::begin(p->getContactFaces()), std::end(p->getContactFaces()));
      return ret;
   };

   std::unordered_set<networkPoint *> getContactPointsOfPoints() const {
      std::unordered_set<networkPoint *> ret;
      for (const auto &p : this->getPoints())
         ret.insert(begin(p->getContactPoints()), end(p->getContactPoints()));
      return ret;
   };
   std::unordered_set<networkPoint *> getContactPointsOfPoints(const std::vector<Network *> &nets) const {
      /*
      getContactPointsOfPointsは，自身の保有するPointsが接した点を返す．
      Pointsの保有するmap_Net_ContactPointsから指定されたNetworkのPointsを抽出している．
      */
      std::unordered_set<networkPoint *> ret;
      for (const auto &p : this->getPoints())
         ret.insert(begin(p->getContactPoints(nets)), end(p->getContactPoints(nets)));
      return ret;
   };

  public:
   //% ------------------------------------------------------ */
   void makeMirroredPoints(const Buckets<networkFace *> &B_face, const double mirroring_distance) {
      auto points = this->getPoints();
      for (const auto &p : points)
         p->makeMirroredPoints(B_face, mirroring_distance);
   };
   void clearMirroredPoints() {
      auto points = this->getPoints();
      for (const auto &p : points)
         p->clearMirroredPoints();
   };
   //% ------------------------------------------------------ */

   // b$ =================================================================== */
   // b$       ３次元空間分割に関する．　バケツ．Faces,Points,ParametricPoints          */
   // b$ =================================================================== */

  public:
   Buckets<networkFace *> BucketFaces;
   Buckets<networkPoint *> BucketParametricPoints;
   Buckets<networkPoint *> BucketPoints;

   const Buckets<networkFace *> &getBucketFaces() const { return BucketFaces; };
   const Buckets<networkPoint *> &getBucketPoints() const { return BucketPoints; };
   const Buckets<networkPoint *> &getBucketParametricPoints() const { return BucketParametricPoints; };

   const double expand_bounds = 1.5;

   /*
   面は点と違って，複数のバケツ（セル）と接することがある．
    */

   void makeBucketFaces(const double spacing = 1E+20) {
      this->setGeometricProperties();
      this->BucketFaces.clear();  // こうしたら良くなった
      std::cout << this->getName() << "this->scaledBounds(expand_bounds) = " << this->scaledBounds(expand_bounds) << std::endl;
      this->BucketFaces.initialize(this->scaledBounds(expand_bounds), spacing);
      double particlize_spacing;
      for (const auto &f : this->getFaces()) {
         particlize_spacing = Min(Tdd{spacing / 4., Min(extLength(f->getLines())) / 4.}) / 10.;
         for (const auto [xyz, t0t1] : triangleIntoPoints(f->getXVertices(), particlize_spacing))
            this->BucketFaces.add(xyz, f);
         // add extra points to make sure that the face is included in the bucket
         for (const auto X : f->getXVertices())
            this->BucketFaces.add(X, f);
         this->BucketFaces.add(Mean(f->getXVertices()), f);
      }

      this->BucketFaces.setVector();
      std::cout << this->getName() << ", BucketFaces.setVector() done" << std::endl;
   };

   double last_makeBucketPoints_spacing = 0.;
   void makeBucketPoints(const double spacing) {
      this->last_makeBucketPoints_spacing = spacing;
      // this->BucketPoints.hashing_done = false;
      this->setGeometricProperties();
      this->BucketPoints.clear();
      this->BucketPoints.initialize(this->scaledBounds(expand_bounds), spacing);
      std::cout << this->getName() << ", resize done" << std::endl;
      if (!this->BucketPoints.add(this->getPoints()))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "points are not added");
      else {
         std::cout << green << this->getName() << ", all points are added" << colorReset << std::endl;
         std::cout << green << "BucketPoints.all_stored_objects.size() = " << this->BucketPoints.all_stored_objects.size() << colorReset << std::endl;
      }
      this->BucketPoints.setVector();
      std::cout << this->getName() << ", BucketPoints.setVector() done" << std::endl;
   };

   void remakeBucketPoints() {
      this->makeBucketPoints(this->last_makeBucketPoints_spacing);
   };

   void makeBucketParametricPoints(const double spacing) {
      // this->BucketPoints.hashing_done = false;
      this->setGeometricProperties();
      this->BucketParametricPoints.initialize(this->scaledBounds(expand_bounds), spacing);
      this->BucketParametricPoints.add(this->getParametricPoints());
      std::cout << "this->getParametricPoints().size()=" << this->getParametricPoints().size() << std::endl;
   };
   void makeBucketParametricPoints(const T3Tdd &bounds, const double spacing) {
      // this->BucketPoints.hashing_done = false;
      this->BucketParametricPoints.initialize(this->scaledBounds(expand_bounds), spacing);
      this->BucketParametricPoints.add(this->getParametricPoints());
   };
   void clearBucketParametricPoints() {
      // this->BucketPoints.hashing_done = false;
      this->BucketParametricPoints.clear();
   };

   //! セルに含まれるかどうかではなく，ポリゴンの内部にあるかどうかで判定するための関数

   bool isInside_MethodBucket(const Tddd &X) const {
      networkFace *closest_face = nullptr;
      Tddd V_to_closest = {1E+20, 1E+20, 1E+20}, V = {0., 0., 0.};
      double min_d = 1E+20, d;
      this->BucketFaces.apply_to_the_nearest_bound(X, [&](const int i, const int j, const int k) {
         if (closest_face == nullptr)
            for (const auto &f : this->BucketFaces.data[i][j][k]) {
               V = Nearest(X, f) - X;
               d = Norm(V);
               if (d < min_d) {
                  min_d = d;
                  closest_face = f;
                  V_to_closest = V;
               }
            }
      });
      if (closest_face == nullptr)
         return false;
      else
         return Dot(closest_face->normal, V_to_closest) >= 0.;
   };

   // b$ ------------------------------------------------------ */

   // @ ============================================================== */
   // @                 RigidBodyDynamicsを使ったメソッド                 */
   // @ ============================================================== */

   void RigidBodyUpdatePoints() {
      for (const auto &p : this->getPoints())
         p->setXSingle(rigidTransformation(this->ICOM, this->center_of_mass, this->Q.R(), p->initialX));
      // p->setXSingle(this->quaternion.Rv(p->initialX - this->ICOM) - (p->initialX - this->ICOM) + trans + p->initialX);
      this->setGeometricProperties();
   };

   void calcPhysicalProperties() {
      this->force.fill(0.);
      this->inertia.fill(0.);
      this->mass = 0.;
      for (const auto &p : this->Points) {
         // for (auto i = 0; i < 3; ++i)
         // 	this->force[i] += p->force[i];
         std::get<0>(this->force) += std::get<0>(p->force);
         std::get<1>(this->force) += std::get<1>(p->force);
         std::get<2>(this->force) += std::get<2>(p->force);
         this->mass += p->mass;  // total mass
      }
      std::get<0>(this->inertia) = std::get<1>(this->inertia) = std::get<2>(this->inertia) = this->mass;  // total mass

      // this->center_of_mass = {0, 0, 0};
      // //! -------------------------------------- */
      // for (const auto &p : this->Points)
      // 	this->center_of_mass += p->X;
      // this->center_of_mass /= this->mass;
      // //! -------------------------------------- */
      Tddd tmp = {0, 0, 0};
      for (const auto &p : this->Points) {
         tmp = p->X - this->center_of_mass;
         std::get<3>(this->force) += std::get<0>(p->force) * std::get<0>(tmp);
         std::get<4>(this->force) += std::get<1>(p->force) * std::get<1>(tmp);
         std::get<5>(this->force) += std::get<2>(p->force) * std::get<2>(tmp);
         //
         std::get<3>(this->inertia) += p->mass * std::pow(std::get<0>(tmp), 2.);
         std::get<4>(this->inertia) += p->mass * std::pow(std::get<1>(tmp), 2.);
         std::get<5>(this->inertia) += p->mass * std::pow(std::get<2>(tmp), 2.);
      }
   };

   //* ------------------- 物体全体にかかる力の計算方法２つ ------------------- */

   // まずは，点にかかる力を計算してから使うこと．
   void sumForceOfPoints() {
      //%this->center_of_massが設定してある必要がある
      //%ネットワークの各点にかかっている力の総和を計算する
      this->force.fill(0.);
      for (const auto &p : this->Points) {
         auto F = p->force;
         auto r = (p->X - this->center_of_mass);
         auto [Fx, Fy, Fz, Fx_, Fy_, Fz_] = F;
         auto [Nx, Ny, Nz] = Cross(r, Tddd{Fx, Fy, Fz});  // モーメント
         std::get<0>(this->force) += Fx;
         std::get<1>(this->force) += Fy;
         std::get<2>(this->force) += Fz;
         std::get<3>(this->force) += Nx;
         std::get<4>(this->force) += Ny;
         std::get<5>(this->force) += Nz;
      }
   };
   void integrateForceOnFace() {
      /*(*checked by mathematica*)
      f0 = {f0x, f0y, f0z}
      f1 = {f1x, f1y, f1z}
      f2 = {f2x, f2y, f2z}
      x0 = {x0x, x0y, x0z}
      x1 = {x1x, x1y, x1z}
      x2 = {x2x, x2y, x2z}
      f[x0_, x1_] := Dot[{f0, f1, f2}, {x0, x1, 1 - x0 - x1}]
      X[z0_, z1_] := Dot[{x0, x1, x2}, {z0, z1, 1 - z0 - z1}]
      CForm@Flatten[
      Join[Simplify@
              Integrate[f[z0, z1], {\[Xi]0, 0, 1}, {z1, 0, 1 - z0}],
      Simplify@
              Integrate[
              Cross[X[z0, z1], f[z0, z1]], {z0, 0, 1}, {z1, 0, 1 - z0}]]]
      */
      this->force.fill(0.);
      for (const auto &f : this->Faces) {
         auto [p0, p1, p2] = f->getPoints();
         auto [f0x, f0y, f0z, f0x_, f0y_, f0z_] = p0->force;
         auto [f1x, f1y, f1z, f1x_, f1y_, f1z_] = p1->force;
         auto [f2x, f2y, f2z, f2x_, f2y_, f2z_] = p2->force;
         auto [x0x, x0y, x0z] = p0->X - this->center_of_mass;
         auto [x1x, x1y, x1z] = p1->X - this->center_of_mass;
         auto [x2x, x2y, x2z] = p2->X - this->center_of_mass;
         this->force += T6d{(f0x + f0y + f0z) / 6.,
                            (f1x + f1y + f1z) / 6.,
                            (f2x + f2y + f2z) / 6.,
                            (f2z * x1x + f2z * x1y + 2 * f2z * x1z + f2x * (2 * x1x + x1y + x1z) + f2y * (x1x + 2 * x1y + x1z) - 2 * f1x * x2x - f1y * x2x - f1z * x2x - f1x * x2y - 2 * f1y * x2y - f1z * x2y -
                             f1x * x2z - f1y * x2z - 2 * f1z * x2z) /
                                24.,
                            (-(f2z * x0x) - f2z * x0y - 2 * f2z * x0z - f2x * (2 * x0x + x0y + x0z) - f2y * (x0x + 2 * x0y + x0z) + 2 * f0x * x2x + f0y * x2x +
                             f0z * x2x + f0x * x2y + 2 * f0y * x2y + f0z * x2y + f0x * x2z + f0y * x2z + 2 * f0z * x2z) /
                                24.,
                            (f1z * x0x + f1z * x0y + 2 * f1z * x0z + f1x * (2 * x0x + x0y + x0z) + f1y * (x0x + 2 * x0y + x0z) - 2 * f0x * x1x - f0y * x1x - f0z * x1x - f0x * x1y - 2 * f0y * x1y - f0z * x1y -
                             f0x * x1z - f0y * x1z - 2 * f0z * x1z) /
                                24.};
      }
   };
   //    void calcAccelFromForce() {
   //       //! 慣性を設定しておく必要がある．
   //       this->acceleration = this->force / this->inertia;
   //    };

   //@ ------------------------------------------------------ */

   // BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM
   /*Network_BEM_detail
   境界要素法用のNetworkには，以下が設定されている．

   * bool BASE;
   * bool CORNER;
   * bool Dirichlet;
   * bool Neumann;

   `networkPoint`も`this->getNetwork()->setB()`などとして，自身の境界条件判断できる．
   Network_BEM_detail*/
   /*Network_BEM_code*/
  public:
   bool isFluid;
   bool isRigidBody;
   bool isSoftBody;
   bool isFloatingBody;
   bool isAbsorber;
#ifdef BEM
  public:
   double _current_time_;
   //
   std::array<bool, 6> isFixed = {false, false, false, false, false, false};
   JSON inputJSON;
   std::tuple<std::string, double> velocity_name_start;

   const double move_amplitude = 0.4;
   T6d velocityPredefined() {
      double t = _current_time_;
      double s = M_PI / 2.;
      // double a = move_amplitude;
      double k = M_PI / 1.;
      /* ------------------------------------------------------ */
      // T6d move_dir = {std::cos(k * t), std::sin(k * t), 0., 0., 0., 0.};
      // T6d ddt_move_dir = {-k * std::sin(k * t), k * std::std::cos(k * t), 0., 0., 0., 0.};
      // // /* |U|*n_p . n_surface = phin <-- given
      // auto tmp = (-move_amplitude * std::exp(-t) * (sin(k * t - s) - std::sin(-s)) + move_amplitude * std::exp(-t) * (std::cos(k * t - s) * k)) * move_dir;
      // tmp += move_amplitude * std::exp(-t) * (sin(k * t - s) - std::sin(-s)) * ddt_move_dir;
      // return tmp;
      /* ------------------------------------------------------ */
      Tddd tmp = Normalize(Tddd{1., 1., 0.});
      T6d move_dir = {std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), 0., 0., 0.};
      return (-move_amplitude * std::exp(-t) * (std::sin(k * t - s) - std::sin(-s)) + move_amplitude * std::exp(-t) * (std::cos(k * t - s) * k)) * move_dir;
   };

   Tddd translationPredefined() {
      double t = _current_time_;
      double s = M_PI / 2.;
      double k = M_PI / 1.;
      /* ------------------------------------------------------ */
      // Tddd move_dir = {std::cos(k * t), std::sin(k * t), 0.};
      // return move_amplitude * std::exp(-t) * (sin(k * t - s) - std::sin(-s)) * move_dir;
      /* ------------------------------------------------------ */
      Tddd move_dir = Normalize(Tddd{1., 1., 0.});
      return move_amplitude * std::exp(-t) * (std::sin(k * t - s) - std::sin(-s)) * move_dir;
   };
   V_d nabla_phi;
   V_d meanX;
#endif
   /*Network_BEM_code*/
   // BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM

  private:
   using V_netLp = std::vector<networkLine *>;
   std::string name = "no_name";
   // 接続している点をLinesをとおして素早く参照できる
  protected:
   std::unordered_set<networkTetra *> Tetras;
   std::unordered_set<networkPoint *> Points;
   std::unordered_set<networkFace *> Faces;
   std::vector<networkPoint *> Points_vector;
   std::vector<networkFace *> Faces_vector;
   V_netPp PointGarbage;
   V_netFp FaceGarbage;

  public:
   std::unordered_set<networkLine *> Lines;
   // 2021/08/10追加
   // bool isMember(const networkPoint *const p_IN) const { return (this->Points.find(p_IN) != this->Points.end()); };
   bool MemberQ(networkPoint *const &p_IN) const { return (this->Points.find(p_IN) != this->Points.end()); };
   // bool isMember(const networkFace *const f_IN) const { return (this->Faces.find(f_IN) != this->Faces.end()); };
   bool MemberQ(networkFace *const &f_IN) const { return (this->Faces.find(f_IN) != this->Faces.end()); };
   //
   void add(netPp const p_IN) { this->Points.emplace(p_IN); };
   void add(netFp const f_IN) { this->Faces.emplace(f_IN); };
   void add(netTp const t_IN) { this->Tetras.emplace(t_IN); };
   void add(const V_netPp &ps_IN) {
      for (const auto &p : ps_IN)
         this->add(p);
   };
   //
   void add(const V_netFp &fs_IN) {
      for (const auto &f : fs_IN)
         this->add(f);
   };
   //
   bool Erase(networkFace *f) {
      auto it = this->Faces.find(f);
      if (it != this->Faces.end()) {
         this->Faces.erase(it);
         return true;
      } else
         return false;
   };
   void Erase(const V_netFp &fs) {
      for (auto &f : fs)
         Erase(f);
   };
   bool Erase(networkPoint *p) {
      auto it = this->Points.find(p);
      if (it != this->Points.end()) {
         this->Points.erase(it);
         return true;
      } else
         return false;
   };

   //% ------------------------- obj ------------------------ */
   //% NetwrokObjは，Networkとしてコンストラクトするために，下のコンストラクタように修正した．*/
   // \label{Network::constructor}
   // Netwrokコンストラク
   Network(const std::string &filename = "file_name_is_not_given", const std::string &name_IN = "no_name")
       : CoordinateBounds(Tddd{{0., 0., 0.}}),
         RigidBodyDynamics(),
         name(name_IN),
         BucketFaces(CoordinateBounds(Tddd{{0., 0., 0.}}), 1.),
         BucketPoints(CoordinateBounds(Tddd{{0., 0., 0.}}), 1.),
         BucketParametricPoints(CoordinateBounds(Tddd{{0., 0., 0.}}), 1.),
         IGNORE(false),
         grid_pull_depth(0),
         velocity_name_start({"fixed", 0.}),
         inputJSON(),
         octreeOfFaces(nullptr),
         octreeOfPoints(nullptr),
         surfaceNet(nullptr) {
      if (filename.contains(".obj") || filename.contains(".off")) {
         Load3DFile objLoader(filename);
         if (!objLoader.f_v.empty())
            this->setFaces(objLoader.f_v, this->setPoints(objLoader.v));  // indexの書き換えも可能だがする必要は今のところない
         if (!objLoader.l_v.empty())
            this->setLines(objLoader.l_v, this->setPoints(objLoader.v));  // indexの書き換えも可能だがする必要は今のところない
         this->displayStates();
      } else {
         std::cout << filename << std::endl;
      }
   };
   //% ------------------------------------------------------ */
   const std::string &getName() const { return this->name; };
   void setName(const std::string nameIN) { this->name = nameIN; };
   V_netPp linkXPoints(Network &water, Network &obj);
   /* ------------------------------------------------------ */
   void displayStates() {
      std::cout << "/* -------------------------------------------------------------------------- */" << colorReset << std::endl;
      const int size = 15;
      V_i num(size, 0), pointsLines(size, 0);
      // 点の持つ線の数をカウント
      for (const auto &p : this->Points)
         for (auto i = 0; i < size; i++)
            if (p->Lines.size() == (size_t)i)
               pointsLines[i]++;
      V_i connection(3);
      for (const auto &l : this->getLines()) {
         if ((l->getFaces()).size() == 0)
            connection[0]++;
         else if ((l->getFaces()).size() == 1)
            connection[1]++;
         else if ((l->getFaces()).size() == 2)
            connection[2]++;
      }
      std::cout << Magenta << std::setw(25) << "this->getName(): " << this->getName() << ", " << this << std::endl;
      std::cout << Magenta << std::setw(25) << "   Lines.size(): " << getLines().size() << std::endl;
      std::cout << Magenta << std::setw(25) << "  Points.size(): " << Points.size() << std::endl;
      std::cout << Magenta << std::setw(25) << "   Faces.size(): " << Faces.size() << std::endl;
      std::cout << Magenta << std::setw(25) << "   this->bounds: " << this->bounds << std::endl;
      std::cout << Green << std::setw(25) << "                  " << GridVector(Subdivide(0., (double)size - 1., size - 1), 6) << std::endl;
      std::cout << Magenta << std::setw(25) << "Lines of points : " << GridVector(pointsLines, 6) << std::endl;
      std::cout << Magenta << std::setw(25) << " Connection of faces : " << GridVector(connection, 6);
      if (connection[0] == 0 && connection[1] == 0 && connection[2] != 0)
         std::cout << Blue << " 全ての線が2つの面と接続しているため，格子は閉じた面を形成している" << colorReset;
      std::cout << colorReset << std::endl;
      std::cout << "/* -------------------------------------------------------------------------- */" << colorReset << std::endl;
   };
   double rotation_offset;
   V_d translation_offset;
   double rotation;
   V_d translation;
   void rotate(const double theta, const Tddd &axis) {
      for (auto &p : this->getPoints())
         p->X = Dot(RotationMatrix(theta, axis), p->X);
      this->setGeometricProperties();
   };
   /* ------------------------------------------------------ */
   /*             ネットワークの姿勢に関する変数とメソッド          */
   /* ------------------------------------------------------ */
   void rotate(const Quaternion &Q, const Tddd &barycenter = {0., 0., 0.}) {
      for (auto &p : this->getPoints())
         p->X = Q.Rv(p->X - barycenter) + barycenter;
      setGeometricProperties();
   };
   void scale(const double ratio, const Tddd &center = {0, 0, 0}) {
      for (auto &p : this->getPoints())
         p->setX(center + (p->X - center) * ratio);
      this->setGeometricProperties();
   };
   void scale(const Tddd &ratio, const Tddd &center = {0, 0, 0}) {
      for (auto &p : this->getPoints())
         p->setX(center + (p->X - center) * ratio);
      this->setGeometricProperties();
   };
   void translate(const Tddd &translation) {
      for (auto &p : this->getPoints())
         p->setXSingle(p->X + translation);
      this->setGeometricProperties();
   };
   void translateFromInitialX(const Tddd &translation) {
      for (auto &p : this->getPoints())
         p->setXSingle(p->initialX + translation);
      this->setGeometricProperties();
   };
   void resetInitialX() {
      for (auto &p : this->getPoints())
         p->initialX = p->X;
   };
   void reverseNormal() {
      for (auto &f : this->Faces)
         f->reverseNormal();
   };
   /* -------------------------------------------------------------------------- */
   void applyTransformations(const JSON &json) {
      if (json.find("center_of_mass")) {
         std::get<0>(this->center_of_mass) = stob(json()["center_of_mass"])[0];
         std::get<1>(this->center_of_mass) = stob(json()["center_of_mass"])[1];
         std::get<2>(this->center_of_mass) = stob(json()["center_of_mass"])[2];
      }
      if (json.find("ignore")) {
         this->IGNORE = stob(json.at("ignore"))[0];
      }
      if (json.find("rotate")) {
         auto rotate = stod(json()["rotate"]);
         if (rotate.size() > 1)
            this->rotate(rotate[0], Tddd{rotate[1], rotate[2], rotate[3]});
      }
      if (json.find("scale")) {
         auto scale = stod(json()["scale"]);
         if (scale.size() > 1)
            this->scale({scale[0], scale[1], scale[2]});
         else
            this->scale(scale[0]);
      }
      if (json.find("translate")) {
         auto translate = stod(json.at("translate"));
         if (translate.size() > 1)
            this->translate({translate[0], translate[1], translate[2]});
         resetInitialX();
      }

      if (json.find("radius"))
         for (auto &p : this->getPoints())
            p->radius = stod(json.at("radius"))[0];

      if (json.find("mass")) {
         if (json.at("mass").size() == 1) {
            std::get<2>(this->inertia) = std::get<1>(this->inertia) = std::get<0>(this->inertia) = this->mass = stod(json.at("mass"))[0];
         }
         if (json.at("mass").size() == 3) {
            std::get<0>(this->inertia) = stod(json.at("mass"))[0];
            std::get<1>(this->inertia) = stod(json.at("mass"))[1];
            std::get<2>(this->inertia) = stod(json.at("mass"))[2];
         }
      }
      if (json.find("MOI")) {
         std::get<3>(this->inertia) = stod(json.at("MOI"))[0];
         std::get<4>(this->inertia) = stod(json.at("MOI"))[1];
         std::get<5>(this->inertia) = stod(json.at("MOI"))[2];
      }
      if (json.find("COM"))
         this->COM = this->initial_center_of_mass = Tddd{stod(json.at("COM"))[0], stod(json.at("COM"))[1], stod(json.at("COM"))[2]};

      // if (json.find("remesh"))
      // {
      // 	auto minlen = stod(json()["remesh"]);
      // 	if (minlen.size() > 0)
      // 		remesh(&this, minlen[0]);
      // }
      // if (json.find("coarsen"))
      // {
      // 	auto minlen = stod(json()["coarsen"]);
      // 	if (minlen.size() > 0)
      // 		coarsen(&this, minlen[0]);
      // }
      if (json.find("reverseNormal")) {
         std::string TorF = json()["reverseNormal"][0];
         if (TorF.compare("True") == 0 || TorF.compare("true") == 0 || TorF.compare("1") == 0) {
            this->reverseNormal();
            std::cout << "reverse done" << std::endl;
         }
      }
   };
   /* -------------------------------------------------------------------------- */
   networkFace *getNearestFace(const Tddd &xyz) const {
      networkFace *ret = nullptr;
      double nearest_dist = 1E+30, dist;
      for (const auto &f : this->Faces) {
         if (nearest_dist > (dist = Norm(f->X - xyz))) {
            nearest_dist = dist;
            ret = f;
         }
      }
      return ret;
   };
   // これを使おう．
   networkPoint *getNearestPoint(const Tddd &xyz) const {
      networkPoint *ret = nullptr;
      double nearest_dist = 1E+30, dist;
      for (const auto &p : this->Points) {
         if (nearest_dist > (dist = Norm(p->X - xyz))) {
            nearest_dist = dist;
            ret = p;
         }
      }
      return ret;
   };

   const std::unordered_set<networkPoint *> &getPoints() const { return this->Points; };

   /* ------------------------------ 境界条件に関する設定関数 ------------------------------ */

   std::unordered_set<networkPoint *> getPointsCORNER() const {
      std::unordered_set<networkPoint *> ret;
      ret.reserve(this->Points.size());
      for (const auto &p : this->Points)
         if (p->CORNER)
            ret.emplace(p);
      return ret;
   };
   std::unordered_set<networkPoint *> getPointsNeumann() const {
      std::unordered_set<networkPoint *> ret;
      ret.reserve(this->Points.size());
      for (const auto &p : this->Points)
         if (p->Neumann)
            ret.emplace(p);
      return ret;
   };
   std::unordered_set<networkPoint *> getPointsDirichlet() const {
      std::unordered_set<networkPoint *> ret;
      ret.reserve(this->Points.size());
      for (const auto &p : this->Points)
         if (p->Dirichlet)
            ret.emplace(p);
      return ret;
   };

   void setMinDepthFromCORNER() {
      std::vector<networkPoint *> points = ToVector(this->getPoints());
      for (const auto &p : points) {
         if (p->CORNER)
            p->minDepthFromCORNER_ = p->minDepthFromCORNER = 0;
         else
            p->minDepthFromCORNER_ = p->minDepthFromCORNER = 100000000;
         //
         if (p->isMultipleNode)
            p->minDepthFromMultipleNode_ = p->minDepthFromMultipleNode = 0;
         else
            p->minDepthFromMultipleNode_ = p->minDepthFromMultipleNode = 100000000;
      }

      for (auto i = 0; i < 100; i++) {
#pragma omp parallel
         for (auto &p : points)
#pragma omp single nowait
            for (auto &q : p->getNeighbors()) {
               p->minDepthFromCORNER_ = std::min(p->minDepthFromCORNER_, q->minDepthFromCORNER + 1);
               p->minDepthFromMultipleNode_ = std::min(p->minDepthFromMultipleNode_, q->minDepthFromMultipleNode + 1);
            }
         // apply
#pragma omp parallel
         for (const auto &p : points)
#pragma omp single nowait
         {
            p->minDepthFromCORNER = p->minDepthFromCORNER_;
            p->minDepthFromMultipleNode = p->minDepthFromMultipleNode_;
         }
      }
   };

   /* -------------------------------------------------------------------------- */

   const std::unordered_set<networkFace *> &getFaces() const { return this->Faces; };
   //
   const std::vector<networkFace *> &getFacesVector() const { return this->Faces_vector; }
   const std::vector<networkPoint *> &getPointsVector() const { return this->Points_vector; }
   // Update Faces_vector based on the contents of Faces
   void setFacesVector() {
      this->Faces_vector.clear();
      this->Faces_vector.reserve(this->Faces.size());
      this->Faces_vector.assign(this->Faces.begin(), this->Faces.end());
   }
   // Update Points_vector based on the contents of Points
   void setPointsVector() {
      this->Points_vector.clear();
      this->Points_vector.reserve(this->Points.size());
      this->Points_vector.assign(this->Points.begin(), this->Points.end());
   }

   const std::unordered_set<networkTetra *> &getTetras() const { return this->Tetras; };
   std::unordered_set<networkPoint *> getParametricPoints() const {
      std::unordered_set<networkPoint *> ret, tmp;
      for (const auto &f : this->getFaces()) {
         tmp = f->getParametricPoints();
         ret.insert(tmp.begin(), tmp.end());
      }
      return ret;
   };
   Tddd getMeanX() const {
      Tddd ret = {0, 0, 0};
      for (const auto &p : this->Points)
         ret += p->X;
      return ret / this->Points.size();
   };
   // VV_d getLocations() const
   // {
   // 	VV_d ret;
   // 	getLocations(ret);
   // 	return ret;
   // };
   void getLocations(VV_d &ret) const {
      ret.resize(Points.size());
      int i(0);
      // for (const auto &p : Points)
      // 	ret[i++] = p->xyz;
      for (const auto &p : Points)
         ret[i++] = {std::get<0>(p->X), std::get<1>(p->X), std::get<2>(p->X)};
   };
   // std::vector<Tddd> getLocationsTuple() const
   // {
   // 	std::vector<Tddd> ret(this->Points.size());
   // 	int i = 0;
   // 	for (const auto &p : this->Points)
   // 		ret[i++] = p->X;
   // 	return ret;
   // };
   //-------------------------
   // get lines
   void setLinesStatus(bool TorF) {
      for (const auto &p : Points)
         p->setLinesStatus(TorF);
   };

   std::unordered_set<networkLine *> getLines() const {
      std::unordered_set<networkLine *> ret;
      ret.reserve(this->Points.size() * 3);
      for (auto p : this->Points)
         for (auto l : p->getLines())
            ret.emplace(l);
      return ret;
   };

   std::unordered_set<networkLine *> getLinesUO() const {
      std::unordered_set<networkLine *> ret;
      ret.reserve(3 * this->Faces.size());
      for (const auto &f : this->Faces) {
         auto [l0, l1, l2] = f->getLines();
         ret.emplace(l0);
         ret.emplace(l1);
         ret.emplace(l2);
      }
      return ret;
   };

   V_netLp getLinesIntxn() const {
      V_netLp ret(0);
      for (const auto &p : this->getPoints())
         for (const auto &l : p->Lines)
            if (l->isIntxn())
               ret.emplace_back(l);
      return DeleteDuplicates(ret);
   };
   //======================================
   /*setXPoints_detail
   `setXPoints`は，thisとtargetの`networkFace`と`networkLine`の干渉を検知して，干渉点として`networkPoint`を新たに作成する．

   - 干渉点`networkPoint`は，干渉した`networkFace`と`networkLine`の`XPoints`に`push`される．
   - 干渉点として作成される`networkPoint`には，`xline`と`xface`が与えられる．これらは，デフォルトで`xline=nullptr`，`xface=nullptr`となっている．

   setXPoints_detail*/

   void setXPoints(Network &target, Network *interactionNet);

   /**
   p->setFaces()もf->setFaces()も,自身のp->Lines,f->Linesを元に，p->Faces,f->Facesを決定し保存する．
   p->Linesとf->Linesが正しく設定してあるかチェックする．
   特に，flipやdivideの後には，p->Lines,f->Linesが正しく設定されていない可能性があるので，要注意．
   */

   void setGeometricProperties() {
      if (!this->getPoints().empty()) {
         for (const auto &p : this->getPoints())
            p->setFaces();  // % point->Lines must have been determined
         for (const auto &l : this->getLines())
            l->setBoundsSingle();
         for (const auto &f : this->getFaces())
            f->setGeometricProperties(ToX(f->setPoints()));  // @ face->Lines must have been determined
         for (const auto &f : this->getFaces())
            f->setDodecaPoints();
         CoordinateBounds::setBounds(ToX(this->getPoints()));
      } else {
         CoordinateBounds::setBounds(Tddd{{0., 0., 0.}});
      }
   };

   /**
   pointのFacesとLinesの関係は整合性があるか？faceのFacesとLinesの関係は整合性があるか？をチェック．
   setGeometricProperties()を実行していれば，この整合性は保たれるはずではある．
   */

   std::array<bool, 2> validateConsistency() {

      bool consistent_Line_Face_in_Points = true;
      bool consistent_Line_Face_in_Faces = true;

      for (const auto &p : this->getPoints()) {
         std::unordered_set<networkFace *> tmp1, tmp2;
         for (const auto &f : p->getFacesFromLines())
            tmp1.emplace(f);
         for (const auto &f : p->getFaces())
            tmp2.emplace(f);
         auto s = Intersection(tmp1, tmp2).size();
         if (s != tmp1.size() || s != tmp2.size())
            consistent_Line_Face_in_Points = false;
      }

      for (const auto &f : this->getFaces()) {
         std::unordered_set<networkPoint *> tmp1, tmp2;
         for (const auto &p : f->getPoints())
            tmp1.emplace(p);
         for (const auto &p : f->getPoints())
            tmp2.emplace(p);
         auto s = Intersection(tmp1, tmp2).size();
         if (s != tmp1.size() || s != tmp2.size())
            consistent_Line_Face_in_Faces = false;
      }

      return {consistent_Line_Face_in_Points,
              consistent_Line_Face_in_Faces};
   };

   /* -------------------------------------------------------------------------- */

  public:
   V_netPp setPoints(const std::vector<Tddd> &v_IN) {
      std::unordered_map<std::shared_ptr<Tddd>, int> map_Tddd_Int;
      std::vector<Tddd> vecTddd(v_IN.size());
      map_Tddd_Int.reserve(v_IN.size());
      int n = 0;
      for (const auto &v : v_IN) {
         std::shared_ptr<Tddd> x(new Tddd(v));
         vecTddd[n] = *x;
         map_Tddd_Int[x] = n++;
      }
      CoordinateBounds bound(MinMaxColumns(vecTddd));
      Buckets<std::shared_ptr<Tddd>> bucket(bound, bound.getScale() / 100.);
      std::vector<networkPoint *> ret(v_IN.size());
      int overlap_index = 0, overlaps = 0;
      for (const auto &[x, n] : map_Tddd_Int) {
         auto [i, j, k] = bucket.indices(*x);
         overlap_index = 0;
         for (const auto &tdd : bucket.data[i][j][k])
            if (Norm(*tdd - *x) < 1E-10) {
               overlap_index = map_Tddd_Int[tdd];
               break;
            }
         if (!overlap_index) {
            bucket.add(*x, x);
            ret[n] = new networkPoint(this, *x);
         } else {
            ret[n] = ret[overlap_index];
            overlaps++;
         }
      }
      return ret;
   };

   V_netPp setPoints(const VV_d &v_IN) {
      std::unordered_map<std::shared_ptr<Tddd>, int> map_Tddd_Int;
      std::vector<Tddd> vecTddd(v_IN.size());
      map_Tddd_Int.reserve(v_IN.size());
      int n = 0;
      for (const auto &v : v_IN) {
         if (v.size() != 3) throw std::runtime_error("v.size() != 3");
         std::shared_ptr<Tddd> x(new Tddd({v[0], v[1], v[2]}));
         vecTddd[n] = *x;
         map_Tddd_Int[x] = n++;
      }
      CoordinateBounds bound(MinMaxColumns(vecTddd));
      Buckets<std::shared_ptr<Tddd>> bucket(bound, bound.getScale() / 100.);
      std::vector<networkPoint *> ret(v_IN.size());
      int overlap_index = 0, overlaps = 0;
      for (const auto &[x, n] : map_Tddd_Int) {
         auto [i, j, k] = bucket.indices(*x);
         overlap_index = 0;
         for (const auto &tdd : bucket.data[i][j][k])
            if (Norm(*tdd - *x) < 1E-10) {
               overlap_index = map_Tddd_Int[tdd];
               break;
            }
         if (!overlap_index) {
            bucket.add(*x, x);
            ret[n] = new networkPoint(this, *x);
         } else {
            ret[n] = ret[overlap_index];
            overlaps++;
         }
      }
      // std::cout << "overlaped nodes : " << overlaps << std::endl;
      // std::cout << "setPoints elapsed time :" << timer() << std::endl;
      return ret;
   };
   /* -------------------------------------------------------------------------- */
  public:
   void setFaces(const std::vector<T3Tddd> &v_IN, const double resolution = 1E-10) {
      std::vector<std::array<networkPoint *, 3>> v_Points(v_IN.size());

      CoordinateBounds bound(v_IN);
      Buckets<networkPoint *> bucket;
      bucket.initialize(bound, bound.getScale() / 20.);

      auto findOverlap = [&](const Tddd &X) -> networkPoint * {
         for (const auto &p : bucket.getData(X))
            if (Norm(p->X - X) < resolution)
               return p;
         return nullptr;
      };

      networkPoint *p = nullptr;
      int i = 0;
      for (const auto &position : v_IN) {
         int j = 0;
         for (const auto &X : position) {
            auto p = findOverlap(X);
            if (p == nullptr) {
               p = v_Points[i][j] = new networkPoint(this, X);
               bucket.add(X, p);
            } else
               v_Points[i][j] = p;
            j++;
         }
         i++;
      }

      std::cout << "this->Points.size() : " << this->Points.size() << std::endl;

      for (const auto &points : v_Points)
         new networkFace(this, points);

      this->setGeometricProperties();

      // check size
      std::cout << "v_IN.size() : " << v_IN.size() << std::endl;
      std::cout << "this->Faces.size() : " << this->Faces.size() << std::endl;
   };

   /* -------------------------------------------------------------------------- */

   void setFaces(const VV_i &f_v, const V_netPp &points) {
      try {
         int count = 0;
         for (const auto &index : f_v) {
            if (points[index[0]] && points[index[1]] && points[index[2]])
               new networkFace(this, {link(points[index[0]], points[index[1]], this), link(points[index[1]], points[index[2]], this), link(points[index[2]], points[index[0]], this)}, {points[index[0]], points[index[1]], points[index[2]]});
            else {
               std::stringstream ss;
               ss << "index = " << index << std::endl;
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
            }
         }
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
      this->setGeometricProperties();
   };

   /* -------------------------------------------------------------------------- */

   void setLines(const VV_i &l_v, const V_netPp &points) {
      std::cout << "l_v " << l_v << std::endl;
      try {
         for (const auto &index : l_v) {
            for (auto i = 0; i < index.size() - 1; i++) {
               if (points[index[i]] && points[index[i + 1]])
                  new networkLine(this, points[index[i]], points[index[i + 1]]);
               else {
                  std::stringstream ss;
                  ss << "index = " << index << std::endl;
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
               }
            }
         }
      } catch (std::exception &e) {
         std::cerr << e.what() << colorReset << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
      this->setGeometricProperties();
   };

   /* -------------------------------------------------------------------------- */

  public:
   V_netPp getXPoints() {
      V_netPp ret({});
      for (const auto &line : this->getLines())
         for (const auto &info : line->XPoints)
            ret.emplace_back(info);
      return ret;
   };
   void clearXPoints() {
      for (const auto &p : this->Points)
         for (const auto &line : p->Lines)
            line->XPoints.clear();

      for (const auto &face : this->Faces) {
         face->clearXPoints();
         // for (const auto &line : face->getLines())
         // 	line->XPoints.clear();

         std::ranges::for_each(face->getLines(), [&](const auto &l) { l->XPoints.clear(); });
      }
   };
   /* ------------------------------------------------------ */

   void DeleteParticles() {
      auto points = this->Points;
      for (const auto &p : points)
         if (std::get<0>(p->particlize_info))
            p->Delete();
   };

   /* ------------------------------------------------------ */
  public:
   std::vector<MooringLine *> mooringLines;
};
//% ========================================================================== */
//% ========================================================================== */
//% ========================================================================== */
networkFace *isFace(const networkPoint *const a, const networkPoint *const b, const networkPoint *const c) {
   for (const auto &p : {a, b, c})
      for (const auto &f : p->Faces)
         if (f->MemberQ(a) && f->MemberQ(b) && f->MemberQ(c))
            return f;
   return nullptr;
};
networkLine *ConnectedQ(const networkFace *const a, const networkFace *const b) {
   auto [l0, l1, l2] = a->getLines();
   if ((*l0)(a) == b)
      return l0;
   else if ((*l1)(a) == b)
      return l1;
   else if ((*l2)(a) == b)
      return l2;
   return nullptr;
};
networkLine *ConnectedQ(const networkLine *const l, const networkPoint *const a) {
   // 両方が参照し合う状態かどうか
   /*
    __line__
    |      |
    |     *|-><--* point
    --------
   */
   for (const auto &L : a->getLines())
      if (L == l) {
         auto [p0, p1] = L->getPoints();
         if (p0 == a || p1 == a)
            return L;
      }
   return nullptr;
};
networkLine *ConnectedQ(const networkPoint *const a, const networkLine *const l) { return ConnectedQ(l, a); };

networkLine *ConnectedQ(const networkPoint *const a, const networkPoint *const b) {
   // 両方が参照し合う状態かどうか
   /*
    point  __line__
    *--><--|*     |
           |     *|-><--* point
           --------
   */
   networkLine *L;
   for (const auto &l : a->getLines())
      if (L = ConnectedQ(l, a))
         if (ConnectedQ(L, b))
            return L;
   return nullptr;
};

T_3L ConnectedQ(const networkPoint *const a, const networkPoint *const b, const networkPoint *const c) {
   return T_3L{ConnectedQ(a, b), ConnectedQ(b, c), ConnectedQ(c, a)};

   // for (const auto &la : a->getLines())
   //    if ((*la)(a) == b)
   //       for (const auto &lb : b->getLines())
   //          if ((*lb)(b) == c)
   //             for (const auto &lc : c->getLines())
   //                if ((*lc)(c) == a)
   //                   return T_3L{la, lb, lc};
   // return T_3L{nullptr, nullptr, nullptr};
};

T_3L ConnectedQ(const T_3P &abc) { return ConnectedQ(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc)); };

#include "my_vtk.hpp"
#include "networkPoint.hpp"

// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM
#ifdef BEM
V_Netp takeNetworks(const V_netPp &ps) {
   V_Netp ret({});
   for (auto &p : ps)
      network::add(ret, p->getNetwork());
   return ret;
};

inline V_d networkPoint::nabla_phi() const { return this->getNetwork()->nabla_phi; };

#endif
// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM

////////////////////////////////////////////////////////////////////////
inline V_netFp networkLine::getFaces(const bool TorF) const {
   V_netFp ret;
   for (const auto &f : this->getFaces())
      if (f->isStatus(TorF))
         ret.emplace_back(f);
   return ret;
};
////////////////////////////////////////////////////////
// int intersectQ(networkFace *f, networkLine *l)
// {
// 	return network::find(f->getLinesPenetrate(), l);
// };

inline V_netPp networkFace::getXPoints() const {  // 特別な場合のみしか使えないので改善が必要
   bool TorF = true;
   network::setStatus(this->network->getLines(), !TorF);
   auto Ps = getPointsPenetrate();
   V_netPp object{Ps[0]}, nextobject, ret{Ps[0]};
   netP *P;
   while (true) {
      for (const auto &p : object)
         for (const auto &l : p->getLines_toggle(!TorF)) {
            if ((P = (*l)(p))->getStatus() != TorF && P->getXFace() == this) {
               ret.emplace_back(P);
               P->setStatus(TorF);
               nextobject.emplace_back(P);
            }
         }
      if (nextobject.empty())
         break;
      object = nextobject;
      nextobject.clear();
   }
   ret.emplace_back(Ps[1]);
   return ret;
};
//@ ------------------------------------------------------ */
//@                          外部関数                       */
//@ ------------------------------------------------------ */
/*networkFace::Delete_detail
##### networkFace::Delete()
1. `this->Lines`から`this`を削除する
2. `this->Lines`の内，面を全て失った線は，点の接続を解除する
> **NOTE:**
`networkLine`から`networkFace`からへの一方的なポイントは，絶対に起きないようにしているつもりである．
**`networkFace`のコンストラクタでは，自動で`networkLine`に`networkFace`を格納するようにしている．**
なので，`networkFace`を削除する際に，空になった`networkLine`は削除して構わない，
これが可能なのは，一方向のみのポイントが存在しないからである．

networkFace::Delete_detail*/
/*networkFace::Delete_code*/
// #define debug_Delete

inline void networkFace::Delete() {
   /*!
   networkに上に穴が空いている状態もあり得るため，
   faceのDelete()は，pointのDelete()を伴わない．
   */
   try {
      if (!this->network->Erase(this))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "is not a member of " + this->network->getName());
      else
         std::ranges::for_each(this->Lines, [&](const auto &l) { l->Erase(this); });
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

/* -------------------------------------------------------------------------- */

/*!
pointのDelete()と同時に，pointが接続している辺のDelete()を必ず伴うべきだ．
*/
inline void networkPoint::Delete() {
   // std::cout << "networkPoint::Delete()\n";
   /*
   節点を消去する際に，隣接する面も消去すると，
   フリップや分割などの過程で面が意図せず消去される危険があるのでそのようなことはしない．
   */
   try {
      if (!this->network->Erase(this)) {
         std::stringstream ss;
         ss << "この点は，networkに含まれていない" << std::endl;
         ss << "network name : " << this->network->getName() << std::endl;
         if (this->network->MemberQ(this))
            ss << this->network->getName() << "に含まれている" << std::endl;
         else
            ss << this->network->getName() << "にも含まれていない" << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      } else {
         // std::cout << "networkPoint::Delete() 0\n";
         auto ps = this->getNeighbors();
         // std::cout << "networkPoint::Delete() 0.1\n";
         V_netLp delL = {};
         netLp line;
         for (const auto &p : ps) {
            // std::cout << "networkPoint::Delete() 1\n";
            if ((line = unlink(p, this)))
               delL.emplace_back(line);
         }
         // 完全にリンクしているものだけを削除する
         //
         // for (const auto &l : delL)
         //    if (l)
         //       for (const auto &f : l->getFaces())
         //          if (f)
         //             delete f;

         for (auto i = 0; i < delL.size(); i++) {
            // std::cout << "networkPoint::Delete() 2\n";
            if (delL[i]->getFaces().empty()) {
               auto tmp = delL[i];
               delL[i] = nullptr;
               delete tmp;
            }
         }
         // std::cout << "networkPoint::Delete() done\n";
         // this->resetXinfo();
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

/*!
lineのDelete()と同時に，lineが接続している面のDelete()を必ず伴うべきだ．
*/
inline void networkLine::Delete() {
   if (this->network->Lines.erase(this)) {
      if (this->Point_A != nullptr) {
         (this->Point_A)->Erase(this);
         this->Point_A = nullptr;
      }
      if (this->Point_B != nullptr) {
         (this->Point_B)->Erase(this);
         this->Point_B = nullptr;
      }
      // this->Faces.clear();
      // std::cout << "delete line " << this << std::endl;
   }
};
// inline void networkLine::Delete() {
//    this->network->Lines.erase(this);
//    try {
//       int maxloop = 1000 + this->Faces.size();
//       int counter = 0;
//       while (!this->Faces.empty() && counter++ < maxloop) {
//          auto tmp = this->Faces;
//          tmp[0]->Delete();
//       }

//       if (this->Point_A != nullptr) {
//          (this->Point_A)->Erase(this);
//          this->Point_A = nullptr;
//       }
//       if (this->Point_B != nullptr) {
//          (this->Point_B)->Erase(this);
//          this->Point_B = nullptr;
//       }

//       auto tmpP = this->getXPoints();
//       for (auto &p : tmpP) {
//          p->resetXinfo();
//          network::erase(this->XPoints, p);
//       }
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorReset << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };

///////////////////////////////////////////////////////////////
inline networkPoint::~networkPoint() { this->Delete(); };
inline networkFace::~networkFace() { this->Delete(); };
inline networkLine::~networkLine() { this->Delete(); };
///////////////////////////////////////////////////////////////
inline Network::~Network() {
   if (octreeOfFaces)
      delete this->octreeOfFaces;
   //
   try {
      std::cout << this->getName() << "->Points.size() = " << this->Points.size() << std::endl;
      std::cout << this->getName() << "->Faces.size() = " << this->Faces.size() << std::endl;
      std::cout << this->getName() << "->Lines.size() = " << this->Lines.size() << std::endl;
      std::cout << "destroying Points in " << Red << this->getName() << colorReset << std::endl;
      {
         auto tmp = this->Points;
         for (const auto &p : tmp)
            delete p;
      }
      std::cout << "destroying Faces in " << Red << this->getName() << colorReset << std::endl;
      {
         auto tmp = this->Faces;
         for (const auto &f : tmp)
            delete f;
      }
      std::cout << "destroying Lines in " << Red << this->getName() << colorReset << std::endl;
      {
         auto tmp = this->Lines;
         for (const auto &l : tmp)
            delete l;
      }
      std::cout << this->getName() << "->Points.size() = " << this->Points.size() << std::endl;
      std::cout << this->getName() << "->Faces.size() = " << this->Faces.size() << std::endl;
      std::cout << this->getName() << "->Lines.size() = " << this->Lines.size() << std::endl;
      // try {
      //    int maxloop = 100000 + this->Faces.size();
      //    int counter = 0;
      //    while (!this->Faces.empty() && counter++ < maxloop) {
      //       auto tmp = this->Faces;
      //       (*tmp.begin())->Delete();
      //    }
      //    if (counter > maxloop) {
      //       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      //    }
      // } catch (std::exception &e) {
      //    std::cerr << e.what() << colorReset << std::endl;
      //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // };

      // std::cout << "|        |"
      //           << "Faces" << std::endl;
      // std::cout << "|        |" << this->Faces.size() << std::endl;
      // std::cout << "|________|" << std::endl;

      // //--------------------
      // try {
      //    int maxloop = 100000 + this->Points.size();
      //    int counter = 0;
      //    while (!this->Points.empty() && counter++ < maxloop) {
      //       //
      //       // auto tmp = this->Points;
      //       // tmp[0]->Delete();
      //       //
      //       // std::get<0>(this->Points)->Delete();
      //       (*this->Points.begin())->Delete();
      //    }
      //    if (counter > maxloop) {
      //       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      //    }
      // } catch (const error_message &e) {
      //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // }

      // std::cout << "|        |"
      //           << "Points" << std::endl;
      // std::cout << "|        |" << this->Points.size() << std::endl;
      // std::cout << "|________|" << std::endl;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorReset << std::endl;
   };
};

// Networkコンストラクタ2
// inline Network::Network(const V_Netp &base_nets, const std::string &name_IN = "no_name_xnetwork")
//     : CoordinateBounds({0, 0, 0}),
//       name(name_IN),
//       /* ------------------------------------------------------ */
//       force({0., 0., 0., 0., 0., 0.}),
//       inertia({0., 0., 0., 0., 0., 0.}),
//       forced_acceleration({0., 0., 0., 0., 0., 0.}),
//       forced_velocity({0., 0., 0., 0., 0., 0.}),
//       acceleration({0., 0., 0., 0., 0., 0.}),
//       velocity({0., 0., 0., 0., 0., 0.}),
//       mass(0.),
//       center_of_mass({0., 0., 0.}),
//       initial_center_of_mass({0., 0., 0.}),
//       BucketFaces(CoordinateBounds({0., 0., 0.}), 1.),
//       BucketPoints(CoordinateBounds({0., 0., 0.}), 1.),
//       BucketParametricPoints(CoordinateBounds({0., 0., 0.}), 1.),
//       octreeOfFaces(nullptr),
//       octreeOfPoints(nullptr) {
//    Timer timer;
// #ifdef BEM
//    isRigidBody = false;
//    isSoftBody = false;
// #endif
//    /* ----------------------------------- */
//    // 物理
//    // this->force = V_d(3, 0.);
//    /* ------------------------------------ */
//    // オブジェクトが交差している点を計算，保存
//    for (auto i = 0; i < base_nets.size(); i++)
//       for (auto j = 0; j < base_nets.size(); j++)
//          if (i != j) {
//             Print(base_nets[i]->getName() + base_nets[j]->getName());
//             base_nets[i]->setXPoints(*base_nets[j], this);
//          }

//    V_netPp accumXPoints({});
//    for (auto i = 0; i < base_nets.size(); i++)
//       for (auto j = i + 1; j < base_nets.size(); j++) {
//          auto tmp = linkXPoints(*base_nets[i], *base_nets[j]);
//          accumXPoints.insert(accumXPoints.end(), tmp.begin(), tmp.end());
//       }

//    /// CHECK IF INTERSECTION EXISTS
//    if (accumXPoints.empty() /*no intersection*/) {
//       message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Intersection Network has been generated but no intersection point exists");
//    } else {
//       this->add(accumXPoints);
//       this->setGeometricProperties();
//       Print("Intersection Network has been generated", Magenta);
//    }
//    displayStates();
//    std::cout << Green << "Elapsed time : " << timer() << colorReset << std::endl;
// };
////////////////////////////////////////////////////////////////////////
inline V_netPp Network::linkXPoints(Network &water, Network &obj) {
   V_netPp accumXPoints;
   std::cout << __PRETTY_FUNCTION__ << "Xpointをset " << water.getName() << " " << obj.getName() << std::endl;
   bool findsameface = false;
   std::vector<Network *> nets({&water, &obj});
   for (const auto &net : nets)
      for (const auto &f : net->Faces) {  // thisの干渉点だけを対象にする！！！
         auto xps = f->getPointsPenetrate(this /*干渉点の制限！*/);
         for (size_t i = 0; i < xps.size(); i++, findsameface = false)  // あるcrossinfoを持つLineにたいして
         {
            for (size_t j = i + 1; j < xps.size(); j++)
               if (xps[i]->getXFace() == xps[j]->getXFace()) {
                  // ２本のLineが，同じ面を貫通しているLineである場合
                  auto newline = link(xps[i], xps[j], this);
                  // newline->setIntxn(true);
                  accumXPoints.emplace_back(xps[i]);
                  accumXPoints.emplace_back(xps[j]);
                  findsameface = true;
                  //                break;
               }
            // ２本のLineが，異なる面を貫通しているLineである場合
            // それぞれが貫通している面のcrossinfoを調べる．その中に元の面を貫通するLineがあった場合，元の貫通点と接続する
            if (!findsameface)
               for (const auto &XPOfL : xps)
                  for (const auto &XPOfL2 : XPOfL->getXFace()->getPointsPenetrate(this /*干渉点の制限！*/))
                     if (XPOfL2->getXFace() == f) {
                        auto newline = link(XPOfL, XPOfL2, this);
                        // newline->setIntxn(true);
                        accumXPoints.emplace_back(XPOfL);
                        accumXPoints.emplace_back(XPOfL2);
                        //                    break;
                     }
         }
      }
   for (const auto &f : obj.Faces)  // thisの干渉点だけを対象にする！！！
   {
      auto xps = f->getPointsPenetrate(this /*干渉点の制限！*/);
      for (size_t i = 0; i < xps.size(); i++, findsameface = false)
         for (size_t j = i + 1; j < xps.size(); j++)
            if (xps[i]->getXFace() == xps[j]->getXFace()) {
               auto newline = link(xps[i], xps[j], this);
               // newline->setIntxn(true);
               accumXPoints.emplace_back(xps[i]);
               accumXPoints.emplace_back(xps[j]);
               findsameface = true;
               // break;
            }
   }
   return accumXPoints;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*XNetwork_code*/

/*networkFace::divide_detail
divideの操作を整理して理解するために，networkFaceにdivideの一部操作を任せている．
ここでnetworkFaceは，与えられた線と自身がすでに所持していた線を，適切に接続し直す．
点と線の接続はnetworkFaceに行わせず，netowrkLineに行わせる．
networkFace::divide_detail*/
//
/*networkFace::divide_code*/
inline netFp networkFace::divide(netLp DivL /*this*/, netLp newDivL, netLp newMidL, int type) {
   if (!MemberQ_(this->Lines, DivL))
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "divL is not found in this->Lines"));
   auto oldF = this;
   auto newF = new networkFace(this);
   auto fL = oldF->getLineFront(DivL);
   auto bL = oldF->getLineBack(DivL);
   //
   //          /     \                     /     \
   //         /       \                   /       \
   //        /         \                 /         \
   //       /    / \    \               /    / \    \
   // bL<--/--- /   \----\-->fL   bL<--/----/   \----\-->fL
   //   --/--->/oldF \<---\--         /    /newF \    \
   //    /     --|A---     \         /     --|----     \
   //   /        V|         \       /        V          \
   //    -------DivL-------         --------DivL---------
   //
   if (oldF->replace(fL, newMidL)) {  // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
      std::stringstream ss;
      // ss << "oldF->getLines() = " << oldF->getLines();
      ss << ", oldF = " << oldF << ", fL = " << fL << ", newMidL = " << newMidL;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   }

   // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
   if (newF->replace(bL, newMidL)) {
      std::stringstream ss;
      ss << "newF = " << newF << ", bL = " << bL << ", newMidL = " << newMidL;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   }
   // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
   if (fL->replace(oldF, newF))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

   //          /     \                     /     \
	//         /       \                   /       \
	//        /         \                 /         \
	//       /    / \    \               /    / \    \
	// bL<--/----/   \<---\-- newMidL<--/--  /   \----\-->fL
   //   --/--->/oldF \----\->         /    /newF \<---\---
   //    /     --|A---     \         /     --|----     \
	//   /        V|         \       /        V          \
	//    -------DivL-------         --------DivL---------
   //
   newMidL->set(oldF, newF);

   //          /     \                     /     \
	//         /       \                   /       \
	//        /         \                 /         \
	//       /    / \    \               /    / \    \
	// bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
   //  ---/--->/oldF \<---\--       --/--->/newF \<---\---
   //    /     --|A---     \         /     --|----     \
	//   /        V|         \       /        V          \
	//    -------DivL-------         --------DivL---------
   if (type == 1) {
      // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
      if (newF->replace(DivL, newDivL))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      newDivL->set(newF);
      //          /     \                     /     \
		//         /       \                   /       \
		//        /         \                 /         \
		//       /    / \    \               /    / \    \
		// bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
      //  ---/--->/oldF \<---\--       --/--->/newF \<---\---
      //    /     --|A---     \         /     --|A---     \
		//   /        V|         \       /        V|         \
		//    -------DivL/*this*/---     -------newDivL------
   } else if (type == 2) {
      // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
      if (oldF->Replace(DivL, newDivL))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      DivL->Add(newF);
      //          /     \                     /     \
		//         /       \                   /       \
		//        /         \                 /         \
		//       /    / \    \               /    / \    \
		// bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
      //  ---/--->/oldF \<---\--       --/--->/newF \<---\---
      //    /     --|A---     \         /     --|A---     \
		//   /        V|         \       /        V|         \
		//    ------newDivL------         ----DivL/*this*/----
   }
   return newF;
};

std::unordered_map<netP *, Tdd> &operator+=(std::unordered_map<netP *, Tdd> &ii_dd, const std::unordered_map<netP *, Tdd> &jj_dd) {
   std::unordered_map<netP *, Tdd>::iterator it;
   for (auto &[jj, dd] : jj_dd)
      if ((it = ii_dd.find(jj)) != ii_dd.end())
         it->second += dd;
      else
         ii_dd[jj] = dd;
   return ii_dd;
};

std::map<netP *, Tdd> &operator+=(std::map<netP *, Tdd> &ii_dd, const std::map<netP *, Tdd> &jj_dd) {
   std::map<netP *, Tdd>::iterator it;
   for (auto &[jj, dd] : jj_dd)
      if ((it = ii_dd.find(jj)) != ii_dd.end())
         it->second += dd;
      else
         ii_dd[jj] = dd;
   return ii_dd;
};

std::map<netP *, Tdd> &operator+=(std::map<netP *, Tdd> &ii_dd, const std::tuple<netP *, Tdd> &jj_dd) {
   std::map<netP *, Tdd>::iterator it;
   if ((it = ii_dd.find(std::get<0>(jj_dd))) != ii_dd.end())
      it->second += std::get<1>(jj_dd);
   else
      ii_dd[std::get<0>(jj_dd)] = std::get<1>(jj_dd);
   return ii_dd;
};

/* -------------------------------------------------------------------------- */

T4T3P ToT4T3P(const T_4P &p0123) {
   auto [p0, p1, p2, p3] = p0123;
   return {{{p0, p1, p2}, {p0, p1, p3}, {p0, p2, p3}, {p1, p2, p3}}};
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

// struct QuadPoints {
//    /*this nodes forms a quadratic element like below
//    0
//    | \
//    |   \
//    |     \
//    3(p2)-l1-5(p1)
//    | \    |   \
//    |  l2  l0   \
//    |   \  |     \
//    1 ----4(p0 --- 2
//    */
//    const std::array<networkPoint *, 3> p_of_face;
//    const std::array<networkLine *, 3> l_of_face;
//    const std::array<bool, 3> ignore;
//    const std::array<networkPoint *, 6> points;
//    QuadPoints(const networkPoint *origin, const networkFace *f)
//        : p_of_face(f->getPoints(origin)),
//          l_of_face(f->getLinesTupleFrom(std::get<0>(p_of_face))),
//          ignore({!std::get<1>(l_of_face)->CORNER,
//                  !std::get<2>(l_of_face)->CORNER,
//                  !std::get<0>(l_of_face)->CORNER}),
//          points({f->getPointOpposite(std::get<1>(l_of_face)) /*0*/,
//                  f->getPointOpposite(std::get<2>(l_of_face)) /*1*/,
//                  f->getPointOpposite(std::get<0>(l_of_face)) /*2*/,
//                  std::get<2>(p_of_face) /*3*/,
//                  std::get<0>(p_of_face) /*4*/,
//                  std::get<1>(p_of_face) /*5*/}){};
// };

struct QuadPoints {
   networkFace *f;
   const std::tuple<networkPoint *, networkLine *, networkPoint *, networkLine *, networkPoint *, networkLine *> PLPLPL;
   std::array<bool, 3> available;  // availability of the lines (face)
   std::array<networkFace *, 3> faces;
   std::array<networkPoint *, 6> points;
   std::function<bool(const networkLine *)> conditionChecker;
   std::array<std::tuple<networkPoint *, networkFace *>, 6> points_faces;

   /*       0
          /   \
         / f1  \
        /       \
    3(p2)--(l1)--5(p1)
      /  \       / \
     /  (l2) fc l0  \
    / f2   \   /  f0 \
   1------ 4(p0)------2

   points_facesは，以下のように，点と点が所属する面を保持する（多重節点を使う際に役立つ）．

   {{point 0, face 1},
   {point 1, face 2},
   {point 2, face 0},
   {point 3, face center},
   {point 4, face center},
   {point 5, face center}}

   与えられた点または線が，面のPLPLPLのメンバーなら，
   その点または線を先頭にするように，PLPLPLを設定し，保持する．

   */

   // QuadPoints(const networkPoint *origin_point, const networkFace *f) : f(f), PLPLPL(f->getPointsAndLinesFromNearest(origin_point)) { initialize(); };
   // QuadPoints(const networkLine *origin_line, const networkFace *f) : f(f), PLPLPL(f->getPointsAndLines(origin_line)) { initialize(); };

   QuadPoints(
       networkFace *f,
       const networkPoint *origin_point = nullptr,
       std::function<bool(const networkLine *)> checker = [](const networkLine *line) { return !line->CORNER && line->getFaces().size() == 2; })
       : f(f), PLPLPL(f->getPointsAndLinesFromNearest(origin_point)), conditionChecker(checker) { initialize(); };

   QuadPoints(
       networkFace *f,
       const networkLine *origin_line = nullptr,
       std::function<bool(const networkLine *)> checker = [](const networkLine *line) { return !line->CORNER && line->getFaces().size() == 2; })
       : f(f), PLPLPL(f->getPointsAndLines(origin_line)), conditionChecker(checker) { initialize(); };

   void initialize() {
      auto [p0, l0, p1, l1, p2, l2] = this->PLPLPL;

      this->available = {conditionChecker(l1), conditionChecker(l2), conditionChecker(l0)};

      auto fs0 = l0->getFaces();
      auto fs1 = l1->getFaces();
      auto fs2 = l2->getFaces();

      networkFace *f0, *f1, *f2;

      // if (fs0.size() == 2 && fs1.size() == 2 && fs2.size() == 2)

      if (fs0.empty() || fs1.empty() || fs2.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "fs0.empty() || fs1.empty() || fs2.empty()");

      {
         f0 = fs0.size() == 1 ? fs0[0] : (fs0[0] == f ? fs0[1] : fs0[0]);
         f1 = fs1.size() == 1 ? fs1[0] : (fs1[0] == f ? fs1[1] : fs1[0]);
         f2 = fs2.size() == 1 ? fs2[0] : (fs2[0] == f ? fs2[1] : fs2[0]);
      }
      // else throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "fs0.size() != 2 || fs1.size() != 2 || fs2.size() != 2");

      this->points = {f1->getPointOpposite(l1) /*0*/,
                      f2->getPointOpposite(l2) /*1*/,
                      f0->getPointOpposite(l0) /*2*/,
                      p2 /*3*/,
                      p0 /*4*/,
                      p1 /*5*/};

      this->faces = {f1 /*0*/, f2 /*1*/, f0};

      this->points_faces[0] = {this->points[0], f1};
      this->points_faces[1] = {this->points[1], f2};
      this->points_faces[2] = {this->points[2], f0};
      this->points_faces[3] = {this->points[3], f};
      this->points_faces[4] = {this->points[4], f};
      this->points_faces[5] = {this->points[5], f};
   };
};

/*
線形補間するとして，
pの重みがわかる必要がある．

線形三角補間の場合，割り振りは，{t0,t1,1-t0-t1}である．
2次三角補間の場合，割り振りは，{t0,t1,1-t0-t1,t0*t1,t0*(1-t0-t1),t1*(1-t0-t1)}である．


最後は，
{t0,t1,1-t0-t1}

{t0,t1,1-t0-t1,t0*t1,t0*(1-t0-t1),t1*(1-t0-t1)}

(t0)*
(t1)*
(1-t0-t1) ->

p2 = (Dot({-1, -1, 0, 2, 4, 4},{p0,p1,p2,p3,p4,p5}) + Dot({-1, -1, 0, 2, 4, 4},{p0,p1,p2,p3,p4,p5}))/16


今ある関数でうまくできる？

N = {-1, -1, 0, 2, 4, 4}/8

2次補間の重みを，さらに節点上の補間の重みを通して，節点上の補間の重みを求める．

t0 * p0の変数
t1 * p1の変数
(1-t0-t1) * p2の変数
[t0*t1] * (Dot(N,{p0,p1,p2,p3,p4,p5}) + Dot(N,{p0,p1,p2,p3,p4,p5}))/2
[t0*(1-t0-t1)] * (Dot(N,{p0,p1,p2,p3,p4,p5}) + Dot(N,{p0,p1,p2,p3,p4,p5}))/2
[t1*(1-t0-t1)] * (Dot(N,{p0,p1,p2,p3,p4,p5}) + Dot(N,{p0,p1,p2,p3,p4,p5}))/2

これをみると，
{t0, t1, 1-t0-t1, N * t0*t1, N * t0*(1-t0-t1), N * t1*(1-t0-t1)}
となっていることがわかる．N = {-1, -1, 0, 2, 4, 4}/8．

言い換えると，
t0は，p0の重み
t1は，p1の重み
1-t0-t1は，p2の重み
N * t0*t1は，p3の重みであるが，p3は存在しないので，これはさらに分配され，f_l0とfの節点に分配される．
N * t0*(1-t0-t1)は，p4の重みであるが，p4は存在しないので，これはさらに分配され，f_l1とfの節点に分配される．
N * t1*(1-t0-t1)は，p5の重みであるが，p5は存在しないので，これはさらに分配され，f_l2とfの節点に分配される．

*/
#include "minMaxOfFunctions.hpp"
struct DodecaPoints {

   networkFace *f;  // given center face
   const std::function<bool(const networkLine *)> conditionChecker;
   const QuadPoints quadpoint;
   const std::tuple<networkPoint *, networkLine *, networkPoint *, networkLine *, networkPoint *, networkLine *> PLPLPL;
   const networkLine *l0, *l1, *l2;
   const std::array<bool, 3> available;  // availability of the lines (face)
   networkFace *f0, *f1, *f2;

   const QuadPoints quadpoint_l0;
   const QuadPoints quadpoint_l1;
   const QuadPoints quadpoint_l2;

   const T3Tdd shape_map_center = {{{0., 0.5} /*quad 4 -> linear 0*/, {0.5, 0.} /*quad 5 -> linear 1*/, {0.5, 0.5} /*quad 3 -> linear 2*/}};
   //! これとN3の内積を新なパラメタとして利用すると，(t0,t1)=(1,0)で{0., 0.5}に，(t0,t1)=(0,1)で{0.5, 0.}に，t0=t1=0で{0.5, 0.5}になる．
   const T3Tdd shape_map_l0_face = {{{0.5, 0.} /*0*/, {0., 0.5} /*1*/, {0., 0.} /*2*/}};
   const T3Tdd shape_map_l1_face = {{{0., 0.} /*0*/, {0.5, 0.} /*1*/, {0., 0.5} /*2*/}};
   const T3Tdd shape_map_l2_face = {{{0., 0.5} /*0*/, {0., 0.} /*1*/, {0.5, 0.} /*2*/}};

   DodecaPoints(
       networkFace *f,
       const networkPoint *origin_point = nullptr,
       std::function<bool(const networkLine *)> checker = [](const networkLine *line) { return !line->CORNER && line->getFaces().size() == 2; })
       : f(f),
         conditionChecker(checker),
         quadpoint(f, origin_point, checker),
         PLPLPL(quadpoint.PLPLPL),
         l0(std::get<1>(PLPLPL)),
         l1(std::get<3>(PLPLPL)),
         l2(std::get<5>(PLPLPL)),
         available(quadpoint.available),
         f0(l0->getFaces().size() == 1 ? l0->getFaces()[0] : (l0->getFaces()[0] == f ? l0->getFaces()[1] : l0->getFaces()[0])),
         f1(l1->getFaces().size() == 1 ? l1->getFaces()[0] : (l1->getFaces()[0] == f ? l1->getFaces()[1] : l1->getFaces()[0])),
         f2(l2->getFaces().size() == 1 ? l2->getFaces()[0] : (l2->getFaces()[0] == f ? l2->getFaces()[1] : l2->getFaces()[0])),
         quadpoint_l0(f0, l0, checker),
         quadpoint_l1(f1, l1, checker),
         quadpoint_l2(f2, l2, checker) {};
   /*
   これをみると，
   {t0, t1, 1-t0-t1, N * t0*t1, N * t0*(1-t0-t1), N * t1*(1-t0-t1)}
   となっていることがわかる．N = {-1, -1, 0, 2, 4, 4}/8．
   言い換えると，
   t0は，p0の重み
   t1は，p1の重み
   1-t0-t1は，p2の重み
   N * t0*t1は，p3の重みであるが，p3は存在しないので，これはさらに分配され，f_l0とfの節点に分配される．
   N * t0*(1-t0-t1)は，p4の重みであるが，p4は存在しないので，これはさらに分配され，f_l1とfの節点に分配される．
   N * t1*(1-t0-t1)は，p5の重みであるが，p5は存在しないので，これはさらに分配され，f_l2とfの節点に分配される．
   */

   /* -------------------------------------------------------------------------- */

   // template <typename T>
   // std::tuple<Tdd, T> findMinimumCommon(const std::function<T(networkPoint *)> &conversion) {

   //    Tdd first_guess = {0.25, 0.25};
   //    T first_point = this->interpolate(first_guess, conversion), tmp;

   //    if (first_point > (tmp = this->interpolate(0., 1., conversion))) {
   //       first_guess = {0., 1.};
   //       first_point = tmp;
   //    }
   //    if (first_point > (tmp = this->interpolate(1., 0., conversion))) {
   //       first_guess = {1., 0.};
   //       first_point = tmp;
   //    }
   //    if (first_point > (tmp = this->interpolate(0., 0., conversion))) {
   //       first_guess = {0., 0.};
   //       first_point = tmp;
   //    }

   //    const Tdd dt0t1 = {0., 1E-5};

   //    BroydenMethod<Tdd> BM(first_guess, first_guess + dt0t1);

   //    Tddd F, dFdt0, dFdt1;
   //    auto grad_FF = [&](const double t0, const double t1) -> Tdd {
   //       F = this->interpolate(t0, t1, conversion);
   //       dFdt0 = this->D_interpolate<1, 0>(t0, t1, conversion);
   //       dFdt1 = this->D_interpolate<0, 1>(t0, t1, conversion);
   //       return Tdd{Dot(dFdt0, F), Dot(dFdt1, F)};
   //    };

   //    for (auto i = 0; i < 50; i++) {
   //       BM.updateGoodBroyden(grad_FF(BM.X[0], BM.X[1]), grad_FF(BM.X[0] - BM.dX[0], BM.X[1] - BM.dX[1]));
   //       BM.X[0] = std::clamp(BM.X[0], 0., 1.);
   //       BM.X[1] = std::clamp(BM.X[1], 0., 1. - BM.X[0]);
   //       if (Norm(BM.dX) < 1e-13 && i > 3) {
   //          std::cout << "iteration : " << i << "\n";
   //          break;
   //       }
   //    }
   //    return {BM.X, this->interpolate(BM.X[0], BM.X[1], conversion)};
   // };

   /* -------------------------------------------------------------------------- */

   /*

   `findMininum`は，節点のデータから補間された関数上での最小値を求める関数である．
   なので，滑らかな関数として補間できない場合は，この関数は使えない．

   節点の空間座標を補間した表面上で，ある別の関数の最小値を探したい場合，
   対象となる'別の関数'の微分も用意しておく必要がある．

   */

   // #define USE_BROYDEN_METHOD_FOR_FIND_MINIMUM
   // #define USE_NEWTON_METHOD_FOR_FIND_MINIMUM
#define USE_SUBDIVISION_METHOD_FOR_FIND_MINIMUM

   template <std::size_t N>
   std::tuple<Tdd, std::array<double, N>> findMinimum(const std::function<std::array<double, N>(networkPoint *)> &conversion,
                                                      const auto &minimizing_function,
                                                      const std::array<double, 2> t0_range = {0., 1.}) {

      Tdd first_guess = {0.25, 0.25};
      std::array<double, N> first_point = this->interpolate(first_guess, conversion), tmp;

      // if (first_point > (tmp = this->interpolate(0., 1., conversion))) {
      //    first_guess = {0., 1.};
      //    first_point = tmp;
      // }
      // if (first_point > (tmp = this->interpolate(1., 0., conversion))) {
      //    first_guess = {1., 0.};
      //    first_point = tmp;
      // }
      // if (first_point > (tmp = this->interpolate(0., 0., conversion))) {
      //    first_guess = {0., 0.};
      //    first_point = tmp;
      // }

      std::array<double, N> F, dFdt0, dFdt1, dFdt0t0, dFdt0t1, dFdt1t1;

#ifdef USE_BROYDEN_METHOD_FOR_FIND_MINIMUM
      const Tdd dt0t1 = {0., 1E-10};
      auto grad_FF = [&](const double t0, const double t1) -> Tdd {
         F = this->interpolate(t0, t1, conversion);
         dFdt0 = this->D_interpolate<1, 0>(t0, t1, conversion);
         dFdt1 = this->D_interpolate<0, 1>(t0, t1, conversion);
         return Tdd{Dot(dFdt0, F), Dot(dFdt1, F)};
      };
      BroydenMethod<Tdd> BM(first_guess, first_guess + dt0t1);
      for (auto i = 0; i < 50; i++) {
         BM.updateGoodBroyden(grad_FF(BM.X[0], BM.X[1]), grad_FF(BM.X[0] - BM.dX[0], BM.X[1] - BM.dX[1]));
         BM.X[0] = std::clamp(BM.X[0], t0_range[0], t0_range[1]);
         BM.X[1] = std::clamp(BM.X[1], 0., 1. - BM.X[0]);
         if (Norm(BM.dX) < 1e-13 && i > 3) {
            // std::cout << "iteration : " << i << "\n";
            break;
         }
      }
      return {BM.X, this->interpolate(BM.X[0], BM.X[1], conversion)};
#elif defined(USE_NEWTON_METHOD_FOR_FIND_MINIMUM)

      NewtonRaphson<Tdd> NR(first_guess);
      for (auto i = 0; i < 50; i++) {
         auto [t0, t1] = NR.X;

         F = this->interpolate(t0, t1, conversion);
         dFdt0 = this->D_interpolate<1, 0>(t0, t1, conversion);
         dFdt1 = this->D_interpolate<0, 1>(t0, t1, conversion);
         Tdd grad = {Dot(dFdt0, F), Dot(dFdt1, F)};

         dFdt0t1 = this->D_interpolate<1, 1>(t0, t1, conversion);
         dFdt0t0 = this->D_interpolate<2, 0>(t0, t1, conversion);
         dFdt1t1 = this->D_interpolate<0, 2>(t0, t1, conversion);
         T2Tdd Hessian = {Tdd{Dot(dFdt0t0, F) + Dot(dFdt0, dFdt0), Dot(dFdt0t1, F) + Dot(dFdt1, dFdt0)},
                          Tdd{Dot(dFdt0t1, F) + Dot(dFdt0, dFdt1), Dot(dFdt1t1, F) + Dot(dFdt1, dFdt1)}};

         NR.update(grad, Hessian);

         NR.X[0] = std::clamp(NR.X[0], t0_range[0], t0_range[1]);
         NR.X[1] = std::clamp(NR.X[1], 0., 1. - NR.X[0]);
         if (Norm(NR.dX) < 1e-13 && i > 2)
            // std::cout << "iteration : " << i << "\n";
            break;
      }
      return {NR.X, this->interpolate(NR.X[0], NR.X[1], conversion)};

         /* -------------------------------------------------------------------------- */
#elif defined(USE_SUBDIVISION_METHOD_FOR_FIND_MINIMUM)
      double min = 1E10;
      Tdd t0t1_optimum = {0.25, 0.25};
      double t0_min = 0;
      double t0_max = 1;
      double t1_min = 0;
      double t1_max = 1;
      double div = 3, TMP;
      auto samples = convertions(conversion);
      double t0, t1;
      int count = 0;
      double eps = 1E-10;
      for (auto i = 0; i < 10; i++) {
         double dt0 = (t0_max - t0_min) / div;
         double dt1 = (t1_max - t1_min) / div;
         for (int j = 0; j <= div; j++) {
            t0 = t0_min + j * dt0;
            for (int k = 0; k <= div; k++) {  // Adjust based on your constraints
               t1 = t1_min + k * dt1;
               if (t0 + t1 > 1) break;
               TMP = minimizing_function(this->Dot4(N6_new(t0, t1), samples));
               if (min > TMP) {
                  min = TMP;
                  t0t1_optimum = {t0, t1};
               }
            }
         }

         t0_min = std::clamp(t0t1_optimum[0] - 0.5 * dt0, t0_range[0], t0_range[1] - eps);
         t0_max = std::clamp(t0t1_optimum[0] + 0.5 * dt0, t0_min, t0_range[1]);
         t1_min = std::clamp(t0t1_optimum[1] - 0.5 * dt1, 0., 1. - eps);
         t1_max = std::clamp(t0t1_optimum[1] + 0.5 * dt1, t1_min, 1.);
      }

      // std::cout << "count : " << count << "\n";
      return {t0t1_optimum, this->interpolate(t0t1_optimum[0], t0t1_optimum[1], conversion)};

#endif
   };

   std::tuple<Tdd, Tddd> Nearest(const Tddd &A, const Tdd &t0_range = {0., 1.}) {

      Tdd first_guess = {0.25, 0.25};
      Tddd first_point = this->X(first_guess), tmp;
      Tddd F, dFdt0, dFdt1, dFdt0t0, dFdt0t1, dFdt1t1;

      NewtonRaphson<Tdd> NR(first_guess);
      for (auto i = 0; i < 50; i++) {
         auto [t0, t1] = NR.X;
         F = this->X(t0, t1) - A;
         dFdt0 = this->D_X<1, 0>(t0, t1);
         dFdt1 = this->D_X<0, 1>(t0, t1);
         Tdd grad = {Dot(dFdt0, F), Dot(dFdt1, F)};

         dFdt0t1 = this->D_X<1, 1>(t0, t1);
         dFdt0t0 = this->D_X<2, 0>(t0, t1);
         dFdt1t1 = this->D_X<0, 2>(t0, t1);
         T2Tdd Hessian = {Tdd{Dot(dFdt0t0, F) + Dot(dFdt0, dFdt0), Dot(dFdt0t1, F) + Dot(dFdt1, dFdt0)},
                          Tdd{Dot(dFdt0t1, F) + Dot(dFdt0, dFdt1), Dot(dFdt1t1, F) + Dot(dFdt1, dFdt1)}};

         NR.update(grad, Hessian);

         NR.X[0] = std::clamp(NR.X[0], t0_range[0], t0_range[1]);
         NR.X[1] = std::clamp(NR.X[1], 0., 1. - NR.X[0]);
         if (Norm(NR.dX) < 1e-13 && i > 2)
            // std::cout << "iteration : " << i << "\n";
            break;
      }
      return {NR.X, this->X(NR.X[0], NR.X[1])};
   };

   /* -------------------------------------------------------------------------- */

   std::array<T6d, 4> N6_new(const double t0, const double t1) const {
      // auto N3 = TriShape<3>(t0, t1);
      auto N6_interface = TriShape<6>(t0, t1);

      T6d Nc = {0., 0., 0., N6_interface[2], N6_interface[0], N6_interface[1]};
      T6d Nl0, Nl1, Nl2;
      if (this->conditionChecker(this->l0)) {
         Nc += 0.5 * TriShape<6>(0.25, 0.25, this->quadpoint.available) * N6_interface[3];
         Nl0 = 0.5 * TriShape<6>(0.25, 0.25, this->quadpoint_l0.available) * N6_interface[3];
      } else {
         Nc += 0.5 * T6d{0., 0., 0., 0., 1., 1.} * N6_interface[3];
         Nl0.fill(0.);
      }

      if (this->conditionChecker(this->l1)) {
         Nc += 0.5 * TriShape<6>(0.5, 0.25, this->quadpoint.available) * N6_interface[4];
         Nl1 = 0.5 * TriShape<6>(0.25, 0.25, this->quadpoint_l1.available) * N6_interface[4];
      } else {
         Nc += 0.5 * T6d{0., 0., 0., 1., 0., 1.} * N6_interface[4];
         Nl1.fill(0.);
      }

      if (this->conditionChecker(this->l2)) {
         Nc += 0.5 * TriShape<6>(0.25, 0.5, this->quadpoint.available) * N6_interface[5];
         Nl2 = 0.5 * TriShape<6>(0.25, 0.25, this->quadpoint_l2.available) * N6_interface[5];
      } else {
         Nc += 0.5 * T6d{0., 0., 0., 1., 1., 0.} * N6_interface[5];
         Nl2.fill(0.);
      }

      // falseの倍は線形にしてはどうか？

      return {Nc, Nl0, Nl1, Nl2};
   };

   template <int i, int j>
   std::array<T6d, 4> N(const double t0, const double t1) const {
      return N6_new<i, j>(t0, t1);
   };

   template <int i, int j>
   std::array<T6d, 4> D_N6_new(const double t0, const double t1) const {
      // auto N3 = TriShape<3>(t0, t1);
      auto N6_interface = D_TriShape_Quadratic<i, j>(t0, t1);

      T6d Nc = {0., 0., 0., N6_interface[2], N6_interface[0], N6_interface[1]};
      T6d Nl0, Nl1, Nl2;
      if (this->conditionChecker(this->l0)) {
         Nc += 0.5 * TriShape<6>(0.25, 0.25, this->quadpoint.available) * N6_interface[3];
         Nl0 = 0.5 * TriShape<6>(0.25, 0.25, this->quadpoint_l0.available) * N6_interface[3];
      } else {
         Nc += T6d{0., 0., 0., 0., 0.5, 0.5} * N6_interface[3];
         Nl0.fill(0.);
      }

      if (this->conditionChecker(this->l1)) {
         Nc += 0.5 * TriShape<6>(0.5, 0.25, this->quadpoint.available) * N6_interface[4];
         Nl1 = 0.5 * TriShape<6>(0.25, 0.25, this->quadpoint_l1.available) * N6_interface[4];
      } else {
         Nc += T6d{0., 0., 0., 0.5, 0., 0.5} * N6_interface[4];
         Nl1.fill(0.);
      }

      if (this->conditionChecker(this->l2)) {
         Nc += 0.5 * TriShape<6>(0.25, 0.5, this->quadpoint.available) * N6_interface[5];
         Nl2 = 0.5 * TriShape<6>(0.25, 0.25, this->quadpoint_l2.available) * N6_interface[5];
      } else {
         Nc += T6d{0., 0., 0., 0.5, 0.5, 0.} * N6_interface[5];
         Nl2.fill(0.);
      }
      return {Nc, Nl0, Nl1, Nl2};
   };

   template <int i, int j>
   std::array<T6d, 4> D_N(const double t0, const double t1) const {
      return D_N6_new<i, j>(t0, t1);
   };

   std::array<double, 3> X(const double t0, const double t1) const {
      auto [Nc, N0, N1, N2] = N6_new(t0, t1);
      return Dot(Nc, ToX(this->quadpoint.points)) + Dot(N0, ToX(this->quadpoint_l0.points)) + Dot(N1, ToX(this->quadpoint_l1.points)) + Dot(N2, ToX(this->quadpoint_l2.points));
   };

   std::array<double, 3> X(const std::array<double, 2> &t0t1) const { return X(std::get<0>(t0t1), std::get<1>(t0t1)); };

   template <int i, int j>
   std::array<double, 3> D_X(const double t0, const double t1) const {
      auto [Nc, N0, N1, N2] = D_N6_new<i, j>(t0, t1);
      return Dot(Nc, ToX(this->quadpoint.points)) + Dot(N0, ToX(this->quadpoint_l0.points)) + Dot(N1, ToX(this->quadpoint_l1.points)) + Dot(N2, ToX(this->quadpoint_l2.points));
   };

   std::array<double, 3> cross(const double t0, const double t1) const { return Cross(D_X<1, 0>(t0, t1), D_X<0, 1>(t0, t1)); };

   /* -------------------------------------------------------------------------- */

   template <typename T>
   T Dot4(const std::array<T6d, 4> &N,
          const std::array<std::array<T, 6>, 4> &samples) const {
      return Dot(N[0], samples[0]) + Dot(N[1], samples[1]) + Dot(N[2], samples[2]) + Dot(N[3], samples[3]);
   };

   template <typename T>
   std::array<std::array<T, 6>, 4> convertions(const std::function<T(networkPoint *)> &conversion) const {
      std::array<std::array<T, 6>, 4> samples;
      for (auto i = 0; i < 6; ++i) {
         samples[0][i] = conversion(this->quadpoint.points[i]);
         samples[1][i] = conversion(this->quadpoint_l0.points[i]);
         samples[2][i] = conversion(this->quadpoint_l1.points[i]);
         samples[3][i] = conversion(this->quadpoint_l2.points[i]);
      }
      return samples;
   };

   Tddd interpolate(const double t0, const double t1, const std::function<Tddd(networkPoint *)> &conversion) const {
      return this->Dot4(N6_new(t0, t1), convertions(conversion));
   };

   Tddd interpolate(const Tdd t0t1, std::function<Tddd(networkPoint *)> conversion) const {
      return interpolate(std::get<0>(t0t1), std::get<1>(t0t1), conversion);
   };

   double interpolate(const double t0, const double t1, std::function<double(networkPoint *)> conversion) const {
      return this->Dot4(N6_new(t0, t1), convertions(conversion));
   };

   double interpolate(const Tdd t0t1, std::function<double(networkPoint *)> conversion) const {
      return interpolate(std::get<0>(t0t1), std::get<1>(t0t1), conversion);
   };

   template <std::size_t N, typename T>
   std::array<T, N> interpolate(const std::array<Tdd, N> &t0t1, std::function<T(networkPoint *)> conversion) const {
      std::array<T, N> ret;
      auto samples = convertions(conversion);
      for (auto i = 0; i < t0t1.size(); ++i)
         ret[i] = this->Dot4(N6_new(std::get<0>(t0t1[i]), std::get<1>(t0t1[i])), samples);
      return ret;
   };

   template <std::size_t N>
   std::vector<std::array<Tddd, N>> interpolate(const std::vector<std::array<Tdd, N>> &t0t1, std::function<Tddd(networkPoint *)> conversion) const {
      std::vector<std::array<Tddd, N>> ret(t0t1.size());
      auto samples = convertions(conversion);
      for (auto i = 0; i < t0t1.size(); ++i)
         for (auto j = 0; j < N; ++j)
            ret[i][j] = this->Dot4(N6_new(std::get<0>(t0t1[i][j]), std::get<1>(t0t1[i][j])), samples);
      return ret;
   };

   template <std::size_t N, typename T>
   std::vector<T> interpolate(const std::vector<Tdd> &t0t1, std::function<T(networkPoint *)> conversion) const {
      std::vector<T> ret(t0t1.size());
      auto samples = convertions(conversion);
      for (auto i = 0; i < N; ++i)
         ret[i] = this->Dot4(N6_new(std::get<0>(t0t1[i]), std::get<1>(t0t1[i])), samples);
      return ret;
   };

   /* -------------------------------------------------------------------------- */

   template <int i, int j>
   double D_interpolate(const double t0, const double t1, std::function<double(networkPoint *)> conversion) const {
      return Dot4(D_N6_new<i, j>(t0, t1), convertions(conversion));
   };

   template <int i, int j>
   double D_interpolate(const Tdd t0t1, std::function<double(networkPoint *)> conversion) const {
      return D_interpolate<i, j>(std::get<0>(t0t1), std::get<1>(t0t1), conversion);
   };

   template <int i, int j>
   std::vector<double> D_interpolate(const std::vector<Tdd> t0t1, std::function<double(networkPoint *)> conversion) const {
      std::vector<double> ret(t0t1.size());
      auto samples = convertions(conversion);
      for (auto k = 0; k < t0t1.size(); ++k)
         ret[k] = this->Dot4(D_N6_new<i, j>(std::get<0>(t0t1[k]), std::get<1>(t0t1[k])), samples);
      return ret;
   };

   template <int i, int j>
   Tddd D_interpolate(const double t0, const double t1, std::function<Tddd(networkPoint *)> conversion) const {
      return this->Dot4(D_N6_new<i, j>(t0, t1), convertions(conversion));
   };

   template <int i, int j>
   Tddd D_interpolate(const Tdd t0t1, std::function<Tddd(networkPoint *)> conversion) const {
      return this->Dot4(D_N6_new<i, j>(std::get<0>(t0t1), std::get<1>(t0t1)), convertions(conversion));
   };

   template <int i, int j>
   std::vector<Tddd> D_interpolate(const std::vector<Tdd> t0t1, std::function<Tddd(networkPoint *)> conversion) const {
      std::vector<Tddd> ret(t0t1.size());
      auto samples = convertions(conversion);
      for (auto k = 0; k < t0t1.size(); ++k)
         ret[k] = this->Dot4(D_N6_new<i, j>(std::get<0>(t0t1[k]), std::get<1>(t0t1[k])), samples);
      return ret;
   };
};

#include "MooringLine.hpp"
#include "NetworkUtility.hpp"
#include "networkFace.hpp"
#include "networkLine.hpp"
#include "networkTetra.hpp"

#endif
