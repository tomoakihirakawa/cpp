#ifndef RigidBodyDynamics_H
#define RigidBodyDynamics_H

#include "basic.hpp"
#include "interpolations.hpp"

// @ ============================================================== */
// @                        剛体の力学に関する                         */
// @ ============================================================== */
//@ メッシュが関わらない剛体の運動を表すクラス
struct RigidBodyDynamics {

   RigidBodyDynamics() : quaternion(Quaternion()){};
   virtual ~RigidBodyDynamics() = default;

   T6d forced_velocity = {0., 0., 0., 0., 0., 0.};
   T6d forced_acceleration = {0., 0., 0., 0., 0., 0.};

   InterpolationBspline<std::array<double, 6>> intpMotionRigidBody;
   //* ------------------------------------------------------ */
   //*                     運動を表す量                         */
   //* ------------------------------------------------------ */
   T6d force = {0., 0., 0., 0., 0., 0.};
   T6d velocity = {0., 0., 0., 0., 0., 0.};  // = {velocity,angular velocity}
   T6d acceleration = {0., 0., 0., 0., 0., 0.};
   InterpolationLagrange<T6d> interp_accel;
   Tddd velocityTranslational() const { return {std::get<0>(velocity), std::get<1>(velocity), std::get<2>(velocity)}; };
   Tddd velocityRotational() const { return {std::get<3>(velocity), std::get<4>(velocity), std::get<5>(velocity)}; };
   Tddd accelTranslational() const { return {std::get<0>(acceleration), std::get<1>(acceleration), std::get<2>(acceleration)}; };
   Tddd accelRotational() const { return {std::get<3>(acceleration), std::get<4>(acceleration), std::get<5>(acceleration)}; };

   //! 固定された空間座標におけるベクトルであることを頭に入れておくこと．
   //! 回転，移動をする物体の座標系ではないので，固定座標にとって，回転前後でinertiaは書き換える必要がある．
   //! inertiaの慣性モーメントはそのまま固定座標における回転行列をかけて，更新すればいい
   T6d &F = this->force;
   T6d &A = this->acceleration;
   T6d &V = this->velocity;
   //* ------------------------------------------------------ */
   //*                   位置や姿勢を表す量                      */
   //* ------------------------------------------------------ */
   Tddd center_of_mass = {0., 0., 0.};  // 現在の座標を表す
   Quaternion quaternion;               // 現在の姿勢を表す，初期のクォータニオンは固定で{1,0,0,0}
   Tddd &COM = this->center_of_mass;
   Quaternion &Q = this->quaternion;
   //* ------------------------------------------------------ */
   //*                   　　 不変の量                          */
   //* ------------------------------------------------------ */
   Tddd initial_center_of_mass = {0., 0., 0.};  // 初期のの座標を表す
   Tddd &ICOM = this->initial_center_of_mass;
   double mass = 1E+20;
   T6d inertia = {1E+20, 1E+20, 1E+20, 1E+20, 1E+20, 1E+20};
   T6d &I = this->inertia;
   /* ------------------------------------------------------ */
   double mass_tmp = 1E+20;
   T6d inertia_tmp = {1E+20, 1E+20, 1E+20, 1E+20, 1E+20, 1E+20};
   // make it very heavy to stop moving
   void makeImmovable() {
      this->mass_tmp = this->mass;
      this->inertia_tmp = this->inertia;
      this->mass = 1e+20;
      this->inertia = {1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20};
   };
   // release
   void makeMovable() {
      this->mass = this->mass_tmp;
      this->inertia = this->inertia_tmp;
   };

   std::tuple<double, double, double, T3Tddd, T3Tddd> getInertiaGC()  // Global coordinate
   {
      auto [mx, my, mz, Ix, Iy, Iz] = this->inertia;
      auto R = this->quaternion.Rv();
      auto RT = Transpose(R);
      T3Tddd I_accounting_float_attitude = Dot(RT, Dot(T3Tddd{{{Ix, 0., 0.},
                                                               {0., Iy, 0.},
                                                               {0., 0., Iz}}},
                                                       R));
      // T3Tddd IG_inv = Dot(Transpose(R), Dot(T3Tddd{{{1. / Ix, 0., 0.}, {0., 1. / Iy, 0.}, {0., 0., 1. / Iz}}}, R));
      T3Tddd IG_inv = Dot(RT, T3Tddd{{R[0] / Ix, R[1] / Iy, R[2] / Iz}});

      return {mx, my, mz, I_accounting_float_attitude, IG_inv};
      // return {mx, my, mz, T3Tddd{{{Ix, 0., 0.}, {0., Iy, 0.}, {0., 0., Iz}}}};
   };

   T6d getInertiaBC() { return this->inertia; };  // Body coordinate
   Tddd getMass3D() { return {std::get<0>(this->inertia), std::get<1>(this->inertia), std::get<2>(this->inertia)}; };

   /* ------------------------------------------------------ */

   Tddd velocityRigidBody(const Tddd &X) const { return velocityTranslational() + Cross(velocityRotational(), X - this->COM); };  //! \label{velocityRigidBody}
   Tddd accelRigidBody(const Tddd &X) const { return accelTranslational() + Cross(accelRotational(), X - this->COM); };           //! \label{accelRigidBody}
   void calcAccelFromForce() {
      //! 慣性を設定しておく必要がある．
      this->acceleration = this->force / this->inertia;
   };
};

#endif