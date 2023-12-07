#ifndef MooringLine_H
#define MooringLine_H
#include "Network.hpp"
#include "integrationOfODE.hpp"

/*DOC_EXTRACT 1_MooringLine

## 浮体係留用に`Network`の派生クラスを作成

`networkLine`には，`natural_length`，`stiffness`，`damping`，`weight_per_unit_length`の4つのパラメータを持たせる．

`natural_length`は，`moorinLine`の`total_length`と`MooringLine`の`getPoints().size()`から決まる．

それをまとめる`MooringLine`は，`total_length`を持つ．

TODO: 内部で，自動でタイムステップを細かく取り，与えられた時間までシミュレーションするようにしたい．

*/

/* -------------------------------------------------------------------------- */

class MooringLine : public Network {
  public:
   // default destructor
   ~MooringLine() = default;

   double total_length;
   networkPoint* firstPoint = nullptr;
   networkPoint* lastPoint = nullptr;

   MooringLine(const std::array<double, 3> X0, const std::array<double, 3> X1, const double total_length, const int n) : Network(), total_length(total_length) {
      std::vector<networkPoint*> points;
      for (int i = 0; i < n; ++i) {
         auto p = new networkPoint(this, X0 + (X1 - X0) * i / (n - 1.));
         points.push_back(p);
         if (i == 0) firstPoint = p;
         if (i == n - 1) lastPoint = p;
      }
      for (auto i = 0; i < points.size() - 1; ++i) new networkLine(this, points[i], points[i + 1]);

      setNaturalLength();
   }
   double natural_length() const { return total_length / (this->getPoints().size() - 1); };

   //! total_length [m]
   //! weight_per_unit_length [kg/m]
   //! density [kg/m^3]
   //! natural_length [m]

   void setNaturalLength() {
      for (auto& l : this->getLines()) l->natural_length = natural_length();
   };

   void setDensityStiffnessDampingDiameter(const double density, const double stiffness, const double damp, const double diam) {
      for (auto& l : this->getLines()) {
         l->weight_per_unit_length = density * l->natural_length;
         l->stiffness = stiffness;
         l->damping = damp;
         l->diameter = diam;
      }
      for (auto& p : this->getPoints()) {
         p->mass = 0;
         for (auto& l : p->getLines())
            p->mass += l->natural_length * l->weight_per_unit_length / 2.;
      }
   };

   /* ------------------------------- FOR SIMULATIONS ------------------------------- */

   // std::array<double, 3> getTension(const networkPoint* p) {
   //    std::array<double, 3> force, v, relative_velocity;
   //    force.fill(0);
   //    double disp, strain;
   //    for (const auto& l : p->getLines()) {
   //       v = (*l)(p)->RK_X.getX() - p->RK_X.getX();
   //       disp = Norm(v) - l->natural_length;
   //       strain = disp / l->natural_length;
   //       if (disp > 0.)
   //          force += l->stiffness * strain * Normalize(v);
   //       relative_velocity = (*l)(p)->RK_velocity.getX() - p->RK_velocity.getX();
   //       force += l->damping * relative_velocity / p->RK_velocity.dt_fixed;
   //    }
   //    return force;
   // }

   // std::array<double, 3> getDragForce(const networkPoint* p) {
   //    std::array<double, 3> drag_force, relative_velocity, mean_v, fluid_velocity, normalized_relative_velocity, A2B;
   //    drag_force.fill(0);
   //    fluid_velocity.fill(0);
   //    double Cd = 0.3, A;
   //    for (const auto& l : p->getLines()) {
   //       mean_v = 0.5 * ((*l)(p)->RK_velocity.getX() + p->RK_velocity.getX());
   //       relative_velocity = fluid_velocity - mean_v;
   //       // A = M_PI * std::pow(l->diameter, 2);
   //       A2B = (*l)(p)->RK_X.getX() - p->RK_X.getX();
   //       normalized_relative_velocity = Normalize(relative_velocity);
   //       A = l->diameter * Norm(A2B) * 0.5;                                                            //! area
   //       A *= Norm(normalized_relative_velocity - Dot(normalized_relative_velocity, Normalize(A2B)));  //! projected area
   //       drag_force += 0.5 * _WATER_DENSITY_ * Dot(relative_velocity, relative_velocity) * Cd * A * normalized_relative_velocity;
   //    }
   //    return drag_force;
   // }

   // // std::array<double, 3> getAcceleration(const networkPoint* p) {
   // //    return getTension(p) / p->mass + getDragForce(p) / p->mass + _GRAVITY3_;
   // // }

   // std::array<double, 3> getForce(const networkPoint* p) {
   //    return getTension(p) + getDragForce(p) + p->mass * _GRAVITY3_;
   // }

   void update() {
      std::array<double, 3> v;
      for (auto& p : this->getPoints()) {
         p->setX(p->RK_X.get_x());
         v = p->RK_velocity.get_x();
         p->velocity[0] = v[0];
         p->velocity[1] = v[1];
         p->velocity[2] = v[2];
      }
   }

   void simulate(const double current_time,
                 const double dt,
                 const std::function<void(networkPoint*)> setBoundaryCondition) {

      double dt_acum = 0;
      const double C_VELOCITY = 0.05, C_ACCEL = 0.05;

      while (dt_acum < dt) {
         double dt_cfl = dt;
         double dt_tmp;

         // CFL based on velocity
         for (const auto& p : this->getPoints())
            for (const auto& l : p->getLines()) {
               double norm = (l->stiffness / p->mass * 0.001) * Norm(Tddd{p->velocity[0], p->velocity[1], p->velocity[2]});
               if (norm == 0)
                  continue;
               dt_tmp = C_VELOCITY * l->length() / norm;
               if (dt_cfl > dt_tmp)
                  dt_cfl = dt_tmp;
            }

         // CFL based on acceleration
         for (const auto& p : this->getPoints())
            for (const auto& l : p->getLines()) {
               double norm = (l->stiffness / p->mass * 0.001) * Norm(Tddd{p->acceleration[0], p->acceleration[1], p->acceleration[2]});
               if (norm != 0.) {
                  dt_tmp = C_ACCEL * std::sqrt(l->length() / norm);
                  if (dt_cfl > dt_tmp)
                     dt_cfl = dt_tmp;
               }
               norm = (l->stiffness / p->mass * 0.001) * _GRAVITY_;
               if (norm != 0) {
                  dt_tmp = C_ACCEL * std::sqrt(l->length() / norm);
                  if (dt_cfl > dt_tmp)
                     dt_cfl = dt_tmp;
               }
            }

         if (dt_acum + dt_cfl >= dt)
            dt_cfl = dt - dt_acum;

         //! initialize RK
         for (auto& p : this->getPoints()) {
            p->RK_velocity.initialize(dt_cfl, current_time, p->velocityTranslational(), 4);
            p->RK_X.initialize(dt_cfl, current_time, p->X, 4);
            p->RK_force.initialize(dt_cfl, current_time, p->forceTranslational(), 4);
         }

         // Print("simulate", dt, current_time);
         std::array<double, 3> a, force;
         while (1) {
            for (auto& p : this->getPoints()) {
               force = p->getForce();
               p->force[0] = force[0];
               p->force[1] = force[1];
               p->force[2] = force[2];
               a = force / p->mass;
               p->acceleration[0] = a[0];
               p->acceleration[1] = a[1];
               p->acceleration[2] = a[2];
               setBoundaryCondition(p);
            }

            for (auto& p : this->getPoints()) {
               p->RK_X.push(p->velocityTranslational());
               p->RK_velocity.push(p->accelTranslational());
               p->RK_force.push(p->forceTranslational());
               // std::cout << "a = " << p->accelTranslational() << " v = " << p->velocityTranslational() << std::endl;
            }

            if ((*this->getPoints().begin())->RK_X.finished)
               break;
         }

         this->update();

         dt_acum += dt_cfl;
         if (dt_acum >= dt)
            return;
      }
   }
};

#endif