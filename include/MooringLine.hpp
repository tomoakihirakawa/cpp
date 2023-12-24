#ifndef MooringLine_H
#define MooringLine_H
#include "Network.hpp"
#include "integrationOfODE.hpp"

/*DOC_EXTRACT 1_MooringLine

## 浮体係留用に`Network`の派生クラスを作成

`networkLine`には，`natural_length`，`stiffness`，`damping`，`weight_per_unit_length`の4つのパラメータを持たせる．
`natural_length`は，`moorinLine`の`total_length`と`MooringLine`の`getPoints().size()`から決まる．
それをまとめる`MooringLine`は，`total_length`を持つ．

```cpp
const int n_points = 100;
const double total_length = 522.; //! [m]
std::array<double,3> X_anchor = {r * std::cos(0.), r * std::sin(0.), h}
std::array<double,3> X_ = {0., 0., 0.}
auto mooring = new MooringLine(X_anchor, X_fairlead, total_length, n_points);　//mooringオブジェクトの作成
```

次のように物性の設定を行う．

```cpp
mooring->setDensityStiffnessDampingDiameter(density, stiffness, damp, diam);
```

考えられる係留索クラスの使い道は，
単点（フェアリード）の位置がわかっていて，そこでの係留索が及ぼす力を計算するというもの．
前時刻の係留索の曲線から，少しずつ係留索を時間発展させて，現在の係留索の曲線を求め，その節点の配置から，力を計算できる．

移動速度がわかっている節点にはその条件を与え，わからないものに関しては，各時刻の張力と重力と抗力から，運動方程式を解いていく．


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

   /* -------------------------------------------------------------------------- */
   /*                                     CFL                                    */
   /* -------------------------------------------------------------------------- */
   //! function to get dt automatically
   double get_dt(const double max_dt) {
      const double C_VELOCITY = 0.01, C_ACCEL = 0.01;
      double dt_cfl = max_dt, dt_tmp;
      // CFL based on velocity
      for (const auto& p : this->getPoints())
         for (const auto& l : p->getLines()) {
            double norm = (l->stiffness / p->mass * 0.00001) * Norm(Tddd{p->velocity[0], p->velocity[1], p->velocity[2]});
            if (norm == 0)
               continue;
            dt_tmp = C_VELOCITY * l->length() / norm;
            if (dt_cfl > dt_tmp)
               dt_cfl = dt_tmp;
         }

      // CFL based on acceleration
      for (const auto& p : this->getPoints())
         for (const auto& l : p->getLines()) {
            double norm = (l->stiffness / p->mass * 0.0001) * Norm(Tddd{p->acceleration[0], p->acceleration[1], p->acceleration[2]});
            if (norm != 0.) {
               dt_tmp = C_ACCEL * std::sqrt(l->length() / norm);
               if (dt_cfl > dt_tmp)
                  dt_cfl = dt_tmp;
            }
            norm = (l->stiffness / p->mass * 0.0001) * _GRAVITY_;
            if (norm != 0) {
               dt_tmp = C_ACCEL * std::sqrt(l->length() / norm);
               if (dt_cfl > dt_tmp)
                  dt_cfl = dt_tmp;
            }
         }

      return dt_cfl;
   }

   /*DOC_EXTRACT 2_MooringLine

   `simulate`関数は，`netwrokPoint`が持つ`getForce`関数を用いて力を計算する．
   他の関数では使われない４次のルンゲクッタクラスを使っている．`RK_velocity_sub`，`RK_X_sub`，`RK_force_sub`

   `RK_velocity_sub`，`RK_X_sub`，`RK_force_sub`だけにとどまり，実際には節点は移動しないし，速度も変わらないので，
   次時刻の速度や位置は，`RK_velocity_sub`，`RK_X_sub`の`get_x`関数で取得する．



   */

   double DragForceCoefficient = 0.3;

   void simulate(const double current_time,
                 const double dt,
                 const std::function<void(networkPoint*)> setBoundaryCondition) {

      double dt_acum = 0;
      bool first = true;

      //! 実際の値は変更されない
      while (dt_acum < dt) {
         double dt_cfl = get_dt(dt);
         if (dt_acum + dt_cfl >= dt)
            dt_cfl = dt - dt_acum;
         /* -------------------------------------------------------------------------- */
         //! initialize RK
         for (auto& p : this->getPoints()) {
            p->RK_velocity_sub.initialize(dt_cfl, current_time, first ? p->velocityTranslational() : p->RK_velocity_sub.get_x(), 4);
            p->RK_X_sub.initialize(dt_cfl, current_time, first ? p->X : p->RK_X_sub.get_x(), 4);
         }
         first = false;
         /* -------------------------------------------------------------------------- */
         // Print("simulate ", dt_cfl, ",", dt_acum, ",", dt);
         std::array<double, 3> a;

         while (1) {

            for (auto& p : this->getPoints()) {
               a = (p->getTension() + p->getDragForce(DragForceCoefficient) + p->getGravitationalForce()) / p->mass;
               std::get<0>(p->acceleration) = std::get<0>(a);  // accelは変更しても構わない
               std::get<1>(p->acceleration) = std::get<1>(a);  // accelは変更しても構わない
               std::get<2>(p->acceleration) = std::get<2>(a);  // accelは変更しても構わない
               setBoundaryCondition(p);
            }

            for (auto& p : this->getPoints()) {
               p->RK_X_sub.push(p->RK_velocity_sub.get_x());
               p->RK_velocity_sub.push(p->accelTranslational());
            }

            if ((*this->getPoints().begin())->RK_X_sub.finished)
               break;
         }

         dt_acum += dt_cfl;
         if (dt_acum >= dt)
            return;
      }
   }

   void setEquilibriumState(const std::function<void(networkPoint*)> setBoundaryCondition) {
      double norm_total_velocity = 0;
      double current_DragForceCoefficient = this->DragForceCoefficient;
      this->DragForceCoefficient = 1000.;
      for (auto i = 0; i < 1000; ++i) {
         this->simulate(0, 5., setBoundaryCondition);
         norm_total_velocity = 0;
         for (const auto& p : this->getPoints())
            norm_total_velocity += Norm(p->RK_velocity_sub.get_x());
         norm_total_velocity /= this->getPoints().size();
         if (norm_total_velocity < 1e-3)
            break;
         std::cout << norm_total_velocity << std::endl;
         applyMooringSimulationResult();
      }
      this->DragForceCoefficient = current_DragForceCoefficient;
   }

   void applyMooringSimulationResult() {
      for (auto& p : this->getPoints()) {
         p->setX(p->RK_X_sub.get_x());
         auto v = p->RK_velocity_sub.get_x();
         p->velocity[0] = v[0];
         p->velocity[1] = v[1];
         p->velocity[2] = v[2];
      }
   }
};

#endif