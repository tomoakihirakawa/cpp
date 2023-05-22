#ifndef integrationOfODE_H
#define integrationOfODE_H

#pragma once

#include "basic.hpp"
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

// \label{ODE:RungeKutta}
template <typename T>
struct RungeKuttaCommon {
   double dt_fixed;  // 基本となるステップ間隔,固定
   /* ----------------------- 毎回変わる値 ----------------------- */
   double dt;           // 次回計算する必要がある値になる増分,回数ごとに変化
   std::vector<T> _dX;  // これにプッシュしていく，c1k1,c2k2,c3k3,c4k4,,,//このベクトルの長さで何回目かを判断数し，必要な回数まできたら，finished=trueとなる
   T dX;                //_dXやdtを基に計算される
   /* ------------------------- 初期値 ------------------------ */
   double t_init;  // 初期値,固定
   T Xinit;
   /* ------------------------ 積分結果 ------------------------ */
   bool finished;  // 全てそろって積分が評価できるようになったらtrue
                   // 最終結果はdXに代入される;
   void displayStatus() {
      std::cout << "dt : " << red << this->dt << colorOff << std::endl;
      // std::cout << "_dX : " << red << this->_dX << colorOff << std::endl;
      std::cout << "dX : " << red << this->dX << colorOff << std::endl;
   };
   int steps;
   int current_step;
   RungeKuttaCommon(){};
   RungeKuttaCommon(const double dt_IN, const double t0, const T &X0, int stepsIN)
       : dt_fixed(dt_IN), dt(0), t_init(t0), Xinit(X0), _dX({}), dX(X0), steps(stepsIN), current_step(0){};
   ~RungeKuttaCommon(){};

   void initialize(const double dt_IN, const double t0, const T &X0, int stepsIN) {
      this->dt_fixed = dt_IN;
      this->dt = 0;
      this->t_init = t0;
      this->Xinit = X0;
      this->_dX = {};
      this->dX = X0;
      this->steps = stepsIN;
      if (stepsIN == 0 && stepsIN > 4)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "stepsINが1~4の整数でなければならない");
      this->current_step = 0;
   };

   T getXinit() const { return this->Xinit; };

   T getX() const {
      if (this->current_step == 0)
         return this->Xinit;
      else
         return this->Xinit + this->dX;
   };
   T get_x() const { return getX(); };

   T getdX() const { return this->dX; };
   double gett() const { return this->t_init + this->dt; };
   T get_t() const { return gett(); };

   double getTime() const { return this->t_init + this->dt; };
   double getNextTime() const { return this->t_init + this->getdt(); };
   double getdt() const {
      // 次の計算は，t+dtを狙って計算することになる．
      // ここでのdtを返す
      if (this->steps == 1)
         return dt_fixed;  //! 何もプッシュしていない初期状態
      else if (this->steps == 2) {
         switch (_dX.size()) {
            case 0:
               return dt_fixed;
            case 1:
               return dt_fixed / 2.;
            default:
               return dt_fixed;
         }
      } else if (this->steps == 3) {
         switch (_dX.size()) {
            case 0:  //! 何もプッシュしていない初期状態
               return dt_fixed / 2.;
            case 1:
               return dt_fixed;
            case 2:
               return dt_fixed / 6.;
            default:
               return dt_fixed;
         }
      } else if (this->steps == 4) {
         switch (_dX.size()) {
            case 0:  //! 何もプッシュしていない初期状態
               return dt_fixed / 2.;
            case 1:  // -> {f2, f1(t,v(t))}
               return dt_fixed / 2.;
            case 2:  // -> {f2, f1(t,v(t))}
               return dt_fixed;
            case 3:  // -> {f2, f1(t,v(t))}
               return dt_fixed / 6.;
            default:
               return dt_fixed;
         }
      } else
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
   };

   bool repush(const T &dXdt_IN) {
      if (!_dX.empty())
         _dX.pop_back();
      return this->push(dXdt_IN);
   };

   bool push(const T &dXdt_IN) {
      // _dX.insert(_dX.begin(), dXdt_IN * dt_fixed);
      _dX.push_back(dXdt_IN * dt_fixed);
      current_step++;
      if (this->steps == 1) {
         switch (_dX.size()) {
            case 1:  // -> {f1(t,v(t))}
               dX = dXdt_IN * (dt = dt_fixed);
               return this->finished = true;
            default:
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
         }
      } else if (this->steps == 2) {
         switch (_dX.size()) {
            case 1:  // -> {f1(t,v(t))}
               dX = dXdt_IN * (dt = dt_fixed);
               return this->finished = false;
            case 2:  // -> {f3, f2, f1(t,v(t))}
               dX = (_dX[0] + _dX[1]) / 2.;
               dt = dt_fixed;
               return this->finished = true;
            default:
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
         }
      } else if (this->steps == 3) {
         switch (_dX.size()) {
            case 1:  // -> {f1(t,v(t))}
               dX = dXdt_IN * (dt = dt_fixed / 2.);
               return finished = false;
            case 2:  // -> {f2, f1(t,v(t))}
               dX = dXdt_IN * (dt = dt_fixed);
               return finished = false;
            case 3:  // -> {f3, f2, f1(t,v(t))}
               dX = (_dX[0] + 4. * _dX[1] + _dX[2]) / 6.;
               dt = dt_fixed;
               return finished = true;
            default:
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
         }
      } else if (this->steps == 4) {
         switch (_dX.size()) {
            case 1:  // -> {f1(t,v(t))}
               dX = dXdt_IN * (dt = dt_fixed / 2.);
               return finished = false;
            case 2:  // -> {f2, f1(t,v(t))}
               dX = dXdt_IN * (dt = dt_fixed / 2.);
               return finished = false;
            case 3:  // -> {f3, f2, f1(t,v(t))}
               dX = dXdt_IN * (dt = dt_fixed);
               return finished = false;
            case 4:  // -> dXdt[1][] = {f4, f3, f2, f1(t,v(t))}
               dX = (_dX[0] + 2. * _dX[1] + 2. * _dX[2] + _dX[3]) / 6.;
               dt = dt_fixed;
               return finished = true;
            default:
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
         }
      } else {
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
      }
   };

   /* ------------------------------------------------------ */
   /*            微分を与えて，次の時刻に使うべきXを返す            */
   /* ------------------------------------------------------ */
   T getX(const T &dXdt_IN) const {
      if (this->steps == 1)
         return this->Xinit + dXdt_IN * dt_fixed;  //! 何もプッシュしていない初期状態
      else if (this->steps == 2) {
         switch (_dX.size()) {
            case 0:  //! 何もプッシュしていない初期状態
               return this->Xinit + dXdt_IN * dt_fixed;
            default:
               return this->Xinit + (_dX[0] + dXdt_IN * dt_fixed) / 2.;
         }
      } else if (this->steps == 3) {
         switch (_dX.size()) {
            case 0:  //! 何もプッシュしていない初期状態
               return this->Xinit + dXdt_IN * dt_fixed / 2.;
            case 1:  // -> {f2, f1(t,v(t))}
               return this->Xinit + dXdt_IN * dt_fixed;
            default:
               return this->Xinit + (_dX[0] + 4. * _dX[1] + dXdt_IN * dt_fixed) / 6.;
         }
      } else if (this->steps == 4) {
         switch (_dX.size()) {
            case 0:  //! 何もプッシュしていない初期状態
               return this->Xinit + dXdt_IN * dt_fixed / 2.;
            case 1:  // -> {f2, f1(t,v(t))}
               return this->Xinit + dXdt_IN * dt_fixed / 2.;
            case 2:  // -> {f3, f2, f1(t,v(t))}
               return this->Xinit + dXdt_IN * dt_fixed;
            default:  // -> {f3, f2, f1(t,v(t))}
               return this->Xinit + (_dX[0] + 2. * _dX[1] + 2. * _dX[2] + dXdt_IN * dt_fixed) / 6.;
         }
      } else {
         std::stringstream ss;
         ss << std::to_string(this->steps);
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．steps=" + ss.str());
      }
   };

   T getdXdt(const T &dXdt_IN) const {
      if (this->steps == 1)
         return dXdt_IN;  //! 何もプッシュしていない初期状態
      else if (this->steps == 2) {
         switch (_dX.size()) {
            case 0:  //! 何もプッシュしていない初期状態
               return dXdt_IN * dt_fixed;
            default:
               return (_dX[0] / dt_fixed + dXdt_IN);
         }
      } else if (this->steps == 3) {
         switch (_dX.size()) {
            case 0:  //! 何もプッシュしていない初期状態
               return dXdt_IN;
            case 1:  // -> {f2, f1(t,v(t))}
               return dXdt_IN;
            default:
               return (_dX[0] / dt_fixed + 4. * _dX[1] / dt_fixed + dXdt_IN) / 6.;
         }
      } else if (this->steps == 4) {
         switch (_dX.size()) {
            case 0:  //! 何もプッシュしていない初期状態
               return dXdt_IN;
            case 1:  // -> {f2, f1(t,v(t))}
               return dXdt_IN;
            case 2:  // -> {f3, f2, f1(t,v(t))}
               return dXdt_IN;
            default:  // -> {f3, f2, f1(t,v(t))}
               return (_dX[0] / dt_fixed + 2. * _dX[1] / dt_fixed + 2. * _dX[2] / dt_fixed + dXdt_IN) / 6.;
         }
      } else {
         std::stringstream ss;
         ss << std::to_string(this->steps);
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．steps=" + ss.str());
      }
   };
};

template <typename T>
struct RungeKutta : public RungeKuttaCommon<T> {
   RungeKutta(const double dt_IN, const double t0, const T &X0, int stepsIN) : RungeKuttaCommon<T>(dt_IN, t0, X0, stepsIN){};
   RungeKutta() : RungeKuttaCommon<T>(){};
};
template <>
struct RungeKutta<std::vector<double>> : public RungeKuttaCommon<std::vector<double>> {
   RungeKutta(const double dt_IN, const double t0, const std::vector<double> &X0, int stepsIN) : RungeKuttaCommon<std::vector<double>>(dt_IN, t0, X0, stepsIN) {
      this->dX = V_d(X0.size(), 0.);
   };
   RungeKutta() : RungeKuttaCommon<std::vector<double>>(){};
};
template <>
struct RungeKutta<double> : public RungeKuttaCommon<double> {
   RungeKutta(const double dt_IN, const double t0, const double &X0, int stepsIN) : RungeKuttaCommon<double>(dt_IN, t0, X0, stepsIN) {
      this->dX = 0.;
   };
   RungeKutta() : RungeKuttaCommon<double>(){};
};
template <std::size_t N>
struct RungeKutta<std::array<double, N>> : public RungeKuttaCommon<std::array<double, N>> {
   RungeKutta(const double dt_IN, const double t0, const std::array<double, N> &X0, int stepsIN) : RungeKuttaCommon<std::array<double, N>>(dt_IN, t0, X0, stepsIN) { this->dX.fill(0.); };
   RungeKutta() : RungeKuttaCommon<std::array<double, N>>(){};
};
/* -------------------------------------------------------------------------- */

// \label{ODE:LeapFrog}

/*DOC_EXTRACT ODE

## LeapFrog

時間間隔$\Delta t$が変化する場合でも使える形でプログラムしている．
$\Delta t$が変化する場合，"半分蹴って-移動-半分蹴って"，"半分蹴って-移動-半分蹴って"の手順を繰り返す．

*/

template <typename T>
class LeapFrog {
  public:
   LeapFrog() : is_first(true), finished(false){};
   LeapFrog(double dt, double t0, const T &x0, const T &v0) : dt(dt), t(t0), x(x0), v(v0), is_first(true), finished(false){};

   void initialize(double dt, double t0, const T &x0, const T &v0) {
      this->dt = dt;
      this->t = t0;
      this->x = x0;
      this->v = v0;
      this->is_first = true;
      this->finished = false;
   }

   bool is_first, finished;
   T v_full;

   void push(const T &a) {
      double half_dt = 0.5 * dt;
      if (is_first) {
         // give a = a(get_t(),get_x())
         v_old = v;
         a_old = a;
         v += half_dt * a;  // half-step update of v
         x += dt * v;       // full-step update of x
         t += dt;
         is_first = false;
         finished = false;
      } else {
         // give a = a(get_t(),get_x())
         v_old = v;
         a_old = a;
         v += half_dt * a;  // half-step update of v
         is_first = true;
         finished = true;
      }
   }

   // The Leapfrog method (also known as the midpoint method)
   void repush(const T &a) {
      double half_dt = 0.5 * dt;
      if (!is_first) {
         // give a = a(get_t(),get_x())
         v_old = v;
         a_old = a;
         v += half_dt * (a - a_old);  // half-step update of v
         x += dt * (v - v_old);       // full-step update of x
         is_first = false;
         finished = false;
      } else {
         // give a = a(get_t(),get_x())
         v_old = v;
         a_old = a;
         v += half_dt * (a - a_old);  // half-step update of v
         is_first = true;
         finished = true;
      }
   }

   const double get_dt() const { return dt; }
   const double get_t() const { return t; }
   const T &get_x() const { return x; }
   const T &get_v() const { return v; }

  private:
   double dt, t;
   T x, v, a_old, v_old;
};

/* -------------------------------------------------------------------------- */

template <typename T>
class VelocityVerlet {
  public:
   VelocityVerlet(double dt, double t0, const T &x0, const T &v0, const T &a0) : dt(dt), t(t0), x(x0), v(v0), a(a0){};

   void push(const T &new_a) {
      v += dt * (a + new_a) * 0.5;  // this should be perfomed after x is updated though this sometimes stabilizes the calculation
      x += dt * v + 0.5 * dt * dt * new_a;
      t += dt;
      a = new_a;
   }

   const double get_t() const { return t; }
   const T &get_x() const { return x; }
   const T &get_v() const { return v; }
   const T &get_a() const { return a; }

  private:
   double dt, t;
   T x, v, a;
};

/* -------------------------------------------------------------------------- */

template <typename T>
class Beeman {
  public:
   Beeman(double dt, double t0, const T &x0, const T &v0, const T &a0, const T &old_a) : dt(dt), t(t0), x(x0), v(v0), a(a0), old_a(old_a){};
   Beeman(double dt, double t0, const T &x0, const T &v0, const T &a0) : dt(dt), t(t0), x(x0), v(v0), a(a0), old_a(a0){};

   void push(const T &new_a) {
      x += v * dt + (2.0 / 3) * a * dt * dt - (1.0 / 6) * old_a * dt * dt;
      T old_v = v;
      v += (1.0 / 3) * new_a * dt + (5.0 / 6) * a * dt - (1.0 / 6) * old_a * dt;
      t += dt;
      old_a = a;
      a = new_a;
   }

   double get_t() const { return t; }
   const T &get_x() const { return x; }
   const T &get_v() const { return v; }
   const T &get_a() const { return a; }

  private:
   double dt, t;
   T x, v, a, old_a;
};

/* -------------------------------------------------------------------------- */

template <typename T>
class Gear {
  public:
   Gear(double dt, double t0, const T &x0, const T &v0, const T &a0) : dt(dt), t(t0), x(x0), v(v0), a(a0){};

   void push(const T &new_a) {
      T pred_x = x + v * dt + 0.5 * a * dt * dt;
      T pred_v = v + a * dt;
      T error = new_a - a;
      x = pred_x + (1.0 / 6) * error * dt * dt;
      v = pred_v + 0.5 * error * dt;
      t += dt;
      a = new_a;
   }

   double get_t() const { return t; }
   const T &get_x() const { return x; }
   const T &get_v() const { return v; }
   const T &get_a() const { return a; }

  private:
   double dt, t;
   T x, v, a;
};

/* -------------------------------------------------------------------------- */

#endif