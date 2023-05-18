#include "GNUPLOT.hpp"
#include "integrationOfODE.hpp"

/*DOC_EXTRACT ODE

# ODEの初期値問題

## 減衰調和振動子/Damped Harmonic Oscillatorの例

減衰調和振動子の式から，
次のような加速度$a(x,v)=\frac{d^2x}{dt^2}$を
\ref{DampedHrmonicOscillator:declOfAccel}{プログラム中で宣言}し，

$$
\begin{align*}
m \frac{d^2x}{dt^2} + b \frac{dx}{dt} + k x &= 0\\
\rightarrow a(x,v) &= -\gamma v - \omega^2 x, \quad v=\frac{dx}{dt},\quad \gamma=\frac{b}{m}, \quad \omega^2=\frac{k}{m}
\end{align*}
$$

$\gamma = 1, \omega = 10$として，初期値問題をといてみる．
加速度の評価回数$N$を合わせて比較した例：

| ![](figN25.png) | ![](figN50.png) |  ![](figError.png) |
|:---:|:---:|:---:|
|$N=25$ evaluations|$N=50$ evaluations|the sum of differences|

   * \ref{DampedHrmonicOscillator:BackwardEuler}{後退オイラー}の１回の計算で溜まる誤差は$O(\Delta t^2)$．次時刻における速度と加速度が正確に計算できなければ使えない．
   * \ref{DampedHrmonicOscillator:LeapFrog}{リープフロッグ}の１回の計算で溜まる誤差は$O({\Delta t}^3)$となる．\ref{ODE:RungeKutta}{LeapFrogのクラス}
   * \ref{DampedHrmonicOscillator:RungeKutta}{4次のルンゲクッタ}の１回の計算で溜まる誤差は$O({\Delta t}^5)$となる．しかし，加速度を4階も計算する必要がある．\ref{ODE:LeapFrog}{RungeKuttaのクラス}

*/

const double m = 1;
const double gamma = 1;
const double omega = 10;
const double b = gamma * m;
const double k = omega * omega * m;

// \label{DampedHrmonicOscillator:declOfAccel}
double acceleration(double x, double v) {
   return -(b * v + k * x) / m;
}

double solution_x(double t) {
   return exp(-gamma / 2 * t) * cos(sqrt(omega * omega - gamma * gamma / 4.) * t);
}

double solution_v(double t) {
   return -gamma / 2 * exp(-gamma / 2 * t) * cos(sqrt(omega * omega - gamma * gamma / 4.) * t) - exp(-gamma / 2 * t) * sin(sqrt(omega * omega - gamma * gamma / 4.) * t) * sqrt(omega * omega - gamma * gamma / 4.);
}

double error_x(const auto& result_t_x) {
   double acum_error = 0;
   for (const auto& tv : result_t_x)
      acum_error += std::abs(tv[1] - solution_x(tv[0]));
   return acum_error;
}

int main() {
   int N = 50;                  // number of steps
   double dt = 0.02;            // time step
   double t0 = 0;               // initial time
   double x0 = 1;               // initial position
   double v0 = solution_v(t0);  // initial velocity
   double a0 = 0;               // initial acceleration

   /* -------------------------------------------------------------------------- */
   // Backward Euler
   // \label{DampedHrmonicOscillator:BackwardEuler}
   auto result_BKE = [&](auto N) {
      double dt = 1. / N;
      std::vector<std::vector<double>> BKE_t_x({{t0, x0}}), BKE_t_v({{t0, v0}});
      {
         double v = v0, x = x0, t = t0;
         for (auto j = 0; j < N; j++) {
            //  double x_next_appx = x + dt * v;x
            //  double v_next_appx = v + dt * acceleration(t + dt, x_next_appx);
            double v_next_appx = solution_v(t + dt);
            x += dt * v_next_appx;
            v += dt * acceleration(t + dt, x);
            t += dt;
            BKE_t_x.push_back({t, x});
            BKE_t_v.push_back({t, v});
         }
      }
      return BKE_t_x;
   };
   /* -------------------------------------------------------------------------- */
   // LeapFrog
   // \label{DampedHrmonicOscillator:LeapFrog}
   auto result_LPFG = [&](auto N_IN) {
      auto N = (int)(N_IN / 2);
      double dt = 1. / N;
      std::vector<std::vector<double>> LPFG_t_x({{t0, x0}}), LPFG_t_v({{t0, v0}});
      LeapFrog LPFG(dt, t0, x0, v0);
      for (auto j = 0; j < 2 * N; j++) {
         double t = LPFG.get_t();
         double x = LPFG.get_x();
         double v = LPFG.get_v();
         LPFG_t_x.push_back({t, x});
         LPFG_t_v.push_back({t, v});
         LPFG.push(acceleration(x, v));
      }
      return LPFG_t_x;
   };
   /* -------------------------------------------------------------------------- */
   // Runge-Kutta
   // \label{DampedHrmonicOscillator:RungeKutta}
   const int order = 4;
   auto result_RK = [&](auto N_IN) {
      auto N = (int)(N_IN / order);
      double dt = 1. / N;
      std::vector<std::vector<double>> RK_t_v({{t0, v0}}), RK_t_x({{t0, x0}});
      double t = t0, x = x0, v = v0;
      for (auto j = 0; j < N; j++) {
         RungeKutta RK_x(dt, t, x, order);
         RungeKutta RK_v(dt, t, v, order);
         do {
            t = RK_x.get_t();
            x = RK_x.get_x();
            v = RK_v.get_x();
            RK_x.push(v);
            RK_v.push(acceleration(x, v));
         } while (!RK_x.finished);
         t = RK_x.get_t();
         x = RK_x.get_x();
         v = RK_v.get_x();
         RK_t_x.push_back({t, x});
         RK_t_v.push_back({t, v});
      }
      return RK_t_x;
   };
   /* -------------------------------------------------------------------------- */
   // Exact solution
   std::vector<std::vector<double>> exact_x, exact_v;
   for (auto j = 0; j < 10 * N; j++) {
      double t = j * 1. / (10. * N);
      exact_x.push_back({t, solution_x(t)});
      exact_v.push_back({t, solution_v(t)});
   }
   /* -------------------------------------------------------------------------- */
   {
      // plot time vs position
      GNUPLOT plot;
      plot.Set({{"key", "left"}});
      plot.SaveData(exact_x, {{"lc", plot.rgb({255, 0, 0})}, {"w", "l"}, {"lw", "2"}, {"title", "exact x"}});
      plot.SaveData(result_BKE(N), {{"lc", plot.rgb({100, 0, 255})}, {"w", "lp"}, {"lw", ".5"}, {"title", "Backward Euler"}});
      plot.SaveData(result_LPFG(N), {{"lc", plot.rgb({0, 100, 255})}, {"w", "lp"}, {"lw", ".5"}, {"title", "LeapFrog"}});
      plot.SaveData(result_RK(N), {{"lc", plot.rgb({255, 100, 0})}, {"w", "lp"}, {"lw", ".5"}, {"title", "Runge-Kutta" + std::to_string(order)}});
      plot.Plot2D();
      std::cin.ignore();
   }
   /* -------------------------------------------------------------------------- */
   // plot error
   {
      std::vector<std::vector<double>> error_BKE, error_LPFG, error_RK;
      for (auto N = 10; N < 100; ++N) {
         auto tmp = result_BKE(N);
         error_BKE.push_back({(double)N, error_x(tmp) / tmp.size()});
         tmp = result_LPFG(N);
         error_LPFG.push_back({(double)N, error_x(tmp) / tmp.size()});
         tmp = result_RK(N);
         error_RK.push_back({(double)N, error_x(tmp) / tmp.size()});
      }

      GNUPLOT plot;
      plot.Set({{"key", "left"}, {"logscale", "y"}});
      plot.SaveData(error_BKE, {{"lc", plot.rgb({100, 0, 255})}, {"w", "lp"}, {"lw", ".5"}, {"title", "Backward Euler"}});
      plot.SaveData(error_LPFG, {{"lc", plot.rgb({0, 100, 255})}, {"w", "lp"}, {"lw", ".5"}, {"title", "LeapFrog"}});
      plot.SaveData(error_RK, {{"lc", plot.rgb({255, 100, 0})}, {"w", "lp"}, {"lw", ".5"}, {"title", "Runge-Kutta" + std::to_string(order)}});
      plot.Plot2D();
      std::cin.ignore();
   }
};
