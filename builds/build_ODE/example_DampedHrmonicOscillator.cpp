#include "GNUPLOT.hpp"
#include "integrationOfODE.hpp"

/*DOC_EXTRACT ODE

## 減衰調和振動子/Damped Harmonic Oscillator

$m * \frac{d^2x}{dt^2} + b * \frac{dx}{dt} + k * x = 0$

| ![](example_DampedHrmonicOscillator.png) | ![](example_DampedHrmonicOscillator_last.png) |
|:---:|:---:|

*/

const double m = 1;
const double gamma = 1;
const double omega = 10;
const double b = gamma * m;
const double k = omega * omega * m;

double acceleration(double x, double v) {
   return -(b * v + k * x) / m;
}

double solution_x(double t) {
   return exp(-gamma / 2 * t) * cos(sqrt(omega * omega - gamma * gamma / 4.) * t);
}

double solution_v(double t) {
   return -gamma / 2 * exp(-gamma / 2 * t) * cos(sqrt(omega * omega - gamma * gamma / 4.) * t) - exp(-gamma / 2 * t) * sin(sqrt(omega * omega - gamma * gamma / 4.) * t) * sqrt(omega * omega - gamma * gamma / 4.);
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
   std::vector<std::vector<double>> BKE_t_v, BKE_t_x;
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
   /* -------------------------------------------------------------------------- */
   // LeapFrog
   std::vector<std::vector<double>> LPFG_t_v, LPFG_t_x;
   LeapFrog LPFG(dt, t0, x0, v0);
   for (auto j = 0; j < 2 * N; j++) {
      double t = LPFG.get_t();
      double x = LPFG.get_x();
      double v = LPFG.get_v();
      LPFG_t_x.push_back({t, x});
      LPFG_t_v.push_back({t, v});
      LPFG.push(acceleration(x, v));
   }
   /* -------------------------------------------------------------------------- */
   // Runge-Kutta
   const int order = 4;
   std::vector<std::vector<double>> RK_t_v, RK_t_x;
   double t = t0;
   double x = x0;
   double v = v0;
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
   /* -------------------------------------------------------------------------- */
   // Exact solution
   std::vector<std::vector<double>> exact_x, exact_v;
   for (auto j = 0; j < 10 * N; j++) {
      double t = j * dt / 10;
      exact_x.push_back({t, solution_x(t)});
      exact_v.push_back({t, solution_v(t)});
   };
   /* -------------------------------------------------------------------------- */
   GNUPLOT plot;
   plot.Set({{"key", "left"}});
   plot.SaveData(exact_x, {{"lc", plot.rgb({255, 0, 0})}, {"w", "l"}, {"lw", "2"}, {"title", "exact x"}});
   plot.SaveData(BKE_t_x, {{"lc", plot.rgb({100, 0, 255})}, {"w", "lp"}, {"lw", ".5"}, {"title", "Backward Euler"}});
   plot.SaveData(LPFG_t_x, {{"lc", plot.rgb({0, 100, 255})}, {"w", "lp"}, {"lw", ".5"}, {"title", "LeapFrog"}});
   plot.SaveData(RK_t_x, {{"lc", plot.rgb({255, 100, 0})}, {"w", "lp"}, {"lw", ".5"}, {"title", "Runge-Kutta" + std::to_string(order)}});
   plot.Plot2D();
   std::cin.ignore();
};
