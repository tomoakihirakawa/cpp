/*DOC_EXTRACT ODE
# Runge-Kutta Integration of ODE
This C++ program demonstrates the application of various Runge-Kutta methods (first to fourth order) for solving a first-order ordinary differential equation (ODE).
![](res.png)
*/

#include "GNUPLOT.hpp"
#include "integrationOfODE.hpp"

double dydt(double y, double t) { return sin(t) * sin(t) * y; };
double solution(double y0, double t) { return y0 * exp((2 * t - sin(2. * t)) / 4.); };

int main() {
   int n = 4;
   std::vector<std::vector<std::vector<double>>> ansRK(n);
   std::array<double, 1> y0 = {2};
   double t0 = 0.;    // 初期値
   double dt = 1.;    // 時間ステップ
   double t_end = 5;  // 終了時刻
   for (auto i = 1; i <= 4; ++i) {
      std::cout << i << "次のルンゲクッタ" << std::endl;
      std::array<double, 1> y = y0;
      double t = t0;
      for (auto j = 0; j < 100; j++) {
         RungeKutta rk(dt, t, y, i);
         ansRK[i - 1].push_back({t, y[0]});
         while (true) {
            rk.displayStatus();
            if (rk.push({dydt(rk.getX()[0], rk.gett())}))
               break;
         }
         y = rk.getX();
         t += dt;
         Print(y[0], Magenta);
         if (t > t_end)
            break;
      }
   }

   std::vector<std::vector<double>> exact;
   for (auto j = 0; j < 1000; j++) {
      double t = j * 0.05;
      if (t > t_end)
         break;
      exact.push_back({t, solution(y0[0], t)});
   };

   GNUPLOT plot;
   plot.Set({{"key", "left"}});
   plot.SaveData(exact, {{"lc", plot.rgb({255, 0, 0})}, {"w", "l"}, {"lw", "4"}, {"title", "exact"}});
   plot.SaveData(ansRK[3], {{"lc", plot.rgb({205, 0, 205})}, {"w", "lp"}, {"lw", "2"}, {"title", "RK4"}});
   plot.SaveData(ansRK[2], {{"lc", plot.rgb({0, 205, 205})}, {"w", "lp"}, {"lw", "2"}, {"title", "RK3"}});
   plot.SaveData(ansRK[1], {{"lc", plot.rgb({205, 205, 0})}, {"w", "lp"}, {"lw", "2"}, {"title", "RK2"}});
   plot.SaveData(ansRK[0], {{"lc", plot.rgb({105, 205, 0})}, {"w", "lp"}, {"lw", "2"}, {"title", "RK1"}});
   plot.Plot2D();
   std::cin.ignore();
};
