#include "GNUPLOT.hpp"
#include "fundamental.hpp"

double dydt(double y, double t)
{
  return sin(t) * sin(t) * y;
};
double solution(double y0, double t)
{
  return y0 * exp((2 * t - sin(2. * t)) / 4.);
};

using VV_d = std::vector<std::vector<double>>;

int main()
{

  VV_d ansRK2, ansRK3, ansRK4, exact;

  double y0 = 2, t0 = 0.; //初期値
  double dt = 1.;         //時間ステップ
  double t_end = 15;      //終了時刻

  int max_step = 150;
  { //2次のルンゲクッタ
    RK2 rk(dt);
    double y = y0;
    double t = t0;
    for (auto j = 0; j < max_step; j++)
    {
      ansRK2.push_back({t, y});
      while (true)
      {
        bool finish = rk.improve(dydt(y + rk.dx, t + rk.dt));
        rk.displayStatus();
        if (finish)
          break;
      }
      y += rk.getImproved();
      t += dt;
      rk.Clear();
      Print(y, Magenta);
      if (t > t_end)
        break;
    }
  }
  { //3次のルンゲクッタ
    RK3 rk(dt);
    double y = y0;
    double t = t0;
    for (auto j = 0; j < max_step; j++)
    {
      ansRK3.push_back({t, y});
      while (true)
      {
        bool finish = rk.improve(dydt(y + rk.dx, t + rk.dt));
        rk.displayStatus();
        if (finish)
          break;
      }
      y += rk.getImproved();
      t += dt;
      rk.Clear();
      Print(y, Magenta);
      if (t > t_end)
        break;
    }
  }
  { //4次のルンゲクッタ
    RK4 rk(dt);
    double y = y0;
    double t = t0;
    for (auto j = 0; j < max_step; j++)
    {
      ansRK4.push_back({t, y});
      while (true)
      {
        bool finish = rk.improve(dydt(y + rk.dx, t + rk.dt));
        rk.displayStatus();
        if (finish)
          break;
      }
      y += rk.getImproved();
      t += dt;
      rk.Clear();
      Print(y, Magenta);
      if (t > t_end)
        break;
    }
  }

  ////////////////////////////////////////////////////////////////////
  VV_d ansTayler1, ansTayler2, ansTayler3;
  { //テイラー展開
    double y = y0;
    double t = t0;
    for (auto j = 0; j < max_step; j++)
    {
      ansTayler1.push_back({t, y});
      y += dt * dydt(y, t);
      t += dt;
      if (t > t_end)
        break;
    }
  }

  auto dy2dt2 = [](const double y, const double t) {
    return 2. * sin(t) * cos(t) * y + sin(t) * sin(t) * dydt(y, t);
  };
  auto dy3dt3 = [&dy2dt2](const double y, const double t) {
    return -sin(2. * t) * 2. * y + sin(2. * t) * dydt(y, t) + 2. * cos(t) * sin(t) * dydt(y, t) + sin(t) * sin(t) * dy2dt2(y, t);
  };

  { //テイラー展開
    double y = y0;
    double t = t0;
    for (auto j = 0; j < max_step; j++)
    {
      ansTayler2.push_back({t, y});
      y += dt * dydt(y, t) + dt * dt / 2. * dy2dt2(y, t);
      t += dt;
      if (t > t_end)
        break;
    }
  }

  { //テイラー展開
    double y = y0;
    double t = t0;
    for (auto j = 0; j < max_step; j++)
    {
      ansTayler3.push_back({t, y});
      y += dt * dydt(y, t) + dt * dt / 2. * dy2dt2(y, t) + dt * dt * dt / 6. * dy3dt3(y, t);
      t += dt;
      if (t > t_end)
        break;
    }
  }

  ////////////////////////////////////////////////////////////////////////////

  VV_d ansAB2, ansAB3;

  { //2次のアダムスバッシュフォース
    double y = y0;
    double t = t0;
    double dydt_ = 0.;
    double pre_y, pre_t;
    for (auto j = 0; j < max_step; j++)
    {
      ansAB2.push_back({t, y});

      auto tmp_y = y, tmp_t = t;

      if (j == 0)
        y += dt * dydt(y, t);
      else
        y += dt / 2. * (3. * dydt(y, t) - dydt(pre_y, pre_t));

      t += dt;
      if (t > t_end)
        break;

      pre_y = tmp_y;
      pre_t = tmp_t;
    }
  }

  { //3次のアダムスバッシュフォース
    double y = y0;
    double t = t0;
    double dydt_ = 0.;
    double pre_y, pre_t;
    double pre2_y, pre2_t;
    for (auto j = 0; j < max_step; j++)
    {
      ansAB3.push_back({t, y});

      auto tmp2_y = pre_y, tmp2_t = pre_t;
      auto tmp_y = y, tmp_t = t;

      if (j == 0)
        y += dt * dydt(y, t);
      else if (j == 1)
        y += dt / 2. * (3. * dydt(y, t) - dydt(pre_y, pre_t));
      else
        y += dt / 12. * (23. * dydt(y, t) - 16. * dydt(pre_y, pre_t) + 5. * dydt(pre2_y, pre2_t));

      t += dt;
      if (t > t_end)
        break;

      pre2_y = tmp2_y;
      pre2_t = tmp2_t;
      pre_y = tmp_y;
      pre_t = tmp_t;
    }
  }

  for (auto j = 0; j < max_step * 2; j++)
  {
    double t = j * 0.05;
    if (t > t_end)
      break;
    exact.push_back({t, solution(y0, t)});
  };

  GNUPLOT plot;
  plot.Set({{"key", "left"}});
  plot.SaveData(exact, {{"lc", plot.rgb({255, 0, 0})}, {"w", "l"}, {"lw", "4"}, {"title", "exact"}});
  plot.SaveData(ansRK4, {{"lc", plot.rgb({205, 0, 205})}, {"w", "lp"}, {"lw", "2"}, {"title", "RK4"}});
  plot.SaveData(ansRK3, {{"lc", plot.rgb({0, 205, 205})}, {"w", "lp"}, {"lw", "2"}, {"title", "RK3"}});
  plot.SaveData(ansRK2, {{"lc", plot.rgb({205, 205, 0})}, {"w", "lp"}, {"lw", "2"}, {"title", "RK2"}});
  plot.SaveData(ansTayler1, {{"lc", plot.rgb({105, 105, 0})}, {"w", "lp"}, {"lw", "2"}, {"title", "Tayler1"}});
  plot.SaveData(ansTayler2, {{"lc", plot.rgb({205, 105, 80})}, {"w", "lp"}, {"lw", "2"}, {"title", "Tayler2"}});
  plot.SaveData(ansTayler3, {{"lc", plot.rgb({155, 105, 180})}, {"w", "lp"}, {"lw", "2"}, {"title", "Tayler3"}});
  plot.SaveData(ansAB2, {{"lc", plot.rgb({55, 205, 230})}, {"w", "lp"}, {"lw", "2"}, {"title", "Adams-Bashforth2"}});
  plot.SaveData(ansAB3, {{"lc", plot.rgb({155, 205, 230})}, {"w", "lp"}, {"lw", "2"}, {"title", "Adams-Bashforth3"}});
  plot.Plot2D();
  std::cin.ignore();
};
