//
//
#include "GNUPLOT.hpp"
#include "integrationOfODE.hpp"

double dydt(double y, double t)
{
	return sin(t) * sin(t) * y;
};
double solution(double y0, double t)
{
	return y0 * exp((2 * t - sin(2. * t)) / 4.);
};

int main()
{
	std::vector<std::vector<double>> ansRK2, ansRK3, ansRK4, exact;

	double y0 = 2, t0 = 0.; //初期値
	double dt = 1.;			//時間ステップ
	double t_end = 5;		//終了時刻

	{
		Print("2次のルンゲクッタ");
		double y = y0;
		double t = t0;
		for (auto j = 0; j < 100; j++)
		{
			// タプルでもできることのチェック
			Tddd Y = {y, y, y};
			RungeKutta_ rk(dt, t, Y, 2);
			ansRK2.push_back({t, std::get<0>(Y)});
			while (true)
			{
				rk.displayStatus();
				//!ここが境界値問題を解いて求めるところ
				/* ------------------------ 微分を評価 ----------------------- */
				auto tmp = dydt(std::get<0>(rk.getX()), rk.gett());
				Tddd dYdT = {tmp, tmp, tmp};
				if (rk.push(dYdT))
					break;
			}
			y = std::get<0>(rk.getX());
			t += dt;
			Print(y, Magenta);
			if (t > t_end)
				break;
		}
	}
	{
		Print("3次のルンゲクッタ");
		double y = y0;
		double t = t0;
		for (auto j = 0; j < 100; j++)
		{
			RungeKutta3 rk(dt, t, y);
			ansRK3.push_back({t, y});
			while (true)
			{
				rk.displayStatus();
				if (rk.push(dydt(rk.getX()[0], rk.gett())))
					break;
			}
			y = rk.getX()[0];
			t += dt;
			Print(y, Magenta);
			if (t > t_end)
				break;
		}
	}
	{
		Print("4次のルンゲクッタ");
		double y = y0;
		double t = t0;
		for (auto j = 0; j < 100; j++)
		{
			RungeKutta4 rk(dt, t, y);
			ansRK4.push_back({t, y});
			while (true)
			{
				rk.displayStatus();
				if (rk.push(dydt(rk.getX()[0], rk.gett())))
					break;
			}
			y = rk.getX()[0];
			t += dt;
			Print(y, Magenta);
			if (t > t_end)
				break;
		}
	}

	for (auto j = 0; j < 1000; j++)
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
	plot.Plot2D();
	std::cin.ignore();
};
