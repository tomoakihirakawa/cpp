#ifndef integrationOfODE_H
#define integrationOfODE_H
#pragma once

#include "fundamental.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

//* ----------------- 新しく作ったルンゲクッタのプログラム2021/05/01 ----------------- */
//* ベクトルに拡張した
//* 注意：入力はスカラーでもいいが，getXはベクトルしか返せない

template <typename T>
struct RungeKuttaCommon
{
public:
	double dt_fixed; //基本となるステップ間隔,固定
	/* ----------------------- 毎回変わる値 ----------------------- */
	double dt;			//次回計算する必要がある値になる増分,回数ごとに変化
	std::vector<T> _dX; //これにプッシュしていく，c1k1,c2k2,c3k3,c4k4,,,//このベクトルの長さで何回目かを判断数し，必要な回数まできたら，finished=trueとなる
	T dX;				//_dXやdtを基に計算される
	/* ------------------------- 初期値 ------------------------ */
	double t_init; //初期値,固定

public:
	T Xinit;
	/* ------------------------ 積分結果 ------------------------ */
	bool finished; //全てそろって積分が評価できるようになったらtrue
				   //最終結果はdXに代入される;

public:
	void displayStatus()
	{
		std::cout << "dt : " << red << this->dt << reset << std::endl;
		std::cout << "_dX : " << red << this->_dX << reset << std::endl;
		std::cout << "dX : " << red << this->dX << reset << std::endl;
	};
	int steps;
	RungeKuttaCommon(const double dt_IN, const double t0, const T &X0 /*initial value*/, int stepsIN)
		: dt_fixed(dt_IN), dt(0), t_init(t0), Xinit(X0), _dX({}), dX(X0), steps(stepsIN){};
	~RungeKuttaCommon(){};

	T getXinit() { return this->Xinit; };

	T getX() { return this->Xinit + this->dX; };

	T getdX() { return this->dX; };
	double gett() { return this->t_init + this->dt; };
	double getdt() { return this->dt; };
	bool push(const T &dXdt_IN)
	{
		_dX.insert(_dX.begin(), dXdt_IN * dt_fixed);
		if (this->steps == 2)
		{
			switch (_dX.size())
			{
			case 1: // -> {f1(t,v(t))}
				dX = dXdt_IN * (dt = dt_fixed);
				return this->finished = false;
			case 2: // -> {f3, f2, f1(t,v(t))}
				dX = (_dX[0] + _dX[1]) / 2.;
				dt = dt_fixed;
				return this->finished = true;
			default:
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
			}
		}
		else if (this->steps == 3)
		{
			switch (_dX.size())
			{
			case 1: // -> {f1(t,v(t))}
				dX = dXdt_IN * (dt = dt_fixed / 2.);
				return finished = false;
			case 2: // -> {f2, f1(t,v(t))}
				dX = dXdt_IN * (dt = dt_fixed);
				return finished = false;
			case 3: // -> {f3, f2, f1(t,v(t))}
				dX = (_dX[0] + 4. * _dX[1] + _dX[2]) / 6.;
				dt = dt_fixed;
				return finished = true;
			default:
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
			}
		}
		else if (this->steps == 4)
		{
			switch (_dX.size())
			{
			case 1: // -> {f1(t,v(t))}
				dX = dXdt_IN * (dt = dt_fixed / 2.);
				return finished = false;
			case 2: // -> {f2, f1(t,v(t))}
				dX = dXdt_IN * (dt = dt_fixed / 2.);
				return finished = false;
			case 3: // -> {f3, f2, f1(t,v(t))}
				dX = dXdt_IN * (dt = dt_fixed);
				return finished = false;
			case 4: // -> dXdt[1][] = {f4, f3, f2, f1(t,v(t))}
				dX = (_dX[0] + 2. * _dX[1] + 2. * _dX[2] + _dX[3]) / 6.;
				dt = dt_fixed;
				return finished = true;
			default:
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
			}
		}
		else
		{
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
		}
	};
};

template <typename T>
struct RungeKutta_ : public RungeKuttaCommon<T>
{
	RungeKutta_(const double dt_IN, const double t0, const T &X0 /*initial value*/, int stepsIN)
		: RungeKuttaCommon<T>(dt_IN, t0, X0, stepsIN){};
};
template <>
struct RungeKutta_<std::vector<double>> : public RungeKuttaCommon<std::vector<double>>
{
	RungeKutta_<std::vector<double>>(const double dt_IN, const double t0, const std::vector<double> &X0 /*initial value*/, int stepsIN)
		: RungeKuttaCommon<std::vector<double>>(dt_IN, t0, X0, stepsIN)
	{
		this->dX = V_d(X0.size(), 0.);
	};
};
template <>
struct RungeKutta_<double> : public RungeKuttaCommon<double>
{
	RungeKutta_<double>(const double dt_IN, const double t0, const double &X0 /*initial value*/, int stepsIN)
		: RungeKuttaCommon<double>(dt_IN, t0, X0, stepsIN)
	{
		this->dX = 0.;
	};
};
template <>
struct RungeKutta_<Tddd> : public RungeKuttaCommon<Tddd>
{
	RungeKutta_<Tddd>(const double dt_IN, const double t0, const Tddd &X0 /*initial value*/, int stepsIN)
		: RungeKuttaCommon<Tddd>(dt_IN, t0, X0, stepsIN)
	{
		this->dX = {0, 0, 0};
	};
};
/* ------------------------------------------------------ */
class RungeKutta
{
protected:
	double dt_fixed; //基本となるステップ間隔,固定
	/* ----------------------- 毎回変わる値 ----------------------- */
	double dt; //次回計算する必要がある値になる増分,回数ごとに変化
	VV_d _dX;  //これにプッシュしていく，c1k1,c2k2,c3k3,c4k4,,,//このベクトルの長さで何回目かを判断数し，必要な回数まできたら，finished=trueとなる
	V_d dX;	   //_dXやdtを基に計算される
	/* ------------------------- 初期値 ------------------------ */
	double t_init; //初期値,固定

public:
	V_d Xinit;
	/* ------------------------ 積分結果 ------------------------ */
	bool finished; //全てそろって積分が評価できるようになったらtrue
				   //最終結果はdXに代入される;

public:
	void displayStatus()
	{
		std::cout << "dt : " << red << this->dt << reset << std::endl;
		std::cout << "_dX : " << red << this->_dX << reset << std::endl;
		std::cout << "dX : " << red << this->dX << reset << std::endl;
	};
	RungeKutta(const double dt_IN, const double t0, const V_d &X0 /*initial value*/)
		: dt_fixed(dt_IN), dt(0), t_init(t0), Xinit(X0), _dX({}), dX(X0.size(), 0.){};
	RungeKutta(const double dt_IN, const double t0, const double X0 /*initial value*/)
		: dt_fixed(dt_IN), dt(0), t_init(t0), Xinit({X0}), _dX({}), dX(1, 0.){};
	~RungeKutta(){};
	V_d getX() { return this->Xinit + this->dX; };
	V_d getdX() { return this->dX; };
	double gett() { return this->t_init + this->dt; };
	double getdt() { return this->dt; };
	virtual bool push(const V_d &dXdt_IN) = 0;
	virtual bool push(const double dXdt_IN) = 0;
};
class RungeKutta2 : public RungeKutta
{
public:
	RungeKutta2(const double dt_IN, const double t0, const V_d &x0) : RungeKutta(dt_IN, t0, x0){};
	RungeKutta2(const double dt_IN, const double t0, const double x0) : RungeKutta(dt_IN, t0, x0){};
	bool push(const V_d &dXdt_IN) override
	{
		_dX.insert(_dX.begin(), dXdt_IN * dt_fixed);
		switch (_dX.size())
		{
		case 1: // -> {f1(t,v(t))}
			dX = dXdt_IN * (dt = dt_fixed);
			return this->finished = false;
		case 2: // -> {f3, f2, f1(t,v(t))}
			dX = (_dX[0] + _dX[1]) / 2.;
			dt = dt_fixed;
			return this->finished = true;
		default:
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
		}
	};
	bool push(const double dXdt_IN) override { return push(V_d{dXdt_IN}); /*スカラーとして微分をプッシュしようとした場合ベクトルに直しプッシュ*/ };
};
class RungeKutta3 : public RungeKutta
{
public:
	RungeKutta3(const double dt_IN, const double t0, const V_d &x0) : RungeKutta(dt_IN, t0, x0){};
	RungeKutta3(const double dt_IN, const double t0, const double x0) : RungeKutta(dt_IN, t0, x0){};
	bool push(const V_d &dXdt_IN) override
	{
		_dX.insert(_dX.begin(), dXdt_IN * dt_fixed);
		switch (_dX.size())
		{
		case 1: // -> {f1(t,v(t))}
			dX = dXdt_IN * (dt = dt_fixed / 2.);
			return finished = false;
		case 2: // -> {f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = dt_fixed);
			return finished = false;
		case 3: // -> {f3, f2, f1(t,v(t))}
			dX = (_dX[0] + 4. * _dX[1] + _dX[2]) / 6.;
			dt = dt_fixed;
			return finished = true;
		default:
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
		}
	};
	bool push(const double dXdt_IN) override { return push(V_d{dXdt_IN}); /*スカラーとして微分をプッシュしようとした場合ベクトルに直しプッシュ*/ };
};
class RungeKutta4 : public RungeKutta
{
public:
	RungeKutta4(const double dt_IN, const double t0, const V_d &x0) : RungeKutta(dt_IN, t0, x0){};
	RungeKutta4(const double dt_IN, const double t0, const double x0) : RungeKutta(dt_IN, t0, x0){};
	bool push(const V_d &dXdt_IN) override
	{
		_dX.insert(_dX.begin(), dXdt_IN * dt_fixed);
		switch (_dX.size())
		{
		case 1: // -> {f1(t,v(t))}
			dX = dXdt_IN * (dt = dt_fixed / 2.);
			return finished = false;
		case 2: // -> {f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = dt_fixed / 2.);
			return finished = false;
		case 3: // -> {f3, f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = dt_fixed);
			return finished = false;
		case 4: // -> dXdt[1][] = {f4, f3, f2, f1(t,v(t))}
			dX = (_dX[0] + 2. * _dX[1] + 2. * _dX[2] + _dX[3]) / 6.;
			dt = dt_fixed;
			return finished = true;
		default:
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．");
		}
	};
	bool push(const double dXdt_IN) override { return push(V_d{dXdt_IN}); /*スカラーとして微分をプッシュしようとした場合ベクトルに直しプッシュ*/ };
};

/* ------------------- 古いるんげくったのプログラム ------------------- */

class derivativeImprover
{
protected:
	double _dt;
	V_d _dx;
	VV_d _dX;

public:
	double improved;
	bool finished;
	double dt, dx;
	double t_init, x_init;
	V_d dX, X_init, improved_dX;

	derivativeImprover(const double &dt_IN)
		: _dx(0), _dt(dt_IN), dt(0), dx(0), t_init(0), x_init(0), _dX(0), dX(0){};
	derivativeImprover(const double &dt_IN, const double &t_init_IN, const double &x_init_IN)
		: _dx(0), _dt(dt_IN), dt(0), dx(0), t_init(t_init_IN), x_init(x_init_IN), _dX(0), dX(0){};
	derivativeImprover(const double &dt_IN, const double &t_init_IN, const V_d &X_init_IN)
		: _dx(0), _dt(dt_IN), dt(0), dx(0), t_init(t_init_IN), X_init(X_init_IN), _dX(0), dX(X_init_IN.size(), 0.), improved_dX(X_init_IN.size(), 0.){};

	~derivativeImprover(){
		// Print("destruct derivativeImprover", red);
	};

	void Clear()
	{
		finished = false;
		_dx.clear();
	};

	//!dxです
	double getImproved()
	{
		if (finished)
			return improved;
		else
			return 0.;
	};

	V_d getImproved_dX()
	{
		if (finished)
			return improved_dX;
		else
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};

	double getImproved_x()
	{
		if (finished)
			return this->x_init + this->improved;
		else
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};

	V_d getImproved_X()
	{
		if (finished)
			return this->X_init + this->improved_dX;
		else
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};

	virtual bool improve(const double &dxdt_IN) = 0;
	virtual bool improve(const V_d &dXdt_IN) = 0;

	void displayStatus()
	{
		std::cout << "step : " << Red << this->_dx.size() << reset << std::endl;
		std::cout << "dt : " << red << this->dt << reset << std::endl;
		std::cout << "dx : " << red << this->dx << reset << std::endl;
	};

	void pushFront_dx(const double &v)
	{
		_dx.insert(_dx.begin(), v);
	};

	void pushFront_dX(const V_d &V)
	{
		_dX.insert(_dX.begin(), V);
	};

	double get_t() { return this->t_init + this->dt; };
	double get_x() { return this->x_init + this->dx; };
	V_d get_X() { return this->X_init + this->dX /*初期は{0,0,0..}*/; };
};
class RK2 : public derivativeImprover
{
public:
	RK2(const double &dt_IN) : derivativeImprover(dt_IN){};
	//* ------------------------ スカラー ------------------------ */
	RK2(const double &dt_IN, const double &t_init, const double &x_init)
		: derivativeImprover(dt_IN, t_init, x_init){};

	bool improve(const double &dxdt_IN) override
	{
		pushFront_dx(dxdt_IN * _dt);
		switch (_dx.size())
		{
		case 1: // -> {f1(t,v(t))}
			dx = dxdt_IN * (dt = _dt);
			return this->finished = false;
		case 2: // -> {f3, f2, f1(t,v(t))}
			dx = (dt = 0.);
			improved = (_dx[0] + _dx[1]) / 2.;
			return this->finished = true;
		default:
			return false;
		}
	};
	//! ------------------------ ベクトル用 ----------------------- */
	RK2(const double &dt_IN, const double &t_init, const V_d &x_init)
		: derivativeImprover(dt_IN, t_init, x_init){};
	bool improve(const V_d &dXdt_IN) override
	{
		pushFront_dX(dXdt_IN * _dt);
		switch (_dX.size())
		{
		case 1: // -> {f1(t,v(t))}
			dX = dXdt_IN * (dt = _dt);
			return this->finished = false;
		case 2: // -> {f3, f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = 0.);
			improved_dX = (_dX[0] + _dX[1]) / 2.;
			return this->finished = true;
		default:
			return false;
		}
	};
};
class RK3 : public derivativeImprover
{
public:
	RK3(const double &dt_IN) : derivativeImprover(dt_IN){};
	//* ------------------------ スカラー ------------------------ */
	RK3(const double &dt_IN, const double &t_init, const double &x_init)
		: derivativeImprover(dt_IN, t_init, x_init){};
	bool improve(const double &dxdt_IN) override
	{
		pushFront_dx(dxdt_IN * _dt);
		switch (_dx.size())
		{
		case 1: // -> {f1(t,v(t))}
			dx = dxdt_IN * (dt = _dt / 2.);
			return finished = false;
		case 2: // -> {f2, f1(t,v(t))}
			dx = dxdt_IN * (dt = _dt);
			return finished = false;
		case 3: // -> {f3, f2, f1(t,v(t))}
			dx = (dt = 0.);
			improved = (_dx[0] + 4. * _dx[1] + _dx[2]) / 6.;
			return finished = true;
		default:
			return false;
		}
	};
	//! ------------------------ ベクトル用 ----------------------- */
	RK3(const double &dt_IN, const double &t_init, const V_d &x_init)
		: derivativeImprover(dt_IN, t_init, x_init){};
	bool improve(const V_d &dXdt_IN) override
	{
		pushFront_dX(dXdt_IN * _dt);
		switch (_dX.size())
		{
		case 1: // -> {f1(t,v(t))}
			dX = dXdt_IN * (dt = _dt / 2.);
			return finished = false;
		case 2: // -> {f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = _dt);
			return finished = false;
		case 3: // -> {f3, f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = 0.);
			improved_dX = (_dX[0] + 4. * _dX[1] + _dX[2]) / 6.;
			return finished = true;
		default:
			return false;
		}
	};
};
class RK4 : public derivativeImprover
{
public:
	RK4(const double &dt_IN) : derivativeImprover(dt_IN){};
	//* ------------------------ スカラー ------------------------ */
	RK4(const double &dt_IN, const double &t_init, const double &x_init)
		: derivativeImprover(dt_IN, t_init, x_init){};
	bool improve(const double &dxdt_IN) override
	{
		pushFront_dx(dxdt_IN * _dt);
		switch (_dx.size())
		{
		case 1: // -> {f1(t,v(t))}
			dx = dxdt_IN * (dt = _dt / 2.);
			return finished = false;
		case 2: // -> {f2, f1(t,v(t))}
			dx = dxdt_IN * (dt = _dt / 2.);
			return finished = false;
		case 3: // -> {f3, f2, f1(t,v(t))}
			dx = dxdt_IN * (dt = _dt);
			return finished = false;
		case 4: // -> dxdt[1][] = {f4, f3, f2, f1(t,v(t))}
			dx = (dt = 0.);
			improved = (_dx[0] + 2. * _dx[1] + 2. * _dx[2] + _dx[3]) / 6.;
			return finished = true;
		default:
			return false;
		}
	};
	//! ------------------------ ベクトル用 ----------------------- */
	RK4(const double &dt_IN, const double &t_init, const V_d &X_init)
		: derivativeImprover(dt_IN, t_init, X_init){};
	bool improve(const V_d &dXdt_IN) override
	{
		pushFront_dX(dXdt_IN * _dt);
		switch (_dX.size())
		{
		case 1: // -> {f1(t,v(t))}
			dX = dXdt_IN * (dt = _dt / 2.);
			return finished = false;
		case 2: // -> {f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = _dt / 2.);
			return finished = false;
		case 3: // -> {f3, f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = _dt);
			return finished = false;
		case 4: // -> dXdt[1][] = {f4, f3, f2, f1(t,v(t))}
			dX = dXdt_IN * (dt = 0.);
			improved_dX = (_dX[0] + 2. * _dX[1] + 2. * _dX[2] + _dX[3]) / 6.;
			return finished = true;
		default:
			return false;
		}
	};
};

class derivativeImprovers
{
public:
	double dt;
	int degree;
	V_d improved;
	std::vector<bool> finished;

	std::vector<derivativeImprover *> derImps;
	derivativeImprovers(const double dt_IN, const int degree_IN) : dt(dt_IN), degree(degree_IN), derImps({}), improved(degree_IN), finished(degree_IN){};

	~derivativeImprovers()
	{
		auto tmp = this->derImps;
		for (const auto d : tmp)
			delete d;
	};

	V_d dx(const std::vector<int> &vi) const
	{
		V_d ret({});
		for (const auto i : vi)
			ret.emplace_back(dx(i));
		return ret;
	};

	double dx(const int i) const
	{
		if (i > degree)
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

		return this->derImps[i]->dx;
	};

	std::vector<bool> improve(const V_d &dxdt)
	{
		std::vector<bool> ret(dxdt.size());
		if (dxdt.size() != degree)
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		else
		{
			for (auto i = 0; i < dxdt.size(); i++)
			{
				ret[i] = derImps[i]->improve(dxdt[i]);
				improved[i] = derImps[i]->getImproved();
				finished[i] = derImps[i]->finished;
			}
		}
		return ret;
	};

	bool isFinished() const
	{
		for (const auto f : finished)
			if (!f)
			{
				return false;
			}
		return true;
	};

	V_d getImproved(const std::vector<int> &vi) const
	{
		for (const auto f : finished)
			if (!f)
			{
				std::cout << "improver has not finished" << std::endl;
				return {};
			}

		V_d ret({});
		for (const auto i : vi)
			ret.emplace_back(this->improved[i]);
		return ret;
	};

	double getImproved(const int i) const
	{
		for (const auto f : finished)
			if (!f)
			{
				std::cout << "improver has not finished" << std::endl;
				return {};
			}

		return this->improved[i];
	};

	V_d getImproved()
	{
		for (const auto f : finished)
			if (!f)
			{
				std::cout << "improver has not finished" << std::endl;
				return {};
			}
		return improved;
	};
};

////////////////////////////////
class RK2s : public derivativeImprovers
{
public:
	RK2s(const double dt_IN, const int degree_IN) : derivativeImprovers(dt_IN, degree_IN)
	{
		for (auto i = 0; i < degree; i++)
			this->derImps.emplace_back(new RK2(dt));
	};
};
class RK3s : public derivativeImprovers
{
public:
	RK3s(const double dt_IN, const int degree_IN) : derivativeImprovers(dt_IN, degree_IN)
	{
		for (auto i = 0; i < degree; i++)
			this->derImps.emplace_back(new RK3(dt));
	};
};
class RK4s : public derivativeImprovers
{
public:
	RK4s(const double dt_IN, const int degree_IN) : derivativeImprovers(dt_IN, degree_IN)
	{
		for (auto i = 0; i < degree; i++)
			this->derImps.emplace_back(new RK4(dt));
	};
};

/*RungeKutta_code*/

#endif