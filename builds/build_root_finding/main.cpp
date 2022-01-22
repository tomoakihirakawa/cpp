#include "GNUPLOT.hpp"
#include "fundamental.hpp"
#include "rootFinding.hpp"

// #define simple_test
#ifdef simple_test
double f(const double x)
{
	return sin(x);
}
double dfdx(const double x)
{
	return cos(x);
}

V_d F(const V_d &X)
{
	return {sin(X[0]), cos(X[1]), sin(X[2])};
}
VV_d dFdx(const V_d &X)
{
	return {{cos(X[0]), 0, 0},
			{0, -sin(X[1]), 0},
			{0, 0, cos(X[2])}};
}

int main()
{
	/* ----------------------- 一次元の場合 ----------------------- */
	double xn = 1.;
	double xn1 = 0.;
	for (auto i = 0; i < 10; i++)
	{
		Print(xn);
		xn1 = xn - f(xn) / dfdx(xn);
		xn = xn1;
	}
	/* ----------------------- 多次元の場合 ----------------------- */
	V_d Xn = {1., 1., 1.};
	V_d Xn1(3), dX(3);
	for (auto i = 0; i < 10; i++)
	{
		Print(Xn);
		ludcmp lu(dFdx(Xn));
		lu.solve(-F(Xn), dX);
		Xn1 = Xn + dX;
		Xn = Xn1;
	}
};
#else

// class NewtonRaphson
// {
// public:
// 	V_d X;
// 	V_d dX; // tmp
// 	NewtonRaphson(const V_d &Xinit) : X(Xinit), dX(Xinit){};
// 	void update(const V_d &F, const VV_d &dFdx)
// 	{
// 		ludcmp lu(dFdx);
// 		lu.solve(-F, dX);
// 		X += dX;
// 	};
// };

double f(const double x)
{
	return sin(x);
}
double dfdx(const double x)
{
	return cos(x);
}

Tddd F(const Tddd &X)
{
	return {sin(std::get<0>(X)), cos(std::get<1>(X)), sin(std::get<2>(X))};
}
T3Tddd dFdx(const Tddd &X)
{
	return {{cos(std::get<0>(X)), 0, 0},
			{0, -sin(std::get<1>(X)), 0},
			{0, 0, cos(std::get<2>(X))}};
}

int main()
{
	/* ----------------------- 1次元の場合 ----------------------- */
	std::cout << "1次元の場合" << std::endl;
	{
		NewtonRaphson nr(4.2);
		for (auto i = 0; i < 10; i++)
		{
			nr.update(f(nr.X), dfdx(nr.X));
			Print(Norm(nr.dX), red);
			Print(Norm(nr.X), red);
		}
	}
	/* ----------------------- 多次元の場合 ----------------------- */
	std::cout << "多次元の場合" << std::endl;
	{
		NewtonRaphson nr(Tddd{1., 1., 1.});
		for (auto i = 0; i < 10; i++)
		{
			nr.update(F(nr.X), dFdx(nr.X));
			Print(Norm(nr.dX), red);
			Print(Norm(nr.X), red);
		}
	}
};
#endif