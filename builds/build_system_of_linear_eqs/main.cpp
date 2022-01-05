#include "GNUPLOT.hpp"
#include "minMaxOfFunctions.hpp"

int main()
{
	VV_d A = {{4., 1.}, {1., 3.}};
	V_d x = {1., 2.};
	V_d b = {1., 2.};
	GradientMethod gd(A);
	{
		TimeWatch tm;
		auto x = gd.solve(b);
		Print(x);
		Print(gd.count);
		std::cout << "time = " << tm.get() << std::endl;
	}

	{
		TimeWatch tm;
		auto x = gd.solveCG(b);
		Print(x);
		Print(gd.count);
		std::cout << "time = " << tm.get() << std::endl;
	}

	{
		TimeWatch tm;
		ludcmp lu(A);
		lu.solve(b, x);
		std::cout << x << std::endl;
		std::cout << "time = " << tm.get() << std::endl;
	}
};
