#include "GNUPLOT.hpp"
#include "minMaxOfFunctions.hpp"

int main() {
	VV_d A = {{4., 1.}, {1., 3.}};
	GradientDescent gd(A);
	{
		TimeWatch tm;
		auto x = gd.solve({1., 2.});
		Print(x);
		Print(gd.count);
		std::cout << tm.get()[0].count() << std::endl;
	}

	{
		TimeWatch tm;
		auto x = gd.solveCG({1., 2.});
		Print(x);
		Print(gd.count);
		std::cout << tm.get()[0].count() << std::endl;
	}
};
