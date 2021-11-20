
#include "fundamental.hpp"
#include "./matrix_to_be_solved.hpp"

int main()
{
	TimeWatch watch;
	std::cout << Magenta << watch.get() << reset << std::endl;
	ludcmp lu(mat);
	std::cout << Magenta << watch.get() << reset << std::endl;
	double total = 0.;
	for (const auto &n : (lu.Inverse() - inv))
		total += Norm(n);

	std::cout << total << std::endl;
	std::cout << Magenta << watch.get() << reset << std::endl;
}
