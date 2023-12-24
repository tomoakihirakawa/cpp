#include <iostream>
#include <unordered_map>
#include <map>
#include <fundamental.hpp>
#include <utility>

int main()
{
	TimeWatch watch;
	std::cout << Green << "Elapsed time: " << Red << watch.get()[0].count() << reset << " s\n";
	int s = 1000000;
	/* ------------------------------------------------------ */
	std::vector<std::pair<int, double>> v(s);
	for (auto i = 0; i < s; i++)
		v.emplace_back(std::make_pair(i, (double)i));
	for (auto i = 0; i < s; i++)
		v[i];
	std::cout << Green << "Elapsed time: " << Red << watch.get()[0].count() << reset << " s\n";
	/* ------------------------------------------------------ */
	std::unordered_map<int, double> u_m;
	u_m.reserve(s);
	for (auto i = 0; i < s; i++)
		u_m[i] = (double)i;
	for (auto i = 0; i < s; i++)
		u_m[i];
	std::cout << Green << "Elapsed time: " << Red << watch.get()[0].count() << reset << " s\n";
	/* ------------------------------------------------------ */
	std::map<int, double> m;
	for (auto i = 0; i < s; i++)
		m[i] = (double)i;
	for (auto i = 0; i < s; i++)
		m[i];
	std::cout << Green << "Elapsed time: " << Red << watch.get()[0].count() << reset << " s\n";
};
