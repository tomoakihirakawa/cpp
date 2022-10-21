#include <fundamental.hpp>

// std::unordered_map<Tii, Tdd> operator*(const double din, const std::unordered_map<Tii, Tdd> &v)
// {
// 	auto ret = v;
// 	for (auto &[ii, dd] : ret)
// 		dd *= din;
// 	return ret;
// };
// std::unordered_map<Tii, Tdd> operator*(const std::unordered_map<Tii, Tdd> &v, const double din)
// {
// 	auto ret = v;
// 	for (auto &[ii, dd] : ret)
// 		dd *= din;
// 	return ret;
// };
// std::unordered_map<Tii, Tdd> &operator+=(std::unordered_map<Tii, Tdd> &ii_dd, const std::unordered_map<Tii, Tdd> &jj_dd)
// {
// 	// for (auto &[ii, dd] : ii_dd)
// 	// 	dd += jj_dd[ii];
// 	return ii_dd;
// };

int main()
{

	std::unordered_map<Tii, Tdd> ii_dd, jj_dd;
	for (auto k = 0; k <= 10; ++k)
		for (auto m = -k; m <= k; ++m)
			ii_dd[Tii{k, m}] = {200, 200};

	for (auto k = 0; k <= 10; ++k)
		for (auto m = -k; m <= k; ++m)
			jj_dd[Tii{k, m}] = {(double)k, (double)m};

	jj_dd += ii_dd;
	jj_dd = jj_dd * 10.;

	for (const auto &[jj, dd] : jj_dd)
		std::cout << "jj = " << jj << ", dd = " << dd << std::endl;

	// int i = 2;
	// for (auto l = 0; l <= i; ++l)
	// 	for (auto m = -i; m <= i; ++m)
	// 		for (auto x = 0; x <= 10; ++x)
	// 		{
	// 			double q = (2. * M_PI * x / 10.);
	// 			std::cout << "std::sph_legendre ("
	// 					  << l << ","
	// 					  << m << ","
	// 					  << q << ") = "
	// 					  << std::sph_legendre(l, m, q) << std::endl;

	// 			std::cout << "std::assoc_legendre ("
	// 					  << l << ","
	// 					  << m << ","
	// 					  << q << ") = "
	// 					  << std::assoc_legendre(l, m, q) << std::endl;
	// 		}
	// return 0;
};
