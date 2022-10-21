#include "InterpolationRBF.hpp"
#include "Network.hpp"

#if defined(__APPLE__)
std::string home_dir = "/Users/tomoaki";
#else
std::string home_dir = "/home/tomoaki";
#endif

//* ------------------------------------------------------ */
//*                           メイン                        */
//* ------------------------------------------------------ */
int main()
{
	JSON json(std::ifstream("./input.json"));
	Network net(json["objfile"][0], "example_file");
	net.rotate(T4d{-M_PI / 2., 1, 0, 0});
	/* ------------------------------------------------------ */
	double r = stod(json["r"][0]);
	VV_d polared_points;
	uomap_P_Tddd P_normal, P_grad_total_xyz;
	for (const auto &p : net.getPoints())
	{
		std::cout << "make intp" << std::endl;
		auto intp = parametricPolarInterpolation_(p);
		auto intp2 = parametricPolarInterpolation2_(p);
		std::cout << "intp done" << std::endl;
		for (const auto &t0 : Subdivide(-r, r, 20))
			for (const auto &t1 : Subdivide(-r, r, 20))
			{
				std::cout << intp({t0, t1}) << std::endl;
				polared_points.emplace_back(ToVector(intp({t0, t1})));
			}
		P_normal[p] = p->getNormalTuple();
		auto [dxdt0t1, dydt0t1, dzdt0t1] = intp2.grad({0, 0});
		auto [dxdt0, dxdt1] = dxdt0t1;
		auto [dydt0, dydt1] = dydt0t1;
		auto [dzdt0, dzdt1] = dzdt0t1;
		P_grad_total_xyz[p] = {Total(dxdt0t1), Total(dydt0t1), Total(dzdt0t1)};
	}
	// あとはグラッド
	VV_VarForOutput data = {{"normal", P_normal}, {"P_grad_total_xyz", P_grad_total_xyz}};
	mk_vtu(home_dir + "/vtu/polared_points.vtu", {polared_points});
	mk_vtu(home_dir + "/vtu/polared_points_base.vtu", {net.getPoints()}, data);
}

// int main()
// {
// 	JSON json(std::ifstream("./input.json"));
// 	Network net(json["objfile"][0], "example_file");
// 	/* ------------------------------------------------------ */
// 	std::vector<Tddd> polared_points;
// 	auto lines = net.getLines();
// 	for (auto &l : lines)
// 		l->CORNER = true;
// 	for (const auto &f : net.getFaces())
// 	{
// 		auto intp = polarInterpolation2(f);
// 		for (const auto &t0 : Subdivide(0., 1., 20))
// 		{
// 			for (const auto &t1 : Subdivide(0., 1., 20))
// 				polared_points.emplace_back(intp(t0, t1));
// 		}
// 	}
// 	mk_vtu(home_dir + "/vtu/polared_points_original.vtu", net.getFaces());
// 	mk_vtu(home_dir + "/vtu/polared_points.vtu", {polared_points});
// }
