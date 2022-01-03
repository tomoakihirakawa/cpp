#include "InterpolationRBF.hpp"
#include "Network.hpp"

std::string home_dir = std::getenv("HOME");

bool refine(netLp l, double len)
{
	if (l->length() > len)
	{
		l->divide();
		return true;
	}
	return false;
};
/* ------------------------------------------------------ */
auto creteOBJ = [](std::ofstream &ofs, Network &net)
{
	std::map<netPp, int> P_i;
	int i = 0;
	for (const auto &p : net.getPoints())
	{
		P_i[p] = ++i;
		auto X = p->getX();
		ofs << "v "
			<< X[0] << " "
			<< X[2] << " "
			<< X[1] << std::endl;
	}
	ofs << std::endl;

	for (const auto &f : net.getFaces())
	{
		auto ps = f->getPoints();
		ofs << "f "
			<< P_i[ps[0]] << " "
			<< P_i[ps[2]] << " "
			<< P_i[ps[1]] << std::endl;
	}
};
//* ------------------------------------------------------ */
//*                           メイン                        */
//* ------------------------------------------------------ */
#define bem
#if defined(bem)
int main()
{
	JSON json(std::ifstream("./input.json"));
	Network net(json["objfile"][0], json["name"][0]);
	auto name = json["name"][0];
	net.displayStates();

	double lim_len = .2;
	auto extractLength = [](const V_netLp &lines)
	{
		V_d ret;
		for (auto l : lines)
			ret.emplace_back(l->length());
		return ret;
	};
	//* ------------------------------------------------------ */
	//*                          線の分割                        */
	//* ------------------------------------------------------ */
	for (auto count = 0; count < 25; ++count)
	{
		//! ------------------------------------------------------ */
		mk_vtu(json["outputdir"][0] + name + std::to_string(count) + ".vtu", net.getFaces());
		std::ofstream ofs(json["outputdir"][0] + name + std::to_string(count) + ".obj");
		creteOBJ(ofs, net);
		ofs.close();

		//! ------------------------------------------------------ */
		Histogram h(extLength(net.getLines()));
		std::cout << Grid({"count", h.count}, 50) << std::endl;
		std::cout << Grid({"cumulative_count", h.cumulative_count}, 50) << std::endl;
		std::cout << Grid({"diff", h.diff}, 50) << std::endl;
		std::cout << Grid({"interval", h.interval}, 50) << std::endl;

		int num = 0;
		/* ------------------------*/
		//* | 0 | 1 | 2 |[3]| 4 |   diff, mid_interval
		//! 0   1   2  [3]  4   5   interval, cumulative_count
		/* ------------------------*/
		for (auto j = 0; j < h.diff.size(); ++j)
			if (h.cumulative_count[j + 1] / h.data_size >= 0.6 && h.diff[j] >= h.diff[num == 0 ? j : num])
				num = j;
		std::cout << "num " << num << std::endl;
		std::cout << h.interval[num] << std::endl;
		if (!num == 0)
		{
			bool found = false;
			do
			{
				found = false;
				auto lines = net.getLines();
				sortByLength(lines);
				for (const auto &l : Reverse(lines))
				{
					std::cout << l->length() << ", " << h.interval[num + 1] << std::endl;
					if (l->length() >= h.interval[num])
					{
						auto p = l->divide();
						found = true;
						for (auto &ps : BFS(p, 3))
							for (auto &q : ps)
							{
								for (auto &l : q->getLines())
								{
									double lim_degree = 1.; // degree
									auto fs = l->getFaces();
									// if (MyVectorAngle(fs[0]->getNormal(), fs[1]->getNormal()) / M_PI * 180. < lim_degree)
									// 	l->flipIfIllegal();
									/* ------------------------------------------------------ */
									found = l->flipIfBetter();
									// if (MyVectorAngle(fs[0]->getNormalTuple(), fs[1]->getNormalTuple()) / M_PI * 180. < lim_degree)
									// {
									// 	if (!l->islegal())
									// 	{
									// 		found = true;
									// 		// 平らな三角形通しの場合のみflipで対応すべき
									// 		// auto [Pa, Pb] = l->getPointsTuple();
									// 		// auto diff_size_points = std::abs((int)(Pa->getLines().size() - Pb->getLines().size()));
									// 		// auto fs = l->getFaces();
									// 		auto diff_area = std::abs(fs[0]->getArea() - fs[1]->getArea());
									// 		l->flip();
									// 		// まだ不正な線のままで，点の数が違いすぎる場合戻す．
									// 		if (!l->islegal())
									// 		{
									// 			// auto [Pa, Pb] = l->getPointsTuple();
									// 			// if (std::abs((int)(Pa->getLines().size() - Pb->getLines().size())) > diff_size_points)
									// 			// 	l->flip();
									// 			auto fs = l->getFaces();
									// 			if (std::abs(std::abs(fs[0]->getArea() - fs[1]->getArea())) > diff_area)
									// 				l->flip();
									// 		}
									// 	}
									// }
									if (found)
										break;
								}
								if (found)
									break;
							}
						if (found)
							break;
					}
				}
				for (auto &p : net.getPoints())
					LaplacianSmoothingIfFlat(p);
				for (auto &p : net.getPoints())
					LaplacianSmoothingIfFlat(p);
			} while (found);
		}
	}
}
#elif defined(sphare)
int main()
{
	// https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/
	const double X = .525731112119133606f;
	const double Z = .850650808352039932f;
	const double N = 0.;

	const VV_d vertices =
		{{-X, N, Z},
		 {X, N, Z},
		 {-X, N, -Z},
		 {X, N, -Z},
		 {N, Z, X},
		 {N, Z, -X},
		 {N, -Z, X},
		 {N, -Z, -X},
		 {Z, X, N},
		 {-Z, X, N},
		 {Z, -X, N},
		 {-Z, -X, N}};

	const VV_i connection =
		{{0, 4, 1},
		 {0, 9, 4},
		 {9, 5, 4},
		 {4, 5, 8},
		 {4, 8, 1},
		 {8, 10, 1},
		 {8, 3, 10},
		 {5, 3, 8},
		 {5, 2, 3},
		 {2, 7, 3},
		 {7, 10, 3},
		 {7, 6, 10},
		 {7, 11, 6},
		 {11, 0, 6},
		 {0, 1, 6},
		 {6, 1, 10},
		 {9, 0, 11},
		 {9, 11, 2},
		 {9, 2, 5},
		 {7, 2, 11}};

	Network net;
	net.generateNetwork(vertices, connection);
	net.setBounds();
	net.displayStates();

	double lim_len = .2;
	auto extractLength = [](const V_netLp &lines)
	{
		V_d ret;
		for (auto l : lines)
			ret.emplace_back(l->length());
		return ret;
	};
	//* ------------------------------------------------------ */
	//*                          線の分割                        */
	//* ------------------------------------------------------ */
	for (auto count = 0; count < 10; ++count)
	{
		//! ------------------------------------------------------ */
		mk_vtu("./vtu/sphere_divide" + std::to_string(count) + ".vtu", net.getFaces());
		std::ofstream ofs(home_dir + "/vtu/sphere_divide" + std::to_string(count) + ".obj");
		creteOBJ(ofs, net);
		ofs.close();

		//! ------------------------------------------------------ */
		Histogram h(extLength(net.getLines()));
		std::cout << Grid({"count", h.count}, 50) << std::endl;
		std::cout << Grid({"cumulative_count", h.cumulative_count}, 50) << std::endl;
		std::cout << Grid({"diff", h.diff}, 50) << std::endl;
		std::cout << Grid({"interval", h.interval}, 50) << std::endl;

		int num = 0;
		/* ------------------------*/
		//* | 0 | 1 | 2 |[3]| 4 |   diff, mid_interval
		//! 0   1   2  [3]  4   5   interval, cumulative_count
		/* ------------------------*/
		for (auto j = 0; j < h.diff.size(); ++j)
			if (h.cumulative_count[j + 1] / h.data_size >= 0.8 && h.diff[j] >= h.diff[num == 0 ? j : num])
				num = j;
		std::cout << "num " << num << std::endl;
		std::cout << h.interval[num] << std::endl;
		if (!num == 0)
		{

			bool found = false;
			do
			{
				found = false;
				for (const auto &l : net.getLines())
				{
					std::cout << l->length() << ", " << h.interval[num + 1] << std::endl;
					if (l->length() >= h.interval[num])
					{
						std::cout << 1 << std::endl;
						auto p = l->divide();
						std::cout << 2 << std::endl;
						p->setX(p->getX() / Norm(p->getX()));
						std::cout << 3 << std::endl;
						found = true;
						std::cout << 4 << std::endl;
						for (auto &p : Flatten(BFS(p, 3)))
							for (auto &L : p->getLines())
								L->flipIfIllegal();
						std::cout << 5 << std::endl;
					}
				}
			} while (found);
		}
	}
	//* ------------------------------------------------------ */
	//*                          面のマージ                      */
	//* ------------------------------------------------------ */
	int cc = 0;
	for (auto count = 0; count < 10; ++count)
	{
		//! ------------------------------------------------------ */
		Histogram h(extLength(net.getLines()));
		std::cout << Grid({"count", h.count}, 50) << std::endl;
		std::cout << Grid({"cumulative_count", h.cumulative_count}, 50) << std::endl;
		std::cout << Grid({"diff", h.diff}, 50) << std::endl;
		std::cout << Grid({"interval", h.interval}, 50) << std::endl;

		int num = 0;
		/* ------------------------*/
		//* | 0 | 1 | 2 |[3]| 4 |   diff, mid_interval
		//! 0   1   2  [3]  4   5   interval, cumulative_count
		/* ------------------------*/
		for (auto j = 0; j < h.diff.size(); ++j)
			if (h.cumulative_count[j] / h.data_size <= 0.8 && h.diff[j] >= std::abs(h.diff[num]))
				num = j;
		if ((h.diff[num] && !num == h.interval.size()) > 0 || num == 0)
			num++;
		std::cout << "num " << num << std::endl;
		std::cout << h.interval[num] << std::endl;
		try
		{
			if (num != h.interval.size())
			{
				bool found = false;
				do
				{
					found = false;
					std::cout << "size = " << net.getLines().size() << std::endl;
					net.displayStates();
					/* ------------------------------------------------------ */
					// デバッグのために，適当に削除していく　
#define good_merge
#if defined(bad_merge)
					for (const auto &l : net.getLines())
					{
						auto p = l->merge();
						p->setX(p->getX() / Norm(p->getX()));
						for (auto &p : Flatten(BFS(p, 3)))
							for (auto &L : p->getLines())
								L->flipIfIllegal(); //フリップがないと歪な三角形ができる．
						found = true;
						mk_vtu("./vtu/sphere_merge_0" + std::to_string(cc++) + ".vtu", net.getFaces());
						break;
					}
#elif defined(good_merge)
					/* ------------------------------------------------------ */
					// 比較的大きい方から
					for (const auto &l : net.getLines())
					{
						if (!l)
						{
							std::cout << "l = " << l << std::endl;
							std::cout << l->length() << ", " << h.interval[num + 1] << std::endl;
						}
						if (l->length() <= h.interval[num])
						{
							auto p = l->merge();
							p->setX(p->getX() / Norm(p->getX()));
							found = true;
							for (auto &p : Flatten(BFS(p, 3)))
								for (auto &L : p->getLines())
									L->flipIfIllegal();
							//! ------------------------------------------------------ */
							mk_vtu(home_dir + "/vtu/sphere_merge_0" + std::to_string(cc++) + ".vtu", net.getFaces());
							// std::ofstream ofs("./vtu/sphere_merge" + std::to_string(count) + ".obj");
							// creteOBJ(ofs, net);
							// ofs.close();
							break;
						}
						if (found)
							break;
					}
#endif
				} while (found);
			}
		}
		catch (std::exception e)
		{
			std::cout << e.what() << "\n";
		}
	}
}
#endif