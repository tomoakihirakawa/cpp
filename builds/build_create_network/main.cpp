#include "InterpolationRBF.hpp"
#include "Network.hpp"

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
			<< X[1] << " "
			<< X[2] << std::endl;
	}
	ofs << std::endl;

	for (const auto &f : net.getFaces())
	{
		auto ps = f->getPoints();
		ofs << "f "
			<< P_i[ps[0]] << " "
			<< P_i[ps[1]] << " "
			<< P_i[ps[2]] << std::endl;
	}
};

#define make_cube
#if defined(make_cube)
int main()
{
	/* ----------------------- 壁面粒子のもと ---------------------- */
	std::string input_name = "water";
	std::string name = "water";
	Network obj("../../obj/" + input_name + ".obj", name);
	// obj.translate(-Mean(extractX(obj.getPoints())));
	//
	geometry::CoordinateBounds bounds(extractXtuple(obj.getPoints()));
	auto scale = bounds.getScale();
	std::cout << scale << std::endl;
	// obj.scale({1 / scale, 1 / scale, 1 / scale});
	//
	mk_vtu("./vtu/" + name + ".vtu", obj.getFaces());
	//
	double meanlength = Mean(extLength(obj.getLines()));
	//
	for (auto i = 0; i < 15; ++i)
	{
		Histogram h(extLength(obj.getLines()));
		std::cout << Grid({"count", "cumulative_count", "diff", "interval"}, 50) << std::endl;
		std::cout << Grid({h.count, h.cumulative_count, h.diff, h.interval}, 50) << std::endl;

		int num = 0;
		for (auto j = 0; j < h.diff.size(); ++j)
			if (h.cumulative_count[j + 1] / h.data_size >= 0.5 && h.diff[j] >= h.diff[num == 0 ? j : num])
				num = j;

		bool cut = false;
		do
		{
			cut = false;
			for (const auto &l : obj.getLines())
				if (l->length() >= h.interval[num])
				{
					auto p = l->divide();
					cut = true;
					for (auto &p : Flatten(BFS(p, 3)))
						for (auto &L : p->getLines())
						{
							auto fs = L->getFaces();
							if (fs.size() > 1 && 10 > 180. / M_PI * MyVectorAngle(fs[0]->getNormal(), fs[1]->getNormal()))
								L->flipIfIllegal();
						}
				}
		} while (cut);

		mk_vtu("./vtu/" + name + std::to_string(i) + ".vtu", obj.getFaces());
		std::ofstream ofs("./vtu/" + name + std::to_string(i) + ".obj");
		creteOBJ(ofs, obj);
		ofs.close();

		/* ------------------------------------------------------ */
		int count = 0;
		//
		mk_vtu("./vtu/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.01, {0.})});
		mk_vtu("./vtu/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.01, {0.05})});
		mk_vtu("./vtu/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.01, {0.1})});
		mk_vtu("./vtu/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.01, {0.15})});
		mk_vtu("./vtu/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.01, {0.2})});
		mk_vtu("./vtu/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.01, {0.25})});
		mk_vtu("./vtu/particles" + std::to_string(count++) + ".vtu", {particlize(*obj.getFaces().begin(), 0.01, {0., 0.1, 0.2, 0.3, 0.4, 0.5})});
	}
}
#elif defined(check_interpolation)
int main()
{
	/* ----------------------- 壁面粒子のもと ---------------------- */
	NetworkObj obj("../../obj/tank2.obj");
	obj.scale({1., 1., 1.});
	mk_vtu("./vtu/obj_org.vtu", obj.getFaces());

	double dx = .01;
	auto tmp = new Network;
	for (const auto &X : Flatten(InterpolateFacesLinear(obj.getFaces(), dx, {0})))
		new networkPoint(tmp, tmp, {X[0], X[1], X[2]});
	mk_vtu("./vtu/obj.vtu", {tmp->getPoints()});
}
#else

int main()
{
	// https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/

	Network net;

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
	int count = 0;

	while (true)
	{
		/* ------------ 格子点を中心にずらすことで，重複をさけ粒子点を配置できる． ----------- */

		int l = 0;
		// for (const auto &particles : InterpolateFacesLinear(net.getFaces(), 0.02, {0, 1, 2, 3}))
		// {
		// 	mk_vtu("./vtu/particles" + std::to_string(count) + "_" + std::to_string(l++) + ".vtu", {particles});
		// }
		// }
		// mk_vtu("./vtu/particles" + std::to_string(count++) + ".vtu", {particlesXs});
		//
		mk_vtu("./vtu/sphere" + std::to_string(count++) + ".vtu", net.getFaces());
		std::ofstream ofs("./vtu/sphere" + std::to_string(count) + ".obj");
		creteOBJ(ofs, net);
		ofs.close();
		// auto lines = net.getLines();
		// /* ---------------------- 最大最小のチェック --------------------- */
		// double maxlength = -1E+10, minlength = 1E+10;
		// for (auto l : lines)
		// {
		// 	maxlength = (l->length() > maxlength) ? l->length() : maxlength;
		// 	minlength = (l->length() < minlength) ? l->length() : minlength;
		// }
		// std::cout << "maxlength = " << maxlength << std::endl;
		// std::cout << "minlength = " << minlength << std::endl;
		// /* --------------------- 上位10％を含む長さ --------------------- */
		// double select_len = 0.;
		// auto div = Reverse(Subdivide(minlength, maxlength, 100)); /*大きい場合*/
		// std::cout << "div = " << div << std::endl;
		// for (auto i = 0; i < div.size(); i++)
		// {
		// 	int c = 0;
		// 	for (auto l : lines)
		// 		if (l->length() > /*大きい場合*/ div[i])
		// 			c++;
		// 	if (c / lines.size() > 0.01)
		// 	{
		// 		select_len = div[i];
		// 		break;
		// 	}
		// }

		// std::cout << "select_len = " << select_len << std::endl;
		// double var = Variance(extractLength(lines));
		// std::cout << "variance = " << var << std::endl;
		// if (maxlength > lim_len && minlength > lim_len)
		// 	select_len = lim_len;
		// if (maxlength < lim_len && var < 0.001)
		// 	break;

		//! ------------------------------------------------------ */
		Histogram h(extLength(net.getLines()));
		std::cout << Grid({"count", h.count}, 50) << std::endl;
		std::cout << Grid({"cumulative_count", h.cumulative_count}, 50) << std::endl;
		std::cout << Grid({"diff", h.diff}, 50) << std::endl;
		std::cout << Grid({"interval", h.interval}, 50) << std::endl;

		int num = 0;
		for (auto j = 0; j < h.diff.size(); ++j)
			if (h.cumulative_count[j + 1] / h.data_size >= 0.5 && h.diff[j] >= h.diff[num == 0 ? j : num])
				num = j;
		std::cout << "num " << num << std::endl;
		std::cout << h.interval[num] << std::endl;
		//
		if (!num == 0)
		{
			bool cut = false;
			do
			{
				cut = false;

				auto lines = net.getLines();

				std::sort(lines.begin(), lines.end(), [](auto l0, auto l1)
						  { return l0->length() > l1->length(); });

				for (const auto &l : lines)
				{
					// std::cout << l->length() << ", " << h.interval[num + 1] << std::endl;
					if (l->length() >= h.interval[num])
					{
						auto p = l->divide();
						p->setX(p->getX() / Norm(p->getX()));
						cut = true;
						for (auto &p : Flatten(BFS(p, 3)))
							for (auto &L : p->getLines())
								L->flipIfIllegal();
					}
				}

			} while (cut);
		}
		//! ------------------------------------------------------ */

		// for (auto &l : lines)
		// {
		// 	if (l->length() > select_len)
		// 	{
		// 		auto p = l->divide();
		// 		p->setX(p->getX() / Norm(p->getX()));
		// 		// std::cout << "divide" << std::endl;
		// 	}
		// }

		for (auto &l : net.getLines())
		{
			l->flipIfIllegal();
		}
		// if (net.getPoints().size() > 100)
		// 	for (auto q : net.getPoints())
		// 	{
		// 		LaplacianSmoothing(q);
		// 		q->setX(q->getX() / Norm(q->getX()));
		// 	}
	}
	// mk_vtu("./vtu/sphere" + std::to_string(count++) + ".vtu", net.getFaces());
	// std::ofstream ofs("./vtu/sphere" + std::to_string(count) + ".obj");
	// creteOBJ(ofs, net);
	// ofs.close();
}
#endif