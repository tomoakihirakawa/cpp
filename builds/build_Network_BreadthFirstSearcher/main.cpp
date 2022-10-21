#include "InterpolationRBF.hpp"
#include "Network.hpp"

#if defined(__APPLE__)
std::string home_dir = "/Users/tomoaki";
#else
std::string home_dir = "/home/tomoaki";
#endif

bool refine(netLp l, double len)
{
	if (l->length() > len)
	{
		l->divide();
		return true;
	}
	return false;
};

int main()
{
	/* ----------------------- 壁面粒子のもと ---------------------- */
	std::string input_name = "bunny";
	std::string name = "bunny";
	NetworkObj obj("../../obj/" + input_name + ".obj");
	obj.translate(-Mean(extractX(obj.getPoints())));
	//
	geometry::CoordinateBounds bounds(extractXtuple(obj.getPoints()));
	auto scale = bounds.getScale();
	std::cout << scale << std::endl;
	obj.scale({1 / scale, 1 / scale, 1 / scale});
	mk_vtu(home_dir + "/vtu/" + name + ".vtu", obj.getFaces());
	//
	int count = 0;
	for (const auto &q : BFSUO(*obj.getPoints().begin(), 10))
		mk_vtu(home_dir + "/vtu/" + name + "BFS_Points" + std::to_string(count++) + ".vtu", {q});
	count = 0;
	for (const auto &q : BFSUO(*obj.getFaces().begin(), 10))
		mk_vtu(home_dir + "/vtu/" + name + "BFS_Faces" + std::to_string(count++) + ".vtu", q);
	// mk_vtu(home_dir + "/vtu/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.01, {0.05})});
}
