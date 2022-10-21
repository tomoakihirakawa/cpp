#include "InterpolationRBF.hpp"
#include "Network.hpp"

std::string home_dir = std::getenv("HOME");
std::string output_dir = home_dir + "/Network";

int main()
{
	/* ----------------------- 壁面粒子のもと ---------------------- */
	std::string input_name = "bunny";
	Network obj("../../obj/" + input_name + ".obj", "noname");
	std::string name = "bunny";
	obj.translate(-Mean(extractX(obj.getPoints())));
	//
	geometry::CoordinateBounds bounds(extractXtuple(obj.getPoints()));
	auto scale = bounds.getScale();
	std::cout << scale << std::endl;
	obj.scale({1 / scale, 1 / scale, 1 / scale});
	//
	mk_vtu(output_dir + "/" + name + ".vtu", obj.getFaces());
	//
	double meanlength = Mean(extLength(obj.getLines()));
	//
	/* ------------------------------------------------------ */
	int count = 0;
	for (const auto &depth : Subdivide(-0.01, 0.01, 10))
	{
		for (const auto &f : obj.getFaces())
			f->particlize(0.005, {depth});
		mk_vtu(output_dir + "/particles" + std::to_string(count++) + ".vtu", {obj.getPoints()});
		obj.DeleteParticles();
	}
	// mk_vtu(output_dir + "/particles" + std::to_string(count++) + ".vtu", {obj.getPoints()});
	// mk_vtu(output_dir + "/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.005, {0.05})});
	// mk_vtu(output_dir + "/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.005, {0.1})});
	// mk_vtu(output_dir + "/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.005, {0.15})});
	// mk_vtu(output_dir + "/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.005, {0.2})});
	// mk_vtu(output_dir + "/particles" + std::to_string(count++) + ".vtu", {particlize(obj.getFaces(), 0.01, {0.25})});
	// mk_vtu(output_dir + "/particles" + std::to_string(count++) + ".vtu", {particlize(*obj.getFaces().begin(), 0.01, {0., 0.1, 0.2, 0.3, 0.4, 0.5})});
}