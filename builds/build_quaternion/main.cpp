#include "basic.hpp"
#include "Network.hpp"
#include "vtkWriter.hpp"
/**
 * クォータニオン
 *
 */

void translate(Network* const net, const Tddd& shift)
{
	for (auto& p : net->getPoints())
		p->setXSingle(p->initialX + shift);
	net->setGeometricProperties();
};

void rotate(Network* const net, const Quaternion& Q)
{
	Tddd c = { 0, 0, 0 };
	for (auto& p : net->getPoints())
		p->setXSingle(Q.Rv(p->initialX - c));
	net->setGeometricProperties();
};

int main()
{
	{
		auto net = new Network("bunny.obj", "bunny");
		auto mean = Mean(ToX(net->getPoints()));
		translate(net, -mean);
		net->resetInitialX();
		int l = 0;
		for (const auto& axis : std::vector<Tddd>{ {1, 0, 0}, {0, 1, 0}, {0, 0, 1} })
		{
			for (auto i = 0; i < 50; ++i)
			{
				auto Q = Quaternion(axis, 2 * M_PI / 50. * i);
				vtkPolygonWriter<networkPoint*> vtp;
				rotate(net, Q);
				for (const auto& f : net->getFaces())
				{
					auto [p0, p1, p2] = f->getPointsTuple();
					vtp.add({ p0, p1, p2 });
					vtp.addPolygon({ p0, p1, p2 });
				}
				std::cout << l << ", {yaw,pitch,roll} = " << Q.YPR() << std::endl;
				std::ofstream ofs("./output/bunny" + std::to_string(l++) + ".vtp");
				vtp.write(ofs);
			}
		}
	}

	{
		auto net = new Network("cow.obj", "cow");
		auto mean = Mean(ToX(net->getPoints()));
		translate(net, -mean);
		net->resetInitialX();
		int l = 0;
		Quaternion Q;
		for (auto i = 0; i < 50; ++i)
		{
			vtkPolygonWriter<networkPoint*> vtp;
			// Q *= Quaternion({ 1, 0, 0 }, 2 * M_PI / 50.);
			Q *= Quaternion({ 0, 0, 1 }, 2 * M_PI / 50.);
			rotate(net, Q);
			for (const auto& f : net->getFaces())
			{
				auto [p0, p1, p2] = f->getPointsTuple();
				vtp.add({ p0, p1, p2 });
				vtp.addPolygon({ p0, p1, p2 });
			}
			std::cout << l << ", {yaw,pitch,roll} = " << Q.YPR() << std::endl;
			std::ofstream ofs("./output/armadillo" + std::to_string(l++) + ".vtp");
			vtp.write(ofs);
		}
	}

	{
		auto net = new Network("camel.obj", "camel");
		auto mean = Mean(ToX(net->getPoints()));
		translate(net, -mean);
		net->resetInitialX();
		int l = 0;
		Quaternion Qyaw, Qpitch, Qroll;
		for (auto i = 0; i < 50; ++i)
		{
			vtkPolygonWriter<networkPoint*> vtp;
			Qyaw *= Quaternion({ 0, 0, 1 }, 2 * M_PI / 50.);
			Qpitch *= Quaternion({ 0, 1, 0 }, 2 * M_PI / 50.);
			Qroll *= Quaternion({ 1, 0, 0 }, 2 * M_PI / 50.);
			auto Q = Quaternion({ 0,0,1 }, Qyaw.yaw());
			Q *= Quaternion({ 0,1,0 }, Qpitch.pitch());
			Q *= Quaternion({ 1,0,0 }, Qroll.roll());
			rotate(net, Q);
			for (const auto& f : net->getFaces())
			{
				auto [p0, p1, p2] = f->getPointsTuple();
				vtp.add({ p0, p1, p2 });
				vtp.addPolygon({ p0, p1, p2 });
			}
			std::cout << l << ", {yaw,pitch,roll} = " << Q.YPR() << std::endl;
			std::ofstream ofs("./output/camel" + std::to_string(l++) + ".vtp");
			vtp.write(ofs);
		}
	}
};
