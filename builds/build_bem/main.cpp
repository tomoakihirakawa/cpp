int time_step;

double real_time = 0;

#define simulation

#include "GNUPLOT.hpp"
#include "Network.hpp"
#include "integrationOfODE.hpp"
#include <filesystem>

pvd cpg_pvd("./vtu/bem.pvd");

//#define debug_bem
#include "bem.hpp"
#include "svd.hpp"

using V_i = std::vector<int>;
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
using V_Netp = std::vector<Network *>;
using V_netFp = std::vector<networkFace *>;
using VV_netFp = std::vector<V_netFp>;

networkLine *longerLine(networkLine *line_IN, double ratio = 1.01)
{
	networkLine *ret = line_IN;
	double len = line_IN->length(), v;

	for (const auto &p : line_IN->getPoints())
		for (const auto &l : p->getLines())
		{
			if (l != line_IN)
			{ /*omit comparison with line_IN line self*/
				v = l->length();
				if (v > ratio * len)
				{
					len = v;
					ret = l;
				}
			}
		}
	return ret;
};

std::string home_dir = std::getenv("HOME");
/* ------------------------------------------------------ */
double move_amplitude = 0.35;
Tddd translation(const double t)
{
	double s = M_PI / 2.;
	double k = M_PI / 1.;
	/* ------------------------------------------------------ */
	// Tddd move_dir = {cos(k * t), sin(k * t), 0.};
	// return move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) * move_dir;
	/* ------------------------------------------------------ */
	Tddd move_dir = Normalize(Tddd{1., 1., 0.});
	return move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) * move_dir;
};

T6d forced_velocity(const double t)
{
	double s = M_PI / 2.;
	double a = move_amplitude;
	double k = M_PI / 1.;
	/* ------------------------------------------------------ */
	// T6d move_dir = {cos(k * t), sin(k * t), 0., 0., 0., 0.};
	// T6d ddt_move_dir = {-k * sin(k * t), k * cos(k * t), 0., 0., 0., 0.};
	// // /* |U|*n_p . n_surface = phin <-- given
	// auto tmp = (-move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) + move_amplitude * exp(-t) * (cos(k * t - s) * k)) * move_dir;
	// tmp += move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) * ddt_move_dir;
	// return tmp;
	/* ------------------------------------------------------ */
	Tddd tmp = Normalize(Tddd{1., 1., 0.});
	T6d move_dir = {std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), 0., 0., 0.};
	return (-move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) + move_amplitude * exp(-t) * (cos(k * t - s) * k)) * move_dir;
};
/* ------------------------------------------------------ */
double phin(networkFace *f)
{
	auto normal = f->getNormalTuple();
	return Dot(f->getNetwork()->velocity,
			   T6d{
				   std::get<0>(normal),
				   std::get<1>(normal),
				   std::get<2>(normal),
				   0, 0, 0});
};

std::unordered_set<networkFace *> DeleteDuplicates_FacingContanctFaces(networkPoint *p)
{
	/* -------------------- contactfaces -------------------- */
	std::unordered_set<networkFace *> contactfaces;
	for (auto &f : p->getContactFaces())
	{
		bool duplicate = false;
		for (auto &cface : contactfaces)
		{
			auto angle = MyVectorAngle(cface->getNormalTuple(), f->getNormalTuple());
			if (2. * M_PI / 180. > angle)
			{
				duplicate = true;
				break;
			}
		}
		if (!duplicate)
			contactfaces.emplace(f);
	};
	/* ------------------------------------------------------ */
	return contactfaces;
};
/* ------------------------------------------------------ */
// ある点の接する面のphinを抜き出す．
std::vector<std::tuple<Tddd, double>> phins_contact(networkPoint *p)
{
	std::vector<std::tuple<Tddd, double>> ret;
	for (auto &f : DeleteDuplicates_FacingContanctFaces(p))
	{
		auto n = f->getNormalTuple();
		ret.emplace_back(std::tuple<Tddd, double>{n, Dot(f->getNetwork()->velocity, ToT6d(n))});
	}
	return ret;
}

double phin_contact(networkPoint *const p)
{
	auto contactfaces = DeleteDuplicates_FacingContanctFaces(p);
	double ret = 0.;
	int count = 0;
	for (auto &f : contactfaces)
	{
		auto n_face = f->getNormalTuple();
		auto n_point = p->getNormalAreaAveraged();
		auto U_normal_of_face = Dot(f->getNetwork()->velocity, ToT6d(n_face)) * n_face;
		auto divide = Dot(n_face, n_point);
		ret += Dot(U_normal_of_face, n_point);
		count++;
	}
	return ret / (double)contactfaces.size();
}

/* ------------------------------------------------------ */

std::unordered_set<networkFace *> facingFace(const networkFace *const f)
{
	//! getContactFacesは，点の周辺の点の影となった面を含まない．
	//! しかし，面の接触面を判断するには，そのような除外だけでは十分ではなく，
	//! さらに，この面と向き合っているかどうかの判定し除外する必要がある．それが以下．
	std::unordered_set<networkFace *> ret;
	for (const auto &q : f->getPoints())
		for (const auto &F : q->getContactFaces())
		{
			auto angle = MyVectorAngle(f->getNormalTuple(), F->getNormalTuple());
			if ((M_PI - angle) / M_PI * 180. < 30.)
				ret.emplace(F);
		}
	return ret;
};
/* ------------------------------------------------------ */
netFp closest_facingFace(const networkFace *const f_IN)
{
	networkFace *closest_face = nullptr;
	double min_distance = 1E+50;
	for (const auto &F : facingFace(f_IN))
	{
		auto dist = Norm(vectorToTriangle(F, f_IN->getXtuple()));
		if (dist < min_distance)
		{
			closest_face = F;
			min_distance = dist;
		}
	}
	return closest_face;
};
/* ------------------------------------------------------ */
double phin_from_Neumann_surface(networkPoint *const p)
{
	std::unordered_set<networkFace *> suraceNeumann;
	/*
	@ 周辺面のうちNeumann面を抜き出し，それぞれに対してclosest_facingFaceを取得する．
	@ 次に，面それぞれの速度を計算し，点の法線方向速度を計算する．
	*/
	for (const auto &F : p->getFaces())
		if (F->Neumann)
		{
			auto cf = closest_facingFace(F);
			if (cf)
				suraceNeumann.emplace(cf);
			else
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Fはノイマンなのになぜfacing faceがないのか?");
		}
	if (suraceNeumann.empty())
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

	double ret = 0.;
	for (auto &f : suraceNeumann)
	{
		auto n_face = f->getNormalTuple();
		auto n_point = p->getNormalAreaAveraged();
		auto U_normal_of_face = Dot(f->getNetwork()->velocity, ToT6d(n_face)) * n_face;
		// auto divide = Dot(n_face, n_point);
		ret += Dot(U_normal_of_face, n_point);
	}
	return ret / (double)suraceNeumann.size();
}
/* ------------------------------------------------------ */
double phin_contact(const netFp f_IN)
{
	auto closest_face = closest_facingFace(f_IN);
	if (closest_face)
	{
		// auto n_point = p->getNormalTuple();
		// auto U_normal_of_face = Dot(f->getNetwork()->velocity, ToT6d(n_face)) * n_face;
		return Dot(closest_face->getNetwork()->velocity,
				   ToT6d(f_IN->getXtuple()));
	}
	else
	{
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		return 0;
	}
}
/* ------------------------------------------------------ */
V_netLp getLinesAround(netPp p)
{
	V_netLp ret = {};
	for (const auto &f : p->getFaces())
		for (const auto &line : f->getLines())
			ret.emplace_back(line);
	return DeleteDuplicates(ret);
};
V_netLp getLinesAround(netLp l)
{
	V_netLp ret = {l};
	for (const auto &p : l->getPoints())
		for (const auto &f : p->getFaces())
			for (const auto &line : f->getLines())
				ret.emplace_back(line);
	return DeleteDuplicates(ret);
};
/* ------------------------------------------------------ */
// using V_RKRK = std::vector<std::shared_ptr<derivativeImprover>>;
using V_RKRK = std::vector<derivativeImprover *>;
using map_P_d = std::map<netP *, double>;
using map_P_Vd = std::map<netP *, V_d>;
using map_P_VVd = std::map<netP *, VV_d>;
using map_F_P_Vd = std::map<netF *, map_P_Vd>;
using map_P_P_Vd = std::map<netP *, map_P_Vd>;
using map_pairPF_P_Vd = std::map<std::pair<netP *, netF *> /*タプル*/, std::map<netP *, V_d>>;
using map_pairPF_P_Tdd = std::map<std::pair<netP *, netF *> /*タプル*/, std::map<netP *, Tdd>>;
using map_pairPF_pairPF_Tdd = std::map<std::pair<netP *, netF *> /*タプ*/, std::map<std::pair<netP *, netF *> /*タプル*/, Tdd>>;
using map_P_F_P_Vd = std::map<netP *, map_F_P_Vd>;

using VV_SorIorMap = std::vector<std::vector<std::variant<std::string, int, map_P_Vd>>>;

V_netFp takeFaces(const V_Netp &nets)
{
	V_netFp ret({});
	for (const auto &n : nets)
		ret.insert(ret.end(), n->getFaces().begin(), n->getFaces().end());
	return DeleteDuplicates(ret);
};
/* ------------------------------------------------------ */
V_d volume({});

auto modify = [](Network &water, const JSON &json1)
{
	auto isExsist = [&json1](std::string str)
	{ return (json1().find(str) != json1().end() && !json1()[str].empty()); };
	if (isExsist("rotate"))
	{
		auto rotate = stod(json1()["rotate"]);
		if (rotate.size() > 1)
			water.rotate(rotate[0], Tddd{rotate[1], rotate[2], rotate[3]});
	}
	if (isExsist("scale"))
	{
		auto scale = stod(json1()["scale"]);
		if (scale.size() > 1)
			water.scale({scale[0], scale[1], scale[2]});
	}
	if (isExsist("translate"))
	{
		auto translate = stod(json1()["translate"]);
		if (translate.size() > 1)
			water.translate({translate[0], translate[1], translate[2]});
	}
	// if (isExsist("remesh"))
	// {
	// 	auto minlen = stod(json1()["remesh"]);
	// 	if (minlen.size() > 0)
	// 		remesh(&water, minlen[0]);
	// }
	// if (isExsist("coarsen"))
	// {
	// 	auto minlen = stod(json1()["coarsen"]);
	// 	if (minlen.size() > 0)
	// 		coarsen(&water, minlen[0]);
	// }
	if (isExsist("reverseNormal"))
	{
		std::string TorF = json1()["reverseNormal"][0];
		if (TorF.compare("True") == 0 || TorF.compare("true") == 0 || TorF.compare("1") == 0)
		{
			water.reverseNormal();
			std::cout << "reverse done" << std::endl;
		}
	}
};

VV_d extVelocities(const V_netPp &ps)
{
	VV_d ret(0);
	for (const auto &p : ps)
		ret.emplace_back(ToVector(p->U_BEM));
	return ret;
};

// double dt = 0.01;
/* ------------------------------------------------------ */
struct values_for_overdetermined_interpolation3
{
	values_for_overdetermined_interpolation3(const std::unordered_set<networkPoint *> &Points,
											 networkPoint *centerp)
	{
		set(ToVector(Points), centerp);
	}
	values_for_overdetermined_interpolation3(const V_netPp &Points,
											 networkPoint *centerp)
	{
		set(Points, centerp);
	}
	/* ------------------------------------------------------ */
	std::vector<std::tuple<Tddd, double>> kernel_X_s;
	std::vector<std::tuple<Tddd, double>> eq;
	//! ------------------------------------------------------ */
	void set(const V_netPp &ps, networkPoint *const centerp)
	{
		try
		{
			this->kernel_X_s.clear();
			this->eq.clear();
			auto eps = 1E-2;
			double alpha = 1.;
			double scale;
			// double fixed_scale = alpha * std::sqrt(Total(extractAreas(centerp->getFaces()) / M_PI));//!
			double fixed_scale = alpha * Mean(extLength(centerp->getLines()));
			/* ------------------------------------------------------ */
			for (const auto &q : ps)
			{
				// scale = alpha * Mean(extLength(q->getLines()));
				// scale = alpha * std::sqrt(Total(extractAreas(q->getFaces()) / M_PI));
				scale = fixed_scale; //!
				/* ------------------------------------------------------ */
				auto [phi, phi_n] = q->phiphin;
				auto X = q->getXtuple();
				this->kernel_X_s.emplace_back(std::tuple<Tddd, double>{X, scale});
				this->eq.emplace_back(std::tuple<Tddd, double>{X, phi});
				/* ------------------------------------------------------ */
				//プラスマイナスのミスが多かった
				auto N = q->getNormalAreaAveraged();
				// if (q->multiple_phiphin.empty())
				// {
				// 	for (const auto& n : std::vector<Tddd>{ N * eps * scale, -N * eps * scale })
				// 	{
				// 		this->kernel_X_s.emplace_back(std::tuple<Tddd, double>{X + n, scale});
				// 		this->eq.emplace_back(std::tuple<Tddd, double>{ X + n, phi + phi_n * Dot(N, n)});
				// 	}
				// }
				// else
				// {
				// 	for (const auto& F_phiphin : q->multiple_phiphin) {
				// 		auto F = std::get<0>(F_phiphin);
				// 		auto [phi, phin] = std::get<1>(F_phiphin);
				// 		N = F->getNormalTuple();
				// 		for (const auto& n : std::vector<Tddd>{ N * eps * scale, -N * eps * scale })
				// 		{
				// 			this->kernel_X_s.emplace_back(std::tuple<Tddd, double>{X + n, scale});
				// 			this->eq.emplace_back(std::tuple<Tddd, double>{ X + n, phi + phi_n * Dot(N, n)});
				// 		}
				// 	}
				// }
			}
			// for (const auto& q : ps)
			// 	if (q->CORNER)
			// 		for (const auto& dir_phin/*面の法線ベクトルと，その方向の流速（面のvelocityから計算）*/ : phins_contact(q))
			// 		{
			// 			//scale = alpha * Mean(extLength(q->getLines()));
			// 			// scale = std::sqrt(Total(extractAreas(q->getFaces()) / M_PI));
			// 			scale = fixed_scale;//!
			// 			auto [phi, tmp] = q->phiphin;
			// 			auto X = q->getXtuple();
			// 			auto N = std::get<0>(dir_phin);
			// 			auto phi_n = std::get<1>(dir_phin);
			// 			for (const auto& n : std::vector<Tddd>{ -N * eps * scale })
			// 			{
			// 				this->kernel_X_s.emplace_back(std::tuple<Tddd, double>{X + n, scale});
			// 				this->eq.emplace_back(std::tuple<Tddd, double>{ X + n, phi + phi_n * Dot(N, n)});
			// 			}
			// 		}
			//! ---------------------- 面を加える -------------------------- */
			// std::unordered_set<networkFace*> faces;
			// for (const auto& q : ps)
			// 	for (const auto& f : q->getFaces())
			// 		faces.emplace(f);

			// for (const auto& f : faces)
			// if (!centerp->CORNER)
			// 	for (const auto& f : centerp->getFaces())
			// 	{
			// 		double t0 = 1 / 4., t1 = 1 / 4.;
			// 		double t2 = 1 - t0 - t1;
			// 		V_d xi = { t0 * (2 * t0 - 1),
			// 					t1 * (2 * t1 - 1),
			// 					t2 * (2 * t2 - 1),
			// 					4 * t0 * t1,
			// 					4 * t1 * t2,
			// 					4 * t0 * t2 };
			// 		auto [p0, p1, p2, p3, p4, p5] = f->get6PointsTuple(centerp);
			// 		auto intpX = Dot(xi, VV_d{ p0->getX(),
			// 									p1->getX(),
			// 									p2->getX(),
			// 									p3->getX(),
			// 									p4->getX(),
			// 									p5->getX() });
			// 		double phi = Dot(xi, V_d{ std::get<0>(p0->phiphin),
			// 							std::get<0>(p1->phiphin),
			// 							std::get<0>(p2->phiphin),
			// 							std::get<0>(p3->phiphin),
			// 							std::get<0>(p4->phiphin),
			// 							std::get<0>(p5->phiphin) });
			// 		scale = fixed_scale;//!
			// 		/* ------------------------------------------------------ */
			// 		// Tddd X = Tddd{ intpX[0],intpX[1],intpX[2] };
			// 		auto X = f->getXtuple();
			// 		this->kernel_X_s.emplace_back(std::tuple<Tddd, double>{X, scale});
			// 		this->eq.emplace_back(std::tuple<Tddd, double>{ X, phi});
			// 	}
			/* ------------------------------------------------------ */
			// for (const auto& f : centerp->getFaces())
			// {
			// 	// double t0 = 1 / 4., t1 = 1 / 4.;
			// 	// double t2 = 1 - t0 - t1;
			// 	// V_d xi = { t0 * (2 * t0 - 1),
			// 	// 			t1 * (2 * t1 - 1),
			// 	// 			t2 * (2 * t2 - 1),
			// 	// 			4 * t0 * t1,
			// 	// 			4 * t1 * t2,
			// 	// 			4 * t0 * t2 };
			// 	// auto [p0, p1, p2, p3, p4, p5] = f->get6PointsTuple(centerp);
			// 	// auto intpX = interpolationTriangleQuad({ p0->getX(),p1->getX(),p2->getX(),p3->getX(),p4->getX(),p5->getX() });
			// 	// auto intpPhi = Dot(V_d{ std::ge<0>(p0->phiphin),
			// 	// 					std::ge<1>(p0->phiphin),
			// 	// 					std::ge<2>(p0->phiphin),
			// 	// 					std::ge<3>(p0->phiphin),
			// 	// 					std::ge<4>(p0->phiphin),
			// 	// 					std::ge<5>(p0->phiphin) }, xi);
			// 	// auto intp1 = interpolationTriangleQuad(extractX(fs[1]->get6Points(l)));
			// 	// auto estimate_line_mid = ((intp0(1 / 2., 1 / 4.) + intp1(1 / 2., 1 / 4.)) / 2.);
			// 	/* ------------------------------------------------------ */
			// 	// scale = alpha * Mean(extLength(q->getLines()));
			// 	// scale = alpha * std::sqrt(f->getArea() / M_PI);
			// 	scale = fixed_scale;//!
			// 	/* ------------------------------------------------------ */
			// 	double phi = 0, phi_n = 0;
			// 	for (const auto& q : f->getPoints())
			// 	{
			// 		phi += std::get<0>(q->phiphin) / 3.;
			// 		phi_n += std::get<1>(q->phiphin) / 3.;//正確ではない．
			// 	}
			// 	auto X = f->getXtuple();
			// 	this->kernel_X_s.emplace_back(std::tuple<Tddd, double>{X, scale});
			// 	this->eq.emplace_back(std::tuple<Tddd, double>{ X, phi});
			// 	/* ------------------------------------------------------ */
			// 	//プラスマイナスのミスが多かった
			// 	//phi_nは正確ではない．
			// 	// auto N = f->getNormalTuple();
			// 	// for (const auto& n : std::vector<Tddd>{ N * eps * scale, -N * eps * scale })
			// 	// {
			// 	// 	this->kernel_X_s.emplace_back(std::tuple<Tddd, double>{X + n, scale});
			// 	// 	this->eq.emplace_back(std::tuple<Tddd, double>{ X + n, phi + phi_n * Dot(N, n)});
			// 	// }
			// }
			//! ------------------------------------------------------ */
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	}
};

/* ------------------------------------------------------ */
struct values_for_overdetermined_interpolation2
{
	values_for_overdetermined_interpolation2(const std::unordered_set<networkPoint *> &Points,
											 networkPoint *centerp)
	{
		set(ToVector(Points), centerp);
	}
	values_for_overdetermined_interpolation2(const V_netPp &Points,
											 networkPoint *centerp)
	{
		set(Points, centerp);
	}
	/* ------------------------------------------------------ */
	std::vector<Tddd> kernel_locations;
	std::vector<std::tuple<Tddd, double, double>> eq;
	std::vector<std::tuple<Tddd, Tddd, double, double>> eq_der_dot;
	//
	void set(const V_netPp &ps, networkPoint *centerp)
	{
		this->kernel_locations.clear();
		this->eq.clear();
		this->eq_der_dot.clear();
		for (const auto &q : ps)
		{
			auto m = 2. * Mean(extLength(q->getLines()));
			/* ------------------------------------------------------ */
			auto [phi, phin] = q->phiphin;
			auto X = q->getXtuple();
			this->kernel_locations.emplace_back(X);
			this->eq.emplace_back(std::tuple<Tddd, double, double>{X, phi, m});
			/* ------------------------------------------------------ */
			// if (!q->CORNER) {
			auto n = q->getNormalTuple();
			this->kernel_locations.emplace_back(X + n);
			auto eps = 0.1;
			this->eq.emplace_back(std::tuple<Tddd, double, double>{X + n, phi + std::get<1>(q->phiphin) * Dot(n, n) * eps, m});
			/* ------------------------------------------------------ */
			this->kernel_locations.emplace_back(X - n);
			this->eq.emplace_back(std::tuple<Tddd, double, double>{X - n, phi + std::get<1>(q->phiphin) * Dot(n, -n) * eps, m});
			// }
		}

		/* ------------------------------------------------------ */
		// for (const auto& q : Join(centerp->getNeighbors(), { centerp }))
		// {
		// 	//もし角点が大丈夫
		// 	{
		// 		// Xは同じでないといけないようだ
		// 		auto X = q->getXtuple();
		// 		auto N = q->getNormalTuple();
		// 		// auto Xmod = X + N * Mean(extLength(q->getLines()));
		// 		this->kernel_locations.emplace_back(Xmod);
		// 		this->eq_der_dot.emplace_back(std::tuple<Tddd, Tddd, double>{X, N, std::get<1>(q->phiphin) });
		// 	}
		// 	// {
		// 	// 	auto [phi, phin] = q->phiphin;
		// 	// 	auto X = q->getXtuple();
		// 	// 	auto N = q->getNormalTuple();
		// 	// 	X = X - N * Mean(extLength(q->getLines()));
		// 	// 	this->kernel_locations.emplace_back(X);
		// 	// 	this->eq.emplace_back(std::tuple<Tddd, double>{ X, phi });
		// 	// }
		// }
		/* ------------------------------------------------------ */
		// std::unordered_set<networkFace*> fs;
		// for (const auto& q : ps)
		// 	for (const auto& f : q->getFaces())
		// 		fs.emplace(f);
		// for (const auto& f : fs)
		// {
		// 	double mean_phin = 0.;
		// 	for (const auto& q : f->getPoints())
		// 	{
		// 		auto [phi, phin] = q->phiphin;
		// 		mean_phin += phin / 3.;
		// 	}
		// 	auto X = f->getXtuple();
		// 	this->kernel_locations.emplace_back(X);
		// 	this->eq_der_dot.emplace_back(std::tuple<Tddd, Tddd, double>{ X, f->getNormalTuple(), mean_phin });
		// }
	}
};
//
struct values_for_overdetermined_interpolation
{
	VV_d X;
	V_d Phi;
	V_d Phi_n_of_X;
	//
	VV_d Y;
	V_d Phi_n;
	VV_d normals;
	//
	values_for_overdetermined_interpolation(const std::unordered_set<networkPoint *> &Points)
	{
		set(ToVector(Points));
	}
	values_for_overdetermined_interpolation(const V_netPp &Points)
	{
		set(Points);
	}
	void set(const V_netPp &ps)
	{
		// this->X = obj3D::extractX(ps);
		// this->Phi = extPhi(ps);
		// this->Phi_n_of_X = extPhin(ps);
		this->X.clear();
		this->Phi.clear();
		this->Phi_n_of_X.clear();
		//
		std::unordered_set<networkFace *> contactfaces;
		for (const auto &q : ps)
		{
			this->X.emplace_back(q->getX());
			this->Phi.emplace_back(std::get<0>(q->phiphin));
			this->Phi_n_of_X.emplace_back(std::get<1>(q->phiphin));
			for (auto &f : q->getContactFaces())
			{
				bool duplicate = false;
				for (auto &cface : contactfaces)
					if (M_PI / 180. > MyVectorAngle(cface->getNormalTuple(), f->getNormalTuple()))
					{
						duplicate = true;
						break;
					}
				if (!duplicate)
					contactfaces.emplace(f);
			}
		}
		/* ------------------------ 鏡像の追加 ----------------------- */
		// for (const auto& q : ps)
		// 	for (const auto& f : contactfaces)
		// 	{
		// 		auto d = 1E-5;
		// 		auto X = f->mirrorPosition(q, 2.) + Tddd{ RandomReal({ -d,d }), RandomReal({ -d,d }), RandomReal({ -d,d }) };
		// 		Xs.emplace_back(ToVector(X));
		// 		Phis.emplace_back(std::get<0>(q->phiphin) + Norm((X - q->getXtuple())) * std::get<1>(p->phiphin));
		// 	}
		// for (const auto& q : Flatten(BFS(p, 2, { p->getNetwork() })))
		// 	for (const auto& f : q->getContactFaces())
		// 	{
		// 		Tddd move_dir = { 1.,0.,0. };
		// 		auto d = 1E-14;
		// 		auto r = Tddd{ RandomReal({ -d,d }), RandomReal({ -d,d }), RandomReal({ -d,d }) };
		// 		// auto X = f->mirrorPosition(q, 2.) + Tddd{ RandomReal({ -d,d }), RandomReal({ -d,d }), RandomReal({ -d,d }) };
		// 		auto x = vectorToTriangle(geometry::Triangle(f->getLocationsTuple()), q->getXtuple()) + q->getXtuple();
		// 		Xs.emplace_back(ToVector(x + r));
		// 		Phis.emplace_back(std::get<0>(q->phiphin) - Dot((x - q->getXtuple()), Normalize(move_dir)) * std::get<1>(p->phiphin));
		// 	}
		/* ------------------------------------------------------ */
		//! 接触面の境界条件を満たすようにする
		// for (const auto& q : ps) {
		// 	for (const auto& f : q->getContactFaces()) {
		// 		auto R = 1E-14 * RandomReal({ -1.,1. });
		// 		auto r = Tddd{ R,R,R };
		// 		// auto r = Tddd{ 1E-14,1E-14,1E-14 };
		// 		// auto r = Tddd{ RandomReal({ -d,d }), RandomReal({ -d,d }), RandomReal({ -d,d }) };
		// 		// auto x = vectorToTriangle(geometry::Triangle(f->getLocationsTuple()), q->getXtuple()) + q->getXtuple();
		// 		auto x = q->getXtuple();
		// 		Y.emplace_back(ToVector(x + r));
		// 		auto n = f->getNormalTuple();
		// 		Phi_n.emplace_back(phin(f));//!?
		// 		normals.emplace_back(ToVector(n));
		// 	}
		// 	auto R = 1E-14 * RandomReal({ -1.,1. });
		// 	auto r = Tddd{ R,R,R };
		// 	auto x = q->getXtuple();
		// 	Y.emplace_back(ToVector(x + r));
		// 	Phi_n.emplace_back(std::get<1>(q->phiphin));
		// 	normals.emplace_back(q->getNormal());
		// }
	}
};
#define ENHANCE_CORNER_ACCURACY
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
Tddd grad_phi_tangential(const networkFace *const f)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	auto A = f->getArea();
	auto n = f->getNormalTuple();
	auto X0 = p0->getXtuple();
	auto X1 = p1->getXtuple();
	auto X2 = p2->getXtuple();
	auto p0_phi = std::get<0>(p0->phiphin);
	auto p1_phi = std::get<0>(p1->phiphin);
	auto p2_phi = std::get<0>(p2->phiphin);
	return (
			   (p2_phi - p1_phi) * Cross(n, X1 - X0) +
			   (p0_phi - p1_phi) * Cross(n, X2 - X1)) /
		   (2 * A);
};
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
struct derivatives
{
	// #define derivatives_debug
	// public:
	uomap_P_Tddd P_tension, P_gradPhi, P_gradPhi_tangential, P_phin_vector, P_dxdt, P_dxdt_mod, P_laplacian;
	uomap_P_d P_DphiDt, P_kappa, P_pressure, P_aphiat;
	~derivatives(){
		// std::cout << "derivatives 破棄" << std::endl;
	};
	derivatives(std::unordered_set<networkPoint *> Points, const double dt)
	{
		// 	this->set(Points);
		// };
		// void set(const std::unordered_set<networkPoint*>& Points)
		// {
#ifdef derivatives_debug
		std::cout << Red << "initialize for parallelization" << reset << std::endl;
#endif
		//! initialize for parallelization
		for (const auto &p : Points)
		{
			this->P_gradPhi[p] = {0, 0, 0};
			this->P_gradPhi_tangential[p] = {0, 0, 0};
			this->P_phin_vector[p] = {0, 0, 0};
			this->P_dxdt[p] = {0, 0, 0};	 //流速
			this->P_dxdt_mod[p] = {0, 0, 0}; //流速
			this->P_DphiDt[p] = 0.;
			this->P_aphiat[p] = 0.;
			this->P_pressure[p] = 0.;
			this->P_kappa[p] = 0.;
			this->P_laplacian[p] = {0, 0, 0};
			this->P_tension[p] = {0, 0, 0};
		}
		auto pointsbegin = Points.begin();
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (auto it = Points.begin(); it != Points.end(); ++it)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			// std::cout << "曲率の計算" << std::endl;
			auto p = *it;

			if (!isFinite(p->phiphin))
			{
				std::cout << "p->phiphinはfiniteではない！！" << std::endl;
				std::cout << "p->phiphin = " << p->phiphin << std::endl;
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			}

			// auto fs = Flatten(BFS(p->getFaces(), 2, { p->getNetwork() }));
			// auto ps = Flatten(BFS(p, 2, { p->getNetwork() }));
			// auto interpNormals = InterpolationVectorRBF(
			// 	Join(extractX(fs), extractX(ps)),
			// 	Join(extNormals(fs), extNormals(ps)),
			// 	p->getX());
			// auto n = Normalize(interpNormals(p->getX()));

			V_netPp ps = Flatten(BFS(p, 4));
			auto interpNormals = InterpolationVectorRBF(ToVector(extX(ps)), ToVector(extNormals(ps)), p->getX());

			//!衝突
			// auto n = p->normalContanctSurface(1., 2.);
			auto n = p->getNormalTuple();
			p->normal_BEM = n;
			// 法線方向も怪しい．．．
			// p->normal_BEM = p->getNormal();
			// auto n = p->getNormal();
			p->kappa_BEM = interpNormals.div(p->getX()) / 2.; //中心方向法線ベクトルの場合，マイナスをつける．
			if (!isFinite(p->kappa_BEM))
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
				//! https://en.wikipedia.org/wiki/Mean_curvature
#ifdef derivatives_debug
			std::cout << "流速の計算" << std::endl;
#endif
			{
				// #define SIMPLE_RBF
// #define ENHANCED_RBF
// #define AngleWeighted
#define Dombre2019
				// #define Dombre2019_2
				/* ------------------------------------------------------ */
				Tddd grad = {0, 0, 0};
				/* ------------------------------------------------------ */

#if defined(SIMPLE_RBF)
				auto values = values_for_overdetermined_interpolation(ps);
				if (it == pointsbegin)
					std::cout << Red << "シンプルなRBFを使って，流速を計算する" << reset << std::endl;
				auto interp = InterpolationRBF(values.X, values.Phi, 2. * Mean(extLength(p->getLines() /*分布間隔よりも大きく，全体よりも小さく*/)));
				auto tmp = interp.grad(p->getX());
				grad = {tmp[0], tmp[1], tmp[2]};
				// phinは正しくないので補正
				grad -= Dot(grad, p->getNormalAreaAveraged()) * p->getNormalAreaAveraged();
				grad += std::get<1>(p->phiphin) * p->getNormalAreaAveraged();
#elif defined(ENHANCED_RBF)

				// bool all_Dirichlet = std::all_of(ps.begin(), ps.end(), [](const auto& p) {return p->Dirichlet;});
				// bool any_CORNER = std::any_of(ps.begin(), ps.end(), [](const auto& p) {return p->CORNER;});
				auto values = values_for_overdetermined_interpolation3(ps, p);
				if (it == pointsbegin)
					std::cout << Red << "工夫したなRBFを使って，流速を計算する" << reset << std::endl;
				// std::cout << Red << "過剰決定のRBFを使って，流速を計算する" << reset << std::endl;
				// auto interp = InterpolationRBF(
				// 	values.kernel_locations,
				// 	values.eq
				// );
				auto interp = InterpolationRBF_(
					values.kernel_X_s,
					values.eq);
				// phinは一応考慮されている．
				grad = interp.grad(p->getXtuple());
				grad -= Dot(grad, p->getNormalAreaAveraged()) * p->getNormalAreaAveraged();
				grad += std::get<1>(p->phiphin) * p->getNormalAreaAveraged();
				/* ------------------------------------------------------ */
				if (!isFinite(grad))
				{
					// std::cout << "interp.V = " << interp.V << std::endl;
					// std::cout << "interp.w = " << interp.w << std::endl;
					// std::cout << "interp.F = " << interp.F << std::endl;
					// std::cout << "Phi = " << values.Phi << std::endl;
					// std::cout << "X = " << values.X << std::endl;
					// std::cout << "interp =" << interp(p->getX()) << std::endl;
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not finite!");
				}
				// grad -= Dot(grad, p->getNormalAreaAveraged()) * p->getNormalAreaAveraged();
				// grad += std::get<1>(p->phiphin) * p->getNormalAreaAveraged();

#elif defined(AngleWeighted) || defined(AngleWeighted2) || defined(Dombre2019_2) || defined(Dombre2019)
				// Meyer et al. (2002).
				// 	線形の流速近似と比較
				// 	ext namespaceを作る．
				// 	Meyerの式を実装
				double TotalArea = 0., TotalArea_for_phi = 0., TotalArea_for_phin = 0.;
				double Theta = 0., TotalTheta = 0.;
				Tddd grad_tangential = {0, 0, 0}, grad_normal = {0, 0, 0};
				for (const auto &f : p->getFaces())
				{
					auto A = f->getArea();
					auto n = f->getNormalTuple();
					TotalArea += A;
					auto [p0, p1, p2] = f->getPointsTuple();
#if defined(Dombre2019)
#if defined(ENHANCE_CORNER_ACCURACY)
					// if (!p->CORNER || (p->CORNER && f->Dirichlet)) {
					// }
					// if (p->CORNER) {
					// 	if (f->Dirichlet) {
					// 		//@ 角点付近の場合，角の法線流速は怪しいので使わない．
					// 		//@ 角点なので，ノイマンを尊重すべきでは？
					// 		TotalArea_for_phi += f->getArea();
					// 		grad_tangential += grad_phi_tangential(f) * f->getArea();
					// 	}
					// }
					// else {
					//@角の場合は，Neumannは考慮しない
					TotalArea_for_phi += A;
					grad_tangential += grad_phi_tangential(f) * A;
					// }

					//@ 角点付近の場合，角の法線流速は怪しいので使わない．
					//@ 角点なので，ノイマンを尊重すべきでは？

					// case 1
					//  {
					//  	auto [p0, p1, p2] = f->getPointsTuple();
					//  	auto f_phin = (std::get<1>(p0->phiphin) + std::get<1>(p1->phiphin) + std::get<1>(p2->phiphin)) / 3.;
					//  	TotalArea_for_phin += f->getArea();
					//  	grad_normal += f_phin * f->getArea() * f->getNormalTuple();
					//  }
					//  case2
					if (p->CORNER)
					{
						if (f->Dirichlet)
						{
							TotalArea_for_phin += A;
							grad_normal += std::get<1>(p->phiphin) * A * n;
						}
					}
					else
					{
						TotalArea_for_phin += A;
						grad_normal += std::get<1>(p->phiphin) * A * n;
					}
#else
					grad += grad_phi_tangential(f) * A + std::get<1>(p->phiphin) * A * n;
#endif
#elif defined(Dombre2019_2)
					grad += grad_phi_tangential(f) * A;
#elif defined(AngleWeighted)
					Theta = MyVectorAngle(p2->getXtuple() - p0->getXtuple(), p1->getXtuple() - p0->getXtuple());
					TotalTheta += Theta;
					grad += grad_phi_tangential(f) * Theta;
#elif defined(AngleWeighted2)
					Theta = MyVectorAngle(p2->getXtuple() - p0->getXtuple(), p1->getXtuple() - p0->getXtuple());
					TotalTheta += Theta;
					grad += grad_phi_tangential(f) * Theta + std::get<1>(p->phiphin) * Theta * n;
					// grad += grad_phi * Theta;
#endif
				}
#if defined(Dombre2019)
#if defined(ENHANCE_CORNER_ACCURACY)
				if (it == pointsbegin)
					std::cout << Red << "Dombre2019と同じ方法で，流速を計算する：面の∇Φの面積重み平均" << reset << std::endl;
				grad_tangential /= TotalArea_for_phi;
				grad_normal /= TotalArea_for_phin;
				grad = grad_tangential + grad_normal;
#else
				/*
				Dombre2019の方法は周辺の節点のΦnを用いる．
				*/
				if (it == pointsbegin)
					std::cout << Red << "Dombre2019と同じ方法で，流速を計算する：面の∇Φの面積重み平均" << reset << std::endl;
				grad /= TotalArea;
#endif

#elif defined(Dombre2019_2)
				if (it == pointsbegin)
					std::cout << Red << "Dombre2019と同じ方法で，流速を計算する：面の∇Φの面積重み平均" << reset << std::endl;
				grad /= TotalArea;
				grad += std::get<1>(p->phiphin) * p->getNormalAreaAveraged();
#elif defined(AngleWeighted)
				/*
				各節点のΦnは境界値問題を解くことで決まっているので，
				Dombre2019の方法とは違って，周辺の点のΦnは用いなくても∇Φは計算できる．
				*/
				if (it == pointsbegin)
					std::cout << Red << "Dombre2019とは違うが，似た方法で流速を計算する：面の∇Φの角度重み平均" << reset << std::endl;
				grad /= TotalTheta;
				grad -= Dot(grad, p->getNormalTuple()) * p->getNormalTuple();
				grad += std::get<1>(p->phiphin) * p->getNormalTuple();
#elif defined(AngleWeighted2)
				/*
				各節点のΦnは境界値問題を解くことで決まっているので，
				Dombre2019の方法とは違って，周辺の点のΦnは用いなくても∇Φは計算できる．
				*/
				if (it == pointsbegin)
					std::cout << Red << "Dombre2019とは違うが，似た方法2で流速を計算する：面の∇Φの角度重み平均" << reset << std::endl;
				grad /= TotalTheta;
#endif
#endif
				/* ------------------------------------------------------ */
				if (!isFinite(grad))
				{
					std::cout << "gradはfiniteではない！！" << std::endl;
					std::cout << "p->phiphin = " << p->phiphin << std::endl;
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
				}
				p->grad_phi_BEM = grad;
				p->U_BEM = p->grad_phi_BEM;
				//! Uの修正
				if ((p->Neumann || p->CORNER) && !p->getContactFaces().empty())
				{
					for (const auto &f : p->getContactFaces())
					{
						auto n = f->getNormalTuple();
						p->U_BEM -= Dot(p->U_BEM, n) * n;
						p->U_BEM += phin(f) * n;
					}
				}
				p->U_tangential_BEM = p->U_BEM - Dot(p->U_BEM, n) * n;
				p->U_normal_BEM = Dot(p->U_BEM, n) * n;
			}
		}
		/* ------------------------------------------------------ */
#ifdef derivatives_debug
		std::cout << "ラプラシアンを計算" << std::endl;
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			auto ps = Flatten(BFS(p, 2, {p->getNetwork()}));
			auto intp = InterpolationVectorRBF(obj3D::extractX(ps), extVelocities(ps), p->getX());
			auto tmp = intp.laplacian(p->getX());
			p->laplacian_U_BEM = {tmp[0], tmp[1], tmp[2]};
		}
		/* ------------------------------------------------------ */
		// Tddd mean_position = Points[0]->getNetwork()->getMeanX();
		// Tddd mean_position = {0, 0, 0.5}; // Points[0]->getNetwork()->getMeanX();
#ifdef derivatives_debug
		std::cout << "DphiDtを計算" << std::endl;
#endif
		double gamma = 72.75 * 1E-3; //[N/m] 水20度
		double gravity = 9.80665;	 // [m/s2]
		double density = 997.;		 // [kg/m3]
		double nu = 0.01005 / density;
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (auto it = Points.begin(); it != Points.end(); ++it)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			auto p = *it;
			this->P_kappa[p] = p->kappa_BEM;
			this->P_gradPhi[p] = p->U_BEM;
			this->P_gradPhi_tangential[p] = p->U_tangential_BEM;
			this->P_phin_vector[p] = p->U_normal_BEM;
			this->P_laplacian[p] = p->laplacian_U_BEM;

			//!この微分は何のために評価したのか？
			//! ルンゲクッタは，この微分の値，grad phiとDphiDtを複数使って，DxDtとDphiDtの積分精度を向上させる
			// this->P_dxdt[p] = p->U_BEM; //流速
			// double fact = 1 + tanh(time_step / 2. - 3.);
			// auto nondim gamma* = gamma/g/rhp/L^2 = 0.007275 * 1E-3
			// auto c = (1 / 2. + tanh(real_time / 3. - 3.) / 2.);

			auto lap_tang = nu * p->laplacian_U_BEM - Dot(nu * p->laplacian_U_BEM, p->normal_BEM) * p->normal_BEM;
			// auto lap_tang = 100. * nu * p->laplacian_U_BEM;
			p->U_mod_BEM = p->U_BEM - lap_tang;
			// p->U_BEM = p->U_mod_BEM;
			this->P_dxdt_mod[p] = p->U_mod_BEM; //流速
			// 2021/08/08

			if ((p->Neumann || p->CORNER) && !p->getContactFaces().empty())
			{
				for (const auto &f : p->getContactFaces())
				{
					auto n = f->getNormalTuple();
					p->U_BEM -= Dot(p->U_BEM, n) * n;
					p->U_BEM += phin(f) * n;
				}
			}

			Tddd tension = {0, 0, 0};
			double area = 0.;
			double Area_tot = 0.;
			auto p_next_X = p->getXtuple() + p->U_BEM * dt;
			auto p_current_X = p->getXtuple();
			// なぜ面に依存した張力によって補正すべきなのか？
			// 線ではだめ，
			for (const auto &f : p->getFaces())
			{
				if ((p->CORNER && f->Dirichlet) || (!p->CORNER))
				{
					Tddd f_next_X = {0, 0, 0};
					Tddd f_current_X = {0, 0, 0};
					auto [fp0, fp1, fp2] = f->getPointsTuple();
					auto fp0_next_X = fp0->getXtuple() + fp0->U_BEM * dt;
					auto fp1_next_X = fp1->getXtuple() + fp1->U_BEM * dt;
					auto fp2_next_X = fp2->getXtuple() + fp2->U_BEM * dt;
					f_next_X = fp0_next_X / 3. + fp1_next_X / 3. + fp2_next_X / 3.;
					area = TriangleArea(T3Tddd{fp0_next_X, fp1_next_X, fp2_next_X});
					Area_tot += area;
					f_current_X += fp0->getXtuple();
					f_current_X += fp1->getXtuple();
					f_current_X += fp2->getXtuple();
					f_current_X /= 3.;
					{
						// auto v = (f_next_X + f_current_X) / 2 - (p_next_X + p_current_X) / 2;
						// auto v = f_next_X - p_next_X;
						// auto dist = Norm(v);
						// auto k = dist / dt;
						// // k += k / 5. * k / 5. * k / 5.;
						// tension += k * Normalize(v) * f->getAngle(p) / Theta_tot;
						/* ------------------------------------------------------ */
						// tension += (area / dt) * Normalize(f_next_X - p_next_X);
						tension += (area / dt) * (f_next_X - p_next_X);
					}
				}
			}
			tension /= Area_tot;
			if (!p->CORNER)
				tension -= Dot(tension, p->getNormalAreaAveraged()) * p->getNormalAreaAveraged();

			Tddd dumping;
			if (!p->CORNER)
				dumping -= p->U_tangential_BEM;
			P_tension[p] = tension;

			p->U_update_BEM = p->U_BEM;
			//! tensionが強いと角でとがるようなエラーがでる．
			double K = 0.1;
			p->U_update_BEM = p->U_BEM + K * tension; // p->U_normal_BEM + tangent; //要修正．EMTをプログラム
			/* ------------------------------------------------------ */
			if ((p->Neumann || p->CORNER) && !p->getContactFaces().empty())
			{
				for (const auto &f : p->getContactFaces())
				{
					auto n = f->getNormalTuple();
					p->U_BEM -= Dot(p->U_BEM, n) * n;
					p->U_BEM += phin(f) * n;
					//
					p->U_update_BEM -= Dot(p->U_update_BEM, n) * n;
					p->U_update_BEM += phin(f) * n;
				}
			}

			this->P_dxdt[p] = p->U_update_BEM; //流速
			// //!この場合マイナスでないと，上の部分が半たんする
			auto DphiDt = p->DphiDt_EMT(p->U_update_BEM, 0.);
			auto DphiDt_true = p->DphiDt(0.);

			//!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！
			this->P_DphiDt[p] = DphiDt; //全ての場所でこっちがまず求まる
			double Vphi_Vphi = Dot(p->U_BEM, p->U_BEM);
			auto aphiat = DphiDt_true - Vphi_Vphi;
			this->P_aphiat[p] = aphiat;
			this->P_pressure[p] = 10000. * (-1 / 2. * Vphi_Vphi - gravity * (std::get<2>(p->getXtuple())) - aphiat);
			//ここで圧力の項が抜けているが，これは全く流速に関係ないことに気づく．
			//なぜなら，表面上のどこでも同じだけ増加に寄与する大気圧は，
			// grad phiの計算によって，定数のため相殺されるからだ．
		}
#ifdef derivatives_debug
		std::cout << "derivatives終了" << std::endl;
#endif
	}
};
//* ------------------------------------------------------ */
//*                        境界値問題を解く                   */
//* ------------------------------------------------------ */
// #define solve_equations_on_all_points
#define solve_equations_on_all_points_rigid_mode
#define solveBVP_debug

void addIG(std::map<std::pair<netP *, netF *>, Tdd> &PF_phiphin, const netPp P_IN, const Tdd &igign)
{
	// #ifdef ENHANCE_CORNER_ACCURACY
	// 	if (P_IN->CORNER) {
	// 		//@ 角点は使えないので，係数を散らす．これはやめる
	// 		// auto Atot = Total(P_IN->getFaceAreas());
	// 		// for (const auto& F : P_IN->getFaces())
	// 		// {//! 分配
	// 		// 	auto A = F->getArea();
	// 		// 	if (PF_phiphin.count({ nullptr, F }))
	// 		// 		std::get<0>(PF_phiphin[{nullptr, F}]) += std::get<0>(igign) * A / Atot; //phiは忘れずに計算
	// 		// 	else
	// 		// 		PF_phiphin[{nullptr, F}] = { std::get<0>(igign) * A / Atot, 0. }; //phiは忘れずに計算
	// 		// }
	// 		if (Norm(igign) > 1E-10) {
	// 			double w = 0.;
	// 			for (const auto& F : P_IN->getFaces())
	// 				w += 1. / Norm(F->getXtuple() - P_IN->getXtuple());
	// 			for (const auto& F : P_IN->getNeighbors())
	// 				w += 1. / Norm(F->getXtuple() - P_IN->getXtuple());
	// 			for (const auto& F : P_IN->getFaces())
	// 			{//! 分配
	// 				auto R = 1. / Norm(F->getXtuple() - P_IN->getXtuple());
	// 				if (PF_phiphin.count({ nullptr, F }))
	// 					std::get<0>(PF_phiphin[{nullptr, F}]) += std::get<0>(igign) * R / w; //phiは忘れずに計算
	// 				else
	// 					PF_phiphin[{nullptr, F}] = { std::get<0>(igign) * R / w,0. }; //phiは忘れずに計算
	// 			}
	// 			for (const auto& F : P_IN->getNeighbors())
	// 			{//! 分配
	// 				auto R = 1. / Norm(F->getXtuple() - P_IN->getXtuple());
	// 				addIG(PF_phiphin, F, igign * R / w);
	// 				// if (PF_phiphin.count({ F,nullptr }))
	// 				// 	std::get<0>(PF_phiphin[{F, nullptr}]) += std::get<0>(igign) * R / w; //phiは忘れずに計算
	// 				// else
	// 				// 	PF_phiphin[{F, nullptr}] = { std::get<0>(igign) * R / w,0. }; //phiは忘れずに計算
	// 			}
	// 		}
	// 	}
	// 	else
	// #endif
	{
		if (PF_phiphin.count({P_IN, nullptr}))
			std::get<0>(PF_phiphin[{P_IN, nullptr}]) += std::get<0>(igign); // phiは忘れずに計算
		else
			PF_phiphin[{P_IN, nullptr}] = {std::get<0>(igign), 0.}; // phiは忘れずに計算
	}
}

void addIGn(std::map<std::pair<netP *, netF *>, Tdd> &PF_phiphin, const netPp P_IN, const Tdd &igign)
{
	// #ifdef ENHANCE_CORNER_ACCURACY
	// if (P_IN->CORNER) {
	// 	//@ 角点は使えないので，係数を散らす
	// 	auto Atot = Total(P_IN->getFaceAreas());
	// 	// for (const auto& F : P_IN->getFaces())
	// 	// {//! 分配
	// 	// 	auto A = F->getArea();
	// 	// 	if (PF_phiphin.count({ nullptr, F }))
	// 	// 		std::get<1>(PF_phiphin[{nullptr, F}]) += std::get<1>(igign) * A / Atot; //phiは忘れずに計算
	// 	// 	else
	// 	// 		PF_phiphin[{nullptr, F}] = { 0.,std::get<1>(igign) * A / Atot }; //phiは忘れずに計算
	// 	// }
	// 	if (Norm(igign) > 1E-10) {
	// 		double w = 0.;
	// 		for (const auto& F : P_IN->getFaces())
	// 			w += 1. / Norm(F->getXtuple() - P_IN->getXtuple());
	// 		for (const auto& F : P_IN->getNeighbors())
	// 			w += 1. / Norm(F->getXtuple() - P_IN->getXtuple());
	// 		for (const auto& F : P_IN->getFaces())
	// 		{//! 分配
	// 			auto R = 1. / Norm(F->getXtuple() - P_IN->getXtuple());
	// 			if (PF_phiphin.count({ nullptr, F }))
	// 				std::get<1>(PF_phiphin[{nullptr, F}]) += std::get<1>(igign) * R / w; //phiは忘れずに計算
	// 			else
	// 				PF_phiphin[{nullptr, F}] = { 0.,std::get<1>(igign) * R / w }; //phiは忘れずに計算
	// 		}
	// 		for (const auto& F : P_IN->getNeighbors())
	// 		{//! 分配
	// 			auto R = 1. / Norm(F->getXtuple() - P_IN->getXtuple());
	// 			addIGn(PF_phiphin, F, igign * R / w);
	// 			// if (PF_phiphin.count({ F,nullptr }))
	// 			// 	addIGn(PF_igign, p_IN, igign * R / w);
	// 			// // std::get<1>(PF_phiphin[{F, nullptr}]) += std::get<1>(igign) * R / w; //phiは忘れずに計算
	// 			// else
	// 			// 	PF_phiphin[{F, nullptr}] = { 0.,std::get<1>(igign) * R / w }; //phiは忘れずに計算
	// 		}
	// 	}
	// }
	// else
	// #endif
	{
		if (PF_phiphin.count({P_IN, nullptr}))
			std::get<1>(PF_phiphin[{P_IN, nullptr}]) += std::get<1>(igign); // phiは忘れずに計算
		else
			PF_phiphin[{P_IN, nullptr}] = {0., std::get<1>(igign)}; // phiは忘れずに計算
	}
}

void addIGn(std::map<std::pair<netP *, netF *>, Tdd> &PF_phiphin, const netFp F_IN, const Tdd &igign)
{
	{
		// auto [p0, p1, p2] = F_IN->getPointsTuple();
		// addIGn(PF_phiphin, p0, igign / 3.);
		// addIGn(PF_phiphin, p1, igign / 3.);
		// addIGn(PF_phiphin, p2, igign / 3.);
		if (PF_phiphin.count({nullptr, F_IN}))
			std::get<1>(PF_phiphin[{nullptr, F_IN}]) += std::get<1>(igign); // phiは忘れずに計算
		else
			PF_phiphin[{nullptr, F_IN}] = {0., std::get<1>(igign)}; // phiは忘れずに計算
	}
}

void solveBVP(const std::unordered_set<networkPoint *> &PointsIN, const std::unordered_set<networkFace *> &FacesIN)
{
	auto Points = ToVector(PointsIN);
	auto Faces = ToVector(FacesIN);
	using map_P_Vd = std::map<netP *, V_d>;

	// for (auto it = PointsIN.begin(); it != PointsIN.end();++it)
	// {
	// 	// std::cout << "曲率の計算" << std::endl;
	// 	auto p = *it;
	// 	if (!isFinite(p->phiphin))
	// 	{
	// 		std::cout << "p->phiphinはfiniteではない！！" << std::endl;
	// 		std::cout << "p->phiphin = " << p->phiphin << std::endl;
	// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	// 	}
	// }

	// #ifdef solve_equations_on_all_points_rigid_mode
	//* ------------------------------------------------------ */
	//%                     各点で方程式を作る場合                */
	//* ------------------------------------------------------ */
	std::cout << "各点で方程式を作る場合" << std::endl;
	map_pairPF_pairPF_Tdd P_P_IGIGn; //途中でreverseしているので注意
	int ii = 0;
	for (const auto &pO : Points)
	{
#ifdef ENHANCE_CORNER_ACCURACY
		if (false)
			if (pO->CORNER)
				for (const auto &fO : pO->getFaces())
					P_P_IGIGn[{nullptr, fO}] = {}; //! P_P_IGIGn[row][col] = {IG,IGn}
		P_P_IGIGn[{pO, nullptr}] = {};			   //! P_P_IGIGn[row][col] = {IG,IGn}
#else
		P_P_IGIGn[{pO, nullptr}] = {}; //! P_P_IGIGn[row][col] = {IG,IGn}
#endif
	}
#ifdef _OPENMP
	std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;
#pragma omp parallel
#endif
	// for (auto it = P_P_IGIGn.begin(); it != P_P_IGIGn.end();++it)

	for (auto &[pairPF, P_igign] : P_P_IGIGn)
	{
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			auto [pO, fO] = pairPF;
			// #ifdef ENHANCE_CORNER_ACCURACY
			// 			if (fO)
			// 			{
			// 				//! 角点が和となることはやめようか．
			// 				//@ 角点の面が該当するが，角点自身の行の方程式は該当しない．すでに方程式は与えている．
			// 				//* ------------------------------------------------------ */
			// 				//*                    面上で方程式を作る場合                  */
			// 				//* ------------------------------------------------------ */
			// 				for (const auto &integrating_f : Faces)
			// 				{
			// 					for (const auto &[p, igign] : BEM::calc_P_IGIGnLinearTuple(integrating_f, fO))
			// 					{
			// 						//@ 面の頂点pOとpが一致した場合，pOはpOの周辺面のΦの面積平均で表されるので，igignに面積比をかける．
			// 						//@ { p, nullptr }はないので，分配する．Φ[{ p, nullptr }]= (Af0*Φ[{ p, f0 }]+Af1*Φ[{ p, f1 }]+Af2*Φ[{ p, f2 }])/A_totのように考える
			// 						addIG(P_igign, p, igign);
			// 						if (integrating_f != fO)
			// 						{
			// 							addIGn(P_igign, p, igign);
			// 							//@ 対角成分は，リジッドモード法を使って計算できる
			// 							// addIGn(P_igign, fO, -igign);
			// 							// やはり， nullptr,Fが必要なのでは？
			// 						}
			// 					}
			// 				}
			// 			}
			// 			else
			// #endif
			{
				//* ------------------------------------------------------ */
				//*                     点で方程式を作る場合                  */
				//* ------------------------------------------------------ */
				for (const auto &integrating_f : Faces)
				{
					for (const auto &[p, igign] : BEM::calc_P_IGIGnLinearTuple(integrating_f, pO))
					{
						addIG(P_igign, p, igign);
						if (p != pO)
						{
							addIGn(P_igign, p, igign);
							//@ 対角成分は，リジッドモード法を使って計算できる
							addIGn(P_igign, pO, -igign);
						}
					}
				}
			}
		}
	}

	// 面に与えるので，
	// 	multiple_phiphinはには面の値をいれるといい．
	// 	1 / 3の値
	// 	φnにはすなおにノイマン条件
	// 角点のΦnは周囲のディリクレ面から計算する．

#ifdef ENHANCE_CORNER_ACCURACY
	for (auto &[pairPF, P_igign] : P_P_IGIGn)
	{
		auto [pO, fO] = pairPF;
		if (fO)
		{
			for (const auto &P : fO->getPoints())
				if (P_igign.count({P, nullptr}))
					std::get<1>(P_igign[{P, nullptr}]) -= -2. * M_PI / 6.;
				else
					P_igign[{P, nullptr}] = {0, 2. * M_PI / 6.};
			//
			if (P_igign.count({nullptr, fO}))
				std::get<1>(P_igign[{nullptr, fO}]) -= -2. * M_PI / 2.;
			else
				P_igign[{nullptr, fO}] = {0, 2. * M_PI / 2.};
		}
	}
#endif

	// for (const auto& p : corner_points) {
	// 	P_P_IGIGn[{p, nullptr}].clear();
	// 	auto& igign = P_P_IGIGn[{p, nullptr}][{p, nullptr}];
	// 	auto Atot = Total(extractAreas(p->getFaces()));
	// 	for (const auto& f : p->getFaces())
	// 	{
	// 		P_P_IGIGn[{p, nullptr}][{p, f}] = { 0,f->getArea() / Atot };
	// 	}
	// 	P_P_IGIGn[{p, nullptr}][{p, nullptr}] = { 0,-1 };
	// }
	// for (auto it = P_P_IGIGn.begin(); it != P_P_IGIGn.end();++it)
	// {
	// 	auto neighbor_f = std::get<1>(it->first);
	// 	if (neighbor_f) {
	// 		auto origin = std::get<0>(it->first);
	// 		auto& P_igign = it->second;/*P_Tdd*/
	// 		auto [p0, p1, p2] = neighbor_f->getPointsTuple();
	// 		std::get<1>(P_igign[{p0, nullptr}]) -= -2. * M_PI / 3;
	// 		std::get<1>(P_igign[{p1, nullptr}]) -= -2. * M_PI / 3;
	// 		std::get<1>(P_igign[{p2, nullptr}]) -= -2. * M_PI / 3;
	// 		std::get<1>(P_igign[{origin, neighbor_f}]) -= -2. * M_PI / 2;
	// 	}
	// }
	std::cout << "並列化 DONE" << std::endl;

	std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;
	std::cout << "未知変数の係数行列は左，既知変数の係数行列は右" << std::endl;

	for (auto &[_, p_igign] : P_P_IGIGn)
	{ // each row
		for (auto &[pairPF, igign] : p_igign)
		{ // each column
			// この入れ替えが，左辺と右辺の入れ替えを意味する
			// igign[0]が左辺，未知となる
			// igign[1]が右辺，既知となる
			// CORNERの点はDとして計算している
			auto [pO, fO] = pairPF;
			if (fO)
			{
				if (fO->Neumann)
				{ /*Neumann*/
					igign = {-std::get<1>(igign), -std::get<0>(igign)};
				}
			}
			else
			{
				if (pO->Neumann)
				{ /*Neumann*/
					igign = {-std::get<1>(igign), -std::get<0>(igign)};
				}
			}
		}
	}

	/**
	 * 変数ベクトル knowns について:
	 * 境界条件をisD(),isN()で確認し，
	 * isD()=true ->　P_igign[p][0]がknown
	 * isN()=true ->　P_igign[p][1]がknown
	 * となる
	 */
	//! Pointsの順番と合わせてとるように注意
	/* ------------------------------------------------------ */
	// <-- Points size -->
	// {{},{},{},{},{}} . {<--  face size -->} = {<-- point size -->}

	// V_d knowns = V_d(Points.size()); //!ok
	//方程式の数は，pointsだけある
	VV_d mat_ukn(P_P_IGIGn.size(), V_d(P_P_IGIGn.size(), 0.)); //! ok
	VV_d mat_kn(P_P_IGIGn.size(), V_d(P_P_IGIGn.size(), 0.));  //! ok

	// 順番を間違えないようにベクトルを作成
	std::vector<std::pair<networkPoint *, networkFace *>> pair_PF(P_P_IGIGn.size());
	auto i = 0;
	for (auto it = P_P_IGIGn.begin(); it != P_P_IGIGn.end(); ++it)
		pair_PF[i++] = it->first;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i = 0; i < pair_PF.size(); i++)
	{
		auto tmp = P_P_IGIGn.at(pair_PF[i]);
		for (auto j = 0; j < pair_PF.size(); j++)
		{
			if (tmp.find(pair_PF[j]) != tmp.end())
			{
				auto igign = tmp.at(pair_PF[j]);
				mat_ukn[i][j] = std::get<0>(igign); //左辺，phinがwaterのunknown->IG
				mat_kn[i][j] = std::get<1>(igign);	//右辺
			}
		}
	}

	if (!isFinite(mat_ukn))
	{
		std::cout << "mat_ukn = " << mat_ukn << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "解なし");
	}

	if (!isFinite(mat_kn))
	{
		std::cout << "mat_kn = " << mat_ukn << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "解なし");
	}

	/**
	 * Dot(mat_ukn,phiORphin) = Dot(mat_kn,knowns)
	 * => phiORphin = Dot(mat_ukn^-1, Dot(mat_kn,knowns))
	 */
	V_d knowns = V_d(pair_PF.size()); //! ok
	//! Pointsの順番と合わせてとるように注意
	// for (auto i = 0; i < Points.size(); i++)
	// 	knowns[i] = Points[i]->phiphin[Points[i]->Neumann ? 1 /*phin*/ : 0 /*phi*/];

	for (auto i = 0; i < pair_PF.size(); i++)
	{
		auto [p, f] = pair_PF[i];
		if (f)
		{
			if (f->Neumann)
				knowns[i] = std::get<1>(f->phiphin);
			else
				knowns[i] = std::get<0>(f->phiphin);
		}
		else
		{
			if (p->Neumann)
				knowns[i] = std::get<1>(p->phiphin);
			else
				knowns[i] = std::get<0>(p->phiphin);
		}
	}
	std::cout << "--------------------- 境界積分方程式を解く ---------------------" << std::endl;

	/* ------------------------------------------------------ */
	V_d phiORphin(knowns.size());
	std::cout << "phiORphin.size()= " << phiORphin.size() << std::endl;
	std::cout << "  mat_kn.size() = " << mat_kn.size() << std::endl;
#if defined(solve_equations_on_all_points) || defined(solve_equations_on_all_points_RBF) || defined(solve_equations_on_all_points_rigid_mode)
	//* 未知変数の計算
	std::cout << "LU decomposition" << std::endl;
	ludcmp_parallel lu(mat_ukn /*未知の行列係数（左辺）*/);
	lu.solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#else
	//* 未知変数の計算
	std::cout << "SVD decomposition" << std::endl;
	SVD svd(mat_ukn);
	svd.solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#endif
	i = 0;
	for (auto i = 0; i < pair_PF.size(); i++)
	{
		auto p = std::get<0>(pair_PF[i]);
		auto f = std::get<1>(pair_PF[i]);
		if (f)
		{
			if (f->Neumann)
				std::get<0>(f->phiphin) = phiORphin[i];
			else
				std::get<1>(f->phiphin) = phiORphin[i];
		}
		else
		{
			if (p->Neumann)
				std::get<0>(p->phiphin) = phiORphin[i];
			else
				std::get<1>(p->phiphin) = phiORphin[i];
		}
	}

	if (!isFinite(phiORphin))
	{
		std::cout << "phiORphin = " << phiORphin << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "解なし");
	}

	// #ifdef ENHANCE_CORNER_ACCURACY
	// 	for (const auto& p : Points)
	// 		if (p->CORNER && p->Neumann) {
	// 			Tdd phiphin = { 0,0 };
	// 			for (const auto& f : p->getFaces())
	// 				phiphin += f->getArea() * f->phiphin;
	// 			phiphin /= Total(p->getFaceAreas());
	// 			std::get<0>(p->phiphin) = std::get<0>(phiphin);
	// 		}
	// 		else if (p->CORNER && p->Dirichlet) {
	// 			Tdd phiphin = { 0,0 };
	// 			for (const auto& f : p->getFaces())
	// 				phiphin += f->getArea() * f->phiphin;
	// 			phiphin /= Total(p->getFaceAreas());
	// 			std::get<1>(p->phiphin) = std::get<1>(phiphin);
	// 		}
	// #endif

	// とりあえず平均とする．
	// for (const auto& p : Points)
	// 	if (!p->multiple_phiphin.empty()) {
	// 		double phi = 0, phin = 0;
	// 		for (const auto& [F, phiphin] : p->multiple_phiphin)
	// 			phi += std::get<0>(phiphin);
	// 		std::get<0>(p->phiphin) = phi / (double)(p->multiple_phiphin.size());
	// 		double A_tot = 0;
	// 		for (const auto& [F, phiphin] : p->multiple_phiphin) {
	// 			auto A = F->getArea();
	// 			A_tot += A;
	// 			phin += std::get<1>(phiphin) * A;
	// 		}
	// 		std::get<1>(p->phiphin) = phin / A_tot;
	// 	}

	// for (auto i = 0; i < Points.size(); i++)
	// 	Points[i]->phiphin[Points[i]->Neumann ? 0 /*phi*/ : 1 /*Dの場合phinが解として得られている*/] = phiORphin[i] /*解*/;
};

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

auto Pd2PVd = [](map_P_d &P_kappa)
{
	map_P_Vd ret;
	for (auto &[p, kappa] : P_kappa)
		ret[p] = {kappa};
	return ret;
};

V_d laplacianV(const netPp p)
{
	auto ps = Flatten(BFS(p, 2, {p->getNetwork()}));
	auto intp = InterpolationVectorRBF(obj3D::extractX(ps), extVelocities(ps), p->getX());
	return intp.laplacian(p->getX());
};

V_d meanX(const V_netPp ps)
{
	V_d ret(3, 0.);
	for (const auto &p : ps)
		ret += p->getX();
	return ret / ((double)(ps.size()));
};
void normalizePhi(const V_netPp ps)
{
	double ret = 0.;
	for (const auto &p : ps)
		ret += std::get<0>(p->phiphin);
	ret /= ((double)(ps.size()));
	for (const auto &p : ps)
		std::get<0>(p->phiphin) -= ret;
};
/* ------------------------------------------------------ */
// V_d moveToOptimumX(networkPoint * p)
// {
// 	// 自分あり
// 	VV_d Xs = { p->getX() };
// 	VV_d PhiPhins = { p->phiphin };
// 	VV_d params = { {0., 0.} };
// 	// 角度の考慮なし
// 	for (const auto& tup : p->getNeighbors_Depth2_OnPolarAsTuple())
// 	{
// 		auto [t0, t1, X, q] = tup;
// 		Xs.emplace_back(X);
// 		PhiPhins.emplace_back(q->phiphin);
// 		params.emplace_back(V_d{ t0, t1 });
// 	}
// 	auto intpX = InterpolationVectorRBF(params, Xs);

// 	double smallest_angle_of_best = 0;
// 	V_d selected_t0t1 = { 0., 0. };
// 	V_d min_angles = {};
// 	V_d T0T1 = Subdivide(-0.05, 0.05, 10);
// 	for (const auto& t0 : T0T1)
// 		for (const auto& t1 : T0T1)
// 		{
// 			p->setX(intpX({ t0, t1 }));
// 			min_angles.clear();
// 			for (const auto& f : p->getFaces())
// 				min_angles.emplace_back(Min(f->getAngles()));
// 			if (Min(min_angles) > smallest_angle_of_best /*最小の角度がより大きい場合は交換*/)
// 			{
// 				smallest_angle_of_best = Min(min_angles);
// 				selected_t0t1 = { t0, t1 };
// 			}
// 		}
// 	auto ret = intpX(selected_t0t1);
// 	p->setX(ret); //最適な位置へ移動
// 	return ret;
// };

// b! ------------------------------------------------------ */
// b!           格子のdivide, merge．それに伴うΦ，Φnの付与       */
// b! ------------------------------------------------------ */
Tdd phiphin_from_faces_IDW(const networkLine *const l)
{
	double TotalWphi = 0., TotalWphin = 0., phi = 0., phin = 0.;
	std::unordered_set<networkPoint *> points;
	auto [p0, p1] = l->getPointsTuple();
	Tddd Xofl = (p0->getXtuple() + p1->getXtuple()) / 2.;
	auto fs = l->getFaces();
	for (const auto &f : fs)
		for (const auto &p : f->getPoints())
			points.emplace(p);

	for (const auto &q : points)
	{
		auto w = 1. / Distance(q, Xofl);
		phi += w * std::get<0>(q->phiphin);
		TotalWphi += w;
		phin += w * std::get<1>(q->phiphin);
		TotalWphin += w;
	}
	return {phi / TotalWphi, phin / TotalWphin};
};
Tdd phiphin_from_faces(const networkLine *const l)
{
	Tdd phiphin = {0., 0.};
	for (const auto &f : l->getFaces())
	{
		auto [p0, p1, p2] = f->getPointsTuple();
		phiphin += (p0->phiphin + p1->phiphin + p2->phiphin) / 6.;
	}
	return phiphin;
};
Tdd phiphin_from_points(const networkLine *const l)
{
	Tdd phiphin = {0., 0.};
	for (const auto &p : l->getPoints())
		phiphin += p->phiphin / 2.;
	return phiphin;
};

Tdd phiphin_from_points_faces(const networkLine *const l)
{
	return phiphin_from_points(l) / 2. + phiphin_from_faces(l) / 2.;
};

void remesh(Network &water)
{
	std::cout << "remesh" << std::endl;
	water.setBounds();
	double mean_length = Mean(extLength(water.getLines()));
	bool isfound = false, ismerged = false;
	;
	int count = 0;
	double lim_degree = 5.;
	do
	{
		// なくなるまでやるか？
		isfound = false;
		ismerged = false;
		auto points = water.getPoints();
		V_netPp ps;
		for (const auto &p : RandomSample(ToVector(points)))
		{
			//* ------------------------------------------------------ */
			//*                立体角が小さすぎる場合merge                 */
			//* ------------------------------------------------------ */
			if (p->getSolidAngle() < 0.05)
			{
				p->sortLinesByLength();
				auto l = p->getLines()[0];
				//@ case 1
				// Tdd phiphin = phiphin_from_faces_IDW(l);
				//@ case2 lの面全体を考慮
				// Tdd phiphin = phiphin_from_faces(l);
				//@ case3 lの点だけを考慮
				// Tdd phiphin = phiphin_from_points(l);
				//@ case4 lの面全体を考慮
				Tdd phiphin = phiphin_from_points_faces(l);
				/* ------------------------------------------------------ */
				auto q = l->merge();
				q->phiphin = phiphin;
				ismerged = true;
				break;
			}
			//! ------------------------------------------------------ */
			//!             辺の長さが長すぎるまたは短すぎる場合             */
			//! ------------------------------------------------------ */
			for (const auto &l : RandomSample(p->getLines()))
			{
				auto [p0, p1] = l->getPointsTuple();
				if (!((p0->Neumann && p1->Dirichlet) || (p0->Dirichlet && p1->Neumann)))
				{
					l->flipIfBetter(5.);
					//@ ------------------------------------------------------ */
					if (l->length() > mean_length * 3. / 2. /*長すぎる*/)
					{
						//@ case 1
						// Tdd phiphin = phiphin_from_faces_IDW(l);
						//@ case2 lの面全体を考慮
						// Tdd phiphin = phiphin_from_faces(l);
						//@ case3 lの点だけを考慮
						// Tdd phiphin = phiphin_from_points(l);
						//@ case4 lの面全体を考慮
						Tdd phiphin = phiphin_from_points_faces(l);
						/* ------------------------------------------------------ */
						auto q = l->divide();
						q->phiphin = phiphin;
						isfound = true;
						break;
					}
					//@ ------------------------------------------------------ */
					if (l->length() < mean_length / 4.)
					{
						//@ case1
						// Tdd phiphin = phiphin_from_faces_IDW(l);
						//@ case2 lの面全体を考慮
						// Tdd phiphin = phiphin_from_faces(l);
						//@ case3 lの点だけを考慮
						// Tdd phiphin = phiphin_from_points(l);
						//@ case4 lの面全体を考慮平均
						Tdd phiphin = phiphin_from_points_faces(l);
						/* ------------------------------------------------------ */
						auto q = l->merge();
						q->phiphin = phiphin;
						ismerged = true;
						break;
					}
					if (ismerged || isfound)
						break;
				}
			}
			if (ismerged || isfound)
				break;
		}
	} while ((ismerged || isfound) && count++ < 500);
	if (count == 10)
		std::cout << "remesh too much" << std::endl;
};
//////////////////////////////////////////////////////////////////////////////////////////////////
void setBoundaryConditions(Network &water, const Buckets<networkFace> &B)
{
	auto radius = Mean(extLength(water.getLines()));
	// b% -------------------------------------------------------- */
	// b%            境界条件（角点・ディリクレ・ノイマン）の決定         */
	// b% -------------------------------------------------------- */
	// b% step1 衝突の判定
	//!!! 衝突の判定がよくエラーが出る箇所
	for (const auto &p : water.getPoints())
	{
		//!ここも重要：点と面の衝突をどのようにすれば矛盾なく判定できるか．
		p->clearContactFaces();
		p->radius = 1. * radius;	  // Mean(extLength(p->getLines()));
		p->addContactFaces(B, false); /**shadowあり*/
	}
	// b% step2 面の境界条件を決定
	/*面Aの点が接触している面Bを取得．A,B面が向き合っていればノイマン*/
	for (const auto &f : water.getFaces())
	{
		// auto isDirichlet = facingFace(f).empty();
		auto isDirichlet = !closest_facingFace(f); // nullptrならDirichlet
		f->Neumann = !isDirichlet;
		f->Dirichlet = isDirichlet;
	}
	// b% step3 線の境界条件を決定
	for (const auto &l : water.getLines())
	{
		auto faces = l->getFaces();
		l->Neumann = std::all_of(faces.begin(), faces.end(), [](const auto &f)
								 { return f->Neumann; });
		l->Dirichlet = std::all_of(faces.begin(), faces.end(), [](const auto &f)
								   { return f->Dirichlet; });
		if (!l->Neumann && !l->Dirichlet)
			l->CORNER = true;
		else
			l->CORNER = false;
	}
	// b% step4 点の境界条件を決定
	/*周りの面が全てノイマンなら点はノイマン．
		周りの面がディリクレとノイマン両方を持っていれば角点
		　それ以外は，ディリクレ．*/
	for (const auto &p : water.getPoints())
	{
		auto faces = p->getFaces();
		bool isNeumann = std::all_of(faces.begin(), faces.end(), [](const auto &f)
									 { return f->Neumann; });
		bool isDirichlet = std::all_of(faces.begin(), faces.end(), [](const auto &f)
									   { return f->Dirichlet; });
		if (isNeumann)
			p->setN();
		else if (isDirichlet)
			p->setD();
		else
			p->setC();
	}
	//% ------------------------------------------------------ */
	// b* ------------------------------------------------------ */
	// b*       　    ノイマン境界の点や面にはΦnを与える               */
	// b* ------------------------------------------------------ */
	std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << reset << std::endl;
	for (const auto &p : water.getPoints())
	{
		if (p->Neumann || p->CORNER)
		{
			// TODO 多重に対応させるため，{p,ノイマンf}に対してΦnを与える
			std::get<1>(p->phiphin) = phin_from_Neumann_surface(p);
			if (!isFinite(p->phiphin))
			{
				std::cout << "p->phiphinはfiniteではない！！" << std::endl;
				if (p->Neumann)
					std::cout << "Neumann" << std::endl;
				if (p->Dirichlet)
					std::cout << "Dirichlet" << std::endl;
				if (p->CORNER)
					std::cout << "CORNER" << std::endl;
				std::cout << "forced_velocity(t) = " << forced_velocity(real_time) << std::endl;
				std::cout << "p->phiphin = " << p->phiphin << std::endl;
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			}
		}
	}
	std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << reset << std::endl;
	for (const auto &f : water.getFaces())
	{
		auto [p0, p1, p2] = f->getPointsTuple();
		auto phiphin = (p0->phiphin + p1->phiphin + p2->phiphin) / 3.;
		if (f->Neumann)
			std::get<1>(f->phiphin) = phin_contact(f);
		else
			std::get<0>(f->phiphin) = std::get<0>(phiphin);
	}
};
//////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
	try
	{
		// b* -------------------- メッシュ・設定の読み込み -------------------- */
		//  b! ------------------------------------------------------ */
		Print("------------------- ./tank.json -----------------");
		std::ifstream obj0_str("./tank.json");
		JSON json0(obj0_str);
		Print("------------------------------------");
		Network obj(json0["objfile"][0], json0["name"][0]);
		for (auto &p : obj.getPoints())
			p->radius = stod(json0["radius"])[0];
		Print("------------------------------------");
		modify(obj, json0);
		Print("------------------------------------");
		std::string tank_output_directory = json0["output_directory"][0];
		std::filesystem::create_directories(tank_output_directory);
		std::string tank_output_pvd_file_name = json0["output_pvd_file_name"][0];
		PVDWriter tankPVD(tank_output_directory + "/" + tank_output_pvd_file_name + ".pvd");
		std::string tank_output_vtu_file_name = json0["output_vtu_file_name"][0];
		// std::system("ls -l" + tank_output_directory);
		std::filesystem::copy_file("tank.json", tank_output_directory + "/tank.json", std::filesystem::copy_options::overwrite_existing);
		// b! ------------------------------------------------------ */
		Print("------------------- ./water.json -----------------");
		std::ifstream obj1_str("./water.json");
		JSON json1(obj1_str);
		Network water(json1["objfile"][0], json1["name"][0]);
		for (auto &p : water.getPoints())
			p->radius = stod(json1["radius"])[0];
		modify(water, json1);
		std::string water_output_directory = json1["output_directory"][0];
		std::filesystem::create_directories(water_output_directory);
		std::string water_output_pvd_file_name = json1["output_pvd_file_name"][0];
		PVDWriter waterPVD(water_output_directory + "/" + water_output_pvd_file_name + ".pvd");
		std::string water_output_vtu_file_name = json1["output_vtu_file_name"][0];
		// std::system("ls -l" + water_output_directory);
		std::filesystem::copy_file("water.json", water_output_directory + "/water.json", std::filesystem::copy_options::overwrite_existing);
		// b! ------------------------------------------------------ */
		//  b* ------------------------------------------------------ */
		//  b*                         メインループ                     */
		//  b* ------------------------------------------------------ */
		for (time_step = 0; time_step < 10000; time_step++)
		{
			//!体積を保存するようにリメッシュする必要があるだろう．
			auto radius = Mean(extLength(water.getLines()));
			Print("makeBucketFaces", Green);
			obj.makeBucketFaces(radius);
			Print("setBoundaryConditions", Green);
			setBoundaryConditions(water, obj.BucketFaces);
			remesh(water);
			// b# ------------------------------------------------------ */
			// b#                       刻み時間の決定                     */
			// b# ------------------------------------------------------ */
			const auto Points = water.getPoints();
			// mk_vtu(output_directory + "/getFacesSort2.vtu", {water.getPoints()[0]->getNeighborsSort2()});
			// mk_vtu(output_directory + "/Points.vtu", {water.getPoints()});
			const auto Faces = water.getFaces();
			// mk_vtu(water_output_directory + "/Faces.vtu", { water.getFaces() });
			auto p = (*std::min_element(Points.begin(), Points.end(),
										[](auto a, auto b)
										{
											auto A = Min(extLength(a->getLines())) / Norm(a->U_BEM);
											auto B = Min(extLength(b->getLines())) / Norm(b->U_BEM);
											if (!isFinite(A))
												return false;
											return A < B;
										}));

			double min_dt = Min(extLength(p->getLines())) / Norm(p->U_BEM) / 2.;
			double dt;
			if (time_step < 2)
				dt = 0.001;
			else if (min_dt < 0.01)
				dt = min_dt;
			else
				dt = 0.01;
			Print("===========================================================================");
			Print("time_step :" + Red + std::to_string(time_step) + reset);
			Print("real time :" + Red + std::to_string(real_time) + reset);
			Print("---------------------------------------------------------------------------");
			// b@ ------------------------------------------------------ */
			// b@        初期値問題を解く（時間微分方程式を数値積分する）        */
			// b@ ------------------------------------------------------ */
			std::map<netPp, RungeKutta_<double> *> P_RK_phi;
			std::map<netPp, RungeKutta_<Tddd> *> P_RK_X;
			// 面ベースで境界条件を判断する
			/* ------------------------------------------------------ */
			// for (const auto& p : Points)
			// 	p->multiple_phiphin.clear();
			for (const auto &p : Points)
			{
				P_RK_phi[p] = (new RungeKutta_(dt, real_time, std::get<0>(p->phiphin), 4));
				P_RK_X[p] = (new RungeKutta_(dt, real_time, p->getXtuple(), 4));
			}
			do
			{
				//! 壁面の動きは，マイステップ更新することにした．この結果はphin()で参照される
				auto RK_time = P_RK_phi[*Points.begin()]->gett(); //%各ルンゲクッタの時刻を使う
				std::cout << "RK_time = " << RK_time << ", real_time = " << real_time << std::endl;
				obj.velocity = forced_velocity(RK_time); // T6d //@ Φnを計算するために，物体表面の速度forced_velocityは，保存しておく必要がある
				// b% -------------------------------------------------------- */
				// b%            境界条件（角点・ディリクレ・ノイマン）の決定         */
				// b% -------------------------------------------------------- */
				auto radius = Mean(extLength(water.getLines()));
				Buckets<networkFace> B(obj.getBounds(), radius);
				for (const auto &f : obj.getFaces())
				{
					f->clearParametricPoints();
					f->particlize(radius / 3., {-radius / 2});
				}
				for (const auto &f : obj.getFaces())
					for (const auto &p : f->getParametricPoints())
						B.add(p->getXtuple(), f);
				setBoundaryConditions(water, B);
				// b* ------------------------------------------------------ */
				// b*           　境界値問題を解く-> {Φ,Φn}が決まる              */
				// b* ------------------------------------------------------ */
				std::cout << Green << "境界値問題を解く-> {Φ,Φn}が決まる" << reset << std::endl;
				solveBVP(Points, Faces);
				// b* ------------------------------------------------------ */
				// b*                    微分∇ΦやDUDtを計算                    */
				// b* ------------------------------------------------------ */
				std::cout << Green << "微分∇ΦやDUDtを計算" << reset << std::endl;
				derivatives ders(Points, dt);
				// b* ------------------------------------------------------ */
				// b*                 ディリクレ境界ではΦを時間積分               */
				// b* ------------------------------------------------------ */
				std::cout << Green << "ディリクレ境界ではΦを時間積分，ノイマン境界ではΦnを陽に与える" << reset << std::endl;
				for (const auto &p : Points)
				{
					//@ Φの時間発展，Φnの時間発展はない
					{
						auto rk = P_RK_phi[p];
						rk->push(ders.P_DphiDt[p]);
						if (p->CORNER || p->Dirichlet)
							std::get<0>(p->phiphin) = rk->getX(); // 角点の法線方向はわからないので，ノイマンの境界条件phinを与えることができない．
					}
					//@ 位置xの時間発展
					{
						auto rk = P_RK_X[p];
						rk->push(ders.P_dxdt[p]);
						p->setXSingle(rk->getX());
					}
				}
				std::cout << Green << "setBounds" << reset << std::endl;
				water.setBounds();
				std::cout << Green << "real_timeを取得" << reset << std::endl;
				real_time = P_RK_phi[*Points.begin()]->gett();
				std::cout << Green << "translation(real_time) = " << translation(real_time) << reset << std::endl;
				obj.translateFromInitialX(translation(real_time));
			} while (!(P_RK_phi[*Points.begin()]->finished));

			for (const auto &p : Points)
			{
				delete P_RK_phi[p];
				delete P_RK_X[p];
			}

			// b* ------------------------- 出力 ------------------------- */
			{

				auto hist = Histogram(extLength(water.getLines()));
				std::stringstream ss;
				ss << "\"cumulative_count\":" << hist.cumulative_count << ","
				   << "\"diff\":" << hist.diff << ","
				   << "\"count\":" << hist.count << ","
				   << "\"interval\":" << hist.interval << ","
				   << "\"bin_width\":" << hist.bin_width << ","
				   << "\"mid_interval\":" << hist.mid_interval << ",";
				std::string s = ss.str();
				std::replace(s.begin(), s.end(), '{', '[');
				std::replace(s.begin(), s.end(), '}', ']');
				std::cout << "{" << s << "}" << std::endl;
				if (*hist.data.rbegin() > 5)
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "length > 5");

				int ii = 0;
				derivatives ders(Points, dt);
				//-------------------------------------------
				map_P_Vd P_phiphin;
				// V_d lim_len = Subdivide(0.3, 0.2, 5 - 1);
				/* ------------------------------------------------------ */
				std::cout << "-------------------- 次の時刻の変数の値を得る -------------------- " << std::endl;
				std::cout << "------------------ getImprovedの後に微分を評価 ----------------------- " << std::endl;

				auto P_gradPhi = ders.P_gradPhi;
				auto P_gradPhi_tangential = ders.P_gradPhi_tangential;
				auto P_phin_vector = ders.P_phin_vector;
				auto P_dxdt = ders.P_dxdt;
				auto P_dxdt_mod = ders.P_dxdt_mod;
				auto P_DphiDt = ders.P_DphiDt;
				auto P_aphiat = ders.P_aphiat;
				auto P_pressure = ders.P_pressure;
				auto P_kappa = ders.P_kappa;
				auto P_laplacian = ders.P_laplacian;
				/* ------------------------------------------------------ */
				Print("圧力の計算", Magenta);
				uomap_P_Vd P_lines_length, P_Intxn_length;
				uomap_P_d P_state, P_solidangleBIE, P_height, P_phi, P_phin, P_face_size, P_radius, P_lines_size, P_ishit;
				uomap_P_d P_BC, P_Intxn_size, P_ContactFaces, P_is_multiple_phiphin;
				uomap_P_Tddd P_normal, P_normal_BEM, P_mirrorPosition, P_U_normal_BEM, P_U_tangential_BEM, P_U_BEM;
				Print("ders.P_phiphin_InnerOuterCornerPを出力");
				try
				{
					for (const auto &p : water.getPoints())
					{
						P_U_BEM[p] = p->U_BEM;
						P_solidangleBIE[p] = p->getSolidAngle();
						P_U_normal_BEM[p] = p->U_normal_BEM;
						P_U_tangential_BEM[p] = p->U_tangential_BEM;
						P_state[p] = p->getStatus();
						P_height[p] = p->getX()[2];
						P_phi[p] = std::get<0>(p->phiphin);
						P_phin[p] = std::get<1>(p->phiphin);
						P_normal[p] = p->getNormalTuple();
						P_normal_BEM[p] = p->normal_BEM;
						P_face_size[p] = (double)p->getFaces().size();
						P_lines_size[p] = (double)p->getLines().size();
						P_lines_length[p] = extLength(p->getLines());
						P_Intxn_size[p] = (double)takeIntxn(p->getLines()).size();
						P_ContactFaces[p] = (double)p->getContactFaces().size();
						P_Intxn_length[p] = extLength(takeIntxn(p->getLines()));
						P_BC[p] = p->Dirichlet ? 0. : (p->Neumann ? 1. : (p->CORNER ? 2. : 1 / 0.));
						if (!p->getContactFaces().empty())
							P_mirrorPosition[p] = (*p->getContactFaces().begin())->mirrorPosition(p) - p->getXtuple();
						P_radius[p] = p->radius;
						// P_ishit[p] = (double)(!p->getContactFaces().empty());
						P_ishit[p] = (double)(p->getStatus());
						P_is_multiple_phiphin[p] = (double)(p->multiple_phiphin.size());
					}
				}
				catch (std::exception &e)
				{
					std::cerr << e.what() << reset << std::endl;
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
				};
				try
				{
					// phiのアップデートがされていない
					VV_VarForOutput data = {
						{"U_BEM", P_U_BEM},
						{"U_normal_BEM", P_U_normal_BEM},
						{"U_tangential_BEM", P_U_tangential_BEM},
						{"kappa", P_kappa},
						{"ContactFaces", P_ContactFaces},
						// {"dxdt_mod", P_dxdt_mod},
						{"grad_phi", P_gradPhi},
						{"phin_vector", P_phin_vector},
						{"grad_phi_tangential", P_gradPhi_tangential},
						{"height", P_height},
						{"phi", P_phi},
						{"phin", P_phin},
						{"solidangle", P_solidangleBIE},
						{"normal", P_normal},
						{"normal_BEM", P_normal_BEM},
						{"laplacian_phi", P_laplacian},
						{"face_size", P_face_size},
						{"line_length", 10, P_lines_length},
						{"line_size", P_lines_size},
						{"Intxn_size", P_Intxn_size},
						{"Intxn_length", 10, P_Intxn_length},
						{"boundary condition", P_BC},
						{"radius", P_radius},
						{"state", P_state},
						{"tension", ders.P_tension},
						{"DphiDt", P_DphiDt},
						{"aphiat", P_aphiat},
						{"pressure", P_pressure},
						{"vector to mirrorPosition", P_mirrorPosition},
						{"is hit", P_ishit},
						{"is_multiple_phiphin", P_is_multiple_phiphin}};
					//流体
					{
						auto filename = water_output_vtu_file_name + std::to_string(time_step) + ".vtu";
						mk_vtu(water_output_directory + "/" + filename, water.getFaces(), data);
						waterPVD.push(filename, real_time);
						waterPVD.output();
					}
					//壁面
					{
						auto filename = tank_output_vtu_file_name + std::to_string(time_step) + ".vtu";
						mk_vtu(tank_output_directory + "/" + filename, obj.getFaces());
						tankPVD.push(filename, real_time);
						tankPVD.output();
					}
				}
				catch (std::exception &e)
				{
					std::cerr << e.what() << reset << std::endl;
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
				};
				// mk_vtu(output_directory + "/" + obj.getName() + std::to_string(time_step) + ".vtu", obj.getFaces(), datacpg);
				// mk_vtu(output_directory + "/" + name + std::to_string(time_step) + ".vtu", cpg.well_Faces, datacpg);
				// cpg_pvd.push(name, name + std::to_string(time_step) + ".vtu", time_step * t_rep * dt);
				// cpg_pvd.output_();
			}
		}
	}
	catch (error_message &e)
	{
		e.print();
	};
	return 0;
};
