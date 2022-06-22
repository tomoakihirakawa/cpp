#define BEM
#define use_lapack

int time_step;

double real_time = 0;

#define simulation

#include "GNUPLOT.hpp"
#include "Network.hpp"
#include "minMaxOfFunctions.hpp"
#include "integrationOfODE.hpp"
#include <filesystem>
#include "rootFinding.hpp"
#include "kernelFunctions.hpp"
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
// #define rotation_test
// #define up_down_test
// #define testCase1
// #define experiment_Chaplin1999_Retzler2000
#define experiment_Retzler2000
// #define experiment_Li2002
// #define experiment_Goring1979
// #define Hu2002
// #define WenWenLi2002_MachReflection
// #define experiment_sawai
// #define XueAndLin2011_large_amplitude

/* ------------------------------------------------------ */

namespace forced_motion
{
	double start = .1;
	double move_amplitude = .005;
	double s = M_PI / 2.;
	double a = move_amplitude;
	double k = M_PI / 1.;
	auto g = _GRAVITY_;
	T6d move_dir = {1., 0., 0., 0, 0, 0};
#if defined(experiment_sawai)
	auto h = 0.1;
	auto L = 0.25;
#elif defined(experiment_Li2002)
	double h = 0.3048;
	double H = 0.3 * h; //造波する初期入射波の波高
	double x = 0;
	double c = std::sqrt(g * (H + h));
#elif defined(experiment_Goring1979)
	double h = 0.25;
	double H = 0.1 * h; //造波する初期入射波の波高
	double x = 0;
	double c = std::sqrt(g * (H + h));
#elif defined(WenWenLi2002_MachReflection)
	double h = 0.06;
	double H = 0.5 * h;
	double x = 0;
	double c = std::sqrt(g * (H + h));
#elif defined(experiment_)
	double h = 0.3048;
	double H = 0.6 * h;
	double x = 0;
	double c = std::sqrt(g * (H + h));
#elif defined(Hu2002)
	double h = 0.5;
	double w = 1.257 / std::sqrt(h / g); /*5.57065*/
	double A = 0.046 * h;
#elif defined(experiment_Retzler2000)
	const std::vector<Tdd> sample = {
		{-0.15000000000000002, 0.},
		{-0.11, 0.},
		{-0.1, 0.},
		{-0.09, 0.},
		{-0.08, 0.},
		{-0.07, 0.},
		{-0.06, 0.},
		{-0.050409836065573754, 0.007920792079208039},
		{-0.02397540983606558, 0.06930693069306948},
		{-0.0006147540983606481, 0.13762376237623775},
		{0.0282786885245902, 0.24851485148514862},
		{0.05040983606557381, 0.3594059405940595},
		{0.05901639344262294, 0.40297029702970305},
		{0.06946721311475412, 0.4425742574257427},
		{0.08299180327868855, 0.4792079207920793},
		{0.10020491803278692, 0.516831683168317},
		{0.11065573770491804, 0.5415841584158416},
		{0.12110655737704923, 0.5663366336633664},
		{0.132172131147541, 0.5910891089108912},
		{0.15000000000000008, 0.6198019801980199},
		{0.1616803278688525, 0.6306930693069308},
		{0.17643442622950822, 0.6445544554455447},
		{0.18872950819672135, 0.6603960396039605},
		{0.20102459016393448, 0.6792079207920793},
		{0.21823770491803285, 0.6801980198019802},
		{0.23053278688524592, 0.6445544554455447},
		{0.25081967213114753, 0.5742574257425743},
		{0.27848360655737703, 0.40297029702970305},
		{0.30000000000000004, 0.23564356435643574},
		{0.319672131147541, 0.10396039603960416},
		{0.3319672131147542, 0.02772277227722786},
		{0.35040983606557385, -0.03960396039603942},
		{0.37069672131147546, -0.085148514851485},
		{0.38913934426229513, -0.08316831683168302},
		{0.4002049180327869, -0.06831683168316827},
		{0.4254098360655738, -0.03069306930693061},
		{0.45000000000000007, -0.009900990099009799},
		{0.5, 0.},
		{0.55, 0.},
		{0.6, 0.},
		{0.65, 0.},
		{0.7, 0.}};
	const auto intp = InterpolationBspline(3, sample);
#endif
	Tddd displacement(const double t)
	{
#ifdef XueAndLin2011_large_amplitude
		auto w = 3.5317;
		Tddd move_dir = {1., 0., 0.};
		move_amplitude = 0.1;
		if (t > start)
			return -move_amplitude * cos(w * (t - start)) * move_dir;
		else
			return {0., 0., 0.};
#elif defined(experiment_sawai)
		auto w = std::sqrt(M_PI * g / L * tanh(M_PI * h / L));
		Tddd move_dir = {1., 0., 0.};
		return move_amplitude * sin(w * t) * move_dir;
#elif defined(experiment_Li2002)
		//使わない
		Tddd move_dir = {1., 0., 0.};
		return move_dir;
#elif defined(Hu2002)
		if (t >= start)
		{
			auto tmp = -A * cos(w * (t - start));
			tmp -= -A;
			return {
				tmp,
				0.,
				0.,
			};
		}
		else
			return {0., 0., 0.};
#elif defined(experiment_Goring1979)
		double k = std::sqrt(3. * H / (4. * h * h * h));
		Tddd move_dir = {1., 0., 0.};
		// auto xi = atanh((tanh(c * k * (t - start))) / sqrt(1 + h)) / (sqrt(1 + h) * k);
		// auto xi0 = atanh((tanh(c * k * (-start))) / sqrt(1 + h)) / (sqrt(1 + h) * k);
		auto tmp = (sqrt(H) * atanh((sqrt(H) * tanh(c * k * (-start + t))) / sqrt(h + H))) / (sqrt(h + H) * k);
		tmp -= (sqrt(H) * atanh((sqrt(H) * tanh(c * k * (-start))) / sqrt(h + H))) / (sqrt(h + H) * k);
		return tmp * move_dir;
#elif defined(WenWenLi2002_MachReflection)
		Tddd move_dir = {1., 0., 0.};
		return move_dir;
#else
		/* ------------------------------------------------------ */
		Tddd move_dir = {cos(k * t), sin(k * t), 0.};
		return move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) * move_dir;
		/* ------------------------------------------------------ */
		// Tddd move_dir = Normalize(Tddd{1., 1., 0.});
		// return move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) * move_dir;
#endif
	};

	T6d velocity(double t)
	{
#ifdef XueAndLin2011_large_amplitude
		auto w = 3.5317;
		move_amplitude = 0.1;
		if (t > start)
			return move_amplitude * w * sin(w * (t - start)) * move_dir;
		else
			return {0., 0., 0., 0., 0., 0.};
#elif defined(experiment_sawai)
		auto w = std::sqrt(M_PI * g / L * tanh(M_PI * h / L));
		if (t > start)
			return -move_amplitude * w * sin(w * (t - start)) * move_dir;
		else
			return {0., 0., 0., 0., 0., 0.};
#elif defined(experiment_Li2002)
		double k = std::sqrt(3. * H / (4. * h * h * h));
		double tmp = cosh(k * (x - c * (t - start)));
		double eta = H * std::pow(1. / tmp, 2.);
		double u = c * eta / (h + eta);
		return u * move_dir;
#elif defined(experiment_Goring1979)
		double kappa = std::sqrt(3. * H / (4. * h * h * h));
		double eta = H * std::pow(1. / cosh(kappa * (-c * (t - start))), 2.);
		double u = c * eta / (h + eta);
		return u * move_dir;
#elif defined(rotation_test)
		double T = .25;
		// if (t > start)
		// 	return {0., 0., 0., M_PI / 8., 0., 0.};
		// else
		// 	return {0., 0., 0., 0., 0., 0.};

		if (t > start)
			return {0., 0., 0., 40. * M_PI / 180. * sin(M_PI / T * (t - start)), 0., 0.};
		else
			return {0., 0., 0., 0., 0., 0.};
#elif defined(up_down_test)
		double T = .25;
		if (t > start)
			return {0., 0., -0.07 * sin(M_PI / T * (t - start)), 0., 0., 0.};
		else
			return {0., 0., 0., 0., 0., 0.};
#elif defined(Hu2002)
		if (t >= start)
			return {A * w * sin(w * (t - start)), 0., 0., 0., 0., 0.};
		else
			return {0., 0., 0., 0., 0., 0.};
#elif defined(testCase1)
		double h = 0.5;
		double w = 1. / std::sqrt(h / g); /*5.57065*/
		double A = 0.1 * h;
		if (t >= start)
			return {A * w * sin(w * (t - start)),
					A * w * cos(w * (t - start)), 0., 0., 0., 0.};
		else
			return {0., 0., 0., 0., 0., 0.};
			// double T = .01;
			// double b = 0.001 * M_PI / T;
			// if (t > start)
			// 	return {b * sin(2 * M_PI / T * (t - start)),
			// 			b * cos(2 * M_PI / T * (t - start)),
			// 			0., 0., 0., 0.};
			// else
			// 	return {0., 0., 0., 0., 0., 0.};
#elif defined(experiment_Chaplin1999_Retzler2000)
		double T = .1;
		if (t > start)
			return {0.3 * sin(M_PI / T * (t - start)), 0., 0., 0., 0., 0.};
		else
			return {0., 0., 0., 0., 0., 0.};
#elif defined(experiment_Retzler2000)
		double T = .1;
		t -= start;
		return {intp(t), 0., 0., 0., 0., 0.};
		// if (t > 0.633918 || t < 0)
		// 	return {0., 0., 0., 0., 0., 0.};
		// else if (0.121637 < t && t <= 0.512281)
		// 	return {0.292167, 0., 0., 0., 0., 0.};
		// else if (0.512281 < t && t <= 0.633918)
		// 	return {1.52264 - 2.40195 * t, 0., 0., 0., 0., 0.};
		// else
		// 	return {2.40195 * t, 0., 0., 0., 0., 0.};
#elif defined(WenWenLi2002_MachReflection)
		double eta = H * std::pow(1. / cosh(std::sqrt(3. * H / (4. * std::pow(h, 3))) * (x - c * (t - start))), 2.);
		double u = c * eta / (h + eta);
		if (t > start)
			return u * move_dir;
		else
			return {0., 0., 0., 0., 0., 0.};
#else
		/* ------------------------------------------------------ */
		T6d move_dir = {cos(k * t), sin(k * t), 0., 0., 0., 0.};
		T6d ddt_move_dir = {-k * sin(k * t), k * cos(k * t), 0., 0., 0., 0.};
		// /* |U|*n_p . n_surface = phin <-- given
		auto tmp = (-move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) + move_amplitude * exp(-t) * (cos(k * t - s) * k)) * move_dir;
		tmp += move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) * ddt_move_dir;
		return tmp;
#endif
	};

	T6d acceleration(const double t)
	{
#ifdef XueAndLin2011_large_amplitude
		auto w = 3.5317;
		move_amplitude = 0.1;
		if (t > start)
			return move_amplitude * w * w * cos(w * (t - start)) * move_dir;
		else
			return {0., 0., 0., 0., 0., 0.};
#elif defined(experiment_sawai)
		auto w = std::sqrt(M_PI * g / L * tanh(M_PI * h / L));
		return -move_amplitude * w * w * cos(w * t) * move_dir;
#elif defined(experiment_Li2002)
		auto u = (2 * std::sqrt(3) * pow(c, 2) * h * H * std::sqrt(H / pow(h, 3)) * std::sinh(std::sqrt(3) * std::sqrt(H / pow(h, 3)) * (-(c * t) + x)));
		u /= pow(h + 2 * H + h * std::cosh(std::sqrt(3) * std::sqrt(H / pow(h, 3)) * (-(c * t) + x)), 2);
		return u * move_dir;
#elif defined(WenWenLi2002_MachReflection)
		auto u = (2 * std::sqrt(3) * pow(c, 2) * h * H * std::sqrt(H / pow(h, 3)) * std::sinh(std::sqrt(3) * std::sqrt(H / pow(h, 3)) * (-(c * (t - start)) + x)));
		u /= pow(h + 2 * H + h * std::cosh(std::sqrt(3) * std::sqrt(H / pow(h, 3)) * (-(c * (t - start)) + x)), 2);
		return u * move_dir;
#else
		return {0., 0., 0., 0., 0., 0.};
#endif
	};
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
			// auto angle = MyVectorAngle(f->getNormalTuple(), -F->getNormalTuple());
			// if ((angle / M_PI * 180. < 60.) || ((M_PI - angle) / M_PI * 180. < 60.))
			ret.emplace(F);
		}
	return ret;
};
/* ------------------------------------------------------ */
netFp closestFacingFace(const networkFace *const f_IN)
{
	//!もしない場合はnullptrを返すので注意
	networkFace *closest_face = nullptr;
	double min_distance = 1E+100;
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
std::unordered_map<networkFace *, networkFace *> getNearestContactFacesOfSurroundedNeumannFace(const networkPoint *const p)
{
	std::unordered_map<networkFace *, networkFace *> structure_face;
	networkFace *tmp;
	//@ 点の隣接面のうちNeumann面を抜き出し，それぞれに対してclosestFacingFaceを取得する．
	for (const auto &F : p->getFaces())
		if (F->Neumann && (tmp = closestFacingFace(F) /*nulltrなら入らない*/))
			structure_face[F] = tmp;
	return structure_face;
};

Tddd getNearestContactFacesX(const Tddd &X,
							 const double radius,
							 const std::vector<T3Tddd> &vertices)
{
	Tddd ret = {1E+100, 1E+100, 1E+100};
	for (const auto &vertex : vertices)
	{
		auto intxn = IntersectionSphereTriangle_(X, radius, vertex);
		if (intxn.isIntersecting)
		{
			auto r = intxn.X - X;
			auto nr = Norm(r);
			if (Norm(ret - X) > nr)
				ret = intxn.X;
		}
	}
	return ret;
};
std::tuple<networkFace *, Tddd> getNearestContactFacesX(const networkPoint *const p)
{
	std::tuple<networkFace *, Tddd> ret = {nullptr, {1E+100, 1E+100, 1E+100}};
	auto pX = p->getXtuple();
	if (!p->getContactFaces().empty())
		for (const auto &[f, X] : p->getContactFacesX())
		{
			if (isFinite(X - pX))
				if (Norm(std::get<1>(ret) - pX) > Norm(X - pX))
					ret = {f, X};
		}
	return ret;
};

Tddd uNeumann(const networkPoint *const p)
{
	auto tmpContactFaces = p->getContactFacesXCloser();
	if (tmpContactFaces.empty())
		return {0., 0., 0.};
	std::vector<Tddd> V;
	std::vector<double> W;
	Tddd v;
	double w;
	for (const auto &[f, X] : tmpContactFaces)
	{
		v = f->getNetwork()->velocityRigidBody(X);
		w = kernel_Bspline3(Norm(X - p->getXtuple()), p->radius);
		V.emplace_back(v);
		W.emplace_back(w);
	}
	auto [f, X] = tmpContactFaces[0];
	return optimumVector_(V, {0., 0., 0.}, W);
};

/* ------------------------------------------------------ */

// double phi_stokes2nd(double x, double z, const double t)
// {
// 	const double g = 9.8;
// 	auto tmp = a * w / k / sinh(k * h);
// 	tmp *= (cosh(k * (z + h)) * sin(q) + k * a * 3 * cosh(2 * k * (z + h)) / (8 * pow(sinh(k * h), 3)) * sin(2 * q));
// 	tmp -= pow(k * a, 2) / (2 * sinh(2 * k * h) * g * t / k);
// 	return tmp;
// };
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
double accel_normal_from_Neumann_surface(const networkPoint *const p)
{
	/*
	 * 面積平均に変更した
	 */
	//@ 点の隣接面のうちNeumann面を抜き出し，それぞれに対してclosestFacingFaceを取得する．
	//@ 次に，面それぞれの速度を計算し，点の法線方向速度を計算する．
	double ret = 0., A = 0, Atot = 0;
	int i = 0;
	for (auto &[sF, cF] : getNearestContactFacesOfSurroundedNeumannFace(p) /*この点に隣接する流体面それぞれが最も近くで接している構造物の面*/)
	{
		A = sF->getArea();
		Atot += A;
		ret += A * Dot(cF->getNetwork()->acceleration, ToT6d(sF->getNormalTuple()));
		i++;
	}
	return i == 0 ? ret : ret / Atot;
}
/* ------------------------------------------------------ */
double velocity_normal_from_Neumann_surface(const networkPoint *const p)
{
	/*
	 * 面積平均に変更した
	 */
	//@ 点の隣接面のうちNeumann面を抜き出し，それぞれに対してclosestFacingFaceを取得する．
	//@ 次に，面それぞれの速度を計算し，点の法線方向速度を計算する．
	double ret = 0., A = 0, Atot = 0;
	int i = 0;
	for (auto &[sF, cF] : getNearestContactFacesOfSurroundedNeumannFace(p) /*この点に隣接する流体面それぞれが最も近くで接している構造物の面*/)
	{
		A = sF->getArea();
		Atot += A;
		ret += A * Dot(cF->getNetwork()->velocity, ToT6d(sF->getNormalTuple()));
		i++;
	}
	return i == 0 ? ret : ret / Atot;
}
/* ------------------------------------------------------ */
T6d velocity_from_Neumann_surface(const networkPoint *const p)
{
	T6d ret = {0., 0., 0., 0., 0., 0.};
	double A = 0, Atot = 0;
	int i = 0;
	for (auto &[sF, cF] : getNearestContactFacesOfSurroundedNeumannFace(p) /*この点に隣接する流体面それぞれが最も近くで接している構造物の面*/)
	{
		auto n_cF = cF->getNormalTuple();
		auto n_sF = sF->getNormalTuple();
		A = sF->getArea();
		Atot += A;
		Tddd V = cF->getNetwork()->velocityRigidBody(p->getXtuple());
		ret += A * Dot(V, p->getNormalTuple());
		i++;
	}
	return i == 0 ? ret : ret / Atot;
}
Tddd local_velocity_from_Neumann_surface(const networkPoint *const p)
{
	Tddd ret = {0., 0., 0.};
	double A = 0, Atot = 0;
	int i = 0;
	for (auto &[sF, cF] : getNearestContactFacesOfSurroundedNeumannFace(p) /*この点に隣接する流体面それぞれが最も近くで接している構造物の面*/)
	{
		auto n_cF = cF->getNormalTuple();
		auto n_sF = sF->getNormalTuple();
		A = sF->getArea();
		Atot += A;
		Tddd V = cF->getNetwork()->velocityRigidBody(p->getXtuple());
		ret += A * V;
		i++;
	}
	return i == 0 ? ret : ret / Atot;
}
/* ------------------------------------------------------ */
double phin_contact(const netFp f_IN)
{
	auto closest_face = closestFacingFace(f_IN);
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
// using map_P_TiiiTdd = std::map<netP *, std::map<Tiii, Tdd>>;
using map_P_VVd = std::map<netP *, VV_d>;
using map_F_P_Vd = std::map<netF *, map_P_Vd>;
using map_P_P_Vd = std::map<netP *, map_P_Vd>;
using pair_PB = std::pair<netP *, bool>;
using map_pairPB_Tdd = std::map<pair_PB, Tdd>;
using map_pairPB_pairPB_Tdd = std::map<pair_PB /*タプル*/, std::map<pair_PB, Tdd>>;
using map_P_P_Tdd = std::map<netP *, std::map<netP *, Tdd>>;
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
	if (isExsist("center_of_mass"))
	{
		std::get<0>(water.center_of_mass) = stob(json1()["center_of_mass"])[0];
		std::get<1>(water.center_of_mass) = stob(json1()["center_of_mass"])[1];
		std::get<2>(water.center_of_mass) = stob(json1()["center_of_mass"])[2];
	}
	if (isExsist("ignore"))
	{
		water.IGNORE = stob(json1()["ignore"])[0];
	}
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
		else
			water.scale(scale[0]);
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

//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
Tddd gradTangential_LinearElement(const Tddd &V, const T3Tddd &X012)
{
	auto [X0, X1, X2] = X012;
	return Cross(TriangleNormal(X012),
				 std::get<0>(V) * (X2 - X1) +
					 std::get<1>(V) * (X0 - X2) +
					 std::get<2>(V) * (X1 - X0)) /
		   (2 * TriangleArea(X012));
};

// Tddd gradTangential_LinearElement(const Tddd &V, const networkFace *const f)
// {
// 	return gradTangential_LinearElement(V, f->getXVertices());
// 	// auto [X0, X1, X2] = f->getXVertices();
// 	// return Cross(f->getNormalTuple(),
// 	// 			 std::get<0>(V) * (X2 - X1) +
// 	// 				 std::get<1>(V) * (X0 - X2) +
// 	// 				 std::get<2>(V) * (X1 - X0)) /
// 	// 	   (2 * f->getArea());
// };

// using T3T2ddd = std::tuple<T2Tddd, T2Tddd, T2Tddd>;

// Tddd gradTangential_LinearElement(const T3T2ddd &V)
// {
// 	auto [XV0, XV1, XV2] = V;
// 	return gradTangential_LinearElement({std::get<0>(XV0), std::get<0>(XV1), std::get<0>(XV2)},
// 										{std::get<1>(XV0), std::get<1>(XV1), std::get<1>(XV2)});
// 	// auto [X0, X1, X2] = f->getXVertices();
// 	// return Cross(f->getNormalTuple(),
// 	// 			 std::get<0>(V) * (X2 - X1) +
// 	// 				 std::get<1>(V) * (X0 - X2) +
// 	// 				 std::get<2>(V) * (X1 - X0)) /
// 	// 	   (2 * f->getArea());
// };

Tddd gradPhiTangential(const networkFace *const f)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
										{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
}

Tddd gradPhiTangential(const networkFace *const f, const networkPoint *p)
{
	auto [p0, p1, p2] = f->getPointsTuple(p);
	return gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
										{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
};

Tddd gradPhi(const networkFace *const f)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	Tddd ret = gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
											{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
	if (f->Neumann)
		ret += f->getNormalTuple() * (p0->phin_Neumann + p1->phin_Neumann + p2->phin_Neumann) / 3.;
	else if (f->Dirichlet)
		ret += f->getNormalTuple() * (p0->phin_Dirichlet + p1->phin_Dirichlet + p2->phin_Dirichlet) / 3.;
	else
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	return ret;
}

Tddd gradPhi(const networkFace *const f,
			 const networkPoint *const p,
			 const double t0 = 1. /*phinにのみ関するもの*/)
{
	Tddd N = {t0, (1 - t0) / 2., (1 - t0) / 2.};
	auto [p0, p1, p2] = f->getPointsTuple(p);
	Tddd ret = gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
											{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
	if (f->Neumann)
		ret += f->getNormalTuple() * Dot(Tddd{p0->phin_Neumann, p1->phin_Neumann, p2->phin_Neumann}, N);
	else if (f->Dirichlet)
		ret += f->getNormalTuple() * Dot(Tddd{p0->phin_Dirichlet, p1->phin_Dirichlet, p2->phin_Dirichlet}, N);
	else
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	return ret;
}

Tddd gradPhiTangential(const networkFace *const f,
					   const networkPoint *const p,
					   const double t0 = 1. /*phinにのみ関するもの*/)
{
	Tddd N = {t0, (1 - t0) / 2., (1 - t0) / 2.};
	auto [p0, p1, p2] = f->getPointsTuple(p);
	Tddd ret = gradTangential_LinearElement({std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)},
											{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
	return ret;
}

double getPhi(const networkLine *const l)
{
	auto fs = l->getFaces();
	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_0(fs[0], l);
	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_1(fs[1], l);
	return Dot(intp_l0_0.N(.5, .5), ToPhi(intp_l0_0.Points)) * 0.5 +
		   Dot(intp_l0_1.N(.5, .5), ToPhi(intp_l0_1.Points)) * 0.5;
}

T3Tddd grad_U_tangential_LinearElement(const networkFace *const f)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	/*
		{grad(Ux),
		grad(Uy),
		grad(Uz)}
	*/
	return {gradTangential_LinearElement({std::get<0>(p0->U_BEM), std::get<0>(p1->U_BEM), std::get<0>(p2->U_BEM)}, {p0->getXtuple(), p1->getXtuple(), p2->getXtuple()}),
			gradTangential_LinearElement({std::get<1>(p0->U_BEM), std::get<1>(p1->U_BEM), std::get<1>(p2->U_BEM)}, {p0->getXtuple(), p1->getXtuple(), p2->getXtuple()}),
			gradTangential_LinearElement({std::get<2>(p0->U_BEM), std::get<2>(p1->U_BEM), std::get<2>(p2->U_BEM)}, {p0->getXtuple(), p1->getXtuple(), p2->getXtuple()})};
};

// T3Tddd grad_U_LinearElement2(const networkPoint *const p)
// {
// 	T3Tddd grad_U_tangential = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
// 	double Atot = 0;
// 	for (const auto &f : p->getFaces())
// 	{
// 		grad_U_tangential += f->getArea() * grad_U_tangential_LinearElement(f);
// 		Atot += f->getArea();
// 	}
// 	/*
// 		grad_U_tangential = {grad(Ux),grad(Uy),grad(Uz)}
// 	*/
// 	grad_U_tangential /= Atot;
// 	Tddd x = {1, 0, 0}, y = {0, 1, 0};
// 	auto V00 = Dot(std::get<0>(tmp), x);
// 	auto V10 = Dot(std::get<1>(tmp), x);
// 	auto V20 = Dot(std::get<2>(tmp), x);
// 	auto V01 = Dot(std::get<0>(tmp), y);
// 	auto V11 = Dot(std::get<1>(tmp), y);
// 	auto V21 = Dot(std::get<2>(tmp), y);
// 	return T3Tddd{Tddd{V00, V01, V20},
// 				  Tddd{V10, V11, V21},
// 				  Tddd{V20, V21, -V00 - V11}};
// };

T3Tddd grad_U_LinearElement(const networkPoint *const p)
{
	T3Tddd gradUtang = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	/*
	スカラー量の接線方向勾配を計算することはできるが，法線方向はわからない．
	しかし，連続の式を使えば，phiの法線方向の勾配は，接線方向の勾配から計算することができる．
	*/
	double Atot = 0;
	for (const auto &f : p->getFaces())
	{
		gradUtang += f->getArea() * grad_U_tangential_LinearElement(f);
		Atot += f->getArea();
	}
	gradUtang /= Atot; //接線方向勾配
	auto q = (p->getNeighbors())[0];
	auto V = q->getXtuple() - p->getXtuple();
	auto n = p->getNormal_BEM();
	Tddd s0 = Normalize(V - n * Dot(V, n));
	Tddd s1 = Normalize(Cross(n, s0));
	Tddd V0 = Dot(gradUtang, s0);
	Tddd V1 = Dot(gradUtang, s1);
	Tddd V2 = {std::get<2>(V0), std::get<2>(V1), -std::get<0>(V0) - std::get<1>(V1)};
	return T3Tddd{V0, V1, V2};
};
/* ------------------------------------------------------ */

T4Tddd gradPhi(const networkPoint *const p)
{
	const double t0 = 1;
	// const double t0 = 2 / 3.;
	// const double t0 = 1 / 3.;
	// const double t0 = 11. / 18.;
	// const double t0 = 5/6;
	// const double t0 = 1 - 1 / 3;
	Tddd grad, grad_normal, grad_tangential, n = p->getNormal_BEM();
	grad_normal = std::get<1>(p->phiphin) * n;
	V_Tddd U, Us, V, N;
	V_d weights;
	double w;
	for (const auto &f : p->getFaces())
	{
#if defined(Zhang2001_InverseAreaWeighted)
		w = 1 / f->getArea();
		U.emplace_back(gradPhi(f, p, t0) * w);
		weights.emplace_back(w);
#elif defined(use_inscribedcircle_area_weigted_normal)
		w = f->getInscribedCircleArea();
		U.emplace_back(gradPhi(f, p, t0) * w);
		weights.emplace_back(w);
#elif defined(use_angle_weigted_normal)
		w = std::get<0>(f->getAngles(p));
		U.emplace_back(gradPhi(f, p, t0) * w);
		weights.emplace_back(w);
#elif defined(use_area_weigted_normal)
		// if (p->CORNER)
		// {
		// 	if (f->Dirichlet)
		// 	{
		// 		w = f->getArea();
		// 		U.emplace_back(gradPhi(f, p, t0) * w);
		// 		weights.emplace_back(w);
		// 	}
		// }
		// else
		{
			w = f->getArea();
			U.emplace_back(gradPhi(f, p, t0) * w);
			weights.emplace_back(w);
		}
#elif defined(use_subarea_weigted_normal)
		auto [p0, p1, p2] = f->getPointsTuple(p);
		double phin0, phin1, phin2;
		if (f->Neumann)
		{
			phin0 = p0->phin_Neumann;
			phin1 = p1->phin_Neumann;
			phin2 = p2->phin_Neumann;
		}
		else if (f->Dirichlet)
		{
			phin0 = p0->phin_Dirichlet;
			phin1 = p1->phin_Dirichlet;
			phin2 = p2->phin_Dirichlet;
		}
		double A = TriangleArea(p0->getXtuple(), (p0->getXtuple() + p1->getXtuple()) / 2, (p0->getXtuple() + p1->getXtuple() + p2->getXtuple()) / 3);
		Tddd u = (phin0 + (phin0 + phin1) / 2 + (phin0 + phin1 + phin2) / 3) / 6 * f->getNormalTuple() + gradPhiTangential(f);
		U.emplace_back(u * A);
		weights.emplace_back(A);
		A = TriangleArea(p0->getXtuple(), (p0->getXtuple() + p2->getXtuple()) / 2, (p0->getXtuple() + p1->getXtuple() + p2->getXtuple()) / 3);
		u = (phin0 + (phin0 + phin2) / 2 + (phin0 + phin1 + phin2) / 3) / 6 * f->getNormalTuple() + gradPhiTangential(f);
		U.emplace_back(u * A);
		weights.emplace_back(A);
#elif defined(use_arithmetic_averaged_normal)
		w = 1.;
		U.emplace_back(gradPhi(f, p, t0) * w);
		weights.emplace_back(w);
#endif
	}

	grad = Total(U) / Total(weights);
	grad_tangential = grad - grad_normal;

	/* ------------------------------------------------------ */
	// if (p->CORNER)
	// {
	// 	auto n_D = p->getNormalDirichlet_BEM();
	// 	auto n_N = p->getNormalNeumann_BEM();
	// 	auto S = Normalize(Cross(n_D, n_N));
	// 	auto m = Dot(n_D, n_N);
	// 	grad = Dot(Total(Us) / Total(weights), S) * S;
	// 	if (std::abs(m * m - 1) < 1E-10)
	// 		grad += (p->phin_Dirichlet * n_D + p->phin_Neumann * n_N) / 2.;
	// 	else
	// 	{
	// 		grad += p->phin_Dirichlet * n_D;
	// 		grad += (p->phin_Neumann - m * p->phin_Dirichlet) * (m * n_D - n_N) / ((m - 1.) * (m + 1.));
	// 	}
	// 	if (!isFinite(grad))
	// 		grad = {0., 0., 0.};

	// 	grad_normal = Dot(grad, n) * n;
	// 	grad_tangential = grad - grad_normal;
	// }
	/* ------------------------------------------------------ */

	return {grad_tangential, grad_normal, grad, n};
};

/* ------------------------------------------------------ */

Tddd gradPhiWithoutTangential(const networkFace *const f)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	Tddd ret = {0., 0., 0.};
	if (f->Neumann)
		ret += f->getNormalTuple() * (p0->phin_Neumann + p1->phin_Neumann + p2->phin_Neumann) / 3.;
	else if (f->Dirichlet)
		ret += f->getNormalTuple() * (p0->phin_Dirichlet + p1->phin_Dirichlet + p2->phin_Dirichlet) / 3.;
	else
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	return ret;
}

Tddd gradPhiWithoutTangential(networkPoint *const p)
{
	V_Tddd U;
	for (const auto &f : p->getFaces())
		U.emplace_back(gradPhiWithoutTangential(f));
	return optimumVector_(U, std::get<2>(gradPhi(p)));
}

// b@ ------------------------------------------------------ */

Tddd nextX_U(const networkPoint *const p, const double dt)
{
	return p->getXBuffer();
};
std::vector<Tddd> nextX_U(const std::vector<networkPoint *> &ps, const double dt)
{
	std::vector<Tddd> ret;
	for (const auto &p : ps)
		ret.emplace_back(nextX_U(p, dt));
	return ret;
};
Tddd nextX_U_Ua(const networkPoint *const p, const double dt)
{
	return p->getXBuffer() + p->U_BUFFER * dt;
};
std::vector<Tddd> nextX_U_Ua(const std::vector<networkPoint *> &ps, const double dt)
{
	std::vector<Tddd> ret;
	for (const auto &p : ps)
		ret.emplace_back(nextX_U_Ua(p, dt));
	return ret;
};
std::tuple<Tddd, Tddd> nextX_U_Ua(const std::tuple<networkPoint *, networkPoint *> ps, const double dt)
{
	return {std::get<0>(ps)->getXBuffer() + std::get<0>(ps)->U_BUFFER * dt,
			std::get<1>(ps)->getXBuffer() + std::get<1>(ps)->U_BUFFER * dt};
};
Tddd next_normal_U(const networkFace *const f, const double dt)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return TriangleNormal(f->getXVertices() + dt * T3Tddd{p0->getXBuffer(), p1->getXBuffer(), p2->getXBuffer()});
};

Tddd next_area_weighted_normal(const networkPoint *p, const double dt)
{
	T3Tddd VVV;
	Tddd ret = {0., 0., 0.};
	for (const auto &f : p->getFaces())
	{
		auto [p0, p1, p2] = f->getPointsTuple();
		VVV = f->getXVertices() + dt * T3Tddd{p0->getXBuffer(), p1->getXBuffer(), p2->getXBuffer()};
		ret += TriangleNormal(VVV) * TriangleArea(VVV);
	}
	return ret;
};

Tddd next_normal_U_Ua(const networkFace *const f, const double dt)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return TriangleNormal(f->getXVertices() + dt * T3Tddd{nextX_U_Ua(p0, dt), nextX_U_Ua(p1, dt), nextX_U_Ua(p2, dt)});
};
/* ------------------------------------------------------ */
bool isFacing(const Tddd &F, const Tddd &f, double rad = 1E-10)
{
	return (Dot(F, -f) >= cos(rad));
};
Tddd closestXFacing(const Tddd &p_next_X, const double radius, const std::vector<T3Tddd> &vertices, const Tddd &n)
{
	Tddd r = {1E+100, 1E+100, 1E+100};
	for (const auto &vertex : vertices)
	{
		if (isFacing(TriangleNormal(vertex), n, M_PI / 180 * 20))
		{
			auto intxn = IntersectionSphereTriangle_(p_next_X, radius, vertex);
			if (intxn.isIntersecting)
				if (Norm(r) >= Norm(intxn.X - p_next_X))
					r = intxn.X - p_next_X;
		}
	}
	return r;
};
Tddd closestX(const Tddd &p_next_X, const double radius, const std::vector<T3Tddd> &vertices)
{
	Tddd r = {1E+100, 1E+100, 1E+100};
	for (const auto &vertex : vertices)
	{
		auto intxn = IntersectionSphereTriangle_(p_next_X, radius, vertex);
		if (intxn.isIntersecting)
			if (Norm(r) >= Norm(intxn.X - p_next_X))
				r = intxn.X - p_next_X;
	}
	return r;
};

Tddd nextFaceNormal_U_Ua(const networkFace *const f, const double dt)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return TriangleNormal(nextX_U_Ua(p0, dt),
						  nextX_U_Ua(p1, dt),
						  nextX_U_Ua(p2, dt));
};
Tddd nextX_U_Ua(const networkFace *const f, const double dt)
{
	// auto [p0, p1, p2] = f->getPointsTuple();
	return (nextX_U_Ua(std::get<0>(f->getPointsTuple()), dt) +
			nextX_U_Ua(std::get<1>(f->getPointsTuple()), dt) +
			nextX_U_Ua(std::get<2>(f->getPointsTuple()), dt)) /
		   3.;
};
/* ------------------------------------------------------ */
Tddd nextX_U_Ua_Uc(const networkPoint *const p, const double dt)
{
	return p->getXBuffer() + (p->U_BUFFER + p->U_cling_to_Neumann) * dt;
};

Tddd nextFaceNormal_U_Ua_Uc(const networkFace *const f, const double dt)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return TriangleNormal(nextX_U_Ua_Uc(p0, dt),
						  nextX_U_Ua_Uc(p1, dt),
						  nextX_U_Ua_Uc(p2, dt));
};
Tddd nextX_U_Ua_Uc(const networkFace *const f, const double dt)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return (nextX_U_Ua_Uc(p0, dt) + nextX_U_Ua_Uc(p1, dt) + nextX_U_Ua_Uc(p2, dt)) / 3.;
};
/* ------------------------------------------------------ */
std::vector<Tddd> vectorsClingToNeumann_old(const networkPoint *p, const double dt)
{
	if (p->Neumann || p->CORNER)
	{
		//接触面候補を多く取得しておく
		std::unordered_set<networkFace *> Fs;
		for (auto &[f, X] : p->getContactFacesX())
		{
			auto [p0, p1, p2] = f->getPointsTuple();
			for (auto &F : p0->getFaces())
				Fs.emplace(F);
			for (auto &F : p1->getFaces())
				Fs.emplace(F);
			for (auto &F : p2->getFaces())
				Fs.emplace(F);
		}

		//接触面候補の次の時刻の位置を予測
		std::vector<T3Tddd> next_Vrtx(Fs.size());
		int i = 0;
		for (auto &f : Fs)
		{
			auto [p0, p1, p2] = f->getPointsTuple();
			auto X0 = p0->getXtuple() + f->getNetwork()->velocityRigidBody(p0->getXtuple()) * dt;
			auto X1 = p1->getXtuple() + f->getNetwork()->velocityRigidBody(p1->getXtuple()) * dt;
			auto X2 = p2->getXtuple() + f->getNetwork()->velocityRigidBody(p2->getXtuple()) * dt;
			next_Vrtx[i++] = {X0, X1, X2};
		}

		if (next_Vrtx.empty())
			return {{0., 0., 0.}};

		//各面の次の時刻の配置を予測
		//最短距離でこの面と向き合う面をnext_Vrtxから探しF_clingsに保存
		auto Xp = nextX_U_Ua_Uc(p, dt);
		std::vector<Tddd> F_clings;
		for (const auto &f : p->getFacesNeumann())
		{
			auto n = nextFaceNormal_U_Ua_Uc(f, dt);
			auto Xf = nextX_U_Ua_Uc(f, dt);
			auto V = Xf - Xp;
			auto X = Xp;
			auto to_closest_X = closestXFacing(X, p->radius, next_Vrtx, n);
			if (isFinite(to_closest_X))
			{
				auto r = Dot(to_closest_X, n) * n;
				if (isFinite(r))
					F_clings.push_back(r);
			}
		}

		Tddd r = Mean(F_clings);
		if (isFinite(r))
			return F_clings;
		else
			return {{0., 0., 0.}};
	}
	else
		return {{0., 0., 0.}};
};
//$ ------------------------------------------------------ */
//$ ------------------------------------------------------ */
//$ ------------------------------------------------------ */
double next_length(const networkLine *const l, const double dt)
{
	auto [p0, p1] = l->getPointsTuple();
	return Norm(nextX_U_Ua(p0, dt) - nextX_U_Ua(p1, dt));
};
V_d next_length(const std::unordered_set<networkLine *> &ls,
				const double dt)
{
	V_d ret;
	for (const auto &l : ls)
		ret.emplace_back(next_length(l, dt));
	return ret;
};

double next_max_length(const networkPoint *const p, const double dt)
{
	auto p_next_X = nextX_U_Ua(p, dt);
	int size = 0;
	double max_len = -1;
	for (const auto &q : p->getNeighbors())
	{
		auto tmp = Norm(nextX_U_Ua(q, dt) - p_next_X);
		if (max_len < 0 || max_len < tmp)
			max_len = Norm(nextX_U_Ua(q, dt) - p_next_X);
	}
	return max_len;
};

double next_min_length(const networkPoint *const p, const double dt)
{
	auto p_next_X = nextX_U_Ua(p, dt);
	int size = 0;
	double min_len = -1, tmp;
	for (const auto &q : p->getNeighbors())
	{
		tmp = Norm(nextX_U_Ua(q, dt) - p_next_X);
		if (min_len < 0 || min_len > tmp)
			min_len = Norm(nextX_U_Ua(q, dt) - p_next_X);
	}
	return min_len;
};

double next_mean_length(const networkPoint *const p, const double dt)
{
	auto p_next_X = nextX_U_Ua(p, dt);
	int size = 0;
	double mean_len = 0.;
	for (const auto &q : p->getNeighbors())
	{
		mean_len += Norm(nextX_U_Ua(q, dt) - p_next_X);
		size++;
	}
	return mean_len / (double)size;
};

V_d next_length_shorter(const networkPoint *const p, const double dt)
{
	V_d ret(p->getLines().size());
	int i = 0;
	for (const auto &l : p->getLines())
		ret[i++] = next_length(l, dt);
	std::sort(ret.begin(), ret.end(), [](const auto &a, const auto &b)
			  { return a < b; });
	return ret;
};

std::vector<networkLine *> next_lines_shorter(const networkPoint *const p, const double dt)
{
	std::vector<networkLine *> ret = p->getLines();
	std::sort(ret.begin(), ret.end(), [dt](const auto &a, const auto &b)
			  { return next_length(a, dt) < next_length(b, dt); });
	return ret;
};

double next_mean_length_except(const networkPoint *const p, const double dt, networkLine *const L)
{
	auto p_next_X = nextX_U_Ua(p, dt);
	int size = 0;
	double mean_len = 0.;
	for (const auto &l : p->getLines())
	{
		if (L != l)
		{
			mean_len += Norm(nextX_U_Ua((*l)(p), dt) - p_next_X);
			size++;
		}
	}
	return mean_len / (double)size;
};

T3Tddd next_vertices(const networkFace *const f, const double dt)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return f->getXVertices() + dt * T3Tddd{nextX_U_Ua(p0, dt), nextX_U_Ua(p1, dt), nextX_U_Ua(p2, dt)};
};

double next_area(const networkFace *const f, const double dt)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return TriangleArea(f->getXVertices() + dt * T3Tddd{nextX_U_Ua(p0, dt), nextX_U_Ua(p1, dt), nextX_U_Ua(p2, dt)});
};

Tddd next_face_center(const networkFace *const f, const double dt)
{
	return Mean(next_vertices(f, dt));
};

double next_min_length_to_face(const networkPoint *const p, const double dt)
{
	auto p_next_X = nextX_U_Ua(p, dt);
	int size = 0;
	double min_len = -1, tmp;
	for (const auto &f : p->getFaces())
	{
		tmp = Norm(next_face_center(f, dt) - p_next_X);
		if (min_len < 0 || min_len > tmp)
			min_len = tmp;
	}
	return min_len;
};

double next_neighbors_mean_Area(const networkPoint *const p, const double dt)
{
	double ret = 0.;
	int i = 0;
	for (const auto &F : p->getFaces())
	{
		ret += next_area(F, dt);
		i++;
	}
	return ret / i;
};

double next_neighbors_mean_Area(const networkFace *const f, const double dt)
{
	double ret = 0.;
	int i = 0;
	for (const auto &F : f->getNeighbors())
	{
		ret += next_area(F, dt);
		i++;
	}
	return ret / i;
};

Tddd next_neighbors_center_AreaAveraged(const networkPoint *const p, const double dt)
{
	int i = 0;
	Tddd ret = {0., 0., 0.};
	if (p->CORNER)
	{
		for (const auto &l : p->getLines())
		{
			if (l->CORNER)
			{
				ret += nextX_U_Ua(((*l)(p)), dt);
				i++;
			}
		}
		return ret / (double)i;
	}
	else
	{
		double A = 0;
		for (const auto &f : p->getFaces())
		{
			double a = next_area(f, dt);
			A += a;
			ret += next_face_center(f, dt) * a;
			i++;
		}
		return ret / A;
	}
};

Tddd next_neighbors_center(const networkPoint *const p, const double dt)
{
	double s = 0;
	Tddd ret = {0., 0., 0.};
	for (const auto &l : p->getLines())
	{
		if (p->CORNER && l->CORNER)
		{
			ret += nextX_U_Ua(((*l)(p)), dt);
			s += 1.;
		}
		if (!p->CORNER)
		{
			for (const auto &l : p->getLines())
			{
				ret += nextX_U_Ua(((*l)(p)), dt);
				s += 1.;
			}
		}
	}
	return ret / s;
};

int getMinDepth(const networkPoint *const p)
{
	int min_depth = 999999;
	for (const auto &[net, depth] : p->net_depth)
		if (min_depth >= depth)
			min_depth = depth;
	return min_depth;
};

std::tuple<Network *, int> getMinDepthAndNetwork(const networkPoint *const p)
{
	std::tuple<Network *, int> ret = {nullptr, 999999};
	for (const auto &[net, depth] : p->net_depth)
		if (std::get<1>(ret) >= depth)
		{
			std::get<0>(ret) = net;
			std::get<1>(ret) = depth;
		}
	return ret;
};

bool isWellConditon(const networkPoint *p)
{
	auto [minA, maxA] = MinMax(extAreas(p->getFaces()));
	auto [minL, maxL] = MinMax(extLength(p->getLines()));
	return (minL / maxL > 0.8 /*差は0.2まで*/ || minA / maxA > 0.8 /*差は0.2まで*/);
};

Tddd next_neighbors_center_modified_for_corner(const networkPoint *const p, const double dt)
{
	bool corner_found = false, too_large = false;
	int i = 0, min_depth_from_all = 999999;
	Network *closest_net = nullptr;
	V_i min_depths(p->getNeighbors().size(), 999999);
	/* ------------------------------------------------------ */
	for (const auto &q : p->getNeighbors())
	{
		for (const auto &[net, depth] : q->net_depth)
			if (min_depth_from_all >= depth)
			{
				if (std::get<1>(net->grid_pull_factor) > 1E-10)
				{
					min_depth_from_all = depth;
					closest_net = net;
				}
			}
	}
	i = 0;
	/* ------------------------------------------------------ */
	Tddd ret = {0., 0., 0.};
	int max_limit = 20;
	int my_depth = 99999;
	if (p->net_depth.find(closest_net) != p->net_depth.end())
	{
		my_depth = p->net_depth.at(closest_net);
	}
	// else
	// {
	// 	i = 0;
	// 	for (const auto &q : p->getNeighbors())
	// 	{
	// 		ret += nextX_U_Ua(q, dt);
	// 		i++;
	// 	}
	// 	return ret / (double)i;
	// }

	// if (my_depth <= max_limit && isWellConditon(p))
	// {
	double I = 0;
	double mean_lines_size = 0;
	bool can_use_closest_net = true;
	i = 0;
	for (const auto &q : p->getNeighbors())
	{
		mean_lines_size += (double)q->getLines().size();
		if (q->net_depth.find(closest_net) == q->net_depth.end())
			can_use_closest_net = false;
		i++;
	}
	mean_lines_size /= i;
	double w;
	for (const auto &q : p->getNeighbors())
	{
		// if (q->CORNER)
		// {
		// 	if (std::get<1>(closest_net->grid_pull_factor) > 1E-10)
		// 	{
		// 		ret += nextX_U_Ua(q, dt);
		// 		I += 1.;
		// 	}
		// }

		w = 1; //(double)q->getLines().size() / mean_lines_size;
		if (can_use_closest_net && closest_net && isWellConditon(p))
		{
			auto c = 0.1 * (p->net_depth.at(closest_net) - q->net_depth.at(closest_net));
			w += c;
		}

		if (p->CORNER)
		{
			if (q->CORNER)
			{
				I += w;
				ret += nextX_U_Ua(q, dt) * w;
			}
		}
		else
		{
			I += w;
			ret += nextX_U_Ua(q, dt) * w;
		}
	}
	return ret / I;
};

// b@ ------------------------------------------------------ */
// b@ ------------------------------------------------------ */

Tddd fitToNeumannVelocity(Tddd VECTOR, const networkPoint *const p)
{
	if (p->Neumann || p->CORNER)
	{
		Tddd nf, vn0, vn1;
		double max_w = kernel_Bspline3(0., p->radius);
		for (const auto &[f, hit_X] : Reverse(p->getContactFacesXCloser()) /*遠い方から*/)
		{
			nf = f->getNormalTuple();
			vn0 = Dot(VECTOR, nf) * nf;
			vn1 = Dot(f->getNetwork()->velocityRigidBody(hit_X), nf) * nf;
			VECTOR += kernel_Bspline3(Norm(hit_X - p->getXtuple()), p->radius) / max_w * (vn1 - vn0);
		}
	}
	return VECTOR;
};

Tddd condition_Ua(Tddd VECTOR, const networkPoint *const p)
{
	if (p->CORNER)
	{
		auto cross = Normalize(Cross(p->getNormalNeumann_BEM(), p->getNormalDirichlet_BEM()));
		VECTOR = Dot(VECTOR, cross) * cross;
		// for (const auto &l : p->getLines())
		// 	if (l->CORNER)
		// 	{
		// 		Tddd dir = Normalize((*l)(p)->getXBuffer() - p->getXBuffer());
		// 		VECTOR = Dot(VECTOR, dir) * dir;
		// 	}
		for (const auto &f : p->getFaces())
			if (f->Neumann)
				VECTOR -= Dot(VECTOR, f->getNormalTuple()) * f->getNormalTuple();
		return VECTOR;
	}
	else if (p->Dirichlet)
	{
		return VECTOR - Dot(VECTOR, p->getNormal_BEM()) * p->getNormal_BEM();
	}
	else
	{
		for (const auto &f : p->getFaces())
			VECTOR -= Dot(VECTOR, f->getNormalTuple()) * f->getNormalTuple();
		return VECTOR;
	}
};

//! ------------------------------------------------------ */

Tddd vectorsClingToNeumannSurface(const networkPoint *p, const std::vector<T3Tddd> &next_Vrtx)
{
	auto closestXFacing = [](const Tddd &p_next_X, const double radius, const std::vector<T3Tddd> &vertices, const Tddd &n)
	{
		Tddd r = {1E+100, 1E+100, 1E+100};
		for (const auto &vertex : vertices)
		{
			if (isFacing(TriangleNormal(vertex), n, M_PI / 180 * 20))
			{
				auto intxn = IntersectionSphereTriangle_(p_next_X, radius, vertex);
				if (intxn.isIntersecting)
					if (Norm(r) >= Norm(intxn.X - p_next_X))
						r = intxn.X - p_next_X;
			}
		}
		return r;
	};

	if (next_Vrtx.empty())
		return {0., 0., 0.};
	else if (p->Neumann || p->CORNER)
	{
		std::vector<Tddd> F_clings;
		for (const auto &f : p->getFacesNeumann())
		{
			// auto [p0, p1, p2] = f->getPointsTuple();
			// auto n = TriangleNormal(p0->getXBuffer() + p0->U_BUFFER, p1->getXBuffer() + p1->U_BUFFER, p2->getXBuffer() + p2->U_BUFFER);
			auto n = f->getNormalTuple();
			auto to_closest_X = closestXFacing(p->getXBuffer() + p->U_BUFFER, p->radius, next_Vrtx, n);
			if (isFinite(to_closest_X))
				F_clings.push_back(Dot(to_closest_X, n) * n);
		}
		Tddd r = Mean(F_clings);
		if (isFinite(r))
			return r;
		else
			return {0., 0., 0.};
	}
	else
		return {0., 0., 0.};
};

double minViewRatio(const networkPoint *const p)
{
	double a = p->getSolidAngle();
	return (2 * M_PI - Min(Tdd{std::abs(a - 2 * M_PI), std::abs(2 * M_PI - a)})) / (2 * M_PI);
};

double normalVariance(const networkPoint *const p)
{
	auto n = p->getNormalDirichlet_BEM();
	double m = 0, s = 0;
	for (const auto &f : p->getFacesDirichlet())
	{
		m += (M_PI / 2. - VectorAngle(n, f->getNormalTuple())) / (M_PI / 2.);
		s += 1;
	}
	return m / s;
};

Tddd vectorTangentialShift(const networkPoint *p)
{
	auto nextX_U_Ua = [](const networkPoint *p)
	{
		return p->getXBuffer() + p->U_BUFFER;
	};

	auto next_length = [nextX_U_Ua](const networkLine *const l)
	{
		auto [p0, p1] = l->getPointsTuple();
		return Norm(nextX_U_Ua(p0) - nextX_U_Ua(p1));
	};

	auto getBaseLength = [next_length](const networkLine *line)
	{
		auto [p0, p1] = line->getPointsTuple();
		std::unordered_set<networkLine *> lc = Join(extLinesCORNER_(p0->getFaces()), extLinesCORNER_(p1->getFaces()));
		if (!lc.empty())
		{
			V_d ll;
			for (const auto &l : lc)
				ll.emplace_back(next_length(l));
			return Mean(ll);
		}
		else
		{
			// pを引っ張る力は，ノイマン面とディリクレ面で干渉しない
			auto [p0, p1] = line->getPointsTuple();
			double m = 1, s = 0;
			for (const auto &l : Join(p0->getLinesAround(), p1->getLinesAround()))
				if (line != l)
					if ((line->Dirichlet && (l->Dirichlet || l->CORNER)) || (line->Neumann && (l->Neumann || l->CORNER)) || (line->CORNER && l->CORNER))
					{
						m *= next_length(l);
						s += 1;
					}
			return std::pow(m, 1. / s);
		}
	};

	auto vectorToNextNeighborsCenter = [nextX_U_Ua](const networkPoint *const p)
	{
		double s = 0;
		Tddd ret = {0., 0., 0.};
		Tddd pX = nextX_U_Ua(p);
		/* ------------------------------------------------------ */
		if (p->CORNER)
		{
			for (const auto &l : p->getLines())
				if (l->CORNER)
				{
					ret += (nextX_U_Ua(((*l)(p))) - pX);
					s += 1;
				}
		}
		else
		{
			for (const auto &l : p->getLines())
			{
				ret += (nextX_U_Ua(((*l)(p))) - pX);
				s += 1;
			}
		}
		return ret / s;
	};

	Tddd pX = nextX_U_Ua(p);
	// EMTは節点をおしすぎることがあるようだ
	// やはり接線方向でないといけないようだ
	double c_LS = 0.2 /*0.1~0.5*/, c_EMT = 0.1;
	Tddd V_EMT = {0., 0., 0.};
	auto V_LS = vectorToNextNeighborsCenter(p) * c_LS;
	auto a = minViewRatio(p);
	if (p->CORNER)
	{
		if (a < 1. / 3.)
			return {0., 0., 0.};
		else
			return condition_Ua(V_LS, p);
	}
	else if (p->Dirichlet)
	{
		for (const auto &l : p->getLines())
			V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
		return condition_Ua(V_EMT * c_EMT + V_LS, p);
	}
	else
	{
		if (a > 2. / 3.) //比較的滑らかなノイマン面
		{
			for (const auto &l : p->getLines())
				V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
			V_EMT = condition_Ua(V_EMT * c_EMT, p);
			return condition_Ua(V_EMT + V_LS, p);
		}
		else //比較的滑らかでないノイマン面
		{
			V_netFp f;
			for (const auto &l : p->getLines())
			{
				f = l->getFaces();
				a = VectorAngle(f[0]->getNormalTuple(), f[1]->getNormalTuple()) / (2 * M_PI);
				V_EMT += a * (nextX_U_Ua((*l)(p)) - pX);
			}
			return condition_Ua(V_EMT, p);
		}
	}
};

std::tuple<Tddd, double, double, networkLine *, networkFace *> vectorToNearestAjacentSurface(const networkPoint *p)
{
	/*
	UartificialClingは，完璧にOmega(t+\delta t)に張り付くようにしなければ，
	面からはなれることで計算の破綻を招く可能性がある．
	*/
	std::tuple<Tddd, double, double, networkLine *, networkFace *> ret = {{1E+20, 1E+20, 1E+20}, 0., 0., nullptr, nullptr};
	Tddd pX = p->getXBuffer() + p->U_BUFFER;
	if (p->CORNER)
	{
		for (const auto &l : p->getLines())
			if (l->CORNER)
			{
				auto intxn = IntersectionSphereLine(pX, 1E+20, T2Tddd{p->getXBuffer(), (*l)(p)->getXBuffer()});
				if (Norm(std::get<0>(ret) - pX) >= Norm(intxn.X - pX))
				{
					std::get<0>(ret) = intxn.X;
					std::get<1>(ret) = intxn.t;
					std::get<2>(ret) = 1 - intxn.t;
					std::get<3>(ret) = l;
					std::get<4>(ret) = nullptr;
				}
			}
	}
	else
	{
		for (const auto &f : p->getFaces())
		{
			if ((p->Dirichlet && f->Dirichlet) || (p->Neumann && f->Neumann))
			{
				auto [p0, p1, p2] = f->getPointsTuple(p);
				auto intxn = IntersectionSphereTriangle_(pX, 1E+20, T3Tddd{p0->getXBuffer(), p1->getXBuffer(), p2->getXBuffer()});
				if (Norm(std::get<0>(ret) - pX) >= Norm(intxn.X - pX))
				{
					std::get<0>(ret) = intxn.X;
					std::get<1>(ret) = intxn.t0;
					std::get<2>(ret) = intxn.t1;
					std::get<3>(ret) = nullptr;
					std::get<4>(ret) = f;
				}
			}
		}
	}
	return ret;
};

void calculateVectorToSurfaceInBuffer(const Network &net, double dt)
{
	/*
	@ この方法なら，次の時刻における任意の場所でのポテンシャルを見積もることができる．
	@ このことは，任意のノイマン面上に節点を維持する上で便利である．
	@ Ω(t+δt)をまず見積もり，その面上で最適な格子配置となるように流速を修正する．
	*/
	auto Points = ToVector(net.getPoints());
	for (const auto &p : Points)
		p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};

	// for (const auto &p : Points)
	// {
	// 	Tddd v = {0., 0., 0.};
	// 	int i = 0;
	// 	// 角に限った方がいい時もある
	// 	for (const auto &q : extPointsCORNER_(Flatten(BFS(p, 1)) /*p->getNeighbors()*/))
	// 		if (q->CORNER)
	// 		{
	// 			v += q->U_update_BEM * dt;
	// 			i++;
	// 		}
	// 	if (i)
	// 		p->U_BUFFER = condition_Ua(v / i, p);
	// }

	// for (const auto &p : Points)
	// {
	// 	if (p->Neumann)
	// 	{
	// 		auto u = uNeumann(p) * dt;
	// 		auto n = p->getNormal_BEM();
	// 		p->U_BUFFER = condition_Ua(u - Dot(u, n) * n, p);
	// 	}
	// }

	//@ ------------------------------------------------------ */
	//@           次の時刻で最適な格子を目指す修正流速を計算          */
	//@ ------------------------------------------------------ */
	for (auto kk = 0; kk < 20; ++kk)
	{
		//% ------------------------------------------------------ */
		//%           　　　　 vectorTangentialShift   　 　         */
		//%          ラプラス平滑化と引っ張り合わせた接線方向にシフト      */
		//% ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			p->U_BUFFER_BUFFER = vectorTangentialShift(p);
		}

		for (const auto &p : Points)
		{
			if (isFinite(p->U_BUFFER_BUFFER))
				p->U_BUFFER += p->U_BUFFER_BUFFER;
			p->U_BUFFER_BUFFER = {0., 0., 0.};
		}

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
			if (p->CORNER || p->Dirichlet)
			{
				//% ------------------------------------------------------ */
				//%              vectorToNearestAjacentSurface             */
				//%           　　　   周辺ディリクレ面に移動        　　　　    */
				//% ------------------------------------------------------ */
				if (isFinite(p->U_BUFFER_BUFFER))
				{
					p->clungSurface = vectorToNearestAjacentSurface(p);
					p->U_BUFFER_BUFFER = (std::get<0>(p->clungSurface) - (p->getXBuffer() + p->U_BUFFER));
					//! 角のノイマン面から離れるのを防ぐ
					if (p->CORNER)
						for (const auto &f : p->getFacesNeumann())
							p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, f->getNormalTuple()) * f->getNormalTuple();
				}
			}

		for (const auto &p : Points)
			if (p->CORNER || p->Dirichlet)
			{
				if (isFinite(p->U_BUFFER_BUFFER))
					p->U_BUFFER += p->U_BUFFER_BUFFER;
				p->U_BUFFER_BUFFER = {0., 0., 0.};
			}

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			//% ------------------------------------------------------ */
			//%              vectorsClingToNeumannSurface              */
			//%              　　　近傍のノイマン面へ移動           　　　   */
			//% ------------------------------------------------------ */
			if (p->CORNER || p->Neumann)
			{
				//接触面候補の次の時刻の位置を予測
				std::unordered_set<networkFace *> Fs = p->getContactFaces();
				for (auto &f : p->getContactFaces())
					for (auto &q : f->getPoints())
						for (auto &F : q->getFaces())
							Fs.emplace(F);

				std::vector<T3Tddd> nextBodyVertices;

				for (auto &f : Fs)
				{
					auto [p0, p1, p2] = f->getPointsTuple();
					auto X0 = p0->getXtuple() + f->getNetwork()->velocityRigidBody(p0->getXtuple()) * dt;
					auto X1 = p1->getXtuple() + f->getNetwork()->velocityRigidBody(p1->getXtuple()) * dt;
					auto X2 = p2->getXtuple() + f->getNetwork()->velocityRigidBody(p2->getXtuple()) * dt;
					nextBodyVertices.emplace_back(T3Tddd{X0, X1, X2});
				}
				p->U_BUFFER_BUFFER = vectorsClingToNeumannSurface(p, nextBodyVertices);
				//! 角のディリクレ面へのめり込みを防止
				if (p->CORNER)
					p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, p->getNormalDirichlet_BEM()) * p->getNormalDirichlet_BEM();
			}
		}

		for (const auto &p : Points)
		{
			if (isFinite(p->U_BUFFER_BUFFER))
				p->U_BUFFER += p->U_BUFFER_BUFFER;
			p->U_BUFFER_BUFFER = {0., 0., 0.};
		}
	}

	std::cout << "p->U_BUFFER finished" << std::endl;
};

void calculateVectorToClingNeumann(const Network &net, const double dt)
{
	/*
	@ ノイマン面に貼り付けるための必要な調整
	*/
	auto Points = ToVector(net.getPoints());
	// for (const auto &p : Points)
	// 	p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};
	for (auto kk = 0; kk < 50; ++kk)
	{
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			//% ------------------------------------------------------ */
			//%              vectorsClingToNeumannSurface              */
			//%              　　　近傍のノイマン面へ移動           　　　   */
			//% ------------------------------------------------------ */
			if (p->CORNER || p->Neumann)
			{
				//接触面候補の次の時刻の位置を予測
				std::unordered_set<networkFace *> Fs = p->getContactFaces();
				for (auto &f : p->getContactFaces())
					for (auto &q : f->getPoints())
						for (auto &F : q->getFaces())
							Fs.emplace(F);

				std::vector<T3Tddd> nextBodyVertices;
				for (auto &f : Fs)
				{
					auto [p0, p1, p2] = f->getPointsTuple();
					auto X0 = p0->getXtuple() + f->getNetwork()->velocityRigidBody(p0->getXtuple()) * dt;
					auto X1 = p1->getXtuple() + f->getNetwork()->velocityRigidBody(p1->getXtuple()) * dt;
					auto X2 = p2->getXtuple() + f->getNetwork()->velocityRigidBody(p2->getXtuple()) * dt;
					nextBodyVertices.emplace_back(T3Tddd{X0, X1, X2});
				}
				p->U_BUFFER_BUFFER = vectorsClingToNeumannSurface(p, nextBodyVertices);
				//! 角のディリクレ面へのめり込みを防止
				if (p->CORNER)
					p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, p->getNormalDirichlet_BEM()) * p->getNormalDirichlet_BEM();
			}
		}

		for (const auto &p : Points)
		{
			if (isFinite(p->U_BUFFER_BUFFER))
				p->U_BUFFER += 0.1 * p->U_BUFFER_BUFFER;
			p->U_BUFFER_BUFFER = {0., 0., 0.};
		}
	}
	std::cout << "p->U_BUFFER finished" << std::endl;
};

	//! ------------------------------------------------------ */

#define derivatives_debug
struct derivatives
{
	// public:
	uomap_P_Tddd P_tension, P_gradPhi, P_gradPhi_tangential, P_phin_vector, P_dxdt, P_dxdt_mod, P_laplacian, P_U_dot_gradgrad_U;
	uomap_P_d P_DphiDt, P_kappa, P_pressure, P_aphiat, P_aphiant;
	double mean_surface_height_from_zero = 0;
	~derivatives(){
		// std::cout << "derivatives 破棄" << std::endl;
	};
	derivatives(const Network &net, const double dt, bool adjust_mesh = false)
	{
		std::unordered_set<networkPoint *> Points = net.getPoints();
		P_tension.reserve(Points.size());
		P_gradPhi.reserve(Points.size());
		P_gradPhi_tangential.reserve(Points.size());
		P_phin_vector.reserve(Points.size());
		P_dxdt.reserve(Points.size());
		P_dxdt_mod.reserve(Points.size());
		P_laplacian.reserve(Points.size());
		P_U_dot_gradgrad_U.reserve(Points.size());
		P_DphiDt.reserve(Points.size());
		P_kappa.reserve(Points.size());
		P_pressure.reserve(Points.size());
		P_aphiat.reserve(Points.size());
		P_aphiant.reserve(Points.size());
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
			this->P_DphiDt[p] = 0.;
		}
		this->P_tension = this->P_laplacian = this->P_U_dot_gradgrad_U = this->P_dxdt_mod = this->P_dxdt = this->P_phin_vector = this->P_gradPhi_tangential = this->P_gradPhi;
		this->P_kappa = this->P_pressure = this->P_aphiant = this->P_aphiat = this->P_DphiDt;

		int c = 0;
		for (const auto &p : Points)
		{
			if (p->Dirichlet)
			{
				mean_surface_height_from_zero += p->height();
				c++;
			}
		}

		int depthlimit = 9;
		auto Faces = net.getFaces();

		mean_surface_height_from_zero /= c;
		auto pointsbegin = Points.begin();

		/* ------------------------------------------------------ */
		// //@ pullしたい構造物と，接触を感知した流体格子は角に最大の値を与える．
		std::unordered_set<Network *> nets;
		for (const auto &p : Points)
		{
			auto contactFs = p->getContactFaces();
			if (!contactFs.empty())
				for (const auto &cf : contactFs)
				{
					auto contact_net = cf->getNetwork();
					if (p->getNetwork() != contact_net)
						nets.emplace(contact_net);
				}
		}

#ifdef derivatives_debug
		std::cout << "流速の計算" << std::endl;

#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			// std::cout << "曲率の計算" << std::endl;
			if (!isFinite(p->phiphin))
			{
				std::cout << "p->phiphinはfiniteではない！！" << std::endl;
				std::cout << "p->phiphin = " << p->phiphin << std::endl;
				if (p->Neumann)
					std::cout << "p->Neumann" << std::endl;
				if (p->CORNER)
					std::cout << "p->CORNER" << std::endl;
				if (p->Dirichlet)
					std::cout << "p->CORNER" << std::endl;
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			}
			V_netPp ps = Flatten(BFS(p, 3));
			auto interpNormals = InterpolationVectorRBF(ToVector(extX(ps)), ToVector(extNormals(ps)), p->getX());
			p->kappa_BEM = interpNormals.div(p->getX()) / 2.; //中心方向法線ベクトルの場合，マイナスをつける．
			//! https://en.wikipedia.org/wiki/Mean_curvature
			auto grad_N_S_full = gradPhi(p);
			p->U_BEM = std::get<2>(grad_N_S_full);
			p->U_tangential_BEM = std::get<1>(grad_N_S_full);
			p->U_normal_BEM = std::get<0>(grad_N_S_full);
			/* -------------------- おおよそのアップデート流速 ------------------- */
			//@ U_update_BEM は first guess
			if (p->CORNER)
			{
				Tddd tang = Normalize(Cross(p->getNormalDirichlet_BEM(), p->getNormalNeumann_BEM()));
				p->U_update_BEM = p->U_BEM - Dot(p->U_BEM, tang) * tang;
			}
			else if (p->Neumann)
				p->U_update_BEM = uNeumann(p);
			else
				p->U_update_BEM = p->U_BEM;
			/* ------------------------------------------------------ */
			for (const auto &q : Points)
				q->X_BUFFER = q->getXtuple() + q->U_update_BEM * dt;
			//@ この後U_update_BEMをclingなどを使って修正する
		}
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
#ifdef derivatives_debug
		std::cout << "DphiDtを計算" << std::endl;
#endif
		// double gamma = 72.75 * 1E-3;	  //[N/m] 水20度
		// double gravity = _GRAVITY_;		  //[m/s2]
		// double density = _WATER_DENSITY_; //[kg/m3]
		// double nu = 0.01005 / density;

		for (auto &[p, v] : this->P_kappa)
			v = p->kappa_BEM;

		for (auto &[p, v] : this->P_gradPhi)
			v = p->U_BEM;

		for (auto &[p, v] : this->P_laplacian)
			v = p->laplacian_U_BEM;

		for (auto &[p, v] : this->P_dxdt_mod)
			v = p->U_update_BEM;

		for (auto &[p, v] : this->P_gradPhi_tangential)
			v = p->U_tangential_BEM;

		for (auto &[p, v] : this->P_phin_vector)
			v = p->U_normal_BEM;

		/*
			この方法は，ノイマン境界条件のclingにおいて，とても自然に無理なく応用できる．
			ディリクレ境界条件に関しては，計算後に補間によってリグリッドしてもいいかもしれない．
		*/
		/* ------------------------------------------------------ */
		/*
			X_BUFFERには，U_update_BEMで単純に予測した節点が保存されている．
			予測したディリクレ面は正しいと考える．
			予測したノイマン面よりも，実際に物体を移動させて作った面の方が正しい．
			そこで，ノイマン面と角点に関しては，物体を移動させて作った面に貼り付ける．
		*/
		for (const auto &p : Points)
			p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};
		if (adjust_mesh)
			calculateVectorToSurfaceInBuffer(net, dt);
		/* ------------------------------------------------------ */
		//@ ノイマン面に張り付くようにN，Cでの節点のU_BUFFERを修正
		calculateVectorToClingNeumann(net, dt); //@必須
		for (const auto &p : Points)
			p->U_update_BEM += p->U_BUFFER / dt;
		/* ------------------------------------------------------ */
		for (const auto &p : Points)
		{
			//! p->U_update_BEM
			//!  (1) fitToNeumannVelocity　を満たしている
			//
			//! p->U_BUFFER
			//!  (1) fitToNeumannVelocity
			//!  (2) delete_normal
			//!  (3) along_surface_CORNER を満たしている
			P_tension[p] = p->U_BUFFER;
		}

#ifdef derivatives_debug
		std::cout << " DONE" << std::endl;
#endif

		/* ------------------------------------------------------ */

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			// b@ ---------------------- 壁にくっつける修正 --------------------- */
			// p->U_update_BEM += p->U_cling_to_Neumann;
			//@ ------------------------------------------------------------ */
			this->P_dxdt[p] = p->U_update_BEM; //流速
			// //!この場合マイナスでないと，上の部分が半たんする
			auto DphiDt = p->DphiDt(p->U_update_BEM, 0.);
			//!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！
			this->P_DphiDt[p] = DphiDt; // update用
			if (p->Neumann || p->CORNER)
			{
				auto n = p->getNormalNeumann_BEM();
				T6d VW = velocity_from_Neumann_surface(p);
				Tddd angular_velocity = {std::get<3>(VW), std::get<4>(VW), std::get<5>(VW)};
				auto Q = Quaternion();
				auto dQdt = Q.d_dt(angular_velocity);
				std::get<1>(p->phiphin_t) = accel_normal_from_Neumann_surface(p) - Dot(n, Dot(p->U_BEM, grad_U_LinearElement(p)));
				std::get<1>(p->phiphin_t) += Dot(n, Dot(velocity_normal_from_Neumann_surface(p) - p->U_BEM, dQdt.Rv()));
				/*
				∇U=∇∇f={{fxx, fyx, fzx},
						{fxy, fyy, fzy},
						{fxz, fyz, fzz}}
				なので，∇∇f=∇∇f^T
				*/
			}
			if (p->Dirichlet || p->CORNER)
				std::get<0>(p->phiphin_t) = p->DphiDt(0.) - Dot(p->U_BEM, p->U_BEM); //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！

			this->P_aphiat[p] = std::get<0>(p->phiphin_t);	//!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！
			this->P_aphiant[p] = std::get<1>(p->phiphin_t); //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！

			this->P_U_dot_gradgrad_U[p] = Dot(p->U_BEM, grad_U_LinearElement(p));
			this->P_pressure[p] = p->pressure_BEM;

			// 10000. * (-1 / 2. * Vphi_Vphi - gravity * (std::get<2>(p->getXtuple())) - aphiat);
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

void addIG(std::unordered_map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign)
{
	std::unordered_map<netP *, Tdd>::iterator it;
	if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
		std::get<0>(it->second) += std::get<0>(igign); // phiは忘れずに計算
	else
		P_phiphin[P_IN] = {std::get<0>(igign), 0.}; // phiは忘れずに計算
}

void addIGn(std::unordered_map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign)
{
	std::unordered_map<netP *, Tdd>::iterator it;
	if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
		std::get<1>(it->second) += std::get<1>(igign); // phiは忘れずに計算
	else
		P_phiphin[P_IN] = {0., std::get<1>(igign)}; // phiは忘れずに計算
}

void addIG(std::map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign)
{
	std::map<netP *, Tdd>::iterator it;
	if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
		std::get<0>(it->second) += std::get<0>(igign); // phiは忘れずに計算
	else
		P_phiphin[P_IN] = {std::get<0>(igign), 0.}; // phiは忘れずに計算
}

void addIGn(std::map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign)
{
	std::map<netP *, Tdd>::iterator it;
	if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
		std::get<1>(it->second) += std::get<1>(igign); // phiは忘れずに計算
	else
		P_phiphin[P_IN] = {0., std::get<1>(igign)}; // phiは忘れずに計算
}

struct BEM_BVP
{
	std::vector<networkPoint *> Points;
	std::vector<networkFace *> Faces;
	const bool Neumann = false;
	const bool Dirichlet = true;

#if defined(use_lapack)
	lapack_lu *lu;
#else
	ludcmp_parallel *lu;
#endif
	VV_d mat_ukn;
	VV_d mat_kn;
	V_d knowns;
	// V_d unknowns;
	// std::vector<networkPoint *> vec_P;

	std::vector<pair_PB> vec_P;
	BEM_BVP(){};
	~BEM_BVP() { delete this->lu; };

	void solveForPhiPhin_t() const
	{
		// b* ------------------------------------------------------ */
		// b*                         phiphin_t                      */
		// b* ------------------------------------------------------ */
		V_d knowns(vec_P.size());
		V_d phiORphin(vec_P.size());
		for (auto i = 0; i < this->vec_P.size(); i++)
		{
			auto [p, DorN] = vec_P[i];
			if (DorN == Dirichlet)
				knowns[i] = std::get<0>(p->phiphin_t);
			else
				knowns[i] = std::get<1>(p->phiphin_t);
		}
		this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
		/* ------------------------------------------------------ */
		for (auto i = 0; i < this->vec_P.size(); i++)
		{
			// auto p = this->vec_P[i];
			// if (p->Neumann)
			auto [p, DorN] = vec_P[i];
			if (DorN == Dirichlet)
				std::get<1>(p->phiphin_t) = phiORphin[i];
			else
				std::get<0>(p->phiphin_t) = phiORphin[i];
			p->pressure_BEM = -std::get<0>(p->phiphin_t) - _GRAVITY_ * p->height() - Dot(p->U_BEM, p->U_BEM) / 2.;
			p->pressure_BEM *= _WATER_DENSITY_;
			p->pressure = p->pressure_BEM;
		}
	};

	void solve(const Network &water, const Buckets<networkPoint> &FMM_BucketsPoints, const Buckets<networkFace> &FMM_BucketsFaces)
	{
		this->Points = ToVector(water.getPoints());
		this->Faces = ToVector(water.getFaces());
		using map_P_Vd = std::map<netP *, V_d>;
		//* ------------------------------------------------------ */
		//%                     各点で方程式を作る場合                */
		//* ------------------------------------------------------ */
		std::cout << "各点で方程式を作る場合" << std::endl;
		map_P_P_Tdd P_P_IGIGn;

		//@ 各バケツでのモーメントを次数別に保存する．(ユニーク) p->{k,m,Yn,Y}ベクトル
		using uo_P_uoTiiTdd = std::unordered_map<networkPoint *, std::unordered_map<Tii /*k,m*/, Tdd /*YYn*/>>;
		using V_uo_P_uoTiiTdd = std::vector<uo_P_uoTiiTdd>;
		using VV_uo_P_uoTiiTdd = std::vector<V_uo_P_uoTiiTdd>;
		using VVV_uo_P_uoTiiTdd = std::vector<VV_uo_P_uoTiiTdd>;

		VVV_uo_P_uoTiiTdd VVV_P_km_YnY(FMM_BucketsPoints.xsize, VV_uo_P_uoTiiTdd(FMM_BucketsPoints.ysize, V_uo_P_uoTiiTdd(FMM_BucketsPoints.zsize))); //途中でreverseしているので注意
		int ii = 0;
		//
		map_pairPB_pairPB_Tdd PDN_PDN_IGIGn;
		map_pairPB_Tdd init_PDN_IGIGn;
		/* ------------------------------------------------------ */
		/*                     init_PDN_IGIGn                     */
		/* ------------------------------------------------------ */
		for (const auto &org : Points)
		{
			if (org->CORNER)
			{
				init_PDN_IGIGn[{org, Dirichlet}] = {0., 0.};
				init_PDN_IGIGn[{org, Neumann}] = {0., 0.}; // b!
			}
			else if (org->Dirichlet)
				init_PDN_IGIGn[{org, Dirichlet}] = {0., 0.};
			else if (org->Neumann)
				init_PDN_IGIGn[{org, Neumann}] = {0., 0.};
		}
		/* ------------------------------------------------------ */
		/*                        P_P_IGIGn                       */
		/* ------------------------------------------------------ */
		for (const auto &org : Points)
		{
			P_P_IGIGn[org] = {}; //! P_P_IGIGn[row][col] = {IG,IGn}

			//! これで，PDN_PDN_IGIGnは，行の数がCORNERの分だけ倍になった．
			if (org->CORNER)
			{
				PDN_PDN_IGIGn[{org, Dirichlet}] = init_PDN_IGIGn;
				PDN_PDN_IGIGn[{org, Neumann}] = init_PDN_IGIGn; // b!
			}
			else if (org->Dirichlet)
				PDN_PDN_IGIGn[{org, Dirichlet}] = init_PDN_IGIGn;
			else if (org->Neumann)
				PDN_PDN_IGIGn[{org, Neumann}] = init_PDN_IGIGn;
		}
// #define single_layer_FMM
#if defined(single_layer_FMM)
		// double spacing = Mean(extLength(water.getLines())) * 10;
		// Buckets<networkFace> FMM_BucketsFaces(water.bounds(), spacing);
		// Buckets<networkPoint> FMM_BucketsPoints(water.bounds(), spacing);
		// for (const auto &f : water.getFaces())
		// 	FMM_BucketsFaces.add(f->getXtuple(), f);
		// for (const auto &p : water.getPoints())
		// 	FMM_BucketsPoints.add(p->getXtuple(), p);
		std::cout << "FMMのモーメントを計算" << std::endl;
		std::cout << "xsize = " << FMM_BucketsPoints.xsize << std::endl;
		std::cout << "ysize = " << FMM_BucketsPoints.ysize << std::endl;
		std::cout << "zsize = " << FMM_BucketsPoints.zsize << std::endl;

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (auto i = 0; i < FMM_BucketsPoints.buckets.size(); ++i)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			for (auto j = 0; j < FMM_BucketsPoints.buckets[i].size(); ++j)
			{
				for (auto k = 0; k < FMM_BucketsPoints.buckets[i][j].size(); ++k)
				{
					// auto &P_km_YnY = VVV_P_km_YnY[i][j][k]; //@このバケツでのモーメント．次数別に保存m,k
					for (const auto &F : FMM_BucketsFaces.buckets[i][j][k])
					{
						for (const auto &P : FMM_BucketsPoints.buckets[i][j][k])
						{
							Tddd center = FMM_BucketsFaces.indices2X_center(i, j, k);
							for (const auto &[p, km_YYn] : BEM::calc_P_MomentQuadTuple(F, center)) //@このバケツでのモーメント．次数別に保存m,kに計算する
							{
								/*uoTiiTdd*/
								VVV_P_km_YnY[i][j][k][p] += km_YYn;
								// P_km_YnY[p] += km_YYn;
							}
						}
					}
				}
			}
		}
		std::cout << "FMM実行" << std::endl;
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (auto i = 0; i < FMM_BucketsPoints.buckets.size(); ++i)
			for (auto j = 0; j < FMM_BucketsPoints.buckets[i].size(); ++j)
				for (auto k = 0; k < FMM_BucketsPoints.buckets[i][j].size(); ++k)
				{
					for (const auto &P : FMM_BucketsPoints.buckets[i][j][k])
#ifdef _OPENMP
#pragma omp single nowait
#endif
					{
						// std::cout << "(i,j,k) = (" << i << "," << j << "," << k << ")" << std::endl;
						auto &P_igign = P_P_IGIGn[P];
						for (auto I = 0; I < FMM_BucketsFaces.buckets.size(); ++I)
							for (auto J = 0; J < FMM_BucketsFaces.buckets[I].size(); ++J)
								for (auto K = 0; K < FMM_BucketsFaces.buckets[I][J].size(); ++K)
								{
									if ((i - 1 <= I && I <= i + 1) || (j - 1 <= J && J <= j + 1) || (k - 1 <= K && K <= k + 1))
									{
										//@ この点に，比較的近い面の積分
										// std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;
										for (const auto &F : FMM_BucketsFaces.buckets[I][J][K])
										{
											for (const auto &[p, igign] : BEM::calc_P_IGIGnQuadTuple_mod(F, P))
											{
												addIG(P_igign, p, igign);
												if (p != P)
												{
													addIGn(P_igign, p, igign);
													//@ 対角成分は，リジッドモード法を使って計算できる
													addIGn(P_igign, P, -igign);
												}
											}
										}
									}
									else
									{ //@ この点に，比較的遠い面の積分
										Tddd O = FMM_BucketsFaces.indices2X_center(i, j, k);
										for (const auto &[p, kmYYn] : VVV_P_km_YnY[i][j][k])
										{
											Tdd igign = {0., 0.};
											for (const auto &[km, YYn] : kmYYn)
											{
												auto [k, m] = km;
												auto [nr, theta, psi] = ToSphericalCoodrinates(p->getXtuple() - O);
												igign += YYn * real_sph_scale_ommited(k, m, theta, psi) / std::pow(nr, k + 1);
											}
											addIG(P_igign, p, igign);
											if (p != P)
											{
												addIGn(P_igign, p, igign);
												//@ 対角成分は，リジッドモード法を使って計算できる
												addIGn(P_igign, P, -igign);
											}
										}
									}
								}
					}
				}
#else

		/* ------------------------------------------------------ */

		// #define quad_element
#define linear_element
		// #define liear_and_quad_element

#ifdef _OPENMP
		std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;
#ifdef quad_element
		std::cout << "quad" << std::endl;
#elif defined(linear_element)
		std::cout << "linear" << std::endl;
#endif

#pragma omp parallel
#endif
		for (auto &[p0_bool, P_igign] : PDN_PDN_IGIGn)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			auto p0 = std::get<0>(p0_bool);
			std::map<pair_PB, Tdd>::iterator it;
			double p0_ign = 0.;
			pair_PB p_bool;

			for (const auto &integrating_f : Faces)
			{

				/*         *   *       */
				/*        / \ / \      */
				/*       *---*---*     */
				/*      / \ /F\ / \    */
				/*     *---*---*---*   */
				/*        / \ / \      */
				/*       *---*---*     */
				/*
				積分は各節点の変数と重みによって表される．
				その重みを与える．
				*/
#if defined(quad_element)
				for (const auto &[p, igign] : BEM::calc_P_IGIGnQuadTuple_mod(integrating_f, p0))
				{
					if (p->CORNER)
					{
						if (integrating_f->Dirichlet)
							p_bool = {p, Dirichlet};
						else if (integrating_f->Neumann)
							p_bool = {p, Neumann};
					}
					else
						p_bool = {p, p->Dirichlet};

					P_igign[p_bool] += igign;

					if (p != p0)
						p0_ign -= std::get<1>(igign);
				}
#elif defined(linear_element)
				auto func = [&](const std::tuple<netP *, Tdd> &p_igign)
				{
					const auto [p, igign] = p_igign;
					if (p->CORNER)
					{
						if (integrating_f->Dirichlet)
							P_igign[{p, Dirichlet}] += igign;
						else // if (integrating_f->Neumann)
							P_igign[{p, Neumann}] += igign;
					}
					else
						P_igign[{p, p->Dirichlet}] += igign;
					if (p != p0)
						p0_ign -= std::get<1>(igign);
				};
				const auto [q0igign, q1igign, q2igign] = BEM::calc_P_IGIGnLinear3Tuples(integrating_f, p0);
				func(q0igign);
				func(q1igign);
				func(q2igign);
#elif defined(linear_element2) for (const auto &[p, igign] : BEM::calc_P_IGIGnLinearTuple(integrating_f, p0))
				{
					if (p->CORNER)
					{
						if (integrating_f->Dirichlet)
							p_bool = {p, Dirichlet};
						else if (integrating_f->Neumann)
							p_bool = {p, Neumann};
					}
					else
						p_bool = {p, p->Dirichlet};

					P_igign[p_bool] += igign;

					if (p != p0)
						p0_ign -= std::get<1>(igign);
				}
#endif
			}
			std::get<1>(P_igign[p0_bool]) = p0_ign;
		}
#endif
		std::cout << "並列化 DONE" << std::endl;
		std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;

		mat_ukn = VV_d(PDN_PDN_IGIGn.size(), V_d(PDN_PDN_IGIGn.size(), 0.)); //! ok
		mat_kn = VV_d(PDN_PDN_IGIGn.size(), V_d(PDN_PDN_IGIGn.size(), 0.));	 //! ok

		// 順番を間違えないようにベクトルを作成
		vec_P = std::vector<pair_PB>(PDN_PDN_IGIGn.size());
		auto i = 0;
		for (auto it = PDN_PDN_IGIGn.begin(); it != PDN_PDN_IGIGn.end(); ++it)
			vec_P[i++] = it->first;

		{
			// b@ ------------------------------------------------------ */
			// b@                 系数行列mat_ukn．mat_knの計算             */
			// b@ ------------------------------------------------------ */

#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (auto i = 0; i < vec_P.size(); ++i)
			{
				/*@ mat_ukn, mat_kn, vec_P
					+--+--+--+
				|	|--+--+--|
				V	|--+--+--|
					+--+--+--+
				*/
				auto PDN_IGIGn = PDN_PDN_IGIGn.at(vec_P[i]);
				auto p = std::get<0>(vec_P[i]);
				auto DorN = std::get<1>(vec_P[i]);
				/* mat_ukn, mat_kn, vec_P
					   ====>
					+--+--+--+
					|--+--+--|
					|--+--+--|
					+--+--+--+
				*/
				if (!(p->CORNER && DorN == Neumann /*変更する対象の行*/)) //! OK
				{
					for (auto j = 0; j < vec_P.size(); ++j) //! OK
					{
						auto igign = PDN_IGIGn.at(vec_P[j]);
						/* --------------------------------------- */
						// 未知変数の係数行列は左，既知変数の係数行列は右
						if (std::get<1>(vec_P[j]) == Neumann)
							igign = {-std::get<1>(igign), -std::get<0>(igign)};
						//% IGIGn は 左辺に IG*φn が右辺に IGn*φ が来るように計算しているため，移項する場合，符号を変える必要がある．
						/* --------------------------------------- */
						mat_ukn[i][j] = std::get<0>(igign);
						mat_kn[i][j] = std::get<1>(igign);
					}
				}
			}

			/* ------------------------------------------------------ */
			double maxpp = 0;
			for (auto i = 0; i < mat_ukn.size(); i++)
			{
				// LUするのはmat_uknだけなので，mat_knの最大値を使う必要はない
				//  if (maxpp < std::abs(mat_kn[i][i]))
				//  	maxpp = std::abs(mat_kn[i][i]);
				if (maxpp < std::abs(mat_ukn[i][i]))
					maxpp = std::abs(mat_ukn[i][i]);
			}
			/* ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (auto i = 0; i < vec_P.size(); ++i)
			{
				auto PDN_IGIGn = PDN_PDN_IGIGn.at(vec_P[i]);
				auto p = std::get<0>(vec_P[i]);
				auto DorN = std::get<1>(vec_P[i]);
				if (p->CORNER && DorN == Neumann /*変更する対象の行*/) //! OK
				{
					for (auto j = 0; j < vec_P.size(); ++j)
					{
						auto q = std::get<0>(vec_P[j]);
						if (p == q)
						{
							if (std::get<1>(vec_P[j]) == Neumann)
							{
								mat_ukn[i][j] = maxpp; //φの系数
								mat_kn[i][j] = 0;	   //φnの系数
							}
							else
							{
								mat_ukn[i][j] = 0;	  //φnの系数
								mat_kn[i][j] = maxpp; //φの系数移行したからマイナス？　いいえ，移項を考慮した上でこれでいい．
							}
						}
						else
						{
							mat_ukn[i][j] = 0;
							mat_kn[i][j] = 0;
						}
					}
				}
			}

			if (!isFinite(mat_ukn))
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mat_ukn is not finite");
			if (!isFinite(mat_kn))
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mat_kn is not finite");

			// b@ ------------------------------------------------------ */
		}

		{
			// b$ ------------------------------------------------------ */
			// b$　　　                   knownsの計算　                    */
			// b$ ------------------------------------------------------ */
			/**
			 * Dot(mat_ukn,phiORphin) = Dot(mat_kn,knowns)
			 * => phiORphin = Dot(mat_ukn^-1, Dot(mat_kn,knowns))
			 */
			knowns = V_d(vec_P.size()); //! ok
										//! Pointsの順番と合わせてとるように注意
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (auto i = 0; i < vec_P.size(); i++)
			{
				auto p = std::get<0>(vec_P[i]);
				auto DorN = std::get<1>(vec_P[i]);
				if (DorN == Dirichlet) //! OK
					knowns[i] = p->phi_Dirichlet = std::get<0>(p->phiphin);
				else
					knowns[i] = p->phin_Neumann; // std::get<1>(p->phiphin);
			}

			// b$ ------------------------------------------------------ */
		}

		{
			// b% ------------------------------------------------------ */
			// b%                       境界積分方程式を解く                 */
			// b% ------------------------------------------------------ */
			std::cout << "--------------------- 境界積分方程式を解く ---------------------" << std::endl;
			/* ------------------------------------------------------ */
			V_d phiORphin(knowns.size());
			std::cout << "phiORphin.size()= " << phiORphin.size() << std::endl;
			std::cout << "  mat_kn.size() = " << mat_kn.size() << std::endl;
#if defined(solve_equations_on_all_points) || defined(solve_equations_on_all_points_RBF) || defined(solve_equations_on_all_points_rigid_mode)
			//* 未知変数の計算
#if defined(use_lapack)
			std::cout << "lapack lu decomposition" << std::endl;
			this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/);
#else
			std::cout << "parallel lu decomposition" << std::endl;
			this->lu = new ludcmp_parallel(mat_ukn /*未知の行列係数（左辺）*/);
#endif
			std::cout << "try to solve" << std::endl;
			this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#else
			//* 未知変数の計算
			std::cout << "SVD decomposition" << std::endl;
			SVD svd(mat_ukn);
			svd.solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (auto i = 0; i < vec_P.size(); ++i)
			{
				auto p = std::get<0>(vec_P[i]);
				auto DorN = std::get<1>(vec_P[i]);
				if (DorN == Dirichlet)
					std::get<1>(p->phiphin) = p->phin_Dirichlet = phiORphin[i];
				else
					std::get<0>(p->phiphin) = p->phi_Neumann = phiORphin[i];
			}
			if (!isFinite(phiORphin))
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "phiORphin is not finite");
		}
	};
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

V_d meanX(const std::unordered_set<networkPoint *> &ps)
{
	V_d ret(3, 0.);
	for (const auto &p : ps)
		ret += p->getX();
	return ret / ((double)(ps.size()));
};
void normalizePhi(const std::unordered_set<networkPoint *> &ps)
{
	double ret = 0.;
	for (const auto &p : ps)
		ret += std::get<0>(p->phiphin);
	ret /= ((double)(ps.size()));
	for (const auto &p : ps)
		std::get<0>(p->phiphin) -= ret;
};
/* ------------------------------------------------------ */

// b! ------------------------------------------------------ */
// b!           格子のdivide, merge．それに伴うΦ，Φnの付与       */
// b! ------------------------------------------------------ */
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

Tdd estimate_phiphin(const networkLine *const l)
{
	// auto fs = l->getFaces();
	// interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_0(fs[0], l);
	// interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_1(fs[1], l);
	// auto phi0 = Dot(intp_l0_0.N(.5, .5), ToPhi(intp_l0_0.Points));
	// auto phi1 = Dot(intp_l0_1.N(.5, .5), ToPhi(intp_l0_1.Points));
	// auto phin0 = Dot(intp_l0_0.N(.5, .5), ToPhin(intp_l0_0.Points));
	// auto phin1 = Dot(intp_l0_1.N(.5, .5), ToPhin(intp_l0_1.Points));
	// return {(phi0 + phi1) / 2., (phin0 + phin1) / 2.};
	//
	auto [a, b] = l->getPointsTuple();
	return (a->phiphin + b->phiphin) / 2.;
};

/* ------------------------------------------------------ */

void remesh(Network &water,
			const Tdd &limit_angle_D,
			const Tdd &limit_angle_N,
			bool force = false,
			int max_count = 100)
{
	std::cout << "remeshing" << std::endl;
	water.setBounds();
	double mean_length = Mean(extLength(water.getLines()));
	bool isfound = false, ismerged = false;
	int count = 0;
	// double lim_degree_Neumann = limit_angle;
	// double lim_degree = limit_angle;

	networkLine *l;
	Tddd X, V;
	networkPoint *q;
	double meanArea;
	V_netFp Fs;
	Tdd phiphin;
	V_netLp lines;
	double local_mean_length;
	do
	{
		// なくなるまでやるか？
		isfound = false;
		ismerged = false;
		/* ------------------------------------------------------ */
		for (const auto &p : RandomSample(ToVector(water.getPoints())))
		{
			/* ------------------------------------------------------ */
			meanArea = Mean(p->getFaceAreas());
			if (false)
				for (const auto &f : p->getFaces())
				{
					if (f->getArea() / meanArea < 1E-3)
					{
						p->sortLinesByLength();
						l = *(p->getLines().rbegin());
						phiphin = estimate_phiphin(l);
						auto [a, b] = l->getPointsTuple();
						X = (a->getXtuple() + b->getXtuple()) / 2.;
						q = l->divide();
						q->phiphin = phiphin;
						q->setX(X);
						isfound = true;
						break;
					}
				}
			if (isfound)
				break;
			/* ------------------------------------------------------ */
			//* ------------------------------------------------------ */
			//*                立体角が小さすぎる場合merge                 */
			//* ------------------------------------------------------ */
			if (false)
				if (p->getSolidAngle() < 4 * M_PI / 100. || (4 * M_PI - p->getSolidAngle()) < 4 * M_PI / 100.)
				{
					p->sortLinesByLength();
					l = p->getLines()[0];
					//@ case2 lの面全体を考慮
					// Tdd phiphin = phiphin_from_faces(l);
					//@ case3 lの点だけを考慮
					// Tdd phiphin = phiphin_from_points(l);
					//@ case4 lの面全体を考慮
					phiphin = estimate_phiphin(l);
					/* ------------------------------------------------------ */
					// auto [a, b] = l->getPointsTuple();
					// if (a->CORNER)
					// 	b = a;
					// else if (b->CORNER)
					// 	a = b;
					// Tddd V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
					// for (const auto &f : a->getContactFaces())
					// 	V -= f->getNormalTuple() * Dot(V, f->getNormalTuple());
					// for (const auto &f : b->getContactFaces())
					// 	V -= f->getNormalTuple() * Dot(V, f->getNormalTuple());
					// Tddd X = V + a->getXtuple();
					// auto q = l->merge();
					// q->phiphin = phiphin;
					// q->setX(X);
					// ismerged = true;
					// break;
					/* ------------------------------------------------------ */
					auto [a, b] = l->getPointsTuple(p);
					X = b->getXtuple();
					q = l->merge();
					q->phiphin = phiphin;
					q->setX(X);
					ismerged = true;
					break;
				}
			//! ------------------------------------------------------ */
			//!             辺の長さが長すぎるまたは短すぎる場合             */
			//! ------------------------------------------------------ */
			local_mean_length = Mean(extLength(extractLines(Flatten(BFS(p, 2)))));
			lines = p->getLines();
			sortByLength(lines);
			for (const auto &l : Reverse(lines))
			{
				auto [p0, p1] = l->getPointsTuple();
				Fs = l->getFaces();
				//@ ------------------------------------------------------ */
				// if (l->length() > mean_length * 1.75 /*長すぎる*/)
				if (l->length() > local_mean_length * 1.5 /*長すぎる*/)
				{
					//@ case2 lの面全体を考慮
					// Tdd phiphin = phiphin_from_faces(l);
					//@ case3 lの点だけを考慮
					// Tdd phiphin = phiphin_from_points(l);
					//@ case4 lの面全体を考慮
					phiphin = estimate_phiphin(l);
					/* ------------------------------------------------------ */
					q = l->divide();
					q->phiphin = phiphin;
					isfound = true;
					break;
				}
				//@ ------------------------------------------------------ */
				if (l->length() < local_mean_length / 20.)
				{
					/*
					b!マージの原則：マージによってノイマン面は変形してはいけない．
					*/
					//@ case5 ２点の平均，移動位置は，ノイマンを崩さない方向：ノイマン面の法線方向成分には移動しない．
					auto [a, b] = l->getPointsTuple();
					if (!(a->CORNER && b->CORNER))
					{
						if (a->CORNER)
							b = a;
						else if (b->CORNER)
							a = b;
					}
					phiphin = (a->phiphin + b->phiphin) / 2.;
					V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
					V -= a->getNormalTuple() * Dot(V, a->getNormalTuple());
					V -= b->getNormalTuple() * Dot(V, b->getNormalTuple());
					// for (const auto &f : a->getContactFaces())
					// 	V -= f->getNormalTuple() * Dot(V, f->getNormalTuple());
					// for (const auto &f : b->getContactFaces())
					// 	V -= f->getNormalTuple() * Dot(V, f->getNormalTuple());
					X = V + a->getXtuple();
					q = l->merge();
					q->phiphin = phiphin;
					q->setX(X);
					/* ------------------------------------------------------ */
					ismerged = true;
					break;
				}
				//@ ------------------------------------------------------ */
				auto [min, max] = MinMax(extAreas(Fs));
				if (false && min / max < 1 / 100.)
				{
					//@ case5 ２点の平均，移動位置は，ノイマンを崩さない方向：ノイマン面の法線方向成分には移動しない．
					auto [a, b] = l->getPointsTuple();
					if (!(a->CORNER && b->CORNER))
					{
						if (a->CORNER)
							b = a;
						else if (b->CORNER)
							a = b;
					}
					phiphin = (a->phiphin + b->phiphin) / 2.;
					V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
					V -= a->getNormalTuple() * Dot(V, a->getNormalTuple());
					V -= b->getNormalTuple() * Dot(V, b->getNormalTuple());
					// for (const auto &f : a->getContactFaces())
					// 	V -= f->getNormalTuple() * Dot(V, f->getNormalTuple());
					// for (const auto &f : b->getContactFaces())
					// 	V -= f->getNormalTuple() * Dot(V, f->getNormalTuple());
					X = V + a->getXtuple();
					q = l->merge();
					q->phiphin = phiphin;
					q->setX(X);
					/* ------------------------------------------------------ */
					ismerged = true;
					break;
				}
				if (ismerged || isfound)
					break;
			}

			// for (const auto &l : water.getLines())
			// {
			// 	auto [p0, p1] = l->getPointsTuple();
			// 	Fs = l->getFaces();
			// 	if (!l->CORNER)
			// 		if (force)
			// 		{
			// 			isfound = l->flipIfTopologicalyBetter((l->Neumann ? lim_degree_Neumann : lim_degree));
			// 			break;
			// 		}
			// 		else
			// 		{
			// 			isfound = l->flipIfBetter((l->Neumann ? lim_degree_Neumann : lim_degree));
			// 		}
			// }

			if (ismerged || isfound)
				break;
		}
	} while ((ismerged || isfound) && count++ < max_count);

	flipIf(water, limit_angle_D, limit_angle_N, force);
};
// b! ------------------------------------------------------ */
// b! ------------------------------------------------------ */
// b! ------------------------------------------------------ */

bool isNeumann(const networkFace *const f)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	auto faces_p0 = p0->getContactFaces();
	auto faces_p1 = p1->getContactFaces();
	auto faces_p2 = p2->getContactFaces();
	bool isNeumann = f->isThereAnyFacingFace(faces_p0, M_PI / 9.) &&
					 f->isThereAnyFacingFace(faces_p1, M_PI / 9.) &&
					 f->isThereAnyFacingFace(faces_p2, M_PI / 9.);
	return isNeumann;
};

void setBoundaryConditions(Network &water, const std::vector<Network *> &objects)
{
	auto radius = Mean(extLength(water.getLines()));
	auto Points = ToVector(water.getPoints());
	Print("makeBucketFaces", Green);
	for (const auto &obj : objects)
		obj->makeBucketFaces(radius);

	// b% -------------------------------------------------------- */
	// b%            境界条件（角点・ディリクレ・ノイマン）の決定         */
	// b% -------------------------------------------------------- */
	// b% step1 点の衝突の判定
	std::cout << "step1 点の衝突の判定" << std::endl;
	for (const auto &p : Points)
		p->clearContactFaces();
	//!!! 衝突の判定がよくエラーが出る箇所
	for (const auto &obj : objects)
	{
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
		{
			//!ここも重要：点と面の衝突をどのようにすれば矛盾なく判定できるか．
			// p->radius = radius / 2.5; // Mean(extLength(p->getLines()));
			auto toF = extXtuple(p->getFaces()) - p->getXtuple();
			auto toP = extXtuple(p->getNeighbors()) - p->getXtuple();
			double a = Norm(*std::min_element(toP.begin(), toP.end(), [](const auto &a, const auto &b)
											  { return Norm(a) < Norm(b); }));
			double b = Norm(*std::min_element(toF.begin(), toF.end(), [](const auto &a, const auto &b)
											  { return Norm(a) < Norm(b); }));
			// p->radius = Norm(Mean(toF);
			// p->radius = Mean(extLength(p->getLines())) / 2.;
			p->radius = Mean(extLength(extractLines(Flatten(BFS(p, 2))))) / 5.;
			// p->radius = radius / 2;
			p->addContactFaces(obj->getBucketFaces(), false); /**shadowあり*/
		}
	}
	// b% step2 面の境界条件を判定
	std::cout << "step2 面の境界条件を判定" << std::endl;
	/*面Aの点が接触している面Bを取得．A,B面が向き合っていればノイマン*/
	for (const auto &f : water.getFaces())
	{
		f->Neumann = isNeumann(f); //(!faces_p0.empty()) && (!faces_p1.empty()) && (!faces_p2.empty());
		f->Dirichlet = !f->Neumann;
	}
	// b% step3 線の境界条件を決定
	std::cout << "step3 線の境界条件を決定" << std::endl;
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
	std::cout << "step4 点の境界条件を決定" << std::endl;
	/*
	  周りの面が全てノイマンなら点はノイマン．
	  周りの面がディリクレとノイマン両方を持っていれば角点
	　それ以外は，ディリクレ．
	*/
	for (const auto &p : Points)
	{
		auto faces = p->getFacesUO();
		bool isNeumann = std::all_of(faces.begin(), faces.end(), [](const auto &f)
									 { return f->Neumann; }) &&
						 !p->getContactFacesX().empty();
		// getContactFacesXがないと，周りの面がノイマンでも，phinを計算できないことがある．
		bool isDirichlet = std::all_of(faces.begin(), faces.end(), [](const auto &f)
									   { return f->Dirichlet; });

		if (isNeumann)
			p->setN();
		else if (isDirichlet)
			p->setD();
		else
			p->setC();
	}
	// b! ------------------------------------------------------ */
	// b!       　    ノイマン境界の点や面にはΦnを与える               */
	// b! ------------------------------------------------------ */
	// b! 点
	std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << reset << std::endl;
	for (const auto &p : Points)
	{
		if (p->Neumann || p->CORNER)
		{
			std::get<1>(p->phiphin) = p->phin_Neumann = Dot(uNeumann(p), p->getNormalNeumann_BEM());
			if (!isFinite(p->phiphin))
			{
				std::cout << "p->phiphinはfiniteではない！！" << std::endl;
				if (p->Neumann)
					std::cout << "Neumann" << std::endl;
				if (p->Dirichlet)
					std::cout << "Dirichlet" << std::endl;
				if (p->CORNER)
					std::cout << "CORNER" << std::endl;
				std::cout << "forced_velocity(t) = " << forced_motion::velocity(real_time) << std::endl;
				std::cout << "p->phiphin = " << p->phiphin << std::endl;
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			}
		}
	}
	// b! 面
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

// b* ------------------------- 出力 ------------------------- */
VV_VarForOutput dataForOutput(const Network &water, const double dt)
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
	// if (*hist.data.rbegin() > 5)
	// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "length > 5");

	int ii = 0;
	derivatives ders(water, dt, true);
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
	//
	// auto P_aphiat = ders.P_aphiat;
	// auto P_aphiant = ders.P_aphiant;
	// auto P_pressure = ders.P_pressure;
	auto P_U_dot_gradgrad_U = ders.P_U_dot_gradgrad_U;
	auto P_kappa = ders.P_kappa;
	auto P_laplacian = ders.P_laplacian;

	uomap_P_d P_volume, P_phi_Neumann, P_phi_Dirichlet;
	uomap_P_Tddd P_accel_body, P_NearestContactFacesX, P_phin_Neumann, P_phin_Dirichlet;
	uomap_P_Tddd P_velocity_body;
	/* ------------------------------------------------------ */
	uomap_P_Vd P_lines_length, P_Intxn_length;
	uomap_P_d P_state, P_solidangleBIE, P_height, P_phi, P_phin, P_face_size, P_radius, P_lines_size, P_ishit;
	uomap_P_d P_BC, P_Intxn_size, P_ContactFaces, P_is_multiple_phiphin, P_min_depth;
	uomap_P_Tddd P_uNeumann, P_normal, P_normal_BEM, P_mirrorPosition, P_U_normal_BEM, P_U_tangential_BEM, P_U_BEM, P_U_update_BEM, P_U_cling_to_Neumann;
	uomap_P_d P_IG, P_IGn, P_isGoodForQuad, P_smin_min, P_s_m, P_minViewRatio;
	uomap_P_d P_isWellConditon, P_aphiat, P_aphiant, P_pressure, P_update_vs_cling, P_normalVariance;
	Print("ders.P_phiphin_InnerOuterCornerPを出力");
	try
	{
		for (const auto &p : water.getPoints())
		{
			if (p->Neumann || p->CORNER)
			{
				P_accel_body[p] = p->getNormal_BEM() * accel_normal_from_Neumann_surface(p); //!
				P_velocity_body[p] = local_velocity_from_Neumann_surface(p);				 //!
			}
			// auto [m0, s0, min0, smin0, max0] = distorsion(p, dt);
			// P_s_m[p] = s0 / m0;
			// P_smin_min[p] = smin0 / min0;
			P_volume[p] = water.getVolume();
			P_phin_Neumann[p] = p->getNormalNeumann_BEM() * p->phin_Neumann;
			P_phin_Dirichlet[p] = p->getNormalDirichlet_BEM() * p->phin_Dirichlet;
			P_phi_Neumann[p] = p->phi_Neumann;
			P_phi_Dirichlet[p] = p->phi_Dirichlet;
			if (p->Neumann || p->CORNER)
				if (std::get<0>(getNearestContactFacesX(p)))
				{
					auto tmp = std::get<1>(getNearestContactFacesX(p)) - p->getXtuple();
					if (Norm(tmp) < p->radius)
						P_NearestContactFacesX[p] = tmp;
				}
			P_U_cling_to_Neumann[p] = p->U_cling_to_Neumann;
			P_update_vs_cling[p] = Norm(p->U_cling_to_Neumann) / Norm(p->U_update_BEM);
			P_U_BEM[p] = p->U_BEM;
			P_U_update_BEM[p] = p->U_update_BEM;
			P_solidangleBIE[p] = p->getSolidAngle();
			P_minViewRatio[p] = minViewRatio(p);
			P_normalVariance[p] = normalVariance(p);
			P_U_normal_BEM[p] = p->U_normal_BEM;
			P_U_tangential_BEM[p] = p->U_tangential_BEM;
			// P_state[p] = p->getStatus();
			P_height[p] = p->getX()[2];
			P_phi[p] = std::get<0>(p->phiphin);
			P_phin[p] = std::get<1>(p->phiphin);
			// P_normal[p] = p->getNormalTuple();
			P_normal_BEM[p] = p->getNormal_BEM();
			// P_face_size[p] = (double)p->getFaces().size();
			// P_lines_size[p] = (double)p->getLines().size();
			// P_lines_length[p] = extLength(p->getLines());
			// P_Intxn_size[p] = (double)takeIntxn(p->getLines()).size();
			P_ContactFaces[p] = (double)p->getContactFaces().size();
			// P_Intxn_length[p] = extLength(takeIntxn(p->getLines()));
			P_BC[p] = p->Dirichlet ? 0. : (p->Neumann ? 1. : (p->CORNER ? 2. : 1 / 0.));
			// if (!p->getContactFaces().empty())
			// 	P_mirrorPosition[p] = 2. * ((*p->getContactFaces().begin()).second) - p->getXtuple();
			P_radius[p] = p->radius;
			// P_ishit[p] = (double)(!p->getContactFaces().empty());
			// P_ishit[p] = (double)(p->getStatus());
			// P_is_multiple_phiphin[p] = (double)(p->multiple_phiphin.size());
			// P_min_depth[p] = p->minDepthFromCORNER; //;getMinDepth(p);
			bool isgood = true;
			for (const auto &l : p->getLines())
				isgood = isgood && l->isGoodForQuadInterp();
			P_isGoodForQuad[p] = isgood;
			P_isWellConditon[p] = isWellConditon(p);
			P_aphiat[p] = std::get<0>(p->phiphin_t);
			P_aphiant[p] = std::get<1>(p->phiphin_t);
			P_pressure[p] = p->pressure_BEM;
			P_uNeumann[p] = uNeumann(p);
		}
		/* ------------------------------------------------------ */
		// Tddd XYZ = {1E+100, 1E+100, 1E+100};
		// networkPoint *originP;
		// for (const auto &p : water.getPoints())
		// {
		// 	double tmp = Norm(p->getXtuple());
		// 	if (Norm(XYZ) > tmp)
		// 	{
		// 		XYZ = p->getXtuple();
		// 		originP = p;
		// 	}
		// 	P_IG[p] = 0.;
		// 	P_IGn[p] = 0.;
		// }
		// for (const auto &f : water.getFaces())
		// 	for (const auto &[p, igign] : BEM::calc_P_IGIGnQuadTuple_mod(f, originP))
		// 	{
		// 		auto [p0, p1, p2] = f->getPointsTuple();
		// 		auto [IG, IGn] = igign;
		// 		P_IG[p0] += IG;
		// 		P_IG[p1] += IG;
		// 		P_IG[p2] += IG;
		// 		P_IGn[p0] += IGn;
		// 		P_IGn[p1] += IGn;
		// 		P_IGn[p2] += IGn;
		// 	}
		/* ------------------------------------------------------ */
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
			// {"s/m", P_s_m},
			// {"smin/min", P_smin_min},
			{"accel_body", P_accel_body},
			// {"volume", P_volume},
			{"U_cling_to_Neumann", P_U_cling_to_Neumann},
			{"absU_cling_to_Neumann/absU_update_BEM", P_update_vs_cling},
			{"velocity_body", P_velocity_body},
			{"NearestContactFacesX", P_NearestContactFacesX},
			{"φn_Neumann", P_phin_Neumann},
			{"φn_Dirichlet", P_phin_Dirichlet},
			{"φ_Neumann", P_phi_Neumann},
			{"φ_Dirichlet", P_phi_Dirichlet},
			{"min_depth", P_min_depth},
			// P_NearestContactFacesXで確かに最寄の構造物までをさすことができている．！
			// これで改善できない？
			{"U_update_BEM", P_U_update_BEM},
			{"U_BEM", P_U_BEM},
			{"U_normal_BEM", P_U_normal_BEM},
			{"U_tangential_BEM", P_U_tangential_BEM},
			{"kappa", P_kappa},
			{"ContactFaces", P_ContactFaces},
			// {"dxdt_mod", P_dxdt_mod},
			{"grad_phi", P_gradPhi},
			{"phin_vector", P_phin_vector},
			{"gradPhiTangential", P_gradPhi_tangential},
			{"height", P_height},
			{"φ", P_phi},
			{"φn", P_phin},
			{"solidangle", P_solidangleBIE},
			{"minViewRatio", P_minViewRatio},
			{"normalVariance", P_normalVariance},
			{"normal", P_normal},
			{"uNeumann", P_uNeumann},
			{"isWellConditon", P_isWellConditon},
			{"normal_BEM", P_normal_BEM},
			{"all_lines_are_GoodForQuad", P_isGoodForQuad},
			// {"laplacian_phi", P_laplacian},
			// {"face_size", P_face_size},
			// {"line_length", 10, P_lines_length},
			// {"line_size", P_lines_size},
			// {"Intxn_size", P_Intxn_size},
			// {"Intxn_length", 10, P_Intxn_length},
			{"boundary condition", P_BC},
			{"radius", P_radius},
			// {"state", P_state},
			{"tension", ders.P_tension},
			{"DφDt", P_DphiDt},
			{"φt", P_aphiat},
			{"φnt", P_aphiant},
			{"pressure", P_pressure},
			{"U.∇U", P_U_dot_gradgrad_U},
			{"IG", P_IG},
			{"IGn", P_IGn}
			// {"vector to mirrorPosition", P_mirrorPosition},
			// {"is hit", P_ishit}
		};
		return data;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
	std::cout << __PRETTY_FUNCTION__ << " done" << std::endl;
	// mk_vtu(output_directory + "/" + obj.getName() + std::to_string(time_step) + ".vtu", obj.getFaces(), datacpg);
	// mk_vtu(output_directory + "/" + name + std::to_string(time_step) + ".vtu", cpg.well_Faces, datacpg);
	// cpg_pvd.push(name, name + std::to_string(time_step) + ".vtu", time_step * t_rep * dt);
	// cpg_pvd.output_();
};
/* ------------------------------------------------------ */

double dt_CFL(const Network &water, double min_dt, const Tdd coeff = {0.4, 1})
{
	auto [c0, c1] = coeff;
	for (const auto &p : water.getPoints())
	{
		for (const auto &q : p->getNeighbors())
		{
			if (min_dt > c0 * Norm(p->getXtuple() - q->getXtuple()) / Norm(p->U_update_BEM - q->U_update_BEM))
			{
				min_dt = c0 * Norm(p->getXtuple() - q->getXtuple()) / Norm(p->U_update_BEM - q->U_update_BEM);
			}
			if (min_dt > c1 * Norm(p->getXtuple() - q->getXtuple()) / Norm(p->U_update_BEM))
			{
				min_dt = c1 * Norm(p->getXtuple() - q->getXtuple()) / Norm(p->U_update_BEM);
			}
		}
	}
	return min_dt;
};

void show_info(const Network &net)
{
	int total = 0, total_c_face = 0, c = 0, n = 0, d = 0;
	for (const auto &p : net.getPoints())
	{
		total++;
		if (p->CORNER)
		{
			c++;
			total_c_face += p->getFaces().size();
		}
		else if (p->Neumann)
			n++;
		else if (p->Dirichlet)
			d++;
	}
	std::cout << "Total : " << total << std::endl;
	int doublenode = total - c + total_c_face;
	std::cout << "Total case double-node : " << doublenode << std::endl;
	std::cout << "node reduction : " << (double)(doublenode - total) / (double)doublenode << std::endl;
	std::cout << "CORNER : " << c << std::endl;
	std::cout << "Total CORNER faces : " << total_c_face << std::endl;
	std::cout << "Neumann : " << n << std::endl;
	std::cout << "Dirichlet : " << d << std::endl;
};
JSONoutput jsonout;

// b! ------------------------------------------------------ */

struct calculateForceFromPressure
{
	Tddd force, torque;
	double area;
	T6d acceleration;
	std::vector<std::tuple<Tddd, T3Tddd>> PressureVeticies;
	calculateForceFromPressure(const std::unordered_set<networkFace *> faces, const Network *PasObj)
		: force({0., 0., 0.}), torque({0., 0., 0.}), area(0.), PressureVeticies({}), acceleration({0., 0., 0., 0., 0., 0.})
	{
		// PasObjと接したfaceの頂点にpressureが設定されている前提
		Tddd pressure012 = {1E+60, 1E+60, 1E+60};
		for (const auto &f : faces)
		{
			pressure012 = {1E+60, 1E+60, 1E+60};
			auto [p0, p1, p2] = f->getPointsTuple();
			for (const auto &f0 : p0->getContactFaces())
				if (f0->getNetwork() == PasObj)
					std::get<0>(pressure012) = p0->pressure;
			for (const auto &f1 : p1->getContactFaces())
				if (f1->getNetwork() == PasObj)
					std::get<1>(pressure012) = p1->pressure;
			for (const auto &f2 : p2->getContactFaces())
				if (f2->getNetwork() == PasObj)
					std::get<2>(pressure012) = p2->pressure;
			if (isFinite(pressure012))
				PressureVeticies.emplace_back(std::tuple<Tddd, T3Tddd>{pressure012, T3Tddd{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()}});
		}
	};

	Tddd getTorque(const Tddd &COM)
	{
		Tddd ret = {0., 0., 0.};
		for (const auto &[P012, X012] : PressureVeticies)
		{
			auto intpP = interpolationTriangleLinear0101(P012);
			auto intpX = interpolationTriangleLinear0101(X012);
			for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
			{
				auto P = TriangleNormal(X012) * intpP(x0, x1);
				auto X = intpX(x0, x1);
				ret += Cross(P, X - COM) * intpX.J(x0, x1) * w0w1;
			}
		}
		this->torque = ret;
		return ret;
	};

	Tddd getForce()
	{
		Tddd ret = {0., 0., 0.};
		for (const auto &[P012, X012] : PressureVeticies)
		{
			auto intpP = interpolationTriangleLinear0101(P012);
			auto intpX = interpolationTriangleLinear0101(X012);
			for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
				ret += TriangleNormal(X012) * intpP(x0, x1) * intpX.J(x0, x1) * w0w1;
		}
		this->force = ret;
		return ret;
	};

	T6d getAcceleration(const Tddd &COM, const T6d &I)
	{
		auto [mx, my, mz, Ix, Iy, Iz] = I;
		auto [Tx, Ty, Tz] = getTorque(COM);
		auto [Fx, Fy, Fz] = getForce();
		return {Fx / mx, Fy / my, Fz / mz, Tx / Ix, Ty / Iy, Tz / Iz};
	};
};

// b! ------------------------------------------------------ */

struct outputInfo
{
	std::string pvd_file_name;
	std::string vtu_file_name;
	PVDWriter *PVD;
	outputInfo(){};
};

// b! ------------------------------------------------------ */

int main()
{
	try
	{
		//* ------------------------------------------------------ */
		//*                         setting                        */
		//* ------------------------------------------------------ */
		JSON settingJSON(std::ifstream("./setting.json"));
		std::map<Network *, outputInfo> NetOutputInfo;
		std::vector<Network *> RigidBodyObject;
		std::vector<Network *> FluidObject;
		std::string output_directory = settingJSON["output_directory"][0];
		double max_dt = stod(settingJSON["max_dt"])[0];
		std::filesystem::create_directories(output_directory);
		/* ------------------------------------------------------ */
		double stop_remesh_time = 1E+10;
		double force_remesh_time = 0;
		double preparation_time = stod(settingJSON["preparation_time"])[0];
		double preparation_max_dt = stod(settingJSON["preparation_max_dt"])[0];
		if (settingJSON.find("stop_remesh_time"))
			stop_remesh_time = stod(settingJSON["stop_remesh_time"][0]);

		if (settingJSON.find("force_remesh_time"))
			force_remesh_time = stod(settingJSON["force_remesh_time"][0]);

		int grid_refinement = 0;
		if (settingJSON.find("grid_refinement"))
			grid_refinement = stoi(settingJSON["grid_refinement"][0]);
		/* ------------------------------------------------------ */

		for (const auto &FileName : settingJSON["inputfiles"])
		{
			JSON J(std::ifstream("./" + FileName));
			Print("---------------------------------------------------");
			Print("-----------------" + FileName + "------------------");
			Print("---------------------------------------------------");
			for (auto &[key, value] : J())
				std::cout << key << ": " << value << std::endl;
			if (!J.find("ignore") || !stob(J["ignore"])[0])
			{
				auto net = new Network(J["objfile"][0], J["name"][0]);
				for (auto &p : net->getPoints())
					p->radius = stod(J["radius"])[0];
				if (J.find("isMovementPredefined"))
					net->isMovementPredefined = stob(J["isMovementPredefined"])[0];
				modify(*net, J);
				NetOutputInfo[net].pvd_file_name = J["output_pvd_file_name"][0];
				NetOutputInfo[net].vtu_file_name = J["output_vtu_file_name"][0];
				NetOutputInfo[net].PVD = new PVDWriter(output_directory + "/" + J["output_pvd_file_name"][0] + ".pvd");
				std::filesystem::copy_file(FileName, output_directory + "/" + FileName, std::filesystem::copy_options::overwrite_existing);
				if (J["type"][0] == "RigidBody")
					RigidBodyObject.emplace_back(net);
				else if (J["type"][0] == "Fluid")
					FluidObject.emplace_back(net);
			}
			else
			{
				Print("skipped");
			}
			Print("---------------------------------------------------");
		}
		std::filesystem::copy_file("./main.cpp", output_directory + "/main.cpp", std::filesystem::copy_options::overwrite_existing);
		auto water = FluidObject[0];
		/* ------------------------------------------------------ */
		PVDWriter cornerPointsPVD(output_directory + "/cornerPointsPVD.pvd");
		PVDWriter cornerPVD(output_directory + "/corner.pvd");
		//  b@ ------------------------------------------------------ */
		//  b@                     グリッドの調整                         */
		//  b@ ------------------------------------------------------ */
		// if (true)
		// {
		// 	setBoundaryConditions(water, RigidBodyObject);
		// 	for (const auto &p : water.getPoints())
		// 	{
		// 		double h = 0.5;
		// 		if (!p->Neumann)
		// 		{
		// 			auto tmp = p->getXtuple();
		// 			p->setXSingle(Tddd{std::get<0>(tmp), std::get<1>(tmp), h});
		// 		}
		// 	}
		// 	water.setBounds();
		// }
		//
		{
			PVDWriter preRefinementPVD(output_directory + "/preRefinement.pvd");
			double dt = 1E-6;
			double min_dt = 1E+10;
			for (auto i = 0; i < grid_refinement; i++)
			{
				std::cout << "refinement i = " << i << std::endl;
				// b% -------------------------------------------------------- */
				// b%            境界条件（角点・ディリクレ・ノイマン）の決定               */
				// b% -------------------------------------------------------- */
				setBoundaryConditions(*water, RigidBodyObject);
				auto data = dataForOutput(*water, dt);
				auto filename = output_directory + "/preRefinement" + std::to_string(i) + ".vtu";
				mk_vtu(filename, {water->getFaces()}, data);
				preRefinementPVD.push(filename, (double)i);
				preRefinementPVD.output();
				std::string name = output_directory + "/water100_refined_" + std::to_string(i) + ".obj";
				// std::string name = output_directory + "/water30_refined_laplacian_0d05_" + std::to_string(i) + ".obj";
				std::ofstream ofs(name);
				std::cout << "name = " << name << std::endl;
				creteOBJ(ofs, *water);
				ofs.close();
				// b* ------------------------------------------------------ */
				// b*                    微分∇ΦやDUDtを計算                    */
				// b* ------------------------------------------------------ */
				std::cout << Green << "微分∇ΦやDUDtを計算" << reset << std::endl;
				int loop = 2;
				double rad = M_PI / 180;
				for (auto ii = 0; ii <= loop; ++ii)
				{
					min_dt = dt_CFL(*water, min_dt);
					if (min_dt < preparation_max_dt)
						dt = min_dt;
					else
						dt = preparation_max_dt;

					int RK_order = 3;
					std::map<netPp, RungeKutta_<Tddd> *> P_RK_X;
					auto Points = water->getPoints();
					for (const auto &p : Points)
						P_RK_X[p] = (new RungeKutta_(dt, real_time, p->getXtuple(), RK_order));

					do
					{
						derivatives ders(*water, dt, true);
						for (const auto &p : Points)
						{
							auto rk = P_RK_X[p];
							rk->push(ders.P_dxdt[p] * (i < 30 ? 0.1 : 1.));
							// rk->push(ders.P_dxdt[p]);
							p->setXSingle(rk->getX());
						}
						water->setBounds();
					} while (!(P_RK_X[*Points.begin()]->finished));

					for (const auto &p : water->getPoints())
						delete P_RK_X[p];

					flipIf(*water, {15 * rad, rad}, {10 * rad /*結構小さく*/, rad}, false);
					// if (time_step > 1 && time_step < 10)
					// {
					// 	flipIf(water, {8 * rad, 30 * rad}, {8 * rad /*結構小さく*/, 10 * rad}, true, 5 /*やはり制限はひつようか？*/);
					// 	flipIf(water, {8 * rad, 10 * rad}, {8 * rad /*結構小さく*/, 30 * rad}, true, 5 /*やはり制限はひつようか？*/);
					// }

					// flipIf(water, {rad, rad}, {10 * rad /*結構小さく*/, rad}, false);
					// if (i > 5)
					// 	if (ii == 0)
					// 	{
					// 		flipIf(water, {10 * rad, 30 * rad}, {10 * rad /*結構小さく*/, 10 * rad}, true, 10 /*やはり制限はひつようか？*/);
					// 		flipIf(water, {10 * rad, rad}, {10 * rad /*結構小さく*/, 30 * rad}, true, 10 /*やはり制限はひつようか？*/);
					// 		// arrangeCORNER(water, M_PI / 180, 10);
					// 		// 	arrangeCORNER(water, rad / 100., 10); //このアレンジで崩れることもあった．
					// 		// }else if (ii == 2)
					// 		// {
					// 		// 	flipIf(water, {rad, 30 * rad}, {5 * rad /*結構小さく*/, 5 * rad}, true, 2 /*やはり制限はひつようか？*/);
					// 		// 	flipIf(water, {rad, 5 * rad}, {5 * rad /*結構小さく*/, 30 * rad}, true, 2 /*やはり制限はひつようか？*/);
					// 	}
					// flipIf(water, {10 * rad, rad}, {10 * rad /*結構小さく*/, rad}, false);
				}
				// if (i % 2 == 0)
				// 	remesh(water, {10 * rad, 10 * rad}, {10 * rad /*結構小さく*/, 10 * rad}, false);
			};
			// 内部の角点に常にいるようにするには？？？？？？？？？？？？？・
			// 面が吸い付く必要がある
			std::ofstream ofs(output_directory + "/water_refined.obj");
			creteOBJ(ofs, *water);
			ofs.close();
		}
		//  b* ------------------------------------------------------ */
		//  b*                         メインループ                     */
		//  b* ------------------------------------------------------ */
		TimeWatch watch;
		for (time_step = 0; time_step < 2000; time_step++)
		{
			show_info(*water);
			//!体積を保存するようにリメッシュする必要があるだろう．
			// auto radius = Mean(extLength(water->getLines()));
			setBoundaryConditions(*water, RigidBodyObject);
			// b* ------------------------------------------------------ */
			if (real_time < stop_remesh_time && time_step > 0)
			{
				double rad = M_PI / 180;
				flipIf(*water, {15 * rad, rad}, {15 * rad /*結構小さく*/, rad}, false);
				if (time_step > 1 && time_step < 10)
				{
					flipIf(*water, {10 * rad, 30 * rad}, {10 * rad /*結構小さく*/, 20 * rad}, true, 20 /*やはり制限はひつようか？*/);
					flipIf(*water, {10 * rad, 10 * rad}, {10 * rad /*結構小さく*/, 30 * rad}, true, 20 /*やはり制限はひつようか？*/);
				}
			}
			// b# ------------------------------------------------------ */
			// b#                       刻み時間の決定                     */
			// b# ------------------------------------------------------ */
			const auto Points = water->getPoints();
			const auto Faces = water->getFaces();

			// double min_dt = 1E+10;
			// double dt = dt_CFL(*water, max_dt, {0.2, 1});
			double dt = max_dt;
			// if (2. < real_time && real_time < 4.)
			// 	dt = 0.005;

			// if (real_time <= preparation_time)
			// 	dt = preparation_max_dt;
			// else if (min_dt < max_dt)
			// 	dt = min_dt;
			// else
			// 	dt = max_dt;
			Print("===========================================================================");
			Print("       dt :" + Red + std::to_string(dt) + reset);
			Print("time_step :" + Red + std::to_string(time_step) + reset);
			Print("real time :" + Red + std::to_string(real_time) + reset);
			Print("---------------------------------------------------------------------------");

			// b* ------------------------- 出力 ------------------------- */
			auto data = dataForOutput(*water, dt);
			//流体
			for (const auto &FileName : FluidObject)
			{
				auto filename = NetOutputInfo[FileName].vtu_file_name + std::to_string(time_step) + ".vtu";
				mk_vtu(output_directory + "/" + filename, FileName->getFaces(), data);
				NetOutputInfo[FileName].PVD->push(filename, real_time);
				NetOutputInfo[FileName].PVD->output();
			}
			for (const auto &FileName : RigidBodyObject)
			{
				auto filename = NetOutputInfo[FileName].vtu_file_name + std::to_string(time_step) + ".vtu";
				mk_vtu(output_directory + "/" + filename, FileName->getFaces());
				NetOutputInfo[FileName].PVD->push(filename, real_time);
				NetOutputInfo[FileName].PVD->output();
			}
			//流体
			{
				std::vector<Tddd> points;
				for (const auto &p : water->getPoints())
					if (p->CORNER)
						points.emplace_back(p->getXtuple());

				auto filename = "cornerPoints" + std::to_string(time_step) + ".vtu";
				mk_vtu(output_directory + "/" + filename, points);
				cornerPointsPVD.push(filename, real_time);
				cornerPointsPVD.output();
			}
			{
				std::unordered_set<networkFace *> faces;
				for (const auto &f : water->getFaces())
					for (const auto p : f->getPoints())
						if (p->CORNER)
						{
							faces.emplace(f);
							break;
						}
				auto filename = "corner" + std::to_string(time_step) + ".vtu";
				mk_vtu(output_directory + "/" + filename, faces, data);
				cornerPVD.push(filename, real_time);
				cornerPVD.output();
			}
			// mk_vtu(output_directory + "/" + obj.getName() + std::to_string(time_step) + ".vtu", obj.getFaces(), datacpg);
			// mk_vtu(output_directory + "/" + name + std::to_string(time_step) + ".vtu", cpg.well_Faces, datacpg);
			// cpg_pvd.push(name, name + std::to_string(time_step) + ".vtu", time_step * t_rep * dt);
			// cpg_pvd.output_();
			// b* ------------------------------------------------------ */
			jsonout.push("time", real_time);
			jsonout.push("volume", water->getVolume());

			double maxUcling = 0;
			for (const auto &p : water->getPoints())
				if (p->Neumann || p->CORNER)
				{
					auto tmp = Norm(p->U_cling_to_Neumann) / Norm(p->U_update_BEM);
					if (maxUcling < tmp && isFinite(tmp))
						maxUcling = tmp;
				}
			jsonout.push("maxUcling", maxUcling);
			double maxUcling_minlength = 0;
			for (const auto &p : water->getPoints())
				if (p->Neumann || p->CORNER)
				{
					double min = Min(extLength(p->getLines()));
					auto tmp = Norm(p->U_cling_to_Neumann) * dt / min;
					if (maxUcling_minlength < tmp && isFinite(tmp))
						maxUcling_minlength = tmp;
				}
			jsonout.push("maxUcling_minlength", maxUcling_minlength);

			double maxUcling_meanlength = 0;
			for (const auto &p : water->getPoints())
				if (p->Neumann || p->CORNER)
				{
					double mean = Mean(extLength(p->getLines()));
					auto tmp = Norm(p->U_cling_to_Neumann) * dt / mean;
					if (maxUcling_meanlength < tmp && isFinite(tmp))
						maxUcling_meanlength = tmp;
				}
			jsonout.push("maxUcling_meanlength", maxUcling_meanlength);

			std::ofstream os(output_directory + "/result.json");
			jsonout.output(os);
			os.close();
			/* ------------------------------------------------------ */
			int depthlimit = 9;

			double spacing = Mean(extLength(water->getLines())) * 3;
			Buckets<networkFace> FMM_BucketsFaces(water->bounds(), spacing);
			Buckets<networkPoint> FMM_BucketsPoints(water->bounds(), spacing);
			for (const auto &f : water->getFaces())
				FMM_BucketsFaces.add(f->getXtuple(), f);
			for (const auto &p : water->getPoints())
				FMM_BucketsPoints.add(p->getXtuple(), p);

			int RK_order = 4;

			// b@ ------------------------------------------------------ */
			// b@        初期値問題を解く（時間微分方程式を数値積分する）           */
			// b@ ------------------------------------------------------ */
			std::map<netPp, RungeKutta_<double> *> P_RK_phi;
			std::map<netPp, RungeKutta_<Tddd> *> P_RK_X;
			std::map<Network *, RungeKutta_<Tddd> *> Net_RK_COM;
			std::map<Network *, RungeKutta_<T4d> *> Net_RK_Q;
			// 面ベースで境界条件を判断する
			/* ------------------------------------------------------ */
			// for (const auto& p : Points)
			// 	p->multiple_phiphin.clear();
			for (const auto &p : Points)
			{
				P_RK_phi[p] = (new RungeKutta_(dt, real_time, std::get<0>(p->phiphin), RK_order));
				P_RK_X[p] = (new RungeKutta_(dt, real_time, p->getXtuple(), RK_order));
			}

			// for (const auto &s : RigidBodyObject)
			// 	for (const auto &p : s->getPoints())
			// 		P_RK_X_body[p] = (new RungeKutta_(dt, real_time, p->getXtuple(), RK_order));
			for (const auto &net : RigidBodyObject)
			{
				Net_RK_COM[net] = (new RungeKutta_(dt, real_time, net->COM, RK_order));
				Net_RK_Q[net] = (new RungeKutta_(dt, real_time, net->Q(), RK_order));
			}

			int RK_step = 1;
			do
			{
				//! 壁面の動きは，マイステップ更新することにした．この結果はphin()で参照される
				auto RK_time = P_RK_phi[*Points.begin()]->gett(); //%各ルンゲクッタの時刻を使う
				std::cout << "RK_step = " << RK_step++ << "/" << RK_order << ", RK_time = " << RK_time << ", real_time = " << real_time << std::endl;
				// * ------------------------------------------------------ */
				// *                　物体のノイマン境界の設置               　　*/
				// * ------------------------------------------------------ */
				for (const auto &obj : RigidBodyObject)
					if (obj->isMovementPredefined)
					{
						obj->velocity = forced_motion::velocity(RK_time);		  // T6d //@ Φnを計算するために，物体表面の速度forced_velocityは，保存しておく必要がある
						obj->acceleration = forced_motion::acceleration(RK_time); // T6d //@ 圧力を計算するために，物体表面の加速度は，保存しておく必要がある
#ifndef experiment_Retzler2000
						obj->COM = forced_motion::displacement(RK_time);
						obj->RigidBodyMovePoints();
#endif
					}
				// // b# ------------------------------------------------------ */
				// // b#                　ノイマン境界にフィットさせる               */
				// // b# ------------------------------------------------------ */
				// for (const auto &p : Points)
				// 	p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};
				// calculateVectorToClingNeumann(*water, /*dt=*/0); //@必須
				// for (const auto &p : Points)
				// 	if (p->CORNER || p->Neumann)
				// 		p->setXSingle(p->getXtuple() + p->U_BUFFER);
				// water->setBounds();
				// b% -------------------------------------------------------- */
				// b%            境界条件（角点・ディリクレ・ノイマン）の決定          */
				// b% -------------------------------------------------------- */
				setBoundaryConditions(*water, RigidBodyObject);
				std::cout << Blue << "Elapsed time: " << Red << watch() << reset << " s\n";
				// b! ------------------------------------------------------ */
				// b!           　境界値問題を解く-> {Φ,Φn}が決まる                */
				// b! ------------------------------------------------------ */
				std::cout << Green << "境界値問題を解く-> {Φ,Φn}が決まる" << reset << std::endl;
				BEM_BVP BVP;
				BVP.solve(*water, FMM_BucketsPoints, FMM_BucketsFaces);
				std::cout << Blue << "Elapsed time: " << Red << watch() << reset << " s\n";
				// b* ------------------------------------------------------ */
				// b*                    微分∇ΦやDUDtを計算                    */
				// b* ------------------------------------------------------ */
				std::cout << Green << "微分∇ΦやDUDtを計算" << reset << std::endl;
				derivatives ders(*water, P_RK_phi[*Points.begin()]->getdt(), true);
				std::cout << Blue << "Elapsed time: " << Red << watch() << reset << " s\n";
				// b* ------------------------------------------------------ */
				// b*           　境界値問題を解く-> {Φt,Φtn}が決まる              */
				// b* ------------------------------------------------------ */
				std::cout << Green << "境界値問題を解く-> {Φt,Φtn}が決まる" << reset << std::endl;
				BVP.solveForPhiPhin_t();
				std::cout << Blue << "Elapsed time: " << Red << watch() << reset << " s\n";
				for (const auto obj : RigidBodyObject)
				{
					auto tmp = calculateForceFromPressure(water->getFaces(), obj);
					std::cout << red << "force = " << tmp.force << ", area = " << tmp.area << reset << std::endl;
				}
				// b* ------------------------------------------------------ */
				// b*                 ディリクレ境界ではΦを時間積分                 */
				// b* ------------------------------------------------------ */
				std::cout << Green << "ディリクレ境界ではΦを時間積分，ノイマン境界ではΦnを陽に与える" << reset << std::endl;
				for (const auto &p : Points)
					if (p->Dirichlet || p->CORNER)
					{
						// b! ここでノイマン面が変な場所に移動させられてはたまらない
						//@ Φの時間発展，Φnの時間発展はない
						{
							auto rk = P_RK_phi[p];
							rk->push(ders.P_DphiDt[p]);
							if (p->CORNER || p->Dirichlet)
								std::get<0>(p->phiphin) = p->phi_Dirichlet = rk->getX(); // 角点の法線方向はわからないので，ノイマンの境界条件phinを与えることができない．
						}
						//@ 位置xの時間発展
						{
							auto rk = P_RK_X[p];
							rk->push(ders.P_dxdt[p]);
							p->setXSingle(rk->getX());
						}
					}

#ifdef experiment_Retzler2000
				for (const auto &net : RigidBodyObject)
				{
					{
						auto rk = Net_RK_COM[net];
						rk->push(net->velocityTranslational());
						net->COM = rk->getX();
					}
					{
						auto rk = Net_RK_Q[net];
						Quaternion q;
						q = q.d_dt(net->velocityRotational()); // w->クォータニオン
						rk->push(q());						   //クォータニオン->T4dとしてプッシュ
						net->Q = rk->getX();
					}
					net->RigidBodyMovePoints();
				}
#endif
				/* ------------------------------------------------------ */
				std::cout << Green << "setBounds" << reset << std::endl;
				water->setBounds();
				std::cout << Green << "real_timeを取得" << reset << std::endl;
				real_time = P_RK_phi[*Points.begin()]->gett();
				std::cout << Blue << "Elapsed time: " << Red << watch() << reset << " s\n";
			} while (!(P_RK_phi[*Points.begin()]->finished));

			for (const auto &p : Points)
			{
				delete P_RK_phi[p];
				delete P_RK_X[p];
			}

			for (const auto &s : RigidBodyObject)
			{
				delete Net_RK_COM[s];
				delete Net_RK_Q[s];
			}
			//

			for (const auto &p : Points)
			{
				p->U_BEM_last = p->U_BEM;
				p->U_tangential_BEM_last = p->U_tangential_BEM;
			}

			/* ------------------------------------------------------ */
			/* ------------------------------------------------------ */
			// for (const auto &p : Points)
			// 	p->X_BUFFER = p->getXtuple();

			// calculateVectorToSurfaceInBuffer(*water, 0.);

			// for (const auto &p : Points)
			// 	p->calculateBufferPotentialsOnClungSurface();

			// for (const auto &p : Points)
			// 	p->setX(p->getXtuple() + p->U_BUFFER);

			// for (const auto &p : Points)
			// 	p->copyPotentialsBuffer();

			// water->setBounds();
			/* ------------------------------------------------------ */
			/* ------------------------------------------------------ */

			// normalizePhi(wavemaker->getPoints());
			//
			std::ofstream ofs(output_directory + "/water_current.obj");
			creteOBJ(ofs, *water);
			ofs.close();
		}
	}
	catch (error_message &e)
	{
		e.print();
	};
	return 0;
};
