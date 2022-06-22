#ifndef NetworkUtility_H
#define NetworkUtility_H
#pragma once

#include "Network.hpp"
#include "InterpolationRBF.hpp"

// using T_2P = std::tuple<networkPoint *, networkPoint *>;
// using T_3P = std::tuple<networkPoint *, networkPoint *, networkPoint *>;
// using T_4P = std::tuple<networkPoint *, networkPoint *, networkPoint *, networkPoint *>;
// using T_5P = std::tuple<networkPoint *, networkPoint *, networkPoint *, networkPoint *, networkPoint *>;
// using T_6P = std::tuple<networkPoint *, networkPoint *, networkPoint *, networkPoint *, networkPoint *, networkPoint *>;

// T2Tddd ToX(const T_2P &ps) { return {std::get<0>(ps)->getXtuple(), std::get<1>(ps)->getXtuple()}; };
// T3Tddd ToX(const T_3P &ps) { return {std::get<0>(ps)->getXtuple(), std::get<1>(ps)->getXtuple(), std::get<2>(ps)->getXtuple()}; };
// T4Tddd ToX(const T_4P &ps) { return {std::get<0>(ps)->getXtuple(), std::get<1>(ps)->getXtuple(), std::get<2>(ps)->getXtuple(), std::get<3>(ps)->getXtuple()}; };
// T5Tddd ToX(const T_5P &ps) { return {std::get<0>(ps)->getXtuple(), std::get<1>(ps)->getXtuple(), std::get<2>(ps)->getXtuple(), std::get<3>(ps)->getXtuple(), std::get<4>(ps)->getXtuple()}; };
// T6Tddd ToX(const T_6P &ps) { return {std::get<0>(ps)->getXtuple(), std::get<1>(ps)->getXtuple(), std::get<2>(ps)->getXtuple(), std::get<3>(ps)->getXtuple(), std::get<4>(ps)->getXtuple(), std::get<5>(ps)->getXtuple()}; };

/* ------------------------------------------------------ */
void creteOBJ(std::ofstream &ofs, Network &net)
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
//% ------------------------------------------------------ */
//%                         反射点の計算                     */
//% ------------------------------------------------------ */
Tddd getClosestFacePoint(const networkPoint *p,
						 const std::unordered_set<networkFace *> &boundary_face,
						 const double mirroring_distance)
{
	Tddd ret = {1E+40, 1E+40, 1E+40}; //反射点を返す
	auto sphere = geometry::Sphere(p->getXtuple(), mirroring_distance);
	for (const auto &f : boundary_face)
	{
		auto ITX = intersection(sphere, geometry::Triangle(f->getXVertices()));
		if (ITX.isIntersecting)
			if (Norm(ret) > Norm(ITX.X - p->getXtuple()))
				ret = ITX.X - p->getXtuple();
	}
	return ret;
};
std::unordered_map<networkPoint *, Tddd> getMap_PointToClosestFacePoint(const std::unordered_set<networkPoint *> &dummy_points,
																		const std::unordered_set<networkFace *> &boundary_face,
																		const double mirroring_distance)
{
	std::unordered_map<networkPoint *, Tddd> ret; //反射点を返す
	for (const auto &p : dummy_points)
	{
		Tddd pToX = getClosestFacePoint(p, boundary_face, mirroring_distance);
		if (Norm(pToX) < 1E+20)
		{
			auto it = ret.find(p);
			if (it == ret.end())
				ret[p] = pToX;
			else if (Norm(it->second) > Norm(pToX))
				ret[p] = pToX;
		}
	}
	return ret;
};
/* ------------------------------------------------------ */
class polarInterpolation_ : public InterpolationRBF_<Tdd, Tddd>
{
public:
	polarInterpolation_(networkPoint *const origin, networkLine *base_line = nullptr)
		: InterpolationRBF_<Tdd, Tddd>()
	{
		try
		{
			double scale = Mean(extLength(origin->getLines()));
			std::vector<std::tuple<Tdd, double>> kernel_X_s_IN;
			std::vector<std::tuple<Tdd, Tddd>> EQs_IN;
			Tddd X_;
			Tdd t0t1;
			for (const auto &tup : origin->getNeighbors_Depth2_OnPolarAsTuple(base_line))
			{
				auto [t0, t1, X, p] = tup;
				X_ = Tddd{X[0], X[1], X[2]};
				t0t1 = {t0, t1};
				kernel_X_s_IN.emplace_back(std::tuple<Tdd, double>{t0t1, scale});
				EQs_IN.emplace_back(std::tuple<Tdd, Tddd>{t0t1, X_});
			}
			t0t1 = Tdd{0, 0};
			kernel_X_s_IN.emplace_back(std::tuple<Tdd, double>{t0t1, scale /*パラメトリックなので*/});
			EQs_IN.emplace_back(std::tuple<Tdd, Tddd>{t0t1, origin->getXtuple()});

			this->initialize(kernel_X_s_IN, EQs_IN);
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
};

class parametricPolarInterpolation_ : public InterpolationRBF_<Tdd, Tddd>
{
public:
	parametricPolarInterpolation_(networkPoint *const origin, networkLine *base_line = nullptr)
		: InterpolationRBF_<Tdd, Tddd>()
	{
		try
		{
			double scale = 1.25;
			std::vector<std::tuple<Tdd, double>> kernel_X_s_IN;
			std::vector<std::tuple<Tdd, Tddd>> EQs_IN;
			Tddd X_;
			Tdd t0t1;
			for (const auto &tup : origin->getNeighbors_Depth2_OnPolarAsTuple_parametric(base_line))
			{
				auto [t0, t1, X, p] = tup;
				X_ = {X[0], X[1], X[2]};
				t0t1 = {t0, t1};
				kernel_X_s_IN.emplace_back(std::tuple<Tdd, double>{t0t1, scale /*パラメトリックなので*/});
				EQs_IN.emplace_back(std::tuple<Tdd, Tddd>{t0t1, X_});
			}
			t0t1 = {0, 0};
			kernel_X_s_IN.emplace_back(std::tuple<Tdd, double>{t0t1, scale /*パラメトリックなので*/});
			EQs_IN.emplace_back(std::tuple<Tdd, Tddd>{t0t1, origin->getXtuple()});
			this->initialize(kernel_X_s_IN, EQs_IN);
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
};

class parametricPolarInterpolation2_ : public InterpolationRBF_<Tdd, Tddd>
{
public:
	parametricPolarInterpolation2_(networkPoint *const origin, networkLine *base_line = nullptr)
		: InterpolationRBF_<Tdd, Tddd>()
	{
		try
		{
			double scale = 1.25;
			std::vector<std::tuple<Tdd, double>> kernel_X_s_IN;
			std::vector<std::tuple<Tdd, Tddd>> EQs_IN;
			Tddd X_;
			Tdd t0t1;
			for (const auto &tup : origin->getNeighbors_Depth2_OnPolarAsTuple_parametric(base_line))
			{
				auto [t0, t1, X, p] = tup;
				X_ = {X[0], X[1], X[2]};
				t0t1 = {t0, t1};
				kernel_X_s_IN.emplace_back(std::tuple<Tdd, double>{t0t1, scale /*パラメトリックなので*/});
				EQs_IN.emplace_back(std::tuple<Tdd, Tddd>{t0t1, X_ - origin->getXtuple()});
			}
			t0t1 = {0, 0};
			kernel_X_s_IN.emplace_back(std::tuple<Tdd, double>{t0t1, scale /*パラメトリックなので*/});
			EQs_IN.emplace_back(std::tuple<Tdd, Tddd>{t0t1, {0, 0, 0}});
			this->initialize(kernel_X_s_IN, EQs_IN);
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
};

/* ------------------------------------------------------ */
class polarInterpolation : public InterpolationVectorRBF
{
	/*
	 * RBFを使って，平面化したパラメタ空間で，実際の3次元の座標を補間する．
	 * 補間した3次元の値を取得する際に，平面パラメタを与える必要がでてくる．
	 * 三角形の線形補間パラメタt0,t1と各頂点で作られたRBF補間関数の平面パラメタXYの関係を計算する．
	 * 1) 面の周りの点と線を順番をそろえて取得する
	 * 2) 次に，各頂点において，RBF補間関数をそれぞれ作成する．合計３つ．
	 * 3) t0,t1を与えれば，平面パラメタは，それぞれの補間関数に対して計算することができ，3つのRBF補間関数を利用できる．
	 TODO: 微分を計算できるようにする．
	 */
public:
	polarInterpolation(networkPoint *const origin, networkLine *base_line = nullptr) : InterpolationVectorRBF()
	{
		try
		{
			VV_d xyz;
			VV_d param;
			V_netPp points;
			// for (const auto &tup : origin->getNeighbors_Depth2_OnPolarAsTuple(base_line))
			for (const auto &tup : origin->getNeighbors_Depth2_OnPolarAsTuple(base_line))
			{
				auto [t0, t1, X, p] = tup;
				points.emplace_back(p);
				param.emplace_back(V_d{t0, t1});
				xyz.emplace_back(X);
			}
			points.emplace_back(origin);
			param.emplace_back(V_d{0, 0});
			xyz.emplace_back(origin->getX());
			this->set(param, xyz);
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
};
/* ------------------------------------------------------ */
class polarInterpolation2
{
	/*
	 * 頂点に対して補間関数を作るpolarInterpolationを利用し，面に対する補間関数を作成する．
	 */
	polarInterpolation *intp0;
	polarInterpolation *intp1;
	polarInterpolation *intp2;
	std::tuple<networkPoint *, networkPoint *, networkPoint *> P;
	std::tuple<networkLine *, networkLine *, networkLine *> L;

public:
	~polarInterpolation2()
	{
		delete this->intp0;
		delete this->intp1;
		delete this->intp2;
	};
	polarInterpolation2(networkFace *const f)
		: intp0(nullptr), intp1(nullptr), intp2(nullptr)
	{
		/*
			  p2
			  /\
		 l2  /  \ l1
			/    \
		   --------
		p0    l0    p1
		*/
		try
		{
			this->P = f->getPointsTuple();
			this->L = f->getLinesTuple();
			this->intp0 = new polarInterpolation(std::get<0>(P), std::get<0>(L));
			this->intp1 = new polarInterpolation(std::get<1>(P), std::get<1>(L));
			this->intp2 = new polarInterpolation(std::get<2>(P), std::get<2>(L));
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	/* ------------------------------------------------------ */
	Tddd operator()(double t0, double t1)
	{
		//面上のある点
		Tddd shape = {t0, t1 * (1. - t0), (-1. + t0) * (-1. + t1)};
		Tddd X = Dot(shape, T3Tddd{std::get<0>(P)->getXtuple(), std::get<1>(P)->getXtuple(), std::get<2>(P)->getXtuple()});
		//この点が各RBFのパラメタにおいて何になるか計算
		auto A = std::get<0>(P)->getXtuple();
		auto B = std::get<1>(P)->getXtuple();
		auto C = std::get<2>(P)->getXtuple();
		//
		auto q0 = MyVectorAngle(B - A, X - A) * (2. * M_PI / Total(std::get<0>(P)->getAngles()));
		auto q1 = MyVectorAngle(C - B, X - B) * (2. * M_PI / Total(std::get<1>(P)->getAngles()));
		auto q2 = MyVectorAngle(A - C, X - C) * (2. * M_PI / Total(std::get<2>(P)->getAngles()));
		auto r0 = Norm(X - A);
		auto r1 = Norm(X - B);
		auto r2 = Norm(X - C);
		if (!isFinite(q0))
			q0 = 0;
		if (!isFinite(q1))
			q1 = 0;
		if (!isFinite(q2))
			q2 = 0;
		auto tmp0 = (*(this->intp0))({r0 * cos(q0), r0 * sin(q0)});
		auto tmp1 = (*(this->intp1))({r1 * cos(q1), r1 * sin(q1)});
		auto tmp2 = (*(this->intp2))({r2 * cos(q2), r2 * sin(q2)});
		return Dot(shape, T3Tddd{{tmp0[0], tmp0[1], tmp0[2]}, {tmp1[0], tmp1[1], tmp1[2]}, {tmp2[0], tmp2[1], tmp2[2]}});
		// auto ret = (*(this->intp0))({r0, q0});
		// auto ret = (*(this->intp1))({r1*cos(q1), r1*sin(q1)});
		// return {ret[0], ret[1], ret[2]};
	};
	/* ------------------------------------------------------ */
	double Sqrt(const double x) const { return std::sqrt(x); };
	double Power(const double x, const double n) const { return std::pow(x, n); };
	T3Tdd Xi(const Tddd &P0, const Tddd &P1, const Tddd &P2, const double t0, const double t1)
	{
		// P0のパラメタ-Xi0,Xi1を返す．
		Tddd shape = {t0, t1 * (1. - t0), (-1. + t0) * (-1. + t1)};
		Tddd X = Dot(shape, T3Tddd{P0, P1, P2});
		auto q0 = MyVectorAngle(P1 - P0, X - P0);
		auto q1 = MyVectorAngle(P2 - P1, X - P1);
		auto q2 = MyVectorAngle(P0 - P2, X - P2);
		auto r0 = Norm(X - P0);
		auto r1 = Norm(X - P1);
		auto r2 = Norm(X - P2);
		if (!isFinite(q0))
			q0 = 0;
		if (!isFinite(q1))
			q1 = 0;
		if (!isFinite(q2))
			q2 = 0;
		return {{r0 * cos(q0), r0 * sin(q0)}, {r1 * cos(q1), r1 * sin(q1)}, {r2 * cos(q2), r2 * sin(q2)}};
	};
	// T3Tdd DXiDt1(const Tddd &P0, const Tddd &P1, const Tddd &P2, const double t0, const double t1){

	// };
	// Tddd grad_Xi0Xi1()
	// {
	// 	auto X0 = (-((P1x * P2x + P1y * P2y + P1z * P2z) * (-1 + t0)) + (Power(P1x, 2) + Power(P1y, 2) - P1x * P2x - P1y * P2y + P1z * (P1z - P2z)) * t1) / Sqrt(Power(P1x, 2) + Power(P1y, 2) + Power(P1z, 2));
	// 	return {X0, std::sqrt(1. - X0 * X0)};
	// };
	// //ToDO
	// V_d N(const double t0, const double t1) const
	// {
	// 	for (auto i = 0; i < P.size(); i++)
	// 		PHI[i] = phi(x, P[i]);
	// 	return Dot(PHI, this->invF);
	// };

	// VV_d gradN(const V_d &x) const
	// {
	// 	if (this->parameter_dim != x.size())
	// 	{
	// 		std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
	// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
	// 	}
	// 	VV_d tmp(P.size());
	// 	for (auto i = 0; i < P.size(); i++)
	// 		for (const auto &gp : grad_phi(x, P[i]))
	// 			tmp[i].emplace_back(gp);
	// 	return Dot(Transpose(tmp), this->invF);
	// };

	// VV_d gradN(const double t0, const double t1) const
	// {
	// 	return this->gradN({t0, t1});
	// };

	// V_d cross(const double t0, const double t1) const
	// {
	// 	auto dxdt = Dot(this->gradN(t0, t1), this->VV_values);
	// 	return Cross(dxdt[0], dxdt[1]);
	// };

	// double J(const double t0, const double t1) const
	// {
	// 	auto dxdt = Dot(this->gradN(t0, t1), this->VV_values);
	// 	return std::sqrt(Dot(dxdt[0], dxdt[0]) * Dot(dxdt[1], dxdt[1]) - pow(Dot(dxdt[0], dxdt[1]), 2.));
	// };
};
/* ------------------------------------------------------ */
void isConsistent(const std::vector<Network *> &all_obj)
{
	for (const auto &net : all_obj)
	{
		for (auto l : net->getLines())
			isConsistent(l);
	}
};
///////////////////////////////////////////////////////////
// bool isConvexPolygon(const V_netPp &ps, const V_d &normal) { return geometry::isConvexPolygon(obj3D::extractX(ps), normal); };
// bool isConcavePolygon(const V_netPp &ps, const V_d &normal) { return geometry::isConcavePolygon(obj3D::extractX(ps), normal); };
///////////////////////////////////////////////////////////
bool isInConvexPolygon(const networkPoint *const p)
{
	try
	{
		auto ps = p->getNeighborsSort(); //ソートできない場合エラーとなる
		if (ps.size() < 2)
			return false;
		if (ps.size() == 3)
			return true;
		if (ps.size() != p->getNeighbors().size())
			return false; //面が重複する場合を覗くために設置

		return geometry::isConvexPolygon(extX(ps), p->getNormalTuple());
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
///////////////////////////////////////////////////////////
bool isStraight(const V_d &v0, const V_d &v1, const double angle = 1E-2 /*閾値*/)
{
	//        ^
	//      v1|angle (positive)
	// ---v0-->----
	//      v1|angle (negative)
	//        v
	if (std::abs(MyVectorAngle(v0, v1)) < angle)
		return true;
	return false;
};
///////////////////////////////////////////////////////////
bool isStraight(const netPp p0, const netPp p1, const netPp p2, const double angle = 1E-2 /*閾値*/)
{
	return (isStraight(p0->getX() - p1->getX(), p1->getX() - p2->getX()) < angle);
};
///////////////////////////////////////////////////////////
bool isFlat(const netPp p, double minangle = M_PI / 180.)
{
	if (p->getLines().empty())
	{
		// mk_vtu("./vtu/p->getLines().empty().vtu", {{p}});
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty()");
	}

	std::vector<Tddd> normals;
	for (const auto &f : p->getFaces())
		normals.emplace_back(f->getNormalTuple());
	//面の法線方向が作る内角が全て1E-2より小さい場合，flat
	for (auto i = 0; i < normals.size(); i++)
		for (auto j = i + 1; j < normals.size(); j++)
			if (Dot(normals[i], normals[j]) < cos(minangle))
				return false; // not flat
	return true;
};
///////////
bool isFlat(const netLp line, double minangle = M_PI / 180.)
{
	auto fs = line->getFaces();
	if (fs.size() != 2)
		return false;
	else if (Dot(fs[0]->getNormalTuple(), fs[1]->getNormalTuple()) < cos(minangle)) // dotはflatなら大きくなる
		return false;
	else
		return true;
};
///////////////////////////////////////////////////////////
Tddd AreaWeightedSmoothingVector(netPp p)
{
	try
	{
		if (!isEdgePoint(p))
		{
			Tddd ret = {0., 0., 0.}, p_next = p->getXtuple(), V = {0., 0, 0.};
			T3Tddd F;
			double scale = Mean(extLength(p->getLines()));
			double Wtot, W;
			for (auto k = 0; k < 1000; ++k)
			{
				Wtot = 0.;
				V *= 0.;
				for (const auto &f : p->getFaces())
				{
					auto [p0, p1, p2] = f->getPointsTuple(p);
					F = {p_next, p1->getXtuple(), p2->getXtuple()};
					W = TriangleArea(F);
					Wtot += W;
					V += scale * W * Normalize(Mean(F) - p_next);
				}
				V /= Wtot;
				ret = (ret + V) / 2.;
				if (Norm(ret - V) < scale * 1E-10)
					return ret;
				p_next = p->getXtuple() + ret; //このp_nextの位置で引っ張られるかチェック
			}
			return ret;
		}
		return {0., 0., 0.};
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
/* ------------------------------------------------------ */
void LaplacianSmoothing(netPp p)
{
	try
	{
		if (!isEdgePoint(p) && isInConvexPolygon(p))
		{ /*端の点はsmoothingしない*/
			auto ps = p->getNeighbors();
			if (ps.size() > 2) // 2点の場合は2点の中点に動いてしまうので，実行しない
				p->setX(Mean(obj3D::extractX(ps)));
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
void LaplacianSmoothing(const V_netPp &ps)
{
	try
	{
		for (const auto &p : ps)
			LaplacianSmoothing(p);
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
///////////////////////////////////////////////////////////
Tddd exact_along_surface(const networkPoint *const p, Tddd VECTOR)
{
	std::vector<std::tuple<Tddd, Tddd>> normals;
	for (const auto &f : p->getFaces())
		for (const auto &F : p->getFaces())
			if (!isFlat(f->getNormalTuple(), F->getNormalTuple(), 1E-2 * M_PI / 180.))
				normals.push_back({f->getNormalTuple(), F->getNormalTuple()});
	Tddd c;
	for (const auto &[N0, N1] : normals)
	{
		c = Normalize(Cross(N0, N1));
		VECTOR = Dot(VECTOR, c) * c;
	}
	return VECTOR;
};
/* ------------------------------------------------------ */
void AreaWeightedSmoothingPreserveShape(netPp p, const double lim_rad = 1E-10)
{
	/*
	AreaWeightedSmoothingVectorの設定に伴って移動が変わる
	*/
	try
	{
		if (!isEdgePoint(p))
		{
			auto ps = p->getNeighbors();
			bool isflat = true;
			auto V = AreaWeightedSmoothingVector(p);
			Tddd N;
			bool allflat = true;
			for (const auto &f : p->getFaces())
				for (const auto &F : p->getFaces())
					if (!isFlat(f->getNormalTuple(), F->getNormalTuple(), 1E-2 * M_PI / 180.))
						allflat = false;
			//もし面が全てフラットなら調べる必要はない
			for (auto kk = 0; kk < 5; ++kk)
			{
				if (ps.size() > 2)
				{
					V = exact_along_surface(p, V);
					for (const auto &f : p->getFaces())
					{
						auto [p0, p1, p2] = f->getPointsTuple(p);
						bool isfiniteangles = isFinite(TriangleAngles(T3Tddd{p0->getXtuple() + V, p1->getXtuple(), p2->getXtuple()}));
						bool isfiniteareas = isFinite(TriangleArea(T3Tddd{p0->getXtuple() + V, p1->getXtuple(), p2->getXtuple()}));
						N = TriangleNormal(T3Tddd{p0->getXtuple() + V, p1->getXtuple(), p2->getXtuple()});
						bool isfinitenormal = isFinite(N);
						if (allflat)
						{
							//全てフラットなので，反転しなければいい．
							if (!isfiniteangles || !isfiniteareas || !isfinitenormal || !isFlat(N, f->getNormalTuple(), M_PI / 180.))
							{
								isflat = false;
								break;
							}
						}
						else
						{
							if (!isfiniteangles || !isfiniteareas || !isfinitenormal || !isFlat(N, f->getNormalTuple(), lim_rad))
							{
								isflat = false;
								break;
							}
						}
					}
					if (isflat)
					{
						p->setX(p->getXtuple() + V);
						return;
					}
				}
				V /= 10.;
			}
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
void AreaWeightedSmoothingPreserveShape(const V_netPp &ps /*copy*/, const double lim_rad = 1E-10)
{
	for (const auto &p : ps)
		AreaWeightedSmoothingPreserveShape(p, lim_rad);
};
void AreaWeightedSmoothingPreserveShape(const std::unordered_set<networkPoint *> &ps /*copy*/, const double lim_rad = 1E-10)
{
	for (const auto &p : ps)
		AreaWeightedSmoothingPreserveShape(p, lim_rad);
};
/* ------------------------------------------------------ */
void LaplacianSmoothingPreserveShape(netPp p, const double lim_rad = 1E-10)
{
	try
	{
		if (!isEdgePoint(p))
		{ /*端の点はsmoothingしない*/
			auto ps = p->getNeighbors();
			bool isflat = true;
			bool allflat = true;
			Tddd N;
			for (const auto &f : p->getFaces())
				for (const auto &F : p->getFaces())
					if (!isFlat(f->getNormalTuple(), F->getNormalTuple(), 1E-2 * M_PI / 180.))
					{
						allflat = false;
						break;
					}
			if (ps.size() > 2) // 2点の場合は2点の中点に動いてしまうので，実行しない
			{
				Tddd V = Mean(extractXtuple(ps)) - p->getXtuple();
				V = exact_along_surface(p, V);

				for (const auto &f : p->getFaces())
				{
					auto [p0, p1, p2] = f->getPointsTuple(p);
					bool isfiniteangles = isFinite(TriangleAngles(T3Tddd{p0->getXtuple() + V, p1->getXtuple(), p2->getXtuple()}));
					bool isfiniteareas = isFinite(TriangleArea(T3Tddd{p0->getXtuple() + V, p1->getXtuple(), p2->getXtuple()}));
					bool isfinitenormal = isFinite(TriangleNormal(T3Tddd{p0->getXtuple() + V, p1->getXtuple(), p2->getXtuple()}));
					// if (!isfiniteangles || !isfiniteareas || !isfinitenormal)
					// {
					// 	isflat = false;
					// 	break;
					// }
					N = TriangleNormal(T3Tddd{p0->getXtuple() + V, p1->getXtuple(), p2->getXtuple()});
					if (allflat)
					{
						//全てフラットなので，反転しなければいい．
						if (!isfiniteangles || !isfiniteareas || !isfinitenormal || !isFlat(N, f->getNormalTuple(), M_PI / 180.))
						{
							isflat = false;
							break;
						}
					}
					else
					{
						if (!isfiniteangles || !isfiniteareas || !isfinitenormal || !isFlat(N, f->getNormalTuple(), lim_rad))
						{
							isflat = false;
							break;
						}
					}
				}
				if (isflat)
					p->setX(p->getXtuple() + V);
			}
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
void LaplacianSmoothingPreserveShape(const V_netPp &ps /*copy*/, const double lim_rad = 1E-10)
{
	for (const auto &p : ps)
		LaplacianSmoothingPreserveShape(p, lim_rad);
};
void LaplacianSmoothingPreserveShape(const std::unordered_set<networkPoint *> &ps /*copy*/, const double lim_rad = 1E-10)
{
	for (const auto &p : ps)
		LaplacianSmoothingPreserveShape(p, lim_rad);
};
/* ------------------------------------------------------ */
void arrangeCORNER(Network &water, const double limit_inner_angle = M_PI / 180 * 5, const int times = 1000)
{
	int count = 0;
	std::cout << "remeshing" << std::endl;
	water.setBounds();
	for (const auto &l : water.getLines())
	{
		auto [p0, p1] = l->getPointsTuple();
		if (!l->CORNER)
		{
			if (!p0->CORNER && !p1->CORNER)
			{
				for (auto &f : l->getFaces())
				{
					auto oppo0 = f->getPointOpposite(l);
					if (oppo0->CORNER)
					{
						if (oppo0->getLinesDirichlet().size() <= 2)
						{
							for (const auto &L : oppo0->getLinesCORNER())
							{
								if ((*L)(oppo0)->getLinesDirichlet().size() == 1)
								{
									if (l->canflip(limit_inner_angle))
									{
										l->flip();
										count++;
										break;
									}
								}
							}
						}
						else if (oppo0->getLinesNeumann().size() <= 2)
						{
							for (const auto &L : oppo0->getLinesCORNER())
							{
								if ((*L)(oppo0)->getLinesNeumann().size() == 1)
								{
									if (l->canflip(limit_inner_angle))
									{
										l->flip();
										count++;
										break;
									}
								}
							}
						}
					}
				}
			}
			else
			{
				for (auto &q : l->getPoints())
				{
					if (q->CORNER && q->getLinesDirichlet().size() >= 3)
					{
						for (auto &f : l->getFaces())
						{
							auto oppo0 = f->getPointOpposite(l);
							if (oppo0->CORNER && oppo0->getLinesDirichlet().size() == 1)
							{
								if (l->canflip(limit_inner_angle))
								{
									l->flip();
									count++;
									break;
								}
							}
						}
					}
					else if (q->CORNER && q->getLinesNeumann().size() >= 3)
					{
						for (auto &f : l->getFaces())
						{
							auto oppo0 = f->getPointOpposite(l);
							if (oppo0->CORNER && oppo0->getLinesNeumann().size() == 1)
							{
								if (l->canflip(limit_inner_angle))
								{
									l->flip();
									count++;
									break;
								}
							}
						}
					}
				}
			}
		}
		/* ------------------------------------------------------ */
		if (count > times)
			break;
	}
};

void flipIf(Network &water, double limit_angle = M_PI / 180., bool force = false, int times = 0)
{
	// 2022/04/13こっちにBEMのmainから持ってきた
	std::cout << "remeshing" << std::endl;
	water.setBounds();
	double mean_length = Mean(extLength(water.getLines()));
	bool isfound = false, ismerged = false;
	int count = 0;
	for (const auto &l : water.getLines())
	{
		auto [p0, p1] = l->getPointsTuple();
		if (!l->CORNER)
		// if (!p0->CORNER && !p1->CORNER)
		{
			if (force && (times == 0 || count < times))
			{
				// p0とp1が角の場合，6という数にこだわる必要がない．
				// つまり6がトポロジカルにベターではない．
				isfound = l->flipIfTopologicalyBetter(limit_angle, M_PI / 180. * 10.);
				if (isfound)
					count++;
			}
			else
			{
				isfound = l->flipIfBetter(limit_angle);
			}
		}
	}
};
void flipIf(Network &water, const Tdd &limit, bool force = false, int times = 0)
{
	/*
	* フリップによって表面の形状が大きく変わってしまうのはよくない．
	* フリップによって三角形の内角がとても小さくなるのもよくない．
	flipIfでは，このようなフリップをどこまで許容するかを与えることができる．
	* limit_angle：フリップを許容する，辺に隣接する面の法線方向の内角
	* limit_inner_angle：フリップを許容する，フリップによってできる三角形の最小内角の最小
	*/
	auto [limit_angle, limit_inner_angle] = limit;
	// 2022/04/13こっちにBEMのmainから持ってきた
	std::cout << "remeshing" << std::endl;
	water.setBounds();
	double mean_length = Mean(extLength(water.getLines()));
	bool isfound = false, ismerged = false;
	int count = 0;
	for (const auto &l : water.getLines())
	{
		auto [p0, p1] = l->getPointsTuple();
		if (!l->CORNER)
		// if (!p0->CORNER && !p1->CORNER)
		{
			if (force && (times == 0 || count < times))
			{
				isfound = l->flipIfTopologicalyBetter(limit_angle, limit_inner_angle);
				if (isfound)
					count++;
			}
			else
			{
				isfound = l->flipIfBetter(limit_angle, limit_inner_angle);
				if (isfound)
					count++;
			}
		}
	}
};
void flipIf(Network &water, const Tdd &limit_Dirichlet, const Tdd &limit_Neumann, bool force = false, int times = 0)
{
	/*
	* フリップによって表面の形状が大きく変わってしまうのはよくない．
	* フリップによって三角形の内角がとても小さくなるのもよくない．
	flipIfでは，このようなフリップをどこまで許容するかを与えることができる．
	* limit_angle：フリップを許容する，辺に隣接する面の法線方向の内角
	* limit_inner_angle：フリップを許容する，フリップによってできる三角形の最小内角の最小
	*/
	auto [limit_angle_D, limit_inner_angle_D] = limit_Dirichlet;
	auto [limit_angle_N, limit_inner_angle_N] = limit_Neumann;
	// 2022/04/13こっちにBEMのmainから持ってきた
	std::cout << "remeshing" << std::endl;
	water.setBounds();
	double mean_length = Mean(extLength(water.getLines()));
	bool isfound = false, ismerged = false;
	int count = 0;
	for (const auto &l : water.getLines())
	{
		auto [p0, p1] = l->getPointsTuple();
		if (!l->CORNER)
		// if (!p0->CORNER && !p1->CORNER)
		{
			if (force && (times == 0 || count < times))
			{
				// p0とp1が角の場合，6という数にこだわる必要がない．
				// つまり6がトポロジカルにベターではない．
				if (l->Dirichlet)
				{
					isfound = l->flipIfTopologicalyBetter(limit_angle_D, limit_inner_angle_D);
					if (isfound)
						count++;
				}
				else
				{
					isfound = l->flipIfTopologicalyBetter(limit_angle_N, limit_inner_angle_N);
					if (isfound)
						count++;
				}
			}
			else
			{
				if (l->Dirichlet)
				{
					isfound = l->flipIfBetter(limit_angle_D, limit_inner_angle_D);
					if (isfound)
						count++;
				}
				else
				{
					isfound = l->flipIfBetter(limit_angle_N, limit_inner_angle_N);
					if (isfound)
						count++;
				}
			}
		}
	}
};

/* ------------------------------------------------------ */
void LaplacianSmoothingIfFlat(netPp p)
{
	try
	{
		if (!isEdgePoint(p) && isFlat(p) && isInConvexPolygon(p))
		{ /*端の点はsmoothingしない*/
			auto ps = p->getNeighbors();
			if (ps.size() > 2) // 2点の場合は2点の中点に動いてしまうので，実行しない
				p->setX(Mean(obj3D::extractX(ps)));
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
//////////////////////////////////////////////////////////
void LaplacianSmoothingIfFlat(V_netPp ps /*copy*/)
{
	for (const auto &p : ps)
		LaplacianSmoothingIfFlat(p);
};
///////////////////////////////////////////////////////////
void LaplacianSmoothingIfFlat_ExceptIntX(V_netPp ps /*copy*/)
{
	try
	{
		std::shuffle(std::begin(ps), std::end(ps), std::default_random_engine());
		for (const auto &p : ps)
			if (!isEdgePoint(p) && isFlat(p) && isInConvexPolygon(p))
			{ //端の点はsmoothingしない
				auto countintx = 0;
				for (const auto &l : p->getLines())
					if (l->isIntxn())
						countintx++;
				if (!(countintx > 1))
				{
					auto tmp = Mean(obj3D::extractX(p->getNeighbors()));
					p->setX(tmp);
				}
			}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
///////////////////////////////////////////////////////////
void LaplacianSmoothingIfOnStraightLine(V_netPp ps)
{
	try
	{
		std::shuffle(std::begin(ps), std::end(ps), std::default_random_engine());
		for (auto p : ps)
			if (!isEdgePoint(p) && !isFlat(p)) //端の点はsmoothingしない
			{
				std::map<netPp, V_d> mapP_vec;
				for (auto q : p->getNeighbors())
					mapP_vec[q] = (q->getX() - p->getX()) / Norm(q->getX() - p->getX());

				netPp p0, p1;
				bool isStraight = false;
				auto Ps = TakeFirst(mapP_vec);
				double minarea = 1E-2;
				for (auto i = 0; i < Ps.size() - 1; i++)
					for (auto j = i + 1; j < Ps.size(); j++)
					{
						// (Ps[i]) <--------- (p) <-------- (Ps[j]) : angle = 0
						double abs_angle = std::abs(MyVectorAngle(Ps[i]->getX() - p->getX(), p->getX() - Ps[j]->getX()));
						if (abs_angle < minarea)
						{
							p0 = Ps[i];
							p1 = Ps[j];
							isStraight = true;
							minarea = abs_angle;
						}
					}
				if (isStraight)
					p->setX((p0->getX() + p1->getX()) / 2.);
			}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
///////////////////////////////////////////////////////////
bool LaplacianSmoothing_IfOnStraight_IntXLine(netPp p)
{
	try
	{
		// V_d A;
		// for(const auto& f:p->getFaces())
		// {
		//   A = f->getAngles(p);//A[0]はこのpの内角
		//   if(A[1]>M_PI*0.8 || A[2]>M_PI*0.8)
		//     return false;
		// }

		// Print("平滑化可能かチェック");
		///////////////////////////////
		auto init_x = p->getX();
		std::map<networkFace *, V_d> Map_F_normal;
		for (const auto &f : p->getFaces())
			Map_F_normal[f] = f->getNormal();
		////////////////////////////////
		V_netPp ps({});
		V_netLp intxlines({});
		netPp q;
		for (const auto &l : p->getLines())
			if (l->isIntxn())
				if ((q = (*l)(p)) && q->isXPoint() && p->isXPoint())
				{
					network::add(ps, q);
					network::add(intxlines, l);
				}

		if (ps.size() == 2)
		{
			// Print("平滑化可能");
			auto v0 = p->getX() - ps[0]->getX();
			auto v1 = ps[1]->getX() - p->getX();
			auto angle = MyVectorAngle(v0, v1);
			if (isFinite(angle) && std::abs(angle) < 1E-2)
			{
				// setできない場合元に戻すようにした
				// setできない理由は，getNormalが計算できない場合である
				auto success = p->setXcarefully((ps[0]->getX() + ps[1]->getX()) / 2.);
				if (!success)
					return false; //問答無用終了
				//ここでは，setできたとしても，変化がある場合元に戻すので，さらにそのチェックを以下で行う．
				for (const auto &f : p->getFaces())
				{
					// Print("修正前後の面の法線方向に大きな違いがある場合，元に戻しfalseを返し終了");
					// Print("finiteの場合のみ比較しよう");
					auto A = Map_F_normal[f];
					auto B = f->getNormal();
					auto AB = Dot(A, B) / (Norm(A) * Norm(B));
					auto dif_deg = acos(AB) * 180. / M_PI;
					if (4. /*[deg]*/ < dif_deg /*[0,180] deg*/)
					{
						// std::cout << Red << "dif_deg = " << dif_deg << reset << std::endl;
						p->setX(init_x);
						// p->setXcarefully(init_x);
						return false;
					}
				}
			}
		}
		return false;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
///////////////////////////////////////////////////////////
std::vector<bool> LaplacianSmoothing_IfOnStraight_IntXLine(const V_netPp &ps)
{
	std::vector<bool> ret;
	for (auto &p : ps)
		ret.emplace_back(LaplacianSmoothing_IfOnStraight_IntXLine(p));
	return ret;
};
///////////////////////////////////////////////////////////
void LaplacianSmoothingWhileKeepingShape(const V_netPp &ps)
{
	LaplacianSmoothingIfFlat_ExceptIntX(ps);
	LaplacianSmoothing_IfOnStraight_IntXLine(ps); // only intx line
};

/////////////////////////////////////////////////////////////////
double calcVolume(const netFp f)
{
	auto p012 = f->getPoints();
	return (f->getArea() * Dot(f->getNormal(), V_d{0., 0., 1.})) * (p012[0]->getX()[2] + p012[1]->getX()[2] + p012[2]->getX()[2]) / 3.;
};
double calcVolume(const V_netFp &fs)
{
	double tmp(0.);
	for (const auto &f : fs)
	{
		double v = calcVolume(f);
		if (std::abs(v) < 1E+30)
			tmp += v;
	}
	return tmp;
};
/////////////////////////////////////////////////////////////////
double Distance(const networkPoint *const p, const Tddd &X)
{
	return Norm(p->getXtuple() - X);
};
double Distance(const networkPoint *const p, const networkPoint *const q)
{
	return Norm(p->getXtuple() - q->getXtuple());
};
double Distance(const networkLine *const l, const networkPoint *const p)
{
	return Norm(p->getXtuple() - l->getXtuple());
};
double Distance(const networkPoint *const p, const networkLine *const l)
{
	return Norm(p->getXtuple() - l->getXtuple());
};
V_d Distance(const networkPoint *const p, const std::vector<networkPoint *> &ps)
{
	V_d ret(ps.size());
	int i = 0;
	for (const auto &q : ps)
		ret[i++] = Distance(q, p);
	return ret;
};
V_d Distance(const std::vector<networkPoint *> &ps, const networkPoint *p)
{
	return Distance(p, ps);
};
void sortByLength(V_netLp &lines)
{
	std::sort(lines.begin(), lines.end(), [](const auto &v, const auto &w)
			  { return v->length() < w->length(); });
};
void sortByDistance(V_netPp &points, const netPp origin)
{
	std::sort(points.begin(), points.end(), [&origin](const netPp &v, const netPp &w)
			  { return Norm(v->getXtuple() - origin->getXtuple()) < Norm(w->getXtuple() - origin->getXtuple()); });
};

////////////////////////////////////////////////////////////////

std::vector<V_netPp> triangulate(const V_netPp &objects, const V_d &normal, const double smallangle = 0.)
{
	// std::cout << __PRETTY_FUNCTION__ << std::endl;
	try
	{
		if (!DuplicateFreeQ(objects))
		{
			std::cout << "objects = " << objects << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Duplication is found !");
		}
		else if (objects.size() < 3)
		{
			Print(objects);
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "objects.size()<3");
		}
		else
		{
			std::vector<V_netPp> ret({});
			geometry::polygon poly(obj3D::extractX(objects));
			for (const auto &ind : poly.triangulate(normal, smallangle))
			{
				V_netPp points = {objects[ind[0]], objects[ind[1]], objects[ind[2]]};
				TriangleNormal(points[0]->getX(), points[1]->getX(), points[2]->getX());
				ret.emplace_back(points);
			}
			return ret;
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
//////////////////////////////////////////////////////////////
//////////////////////////// 外部関数 /////////////////////////
//////////////////////////////////////////////////////////////
std::vector<netPp> getPointsWithoutLines(const std::vector<netPp> &obj)
{
	std::vector<netPp> ret;
	for (const auto &o : obj)
		if (o->getLines().empty())
			ret.emplace_back(o);
	return ret;
};
///////
std::vector<netPp> getPointsWithoutFaces(const std::vector<netPp> &obj)
{
	std::vector<netPp> ret;
	V_netLp ls;
	bool facefound;
	for (const auto &o : obj)
	{
		ls = o->getLines();
		if (ls.empty())
			ret.emplace_back(o);
		else
		{
			facefound = false;
			for (const auto &l : ls)
				if (!l->getFaces().empty())
				{
					facefound = true;
					break;
				}
			if (!facefound)
				ret.emplace_back(o);
		}
	}
	return ret;
};

V_d areas(const V_netFp &fs)
{
	V_d ret(fs.size());
	for (auto i = 0; i < fs.size(); i++)
		ret[i] = fs[i]->getArea();
	return ret;
};
/* ------------------------------------------------------ */
std::vector<std::string> getNames(const V_netFp &fs)
{
	std::vector<std::string> ret;
	for (const auto &f : fs)
		ret.emplace_back(f->getNetwork()->getName());
	return ret;
};

std::vector<std::string> getNames(const V_netPp &ps)
{
	std::vector<std::string> ret;
	for (const auto &p : ps)
		ret.emplace_back(p->getNetwork()->getName());
	return ret;
};

std::vector<std::string> getNames(const V_Netp &nets)
{
	std::vector<std::string> ret;
	for (const auto &n : nets)
		ret.emplace_back(n->getName());
	return ret;
};

V_netLp takeIntxn(const V_netLp &ls)
{
	V_netLp ret({});
	for (const auto &l : ls)
		if (l->isIntxn())
			ret.emplace_back(l);
	return ret;
};

template <typename T>
std::vector<T *> takeNetwork(const std::vector<T *> &ps, const std::vector<Network *> net)
{
	std::vector<T *> ret = {};
	for (const auto &p : ps)
		if (MemberQ(net, p->getNetwork()))
			ret.emplace_back(p);
	return ret;
};

////////////////////////////////////////////////////////////////
////////////////////////// dislay //////////////////////////////
////////////////////////////////////////////////////////////////

template <typename T>
void displayNames(const std::vector<T *> &ps)
{
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	std::map<Network *, int> storages;
	for (const auto &p : ps)
		if (storages.find(p->getNetwork()) != storages.end())
			storages[p->getNetwork()] += 1;
		else
			storages[p->getNetwork()] = 1;
	for (auto [n, i] : storages)
		std::cout << std::setw(8) << n->getName() << "  " << std::setw(8) << i << std::endl;
	if (storages.empty())
		std::cout << Red << "not stored !?" << reset << std::endl;
};

template <typename T>
void displayStorages(const std::vector<T *> &ps)
{
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	std::map<Network *, int> storages;
	for (const auto &p : ps)
		if (storages.find(p->getStorage()) != storages.end())
			storages[p->getStorage()] += 1;
		else
			storages[p->getStorage()] = 1;
	for (auto [n, i] : storages)
		std::cout << std::setw(8) << n->getName() << "  " << std::setw(8) << i << std::endl;
	if (storages.empty())
		std::cout << Red << "not stored !?" << reset << std::endl;
};

void display(const V_netLp &ls)
{
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	std::map<size_t, int> data;
	int intx = 0;
	for (const auto &l : ls)
	{
		auto s = l->getFaces().size();
		if (l->isIntxn())
			intx++;
		if (data.find(s) != data.end())
			data[s] += 1;
		else
			data[s] = 1;
	}
	for (auto [n, i] : data)
		std::cout << std::setw(8) << "Faces size = " << n << "  " << std::setw(8) << i << std::endl;

	std::cout << std::setw(8) << "intx = " << intx << "  " << std::setw(8) << intx << std::endl;
};

void display(Network *net)
{
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	V_i num(15, 0), pointsLines(15, 0), facesPoints(15, 0), facesLines(15, 0), linesFaces(15, 0);
	;
	for (const auto &p : net->getPoints())
		for (auto i = 0; i < 15; i++)
			if (p->getLines().size() == (size_t)i)
				pointsLines[i]++;
	for (const auto &f : net->getFaces())
	{
		auto f_Points = f->getPoints();
		for (auto i = 0; i < 15; i++)
		{
			if (f_Points.size() == (size_t)i)
				facesPoints[i]++;
			if (f->getLines().size() == (size_t)i)
				facesLines[i]++;
		}
	}
	for (const auto &l : net->getLines())
	{
		auto fs = l->getFaces();
		for (auto i = 0; i < 15; i++)
			if (fs.size() == (size_t)i)
				linesFaces[i]++;
	}
	std::cout << "Network name : " << net->getName() << std::endl;
	std::cout << "--------------------Names---------------------------" << std::endl;
	std::cout << Blue << "Points : " << reset;
	// displayNames(net->getPoints());
	std::cout << Magenta << " Faces : " << reset;
	// displayNames(net->getFaces());
	std::cout << "------------------Storages--------------------------" << std::endl;
	std::cout << Blue << "Points : " << reset;
	// displayStorages(net->getPoints());
	std::cout << Magenta << " Faces : " << reset;
	// displayStorages(net->getFaces());
	std::cout << "--------------------Lines--------------------------" << std::endl;
	display(net->getLines());
	std::cout << "-----------------------------------------------" << std::endl;
	std::cout << Magenta << "Points.size() : " << net->getPoints().size() << std::endl;
	std::cout << Magenta << " Faces.size() : " << net->getFaces().size() << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
	int i = 0;
	for (auto &n : num)
	{
		n = i++;
	};
	std::cout << Red << "                  ";
	for (const auto &n : num)
	{
		std::cout << std::setw(6) << n;
	};
	std::cout << reset << std::endl;

	std::cout << magenta << "Lines of points : ";
	for (const auto &n : pointsLines)
	{
		std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
	}
	std::cout << reset << std::endl;

	std::cout << magenta << "Points of faces : ";
	for (const auto &n : facesPoints)
	{
		std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
	};
	std::cout << reset << std::endl;

	std::cout << magenta << " Lines of faces : ";
	for (const auto &n : facesLines)
	{
		std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
	};
	std::cout << reset << std::endl;

	std::cout << magenta << " Faces of Lines : ";
	for (const auto &n : linesFaces)
	{
		std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
	};
	std::cout << reset << std::endl;

	V_i connection(3);
	net->setLinesStatus(true);
	for (const auto &p : net->getPoints())
		for (const auto &line : p->getLines_toggle(true))
		{
			if ((line->getFaces()).size() == 0)
				connection[0]++;
			else if ((line->getFaces()).size() == 1)
				connection[1]++;
			else if ((line->getFaces()).size() == 2)
				connection[2]++;
		}
	std::cout << magenta << " Connection of faces : " << connection;
	if (connection[0] == 0 && connection[1] == 0 && connection[2] != 0)
		std::cout << Blue << " <- Network is a closed surface";
	std::cout << reset << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
};

//* ------------------------------------------------------ */
//*                     粒子法のためのユーティリティ            */
//* ------------------------------------------------------ */
std::vector<Tddd> fullparticlize(const Tddd &X0, const Tddd &X1, const Tddd &X2, const double dx)
{
	/*
	|-o--o--o--o-| n = 4
	x,yはパラメタ
	*/
	interpolationTriangleLinearTuple intp({X0, X1, X2});
	auto y_list = [&intp, &dx](const double x)
	{
		if (1. - x < 1E-10)
			return V_d{0.};
		int n = std::round(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
		if (n > 0)
			return Subdivide(0., 1. - x, n);
		else
			return V_d{};
	};
	/* ------------------------------------------------------ */
	Tddd v = intp(0., 1.) - intp(0., 0.);
	Tddd u = intp(1., 0.) - intp(0., 0.);
	double height = Norm(u - Dot(Normalize(v), u) * Normalize(v));
	int n = std::round(height / dx);
	/* ------------------------------------------------------ */
	std::vector<Tddd> ret = {};
	if (n > 0)
	{
		for (const auto &x : Subdivide(0., 1., n))
			for (const auto &y : y_list(x))
				ret.emplace_back(intp(x, y));
	}
	return ret;
};
std::vector<Tddd> fullparticlize(networkFace const *f, const double dx)
{
	auto ps = f->getPoints();
	return fullparticlize(ps[0]->getXtuple(), ps[1]->getXtuple(), ps[2]->getXtuple(), dx);
};

/*
	% particlizeは座標とパラメタを返す
	% 面をdx間隔で粒子化し，各点の情報をベクトルで返す
*/

// // b% ------------------- タプルからparticlize ------------------ */

// std::vector<std::tuple<Tddd /*実際の座標*/, Tdd /*パラメタt0t1*/>>
// particlize(const T3Tddd &X0X1X2, const double dx)
// {
// 	/*
// 	|-o--o--o--o-| n = 4
// 	x,yはパラメタ
// 	*/
// 	interpolationTriangleLinear3D intp(X0X1X2);
// 	auto y_list = [&intp, &dx](const double x)
// 	{
// 		int n = std::round(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
// 		// int n = std::ceil(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
// 		V_d ret(n);
// 		for (auto i = 0; i < n; ++i)
// 			ret[i] = ((1. - x) / n) * (i + 0.5);
// 		return ret;
// 	};
// 	/* ------------------------------------------------------ */
// 	Tddd v = intp(0., 1.) - intp(0., 0.);
// 	Tddd u = intp(1., 0.) - intp(0., 0.);
// 	double height = Norm(u - Dot(Normalize(v), u) * Normalize(v));
// 	// int n = std::ceil(height / dx);
// 	int n = std::round(height / dx);
// 	// std::cout << Grid({"n", n, "height", height}) << std::endl;
// 	double dH = 1. / n;
// 	V_d x_list(n);
// 	for (auto i = 0; i < n; i++)
// 		x_list[i] = dH * (i + 0.5);
// 	/* ------------------------------------------------------ */
// 	std::vector<std::tuple<Tddd /*実際の座標*/, Tdd /*パラメタt0t1*/>> ret;
// 	ret.reserve(x_list.size() * x_list.size());
// 	for (const auto &x : x_list)
// 		for (const auto &y : y_list(x))
// 			ret.push_back({intp(x, y), Tdd{x, y}});
// 	return ret;
// };

// b% -------------- networkFaceを使ったparticlize ------------- */

std::vector<std::tuple<Tddd, Tdd>> particlize(const networkFace *const f, const double dx) { return particlize(f->getXVertices(), dx); };
std::vector<std::tuple<Tddd, Tdd>> particlize(const V_netFp &faces, const double dx)
{
	std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
	for (const auto &f : faces)
	{
		tmp = particlize(f, dx);
		ret.insert(ret.end(), tmp.begin(), tmp.end());
	}
	return ret;
};
std::vector<std::tuple<Tddd, Tdd>> particlize(const std::unordered_set<networkFace *> &faces, const double dx)
{
	std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
	for (const auto &f : faces)
	{
		tmp = particlize(f, dx);
		ret.insert(ret.end(), tmp.begin(), tmp.end());
	}
	return ret;
};

/* ------------------------ 深さあり ------------------------ */

std::vector<std::tuple<Tddd, Tdd>> particlize(const networkFace *const f, const double dx, const V_d &depth_list /*法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など*/)
{
	std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
	double alpha;
	Tddd X0, X1, X2;
	for (const auto &d : depth_list /*double実際の長さ*/)
	{
		auto [p0, p1, p2] = f->getPointsTuple();
		X0 = p0->getXtuple() + d / Dot(p0->getNormalTuple(), f->getNormalTuple()) * p0->getNormalTuple();
		X1 = p1->getXtuple() + d / Dot(p1->getNormalTuple(), f->getNormalTuple()) * p1->getNormalTuple();
		X2 = p2->getXtuple() + d / Dot(p2->getNormalTuple(), f->getNormalTuple()) * p2->getNormalTuple();
		tmp = particlize({X0, X1, X2}, dx);
		ret.insert(ret.end(), tmp.begin(), tmp.end());
	}
	return ret;
};

Tddd particlize(const networkFace *f, const Tdd &t0t1, const double d)
{
	auto [p0, p1, p2] = f->getPointsTuple();
	return Dot(Tddd{std::get<0>(t0t1), std::get<1>(t0t1), 1. - std::get<0>(t0t1) - std::get<1>(t0t1)},
			   T3Tddd{p0->getXtuple() + d / Dot(p0->getNormalTuple(), f->getNormalTuple()) * p0->getNormalTuple(),
					  p1->getXtuple() + d / Dot(p1->getNormalTuple(), f->getNormalTuple()) * p1->getNormalTuple(),
					  p2->getXtuple() + d / Dot(p2->getNormalTuple(), f->getNormalTuple()) * p2->getNormalTuple()});
};

// b% ------------------------------------------------------ */
using TPPP = std::tuple<networkPoint *, networkPoint *, networkPoint *>;
Tddd oppositeX(const std::tuple<networkFace * /*補間に使った三角形の頂点*/,
								TPPP /*補間に使った三角形の頂点*/,
								Tdd /*パラメタt0,t1*/,
								double /*深さ方向距離*/,
								double /*粒子間隔*/> &particlize_info)
{
	// 粒子化された点の面にたいする反対側の位置を返す．
	// 単純にp+2*Dot(f->p0 - p,n)*nでは，角の面において重なりが生じてしまう．
	// パラメタを使って計算すれば重なりは生じない．
	auto [f, p0p1p2, t0t1, d, dx] = particlize_info;
	if (f)
	{
		auto [p0, p1, p2] = p0p1p2; // f->getPointsTuple();
		Tddd n0 = p0->getNormalTuple();
		Tddd n1 = p1->getNormalTuple();
		Tddd n2 = p2->getNormalTuple();
		return Dot(Tddd{std::get<0>(t0t1), std::get<1>(t0t1), 1. - std::get<0>(t0t1) - std::get<1>(t0t1)},
				   T3Tddd{p0->getXtuple() + d / Dot(p0->getNormalTuple(), f->getNormalTuple()) * p0->getNormalTuple(),
						  p1->getXtuple() + d / Dot(p1->getNormalTuple(), f->getNormalTuple()) * p1->getNormalTuple(),
						  p2->getXtuple() + d / Dot(p2->getNormalTuple(), f->getNormalTuple()) * p2->getNormalTuple()});
	}
	else
	{
		std::stringstream ss;
		// ss << "particlize_info = " << f << ", " << t0 << ", " << t1 << ", " << d << ", " << dx << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
	};
};

Tddd oppositeX(const networkPoint *p)
{
	return oppositeX(p->particlize_info);
};
/* ------------------------------------------------------ */

// std::vector<Tddd> particlize(const V_netFp &faces, const double dx, const V_d &depth_list /*法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など*/)
// {
// 	std::vector<Tddd> ret, tmp;
// 	double alpha;
// 	Tddd X0, X1, X2;
// 	for (const auto &f : faces)
// 	{
// 		tmp = particlize(f, dx, depth_list);
// 		ret.insert(ret.end(), tmp.begin(), tmp.end());
// 	}
// 	return ret;
// };
// std::vector<Tddd> particlize(const std::unordered_set<networkFace *> &faces, const double dx, const V_d &depth_list /*法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など*/)
// {
// 	std::vector<Tddd> ret, tmp;
// 	double alpha;
// 	Tddd X0, X1, X2;
// 	for (const auto &f : faces)
// 	{
// 		tmp = particlize(f, dx, depth_list);
// 		ret.insert(ret.end(), tmp.begin(), tmp.end());
// 	}
// 	return ret;
// };
// std::vector<std::tuple<Tddd, double, double, double /*深さ方向*/, double>> particlizeInfo(const std::unordered_set<networkFace *> &faces, const double dx, const V_d &depth_list /*法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など*/)
// {
// 	std::vector<std::tuple<Tddd, double /*t0*/, double /*t1*/, double /*深さ方向*/, double>> ret, tmp;
// 	double alpha;
// 	Tddd X0, X1, X2;
// 	for (const auto &f : faces)
// 	{
// 		tmp = particlizeInfo(f, dx, depth_list);
// 		ret.insert(ret.end(), tmp.begin(), tmp.end());
// 	}
// 	return ret;
// };
/* ------------------------------------------------------ */
//
VVV_d InterpolateFacesLinear(const netFp f, const double dx, const V_d &depth_list /*法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など*/)
{
	//* およそにdx等分になるように分割する．法線方向にもdxで縮めながらまたは広げてながら，depth_list分だけ線形補間で点を保存していく
	VVV_d particles(depth_list.size());
	try
	{
		int count = 0;
		for (const auto &depth : depth_list)
		{
			V_netPp ps;
			VV_d X(3, V_d(3, 0.));
			V_d c;			   //*center point od a face
			double L0, L1, L2; //*パラメタt0が0->1の実際の長さ
			int m;			   //*範囲が変化するパラメタt1の分割数m
			V_d vec_t1;
			auto &ptcls = particles[count++];
			// for (const auto &f : faces) {
			ps = f->getPoints();
			X = obj3D::extractX(ps);
			for (auto i = 0; i < 3; i++)
				X[i] = X[i] + ps[i]->getNormal() * (depth) / Dot(ps[i]->getNormal(), f->getNormal());
			// X[i] = X[i] + ps[i]->getNormal() * dx * (double)l;  //深さ方向への移動
			//! edge points at this depth
			for (auto &x : X)
				ptcls.emplace_back(x);
			//! move edge points inwards

			interpolationTriangleLinear intp_to_shift(X);
			auto X01 = X[1] - X[0];
			auto X12 = X[2] - X[1];
			auto X20 = X[0] - X[2];
			L0 = std::sqrt(std::pow(Norm(X01) * Norm(X12), 2.) - std::pow(Dot(X01, X12), 2)) / Norm(X12); // X0から対角線までの距離
			auto dt0 = dx / L0 / 2.;
			L1 = std::sqrt(std::pow(Norm(X12) * Norm(X20), 2.) - std::pow(Dot(X12, X20), 2)) / Norm(X20); // X1から対角線までの距離
			auto dt1 = dx / L1 / 2.;
			L2 = std::sqrt(std::pow(Norm(X20) * Norm(X01), 2.) - std::pow(Dot(X20, X01), 2)) / Norm(X01); // X2から対角線までの距離
			auto dt2 = dx / L2 / 2.;
			X[0] = intp_to_shift(dt0, dt1);
			X[1] = intp_to_shift(1. - dt2 - dt1, dt1);
			X[2] = intp_to_shift(dt0, 1. - dt0 - dt2);
			// c = Mean(X);
			// for (auto &x : X)
			// 	x += (c - x) / Norm((c - x)) * dx; //内部に移動させる
			X01 = X[1] - X[0];
			X12 = X[2] - X[1];
			X20 = X[0] - X[2];
			L0 = std::sqrt(std::pow(Norm(X01) * Norm(X12), 2.) - std::pow(Dot(X01, X12), 2)) / Norm(X12); // X0から対角線までの距離
			interpolationTriangleLinear intp(X);
			for (const auto &t0 : Subdivide(0., 1., std::round(L0 / dx)))
			{
				m = std::round(Norm(intp(t0, 0.) - intp(t0, 1. - t0)) / dx);
				vec_t1 = Subdivide(0., 1. - t0, m > 0 ? m : 1);
				for (const auto &t1 : vec_t1)
					ptcls.emplace_back(intp(t0, t1));
			}
			// }
		}
		return particles;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};

VVV_d InterpolateFacesLinear(const V_netFp &faces, const double dx, const V_d &depth_list /*法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など*/)
{
	//* およそにdx等分になるように分割する．法線方向にもdxで縮めながらまたは広げてながら，depth_list分だけ線形補間で点を保存していく
	VVV_d particles(depth_list.size());
	try
	{
		int count = 0;
		for (const auto &depth : depth_list)
		{
			V_netPp ps;
			VV_d X(3, V_d(3, 0.));
			V_d c;			   //*center point od a face
			double L0, L1, L2; //*パラメタt0が0->1の実際の長さ
			int m;			   //*範囲が変化するパラメタt1の分割数m
			V_d vec_t1;
			auto &ptcls = particles[count++];
			for (const auto &f : faces)
			{
				ps = f->getPoints();
				X = obj3D::extractX(ps);
				for (auto i = 0; i < 3; i++)
					X[i] = X[i] + ps[i]->getNormal() * (depth) / Dot(ps[i]->getNormal(), f->getNormal());
				//深さ方向への移動．ただし，節点の法線方向なので，必ずしも，面の法線方向と一致しない．
				//! edge points at this depth
				for (auto &x : X)
					ptcls.emplace_back(x);
				//! move edge points inwards

				interpolationTriangleLinear intp_to_shift(X);
				auto X01 = X[1] - X[0];
				auto X12 = X[2] - X[1];
				auto X20 = X[0] - X[2];
				L0 = std::sqrt(std::pow(Norm(X01) * Norm(X12), 2.) - std::pow(Dot(X01, X12), 2)) / Norm(X12); // X0から対角線までの距離
				auto dt0 = dx / L0 / 2.;
				L1 = std::sqrt(std::pow(Norm(X12) * Norm(X20), 2.) - std::pow(Dot(X12, X20), 2)) / Norm(X20); // X1から対角線までの距離
				auto dt1 = dx / L1 / 2.;
				L2 = std::sqrt(std::pow(Norm(X20) * Norm(X01), 2.) - std::pow(Dot(X20, X01), 2)) / Norm(X01); // X2から対角線までの距離
				auto dt2 = dx / L2 / 2.;
				X[0] = intp_to_shift(dt0, dt1);
				X[1] = intp_to_shift(1. - dt2 - dt1, dt1);
				X[2] = intp_to_shift(dt0, 1. - dt0 - dt2);
				// c = Mean(X);
				// for (auto &x : X)
				// 	x += (c - x) / Norm((c - x)) * dx; //内部に移動させる
				X01 = X[1] - X[0];
				X12 = X[2] - X[1];
				X20 = X[0] - X[2];
				L0 = std::sqrt(std::pow(Norm(X01) * Norm(X12), 2.) - std::pow(Dot(X01, X12), 2)) / Norm(X12); // X0から対角線までの距離
				interpolationTriangleLinear intp(X);
				for (const auto &t0 : Subdivide(0., 1., std::round(L0 / dx)))
				{
					m = std::round(Norm(intp(t0, 0.) - intp(t0, 1. - t0)) / dx);
					vec_t1 = Subdivide(0., 1. - t0, m > 0 ? m : 1);
					for (const auto &t1 : vec_t1)
						ptcls.emplace_back(intp(t0, t1));
				}
			}
		}
		return particles;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};

// /* ------------------------ 近傍点探索 ----------------------- */
// //テンプレートどんなオブジェクトでもできるはず
// template <typename T>
// class Buckets
// {
// public:
// 	//!getX()でxyz座標を取得できるオブジェクトTのための，バケツ
// 	Tdd xbounds;
// 	Tdd ybounds;
// 	Tdd zbounds;
// 	int xsize, ysize, zsize;
// 	T3Tdd bounds;
// 	Tiii dn;
// 	std::vector</*x*/ std::vector</*y*/ std::vector</*z*/ std::vector<T *>>>> buckets;
// 	double dL;

// 	Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : bounds(boundingboxIN), dL(dL_IN), dn({0, 0, 0})
// 	{
// 		set(boundingboxIN, dL_IN);
// 	};
// 	void set(const T3Tdd &boundingboxIN, const double dL_IN)
// 	{
// 		this->bounds = boundingboxIN;
// 		this->dL = dL_IN;
// 		this->xbounds = std::get<0>(this->bounds);
// 		this->ybounds = std::get<1>(this->bounds);
// 		this->zbounds = std::get<2>(this->bounds);
// 		this->xsize = (int)((std::get<1>(this->xbounds) - std::get<0>(this->xbounds)) / this->dL);
// 		this->ysize = (int)((std::get<1>(this->ybounds) - std::get<0>(this->ybounds)) / this->dL);
// 		this->zsize = (int)((std::get<1>(this->zbounds) - std::get<0>(this->zbounds)) / this->dL);
// 		this->dn = {xsize, ysize, zsize};
// 		this->buckets.resize(xsize,
// 							 std::vector</*y*/ std::vector</*z*/ std::vector<T *>>>(ysize,
// 																					std::vector</*z*/ std::vector<T *>>(zsize, std::vector<T *>(0))));
// 	};
// 	//x座標を内包するバケツのインデックスを返す
// 	void indices(const Tddd &x, int &i, int &j, int &k) const
// 	{
// 		i = (int)((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL);
// 		j = (int)((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL);
// 		k = (int)((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL);
// 	};
// 	Tiii indices(const Tddd &x) const
// 	{
// 		return {(int)((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL),
// 				(int)((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL),
// 				(int)((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL)};
// 	};
// 	Tddd indices2X(const int i, const int j, const int k) const
// 	{
// 		return {i * this->dL + std::get<0>(this->xbounds) /*min*/,
// 				j * this->dL + std::get<0>(this->ybounds) /*min*/,
// 				k * this->dL + std::get<0>(this->zbounds) /*min*/};
// 	};
// 	//インデックスがboundsに収まっているかどうか
// 	bool isInside(const int i, const int j, const int k) const
// 	{
// 		return (!(i < 0 || j < 0 || k < 0 || i >= this->xsize || j >= this->ysize || k >= this->zsize));
// 	};
// 	bool isInside(const Tiii &i) const
// 	{
// 		return isInside(std::get<0>(i), std::get<1>(i), std::get<2>(i));
// 	};
// 	//
// 	void add(const Tddd &x, T *p)
// 	{
// 		try
// 		{
// 			int i, j, k;
// 			this->indices(x, i, j, k);
// 			if (isInside(i, j, k))
// 				this->buckets[i][j][k].emplace_back(p);
// 		}
// 		catch (error_message &e)
// 		{
// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 		}
// 	};
// 	//座標を入力し，バケツを指定する．
// 	std::vector<T *> getAllBuckets() const { return this->buckets; };
// 	std::vector<T *> getBucket(const Tddd &x) const
// 	{
// 		Tiii i = this->indices(x);
// 		return this->buckets[std::get<0>(i)][std::get<1>(i)][std::get<2>(i)];
// 	};
// 	/* -------------- デフォルトのバケツは，深さ毎に粒子を保存していく -------------- */
// 	std::vector<std::vector<T *>> getObjects(const Tddd &x, const int depth /*limit depth*/, const int limit_number = 100000) const
// 	{
// 		std::vector<std::vector<T *>> ret(depth);
// 		// std::vector<T *> ret_last(0);
// 		int i0, j0, k0;
// 		this->indices(x, i0, j0, k0);
// 		//* ------------------------ depth=0 ------------------------ */
// 		auto &r = ret[0];
// 		if (isInside(i0, j0, k0))
// 			r.insert(r.end(), this->buckets[i0][j0][k0].begin(), this->buckets[i0][j0][k0].end());
// 		int tot = r.size();
// 		int i, j, k;
// 		for (auto d = 1; d < depth; d++)
// 		{
// 			// ret_last.clear();
// 			// if (d == 0)
// 			// {
// 			// 	//* ------------------------ depth=0 ------------------------ */
// 			// 	if (isInside(i0, j0, k0))
// 			// 		ret_last.insert(ret_last.end(), this->buckets[i0][j0][k0].begin(), this->buckets[i0][j0][k0].end());
// 			// }
// 			// else
// 			// {
// 			//* ------------------------ depth>0の周囲 ------------------------ */
// 			// int i, j, k;
// 			// for (auto j = j0 - d; j <= j0 + d; j++)
// 			// 	for (auto k = k0 - d; k <= k0 + d; k++)
// 			// 	{
// 			// 		i = i0 + d;
// 			// 		if (isInside(i, j, k))
// 			// 			ret_last.insert(ret_last.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 			// 		i = i0 - d;
// 			// 		if (isInside(i, j, k))
// 			// 			ret_last.insert(ret_last.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 			// 	}
// 			// for (auto i = i0 - d + 1; i < i0 + d; i++)
// 			// 	for (auto k = k0 - d; k <= k0 + d; k++)
// 			// 	{
// 			// 		j = j0 + d;
// 			// 		if (isInside(i, j, k))
// 			// 			ret_last.insert(ret_last.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 			// 		j = j0 - d;
// 			// 		if (isInside(i, j, k))
// 			// 			ret_last.insert(ret_last.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 			// 	}
// 			// for (auto i = i0 - d + 1; i < i0 + d; i++)
// 			// 	for (auto j = j0 - d + 1; j < j0 + d; j++)
// 			// 	{
// 			// 		k = k0 + d;
// 			// 		if (isInside(i, j, k))
// 			// 			ret_last.insert(ret_last.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 			// 		k = k0 - d;
// 			// 		if (isInside(i, j, k))
// 			// 			ret_last.insert(ret_last.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 			// 	}
// 			//* ------------------------------------------------------ */
// 			r = ret[d];
// 			for (auto j = (j0 - d < 0 ? 0 : j0 - d); j <= j0 + d && j < this->ysize; j++)
// 				for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; k++)
// 				{
// 					i = i0 + d;
// 					if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].cbegin(), this->buckets[i][j][k].cend());
// 					i = i0 - d;
// 					if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].cbegin(), this->buckets[i][j][k].cend());
// 				}
// 			for (auto i = (i0 - d + 1 < 0 ? 0 : i0 - d + 1); i < i0 + d && i < this->xsize; i++)
// 			{
// 				for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; k++)
// 				{
// 					j = j0 + d;
// 					if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].cbegin(), this->buckets[i][j][k].cend());
// 					j = j0 - d;
// 					if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].cbegin(), this->buckets[i][j][k].cend());
// 				}
// 				for (auto j = (j0 - d + 1 < 0 ? 0 : j0 - d + 1); j < j0 + d && j < this->ysize; j++)
// 				{
// 					k = k0 + d;
// 					if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].cbegin(), this->buckets[i][j][k].cend());
// 					k = k0 - d;
// 					if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].cbegin(), this->buckets[i][j][k].cend());
// 				}
// 			}
// 			tot += r.size();
// 			if (tot >= limit_number)
// 				return ret;
// 		}
// 		return ret;
// 	};
// 	std::vector<std::vector<T *>> getObjects(const Tddd &x, const double smoothing_length, const int limit_number = 100000) const
// 	{
// 		int depth = std::ceil(smoothing_length / this->dL);
// 		return getObjects(x, depth, limit_number);
// 	};
// 	/* ------------------------------------------------------ */
// 	std::vector<T *> getObjectsFlattened(const Tddd &x, const int d /*limit depth*/, const int limit_number = 100000) const
// 	{
// 		std::vector<T *> ret;
// 		int i0, j0, k0;
// 		this->indices(x, i0, j0, k0);
// 		int i_beg = ((i0 - d) >= 0 ? (i0 - d) : 0);
// 		int j_beg = ((j0 - d) >= 0 ? (j0 - d) : 0);
// 		int k_beg = ((k0 - d) >= 0 ? (k0 - d) : 0);
// 		int i_end = ((i0 + d) <= this->xsize ? (i0 + d) : this->xsize);
// 		int j_end = ((j0 + d) <= this->ysize ? (j0 + d) : this->ysize);
// 		int k_end = ((k0 + d) <= this->zsize ? (k0 + d) : this->zsize);
// 		for (auto it = this->buckets.begin() + i_beg; it != this->buckets.begin() + i_end; ++it)
// 			for (auto jt = it->begin() + j_beg; jt != it->begin() + j_end; ++jt)
// 				for (auto kt = jt->begin() + k_beg; kt != jt->begin() + k_end; ++kt)
// 					ret.insert(ret.end(), kt->begin(), kt->end());
// 		return ret;
// 	};
// 	std::vector<T *> getObjectsFlattened(const Tddd &x, const double smoothing_length, const int limit_number = 100000) const
// 	{
// 		int depth = std::ceil(smoothing_length / this->dL);
// 		return getObjectsFlattened(x, depth, limit_number);
// 	};
// 	/* ------------------------------------------------------ */
// 	void apply(const std::function<void(T *)> &fun, const Tddd &x, const int d /*limit depth*/, const int limit_number = 100000)
// 	{
// 		int i0, j0, k0;
// 		this->indices(x, i0, j0, k0);
// 		int i_beg = ((i0 - d) >= 0 ? (i0 - d) : 0);
// 		int j_beg = ((j0 - d) >= 0 ? (j0 - d) : 0);
// 		int k_beg = ((k0 - d) >= 0 ? (k0 - d) : 0);
// 		int i_end = ((i0 + d) <= this->xsize ? (i0 + d) : this->xsize);
// 		int j_end = ((j0 + d) <= this->ysize ? (j0 + d) : this->ysize);
// 		int k_end = ((k0 + d) <= this->zsize ? (k0 + d) : this->zsize);
// 		for (auto it = this->buckets.begin() + i_beg; it != this->buckets.begin() + i_end; ++it)
// 			for (auto jt = it->begin() + j_beg; jt != it->begin() + j_end; ++jt)
// 				for (auto kt = jt->begin() + k_beg; kt != jt->begin() + k_end; ++kt)
// 					for (auto objT = kt->begin(); objT != kt->end(); ++objT)
// 						fun(*objT);
// 	};
// 	/* -------------- デフォルトのバケツは，深さ毎に粒子を保存していく -------------- */
// 	// std::map<V_d, std::vector<T *>> getMapBuckets(const V_d &x, const int depth /*limit depth*/, const int limit_number = 100000) const
// 	// {
// 	// 	V_i ind = this->indices(x);
// 	// 	std::map<V_d, std::vector<T *>> ret;
// 	// 	int tot = 0;
// 	// 	auto i0 = ind[0], j0 = ind[1], k0 = ind[2];
// 	// 	for (auto d = 0; d < depth; d++)
// 	// 	{
// 	// 		if (d == 0)
// 	// 		{
// 	// 			//* ------------------------ depth=0 ------------------------ */
// 	// 			auto &tmp = ret[indices2X(i0, j0, k0)];
// 	// 			tmp = {};
// 	// 			if (isInside(i0, j0, k0))
// 	// 				tmp.insert(tmp.end(), this->buckets[i0][j0][k0].begin(), this->buckets[i0][j0][k0].end());
// 	// 		}
// 	// 		else
// 	// 		{
// 	// 			//* ------------------------ depth>0の周囲 ------------------------ */
// 	// 			int i, j, k;
// 	// 			for (auto j = j0 - d; j <= j0 + d; j++)
// 	// 				for (auto k = k0 - d; k <= k0 + d; k++)
// 	// 				{
// 	// 					i = i0 + d;
// 	// 					auto &tmp = ret[indices2X(i, j, k)];
// 	// 					tmp = {};
// 	// 					if (isInside(i, j, k))
// 	// 					{
// 	// 						tmp.insert(tmp.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 	// 					}
// 	// 					i = i0 - d;
// 	// 					tmp = ret[indices2X(i, j, k)];
// 	// 					tmp = {};
// 	// 					if (isInside(i, j, k))
// 	// 					{
// 	// 						tmp.insert(tmp.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 	// 					}
// 	// 				}
// 	// 			for (auto i = i0 - d + 1; i < i0 + d; i++)
// 	// 				for (auto k = k0 - d; k <= k0 + d; k++)
// 	// 				{
// 	// 					j = j0 + d;
// 	// 					auto &tmp = ret[indices2X(i, j, k)];
// 	// 					tmp = {};
// 	// 					if (isInside(i, j, k))
// 	// 					{
// 	// 						tmp.insert(tmp.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 	// 					}
// 	// 					j = j0 - d;
// 	// 					tmp = ret[indices2X(i, j, k)];
// 	// 					tmp = {};
// 	// 					if (isInside(i, j, k))
// 	// 					{
// 	// 						tmp.insert(tmp.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 	// 					}
// 	// 				}
// 	// 			for (auto i = i0 - d + 1; i < i0 + d; i++)
// 	// 				for (auto j = j0 - d + 1; j < j0 + d; j++)
// 	// 				{
// 	// 					k = k0 + d;
// 	// 					auto &tmp = ret[indices2X(i, j, k)];
// 	// 					tmp = {};
// 	// 					if (isInside(i, j, k))
// 	// 					{
// 	// 						tmp.insert(tmp.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 	// 					}
// 	// 					k = k0 - d;
// 	// 					tmp = ret[indices2X(i, j, k)];
// 	// 					tmp = {};
// 	// 					if (isInside(i, j, k))
// 	// 					{
// 	// 						tmp.insert(tmp.end(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 	// 					}
// 	// 				}
// 	// 			//* ------------------------------------------------------ */
// 	// 		}
// 	// 		tot = 0;
// 	// 		for (const auto &[x, V] : ret)
// 	// 			tot += V.size();
// 	// 		if (tot >= limit_number)
// 	// 			return ret;
// 	// 	}
// 	// 	return ret;
// 	// };
// 	// std::map<V_d, std::vector<T *>> getMapBuckets(const V_d &x, const double smoothing_length, const int limit_number = 100000) const
// 	// {
// 	// 	int depth = std::ceil(smoothing_length / this->dL);
// 	// 	return getMapBuckets(x, depth, limit_number);
// 	// };
// 	// std::vector<std::vector<T *>> get125Buckets(const V_d &x) const
// 	// {
// 	// 	V_i ind = this->indices(x);
// 	// 	std::vector<std::vector<T *>> ret(125, std::vector<T *>(0));
// 	// 	int n = 0;
// 	// 	for (auto i = ind[0] - 2; i <= ind[0] + 2; i++)
// 	// 		for (auto j = ind[1] - 2; j <= ind[1] + 2; j++)
// 	// 			for (auto k = ind[2] - 2; k <= ind[2] + 2; k++)
// 	// 			{
// 	// 				if (isInside(i, j, k))
// 	// 					ret[n] = this->buckets[i][j][k];
// 	// 				n++;
// 	// 			}
// 	// 	return ret;
// 	// };
// 	// std::vector<std::vector<T *>> get27Buckets(const V_d &x) const
// 	// {
// 	// 	V_i ind = this->indices(x);
// 	// 	std::vector<std::vector<T *>> ret(27, std::vector<T *>(0));
// 	// 	int n = 0;
// 	// 	for (auto i = ind[0] - 1; i <= ind[0] + 1; i++)
// 	// 		for (auto j = ind[1] - 1; j <= ind[1] + 1; j++)
// 	// 			for (auto k = ind[2] - 1; k <= ind[2] + 1; k++)
// 	// 			{
// 	// 				if (isInside(i, j, k))
// 	// 					ret[n] = this->buckets[i][j][k];
// 	// 				n++;
// 	// 			}
// 	// 	return ret;
// 	// };
// };

// networkPoint限定のもの
//  template <>
//  class Buckets<networkPoint> : public BaseBuckets<networkPoint>
//  {
//  public:
//  	// Buckets<networkPoint>(const geometry::CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<networkPoint>(c_bounds, dL_IN){};
//  	// Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<networkPoint>(boundingboxIN, dL_IN){};
//  };
/* ------------------------------------------------------ */
//この機能を入れ込む
// class PointsBuckets : public Buckets<networkPoint>
// {
// public:
// 	PointsBuckets(const T3Tdd &boundingboxIN, const double dL_IN)
// 		: Buckets<networkPoint>(boundingboxIN, dL_IN){};
// 	//
// 	void add(const V_netPp &ps)
// 	{
// 		try
// 		{
// 			if (ps.empty())
// 				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "points is empty vector");
// 			for (const auto &p : ps)
// 				Buckets::add(p->getXtuple(), p);
// 		}
// 		catch (std::exception &e)
// 		{
// 			std::cerr << e.what() << reset << std::endl;
// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 		};
// 	};
// 	void clear()
// 	{
// 		this->buckets.clear();
// 	};
// };

std::vector<networkPoint *> operator+(std::vector<networkPoint *> A /*copy and return*/, const std::vector<networkPoint *> &B)
{
	A.insert(A.end(), B.begin(), B.end());
	return A;
};
std::vector<networkPoint *> &operator+=(std::vector<networkPoint *> &v, const std::vector<networkPoint *> &w)
{
	return (v = v + w);
};
#endif
