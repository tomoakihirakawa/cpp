#ifndef networkPoint_H
#define networkPoint_H
#pragma once

#include "Network.hpp"
#include "InterpolationRBF.hpp"

inline Tddd networkPoint::normalDirichletFace() const
{
	//接触面としての法線方向を計算する良い方法が必要．
	Tddd normal = {0., 0., 0.};
	int count = 0;
	for (auto &f : this->getFaces())
	{
		auto [p0, p1, p2] = f->getPointsTuple();
		if ((p0->Dirichlet || p0->CORNER) && (p1->Dirichlet || p1->CORNER) && (p2->Dirichlet || p2->CORNER))
			normal += f->getAngle(this) * f->getNormalTuple();
	}
	if (count == 0)
		error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Dirichlet surface not found");
	return Normalize(normal);
};

inline Tddd networkPoint::normalNeumannFace() const
{
	Tddd normal = {0., 0., 0.};
	int count = 0;
	for (auto &f : this->getFaces())
	{
		auto [p0, p1, p2] = f->getPointsTuple();
		if ((p0->Neumann || p0->CORNER) && (p1->Neumann || p1->CORNER) && (p2->Neumann || p2->CORNER))
		{
			normal += f->getAngle(this) * f->getNormalTuple();
		}
	}
	if (count == 0)
		error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Dirichlet surface not found");
	return Normalize(normal);
};

inline Tddd networkPoint::normalContanctSurface(const double pw0 = 1., const double pw1 = 1.) const
{

	//接触面としての法線方向を計算する良い方法が必要．
	double len = this->radius;
	Tddd n = {0, 0, 0};
	double totlen = 0;
	n = (1. / pow(len, pw0)) * this->getNormalTuple();
	totlen = 1. / pow(len, pw0);
	if (!this->getContactFaces().empty())
	{
		std::unordered_set<networkFace *> contactfaces;
		for (auto &f : this->getContactFaces())
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
		//
		for (const auto &f : contactfaces)
		{
			len = Norm(f->mirrorPosition(this) - this->getXtuple());
			n += (1. / pow(len, pw1)) * sgn(Dot(f->getNormalTuple(), this->getNormalTuple())) * f->getNormalTuple();
			totlen += 1. / pow(len, pw1);
		}
	}
	return n / totlen;
};

inline Tddd networkPoint::reflect(const Tddd &v) const
{
	// particlizeした点だけが使える．
	auto f = std::get<0>(particlize_info);
	if (f)
	{
		auto n = f->getNormalTuple();
		return v - 2. * Dot(v, n) * n;
	}
	else
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
};
/*
@ addContactPointsは，引数として半径が与えられない場合，各点のradiusが成す球として，互いの重なりを調べ，重なった点はContactPointsに保存する．
*/
inline void networkPoint::addContactFaces(const Buckets<networkFace> &B, bool include_self_network = true)
{
	//!まずは，衝突があり得そうな面を多めに保存する．
	std::unordered_set<networkFace *> tmpContactFaces;
	for (const auto &f : DeleteDuplicates(Flatten(B.getObjects(this->getXtuple(), 2. * this->radius))))
	{
		if (!(!include_self_network && (f->getNetwork() == this->getNetwork())))
		{
			auto [p0, p1, p2] = f->getPointsTuple();
			if (intersection(geometry::Sphere(this->getXtuple(), this->radius), geometry::Triangle(p0->getXtuple(), p1->getXtuple(), p2->getXtuple())).isIntersecting)
				tmpContactFaces.emplace(f);
		}
	}
	//!各面について衝突があり得るか詳しく調べて判断する．
	for (const auto &f : tmpContactFaces)
	{
		auto shadowed = false;
		//@ 周りの点に隠れて面が見えない状況の場合，その面とは接し得ないので，考慮しない．
		auto [p0, p1, p2] = f->getPointsTuple();
		// BEMのために導入したもの．周辺の点の方向にある面には衝突し得ないため，無視する．
		// ただし修正面の方向をむくせんを抜き出す必要がある．
		auto n = f->getNormalTuple();
		Tddd thisPointToFace = vectorToTriangle(f, this->getXtuple());
		/*
			 |
		Face |<---thisPointToFace---- thisPoint
			 |
		*/
		/*
		_____________________
		|       \  A
		| shadow \ |
		|      60 \|
		|  q<------*--------q
		|  |
		|  |
		|  *
		*/
		for (const auto &q : this->getNeighbors())
		{
			auto thisPointToNeighbor = q->getXtuple() - this->getXtuple();
			auto angle = MyVectorAngle(thisPointToNeighbor, thisPointToFace);
			if (180. / M_PI * angle < 60. || Norm(thisPointToFace) > this->radius)
			{
				shadowed = true;
				break;
			}
		}
		if (Dot(f->getNormalTuple(), this->getNormalTuple()) < 0. /*少なくとも向き合っていなければいけない*/)
			if (!shadowed)
				this->ContactFaces.emplace(f);
	}
};
//点なのに，sphereでインターセクションをチェックするのは効率的ではない．
inline void networkPoint::addContactPoints(const Buckets<networkPoint> &B,
										   const bool include_self_network = true)
{
	for (const auto &q : B.getObjects_unorderedset(this->getXtuple(), 2. * this->radius /*depth*/))
		if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
			if (this != q)
				if (intersection(geometry::Sphere(this->getXtuple(), this->radius), geometry::Sphere(this->getXtuple(), q->radius)).isIntersecting)
					this->ContactPoints.emplace(q);
};
// radiusのしていがある場合は，その半径内の点を取得する．
inline void networkPoint::addContactPoints(const Buckets<networkPoint> &B,
										   const double radius,
										   const bool include_self_network = true)
{
	for (const auto &q : B.getObjects_unorderedset(this->getXtuple(), 2. * this->radius /*depth*/))
		if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
			if (this != q)
				if (Norm(this->getXtuple() - q->getXtuple()) <= radius)
					this->ContactPoints.emplace(q);
};
inline void networkPoint::addContactPoints(const Buckets<networkPoint> &B,
										   const int limit_depth,
										   const int limit_num,
										   const bool include_self_network = true)
{
	for (const auto &q : B.getObjects_unorderedset(this->getXtuple(), limit_depth, limit_num))
		if (!(!include_self_network && (q->getNetwork() == this->getNetwork())))
			if (this != q)
				if (Norm(this->getXtuple() - q->getXtuple()) <= radius)
					this->ContactPoints.emplace(q);
};
// radiusのしていがある場合は，その半径内の点を取得する．
inline void networkPoint::addContactPoints(const std::vector<Buckets<networkPoint>> &Bs,
										   const double radius,
										   const bool include_self_network = true)
{
	for (const auto &B : Bs)
		addContactPoints(B, radius, include_self_network);
};
inline void networkPoint::addContactPoints(const std::vector<Buckets<networkPoint>> &Bs,
										   const int limit_depth,
										   const int limit_num,
										   const bool include_self_network = true)
{
	for (const auto &B : Bs)
		addContactPoints(B, limit_depth, limit_num, include_self_network);
};
// inline void networkPoint::saveContactPoints(const Buckets<networkPoint> &B, double radius = 1E-40, bool exclude_self = false)
// {
// 	if (radius <= 1E-30)
// 		radius = this->radius;

// 	for (const auto &q : DeleteDuplicates(Flatten(B.getObjects(this->getXtuple(), 2. * radius /*depth*/))))
// 		if (!(exclude_self && (q->getNetwork() == this->getNetwork())))
// 		{
// 			if (this != q)
// 			{
// 				auto ixn = intersection(geometry::Sphere(this->getXtuple(), this->radius),
// 										geometry::Sphere(this->getXtuple(), q->radius));
// 				if (ixn.isIntersecting)
// 					this->ContactPoints.emplace_back(q);
// 			}
// 		}
// };

// inline V_netFp networkPoint::getFacesSort() const
// {
//     V_netFp ret;
//     auto this_ls = this->getLines();
//     auto first_l = this_ls[0];
//     auto next_l = first_l;
//     ///////
//     netFp next_f = nullptr;
//     netLp next_l;
//     // first try
//     for (const auto &current_f : next_l->getFaces())
//     {
//         auto tmp_l = current_f->getLine(next_l, -1);
//         if (MemberQ(this_ls, tmp_l))
//         {
//             ret.emplace_back(current_f);
//             next_l = tmp_l;
//             next_f = current_f;
//             break;
//         }
//     }
//     //////
//     if (next_f)
//     {
//         do
//         {
//             for (const auto &current_f : next_l->getFacesExcept(next_f))
//             {
//                 auto tmp_l = current_f->getLine(next_l, -1);
//                 if (network::erase(this_ls, tmp_l))
//                 {
//                     ret.emplace_back(current_f);
//                     next_l = tmp_l;
//                     next_f = current_f;
//                     break;
//                 }
//             }
//         } while ()
//     }
//     else
//         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "cannot find line at first try");
// };

//追加2021/06/18
V_netFp networkPoint::getFacesSort() const
{
	V_netFp ret = {(this->getFaces())[0]};
	netFp f;
	netLp l;
	int count = 0;
	/*
	 *        *-------*
	 *       / \    /  \
	 *      / f \  /    \
	 *     *-- l -0------*
	 *      \    / \    /
	 *       \  /   \  /
	 *        *-------*
	 */
	do
	{
		f = *ret.rbegin();
		l = f->getLinesFrom(this)[2];
		if ((*l)(f) != ret[0])
			ret.emplace_back((*l)(f));
		else
			break;
	} while (count++ < 1000);
	return ret;
};
/* ------------------------------------------------------ */
//追加2021/09/02
inline V_netFp networkPoint::getFacesSort(networkLine *const line) const
{
	try
	{
		auto p1 = (*line)(this);
		if (!p1)
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		networkFace *begin_face;
		//左回りの面を選ぶ．その面からスタート
		for (const auto &f : line->getFaces())
			if (f->getPointFront(line) == p1)
			{
				begin_face = f;
				break;
			}

		V_netFp ret = {begin_face};
		ret.reserve(10);
		int c = 0;
		networkLine *next_line = line;
		do
		{
			next_line = begin_face->getLineBack(next_line);
			if (next_line != line)
				ret.emplace_back(begin_face = (*next_line)(begin_face));
			else
				break;
		} while (c++ < 1000);
		if (c >= 1000)
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

		return ret;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
/* ------------------------------------------------------ */
inline V_netFp networkPoint::getFaces(networkLine *line) const
{
	return getFacesSort(line);
};
//追加2021/06/18
V_netFp networkPoint::getFacesSort2() const
{
	/*
	 *           *
	 *          /  \
	 *         / f11\
	 *  *-----*-------*-----*
	 *   \f1 /l0 f10/  \f9 /
	 *    \l1 f0\  / f8 \ /
	 *     *--l2-@-------*
	 *    / \f2 /  \f6 /  \
	 *   /f3 \ /f4  \ / f7 \
	 *  *-----*-------*-----*
	 *         \ f5 /
	 *          \  /
	 *            *
	 */
	V_netFp ret;
	netLp l;
	for (const auto &f : this->getFacesSort())
	{
		ret.emplace_back(f);
		l = f->getLinesFrom(this)[1];
		ret.emplace_back((*l)(f)); // nullptrかもしれない
	}
	return ret;
};
//
V_netPp networkPoint::getNeighborsSort2() const
{
	V_netPp ret;
	netPp p;
	netLp l;
	for (const auto &f : getFacesSort())
	{
		p = f->getPoints(this)[1];
		ret.emplace_back(p);
		l = f->getLineOpposite(this);
		p = (*l)(f)->getPointOpposite(l);
		ret.emplace_back(p);
	}
	return ret;
};
//
// using V_Tup_dddPp = std::vector<std::tuple<double, double, double, networkPoint *>>;

inline V_d networkPoint::getAngles() const
{
	//この点の周囲の角度を返す
	V_d ret(this->getFaces().size());
	int i = 0;
	for (const auto &f : this->getFaces())
		ret[i++] = f->getAngle(this);
	return ret;
};

inline V_d networkPoint::getAngles(networkLine *const base_line) const
{
	if (!base_line)
		return this->getAngles();
	//この点の周囲の角度を返す
	V_d ret(this->getFaces().size());
	int i = 0;
	for (const auto &f : this->getFacesSort(base_line))
		ret[i++] = f->getAngle(this);
	return ret;
};

using V_Tup_ddVdPp = std::vector<std::tuple<double, double, V_d, networkPoint *>>;
V_Tup_ddVdPp networkPoint::getNeighbors_Depth1_OnPolarAsTuple(networkLine *base_line = nullptr) const
{
	/**
	 * tupleは静的型しか持つことができない．
	 */
	auto fs = (base_line ? getFacesSort(base_line) : getFacesSort());
	netLp l;
	netPp p;
	double tmp;
	V_d angles = this->getAngles(base_line);
	double total_angles = Total(angles);
	double theta = 0;
	V_Tup_ddVdPp ret(fs.size());
	for (auto i = 0; i < fs.size(); i++)
	{
		p = fs[i]->getPoints(this)[1];
		auto r = Norm(p->getXtuple() - this->getXtuple());
		ret[i] = std::make_tuple(r * cos(theta), r * sin(theta), p->getX(), p);
		theta += 2. * M_PI * angles[i] / total_angles; // angleの最後は不要
	}
	return ret;
};

inline V_Tup_ddVdPp networkPoint::getNeighbors_Depth2_OnPolarAsTuple(networkLine *base_line = nullptr) const
{
	// 2021/10/26以前
	// /**
	//  * tupleは静的型しか持つことができない．
	//  */
	// V_Tup_ddVdPp ret;
	// auto fs = (base_line ? getFacesSort(base_line) : getFacesSort());
	// netLp l;
	// V_d angles = this->getAngles(base_line);
	// double total_angles = Total(angles);
	// //
	// double mean_len = 0.;
	// auto ls = this->getLines();
	// for (const auto &l : ls)
	// 	mean_len += l->length();
	// mean_len /= (double)ls.size();
	// //
	// double theta = 0., dt, r;
	// netPp p;
	// int i = 0;
	// for (const auto &f : fs)
	// {
	// 	p = f->getPoints(this)[1];
	// 	r = Norm3d(p->getX() - this->getX());
	// 	ret.emplace_back(std::make_tuple(r * cos(theta), r * sin(theta), p->getX(), p));
	// 	dt = 2. * M_PI * angles[i] / total_angles;
	// 	l = f->getLineOpposite(this);
	// 	p = (*l)(f)->getPointOpposite(l);
	// 	// 深さ２までの距離を計算　
	// 	auto ab = l->getPoints();
	// 	auto mid = ab[0]->getX() / 2. + ab[1]->getX() / 2.;
	// 	r = Norm3d(mid - this->getX());
	// 	r = r + Norm3d(p->getX() - mid);
	// 	ret.emplace_back(std::make_tuple(r * cos(theta + dt / 2), r * sin(theta + dt / 2), p->getX(), p));
	// 	theta += 2. * M_PI * angles[i] / total_angles; // angleの最後は不要
	// 	i++;
	// }
	// return ret;
	//! ------------------------------------------------------ */
	/**
	 * tupleは静的型しか持つことができない．
	 */
	V_Tup_ddVdPp ret;
	netLp l;
	V_d angles = this->getAngles(base_line);
	double total_angles = Total(angles);
	double theta = 0., dAngle, r;
	netPp p;
	int i = 0;
	for (const auto &f : (base_line ? getFacesSort(base_line) : getFacesSort()))
	{
		p = f->getPoints(this)[1];
		r = Norm(p->getXtuple() - this->getXtuple());
		ret.emplace_back(std::make_tuple(r * cos(theta), r * sin(theta), p->getX(), p));
		dAngle = 2. * M_PI * angles[i] / total_angles;
		l = f->getLineOpposite(this);
		if (!l->CORNER)
		{
			p = (*l)(f)->getPointOpposite(l);
			// 深さ２までの距離を計算
			auto [a, b] = l->getPointsTuple();
			auto mid = (a->getXtuple() + b->getXtuple()) / 2.;
			r = Norm(mid - this->getXtuple()) + Norm(p->getXtuple() - mid);
			ret.emplace_back(std::make_tuple(r * cos(theta + dAngle / 2), r * sin(theta + dAngle / 2), p->getX(), p));
		}
		theta += dAngle; // angleの最後は不要
		i++;
	}
	return ret;
};
inline V_Tup_ddVdPp networkPoint::getNeighbors_Depth2_OnPolarAsTuple_parametric(networkLine *base_line = nullptr) const
{
	double mean_len = 0.;
	auto ls = this->getLines();
	for (const auto &l : ls)
		mean_len += l->length();
	mean_len /= (double)ls.size();
	/* ------------------------- */
	V_Tup_ddVdPp ret;
	netLp l;
	V_d angles = this->getAngles(base_line);
	double total_angles = Total(angles);
	double theta = 0., dAngle, r;
	netPp p;
	int i = 0;
	for (const auto &f : (base_line ? getFacesSort(base_line) : getFacesSort()))
	{
		p = f->getPoints(this)[1];
		// r = Norm(p->getXtuple() - this->getXtuple());
		r = Norm(p->getXtuple() - this->getXtuple()) / mean_len;
		ret.emplace_back(std::make_tuple(r * cos(theta), r * sin(theta), p->getX(), p));
		dAngle = 2. * M_PI * angles[i] / total_angles;
		l = f->getLineOpposite(this);
		if (!l->CORNER)
		{
			p = (*l)(f)->getPointOpposite(l);
			// 深さ２までの距離を計算
			auto [a, b] = l->getPointsTuple();
			auto mid = (a->getXtuple() + b->getXtuple()) / 2.;
			r = Norm(mid - this->getXtuple()) + Norm(p->getXtuple() - mid);
			r /= mean_len;
			// r = 1 + Norm(p->getXtuple() - mid) / Norm(mid - this->getXtuple());
			ret.emplace_back(std::make_tuple(r * cos(theta + dAngle / 2), r * sin(theta + dAngle / 2), p->getX(), p));
		}
		theta += dAngle; // angleの最後は不要
		i++;
	}
	return ret;
};

// V_Tup_ddVdPp networkPoint::getNeighbors_Depth2_OnPolarAsTuple() const
// {
// 	/**
// 	 * tupleは静的型しか持つことができない．
// 	 */
// 	V_Tup_ddVdPp ret;
// 	auto fs = getFacesSort();
// 	netLp l;
// 	V_d angles(fs.size());
// 	double total_angles = 0.;
// 	//
// 	for (auto i = 0; i < fs.size(); i++)
// 	{
// 		angles[i] = fs[i]->getAngles(this)[0];
// 		total_angles += angles[i];
// 	}
// 	//
// 	double theta = 0., dt;
// 	netPp p;
// 	int i = 0;
// 	for (const auto &f : fs)
// 	{
// 		p = f->getPoints(this)[1];
// 		ret.emplace_back(std::make_tuple(cos(theta), sin(theta), p->getX(), p));
// 		theta += 2. * M_PI * angles[i] / total_angles; //angleの最後は不要
// 		dt = (angles[(i + 1) % fs.size()] - angles[i]);
// 		l = f->getLineOpposite(this);
// 		p = (*l)(f)->getPointOpposite(l);
// 		ret.emplace_back(std::make_tuple(2. * cos(theta + dt / 2), 2. * sin(theta + dt / 2), p->getX(), p));
// 		i++;
// 	}
// 	return ret;
// };

// V_Tup_ddVdPp networkPoint::getNeighbors_Depth2_OnPolarAsTuple() const
// {
// 	/**
// 	 * tupleは静的型しか持つことができない．
// 	 */
// 	V_Tup_ddVdPp ret;
// 	auto fs = getFacesSort();
// 	netLp l;
// 	netPp p;
// 	double theta, dt = 2. * M_PI / (fs.size());
// 	int i = 0;
// 	for (const auto &f : fs)
// 	{
// 		theta = dt * (i++);
// 		p = f->getPoints(this)[1];
// 		ret.emplace_back(std::make_tuple(cos(theta), sin(theta), p->getX(), p));
// 		l = f->getLineOpposite(this);
// 		p = (*l)(f)->getPointOpposite(l);
// 		ret.emplace_back(std::make_tuple(2. * cos(theta + dt / 2), 2. * sin(theta + dt / 2), p->getX(), p));
// 	}
// 	return ret;
// };

//!起点を指定できるようにする．
using V_Tup_ddVd = std::vector<std::tuple<double, double, V_d>>;
V_Tup_ddVd networkPoint::getNeighbors_Depth2_OnPolarAsTuple2() const
{
	/**
	 * tupleは静的型しか持つことができない．
	 */
	V_Tup_ddVd ret;
	auto fs = getFacesSort();
	netLp l;
	netPp p;
	double theta, dt = 2. * M_PI / (fs.size());
	for (auto i = 0; i < fs.size(); i++)
	{
		theta = dt * i;
		p = fs[i]->getPoints(this)[1];
		ret.emplace_back(std::make_tuple(cos(theta), sin(theta), p->getX()));

		auto q = fs[(i + 1) % fs.size()]->getPoints(this)[1];
		ret.emplace_back(std::make_tuple(cos(theta + dt / 2), sin(theta + dt / 2), (p->getX() + q->getX()) / 2.));

		l = fs[i]->getLineOpposite(this);
		p = (*l)(fs[i])->getPointOpposite(l);
		ret.emplace_back(std::make_tuple(2. * cos(theta + dt / 2), 2. * sin(theta + dt / 2), p->getX()));
	}

	return ret;
};
//
using V_Var_VVdVPp = std::vector<std::variant<VV_d, V_netPp>>;
V_Var_VVdVPp networkPoint::getNeighbors_Depth1_OnPolarAsVariant() const
{
	/**
	 * variantは，動的な可変長型の変数をもつものとして宣言できる．
	 */
	V_netPp VnetP;
	VV_d VVd;
	auto fs = getFacesSort();
	netLp l;
	netPp p;
	double theta, dt = 2. * M_PI / (fs.size());
	int i = 0;
	for (const auto &f : fs)
	{
		theta = dt * (i++);
		VVd.push_back({cos(theta), sin(theta)});
		VnetP.emplace_back(f->getPoints(this)[1]);
		// l = f->getLineOpposite(this);
		// p = (*l)(f)->getPointOpposite(l);
		// VVd.push_back({2. * cos(theta + dt / 2), 2. * sin(theta + dt / 2)});
		// VnetP.emplace_back(p);
	}
	return {VVd, VnetP};
};
V_Var_VVdVPp networkPoint::getNeighbors_Depth2_OnPolarAsVariant() const
{
	/**
	 * variantは，動的な可変長型の変数をもつものとして宣言できる．
	 */
	V_netPp VnetP;
	VV_d VVd;
	auto fs = getFacesSort();
	netLp l;
	netPp p;
	double theta, dt = 2. * M_PI / (fs.size());
	int i = 0;
	for (const auto &f : fs)
	{
		theta = dt * (i++);
		VVd.push_back({cos(theta), sin(theta)});
		VnetP.emplace_back(f->getPoints(this)[1]);
		l = f->getLineOpposite(this);
		p = (*l)(f)->getPointOpposite(l);
		VVd.push_back({2. * cos(theta + dt / 2), 2. * sin(theta + dt / 2)});
		VnetP.emplace_back(p);
	}
	return {VVd, VnetP};
};

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
inline void networkPoint::resetXinfo()
{
	if (this->xline != nullptr)
		this->xline->clearXPoints(this);
	this->xline = nullptr;
	if (this->xface != nullptr)
		this->xface->clearXPoints(this);
	this->xface = nullptr;
};
inline V_netFp networkPoint::getFaces_intersectQ(const bool TorF = true) const
{
	V_netFp ret;
	for (const auto &l : this->getLines())
		for (const auto &f : l->getFaces())
			if (f->intersectQ() == TorF)
			{
				if (std::find(ret.begin(), ret.end(), f) == ret.end())
					ret.emplace_back(f);
			}
	return ret;
};

// overload
inline void networkPoint::setBounds()
{
	try
	{
		// networkPointは，他のオブジェクトの座標の元となる座標情報を保持している．
		//そのため，networkPointが更新されれば，当然他のオブジェクトも更新されなければならない．
		//例えば，面の面積や法線方向など
		//  auto cbounds = CoordinateBounds({this->xyz});
		//  auto cbounds = CoordinateBounds({this->getX()});
		object3D::setBounds(geometry::CoordinateBounds(this->X));
		// object3D::setBounds(cbounds);
		for (const auto &l : this->getLines())
			l->setBounds();
		for (const auto &f : this->getFaces())
			f->setBounds();
		//! networkのsetBoundsをここで行うと，並列計算などで困ることになるだろう．点が全て動いてからboundsを設定しよう
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};

inline void networkPoint::setX(const Tddd &xyz_IN)
{
	try
	{
		object3D::setBounds(xyz_IN);
		for (const auto &l : this->getLines())
			l->setBounds();
		for (const auto &f : this->getFaces())
			f->setBounds();
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
// setXは，接続するlineやfaceのsetBoundsを行う．
inline void networkPoint::setX(const V_d &xyz_IN)
{
	this->setX(Tddd{xyz_IN[0], xyz_IN[1], xyz_IN[2]});
};
// setできない場合，元に戻すようにした．3月11日(木)
// mergeのために，実行せずチェックするだけの関数も作らなければ
//  inline bool networkPoint::setXcarefully(const V_d &xyz_IN)
//  {
//      auto current_X = this->getX();
//      object3D::setBounds(xyz_IN);
//      for (const auto &l : getLines())
//          if(!(l->setBounds()))
//          {Print("l->setBounds() 失敗した場合");
//              if(this->setXcarefully(current_X))
//                  return false;
//              else
//                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "戻そうとしたが失敗した");
//          }
//      for (const auto &f : getFaces())
//          if(!(f->setBounds()))
//          {Print("f->setBounds() 失敗した場合");
//              if(this->setXcarefully(current_X))
//                  return false;
//              else
//                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "戻そうとしたが失敗した");
//          }
//      //何もなければ正常に終了
//      return true;
//  };
inline bool networkPoint::setXcarefully(const V_d &xyz_IN)
{
	auto current_X = this->getX();
	object3D::setBounds(Tddd{xyz_IN[0], xyz_IN[1], xyz_IN[2]});
	for (const auto &l : getLines())
		if (!(l->setBounds()))
		{
			this->setX(current_X);
			return false;
		}
	for (const auto &f : getFaces())
		if (!(f->setBounds()))
		{
			this->setX(current_X);
			return false;
		}
	//何もなければ正常に終了
	return true;
};
inline Tddd networkPoint::getNormalTuple() const
{
	//角度の重みを掛けた法線ベクトルを足し合わせる
	Tddd normal = {0., 0., 0.};
	for (const auto &f : this->getFaces())
		normal += f->getAngle(this) * f->getNormalTuple();
	return Normalize(normal);

	// std::vector<Tddd> normals({});
	// for (const auto &f : this->getFaces())
	// 	normals.emplace_back(f->getNormalTuple());

	// if (normals.empty())
	// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "normals is empty");
	// else
	// 	return Mean(normals);
};
inline std::vector<double> networkPoint::getFaceAreas() const
{
	auto fs = this->getFaces();
	V_d ret(fs.size());
	for (auto i = 0; i < fs.size(); ++i)
		ret[i] = fs[i]->getArea();
	return ret;
};
inline Tddd networkPoint::getNormalAreaAveraged() const
{
	//角度の重みを掛けた法線ベクトルを足し合わせる
	Tddd normal = {0., 0., 0.};
	for (const auto &f : this->getFaces())
		normal += f->getArea() * f->getNormalTuple();
	return Normalize(normal);
};
inline V_d networkPoint::getNormal() const
{
	//角度の重みを掛けた法線ベクトルを足し合わせる
	auto normal = this->getNormalTuple();
	return ToVector(normal);
	// VV_d normals({});
	// for (const auto &f : this->getFaces())
	// 	normals.emplace_back(f->getNormal());

	// if (normals.empty())
	// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "normals is empty");
	// else
	// 	return Mean(normals);
};

/*SolidAngle_detail

  SolidAngle_detail*/
/*SolidAngle_detail_code*/
/*角点は周囲に点があるので，取得する点を制限しなければならない*/
inline double networkPoint::getSolidAngle(bool TorF)
{
	auto fs = this->getFaces();
	VV_d ns(fs.size(), V_d(3, 0.));
	int i = 0;
	for (const auto &f : fs)
		ns[i++] = f->getNormal();
	V_d mean = ns[0], var = {0., 0., 0.};
	for (auto i = 0; i < ns.size(); i++)
		var += (ns[i] - mean) * (ns[i] - mean);
	double s = std::sqrt(Norm(var));
	if (s < 1E-5)
		return 2. * M_PI;

	return geometry::SolidAngle(this->getX(), obj3D::extractX(this->getNeighborsSort(TorF)));
	;
};
inline double networkPoint::getSolidAngle(const V_netFp &faces /*available faces*/)
{
	V_netFp fs({});
	for (const auto &f : this->getFaces())
		if (MemberQ(faces, f))
			fs.emplace_back(f);
	if (fs.size() < 3)
		return 0;

	VV_d ns(fs.size(), V_d(3, 0.));
	int i = 0;
	for (const auto &f : fs)
		ns[i++] = f->getNormal();
	V_d mean = ns[0], var = {0., 0., 0.};
	for (auto i = 0; i < ns.size(); i++)
		var += (ns[i] - mean) * (ns[i] - mean);
	double s = std::sqrt(Norm(var));
	if (s < 1E-5)
		return 2. * M_PI;

	return geometry::SolidAngle(this->getX(), obj3D::extractX(fs));
	;
};
/*SolidAngle_detail_code*/
/*getNeighborsSort_detail
networkFaceのPointsが反時計周りに並んで保存されていれば，
その並びに沿って，`getNeighborsSort`も必ず反時計周りの多角形をなす点ベクトルを返す．
getNeighborsSort_detail*/
/*getNeighborsSort_sort*/
inline V_netPp networkPoint::getNeighborsSort() const
{
	try
	{
		VV_netPp Vps({});
		V_netPp qs(0);
		for (const auto &f : this->getFaces())
		{
			// qs = f->getPoints();
			// if (qs[0] == this)
			// 	Vps.push_back(V_netPp{qs[1], qs[2]});
			// else if (qs[1] == this)
			// 	Vps.push_back(V_netPp{qs[2], qs[0]});
			// else if (qs[2] == this)
			// 	Vps.push_back(V_netPp{qs[0], qs[1]});
			// else
			// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "ps is empty");
			//これは以下で実現できる．変更2021/06/18
			qs = f->getPoints(this);
			Vps.push_back(V_netPp{qs[1], qs[2]});
		}
		auto ret = FlattenAsChain(Vps);
		if (*ret.begin() == *ret.rbegin())
			ret.pop_back();
		else
		{
			return ret;
			//この場合，
			// * thisが端点であるか
			// * 周りが数珠つなぎとなっていない（エラー）
			// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not chain !");
		}
		return ret;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
//
inline V_netPp networkPoint::getNeighborsSort(bool TorF)
{
	try
	{
		VV_netPp Vps({});
		V_netPp qs(0);
		for (const auto &f : this->getFaces(TorF))
		{
			qs = f->getPoints();
			if (qs[0] == this)
				Vps.push_back(V_netPp{qs[1], qs[2]});
			else if (qs[1] == this)
				Vps.push_back(V_netPp{qs[2], qs[0]});
			else if (qs[2] == this)
				Vps.push_back(V_netPp{qs[0], qs[1]});
			else
				throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "ps is empty"));
		}
		auto ret = FlattenAsChain(Vps);
		if (*ret.begin() == *ret.rbegin())
			ret.pop_back();
		else
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not chain !");
		return ret;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
/*getNeighborsSort_sort*/
inline V_netPp networkPoint::getXNeighbors() const
{
	V_netPp ret({});
	for (const auto &l : this->getLines())
	{
		netPp closestXp = nullptr;
		double minlen = 1E+100;
		for (const auto &p : l->getXPoints())
			if (minlen > distance(this, p))
			{
				closestXp = p;
				minlen = distance(this, p);
			}
		if (closestXp)
			ret.emplace_back(closestXp);
	}
	return ret;
};

///////////// networkPoint //////////////
// inline networkPoint::networkPoint(Network *network_IN, const V_d &xyz_IN, networkLine *line_IN, networkFace *face_IN)
//     : networkObject(network_IN),
//       object3D(xyz_IN),
//       xline(line_IN),
//       xface(face_IN)
// {
//     // std::cout << "new networkPoint ... ";
//     networkObject::storage = network_IN; //なぜか初期化リストに入れれない
//     // std::cout << "done" << std::endl;
// };
// コンストラクタ
inline networkPoint::networkPoint(Network *network_IN,
								  Network *storage_IN,
								  const Tddd &xyz_IN,
								  networkLine *xline_IN,
								  networkFace *xface_IN)
	: networkObject(network_IN),
	  object3D(xyz_IN /*CoordinateBoundsに暗黙に置換される*/),
	  initialX(xyz_IN),
	  xline(xline_IN),
	  xface(xface_IN),
	  /* ------------------------------ */
	  force({0., 0., 0., 0., 0., 0.}),
	  inertia({0., 0., 0., 0., 0., 0.}),
	  acceleration({0., 0., 0., 0., 0., 0.}),
	  velocity({0., 0., 0., 0., 0., 0.}),
	  mass(0.),
	  /* ------------------------------------------------------ */
	  density(0.),
	  volume(0.),
	  radius(1.),
/* ------------------------------------------------------ */
#ifdef BEM
	  phiphin({0., 0.}),
	  grad_phi_BEM({0., 0., 0.}),
	  U_BEM({0., 0., 0.}),
	  U_update_BEM({0., 0., 0.}),
	  U_mod_BEM({0., 0., 0.}),
	  U_tangential_BEM({0., 0., 0.}),
	  U_normal_BEM({0., 0., 0.}),
	  laplacian_U_BEM({0., 0., 0.}),
	  normal_BEM({0., 0., 0.}),
#endif
	  /* ------------------------------------------------------ */
	  particlize_info({nullptr, {nullptr, nullptr, nullptr}, {0., 0.}, 0., 0.})
{

	// std::cout << "new networkPoint ... ";

	if (xline != nullptr)
		xline->addXPoint(this);
	if (xface != nullptr)
		xface->addXPoint(this);

	this->storage = storage_IN; //なぜか初期化リストに入れれない
	this->storage->add(this);
	// std::cout << "done" << std::endl;
	/* ------------------------------------------------------ */

#ifdef DEM
	this->contactP = {};
	this->mass = 1.;
	this->radius = 1.;
	this->W = 0.;
	//! SPH
	this->a_viscosity = 0.;
	this->mu_SPH = 0.001005;
	this->DUDt_SPH = {0., 0., 0.};
	this->DrhoDt_SPH = 0.;
	this->density = 0.;
	this->lap_U = {0., 0., 0.};
	this->div_U = 0;
	this->gradP_SPH = {0., 0., 0.};
	this->U_SPH = {0., 0., 0.};
	this->mu_lap_rho_g_SPH = {0., 0., 0.};
	this->tmp_U_SPH = {0., 0., 0.};
	this->interpolated_normal_SPH = {0., 0., 0.};
	this->pn_is_set = false;
	this->radius_SPH = 0.;
#endif

	this->Dirichlet = false;
	this->Neumann = false;
};
//逆方向の立体角
inline double networkPoint::getSolidAngle()
{
	auto fs = this->getFaces();
	// VV_d ns(fs.size(), V_d(3, 0.));
	// int i = 0;
	// for (const auto &f : fs)
	//   ns[i++] = f->getNormal();
	// V_d mean = ns[0], var = {0., 0., 0.};
	// for (auto i = 0; i < ns.size(); i++)
	//   var += (ns[i] - mean) * (ns[i] - mean);
	// double s = std::sqrt(Norm(var));
	// if (s < 1E-5)
	//   return 2. * M_PI;

#define check_solidangle
#ifdef check_solidangle

	//  mk_vtu("./vtu/getFaces.vtu", fs);
	auto ps = this->getNeighborsSort();
	//  mk_vtu("./vtu/getFaces_sortpoint.vtu", {ps});
	auto xyz = obj3D::extractX(ps);
	VV_d x_on_sphere({});
	for (const auto &x : xyz)
		x_on_sphere.emplace_back((x - this->getX()) / Norm(x - this->getX()));
	auto normal = -this->getNormal() / Norm(this->getNormal());
	double total = 0.;
	int sz = xyz.size();
	for (auto i = 0; i < sz; i++)
	{
		auto tmp = geometry::SolidAngle({0., 0., 0.}, {x_on_sphere[(i + sz) % sz], x_on_sphere[(i + sz + 1) % sz], normal});
		if (Between(tmp, {0., 4. * M_PI}))
			total += tmp;
	}
	return total;
	// std::cout << "total = " << total << std::endl;
	// std::cout << "solid = " << geometry::SolidAngle(this->getX(), obj3D::extractX(this->getNeighborsSort())) << std::endl;
	// std::cin.ignore();
#endif
	return geometry::SolidAngle(this->getX(), obj3D::extractX(this->getNeighborsSort()));
};
/////////////////////////////////////////

std::vector<bool> isIntxns(const V_netLp &ls)
{
	std::vector<bool> ret(ls.size(), false);
	for (auto i = 0; i < ls.size(); i++)
		if (ls[i]->isIntxn())
			ret[i] = true;
	return ret;
};

std::vector<bool> isXPoints(const V_netPp &ps)
{
	std::vector<bool> ret(ps.size(), false);
	for (auto i = 0; i < ps.size(); i++)
		if (ps[i]->isXPoint())
			ret[i] = true;
	return ret;
};

std::vector<std::vector<bool>> isXPoints(const VV_netPp &ps)
{
	std::vector<std::vector<bool>> ret(ps.size());
	for (auto i = 0; i < ps.size(); i++)
		ret[i] = isXPoints(ps[i]);
	return ret;
};

V_netLp getXLines(const V_netPp &ps)
{
	V_netLp ret(ps.size());
	for (auto i = 0; i < ps.size(); i++)
		ret[i] = ps[i]->getXLine();
	return ret;
};

VV_netLp getXLines(const VV_netPp &ps)
{
	VV_netLp ret(ps.size());
	for (auto i = 0; i < ps.size(); i++)
		ret[i] = getXLines(ps[i]);
	return ret;
};

#endif