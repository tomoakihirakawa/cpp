
// できるだけカプセル化することが，今後の変更を最小に抑えてくれる．
// この観点から幾何学関連の操作（面の抽出や三角分割など）は，bemには持ち込まない．
#ifndef BEM_H
#define BEM_H

#include "InterpolationRBF.hpp"
#include "Network.hpp"
#include "svd.hpp"

namespace BEM {
using V_d = std::vector<double>;
using V_i = std::vector<int>;
using VV_d = std::vector<V_d>;
using VVV_d = std::vector<VV_d>;

using netFp = networkFace *;

using netP = networkPoint;
using netPp = networkPoint *;
using V_netPp = std::vector<networkPoint *>;
using VV_netPp = std::vector<V_netPp>;

using netF = networkFace;
using V_netFp = std::vector<netF *>;

using V_Netp = std::vector<Network *>;

/*DN_detail


DNは，`Network *`を以下の種類に分けて保存・整理しておくためのクラス．

* Base
* Corner
* Dirichlet
* Neumann

DN_detail*/
/*DN_code*/
// class DN
// {
// 	using Netp = Network *;
// 	using V_Netp = std::vector<Netp>;
// 	using VV_Netp = std::vector<V_Netp>;

// public:
// 	int Bindex = -2;
// 	int Cindex = -1;
// 	int Dindex = 0;
// 	int Nindex = 1;

// 	V_Netp B, D, N, C;

// 	DN() : B(0), D(0), N(0), C(0){};
// 	~DN(){};

// 	void clear()
// 	{
// 		this->B.clear();
// 		this->D.clear();
// 		this->N.clear();
// 		this->C.clear();
// 	};

// 	void set(const int i, const Netp &net)
// 	{
// 		if (i == Dindex)
// 			network::add(this->D, net);
// 		else if (i == Nindex)
// 			network::add(this->N, net);
// 		else if (i == Cindex)
// 			network::add(this->C, net);
// 		else if (i == Bindex)
// 			network::add(this->B, net);
// 		else
// 			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "index must be -1, 0, or 1"));
// 	};

// 	void set(const int i, const V_Netp &nets)
// 	{
// 		for (const auto &net : nets)
// 			this->set(i, net);
// 	};

// 	void setB(const V_Netp &nets) { set(Bindex, nets); };
// 	void setC(const V_Netp &nets) { set(Cindex, nets); };
// 	void setD(const V_Netp &nets) { set(Dindex, nets); };
// 	void setN(const V_Netp &nets) { set(Nindex, nets); };
// 	void setB(const Netp &nets) { set(Bindex, {nets}); };
// 	void setC(const Netp &nets) { set(Cindex, {nets}); };
// 	void setD(const Netp &nets) { set(Dindex, {nets}); };
// 	void setN(const Netp &nets) { set(Nindex, {nets}); };

// 	V_Netp getAll() const
// 	{
// 		VV_Netp tmp = {this->B, this->C, this->D, this->N};
// 		return DeleteDuplicates(Flatten(tmp));
// 	}

// 	V_Netp get(int i)
// 	{
// 		if (i == Dindex)
// 			return this->D;
// 		else if (i == Nindex)
// 			return this->N;
// 		else if (i == Cindex)
// 			return this->C;
// 		else if (i == Bindex)
// 			return this->B;
// 		else
// 			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "index must be -1, 0, or 1"));
// 	};

// 	int index(const Netp net) const
// 	{
// 		if (MemberQ(this->D, net))
// 			return this->Dindex;
// 		else if (MemberQ(this->N, net))
// 			return this->Nindex;
// 		else if (MemberQ(this->C, net))
// 			return this->Cindex;
// 		else if (MemberQ(this->B, net))
// 			return this->Bindex;
// 		else
// 			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "index must be -1, 0, or 1"));
// 	};

// 	int index(const netPp &p) const { return this->index(p->getNetwork()); };
// 	int isD(const netPp &p) const { return MemberQ(this->D, p->getNetwork()); };
// 	int isN(const netPp &p) const { return MemberQ(this->N, p->getNetwork()); };
// 	int isC(const netPp &p) const { return MemberQ(this->C, p->getNetwork()); };
// 	int isB(const netPp &p) const { return MemberQ(this->B, p->getNetwork()); };
// };
/*DN_code*/

// V_Netp TakeBaseNetwork(V_Netp nets_IN)
// {
// 	V_Netp ret;
// 	for (auto &n : nets_IN)
// 		if (n->isB())
// 			ret.emplace_back(n);
// 	return ret;
// };
//////////////////////////////////////
using map_P_Vd = std::map<netPp, V_d>;
using map_P_d = std::map<netPp, double>;

// using uo_map_P_d = std::unordered_map<netP *, double>;
// using uo_map_P_Vd = std::unordered_map<netP *, V_d>;

// //----------------------------------------------------------------
// V_netPp searchPoints(const int depth, const netPp p, const V_Netp &nets = {})
// {
// 	depth_searcher<networkPoint> S(depth);
// 	S.clear();
// 	S.set(p);
// 	for (const auto &n : nets)
// 		S.addNetwork(n);

// 	//これはsearch(true)であるべき．微分は自身の値があった方が精度がいい
// 	S.search(true);
// 	return ToVector(S.getObjectsUO());
// };
// //----------------------------------------------------------------
// V_netFp searchFaces(const int depths, const netPp p, const V_Netp &nets = {})
// {
// 	V_netFp fs({});
// 	for (const auto &p : searchPoints(depths, p, nets))
// 		for (const auto &f : p->getFaces())
// 			if (MemberQ(nets, f->getNetwork()))
// 				fs.emplace_back(f);
// 	return DeleteDuplicates(fs);
// };
// ////////////////////////////////////////////////////////////////////////////////////////
// V_netPp searchPoints(const V_i &depths, const netPp p, const V_netPp &refp, const V_Netp &nets = {})
// {
// 	V_netPp ret;
// 	depth_searcher<networkPoint> S;
// 	for (int d = depths[0]; d <= depths[1]; d++)
// 	{
// 		// depth_searcher<networkPoint> S(d);
// 		S.clear();
// 		S.setDepth(d);
// 		S.set(p);
// 		for (const auto &n : nets)
// 			S.addNetwork(n);

// 		S.search(true); // must be true

// 		ret = Intersection(refp, S.getObjects());

// 		if (ret.size() > 1)
// 			return ret;

// 		if (d == depths[1])
// 		{
// 			std::cout << message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not search enough points, depth = " + std::to_string(S.depth)) << std::endl;
// 			// Print(depths, red);
// 			// mk_vtu("./vtu/can_not_search_enough_points_getReachedLines.vtu", S.getReachedLines());
// 			// mk_vtu("./vtu/can_not_search_enough_points_getEnteredLines.vtu", S.getEnteredLines());
// 			// mk_vtu("./vtu/can_not_search_enough_points.vtu", {S.getObjects()});
// 			// mk_vtu("./vtu/can_not_search_enough_points_refp.vtu", {refp});
// 			// mk_vtu("./vtu/can_not_search_enough_points_ret.vtu", {ret});
// 			// std::cout << "p->getNeighbors() = " << p->getNeighbors() << std::endl;
// 			// std::cout << "p->getLines() = " << p->getLines() << std::endl;
// 			// std::cout << "S.getReachedLines().size() = " << S.getReachedLines().size() << std::endl;
// 			// std::cout << "S.getEnteredLines().size() = " << S.getEnteredLines().size() << std::endl;
// 			// std::cout << "S.getObjects().size() = " << S.getObjects().size() << std::endl;
// 			// std::cout << "S.getObjects_().size() = " << S.getObjects_().size() << std::endl;
// 			// std::cout << "S.getNetworks().size() = " << S.getNetworks().size() << std::endl;
// 			// std::cout << nets << std::endl;
// 			// for (auto n : nets)
// 			//   std::cout << n->getName() << std::endl;
// 		}
// 	}
// 	// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not search enough points");
// 	return {};
// };
// // //////////////////
// class takeToDepth
// {
// public:
// 	VV_d X;
// 	V_d V;
// 	VV_d VV;
// 	V_netPp points;
// 	netPp origin;

// 	V_netPp takefirstpoints(const V_netPp &ps, const int maxnum = 0)
// 	{
// 		V_netPp ret = {};
// 		int count = 0;
// 		for (const auto p : ps)
// 		{
// 			ret.emplace_back(p);
// 			if (maxnum != 0 && ++count >= maxnum)
// 				break;
// 		}
// 		return ret;
// 	};

// 	takeToDepth(const V_i &depths, const netPp p, const V_netPp &P, const V_Netp &nets = {}, const int maxnum = 0)
// 		: X({}), V({}), VV({}), points({}), origin(p)
// 	{
// 		auto tmp = searchPoints(depths, p, P, nets);
// 		sortByDistance(tmp, origin);
// 		if (maxnum != 0)
// 			this->points = takefirstpoints(tmp, maxnum);
// 		else
// 			this->points = tmp;

// 		for (const auto &q : P)
// 			if (MemberQ(this->points, q))
// 			{
// 				this->X.emplace_back(ToVector(q->getXtuple()));
// 			}
// 	};
// 	takeToDepth(const V_i &depths, const netPp p, const map_P_d &P_d, const V_Netp &nets = {}, const int maxnum = 0)
// 		: X({}), V({}), VV({}), points({}), origin(p)
// 	{
// 		auto tmp = searchPoints(depths, p, TakeFirst(P_d), nets);
// 		sortByDistance(tmp, origin);
// 		if (maxnum != 0)
// 			this->points = takefirstpoints(tmp, maxnum);
// 		else
// 			this->points = tmp;

// 		for (const auto &[q, d] : P_d)
// 			if (MemberQ(this->points, q))
// 			{
// 				this->X.emplace_back(ToVector(q->getXtuple()));
// 				this->V.emplace_back(d);
// 			}
// 	};
// 	takeToDepth(const V_i &depths, const netPp p, const map_P_Vd &P_V, const V_Netp &nets = {}, const int maxnum = 0)
// 		: X({}), V({}), VV({}), points({}), origin(p)
// 	{
// 		auto tmp = searchPoints(depths, p, TakeFirst(P_V), nets);
// 		sortByDistance(tmp, origin);
// 		if (maxnum != 0)
// 			this->points = takefirstpoints(tmp, maxnum);
// 		else
// 			this->points = tmp;

// 		for (const auto &[q, v] : P_V)
// 			if (MemberQ(this->points, q))
// 			{
// 				this->X.emplace_back(ToVector(q->getXtuple()));
// 				this->VV.emplace_back(v);
// 			}
// 	};
// 	// takeToDepth(const V_i &depths, const netPp p, const uo_map_P_Vd &P_V, const V_Netp &nets = {}, const int maxnum = 0)
// 	//     : X({}), V({}), VV({}), points({}), origin(p)
// 	// {
// 	//   auto tmp = searchPoints(depths, p, TakeFirst(P_V), nets);
// 	//   sortByDistance(tmp, origin);
// 	//   if (maxnum != 0)
// 	//     this->points = takefirstpoints(tmp, maxnum);
// 	//   else
// 	//     this->points = tmp;

// 	//   for (const auto &[q, v] : P_V)
// 	//     if (MemberQ(this->points, q))
// 	//     {
// 	//       this->X.emplace_back(q->getX());
// 	//       this->VV.emplace_back(v);
// 	//     }
// 	// };
// };

/////////////////////////////////////////////////////////////////
// #include "CompGrid.hpp"

///////////////////////////////////////////////////////////////////
/// USE LIKE Dot( N(ab), sample )
// V_d N(const V_d &ab) { return {ab[0], ab[1], 1. - (ab[0] + ab[1])}; };

// V_d dNd_(const int i)
// {
// 	switch (i)
// 	{
// 	case 0:
// 		return {1., 0., -1.};
// 	case 1:
// 		return {0., 1., -1.};
// 	default:
// 		throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "0 or 1"));
// 		return {0., 0., 0.};
// 	}
// };

// VV_d dNd_() { return {{1., 0., -1.}, {0., 1., -1.}}; };

//--------------------
// VV_d M(const netFp &f)
// {
// 	auto xyzs = extractX(f->getPoints());
// 	auto dXds0 = Dot(dNd_(0), xyzs);
// 	auto dXds1 = Dot(dNd_(1), xyzs);
// 	return {dXds0 / Norm(dXds0) /*s*/,
// 			dXds1 / Norm(dXds1) /*m*/,
// 			ToVector(Normalize(f->normal)) /*n*/};
// };
//--------------------
/*global velocity*/
// V_d vG(const netFp &f, const V_d &Phi, const double phin)
// {
// 	return Dot(Inverse(M(f)), {Dot(dNd_(0), Phi),
// 							   Dot(dNd_(1), Phi),
// 							   phin} /*local velocity*/);
// };
//--------------------
// V_d meanvG(const netPp p, map_P_Vd &P_phiphin, const V_Netp &nets, const V_netFp &faces = {})
// {
// 	Print("meanvG");
// 	VV_d VvG = {};
// 	V_netPp ps;
// 	V_d Phi;
// 	bool found = false;
// 	V_d phiphin0, phiphin1, phiphin2;
// 	for (const auto &f : p->getFaces())
// 	{ // need to check
// 		ps = f->getPoints();
// 		Print(ps);
// 		if (faces.empty() || MemberQ(faces, f))
// 			if (ps.size() == 3 && AllMemberQ(P_phiphin, ps) && MemberQ(nets, f->getNetwork()))
// 			{
// 				phiphin0 = P_phiphin[ps[0]];
// 				phiphin1 = P_phiphin[ps[1]];
// 				phiphin2 = P_phiphin[ps[2]];

// 				Phi = {phiphin0[0], phiphin1[0], phiphin2[0]};
// 				double phin = (phiphin0[1] + phiphin1[1] + phiphin2[1]) / 3.;
// 				VvG.emplace_back(vG(f, Phi, phin));
// 				found = true;
// 			}
// 		// Print(VvG);
// 	}
// 	if (VvG.empty())
// 		return {};
// 	else
// 		return Mean(VvG);
// };
///////////////////

V_d parameterize(const V_netPp &points, const netPp &p) {
   if (p == points[0])
      return {1., 0., 0.};
   else if (p == points[1])
      return {0., 1., 0.};
   else if (p == points[2])
      return {0., 0., 1.};

   Print(__PRETTY_FUNCTION__);
   abort();

   return {};
};
/////////////////
// class searchPoints_X_phi_phin
// {
// public:
// 	VV_d X;
// 	V_d phi;
// 	V_d phin;
// 	V_netPp points;
// 	//
// 	searchPoints_X_phi_phin(const V_i &depths, const netPp p, const map_P_Vd &sampleP_phiphin, const V_Netp &nets = {})
// 		: X({}), phi({}), phin({}), points({})
// 	{
// 		this->points = searchPoints(depths, p, TakeFirst(sampleP_phiphin), nets);
// 		for (const auto &[q, phiphin] : sampleP_phiphin)
// 			if (MemberQ(this->points, q))
// 			{
// 				this->X.emplace_back(ToVector(q->getXtuple()));
// 				this->phi.emplace_back(phiphin[0]);
// 				this->phin.emplace_back(phiphin[1]);
// 			}
// 	};
// };
//----------------------------------------------------------------
// V_d AverageNearNeighbors_phiphin(const netPp p, const map_P_Vd &sampleP_phiphin, const V_Netp &nets)
// {
// 	if (p->getLines().empty())
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	takeToDepth s({5, 10}, p, sampleP_phiphin, nets, 40);
// 	auto phiphin = Transpose(s.VV);
// 	return {Mean(phiphin[0]), Mean(phiphin[1])};
// };
//----------------------------------------------------------
// V_d InterpolationRBF_phiphin(const netPp p, const map_P_Vd &P_phiphin, const V_Netp &nets, const int num = 40)
// {
// 	if (p->getLines().empty())
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

// 	takeToDepth s({5, 10}, p, P_phiphin, nets, num);
// 	auto phiphin = Transpose(s.VV);
// 	auto interp0 = InterpolationRBF(s.X, phiphin[0], ToVector(p->getXtuple()));
// 	auto interp1 = InterpolationRBF(s.X, phiphin[1], ToVector(p->getXtuple()));
// 	return {interp0(ToVector(p->getXtuple())), interp1(ToVector(p->getXtuple()))};
// };

// double InterpolationRBF_phin(const netPp p, map_P_Vd &P_phiphin, const V_Netp &nets, const int num = 40)
// {
// 	if (p->getLines().empty())
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

// 	// takeToDepth s({4, 10}, p, P_phiphin, nets);
// 	takeToDepth s({5, 10}, p, P_phiphin, nets, num);
// 	auto interp = InterpolationRBF(s.X, Transpose(s.VV)[1] /*phiphin*/, ToVector(p->getXtuple()));
// 	return interp(ToVector(p->getXtuple()));
// };

// double InterpolationRBF_phi(const netPp p, map_P_Vd &P_phiphin, const V_Netp &nets, const int num = 40)
// {
// 	Print(message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
// 	if (p->getLines().empty())
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

// 	// takeToDepth s({4, 10}, p, P_phiphin, nets);
// 	takeToDepth s({5, 10}, p, P_phiphin, nets, num);

// 	auto phiphin = Transpose(s.VV);
// 	auto interp = InterpolationRBF(s.X, phiphin[0], ToVector(p->getXtuple()));
// 	return interp(ToVector(p->getXtuple()));
// };
// double InterpolationRBF_phi(const netPp p, map_P_Vd &P_phiphin, Network *net)
// {
// 	std::vector<Network *> nets = {net};
// 	return InterpolationRBF_phi(p, P_phiphin, nets);
// };
// //----------------------------------------------------------
// V_d IDW_phiphin(const netPp p, const map_P_Vd &sampleP_phiphin, const V_Netp &nets, const double powerIDW = 2.)
// {
// 	if (p->getLines().empty())
// 	{
// 		// return {0., 0.};
// 		mk_vtu("./vtu/p->getLines().empty()).vtu", {{p}});
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty())");
// 	}
// 	takeToDepth s({5, 10}, p, sampleP_phiphin, nets, 30);
// 	auto phiphin = Transpose(s.VV);
// 	auto interp0 = InterpolationIDW(s.X, phiphin[0], powerIDW);
// 	auto interp1 = InterpolationIDW(s.X, phiphin[1], powerIDW);
// 	return {interp0(p->getX()), interp1(p->getX())};
// };
// double IDW_phi(const netPp p, const map_P_Vd &sampleP_phiphin, const V_Netp &nets, const double powerIDW = 2.)
// {
// 	if (p->getLines().empty())
// 	{
// 		// return {0., 0.};
// 		mk_vtu("./vtu/p->getLines().empty()).vtu", {{p}});
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty())");
// 	}
// 	takeToDepth s({5, 10}, p, sampleP_phiphin, nets, 30);
// 	auto phiphin = Transpose(s.VV);
// 	auto interp0 = InterpolationIDW(s.X, phiphin[0], powerIDW);
// 	return interp0(p->getX());
// };
// double IDW_phin(const netPp p, const map_P_Vd &sampleP_phiphin, const V_Netp &nets, const double powerIDW = 2.)
// {
// 	if (p->getLines().empty())
// 	{
// 		// return {0., 0.};
// 		mk_vtu("./vtu/p->getLines().empty()).vtu", {{p}});
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty())");
// 	}
// 	takeToDepth s({5, 10}, p, sampleP_phiphin, nets, 30);
// 	auto phiphin = Transpose(s.VV);
// 	auto interp0 = InterpolationIDW(s.X, phiphin[1], powerIDW);
// 	return interp0(p->getX());
// };
// //-----------------------------------------------------------
// // ここで，phinがつかわれるので，BIEでphinが得られていない外部の点のnablaPhiは，ランダムで正しくない.
// V_d InterpolationRBFvG(const netPp p, map_P_Vd &P_phiphin, const V_Netp &nets)
// {
// 	if (p->getLines().empty())
// 	{
// 		// return {0., 0.};
// 		mk_vtu("./vtu/p->getLines().empty()).vtu", {{p}});
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty())");
// 	}

// 	takeToDepth s({6, 10}, p, P_phiphin, nets, 40);

// 	auto phiphin = Transpose(s.VV);
// 	auto interp = InterpolationRBF(s.X, phiphin[0], p->getX());
// 	auto normal = p->getNormalFromSameBC();
// 	auto v = interp.grad(p->getX());

// 	return v + (-Dot(v, normal) + P_phiphin[p][1]) * normal;
// };
////
// V_d InterpolationRBFvG(const netPp p, uo_map_P_Vd &P_phiphin, const V_Netp &nets)
// {
//   if (p->getLines().empty())
//   {
//     // return {0., 0.};
//     mk_vtu("./vtu/p->getLines().empty()).vtu", {{p}});
//     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty())");
//   }

//   takeToDepth s({6, 10}, p, P_phiphin, nets, 40);

//   auto phiphin = Transpose(s.VV);
//   auto interp = InterpolationRBF(s.X, phiphin[0], p->getX());
//   auto normal = p->getNormalFromSameBC();
//   auto v = interp.nabla(p->getX());

//   return v + (-Dot(v, normal) + P_phiphin[p][1]) * normal;
// };

// V_d InterpolationIDWvG(const netPp p, map_P_Vd &P_phiphin, const V_Netp &nets, const double pp)
// {
// 	if (p->getLines().empty())
// 	{
// 		// return {0., 0.};
// 		mk_vtu("./vtu/p->getLines().empty()).vtu", {{p}});
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty())");
// 	}

// 	takeToDepth s({2, 5}, p, P_phiphin, nets);
// 	auto phiphin = Transpose(s.VV);
// 	auto interp = InterpolationIDW(s.X, phiphin[0], pp);
// 	auto normal = p->getNormalFromSameBC();
// 	auto v = interp.nabla(p->getX());
// 	return v + (-Dot(v, normal) + P_phiphin[p][1]) * normal;
// };

// ///////////////////

// /////---------
// V_d global_velocity(const netPp p, map_P_Vd &P_phiphin, const V_Netp &nets, const V_netFp &faces = {})
// {
// 	auto vg = meanvG(p, P_phiphin, nets, faces);
// 	if (vg.empty())
// 		return InterpolationRBFvG(p, P_phiphin, nets);
// 	else
// 		return vg;
// };

// //////////////////

// map_P_d DphiDtRBF(map_P_Vd &P_phiphin, const V_netFp &Faces, const V_Netp &nets)
// {
// 	Print("DphiDtRBF", red);
// 	map_P_d ret;
// 	V_d u = {0., 0., 0.};
// 	for (const auto &p : TakeFirst(P_phiphin))
// 	{
// 		u = InterpolationRBFvG(p, P_phiphin, {p->getNetwork()});
// 		ret[p] = Dot(u, u) / 2. - p->getX()[2];
// 	}
// 	return ret;
// };

// map_P_Vd nablaPhiRBF(map_P_Vd &P_phiphin, std::vector<Network *> nets = {})
// {
// 	map_P_Vd ret;
// 	for (const auto &p : TakeFirst(P_phiphin))
// 		ret[p] = InterpolationRBFvG(p, P_phiphin, {p->getNetwork()});
// 	return ret;
// };

// map_P_Vd nablaPhiRBF(map_P_Vd &P_phiphin, Network *net)
// {
// 	return nablaPhiRBF(P_phiphin, std::vector<Network *>{net});
// };

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

bool isConverged(std::map<netP *, V_d> &new_igign, std::map<netP *, V_d> &pre_igign, const V_netPp &ps, const double conv) {
   V_d tmp, pre;
   double t;
   for (const auto &p : ps) {
      pre = pre_igign[p];
      tmp = new_igign[p];
      if ((t = std::abs(tmp[0])) > 1E-20)
         if (std::abs(tmp[0] - pre[0]) / t > conv)
            return false;

      if ((t = std::abs(tmp[1])) > 1E-20)
         if (std::abs(tmp[1] - pre[1]) / t > conv)
            return false;
   }
   return true;
};
/* ------------------------------------------------------ */
bool isConverged(std::unordered_map<netP *, Tdd> &new_igign, std::unordered_map<netP *, Tdd> &pre_igign, const V_netPp &ps, const double conv) {
   Tdd tmp, pre;
   double t;
   for (const auto &p : ps) {
      pre = pre_igign[p];
      tmp = new_igign[p];
      if ((t = std::abs(std::get<0>(tmp))) > 1E-20)
         if (std::abs(std::get<0>(tmp) - std::get<0>(pre)) / t > conv)
            return false;

      if ((t = std::abs(std::get<1>(tmp))) > 1E-20)
         if (std::abs(std::get<1>(tmp) - std::get<1>(pre)) / t > conv)
            return false;
   }
   return true;
};

const V_d zeros(2, 0.);
/* ------------------ 三角形面の2次のラグランジュ補間 ------------------ */
VVV_d GWGWfrom0to1 =
    {
        VV_d{},
        GaussianQuadratureWeights(1, 0., 1.),
        GaussianQuadratureWeights(2, 0., 1.),
        GaussianQuadratureWeights(3, 0., 1.),
        GaussianQuadratureWeights(4, 0., 1.),
        GaussianQuadratureWeights(5, 0., 1.),
        GaussianQuadratureWeights(6, 0., 1.),
        GaussianQuadratureWeights(7, 0., 1.),
        GaussianQuadratureWeights(8, 0., 1.),
        GaussianQuadratureWeights(9, 0., 1.),
        GaussianQuadratureWeights(10, 0., 1.),
        GaussianQuadratureWeights(11, 0., 1.),
        GaussianQuadratureWeights(12, 0., 1.),
        GaussianQuadratureWeights(13, 0., 1.),
        GaussianQuadratureWeights(14, 0., 1.),
        GaussianQuadratureWeights(15, 0., 1.),
        GaussianQuadratureWeights(16, 0., 1.),
        GaussianQuadratureWeights(17, 0., 1.),
        GaussianQuadratureWeights(18, 0., 1.),
        GaussianQuadratureWeights(19, 0., 1.),
        GaussianQuadratureWeights(20, 0., 1.),
        GaussianQuadratureWeights(21, 0., 1.),
        GaussianQuadratureWeights(22, 0., 1.),
        GaussianQuadratureWeights(23, 0., 1.),
        GaussianQuadratureWeights(24, 0., 1.),
        GaussianQuadratureWeights(25, 0., 1.),
        GaussianQuadratureWeights(26, 0., 1.)};

// class checkSurfaceInterp
// {
// 	// たぶん範囲を間違っている
// public:
// 	double area;
// 	double solution_intgration_out;
// 	double solution_intgration_in;
// 	double solution_intgration_surface;
// 	double solution_intgration_surface_simplified;
// 	double solution_intgration_surface_simplified_small;
// 	double solution_intgration_surface_simplified_face;
// 	V_d solution_intgration_suface_from_eps_1;
// 	double r_1_intgration_surface;
// 	double r_1_intgration_surface_face;
// 	checkSurfaceInterp(interpolationTriangle *intp,
// 					   const netPp origin,
// 					   const V_netFp &fs,
// 					   int gausspoints,
// 					   int Npoints = 12,
// 					   bool checkSing = false)
// 		: area(0.),
// 		  solution_intgration_out(0.),
// 		  solution_intgration_in(0.),
// 		  solution_intgration_surface(0.),
// 		  solution_intgration_surface_simplified(0.),
// 		  solution_intgration_surface_simplified_face(0.),
// 		  solution_intgration_suface_from_eps_1(100, 0.),
// 		  r_1_intgration_surface(0.),
// 		  r_1_intgration_surface_face(0.),
// 		  solution_intgration_surface_simplified_small(0.)
// 	{
// 		/* ---------------------- 原点の決定 -------------------------------- */
// 		V_netPp ps;
// 		V_d in_point = {0., 0., 0.};
// 		V_d out_point = {1E+3, 1E+3, 1E+3};
// 		/* ------------------------------------------------------ */
// 		V_d a(3, 0.);
// 		double weight = 0;
// 		V_d X(3, 0.);
// 		V_d r(3, 0.);
// 		V_d grad_phi(3, 0.);
// 		double t0, t1, J, w0, w1;
// 		double x_sing_t0 = 1 / 2.;
// 		auto beta = 4.;
// 		/* ------------------------------------------------------ */
// 		// VVV_d check;
// 		auto originface = origin->getFaces()[0];
// 		intp->set(extractX((originface)->get12Points(origin)));
// 		V_d oFace = (*intp)(0.5, 0.5);
// 		auto gwgw = GWGWfrom0to1[gausspoints];
// 		auto sing_gwgw = GaussianQuadratureWeights(gausspoints, InvSg(0., x_sing_t0, beta), InvSg(1., x_sing_t0, beta));
// 		for (const auto &f : fs)
// 		{
// 			if (Npoints == 12)
// 				ps = f->get12Points(origin);
// 			else if (Npoints == 6)
// 				ps = f->get6PointsTuple(origin);

// 			// check.emplace_back(extractX(ps));
// 			intp->set(extractX(ps) /*get12Points()で得られた点*/);
// 			// bool isSingular = (checkSing && ps.size() > 6 && ps[4] == origin);
// 			bool isSingular = (originface == f);
// 			for (const auto &x0w0 : isSingular ? sing_gwgw : gwgw)
// 			{
// 				t0 = isSingular ? Sg(x0w0[0], x_sing_t0, beta) : /**/ x0w0[0];
// 				w0 = isSingular ? x0w0[1] * DSg(x0w0[0], x_sing_t0, beta) : /**/ x0w0[1];
// 				for (const auto &x1w1 : isSingular ? sing_gwgw : gwgw)
// 				{
// 					t1 = isSingular ? Sg(x1w1[0], x_sing_t0, beta) : /**/ x1w1[0];
// 					w1 = isSingular ? x1w1[1] * DSg(x1w1[0], x_sing_t0, beta) : /**/ x1w1[1];
// 					// t1 = x1w1[0];
// 					// w1 = x1w1[1];
// 					// t1 = x1w1[0];
// 					J = intp->J(t0, t1);
// 					X = (*intp)(t0, t1);
// 					weight = w0 * w1;
// 					/* ------------------------------------------------------ */
// 					area += J * weight;
// 					/* ------------------------------------------------------ */
// 					r = X - in_point;
// 					grad_phi = -r / pow(Norm(r), 3);
// 					solution_intgration_in += Dot(grad_phi, Normalize(intp->cross(t0, t1))) * J * weight;
// 					/* ------------------------------------------------------ */
// 					r = X - out_point;
// 					grad_phi = -r / pow(Norm(r), 3);
// 					solution_intgration_out += Dot(grad_phi, Normalize(intp->cross(t0, t1))) * J * weight;
// 					/* ------------------------------------------------------ */
// 					// r = X - originX;
// 					// grad_phi = -r / pow(Norm(r), 3);
// 					// solution_intgration_surface += Dot(grad_phi, Normalize(intp->cross(t0, t1))) * J * weight;
// 					/* ------------------------------------------------------ */
// 					r = X - origin->getX();
// 					solution_intgration_surface_simplified += Dot(-r, intp->cross(t0, t1)) / pow(Norm(r), 3) * weight;
// 					/* ------------------------------------------------------ */
// 					r = X - oFace;
// 					solution_intgration_surface_simplified_face += Dot(-r, intp->cross(t0, t1)) / pow(Norm(r), 3) * weight;
// 					/* ------------------------------------------------------ */
// 					r = X - origin->getX();
// 					r_1_intgration_surface += 1. / Norm(r) * J * weight;
// 					/* ------------------------------------------------------ */
// 					r = X - oFace;
// 					r_1_intgration_surface_face += 1. / Norm(r) * J * weight;
// 				}
// 			}
// 		}

// 		// std::cout << check << std::endl;
// 		// std::cin.ignore();
// 		/* ------------------------------------------------------ */
// 		std::cout << "origin->getLines().size() = " << origin->getLines().size() << std::endl;
// 		std::cout << "gausspoints =" << gausspoints << ", checkSing =" << checkSing << ", beta =" << beta << std::endl;
// 		int k = gausspoints;
// 		std::cout << "面積 " << k << " = " << area
// 				  << ", 誤差=" << area / (4. * M_PI * 0.1 * 0.1) - 1 << std::endl;
// 		//
// 		auto tmp = solution_intgration_in;
// 		std::cout << "ガウス in " << k << " = " << tmp
// 				  << ", 誤差=" << tmp / (-4. * M_PI) - 1 << std::endl;
// 		//
// 		tmp = solution_intgration_out;
// 		std::cout << "ガウス out " << k << " = " << tmp << std::endl;
// 		//
// 		// tmp = solution_intgration_surface;
// 		// std::cout << Red << "ガウス surface " << k << " = " << tmp
// 		// 		  << ", 誤差=" << tmp / (-2. * M_PI) - 1. << reset << std::endl;
// 		//
// 		tmp = solution_intgration_surface_simplified;
// 		std::cout << Magenta << "ガウス surface ある節点を原点とした場合" << k << " = " << tmp
// 				  << ", 誤差=" << tmp / (-2. * M_PI) - 1. << reset << std::endl;

// 		tmp = solution_intgration_surface_simplified_face;
// 		std::cout << Magenta << "ガウス surface ある面上の点を原点とした場合" << k << " = " << tmp
// 				  << ", 誤差=" << tmp / (-2. * M_PI) - 1. << reset << std::endl;

// 		// std::cout << Green << "solution_intgration_surface_simplified_small = " << solution_intgration_surface_simplified_small << reset << std::endl;

// 		// auto Vtmp = solution_intgration_suface_from_eps_1;
// 		// std::cout << Green << "ガウス surface" << k << " = " << Vtmp
// 		//           << ", 誤差=" << Vtmp / (-2. * M_PI) - 1. << reset << std::endl;

// 		// Vtmp = minimazing_value;
// 		// std::cout << green << "minimazing_value" << k << " = " << minimazing_value << std::endl;
// 		//
// 		tmp = r_1_intgration_surface;
// 		std::cout << "ガウス 1/r　" << k << " = " << tmp << ",  / M_PI=" << tmp / M_PI << std::endl;
// 		tmp = r_1_intgration_surface_face;
// 		std::cout << "ガウス 1/r　" << k << " = " << tmp << ",  / M_PI=" << tmp / M_PI << std::endl;
// 	}
// };

// double checkSurfaceAreaQuadIDW(const V_netFp &fs, int gausspoints, bool checkSing = false)
// {
// 	//    \ \            ^
// 	//     7  0  6       |<-(t0=1.5)
// 	//   8  3   5  11    |<-(t0=1)
// 	//     1  4  2       |<-(t0=0)
// 	//      9  10  \     |
// 	//           \  \
	// //         (t1=1)(t1=0)
// 	double ret = 0.;
// 	V_netPp ps;
// 	for (const auto &f : fs)
// 	{
// 		auto origin = f->getPoints()[0];
// 		//!この積分範囲はしゅうせする必要がない
// 		interpolationCenterTriangleQuadIDW12 intp(extractX(f->get12Points(origin)) /*get12Points()で得られた点*/);
// 		// interpolationTriangleQuadByFixedRangeCenter intp(extractX(f->get6PointsTuple(origin)) /*get12Points()で得られた点*/);
// 		// interpolationLagLinear intp(extractX(f->getPoints(origin)) /*get12Points()で得られた点*/);
// 		if (checkSing)
// 		{
// 			//! IG
// 			double x0, x1, w0, w1;
// 			double x_sing_t0 = 0.;
// 			auto beta_t0 = 2.;
// 			auto gwgw_t0 = GaussianQuadratureWeights(gausspoints, InvSg(0., x_sing_t0, beta_t0), InvSg(1., x_sing_t0, beta_t0));
// 			for (const auto &x0w0 : gwgw_t0)
// 			{
// 				x0 = Sg(x0w0[0], x_sing_t0, beta_t0);
// 				w0 = x0w0[1] * DSg(x0w0[0], x_sing_t0, beta_t0);
// 				for (const auto &x1w1 : GWGWfrom0to1[gausspoints])
// 				{
// 					ret += intp.J(x0, x1w1[0]) * w0 * x1w1[1];
// 				}
// 			}
// 		}
// 		else
// 		{
// 			for (const auto &x0w0 : GWGWfrom0to1[gausspoints])
// 			{
// 				for (const auto &x1w1 : GWGWfrom0to1[gausspoints])
// 				{
// 					ret += intp.J(x0w0[0], x1w1[0]) * x0w0[1] * x1w1[1];
// 				}
// 			}
// 		}
// 	}
// 	return ret;
// };

// #define use_old_p_igign_Quad_IDW
#ifdef use_old_p_igign_Quad_IDW
void p_igign_Quad_IDW(const V_netPp &ps, const netPp origin, const int i, std::map<netP *, V_d> &new_igign_OUT) {
   //    \ \            ^
   //     7  0  6       |<-(t0=1.5)
   //   8  3   5  11    |<-(t0=1)
   //     1  4  2       |<-(t0=0)
   //      9  10  \     |
   //           \  \
	//         (t1=1)(t1=0)

   for (auto &[p, igign] : new_igign_OUT)
      igign = zeros;
   //! この積分範囲はしゅうせする必要がない
   interpolationCenterTriangleQuadIDW12 intp(extractX(ps) /*get12Points()で得られた点*/);
   V_d IGIGn(2, 0.), r(3, 0.), N(ps.size(), 0.), A = origin->getX();
   double common, nr;
   // 0から1/2に変えよう
   //  ここの積分の整理から始めよう
   if (ps[4] == origin) {
      //! IG
      double x_sing_t0 = 0.;
      auto beta_t0 = 2.;
      auto gwgw_t0 = GaussianQuadratureWeights(i, InvSg(0., x_sing_t0, beta_t0), InvSg(1., x_sing_t0, beta_t0));
      double x0, x1, w0, w1;
      for (const auto &x0w0 : gwgw_t0) {
         x0 = Sg(x0w0[0], x_sing_t0, beta_t0);
         w0 = x0w0[1] * DSg(x0w0[0], x_sing_t0, beta_t0);
         for (const auto &x1w1 : GWGWfrom0to1[i]) {
            x1 = x1w1[0];
            w1 = x1w1[1];
            {
               common = intp.J(x0, x1) * w0 * w1;
               r = intp(x0, x1) - A;
               nr = Norm(r);
               IGIGn[0] = common / nr;
               IGIGn[1] = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0, x1))) * common;
               N = intp.N(x0, x1);
               for (auto k = 0; k < ps.size(); k++)
                  new_igign_OUT[ps[k]] += IGIGn * N[k];
            }
         }
      }
   } else {
      for (const auto &x0w0 : GWGWfrom0to1[i]) {
         for (const auto &x1w1 : GWGWfrom0to1[i]) {
            common = intp.J(x0w0[0], x1w1[0]) * x0w0[1] * x1w1[1];
            r = intp(x0w0[0], x1w1[0]) - A;
            nr = Norm(r);
            IGIGn[0] = common / nr;
            IGIGn[1] = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0w0[0], x1w1[0]))) * common;
            N = intp.N(x0w0[0], x1w1[0]);
            for (auto k = 0; k < ps.size(); k++)
               new_igign_OUT[ps[k]] += IGIGn * N[k];
         }
      }
   }
};
#else
// void p_igign_Quad_IDW(const V_netPp &ps, const netPp origin, const int i, std::map<netP *, V_d> &new_igign_OUT)
// {
// 	//    \ \            ^
// 	//     7  0  6       |<-(t0=1.5)
// 	//   8  3   5  11    |<-(t0=1)
// 	//     1  4  2       |<-(t0=0)
// 	//      9  10  \     |
// 	//           \  \
	// 	//         (t1=1)(t1=0)
// 	for (auto &[p, igign] : new_igign_OUT)
// 	{
// 		igign[0] = 0.;
// 		igign[1] = 0.;
// 	}
// 	interpolationCenterTriangleQuadIDW12 intp(extractX(ps) /*get12Points()で得られた点*/);
// 	V_d IGIGn(2, 0.), r(3, 0.), N(ps.size(), 0.), A = origin->getX();
// 	double common, nr, weight;
// 	V_d cross;
// 	if (ps[4] == origin)
// 	{
// 		auto beta_IG = 2.;
// 		double x_sing_t0 = 0.;
// 		auto gwgw_IG = GaussianQuadratureWeights(i, InvSg(0., x_sing_t0, beta_IG), InvSg(1., x_sing_t0, beta_IG));
// 		double x0, x1, w0, w1;
// 		for (const auto &x0w0 : gwgw_IG)
// 		{
// 			x0 = Sg(x0w0[0], x_sing_t0, beta_IG);
// 			w0 = x0w0[1] * DSg(x0w0[0], x_sing_t0, beta_IG);
// 			for (const auto &x1w1 : GWGWfrom0to1[i])
// 			{
// 				x1 = x1w1[0];
// 				w1 = x1w1[1];
// 				{
// 					cross = intp.cross(x0, x1);
// 					r = intp(x0, x1) - A;
// 					nr = Norm(r);
// 					IGIGn[0] = 1 / nr * Norm(cross) * w0 * w1;
// 					IGIGn[1] = -Dot(r / pow(nr, 3.), cross) * w0 * w1;
// 					N = intp.N(x0, x1);
// 					for (auto k = 0; k < ps.size(); k++)
// 						new_igign_OUT[ps[k]] += IGIGn * N[k];
// 				}
// 			}
// 		}
// 	}
// 	else
// 	{
// 		for (const auto &x0w0 : GWGWfrom0to1[i])
// 		{
// 			for (const auto &x1w1 : GWGWfrom0to1[i])
// 			{
// 				cross = intp.cross(x0w0[0], x1w1[0]);
// 				r = intp(x0w0[0], x1w1[0]) - A;
// 				nr = Norm(r);
// 				IGIGn[0] = Norm(cross) / nr * x0w0[1] * x1w1[1];
// 				IGIGn[1] = -Dot(r / pow(nr, 3.), cross) * x0w0[1] * x1w1[1];
// 				N = intp.N(x0w0[0], x1w1[0]);
// 				for (auto k = 0; k < ps.size(); k++)
// 					new_igign_OUT[ps[k]] += IGIGn * N[k];
// 			}
// 		}
// 	}
// };

/* -------------------------- 面 ------------------------- */

// 	void p_igign(const netFp &f /*積分する面を構成する点*/,
// 				 const netFp origin_f,
// 				 double t0, double t1, const int i,
// 				 std::map<netP *, V_d> &new_igign_OUT, const V_netPp &ps)
// 	{
// 		for (auto &[p, igign] : new_igign_OUT)
// 			igign = zeros;

// 		// interpolationTriangle *intp;
// 		// if (ps.size() == 3)
// 		// 	intp = new interpolationTriangleLinearByFixedRange();
// 		// else if (ps.size() == 6)
// 		// 	intp = new interpolationTriangleQuadByFixedRangeCenter();
// 		// else if (ps.size() == 12)
// 		// 	intp = new interpolationCenterTriangleQuadIDW12();
// 		// else
// 		// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

// 		auto intp = new interpolationTriangleLinearByFixedRange();

// 		intp->set(extractX(origin_f->getPoints()));
// 		V_d A = (*intp)(t0, t1);
// 		intp->set(extractX(ps));
// 		V_d IGIGn(2, 0.), r(3, 0.), N(ps.size(), 0.);
// 		double nr, weight;
// 		V_d cross;
// // #define quadrature_for_singular
// #ifdef quadrature_for_singular

// 		if (f == origin_f)
// 		{
// 			auto b = 2.;
// 			double t0_sing = 0.5;
// 			auto gwgw_IG = GaussianQuadratureWeights(i, InvSg(0., t0_sing, b), InvSg(1., t0_sing, b));
// 			double x0, x1, w0, w1;
// 			for (const auto &x0w0 : gwgw_IG)
// 			{
// 				x0 = Sg(x0w0[0], t0_sing, b);
// 				w0 = x0w0[1] * DSg(x0w0[0], t0_sing, b);
// 				for (const auto &x1w1 : gwgw_IG)
// 				{
// 					x1 = Sg(x1w1[0], t0_sing, b);
// 					w1 = x1w1[1] * DSg(x1w1[0], t0_sing, b);
// 					{
// 						cross = intp->cross(x0, x1);
// 						r = (*intp)(x0, x1) - A;
// 						nr = Norm(r);
// 						IGIGn[0] = Norm(cross) / nr * w0 * w1;
// 						IGIGn[1] = -Dot(r / pow(nr, 3.), cross) * w0 * w1;
// 						N = intp->N(x0, x1);
// 						for (auto k = 0; k < ps.size(); k++)
// 							new_igign_OUT[ps[k]] += IGIGn * N[k];
// 					}
// 				}
// 			}
// 		}
// 		else
// #endif
// 		{
// 			for (const auto &x0w0 : GWGWfrom0to1[i])
// 				for (const auto &x1w1 : GWGWfrom0to1[i])
// 				{
// 					cross = intp->cross(x0w0[0], x1w1[0]);
// 					r = (*intp)(x0w0[0], x1w1[0]) - A;
// 					nr = Norm(r);
// 					IGIGn[0] = Norm(cross) / nr * x0w0[1] * x1w1[1];
// 					IGIGn[1] = -Dot(r / pow(nr, 3.), cross) * x0w0[1] * x1w1[1];
// 					N = intp->N(x0w0[0], x1w1[0]);
// 					for (auto k = 0; k < ps.size(); k++)
// 						new_igign_OUT[ps[k]] += IGIGn * N[k];
// 				}
// 		}
// 		delete intp;
// 	};

#endif
/* ------------------ 三角形面の2次のラグランジュ補間 ------------------ */
void p_igign_Quad(const V_netPp &ps, const netPp origin, const int i, std::map<netP *, V_d> &new_igign_OUT) {
   //! shape function vector of gaussian points
   //! __GWGW__ = {xi, h * (1. - xi), w_xi * w_h * (1. - xi), w_xi * w_h /*1-xiなし*/}}
   //! __GWGW__ = {xi, h, w_xi * w_h}
   // /Users/tomoaki/Dropbox/markdown/cpp/builds/build_quadrature/main.cpp
   // を見た下さい
   //! あるoriginに対して，各点の影響力igignを計算する．
   //! 影響力係数は全てのFaceそれぞれについて計算するが(origin A v.s. Faces)，
   //! FaceはPoint座標の関数なので，結果はPoints各点の影響力として表され，係数行列として出てくる．
   //! Face達を構成するPoint達は，2次以上の場合，重複していることがある．重複している場合は，影響力を加算していくことで考慮される．
   for (auto &[p, igign] : new_igign_OUT)
      igign = zeros;
   interpolationTriangleQuadByFixedRangeCenter intp(extractX(ps));
   V_d IGIGn(2, 0.), r(3, 0.), N(ps.size(), 0.), A = ToVector(origin->getXtuple());
   double common, nr;
   // 0から1/2に変えよう
   //  ここの積分の整理から始めよう
   if (ps[4] == origin) {
      //! IG
      double x_sing_t0 = 0.;
      auto beta_t0 = 2.;
      auto gwgw_t0 = GaussianQuadratureWeights(i, InvSg(0., x_sing_t0, beta_t0), InvSg(1., x_sing_t0, beta_t0));
      double x0, x1, w0, w1;
      for (const auto &x0w0 : gwgw_t0) {
         x0 = Sg(x0w0[0], x_sing_t0, beta_t0);
         w0 = x0w0[1] * DSg(x0w0[0], x_sing_t0, beta_t0);
         for (const auto &x1w1 : GWGWfrom0to1[i]) {
            x1 = x1w1[0];
            w1 = x1w1[1];
            {
               common = intp.J(x0, x1) * w0 * w1;
               r = intp(x0, x1) - A;
               nr = Norm(r);
               IGIGn[0] = common / nr;
               IGIGn[1] = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0, x1))) * common;
               N = intp.N(x0, x1);
               for (auto k = 0; k < ps.size(); k++)
                  new_igign_OUT[ps[k]] += IGIGn * N[k];
            }
         }
      }
   } else {
      for (const auto &x0w0 : GWGWfrom0to1[i]) {
         for (const auto &x1w1 : GWGWfrom0to1[i]) {
            common = intp.J(x0w0[0], x1w1[0]) * x0w0[1] * x1w1[1];
            r = intp(x0w0[0], x1w1[0]) - A;
            nr = Norm(r);
            IGIGn[0] = common / nr;
            IGIGn[1] = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0w0[0], x1w1[0]))) * common;
            N = intp.N(x0w0[0], x1w1[0]);
            for (auto k = 0; k < ps.size(); k++)
               new_igign_OUT[ps[k]] += IGIGn * N[k];
         }
      }
   }
};

/* ---------------------- 三角形面の線形補間 --------------------- */

void p_igign_Linear(const V_netPp &ps, const netPp origin, const int i, std::map<netP *, V_d> &new_igign_OUT) {
   // for (auto i = 0; i < 3; i++)
   // 	new_igign_OUT[ps[i]] = zeros;
   // V_d r(3, 0.), N(3, 0.), IGIGn(2, 0.), A = origin->getX();
   // double x0, x1, nr, tmp;
   // interpolationLagLinear intp(extractX(ps));
   // for (const auto &a_b_1ab_w : __GWGW__[i]) {
   // 	x0 = a_b_1ab_w[0];
   // 	x1 = a_b_1ab_w[1];
   // 	x1 *= (1 - x0);
   // 	tmp = intp.J(x0, x1) * (1 - x0) * a_b_1ab_w[2] /*変数変換した重みは[2]*/;
   // 	//
   // 	r = intp(x0, x1) - A;
   // 	nr = Norm(r);
   // 	IGIGn[0] = tmp / nr;
   // 	IGIGn[1] = -Dot(r / pow(nr, 3.), intp.normal_cross(x0, x1)) * tmp;
   // 	N = intp.N(x0, x1);
   // 	for (auto i = 0; i < 3; i++)
   // 		new_igign_OUT[ps[i]] += IGIGn * N[i];
   // }
   for (auto &[p, igign] : new_igign_OUT)
      igign = zeros;
   interpolationTriangleLinearByFixedRange intp(extractX(ps));
   V_d IGIGn(2, 0.), r(3, 0.), N(ps.size(), 0.), A = ToVector(origin->getXtuple());
   auto gw = GaussianQuadratureWeights(i, 0., 1.);
   double common, nr;
   double x0, x1, w0, w1;
   bool isSingular = false;  //(ps[0] == origin);  //!linearの場合は0
   double beta = 2.;
   double x_sing = 1.;  //! linearの場合はx_sing = 1.
   for (const auto &x0w0 : isSingular ? GaussianQuadratureWeights(15, InvSg(0., x_sing, beta), InvSg(1., x_sing, beta)) : gw) {
      if (isSingular) {
         x0 = Sg(x0w0[0], x_sing, beta);
         w0 = x0w0[1] * DSg(x0w0[0], x_sing, beta);
      } else {
         x0 = x0w0[0];
         w0 = x0w0[1];
      }
      for (const auto &x1w1 : gw) {
         x1 = x1w1[0];
         w1 = x1w1[1];
         ///////
         common = intp.J(x0, x1) * w0 * w1;
         r = intp(x0, x1) - A;
         nr = Norm(r);
         IGIGn[0] = common / nr;
         IGIGn[1] = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0, x1))) * common;
         N = intp.N(x0, x1);
         for (auto k = 0; k < ps.size(); k++)
            new_igign_OUT[ps[k]] += IGIGn * N[k];
      }
   }
};
//
// void p_igign_Linear(networkPoint *const p0,
// 					networkPoint *const p1,
// 					networkPoint *const p2,
// 					const netPp origin,
// 					const int i,
// 					std::unordered_map<netP *, Tdd> &new_igign_OUT)
// {
// 	// Tdd zeros = {0, 0};
// 	// for (auto &[p, igign] : new_igign_OUT)
// 	// 	igign = zeros;

// 	// タプルに対応させる．
// 	interpolationTriangleLinearByFixedRange3D intp(T3Tddd{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
// 	Tdd IGIGn;
// 	Tddd r, N, A = origin->getXtuple();
// 	//
// 	// auto gw = GaussianQuadratureWeights(i, 0., 1.);

// 	// if (p0 == origin || p1 == origin || p2 == origin)
// 	// 	gw = GaussianQuadratureWeights(15, 0., 1.);
// 	auto gw = __GWGW__[i + 1];
// 	double common, nr;
// 	double x0, x1, w0, w1;
// 	// bool isSingular = false; //(ps[0] == origin);  //!linearの場合は0
// 	bool isSingular = (p0 == origin); //!linearの場合は0
// 	double beta = 2.;
// 	double x_sing = 1.; //!linearの場合はx_sing = 1.
// 	auto &new_igign_OUTp0 = new_igign_OUT[p0];
// 	auto &new_igign_OUTp1 = new_igign_OUT[p1];
// 	auto &new_igign_OUTp2 = new_igign_OUT[p2];
// 	for (const auto &x0w0 : isSingular ? GaussianQuadratureWeights(15, InvSg(0., x_sing, beta), InvSg(1., x_sing, beta)) : gw)
// 	{
// 		if (isSingular)
// 		{
// 			x0 = Sg(x0w0[0], x_sing, beta);
// 			w0 = x0w0[1] * DSg(x0w0[0], x_sing, beta);
// 		}
// 		else
// 		{
// 			x0 = x0w0[0];
// 			w0 = x0w0[1];
// 		}
// 		for (const auto &x1w1 : gw)
// 		{
// 			x1 = x1w1[0];
// 			w1 = x1w1[1];
// 			///////
// 			common = intp.J(x0, x1) * w0 * w1;
// 			r = intp(x0, x1) - A;
// 			nr = Norm(r);
// 			std::get<0>(IGIGn) = common / nr;
// 			std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0, x1))) * common;
// 			N = intp.N(x0, x1);
// 			new_igign_OUTp0 += IGIGn * std::get<0>(N);
// 			new_igign_OUTp1 += IGIGn * std::get<1>(N);
// 			new_igign_OUTp2 += IGIGn * std::get<2>(N);
// 		}
// 	}
// };
// 昔の
void p_igign_Linear(networkPoint *p0,
                    networkPoint *p1,
                    networkPoint *p2,
                    const netPp origin,
                    const int i,
                    std::unordered_map<netP *, Tdd> &new_igign_OUT) {
   // Tdd zeros = {0, 0};
   // for (auto &[p, igign] : new_igign_OUT)
   // 	igign = zeros;

   // タプルに対応させる．
   interpolationTriangleLinearByFixedRange3D intp(T3Tddd{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
   Tdd IGIGn;
   Tddd r, N, A = origin->getXtuple();
   //
   auto gw = GaussianQuadratureWeights(i, 0., 1.);

   // if (p0 == origin || p1 == origin || p2 == origin)
   // 	gw = GaussianQuadratureWeights(15, 0., 1.);
   double common, nr;
   double x0, x1, w0, w1;
   // bool isSingular = false; //(ps[0] == origin);  //!linearの場合は0
   bool isSingular = (p0 == origin);  //! linearの場合は0
   double beta = 2.;
   double x_sing = 1.;  //! linearの場合はx_sing = 1.
   auto &new_igign_OUTp0 = new_igign_OUT[p0];
   auto &new_igign_OUTp1 = new_igign_OUT[p1];
   auto &new_igign_OUTp2 = new_igign_OUT[p2];
   for (const auto &x0w0 : isSingular ? GaussianQuadratureWeights(15, InvSg(0., x_sing, beta), InvSg(1., x_sing, beta)) : gw) {
      if (isSingular) {
         x0 = Sg(x0w0[0], x_sing, beta);
         w0 = x0w0[1] * DSg(x0w0[0], x_sing, beta);
      } else {
         x0 = x0w0[0];
         w0 = x0w0[1];
      }
      for (const auto &x1w1 : gw) {
         x1 = x1w1[0];
         w1 = x1w1[1];
         ///////
         common = intp.J(x0, x1) * w0 * w1;
         r = intp(x0, x1) - A;
         nr = Norm(r);
         std::get<0>(IGIGn) = common / nr;
         std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0, x1))) * common;
         N = intp.N(x0, x1);
         new_igign_OUTp0 += IGIGn * std::get<0>(N);
         new_igign_OUTp1 += IGIGn * std::get<1>(N);
         new_igign_OUTp2 += IGIGn * std::get<2>(N);
      }
   }
};
/* ------------------------------------------------------ */
//* ------------------------------------------------------ */
// #          RBF  面における方程式を点上で作成                  */
// //* ------------------------------------------------------ */
// void p_igign_Quad_IDW_RBF(const V_netPp &ps, const netPp origin, const int i, std::map<netP *, V_d> &new_igign_OUT)
// {
// 	//    \ \            ^
// 	//     7  0  6       |<-(t0=1.5)
// 	//   8  3   5  11    |<-(t0=1)
// 	//     1  4  2       |<-(t0=0)
// 	//      9  10  \     |
// 	//           \  \
	// 	//         (t1=1)(t1=0)
// 	for (auto &[p, igign] : new_igign_OUT)
// 		igign = zeros;
// 	//!この積分範囲はしゅうせする必要がない
// 	interpolationCenterTriangleQuadIDW12 intp(extractX(ps) /*get12Points()で得られた点*/);
// 	V_d IGIGn(2, 0.), r(3, 0.), N(ps.size(), 0.), A = ToVector(origin->getXtuple());
// 	double common, nr, weight;
// 	V_d cross;
// 	if (ps[4] == origin)
// 	{
// 		double s = (double)(origin->getFaces()).size();
// 		// こういう，一度だけ！見たいな場合の処理が面倒．
// 		//! ------------------------------------------------------ */
// 		//!                 RBFを使って補間関数を作成                 */
// 		//! ------------------------------------------------------ */
// 		VV_d xyz = {};
// 		VV_d param = {};
// 		V_netPp points = {};
// 		// for (const auto &tup : origin->getNeighbors_Depth2_OnPolarAsTuple2())
// 		for (const auto &tup : origin->getNeighbors_Depth2_OnPolarAsTuple())
// 		{
// 			auto [t0, t1, X, p] = tup;
// 			points.emplace_back(p);
// 			param.emplace_back(V_d{t0, t1});
// 			xyz.emplace_back(X);
// 		}
// 		points.emplace_back(origin);
// 		param.emplace_back(V_d{0, 0});
// 		xyz.emplace_back(ToVector(origin->getXtuple()));
// 		auto intp = InterpolationVectorRBF(param, xyz);
// 		//! ------------------------------------------------------ */
// 		int gausspoints = 6;
// 		auto gwgw_polar = GaussianQuadratureWeights(gausspoints, 0., M_PI);
// 		auto beta = 1.;
// 		auto x_sing_t0 = 0.;
// 		auto gwgw_r = GaussianQuadratureWeights(gausspoints, InvSg(-1., x_sing_t0, beta), InvSg(1., x_sing_t0, beta));
// 		double x0, x1, w0, w1;
// 		for (const auto &x0w0 : gwgw_r)
// 		{
// 			x0 = Sg(x0w0[0], x_sing_t0, beta); // r
// 			w0 = x0w0[1] * DSg(x0w0[0], x_sing_t0, beta);
// 			for (const auto &x1w1 : gwgw_polar)
// 			{
// 				x1 = x1w1[0];
// 				w1 = x1w1[1];

// 				auto X0 = x0 * cos(x1);
// 				auto X1 = x0 * sin(x1);
// 				{
// 					r = intp({X0, X1}) - A;
// 					nr = Norm(r);
// 					IGIGn[0] = 1. / nr * intp.J(X0, X1) * x0 * (w0 * w1) / s;
// 					IGIGn[1] = -Dot(r / pow(nr, 3.), intp.cross(X0, X1)) * x0 * (w0 * w1) / s;
// 					// このNは未だ決まらない．愛想パラメトリックだからこそ決まっていたもの．
// 					// RBFの場合，サンプルをまず与えて，Fの逆行列をもとめてから初めてわかる．．．．
// 					//!仮にアイソパラメトリックでやってみよう
// 					N = intp.N(X0, X1);
// 					for (auto k = 0; k < points.size(); k++)
// 					{
// 						if (new_igign_OUT.find(points[k]) == new_igign_OUT.end())
// 							new_igign_OUT[points[k]] = {0., 0.};

// 						new_igign_OUT[points[k]] += IGIGn * N[k]; //面の数で割った
// 					}
// 				}
// 			}
// 		}
// 	}
// 	else
// 	{
// 		for (const auto &x0w0 : GWGWfrom0to1[i])
// 		{
// 			for (const auto &x1w1 : GWGWfrom0to1[i])
// 			{
// 				cross = intp.cross(x0w0[0], x1w1[0]);
// 				r = intp(x0w0[0], x1w1[0]) - A;
// 				nr = Norm(r);
// 				IGIGn[0] = Norm(cross) / nr * x0w0[1] * x1w1[1];
// 				IGIGn[1] = -Dot(r / pow(nr, 3.), cross) * x0w0[1] * x1w1[1];
// 				N = intp.N(x0w0[0], x1w1[0]);
// 				for (auto k = 0; k < ps.size(); k++)
// 					new_igign_OUT[ps[k]] += IGIGn * N[k];
// 			}
// 		}
// 	}
// };
// //* ------------------------------------------------------ */
// #                 RBF 各点んで方程式を作る場合              */
//* ------------------------------------------------------ */
// std::map<netP *, V_d> calc_P_IGIGn_RBF(const netFp &f, const netPp origin)
// {
// 	/*
// 	 * 点を原点とした場合の境界積分方程式の一部：IGIGnを計算する
// 	 */
// 	V_netPp ps = f->get12Points(origin); //基準としてoriginを使おうとする．が，面に含まれない点の場合は，適当な面の点が内部で自動で選ばれる．
// 	V_i nGWGW = {6 + 1, 8 + 1, 10 + 1, 12 + 1};
// 	double conv = 1E-2;
// 	V_d zeros(2, 0.);
// 	std::map<netPp, V_d> pre_igign;
// 	for (const auto &p : ps)
// 		pre_igign[p] = zeros;
// 	std::map<netPp, V_d> new_igign = pre_igign;
// 	p_igign_Quad_IDW_RBF(ps, origin, nGWGW[0], /*出力結果*/ pre_igign); //!まず計算
// 	for (auto i = 1; i < nGWGW.size(); i++)
// 	{
// 		p_igign_Quad_IDW_RBF(ps, origin, nGWGW[i], /*出力結果*/ new_igign); //!次に計算
// 		if (isConverged(new_igign, pre_igign, ps, conv))
// 			return new_igign;
// 		else
// 			pre_igign = new_igign; // replace and try again
// 	}
// 	return new_igign;
// };

//* ------------------------------------------------------ */
//%               面における方程式を点上で作成                  */
//* ------------------------------------------------------ */
// std::map<netP *, V_d> calc_P_IGIGn(const netFp &f, const netPp origin)
// {
// 	/* code */
// 	/**
// 	 * 点を原点とした場合の境界積分方程式の一部：IGIGnを計算する
// 	 */
// 	/* --------------------- 線形か2次かを決める --------------------- */
// 	/* ---------------- 引数に与えられた面を，点を使って表現し直す --------------- */
// 	V_netPp ps = f->getPoints(origin);
// 	// V_netPp ps = f->get6PointsTuple(origin);
// 	// V_netPp ps = f->get12Points(origin); //基準としてoriginを使おうとする．が，面に含まれない点の場合は，適当な面の点が内部で自動で選ばれる．
// 	/* ------------------------ 順番を調整 ----------------------- */
// 	// if (ps.size() == 3) {
// 	// 	if (MemberQ(ps, origin)) {
// 	// 		if (ps[1] == origin)
// 	// 			ps = {ps[1] /*Origin*/, ps[2], ps[0]};
// 	// 		if (ps[2] == origin)
// 	// 			ps = {ps[2] /*Origin*/, ps[0], ps[1]};
// 	// 	}
// 	// } else if (ps.size() == 6) {
// 	// 	if (MemberQ(ps, origin)) {
// 	// 		//! {ps[0], ps[1], ps[2], ps[3], ps[4]/*O*//*ここがxi=0の座標のはず*/, ps[5]}
// 	// 		if (ps[5] == origin)
// 	// 			ps = {ps[1], ps[2], ps[0], ps[4], ps[5] /*O ここがxi=0の座標のはず*/, ps[3]};
// 	// 		else if (ps[3] == origin)
// 	// 			ps = {ps[2], ps[0], ps[1], ps[5], ps[3] /*O ここがxi=0の座標のはず*/, ps[4]};
// 	// 	}
// 	// }
// 	/* ------------------------------------------------------ */
// 	// V_i nGWGW = {6, 8, 10, 12};
// 	V_i nGWGW = {5};
// 	double conv = 1E-2;
// 	V_d zeros(2, 0.);
// 	std::map<netPp, V_d> pre_igign;
// 	for (const auto &p : ps)
// 		pre_igign[p] = zeros;
// 	std::map<netPp, V_d> new_igign = pre_igign;
// 	if (ps.size() == 3)
// 	{
// 		p_igign_Linear(ps, origin, nGWGW[0], /*出力結果*/ pre_igign); //!まず計算
// 		new_igign = pre_igign;
// 		for (auto i = 1; i < nGWGW.size(); i++)
// 		{
// 			p_igign_Linear(ps, origin, nGWGW[i], /*出力結果*/ new_igign); //!次に計算
// 			if (isConverged(new_igign, pre_igign, ps, conv))
// 				return new_igign;
// 			else
// 				pre_igign = new_igign; // replace and try again
// 		}
// 	}
// 	else if (ps.size() == 6)
// 	{
// 		p_igign_Quad(ps, origin, nGWGW[0], /*出力結果*/ pre_igign); //!まず計算
// 		new_igign = pre_igign;
// 		for (auto i = 1; i < nGWGW.size(); i++)
// 		{
// 			p_igign_Quad(ps, origin, nGWGW[i], /*出力結果*/ new_igign); //!次に計算
// 			if (isConverged(new_igign, pre_igign, ps, conv))
// 				return new_igign;
// 			else
// 				pre_igign = new_igign; // replace and try again
// 		}
// 	}
// 	else
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	return new_igign;
// };
/* ------------------------------------------------------ */
// std::unordered_map<netP *, Tdd> calc_P_IGIGnTuple(const netFp &f, const netPp origin)
// {
// 	V_netPp ps = f->getPoints(origin);
// 	Tdd zeros = {0, 0};
// 	std::unordered_map<netPp, Tdd> igign = {{ps[0], zeros}, {ps[1], zeros}, {ps[2], zeros}};
// 	p_igign_Linear(ps[0], ps[1], ps[2], origin, 5., /*出力結果*/ igign); //!まず計算
// 	return igign;
// };
/* ------------------------------------------------------ */
std::vector<std::tuple<netP *, Tdd>> calc_P_IGIGnLinearTuple(const networkFace *const f, const networkPoint *const origin) {
   // おそらくIgIgnの計算精度に大きく関わるため，
   const auto [p0, p1, p2] = f->getPoints(origin);
   interpolationTriangleLinearByFixedRange3D intp(T3Tddd{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
   Tdd IGIGn, p0igign = {0., 0.}, p1igign = {0., 0.}, p2igign = {0., 0.};
   Tddd r, N, A = origin->getXtuple(), cross;
   double nr;
   // for (const auto &[x0, x1, w0w1] : (p0 == origin ? __GWGW_beta2_8__ : __GWGW8__Tuple))
   // {
   // 	std::get<0>(IGIGn) = (Norm(cross = intp.cross(x0, x1)) * w0w1) / (nr = Norm(r = intp(x0, x1) - A));
   // 	std::get<1>(IGIGn) = -Dot(r / (nr * nr * nr), cross) * w0w1;
   // 	p0igign += IGIGn * std::get<0>(N = intp.N(x0, x1));
   // 	p1igign += IGIGn * std::get<1>(N);
   // 	p2igign += IGIGn * std::get<2>(N);
   // }
   for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple) {
      std::get<0>(IGIGn) = (Norm(cross = intp.cross(x0, x1)) * w0w1) / (nr = Norm(r = intp(x0, x1) - A));
      std::get<1>(IGIGn) = -Dot(r / (nr * nr * nr), cross) * w0w1;
      p0igign += IGIGn * std::get<0>(N = intp.N(x0, x1));
      p1igign += IGIGn * std::get<1>(N);
      p2igign += IGIGn * std::get<2>(N);
   }
   return {{p0, p0igign}, {p1, p1igign}, {p2, p2igign}};
};

std::tuple<std::tuple<netP *, Tdd>, std::tuple<netP *, Tdd>, std::tuple<netP *, Tdd>>
calc_P_IGIGnLinear3Tuples(const networkFace *const f, const networkPoint *const origin) {
   // おそらくIgIgnの計算精度に大きく関わるため，
   const auto [p0, p1, p2] = f->getPoints(origin);
   Tddd X0 = p0->getXtuple();
   Tddd X1 = p1->getXtuple();
   Tddd X2 = p2->getXtuple();
   interpolationTriangleLinearByFixedRange3D intp(T3Tddd{X0, X1, X2});
   Tdd IGIGn, p0igign = {0., 0.}, p1igign = {0., 0.}, p2igign = {0., 0.};
   Tddd r, N, A = origin->getXtuple(), cross;
   double nr;
   // for (const auto &[x0, x1, w0w1] : (p0 == origin ? __GWGW_beta2_10__ : __GWGW10__Tuple))
   // {
   // 	std::get<0>(IGIGn) = (Norm(cross = intp.cross(x0, x1)) * w0w1) / (nr = Norm(r = (intp(x0, x1) - A)));
   // 	std::get<1>(IGIGn) = -Dot(r, cross) * w0w1 / (nr * nr * nr);
   // 	p0igign += IGIGn * std::get<0>(N = intp.N(x0, x1));
   // 	p1igign += IGIGn * std::get<1>(N);
   // 	p2igign += IGIGn * std::get<2>(N);
   // }
   // for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
   // {
   // 	nr = Norm(r = (intp(x0, x1) - A));
   // 	r = (intp(x0, x1) - A);
   // 	std::get<0>(IGIGn) = (Norm(cross = intp.cross(x0, x1)) * w0w1) / Norm(intp(x0, x1) - A);
   // 	std::get<1>(IGIGn) = -Dot(r, cross) * w0w1 / (nr * nr * nr);
   // 	p0igign += IGIGn * x0;
   // 	p1igign += IGIGn * x1 * (1. - x0);
   // 	p2igign += IGIGn * (-1. + x0) * (-1. + x1);
   // }
   // return {{p0, p0igign}, {p1, p1igign}, {p2, p2igign}};

   // X0 -= X2;
   // X1 -= X2;
   // A -= X2;
   // X2 -= X2;
   // for (const auto &[t0, t1, w0w1] : __GWGW10__Tuple)
   // {
   // 	r = (X0 * t0 + X1 * t1 * (1 - t0) - A);
   // 	nr = Norm(r);
   // 	cross = Cross(X0 - X1 * t1, X1 * (1. - t0));
   // 	std::get<0>(IGIGn) = w0w1 * Norm(cross) / nr;
   // 	std::get<1>(IGIGn) = -w0w1 * Dot(r, cross) / (nr * nr * nr);
   // 	p0igign += IGIGn * t0;
   // 	p1igign += IGIGn * t1 * (1. - t0);
   // 	p2igign += IGIGn * (-1. + t0) * (-1. + t1);
   // }
   // return {{p0, p0igign}, {p1, p1igign}, {p2, p2igign}};

   X0 -= X2;
   X1 -= X2;
   A -= X2;
   X2 -= X2;
   cross = Cross(X0, X1);
   double dot = Dot(-A, cross);
   double nr_cross = Norm(cross);
   for (const auto &[t0, t1, w0w1] : __GWGW10__Tuple) {
      // Dot(r, Cross(X0, X1))
      // Dot((X0 * t0 + X1 * t1 * (1 - t0) - A), Cross(X0, X1))
      // =Dot(X1, Cross((X1 * t1 * (1 - t0) - A),X0))
      // =Dot(X0, Cross(X1,(X1 * t1 * (1 - t0) - A)))
      // =Dot(X0, Cross(X1,- A))
      // =Dot(- A, Cross(X0,X1))
      nr = Norm(X0 * t0 + X1 * t1 * (1 - t0) - A);
      std::get<0>(IGIGn) = w0w1 * nr_cross * (1. - t0) / nr;
      std::get<1>(IGIGn) = -w0w1 * dot * (1. - t0) / (nr * nr * nr);
      p0igign += IGIGn * t0;
      p1igign += IGIGn * t1 * (1. - t0);
      p2igign += IGIGn * (-1. + t0) * (-1. + t1);
   }
   return {{p0, p0igign}, {p1, p1igign}, {p2, p2igign}};
};
/* ------------------------------------------------------ */
std::tuple<std::tuple<netP *, Tdd>, std::tuple<netP *, Tdd>, std::tuple<netP *, Tdd>>
calc_P_IGIGnLinear3Tuples2(const networkFace *const f, const networkPoint *const origin) {
   // おそらくIgIgnの計算精度に大きく関わるため，
   const auto [p0, p1, p2] = f->getPoints(origin);
   std::tuple<std::tuple<netP *, Tdd>, std::tuple<netP *, Tdd>, std::tuple<netP *, Tdd>> ret({{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}});
   Tdd IGIGn;  //, p0igign = {0., 0.}, p1igign = {0., 0.}, p2igign = {0., 0.};
   double nr, tmp;
   Tddd X2 = p2->getXtuple();
   Tddd X0 = p0->getXtuple() - X2;
   Tddd X1 = p1->getXtuple() - X2;
   Tddd A = origin->getXtuple() - X2;
   Tddd cross = Cross(X0, X1);
   Tdd c = {Norm(cross), Dot(A, cross)};
   if (origin == p0 || origin == p1 || origin == p2)
      std::get<1>(c) = 0;
   /* ------------------------------------ */
   // for (const auto &[t0, t1, w0w1] : __GWGW5__Tuple) {
   //    // Dot(r, Cross(X0, X1))
   //    // Dot((X0 * t0 + X1 * t1 * (1 - t0) - A), Cross(X0, X1))
   //    // =Dot(X1, Cross((X1 * t1 * (1 - t0) - A),X0))
   //    // =Dot(X0, Cross(X1,(X1 * t1 * (1 - t0) - A)))
   //    // =Dot(X0, Cross(X1,- A))
   //    // =Dot(- A, Cross(X0,X1))
   //    nr = Norm(X0 * t0 + X1 * t1 * (1. - t0) - A);
   //    tmp = w0w1 * (1. - t0) / nr;
   //    // std::get<0>(IGIGn) = tmp;
   //    // std::get<1>(IGIGn) = tmp / (nr * nr);
   //    IGIGn = {tmp, tmp / (nr * nr)};
   //    // p0igign += IGIGn * t0;
   //    // p1igign += IGIGn * t1 * (1. - t0);
   //    // p2igign += IGIGn * (-1. + t0) * (-1. + t1);
   //    std::get<1>(std::get<0>(ret)) += IGIGn * t0;
   //    std::get<1>(std::get<1>(ret)) += IGIGn * t1 * (1. - t0);
   //    std::get<1>(std::get<2>(ret)) += IGIGn * (-1. + t0) * (-1. + t1);
   // }
   /* ------------------------------------ */
   // double t0, t1, w0w1;
   // for_each(__GW5xGW5__, [&](const auto &GWGW) {
   //    t0 = std::get<0>(GWGW);
   //    t1 = std::get<1>(GWGW);
   //    w0w1 = std::get<2>(GWGW);
   //    tmp = w0w1 * (1. - t0) / (nr = Norm(X0 * t0 + X1 * t1 * (1. - t0) - A));
   //    IGIGn = {tmp, tmp / (nr * nr)};
   //    std::get<1>(std::get<0>(ret)) += IGIGn * t0;
   //    std::get<1>(std::get<1>(ret)) += IGIGn * t1 * (1. - t0);
   //    std::get<1>(std::get<2>(ret)) += IGIGn * (-1. + t0) * (-1. + t1);
   // });
   /* ------------------------------------ */
   for_each(__GW5xGW5__, [&](const auto &GWGW) {
      tmp = std::get<2>(GWGW) * (1. - std::get<0>(GWGW)) / (nr = Norm(X0 * std::get<0>(GWGW) + X1 * std::get<1>(GWGW) * (1. - std::get<0>(GWGW)) - A));
      std::get<1>(std::get<0>(ret)) += (IGIGn = {tmp, tmp / (nr * nr)}) * std::get<0>(GWGW);
      std::get<1>(std::get<1>(ret)) += IGIGn * std::get<1>(GWGW) * (1. - std::get<0>(GWGW));
      std::get<1>(std::get<2>(ret)) += IGIGn * (-1. + std::get<0>(GWGW)) * (-1. + std::get<1>(GWGW));
   });
   /* ------------------------------------ */
   std::get<1>(std::get<0>(ret)) *= c;
   std::get<1>(std::get<1>(ret)) *= c;
   std::get<1>(std::get<2>(ret)) *= c;
   return ret;
};
/* ------------------------------------------------------ */
std::tuple<std::tuple<netP *, Tdd>,
           std::tuple<netP *, Tdd>,
           std::tuple<netP *, Tdd>>
calc_P_IGIGnLinearTupleTuple(const networkFace *const f, const networkPoint *const origin) {
   // おそらくIgIgnの計算精度に大きく関わるため，
   auto [p0, p1, p2] = f->getPoints(origin);
   interpolationTriangleLinearByFixedRange3D intp(T3Tddd{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
   Tdd IGIGn;
   std::tuple<std::tuple<netP *, Tdd>,
              std::tuple<netP *, Tdd>,
              std::tuple<netP *, Tdd>>
       ret = {{p0, {0., 0.}}, {p1, {0., 0.}}, {p2, {0., 0.}}};
   Tddd r, N, A = origin->getXtuple();
   double nr;
   bool isSingular = (p0 == origin);
   for (const auto &[x0, x1, w0w1] : isSingular ? __GWGW_beta2_10__ : __GWGW8__Tuple) {
      std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / (nr = Norm(r = intp(x0, x1) - A));
      std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
      std::get<1>(std::get<0>(ret)) += IGIGn * std::get<0>(N = intp.N(x0, x1));
      std::get<1>(std::get<1>(ret)) += IGIGn * std::get<1>(N);
      std::get<1>(std::get<2>(ret)) += IGIGn * std::get<2>(N);
   }
   return ret;
};

// std::vector<std::tuple<netP *, Tdd>> calc_P_IGIGnLinearTuple(const networkFace *const f, const networkFace *const origin_f)
// {
// 	// auto [P0, P1, P2] = f->getPoints();
// 	// auto [p0, p1, p2] = (P0->CORNER) ? f->getPoints(P0) : ((P1->CORNER) ? f->getPoints(P1) : f->getPoints(P2));
// 	auto [p0, p1, p2] = f->getPoints();
// 	Tdd p0igign = {0, 0};
// 	Tdd p1igign = {0, 0};
// 	Tdd p2igign = {0, 0};
// 	int i = 5;

// 	double mean_len = 5. * Max(extLength(origin_f->getLines()));
// 	if (Norm(p0->getXtuple() - origin_f->getXtuple()) < mean_len ||
// 		Norm(p1->getXtuple() - origin_f->getXtuple()) < mean_len ||
// 		Norm(p2->getXtuple() - origin_f->getXtuple()) < mean_len)
// 		i = 14;
// 	// タプルに対応させる．
// 	interpolationTriangleLinearByFixedRange3D intp(T3Tddd{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()});
// 	// Tddd A = intp(0.25, 0.25);
// 	Tddd A = origin_f->getXtuple(); // p0->getXtuple() + p1->getXtuple() + p2->getXtuple();
// 	// std::cout << "A_ = " << A_ << ", A = " << A << ", intp.N(x0, x1)" << intp.N(0.5, 0.5) << std::endl;
// 	// interpolationTriangleLinearByFixedRange3D intp_normal(T3Tddd{p0->getNormalAreaAveraged(), p1->getNormalAreaAveraged(), p2->getNormalAreaAveraged()});
// 	// interpolationTriangleLinearByFixedRange3D intp_normal(T3Tddd{p0->getNormalTuple(), p1->getNormalTuple(), p2->getNormalTuple()});
// 	Tdd IGIGn;
// 	Tddd r, N; // origin_f->getXtuple();
// 	// auto gw = GaussianQuadratureWeights(i, 0., 1.);
// 	// auto gw = GaussianQuadratureWeights(i, 0., 1.);
// 	double common, nr;
// 	double x0, x1, w0, w1;
// 	// bool isSingular = false; //(ps[0] == origin);  //!linearの場合は0
// 	bool isSingular = false; //(f == origin_f); //! linearの場合は0
// 	double beta = 2.;
// 	double x_sing = 0.5; //! linearの場合はx_sing = 1.
// 	for (const auto &x0w0 : isSingular ? GaussianQuadratureWeightsTuple(i, InvSg(0., x_sing, beta), InvSg(1., x_sing, beta)) : __GW__Tuple[i])
// 	{
// 		if (isSingular)
// 		{
// 			x0 = Sg(std::get<0>(x0w0), x_sing, beta);
// 			w0 = std::get<1>(x0w0) * DSg(std::get<0>(x0w0), x_sing, beta);
// 		}
// 		else
// 		{
// 			x0 = std::get<0>(x0w0);
// 			w0 = std::get<1>(x0w0);
// 		}
// 		for (const auto &x1w1 : __GW__Tuple[i])
// 		{
// 			x1 = std::get<0>(x1w1);
// 			w1 = std::get<1>(x1w1);
// 			/* ------------------------------ */
// 			common = intp.J(x0, x1) * w0 * w1;
// 			r = intp(x0, x1) - A;
// 			nr = Norm(r);
// 			std::get<0>(IGIGn) = common / nr;
// 			//! Dombre2019ならこっち
// 			// std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0, x1))) * common;
// 			std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0 * w1;
// 			// std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), Normalize(intp_normal(x0, x1))) * common;
// 			N = intp.N(x0, x1);
// 			p0igign += IGIGn * std::get<0>(N);
// 			p1igign += IGIGn * std::get<1>(N);
// 			p2igign += IGIGn * std::get<2>(N);
// 		}
// 	}
// 	// if (!isFinite(p0igign) || !isFinite(p1igign) || !isFinite(p2igign))
// 	// {
// 	// 	std::cout << "origin = " << origin << std::endl;
// 	// 	std::cout << "f = " << f << ", f->getPoints() = " << f->getPoints() << std::endl;
// 	// 	std::cout << "origin_f = " << origin_f << "origin_f->getPoints() = " << origin_f->getPoints() << std::endl;
// 	// 	std::cout << "[p0, p1, p2] = " << p0 << ", " << p1 << ", " << p2 << std::endl;
// 	// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	// }
// 	return {{p0, p0igign}, {p1, p1igign}, {p2, p2igign}};
// };

/* ------------------------------------------------------ */
std::vector<std::tuple<networkPoint *, std::unordered_map<Tii /*k,m*/, Tdd /*YYn*/>>> calc_P_MomentQuadTuple(const networkFace *const f, const Tddd &A) {
   /*
   *---*---*
    \ / \ /
     *-l-*
    / \ / \
   *---*---*
   */
   auto [p0, p1, p2] = f->getPoints();
   std::unordered_map<Tii, Tdd> p0igign;
   std::unordered_map<Tii, Tdd> p1igign;
   std::unordered_map<Tii, Tdd> p2igign;
   std::unordered_map<Tii, Tdd> p01igign;
   std::unordered_map<Tii, Tdd> p12igign;
   std::unordered_map<Tii, Tdd> p20igign;
   int i = 10;
   /*
           0*
   f0  / \  f2
      /   \
    1*-----*2
            f1
   */
   auto l0 = p0->getLineBetween(p1);
   auto l1 = p1->getLineBetween(p2);
   auto l2 = p2->getLineBetween(p0);
   //
   auto fs0 = l0->getFaces();
   auto ps6_l0_f00 = fs0[0]->get6PointsTuple(l0);
   auto ps6_l0_f01 = fs0[1]->get6PointsTuple(l0);
   interpolationTriangleQuadByFixedRange3D intp_l0_0(ToX(ps6_l0_f00));
   interpolationTriangleQuadByFixedRange3D intp_l0_1(ToX(ps6_l0_f01));
   //
   auto fs1 = l1->getFaces();
   auto ps6_l1_f10 = fs1[0]->get6PointsTuple(l1);
   auto ps6_l1_f11 = fs1[1]->get6PointsTuple(l1);
   interpolationTriangleQuadByFixedRange3D intp_l1_0(ToX(ps6_l1_f10));
   interpolationTriangleQuadByFixedRange3D intp_l1_1(ToX(ps6_l1_f11));
   //
   auto fs2 = l2->getFaces();
   auto ps6_l2_f20 = fs2[0]->get6PointsTuple(l2);
   auto ps6_l2_f21 = fs2[1]->get6PointsTuple(l2);
   interpolationTriangleQuadByFixedRange3D intp_l2_0(ToX(ps6_l2_f20));
   interpolationTriangleQuadByFixedRange3D intp_l2_1(ToX(ps6_l2_f21));
   //
   auto [f0p0, f0p1, f0p2] = (*l0)(f)->getPoints();
   auto [f1p0, f1p1, f1p2] = (*l1)(f)->getPoints();
   auto [f2p0, f2p1, f2p2] = (*l2)(f)->getPoints();
   interpolationTriangleQuadByFixedRange3D intp(
       T6Tddd{p0->getXtuple(),
              p1->getXtuple(),
              p2->getXtuple(),
              (intp_l0_0(.5, .5) + intp_l0_1(.5, .5)) / 2.,
              (intp_l1_0(.5, .5) + intp_l1_1(.5, .5)) / 2.,
              (intp_l2_0(.5, .5) + intp_l2_1(.5, .5)) / 2.});
   Tdd IGIGn;
   Tddd r;
   T6d N;
   double tmp;
   int max_k = 2;
   // 初期化
   for (auto k = 0; k <= max_k; ++k)
      for (auto m = -k; m <= k; ++m) {
         p0igign[{k, m}] = {0., 0.};
         p1igign[{k, m}] = {0., 0.};
         p2igign[{k, m}] = {0., 0.};
         p01igign[{k, m}] = {0., 0.};
         p12igign[{k, m}] = {0., 0.};
         p20igign[{k, m}] = {0., 0.};
      }
   //
   for (auto k = 0; k <= max_k; ++k)
      for (auto m = -k; m <= k; ++m)
         for (const auto &[x0, x1, w0w1] : __GWGW5__Tuple) {
            // nr = Norm(r = intp(x0, x1) - A);
            // std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
            // std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;

            r = intp(x0, x1) - A;
            auto [nr, theta, psi] = ToSphericalCoodrinates(r);

            std::get<0>(IGIGn) = std::pow(nr, k) * real_sph_scale_ommited(k, -m, theta, psi) * (intp.J(x0, x1) * w0w1);
            tmp = k * std::pow(nr, k - 1) * real_sph_scale_ommited(k, -m, theta, psi);
            // std::get<1>(IGIGn) = 1./nr*Dot(Ox,n)*(intp.J(x0, x1) * w0w1);
            std::get<1>(IGIGn) = Dot(tmp * r / nr, intp.cross(x0, x1)) * w0w1;

            N = intp.N(x0, x1);
            p0igign[{k, m}] += IGIGn * std::get<0>(N);
            p1igign[{k, m}] += IGIGn * std::get<1>(N);
            p2igign[{k, m}] += IGIGn * std::get<2>(N);
            p01igign[{k, m}] += IGIGn * std::get<3>(N);
            p12igign[{k, m}] += IGIGn * std::get<4>(N);
            p20igign[{k, m}] += IGIGn * std::get<5>(N);
         }

   auto N6_l0_0 = intp_l0_0.N(.5, .5) / 2.;
   auto N6_l0_1 = intp_l0_1.N(.5, .5) / 2.;
   auto N6_l1_0 = intp_l1_0.N(.5, .5) / 2.;
   auto N6_l1_1 = intp_l1_1.N(.5, .5) / 2.;
   auto N6_l2_0 = intp_l2_0.N(.5, .5) / 2.;
   auto N6_l2_1 = intp_l2_1.N(.5, .5) / 2.;

   return {
       {p0, p0igign},
       {p1, p1igign},
       {p2, p2igign},
       //
       {std::get<0>(ps6_l0_f00), p01igign * std::get<0>(N6_l0_0) / 2.},
       {std::get<1>(ps6_l0_f00), p01igign * std::get<1>(N6_l0_0) / 2.},
       {std::get<2>(ps6_l0_f00), p01igign * std::get<2>(N6_l0_0) / 2.},
       {std::get<3>(ps6_l0_f00), p01igign * std::get<3>(N6_l0_0) / 2.},
       {std::get<4>(ps6_l0_f00), p01igign * std::get<4>(N6_l0_0) / 2.},
       {std::get<5>(ps6_l0_f00), p01igign * std::get<5>(N6_l0_0) / 2.},
       //
       {std::get<0>(ps6_l0_f01), p01igign * std::get<0>(N6_l0_1) / 2.},
       {std::get<1>(ps6_l0_f01), p01igign * std::get<1>(N6_l0_1) / 2.},
       {std::get<2>(ps6_l0_f01), p01igign * std::get<2>(N6_l0_1) / 2.},
       {std::get<3>(ps6_l0_f01), p01igign * std::get<3>(N6_l0_1) / 2.},
       {std::get<4>(ps6_l0_f01), p01igign * std::get<4>(N6_l0_1) / 2.},
       {std::get<5>(ps6_l0_f01), p01igign * std::get<5>(N6_l0_1) / 2.},
       //
       {std::get<0>(ps6_l1_f10), p12igign * std::get<0>(N6_l1_0) / 2.},
       {std::get<1>(ps6_l1_f10), p12igign * std::get<1>(N6_l1_0) / 2.},
       {std::get<2>(ps6_l1_f10), p12igign * std::get<2>(N6_l1_0) / 2.},
       {std::get<3>(ps6_l1_f10), p12igign * std::get<3>(N6_l1_0) / 2.},
       {std::get<4>(ps6_l1_f10), p12igign * std::get<4>(N6_l1_0) / 2.},
       {std::get<5>(ps6_l1_f10), p12igign * std::get<5>(N6_l1_0) / 2.},
       //
       {std::get<0>(ps6_l1_f11), p12igign * std::get<0>(N6_l1_1) / 2.},
       {std::get<1>(ps6_l1_f11), p12igign * std::get<1>(N6_l1_1) / 2.},
       {std::get<2>(ps6_l1_f11), p12igign * std::get<2>(N6_l1_1) / 2.},
       {std::get<3>(ps6_l1_f11), p12igign * std::get<3>(N6_l1_1) / 2.},
       {std::get<4>(ps6_l1_f11), p12igign * std::get<4>(N6_l1_1) / 2.},
       {std::get<5>(ps6_l1_f11), p12igign * std::get<5>(N6_l1_1) / 2.},
       //
       {std::get<0>(ps6_l2_f20), p20igign * std::get<0>(N6_l2_0) / 2.},
       {std::get<1>(ps6_l2_f20), p20igign * std::get<1>(N6_l2_0) / 2.},
       {std::get<2>(ps6_l2_f20), p20igign * std::get<2>(N6_l2_0) / 2.},
       {std::get<3>(ps6_l2_f20), p20igign * std::get<3>(N6_l2_0) / 2.},
       {std::get<4>(ps6_l2_f20), p20igign * std::get<4>(N6_l2_0) / 2.},
       {std::get<5>(ps6_l2_f20), p20igign * std::get<5>(N6_l2_0) / 2.},
       //
       {std::get<0>(ps6_l2_f21), p20igign * std::get<0>(N6_l2_1) / 2.},
       {std::get<1>(ps6_l2_f21), p20igign * std::get<1>(N6_l2_1) / 2.},
       {std::get<2>(ps6_l2_f21), p20igign * std::get<2>(N6_l2_1) / 2.},
       {std::get<3>(ps6_l2_f21), p20igign * std::get<3>(N6_l2_1) / 2.},
       {std::get<4>(ps6_l2_f21), p20igign * std::get<4>(N6_l2_1) / 2.},
       {std::get<5>(ps6_l2_f21), p20igign * std::get<5>(N6_l2_1) / 2.}};
};
/* ------------------------------------------------------ */
std::vector<std::tuple<netP *, Tdd>> calc_P_IGIGnQuadTuple_mod(const networkFace *const f, const networkPoint *const origin) {
   /*
   *---*---*
    \ / \ /
     *-l-*
    / \ / \
   *---*---*
   */

   auto [p0, p1, p2] = f->getPoints(origin);
   Tdd p0igign = {0, 0};
   Tdd p1igign = {0, 0};
   Tdd p2igign = {0, 0};

   Tdd p01igign = {0, 0};
   Tdd p12igign = {0, 0};
   Tdd p20igign = {0, 0};

   // auto l0 = p0->getLineBetween(p1);
   // auto l1 = p1->getLineBetween(p2);
   // auto l2 = p2->getLineBetween(p0);

   auto [l0, l1, l2] = f->getLinesTupleFrom(p0);
   //
   // 幾何学的補間については，
   // X_surfaceがすでに角点については考慮したものとなっているので．
   // 角点に関して特に修正はいらない．
   //@ φやφnの補間に関しては，角点で滑らかでないことを考慮する必要があるが，φnに関しては，幾何学形状とは違って，同じ場所で2つの異なる値を持っている．
   interpolationTriangleQuadByFixedRange3D intp(
       T6Tddd{p0->getXtuple(),
              p1->getXtuple(),
              p2->getXtuple(),
              l0->X_surface,
              l1->X_surface,
              l2->X_surface});

   Tdd IGIGn;
   Tddd r, A = origin->getXtuple();
   T6d N;
   double nr;
   // bool isSingular = (p0 == origin); //! linearの場合は0

   /*
                   p0*
     f0   /   \
      p01      p20 f2
           /         \
     p1*-- p12 ---*p2
                    f1
    @ まずは，シンプルに上の６点p0,p1,p2,p01,p12,p20に掛かる係数を計算する
   */

   /*
   @ p01,p12,p20は実際には存在しない点
   @ 上で求めたそこでの係数は，p01,p12,p20を存在する点で補間し，その存在する点に分配する．
   @ ここでは，p01,p12,p20は，２つの２次補間の平均で表現する．
   */
   // auto F = !l0->isGoodForQuadInterp() ? f : (*l0)(f);
   auto F = !l0->isGoodForQuadInterp() ? f : (*l0)(f);
   // auto F = !l0->isGoodForQuadInterp_Geo() ? f : (*l0)(f);
   // auto F = (*l0)(f);
   interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_0(f, l0);
   interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_1(F, l0);

   F = !l1->isGoodForQuadInterp() ? f : (*l1)(f);
   // F = (*l1)(f);
   interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l1_0(f, l1);
   interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l1_1(F, l1);

   F = !l2->isGoodForQuadInterp() ? f : (*l2)(f);
   // F = (*l2)(f);
   interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l2_0(f, l2);
   interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l2_1(F, l2);

   if (p0 == origin)
      for (const auto &[x0, x1, w0w1] : __GWGW_beta2_11__) {
         nr = Norm(r = intp(x0, x1) - A);
         std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
         std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
         N = intp.N(x0, x1);
         p0igign += IGIGn * std::get<0>(N);  //@ 係数の計算
         p1igign += IGIGn * std::get<1>(N);
         p2igign += IGIGn * std::get<2>(N);
         p01igign += IGIGn * std::get<3>(N);
         p12igign += IGIGn * std::get<4>(N);
         p20igign += IGIGn * std::get<5>(N);
      }
   else if (MemberQ(intp_l0_0.Points, origin) ||
            MemberQ(intp_l0_1.Points, origin) ||
            MemberQ(intp_l1_0.Points, origin) ||
            MemberQ(intp_l1_1.Points, origin) ||
            MemberQ(intp_l2_0.Points, origin) ||
            MemberQ(intp_l2_1.Points, origin))
      for (const auto &[x0, x1, w0w1] : __GWGW11__Tuple) {
         nr = Norm(r = intp(x0, x1) - A);
         std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
         std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
         N = intp.N(x0, x1);
         p0igign += IGIGn * std::get<0>(N);
         p1igign += IGIGn * std::get<1>(N);
         p2igign += IGIGn * std::get<2>(N);
         p01igign += IGIGn * std::get<3>(N);
         p12igign += IGIGn * std::get<4>(N);
         p20igign += IGIGn * std::get<5>(N);
      }
   else
      for (const auto &[x0, x1, w0w1] : __GWGW8__Tuple) {
         nr = Norm(r = intp(x0, x1) - A);
         std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
         std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
         N = intp.N(x0, x1);
         p0igign += IGIGn * std::get<0>(N);
         p1igign += IGIGn * std::get<1>(N);
         p2igign += IGIGn * std::get<2>(N);
         p01igign += IGIGn * std::get<3>(N);
         p12igign += IGIGn * std::get<4>(N);
         p20igign += IGIGn * std::get<5>(N);
      }

   auto N6_l0_0 = intp_l0_0.N(.5, .5) * 0.5;
   auto N6_l0_1 = intp_l0_1.N(.5, .5) * 0.5;
   auto N6_l1_0 = intp_l1_0.N(.5, .5) * 0.5;
   auto N6_l1_1 = intp_l1_1.N(.5, .5) * 0.5;
   auto N6_l2_0 = intp_l2_0.N(.5, .5) * 0.5;
   auto N6_l2_1 = intp_l2_1.N(.5, .5) * 0.5;

   return {
       {p0, p0igign},
       {p1, p1igign},
       {p2, p2igign},
       /* ------------------------------------------------------ */
       {std::get<0>(intp_l0_0.Points), p01igign * std::get<0>(N6_l0_0)},
       {std::get<1>(intp_l0_0.Points), p01igign * std::get<1>(N6_l0_0)},
       {std::get<2>(intp_l0_0.Points), p01igign * std::get<2>(N6_l0_0)},
       {std::get<3>(intp_l0_0.Points), p01igign * std::get<3>(N6_l0_0)},
       {std::get<4>(intp_l0_0.Points), p01igign * std::get<4>(N6_l0_0)},
       {std::get<5>(intp_l0_0.Points), p01igign * std::get<5>(N6_l0_0)},
       //
       {std::get<0>(intp_l0_1.Points), p01igign * std::get<0>(N6_l0_1)},
       {std::get<1>(intp_l0_1.Points), p01igign * std::get<1>(N6_l0_1)},
       {std::get<2>(intp_l0_1.Points), p01igign * std::get<2>(N6_l0_1)},
       {std::get<3>(intp_l0_1.Points), p01igign * std::get<3>(N6_l0_1)},
       {std::get<4>(intp_l0_1.Points), p01igign * std::get<4>(N6_l0_1)},
       {std::get<5>(intp_l0_1.Points), p01igign * std::get<5>(N6_l0_1)},
       /* ------------------------------------------------------ */
       {std::get<0>(intp_l1_0.Points), p12igign * std::get<0>(N6_l1_0)},
       {std::get<1>(intp_l1_0.Points), p12igign * std::get<1>(N6_l1_0)},
       {std::get<2>(intp_l1_0.Points), p12igign * std::get<2>(N6_l1_0)},
       {std::get<3>(intp_l1_0.Points), p12igign * std::get<3>(N6_l1_0)},
       {std::get<4>(intp_l1_0.Points), p12igign * std::get<4>(N6_l1_0)},
       {std::get<5>(intp_l1_0.Points), p12igign * std::get<5>(N6_l1_0)},
       //
       {std::get<0>(intp_l1_1.Points), p12igign * std::get<0>(N6_l1_1)},
       {std::get<1>(intp_l1_1.Points), p12igign * std::get<1>(N6_l1_1)},
       {std::get<2>(intp_l1_1.Points), p12igign * std::get<2>(N6_l1_1)},
       {std::get<3>(intp_l1_1.Points), p12igign * std::get<3>(N6_l1_1)},
       {std::get<4>(intp_l1_1.Points), p12igign * std::get<4>(N6_l1_1)},
       {std::get<5>(intp_l1_1.Points), p12igign * std::get<5>(N6_l1_1)},
       /* ------------------------------------------------------ */
       {std::get<0>(intp_l2_0.Points), p20igign * std::get<0>(N6_l2_0)},
       {std::get<1>(intp_l2_0.Points), p20igign * std::get<1>(N6_l2_0)},
       {std::get<2>(intp_l2_0.Points), p20igign * std::get<2>(N6_l2_0)},
       {std::get<3>(intp_l2_0.Points), p20igign * std::get<3>(N6_l2_0)},
       {std::get<4>(intp_l2_0.Points), p20igign * std::get<4>(N6_l2_0)},
       {std::get<5>(intp_l2_0.Points), p20igign * std::get<5>(N6_l2_0)},
       //
       {std::get<0>(intp_l2_1.Points), p20igign * std::get<0>(N6_l2_1)},
       {std::get<1>(intp_l2_1.Points), p20igign * std::get<1>(N6_l2_1)},
       {std::get<2>(intp_l2_1.Points), p20igign * std::get<2>(N6_l2_1)},
       {std::get<3>(intp_l2_1.Points), p20igign * std::get<3>(N6_l2_1)},
       {std::get<4>(intp_l2_1.Points), p20igign * std::get<4>(N6_l2_1)},
       {std::get<5>(intp_l2_1.Points), p20igign * std::get<5>(N6_l2_1)}};
};

/* ------------------------------------------------------ */
// std::vector<std::tuple<netP *, Tdd>> calc_P_IGIGnQuadTuple_mod(const networkFace *const f, const networkPoint *const origin)
// {
// 	/*
// 	*---*---*
// 	 \ / \ /
// 	  *-l-*
// 	 / \ / \
	// 	*---*---*
// 	*/

// 	auto [p0, p1, p2] = f->getPoints(origin);
// 	Tdd p0igign = {0, 0};
// 	Tdd p1igign = {0, 0};
// 	Tdd p2igign = {0, 0};

// 	Tdd p01igign = {0, 0};
// 	Tdd p12igign = {0, 0};
// 	Tdd p20igign = {0, 0};

// 	// int i = 10;

// 	/*
// 		0*
// 	f0  / \  f2
// 	   /   \
	// 	 1*-----*2
// 		 f1
// 	*/
// 	auto l0 = p0->getLineBetween(p1);
// 	auto l1 = p1->getLineBetween(p2);
// 	auto l2 = p2->getLineBetween(p0);
// 	//
// 	// 幾何学的補間については，
// 	// X_surfaceがすでに角点については考慮したものとなっているので．
// 	// 角点に関して特に修正はいらない．
// 	//@ φやφnの補間に関しては，角点で滑らかでないことを考慮する必要があるが，φnに関しては，幾何学形状とは違って，同じ場所で2つの異なる値を持っている．
// 	interpolationTriangleQuadByFixedRange3D intp(
// 		T6Tddd{p0->getXtuple(),
// 			   p1->getXtuple(),
// 			   p2->getXtuple(),
// 			   l0->X_surface,
// 			   l1->X_surface,
// 			   l2->X_surface});

// 	Tdd IGIGn;
// 	Tddd r, A = origin->getXtuple();
// 	T6d N;
// 	double common, nr;
// 	double x0, x1, w0, w1;
// 	bool isSingular = (p0 == origin); //! linearの場合は0

// 	if (isSingular)
// 		for (const auto &[x0, x1, w0w1] : __GWGW_beta2_10__)
// 		{
// 			nr = Norm(r = intp(x0, x1) - A);
// 			std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
// 			std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
// 			N = intp.N(x0, x1);
// 			p0igign += IGIGn * std::get<0>(N);
// 			p1igign += IGIGn * std::get<1>(N);
// 			p2igign += IGIGn * std::get<2>(N);
// 			p01igign += IGIGn * std::get<3>(N);
// 			p12igign += IGIGn * std::get<4>(N);
// 			p20igign += IGIGn * std::get<5>(N);
// 		}
// 	else
// 		for (const auto &[x0, x1, w0w1] : __GWGW5__Tuple)
// 		{
// 			nr = Norm(r = intp(x0, x1) - A);
// 			std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
// 			std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
// 			N = intp.N(x0, x1);
// 			p0igign += IGIGn * std::get<0>(N);
// 			p1igign += IGIGn * std::get<1>(N);
// 			p2igign += IGIGn * std::get<2>(N);
// 			p01igign += IGIGn * std::get<3>(N);
// 			p12igign += IGIGn * std::get<4>(N);
// 			p20igign += IGIGn * std::get<5>(N);
// 		}

// 	// 	bool l0_good = !l0->CORNER && l0->islegal() && !(l0->Neumann && !l0->isFlat(M_PI / 2));
// 	// 	bool l1_good = !l1->CORNER && l1->islegal() && !(l1->Neumann && !l1->isFlat(M_PI / 2));
// 	// 	bool l2_good = !l2->CORNER && l2->islegal() && !(l2->Neumann && !l2->isFlat(M_PI / 2));

// 	// 	auto ps6_l0_f00 = f->get6PointsTuple(l0);
// 	// 	auto ps6_l0_f01 = (l0_good ? f : (*l0)(f))->get6PointsTuple(l0);
// 	// 	interpolationTriangleQuadByFixedRange3D intp_l0_0(ToX(ps6_l0_f00));
// 	// 	interpolationTriangleQuadByFixedRange3D intp_l0_1(ToX(ps6_l0_f01));
// 	// 	//
// 	// 	auto ps6_l1_f10 = f->get6PointsTuple(l1);
// 	// 	auto ps6_l1_f11 = (l1_good ? f : (*l1)(f))->get6PointsTuple(l1);
// 	// 	interpolationTriangleQuadByFixedRange3D intp_l1_0(ToX(ps6_l1_f10));
// 	// 	interpolationTriangleQuadByFixedRange3D intp_l1_1(ToX(ps6_l1_f11));
// 	// 	//
// 	// 	auto ps6_l2_f20 = f->get6PointsTuple(l2);
// 	// 	auto ps6_l2_f21 = (l2_good ? f : (*l2)(f))->get6PointsTuple(l2);
// 	// 	interpolationTriangleQuadByFixedRange3D intp_l2_0(ToX(ps6_l2_f20));
// 	// 	interpolationTriangleQuadByFixedRange3D intp_l2_1(ToX(ps6_l2_f21));
// 	// 	//
// 	// 	/*
// 	// 	 approx 0
// 	// 	   Q1 for l2
// 	// 	   Q2 for l1
// 	// 	*--Q0 for l0--*
// 	// 	 \   /   \   /
// 	// 	  \/      \ /
// 	// 	p1*---l0---*p0
// 	// 	  / \  f  / \
	// // 	 /   \   /   \
	// //   Q1*------*------*Q2
// 	// 		   p2
// 	// 	この修正によって，角ではphi，phinが不連続となる.
// 	// 	*/
// 	// 	if (l0_good)
// 	// 	{
// 	// 		intp_l0_0.approxP0();
// 	// 		intp_l0_1.approxP0();
// 	// 		intp_l1_0.approxP2();
// 	// 		intp_l1_1.approxP2();
// 	// 		intp_l2_0.approxP1();
// 	// 		intp_l2_1.approxP1();
// 	// 	}
// 	// 	if (l1_good)
// 	// 	{
// 	// 		intp_l0_0.approxP1();
// 	// 		intp_l0_1.approxP1();
// 	// 		intp_l1_0.approxP0();
// 	// 		intp_l1_1.approxP0();
// 	// 		intp_l2_0.approxP2();
// 	// 		intp_l2_1.approxP2();
// 	// 	}
// 	// 	if (l2_good)
// 	// 	{
// 	// 		intp_l0_0.approxP2();
// 	// 		intp_l0_1.approxP2();
// 	// 		intp_l1_0.approxP1();
// 	// 		intp_l1_1.approxP1();
// 	// 		intp_l2_0.approxP0();
// 	// 		intp_l2_1.approxP0();
// 	// 	}

// 	/* ------------------------------------------------------ */

// 	// auto F = !l0->isGoodForQuadInterp() ? f : (*l0)(f);
// 	auto F = !l0->isGoodForQuadInterp() ? f : (*l0)(f);
// 	// auto F = (*l0)(f);
// 	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_0(f, l0);
// 	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_1(F, l0);

// 	F = !l1->isGoodForQuadInterp() ? f : (*l1)(f);
// 	// F = (*l1)(f);
// 	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l1_0(f, l1);
// 	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l1_1(F, l1);

// 	F = !l2->isGoodForQuadInterp() ? f : (*l2)(f);
// 	// F = (*l2)(f);
// 	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l2_0(f, l2);
// 	interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l2_1(F, l2);

// 	auto N6_l0_0 = intp_l0_0.N(.5, .5) / 2.;
// 	auto N6_l0_1 = intp_l0_1.N(.5, .5) / 2.;
// 	auto N6_l1_0 = intp_l1_0.N(.5, .5) / 2.;
// 	auto N6_l1_1 = intp_l1_1.N(.5, .5) / 2.;
// 	auto N6_l2_0 = intp_l2_0.N(.5, .5) / 2.;
// 	auto N6_l2_1 = intp_l2_1.N(.5, .5) / 2.;

// 	return {
// 		{p0, p0igign},
// 		{p1, p1igign},
// 		{p2, p2igign},
// 		/* ------------------------------------------------------ */
// 		{std::get<0>(intp_l0_0.Points), p01igign * std::get<0>(N6_l0_0) / 2.},
// 		{std::get<1>(intp_l0_0.Points), p01igign * std::get<1>(N6_l0_0) / 2.},
// 		{std::get<2>(intp_l0_0.Points), p01igign * std::get<2>(N6_l0_0) / 2.},
// 		{std::get<3>(intp_l0_0.Points), p01igign * std::get<3>(N6_l0_0) / 2.},
// 		{std::get<4>(intp_l0_0.Points), p01igign * std::get<4>(N6_l0_0) / 2.},
// 		{std::get<5>(intp_l0_0.Points), p01igign * std::get<5>(N6_l0_0) / 2.},
// 		//
// 		{std::get<0>(intp_l0_1.Points), p01igign * std::get<0>(N6_l0_1) / 2.},
// 		{std::get<1>(intp_l0_1.Points), p01igign * std::get<1>(N6_l0_1) / 2.},
// 		{std::get<2>(intp_l0_1.Points), p01igign * std::get<2>(N6_l0_1) / 2.},
// 		{std::get<3>(intp_l0_1.Points), p01igign * std::get<3>(N6_l0_1) / 2.},
// 		{std::get<4>(intp_l0_1.Points), p01igign * std::get<4>(N6_l0_1) / 2.},
// 		{std::get<5>(intp_l0_1.Points), p01igign * std::get<5>(N6_l0_1) / 2.},
// 		/* ------------------------------------------------------ */
// 		{std::get<0>(intp_l1_0.Points), p12igign * std::get<0>(N6_l1_0) / 2.},
// 		{std::get<1>(intp_l1_0.Points), p12igign * std::get<1>(N6_l1_0) / 2.},
// 		{std::get<2>(intp_l1_0.Points), p12igign * std::get<2>(N6_l1_0) / 2.},
// 		{std::get<3>(intp_l1_0.Points), p12igign * std::get<3>(N6_l1_0) / 2.},
// 		{std::get<4>(intp_l1_0.Points), p12igign * std::get<4>(N6_l1_0) / 2.},
// 		{std::get<5>(intp_l1_0.Points), p12igign * std::get<5>(N6_l1_0) / 2.},
// 		//
// 		{std::get<0>(intp_l1_1.Points), p12igign * std::get<0>(N6_l1_1) / 2.},
// 		{std::get<1>(intp_l1_1.Points), p12igign * std::get<1>(N6_l1_1) / 2.},
// 		{std::get<2>(intp_l1_1.Points), p12igign * std::get<2>(N6_l1_1) / 2.},
// 		{std::get<3>(intp_l1_1.Points), p12igign * std::get<3>(N6_l1_1) / 2.},
// 		{std::get<4>(intp_l1_1.Points), p12igign * std::get<4>(N6_l1_1) / 2.},
// 		{std::get<5>(intp_l1_1.Points), p12igign * std::get<5>(N6_l1_1) / 2.},
// 		/* ------------------------------------------------------ */
// 		{std::get<0>(intp_l2_0.Points), p20igign * std::get<0>(N6_l2_0) / 2.},
// 		{std::get<1>(intp_l2_0.Points), p20igign * std::get<1>(N6_l2_0) / 2.},
// 		{std::get<2>(intp_l2_0.Points), p20igign * std::get<2>(N6_l2_0) / 2.},
// 		{std::get<3>(intp_l2_0.Points), p20igign * std::get<3>(N6_l2_0) / 2.},
// 		{std::get<4>(intp_l2_0.Points), p20igign * std::get<4>(N6_l2_0) / 2.},
// 		{std::get<5>(intp_l2_0.Points), p20igign * std::get<5>(N6_l2_0) / 2.},
// 		//
// 		{std::get<0>(intp_l2_1.Points), p20igign * std::get<0>(N6_l2_1) / 2.},
// 		{std::get<1>(intp_l2_1.Points), p20igign * std::get<1>(N6_l2_1) / 2.},
// 		{std::get<2>(intp_l2_1.Points), p20igign * std::get<2>(N6_l2_1) / 2.},
// 		{std::get<3>(intp_l2_1.Points), p20igign * std::get<3>(N6_l2_1) / 2.},
// 		{std::get<4>(intp_l2_1.Points), p20igign * std::get<4>(N6_l2_1) / 2.},
// 		{std::get<5>(intp_l2_1.Points), p20igign * std::get<5>(N6_l2_1) / 2.}};
// };
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
// std::vector<std::tuple<netP *, Tdd>> calc_P_IGIGnQuadTuple_mod(const networkFace *const f, const networkPoint *const origin)
// {
// 	/*
// 	*---*---*
// 	 \ / \ /
// 	  *-l-*
// 	 / \ / \
	// 	*---*---*
// 	*/

// 	auto [p0, p1, p2] = f->getPoints(origin);
// 	Tdd p0igign = {0, 0};
// 	Tdd p1igign = {0, 0};
// 	Tdd p2igign = {0, 0};

// 	Tdd p01igign = {0, 0};
// 	Tdd p12igign = {0, 0};
// 	Tdd p20igign = {0, 0};

// 	// int i = 10;

// 	/*
// 		0*
// 	f0  / \  f2
// 	   /   \
	// 	 1*-----*2
// 		 f1
// 	*/
// 	auto l0 = p0->getLineBetween(p1);
// 	auto l1 = p1->getLineBetween(p2);
// 	auto l2 = p2->getLineBetween(p0);
// 	//
// 	auto fs0 = l0->getFaces();
// 	auto ps6_l0_f00 = fs0[0]->get6PointsTuple(l0);
// 	auto ps6_l0_f01 = fs0[1]->get6PointsTuple(l0);
// 	interpolationTriangleQuadByFixedRange3D intp_l0_0(ToX(ps6_l0_f00));
// 	interpolationTriangleQuadByFixedRange3D intp_l0_1(ToX(ps6_l0_f01));
// 	//
// 	auto fs1 = l1->getFaces();
// 	auto ps6_l1_f10 = fs1[0]->get6PointsTuple(l1);
// 	auto ps6_l1_f11 = fs1[1]->get6PointsTuple(l1);
// 	interpolationTriangleQuadByFixedRange3D intp_l1_0(ToX(ps6_l1_f10));
// 	interpolationTriangleQuadByFixedRange3D intp_l1_1(ToX(ps6_l1_f11));
// 	//
// 	auto fs2 = l2->getFaces();
// 	auto ps6_l2_f20 = fs2[0]->get6PointsTuple(l2);
// 	auto ps6_l2_f21 = fs2[1]->get6PointsTuple(l2);
// 	interpolationTriangleQuadByFixedRange3D intp_l2_0(ToX(ps6_l2_f20));
// 	interpolationTriangleQuadByFixedRange3D intp_l2_1(ToX(ps6_l2_f21));
// 	//
// 	auto [f0p0, f0p1, f0p2] = (*l0)(f)->getPoints();
// 	auto [f1p0, f1p1, f1p2] = (*l1)(f)->getPoints();
// 	auto [f2p0, f2p1, f2p2] = (*l2)(f)->getPoints();
// 	//
// 	// interpolationTriangleQuadByFixedRange3D intp(
// 	// 	T6Tddd{p0->getXtuple(),
// 	// 		   p1->getXtuple(),
// 	// 		   p2->getXtuple(),
// 	// 		   (intp_l0_0(.5, .5) + intp_l0_1(.5, .5)) / 2.,
// 	// 		   (intp_l1_0(.5, .5) + intp_l1_1(.5, .5)) / 2.,
// 	// 		   (intp_l2_0(.5, .5) + intp_l2_1(.5, .5)) / 2.});

// 	// auto [L0, L1, L2] = f->getLinesTupleFrom(p0);

// 	// interpolationTriangleQuadByFixedRange3D intp(
// 	// 	T6Tddd{p0->getXtuple(),
// 	// 		   p1->getXtuple(),
// 	// 		   p2->getXtuple(),
// 	// 		   l0->CORNER ? (p0->getXtuple() + p1->getXtuple()) / 2. : (intp_l0_0(.5, .5) + intp_l0_1(.5, .5)) / 2.,
// 	// 		   l1->CORNER ? (p1->getXtuple() + p2->getXtuple()) / 2. : (intp_l1_0(.5, .5) + intp_l1_1(.5, .5)) / 2.,
// 	// 		   l2->CORNER ? (p2->getXtuple() + p0->getXtuple()) / 2. : (intp_l2_0(.5, .5) + intp_l2_1(.5, .5)) / 2.});

// 	// interpolationTriangleQuadByFixedRange3D intp(
// 	// 	T6Tddd{p0->getXtuple(),
// 	// 		   p1->getXtuple(),
// 	// 		   p2->getXtuple(),
// 	// 		   angle0 > M_PI / 3. ? (p0->getXtuple() + p1->getXtuple()) / 2. : (intp_l0_0(.5, .5) + intp_l0_1(.5, .5)) / 2.,
// 	// 		   angle1 > M_PI / 3. ? (p1->getXtuple() + p2->getXtuple()) / 2. : (intp_l1_0(.5, .5) + intp_l1_1(.5, .5)) / 2.,
// 	// 		   angle2 > M_PI / 3. ? (p2->getXtuple() + p0->getXtuple()) / 2. : (intp_l2_0(.5, .5) + intp_l2_1(.5, .5)) / 2.});

// 	// interpolationTriangleQuadByFixedRange3D intp(
// 	// 	T6Tddd{p0->getXtuple(),
// 	// 		   p1->getXtuple(),
// 	// 		   p2->getXtuple(),
// 	// 		   (p0->CORNER || p1->CORNER || angle0 > M_PI / 3.) ? (p0->getXtuple() + p1->getXtuple()) / 2. : (intp_l0_0(.5, .5) + intp_l0_1(.5, .5)) / 2.,
// 	// 		   (p2->CORNER || p1->CORNER || angle1 > M_PI / 3.) ? (p1->getXtuple() + p2->getXtuple()) / 2. : (intp_l1_0(.5, .5) + intp_l1_1(.5, .5)) / 2.,
// 	// 		   (p0->CORNER || p2->CORNER || angle2 > M_PI / 3.) ? (p2->getXtuple() + p0->getXtuple()) / 2. : (intp_l2_0(.5, .5) + intp_l2_1(.5, .5)) / 2.});

// 	interpolationTriangleQuadByFixedRange3D intp(
// 		T6Tddd{p0->getXtuple(),
// 			   p1->getXtuple(),
// 			   p2->getXtuple(),
// 			   l0->X_surface,
// 			   l1->X_surface,
// 			   l2->X_surface});

// 	// interpolationTriangleQuadByFixedRange3D intp(
// 	// 	T6Tddd{p0->getXtuple(),
// 	// 		   p1->getXtuple(),
// 	// 		   p2->getXtuple(),
// 	// 		   (p0->getXtuple() + p1->getXtuple()) / 2.,
// 	// 		   (p1->getXtuple() + p2->getXtuple()) / 2.,
// 	// 		   (p2->getXtuple() + p0->getXtuple()) / 2.});

// 	// interpolationTriangleQuadByFixedRange3D intp(
// 	// 	T6Tddd{p0->getXtuple(),
// 	// 		   p1->getXtuple(),
// 	// 		   p2->getXtuple(),
// 	// 		   (p0->getXtuple() + p1->getXtuple()) / 2.,
// 	// 		   (p1->getXtuple() + p2->getXtuple()) / 2.,
// 	// 		   (p2->getXtuple() + p0->getXtuple()) / 2.});

// 	// auto [intp_l0_0_,
// 	// 	  intp_l0_1_,
// 	// 	  intp_l1_0_,
// 	// 	  intp_l1_1_,
// 	// 	  intp_l2_0_,
// 	// 	  intp_l2_1_,
// 	// 	  intp_] = f->getQuadInterpolation(origin);

// 	// auto intp_l0_0 = *intp_l0_0_;
// 	// auto intp_l0_1 = *intp_l0_1_;
// 	// auto intp_l1_0 = *intp_l1_0_;
// 	// auto intp_l1_1 = *intp_l1_1_;
// 	// auto intp_l2_0 = *intp_l2_0_;
// 	// auto intp_l2_1 = *intp_l2_1_;
// 	// auto intp = *intp_;

// 	Tdd IGIGn;
// 	Tddd r, A = origin->getXtuple();
// 	T6d N;
// 	// auto gw = GaussianQuadratureWeights(i, 0., 1.);
// 	// auto gw = GaussianQuadratureWeights(i, 0., 1.);
// 	double common, nr;
// 	double x0, x1, w0, w1;
// 	// bool isSingular = false;		  //(ps[0] == origin);  //!linearの場合は0
// 	bool isSingular = (p0 == origin); //! linearの場合は0
// 	// double beta = 2.;
// 	// double x_sing = 1.; //! linearの場合はx_sing = 1.

// 	// for (const auto &x0w0 : isSingular ? GaussianQuadratureWeightsTuple(10, InvSg(0., x_sing, beta), InvSg(1., x_sing, beta)) : __GW5__Tuple)
// 	// {
// 	// 	if (isSingular)
// 	// 	{
// 	// 		x0 = Sg(std::get<0>(x0w0), x_sing, beta);
// 	// 		w0 = std::get<1>(x0w0) * DSg(std::get<0>(x0w0), x_sing, beta);
// 	// 	}
// 	// 	else
// 	// 	{
// 	// 		x0 = std::get<0>(x0w0);
// 	// 		w0 = std::get<1>(x0w0);
// 	// 	}
// 	// 	for (const auto &[x1, w1] : __GW5__Tuple)
// 	// 	{
// 	// 		// x1 = std::get<0>(x1w1);
// 	// 		// w1 = std::get<1>(x1w1);
// 	// 		/* ------------------------------ */
// 	// 		common = intp.J(x0, x1) * w0 * w1;
// 	// 		r = intp(x0, x1) - A;
// 	// 		nr = Norm(r);
// 	// 		std::get<0>(IGIGn) = common / nr;
// 	// 		// std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0, x1))) * common;
// 	// 		std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0 * w1;
// 	// 		// std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), Normalize(intp_normal(x0, x1))) * common;
// 	// 		N = intp.N(x0, x1);
// 	// 		p0igign += IGIGn * std::get<0>(N);
// 	// 		p1igign += IGIGn * std::get<1>(N);
// 	// 		p2igign += IGIGn * std::get<2>(N);
// 	// 		p01igign += IGIGn * std::get<3>(N);
// 	// 		p12igign += IGIGn * std::get<4>(N);
// 	// 		p20igign += IGIGn * std::get<5>(N);
// 	// 	}
// 	// }

// 	if (isSingular)
// 		for (const auto &[x0, x1, w0w1] : __GWGW_beta2_10__)
// 		{
// 			nr = Norm(r = intp(x0, x1) - A);
// 			std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
// 			std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
// 			N = intp.N(x0, x1);
// 			p0igign += IGIGn * std::get<0>(N);
// 			p1igign += IGIGn * std::get<1>(N);
// 			p2igign += IGIGn * std::get<2>(N);
// 			p01igign += IGIGn * std::get<3>(N);
// 			p12igign += IGIGn * std::get<4>(N);
// 			p20igign += IGIGn * std::get<5>(N);
// 		}
// 	else if (MemberQ(std::vector<networkPoint *>{f0p0, f0p1, f0p2, f1p0, f1p1, f1p2, f2p0, f2p1, f2p2}, origin))
// 	{
// 		for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
// 		{
// 			nr = Norm(r = intp(x0, x1) - A);
// 			std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
// 			std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
// 			N = intp.N(x0, x1);
// 			p0igign += IGIGn * std::get<0>(N);
// 			p1igign += IGIGn * std::get<1>(N);
// 			p2igign += IGIGn * std::get<2>(N);
// 			p01igign += IGIGn * std::get<3>(N);
// 			p12igign += IGIGn * std::get<4>(N);
// 			p20igign += IGIGn * std::get<5>(N);
// 		}
// 	}
// 	else
// 		for (const auto &[x0, x1, w0w1] : __GWGW7__Tuple)
// 		{
// 			nr = Norm(r = intp(x0, x1) - A);
// 			std::get<0>(IGIGn) = (intp.J(x0, x1) * w0w1) / nr;
// 			std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0w1;
// 			N = intp.N(x0, x1);
// 			p0igign += IGIGn * std::get<0>(N);
// 			p1igign += IGIGn * std::get<1>(N);
// 			p2igign += IGIGn * std::get<2>(N);
// 			p01igign += IGIGn * std::get<3>(N);
// 			p12igign += IGIGn * std::get<4>(N);
// 			p20igign += IGIGn * std::get<5>(N);
// 		}

// 	auto N6_l0_0 = intp_l0_0.N(.5, .5) / 2.;
// 	auto N6_l0_1 = intp_l0_1.N(.5, .5) / 2.;
// 	auto N6_l1_0 = intp_l1_0.N(.5, .5) / 2.;
// 	auto N6_l1_1 = intp_l1_1.N(.5, .5) / 2.;
// 	auto N6_l2_0 = intp_l2_0.N(.5, .5) / 2.;
// 	auto N6_l2_1 = intp_l2_1.N(.5, .5) / 2.;

// 	return {
// 		{p0, p0igign},
// 		{p1, p1igign},
// 		{p2, p2igign},
// 		//
// 		{std::get<0>(ps6_l0_f00), p01igign * std::get<0>(N6_l0_0) / 2.},
// 		{std::get<1>(ps6_l0_f00), p01igign * std::get<1>(N6_l0_0) / 2.},
// 		{std::get<2>(ps6_l0_f00), p01igign * std::get<2>(N6_l0_0) / 2.},
// 		{std::get<3>(ps6_l0_f00), p01igign * std::get<3>(N6_l0_0) / 2.},
// 		{std::get<4>(ps6_l0_f00), p01igign * std::get<4>(N6_l0_0) / 2.},
// 		{std::get<5>(ps6_l0_f00), p01igign * std::get<5>(N6_l0_0) / 2.},
// 		//
// 		{std::get<0>(ps6_l0_f01), p01igign * std::get<0>(N6_l0_1) / 2.},
// 		{std::get<1>(ps6_l0_f01), p01igign * std::get<1>(N6_l0_1) / 2.},
// 		{std::get<2>(ps6_l0_f01), p01igign * std::get<2>(N6_l0_1) / 2.},
// 		{std::get<3>(ps6_l0_f01), p01igign * std::get<3>(N6_l0_1) / 2.},
// 		{std::get<4>(ps6_l0_f01), p01igign * std::get<4>(N6_l0_1) / 2.},
// 		{std::get<5>(ps6_l0_f01), p01igign * std::get<5>(N6_l0_1) / 2.},
// 		//
// 		{std::get<0>(ps6_l1_f10), p12igign * std::get<0>(N6_l1_0) / 2.},
// 		{std::get<1>(ps6_l1_f10), p12igign * std::get<1>(N6_l1_0) / 2.},
// 		{std::get<2>(ps6_l1_f10), p12igign * std::get<2>(N6_l1_0) / 2.},
// 		{std::get<3>(ps6_l1_f10), p12igign * std::get<3>(N6_l1_0) / 2.},
// 		{std::get<4>(ps6_l1_f10), p12igign * std::get<4>(N6_l1_0) / 2.},
// 		{std::get<5>(ps6_l1_f10), p12igign * std::get<5>(N6_l1_0) / 2.},
// 		//
// 		{std::get<0>(ps6_l1_f11), p12igign * std::get<0>(N6_l1_1) / 2.},
// 		{std::get<1>(ps6_l1_f11), p12igign * std::get<1>(N6_l1_1) / 2.},
// 		{std::get<2>(ps6_l1_f11), p12igign * std::get<2>(N6_l1_1) / 2.},
// 		{std::get<3>(ps6_l1_f11), p12igign * std::get<3>(N6_l1_1) / 2.},
// 		{std::get<4>(ps6_l1_f11), p12igign * std::get<4>(N6_l1_1) / 2.},
// 		{std::get<5>(ps6_l1_f11), p12igign * std::get<5>(N6_l1_1) / 2.},
// 		//
// 		{std::get<0>(ps6_l2_f20), p20igign * std::get<0>(N6_l2_0) / 2.},
// 		{std::get<1>(ps6_l2_f20), p20igign * std::get<1>(N6_l2_0) / 2.},
// 		{std::get<2>(ps6_l2_f20), p20igign * std::get<2>(N6_l2_0) / 2.},
// 		{std::get<3>(ps6_l2_f20), p20igign * std::get<3>(N6_l2_0) / 2.},
// 		{std::get<4>(ps6_l2_f20), p20igign * std::get<4>(N6_l2_0) / 2.},
// 		{std::get<5>(ps6_l2_f20), p20igign * std::get<5>(N6_l2_0) / 2.},
// 		//
// 		{std::get<0>(ps6_l2_f21), p20igign * std::get<0>(N6_l2_1) / 2.},
// 		{std::get<1>(ps6_l2_f21), p20igign * std::get<1>(N6_l2_1) / 2.},
// 		{std::get<2>(ps6_l2_f21), p20igign * std::get<2>(N6_l2_1) / 2.},
// 		{std::get<3>(ps6_l2_f21), p20igign * std::get<3>(N6_l2_1) / 2.},
// 		{std::get<4>(ps6_l2_f21), p20igign * std::get<4>(N6_l2_1) / 2.},
// 		{std::get<5>(ps6_l2_f21), p20igign * std::get<5>(N6_l2_1) / 2.}};
// };
/* ------------------------------------------------------ */
std::vector<std::tuple<netP *, Tdd>> calc_P_IGIGnQuadTuple(const networkFace *const f, const networkPoint *const origin) {

   auto [p0, p1, p2] = f->getPoints(origin);
   Tdd p0igign = {0, 0};
   Tdd p1igign = {0, 0};
   Tdd p2igign = {0, 0};

   Tdd p01igign = {0, 0};
   Tdd p12igign = {0, 0};
   Tdd p20igign = {0, 0};

   int i = 6;

   /*
           0*
   f0  / \  f2
      /   \
    1*-----*2
            f1
   */
   auto l0 = p0->getLineBetween(p1);
   auto l1 = p1->getLineBetween(p2);
   auto l2 = p2->getLineBetween(p0);
   auto [f0p0, f0p1, f0p2] = (*l0)(f)->getPoints();
   auto [f1p0, f1p1, f1p2] = (*l1)(f)->getPoints();
   auto [f2p0, f2p1, f2p2] = (*l2)(f)->getPoints();
   interpolationTriangleQuadByFixedRange3D intp(
       T6Tddd{p0->getXtuple(),
              p1->getXtuple(),
              p2->getXtuple(),
              (p0->getXtuple() + p1->getXtuple()) / 2.,
              (p1->getXtuple() + p2->getXtuple()) / 2.,
              (p2->getXtuple() + p0->getXtuple()) / 2.});

   Tdd IGIGn;
   Tddd r, A = origin->getXtuple();
   T6d N;
   // auto gw = GaussianQuadratureWeights(i, 0., 1.);
   // auto gw = GaussianQuadratureWeights(i, 0., 1.);
   double common, nr;
   double x0, x1, w0, w1;
   // bool isSingular = false;		  //(ps[0] == origin);  //!linearの場合は0
   bool isSingular = (p0 == origin);  //! linearの場合は0
   double beta = 2.;
   double x_sing = 1.;  //! linearの場合はx_sing = 1.
   auto GW = __GW__Tuple[i];
   for (const auto &x0w0 : isSingular ? GaussianQuadratureWeightsTuple(6, InvSg(0., x_sing, beta), InvSg(1., x_sing, beta)) : GW) {
      if (isSingular) {
         x0 = Sg(std::get<0>(x0w0), x_sing, beta);
         w0 = std::get<1>(x0w0) * DSg(std::get<0>(x0w0), x_sing, beta);
      } else {
         x0 = std::get<0>(x0w0);
         w0 = std::get<1>(x0w0);
      }
      for (const auto &x1w1 : GW) {
         x1 = std::get<0>(x1w1);
         w1 = std::get<1>(x1w1);
         /* ------------------------------ */
         common = intp.J(x0, x1) * w0 * w1;
         r = intp(x0, x1) - A;
         nr = Norm(r);
         std::get<0>(IGIGn) = common / nr;
         // std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), Normalize(intp.cross(x0, x1))) * common;
         std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), intp.cross(x0, x1)) * w0 * w1;
         // std::get<1>(IGIGn) = -Dot(r / pow(nr, 3.), Normalize(intp_normal(x0, x1))) * common;
         N = intp.N(x0, x1);
         p0igign += IGIGn * std::get<0>(N);
         p1igign += IGIGn * std::get<1>(N);
         p2igign += IGIGn * std::get<2>(N);
         p01igign += IGIGn * std::get<3>(N);
         p12igign += IGIGn * std::get<4>(N);
         p20igign += IGIGn * std::get<5>(N);
      }
   }

   auto tmp = (p01igign + p12igign + p20igign) / 2. / 3.;

   return {
       // f
       {p0, p0igign + tmp},
       {p1, p1igign + tmp},
       {p2, p2igign + tmp},

       // f0
       {f0p0, p01igign / 2. / 3.},
       {f0p1, p01igign / 2. / 3.},
       {f0p2, p01igign / 2. / 3.},

       // f1
       {f1p0, p12igign / 2. / 3.},
       {f1p1, p12igign / 2. / 3.},
       {f1p2, p12igign / 2. / 3.},

       // f2
       {f2p0, p20igign / 2. / 3.},
       {f2p1, p20igign / 2. / 3.},
       {f2p2, p20igign / 2. / 3.}};

   // auto tmp = (p01igign / 9. + p12igign / 9. + p20igign / 9.);

   // return {
   // 	// f
   // 	{p0, p0igign + p01igign / 6. + p20igign / 6. + tmp},
   // 	{p1, p1igign + p01igign / 6. + p12igign / 6. + tmp},
   // 	{p2, p2igign + p12igign / 6. + p20igign / 6. + tmp},

   // 	// f0
   // 	{f0p0, p01igign / 9.},
   // 	{f0p1, p01igign / 9.},
   // 	{f0p2, p01igign / 9.},

   // 	// f1
   // 	{f1p0, p12igign / 9.},
   // 	{f1p1, p12igign / 9.},
   // 	{f1p2, p12igign / 9.},

   // 	// f2
   // 	{f2p0, p20igign / 9.},
   // 	{f2p1, p20igign / 9.},
   // 	{f2p2, p20igign / 9.}};
};

// /* ------------------------------------------------------ */
// std::unordered_map<netP *, Tdd> calc_P_IGIGnTuple(const netFp &f, const netPp origin)
// {
// 	V_netPp ps = f->getPoints(origin);
// 	V_i nGWGW = {4};
// 	double conv = 1E-2;
// 	Tdd zeros = {0, 0};
// 	std::unordered_map<netPp, Tdd> pre_igign;
// 	pre_igign.reserve(ps.size());
// 	for (const auto &p : ps)
// 		pre_igign[p] = zeros;

// 	std::unordered_map<netPp, Tdd> new_igign = pre_igign;

// 	p_igign_Linear(ps[0], ps[1], ps[2], origin, nGWGW[0], /*出力結果*/ pre_igign); //!まず計算
// 	new_igign = pre_igign;
// 	for (auto i = 1; i < nGWGW.size(); i++)
// 	{
// 		p_igign_Linear(ps[0], ps[1], ps[2], origin, nGWGW[i], /*出力結果*/ new_igign); //!次に計算
// 		if (isConverged(new_igign, pre_igign, ps, conv))
// 			return new_igign;
// 		else
// 			pre_igign = new_igign; //replace and try again
// 	}
// 	return new_igign;
// };
// //* ------------------------------------------------------ */
// //!               面における方程式を面上で作成                  */
// //* ------------------------------------------------------ */
// std::map<netP *, V_d> calc_P_IGIGn(const netFp &f /*積分面*/,
// 								   const netFp origin_f /*原点のある面*/,
// 								   double t0, double t1 /*原点の位置を指定するパラメタ座標*/,
// 								   int points_number)
// {
// 	/**
// 	 * 面を原点とした場合の境界積分方程式の一部：IGIGnを計算する
// 	 */
// 	// 明日はリメッシュについて考察．
// 	// まとめる．論文を読む．
// 	V_i gps = {8};
// 	if (f == origin_f)
// 		gps = {18};
// 	double conv = 1E-3;
// 	V_d zeros(2, 0.);
// 	std::map<netPp, V_d> pre_igign;
// 	V_netPp ps;
// 	/* ------------------------------------------------------ */
// 	if (points_number == 3)
// 		ps = f->getPoints(); //基準としてoriginを使おうとする．が，面に含まれない点の場合は，適当な面の点が内部で自動で選ばれる．
// 	else if (points_number == 6)
// 		ps = f->get6PointsTuple();
// 	// else if (points_number == 12)
// 	// 	ps = f->get12Points();
// 	//
// 	for (const auto &p : ps)
// 		pre_igign[p] = zeros;
// 	//
// 	p_igign(f, origin_f, t0, t1, gps[0], /*出力結果*/ pre_igign, ps); //!まず計算
// 	//
// 	std::map<netPp, V_d> new_igign = pre_igign;
// 	for (auto i = 1; i < gps.size(); i++)
// 	{
// 		p_igign(f, origin_f, t0, t1, gps[i], /*出力結果*/ new_igign, ps); //!次に計算
// 		if (isConverged(new_igign, pre_igign, ps, conv))
// 			return new_igign;
// 		else
// 			pre_igign = new_igign; // replace and try again
// 	}
// 	return new_igign;
// };

///////////////////////////////////

// namespace BEM
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using map_P_d = std::map<netP *, double>;
using map_P_Vd = std::map<netP *, V_d>;
using map_F_P_Vd = std::map<netF *, map_P_Vd>;
using map_P_P_Vd = std::map<netP *, map_P_Vd>;
using map_P_F_P_Vd = std::map<netP *, map_F_P_Vd>;

////////////////////////////////////////////////////////////////////////////////
// template <typename T>
// std::vector<T *> takeB(const std::vector<T *> &nets)
// {
//   std::vector<T *> ret({});
//   for (const auto &n : nets)
//     if (n->isB())
//       ret.emplace_back(n);
//   return DeleteDuplicates(ret);
// };
// template <typename T>
// std::vector<T *> takeC(const std::vector<T *> &nets)
// {
//   std::vector<T *> ret({});
//   for (const auto &n : nets)
//     if (n->isC())
//       ret.emplace_back(n);
//   return DeleteDuplicates(ret);
// };
// template <typename T>
// std::vector<T *> takeN(const std::vector<T *> &nets)
// {
//   std::vector<T *> ret({});
//   for (const auto &n : nets)
//     if (n->isN())
//       ret.emplace_back(n);
//   return DeleteDuplicates(ret);
// };
// template <typename T>
// std::vector<T *> takeD(const std::vector<T *> &nets)
// {
//   std::vector<T *> ret({});
//   for (const auto &n : nets)
//     if (n->isD())
//       ret.emplace_back(n);
//   return DeleteDuplicates(ret);
// };
// V_Netp takeNetworks(const V_netPp &ps)
// {
//   V_Netp ret({});
//   for (auto &p : ps)
//     network::add(ret, p->getNetwork());
//   return ret;
// };
///////

// class networkSampler
// {
// public:
// 	V_Netp nets_for_C_phin;
// 	V_Netp nets_for_N_phin;
// 	V_Netp nets_for_D_phin;
// 	V_Netp nets_for_C_phi;
// 	V_Netp nets_for_N_phi;
// 	V_Netp nets_for_D_phi;
// 	V_Netp nets_for_C_nabla;
// 	V_Netp nets_for_N_nabla;
// 	V_Netp nets_for_D_nabla;

// 	networkSampler(const V_Netp &all_networks)
// 	{
// 		initialize(all_networks);
// 	};

// 	void initialize(const V_Netp &all_networks)
// 	{
// 		nets_for_C_phin = takeN(all_networks);
// 		nets_for_D_phin = takeD(all_networks);
// 		nets_for_N_phin = Join(takeC(all_networks), takeN(all_networks));

// 		nets_for_C_phi = all_networks;
// 		nets_for_D_phi = Join(takeC(all_networks), takeD(all_networks));
// 		nets_for_N_phi = Join(takeC(all_networks), takeN(all_networks));

// 		nets_for_C_nabla = Join(takeC(all_networks), takeN(all_networks));
// 		nets_for_D_nabla = takeD(all_networks);
// 		nets_for_N_nabla = Join(takeC(all_networks), takeN(all_networks));
// 	};

// 	V_Netp net4nabla(const netPp p)
// 	{
// 		if (p->isC())
// 			return nets_for_C_nabla;
// 		else if (p->isD())
// 			return nets_for_D_nabla;
// 		else
// 			return nets_for_N_nabla;
// 	};

// 	V_Netp net4phin(const netPp p)
// 	{
// 		if (p->isC())
// 			return nets_for_C_phin;
// 		else if (p->isD())
// 			return nets_for_D_phin;
// 		else
// 			return nets_for_N_phin;
// 	};

// 	V_Netp net4phi(const netPp p)
// 	{
// 		if (p->isC())
// 			return nets_for_C_phi;
// 		else if (p->isD())
// 			return nets_for_D_phi;
// 		else
// 			return nets_for_N_phi;
// 	};
// };

// V_Netp getReferenceNetwork(const netPp p, const V_Netp &all_networks)
// {
//   V_Netp nets;
//   if (p->isC())
//     nets = takeN(all_networks);
//   else if (p->isD())
//     nets = Join(takeC(all_networks), takeD(all_networks));
//   else
//     nets = Join(takeC(all_networks), takeN(all_networks));
//   return nets;
// };

// V_Netp getReferenceNetwork(const netPp p, const V_netPp &ps)
// {
//   return getReferenceNetwork(p, takeNetworks(ps));
// };

// 	double SolidAngleDN(const netPp p)
// 	{
// #define use_solidangle
// // #define do_not_use_solidangle
// #ifdef use_solidangle
// 		if (isFlat(p))
// 			return 2. * M_PI;
// 		auto s = p->getSolidAngle();
// 		if (p->isD())
// 		{
// 			return 2. * M_PI;
// 			// return (Between(s, {0.5*2.*M_PI, 1.5*2.*M_PI})) ? s : 2. * M_PI;
// 		}
// 		else if (p->isC())
// 			return (Between(s, {0.3 * M_PI, 1.7 * M_PI})) ? s : M_PI;
// 		else
// 			return (Between(s, {M_PI / 4., 2. * M_PI})) ? s : 2. * M_PI;
// #elif defined(do_not_use_solidangle)
// 		if (isFlat(p))
// 			return 2. * M_PI;

// 		if (p->isD())
// 			return 2. * M_PI;
// 		else if (p->isC())
// 			return M_PI;
// 		else
// 			return 2. * M_PI;
// #else
// 		// if (p->isD())
// 		//   return 2. * M_PI;
// 		// else if (p->isC())
// 		//   return M_PI;
// 		// else
// 		//   return 2. * M_PI;

// 		if (isFlat(p))
// 			return 2. * M_PI;

// 		auto s = p->getSolidAngle();
// 		if (p->isD())
// 			return 2. * M_PI;
// 		else if (p->isC())
// 			return M_PI;
// 		else
// 			return (Between(s, {0.7 * M_PI, 1.3 * M_PI})) ? s : 2. * M_PI;
// #endif
// 	};
////////////////////////////////////////////////////
////////////////////////////////////////////////////

// #include "calcDerivatives.hpp"
// 	/////////////////////////////////////////////////////////
// 	/*calc_phiphin_detail
//   毎回のルンゲクッタの初期値となるphiphinを`P_phiphin`として返す．
//   その値は`known_P_phiphin`を元に作成する．
// calc_phiphin_detail*/
// 	/*calc_phiphin_code*/
// 	class calc_phiphin
// 	{
// 	public:
// 		double dt;
// 		calc_phiphin() : dt(0){};
// 		~calc_phiphin(){};
// 		//-------------------------------------
// 		//シミュレーションの初回はこれ
// 		void init(map_P_Vd &P_phiphin, //毎回のルンゲクッタの初期値となるphiphin
// 				  map_P_Vd &known_P_phiphin,
// 				  BEM::CompGrid &cpg)
// 		{
// 			Print("init");
// 			this->dt = 1E-10;
// 			for (const auto &p : cpg.InnerOuterP)
// 			{
// 				auto it = known_P_phiphin.find(p);
// 				if (it != known_P_phiphin.end())
// 					P_phiphin[p] = it->second;
// 				else
// 					P_phiphin[p] = {0., 0.};
// 			}
// 		};
// 		//-------------------------------------
// 		//シミュレーションの初回以外のRKの初めはこれ
// 		void noninit(map_P_Vd &P_phiphin, //毎回のルンゲクッタの初期値となるphiphin
// 					 map_P_Vd &known_P_phiphin,
// 					 BEM::CompGrid &cpg)
// 		{
// 			try
// 			{
// 				Print("noninit");
// 				this->dt = 0.1;
// 				/*角点のnはNeumann境界条件の法線ベクトルとなるよう修正*/
// 				/*内部点の代入*/

// 				// V_Netp all_networks = takeNetworks(TakeFirst(known_P_phiphin));

// 				networkSampler netsamp(takeNetworks(TakeFirst(known_P_phiphin)));

// 				for (const auto &p : cpg.InnerP)
// 				{ //OuterPのphiphinは内部の値から見積もられるので時間発展はいらない
// 					if (/*nets.isD(p)*/ p->isD())
// 					{ /*Dirichlet*/
// 						// Print(__LINE__,Blue);
// 						auto it = known_P_phiphin.find(p);
// 						if (it != known_P_phiphin.end())
// 							P_phiphin[p] = it->second; //{(it->second)[0], 1E+30};
// 						else
// 						{
// 							mk_vtu("./vtu/new_point.vtu", {{p}});
// 							throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Dirichlet条件の水粒子が内部に新たに入ってくることはありえない");
// 							// Print(__LINE__,_Blue);
// 							// P_phiphin[p] = {BEM::InterpolationRBF_phi(p, known_P_phiphin, netsamp.net4phi(p)), 1E+30};
// 							//P_phiphin[p] = BEM::IDW_phiphin(p, known_P_phiphin, all_networks,3);
// 						}
// 					}
// 					else if (/*nets.isN(p)*/ p->isN())
// 					{ /*Neumann*/
// 						// Print(__LINE__,Red);
// 						P_phiphin[p] = {1E+30, p->phi_n()};
// 					}
// 				}
// 				/*外部点の代入*/
// 				for (const auto &p : cpg.OuterP)
// 				{
// 					if (/*nets.isD(p)*/ p->isD())
// 					{ /*Dirichlet*/
// 						P_phiphin[p] = {0., 0.};
// 					}
// 					else if (/*nets.isN(p)*/ p->isN())
// 					{ /*Neumann*/
// 						P_phiphin[p] = {0., p->phi_n()};
// 					}
// 				}
// 				Print("noninit done");
// 			}
// 			catch (const error_message &e)
// 			{
// 				e.print();
// 				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 			};
// 		};
// 	};
/*calc_phiphin_code*/
/////////////////////////////////////////////////////////

}  // namespace BEM

#endif