#ifndef Network_H
#define Network_H
#pragma once

#define BEM // BEMのメンバー変数関数を有効化する

#include "fundamental.hpp"
#include <unordered_set>

// #include <unordered_map>
#include <functional>
#include <typeinfo>

#include "NetworkCommon.hpp"
#include "object3D.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
using V_i = std::vector<int>;
using VV_i = std::vector<std::vector<int>>;

class networkLine;
class networkFace;
class networkPoint;
class Network;
class networkPointX;

using netP = networkPoint;
using netF = networkFace;
using netPp = networkPoint *;
using netFp = networkFace *;
using V_netLp = std::vector<networkLine *>;
using V_netPp = std::vector<networkPoint *>;
using VV_netPp = std::vector<std::vector<networkPoint *>>;
using V_netFp = std::vector<networkFace *>;
using VV_netFp = std::vector<std::vector<networkFace *>>;
using V_Netp = std::vector<Network *>;

/* ------------------------------------------------------ */
template <class T>
class searcher;

/*extractFaces_detail
`getFaces()`を持つクラスのポインターベクトルを引数にとる．
extractFaces_detail*/
template <class T>
V_netFp extractFaces(const std::vector<T *> &Ls)
{
	V_netFp ret;
	for (const auto &l : Ls)
		for (const auto &f : l->getFaces())
			ret.emplace_back(f);
	return DeleteDuplicates(ret);
};
/*extractPoints_detail
`getPoint()`を持つクラスのポインターベクトルを引数にとる．
extractPoints_detail*/
template <class T>
V_netPp extractPoints(const std::vector<T *> &Ls)
{
	V_netPp ret;
	for (const auto &l : Ls)
		for (const auto &f : l->getPoints())
			ret.emplace_back(f);
	return DeleteDuplicates(ret);
};
template <class T>
V_netLp extractLines(const std::vector<T *> &points)
{
	std::unordered_set<networkLine *> ret;
	for (const auto &p : points)
		for (const auto &l : p->getLines())
			ret.emplace(l);
	return ToVector(ret);
};
/* ------------------------------------------------------ */
/*     *     */
/*    / \    */
/*   *===*   */
/*    \ /    */
/*     *     */
/*networkLine_detail

networkLine_detail*/
/*networkLine_code*/
class networkLine : public object3D
{
#ifdef DEM
public:
	double tension;
#endif

	using S_netPp = searcher<networkPoint> *;
	using S_netFp = searcher<networkFace> *;
	// private:
	// bool intxn;

public:
	bool CORNER;
	bool Dirichlet;
	bool Neumann;

	bool isIntxn();
	// void setIntxn(bool TorF) { this->intxn = TorF; };

#if defined(BEM)
	V_d interpoltedX;
#endif

private:
	std::map<S_netPp, bool> mapPointSearcherStatus;
	std::map<S_netFp, bool> mapFaceSearcherStatus;

public:
	bool getStatus(S_netPp s) const
	{
		auto it = this->mapPointSearcherStatus.find(s);
		if (it != this->mapPointSearcherStatus.end())
			return it->second;
		else
			return false;
	};
	void setStatus(S_netPp s, const bool TorF) { this->mapPointSearcherStatus[s] = TorF; };
	bool isStatus(S_netPp s, const bool TorF) { return this->mapPointSearcherStatus[s] == TorF; };

	bool getStatus(S_netFp s) { return this->mapFaceSearcherStatus[s]; };
	void setStatus(S_netFp s, const bool TorF) { this->mapFaceSearcherStatus[s] = TorF; };
	bool isStatus(S_netFp s, const bool TorF) { return this->mapFaceSearcherStatus[s] == TorF; };

	bool doUKM(const S_netPp s) const { return this->mapPointSearcherStatus.find(s) != this->mapPointSearcherStatus.end() ? true : false; };
	bool doUKM(const S_netFp s) const { return this->mapFaceSearcherStatus.find(s) != this->mapFaceSearcherStatus.end() ? true : false; };

	template <class T>
	void remember(searcher<T> *s) { this->setStatus(s, true); };
	bool forget(S_netPp s) { return network::erase(this->mapPointSearcherStatus, s); };
	bool forget(S_netFp s) { return network::erase(this->mapFaceSearcherStatus, s); };

	//---------------------------------
	void clearSearcherStatus()
	{
		this->mapPointSearcherStatus.clear();
		this->mapFaceSearcherStatus.clear();
	};
	//---------------------------------
	Network *network;
	Network *getNetwork() const { return this->network; };

	bool status;

	bool getStatus() const { return this->status; };
	void setStatus(const bool TorF) { this->status = TorF; };
	bool isStatus(const bool TorF) const { return this->status == TorF; };
	//
	V_netPp XPoints;
	V_netPp getXPoints() const { return this->XPoints; };
	void clearXPoints() { return this->XPoints.clear(); };
	bool clearXPoints(const netP *p) { return network::erase(this->XPoints, p); };
	bool addXPoint(netP *p) { return network::add(this->XPoints, p); };
	bool penetrateQ() const { return !XPoints.empty(); };

	V_netFp getFacesPenetrating() const;
	//---------------------------------
protected:
	// V_netPp Points;
	networkPoint *Point_A;
	networkPoint *Point_B;
	V_netFp Faces;
	//---------------------------------
public:
	networkLine(Network *network_IN, netP *sPoint_IN, netP *ePoint_IN);
	// コピーコンストラクタ
	networkLine(const networkLine *l) : object3D(l->getBounds()), status(false), Point_A(nullptr), Point_B(nullptr)
	{
		std::cout << "copying networkLine ...";
		this->network = l->getNetwork();
		this->Faces = l->getFaces();

		this->Dirichlet = l->Dirichlet;
		this->Neumann = l->Neumann;
		this->CORNER = l->CORNER;

		auto tmpL = l->getPoints();
		set(tmpL[0], tmpL[1]);

		setBounds();
		std::cout << " done" << std::endl;
	};
	~networkLine();
	//---------------------------------
	bool setBounds();
	void setBoundsSingle();
	//---------------------------------
	void set(networkPoint *const sPoint_IN,
			 networkPoint *const ePoint_IN)
	{
		this->Point_A = sPoint_IN;
		this->Point_B = ePoint_IN;
	};
	void set(networkFace *const Face0_IN,
			 networkFace *const Face1_IN)
	{
		this->Faces = {Face0_IN, Face1_IN};
	};
	void set(networkFace *const Face0_IN)
	{
		this->Faces = {Face0_IN};
	};
	//---------------------------------
	netPp divide(const V_d &midX = {});
	bool flip();
	bool flipIfIllegal();
	bool flipIfBetter(const double);
	void divideIfIllegal();
	bool isFlat(const double) const;
	bool islegal() const;
	bool isMergeable() const;
	netPp merge();			  // deleteしていない方のpointを返す
	netPp mergeIfMergeable(); // deleteしていない方のpointを返す
	//---------------------------------
	//  netP* operator()(netP* a){return this->Point2Point[a];};
	template <class T>
	T *getTheOther(const std::vector<T *> &PorF, const T *const a) const
	{
		if (PorF.size() < 2)
			return nullptr;
		else if (a == PorF[0])
			return PorF[1];
		else if (a == PorF[1])
			return PorF[0];
		else
			return nullptr;
	};

	netP *operator()(const netP *const a) const
	{
		if (this->Point_A == a)
			return this->Point_B;
		else if (this->Point_B == a)
			return this->Point_A;
		else
			return nullptr;
	};

	netF *operator()(const netF *const a) const
	{
		if (this->Faces.size() < 2)
			return nullptr;
		else if (a == this->Faces[0])
			return this->Faces[1];
		else if (a == this->Faces[1])
			return this->Faces[0];
		else
			return nullptr;
	};

	VV_d getLocations() const { return obj3D::extractX(this->getPoints()); };
	T2Tddd getLocationsTuple() const;

	V_d getNormal() const;
	Tddd getNormalTuple() const;
	//---------------------------------
	double length() const;
	//------------------
	V_netPp getPoints() const { return {this->Point_A, this->Point_B}; };
	std::tuple<networkPoint *, networkPoint *> getPointsTuple(const networkPoint *const p) const
	{
		if (p == this->Point_A)
			return {this->Point_A, this->Point_B};
		else
			return {this->Point_B, this->Point_A};
	};
	std::tuple<networkPoint *, networkPoint *> getPointsTuple() const { return {this->Point_A, this->Point_B}; };
	const V_netFp &getFaces() const { return this->Faces; };
	V_netFp getFacesExcept(const netFp f) const { return TakeExcept(this->Faces, f); };
	V_netFp getFacesExcept(const V_netFp &fs) const { return TakeExcept(this->Faces, fs); };
	V_netFp getFaces(const bool TorF) const;
	//------------------
	// 派生クラスのクラス名で，選ばれる関数
	int Find(netP *p_IN) const { return network::find(this->getPoints(), p_IN); };
	bool Switch(const netP *oldP, netP *newP)
	{
		if (this->Point_A == oldP)
		{
			this->Point_A = newP;
			return true;
		}
		else if (this->Point_B == oldP)
		{
			this->Point_B = newP;
			return true;
		}
		else
			return false;
	};
	bool Erase(const netP *p_IN) { return Switch(p_IN, nullptr); };
	bool Replace(netP *oldP, netP *newP);
	//------------------
	int Find(netF *f_IN) const { return network::find(this->Faces, f_IN); };
	bool Erase(netF *const f_IN) { return network::erase(this->Faces, f_IN); };
	bool Add(netF *f_IN) { return network::add(this->Faces, f_IN); };
	bool Switch(netF *oldF, netF *newF) { return network::myswitch(this->Faces, oldF, newF); };
	bool Replace(netF *oldF, netF *newF, networkLine *newL = nullptr);
	//------------------
	void Delete();
};

using netL = networkLine;
using netLp = networkLine *;
using V_netLp = std::vector<networkLine *>;
using VV_netLp = std::vector<std::vector<networkLine *>>;

/*networkLine_code*/
/////////////////////////////////////////////////////////////////////
/*networkObject_detail
`networkObject`は，networkPointとnetworkFaceの共通部分をまとめたクラスであり，そのためテンプレートクラスとしている．
例えば，`networkPoint`クラスは，networkObject<networkPoint>をベースとして派生させ作成されている．
この方法で，`networkObject`内のテンプレート関数または変数を派生クラスの関数または変数として使うことができる．
networkObject_detail*/
/*networkObject_code*/
/*  @---@    */
/*   \ / \   */
/*    @---@  */
template <class T>
class networkObject
{
// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM
#ifdef BEM
public:
	bool isB() const;
	bool isC() const;
	bool isD() const;
	bool isN() const;
#endif
	// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM

private:
	std::vector<searcher<T> *> searchers;

protected:
	bool status;
	Network *storage;
	V_netLp Lines;

public:
	void sortLinesByLength()
	{
		std::sort(this->Lines.begin(), this->Lines.end(),
				  [](const netLp l0, const netLp l1)
				  {
					  return (l0->length() < l1->length());
				  });
	};
	//-----------v seracher ------------
	bool doUKM(searcher<T> *s) const
	{
		if (this->searchers.empty())
			return false;
		return (std::find(searchers.cbegin(), searchers.cend(), s) != this->searchers.end()) ? true : false;
	};
	bool remember(searcher<T> *s) { return network::add(this->searchers, s); };
	bool forget(searcher<T> *s) { return network::erase(this->searchers, s); };
	void clearSearchers() { this->searchers.clear(); };
	//---------------------------------
	//------------ status -------------
	bool getStatus() const { return this->status; };
	void setStatus(bool TorF) { this->status = TorF; };
	bool isStatus(bool TorF) const { return this->status == TorF; };
	/* ------------------------------------------------------ */
	Network *network;
	Network *getNetwork() const { return this->network; };

	Network *getStorage() const { return this->storage; };

	networkObject(networkObject *obj) : network(obj->network), status(obj->status), Lines(obj->Lines)
	{
		this->Lines.reserve(100);
	};
	// networkObject(const networkObject &obj) : network(obj.network), status(obj.status), Lines(obj.Lines){};
	networkObject(Network *network_IN, const V_netLp &Lines_IN = {}) : network(network_IN), status(false), Lines(Lines_IN)
	{
		this->Lines.reserve(100);
	};
	virtual ~networkObject(){};
	/* ------------------------------------------------------ */
	int Find(netL *l_IN) { return network::find(this->Lines, l_IN); };
	bool Erase(netL *l_IN) { return network::erase(this->Lines, l_IN); };
	bool Add(netL *l_IN)
	{
		if (!MemberQ(this->Lines, l_IN))
		{
			return this->Lines.emplace_back(l_IN);
			return true;
		}
		else
			return false;
	};
	bool Switch(const netL *oldL, netL *newL)
	{
		for (auto &l : this->Lines)
			if (l == oldL)
			{
				l = newL;
				return true;
			}
		return false;
	};

	/*getLines_detail
getLines_detail*/
	/*getLines_code*/
	/* ------------------------------------------------------ */
	void setLines(const V_netLp &ls) { this->Lines = ls; };
	const V_netLp &getLines() const { return this->Lines; };
	V_netLp getLines_toggle(const bool TorF) const
	{
		V_netLp ret(0);
		for (const auto &line : this->Lines)
			if (line->getStatus() == TorF)
			{
				line->setStatus(!TorF);
				ret.emplace_back(line);
			}
		return ret;
	};
	/*getLines_code*/

	/*getNeighbors_detail
`getNeighbors`は，テンプレート関数になっており，また引数として派生クラス自身を指すポインタ`this`を与えるようにしてある．
これは，Neighborsを`networkLine`を使って取得する際に，**派生クラスの型**によって`networkLine`の動作が異なるためである．
`getNeighbors`は，**派生クラスの型**の情報を必ず求めるため，**派生クラス**においてオーバーライドしなければならない．
getNeighbors_detail*/

	/*getNeighbors_code*/
	std::vector<T *> getNeighbors(const T *obj) const
	{
		std::vector<T *> ret(this->Lines.size());
		int i = 0;
		for (const auto &l : this->Lines)
			ret[i++] = (*l)(obj);
		return ret;
	};
	/*getNeighbors_code*/

	netL *getLineBetween(const T *obj, const T *p) const
	{
		for (const auto &l : obj->Lines)
			if (p == (*l)(obj))
				return l;
		return nullptr;
	};
	//
	int Find(const T *obj, const T *lookforobj) const
	{
		for (int i = 0; i < this->Lines.size(); i++)
			if ((*(this->Lines[i]))(obj) == lookforobj)
				return i;
		return -1;
	};

	virtual V_d getNormal() const = 0;
};
/*networkObject_code*/
/* ------------------------------------------------------ */
// networkPoint限定のもの
template <>
struct Buckets<networkPoint> : public BaseBuckets<networkPoint>
{
	Buckets<networkPoint>(const geometry::CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<networkPoint>(c_bounds, dL_IN){};
	Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<networkPoint>(boundingboxIN, dL_IN){};
	void add(const V_netPp &ps);
	void add(const std::unordered_set<networkPoint *> &ps);
	void add(const Tddd &x, networkPoint *const p)
	{
		auto [i, j, k] = this->indices(x);
		if (isInside(i, j, k))
		{
			// std::cout << i << ", " << j << ", " << k << ", isInside(i, j, k) = " << isInside(i, j, k) << std::endl;
			// std::cout << "this->buckets = " << this->buckets.size() << std::endl;
			// std::cout << "this->buckets[i] = " << this->buckets[i].size() << std::endl;
			// std::cout << "this->buckets[i][j] = " << this->buckets[i][j].size() << std::endl;
			// std::cout << "this->buckets[i][j][k] = " << this->buckets[i][j][k].size() << std::endl;
			this->buckets[i][j][k].emplace(p);
		}
	};
};
/* ------------------------------------------------------ */
template <>
struct Buckets<networkFace> : BaseBuckets<networkFace>
{
	Buckets<networkFace>(const geometry::CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<networkFace>(c_bounds, dL_IN){};
	Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<networkFace>(boundingboxIN, dL_IN){};
};
/* ------------------------------------------------------ */
/*networkPoint_detail
networkPoint_detail*/
/*networkPoint_code*/
/*    \ /    */
/*   --@--   */
/*    / \    */
class networkPoint : public networkObject<networkPoint>, public object3D
{
public:
	Tddd initialX; //必ず設定される
	//
	T6d force;
	T6d inertia;
	T6d velocity;
	T6d acceleration;
	double mass;
	/* ------------------------------------------------------ */
	T6d &F = this->force;
	T6d &I = this->inertia;
	T6d &A = this->acceleration;
	T6d &V = this->velocity;
	/* ------------------------------------------------------ */
	double density;
	double volume;
	double radius;
	/* ------------------------------------------------------ */
	Tddd Fxyz() const { return {std::get<0>(this->force), std::get<1>(this->force), std::get<2>(this->force)}; };
	Tddd Txyz() const { return {std::get<3>(this->force), std::get<4>(this->force), std::get<5>(this->force)}; };
	Tddd Vxyz() const { return {std::get<0>(this->velocity), std::get<1>(this->velocity), std::get<2>(this->velocity)}; };
	Tddd Wxyz() const { return {std::get<3>(this->velocity), std::get<4>(this->velocity), std::get<5>(this->velocity)}; };
		/* ------------------------------------------------------ */
#ifdef DEM
public:
	V_netPp contactP;
	V_netPp neighborP; //! for sph
	void setParticle(double volume_IN, double densityIN)
	{
		this->density = densityIN;
		this->volume = volume_IN;
		this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
		this->mass = this->volume * this->density;
	};

	double W;
	/////////////////////////
	//物性
	double mu_SPH;
	//
	Tddd normal_SPH;
	double d_empty_center;
	double pn_SPH;
	bool pn_is_set;
	bool isSurface;
	bool isFreeFalling;
	double radius_SPH;
	double number_density_SPH;
	double density_interpolated_SPH;
	double pressure_SPH, pressure_SPH_;
	double DPDt_SPH;
	double pressure_Tait(const double rho, double C0 = 1466.) const
	{
		// double C0 = 1466.; //[m/s]
		// C0 /= 5;
		// double C0 = Norm(this->U_SPH); //[m/s]
		// double r = 7.15;
		double r = 7.;
		double rho_w = 1000.;
		double B = rho_w * C0 * C0 / r;
		return B * (std::pow(rho / rho_w, r) - 1.);
	};
	double pressure_Tait() { return pressure_Tait(this->density); };
	double getStaticPressure()
	{
		double rho_w = 1000.;
		double g = 9.81;
		return -rho_w * g * getX()[2];
	};
	/////////////////////////
	Tddd lap_U;
	// void setLapU(const V_d &lap_U_IN) { this->lap_U = lap_U_IN; };
	void setLapU(const Tddd &lap_U_IN) { this->lap_U = lap_U_IN; };
	void setDensityVolume(const double &den, const double &v)
	{
		//質量は保存
		this->density = den;
		this->volume = v;
		this->mass = v * den;
		this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
	};
	void setDensity(const double &den)
	{
		//質量は保存
		this->density = den;
		this->volume = this->mass / this->density;
		this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
	};
	void setDensity_ConstantVolume(const double &den)
	{
		//質量は保存
		this->density = den;
		this->mass = this->volume * this->density;
		this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
	};
	void setVolume(const double &v)
	{
		//質量は保存
		this->volume = v;
		this->density = this->mass / this->volume;
		this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
	};
	double div_U;
	Tddd gradP_SPH;
	//////////////////////////
	netFp face_org;
	double a_viscosity;
	Tddd viscosity_term; // nu*laplacian(U)
	Tddd U_SPH;
	Tddd tmp_U_SPH;
	double tmp_density;
	Tddd pre_U_SPH;
	Tddd mu_lap_rho_g_SPH;
	Tddd interpolated_normal_SPH;
	// ダミー粒子としての情報
	/* ------------------- 多段の時間発展スキームのため ------------------- */
	Tddd DUDt_SPH;
	double DrhoDt_SPH;
	//
	int index;
#endif
	// std::tuple<networkFace * /*補間に使った三角形の頂点*/,
	// 		   double /*パラメタt0*/,
	// 		   double /*パラメタt1*/,
	// 		   double /*深さ方向距離*/,
	// 		   double /*粒子間隔*/>
	// 	particlize_info;
	std::tuple<networkFace * /*補間に使った三角形の頂点*/,
			   std::tuple<networkPoint *, networkPoint *, networkPoint *>, /*補間に使った三角形の頂点*/
			   Tdd /*パラメタt0,t1*/,
			   double /*深さ方向距離*/,
			   double /*粒子間隔*/>
		particlize_info;
	Tddd reflect(const Tddd &v) const;
	/* ------------------------ 境界条件 ------------------------ */
public:
	bool CORNER;
	bool Dirichlet;
	bool Neumann;
	//
	bool isCorner() const { return this->CORNER; };
	bool isDirichlet() const { return this->Dirichlet; };
	bool isNeumann() const { return this->Neumann; };
	void setC()
	{
		this->Dirichlet = false;
		this->Neumann = false;
		this->CORNER = true;
	};
	void unsetC()
	{
		this->CORNER = false;
	};
	//
	void setDirichlet()
	{
		this->Dirichlet = true;
		this->Neumann = false;
		this->CORNER = false;
	};
	void unsetDirichlet() { this->Dirichlet = false; };
	//
	void setNeumann()
	{
		this->Dirichlet = false;
		this->Neumann = true;
		this->CORNER = false;
	};
	void unsetNeumann() { this->Neumann = false; };
	//
	void setD() { this->setDirichlet(); };
	void unsetD() { this->unsetDirichlet(); };
	void setN() { this->setNeumann(); };
	void unsetN() { this->unsetNeumann(); };
	//
	bool isCompleteNeumann()
	{
		//近傍が全てノイマンであればノイマンです．
		//導入理由は？
		for (const auto &p : getNeighbors())
			if (!p->Neumann)
				return false;
		return true;
	};

#ifdef BEM
public:
	double phi_n();
	Tdd phiphin;
	std::unordered_map<networkFace *, Tdd> multiple_phiphin;
	V_d nabla_phi() const;
	V_d getNormalFromSameBC();

	Tddd grad_phi_BEM;
	Tddd U_BEM;
	Tddd U_update_BEM; // EMTによってシフトさせた流速，法線方向成分はU_BEMと一致させる必要がある．
	Tddd U_mod_BEM;
	Tddd U_tangential_BEM;
	Tddd U_normal_BEM;
	Tddd laplacian_U_BEM;
	Tddd normal_BEM;
	Tddd pressure_BEM;
	double kappa_BEM;
	// double DphiDt;
	//
	Tddd normalDirichletFace() const;
	Tddd normalNeumannFace() const;
	Tddd normalContanctSurface(const double pw0, const double pw1) const;

	double height(const Tddd &offset = {0., 0., 0.}, const Tddd &g_center = {0., 0., -1E+20}) const
	{
		//重力中心から，この点までの方向を高さ方向として，offsetから測った高さを返す
		return Dot(Normalize(this->getXtuple() - g_center),
				   this->getXtuple() - offset);
	};

	double DphiDt(const double pressure /*zero if on atmosfere*/,
				  const Tddd &offset = {0., 0., 0.},
				  const Tddd &g_center = {0., 0., -1E+20}) const
	{
		// offsetは平均の水面位置
		double gamma = 72.75 * 1E-3; //[N/m] 水20度
		double gravity = 9.80665;	 // [m/s2]
		double density = 997.;		 // [kg/m3]
		double dphidt = -Dot(this->U_BEM, this->U_BEM) / 2.;
		dphidt -= gamma / density * this->kappa_BEM;
		dphidt -= gravity * height(offset, g_center);
		dphidt -= pressure / density;
		return dphidt + Dot(this->U_BEM, this->U_BEM);
	};

	double DphiDt_EMT(const Tddd &U_modified,
					  const double pressure /*zero if on atmosfere*/,
					  const Tddd &offset = {0., 0., 0.},
					  const Tddd &g_center = {0., 0., -1E+20}) const
	{
		// [1] C. Wang, B. C. Khoo, and K. S. Yeo, “Elastic mesh technique for 3D BIM simulation with an application to underwater explosion bubble dynamics,” Comput. Fluids, vol. 32, no. 9, pp. 1195–1212, Nov. 2003. equaiton (11)
		return this->DphiDt(pressure, offset, g_center) + Dot(U_modified - this->U_BEM, this->U_BEM);
	};

	// double phi_t(const double pressure = 0. /*zero if on atmosfere*/,
	// 			 const Tddd &offset = {0., 0., 0.},
	// 			 const Tddd &g_center = {0., 0., -1E+20}) const
	// {
	// 	return this->DphiDt(pressure, offset, g_center) - Dot(this->U_BEM, this->U_BEM);
	// };

	double pressure(const double DphiDt_IN,
					const Tddd &offset = {0., 0., 0.},
					const Tddd &g_center = {0., 0., -1E+20}) const
	{
		// offsetは平均の水面位置
		double gravity = 9.80665; // [m/s2]
		double density = 997.;	  // [kg/m3]
		double tmp = Dot(this->U_BEM, this->U_BEM) / 2. - gravity * height(offset, g_center) - DphiDt_IN;
		return density * tmp;
	};

	double pressure_EMT(const Tddd &U_modified,
						const double DphiDt_IN,
						const Tddd &offset = {0., 0., 0.},
						const Tddd &g_center = {0., 0., -1E+20}) const
	{
		// offsetは平均の水面位置
		double gravity = 9.80665; // [m/s2]
		double density = 997.;	  // [kg/m3]
		double tmp = Dot(this->U_BEM, this->U_BEM) / 2. - gravity * height(offset, g_center) - DphiDt_IN + Dot(this->U_BEM, U_modified);
		return density * tmp;
	};
#endif

private:
	networkLine *xline;
	networkFace *xface;
	// b% ------------------------------------------------------ */
	// b%       面に対する鏡像位置の粒子．一つの面に対して一点決まる　        */
	// b% ------------------------------------------------------ */
	std::unordered_map<networkFace *, networkPoint *> map_Face_MirrorPoint;

public:
	void makeMirroredPoints(const Buckets<networkFace> &B_face, const double mirroring_distance);
	void clearMirroredPoints()
	{
		for (const auto &[f, p] : this->map_Face_MirrorPoint)
			delete p;
		this->map_Face_MirrorPoint.clear();
	};
	//% ------------------------------------------------------ */
	//%                      接触の判別用　                      */
	//% ------------------------------------------------------ */
private:
	std::unordered_set<networkFace *> ContactFaces;
	std::unordered_set<networkPoint *> ContactPoints;
	std::unordered_map<Network *, std::unordered_set<networkPoint *>> map_Net_ContactPoints;

public:
	//% ------------------------------------------------------ */
	//%                          接触の判別                      */
	//% ------------------------------------------------------ */
	std::unordered_set<networkFace *> getContactFaces() const { return this->ContactFaces; };
	void clearContactFaces() { this->ContactFaces.clear(); };
	void addContactFaces(const Buckets<networkFace> &B, bool); //自身と同じfaceを含まない
	//
	void addContactPoints(const Buckets<networkPoint> &B, const bool);															//自身と同じnetのpointを含み得る
	void addContactPoints(const Buckets<networkPoint> &B, const double radius, const bool);										//自身と同じnetのpointを含み得る
	void addContactPoints(const Buckets<networkPoint> &B, const int limit_depth, const int limit_num, const bool);				//自身と同じnetのpointを含み得る
	void addContactPoints(const std::vector<Buckets<networkPoint>> &B, const double radius, const bool);						//自身と同じnetのpointを含み得る
	void addContactPoints(const std::vector<Buckets<networkPoint>> &B, const int limit_depth, const int limit_num, const bool); //自身と同じnetのpointを含み得る

	//
	std::unordered_set<networkPoint *> getContactPoints() const { return this->ContactPoints; };
	const std::unordered_set<networkPoint *> &getContactPoints(Network *const net) const
	{
		auto it = this->map_Net_ContactPoints.find(net);
		if (it != this->map_Net_ContactPoints.end())
			return it->second;
		else
			return this->map_Net_ContactPoints.at(nullptr);
	};
	std::unordered_set<networkPoint *> getContactPoints(const std::vector<Network *> &net) const
	{
		std::unordered_set<networkPoint *> ret;
		for (const auto &n : net)
		{
			auto tmp = this->getContactPoints(n);
			ret.insert(tmp.begin(), tmp.end());
		}
		return ret;
	};
	void clearContactPoints()
	{
		this->ContactPoints.clear();
		this->map_Net_ContactPoints.clear();
		this->map_Net_ContactPoints[nullptr] = {};
	};
	void addContactPoints(networkPoint *const p)
	{
		this->ContactPoints.emplace(p);
		this->map_Net_ContactPoints[p->getNetwork()].emplace(p);
	};
	void addContactPoints(const std::vector<networkPoint *> &P)
	{
		for (const auto &p : P)
			this->addContactPoints(p);
	};
	void addContactPoints(const std::vector<std::vector<networkPoint *>> &VVP)
	{
		for (const auto &VP : VVP)
			this->addContactPoints(VP);
	};
	void addContactPoints(const std::unordered_set<networkPoint *> &P)
	{
		for (const auto &p : P)
		{
			this->ContactPoints.emplace(p);
			this->map_Net_ContactPoints[p->getNetwork()].emplace(p);
		}
	};
	//@ 制限あり
	void addContactPoints(const Buckets<networkPoint> &B,
						  const Tddd &x,
						  const double radius_SPH,
						  const int limit_number)
	{
		int limit_depth = std::ceil(radius_SPH / B.dL);
		int i0, j0, k0;
		B.indices(x, i0, j0, k0);
		//* ------------------------ depth=0 ------------------------ */
		if (B.isInside(i0, j0, k0))
			this->addContactPoints(B.buckets[i0][j0][k0]);
		if (this->getContactPoints().size() >= limit_number)
			return;
		/* ------------------------------------------------------ */
		int i, j, k;
		auto it_begin = B.buckets.begin();
		for (auto d = 1; d <= limit_depth; ++d)
		{

			if (!(i0 - d < 0 || i0 - d >= B.xsize) && !(i0 + d < 0 || i0 + d >= B.xsize))
			{
				auto &a = B.buckets[i0 + d];
				auto &b = B.buckets[i0 - d];
				for (auto j = (j0 - d < 0 ? 0 : j0 - d); j <= j0 + d && j < B.ysize; ++j)
					for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < B.zsize; ++k)
					{
						this->addContactPoints(a[j][k]);
						this->addContactPoints(b[j][k]);
					}
			}
			else
				for (auto j = (j0 - d < 0 ? 0 : j0 - d); j <= j0 + d && j < B.ysize; ++j)
					for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < B.zsize; ++k)
					{
						if (!((i0 + d) < 0 || (i0 + d) >= B.xsize))
							this->addContactPoints(B.buckets[i0 + d][j][k]);
						if (!((i0 - d) < 0 || (i0 - d) >= B.xsize))
							this->addContactPoints(B.buckets[i0 - d][j][k]);
					}

			auto it_start = std::next(it_begin, i0 - d + 1 < 0 ? 0 : i0 - d + 1);
			auto it_end = std::next(it_begin, i0 + d < B.xsize ? i0 + d : B.xsize);
			if (!(j0 + d < 0 || j0 + d >= B.ysize) && !(j0 - d < 0 || j0 - d >= B.ysize) &&
				!(k0 + d < 0 || k0 + d >= B.zsize) && !(k0 - d < 0 || k0 - d >= B.zsize))
				for (auto it = it_start; it != it_end; ++it)
				{
					for (auto k = k0 - d; k <= k0 + d; ++k)
					{
						this->addContactPoints((*it)[j0 + d][k]);
						this->addContactPoints((*it)[j0 - d][k]);
					}
					for (auto j = j0 - d + 1; j < j0 + d; ++j)
					{
						this->addContactPoints((*it)[j][k0 + d]);
						this->addContactPoints((*it)[j][k0 - d]);
					}
				}
			else
				for (auto it = it_start; it != it_end; ++it)
				{
					if (!(j0 + d < 0 || j0 + d >= B.ysize) && !(j0 - d < 0 || j0 - d >= B.ysize))
						for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < B.zsize; ++k)
						{
							this->addContactPoints((*it)[j0 + d][k]);
							this->addContactPoints((*it)[j0 - d][k]);
						}
					else if (!(j0 + d < 0 || j0 + d >= B.ysize))
						for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < B.zsize; ++k)
						{
							this->addContactPoints((*it)[j0 + d][k]);
						}
					else if (!(j0 - d < 0 || j0 - d >= B.ysize))
						for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < B.zsize; ++k)
						{
							this->addContactPoints((*it)[j0 - d][k]);
						}
					/* ------------------------------------------------------ */
					if (!(k0 + d < 0 || k0 + d >= B.zsize) && !(k0 - d < 0 || k0 - d >= B.zsize))
						for (auto j = (j0 - d + 1 < 0 ? 0 : j0 - d + 1); j < j0 + d && j < B.ysize; ++j)
						{
							this->addContactPoints((*it)[j][k0 + d]);
							this->addContactPoints((*it)[j][k0 - d]);
						}
					else if (!(k0 + d < 0 || k0 + d >= B.zsize))
						for (auto j = (j0 - d + 1 < 0 ? 0 : j0 - d + 1); j < j0 + d && j < B.ysize; ++j)
						{
							this->addContactPoints((*it)[j][k0 + d]);
						}
					else if (!(k0 - d < 0 || k0 - d >= B.zsize))
						for (auto j = (j0 - d + 1 < 0 ? 0 : j0 - d + 1); j < j0 + d && j < B.ysize; ++j)
						{
							this->addContactPoints((*it)[j][k0 - d]);
						}
				}

			if (this->getContactPoints().size() >= limit_number)
				return;
		}
	};

	//@ 無制限に点を取得する
	void addContactPoints(const Buckets<networkPoint> &B,
						  const Tddd &x,
						  const double radius_SPH)
	{
		int limit_depth = std::ceil(radius_SPH / B.dL);
		auto [i0, j0, k0] = B.indices(x);
		if (B.isInside(i0, j0, k0))
		{
			int i_start = i0 - limit_depth;
			int j_start = j0 - limit_depth;
			int k_start = k0 - limit_depth;
			if (i_start < 0)
				i_start = 0;
			if (j_start < 0)
				j_start = 0;
			if (k_start < 0)
				k_start = 0;
			int i_end = i0 + limit_depth;
			int j_end = j0 + limit_depth;
			int k_end = k0 + limit_depth;
			if (i_end > B.xsize)
				i_end = B.xsize;
			if (j_end > B.ysize)
				j_end = B.ysize;
			if (k_end > B.zsize)
				k_end = B.zsize;
			//この方法は，B.Bucketsの長さがそれぞれ一定であることを想定している．
			for (i0 = i_start; i0 < i_end; ++i0)
				for (j0 = j_start; j0 < j_end; ++j0)
					for (k0 = k_start; k0 < k_end; ++k0)
						this->addContactPoints(B.buckets[i0][j0][k0]);
		}
	};
	//% ------------------------------------------------------ */
	// インライン化する．
	// コンタクトポイントをpointとface上で作成，完成させる．
	// 次にRBFを使った力の計算をその点を使って行えるようにする．
	//--------------------
	//コンストラクタ
	// networkPoint(Network *network_IN, const V_d &xyz_IN, networkLine *line_IN = nullptr, networkFace *face_IN = nullptr);
	networkPoint(Network *network_IN,
				 Network *storage_IN,
				 const Tddd &xyz_IN,
				 networkLine *xline_IN = nullptr,
				 networkFace *xface_IN = nullptr);
	~networkPoint();
	/* ------------------------------------------------------ */
	V_d getFaceAreas() const;
	//---------------------
	networkLine *getXLine() const { return this->xline; };
	networkFace *getXFace() const { return this->xface; };
	bool isXPoint() const { return (xline != nullptr && xface != nullptr); };
	bool intersectQ() const { return this->isXPoint(); };
	/*for a cross point*/
	//--------------
	// void set(const V_d& xyz_IN){object3D::setBounds(xyz_IN);};
	void setX(const V_d &xyz_IN);
	void setX(const Tddd &xyz_IN);
	void setXSingle(const Tddd &xyz_IN) { object3D::setBounds(xyz_IN); };
	bool setXcarefully(const V_d &xyz_IN);
	void resetXinfo();
	void setLinesStatus(const bool TorF)
	{
		for (const auto &l : Lines)
			l->setStatus(TorF);
	};
	void setBounds();
	void setBoundsSingle() { object3D::setBounds(geometry::CoordinateBounds(this->X)); };
	V_netPp getXNeighbors() const;
	//--------------
	V_netPp getNeighbors() const { return networkObject::getNeighbors(this); };
	// V_netPp getNeighborsExcept(const netP *p_IN) const { return networkObject::getNeighborsExcept(this, p_IN); };
	int Find(const netP *p_IN) const { return networkObject::Find(this, p_IN); };
	netL *getLineBetween(const netP *p) const { return networkObject::getLineBetween(this, p); };

	bool isOpen() const
	{
		for (auto const &l : this->networkObject::getLines())
			if (l->getFaces().size() < 2)
				return true;
		return false;
	};
	bool isClosed() const
	{
		for (auto const &l : this->networkObject::getLines())
			if (l->getFaces().size() < 2)
				return false;
		return true;
	};

	V_netPp getNeighborsSort() const;
	V_netPp getNeighborsSort2() const;
	V_netPp getNeighborsSort(bool TorF);
	//
	// using V_Tup_ddPp = std::vector<std::tuple<double, double, networkPoint *>>;
	using V_Tup_ddVdPp = std::vector<std::tuple<double, double, V_d, networkPoint *>>;
	V_Tup_ddVdPp getNeighbors_Depth1_OnPolarAsTuple(networkLine *base_line) const;
	V_Tup_ddVdPp getNeighbors_Depth2_OnPolarAsTuple(networkLine *base_line) const;
	V_Tup_ddVdPp getNeighbors_Depth2_OnPolarAsTuple_parametric(networkLine *base_line) const;
	using V_Tup_ddVd = std::vector<std::tuple<double, double, V_d>>;
	V_Tup_ddVd getNeighbors_Depth2_OnPolarAsTuple2() const;
	// V_Tup_dddPp getNeighborsOnPolar() const;

	using V_Var_VVdVPp = std::vector<std::variant<VV_d, V_netPp>>;
	V_Var_VVdVPp getNeighbors_Depth1_OnPolarAsVariant() const;
	V_Var_VVdVPp getNeighbors_Depth2_OnPolarAsVariant() const;
	/* ------------------------------------------------------ */
	V_netFp getFaces_intersectQ(const bool TorF) const;
	/* ------------------------------------------------------ */
	V_netFp getFaces(networkLine *line) const;
	/* ------------------------------------------------------ */
	V_netFp getFaces(const bool TorF) const
	{
		V_netFp ret;
		for (const auto &l : this->networkObject::getLines())
			for (const auto &f : l->getFaces(TorF))
				if (std::find(ret.begin(), ret.end(), f) == ret.end())
					ret.emplace_back(f);
		return ret;
	};

	V_d getAngles() const;
	V_d getAngles(networkLine *const base_line) const;
#define use_binary_search
	// 	V_netFp getFaces() const
	// 	{
	// 		V_netFp ret;
	// 		V_netFp::iterator it;
	// 		for (const auto &l : this->networkObject::getLines())
	// 			for (const auto &f : l->getFaces())
	// 			{
	// #if defined(use_binary_search)
	// 				if (!std::binary_search(ret.begin(), ret.end(), f))
	// 					ret.insert(std::lower_bound(ret.begin(), ret.end(), f), f);
	// #else
	// 				if (std::find(ret.begin(), ret.end(), f) == ret.end())
	// 					ret.emplace_back(f);
	// #endif
	// 			}
	// 		return ret;
	// 	};
	//
	/**
	2021/09/02unordered_setを使うよう修正した
	将来的にはunordered setを返す関数に修正すべき
	 **/
	V_netFp getFaces() const
	{
		V_netFp ret;
		std::unordered_set<networkFace *> tmp;
		for (const auto &l : this->networkObject::getLines())
			for (const auto &f : l->getFaces())
				tmp.insert(f);
		ret.assign(tmp.begin(), tmp.end());
		return ret;
	};
	std::unordered_set<networkFace *> getFacesUO() const
	{
		std::unordered_set<networkFace *> ret;
		for (const auto &l : this->networkObject::getLines())
			for (const auto &f : l->getFaces())
				ret.emplace(f);
		return ret;
	};

	V_netFp getFacesSort(networkLine *const l) const;
	V_netFp getFacesSort() const;
	V_netFp getFacesSort2() const;
	// V_netFp getFacesSort() const;
	//--------------
	// bool Replace(netL *oldL, netL *newL)
	// {
	//   this->Switch(oldL, newL); //1
	//   oldL->Erase(this);      //2
	//   newL->Add(this);        //3
	//   return true;
	// };
	//--------------
	V_d getNormal() const override;
	Tddd getNormalTuple() const;
	Tddd getNormalAreaAveraged() const;
	void Delete();
	//--------------

	/*SolidAngle_detail

	SolidAngle_detail*/
	/*SolidAngle_detail_code*/
	double getSolidAngle();
	double getSolidAngle(bool TorF);
	double getSolidAngle(const V_netFp &faces);
	/*SolidAngle_detail_code*/
};

//@ ------------------------ 抽出用関数 ----------------------- */
//@ --------------------------------------------------------- */
// 注意：unordered_setからデータを取る場合．順番は保障されない
std::vector<Tddd> extX(const V_netPp &ps)
{
	std::vector<Tddd> ret;
	for (const auto &p : ps)
		ret.emplace_back(p->getXtuple());
	return ret;
};
std::vector<Tddd> extX(const std::unordered_set<networkPoint *> &ps)
{
	std::vector<Tddd> ret;
	for (const auto &p : ps)
		ret.emplace_back(p->getXtuple());
	return ret;
};
std::vector<Tddd> extNormals(const V_netPp &ps)
{
	std::vector<Tddd> ret;
	for (const auto &p : ps)
		ret.emplace_back(p->getNormalTuple());
	return ret;
};

#ifdef BEM
double extPhi(const networkPoint *const p)
{
	return std::get<0>(p->phiphin);
};

double extPhin(const networkPoint *const p)
{
	return std::get<1>(p->phiphin);
};

V_d extPhi(const V_netPp &ps)
{
	V_d ret(ps.size());
	int i = 0;
	for (const auto &p : ps)
		ret[i++] = std::get<0>(p->phiphin);
	return ret;
};

V_d extPhin(const V_netPp &ps)
{
	V_d ret(ps.size());
	int i = 0;
	for (const auto &p : ps)
		ret[i++] = std::get<1>(p->phiphin);
	return ret;
};
std::vector<Tdd> extPhiphin(const V_netPp &ps)
{
	std::vector<Tdd> ret(ps.size());
	int i = 0;
	for (const auto &p : ps)
		ret[i++] = p->phiphin;
	return ret;
};

#endif
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
/*networkPoint_code*/
double distance(const networkPoint *p0, const networkPoint *p1) { return Norm(p0->getX() - p1->getX()); };
netL *getLineConnected(const networkPoint *obj, const networkPoint *p)
{
	for (const auto &l : obj->networkObject::getLines())
		if (p == (*l)(obj))
			return l;
	return nullptr;
};
netL *link(netP *obj, netP *obj_, Network *net)
{
	try
	{
		// std::cout << "try to link" << std::endl;
		if (!obj && !obj_)
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "both arguments are nullptr");
		if (!obj)
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "first argument is nullptr");
		else if (!obj_)
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "second argument is nullptr");
		else if (obj_ == obj)
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "a point is trying to link itself!");

		auto line = getLineConnected(obj /*このLineリストにobj_がある場合，そのlineオブジェクトを返す*/, obj_);
		auto line_ = getLineConnected(obj_, obj);

		if (line && line_)
		{
			if (line == line_)
				return line;
			else
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Points are linked but by differrent lines");
		}
		else if (!line && line_)
		{
			obj->Add(line_);
			return line_; // only the other side was linked
		}
		else if (line && !line_)
		{
			obj_->Add(line);
			return line; // only the other side was linked
		}
		else
		{
			// std::cout << "lineはコンストラクタで引数のオブジェクトを自身のオブジェクトリスト(Points,Faces)に保存し，さらにオブジェクトのLineリストに自身thisを保存する" << std::endl;
			return new networkLine(net, obj, obj_);
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
/////////////////////////////////////////////////////////////////////
class pathInfo
{
public:
	V_netFp face;
	VV_d xyz; /* {start_xyz, cwp_xyz, crosspoint_xyz} */
	double incidentAngle;
	double r;
	pathInfo(V_netFp face_IN,
			 const VV_d &xyz_IN,
			 const double incidentAngle_IN,
			 const double r_IN) : face(face_IN), xyz(xyz_IN), incidentAngle(incidentAngle_IN), r(r_IN){};
	~pathInfo()
	{
		//    std::cout << Red << "destructed" << reset << std::endl;
	}
	void info()
	{
		std::cout << Blue << "         face :" << face << std::endl;
		std::cout << Blue << "          xyz :" << xyz << std::endl;
		std::cout << Blue << "incidentAngle :" << incidentAngle << std::endl;
		std::cout << Blue << "            r :" << r << reset << std::endl;
	};
};
/////////////////////////////////////////////////////////////////////

/* ------------------------------------------------------ */
/*networkFace_detail
networkFace_detail*/
/*networkFace_code*/
/*     @     */
/*    / \    */
/*   @-->@   */
class networkFace : public networkObject<networkFace>, public object3D
{
public:
	/* ------------------------------------------------------ */
	bool Dirichlet;
	bool Neumann;
	bool isDirichlet() const { return this->Dirichlet; };
	bool isNeumann() const { return this->Neumann; };
	/* ------------------------------------------------------ */
	T6d force;
	V_d inertia;  //={mass,mass,mass,rotational_of_inertia}
	V_d velocity; //={vx,vy,vz,angular_vx,angular_vy,angular_vz}
	V_d acceleration;
	double mass;
	/* ------------------------------------------------------ */
	T6d &F = this->force;
	V_d &I = this->inertia;
	V_d &A = this->acceleration;
	V_d &V = this->velocity;
	/* ------------------------------------------------------ */
	V_d center_of_mass;
#ifdef BEM
	Tdd phiphin;
#endif

#ifdef DEM
public:
	V_netPp contactP;
#endif

private:
	VV_d xyzInverse;
	Tddd normal;
	V_d angles;
	double area;
	V_netPp XPoints;
	V_netPp Points;
	std::tuple<networkPoint *, networkPoint *, networkPoint *> PointsTuple;
	//% ------------------------------------------------------ */
	std::unordered_set<networkPoint *> ParametricPoints;
	//! ------------------------------------------------------ */
	std::unordered_set<networkPoint *> ContactPoints;
	std::unordered_map<Network *, std::unordered_set<networkPoint *>> map_Net_ContactPoints;

public:
	//% ------------------------------------------------------ */
	//%           パーティクライズで作られるParametricPoints        */
	//% ------------------------------------------------------ */
	const std::unordered_set<networkPoint *> &getParametricPoints() const { return this->ParametricPoints; };
	void clearParametricPoints()
	{
		for (const auto &p : this->ParametricPoints)
			delete p;
		this->ParametricPoints.clear();
	};
	void addParametricPoints(networkPoint *const p) { this->ParametricPoints.emplace(p); };
	void addParametricPoints(const V_netPp &P) { this->ParametricPoints.insert(P.begin(), P.end()); };
	void addParametricPoints(const std::unordered_set<networkPoint *> &P) { this->ParametricPoints.insert(P.begin(), P.end()); };
	//! ------------------------------------------------------ */
	//!                          接触の判別                      */
	//! ------------------------------------------------------ */
	const std::unordered_set<networkPoint *> &getContactPoints() const { return this->ContactPoints; };
	const std::unordered_set<networkPoint *> &getContactPoints(Network *const net) const
	{
		auto it = this->map_Net_ContactPoints.find(net);
		if (it != this->map_Net_ContactPoints.end())
			return it->second;
		else
			return this->map_Net_ContactPoints.at(nullptr);
	};
	const std::unordered_set<networkPoint *> getContactPoints(const std::vector<Network *> &nets) const
	{
		std::unordered_set<networkPoint *> ret;
		for (const auto &n : nets)
		{
			auto points = getContactPoints(n);
			ret.insert(points.begin(), points.end());
		}
		return ret;
	};
	void clearContactPoints()
	{
		this->ContactPoints.clear();
		this->map_Net_ContactPoints.clear();
		this->map_Net_ContactPoints[nullptr] = {};
	};
	void addContactPoints(networkPoint *const p)
	{
		this->ContactPoints.emplace(p);
		this->map_Net_ContactPoints[p->getNetwork()].emplace(p);
	};
	void addContactPoints(const std::vector<networkPoint *> &P)
	{
		for (const auto &p : P)
		{
			this->ContactPoints.emplace(p);
			this->map_Net_ContactPoints[p->getNetwork()].emplace(p);
		}
	};
	void addContactPoints(const std::vector<std::vector<networkPoint *>> &VVP)
	{
		for (const auto &VP : VVP)
			for (const auto &p : VP)
			{
				this->ContactPoints.emplace(p);
				this->map_Net_ContactPoints[p->getNetwork()].emplace(p);
			}
	};
	//! ------------------------------------------------------ */
	void reverseNormal()
	{
		std::reverse(this->Lines.begin(), this->Lines.end());
		this->setBounds();
	};
	/* ------------------------------------------------------ */
	Tddd mirror(const Tddd &v) const
	{
		//与えられたvの鏡像を返す
		return v - 2. * Dot(v, this->normal) * this->normal;
	};
	Tddd mirrorPosition(const Tddd &v) const
	{
		//与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
		auto u = std::get<0>(this->PointsTuple)->getXtuple();
		return v + 2. * Dot(u - v, this->normal) * this->normal;
	};
	Tddd mirrorPosition(const Tddd &v, const double scale) const
	{
		//与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
		return v + scale * Dot(this->getXtuple() - v, this->normal) * this->normal;
	};
	Tddd mirrorPosition(const networkPoint *p, const double scale = 2.) const
	{
		//与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
		return p->getXtuple() + scale * Dot(this->getXtuple() - p->getXtuple(), this->normal) * this->normal;
	};
	//-------------------------
	//  networkFace(networkFace* f):object3D(f->getBounds()),networkObject(f),XPoints(0){};
	/*networkFace_constructor_detail
	##### コンストラクタ
	`networkFace`コンストラクタでは，引数の`Line`ベクトルに，
	コンストラクトしている`networkFace`を`set`するが，
	`networkLine`の`set`関数は，`networkFace`の保存を2つまでしか許さないので，setされない`危険性`がある．
	networkFace_constructor_detail*/
	// networkFace(Network *network_IN, const V_netLp &Lines_IN) : object3D(), networkObject(network_IN, Lines_IN), XPoints(0)
	// {
	//   for (const auto &l : Lines_IN)
	//     l->set(this);
	//   setBounds();
	// };
	networkFace(Network *network_IN, Network *storage_IN, const V_netLp &Lines_IN);
	/*コピーコンストラクタ*/
	networkFace(const netFp f);
	~networkFace();
	/* ------------------------------------------------------ */
	std::unordered_set<networkPoint *> particlize(const double dx, const V_d &depth_list);
	/* ------------------------------------------------------ */
	bool intersectQ() const { return (this->penetratedQ() || this->penetrateQ()); };
	bool intersectQ(const netP *p) const
	{
		if (this == p->getXFace())
			return true;
		if (network::find(this->Lines, p->getXLine()) != -1)
			return true;
		return false;
	};
	bool penetrateQ() const
	{
		for (const auto &l : this->Lines)
			if (l->penetrateQ())
				return true;
		return false;
	};
	bool penetratedQ() const { return (!this->XPoints.empty()); /*貫かれているかどうか*/ };
	bool penetratedQ(netL *l) const
	{ /*lによって貫かれているかどうか*/
		for (const auto &p : (this->XPoints))
			if (p->getXLine() == l)
				return true;
		return false;
	};
	//------------- Cross ---------------
	// XPoints
	V_netPp getXPoints() const;
	bool addXPoint(netP *p) { return network::add(this->XPoints, p); };
	void clearXPoints() { this->XPoints.clear(); };
	bool clearXPoints(netP *p) { return network::erase(this->XPoints, p); };
	V_netPp getPointsPenetrated() const { return this->XPoints; };
	V_netPp getPointsPenetrate() const
	{
		V_netPp XPoints;
		for (const auto &line : this->networkObject::getLines())
			for (const auto &p : line->XPoints)
				XPoints.emplace_back(p);
		return XPoints;
	};
	V_netPp getPointsPenetrated(const Network *net) const
	{
		V_netPp ret({});
		for (const auto &p : this->getPointsPenetrated())
			if (p->getNetwork() == net)
				ret.emplace_back(p);
		return ret;
	};
	V_netPp getPointsPenetrate(const Network *net) const
	{
		V_netPp ret({});
		for (const auto &p : this->getPointsPenetrate())
			if (p->getNetwork() == net)
				ret.emplace_back(p);
		return ret;
	};

	V_netPp getPointsIntersect() const
	{
		V_netPp ret = this->getPointsPenetrate();
		ret.insert(ret.end(), this->XPoints.begin(), this->XPoints.end());
		return ret;
	};
	// Lines
	V_netLp getLinesPenetrate() const
	{
		V_netLp ret;
		for (const auto &l : this->Lines)
			if (l->penetrateQ())
				ret.emplace_back(l);
		return ret;
	};
	V_netLp getLinesPenetrated() const
	{
		V_netLp ret;
		for (const auto &p : (this->XPoints))
			if (p->getXLine()->penetrateQ())
				ret.emplace_back(p->getXLine());
		return ret;
	};

	/*getPointsOnLines_detail
`Intersection`（Mathematica like）で面の頂点を順に取得しながら，
もし線上に`xpoint`があればxpointsを`network::sortByDistance`でソートして同時に取得いく．
結果として，this->面の隣接する線に関係する全ての点をccw周りにソートして返す．
getPointsOnLines_detail*/
	/*getPointsOnLines_code*/
	V_netPp getPointsOnLines() const;
	VV_netPp getPointsCutLines() const;
	VV_netPp getPointsOnLinesDivided() const;
	VV_netPp getPointsCutFaces() const;

	/*getPointsCutFaces_code*/

	// use networkObject functions
	V_netFp getNeighbors() { return networkObject::getNeighbors(this); };
	// V_netFp getNeighborsExcept(networkFace *f_IN) { return networkObject::getNeighborsExcept(this, f_IN); };
	int Find(networkFace *f_IN) { return networkObject::Find(this, f_IN); };
	netL *getLineBetween(networkFace *f) { return networkObject::getLineBetween(this, f); };

private:
	/*
	Pointsから，vectorToXPointとなるPointを探し，線の貫いている面の法線方向と比較する．比較した結果が，正なら与えられたPointsは面の後方に位置するとし，trueを返す．
	注意：`this`ネットワークと異なるネットワークに属する点は，干渉点XPointと判断する．
  */
	bool isXPointsBehind(const V_netPp &p) const
	{
		int s = p.size();
		V_d vectorToXPoint;
		for (auto i = 0; i < s; i++)
		{
			if (p[i]->getNetwork() != this->network /* means p[i] is interaction netwrok*/)
			{
				if (p[(i + 1) % s]->getNetwork() == this->network)
				{
					vectorToXPoint = p[i]->getX() - p[(i + 1) % s]->getX();
					if (Dot(p[i]->getXFace()->getNormal(), vectorToXPoint) > 1E-12)
						return true;
				}
				if (p[(i - 1 + s) % s]->getNetwork() == this->network)
				{
					vectorToXPoint = p[i]->getX() - p[(i - 1 + s) % s]->getX();
					if (Dot(p[i]->getXFace()->getNormal(), p[i]->getX() - p[(i - 1 + s) % s]->getX()) > 1E-12)
						return true;
				}
			}
		}
		//変更20200921
		// for(auto i=0; i<s ;i++){
		//   if(p[i]->isXPoint())
		// 	if((p[i]->xface)->getNetwork() != this->network)
		// 	  if(p[i]->getNetwork() == this->network && p[(i+1)%s]->getNetwork() == this->network)
		// 	    if(Dot(p[i]->xface->getNormal(), p[i]->getX() - p[(i+1)%s]->getX()) > 1E-12)
		// 	      return true;
		//   if(p[(i+1)%s]->isXPoint())
		// 	if((p[(i+1)%s]->xface)->getNetwork() != this->network)
		// 	  if(p[(i+1)%s]->getNetwork() == this->network && p[i]->getNetwork() == this->network)
		// 	    if(Dot(p[(i+1)%s]->xface->getNormal(), p[(i+1)%s]->getX() - p[i]->getX()) > 1E-12)
		// 	      return true;
		//}
		return false;
	};

public:
	VV_netPp getPointsCutFacesBehind() const
	{
		VV_netPp ret;
		for (const auto &pSet : this->getPointsCutFaces())
		{
			if (isXPointsBehind(pSet))
				ret.emplace_back(pSet);
		}
		return ret;
	};
	VV_netPp getPointsCutFacesFront() const
	{
		VV_netPp ret;
		for (const auto &pSet : this->getPointsCutFaces())
			if (!isXPointsBehind(pSet))
				ret.emplace_back(pSet);
		return ret;
	};
	//=================================
	VV_d getXInverse() const
	{
		// this->setBounds();
		return this->xyzInverse;
	};
	V_d parameterize(const V_d &xyz_IN) { return Dot(getXInverse(), xyz_IN); };
	//--------------------------
	netL *longestLineAround()
	{
		netL *ret, *line;
		double len = 0, v;
		for (const auto &l : Lines)
			for (const auto &f : l->getFaces())
			{
				if (f != this)
				{
					line = f->longestLine();
					v = Norm(line->getPoints()[1]->X - line->getPoints()[0]->X);
					if (v > len)
					{
						len = v;
						ret = line;
					}
				}
			}
		return ret;
	};
	netL *longestLine()
	{
		netL *ret;
		double len = 0, v;
		for (const auto &l : Lines)
		{
			v = Norm(l->getPoints()[1]->X - l->getPoints()[0]->X);
			if (v > len)
			{
				len = v;
				ret = l;
			}
		}
		return ret;
	};

	// void setBounds()
	// {
	//   object3D::setBounds(CoordinateBounds(getLocations()));
	//   //this->xyzInverse=Inverse(this->getLocations());

	//   V_netPp Points = getPoints();
	//   V_d a = Points[1]->xyz - Points[0]->xyz;
	//   V_d b = Points[2]->xyz - Points[0]->xyz;
	//   V_d c = Cross(a, b);
	//   this->normal = c / Norm(c);
	//   if (!isfinite(this->normal))
	//     throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "this->normal is not finite"));

	//   this->area = std::abs(Norm(c) / Dot(a, b));

	//   int s = 3;
	//   this->angles.resize(s);
	//   for (int i = 0; i < s; i++)
	//   {
	//     a = Points[(i + s + 1) % s]->xyz - Points[i]->xyz;
	//     b = Points[(i + s - 1) % s]->xyz - Points[i]->xyz;
	//     c = Cross(a, b);
	//     double pORn = sgn(Dot(c, this->normal));
	//     if (pORn < 0)
	//       throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "what"));
	//     this->angles[i] = atan2(Norm(c), Dot(a, b)) * pORn;
	//   }
	// };

	////////////////////////////////////////////
	void setPoints(networkPoint *const p0, networkPoint *const p1, networkPoint *const p2)
	{
		this->Points = {p0, p1, p2};
		this->PointsTuple = {p0, p1, p2};
	};
	bool setBounds()
	{
		//!再計算を行う
		//線の持つ点のインターセクションをチェックすることで，面の持つ点を間接的に取得する．
		//* default alignment
		/*          angle[2]
					/\
				   /  \
			 l[2] /    \ l[1]
				 /______\
		angle[0]  l[0]  angle[1]
		*/
		try
		{
			//!点の再取得

			setPointsFromCurrentLines();
			T3Tddd abc = {std::get<0>(this->PointsTuple)->getXtuple(),
						  std::get<1>(this->PointsTuple)->getXtuple(),
						  std::get<2>(this->PointsTuple)->getXtuple()};
			object3D::setBounds(abc);

			if (!isFinite(this->area = TriangleArea(abc)))
			{
				std::cout << "this->PointsTuple = " << this->Points << std::endl;
				std::cout << "abc = " << abc << std::endl;
				std::cout << Grid({"area", this->area}, 40) << std::endl;
				std::cout << Grid({"angle", this->angles}, 40) << std::endl;
				std::cout << "Points X" << extX(this->Points) << std::endl;
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not finite");
			}
			//
			Tddd normal = TriangleNormal(abc);
			if (!isFinite(normal))
			{
				std::cout << "Lines = " << this->Lines << std::endl;
				std::cout << "線が更新されておらずエラーになる可能性がある．normal = " << normal << std::endl;
				std::cout << "全ての線が更新されている必要があるため" << std::endl;

				/*
				@点の座標変更
					V
				!  PointのsetBounds()
					|
				@	|        線のつなぎ変え
					V             V
				!   LineのsetBounds()
					|             |
					V             V
				!   FaceのsetBounds()
				*/

				std::cout << "Lines[0]->getPoints() = " << Lines[0]->getPoints() << std::endl;
				std::cout << "Lines[1]->getPoints() = " << Lines[1]->getPoints() << std::endl;
				std::cout << "Lines[2]->getPoints() = " << Lines[2]->getPoints() << std::endl;

				std::cout << "Lines[0]->getX() = " << Lines[0]->getXtuple() << std::endl;
				std::cout << "Lines[1]->getX() = " << Lines[1]->getXtuple() << std::endl;
				std::cout << "Lines[2]->getX() = " << Lines[2]->getXtuple() << std::endl;

				std::cout << "Points = " << this->Points << std::endl;
				std::cout << "abc = " << abc << std::endl;
				std::cout << "normal = " << normal << std::endl;
				/*
				TriangleNormalは，内部で，
				Normalize(Cross((b - a), (c - a)));の計算をしている．
				Cross((b - a), (c - a))
				Cross(A, X)
				X = {0,0,0}なら
				return {std::get<1>(A) * std::get<2>(X) - std::get<2>(A) * std::get<1>(X),
						std::get<2>(A) * std::get<0>(X) - std::get<0>(A) * std::get<2>(X),
						std::get<0>(A) * std::get<1>(X) - std::get<1>(A) * std::get<0>(X)};
				*/
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not finite");
			}
			else
			{
				// this->normal = {std::get<0>(normal), std::get<1>(normal), std::get<2>(normal)};
				this->normal = normal;
			}
			//
			Tddd angles = TriangleAngles(abc);
			if (!isFinite(angles))
			{
				std::cout << "angles = " << angles << std::endl;
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "not finite");
			}
			else
				this->angles = {std::get<0>(angles), std::get<1>(angles), std::get<2>(angles)};
			//
			return true;
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	////////////////////////////////////////
	//------------------------
	void Delete();
	//-------------------------
	// bool Replace(netL *oldL, netL *newL, netF *newF = nullptr) { return network::replace(this, this->Lines, oldL, newL, newF); };
	bool Replace(netL *oldL, netL *newL, netF *newF = nullptr)
	{
		// switchでないと，順番に意味のあるFaceではおかしくなるので注意
		if (this->Switch(oldL, newL) && oldL->Erase(this) && newL->Add(this))
		{
			if (newF != nullptr)
			{
				oldL->Add(newF);
				newF->Add(oldL);
			}
			return true;
		}
		return false;
	};
	//-------------------------
	V_d getNormal() const override { return {std::get<0>(normal), std::get<1>(normal), std::get<2>(normal)}; };
	const Tddd &getNormalTuple() const { return this->normal; };
	const double &getArea() const { return this->area; };
	//////////////////////
	// このfaceの保存状況に従って，lに対する前後のpointへのインデックスを取得できる
	Tiii point_indicies(const netL *l) const
	{
		try
		{
			for (auto i = 0; i < 3; i++)
				if ((l->getPoints()[0] == this->Points[i] /*back*/ && l->getPoints()[1] == this->Points[(i + 1) % 3] /*front*/) ||
					(l->getPoints()[1] == this->Points[i] /*back*/ && l->getPoints()[0] == this->Points[(i + 1) % 3] /*front*/))
					return {i, (i + 1) % 3, (i + 2) % 3};
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			std::stringstream ss;
			ss << "このlineを基準としてインデックスをつくれない：この線はこの面のいっぺんではない" << std::endl;
			ss << "setBoundsを忘れていませんか？ 面の線を変更した際などは，setBoundsを忘れないように" << std::endl;
			ss << "FaceのgetPoints()は，lineから毎回間接的に取得することはやめて，setBoundsの際に保存するように変更しました" << std::endl;
			ss << "input line l = " << l << std::endl;
			ss << "l->Points = " << l->getPoints() << std::endl;
			ss << "this face->Points = " << this->Points << std::endl;
			ss << "this face->Lines = " << this->Lines << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
		};
	};
	//与えられたpをindex[0]として，this->Pointsのインデックスを返す
	Tiii point_indicies(const networkPoint *p) const
	{
		try
		{
			if (p == this->Points[0])
				return {0, 1, 2};
			else if (p == this->Points[1])
				return {1, 2, 0};
			else if (p == this->Points[2])
				return {2, 0, 1};
			// for (auto i = 0; i < 3; i++)
			// 	if (p == this->Points[i])
			// 		return {i, (i + 1) % 3, (i + 2) % 3};
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			std::stringstream ss;
			ss << "このpointを基準としてインデックスをつくれない：この線はこの面のいっぺんではない" << std::endl;
			ss << "setBoundsを忘れていませんか？ 面の線を変更した際などは，setBoundsを忘れないように" << std::endl;
			ss << "FaceのgetPoints()は，lineから毎回間接的に取得することはやめて，setBoundsの際に保存するように変更しました" << std::endl;
			ss << "input point p = " << p << std::endl;
			ss << "this face->Points = " << this->Points << std::endl;
			ss << "this face->Lines = " << this->Lines << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
		};
	};

	Tddd getAnglesTuple() const { return {this->angles[0], this->angles[1], this->angles[2]}; };
	V_d getAngles() const { return this->angles; };
	V_d getAngles(const netL *l) const
	{
		try
		{
			if (l == this->Lines[0])
				return this->angles;
			else if (l == this->Lines[1])
				return RotateLeft(this->angles, 1);
			else
				return RotateLeft(this->angles, 2);
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	V_d getAngles(const networkPoint *p) const
	{
		try
		{
			auto [i, j, k] = point_indicies(p);
			return {this->angles[i], /*ここにこの線lが位置する*/ this->angles[j], this->angles[k]};
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	double getAngle(const networkPoint *p) const
	{
		if (p == this->Points[0])
			return this->angles[0];
		else if (p == this->Points[1])
			return this->angles[1];
		else if (p == this->Points[2])
			return this->angles[2];
		else
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point is not found");
	};
	/* ------------------------------------------------------ */
	// V_netPp getPoints(const netPp origin = nullptr) const;
	void setPointsFromCurrentLines()
	{
		//!再計算を行う
		//線の持つ点のインターセクションをチェックすることで，面の持つ点を間接的に取得する．
		//* default alignment
		/*         ret[2]
				 /\
				/  \
		  l[2] /    \ l[1]
			  /______\
		ret[0]  l[0]  ret[1]
	 */
		try
		{
			auto l = this->getLines();
			if (l.size() != 3)
			{
				std::string message = ERROR + "このFaceの持つ点の数は，" + std::to_string(l.size()) + "でfaceを形成していません";
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
			}

			this->setPoints(Intersection(l[0]->getPoints(), l[2]->getPoints())[0],
							Intersection(l[1]->getPoints(), l[0]->getPoints())[0],
							Intersection(l[2]->getPoints(), l[1]->getPoints())[0]);
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	/* ------------------------------------------------------ */
	/*
	PointsはsetBoundsと同時にかならず，更新される．
	そのため，問題なくりようできる．
	ただし，位置関係を決めておかないと問題が生じる．
	*/
	const V_netPp &getPoints() const { return this->Points; };
	V_netPp getPoints(const networkLine *const l) const
	{
		if (this->Lines[0] == l)
			return this->Points;
		else if (this->Lines[1] == l)
			return {this->Points[1], this->Points[2], this->Points[0]};
		else
			return {this->Points[2], this->Points[0], this->Points[1]};
	};
	/* ------------------------------------------------------ */
	//% タプル
	using T_LLL = std::tuple<networkLine *, networkLine *, networkLine *>;
	using T_PPP = std::tuple<networkPoint *, networkPoint *, networkPoint *>;
	const T_PPP &getPointsTuple() const { return this->PointsTuple; };
	T_PPP getPointsTuple(const networkPoint *p) const
	{
		if (std::get<1>(this->PointsTuple) == p)
			return T_PPP{std::get<1>(this->PointsTuple), std::get<2>(this->PointsTuple), std::get<0>(this->PointsTuple)};
		else if (std::get<2>(this->PointsTuple) == p)
			return T_PPP{std::get<2>(this->PointsTuple), std::get<0>(this->PointsTuple), std::get<1>(this->PointsTuple)};
		else
			return this->PointsTuple;
	};
	T_PPP getPointsTuple(const networkLine *const l) const
	{
		if (this->Lines[0] == l)
			return this->PointsTuple;
		else if (this->Lines[1] == l)
			return {std::get<1>(this->PointsTuple), std::get<2>(this->PointsTuple), std::get<0>(this->PointsTuple)};
		else
			return {std::get<2>(this->PointsTuple), std::get<0>(this->PointsTuple), std::get<1>(this->PointsTuple)};
	};
	T_LLL getLinesTuple() const { return {this->Lines[0], this->Lines[1], this->Lines[2]}; };
	T_LLL getLinesTuple(const networkLine *const l) const
	{
		if (l == this->Lines[0])
			return {this->Lines[0], this->Lines[1], this->Lines[2]};
		else if (l == this->Lines[1])
			return {this->Lines[1], this->Lines[2], this->Lines[0]};
		else
			return {this->Lines[2], this->Lines[0], this->Lines[1]};
	};
	/* ------------------------------------------------------ */
	int Find(netL *const l_IN) const
	{
		if (this->Lines[0] == l_IN)
			return 0;
		else if (this->Lines[1] == l_IN)
			return 1;
		else if (this->Lines[2] == l_IN)
			return 2;
		else
			return -1;
	};
	bool Switch(netL *const oldL, netL *const newL)
	{
		if (this->Lines[0] == oldL)
		{
			this->Lines[0] = newL;
			return true;
		}
		else if (this->Lines[1] == oldL)
		{
			this->Lines[1] = newL;
			return true;
		}
		else if (this->Lines[2] == oldL)
		{
			this->Lines[2] = newL;
			return true;
		}
		else
			return false;
	};
	std::vector<networkFace *> getNeighbors(const networkFace *const obj) const
	{
		return {(*this->Lines[0])(obj), (*this->Lines[1])(obj), (*this->Lines[2])(obj)};
	};
	/* ------------------------------------------------------ */
	std::tuple<T_PPP, T_LLL> getPointsLinesTuple(const networkLine *const l) const
	{
		return {getPointsTuple(l), getLinesTuple(l)};
	};
	/* ------------------------------------------------------ */
	netP *getPointFront(const netL *l) const { return getPoints(l)[1]; };
	netP *getPointBack(const netL *l) const { return getPoints(l)[0]; };
	netP *getPointOpposite(const netL *l) const { return getPoints(l)[2]; };
	/* ------------------------------------------------------ */
	netL *getLineFront(const networkPoint *p) const { return getLinesFrom(p)[0]; };
	netL *getLineBack(const networkPoint *p) const { return getLinesFrom(p)[2]; };
	netL *getLineOpposite(const networkPoint *p) const { return getLinesFrom(p)[1]; };
	/* ------------------------------------------------------ */
	netL *getLine(int i0, int i1) const
	{
		V_netPp Points = getPoints();
		for (const auto &line : Points[i0]->networkObject::getLines())
			if (Points[i1] == (*line)(Points[i0]))
				return line;
		return nullptr;
	};
	netL *getLine(const netL *l, int j = 0) const
	{
		int s = this->Lines.size();
		if (s != 3)
		{
			std::stringstream ss;
			ss << "辺の数が３ではない：this->Lines" << this->Lines;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
		}

		if (this->Lines[0] == l)
			return this->Lines[(0 + j + s) % 3];
		else if (this->Lines[1] == l)
			return this->Lines[(1 + j + s) % 3];
		else if (this->Lines[2] == l)
			return this->Lines[(2 + j + s) % 3];
		else
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
	netL *getLineFront(const netL *l) const { return getLine(l, 1); };
	netL *getLineBack(const netL *l) const { return getLine(l, -1); };
	netF *getFaceFront(const netL *l) { return (*getLine(l, 1))(this); };
	netF *getFaceBack(const netL *l) { return (*getLine(l, -1))(this); };
	//与えられたpを0番とするpointsベクトル
	//これがないと，objファイルを生成するときにエラーとなる．
	//というか，引数がないgetPointsを作り，this->Poitnsを出力するようにするとエラーになる．
	//つまり，getPoints()は，getPoints(nullptr)としてデフォルト引数でもう一つのgetPointsが実行されないといけない．
	//初めの段階では，this->Pointsが設定されていないことがあるからだ．
	//線を通して点を取得するかどうか選べるように，getPointsThroughLines関数を作るべきだろう
	V_netPp getPoints(const networkPoint *p) const
	{
		if (p == this->Points[0])
			return this->Points;
		if (p == this->Points[1])
			return {this->Points[1], this->Points[2], this->Points[0]};
		if (p == this->Points[2])
			return {this->Points[2], this->Points[0], this->Points[1]};
		else
			return this->Points;
	};

	//与えられたpから出発するlinesベクトル
	V_netLp getLinesFrom(const networkPoint *p) const
	{
		try
		{
			if (p == this->Points[0])
				return this->Lines;
			else if (p == this->Points[1])
				return RotateLeft(this->Lines, 1);
			else if (p == this->Points[2])
				return RotateLeft(this->Lines, 2);
			else
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	//------------------------
	netF *getFaceFront(const networkPoint *p) const { return (*(getLinesFrom(p)[0]))(this); };
	netF *getFaceBack(const networkPoint *p) const { return (*(getLinesFrom(p)[0]))(this); };
	netF *getFaceOpposite(const networkPoint *p) const { return (*(getLinesFrom(p)[0]))(this); };
	//---------------------
	V_netPp get6Points(netPp origin = nullptr) const
	{
		if (origin == nullptr || !MemberQ(this->Points, origin))
		{
			origin = this->Points[0];
		}

		// pが与えられた時は，この面の点として存在するかをまずチェックする．
		//なければ，適当に選んだ点を基準として選び，６点を返す．
		/*
		 *             0
		 *            / \
		 *           /   \
		 *         3/--l--\5
		 *         / \   / \
		 *        /   \ /   \
		 *      1/-----4-----\2
		 */
		try
		{
			V_netPp ret(6, nullptr);
			auto ps = this->getPoints();
			auto l = this->getLinesFrom(ps[0]);
			ret[2] = (*l[0])(this)->getPointOpposite(l[0]);
			ret[0] = (*l[1])(this)->getPointOpposite(l[1]);
			ret[1] = (*l[2])(this)->getPointOpposite(l[2]);
			ret[4] = ps[0];
			ret[5] = ps[1];
			ret[3] = ps[2];
			if (ret[5] == origin)
				return {ret[1], ret[2], ret[0], ret[4], ret[5], ret[3]};
			else if (ret[3] == origin)
				return {ret[2], ret[0], ret[1], ret[5], ret[3], ret[4]};
			else
				return ret;
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};

	V_netPp get6Points(const netLp l /*基準*/) const
	{
		/*            0
		 *            / \
		 *           /   \
		 *         3/--l--\5
		 *         / \   / \
		 *        /   \ /   \
		 *      1/-----4-----\2
		 */
		try
		{
			//上の図を返す
			return this->get6Points(this->getPointOpposite(l) /*4番目(基準)にしたい点*/);
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};

	VV_netPp getFour6Points(netPp p = nullptr) const
	{
		try
		{
			if (p == nullptr || !MemberQ(this->Points, p))
			{
				p = this->Points[0];
			}
			auto ls = this->getLinesFrom(p);
			/*              0
			 *              0
			 *             /  \
			 *            /    \
			 *          3/--l1--\5
			 *          / l2   l0 \
			 *         /   \  /    \
			 *  1    1/----4(p)------\2   2
			 */
			return {this->get6Points(p) /*上の図の添字の様に点を取得する*/,
					(*ls[1])(this)->get6Points(ls[1]) /* ls[1]は->1番*/,
					(*ls[2])(this)->get6Points(ls[2]) /* ls[2]は->0番*/,
					(*ls[0])(this)->get6Points(ls[0]) /*ls[0]は->2番*/};
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};

	std::tuple<netPp, netPp, netPp, netPp, netPp, netPp> get6PointsTuple(netPp p) const
	{
		try
		{
			if (p == nullptr || !MemberQ(this->Points, p))
			{
				p = this->Points[0];
			}
			auto ls = this->getLinesFrom(p);
			/*              0
			 *              0
			 *            /  \
			 *          3/-l1-\5
			 *         /  \  /  \
			 *  1    1/---4(p)---\2   2
			 */
			auto ps = (*ls[1])(this)->get6Points(ls[1]);
			return {ps[0], ps[1], ps[2], ps[3], ps[4], ps[5]};
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};

	V_netPp get12Points(netPp p = nullptr) const
	{
		try
		{
			if (p == nullptr || !MemberQ(this->Points, p))
			{
				p = this->Points[0];
			}
			auto ps = this->getFour6Points(p);
			/*
			 *       7-------0------6
			 *        \    /  \    /
			 *          \ /    \  /
			 *    8-----3/--l1--\5------11
			 *     \    / l2   l0 \	  /
			 *      \  /   \  /    \	 /
			 *       1/----4(p)-----\2
			 *         \    / \     /
			 *          \  /    \  /
			 *           9       10
			 */
			return Join(ps[0], {ps[1][1], ps[1][2], ps[2][1], ps[2][2], ps[3][1], ps[3][2]});
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};

	//-------------------------
	Tddd getMeanX() const { return Mean(getLocationsTuple()); };
	VV_d getLocations() const
	{
		V_netPp Points = getPoints();
		VV_d ret(Points.size());
		int i(0);
		// for (const auto &p : Points)
		// 	ret[i++] = p->xyz;
		for (const auto &p : Points)
			ret[i++] = {std::get<0>(p->X), std::get<1>(p->X), std::get<2>(p->X)};
		return ret;
	};
	T3Tddd getLocationsTuple() const
	{
		V_netPp ps = getPoints();
		return {ps[0]->getXtuple(),
				ps[1]->getXtuple(),
				ps[2]->getXtuple()};
	};
	T3Tddd getX_Vertices() const
	{
		return {std::get<0>(this->PointsTuple)->getXtuple(),
				std::get<1>(this->PointsTuple)->getXtuple(),
				std::get<2>(this->PointsTuple)->getXtuple()};
	};
	//-------------------------
	//つぎつぎと入社して，線との交点をかえし面の切り取りに使う，最終位置を返す//改善が必要
	netL *getCrossLine(const V_d &s, const V_d &vec)
	{
		V_d n = getNormal(), theta(0);
		for (const auto &v : this->getLocations())
			theta.emplace_back(MyVectorAngle(vec /*基準とするベクトル*/, v - s /*対象となるベクトル*/, n /*回転方向を決める法線ベクトル*/));

		double ccwAngle = 2. * M_PI, cwAngle = -2. * M_PI;
		int ccwIndex = 0, cwIndex = 0;
		std::cout << Red << theta << reset << std::endl;
		// 右左（マイナスとプライス）でそれぞれ近い頂点番号（角度が0に近い）を探す
		for (size_t i = 0; i < theta.size(); i++)
			if (theta[i] >= 0 && ccwAngle > theta[i])
			{
				ccwAngle = theta[i];
				ccwIndex = i;
			}
			else if (theta[i] <= 0 && cwAngle < theta[i])
			{
				cwAngle = theta[i];
				cwIndex = i;
			}
		auto Points = getPoints();
		// V_d ccwp = Points[ccwIndex]->xyz;
		// V_d cwp = Points[cwIndex]->xyz;
		V_d ccwp = {std::get<0>(Points[ccwIndex]->X), std::get<1>(Points[ccwIndex]->X), std::get<2>(Points[ccwIndex]->X)};
		V_d cwp = {std::get<0>(Points[cwIndex]->X), std::get<1>(Points[cwIndex]->X), std::get<2>(Points[cwIndex]->X)};
		double c = MyVectorAngle(vec, ccwp - cwp, n);
		double r = Norm(ccwp - s) * sin(c - ccwAngle) / sin(c);
		std::cout << Red << "c:" << c / M_PI * 180 << reset << std::endl;
		std::cout << Red << "r:" << r << reset << std::endl;
		return getLine(ccwIndex, cwIndex); // LineのFace2Faceを使って，次に進むFaceを得る
	};
	//-------------------------
	//-------------------------
	pathInfo getPathInfo(const V_d &s, const V_d &vec)
	{
		V_d n = getNormal(), theta(0);
		for (const auto &v : this->getLocations())
			theta.emplace_back(MyVectorAngle(vec /*基準とするベクトル*/, v - s /*対象となるベクトル*/, n /*回転方向を決める法線ベクトル*/));

		double ccwAngle = 2. * M_PI, cwAngle = -2. * M_PI;
		int ccwIndex = 0, cwIndex = 0;
		// 右左（マイナスとプライス）でそれぞれ近い頂点番号（角度が0に近い）を探す
		for (size_t i = 0; i < theta.size(); i++)
			if (theta[i] >= 0 && ccwAngle > theta[i])
			{
				ccwAngle = theta[i];
				ccwIndex = i;
			}
			else if (theta[i] <= 0 && cwAngle < theta[i])
			{
				cwAngle = theta[i];
				cwIndex = i;
			}
		/*           /\ccwp
		 *           /| \
		 *          / |a \c
		 *         /  s=r=x==len====> vec
		 *        /    \b  \pi-c=incidentAngle
		 *       /       \  \      a=ccwAngle (positive)
		 *      ------------> cwp  b=cwAngle (negative)
		 */
		auto Points = getPoints();

		V_d ccwp = {std::get<0>(Points[ccwIndex]->X), std::get<1>(Points[ccwIndex]->X), std::get<2>(Points[ccwIndex]->X)};
		V_d cwp = {std::get<0>(Points[cwIndex]->X), std::get<1>(Points[cwIndex]->X), std::get<2>(Points[cwIndex]->X)};

		// auto ccwp = Points[ccwIndex]->xyz;
		// auto cwp = Points[cwIndex]->xyz;
		double c = MyVectorAngle(vec, ccwp - cwp, n);
		double r = Norm(ccwp - s) * sin(c - ccwAngle) / sin(c);
		V_d x = r / Norm(vec) * vec + s; // cross point coordinate

		netL *line = getLine(ccwIndex, cwIndex);
		if (line == NULL)
			std::cout << line << std::endl;
		double remainlen = Norm(vec) - r;
		return pathInfo({this, (remainlen <= 0) ? NULL : (*line)(this)},
						{s, cwp, x},
						M_PI - c,
						remainlen);
	};
	pathInfo getPathInfo(const pathInfo &pathinfo)
	{
		auto X = pathinfo.xyz[1] - pathinfo.xyz[2];
		Tddd incidentVec = Dot(RotationMatrix(pathinfo.incidentAngle, this->getNormalTuple()),
							   Tddd{X[0], X[1], X[2]});
		return getPathInfo(pathinfo.xyz[2],
						   ToVector(pathinfo.r * incidentVec / Norm(incidentVec)));
	};
	//-------------------------
	networkFace *linkedFace(const int i) { return (*Lines[i])(this); };
	//-------------------------

	netFp divide(netLp DivL, netLp newDivL, netLp newMidL, int type);
};
/*networkFace_code*/
//@ ------------------------ 抽出用関数など ----------------------- */
std::vector<Tddd> extNormals(const V_netFp &ps)
{
	std::vector<Tddd> ret;
	for (const auto &p : ps)
		ret.emplace_back(p->getNormalTuple());
	return ret;
};
Tddd vectorToTriangle(const networkFace *f, const Tddd &a)
{
	auto n = f->getNormalTuple();
	return n * Dot(n, f->getXtuple() - a);
};
//@ ------------------------------------------------------ */
/* ------------------------------------------------------ */
inline void Buckets<networkPoint>::add(const V_netPp &ps)
{
	for (const auto &p : ps)
		add(p->getXtuple(), p);
};
inline void Buckets<networkPoint>::add(const std::unordered_set<networkPoint *> &ps)
{
	for (const auto &p : ps)
		add(p->getXtuple(), p);
};
/* ------------------------------------------------------ */
std::vector<netL *> link(const V_netPp &obj, Network *net)
{
	try
	{
		std::vector<netL *> ret;
		int s = obj.size();
		for (int i = 0; i < s; i++)
		{
			auto l = link(obj[i], obj[(i + 1) % s], net);
			ret.emplace_back(l);
		}
		return ret;
	}
	catch (const error_message &e)
	{
		Print(obj);
		if (DuplicateFreeQ(obj))
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "no duplication.....???");
		else
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "duplication found!");
	}
};
netL *unlink(netP *obj, netP *obj_)
{
	if (obj_ == obj)
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "a point is trying to unlink itself!");

	if (!obj || !obj_ /*if NULL*/)
		return nullptr;

	auto line = obj->getLineBetween(obj_);
	auto line_ = obj_->getLineBetween(obj);

	if ((line == line_) && line)
	{
		line->Erase(obj);
		line->Erase(obj_);
		obj->Erase(line);
		obj_->Erase(line);
		return line;
	}
	else
	{
		std::cout << Red << obj << " and " << obj_ << " are not linked" << reset << std::endl;
		return nullptr;
	}
};

/////////////////////////////////////////////////////////////////////
#include "searcher.hpp"
//////////////////////////////////////////////////////////////////////
/* Networkは持っているPointsから情報を取り出すメソッドを提供する*/
//*  @-@-@ */
//*  |\|/| */
//*  @-@-@ */
//*  |/|\| */
//*  @-@-@ */
/*Network_code*/
#define use_binary_search_in_Network
// binary_searchを使えば，重複しないように保存する作業の時間が短縮される．
class Network : public object3D
{
public:
	//! ------------------------------------------------------ */
	//!                          接触の判別                      */
	//! ------------------------------------------------------ */
	std::unordered_set<networkPoint *> getContactPointsOfPoints() const
	{
		std::unordered_set<networkPoint *> ret;
		for (const auto &p : this->getPoints())
			ret.insert(begin(p->getContactPoints()), end(p->getContactPoints()));
		return ret;
	};
	std::unordered_set<networkPoint *> getContactPointsOfPoints(const std::vector<Network *> &nets) const
	{
		/*
		getContactPointsOfPointsは，自身の保有するPointsが接した点を返す．
		Pointsの保有するmap_Net_ContactPointsから指定されたNetworkのPointsを抽出している．
		*/
		std::unordered_set<networkPoint *> ret;
		for (const auto &p : this->getPoints())
			ret.insert(begin(p->getContactPoints(nets)), end(p->getContactPoints(nets)));
		return ret;
	};

public:
	//% ------------------------------------------------------ */
	void makeMirroredPoints(const Buckets<networkFace> &B_face, const double mirroring_distance)
	{
		auto points = this->getPoints();
		for (const auto &p : points)
			p->makeMirroredPoints(B_face, mirroring_distance);
	};
	void clearMirroredPoints()
	{
		auto points = this->getPoints();
		for (const auto &p : points)
			p->clearMirroredPoints();
	};
	//% ------------------------------------------------------ */
private:
	Buckets<networkFace> BucketFaces;
	Buckets<networkPoint> BucketParametricPoints;
	Buckets<networkPoint> BucketPoints;

public:
	const Buckets<networkFace> &getBucketFaces() const { return BucketFaces; };
	const Buckets<networkPoint> &getBucketPoints() const { return BucketPoints; };
	const Buckets<networkPoint> &getBucketParametricPoints() const { return BucketParametricPoints; };
	void makeBucketFaces(const double spacing)
	{
		this->setBounds();
		this->BucketFaces.clear(); //こうしたら良くなった
		this->BucketFaces.set(this->bounds(), spacing);
		for (const auto &f : this->getFaces())
		{
			auto [p0, p1, p2] = f->getPointsTuple();
			for (const auto [xyz, t0t1] : triangleIntoPoints(T3Tddd{p0->getXtuple(), p1->getXtuple(), p2->getXtuple()}, spacing / 3.))
				this->BucketFaces.add(xyz, f);
		}
	};

	void makeBucketPoints(const double spacing)
	{
		this->setBounds();
		this->BucketPoints.resize(this->bounds(), spacing);
		this->BucketPoints.add(this->getPoints());
	};
	void makeBucketParametricPoints(const double spacing)
	{
		this->setBounds();
		this->BucketParametricPoints.resize(this->bounds(), spacing);
		this->BucketParametricPoints.add(this->getParametricPoints());
		std::cout << "this->getParametricPoints().size()=" << this->getParametricPoints().size() << std::endl;
	};
	void makeBucketParametricPoints(const T3Tdd &bounds, const double spacing)
	{
		this->BucketParametricPoints.resize(bounds, spacing);
		this->BucketParametricPoints.add(this->getParametricPoints());
	};
	void clearBucketParametricPoints()
	{
		this->BucketParametricPoints.clear();
	};

public:
	T6d forced_velocity;
	T6d forced_acceleration;

public:
	/* ------------------------------------------------------ */
	T6d force;
	T6d inertia;
	T6d velocity;
	T6d acceleration;
	double mass;
	/* ------------------------------------------------------ */
	//! 固定された空間座標におけるベクトルであることを頭に入れておくこと．
	//! 回転，移動をする物体の座標系ではないので，固定座標にとって，回転前後でinertiaは書き換える必要がある．
	//! inertiaの慣性モーメントはそのまま固定座標における回転行列をかけて，更新すればいい
	T6d &I = this->inertia;
	T6d &F = this->force;
	T6d &A = this->acceleration;
	T6d &V = this->velocity;
	/* ------------------------------------------------------ */
	Tddd center_of_mass;
	/* ------------------------------------------------------ */
	// calcPhysicalProperties試作．
	void calcPhysicalProperties()
	{
		this->force = {0, 0, 0, 0, 0, 0};
		this->inertia = {0, 0, 0, 0, 0, 0};
		this->mass = 0.;
		for (const auto &p : this->Points)
		{
			// for (auto i = 0; i < 3; ++i)
			// 	this->force[i] += p->force[i];
			std::get<0>(this->force) += std::get<0>(p->force);
			std::get<1>(this->force) += std::get<1>(p->force);
			std::get<2>(this->force) += std::get<2>(p->force);
			this->mass += p->mass; // total mass
		}
		std::get<0>(this->inertia) = std::get<1>(this->inertia) = std::get<2>(this->inertia) = this->mass; // total mass

		// this->center_of_mass = {0, 0, 0};
		// //! -------------------------------------- */
		// for (const auto &p : this->Points)
		// 	this->center_of_mass += p->getXtuple();
		// this->center_of_mass /= this->mass;
		// //! -------------------------------------- */
		Tddd tmp = {0, 0, 0};
		for (const auto &p : this->Points)
		{
			tmp = p->getXtuple() - this->center_of_mass;
			std::get<3>(this->force) += std::get<0>(p->force) * std::get<0>(tmp);
			std::get<4>(this->force) += std::get<1>(p->force) * std::get<1>(tmp);
			std::get<5>(this->force) += std::get<2>(p->force) * std::get<2>(tmp);
			//
			std::get<3>(this->inertia) += p->mass * std::pow(std::get<0>(tmp), 2.);
			std::get<4>(this->inertia) += p->mass * std::pow(std::get<1>(tmp), 2.);
			std::get<5>(this->inertia) += p->mass * std::pow(std::get<2>(tmp), 2.);
		}
	};
	//@ ------------------- 物体全体にかかる力の計算方法２つ ------------------- */
	//まずは，点にかかる力を計算してから使うこと．
	void sumForceOfPoints()
	{
		//%this->center_of_massが設定してある必要がある
		//%ネットワークの各点にかかっている力の総和を計算する
		this->force = {0, 0, 0, 0, 0, 0};
		for (const auto &p : this->Points)
		{
			auto F = p->force;
			auto r = (p->getXtuple() - this->center_of_mass);
			auto [Fx, Fy, Fz, Fx_, Fy_, Fz_] = F;
			auto [Nx, Ny, Nz] = Cross(r, Tddd{Fx, Fy, Fz}); //モーメント
			std::get<0>(this->force) += Fx;
			std::get<1>(this->force) += Fy;
			std::get<2>(this->force) += Fz;
			std::get<3>(this->force) += Nx;
			std::get<4>(this->force) += Ny;
			std::get<5>(this->force) += Nz;
		}
	};
	void integrateForceOnFace()
	{
		/*(*checked by mathematica*)
		f0 = {f0x, f0y, f0z}
		f1 = {f1x, f1y, f1z}
		f2 = {f2x, f2y, f2z}
		x0 = {x0x, x0y, x0z}
		x1 = {x1x, x1y, x1z}
		x2 = {x2x, x2y, x2z}
		f[x0_, x1_] := Dot[{f0, f1, f2}, {x0, x1, 1 - x0 - x1}]
		X[z0_, z1_] := Dot[{x0, x1, x2}, {z0, z1, 1 - z0 - z1}]
		CForm@Flatten[
		Join[Simplify@
			Integrate[f[z0, z1], {\[Xi]0, 0, 1}, {z1, 0, 1 - z0}],
		Simplify@
			Integrate[
			Cross[X[z0, z1], f[z0, z1]], {z0, 0, 1}, {z1, 0, 1 - z0}]]]
		*/
		this->force = {0, 0, 0, 0, 0, 0};
		for (const auto &f : this->Faces)
		{
			auto [p0, p1, p2] = f->getPointsTuple();
			auto [f0x, f0y, f0z, f0x_, f0y_, f0z_] = p0->force;
			auto [f1x, f1y, f1z, f1x_, f1y_, f1z_] = p1->force;
			auto [f2x, f2y, f2z, f2x_, f2y_, f2z_] = p2->force;
			auto [x0x, x0y, x0z] = p0->getXtuple() - this->center_of_mass;
			auto [x1x, x1y, x1z] = p1->getXtuple() - this->center_of_mass;
			auto [x2x, x2y, x2z] = p2->getXtuple() - this->center_of_mass;
			this->force += T6d{(f0x + f0y + f0z) / 6.,
							   (f1x + f1y + f1z) / 6.,
							   (f2x + f2y + f2z) / 6.,
							   (f2z * x1x + f2z * x1y + 2 * f2z * x1z + f2x * (2 * x1x + x1y + x1z) + f2y * (x1x + 2 * x1y + x1z) - 2 * f1x * x2x - f1y * x2x - f1z * x2x - f1x * x2y - 2 * f1y * x2y - f1z * x2y -
								f1x * x2z - f1y * x2z - 2 * f1z * x2z) /
								   24.,
							   (-(f2z * x0x) - f2z * x0y - 2 * f2z * x0z - f2x * (2 * x0x + x0y + x0z) - f2y * (x0x + 2 * x0y + x0z) + 2 * f0x * x2x + f0y * x2x +
								f0z * x2x + f0x * x2y + 2 * f0y * x2y + f0z * x2y + f0x * x2z + f0y * x2z + 2 * f0z * x2z) /
								   24.,
							   (f1z * x0x + f1z * x0y + 2 * f1z * x0z + f1x * (2 * x0x + x0y + x0z) + f1y * (x0x + 2 * x0y + x0z) - 2 * f0x * x1x - f0y * x1x - f0z * x1x - f0x * x1y - 2 * f0y * x1y - f0z * x1y -
								f0x * x1z - f0y * x1z - 2 * f0z * x1z) /
								   24.};
		}
	};
	//@ ------------------------------------------------------ */
	void calcAccelFromForce()
	{
		//! 慣性を設定しておく必要がある．
		this->acceleration = this->force / this->inertia;
	};

// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM
/*Network_BEM_detail
境界要素法用のNetworkには，以下が設定されている．

* bool BASE;
* bool CORNER;
* bool Dirichlet;
* bool Neumann;

`networkPoint`も`this->getNetwork()->setB()`などとして，自身の境界条件判断できる．
Network_BEM_detail*/
/*Network_BEM_code*/
#ifdef BEM
public:
	double _current_time_;
	bool isRigidBody;
	bool isMovementPredefined;

	const double move_amplitude = 0.4;
	T6d velocityPredefined()
	{
		double t = _current_time_;
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

	Tddd translationPredefined()
	{
		double t = _current_time_;
		double s = M_PI / 2.;
		double k = M_PI / 1.;
		/* ------------------------------------------------------ */
		// Tddd move_dir = {cos(k * t), sin(k * t), 0.};
		// return move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) * move_dir;
		/* ------------------------------------------------------ */
		Tddd move_dir = Normalize(Tddd{1., 1., 0.});
		return move_amplitude * exp(-t) * (sin(k * t - s) - sin(-s)) * move_dir;
	};

private:
	bool BASE;
	bool CORNER;
	bool Dirichlet;
	bool Neumann;

public:
	void setB(const bool TorF) { this->BASE = TorF; };
	void setC(const bool TorF) { this->CORNER = TorF; };
	void setD(const bool TorF) { this->Dirichlet = TorF; };
	void setN(const bool TorF) { this->Neumann = TorF; };
	//
	bool isB() const { return this->BASE; };
	bool isC() const { return this->CORNER; };
	bool isD() const { return this->Dirichlet; };
	bool isN() const { return this->Neumann; };
	//
	void setCornerNeumann()
	{
		this->setB(false);
		this->setC(true);
		this->setD(false);
		this->setN(true);
	};
	void setBaseDirichlet()
	{
		this->setB(true);
		this->setC(false);
		this->setD(true);
		this->setN(false);
	};
	void setBaseNeumann()
	{
		this->setB(true);
		this->setC(false);
		this->setD(false);
		this->setN(true);
	};
	//
	V_d nabla_phi;

	V_d meanX;

#endif
	/*Network_BEM_code*/
	// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM

private:
	using V_netLp = std::vector<networkLine *>;
	std::string name;
	// 接続している点をLinesをとおして素早く参照できる
protected:
	// protectedであり，Network::addは重複を許さないので，PointsとFacesに重複はない
	//  std::vector<std::shared_ptr<networkPoint>> Points;
	//  V_netPp Points;
	//  V_netFp Faces;
	std::unordered_set<networkPoint *> Points;
	std::unordered_set<networkFace *> Faces;
	V_netPp PointGarbage;
	V_netFp FaceGarbage;

public:
	/* ------------------------------------------------------ */
	// bool checkup()
	// {
	// 	if (!DuplicateFreeQ(this->Points))
	// 		return false;
	// 	if (!DuplicateFreeQ(this->Faces))
	// 		return false;
	// 	return true;
	// };
	/* ------------------------------------------------------ */
	// void sortPoints() { std::sort(this->Points.begin(), this->Points.end()); };
	// void sortFaces() { std::sort(this->Faces.begin(), this->Faces.end()); };

	// 2021/08/10追加
	bool isMember(networkPoint *p_IN) const { return (this->Points.find(p_IN) != this->Points.end()); };
	bool isMember(networkFace *f_IN) const { return (this->Faces.find(f_IN) != this->Faces.end()); };
	//
	void add(netPp const p_IN) { this->Points.emplace(p_IN); };
	void add(netFp const f_IN) { this->Faces.emplace(f_IN); };
	void add(const V_netPp &ps_IN)
	{
		for (const auto &p : ps_IN)
			this->add(p);
	};
	//
	void add(const V_netFp &fs_IN)
	{
		for (const auto &f : fs_IN)
			this->add(f);
	};
	//
	void Erase(networkFace *f) { this->Faces.erase(f); };
	void Erase(const V_netFp &fs)
	{
		for (auto &f : fs)
			Erase(f);
	};
	void Erase(networkPoint *p) { this->Points.erase(p); };

	// netPp takeNearestPoint(const V_d &cardinal) const { return network::takeNearest(this->Points, cardinal); };
	// netFp takeNearestFace(const V_d &cardinal) const { return network::takeNearest(this->Faces, cardinal); };
	//------------------
	// #if defined(use_binary_search_in_Network)
	// 	int getPointIndex(const netPp p) const
	// 	{
	// 		auto it = std::lower_bound(this->Points.begin(), this->Points.end(), p);
	// 		if (it == this->Points.end() || *it != p)
	// 			return -1;
	// 		else
	// 			return std::distance(this->Points.begin(), it);
	// 	};
	// 	V_i getPointIndex(const V_netPp &ps) const
	// 	{
	// 		V_i ret(ps.size(), 0);
	// 		for (auto i = 0; i < ps.size(); i++)
	// 			ret[i] = getPointIndex(ps[i]);
	// 		return ret;
	// 	};
	// 	VV_i getPointIndex(const VV_netPp &ps) const
	// 	{
	// 		VV_i ret(ps.size());
	// 		for (auto i = 0; i < ps.size(); i++)
	// 			ret[i] = getPointIndex(ps[i]);
	// 		return ret;
	// 	};
	// 	int getFaceIndex(const netFp p) const
	// 	{
	// 		auto it = std::lower_bound(this->Faces.begin(), this->Faces.end(), p);
	// 		if (it == this->Faces.end() || *it != p)
	// 			return -1;
	// 		else
	// 			return std::distance(this->Faces.begin(), it);
	// 	};
	// 	V_i getFaceIndex(const V_netFp &fs) const
	// 	{
	// 		V_i ret(fs.size(), 0);
	// 		for (auto i = 0; i < fs.size(); i++)
	// 			ret[i] = getFaceIndex(fs[i]);
	// 		return ret;
	// 	};
	// 	VV_i getFaceIndex(const VV_netFp &fs) const
	// 	{
	// 		VV_i ret(fs.size());
	// 		for (auto i = 0; i < fs.size(); i++)
	// 			ret[i] = getFaceIndex(fs[i]);
	// 		return ret;
	// 	};
	// #endif
	//-------------------------
	//% Networkコンストラクタ1
	Network(const std::string &name_IN = "no_name")
		: object3D(geometry::CoordinateBounds({0, 0, 0})),
		  name(name_IN),
		  /* ------------------------------------------------------ */
		  force({0., 0., 0., 0., 0., 0.}),
		  inertia({0., 0., 0., 0., 0., 0.}),
		  forced_acceleration({0., 0., 0., 0., 0., 0.}),
		  forced_velocity({0., 0., 0., 0., 0., 0.}),
		  acceleration({0., 0., 0., 0., 0., 0.}),
		  velocity({0., 0., 0., 0., 0., 0.}),
		  mass(0.),
		  center_of_mass({0., 0., 0.}),
		  BucketFaces(geometry::CoordinateBounds({0., 0., 0.}), 1.),
		  BucketPoints(geometry::CoordinateBounds({0., 0., 0.}), 1.),
		  BucketParametricPoints(geometry::CoordinateBounds({0., 0., 0.}), 1.){};
	Network(const V_Netp &base_nets, const std::string &name_IN);
	~Network();
	//
	//% ------------------------- obj ------------------------ */
	//% NetwrokObjは，Networkとしてコンストラクトするために，下のコンストラクタように修正した．*/
	Network(const std::string &filename, const std::string &name_IN = "no_name")
		: object3D(geometry::CoordinateBounds({0, 0, 0})),
		  name(name_IN),
		  force({0., 0., 0., 0., 0., 0.}),
		  inertia({0., 0., 0., 0., 0., 0.}),
		  forced_acceleration({0., 0., 0., 0., 0., 0.}),
		  forced_velocity({0., 0., 0., 0., 0., 0.}),
		  acceleration({0., 0., 0., 0., 0., 0.}),
		  velocity({0., 0., 0., 0., 0., 0.}),
		  mass(0.),
		  center_of_mass({0., 0., 0.}),
		  BucketFaces(geometry::CoordinateBounds({0., 0., 0.}), 1.),
		  BucketPoints(geometry::CoordinateBounds({0., 0., 0.}), 1.),
		  BucketParametricPoints(geometry::CoordinateBounds({0., 0., 0.}), 1.)
	{
		std::vector<std::vector<std::string>> read_line;
		Load(filename, read_line, {"    ", "   ", "  ", " "});
		LoadObj objLoader;
		objLoader.load(read_line);

		int i = 0;
		auto xyz = objLoader.v;
		double pn = RandomReal({-1., 1.}) > 0 ? 1. : -1;
		double rand = pn * RandomReal({0.5, 1.}) * 1E-13;
		for (const auto &w : objLoader.v)
			xyz[i++] = {w[0] + rand, w[2] + rand, w[1] + rand};
		Network::generateNetwork(xyz, objLoader.f_v);
	};
	//% ------------------------------------------------------ */

	std::string getName() { return this->name; };
	V_netPp linkXPoints(Network &water, Network &obj);
	/* ------------------------------------------------------ */
	void generateNetwork(const VV_d &v_IN, const VV_i &f_v_IN = {})
	{
		/*
		 * setPointsやsetFacesは，deleteduplicatesやsortを行わない．
		 * そのため，それらを使ってnetworkを作成すると，問題が生じるので使わない（privateにしている）．
		 * ネットワークを作成する際は，このgenerateNetworkを使う．
		 */
		Print("generating network..", Red);
		Print(__PRETTY_FUNCTION__, red);
		auto points_consists_faces = this->setPoints(v_IN);
		this->setFaces(f_v_IN, points_consists_faces); // indexの書き換えも可能だがする必要は今のところない

		// this->Points = DeleteDuplicates(this->Points);

		for (const auto &p : this->Points)
			p->network = this;
		for (const auto &f : this->Faces)
			f->network = this;
		for (const auto &l : getLines())
			l->network = this;

		// this->sortPoints();
		// this->sortFaces();

		this->setBounds();
		Print("Network has been generated", Magenta);
		this->displayStates();
	};
	/* ------------------------------------------------------ */
	void displayStates()
	{
		int size = 20;
		V_i num(size, 0), pointsLines(size, 0), facesPoints(size, 0), facesLines(size, 0);

		// 点の持つ線の数をカウント
		for (const auto &p : this->Points)
			for (auto i = 0; i < size; i++)
				if (p->networkObject::getLines().size() == (size_t)i)
					pointsLines[i]++;
		// 面の持つ線と点の数をカウント
		for (const auto &f : Faces)
		{
			auto f_Points = f->getPoints();
			for (auto i = 0; i < size; i++)
			{
				if (f_Points.size() == (size_t)i)
					facesPoints[i]++;
				if (f->networkObject::getLines().size() == (size_t)i)
					facesLines[i]++;
			}
		}

		V_i connection(3);
		setLinesStatus(true);
		for (const auto &p : Points)
			for (const auto &line : p->getLines_toggle(true))
			{
				if ((line->getFaces()).size() == 0)
					connection[0]++;
				else if ((line->getFaces()).size() == 1)
					connection[1]++;
				else if ((line->getFaces()).size() == 2)
					connection[2]++;
			}

		std::cout << Magenta << std::setw(25) << "this : " << Grid({this->getName()}) << ", " << this << std::endl;
		std::cout << Magenta << std::setw(25) << "Lines.size() : " << getLines().size() << std::endl;
		std::cout << Magenta << std::setw(25) << "Points.size() : " << Points.size() << std::endl;
		std::cout << Magenta << std::setw(25) << " Faces.size() : " << Faces.size() << std::endl;
		std::cout << Red << std::setw(25) << "                  " << GridVector(Subdivide(0., (double)size - 1., size - 1), 6) << std::endl;
		std::cout << magenta << std::setw(25) << "Lines of points : " << GridVector(pointsLines, 6) << std::endl;
		std::cout << magenta << std::setw(25) << "Points of faces : " << GridVector(facesPoints, 6) << std::endl;
		std::cout << magenta << std::setw(25) << " Lines of faces : " << GridVector(facesLines, 6) << std::endl;
		std::cout << magenta << std::setw(25) << " Connection of faces : " << GridVector(connection, 6);
		if (connection[0] == 0 && connection[1] == 0 && connection[2] != 0)
			std::cout << Blue << " <- 全ての線が2つの面と接続しているため，格子は閉じた面を形成している";
		std::cout << reset << std::endl;
		std::cout << "CoordinateBounds: " << this->bounds << std::endl;
	};
	/* ------------------------------------------------------ */
	double rotation_offset;
	V_d translation_offset;
	double rotation;
	V_d translation;
	void rotate(const double theta, const Tddd &axis)
	{
		for (auto &p : this->getPoints())
			p->X = Dot(RotationMatrix(theta, axis), p->getXtuple());
		setBounds();
	};
	/* ------------------------------------------------------ */
	/*             ネットワークの姿勢に関する変数とメソッド          */
	/* ------------------------------------------------------ */
	void rotate(const Quaternion &Q, const Tddd &barycenter = {0., 0., 0.})
	{
		try
		{
			auto tmp = this->getPoints();
			for (auto &p : tmp)
				p->X = Q.Rv(p->X - barycenter) + barycenter;
			setBounds();
		}
		catch (std::exception &e)
		{
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		}
	};
	void scale(const double ratio, const Tddd &center = {0, 0, 0})
	{

		try
		{
			for (auto &p : this->getPoints())
				p->setX(center + (p->X - center) * ratio);
			this->setBounds();
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};

	void scale(const Tddd &ratio, const Tddd &center = {0, 0, 0})
	{
		try
		{
			for (auto &p : this->getPoints())
				p->setX(center + (p->X - center) * ratio);
			this->setBounds();
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	void translate(const V_d &direction)
	{
		try
		{
			for (auto &p : this->getPoints())
			{
				std::get<0>(p->X) += direction[0];
				std::get<1>(p->X) += direction[1];
				std::get<2>(p->X) += direction[2];
			}
			setBounds();
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	void translateFromInitialX(const Tddd &translation)
	{
		try
		{
			for (auto &p : this->getPoints())
				p->X = p->initialX + translation;
			setBounds();
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	void resetInitialX()
	{
		for (auto &p : this->getPoints())
			p->initialX = p->X;
	};

	void reverseNormal()
	{
		for (auto f : this->Faces)
		{
			f->reverseNormal();
		}
	};
	/* ------------------------------------------------------ */
	netF *getNearestFace(const V_d &xyz)
	{
		netF *ret;
		double distance = 1E+30;
		for (const auto &f : this->Faces)
			if (distance > Norm(f->getX() - xyz))
			{
				distance = Norm(f->getX() - xyz);
				ret = f;
			}
		return ret;
	};

	const std::unordered_set<networkPoint *> &getPoints() const { return this->Points; };
	const std::unordered_set<networkFace *> &getFaces() const { return this->Faces; };
	std::unordered_set<networkPoint *> getParametricPoints() const
	{
		std::unordered_set<networkPoint *> ret, tmp;
		for (const auto &f : this->getFaces())
		{
			tmp = f->getParametricPoints();
			ret.insert(tmp.begin(), tmp.end());
		}
		return ret;
	};
	Tddd getMeanX() const
	{
		Tddd ret = {0, 0, 0};
		for (const auto &p : this->Points)
			ret += p->getXtuple();
		return ret / this->Points.size();
	};
	VV_d getLocations() const
	{
		VV_d ret;
		getLocations(ret);
		return ret;
	};
	void getLocations(VV_d &ret) const
	{
		ret.resize(Points.size());
		int i(0);
		// for (const auto &p : Points)
		// 	ret[i++] = p->xyz;
		for (const auto &p : Points)
			ret[i++] = {std::get<0>(p->X), std::get<1>(p->X), std::get<2>(p->X)};
	};
	std::vector<Tddd> getLocationsTuple() const
	{
		std::vector<Tddd> ret(this->Points.size());
		int i = 0;
		for (const auto &p : this->Points)
			ret[i++] = p->X;
		return ret;
	};
	//-------------------------
	// get lines
	void setLinesStatus(bool TorF)
	{
		for (const auto &p : Points)
			p->setLinesStatus(TorF);
	};
#define debug_getLines
	V_netLp getLines() const
	{
		V_netLp ret;
		try
		{
			std::unordered_set<networkLine *> tmp;
#ifdef debug_getLines
			for (const auto &p : this->Points)
			{
				auto lines = p->getLines();
				// if (lines.empty())
				// {
				// 	std::stringstream ss;
				// 	ss << "networkPoint ポインター " << p << "は，networkLineを持っていません．"
				// 	   << "このようなpointは削除するべきです";
				// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
				// }
				for (const auto &l : lines)
					tmp.insert(l);
			}
#else
			for (const auto &p : this->Points)
				for (const auto &l : p->getLines())
					tmp.insert(l);
#endif
			ret.assign(tmp.begin(), tmp.end());
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
		return ret;
	};

	std::unordered_set<networkLine *> getLinesUO() const
	{
		std::unordered_set<networkLine *> ret;
		for (const auto &p : this->Points)
			for (const auto &l : p->getLines())
				ret.emplace(l);
		return ret;
	};

	V_netLp getLinesIntxn() const
	{
		V_netLp ret(0);
		for (const auto &ps : this->getPoints())
			for (const auto &l : ps->networkObject::getLines())
				if (l->isIntxn())
					ret.emplace_back(l);
		return DeleteDuplicates(ret);
	};
	//======================================
	/*setXPoints_detail
`setXPoints`は，thisとtargetの`networkFace`と`networkLine`の干渉を検知して，干渉点として`networkPoint`を新たに作成する．

- 干渉点`networkPoint`は，干渉した`networkFace`と`networkLine`の`XPoints`に`push`される．
- 干渉点として作成される`networkPoint`には，`xline`と`xface`が与えられる．これらは，デフォルトで`xline=nullptr`，`xface=nullptr`となっている．

setXPoints_detail*/
	/*setXPoints_code*/
	void setXPoints(Network &target, Network *interactionNet);
	/*setXPoints_code*/
	//======================================
	void setBounds()
	{
		try
		{
			// network::setStatus(this->getPoints(), true);
			for (const auto &p : this->getPoints())
				p->setBoundsSingle();
			for (const auto &l : this->getLines())
				l->setBoundsSingle();
			for (const auto &f : this->getFaces())
				f->setBounds();
			object3D::setBounds(geometry::CoordinateBounds(getLocationsTuple()));
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	// void setBounds(VV_d coordinatebounds) {
	// 	object3D::setBounds(coordinatebounds);
	// };
	/* ------------------------------------------------------ */
private:
	V_netPp setPoints(const VV_d &v_IN)
	{
		Print("setting points...", Red);
		Print("finding overlaps", Red);
		double eps = 1E-10;
		int overlap_size(0), i(0), undefined = -1, s = v_IN.size();
		V_i no_overlap_indices(v_IN.size(), undefined);
		for (const auto &v : v_IN)
		{
			for (auto j = i + 1; j < s; j++)
			{
				if (no_overlap_indices[j] == undefined && std::abs(v[0] - v_IN[j][0]) < eps && std::abs(v[1] - v_IN[j][1]) < eps && std::abs(v[2] - v_IN[j][2]) < eps)
				{ // if overlaping
					no_overlap_indices[j] = i;
					overlap_size++;
				}
			}
			if (no_overlap_indices[i] == undefined)
				no_overlap_indices[i] = i;
			i++;
		}

		Print("overlap_size: " + std::to_string(overlap_size), Red);
		i = 0;
		V_netPp points_consists_faces(v_IN.size());
		for (auto const &j : no_overlap_indices)
		{
			//インデックスが実際の配列番号よりも小さいということは，オバーラップしてることを意味する
			if (j == -1 /*undefined*/)
				throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "index is not defined"));
			points_consists_faces[i] = (i > j) ? points_consists_faces[j] : new networkPoint(this, this, {v_IN[i][0], v_IN[i][1], v_IN[i][2]});
			i++;
		}
		Print("setting points done", Red);
		return points_consists_faces;
	};
	//-------------------------------------
private:
	void setFaces(const VV_i &f_v, const V_netPp &points_consists_faces)
	{
		try
		{
			Print("setting faces...", Red);
			// connect points by lines and create faces, then store those faces in lines
			int count = 0;
			for (const auto &index : f_v)
			{
				if (count++ % 100 == 0)
					std::cout << "|";
				//      Print(index,Red);
				/*{0,1,4}*/
				// V_netPp Ps = points_consists_faces(index);
				// networkLine *line;
				int f_size = index.size();
				if (f_size != 3)
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "f_size != 3");

				if (points_consists_faces[index[0]] && points_consists_faces[index[1]] && points_consists_faces[index[2]])
				{
					V_netLp Ls = {link(points_consists_faces[index[0]], points_consists_faces[index[1]], this),
								  link(points_consists_faces[index[1]], points_consists_faces[index[2]], this),
								  link(points_consists_faces[index[2]], points_consists_faces[index[0]], this)};
					new networkFace(this, this, Ls);
				}
				else
				{
					std::stringstream ss;
					ss << "index = " << index << std::endl;
					// ss << "getPoints(index) = " << Ps << std::endl;
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
				}
			}
			std::cout << std::endl;
			Print("setting faces done", Red);
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};

public:
	V_netPp getXPoints()
	{
		V_netPp ret({});
		for (const auto &line : this->getLines())
			for (const auto &info : line->XPoints)
				ret.emplace_back(info);
		return ret;
	};
	void clearXPoints()
	{
		for (const auto &p : this->Points)
			for (const auto &line : p->networkObject::getLines())
				line->XPoints.clear();

		for (const auto &face : this->Faces)
		{
			face->clearXPoints();
			for (const auto &line : face->networkObject::getLines())
				line->XPoints.clear();
		}
	};
	/* ------------------------------------------------------ */

	void DeleteParticles()
	{
		auto points = this->Points;
		for (const auto &p : points)
			if (std::get<0>(p->particlize_info))
				p->Delete();
	};
};
/*Network_code*/

///////////////////////////////////////////////////////////////////////

VV_netPp extractPoints(const V_netFp &Fs)
{
	VV_netPp ret;
	for (const auto &f : Fs)
		ret.emplace_back(f->getPoints());
	return ret;
};
VV_netPp extractPoints(const std::unordered_set<networkFace *> &faces)
{
	VV_netPp ret(faces.size(), V_netPp(3));
	int i = 0;
	for (const auto &f : faces)
		ret[i++] = f->getPoints();
	return ret;
};
VV_netPp extPoints(const V_netFp &Fs)
{
	VV_netPp ret;
	for (const auto &f : Fs)
		ret.emplace_back(f->getPoints());
	return ret;
};
#include "my_vtk.hpp"
#include "networkPoint.hpp"

// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM
#ifdef BEM
V_Netp takeNetworks(const V_netPp &ps)
{
	V_Netp ret({});
	for (auto &p : ps)
		network::add(ret, p->getNetwork());
	return ret;
};
////////////////
template <typename T>
std::vector<T *> takeB(const std::vector<T *> &nets)
{
	std::vector<T *> ret({});
	for (const auto &n : nets)
		if (n->isB())
			ret.emplace_back(n);
	return DeleteDuplicates(ret);
};
template <typename T>
std::vector<T *> takeC(const std::vector<T *> &nets)
{
	std::vector<T *> ret({});
	for (const auto &n : nets)
		if (n->isC())
			ret.emplace_back(n);
	return DeleteDuplicates(ret);
};
template <typename T>
std::vector<T *> takeN(const std::vector<T *> &nets)
{
	std::vector<T *> ret({});
	for (const auto &n : nets)
		if (n->isN())
			ret.emplace_back(n);
	return DeleteDuplicates(ret);
};
template <typename T>
std::vector<T *> takeD(const std::vector<T *> &nets)
{
	std::vector<T *> ret({});
	for (const auto &n : nets)
		if (n->isD())
			ret.emplace_back(n);
	return DeleteDuplicates(ret);
};
/////////////////////////////
template <typename T>
inline bool networkObject<T>::isB() const { return this->getNetwork()->isB(); };
template <typename T>
inline bool networkObject<T>::isC() const { return this->getNetwork()->isC(); };
template <typename T>
inline bool networkObject<T>::isD() const { return this->getNetwork()->isD(); };
template <typename T>
inline bool networkObject<T>::isN() const { return this->getNetwork()->isN(); };
//
inline V_d networkPoint::nabla_phi() const { return this->getNetwork()->nabla_phi; };
inline double networkPoint::phi_n()
{
	try
	{
		return Dot(this->getNormalFromSameBC(), this->getNetwork()->nabla_phi);
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
/////////////////////////////
inline V_d networkPoint::getNormalFromSameBC()
{
	VV_d normals({});
	V_netFp fs;
	V_netPp ps;
	for (auto d = 1; d <= 4; d++)
	{
		//------------------
		if (d == 1)
			fs = this->getFaces(); //これだけならgetNormal()と同じ
		else
		{
			depth_searcher<networkPoint> S(d);
			S.set(this);
			S.search();
			ps = S.getObjects();
			fs.clear();
			for (const auto &p : ps)
				for (const auto &f : p->getFaces())
					fs.emplace_back(f);
			fs = DeleteDuplicates(fs);
		}
		//------------------
		if (this->isD())
		{
			for (const auto &f : takeD(fs))
				normals.emplace_back(f->getNormal());
		}
		else
		{
			for (const auto &f : takeN(fs))
				normals.emplace_back(f->getNormal());
		}
		//------------------
		if (!normals.empty())
			return Mean(normals) / Norm(Mean(normals));
	}
	std::stringstream ss;
	ss << "Name = " << this->getNetwork()->getName() << "\n"
	   << "this->getX() = " << this->getX() << "\n"
	   << ", isD = " << this->isD() << "\n"
	   << ", isN = " << this->isN() << "\n"
	   << ", getName() = " << this->getNetwork()->getName()
	   << "\n";
	for (auto &q : ps)
		ss << "q->getX() = " << q->getX() << "\n"
		   << ", isD = " << q->isD() << "\n"
		   << ", isN = " << q->isN() << "\n"
		   << ", getName() = " << q->getNetwork()->getName()
		   << "\n";
	for (auto &f : fs)
		ss << "f->getX() = " << f->getX() << "\n"
		   << ", normal = " << f->getNormal() << "\n"
		   << ", isD = " << f->isD() << "\n"
		   << ", isN = " << f->isN() << "\n"
		   << ", getName() = " << f->getNetwork()->getName()
		   << "\n";

	// mk_vtu("./vtu/new_point.vtu", {{this}});

	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "no face is found\n" + ss.str());
};

#endif
// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM

////////////////////////////////////////////////////////////////////////
inline V_netFp networkLine::getFaces(const bool TorF) const
{
	V_netFp ret;
	for (const auto &f : this->getFaces())
		if (f->networkObject::isStatus(TorF))
			ret.emplace_back(f);
	return ret;
};
////////////////////////////////////////////////////////
int intersectQ(networkFace *f, networkLine *l)
{
	return network::find(f->getLinesPenetrate(), l);
};

inline V_netPp networkFace::getXPoints() const
{ //特別な場合のみしか使えないので改善が必要
	bool TorF = true;
	network::setStatus(this->network->getLines(), !TorF);
	auto Ps = getPointsPenetrate();
	V_netPp object{Ps[0]}, nextobject, ret{Ps[0]};
	netP *P;
	while (true)
	{
		for (const auto &p : object)
			for (const auto &l : p->networkObject::getLines_toggle(!TorF))
			{
				if ((P = (*l)(p))->networkObject::getStatus() != TorF && P->getXFace() == this)
				{
					ret.emplace_back(P);
					P->networkObject::setStatus(TorF);
					nextobject.emplace_back(P);
				}
			}
		if (nextobject.empty())
			break;
		object = nextobject;
		nextobject.clear();
	}
	ret.emplace_back(Ps[1]);
	return ret;
};
//@ ------------------------------------------------------ */
//@                          外部関数                       */
//@ ------------------------------------------------------ */
///////////////////////
/*networkFace::Delete_detail
##### networkFace::Delete()
1. `this->Lines`から`this`を削除する
2. `this->Lines`の内，面を全て失った線は，点の接続を解除する
> **NOTE:**
`networkLine`から`networkFace`からへの一方的なポイントは，絶対に起きないようにしているつもりである．
**`networkFace`のコンストラクタでは，自動で`networkLine`に`networkFace`を格納するようにしている．**
なので，`networkFace`を削除する際に，空になった`networkLine`は削除して構わない，
これが可能なのは，一方向のみのポイントが存在しないからである．

networkFace::Delete_detail*/
/*networkFace::Delete_code*/
// #define debug_Delete

/*!
networkに上に穴が空いている状態もあり得るため，
faceのDelete()は，pointのDelete()を伴わない．
*/
inline void networkFace::Delete()
{
	if (!this->storage->isMember(this))
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "is not member");
	try
	{
#ifdef debug_Delete
		std::cout << "networkFace::Delete() ... ";
		// Print("Delete()", Red);
#endif
		auto tmpLines = this->Lines;
		//面の場合は，線で繋がったもう片方の面は残したままでもいいだろう
#ifdef debug_Delete
		std::cout << "1.\n"
					 "+---------+        +----------+\n"
					 "|       --|------->|   Line   |\n"
					 "|  Face   |<-------|--        |\n"
					 "+---------+        +----------+\n"
				  << std::endl;
#endif
		for (const auto &l : tmpLines)
			l->Erase(this);
#ifdef debug_Delete
		std::cout << "+---------+        +----------+\n"
					 "|       --|------->|   Line   |\n"
					 "|  Face   |        |          |\n"
					 "+---------+        +----------+\n"
				  << std::endl;
#endif
		for (const auto &l : tmpLines)
			this->Erase(l);
#ifdef debug_Delete
		std::cout << "+---------+        +----------+\n"
					 "|         |        |   Line   |\n"
					 "|  Face   |        |          |\n"
					 "+---------+        +----------+\n"
				  << std::endl;
#endif

		///////////////////////

#ifdef debug_Delete
		Print("p->resetXinfo()");
#endif
		for (auto &p : this->XPoints)
			p->resetXinfo();
			//この面が保存されている可能性があるネットワークを全てチェックし，消去

#ifdef debug_Delete
		Print("この面が保存されている可能性があるネットワークを全てチェックし，消去");
#endif
		this->storage->Erase(this);
		// if (this->storage->Erase(this))
		// 	return true;
		// else
		// {
		// 	std::stringstream ss;
		// 	ss << "この面は，sotrageに含まれていない" << std::endl;
		// 	ss << "network name : " << this->network->getName() << std::endl;
		// 	ss << "storage name : " << this->storage->getName() << std::endl;
		// 	if (this->network->isMember(this))
		// 		ss << this->network->getName() << "に含まれている" << std::endl;
		// 	else
		// 		ss << this->network->getName() << "にも含まれていない" << std::endl;

		// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
		// }
#ifdef debug_Delete
		Print("ParametricPointsを消去");
#endif
		auto pp = this->getParametricPoints();
		for (const auto &p : pp)
			delete p;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
/*networkFace::Delete_code*/
////////////////////////////////////////////////////////////////////
/*!
pointのDelete()と同時に，pointが接続している辺のDelete()を必ず伴うべきだ．
*/
inline void networkPoint::Delete()
{
	//
	try
	{
		if (!this->storage->isMember(this))
		{
			// std::cout << Red << this << " is already deleted" << reset << std::endl;
			// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		}
		else
		{
			auto ps = this->getNeighbors();
			V_netLp delL = {};
			netLp line;
			for (const auto &p : ps)
				if ((line = unlink(p, this)))
					delL.emplace_back(line);

			//これで，線から参照されることはない
			// l->Delete()しても，この点は，参照されない．
			for (auto i = 0; i < delL.size(); i++)
				if (delL[i]->getFaces().empty())
				{
					for (const auto &p : delL[i]->getPoints())
						delL[i]->Delete(); // disconnect points from the line
					auto tmp = delL[i];
					delL[i] = nullptr;
					delete tmp;
				}

			this->resetXinfo();
			this->storage->Erase(this);
			// if (this->storage->Erase(this))
			// 	return true;
			// else
			// {
			// 	std::stringstream ss;
			// 	ss << "この点は，sotrageに含まれていない" << std::endl;
			// 	ss << "network name : " << this->network->getName() << std::endl;
			// 	ss << "storage name : " << this->storage->getName() << std::endl;
			// 	if (this->network->isMember(this))
			// 		ss << this->network->getName() << "に含まれている" << std::endl;
			// 	else
			// 		ss << this->network->getName() << "にも含まれていない" << std::endl;

			// 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
			// }
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
///////////////////////////////////////////////////////////////
/*!
lineのDelete()と同時に，lineが接続している面のDelete()を必ず伴うべきだ．
*/
inline void networkLine::Delete()
{
	try
	{
		int maxloop = 1000 + this->Faces.size();
		int counter = 0;
		while (!this->Faces.empty() && counter++ < maxloop)
		{
			auto tmp = this->Faces;
			tmp[0]->Delete();
		}

		if (this->Point_A != nullptr)
		{
			(this->Point_A)->Erase(this);
			this->Point_A = nullptr;
		}
		if (this->Point_B != nullptr)
		{
			(this->Point_B)->Erase(this);
			this->Point_B = nullptr;
		}

		auto tmpP = this->getXPoints();
		for (auto &p : tmpP)
		{
			p->resetXinfo();
			network::erase(this->XPoints, p);
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	};
};
///////////////////////////////////////////////////////////////
inline networkPoint::~networkPoint() { this->Delete(); };
inline networkFace::~networkFace() { this->Delete(); };
inline networkLine::~networkLine() { this->Delete(); };
///////////////////////////////////////////////////////////////
inline Network::~Network()
{
	try
	{
		Print(this->getName(), red);
		try
		{
			int maxloop = 100000 + this->Faces.size();
			int counter = 0;
			while (!this->Faces.empty() && counter++ < maxloop)
			{
				auto tmp = this->Faces;
				(*tmp.begin())->Delete();
			}
			if (counter > maxloop)
			{
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			}
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};

		std::cout << "|        |"
				  << "Faces" << std::endl;
		std::cout << "|        |" << this->Faces.size() << std::endl;
		std::cout << "|________|" << std::endl;

		//--------------------
		try
		{
			int maxloop = 100000 + this->Points.size();
			int counter = 0;
			while (!this->Points.empty() && counter++ < maxloop)
			{
				//
				// auto tmp = this->Points;
				// tmp[0]->Delete();
				//
				// this->Points[0]->Delete();
				(*this->Points.begin())->Delete();
			}
			if (counter > maxloop)
			{
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			}
		}
		catch (const error_message &e)
		{
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		}

		std::cout << "|        |"
				  << "Points" << std::endl;
		std::cout << "|        |" << this->Points.size() << std::endl;
		std::cout << "|________|" << std::endl;
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << reset << std::endl;
	};
};
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
template <class net = Network>
class NetworkX : public net
{
public:
	// NetworkX(const V_i &N, const V_d &scale, const int gridtype, const std::string &name_IN = "no_name")
	//     : net(name_IN)
	// {
	//   int Nx = N[0], Ny = N[1];
	//   VV_d v;
	//   v.reserve((2 * Nx) * (2 * Ny) * 9);
	//   VV_i f_v;

	//   setGrid(gridtype, v, f_v, Nx, Ny);

	//   double scalex = scale[0], scaley = scale[1], x, y, z = scale[2];
	//   for (auto &w : v)
	//   {
	//     w[0] = w[0] * scalex / (2. * Nx + 1);
	//     w[1] = w[1] * scaley / (2. * Ny + 1);
	//     //      w[2] = w[2] + .2*cos(3*M_PI*w[0]/scalex)*sin(4*M_PI*w[1]/scaley);
	//     w[2] = w[2] + .1 * pow(cosh(50. * M_PI * (w[0] * w[0] / scalex + w[1] * w[1] / scaley)), -2.);
	//   };

	//   net::generateNetwork(v, f_v);
	// };

	NetworkX(const V_i &N, const V_d &scale, const int gridtype,
			 const std::function<double(const V_d &)> waveheight,
			 const std::string &name_IN = "no_name") : net(name_IN)
	{
		int Nx = N[0], Ny = N[1];
		VV_d v;
		v.reserve((2 * Nx) * (2 * Ny) * 9);
		VV_i f_v;

		setGrid(gridtype, v, f_v, Nx, Ny);

		double scalex = scale[0], scaley = scale[1], x, y, z = scale[2];
		V_d xy;
		double pn = RandomReal({-1., 1.}) > 0 ? 1. : -1;
		double rand = pn * RandomReal({0.5, 1.}) * 1E-20;
		for (auto &w : v)
		{
			w[0] = w[0] * scalex / (2. * Nx + 1.) + scalex * rand;
			w[1] = w[1] * scaley / (2. * Ny + 1.) + scaley * rand;
			//      w[2] = w[2] + .2*cos(3*M_PI*w[0]/scalex)*sin(4*M_PI*w[1]/scaley);
			w[2] = w[2] + waveheight({w[0], w[1]});
		};

		net::generateNetwork(v, f_v);
	};

	void setGrid(const int gridtype, VV_d &v, VV_i &f_v, const int Nx, const int Ny)
	{
		int k = 0, c = 0;
		double x, y, z = 0;
		switch (gridtype)
		{
		case 0:
		{
			for (auto i = -Nx; i <= Nx; i++)
			{
				x = double(i);
				for (auto j = -Ny; j <= Ny; j++)
				{
					y = double(j);
					v.push_back({x, y, z});
					v.push_back({-.5 + x, -.5 + y, z});
					v.push_back({.5 + x, -.5 + y, z});
					f_v.push_back({c = k++, k++, k++});
					//
					v.push_back({.5 + x, -.5 + y, z});
					v.push_back({.5 + x, .5 + y, z});
					f_v.push_back({c, k++, k++});
					//
					v.push_back({.5 + x, .5 + y, z});
					v.push_back({-.5 + x, .5 + y, z});
					f_v.push_back({c, k++, k++});
					//
					v.push_back({-.5 + x, .5 + y, z});
					v.push_back({-.5 + x, -.5 + y, z});
					f_v.push_back({c, k++, k++});
				}
			};
			break;
		}
		case 1:
		{
			for (auto i = -Nx; i <= Nx; i++)
			{
				x = double(i);
				for (auto j = -Ny; j <= Ny; j++)
				{
					y = double(j);
					v.push_back({.5 + x, .5 + y, z});
					v.push_back({-.5 + x, -.5 + y, z});
					v.push_back({.5 + x, -.5 + y, z});
					f_v.push_back({c = k++, k++, k++});
					//
					v.push_back({-.5 + x, .5 + y, z});
					v.push_back({-.5 + x, -.5 + y, z});
					f_v.push_back({c, k++, k++});
				}
			};
			break;
		}
		case 2:
		{
			double X, Y;
			for (auto i = -Nx; i <= Nx; i++)
			{
				x = double(i);
				for (auto j = -Ny; j <= Ny; j++)
				{
					y = double(j);

					X = .5 + x;
					Y = .5 + y;
					v.push_back({(X - Y) / sqrt(2), (X + Y) / sqrt(6), z});

					X = -.5 + x;
					Y = -.5 + y;
					v.push_back({(X - Y) / sqrt(2), (X + Y) / sqrt(6), z});

					X = .5 + x;
					Y = -.5 + y;
					v.push_back({(X - Y) / sqrt(2), (X + Y) / sqrt(6), z});

					f_v.push_back({c = k++, k++, k++});
					//

					X = -.5 + x;
					Y = .5 + y;
					v.push_back({(X - Y) / sqrt(2), (X + Y) / sqrt(6), z});

					X = -.5 + x;
					Y = -.5 + y;
					v.push_back({(X - Y) / sqrt(2), (X + Y) / sqrt(6), z});

					f_v.push_back({c, k++, k++});
				}
			};
			break;
		}
		}
	};
};
///////////////////////////////////////////////////////////////////////
inline void Network::setXPoints(Network &target, Network *interactionNet)
{
	Print("このオブジェクトの範囲外の面は，比較対象に含めない．線がcrossできる面だけを抜き取る", Red);
	target.setBounds();
	this->setBounds();

	Print("obj3D::takeIfBoundariesOverlap", Red);
	auto targetFs = obj3D::takeIfBoundariesOverlap(target.Faces, this->getBounds());
	std::vector<Tddd> locs(0);
	for (const auto &tfs : targetFs)
		for (const auto &loc : tfs->getLocations())
			locs.emplace_back(Tddd{loc[0], loc[1], loc[2]});

	V_netPp points_check;
	if (locs.empty())
	{
		std::cout << target.getName() << "は，";
		std::cout << this->getName() << "と干渉しないので，";
		std::cout << interactionNet->getName() << "は空です．" << std::endl;
	}
	else
	{
		Print("この範囲に関連する線だけを抜き出す", Red);
		// Print("getLinesInBounds(boundsOfTarget)は休止",red);facesのぼboundsといわずそもそもtargetのbounds内にLineふくまれていないようだ
		// Print(target.getBounds(),Cyan);

		// Print(CoordinateBounds(locs), Red);
		auto Ls = DeleteDuplicates(obj3D::takeIfBoundariesOverlap(this->getLines(), geometry::CoordinateBounds(locs)));
		// Print(Ls.size(), Red);

		for (auto &F : targetFs)
		{
			VV_d fs_locs = F->getLocations();
			for (const auto &line : Ls)
			{
				// //not necessary....
				// bool AlreadeyHaveInfoOfThisFace = false;
				// for(const auto& xp:line->getXPoints())
				//   if(xp->getXFace() == targetFs[i])
				//     AlreadeyHaveInfoOfThisFace = true;

				if (isIntersectingSurface(fs_locs, line->getLocations()) == 3)
				{
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						auto xpoint = new networkPoint(interactionNet /*属性*/,
													   interactionNet /*storage*/,
													   pOnSurfaceTuple(F->getLocationsTuple(), line->getLocationsTuple()),
													   line /*自動的にこのlineのxpointにこの点は保存される*/,
													   F /*自動的にこのfaceのxpointにこの点は保存される*/);
						points_check.emplace_back(xpoint);
					}
				}
			}
		}

		// #ifdef _OPENMP
		// 		Print("並列計算", Red);
		// #pragma omp parallel for
		// #endif
		// 		for (auto i = 0; i < targetFs.size(); i++)
		// 		{
		// 			VV_d fs_locs = targetFs[i]->getLocations();
		// 			for (const auto &line : Ls)
		// 			{
		// 				// //not necessary....
		// 				// bool AlreadeyHaveInfoOfThisFace = false;
		// 				// for(const auto& xp:line->getXPoints())
		// 				//   if(xp->getXFace() == targetFs[i])
		// 				//     AlreadeyHaveInfoOfThisFace = true;

		// 				if (isIntersectingSurface(fs_locs, line->getLocations()) == 3)
		// 				{
		// #ifdef _OPENMP
		// #pragma omp critical
		// #endif
		// 					{
		// 						auto xpoint = new networkPoint(interactionNet /*属性*/,
		// 													   interactionNet /*storage*/,
		// 													   pOnSurfaceTuple(targetFs[i]->getLocationsTuple(), line->getLocationsTuple()),
		// 													   line /*自動的にこのlineのxpointにこの点は保存される*/,
		// 													   targetFs[i] /*自動的にこのfaceのxpointにこの点は保存される*/);
		// 						points_check.emplace_back(xpoint);
		// 					}
		// 				}
		// 			}
		// 		}
	}

#if defined(simulation)
	Print(message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
	mk_vtu("./vtu/xpoint" + this->getName() + "_" + target.getName() + std::to_string(time_step) + ".vtu", {points_check});
	mk_vtu("./vtu/xpoint_checked_face" + this->getName() + "_" + target.getName() + std::to_string(time_step) + ".vtu", targetFs);
#endif

	// std::cin.ignore();
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*XNetwork_detail
`Network(Network& netA, Network& netB)`は，

1. netAとnetB間の`xPoint`を`setXPoints`を使って設定する．
2. netAのFaceがpenetrating line（`xpoints`）を持つ場合，そのlineが同じFaceを貫通しているなら`xpoints`を`link`する．
3. （鎖の状態）FaceとFaceが鎖の状態にある場合，それらFaceを貫くpenetrating lineを探し出し，そのline上`xpoints`を`link`する．

> **NOTE:** `water.setXPoints(obj,this);`としているため，ここで作成される`xPoint`の持つ`network`ポインタは，この作成者である`this`干渉ネットワークを指す．

XNetwork_detail*/
/*XNetwork_code*/
// inline Network::Network(Network &water, Network &obj, const std::string &name_IN = "no_name_xnetwork") : object3D(), name(name_IN)
// {
// #ifdef BEM
//   this->setCornerNeumann();
// #endif
//   // obj.clearXPoints();
//   // water.clearXPoints();
//   Print("Xpointをclear", _Blue);
//   //オブジェクトが交差している点を計算，保存
//   water.setXPoints(obj, this);
//   obj.setXPoints(water, this);
//   Print("Xpointをset", _Blue);
//   VV_d vertices;
//   V_netPp accumXPoints;
//   bool findsameface = false;
//   //----------
//   Print("water.Faces " + std::to_string(water.Faces.size()), _Blue);
//   Print("water.Points " + std::to_string(water.Points.size()), _Blue);
//   Print("water.Lines " + std::to_string(water.getLines().size()), _Blue);
//   int n = 0, c = 0;
//   for (auto l : water.getLines())
//     if (l->getNetwork() == &water)
//       c++;
//     else
//       n++;

//   Print(c);
//   Print(n);
//   Print("water.PointGarbage " + std::to_string(water.PointGarbage.size()), _Blue);
//   Print("wate.FaceGarbage " + std::to_string(water.FaceGarbage.size()), _Blue);

//   Print("obj.Faces " + std::to_string(obj.Faces.size()), _Blue);
//   Print("obj.Points" + std::to_string(obj.Points.size()), _Blue);
//   Print("obj.Lines " + std::to_string(obj.getLines().size()), _Blue);
//   Print("obj.PointGarbage " + std::to_string(obj.PointGarbage.size()), _Blue);
//   Print("obj.FaceGarbage " + std::to_string(obj.FaceGarbage.size()), _Blue);

//   /// CHECK IF INTERSECTION EXISTS
//   if (linkXPoints(water, obj).empty() /*no intersection*/)
//   {
//     message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Intersection Network has been generated but no intersection point exists");
//   }
//   else
//   {
//     this->Points = DeleteDuplicates(accumXPoints);
//     setBounds();
//     Print("Intersection Network has been generated", Magenta);
//   }
//   displayStates();
// };
////////////////////////////////////////////////////////////////////////
// Networkコンストラクタ2
inline Network::Network(const V_Netp &base_nets, const std::string &name_IN = "no_name_xnetwork")
	: object3D(geometry::CoordinateBounds({0, 0, 0})),
	  name(name_IN),
	  /* ------------------------------------------------------ */
	  force({0., 0., 0., 0., 0., 0.}),
	  inertia({0., 0., 0., 0., 0., 0.}),
	  forced_acceleration({0., 0., 0., 0., 0., 0.}),
	  forced_velocity({0., 0., 0., 0., 0., 0.}),
	  acceleration({0., 0., 0., 0., 0., 0.}),
	  velocity({0., 0., 0., 0., 0., 0.}),
	  mass(0.),
	  center_of_mass({0., 0., 0.}),
	  BucketFaces(geometry::CoordinateBounds({0., 0., 0.}), 1.),
	  BucketPoints(geometry::CoordinateBounds({0., 0., 0.}), 1.),
	  BucketParametricPoints(geometry::CoordinateBounds({0., 0., 0.}), 1.)
{
#ifdef BEM
	this->setCornerNeumann();
#endif
	/* ----------------------------------- */
	// 物理
	// this->force = V_d(3, 0.);
	/* ------------------------------------ */
	//オブジェクトが交差している点を計算，保存
	for (auto i = 0; i < base_nets.size(); i++)
		for (auto j = 0; j < base_nets.size(); j++)
			if (i != j)
			{
				Print(base_nets[i]->getName() + base_nets[j]->getName());
				base_nets[i]->setXPoints(*base_nets[j], this);
			}

	V_netPp accumXPoints({});
	for (auto i = 0; i < base_nets.size(); i++)
		for (auto j = i + 1; j < base_nets.size(); j++)
		{
			auto tmp = linkXPoints(*base_nets[i], *base_nets[j]);
			accumXPoints.insert(accumXPoints.end(), tmp.begin(), tmp.end());
		}

	/// CHECK IF INTERSECTION EXISTS
	if (accumXPoints.empty() /*no intersection*/)
	{
		message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Intersection Network has been generated but no intersection point exists");
	}
	else
	{
		this->add(accumXPoints);
		setBounds();
		Print("Intersection Network has been generated", Magenta);
	}
	displayStates();
};
////////////////////////////////////////////////////////////////////////
inline V_netPp Network::linkXPoints(Network &water, Network &obj)
{
	V_netPp accumXPoints;
	std::cout << __PRETTY_FUNCTION__ << "Xpointをset " << water.getName() << " " << obj.getName() << std::endl;
	bool findsameface = false;
	std::vector<Network *> nets({&water, &obj});
	for (const auto &net : nets)
		for (const auto &f : net->Faces)
		{ // thisの干渉点だけを対象にする！！！
			auto xps = f->getPointsPenetrate(this /*干渉点の制限！*/);
			for (size_t i = 0; i < xps.size(); i++, findsameface = false) //あるcrossinfoを持つLineにたいして
			{
				for (size_t j = i + 1; j < xps.size(); j++)
					if (xps[i]->getXFace() == xps[j]->getXFace())
					{
						//２本のLineが，同じ面を貫通しているLineである場合
						auto newline = link(xps[i], xps[j], this);
						// newline->setIntxn(true);
						accumXPoints.emplace_back(xps[i]);
						accumXPoints.emplace_back(xps[j]);
						findsameface = true;
						//                break;
					}
				//２本のLineが，異なる面を貫通しているLineである場合
				//それぞれが貫通している面のcrossinfoを調べる．その中に元の面を貫通するLineがあった場合，元の貫通点と接続する
				if (!findsameface)
					for (const auto &XPOfL : xps)
						for (const auto &XPOfL2 : XPOfL->getXFace()->getPointsPenetrate(this /*干渉点の制限！*/))
							if (XPOfL2->getXFace() == f)
							{
								auto newline = link(XPOfL, XPOfL2, this);
								// newline->setIntxn(true);
								accumXPoints.emplace_back(XPOfL);
								accumXPoints.emplace_back(XPOfL2);
								//                    break;
							}
			}
		}
	for (const auto &f : obj.Faces) // thisの干渉点だけを対象にする！！！
	{
		auto xps = f->getPointsPenetrate(this /*干渉点の制限！*/);
		for (size_t i = 0; i < xps.size(); i++, findsameface = false)
			for (size_t j = i + 1; j < xps.size(); j++)
				if (xps[i]->getXFace() == xps[j]->getXFace())
				{
					auto newline = link(xps[i], xps[j], this);
					// newline->setIntxn(true);
					accumXPoints.emplace_back(xps[i]);
					accumXPoints.emplace_back(xps[j]);
					findsameface = true;
					// break;
				}
	}
	return accumXPoints;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*XNetwork_code*/

/*networkFace::divide_detail
divideの操作を整理して理解するために，networkFaceにdivideの一部操作を任せている．
ここでnetworkFaceは，与えられた線と自身がすでに所持していた線を，適切に接続し直す．
点と線の接続はnetworkFaceに行わせず，netowrkLineに行わせる．
networkFace::divide_detail*/
//
/*networkFace::divide_code*/
inline netFp networkFace::divide(netLp DivL /*this*/, netLp newDivL, netLp newMidL, int type)
{
	if (!MemberQ(this->Lines, DivL))
		throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "divL is not found in this->Lines"));
	auto oldF = this;
	auto newF = new networkFace(this);
	auto fL = oldF->getLineFront(DivL);
	auto bL = oldF->getLineBack(DivL);
	//          /     \                     /     \
  //         /       \                   /       \
  //        /         \                 /         \
  //       /    / \    \               /    / \    \
  // bL<--/--- /   \----\-->fL   bL<--/----/   \----\-->fL
	//   --/--->/oldF \<---\--         /    /newF \    \
  //    /     --|A---     \         /     --|----     \
  //   /        V|         \       /        V          \
  //    -------DivL-------         --------DivL---------
	//
	if (oldF->Switch(fL, newMidL))
	{ // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
		std::stringstream ss;
		ss << "oldF->getLines() = " << oldF->networkObject::getLines();
		ss << ", oldF = " << oldF << ", fL = " << fL << ", newMidL = " << newMidL;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
	}

	// oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
	if (newF->Switch(bL, newMidL))
	{
		std::stringstream ss;
		ss << "newF = " << newF << ", bL = " << bL << ", newMidL = " << newMidL;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
	}
	// oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
	if (fL->Switch(oldF, newF))
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

	//          /     \                     /     \
	//         /       \                   /       \
	//        /         \                 /         \
	//       /    / \    \               /    / \    \
	// bL<--/----/   \<---\-- newMidL<--/--  /   \----\-->fL
	//   --/--->/oldF \----\->         /    /newF \<---\---
	//    /     --|A---     \         /     --|----     \
	//   /        V|         \       /        V          \
	//    -------DivL-------         --------DivL---------
	newMidL->set(oldF, newF);
	//          /     \                     /     \
	//         /       \                   /       \
	//        /         \                 /         \
	//       /    / \    \               /    / \    \
	// bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
	//  ---/--->/oldF \<---\--       --/--->/newF \<---\---
	//    /     --|A---     \         /     --|----     \
	//   /        V|         \       /        V          \
	//    -------DivL-------         --------DivL---------
	if (type == 1)
	{
		// oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
		if (newF->Switch(DivL, newDivL))
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		newDivL->set(newF);
		//          /     \                     /     \
		//         /       \                   /       \
		//        /         \                 /         \
		//       /    / \    \               /    / \    \
		// bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
		//  ---/--->/oldF \<---\--       --/--->/newF \<---\---
		//    /     --|A---     \         /     --|A---     \
		//   /        V|         \       /        V|         \
		//    -------DivL/*this*/---     -------newDivL------
	}
	else if (type == 2)
	{
		// oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
		if (oldF->Replace(DivL, newDivL))
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		DivL->Add(newF);
		//          /     \                     /     \
		//         /       \                   /       \
		//        /         \                 /         \
		//       /    / \    \               /    / \    \
		// bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
		//  ---/--->/oldF \<---\--       --/--->/newF \<---\---
		//    /     --|A---     \         /     --|A---     \
		//   /        V|         \       /        V|         \
		//    ------newDivL------         ----DivL/*this*/----
	}
	return newF;
};
/*networkFace::divide_code*/
//
/*networkLine::divide_detail
divideは，**新しい点と面**を作成し，このnetwrokLineが所属する��ットワークに保存する．
もし，このネットワークに**新しい点と面**が保存されてはまずい場合は．問題が生じるだろう．
networkLine::divide_detail*/
//
/*networkLine::divide_code*/
///////////////////////////////////////////////////////////
#ifdef debug_divide
int _I_ = 0;
#endif

#include "networkLine.hpp"

#ifdef GNUPLOT_H
// Points
VV_d getData(V_netPp points)
{
	VV_d vec;
	for (const auto &p : points)
		vec.emplace_back(p->getX());
	return vec;
};

VVV_d getVectorData(V_netPp ps)
{
	VVV_d VVV;
	int s = ps.size();
	// for (auto i = 0; i < s - 1; i++)
	// 	VVV.push_back({ps[i]->xyz, ps[i + 1]->xyz - ps[i]->xyz});
	for (auto i = 0; i < s - 1; i++)
		VVV.push_back({ps[i]->getX(), ps[i + 1]->getX() - ps[i]->getX()});

	return VVV;
};
VVV_d getVectorData(VV_netPp PS)
{
	VVV_d VVV;
	int s;
	for (const auto &ps : PS)
	{
		s = ps.size();
		for (auto i = 0; i < s; i++)
			VVV.push_back({ps[i]->getX(), ps[(i + 1) % s]->getX() - ps[i]->getX()});
	}
	return VVV;
};
// Lines
void pushVectorData(VVV_d &vec, networkLine *line)
{
	auto xyz = line->getLocations();
	vec.push_back({{xyz[0], xyz[1] - xyz[0]}});
	vec.push_back({{xyz[1], xyz[0] - xyz[1]}});
};
void pushVectorData(VVV_d &vec, std::vector<networkLine *> lines)
{
	for (const auto &l : lines)
		pushVectorData(vec, l);
};
VVV_d getVectorData(std::vector<networkLine *> lines)
{
	VVV_d vec;
	pushVectorData(vec, lines);
	return vec;
};
// Faces
void pushVectorData(VVV_d &vec, networkFace *face)
{
	auto xyz = face->getLocations();
	vec.push_back({{xyz[0], xyz[1] - xyz[0]}});
	vec.push_back({{xyz[1], xyz[2] - xyz[1]}});
	vec.push_back({{xyz[2], xyz[0] - xyz[2]}});
};
void pushVectorData(VVV_d &vec, V_netFp faces)
{
	for (const auto &f : faces)
		pushVectorData(vec, f);
};
VVV_d getVectorData(const networkFace *face)
{
	return getVectorData(face->networkObject::getLines());
};
VVV_d getVectorData(const V_netFp &faces)
{
	VVV_d vec;
	pushVectorData(vec, faces);
	return vec;
};
VV_d getData(networkFace *f)
{
	return getData(f->getPoints());
};
VVV_d getData(const V_netFp &fs)
{
	VVV_d ret;
	for (const auto &f : fs)
		ret.push_back(getData(f));
	return ret;
};
#endif

#include "networkFace.hpp"

/*selector_detail
`selector`は，`searcher`を使って，ある基点から干渉ネットワークに至るまでの点や線を取得し，
後に取得した点を`innerP`, `cornerP`, `outerP`に選別する．
（`selector`はファンクターである．ファンクターはクラスなので，関数の結果をメンバ変数として保持できるため便利）
selector_detail*/
/*selector_code*/

class selector
{
	using V_Netp = std::vector<Network *>;

private:
	V_Netp XNets;
	int type;

public:
	V_netPp innerP, cornerP, outerP;
	V_netFp faces;
	searcher9<netP> *S;

	~selector()
	{
		delete (S);
	};
	selector(const std::vector<Network *> &XNets_IN, int type_IN = 1) : type(type_IN), XNets(XNets_IN), innerP({}), cornerP({}), outerP({}), faces({}), S(nullptr)
	{
		Print("selectorを生成, type=" + std::to_string(type), Red);
	};
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	void prepare(Network *target, const V_d &cardinal)
	{
		std::cout << __PRETTY_FUNCTION__ << std::endl;
		std::cout << __PRETTY_FUNCTION__ << std::endl;
		//------------------
		network::setStatus(target->getPoints(), false);
		for (const auto &net : this->XNets)
			network::setStatus(net->getPoints(), false);
		//----- serach -----
		Print("new searcher9<netP>");
		this->S = new searcher9<netP>;
		S->set(network::takeNearest(target->getPoints(), cardinal));
		for (const auto &net : this->XNets)
			S->addNetwork(net);

		try
		{
			S->search(false); // must be false!
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};

		// searcher9の特性
		Print("選別!");
		this->cornerP = S->getObjects();
		this->innerP = S->getObjects_(); //ベースネットワークの内部点，xnetworkと一時的に接続している
		this->outerP.reserve(S->getObjects__().size());
		Print("penetrateQがtrueの点：面を貫通した後にある点がObjects__．貫通した先の点が内部の場合もあるので，これが外部とは限らない！外部の点であることを保証するために，内部点に含まれないObj__だけを抜き出しouterPとする");
		for (const auto &p : S->getObjects__())
			if (!MemberQ(this->innerP, p))
				this->outerP.emplace_back(p); //とても大事
		//
		//準備
		Print("準備");
		network::setStatus(target->getPoints(), false);
		for (const auto &net : this->XNets)
			network::setStatus(net->getPoints(), false);
	};

	///////////////////////////////////////////////
	bool checkConnection(const netFp f)
	{
		auto ls = f->networkObject::getLines();
		int s = ls.size();
		if (s < 3)
			return false;

		if (s > 3)
			return false;

		V_netPp ret(s), Ps, Ps0, Ps1;
		for (auto i = 0; i < s; i++)
		{
			Ps0 = ls[i]->getPoints();
			if (Ps0.size() != 2)
				return false;

			Ps1 = ls[(i + 1) % s]->getPoints();
			if (Ps1.size() != 2)
				return false;

			Ps = Intersection(Ps0, Ps1);
			if (Ps.size() != 1)
				return false;
		}
		return true;
	};

	////////////////////////////////////////////////////////////
	void takePossibleXFaces(Network *target, const V_d &cardinal)
	{
		std::cout << __PRETTY_FUNCTION__ << std::endl;
		prepare(target, cardinal);
		//------------------
		auto nets = this->XNets;
		nets.insert(nets.end(), target);
		this->faces.clear();
		//------------------

		Print("干渉に関連する面だけを取り出す");

		for (const auto &p : this->innerP)
			p->setStatus(true);
		for (const auto &p : this->outerP)
			p->setStatus(true);

		Print("三角分割によって干渉ネットワークに面を挿入する前段階として，三角分割する可能性がある面を抜き出す．");

		// for (const auto &net : nets)
		//   for(auto f:net->Faces)
		//     if(!checkConnection(f))
		//       f->Delete();

		for (const auto &net : nets)
			for (const auto &f : net->getFaces())
				if (target == f->getNetwork() && network::AnyStatusTrue(f->getPoints()))
					this->faces.emplace_back(f);

		for (const auto &l : S->getPenetrateLines())
			for (const auto &f : l->getFaces())
				this->faces.emplace_back(f);

		this->faces = DeleteDuplicates(this->faces);
	};
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	void takeFaces(Network *target, const V_d &cardinal)
	{
		std::cout << __PRETTY_FUNCTION__ << std::endl;
		// #ifdef old_takeFaces
		//     prepare(target, cardinal);
		//     //------------------
		//     auto nets = this->XNets;
		//     nets.insert(nets.end(), target);
		//     this->faces.clear();
		//     //------------------
		//     std::cout << message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "targetと同じ`Network`を持ち，かつ干渉に至る内部の領域の面を取得する．\nただし，干渉ネットワークのFaceの`Face->Netowrk`は，干渉ネットワーク自身ではなく，\n三角分割された面のネットワークを指していることを想定している．") << std::endl;
		//     /*
		//       まず，取得する面の条件は，
		// * innerPと少なくとも１点は接する面
		// * OuterPを含んではいない面
		// * 干渉していない面（面の辺が干渉していないこと）
		// */

		//     network::setStatus(target->Points, false);
		//     for (const auto &net : this->XNets)
		//       network::setStatus(net->Points, false);
		//     for (const auto &p : this->innerP)
		//       p->setStatus(true);

		//     for (const auto &net : nets)
		//     {
		//       auto fs = net->Faces;
		// #ifdef _OPENMP
		//       std::cout << message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "") << std::endl;
		// #pragma omp parallel for
		// #endif
		//       for (auto i = 0; i < fs.size(); i++)
		//         if (target == fs[i]->getNetwork() && !fs[i]->intersectQ() /*３点が内部でも内部でも内部面とは限らない*/)
		//         {
		//           auto ps = fs[i]->getPoints();
		//           if (network::AnyStatusTrue(ps) && !AnyMemberQ(this->outerP, ps))
		//           {
		// #ifdef _OPENMP
		// #pragma omp critical
		// #endif
		//             {
		//               this->faces.emplace_back(fs[i]);
		//             }
		//           }
		//         }
		//     }

		//     Print("しかしこのままでは，内部であっても取得されない面が現れる（innerPに接することなく，cornerPだけに周囲を囲まれている面．外部にもこのような面があり，区別するのが難しい）．\nそこで，取得した内部の面からintxnに到達するまで探査・取得していく．");

		//     network::setStatus(target->Faces, false);
		//     for (const auto &net : this->XNets)
		//       network::setStatus(net->Faces, false);

		//     for (const auto &f : this->faces)
		//       f->setStatus(true);

		//     for (const auto &p : this->outerP)
		//       for (const auto &f : p->getFaces())
		//         f->setStatus(true);

		//     bool found = false;
		//     netFp F;
		//     do
		//     {
		//       found = false;
		//       auto tmp = this->faces;
		//       for (const auto &f : tmp)
		//         if (f /*nullptrではなく*/)
		//           for (const auto &l : f->networkObject::getLines())
		//             if (l /*nullptrではなく*/)
		//               if ((F = (*l)(f)) /*nullptrではなく*/)
		//                 if (/*未取得なら*/ !F->getStatus())
		//                   if (/*交線ではないなら*/ !l->isXLine())
		//                     if (!AnyMemberQ(this->outerP, F->getPoints()) && !F->intersectQ() /*干渉していてはダメ*/)
		//                     {
		//                       // Print("Fは外部点を含んではならない");
		//                       this->faces.emplace_back(F);
		//                       F->setStatus(true);
		//                       found = true;
		//                     }
		//     } while (found);

		//     std::cout << message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "case 1が終了") << std::endl;

		//     this->faces = DeleteDuplicates(this->faces);

		// #else

		prepare(target, cardinal);
		searcherIntxn s;
		s.set(target->getNearestFace(cardinal));
		s.addNetworks(this->XNets);
		s.search();
		this->faces = DeleteDuplicates(s.getObjects());

		// #endif
	};
	///////////////
	void show()
	{
		std::cout << "  innerP : " << innerP.size() << std::endl;
		std::cout << " cornerP : " << cornerP.size() << std::endl;
		std::cout << "  outerP : " << outerP.size() << std::endl;
		std::cout << "   faces : " << faces.size() << std::endl;
		std::cout << "    type : " << type << std::endl;
		for (const auto &n : this->XNets)
			std::cout << n->getName() << ", ";
		std::cout << std::endl;
	};
	/////////////
};
/*selector_code*/

#include "NetworkUtility.hpp"

#endif
