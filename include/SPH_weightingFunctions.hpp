#ifndef SPH_weightingFunctions_H
#define SPH_weightingFunctions_H

#include "kernelFunctions.hpp"

//値にかかる重みをベクトルとして出力する関数
//値を未知変数とする，離散化された方程式を作成する際には，未知変数の重みを抜き出す必要がある
//その際に便利である
V_d W_Bspline(const VV_d &X, V_d dX /*出力する*/, const V_d &a, const double h, const int order = 3)
{
    if (order == 3)
    {
        for (auto i = 0; i < X.size(); i++)
            dX[i] = kernel_Bspline3(X[i], a, h) * dX[i];
    }
    else
    {
        for (auto i = 0; i < X.size(); i++)
            dX[i] = kernel_Bspline5(X[i], a, h) * dX[i];
    }
    return dX;
};
/* ------------------------------------------------------ */
/*                     grad 勾配を計算する際は，             */
/* ------------------------------------------------------ */
VV_d W_grad_Bspine(VV_d X /*出力する*/, const V_d &dX, const V_d &a, const double h, const int order = 3)
{
    V_d zeros(3, 0.);
    if (order == 3)
    {
        for (auto i = 0; i < X.size(); i++)
        {
            if (Norm(X[i] - a) > 1E-10)
            {
                X[i] = grad_kernel_Bspline3(X[i], a, h) * dX[i];
                // X[i] = grad_kernel_Bspline3(X[i], a, e) * dX[i];
            }
            else
                X[i] = zeros;
        }
    }
    else
    {
        for (auto i = 0; i < X.size(); i++)
        {
            if (Norm(X[i] - a) > 1E-10)
            {
                X[i] = grad_kernel_Bspline5(X[i], a, h) * dX[i];
                // X[i] = grad_kernel_Bspline3(X[i], a, e) * dX[i];
            }
            else
                X[i] = zeros;
        }
    }
    return X;
};
//
size_t findIndexOfClosestElement(const VV_d &X, const V_d &a)
{
    double min_r = 1E+100;
    size_t i_a;
    for (auto i = 0; i < X.size(); i++)
        if (Norm(X[i] - a) < min_r)
        {
            i_a = i;
            min_r = Norm(X[i] - a);
        }
    return i_a;
};
//! Dot(値の列,W_grad_Bspine3(X,dX,a,h))
//計算は，
// Dot({v0,v1,v2},{{dwdx0,dwdy0,dwdz0},{dwdx1,dwdy1,dwdz1},{dwdx2,dwdy2,dwdz2}})={dwdx0 v0+dwdx1 v1+dwdx2 v2,dwdy0 v0+dwdy1 v1+dwdy2 v2,dwdz0 v0+dwdz1 v1+dwdz2 v2}
//となり，dwdx0 v0+dwdx1 v1+dwdx2 v2は，確かに正しい計算である．
/* ------------------------------------------------------ */
/*                divergence 発散を計算する際は，  　        */
/* ------------------------------------------------------ */
//! Sum(ElementWiseDot(ベクトル値の列,W_grad_Bspine3(X,dX,a,h)))とする
V_d ElementWiseDot(const VV_d &V, const VV_d &W)
{
    V_d ret(V.size(), 0.);
    for (auto i = 0; i < V.size(); i++)
        ret[i] = Dot(V[i], W[i]);
    return ret;
};
/* ------------------------------------------------------ */
/*                     ラプラシアンの計算方法          　     */
/* ------------------------------------------------------ */
//! Dot(W_laplacian_Bspine3(X, dX, a, h), V);
V_d W_laplacian_Bspine3(const VV_d &X, V_d dX /*出力する*/, const V_d &a, const double e)
{
    V_d r;
    double rr = 0;
    for (auto i = 0; i < X.size(); i++)
    {
        r = X[i] - a;
        rr = Dot(r, r);
        if (rr > 1E-10)
            dX[i] = Dot(r, grad_kernel_Bspline3(X[i], a, e)) / rr * dX[i];
        else
            dX[i] = 0.;
    }
    return dX;
};
V_d W_laplacian_Bspine5(const VV_d &X, V_d dX /*出力する*/, const V_d &a, const double e)
{
    V_d r;
    double rr = 0;
    for (auto i = 0; i < X.size(); i++)
    {
        r = X[i] - a;
        rr = Dot(r, r);
        if (rr > 1E-10)
            dX[i] = Dot(r, grad_kernel_Bspline5(X[i], a, e)) / rr * dX[i];
        else
            dX[i] = 0.;
    }
    return dX;
};
V_d W_laplacian_Bspine(const VV_d &X, V_d dX /*出力する*/, const V_d &a, const double e, const int order = 3)
{
    if (order == 3)
        return W_laplacian_Bspine3(X, dX, a, e);
    else
        return W_laplacian_Bspine5(X, dX, a, e);
};
/* ------------------------------------------------------ */
/*                         改良した核関数                   */
/* ------------------------------------------------------ */
// 2021/07/29以前
// V_d W_Artificial_Viscosity_Bspine_Monaghan1992(const VV_d &X,
//                                                const V_d &mass,
//                                                const V_d &rho,
//                                                const VV_d &U,
//                                                const V_d &a, const double h, const int order = 3)
// {
//     //Xにはaが必ず含まれている，まずは，aをXから探す．
//     int i_a = findIndexOfClosestElement(X, a);

//     V_d ret(3, 0.), Uij(3, 0.), rij(3, 0.), grad(3, 0.);
//     double alpha = 0.02;
//     double cij = 1466.;
//     double muij, PIij;
//     for (auto i = 0; i < X.size(); i++)
//     {
//         if (i_a != i)
//         {
//             Uij = U[i_a] - U[i];
//             rij = X[i_a] - X[i];
//             if (Dot(Uij, rij) < 0)
//             {
//                 muij = h * Dot(Uij, rij) / (Dot(rij, rij) + 0.01 * h * h);
//                 PIij = -alpha * cij * muij / ((rho[i] + rho[i_a]) / 2.);
//                 grad = (order == 3 ? grad_kernel_Bspline3(X[i], a, h) : grad_kernel_Bspline5(X[i], a, h));
//                 ret += (mass[i] * PIij) * grad;
//             }
//         }
//     }
//     return ret;
// };
// 2021/07/29以降
V_d W_Artificial_Viscosity_Bspine_Monaghan1992(const VV_d &X,
                                               const V_d &mass,
                                               const V_d &rho,
                                               const VV_d &U,
                                               const V_d &a, const double h, const int order = 3)
{
    // Xにはaが必ず含まれている，まずは，aをXから探す．
    int i_a = findIndexOfClosestElement(X, a);

    V_d ret(3, 0.), Uij(3, 0.), rij(3, 0.), grad(3, 0.);
    double alpha = 0.02;
    double cij = 1466.;
    double PIij, div;
    for (auto i = 0; i < X.size(); i++)
    {
        if (i_a != i)
        {
            Uij = U[i_a] - U[i];
            rij = X[i_a] - X[i];
            if (Dot(Uij, rij) < 0)
            {
                div = Dot(Uij, rij) / Dot(rij, rij);
                PIij = -alpha * cij * h * div / ((rho[i] + rho[i_a]) / 2.);
                grad = (order == 3 ? grad_kernel_Bspline3(X[i], a, h) : grad_kernel_Bspline5(X[i], a, h));
                ret += (mass[i] * PIij) * grad;
            }
        }
    }
    return ret;
};

//
// VV_d W_grad_Bspine_Monaghan1992_dividedByDensity(const VV_d &X,
// 									const V_d &mass,
// 									const V_d &rho,
// 									const V_d &a, const double e, const int order = 3)
// {
// 	int i_a = findIndexOfClosestElement(X, a);

// 	auto rho_a = rho[i_a];
// 	V_d grad(3, 0.);
// 	VV_d ret(X.size(), V_d(3, 0.));
// 	for (auto i = 0; i < X.size(); i++)
// 	{
// 		if (i_a != i)
// 		{
// 			grad = (order == 3 ? grad_kernel_Bspline3(X[i], a, e) : grad_kernel_Bspline5(X[i], a, e));
// 			ret[i] = mass[i] / std::pow(rho[i], 2.) * grad;
// 			ret[i_a] += mass[i] * grad;
// 		}
// 	}
// 	ret[i_a] /= (rho_a * rho_a);
// 	return ret;
// };
VV_d W_grad_Bspine_Monaghan1992_dividedByDensity(const VV_d &X,
                                                 const V_d &mass,
                                                 const V_d &rho,
                                                 const V_d &a,
                                                 const double h,
                                                 const int order = 3)
{
    int i_a = findIndexOfClosestElement(X, a);
    VV_d W_grad = W_grad_Bspine(X, mass / rho, a, h, order); // return
    W_grad[i_a] = {0., 0., 0.};
    // increment用
    for (auto i = 0; i < X.size(); i++)
        if (i_a != i)
        {
            W_grad[i] = W_grad[i] / rho[i];
            W_grad[i_a] += W_grad[i] / rho[i_a];
        }
    return W_grad;
};
//
VV_d W_grad_Bspine_Monaghan1992(const VV_d &X, const V_d &dX, const V_d &a, const double e, const int order = 3)
{
    VV_d ret(X.size(), V_d(3, 0.));
    V_d tot(a.size(), 0.), tmp(a.size(), 0.);
    int i_a = 0;
    double min_r = 1E+100, nr;
    for (auto i = 0; i < X.size(); i++)
    {
        nr = Norm(X[i] - a);
        if (nr > 1E-10)
        {
            tmp = (order == 3 ? grad_kernel_Bspline3(X[i], a, e) : grad_kernel_Bspline5(X[i], a, e)) * dX[i];
            ret[i] = tmp;
            tot -= tmp;
        }
        if (nr < min_r)
        {
            i_a = i;
            min_r = nr;
        }
    }
    ret[i_a] = tot;
    return -ret;
};

V_d W_laplacian_Bspine_Brookshaw1985(const VV_d &X, const V_d &dX, const V_d &a, const double e, const int order = 3)
{
    V_d ret(X.size(), 0.), r;
    double tot = 0., tmp = 0., nr = 0., min_r = 1E+100;
    int i_a = 0;
    for (auto i = 0; i < X.size(); i++)
    {
        nr = Norm(r = (X[i] - a));
        if (nr > 1E-10)
        {
            tmp = 2. * Dot(r / std::pow(nr, 2.), (order == 3 ? grad_kernel_Bspline3(X[i], a, e) : grad_kernel_Bspline5(X[i], a, e))) * dX[i];
            ret[i] = tmp;
            tot -= tmp;
        }
        if (nr < min_r)
        {
            i_a = i;
            min_r = nr;
        }
    }
    ret[i_a] = tot;
    return -ret;
    // eq.(91) D. J. Price, “Smoothed particle hydrodynamics and magnetohydrodynamics,” J. Comput. Phys., vol. 231, no. 3, pp. 759–794, 2012.
};

    //
    // V_d W_grad_Bspine_Monaghan1992_STD(const VV_d &X,
    //                                    const V_d &mass,
    //                                    const V_d &rho,
    //                                    const V_d &a, const double e, const int order = 3) {
    // 	V_d nr(X.size(), 0.), r(X.size(), 0.);
    // 	double tot = 0., tmp = 0., min_r = 1E+100;
    // 	int i_a = 0;
    // 	for (auto i = 0; i < X.size(); i++) {
    // 		nr[i] = Norm(r[i] = (X[i] - a));
    // 		if (nr[i] < min_r) {
    // 			i_a = i;
    // 			min_r = nr;
    // 		}
    // 	}
    // 	auto rho_a = rho[i_a];
    // 	V_d ret(X.size(), 0.);
    // 	for (auto i = 0; i < X.size(); i++) {
    // 		if (i_a != i) {
    // 			auto grad = (order == 3 ? grad_kernel_Bspline3(X[i], a, e) : grad_kernel_Bspline5(X[i], a, e));
    // 			ret[i] = std::pow(rho[i], 2.) * grad;
    // 			ret[i_a] += mass[i] * grad;
    // 		}
    // 	}
    // 	ret[i_a] /= (rho_a * rho_a);
    // 	return ret;
    // };

    // VV_d ElementWiseProduct(VV_d V, const V_d &W)
    // {
    // 	for (auto i = 0; i < V.size(); i++)
    // 		V[i] = V[i] * W[i];
    // 	return V;
    // };

    // V_d ElementWiseAdd(V_d V, const V_d &W)
    // {
    // 	for (auto i = 0; i < V.size(); i++)
    // 		V[i] = V[i] + W[i];
    // 	return V;
    // };
    //-----------------------------------------------------------------------------------------
    // class InterpolationVectorSPH {
    // 	/* -------------- 各positionでサンプリングしたスカラーの補間 ------------- */
    //    private:
    // 	VV_d X;  //position
    // 	V_d dX;  //volume
    // 	VV_d Values;
    // 	int sample_size;
    // 	std::function<double(V_d, V_d)> kernel;
    // 	std::function<V_d(V_d, V_d)> grad_kernel;
    // 	double sml;

    //    public:
    // 	InterpolationVectorSPH(const VV_d &X_IN, const VV_d &Values_IN, const V_d &dX_IN, double h)
    // 	    : X(X_IN), Values(Values_IN), dX(dX_IN), sample_size(X_IN[0].size()) {
    // 		if (X.size() != Values.size()) throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
    // 		this->sml = h;
    // 		this->kernel = [this](const V_d &x, const V_d &a) { return kernel_Bspline3(x, a, this->sml); };
    // 		this->grad_kernel = [this](const V_d &x, const V_d &a) { return grad_kernel_Bspline3(x, a, this->sml); };
    // 		// std::cout << "this->sml = " << this->sml << std::endl;
    // 	};
    // 	~InterpolationVectorSPH(){};
    // 	/* ------------------------ メソッド ------------------------ */
    // 	V_d operator()(const V_d &a) const {
    // 		return Dot(this->Values, W_SPH(kernel, X, dX, a));
    // 	};
    // 	V_d grad(const V_d &a) const {
    // 		return ElementWiseDot(this->Values, grad_W_SPH(grad_kernel, X, dX, a));
    // 	};
    // 	double div(const V_d &a) const {
    // 		return Sum(this->grad(a));
    // 	};
    // 	V_d laplacian(const V_d &a) const {
    // 		return Dot(this->Values, laplacian_W_SPH(grad_kernel, X, dX, a));
    // 	};
    // 	V_d laplacian_mod(const V_d &a) const {
    // 		//!計算力学の定石p.262 (13.8)
    // 		//!非圧縮性流体に対しては，下のタイプが使われる
    // 		V_d ret(a.size(), 0.), r(a.size(), 0.);
    // 		auto V = (*this)(a);
    // 		double abs_r = 0;
    // 		for (auto j = 0; j < this->X.size(); j++) {
    // 			if (abs_r = Norm(r = a - X[j]) > 1E-14 /*自身の和はとらない*/)
    // 				ret += (V - Values[j]) / pow(abs_r, 2) * Dot(r, grad_kernel(X[j], a)) * dX[j] /*各要素毎の積が計算されるようになっている*/;
    // 		}
    // 		return 2. * ret;
    // 	};
    // };
    // class InterpolationSPH {
    // 	/* -------------- 各positionでサンプリングしたスカラーの補間 ------------- */
    //    private:
    // 	VV_d X;  //position
    // 	V_d dX;  //volume
    // 	V_d Values;
    // 	int sample_size;
    // 	std::function<double(V_d, V_d)> kernel;
    // 	std::function<V_d(V_d, V_d)> grad_kernel;
    // 	double sml;

    //    public:
    // 	InterpolationSPH(const VV_d &X_IN, const V_d &Values_IN, const V_d &dX_IN, double h)
    // 	    : X(X_IN), Values(Values_IN), dX(dX_IN), sample_size(X_IN[0].size()) {
    // 		if (X.size() != Values.size())
    // 			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
    // 		this->sml = h;  //!factorを3.3とした，計算力学の定石p.270
    // 		this->kernel = [this](const V_d &x, const V_d &a) { return kernel_Bspline3(x, a, this->sml / 2.); };
    // 		this->grad_kernel = [this](const V_d &x, const V_d &a) { return grad_kernel_Bspline3(x, a, this->sml / 2.); };
    // 		// std::cout << "this->sml = " << this->sml << std::endl;
    // 	};
    // 	~InterpolationSPH(){};
    // 	/* ------------------------ メソッド ------------------------ */
    // 	double operator()(const V_d &a) const {
    // 		double ret(0.);
    // 		for (auto i = 0; i < this->X.size(); i++) {
    // 			ret += this->Values[i] * kernel(X[i], a) * dX[i];
    // 		}
    // 		return ret;
    // 	};
    // 	V_d grad(const V_d &a) const {
    // 		//!粒子法p.25 (2.89)
    // 		V_d ret(a.size(), 0.);
    // 		for (auto j = 0; j < this->X.size(); j++)
    // 			if (Norm(X[j] - a) > 1E-12 /*自身の和はとらない*/)
    // 				ret += Values[j] * grad_kernel(X[j], a) * dX[j] /*各要素毎の積が計算されるようになっている*/;
    // 		return ret;
    // 	};
    // 	V_d laplacian(const V_d &a) const {
    // 		V_d ret(a.size(), 0.);
    // 		auto V = (*this)(a);
    // 		for (auto i = 0; i < this->X.size(); i++) {
    // 			auto r = a - X[i];
    // 			if (Norm(r) > 1E-14 /*自身の和はとらない*/)
    // 				ret += 2. * (V - this->Values[i]) / pow(Norm(r), 2) * Dot(r, this->grad_kernel(a, X[i]) /*|a-X[i]|がでる*/) * dX[i] /*各要素毎の積が計算されるようになっている*/;
    // 		}
    // 		return ret;
    // 	};
    // };

#endif

// Networkをインクルードしている場合のみ使える
//#if !defined(SPH_weightingFunctions_Network_H) && defined(Network_H)
#define SPH_weightingFunctions_Network_H

#include "kernelFunctions.hpp"

//値にかかる重みをベクトルとして出力する関数
//値を未知変数とする，離散化された方程式を作成する際には，未知変数の重みを抜き出す必要がある
//その際に便利である
V_d W_Bspline(const std::vector<Tddd> &X, V_d dX /*出力する*/, const Tddd &a, const double h, const int order = 3)
{
    if (order == 3)
    {
        for (auto i = 0; i < X.size(); i++)
            dX[i] = kernel_Bspline3(X[i], a, h) * dX[i];
    }
    else
    {
        for (auto i = 0; i < X.size(); i++)
            dX[i] = kernel_Bspline5(X[i], a, h) * dX[i];
    }
    return dX;
};
/* ------------------------------------------------------ */
/*                     grad 勾配を計算する際は，             */
/* ------------------------------------------------------ */
std::vector<Tddd> W_grad_Bspine(std::vector<Tddd> X /*出力する*/, const V_d &dX, const Tddd &a, const double h, const int order = 3)
{
    if (order == 3)
    {
        for (auto i = 0; i < X.size(); i++)
        {
            if (Norm(X[i] - a) > 1E-10)
            {
                X[i] = grad_kernel_Bspline3(X[i], a, h) * dX[i];
                // X[i] = grad_kernel_Bspline3(X[i], a, e) * dX[i];
            }
            else
                X[i] = {0., 0., 0.};
        }
    }
    else
    {
        for (auto i = 0; i < X.size(); i++)
        {
            if (Norm(X[i] - a) > 1E-10)
            {
                X[i] = grad_kernel_Bspline5(X[i], a, h) * dX[i];
                // X[i] = grad_kernel_Bspline3(X[i], a, e) * dX[i];
            }
            else
                X[i] = {0., 0., 0.};
        }
    }
    return X;
};
std::vector<Tddd> W_grad_Bspine(const std::unordered_set<networkPoint *> &ps,
                                const networkPoint *const a,
                                const double h,
                                const int order = 3)
{
    std::vector<Tddd> ret(ps.size());
    int i = 0;
    if (order == 3)
    {
        for (const auto &p : ps)
        {
            if (Norm(p->getXtuple() - a->getXtuple()) > 1E-10)
                ret[i++] = grad_kernel_Bspline3(p->getXtuple(), a->getXtuple(), h) * p->volume;
            else
                ret[i++] = {0., 0., 0.};
        }
    }
    else
    {
        for (const auto &p : ps)
        {
            if (Norm(p->getXtuple() - a->getXtuple()) > 1E-10)
                ret[i++] = grad_kernel_Bspline5(p->getXtuple(), a->getXtuple(), h) * p->volume;
            else
                ret[i++] = {0., 0., 0.};
        }
    }
    return ret;
};
/* ------------------------------------------------------ */
/* ------------------- 仮想マーカーにおける流速を算出 ------------------ */
std::tuple<Tddd, Tddd, double, double> U_tmpU_mass_density_SPH_IDW(const std::unordered_set<networkPoint *> &ps,
                                                                   const Tddd &X,
                                                                   const double power = 2)
{
    Tddd tmp_U = {0, 0, 0}, U = {0, 0, 0};
    double mass = 0., density = 0., total = 0, w = 0.;
    for (const auto &p : ps)
    {
        w = std::pow(1. / Norm(p->getXtuple() - X), power);
        total += w;
        U += w * p->U_SPH;
        tmp_U += w * p->tmp_U_SPH;
        mass += w * p->mass;
        density += w * p->density;
    }
    return {U / total,
            tmp_U / total,
            mass / total,
            density / total};
};
// 2021/11/10
// std::tuple<Tddd, double, double> U_mass_density_SPH_IDW_using_tmpU(const std::unordered_set<networkPoint *> &ps,
//                                                                    const Tddd &X,
//                                                                    const double power = 1)
// {
//     Tddd U = {0, 0, 0};
//     double mass = 0., density = 0.;
//     double total = 0, w = 0.;
//     for (const auto &p : ps)
//     {
//         w = std::pow(1. / Norm(p->getXtuple() - X), power);
//         total += w;
//         U += w * p->tmp_U_SPH;
//         mass += w * p->mass;
//         density += w * p->density;
//     }
//     U /= total;
//     mass /= total;
//     density /= total;
//     return {U, mass, density};
// };
// //
double dummy_pressure_Asai_Modified(const Buckets<networkPoint> &B_water,
                                    const networkPoint *dummy_p /*dummy point*/,
                                    const std::unordered_set<networkFace *> &boundary_face,
                                    const double mirroring_distance,
                                    const double power = 2)
{
    Tddd closestX = getClosestFacePoint(dummy_p, boundary_face, mirroring_distance);
    // ダミー粒子のみ使える関数
    auto [f, t0, t1, d, dx] = dummy_p->particlize_info;
    auto X = dummy_p->getXtuple();
    auto Y = X + 2. * (closestX - X); // oppositeX(dummy_p->particlize_info);
    auto p_at_Y = B_water.getObjects_unorderedset(Y, 8 /*深さ最大*/, 1);
    // double sign = Normalize(Dot(f->getXtuple() - Y, f->getNormalTuple()) * f->getNormalTuple()); // Y->F (+ or -)
    Tddd gravity = {0., 0., -9.81};
    double pressure = 0, total = 0, norm, nu, tmp, pn;
    Tddd accel = {0, 0, 0};
    if (!p_at_Y.empty())
    {
        for (const auto &q : p_at_Y)
        {
            // ここを修正
            // 重み付き平均した反対位置OXと壁との距離を確認し
            // OX->
            nu = q->mu_SPH / q->density;
            pn = q->density * Dot(nu * q->lap_U + gravity - accel, f->getNormalTuple());
            tmp = q->pressure_SPH + /*修正*/ pn * Dot(X - q->getXtuple(), f->getNormalTuple());
            // tmp = (q->pressure_SPH);
            norm = Norm(X - q->getXtuple());
            if (norm < 1E-15)
                return tmp;
            auto w = std::pow(1. / norm, power);
            total += w;
            pressure += w * tmp;
        }
        pressure /= total;
        return pressure;
    }
    return pressure;
};

double dummy_pressure_Asai(const Buckets<networkPoint> &B_water,
                           const networkPoint *dummy_p /*dummy point*/,
                           const double power = 2)
{
    // ダミー粒子のみ使える関数
    auto [f, t0, t1, d, dx] = dummy_p->particlize_info;
    auto Y = oppositeX(dummy_p->particlize_info);
    auto X = dummy_p->getXtuple();
    auto p_at_Y = B_water.getObjects_unorderedset(Y, 8 /*深さ最大*/, 1);
    // double sign = Normalize(Dot(f->getXtuple() - Y, f->getNormalTuple()) * f->getNormalTuple()); // Y->F (+ or -)
    Tddd gravity = {0., 0., -9.81};
    double pressure = 0, total = 0, norm, nu, tmp, pn;
    Tddd accel = {0, 0, 0};
    if (!p_at_Y.empty())
    {
        for (const auto &q : p_at_Y)
        {
            // ここを修正
            // 重み付き平均した反対位置OXと壁との距離を確認し
            // OX->
            nu = q->mu_SPH / q->density;
            pn = q->density * Dot(nu * q->lap_U + gravity - accel, f->getNormalTuple());
            tmp = q->pressure_SPH + /*修正*/ pn * Dot(X - q->getXtuple(), f->getNormalTuple());
            // tmp = (q->pressure_SPH);
            norm = Norm(X - q->getXtuple());
            if (norm < 1E-15)
                return tmp;
            auto w = std::pow(1. / norm, power);
            total += w;
            pressure += w * tmp;
        }
        pressure /= total;
        return pressure;
    }
    return pressure;
};
//
Tddd velocity_SPH_IDW(const std::unordered_set<networkPoint *> &ps,
                      const Tddd &X,
                      const double h,
                      const int order = 3)
{
    Tddd ret = {0, 0, 0};
    double total = 0, norm;
    for (const auto &p : ps)
    {
        norm = Norm(p->getXtuple() - X);
        ret += 1. / norm * p->U_SPH;
        total += 1. / norm;
    }
    return ret / total;
};

Tddd velocity_SPH(const std::unordered_set<networkPoint *> &ps,
                  const Tddd &X,
                  const double h,
                  const int order = 3)
{
    //* これはシンプルに勾配計算の改良版である．
    Tddd ret = {0, 0, 0};
    for (const auto &p : ps)
        ret += p->U_SPH * p->volume * kernel_Bspline5(p->getXtuple(), X, h);
    return ret;
};
double density_SPH(const std::unordered_set<networkPoint *> &ps,
                   const networkPoint *const a,
                   const double h,
                   const int order = 3)
{
    //* これはシンプルに勾配計算の改良版である．
    double ret = a->density * a->volume * kernel_Bspline5(a->getXtuple(), a->getXtuple(), h);
    for (const auto &p : ps)
        if (p != a)
            ret += p->density * p->volume * kernel_Bspline5(p->getXtuple(), a->getXtuple(), h);
    return ret;
};
Tddd normal(const std::unordered_set<networkPoint *> &ps,
            const networkPoint *const a,
            const double h)
{
    Tddd ret = {0, 0, 0};
    for (const auto &p : ps)
        if (p != a)
            ret += grad_kernel_Bspline3(p->getXtuple(), a->getXtuple(), h);
    return Normalize(ret);
};
/* ------------------------------------------------------ */
struct IntersectionsSphereTrianglesLines
{
    std::unordered_set<networkFace *> faces;
    std::vector<std::tuple<networkFace *, networkFace *, T2Tddd>> lines;
    IntersectionsSphereTrianglesLines(const std::unordered_set<networkFace *> &F) : faces({}), lines({})
    {
        addFacesAndLines(F);
    };
    IntersectionsSphereTrianglesLines(const std::unordered_map<networkFace *, Tddd> &F)
    {
        // b! この構造体は，まず与えられた面と，その面同士の交線を作成する．
        // b! 次にgetでは，面または交線が，与えられた点と正対する場合にその位置を返す．
        std::unordered_set<networkFace *> FF;
        FF.reserve(F.size());
        for (const auto &[f, _] : F)
            FF.emplace(f);
        addFacesAndLines(FF);
    };
    void addFacesAndLines(const std::unordered_set<networkFace *> &F)
    {
        // 保存はunordered_setで
        this->faces.reserve(F.size());
        for (const auto &f : F)
            this->faces.emplace(f);

        // ベクトル化
        std::vector<networkFace *> FACES(this->faces.begin(), this->faces.end());

        T3Tddd fi, fj;
        for (auto i = 0; i < FACES.size(); ++i)
        {
            fi = FACES[i]->getX_Vertices();
            fi += 0.1 * (T3Tddd{std::get<0>(fi) - Mean(fi),
                                std::get<1>(fi) - Mean(fi),
                                std::get<2>(fi) - Mean(fi)});
            for (auto j = i; j < FACES.size(); ++j)
            {
                fj = FACES[j]->getX_Vertices();
                fj += 0.1 * (T3Tddd{std::get<0>(fj) - Mean(fj),
                                    std::get<1>(fj) - Mean(fj),
                                    std::get<2>(fj) - Mean(fj)});
                auto intxn = IntersectionTriangles(fi, fj);
                if (intxn.isIntersecting)
                    lines.push_back({FACES[i], FACES[j], intxn.L});
            }
        }
    };

    std::vector<std::tuple<Tddd, Tddd>> X_Y;
    std::vector<std::tuple<networkFace *, networkFace *, Tddd, Tddd>> F_F_X_Y;
    std::vector<std::tuple<Tddd, Tddd>> get(const networkPoint *p, const double h)
    {
        X_Y.clear();
        // b! 次にgetでは，面または交線が，与えられた点と正対する場合にその位置を返す．
        //! 面との干渉
        for (const auto &f : this->faces)
        {
            auto intxn = IntersectionSphereTriangle(p->getXtuple(), h, f->getX_Vertices());
            if (intxn.isIntersecting)
                X_Y.push_back({intxn.X, p->getXtuple() + 2. * (intxn.X - p->getXtuple())});
        }
        //! 線との干渉
        for (const auto &[F0, F1, L] : this->lines)
        {
            auto intxn = IntersectionSphereLine(p->getXtuple(), h, L);
            if (intxn.isIntersecting)
                X_Y.push_back({intxn.X, p->getXtuple() + 2. * (intxn.X - p->getXtuple())});
        }
        return X_Y;
    };
    std::vector<std::tuple<networkFace *, networkFace *, Tddd, Tddd>> getFFXY(const networkPoint *p, const double h)
    {
        F_F_X_Y.clear();
        // b! 次にgetでは，面または交線が，与えられた点と正対する場合にその位置を返す．
        //! 面との干渉
        for (const auto &f : this->faces)
        {
            auto intxn = IntersectionSphereTriangle(p->getXtuple(), h, f->getX_Vertices());
            if (intxn.isIntersecting)
                F_F_X_Y.push_back({f, nullptr, intxn.X, p->getXtuple() + 2. * (intxn.X - p->getXtuple())});
        }
        //! 線との干渉
        for (const auto &[F0, F1, L] : this->lines)
        {
            auto intxn = IntersectionSphereLine(p->getXtuple(), h, L);
            if (intxn.isIntersecting)
                F_F_X_Y.push_back({F0, F1, intxn.X, p->getXtuple() + 2. * (intxn.X - p->getXtuple())});
        }
        return F_F_X_Y;
    };
};
#define free_slip_boundary_condition
// #define boundary_gurd

#if defined(boundary_gurd)
double boundary_gurd_value(const double distance, const double h)
{
    double additional = tanh((2. * M_PI) / (h / 8.) * (-distance + (h / 8.)) - M_PI) + 1.;
    return 1. + additional;
};
#endif

//@ ------------------------------------------------------ */
//@                         速度の発散                       */
//@                          div(U)                        */
//@ ------------------------------------------------------ */
double div_U_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                              const networkPoint *const a,
                              const double h)
{
    auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
    double ret = 0;
    Tddd n, U, velocity;
    for (const auto &p : ps)
    {
        if (p != a)
            ret += p->mass * Dot(p->U_SPH - a->U_SPH, grad_kernel_Bspline5(p->getXtuple(), a->getXtuple(), h));
        // b$ ------------------------ ポリゴン ------------------------ */
        for (const auto &[F0, F1, X, Y] : INTXN.getFFXY(p, h))
        {
            if (F1)
            {
                std::get<0>(velocity) = std::get<0>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
                std::get<1>(velocity) = std::get<1>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
                std::get<2>(velocity) = std::get<2>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
            }
            else
            {
                std::get<0>(velocity) = std::get<0>(F0->getNetwork()->velocity);
                std::get<1>(velocity) = std::get<1>(F0->getNetwork()->velocity);
                std::get<2>(velocity) = std::get<2>(F0->getNetwork()->velocity);
            }
#if defined(free_slip_boundary_condition)
            U = p->U_SPH - 2 * n * Dot(p->U_SPH, n) + Dot(velocity, n) * n;
#elif defined(no_slip_boundary_condition)
            U = -p->U_SPH + Dot(velocity, n) * n;
#else
            U = -n * Dot(p->U_SPH, n) + Dot(velocity, n) * n;
#endif
            ret += p->mass * Dot(U - a->U_SPH, grad_kernel_Bspline5(Y, a->getXtuple(), h))
#if defined(boundary_gurd)
                   * boundary_gurd_value(Norm(Y - p->getXtuple()), h);
#else
                ;
#endif
        }
        // b%$------------------------------------------------------ */
    }
    return ret / a->density;
};

double div_U(const std::unordered_set<networkPoint *> &ps,
             const networkPoint *const a,
             const double h)
{
    //* これはシンプルに勾配計算の改良版である．
    double ret = 0;
    for (const auto &p : ps)
        if (p != a)
            ret += p->mass *
                   Dot(p->U_SPH - a->U_SPH,
                       grad_kernel_Bspline5(p->getXtuple(), a->getXtuple(), h));
    return ret / a->density;
};
//
double div_tmp_U_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                  const networkPoint *const a,
                                  const double h,
                                  const int order = 3)
{
    auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
    //* これはシンプルに勾配計算の改良版である．
    double ret = 0;
    Tddd U, n, velocity;
    for (const auto &p : ps)
    {
        if (p != a)
        {
            ret += p->mass * Dot(p->tmp_U_SPH - a->tmp_U_SPH, grad_kernel_Bspline5(a->getXtuple(), p->getXtuple(), h));
        }
        // b$ ------------------------ ポリゴン ------------------------ */
        for (const auto &[F0, F1, X, Y] : INTXN.getFFXY(p, h))
        {
            if (F1)
            {
                std::get<0>(velocity) = std::get<0>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
                std::get<1>(velocity) = std::get<1>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
                std::get<2>(velocity) = std::get<2>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
            }
            else
            {
                std::get<0>(velocity) = std::get<0>(F0->getNetwork()->velocity);
                std::get<1>(velocity) = std::get<1>(F0->getNetwork()->velocity);
                std::get<2>(velocity) = std::get<2>(F0->getNetwork()->velocity);
            }
            n = Normalize(p->getXtuple() - X);
#if defined(free_slip_boundary_condition)
            U = p->tmp_U_SPH - 2 * n * Dot(p->tmp_U_SPH, n) + Dot(velocity, n) * n;
#elif defined(no_slip_boundary_condition)
            U = -p->tmp_U_SPH + Dot(velocity, n) * n;
#else
            U = -n * Dot(p->tmp_U_SPH, n) + Dot(velocity, n) * n;
#endif
            ret += p->mass * Dot(U - a->tmp_U_SPH, grad_kernel_Bspline5(a->getXtuple(), Y, h))
#if defined(boundary_gurd)
                   * boundary_gurd_value(Norm(Y - p->getXtuple()), h);
#else
                ;
#endif
        }
        // b%$------------------------------------------------------ */
    }

    return ret / a->density;
};
double div_tmp_U(const std::unordered_set<networkPoint *> &ps,
                 const networkPoint *const a,
                 const double h,
                 const int order = 3)
{
    // double ret = 0;
    // for (const auto &p : ps)
    //     if (p != a)
    //         ret += p->volume *
    //                Dot(p->tmp_U_SPH, grad_kernel_Bspline5(p->getXtuple(), a->getXtuple(), h));
    // return ret;
    /* ------------------------------------------------------ */
    //* これはシンプルに勾配計算の改良版である．
    double ret = 0;
    for (const auto &p : ps)
        if (p != a)
            ret += p->mass *
                   Dot(p->tmp_U_SPH - a->tmp_U_SPH,
                       grad_kernel_Bspline5(a->getXtuple(), p->getXtuple(), h));
    return ret / a->density;
};
//# ------------------------------------------------------ */
//#                          圧力勾配                       */
//#                         grad(P)                        */
//# ------------------------------------------------------ */
// Tddd grad_P_Monaghan1992_polygon_boundary(const std::unordered_set<networkPoint *> &ps /*流体だけを渡す*/,
//                                           const networkPoint *const a,
//                                           const double h,
//                                           const int order = 3)
// {
//     //密度の演算子の作用下に入れ込んだMonaghan1992の方法
//     Tddd ret = {0, 0, 0}, W;
//     for (const auto &p : ps)
//     {
//         if (p != a)
//             W = grad_kernel_Bspline5(p->getXtuple(), a->getXtuple(), h);
//         else
//             W = {0, 0, 0};

//         for (const auto &[f, intX] : p->getContactFaces())
//             W += grad_kernel_Bspline5(a->getXtuple() + 2. * (intX - a->getXtuple()), a->getXtuple(), h);
//         ret += p->mass * (p->pressure_SPH / std::pow(p->density, 2) + a->pressure_SPH / std::pow(a->density, 2)) * W;
//     }
//     return ret * a->density;
// };

Tddd grad_P_Monaghan1992_polygon_boundary(const std::unordered_set<networkPoint *> &ps /*流体だけを渡す*/,
                                          const networkPoint *const a,
                                          const double h)
{
    //密度の演算子の作用下に入れ込んだMonaghan1992の方法
    Tddd ret = {0, 0, 0}, W;
    Tddd accel = {0, 0, 0}, n;
    Tddd gravity = {0., 0., -9.81};
    auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
    double pressure;
    for (const auto &p : ps)
    {
        if (p != a)
        {
            W = grad_kernel_Bspline5(p->getXtuple(), a->getXtuple(), h);
            ret += p->mass * (p->pressure_SPH / std::pow(p->density, 2) + a->pressure_SPH / std::pow(a->density, 2)) * W;
        }
        // b$ ------------------------ ポリゴン ------------------------ */
        for (const auto &[F0, F1, X, Y] : INTXN.getFFXY(p, h))
        {
            if (F1)
            {
                std::get<0>(accel) = std::get<0>(F0->getNetwork()->acceleration + F1->getNetwork()->acceleration);
                std::get<1>(accel) = std::get<1>(F0->getNetwork()->acceleration + F1->getNetwork()->acceleration);
                std::get<2>(accel) = std::get<2>(F0->getNetwork()->acceleration + F1->getNetwork()->acceleration);
            }
            else
            {
                std::get<0>(accel) = std::get<0>(F0->getNetwork()->acceleration);
                std::get<1>(accel) = std::get<1>(F0->getNetwork()->acceleration);
                std::get<2>(accel) = std::get<2>(F0->getNetwork()->acceleration);
            }
            pressure = (p->pressure_SPH + /*修正*/ p->density * Dot(p->mu_SPH / p->density * p->lap_U + gravity - accel, Y - p->getXtuple()));
            // pressure = p->pressure_SPH;
            ret += p->mass * (pressure / std::pow(p->density, 2) + a->pressure_SPH / std::pow(a->density, 2)) * grad_kernel_Bspline5(Y, a->getXtuple(), h)
#if defined(boundary_gurd)
                   * boundary_gurd_value(Norm(Y - p->getXtuple()), h);
#else
                ;
#endif
        }
        // b%$------------------------------------------------------ */
    }

    return ret * a->density;
};

Tddd grad_P_Monaghan1992(const std::unordered_set<networkPoint *> &ps,
                         const networkPoint *const a,
                         const double h,
                         const int order = 3)
{
    //密度の演算子の作用下に入れ込んだMonaghan1992の方法
    Tddd ret = {0, 0, 0};
    for (const auto &p : ps)
        if (p != a)
            ret += p->mass *
                   (p->pressure_SPH / std::pow(p->density, 2) + a->pressure_SPH / std::pow(a->density, 2)) *
                   grad_kernel_Bspline5(p->getXtuple(), a->getXtuple(), h);
    return ret * a->density;
};
//% ------------------------------------------------------ */
//%                          粘性項                         */
//%                       laplacian(U)                     */
//% ------------------------------------------------------ */
Tddd laplacian_U_Monaghan1992(const networkPoint *const p,
                              const Tddd &X,
                              const Tddd &U,
                              const networkPoint *const a,
                              const double h,
                              const double alpha = 0.1,
                              const double beta = 0.1)
{
    auto Uij = a->U_SPH - U;
    auto rij = a->getXtuple() - X;
    auto dot_rij_rij = Dot(rij, rij);
    auto dot_Uij_rij = Dot(Uij, rij);
    double Cs = 1466., PIij, div_U, m;
    if (dot_Uij_rij < 0)
    {
        m = h / 3. * dot_Uij_rij / dot_rij_rij /*hを消せばUの発散となっている*/;
        PIij = (-alpha * Cs * m + beta * m * m) / ((p->density + a->density) / 2.);
        return p->mass * PIij * grad_kernel_Bspline5(X, a->getXtuple(), h);
    }
    return {0, 0, 0};
};
Tddd laplacian_U_Monaghan1992_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                               const networkPoint *const a,
                                               const double h,
                                               const double alpha = 0.1,
                                               const double beta = 0.1)
{
    /*
    kernel_Bspline5の定義はwij=w(|p-a|)としている．
    */
    //人工粘性
    //まず，nu*laplacian(U)を計算する
    Tddd ret = {0, 0, 0}, n, U, velocity;
    double Cs = 1466.;
    auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
    for (const auto &p : ps)
    {
        if (p != a)
            ret += laplacian_U_Monaghan1992(p, p->getXtuple(), p->U_SPH, a, h, alpha, beta);
        // b$ ------------------------ ポリゴン ------------------------ */
        for (const auto &[F0, F1, X, Y] : INTXN.getFFXY(p, h))
        {
            if (F1)
            {
                std::get<0>(velocity) = std::get<0>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
                std::get<1>(velocity) = std::get<1>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
                std::get<2>(velocity) = std::get<2>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
            }
            else
            {
                std::get<0>(velocity) = std::get<0>(F0->getNetwork()->velocity);
                std::get<1>(velocity) = std::get<1>(F0->getNetwork()->velocity);
                std::get<2>(velocity) = std::get<2>(F0->getNetwork()->velocity);
            }
#if defined(free_slip_boundary_condition)
            U = p->U_SPH - 2 * n * Dot(p->U_SPH, n) + Dot(velocity, n) * n;
#elif defined(no_slip_boundary_condition)
            U = -p->U_SPH + Dot(velocity, n) * n;
#else
            U = -n * Dot(p->U_SPH, n) + Dot(velocity, n) * n;
#endif
            ret += laplacian_U_Monaghan1992(p, Y, U, a, h, alpha, beta)
#if defined(boundary_gurd)
                   * boundary_gurd_value(Norm(Y - p->getXtuple()), h);
#else
                ;
#endif
        }
        // b$ ------------------------------------------------------ */
    }
    double nu = a->mu_SPH / a->density;
    return -ret / nu;
};

Tddd laplacian_U_Monaghan1992(const std::unordered_set<networkPoint *> &ps,
                              const networkPoint *const a,
                              const double h,
                              const double alpha = 0.1,
                              const double beta = 0.1)
{
    /*
    kernel_Bspline5の定義はwij=w(|p-a|)としている．
    */
    //人工粘性
    //まず，nu*laplacian(U)を計算する
    Tddd ret = {0, 0, 0}, Uij, rij, grad;
    double Cs = 1466., PIij, div_U, dot_rij_rij, dot_Uij_rij, m;
    for (const auto &p : ps)
        if (p != a)
            ret += laplacian_U_Monaghan1992(p, p->getXtuple(), p->U_SPH, a, h, alpha, beta);
    double nu = a->mu_SPH / a->density;
    return -ret / nu;
};

// b!New
Tddd laplacian_U_ShaoAndLo2003(const networkPoint *a, const networkPoint *p, double h, const Tddd &X)
{
    auto Uij = a->U_SPH - p->U_SPH;
    auto Xij = a->getXtuple() - p->getXtuple();
    auto nu_i = p->mu_SPH / p->density;
    auto nu_j = a->mu_SPH / a->density;
    auto W = grad_kernel_Bspline5(a->getXtuple(),
                                  p->getXtuple(), h);
    return (p->mass * 8 * (nu_i + nu_j) * Dot(Uij, Xij) * W) / ((p->density + a->density) * Dot(Xij, Xij));
};
Tddd laplacian_U_ShaoAndLo2003(const networkPoint *a, const networkPoint *p, double h, const Tddd &X, const Tddd &U)
{
    auto Uij = a->U_SPH - U;
    auto Xij = a->getXtuple() - p->getXtuple();
    auto nu_i = p->mu_SPH / p->density;
    auto nu_j = a->mu_SPH / a->density;
    auto W = grad_kernel_Bspline5(a->getXtuple(),
                                  p->getXtuple(), h);
    return (p->mass * 8 * (nu_i + nu_j) * Dot(Uij, Xij) * W) / ((p->density + a->density) * Dot(Xij, Xij));
};
Tddd laplacian_U_ShaoAndLo2003_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                                const networkPoint *const a,
                                                const double h,
                                                const int order = 2)
{
    // nu*laplacian(U)を計算する
    Tddd ret = {0, 0, 0}, tmp = {0, 0, 0}, U, n, velocity;
    double nu_i, nu_j;
    auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
    for (const auto &p : ps)
    {
        if (p != a)
            ret += laplacian_U_ShaoAndLo2003(a, p, h, p->getXtuple());
        // b$ ------------------------ ポリゴン ------------------------ */
        for (const auto &[F0, F1, X, Y] : INTXN.getFFXY(p, h))
        {
            if (F1)
            {
                std::get<0>(velocity) = std::get<0>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
                std::get<1>(velocity) = std::get<1>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
                std::get<2>(velocity) = std::get<2>(F0->getNetwork()->velocity + F1->getNetwork()->velocity);
            }
            else
            {
                std::get<0>(velocity) = std::get<0>(F0->getNetwork()->velocity);
                std::get<1>(velocity) = std::get<1>(F0->getNetwork()->velocity);
                std::get<2>(velocity) = std::get<2>(F0->getNetwork()->velocity);
            }

#if defined(free_slip_boundary_condition)
            U = p->U_SPH - 2 * n * Dot(p->U_SPH, n) + Dot(velocity, n) * n;
#elif defined(no_slip_boundary_condition)
            U = -p->U_SPH + Dot(velocity, n) * n;
#else
            U = -n * Dot(p->U_SPH, n) + Dot(velocity, n) * n;
#endif
            ret += laplacian_U_ShaoAndLo2003(a, p, h, Y, U)
#if defined(boundary_gurd)
                   * boundary_gurd_value(Norm(Y - p->getXtuple()), h);
#else
                ;
#endif
        }
        // b$ ------------------------------------------------------ */
    }
    return ret / (a->mu_SPH / a->density);
};

Tddd laplacian_U_ShaoAndLo2003(const std::unordered_set<networkPoint *> &ps,
                               const networkPoint *const a,
                               const double h,
                               const int order = 2)
{
    // nu*laplacian(U)を計算する
    Tddd ret = {0, 0, 0}, tmp = {0, 0, 0}, Uij, Xij;
    double nu_i, nu_j;
    for (const auto &p : ps)
        if (p != a)
        {
            Uij = a->U_SPH - p->U_SPH;
            Xij = a->getXtuple() - p->getXtuple();
            nu_i = p->mu_SPH / p->density;
            nu_j = a->mu_SPH / a->density;
            tmp = p->mass * 8 * (nu_i + nu_j) * Dot(Uij, Xij) * grad_kernel_Bspline5(a->getXtuple(), p->getXtuple(), h);
            tmp /= ((p->density + a->density) * Dot(Xij, Xij));
            ret += tmp;
        }
    return ret / (a->mu_SPH / a->density);
};

/* ------------------------------------------------------ */
double pressure_EISPH_Hosseini2007_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                                    const networkPoint *const i,
                                                    const double h,
                                                    const double dt)
{
    if (dt < 1E-40)
        return i->pressure_SPH;
    // [1] Nomeritae, E. Daly, S. Grimaldi, and H. H. Bui, “Explicit incompressible SPH algorithm for free-surface flow modelling: A comparison with weakly compressible schemes,” Adv. Water Resour., vol. 97, pp. 156–167, Nov. 2016.
    Tddd Xij;
    double tmp, total_Aij = 0, total_Aij_Pj = 0, Aij;
    auto INTXN = IntersectionsSphereTrianglesLines(i->getContactFaces());
    for (const auto &j : ps)
    {
        if (j != i)
        {
            Xij = i->getXtuple() - j->getXtuple();
            tmp = j->mass * 8. / std::pow(j->density + i->density, 2);
            Aij = tmp * Dot(Xij, grad_kernel_Bspline5(i->getXtuple(), j->getXtuple(), h)) / Dot(Xij, Xij);
            total_Aij += Aij;
            total_Aij_Pj += Aij * j->pressure_SPH;
        }

        // b$ ------------------------ ポリゴン ------------------------ */
        for (const auto &[X, Y] : INTXN.get(j, h))
        {
            Xij = i->getXtuple() - Y;
            tmp = j->mass * 8. / std::pow(j->density + i->density, 2);
            Aij = tmp * Dot(Xij, grad_kernel_Bspline5(i->getXtuple(), Y, h)) / Dot(Xij, Xij)
#if defined(boundary_gurd)
                  * boundary_gurd_value(Norm(Y - j->getXtuple()), h);
#else
                ;
#endif
            total_Aij += Aij;
            total_Aij_Pj += Aij * j->pressure_SPH;
        }
        // b$ ------------------------------------------------------ */
    }
    double Bi = i->div_U / dt;
    return (Bi + total_Aij_Pj) / total_Aij;
};

double pressure_EISPH_Hosseini2007(const std::unordered_set<networkPoint *> &ps,
                                   const networkPoint *const i,
                                   const double h,
                                   const double dt)
{
    if (dt < 1E-40)
        return i->pressure_SPH;
    // [1] Nomeritae, E. Daly, S. Grimaldi, and H. H. Bui, “Explicit incompressible SPH algorithm for free-surface flow modelling: A comparison with weakly compressible schemes,” Adv. Water Resour., vol. 97, pp. 156–167, Nov. 2016.
    Tddd Xij;
    double tmp, total_Aij = 0, total_Aij_Pj = 0, Aij;
    for (const auto &j : ps)
        if (j != i)
        {
            Xij = i->getXtuple() - j->getXtuple();
            tmp = j->mass * 8. / std::pow(j->density + i->density, 2);
            Aij = tmp * Dot(Xij, grad_kernel_Bspline5(i->getXtuple(), j->getXtuple(), h)) / Dot(Xij, Xij);
            total_Aij += Aij;
            total_Aij_Pj += Aij * j->pressure_SPH;
        }
    double Bi = i->div_U / dt;
    return (Bi + total_Aij_Pj) / total_Aij;
};
/* ------------------------------------------------------ */
double DPDt_EISPH_Hosseini2007(const std::unordered_set<networkPoint *> &ps,
                               const networkPoint *const i,
                               const double h,
                               const double dt)
{
    if (dt < 1E-40)
        return 0.;
    else if (i->isSurface)
        return -i->pressure_SPH / dt;
    else
    {
        auto Pn1 = pressure_EISPH_Hosseini2007(ps, i, h, dt);
        return (Pn1 - i->pressure_SPH) / dt;
    }
};
/* ------------------------------------------------------ */

Tddd laplacian_Brookshaw1985(const std::unordered_set<networkPoint *> &ps,
                             networkPoint *a,
                             const double e,
                             const int order = 3)
{
    double tot = 0., tmp = 0.;
    Tddd ret = {0, 0, 0};
    Tddd r;
    for (const auto &p : ps)
        if (p != a)
        {
            r = p->getXtuple() - a->getXtuple();
            tmp = 2. * Dot(r / std::pow(Norm(r), 2), p->volume * grad_kernel_Bspline5(p->getXtuple(), a->getXtuple(), e));
            ret -= tmp * p->U_SPH;
            tot += tmp;
        }
    return ret + tot * a->U_SPH;
    // eq.(91) D. J. Price, “Smoothed particle hydrodynamics and magnetohydrodynamics,” J. Comput. Phys., vol. 231, no. 3, pp. 759–794, 2012.
};

//#endif