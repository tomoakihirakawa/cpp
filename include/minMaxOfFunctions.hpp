#ifndef minMaxOfFunctions_H
#define minMaxOfFunctions_H
#pragma once

#include "fundamental.hpp"
#include <functional>

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

struct GradientMethod
{
public:
	VV_d A;
	int count;
	GradientMethod(const VV_d &A_IN) : A(A_IN){};

	V_d solve(const V_d &b, const V_d &x_init = {}, double eps = 1E-5)
	{
		V_d x(b.size());
		for (auto i = 0; i < b.size(); ++i)
		{
			if (i >= x_init.size())
				x[i] = 0.;
			else if (!isFinite(x_init[i]))
				x[i] = 0.;
			else
				x[i] = x_init[i];
		}
		count = 0;
		double norm;
		V_d p = b - Dot(A, x); //修正ベクトルp
		norm = Norm(p);
		while (!(norm < eps))
		{
			std::cout << "count = " << count << ", norm = " << norm << std::endl;
			p = b - Dot(A, x += Dot(p, p) / Dot(Dot(A, p), p) * p);
			if (!isFinite(norm = Norm(p)))
				return x;
			else if (count++ > 10000)
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
		std::cout << "count = " << count << ", norm = " << norm << std::endl;
		return x;
	};

	/* ------------------------------------------------------ */

	V_d solveCG(const V_d &b, const V_d &x_init = {}, double eps = 1E-5)
	{
		V_d x(b.size());
		for (auto i = 0; i < b.size(); i++)
		{
			if (i >= x_init.size())
				x[i] = 0.;
			else if (!isFinite(x_init[i]))
				x[i] = 0.;
			else
				x[i] = x_init[i];
		}
		count = 0;
		V_d r = b - Dot(A, x); //残差ベクトル
		double norm = Norm(r);
		std::cout << "count = " << count << ", norm = " << norm << std::endl;
		if (norm < eps)
			return x;
		//
		V_d p = r;
		//
		double alpha = Dot(r, p) / Dot(Dot(A, p), p);
		x += alpha * p; // xを修正
		r = b - Dot(A, x /*x=x0+alpha*p0*/);
		// check
		double beta;
		norm = Norm(r);
		while (!(norm < eps))
		{
			std::cout << "count = " << count << ", norm = " << norm << std::endl;
			beta = -Dot(r, Dot(A, p)) / Dot(Dot(A, p), p);
			p = r + beta * p; //修正ベクトルは，最急降下方向を少し変更した物
			alpha = Dot(r, p) / Dot(Dot(A, p), p);
			x += alpha * p;						 // xを修正
			r = b - Dot(A, x /*x=x0+alpha*p0*/); //残差ベクトルr（最急降下方向）
			if (!isFinite(norm = Norm(r)))
				return x;
			else if (count++ > 10000)
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
		std::cout << "count = " << count << ", norm = " << norm << std::endl;
		return x;
	};
	/* ------------------------------------------------------ */
	GradientMethod(){};
	V_d minimize(const std::function<double(V_d)> &F,
				 const std::function<V_d(V_d)> &dF,
				 V_d Xinit,
				 double initial_step = 1E-12,
				 int max_loop = 100)
	{
		V_d X = Xinit, X_last = Xinit, m(X.size(), 0.);
		double s = initial_step;
		for (auto i = 0; i < max_loop; ++i)
		{
			X_last = X;
			double min = 1E+10;
			int size = 50;
			double d, f, S, extend = 2;
			for (auto k = 0; k < size; ++k)
			{
				d = 2. * extend * s / (size - 1);
				S = d * k - extend * s;
				f = F(X_last + S * m);
				if (k == 0 || min >= f)
				{
					min = f;
					s = S;
				}
			}
			X = X_last + s * m;
			/* ------------------------------------------------------ */
			auto dF_X = dF(X);
			auto dF_Xlast = dF(X_last);
			double a = -Dot(dF_X, dF_X - dF_Xlast) / Dot(dF_Xlast, dF_Xlast);
			// double a = -Dot(dF_X, dF_X - dF_Xlast) / Dot(m, dF_X - dF_Xlast);
			/* ------------------------------------------------------ */
			m = dF(X) + a * m;
			// std::cout << "s = " << s << std::endl;
			std::cout << "i = " << i << std::endl;
			std::cout << "X = " << X << std::endl;
			std::cout << "F(X) = " << F(X) << std::endl;
			std::cout << "Norm(m) = " << Norm(m) << std::endl;
			if (Norm(m) < 1E-8)
				break;
		};
		return X;
	};
};

#endif