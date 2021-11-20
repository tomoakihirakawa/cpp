#ifndef minMaxOfFunctions_H
#define minMaxOfFunctions_H
#pragma once

#include "fundamental.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

class GradientMethod
{
public:
	VV_d A;
	int count;
	GradientMethod(const VV_d &A_IN) : A(A_IN){};

	V_d solve(const V_d &b, const V_d &x_init = {}, double eps = 1E-5)
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
		double norm;
		V_d p = b - Dot(A, x); //修正ベクトルp
		while (!((norm = Norm(p)) < eps))
		{
			std::cout << "count = " << count << ", norm = " << norm << std::endl;
			p = b - Dot(A, x += Dot(p, p) / Dot(Dot(A, p), p) * p);
			if (count++ > 10000)
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
		std::cout << "count = " << count << ", norm = " << norm << std::endl;
		return x;
	};

	V_d solveCG(const V_d &b, const V_d &x_init = {}, double eps = 1E-5, int lim_count = 10000)
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
		x += alpha * p; //xを修正
		r = b - Dot(A, x /*x=x0+alpha*p0*/);
		//check
		double beta;
		while (!((norm = Norm(r)) < eps))
		{
			std::cout << "count = " << count << ", norm = " << norm << std::endl;
			beta = -Dot(r, Dot(A, p)) / Dot(Dot(A, p), p);
			p = r + beta * p; //修正ベクトルは，最急降下方向を少し変更した物
			alpha = Dot(r, p) / Dot(Dot(A, p), p);
			x += alpha * p;						 //xを修正
			r = b - Dot(A, x /*x=x0+alpha*p0*/); //残差ベクトルr（最急降下方向）
			if (count++ > lim_count)
				return x;
			// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
		std::cout << "count = " << count << ", norm = " << norm << std::endl;
		return x;
	};
};

#endif