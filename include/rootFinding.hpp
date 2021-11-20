#ifndef rootFinding_H
#define rootFinding_H
#pragma once

#include "fundamental.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

template <typename T>
struct NewtonRaphson_Common
{
	T X;
	T dX; // tmp
	NewtonRaphson_Common(const T &Xinit) : X(Xinit), dX(Xinit){};
};

template <typename T>
struct NewtonRaphson : public NewtonRaphson_Common<T>
{
	NewtonRaphson(const T &Xinit) : NewtonRaphson_Common<T>(Xinit){};
};

template <>
struct NewtonRaphson<V_d> : public NewtonRaphson_Common<V_d>
{
	NewtonRaphson(const V_d &Xinit) : NewtonRaphson_Common<V_d>(Xinit){};
	void update(const V_d &F, const VV_d &dFdx)
	{
		ludcmp lu(dFdx);
		lu.solve(-F, dX);
		X += dX;
	};
};

template <>
struct NewtonRaphson<T4d> : public NewtonRaphson_Common<T4d>
{
	NewtonRaphson(const T4d &Xinit) : NewtonRaphson_Common<T4d>(Xinit)
	{
		this->ans = V_d{0., 0., 0., 0.};
	};
	V_d ans;
	void update(const T4d &F, const T4T4d &dFdx)
	{
		ludcmp lu(ToVector(dFdx));
		lu.solve(ToVector(-F), ans);
		std::get<0>(X) += ans[0];
		std::get<1>(X) += ans[1];
		std::get<2>(X) += ans[2];
		std::get<3>(X) += ans[3];
	};
};

#endif