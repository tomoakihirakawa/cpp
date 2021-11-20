#ifndef rootFinding_H
#define rootFinding_H
#pragma once

#include "fundamental.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

struct NewtonRaphson
{
	V_d X;
	V_d dX; // tmp
	NewtonRaphson(const V_d &Xinit) : X(Xinit), dX(Xinit){};
	void update(const V_d &F, const VV_d &dFdx)
	{
		ludcmp lu(dFdx);
		lu.solve(-F, dX);
		X += dX;
	};
	void update(const T6d &F, const T6dT6d &dFdx)
	{
		this->update(std::vector({std::get<0>(F), std::get<1>(F), std::get<2>(F) std::get<3>(F), std::get<4>(F), std::get<5>(F)}),
					 std::vector({std::get<0>(dFdx), std::get<1>(dFdx), std::get<2>(dFdx) std::get<3>(dFdx), std::get<4>(dFdx), std::get<5>(dFdx)}));
	};
};

#endif