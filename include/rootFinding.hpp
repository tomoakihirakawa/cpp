#ifndef rootFinding_H
#define rootFinding_H
#pragma once

#include "fundamental.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

class NewtonRaphson {
   public:
	V_d X;
	V_d dX;  //tmp
	NewtonRaphson(const V_d &Xinit) : X(Xinit), dX(Xinit){};
	void update(const V_d &F, const VV_d &dFdx) {
		ludcmp lu(dFdx);
		lu.solve(-F, dX);
		X += dX;
	};
};

#endif