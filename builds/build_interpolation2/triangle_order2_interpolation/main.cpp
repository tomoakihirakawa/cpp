#include "GNUPLOT.hpp"
#include "fundamental.hpp"

double func(const double x,
	    const double y,
	    const double z){
  return pow(x,3)*pow(y,2)*pow(z,1);
};

int main(){

  //三角形の場合は，変数変換をして，積分範囲をコンスタントにするとわかりやすい
  auto gws = GaussianQuadratureWeights(6,0,1);
  Print(gws,Blue);
  
  double ans=0;  
  for(const auto& x:gws)
    for(const auto& y:gws)
      ans += x[1]*y[1]*(1-x[0])*func(x[0],(1-x[0])*y[0],(1-x[0])*(1-y[0]));
  
  Print(ans, Red);
   
};

//Mathematica
// ClearAll["Global`*"];
// ijk = {3, 2, 1};

// f[x_, y_, z_, ijk_] := Module[{i, j, k},
//    {i, j, k} = ijk;
//    Return[x^i*y^j*z^k]];

// int = NIntegrate[f[x, y, 1 - x - y, ijk], {x, 0, 1}, {y, 0, 1 - x}]

// << NumericalDifferentialEquationAnalysis`;
// gw = GaussianQuadratureWeights[6, 0, 1]
// p = Transpose[gw][[1]];
// w = Transpose[gw][[2]];
// Sum[w[[j]]*w[[i]]*(1 - p[[i]])*
//   f[p[[i]], (1 - p[[i]]) p[[j]], (1 - p[[i]]) (1 - p[[j]]), ijk],
//  {i, 1, Length[w]},
//  {j, 1, Length[w]}]
