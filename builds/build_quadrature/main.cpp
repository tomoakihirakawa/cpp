#include "GNUPLOT.hpp"
#include "fundamental.hpp"

double func(const double x,
			const double y,
			const double z)
{
	return pow(x, 3) * pow(y, 2) * pow(z, 1);
};

int main()
{

	// /* ----------------------- 普通のガウス点 ---------------------- */

	for (auto i = 0; i < 15; i++)
	{
		VV_d gw = GaussianQuadratureWeights(i, 0., 1.);
		std::cout << Green << "const static std::vector<std::vector<double>> __GW" << std::setprecision(15) << std::to_string(i) << "__ = " << (gw) << ";" << std::endl;
	}

	// 2次元の積分の場合，変数変換に注意
	// for (auto i = 0; i < 15; i++)
	// {
	// 	VV_d gw = GaussianQuadratureWeights(i, 0., 1.);
	// 	VV_d gwgw;
	// 	for (const auto &tw0 : gw)
	// 	{
	// 		double xi = tw0[0];
	// 		double w_xi = tw0[1];
	// 		for (const auto &tw1 : gw)
	// 		{
	// 			double h = tw1[0];
	// 			double w_h = tw1[1];
	// 			gwgw.push_back(
	// 				{xi,
	// 				 h,
	// 				 w_xi * w_h /*1-xiなし*/});
	// 		}
	// 	}
	// 	std::cout << Green << "const static std::vector<std::vector<double>> __GWGW" << std::setprecision(15) << std::to_string(i) << "__ = " << gwgw << ";" << std::endl;
	// }

	// for (auto i = 0; i < 15; i++)
	// {
	// 	VV_d gw = GaussianQuadratureWeights(i, 0., 1.);
	// 	VV_d gwgw;
	// 	for (const auto &tw0 : gw)
	// 	{
	// 		double xi = tw0[0];
	// 		double w_xi = tw0[1];
	// 		for (const auto &tw1 : gw)
	// 		{
	// 			double h = tw1[0];
	// 			double w_h = tw1[1];
	// 			gwgw.push_back({xi /*a*/,
	// 							h * (1. - xi) /*b*/,
	// 							1 - (xi + h * (1. - xi)) /*1-(a+b)*/,
	// 							w_xi * w_h * (1. - xi) /*weight*/,
	// 							w_xi * w_h /*1-xiなし*/});
	// 		}
	// 	}
	// 	std::cout << Red << "const static std::vector<std::vector<double>> __GWGW" << std::setprecision(15) << std::to_string(i) << "__ = " << gwgw << ";" << std::endl;
	// }

	/* ------------------------ 特異積分様 ----------------------- */
	double x_sing = 1.;
	double beta = 2.;
	for (auto i = 0; i < 15; i++)
	{
		//!まず，[0,1]ではない範囲でガウスルジャンドルの点と重みを計算してもらう
		VV_d gw_singular = GaussianQuadratureWeights(i, InvSg(0., x_sing, beta), InvSg(1., x_sing, beta));
		VV_d gw = GaussianQuadratureWeights(i, 0., 1.);
		VV_d gwgw;
		for (const auto &tw0 : gw_singular)
		{
			double xi = Sg(tw0[0], x_sing, beta);			  //! 0,1に戻したガウス点
			double w_xi = tw0[1] * DSg(tw0[0], x_sing, beta); //! 0,1に戻した重み
			for (const auto &tw1 : gw)
			{
				double h = tw1[0];
				double w_h = tw1[1];
				gwgw.push_back({xi,
								h,
								w_xi * w_h});
			}
		}
		std::cout << Green << "const static std::vector<std::vector<double>> __GWGW_beta2_" << std::setprecision(15) << std::to_string(i) << "__ = " << gwgw << ";" << std::endl;
	}

	// Print("GWGWをみてみよう");
	// for(auto i=2; i<10; i++){
	//   auto gws = GaussianQuadratureWeights(i,0,1);
	//   std::cout << "__GWGW" + std::to_string(i) + "__" << gws << std::endl;
	// }

	//三角形の場合は，変数変換をして，積分範囲をコンスタントにするとわかりやすい
	auto gws = GaussianQuadratureWeights(6, 0, 1);
	Print(gws, Blue);

	double ans = 0;
	for (const auto &x : gws)
		for (const auto &y : gws)
			ans += x[1] * y[1] * (1 - x[0]) * func(x[0], (1 - x[0]) * y[0], (1 - x[0]) * (1 - y[0]));

	Print(ans, Red);
};

// Mathematica
//  ClearAll["Global`*"];
//  ijk = {3, 2, 1};

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
