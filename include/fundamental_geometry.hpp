#ifndef fundamental_geometry_H
#define fundamental_geometry_H
#pragma once

#include "fundamental_statistics.hpp"
#include "fundamental_vectors.hpp"

// 面と面の干渉？？？
// 球と面の干渉チェックか．
// そのために，
// タプルを作ろうとしている，

using Tdd = std::tuple<double, double>;
using Tddd = std::tuple<double, double, double>;
using T2Tddd = std::tuple<Tddd, Tddd>;
using T3Tddd = std::tuple<Tddd, Tddd, Tddd>;
using T3Tdd = std::tuple<Tdd, Tdd, Tdd>;

/*
M. Meyer, M. Desbrun, P. Schröder, and A. H. Barr, “Discrete Differential-Geometry Operators for Triangulated 2-Manifolds BT  - Visualization and Mathematics III,” Vis. Math. III, pp. 35–57, 2003.
*/
namespace geometry
{
	/* ------------------------------------------------------ */
	struct Line
	{
		T2Tddd X;
		Line(const T2Tddd &XIN) : X(XIN){};
		Line(const Tddd &X0IN, const Tddd &X1IN) : X({X0IN, X1IN}){};
	};
	/* ------------------------------------------------------ */
	struct Triangle
	{
		double thickness;
		Tddd center;
		T3Tddd X;
		Tddd normal;
		Triangle(const T3Tddd &XIN, double thicknessIN = 0)
			: X(XIN), thickness(thicknessIN), center(Mean(XIN))
		{
			normal = Normalize(Cross(std::get<1>(X) - std::get<0>(X), std::get<2>(X) - std::get<0>(X)));
		};
		Triangle(const Tddd &X0IN, const Tddd &X1IN, const Tddd &X2IN, double thicknessIN = 0)
			: X({X0IN, X1IN, X2IN}), thickness(thicknessIN), center((X0IN + X1IN) / 2.)
		{
			normal = Normalize(Cross(std::get<1>(X) - std::get<0>(X), std::get<2>(X) - std::get<0>(X)));
		};
	};
	/* ------------------------------------------------------ */
	struct Point
	{
		Tddd X;
		Point(const Tddd &XIN) : X(XIN){};
	};
	/* ------------------------------------------------------ */
	struct Sphere
	{
		double radius;
		Tddd X;
		Sphere(const Tddd &XIN, const double radiusIN = 0.) : X(XIN), radius(radiusIN){};
	};
	/* ------------------------------------------------------ */
	// structをわざわざ作るのは，T3Tddではなく，coordinateboundsとして意味を具体的にした状態で持ち回りたいから．それだけ．
	struct CoordinateBounds
	{
		T3Tdd bounds;
		const T3Tdd &operator()() const { return this->bounds; };
		CoordinateBounds() : bounds({{0, 0}, {0, 0}, {0, 0}}){};
		CoordinateBounds(const geometry::CoordinateBounds &bs) : bounds(bs.bounds){};
		CoordinateBounds(const Tddd &X)
		{
			this->bounds = {{std::get<0>(X), std::get<0>(X)},
							{std::get<1>(X), std::get<1>(X)},
							{std::get<2>(X), std::get<2>(X)}};
		};
		CoordinateBounds(const T3Tdd minmax)
		{
			this->bounds = minmax;
		};
		CoordinateBounds(const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ)
		{
			this->bounds = {{minX, maxX}, {minY, maxY}, {minZ, maxZ}};
		};
		CoordinateBounds(const T3Tddd &X)
		{
			// networkPointのsetBoundsで使われる
			auto [Xs, Ys, Zs] = Transpose(X);
			this->bounds = {{Min(Xs), Max(Xs)}, {Min(Ys), Max(Ys)}, {Min(Zs), Max(Zs)}};
		};
		CoordinateBounds(const T2Tddd &X)
		{
			auto [Xs, Ys, Zs] = Transpose(X);
			this->bounds = {{Min(Xs), Max(Xs)}, {Min(Ys), Max(Ys)}, {Min(Zs), Max(Zs)}};
		};
		// CoordinateBounds(const std::unordered_set<Tddd> &X)
		// {
		// 	double xmin = 1E+20, xmax = -1E+20, ymin = 1E+20, ymax = -1E+20, zmin = 1E+20, zmax = -1E+20;
		// 	for (const auto &xyz : X)
		// 	{
		// 		if (xmin >= std::get<0>(xyz))
		// 			xmin = std::get<0>(xyz);
		// 		if (xmax <= std::get<0>(xyz))
		// 			xmax = std::get<0>(xyz);

		// 		if (ymin >= std::get<1>(xyz))
		// 			ymin = std::get<1>(xyz);
		// 		if (ymax <= std::get<1>(xyz))
		// 			ymax = std::get<1>(xyz);

		// 		if (zmin >= std::get<2>(xyz))
		// 			zmin = std::get<2>(xyz);
		// 		if (zmax <= std::get<2>(xyz))
		// 			zmax = std::get<2>(xyz);
		// 	}
		// 	this->bounds = {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}};
		// };
		CoordinateBounds(const std::vector<Tddd> &X)
		{
			double xmin = 1E+20, xmax = -1E+20, ymin = 1E+20, ymax = -1E+20, zmin = 1E+20, zmax = -1E+20;
			for (const auto &xyz : X)
			{
				if (xmin >= std::get<0>(xyz))
					xmin = std::get<0>(xyz);
				if (xmax <= std::get<0>(xyz))
					xmax = std::get<0>(xyz);

				if (ymin >= std::get<1>(xyz))
					ymin = std::get<1>(xyz);
				if (ymax <= std::get<1>(xyz))
					ymax = std::get<1>(xyz);

				if (zmin >= std::get<2>(xyz))
					zmin = std::get<2>(xyz);
				if (zmax <= std::get<2>(xyz))
					zmax = std::get<2>(xyz);
			}
			this->bounds = {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}};
		};
		CoordinateBounds(const geometry::Line &L)
		{
			auto [Xs, Ys, Zs] = Transpose(L.X);
			this->bounds = {{Min(Xs), Max(Xs)}, {Min(Ys), Max(Ys)}, {Min(Zs), Max(Zs)}};
		};
		CoordinateBounds(const geometry::Sphere &S)
		{
			auto [x, y, z] = S.X;
			this->bounds = {{x - S.radius, x + S.radius},
							{y - S.radius, y + S.radius},
							{z - S.radius, z + S.radius}};
		};
		CoordinateBounds(const geometry::Triangle &T)
		{
			auto [Xs, Ys, Zs] = Transpose(T.X);
			this->bounds = {{Min(Xs), Max(Xs)}, {Min(Ys), Max(Ys)}, {Min(Zs), Max(Zs)}};
		};
		double getVolume() const
		{
			return (std::get<1>(std::get<0>(this->bounds)) - std::get<0>(std::get<0>(this->bounds))) *
				   (std::get<1>(std::get<1>(this->bounds)) - std::get<0>(std::get<1>(this->bounds))) *
				   (std::get<1>(std::get<2>(this->bounds)) - std::get<0>(std::get<2>(this->bounds)));
		};
		double getScale() const
		{
			return Norm(Tddd{std::get<1>(std::get<0>(this->bounds)) - std::get<0>(std::get<0>(this->bounds)),
							 std::get<1>(std::get<1>(this->bounds)) - std::get<0>(std::get<1>(this->bounds)),
							 std::get<1>(std::get<2>(this->bounds)) - std::get<0>(std::get<2>(this->bounds))});
		};
		Tddd getCenter() const
		{
			return Tddd{std::get<1>(std::get<0>(this->bounds)) + std::get<0>(std::get<0>(this->bounds)) / 2.,
						std::get<1>(std::get<1>(this->bounds)) + std::get<0>(std::get<1>(this->bounds)) / 2.,
						std::get<1>(std::get<2>(this->bounds)) + std::get<0>(std::get<2>(this->bounds)) / 2.};
		};
		bool isInside(const Tddd X) const
		{
			auto [xbounds, ybounds, zbounds] = this->bounds;
			if (std::get<0>(X) < std::get<0>(xbounds) || std::get<1>(xbounds) < std::get<0>(X))
				return false;
			else if (std::get<1>(X) < std::get<0>(ybounds) || std::get<1>(ybounds) < std::get<1>(X))
				return false;
			else if (std::get<2>(X) < std::get<0>(zbounds) || std::get<1>(zbounds) < std::get<2>(X))
				return false;
			else
				return true;
		};
	};
};

std::ostream &operator<<(std::ostream &stream, const geometry::CoordinateBounds &bounds) { return (stream << bounds.bounds); };
geometry::CoordinateBounds operator+(const geometry::CoordinateBounds &b0, const geometry::CoordinateBounds &b1)
{
	auto [minX0, maxX0] = std::get<0 /*x*/>(b0.bounds);
	auto [minX1, maxX1] = std::get<0 /*x*/>(b1.bounds);
	auto [minY0, maxY0] = std::get<1 /*y*/>(b0.bounds);
	auto [minY1, maxY1] = std::get<1 /*y*/>(b1.bounds);
	auto [minZ0, maxZ0] = std::get<2 /*z*/>(b0.bounds);
	auto [minZ1, maxZ1] = std::get<2 /*z*/>(b1.bounds);
	return geometry::CoordinateBounds((minX0 <= minX1 ? minX0 : minX1),
									  (maxX0 >= maxX1 ? maxX0 : maxX1),
									  (minY0 <= minY1 ? minY0 : minY1),
									  (maxY0 >= maxY1 ? maxY0 : maxY1),
									  (minZ0 <= minZ1 ? minZ0 : minZ1),
									  (maxZ0 >= maxZ1 ? maxZ0 : maxZ1));
};
bool IntersectQ(const geometry::CoordinateBounds &b0, const geometry::CoordinateBounds &b1)
{
	auto [b0x, b0y, b0z] = b0.bounds;
	auto [b1x, b1y, b1z] = b1.bounds;
	return !((std::get<0>(b0x) > std::get<0>(b1x) && std::get<0>(b0x) > std::get<1>(b1x) /*1のxの最大最小が，0のxの最小よりも小さい*/) ||
			 (std::get<0>(b0y) > std::get<0>(b1y) && std::get<0>(b0y) > std::get<1>(b1y) /*1のyの最大最小が，0のyの最小よりも小さい*/) ||
			 (std::get<0>(b0z) > std::get<0>(b1z) && std::get<0>(b0z) > std::get<1>(b1z) /*1のzの最大最小が，0のzの最小よりも小さい*/) ||
			 (std::get<1>(b0x) < std::get<0>(b1x) && std::get<1>(b0x) < std::get<1>(b1x) /*1のxの最大最小が，0のxの最大よりも大きい*/) ||
			 (std::get<1>(b0y) < std::get<0>(b1y) && std::get<1>(b0y) < std::get<1>(b1y) /*1のyの最大最小が，0のyの最大よりも大きい*/) ||
			 (std::get<1>(b0z) < std::get<0>(b1z) && std::get<1>(b0z) < std::get<1>(b1z) /*1のzの最大最小が，0のzの最大よりも大きい*/) /*これがtrueの場合，逆にhitなし*/);
};
struct IntersectionSphereLine
{
	double distance;
	Tddd X;
	bool isIntersecting;
	IntersectionSphereLine(const Tddd &center, double radius, const T2Tddd &AB)
		: isIntersecting(false)
	{
		auto p01 = std::get<0>(AB) - std::get<1>(AB);
		auto t = -Dot(std::get<1>(AB) - center, p01) / Dot(p01, p01);
		this->X = p01 * t + std::get<1>(AB);
		this->distance = Norm(this->X - center);
		if (this->distance <= radius /*may hit*/ && 0. <= t && t <= 1.)
			this->isIntersecting = true; //干渉する最も近い点は，線上にある
	}
};
/* ------------------------------------------------------ */
struct IntersectionTriangles
{
	/*
	チェック
		/Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/三角形と三角形の干渉.nb
	*/
	T3Tddd P0, P1;
	T2Tddd L;
	bool isIntersecting;
	double eps = 1E-10;
	double a0, b0, a1, b1, min0, max0, min1, max1, deno0, deno1;
	double I10InT, I11InT, I00InT, I01InT;
	Tddd I00, I01, I10, I11;
	double lmitMin(const double a0, const double b0) const
	{
		Tddd ret = {0., 0., 0.};
		if (std::abs(1 - a0) > eps)
		{
			if ((1 - a0) > 0)
				std::get<0>(ret) = -b0 / (1 - a0);
			else
				std::get<0>(ret) = (1 - b0) / (1 - a0);
		}

		if (std::abs(-a0) > eps)
		{
			if (-a0 > 0)
				std::get<1>(ret) = b0 / a0;
			else
				std::get<1>(ret) = -(1 - b0) / a0;
		}
		return Max(ret);
	};

	double lmitMax(const double a0, const double b0) const
	{
		Tddd ret = {1., 1., 1.};
		if (std::abs(1 - a0) > eps)
		{
			if ((1 - a0) > 0)
				std::get<0>(ret) = (1 - b0) / (1 - a0);
			else
				std::get<0>(ret) = -b0 / (1 - a0);
		}

		if (std::abs(-a0) > eps)
		{
			if (-a0 > 0)
				std::get<1>(ret) = -(1 - b0) / a0;
			else
				std::get<1>(ret) = b0 / a0;
		}
		return Min(ret); //この内最も小さいものが，上限となる
	};

	IntersectionTriangles(const T3Tddd &P0_IN, const T3Tddd &P1_IN) : P0(P0_IN), P1(P1_IN), isIntersecting(false)
	{
		// if (IntersectQ(geometry::CoordinateBounds(P0), geometry::CoordinateBounds(P0)))
		// 	return;
		auto [p00, p01, p02] = P0;
		auto [p10, p11, p12] = P1;
		auto normalP0 = Normalize(Cross(p01 - p00, p02 - p00));
		auto normalP1 = Normalize(Cross(p11 - p10, p12 - p10));
		if (std::abs(std::abs(Dot(normalP0, normalP1)) - 1) < eps)
		{
			// Print("干渉線を定義することはできない");
			// Print("線を共有するような場合はありえる");
			return;
		}
		deno0 = Dot(normalP1, p01 - p02);
		deno1 = Dot(normalP0, p11 - p12);

		if (std::abs(deno0) < eps)
			for (auto i = 1; i < 3; ++i)
			{
				P0 = RotateLeft(P0, i);
				p00 = std::get<0>(P0);
				p01 = std::get<1>(P0);
				p02 = std::get<2>(P0);
				deno0 = Dot(normalP1, p01 - p02);
				if (!(std::abs(deno0) < eps))
					break;
				if (i == 2)
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			}

		if (std::abs(deno1) < eps)
			for (auto i = 1; i < 3; ++i)
			{
				P1 = RotateLeft(P1, i);
				p10 = std::get<0>(P1);
				p11 = std::get<1>(P1);
				p12 = std::get<2>(P1);
				deno1 = Dot(normalP0, p11 - p12);
				if (!(std::abs(deno1) < eps))
					break;
				if (i == 2)
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
			}

		a0 = Dot(normalP1, p00 - p02) / deno0;
		b0 = Dot(normalP1, p12 - p02) / deno0;
		a1 = Dot(normalP0, p10 - p12) / deno1;
		b1 = Dot(normalP0, p02 - p12) / deno1;

		min0 = lmitMin(a0, b0);
		max0 = lmitMax(a0, b0);
		min1 = lmitMin(a1, b1);
		max1 = lmitMax(a1, b1);

		auto t = min0;
		I00 = Dot({t, b0 - a0 * t, 1 - t - (b0 - a0 * t)}, P0);
		t = max0;
		I01 = Dot({t, b0 - a0 * t, 1 - t - (b0 - a0 * t)}, P0);
		t = min1;
		I10 = Dot({t, b1 - a1 * t, 1 - t - (b1 - a1 * t)}, P1);
		t = max1;
		I11 = Dot({t, b1 - a1 * t, 1 - t - (b1 - a1 * t)}, P1);

		auto tmp = (I10 - I00) / (I01 - I00);
		I10InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
		auto I10IsInside = Between(I10InT, {0, 1});

		tmp = (I11 - I00) / (I01 - I00);
		I11InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
		auto I11IsInside = Between(I11InT, {0, 1});

		tmp = (I00 - I10) / (I11 - I10);
		I00InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
		auto I00IsInside = Between(I00InT, {0, 1});

		tmp = (I01 - I10) / (I11 - I10);
		I01InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
		auto I01IsInside = Between(I01InT, {0, 1});

		if (I10IsInside)
		{
			if (I11IsInside && Norm(I10 - I11) > eps)
			{
				L = {I10, I11};
				isIntersecting = true;
			}
			if (I01IsInside && Norm(I10 - I01) > eps)
			{
				L = {I10, I01};
				isIntersecting = true;
			}
			if (I00IsInside && Norm(I10 - I00) > eps)
			{
				L = {I10, I00};
				isIntersecting = true;
			}
		}
		if (I11IsInside)
		{
			if (I00IsInside && Norm(I11 - I00) > eps)
			{
				L = {I11, I00};
				isIntersecting = true;
			}
			if (I01IsInside && Norm(I11 - I01) > eps)
			{
				L = {I11, I01};
				isIntersecting = true;
			}
		}
		if (I00IsInside && I01IsInside && Norm(I00 - I01) > eps)
		{
			L = {I00, I01};
			isIntersecting = true;
		}
	};
};
/* ------------------------------------------------------ */
struct IntersectionSphereTriangle
{
	double scale, t0, t1;
	Tddd X;
	double eps = 1E-13;
	T3Tddd P;
	bool isIntersecting;
	Tddd center;
	IntersectionSphereTriangle(const Tddd &centerIN, const double radius, const T3Tddd &P_IN)
		: X({0, 0, 0}),
		  scale(0),
		  P(P_IN),
		  isIntersecting(false),
		  center(centerIN)
	{
		auto [p0, p1, p2] = P;
		auto n = Normalize(Cross(p1 - p0, p2 - p0));
		auto [nx, ny, nz] = n;
		auto [p02x, p02y, p02z] = p0 - p2;
		auto [p12x, p12y, p12z] = p1 - p2;
		double determ = -(nz * p02y * p12x) + ny * p02z * p12x + nz * p02x * p12y - nx * p02z * p12y - ny * p02x * p12z + nx * p02y * p12z;
		if (std::abs(determ) < eps)
			return;
		T3Tddd mat = {{nz * p12y - ny * p12z, -(nz * p02y) + ny * p02z, p02z * p12y - p02y * p12z},
					  {-(nz * p12x) + nx * p12z, nz * p02x - nx * p02z, -(p02z * p12x) + p02x * p12z},
					  {ny * p12x - nx * p12y, -(ny * p02x) + nx * p02y, p02y * p12x - p02x * p12y}};
		auto ans = Dot(center - p2, mat / determ);
		this->t0 = std::get<0>(ans);
		this->t1 = std::get<1>(ans);
		this->scale = std::get<2>(ans);
		this->X = scale * n + center;
		if (std::abs(scale) <= radius && (0 <= t0 && t0 <= 1) && (0 <= t1 && t1 <= 1) && (0 <= (1 - t0 - t1) && (1 - t0 - t1) <= 1))
			isIntersecting = true;
	};
	Tddd getNearestX() const
	{
		if (this->isIntersecting)
			return this->X;
		else
		{
			Tddd X_ = {1E+50, 1E+50, 1E+50};
			double mindistance = 1E+50;
			auto intxnL0 = IntersectionSphereLine(center, 1E+50, T2Tddd{std::get<0>(P), std::get<1>(P)});
			auto intxnL1 = IntersectionSphereLine(center, 1E+50, T2Tddd{std::get<1>(P), std::get<2>(P)});
			auto intxnL2 = IntersectionSphereLine(center, 1E+50, T2Tddd{std::get<2>(P), std::get<0>(P)});
			if (intxnL0.isIntersecting)
			{
				X_ = intxnL0.X;
				mindistance = Norm(X_ - center);
			}
			if (intxnL1.isIntersecting)
			{
				if (Norm(intxnL1.X - center) < mindistance)
				{
					X_ = intxnL1.X;
					mindistance = Norm(X_ - center);
				}
			}
			if (intxnL2.isIntersecting)
			{
				if (Norm(intxnL2.X - center) < mindistance)
				{
					X_ = intxnL2.X;
					mindistance = Norm(X_ - center);
				}
			}
			if (Norm(std::get<0>(P) - center) < mindistance)
			{
				X_ = std::get<0>(P);
				mindistance = Norm(X_ - center);
			}
			if (Norm(std::get<1>(P) - center) < mindistance)
			{
				X_ = std::get<1>(P);
				mindistance = Norm(X_ - center);
			}
			if (Norm(std::get<2>(P) - center) < mindistance)
			{
				X_ = std::get<2>(P);
			}
			return X_;
		}
	};
};

/* ------------------------------------------------------ */
const Tddd &Normal(const geometry::Triangle &triangle)
{
	return triangle.normal;
};
double NormalDistance(const geometry::Triangle &T, const Tddd &X) { return Norm(Dot(Normal(T), T.center - X)); };
Tddd vectorToTriangle(const geometry::Triangle &T, const Tddd &a)
{
	// aからTまでの最短ベクトル
	auto n = Normal(T);
	return n * Dot(n, std::get<0>(T.X) - a);
};
int IntersectQ(const geometry::Sphere &sphere, const geometry::Triangle &triangle)
{
	if (NormalDistance(triangle, sphere.X) < sphere.radius)
		return false;
	else
		return true;
};
/* ------------------------------------------------------ */
double VectorAngle(const Tddd &V1, const Tddd &V2) { return std::atan2(Norm(Cross(V1, V2)), Dot(V1, V2)); };
Tddd Angles(const geometry::Triangle &triangle)
{
	auto [a, b, c] = triangle.X;
	return std::make_tuple(VectorAngle(b - a, c - a), VectorAngle(c - b, a - b), VectorAngle(a - c, b - c));
};
double Area(const geometry::Triangle &triangle)
{
	auto [a, b, c] = triangle.X;
	auto A = Norm(a - c);
	auto B = Norm(b - a);
	auto C = Norm(c - b);
	auto s = 0.5 * (A + B + C);
	return std::sqrt(s * (s - A) * (s - B) * (s - C));
};
/* ------------------------------------------------------ */
double scalefactorToReach(const geometry::Line &line, const geometry::Triangle &triangle)
{
	auto [a, b] = line.X;
	auto [p0, p1, p2] = triangle.X;
	//オーダーが匹敵する物を選ぶ
	double log_b_a = log10(Norm(b - a));
	double diff0 = std::abs(log10(Norm(p0 - a) - log_b_a));
	double diff1 = std::abs(log10(Norm(p1 - a) - log_b_a));
	double diff2 = std::abs(log10(Norm(p2 - a) - log_b_a));
	Tddd n = Normal(triangle);
	if (diff0 < diff1 && diff0 < diff2)
		return Dot(p0 - a, n) / Dot(b - a, n);
	else if (diff1 < diff0 && diff1 < diff2)
		return Dot(p1 - a, n) / Dot(b - a, n);
	else
		return Dot(p2 - a, n) / Dot(b - a, n);
};
/* ------------------------------------------------------ */
int IntersectQ(const geometry::Line &line, const geometry::Triangle &triangle)
{
	auto [lA, lB] = line.X;
	auto [tA, tB, tC] = triangle.X;
	auto [lAx, lAy, lAz] = lA;
	auto [lBx, lBy, lBz] = lB;
	auto [tAx, tAy, tAz] = tA;
	auto [tBx, tBy, tBz] = tB;
	auto [tCx, tCy, tCz] = tC;
	if ((tAx > lAx && tBx > lAx && tCx > lAx && tAx > lBx && tBx > lBx && tCx > lBx) ||
		(tAx < lAx && tBx < lAx && tCx < lAx && tAx < lBx && tBx < lBx && tCx < lBx) ||
		(tAy > lAy && tBy > lAy && tCy > lAy && tAy > lBy && tBy > lBy && tCy > lBy) ||
		(tAy < lAy && tBy < lAy && tCy < lAy && tAy < lBy && tBy < lBy && tCy < lBy) ||
		(tAz > lAz && tBz > lAz && tCz > lAz && tAz > lBz && tBz > lBz && tCz > lBz) ||
		(tAz < lAz && tBz < lAz && tCz < lAz && tAz < lBz && tBz < lBz && tCz < lBz))
		return 0;
	double e = 1E-11;
	//これを1E-14とすることで，干渉のチェックが行われるようになる場合があった
	// auto d = factorOfVectorToReachTriangle(tA, tB, tC, a, b);
	auto d = scalefactorToReach(line, triangle);
	if (d < 0. || d > 1.)
		return 0; /*面に到達できていない*/

	auto b_a = lB - lA;
	auto n = Normal(triangle);
	auto ps = lA + (b_a)*d;

	//ポリゴン頂点の最大最小でチェック
	if (((tAx - e > std::get<0>(ps) && tBx - e > std::get<0>(ps) && tCx - e > std::get<0>(ps)) || (tAx + e < std::get<0>(ps) && tBx + e < std::get<0>(ps) && tCx + e < std::get<0>(ps))) ||
		((tAy - e > std::get<1>(ps) && tBy - e > std::get<1>(ps) && tCy - e > std::get<1>(ps)) || (tAy + e < std::get<1>(ps) && tBy + e < std::get<1>(ps) && tCy + e < std::get<1>(ps))) ||
		((tAz - e > std::get<2>(ps) && tBz - e > std::get<2>(ps) && tCz - e > std::get<2>(ps)) || (tAz + e < std::get<2>(ps) && tBz + e < std::get<2>(ps) && tCz + e < std::get<2>(ps))))
		return 1; /*面の最大最小範囲にすら入れていない*/

	auto ps_p0 = tA - ps;
	auto ps_p1 = tB - ps;
	auto ps_p2 = tC - ps;

	ps_p0 = ps_p0 / Norm(ps_p0);
	ps_p1 = ps_p1 / Norm(ps_p1);
	ps_p2 = ps_p2 / Norm(ps_p2);

	if (Dot(Cross(ps_p0, ps_p1), n) >= 0. &&
		Dot(Cross(ps_p1, ps_p2), n) >= 0. &&
		Dot(Cross(ps_p2, ps_p0), n) >= 0.)
		return 3; /*a,bは面と交差*/
	else
		return 2; /*a,bは面と交差していないが，かなり惜しい*/
};
/* ------------------------------------------------------ */
//@ 全接触情報を返す
struct intersection
{
	Tddd X; //接触した物体の最も近い座標
	double distance;
	bool isIntersecting;
	double eps = 1E-10;
	int index_intersection_type;
	/*
	0: not intersecting
	1: intersecting with a sphere
	2: intersecting with a line
	3: intersecting with a triangle
	*/
	/* ------------------------------------------------------ */
	intersection(const geometry::Sphere &sphere0, const geometry::Sphere &sphere1) : X({0, 0, 0}), distance(1E+40), isIntersecting(false)
	{
		if (IntersectQ(geometry::CoordinateBounds(sphere0), geometry::CoordinateBounds(sphere1)))
		{
			double R = Norm(std::get<1>(sphere0.X) - std::get<0>(sphere1.X));
			double overlap = (sphere0.radius + sphere1.radius - R) / 2.;
			this->X = sphere0.X + (sphere0.radius - overlap / 2.) * Normalize(sphere1.X - sphere0.X);
			this->distance = Norm(this->X - sphere0.X);
			this->isIntersecting = (overlap >= 0);
			this->index_intersection_type = 1;
			return;
		}
		else
		{
			this->distance = 1E+40;
			this->isIntersecting = false;
			this->index_intersection_type = 0;
			return;
		}
	};
	/* ------------------------------------------------------ */
	intersection(const geometry::Sphere &sphere, const geometry::Line &line) : X({0, 0, 0}), distance(1E+40), isIntersecting(false)
	{
		if (IntersectQ(geometry::CoordinateBounds(sphere), geometry::CoordinateBounds(line)))
		{
			Tddd X0X1 = std::get<1>(line.X) - std::get<0>(line.X);

			//
			// Dot(C - sphere.X, X0X1) = 0        (1)
			// C = std::get<0>(line.X) + X0X1 * t (2)
			//
			// Dot(std::get<0>(line.X) + X0X1 * t - sphere.X, X0X1) = 0
			// -> Dot(std::get<0>(line.X) - sphere.X, X0X1) + Dot(X0X1 * t , X0X1) = 0
			// -> Dot(X0X1 , X0X1) * t = - Dot(std::get<0>(line.X) - sphere.X, X0X1)
			// -> t = - Dot(std::get<0>(line.X) - sphere.X, X0X1)/Dot(X0X1 , X0X1)
			//

			double t = -Dot(std::get<0>(line.X) - sphere.X, X0X1) / Dot(X0X1, X0X1);
			this->X = std::get<0>(line.X) + X0X1 * t;
			this->distance = Norm(this->X - sphere.X);
			if (this->distance <= sphere.radius /*may hit*/)
			{
				if (-eps <= t && t <= 1 + eps)
				{
					//干渉する最も近い点は，線上にある
					this->isIntersecting = true;
					this->index_intersection_type = 2;
					return;
				}
				else
				{
					//干渉する最も近い点は，端点である可能性をチェック
					intersection intx0(sphere, geometry::Sphere(std::get<0>(line.X)));
					intersection intx1(sphere, geometry::Sphere(std::get<1>(line.X)));

					if (intx0.isIntersecting && intx1.distance > intx0.distance)
					{
						this->isIntersecting = true;
						this->distance = intx0.distance;
						this->X = intx0.X;
						this->index_intersection_type = 1;
						return;
					}
					else if (intx1.isIntersecting && intx0.distance > intx1.distance)
					{
						this->isIntersecting = true;
						this->distance = intx1.distance;
						this->X = intx1.X;
						this->index_intersection_type = 1;
						return;
					}
					else
					{
						this->distance = 1E+40;
						this->isIntersecting = false;
						this->index_intersection_type = 0;
						return;
					}
				}
			}
			else
			{
				this->distance = 1E+40;
				this->isIntersecting = false;
				this->index_intersection_type = 0;
				return;
			}
		}
		else
		{
			this->distance = 1E+40;
			this->isIntersecting = false;
			this->index_intersection_type = 0;
			return;
		}
	};
	/* ------------------------------------------------------ */
	intersection(const geometry::Sphere &sphere, const geometry::Triangle &triangle) : X({0, 0, 0}), distance(1E+40), isIntersecting(false)
	{
		// エラーの原因は初期値を返しているのかもしれない
		// まずは，coordinateboundsをチェックする．
		if (IntersectQ(geometry::CoordinateBounds(sphere), geometry::CoordinateBounds(triangle)))
		{
			auto [x0, y0, z0] = std::get<0>(triangle.X) - sphere.X;
			auto [x1, y1, z1] = std::get<1>(triangle.X) - sphere.X;
			auto [x2, y2, z2] = std::get<2>(triangle.X) - sphere.X;
			Tddd cross = {y2 * (z0 - z1) + y0 * (z1 - z2) + y1 * (-z0 + z2), x2 * (-z0 + z1) + x1 * (z0 - z2) + x0 * (-z1 + z2), x2 * (y0 - y1) + x0 * (y1 - y2) + x1 * (-y0 + y2)};
			double crossSquared = Dot(cross, cross);
			double t0 = ((y2 * z1 - y1 * z2) * (y1 * z0 - y2 * z0 - y0 * z1 + y2 * z1 + y0 * z2 - y1 * z2) + x1 * (x0 * (y1 - y2) * y2 + x0 * (z1 - z2) * z2 + x2 * (-2 * y1 * y2 + y0 * (y1 + y2) - 2 * z1 * z2 + z0 * (z1 + z2))) + pow(x2, 2) * (-(y0 * y1) - z0 * z1 + pow(y1, 2) + pow(z1, 2)) - x0 * x2 * (-(y1 * y2) - z1 * z2 + pow(y1, 2) + pow(z1, 2)) + pow(x1, 2) * (-(y0 * y2) - z0 * z2 + pow(y2, 2) + pow(z2, 2))) / crossSquared;
			double t1 = ((y0 - y2) * y2 * z0 * z1 + (y0 * y1 - 2 * y0 * y2 + y1 * y2) * z0 * z2 + x0 * x2 * (y0 * y1 - 2 * y0 * y2 + y1 * y2 + z0 * z1 - 2 * z0 * z2 + z1 * z2) - z1 * z2 * (-(y0 * y2) + pow(x0, 2) + pow(y0, 2)) + x1 * (x0 * (y0 - y2) * y2 + x0 * (z0 - z2) * z2 - x2 * (-(y0 * y2) + z0 * (z0 - z2) + pow(y0, 2))) + y2 * (-y1 + y2) * (pow(x0, 2) + pow(z0, 2)) + pow(x2, 2) * (-(y0 * y1) - z0 * z1 + pow(y0, 2) + pow(z0, 2)) + (-(y0 * y1) + pow(x0, 2) + pow(y0, 2)) * pow(z2, 2)) / crossSquared;
			double scale = (-(x2 * y1 * z0) + x1 * y2 * z0 + x2 * y0 * z1 - x0 * y2 * z1 - x1 * y0 * z2 + x0 * y1 * z2) / crossSquared;
			//結果
			this->distance = Norm(scale * cross);
			this->X = scale * cross + sphere.X;
			if (this->distance <= sphere.radius /*may hit*/)
			{
				if ((-eps <= t0 && t0 <= 1. + eps) && (-eps <= t1 && t1 <= 1. + eps) && (t0 + t1 <= 1. + eps))
				{ //干渉する最も近い点は，三角形の面内にある．
					this->isIntersecting = true;
					this->index_intersection_type = 3;
					return;
				}
				else
				{
					intersection intx0(sphere, geometry::Line(std::get<0>(triangle.X), std::get<1>(triangle.X)));
					intersection intx1(sphere, geometry::Line(std::get<1>(triangle.X), std::get<2>(triangle.X)));
					intersection intx2(sphere, geometry::Line(std::get<2>(triangle.X), std::get<0>(triangle.X)));
					if (intx1.distance >= intx0.distance && intx2.distance >= intx0.distance && intx0.isIntersecting)
					{
						this->isIntersecting = true;
						this->distance = intx0.distance;
						this->X = intx0.X;
						this->index_intersection_type = intx0.index_intersection_type;
						return;
					}
					else if (intx0.distance >= intx1.distance && intx2.distance >= intx1.distance && intx1.isIntersecting)
					{
						this->isIntersecting = true;
						this->distance = intx1.distance;
						this->X = intx1.X;
						this->index_intersection_type = intx1.index_intersection_type;
						return;
					}
					else if (intx1.distance >= intx2.distance && intx0.distance >= intx2.distance && intx2.isIntersecting)
					{
						this->isIntersecting = true;
						this->distance = intx2.distance;
						this->X = intx2.X;
						this->index_intersection_type = intx2.index_intersection_type;
						return;
					}
					else
					{
						this->distance = 1E+40;
						this->isIntersecting = false;
						this->index_intersection_type = 0;
					}
				}
			}
			else
			{
				this->distance = 1E+40;
				this->isIntersecting = false;
				this->index_intersection_type = 0;
				return;
			}
		}
		else
		{
			this->distance = 1E+40;
			this->isIntersecting = false;
			this->index_intersection_type = 0;
			return;
		}
	};
};
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/*namespace_geometry_detail
  namespace_geometry_detail*/
/*namespace_geometry_code*/
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using V_i = std::vector<int>;
using VV_i = std::vector<std::vector<int>>;
/* ------------------------------------------------------ */
double normalDirDistanceFromTriangle(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a)
{
	return Dot(TriangleNormal(p0, p1, p2), p0 - a);
};
double normalDirDistanceFromTriangle(const T3Tddd &ps, const Tddd &&a)
{
	return Dot(TriangleNormal(std::get<0>(ps), std::get<1>(ps), std::get<2>(ps)), std::get<0>(ps) - a);
};
V_d vectorToTriangle(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a)
{
	auto n = TriangleNormal(p0, p1, p2);
	return n * Dot(n, p0 - a);
};
// line a to b
double factorOfVectorToReachTriangle(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a, const V_d &b)
{
	//オーダーが匹敵する物を選ぶ
	double log_b_a = log10(Norm(b - a));
	double diff0 = std::abs(log10(Norm(p0 - a) - log_b_a));
	double diff1 = std::abs(log10(Norm(p1 - a) - log_b_a));
	double diff2 = std::abs(log10(Norm(p2 - a) - log_b_a));
	V_d n = TriangleNormal(p0, p1, p2);
	if (diff0 < diff1 && diff0 < diff2)
		return Dot(p0 - a, n) / Dot(b - a, n);
	else if (diff1 < diff0 && diff1 < diff2)
		return Dot(p1 - a, n) / Dot(b - a, n);
	else
		return Dot(p2 - a, n) / Dot(b - a, n);
};
///////////////////////////////////////////////////////////
int isPointingTriangle(const V_d &p0, const V_d &p1, const V_d &p2,
					   const V_d &a, const V_d &b)
{
	auto d = factorOfVectorToReachTriangle(p0, p1, p2, a, b);

	if (!std::isfinite(d))
		return false; // nan

	V_d ps = a + (b - a) * d;

	// double e = 1E-20;
	// //ポリゴン頂点の最大最小でチェック
	// if (((p0[0] - e > ps[0] && p1[0] - e > ps[0] && p2[0] - e > ps[0]) || (p0[0] + e < ps[0] && p1[0] + e < ps[0] && p2[0] + e < ps[0])) ||
	//     ((p0[1] - e > ps[1] && p1[1] - e > ps[1] && p2[1] - e > ps[1]) || (p0[1] + e < ps[1] && p1[1] + e < ps[1] && p2[1] + e < ps[1])) ||
	//     ((p0[2] - e > ps[2] && p1[2] - e > ps[2] && p2[2] - e > ps[2]) || (p0[2] + e < ps[2] && p1[2] + e < ps[2] && p2[2] + e < ps[2])))
	//   return false;/*面の最大最小範囲にすら入れていない*/

	auto ps_p0 = p0 - ps;
	auto ps_p1 = p1 - ps;
	auto ps_p2 = p2 - ps;

	ps_p0 = ps_p0 / Norm(ps_p0);
	ps_p1 = ps_p1 / Norm(ps_p1);
	ps_p2 = ps_p2 / Norm(ps_p2);

	V_d n = TriangleNormal(p0, p1, p2);
	if (Dot(Cross(ps_p0, ps_p1), n) >= 0. &&
		Dot(Cross(ps_p1, ps_p2), n) >= 0. &&
		Dot(Cross(ps_p2, ps_p0), n) >= 0.)
		return true;
	else
		return false;
};
/* ------------------------------------------------------ */
int isIntersectingSurface(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a, const V_d &b)
{
	/* 0:頂点の最大最小の範囲の外で，片方にa,bgがある */
	/* 1:拡大した面には入れているが，多角形の頂点の最大最小範囲にすら入れていない */
	/* 2:a,bは多角形の面と交差していないが，かなり惜しい */
	/* 3:a,bは多角形の面と交差 */
	// double e = 1E-11; //これを1E-14とすることで，干渉のチェックが行われるようになる場合があった．
	// if ((p0[0] - e > a[0] && p1[0] - e > a[0] && p2[0] - e > a[0] && p0[0] - e > b[0] && p1[0] - e > b[0] && p2[0] - e > b[0]) ||
	// 	(p0[0] + e < a[0] && p1[0] + e < a[0] && p2[0] + e < a[0] && p0[0] + e < b[0] && p1[0] + e < b[0] && p2[0] + e < b[0]) ||
	// 	(p0[1] - e > a[1] && p1[1] - e > a[1] && p2[1] - e > a[1] && p0[1] - e > b[1] && p1[1] - e > b[1] && p2[1] - e > b[1]) ||
	// 	(p0[1] + e < a[1] && p1[1] + e < a[1] && p2[1] + e < a[1] && p0[1] + e < b[1] && p1[1] + e < b[1] && p2[1] + e < b[1]) ||
	// 	(p0[2] - e > a[2] && p1[2] - e > a[2] && p2[2] - e > a[2] && p0[2] - e > b[2] && p1[2] - e > b[2] && p2[2] - e > b[2]) ||
	// 	(p0[2] + e < a[2] && p1[2] + e < a[2] && p2[2] + e < a[2] && p0[2] + e < b[2] && p1[2] + e < b[2] && p2[2] + e < b[2]))
	// 	return 0;

	if ((p0[0] > a[0] && p1[0] > a[0] && p2[0] > a[0] && p0[0] > b[0] && p1[0] > b[0] && p2[0] > b[0]) ||
		(p0[0] < a[0] && p1[0] < a[0] && p2[0] < a[0] && p0[0] < b[0] && p1[0] < b[0] && p2[0] < b[0]) ||
		(p0[1] > a[1] && p1[1] > a[1] && p2[1] > a[1] && p0[1] > b[1] && p1[1] > b[1] && p2[1] > b[1]) ||
		(p0[1] < a[1] && p1[1] < a[1] && p2[1] < a[1] && p0[1] < b[1] && p1[1] < b[1] && p2[1] < b[1]) ||
		(p0[2] > a[2] && p1[2] > a[2] && p2[2] > a[2] && p0[2] > b[2] && p1[2] > b[2] && p2[2] > b[2]) ||
		(p0[2] < a[2] && p1[2] < a[2] && p2[2] < a[2] && p0[2] < b[2] && p1[2] < b[2] && p2[2] < b[2]))
		return 0;
	double e = 1E-11; //これを1E-14とすることで，干渉のチェックが行われるようになる場合があった．

	auto d = factorOfVectorToReachTriangle(p0, p1, p2, a, b);
	if (d < 0. || d > 1.)
		return 0; /*面に到達できていない*/

	V_d b_a = b - a;
	V_d n = TriangleNormal(p0, p1, p2);
	V_d ps = a + (b_a)*d;

	//ポリゴン頂点の最大最小でチェック
	if (((p0[0] - e > ps[0] && p1[0] - e > ps[0] && p2[0] - e > ps[0]) || (p0[0] + e < ps[0] && p1[0] + e < ps[0] && p2[0] + e < ps[0])) ||
		((p0[1] - e > ps[1] && p1[1] - e > ps[1] && p2[1] - e > ps[1]) || (p0[1] + e < ps[1] && p1[1] + e < ps[1] && p2[1] + e < ps[1])) ||
		((p0[2] - e > ps[2] && p1[2] - e > ps[2] && p2[2] - e > ps[2]) || (p0[2] + e < ps[2] && p1[2] + e < ps[2] && p2[2] + e < ps[2])))
		return 1; /*面の最大最小範囲にすら入れていない*/

	auto ps_p0 = p0 - ps;
	auto ps_p1 = p1 - ps;
	auto ps_p2 = p2 - ps;

	ps_p0 = ps_p0 / Norm(ps_p0);
	ps_p1 = ps_p1 / Norm(ps_p1);
	ps_p2 = ps_p2 / Norm(ps_p2);

	if (Dot(Cross(ps_p0, ps_p1), n) >= 0. &&
		Dot(Cross(ps_p1, ps_p2), n) >= 0. &&
		Dot(Cross(ps_p2, ps_p0), n) >= 0.)
		return 3; /*a,bは面と交差*/
	else
		return 2; /*a,bは面と交差していないが，かなり惜しい*/
};

/////////////////////////////////
// V_d pOnSurface(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a, const V_d &b)
// {
//   V_d n = TriangleNormal(p0, p1, p2);
//   return a + (b - a) * Dot(/*tangential vector*/ Norm(p0 - a) < Norm(p1 - a) ? (p1 - a) : (p0 - a),
//                            /*normal vector*/ n) /
//                  Dot((b - a), n);
//   //return a + (b-a)*Dot(p0-a,n)/Dot(b-a,n);
// }
V_d pOnSurface(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a, const V_d &b)
{
	V_d n = TriangleNormal(p0, p1, p2), b_a = b - a;
	return a + b_a * Dot(p0 - a, n) / Dot(b_a, n); //分母が0の場合はあり得る
}

Tddd pOnSurfaceTuple(const Tddd &p0, const Tddd &p1, const Tddd &p2, const Tddd &a, const Tddd &b)
{
	Tddd n = TriangleNormal(p0, p1, p2), b_a = b - a;
	return a + b_a * Dot(p0 - a, n) / Dot(b_a, n); //分母が0の場合はあり得る
}

V_d pOnSurface(const VV_d &p0p1p2, const VV_d &ab)
{
	if (p0p1p2.size() != 3)
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point size = " + std::to_string(p0p1p2.size()));

	return pOnSurface(p0p1p2[0],
					  p0p1p2[1],
					  p0p1p2[2],
					  ab[0], ab[1]);
}

Tddd pOnSurfaceTuple(const T3Tddd &p0p1p2, const T2Tddd &ab)
{
	return pOnSurfaceTuple(std::get<0>(p0p1p2),
						   std::get<1>(p0p1p2),
						   std::get<2>(p0p1p2),
						   std::get<0>(ab),
						   std::get<1>(ab));
}

int isIntersectingSurface(const VV_d &p0p1p2, const VV_d &ab)
{
	if (p0p1p2.size() != 3)
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point size = " + std::to_string(p0p1p2.size()));

	return isIntersectingSurface(p0p1p2[0],
								 p0p1p2[1],
								 p0p1p2[2],
								 ab[0], ab[1]);
};
/* ------------------------------------------------------ */
class intersectionTriangleLine
{
public:
	V_d X;
	V_d normal;
	V_d vecA2X;
	V_d vecX2B;
	V_d vecX2B_;
	bool isIntersect;
	int indexOfTriangle;
	// Vnewの方向を変更
	// Xonの位置を計算
	intersectionTriangleLine(const V_d &p0,
							 const V_d &p1,
							 const V_d &p2,
							 const V_d &A,
							 const V_d &B)
		: X({}), normal({}), vecA2X({}), vecX2B({}), vecX2B_({}), isIntersect(false), indexOfTriangle(0)
	{
		if (isIntersectingSurface(p0, p1, p2, A, B) == 3 /*三角形と干渉した場合*/)
		{
			this->normal = TriangleNormal(p0, p1, p2);
			this->X = pOnSurface(p0, p1, p2, A, B);
			this->vecA2X = this->X - A;
			this->vecX2B = B - this->X;
			this->vecX2B_ = this->reflect(this->vecX2B, this->normal);
			this->isIntersect = true;
		}
	};
	////////////
	intersectionTriangleLine(const VVV_d &p0p1p2s,
							 const V_d &A,
							 const V_d &B,
							 const V_i &exceptIndices = {})
		: X({}), normal({}), vecA2X({}), vecX2B({}), vecX2B_({}), isIntersect(false), indexOfTriangle(0)
	{
		double closest_distFromWall = 1E+100;
		double normal_distA2X;
		for (auto i = 0; i < p0p1p2s.size(); i++)
		{
			if (!MemberQ(exceptIndices, i))
			{
				intersectionTriangleLine LT(p0p1p2s[i][0], p0p1p2s[i][1], p0p1p2s[i][2], A, B);
				normal_distA2X = Norm(Dot(LT.vecA2X, LT.normal));
				if (LT.isIntersect && isFinite(normal_distA2X) && normal_distA2X < closest_distFromWall)
				{
					closest_distFromWall = normal_distA2X;
					this->X = LT.X;
					this->normal = LT.normal;
					this->vecA2X = LT.vecA2X;
					this->vecX2B = LT.vecX2B;
					this->vecX2B_ = LT.vecX2B_;
					this->isIntersect = LT.isIntersect;
					this->indexOfTriangle = i;
				}
			}
		};
	};
	//////////////
	V_d reflectIfPossible(const V_d &v) const
	{
		if (this->isIntersect)
			return this->reflect(v, this->normal);
		else
			return v;
	};
	V_d reflect(const V_d &v, const V_d &n) const
	{
		return v - 2. * Dot(v, n) * n;
	};
};
/* ------------------------------------------------------ */
namespace geometry
{
	class Point_Line
	{
	public:
		double t; // parameter of v from p_line0
		V_d x;	  // coordinate of point
		V_d v;	  // vector of line
		double d2line;
		double d2line_segment;
		Point_Line(const V_d &p, const V_d &p_line0, const V_d &p_line1)
		{
			this->v = p_line1 - p_line0;
			this->t = -Dot(p_line0 - p, this->v) / Dot(this->v, this->v);
			this->x = p_line0 + this->t * this->v;
			this->d2line = Norm(this->x - p);

			if (this->t < 0)
				this->d2line_segment = Norm(p_line0 - p);
			else if (this->t > 1)
				this->d2line_segment = Norm(p_line1 - p);
			else
				this->d2line_segment = this->d2line;
		};
	};

	//基本的な操作なので，fundamentalに持ってきた
	// どうして，> だとうまくいくのか？
	// まがっているのか？
	// bool isConvexPolygon(const VV_d &ps, const V_d &normal)
	// {
	//   auto s = ps.size();
	//   if (s < 3)
	//     return false;

	//   V_d v0, v1;
	//   for (auto i = 0; i < s; i++)
	//   {
	//     v0 = ps[i] - ps[(s + i - 1) % s];
	//     v1 = ps[(s + i + 1) % s] - ps[i];
	//     if (MyVectorAngle(v0 /*基準*/, v1, normal) < 0. /*ccw*/)
	//       return false;
	//   }
	//   return true;
	// };
	//////////
	bool isConvexPolygon(const VV_d &ps)
	{
		auto s = ps.size();
		if (s < 3)
			return false;

		V_d v0 = *ps.begin() - *ps.rbegin();
		V_d v1 = *(ps.begin() + 1) - *ps.begin();
		V_d normal = Cross(v0, v1);
		normal = normal / Norm(normal);
		double angle = MyVectorAngle(v0 /*基準*/, v1, normal);
		if (!isFinite(normal))
			return false;

		for (auto i = 0; i < s; i++)
		{
			// 0->1 1->2
			angle = angle * MyVectorAngle(ps[i] - ps[(s + i - 1) % s] /*基準*/,
										  ps[(s + i + 1) % s] - ps[i], normal);
			if (!isFinite(angle))
				return false;
			if (angle < 0.)
				return false; //符号が変わったらfalse
		}
		return true;
	};
	/*ccw angle
	 *       *
	 *     / | \
	 *    *  |  *
	 *   /  \|/  \
	 *  *----*----*
	 *   \2 1|1 3/
	 *    \  |  /
	 *     \3|2/
	 *       *
	 */
	// angle_sets = {{a1,a2,a3},{b1,b2,b3},...}
	// bool isInConvexPolygon(const VV_d &angle_sets)
	// {
	//   auto s = angle_sets.size();
	//   for (auto i = 0; i < s + 1; i++)
	//     if (angle_sets[i % s][1] + angle_sets[(i + 1) % s][2] > M_PI)
	//       return false;
	//   return true;
	// };
	/////////////////
	// bool isConcavePolygon(const VV_d &ps, const V_d &normal) { return !geometry::isConvexPolygon(ps, normal); };
	bool isConcavePolygon(const VV_d &ps)
	{
		auto s = ps.size();
		if (s < 3)
			return false;

		V_d v0 = *ps.begin() - *ps.rbegin();
		V_d v1 = *(ps.begin() + 1) - *ps.begin();
		V_d normal = Cross(v0, v1);
		normal = normal / Norm(normal);
		double angle = MyVectorAngle(v0 /*基準*/, v1, normal);
		if (!isFinite(normal))
			return false;

		for (auto i = 0; i < s; i++)
		{
			// 0->1 1->2
			angle = angle * MyVectorAngle(ps[i] - ps[(s + i - 1) % s] /*基準*/,
										  ps[(s + i + 1) % s] - ps[i], normal);
			if (!isFinite(angle))
				return false;
			if (angle < 0.)
				return true; //符号が変わったらtrue
		}
		return false; //符号が変わらなかったのでfalse
	};

	//--------------------------
	class point
	{
	public:
		point(const V_d &xyz, const int i) : active(true), X(xyz), index(i){};
		int index;
		bool active;
		V_d X;
		double area;
		double angle;
	};
	using V_pp = std::vector<point *>;
	using VV_pp = std::vector<std::vector<point *>>;
	//--------------------------
	VV_d extractX(const V_pp &ps)
	{
		VV_d ret(0);
		for (const auto &p : ps)
			ret.emplace_back(p->X);
		return ret;
	};
	V_d extractAreas(const V_pp &ps)
	{
		V_d ret(0);
		for (const auto &p : ps)
			ret.emplace_back(p->area);
		return ret;
	};
	V_d extractAngle(const V_pp &ps)
	{
		V_d ret(0);
		for (const auto &p : ps)
			ret.emplace_back(p->angle);
		return ret;
	};
	V_i extractIndices(const V_pp &ps)
	{
		V_i ret(0);
		for (const auto &p : ps)
			ret.emplace_back(p->index);
		return ret;
	};
	VV_i extractIndices(const VV_pp &pps)
	{
		VV_i ret(0);
		for (const auto &ps : pps)
			ret.emplace_back(extractIndices(ps));
		return ret;
	};
	//-------------------------
	/*polygon_detail
	多角形クラス
	polygon_detail*/
	/*polygon_code*/
	class polygon
	{
	public:
		V_pp points;
		~polygon()
		{
			for (const auto &p : this->points)
				if (p)
					delete p;
		};
		polygon(const VV_d &xyz_IN)
		{
			int s = xyz_IN.size();
			if (s < 3)
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "this is not polygon: size < 3");

			this->points.resize(s);
			for (int i = 0; i < xyz_IN.size(); i++)
				this->points[i] = new point(xyz_IN[i], i);
		};

		void activateAllPoints()
		{
			for (const auto &p : points)
				p->active = true;
		};
		V_pp getActivePoints()
		{
			V_pp ret;
			for (const auto &p : this->points)
				if (p->active)
					ret.emplace_back(p);
			return ret;
		};
		V_pp getAllPoints() { return this->points; };

		void calculateArea(const V_pp &ps, const V_d &normal)
		{
			int s = ps.size();
			for (auto i = 0; i < s; i++)
				ps[i]->area = DirectedArea(ps[i]->X - ps[(s + i - 1) % s]->X, ps[(i + 1) % s]->X - ps[i]->X, normal);
		};

		V_d getExteriorAngles(const V_pp &ps /*ccw*/, const V_d &normal)
		{
			auto s = ps.size();
			V_d ret(s);
			for (auto i = 0; i < s; i++)
				ret[i] = MyVectorAngle(ps[i]->X - ps[(s + i - 1) % s]->X, ps[(i + 1) % s]->X - ps[i]->X, normal);
			/*前後の線からなる多角形の外角*/
			return ret;
		};

		V_d getInteriorAngles(const V_pp &ps /*ccw*/, const V_d &normal)
		{
			return M_PI - getExteriorAngles(ps, normal);
		};

		void calculateAngle(const V_pp &ps)
		{
			auto s = ps.size();
			for (auto i = 0; i < s; i++)
				ps[i]->angle = MyVectorAngle(ps[i]->X - ps[(s + i - 1) % s]->X, ps[(i + 1) % s]->X - ps[i]->X);
			/*前後の線からなる多角形の外角*/
		};

		bool isSmallAngle(const point *p0, const point *p1, const point *p2, const double smallangle)
		{
			if (std::abs(MyVectorAngle(p1->X - p0->X, p2->X - p1->X) /*前後の線からなる多角形の外角*/) < smallangle)
				return true; // too small
			return false;
		};

		/*
		 */

		bool getPointsMeetCondition(const V_pp &ps, const double smallangle_IN, point *&select_p, int &current_index)
		{
			bool found = false;
			int s = ps.size();
			//徐々にsmallangleの制限を弱くしていく
			for (int i = 0; i < s; i++)
				if (ps[i]->area > 0. && ps[i]->angle > smallangle_IN && myIsfinite(ps[i]->area) && myIsfinite(ps[i]->angle))
					if ((select_p == nullptr) /*first time*/ || ps[i]->area < select_p->area /*from 2nd time*/)
					{
						select_p = ps[i];
						current_index = i;
						found = true;
					}
			return found;
		};

		// bool getPointsMeetCondition(const V_pp &ps, const double smallangle_IN, point *&select_p, int &index)
		// {
		//   bool found = false;
		//   int s = ps.size();
		//   double smallangle;
		//   int max_try = 1000; //試行回数
		//   for (int k = 0; k < max_try; k++)
		//   {
		//     //徐々にsmallangleの制限を弱くしていく
		//     for (int i = 0; i < s; i++)
		//     {
		//       smallangle = smallangle_IN * (1. - (double)k / (max_try - 1.)); //smallangle_IN,smallangle_IN*0.999,,smallangle_IN*0.998
		//       if (ps[i]->area > 0. && ps[i]->angle > smallangle)
		//       {
		//         if ((select_p == nullptr) /*first time*/ || ps[i]->area < select_p->area /*from 2nd time*/)
		//         {
		//           select_p = ps[i];
		//           index = i;
		//           found = true;
		//         }
		//       }
		//     }
		//     if (found)
		//       break;
		//   }
		//   return found;
		// };

		bool isConvexPolygon(const V_pp &ps, const V_d &normal) const
		{
			auto s = ps.size();
			for (auto i = 0; i < s; i++)
				if (MyVectorAngle(ps[i]->X - ps[(s + i - 1) % s]->X, ps[(i + 1) % s]->X - ps[i]->X, normal) < 0.)
					return false; //外角がマイナス：時計回り
			return true;
			// for (const auto &angle : this->getExteriorAngles(ps, normal))
			//   if (angle < 0)
			//     return false;
			// return true;
		};

		///////////// 4 点の場合の特別な分割 ///////////////

		// V_i legalConnection4(const V_pp &ps /*4点*/, const V_d &normal) const
		// {
		//   //
		//   //    3 *----* 2
		//   //      |b   |
		//   //      |    |
		//   //      |    |
		//   //      |   a|
		//   //    0 *----* 1
		//   //

		//   auto A = TriangleAngles(ps[0]->X,ps[1]->X,ps[2]->X);
		//   auto B = TriangleAngles(ps[2]->X,ps[3]->X,ps[0]->X);
		//   auto AB = Join(A,B);

		//   auto C = TriangleAngles(ps[1]->X,ps[2]->X,ps[3]->X);
		//   auto D = TriangleAngles(ps[3]->X,ps[0]->X,ps[1]->X);
		//   auto CD = Join(C,D);

		//   if(Max(AB) > Max(CD))
		//     return {{0,1,2},{2,3,2}};
		//   else
		//     return {{1,2,3},{3,0,1}};

		//   double sum;
		//   if (myIsfinite(sum = MyVectorAngle(ps[2]->X - ps[1]->X, ps[0]->X - ps[1]->X) + MyVectorAngle(ps[0]->X - ps[3]->X, ps[2]->X - ps[3]->X)) || M_PI <= sum)
		//     return {3, 1};
		//   else if (myIsfinite(sum = MyVectorAngle(ps[3]->X - ps[2]->X, ps[1]->X - ps[2]->X) + MyVectorAngle(ps[1]->X - ps[0]->X, ps[3]->X - ps[0]->X)) || M_PI <= sum)
		//     return {0, 2};
		//   else
		//     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		// };

		VV_pp legalDivision4(const V_pp &ps /*4点*/, const V_d &normal) const
		{
			// V_i ind = legalConnection4(ps, normal); //与えられたpsのインデックスで，元のindexとは異なる
			try
			{
				// TriangleNormal(ps[ind[0]]->X, ps[ind[1]]->X, ps[(ind[1] + 1) % 4]->X);
				// TriangleNormal(ps[ind[1]]->X, ps[ind[0]]->X, ps[(ind[0] + 1) % 4]->X);

				// return extractIndices(VV_pp{{ps[ind[0]], ps[ind[1]], ps[(ind[1] + 1) % 4]},
				//                             {ps[ind[1]], ps[ind[0]], ps[(ind[0] + 1) % 4]}});

				auto A = TriangleAngles(ps[0]->X, ps[1]->X, ps[2]->X);
				auto B = TriangleAngles(ps[2]->X, ps[3]->X, ps[0]->X);
				auto AB = Join(A, B);

				auto C = TriangleAngles(ps[1]->X, ps[2]->X, ps[3]->X);
				auto D = TriangleAngles(ps[3]->X, ps[0]->X, ps[1]->X);
				auto CD = Join(C, D);

				if (!isfinite(AB) && !isfinite(CD))
					throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "both is not finite");
				else if (!isfinite(AB))
					return {{ps[1], ps[2], ps[3]}, {ps[3], ps[0], ps[1]}}; // CD
				else if (!isfinite(CD))
					return {{ps[0], ps[1], ps[2]}, {ps[2], ps[3], ps[0]}}; // AB
				else if (Max(AB) > Max(CD))
					return {{ps[1], ps[2], ps[3]}, {ps[3], ps[0], ps[1]}}; // CD
				else
					return {{ps[0], ps[1], ps[2]}, {ps[2], ps[3], ps[0]}}; // AB
			}
			catch (const error_message &e)
			{
				e.print();
				std::stringstream ss;
				//
				ss << MyVectorAngle(ps[2]->X - ps[1]->X, ps[0]->X - ps[1]->X) << std::endl;
				ss << MyVectorAngle(ps[0]->X - ps[3]->X, ps[2]->X - ps[3]->X) << std::endl;
				//
				ss << MyVectorAngle(ps[3]->X - ps[2]->X, ps[1]->X - ps[2]->X) << std::endl;
				ss << MyVectorAngle(ps[1]->X - ps[0]->X, ps[3]->X - ps[0]->X) << std::endl;
				//
				ss << MyVectorAngle(ps[1]->X - ps[0]->X, ps[2]->X - ps[0]->X) << std::endl;
				ss << MyVectorAngle(ps[3]->X - ps[1]->X, ps[0]->X - ps[1]->X) << std::endl;
				ss << MyVectorAngle(ps[0]->X - ps[2]->X, ps[1]->X - ps[2]->X) << std::endl;
				ss << MyVectorAngle(ps[1]->X - ps[3]->X, ps[2]->X - ps[3]->X) << std::endl;
				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
			}
		};
		///////////////////////////////////////////////////

		/*polygon::triangulate_detail
	  多角形の三角形分割．外角が`smallangle`よりも狭い三角形は対象にしない．
	  もし該当がなければ，`smallangle`を徐々に小さくしながら該当があるまで何回か繰り返す．
	  polygon::triangulate_detail*/
		/*polygon::triangulate_code*/
		// #define debug_triangle

		VV_i triangulate(const V_d &normal, double smallangle_IN = 0.)
		{
			//必ず1E-10よりも大きくなるように設定している
			// std::cout << __PRETTY_FUNCTION__ << std::endl;
			smallangle_IN += 1E-10;

			if (!isfinite(normal))
				throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "this->normal is not finite"));

			this->activateAllPoints();

			VV_i ret;
			V_pp ps = this->getActivePoints();

			if (ps.size() < 3)
				return {};
			else if (ps.size() == 3)
			{
				if (isfinite(TriangleNormal(ps[0]->X, ps[1]->X, ps[2]->X)))
					return {{ps[0]->index, ps[1]->index, ps[2]->index}};
				else
				{
					std::stringstream ss;
					auto a = ps[0]->X, b = ps[1]->X, c = ps[2]->X;
					ss << "{a,b,c} = " << VV_d{a, b, c} << std::endl;
					ss << "n = " << TriangleNormal(a, b, c) << std::endl;
					ss << "TriangleAngles = " << TriangleAngles(a, b, c) << std::endl;
					ss << "TriangleArea = " << TriangleArea(a, b, c) << std::endl;
					throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str()));
				}
			}
			else if (ps.size() == 4)
			{
				return extractIndices(legalDivision4(ps, normal));
			}
			else if (ps.size() > 4)
			{
				do
				{
#ifdef debug_triangle
					std::cout << "   ps = " << ps << std::endl;
					std::cout << "  ret = " << ret << std::endl;
#endif

					this->calculateArea(ps, normal);
					this->calculateAngle(ps);
					int s = ps.size();

					point *select_p = nullptr;
					int current_index = 0; //このpsのインデックスで元々与えられたベクトルのインデックスではない
					bool found = getPointsMeetCondition(ps, smallangle_IN, select_p, current_index);

#ifdef debug_triangle
					std::cout << "found = " << found << std::endl;
					std::cout << "   ps = " << ps << std::endl;
					std::cout << "  ret = " << ret << std::endl;
#endif

					if (found)
					{
						select_p->active = false;
						ret.push_back({ps[(s + current_index - 1) % s]->index, ps[current_index]->index, ps[(current_index + 1) % s]->index});
					}
					else
					{
						//   break;
						// return ret;
						throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not find triangle that meets conditions");
					}

					ps = this->getActivePoints();
				} while (ps.size() > 4);
			}
			VV_i two_trigs = extractIndices(legalDivision4(ps, normal));
			ret.insert(ret.end(), two_trigs.begin(), two_trigs.end());
			return ret;
		};

		// VV_i triangulate(const V_d &normal, double smallangle_IN = 0.)
		// {
		//   //必ず1E-10よりも大きくなるように設定している
		//   smallangle_IN += 1E-11;
		//   if (!isfinite(normal))
		//     throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "this->normal is not finite"));

		//   this->activateAllPoints();

		//   VV_i ret;
		//   V_pp ps = this->getActivePoints();

		//   if (ps.size() < 3)
		//     return {};
		//   else if (ps.size() == 3)
		//     return {{ps[0]->index, ps[1]->index, ps[2]->index}};
		//   else if (ps.size() > 3)
		//   {
		//     do
		//     {
		//       this->calculateArea(ps, normal);
		//       this->calculateAngle(ps);

		//       point *select_p = nullptr;
		//       int index = 0;
		//       bool found = getPointsMeetCondition(ps, smallangle_IN, select_p, index);
		//       int s = ps.size();
		//       if (found)
		//       {
		//         select_p->active = false;
		//         ret.push_back({ps[(s + index - 1) % s]->index, ps[index]->index, ps[(index + 1) % s]->index});
		//       }
		//       else
		//         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not find triangle that meets conditions");

		//       ps = this->getActivePoints();

		//     } while (ps.size() > 3);
		//   }

		//   ret.emplace_back(std::vector<int>{ps[0]->index, ps[1]->index, ps[2]->index});
		//   return ret;
		// };
		/*polygon::triangulate_code*/
	};

	//-------------------------------------------------

	class spherical_polygon : public geometry::polygon
	{
	public:
		V_pp points;
		V_d org;
		VV_d x_on_sphere;
		~spherical_polygon(){};
		spherical_polygon(const V_d &org_IN, const VV_d &xyz_IN) : geometry::polygon(xyz_IN), org(org_IN)
		{
			VV_d x_on_sphere({});
			for (const auto &x : xyz_IN)
				x_on_sphere.emplace_back((x - org_IN) / Norm(x - org_IN));
		};
		void calculateArea(const V_pp &ps /*, const V_d& normal*/)
		{
			int s = ps.size();
			for (auto i = 0; i < s; i++)
			{
				ps[i]->area = SphericalInteriorArea(this->org, ps[(s + i - 1) % s]->X,
													ps[(s + i) % s]->X,
													ps[(s + i + 1) % s]->X);
			}
		};
		void calculateAngle(const V_pp &ps)
		{
			int s = ps.size();
			for (auto i = 0; i < s; i++)
			{
				ps[i]->angle = SphericalInteriorAngle(this->org, ps[(s + i - 1) % s]->X,
													  ps[(s + i) % s]->X,
													  ps[(s + i + 1) % s]->X);
			}
		};
		/*polygon::triangulate_detail
	  多角形の三角形分割．外角が`smallangle`よりも狭い三角形は対象にしない．
	  もし該当がなければ，`smallangle`を徐々に小さくしながら該当があるまで何回か繰り返す．
	  polygon::triangulate_detail*/
		/*polygon::triangulate_code*/
		VV_i triangulate(/*const V_d& normal, */ double smallangle_IN = 0.)
		{
			//必ず1E-10よりも大きくなるように設定している
			smallangle_IN += 1E-10;
			// if(!isfinite(normal))
			//   throw(error_message(__FILE__,__PRETTY_FUNCTION__,__LINE__, "this->normal is not finite"));

			this->activateAllPoints();

			VV_i ret;
			V_pp ps = this->getActivePoints();

			if (ps.size() < 3)
			{
				return {};
			}
			else if (ps.size() == 3)
			{
				return {{ps[0]->index, ps[1]->index, ps[2]->index}};
			}
			else if (ps.size() > 3)
			{
				do
				{
					this->calculateArea(ps /*, normal*/);
					this->calculateAngle(ps);

					point *select_p = nullptr;
					int index = 0;
					bool found = getPointsMeetCondition(ps, smallangle_IN,
														select_p, index);
					int s = ps.size();

					if (found)
					{
						select_p->active = false;
						ret.push_back({ps[(s + index - 1) % s]->index,
									   ps[index]->index,
									   ps[(index + 1) % s]->index});
					}
					else
					{
						return ret;
					}

					ps = this->getActivePoints();

				} while (ps.size() > 3);
			}

			ret.emplace_back(std::vector<int>{ps[0]->index, ps[1]->index, ps[2]->index});
			return ret;
		};
		/*polygon::triangulate_code*/
	};
	/*polygon_code*/

	// double SolidAngle(const V_d &org, const V_d &p0, const V_d &p1, const V_d &p2)
	// {
	//   //ccw p0 -> p1 -> p2
	//   double cond = Dot(Cross(p1 - p0, p2 - p1), Mean(VV_d{p0- org, p1- org, p2- org})); //Dot(triangle's normal vector, org's view direction vector)
	//   if (cond > 1E-8)
	//   {
	//     return 4. * M_PI - (SphericalInteriorAngle(org, p0, p1, p2) + SphericalInteriorAngle(org, p2, p0, p1) + SphericalInteriorAngle(org, p1, p2, p0) - M_PI);
	//   }
	//   else if (cond < -1E-8)
	//   {
	//     return SphericalInteriorAngle(org, p0, p1, p2) + SphericalInteriorAngle(org, p2, p0, p1) + SphericalInteriorAngle(org, p1, p2, p0) - M_PI;
	//   }
	//   else
	//   {
	//     return 0.5;
	//   }
	// };

	double SolidAngle(const V_d &o, const V_d &A, const V_d &B, const V_d &C)
	{
		//   V_d center = Mean[{A, B, C}];
		V_d oA = (A - o);
		V_d oB = (B - o);
		V_d oC = (C - o);
		double c = MyVectorAngle(oA, oB);
		double a = MyVectorAngle(oB, oC);
		double b = MyVectorAngle(oC, oA);
		double s = (a + b + c) / 2.;
		double eps = 1E-10;
		if (Between(s, {M_PI - eps, M_PI + eps}))
			return 4. * M_PI / 2.;
		else
			return 4. * atan(sqrt(tan(s / 2) * tan((s - a) / 2.) * tan((s - b) / 2) * tan((s - c) / 2)));
	};

	double SolidAngle(const V_d &org, const VV_d &ps)
	{
		// ccw
		double ret(0);
		geometry::spherical_polygon poly(org, ps);
		for (const auto &ind : poly.triangulate(1E-5))
			ret += geometry::SolidAngle(org, ps[ind[0]], ps[ind[1]], ps[ind[2]]);
		return ret;
	};

	//--------------------------
	// //CCW p0p1p2 returns a positive value
	// double SolidAngle(const V_d& p, const V_d& p0, const V_d& p1, const V_d& p2){
	//   V_d r, cross = Cross(p1 - p0, p2 - p0);
	//   V_d unit_n = cross;
	//   double ret(0.);
	//   for(const auto& abw:__GWGW6__){
	//     r = p0*abw[0] + p1*abw[1] + p2*(1.-abw[0]-abw[1]) - p;
	//     ret += Dot(r,unit_n)/pow(Norm(r),3.) * abw[2];
	//   }
	//   return std::abs(ret);
	// };

	// double SolidAngle(const V_d& p, const VV_d& X){
	//   geometry::polygon poly(X);
	//   double ret(0);
	//   auto indices= poly.triangulate(Mean(X)-p);
	//   std::cout << indices << std::endl;
	//   for(const auto& ind:indices){
	//     ret += geometry::SolidAngle(p, X[ind[0]], X[ind[1]], X[ind[2]]);
	//   }
	//   return ret;
	// };

} // namespace geometry
/*namespace_geometry_code*/
#endif
