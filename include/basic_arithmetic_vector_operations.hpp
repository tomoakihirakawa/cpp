#ifndef basic_arithmetic_vector_operations_H
#define basic_arithmetic_vector_operations_H

#include <algorithm> //transformなど
#include <numeric>	 //数値のシーケンスの処理に特化したアルゴリズム
#include <type_traits>
#include <vector>
#include <tuple>

using Tdd = std::tuple<double, double>;
using Tddd = std::tuple<double, double, double>;
using T6d = std::tuple<double, double, double, double, double, double>;
using T7d = std::tuple<double, double, double, double, double, double, double>;
using T2Tddd = std::tuple<Tddd, Tddd>;
using T3Tddd = std::tuple<Tddd, Tddd, Tddd>;
/* ------------------------------------------------------ */
/*                    タプルダブル7                         */
/* ------------------------------------------------------ */
T7d operator-(T7d v)
{
	std::get<0>(v) *= -1;
	std::get<1>(v) *= -1;
	std::get<2>(v) *= -1;
	std::get<3>(v) *= -1;
	std::get<4>(v) *= -1;
	std::get<5>(v) *= -1;
	std::get<6>(v) *= -1;
	return v;
};
/* ------------------------------------------------------ */
T7d &operator-=(T7d &v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	std::get<2>(v) -= u;
	std::get<3>(v) -= u;
	std::get<4>(v) -= u;
	std::get<5>(v) -= u;
	std::get<6>(v) -= u;
	return v;
};
T7d &operator+=(T7d &v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	std::get<4>(v) += u;
	std::get<5>(v) += u;
	std::get<6>(v) += u;
	return v;
};
T7d &operator*=(T7d &v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	std::get<4>(v) *= u;
	std::get<5>(v) *= u;
	std::get<6>(v) *= u;
	return v;
};
T7d &operator/=(T7d &v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	std::get<2>(v) /= u;
	std::get<3>(v) /= u;
	std::get<4>(v) /= u;
	std::get<5>(v) /= u;
	std::get<6>(v) /= u;
	return v;
}; /* ------------------------------------------------------ */
T7d &operator-=(T7d &v, const T7d &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	std::get<3>(v) -= std::get<3>(u);
	std::get<4>(v) -= std::get<4>(u);
	std::get<5>(v) -= std::get<5>(u);
	std::get<6>(v) -= std::get<6>(u);
	return v;
};
T7d &operator+=(T7d &v, const T7d &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	std::get<3>(v) += std::get<3>(u);
	std::get<4>(v) += std::get<4>(u);
	std::get<5>(v) += std::get<5>(u);
	std::get<6>(v) += std::get<6>(u);
	return v;
};
T7d &operator*=(T7d &v, const T7d &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	std::get<3>(v) *= std::get<3>(u);
	std::get<4>(v) *= std::get<4>(u);
	std::get<5>(v) *= std::get<5>(u);
	std::get<6>(v) *= std::get<6>(u);
	return v;
};
T7d &operator/=(T7d &v, const T7d &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	std::get<3>(v) /= std::get<3>(u);
	std::get<4>(v) /= std::get<4>(u);
	std::get<5>(v) /= std::get<5>(u);
	std::get<6>(v) /= std::get<6>(u);
	return v;
};
/* ------------------------------------------------------ */
T7d operator-(T7d v, const T7d &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	std::get<3>(v) -= std::get<3>(u);
	std::get<4>(v) -= std::get<4>(u);
	std::get<5>(v) -= std::get<5>(u);
	std::get<6>(v) -= std::get<6>(u);
	return v;
};
T7d operator+(T7d v, const T7d &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	std::get<3>(v) += std::get<3>(u);
	std::get<4>(v) += std::get<4>(u);
	std::get<5>(v) += std::get<5>(u);
	std::get<6>(v) += std::get<6>(u);
	return v;
};
T7d operator*(T7d v, const T7d &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	std::get<3>(v) *= std::get<3>(u);
	std::get<4>(v) *= std::get<4>(u);
	std::get<5>(v) *= std::get<5>(u);
	std::get<6>(v) *= std::get<6>(u);
	return v;
};
T7d operator/(T7d v, const T7d &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	std::get<3>(v) /= std::get<3>(u);
	std::get<4>(v) /= std::get<4>(u);
	std::get<5>(v) /= std::get<5>(u);
	std::get<6>(v) /= std::get<6>(u);
	return v;
};
/* ------------------------------------------------------ */
T7d operator-(T7d v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	std::get<2>(v) -= u;
	std::get<3>(v) -= u;
	std::get<4>(v) -= u;
	std::get<5>(v) -= u;
	std::get<6>(v) -= u;
	return v;
};
T7d operator+(T7d v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	std::get<4>(v) += u;
	std::get<5>(v) += u;
	std::get<6>(v) += u;
	return v;
};
T7d operator*(T7d v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	std::get<4>(v) *= u;
	std::get<5>(v) *= u;
	std::get<6>(v) *= u;
	return v;
};
T7d operator/(T7d v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	std::get<2>(v) /= u;
	std::get<3>(v) /= u;
	std::get<4>(v) /= u;
	std::get<5>(v) /= u;
	std::get<6>(v) /= u;
	return v;
};
/* ------------------------------------------------------ */
T7d operator-(const double u, T7d v)
{
	std::get<0>(v) = u - std::get<0>(v);
	std::get<1>(v) = u - std::get<1>(v);
	std::get<2>(v) = u - std::get<2>(v);
	std::get<3>(v) = u - std::get<3>(v);
	std::get<4>(v) = u - std::get<4>(v);
	std::get<5>(v) = u - std::get<5>(v);
	std::get<6>(v) = u - std::get<6>(v);
	return v;
};
T7d operator+(const double u, T7d v)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	std::get<4>(v) += u;
	std::get<5>(v) += u;
	std::get<6>(v) += u;
	return v;
};
T7d operator*(const double u, T7d v)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	std::get<4>(v) *= u;
	std::get<5>(v) *= u;
	std::get<6>(v) *= u;
	return v;
};
T7d operator/(const double u, T7d v)
{
	std::get<0>(v) = u / std::get<0>(v);
	std::get<1>(v) = u / std::get<1>(v);
	std::get<2>(v) = u / std::get<2>(v);
	std::get<3>(v) = u / std::get<3>(v);
	std::get<4>(v) = u / std::get<4>(v);
	std::get<5>(v) = u / std::get<5>(v);
	std::get<6>(v) = u / std::get<6>(v);
	return v;
};

/* ------------------------------------------------------ */
/*                    タプルダブル6                         */
/* ------------------------------------------------------ */
T6d operator-(T6d v)
{
	std::get<0>(v) *= -1;
	std::get<1>(v) *= -1;
	std::get<2>(v) *= -1;
	std::get<3>(v) *= -1;
	std::get<4>(v) *= -1;
	std::get<5>(v) *= -1;
	return v;
};
/* ------------------------------------------------------ */
T6d &operator-=(T6d &v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	std::get<2>(v) -= u;
	std::get<3>(v) -= u;
	std::get<4>(v) -= u;
	std::get<5>(v) -= u;
	return v;
};
T6d &operator+=(T6d &v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	std::get<4>(v) += u;
	std::get<5>(v) += u;
	return v;
};
T6d &operator*=(T6d &v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	std::get<4>(v) *= u;
	std::get<5>(v) *= u;
	return v;
};
T6d &operator/=(T6d &v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	std::get<2>(v) /= u;
	std::get<3>(v) /= u;
	std::get<4>(v) /= u;
	std::get<5>(v) /= u;
	return v;
}; /* ------------------------------------------------------ */
T6d &operator-=(T6d &v, const T6d &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	std::get<3>(v) -= std::get<3>(u);
	std::get<4>(v) -= std::get<4>(u);
	std::get<5>(v) -= std::get<5>(u);
	return v;
};
T6d &operator+=(T6d &v, const T6d &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	std::get<3>(v) += std::get<3>(u);
	std::get<4>(v) += std::get<4>(u);
	std::get<5>(v) += std::get<5>(u);
	return v;
};
T6d &operator*=(T6d &v, const T6d &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	std::get<3>(v) *= std::get<3>(u);
	std::get<4>(v) *= std::get<4>(u);
	std::get<5>(v) *= std::get<5>(u);
	return v;
};
T6d &operator/=(T6d &v, const T6d &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	std::get<3>(v) /= std::get<3>(u);
	std::get<4>(v) /= std::get<4>(u);
	std::get<5>(v) /= std::get<5>(u);
	return v;
};
/* ------------------------------------------------------ */
T6d operator-(T6d v, const T6d &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	std::get<3>(v) -= std::get<3>(u);
	std::get<4>(v) -= std::get<4>(u);
	std::get<5>(v) -= std::get<5>(u);
	return v;
};
T6d operator+(T6d v, const T6d &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	std::get<3>(v) += std::get<3>(u);
	std::get<4>(v) += std::get<4>(u);
	std::get<5>(v) += std::get<5>(u);
	return v;
};
T6d operator*(T6d v, const T6d &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	std::get<3>(v) *= std::get<3>(u);
	std::get<4>(v) *= std::get<4>(u);
	std::get<5>(v) *= std::get<5>(u);
	return v;
};
T6d operator/(T6d v, const T6d &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	std::get<3>(v) /= std::get<3>(u);
	std::get<4>(v) /= std::get<4>(u);
	std::get<5>(v) /= std::get<5>(u);
	return v;
};
/* ------------------------------------------------------ */
T6d operator-(T6d v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	std::get<2>(v) -= u;
	std::get<3>(v) -= u;
	std::get<4>(v) -= u;
	std::get<5>(v) -= u;
	return v;
};
T6d operator+(T6d v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	std::get<4>(v) += u;
	std::get<5>(v) += u;
	return v;
};
T6d operator*(T6d v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	std::get<4>(v) *= u;
	std::get<5>(v) *= u;
	return v;
};
T6d operator/(T6d v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	std::get<2>(v) /= u;
	std::get<3>(v) /= u;
	std::get<4>(v) /= u;
	std::get<5>(v) /= u;
	return v;
};
/* ------------------------------------------------------ */
T6d operator-(const double u, T6d v)
{
	std::get<0>(v) = u - std::get<0>(v);
	std::get<1>(v) = u - std::get<1>(v);
	std::get<2>(v) = u - std::get<2>(v);
	std::get<3>(v) = u - std::get<3>(v);
	std::get<4>(v) = u - std::get<4>(v);
	std::get<5>(v) = u - std::get<5>(v);
	return v;
};
T6d operator+(const double u, T6d v)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	std::get<4>(v) += u;
	std::get<5>(v) += u;
	return v;
};
T6d operator*(const double u, T6d v)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	std::get<4>(v) *= u;
	std::get<5>(v) *= u;
	return v;
};
T6d operator/(const double u, T6d v)
{
	std::get<0>(v) = u / std::get<0>(v);
	std::get<1>(v) = u / std::get<1>(v);
	std::get<2>(v) = u / std::get<2>(v);
	std::get<3>(v) = u / std::get<3>(v);
	std::get<4>(v) = u / std::get<4>(v);
	std::get<5>(v) = u / std::get<5>(v);
	return v;
};
/* ------------------------------------------------------ */
/*                    タプルダブル4                         */
/* ------------------------------------------------------ */
T4d operator-(T4d v)
{
	std::get<0>(v) = -std::get<0>(v);
	std::get<1>(v) = -std::get<1>(v);
	std::get<2>(v) = -std::get<2>(v);
	std::get<3>(v) = -std::get<3>(v);
	return v;
};
T4d operator-(T4d v, const T4d &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	std::get<3>(v) -= std::get<3>(u);
	return v;
};
T4d operator+(T4d v, const T4d &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	std::get<3>(v) += std::get<3>(u);
	return v;
};
T4d operator*(T4d v, const T4d &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	std::get<3>(v) *= std::get<3>(u);
	return v;
};
T4d operator/(T4d v, const T4d &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	std::get<3>(v) /= std::get<3>(u);
	return v;
};
/* ------------------------------------------------------ */
T4d operator-(T4d v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	std::get<2>(v) -= u;
	std::get<3>(v) -= u;
	return v;
};
T4d operator+(T4d v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	return v;
};
T4d operator+(const double u, T4d v)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	return v;
};
T4d operator*(T4d v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	return v;
};
T4d operator*(const double u, T4d v)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	return v;
};
T4d operator/(const double u, T4d v)
{
	std::get<0>(v) = u / std::get<0>(v);
	std::get<1>(v) = u / std::get<1>(v);
	std::get<2>(v) = u / std::get<2>(v);
	std::get<3>(v) = u / std::get<3>(v);
	return v;
};
T4d operator/(T4d v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	std::get<2>(v) /= u;
	std::get<3>(v) /= u;
	return v;
};
/* ------------------------------------------------------ */
T4d &operator-=(T4d &v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	std::get<2>(v) -= u;
	std::get<3>(v) -= u;
	return v;
};
T4d &operator+=(T4d &v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	std::get<3>(v) += u;
	return v;
};
T4d &operator*=(T4d &v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	std::get<3>(v) *= u;
	return v;
};
T4d &operator/=(T4d &v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	std::get<2>(v) /= u;
	std::get<3>(v) /= u;
	return v;
}; /* ------------------------------------------------------ */
T4d &operator-=(T4d &v, const T4d &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	std::get<3>(v) -= std::get<3>(u);
	return v;
};
T4d &operator+=(T4d &v, const T4d &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	std::get<3>(v) += std::get<3>(u);
	return v;
};
T4d &operator*=(T4d &v, const T4d &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	std::get<3>(v) *= std::get<3>(u);
	return v;
};
T4d &operator/=(T4d &v, const T4d &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	std::get<3>(v) /= std::get<3>(u);
	return v;
};
/* ------------------------------------------------------ */
/*                    タプルダブル3                         */
/* ------------------------------------------------------ */
Tddd operator-(Tddd v)
{
	std::get<0>(v) = -std::get<0>(v);
	std::get<1>(v) = -std::get<1>(v);
	std::get<2>(v) = -std::get<2>(v);
	return v;
};
/* ------------------------------------------------------ */
Tddd &operator-=(Tddd &v, const Tddd &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	return v;
};
Tddd &operator-=(Tddd &v, const double d)
{
	std::get<0>(v) -= d;
	std::get<1>(v) -= d;
	std::get<2>(v) -= d;
	return v;
};
Tddd &operator+=(Tddd &v, const Tddd &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	return v;
};
Tddd &operator+=(Tddd &v, const double d)
{
	std::get<0>(v) += d;
	std::get<1>(v) += d;
	std::get<2>(v) += d;
	return v;
};
Tddd &operator*=(Tddd &v, const Tddd &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	return v;
};
Tddd &operator*=(Tddd &v, const double d)
{
	std::get<0>(v) *= d;
	std::get<1>(v) *= d;
	std::get<2>(v) *= d;
	return v;
};
Tddd &operator/=(Tddd &v, const Tddd &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	return v;
};
Tddd &operator/=(Tddd &v, const double d)
{
	std::get<0>(v) /= d;
	std::get<1>(v) /= d;
	std::get<2>(v) /= d;
	return v;
};
/* ------------------------------------------------------ */
Tddd operator-(Tddd v, const Tddd &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	return v;
};
Tddd operator+(Tddd v, const Tddd &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	return v;
};
Tddd operator*(Tddd v, const Tddd &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	return v;
};
Tddd operator/(Tddd v, const Tddd &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	return v;
};
/* ------------------------------------------------------ */
Tddd operator-(Tddd v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	std::get<2>(v) -= u;
	return v;
};
Tddd operator+(Tddd v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	return v;
};
Tddd operator*(Tddd v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	return v;
};
Tddd operator/(Tddd v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	std::get<2>(v) /= u;
	return v;
};
/* ------------------------------------------------------ */
Tddd operator-(const double u, Tddd v)
{
	std::get<0>(v) = u - std::get<0>(v);
	std::get<1>(v) = u - std::get<1>(v);
	std::get<2>(v) = u - std::get<2>(v);
	return v;
};
Tddd operator+(const double u, Tddd v)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	return v;
};
Tddd operator*(const double u, Tddd v)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	return v;
};
Tddd operator/(const double u, Tddd v)
{
	std::get<0>(v) = u / std::get<0>(v);
	std::get<1>(v) = u / std::get<1>(v);
	std::get<2>(v) = u / std::get<2>(v);
	return v;
};
/* ------------------------------------------------------ */
/*                    タプルダブル2                         */
/* ------------------------------------------------------ */
Tdd operator-(Tdd v)
{
	std::get<0>(v) *= -1.;
	std::get<1>(v) *= -1.;
	return v;
};
/* ------------------------------------------------------ */
Tdd &operator-=(Tdd &v, const Tdd &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	return v;
};
Tdd &operator-=(Tdd &v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	return v;
};
Tdd &operator+=(Tdd &v, const Tdd &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	return v;
};
Tdd &operator+=(Tdd &v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	return v;
};
Tdd &operator*=(Tdd &v, const Tdd &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	return v;
};
Tdd &operator*=(Tdd &v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	return v;
};
Tdd &operator/=(Tdd &v, const Tdd &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	return v;
};
Tdd &operator/=(Tdd &v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	return v;
};
/* ------------------------------------------------------ */
Tdd operator-(Tdd v, const Tdd &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	return v;
};
Tdd operator+(Tdd v, const Tdd &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	return v;
};
Tdd operator*(Tdd v, const Tdd &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	return v;
};
Tdd operator/(Tdd v, const Tdd &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	return v;
};
/* ------------------------------------------------------ */
Tdd operator-(Tdd v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	return v;
};
Tdd operator+(Tdd v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	return v;
};
Tdd operator*(Tdd v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	return v;
};
Tdd operator/(Tdd v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	return v;
};
/* ------------------------------------------------------ */
Tdd operator-(const double u, Tdd v)
{
	std::get<0>(v) = u - std::get<0>(v);
	std::get<1>(v) = u - std::get<1>(v);
	return v;
};
Tdd operator+(const double u, Tdd v)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	return v;
};
Tdd operator*(const double u, Tdd v)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	return v;
};
Tdd operator/(const double u, Tdd v)
{
	std::get<0>(v) = u / std::get<0>(v);
	std::get<1>(v) = u / std::get<1>(v);
	return v;
};
/* ------------------------------------------------------ */
/*                     行列T3Tddd                          */
/* ------------------------------------------------------ */
// Tddd operator-(Tddd v)
// {
// 	std::get<0>(v) = -std::get<0>(v);
// 	std::get<1>(v) = -std::get<1>(v);
// 	std::get<2>(v) = -std::get<2>(v);
// 	return v;
// };
// /* ------------------------------------------------------ */
T3Tddd &operator-=(T3Tddd &v, const T3Tddd &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	return v;
};
T3Tddd &operator-=(T3Tddd &v, const double d)
{
	std::get<0>(v) -= d;
	std::get<1>(v) -= d;
	std::get<2>(v) -= d;
	return v;
};
T3Tddd &operator+=(T3Tddd &v, const T3Tddd &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	return v;
};
T3Tddd &operator+=(T3Tddd &v, const double d)
{
	std::get<0>(v) += d;
	std::get<1>(v) += d;
	std::get<2>(v) += d;
	return v;
};
T3Tddd &operator*=(T3Tddd &v, const T3Tddd &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	return v;
};
T3Tddd &operator*=(T3Tddd &v, const double d)
{
	std::get<0>(v) *= d;
	std::get<1>(v) *= d;
	std::get<2>(v) *= d;
	return v;
};
T3Tddd &operator/=(T3Tddd &v, const T3Tddd &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	return v;
};
T3Tddd &operator/=(T3Tddd &v, const double d)
{
	std::get<0>(v) /= d;
	std::get<1>(v) /= d;
	std::get<2>(v) /= d;
	return v;
};
/* ------------------------------------------------------ */
T3Tddd operator-(T3Tddd v, const T3Tddd &u)
{
	std::get<0>(v) -= std::get<0>(u);
	std::get<1>(v) -= std::get<1>(u);
	std::get<2>(v) -= std::get<2>(u);
	return v;
};
T3Tddd operator+(T3Tddd v, const T3Tddd &u)
{
	std::get<0>(v) += std::get<0>(u);
	std::get<1>(v) += std::get<1>(u);
	std::get<2>(v) += std::get<2>(u);
	return v;
};
T3Tddd operator*(T3Tddd v, const T3Tddd &u)
{
	std::get<0>(v) *= std::get<0>(u);
	std::get<1>(v) *= std::get<1>(u);
	std::get<2>(v) *= std::get<2>(u);
	return v;
};
T3Tddd operator/(T3Tddd v, const T3Tddd &u)
{
	std::get<0>(v) /= std::get<0>(u);
	std::get<1>(v) /= std::get<1>(u);
	std::get<2>(v) /= std::get<2>(u);
	return v;
};
/* ------------------------------------------------------ */
T3Tddd operator-(T3Tddd v, const double u)
{
	std::get<0>(v) -= u;
	std::get<1>(v) -= u;
	std::get<2>(v) -= u;
	return v;
};
T3Tddd operator+(T3Tddd v, const double u)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	return v;
};
T3Tddd operator*(T3Tddd v, const double u)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	return v;
};
T3Tddd operator/(T3Tddd v, const double u)
{
	std::get<0>(v) /= u;
	std::get<1>(v) /= u;
	std::get<2>(v) /= u;
	return v;
};
/* ------------------------------------------------------ */
T3Tddd operator-(const double u, T3Tddd v)
{
	std::get<0>(v) = u - std::get<0>(v);
	std::get<1>(v) = u - std::get<1>(v);
	std::get<2>(v) = u - std::get<2>(v);
	return v;
};
T3Tddd operator+(const double u, T3Tddd v)
{
	std::get<0>(v) += u;
	std::get<1>(v) += u;
	std::get<2>(v) += u;
	return v;
};
T3Tddd operator*(const double u, T3Tddd v)
{
	std::get<0>(v) *= u;
	std::get<1>(v) *= u;
	std::get<2>(v) *= u;
	return v;
};
T3Tddd operator/(const double u, T3Tddd v)
{
	std::get<0>(v) = u / std::get<0>(v);
	std::get<1>(v) = u / std::get<1>(v);
	std::get<2>(v) = u / std::get<2>(v);
	return v;
};

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(std::vector<T> ret)
{
	std::transform(ret.begin(), ret.end(), ret.begin(), [](const T tmp)
				   { return -tmp; });
	return ret;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> ret)
{
	std::transform(ret.begin(), ret.end(), ret.begin(), [](const std::vector<T> &tmp)
				   { return -tmp; });
	return ret;
};
///////////////////////////////////////////////////
// vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator*(std::vector<T> v, const T din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const T c)
				   { return c * din; });
	return v;
};
//! matrix * scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator*(std::vector<std::vector<T>> v, const T din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &c)
				   { return c * din; });
	return v;
};
//@ scaler * vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator*(const T din, std::vector<T> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const T c)
				   { return c * din; });
	return v;
};
//! scaler * matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator*(const T din, std::vector<std::vector<T>> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &c)
				   { return c * din; });
	return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<std::vector<T>>> operator*(std::vector<std::vector<std::vector<T>>> v, const T &din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<std::vector<T>> &c)
				   { return c * din; });
	return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<std::vector<T>>> operator*(const T &din, std::vector<std::vector<std::vector<T>>> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<std::vector<T>> &c)
				   { return c * din; });
	return v;
};
//@ vector * vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator*(std::vector<T> v, const std::vector<T> &w)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const T a, const T b)
				   { return a * b; });
	return v;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator*=(std::vector<T> &v, const T &w) { return (v = v * w); };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator*=(std::vector<T> &v, const std::vector<T> &w) { return (v = v * w); };
// //vector x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>> operator*(const std::vector<T>& v, std::vector<std::vector<T>> w){
//   std::transform(w.begin(),w.end(),v.cbegin(),w.begin(),[](const T& a, const T& b){ return a*b; });
//   return w;
// };
// //matrix x vector
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator*(std::vector<std::vector<T>> w, const std::vector<T>& v){
//   std::transform(w.begin(),w.end(),v.cbegin(),w.begin(),[](const T& a, const T& b){ return a*b; });
//   return w;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>>& operator*=(std::vector<std::vector<T>>& v, const std::vector<T>& w){
//   return (v = v * w);
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>>& operator*=(std::vector<std::vector<T>>& v, const std::vector<std::vector<T>>& w){
//   return (v = v * w);
// };
/////////////////////////////////////////////////
//@vector / scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator/(std::vector<T> v, const T din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const T &tmp)
				   { return tmp / din; });
	return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator/(std::vector<std::vector<T>> v, const T din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp)
				   { return tmp / din; });
	return v;
};
// vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator/(const T din, std::vector<T> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const T &tmp)
				   { return din / tmp; });
	return v;
};
// matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator/(const T din, std::vector<std::vector<T>> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp)
				   { return din / tmp; });
	return v;
};
// vector x vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator/(std::vector<T> v, const std::vector<T> &w)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](T a, T b)
				   { return a / b; });
	return v;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator/=(std::vector<T> &v, const T &w) { return (v = v / w); };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator/=(std::vector<T> &v, const std::vector<T> &w) { return (v = v / w); };
// //vector x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator/(const std::vector<T>& w, const std::vector< std::vector<T> >& v){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return b/a; });
//   return ret;
// };
// //matrix x vector
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator/(const std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return a/b; });
//   return ret;
// };
// //matrix x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator/(const std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, std::vector<T> b){ return a/b; });
//   return ret;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> >& operator/=(std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   v = v / w;
//   return v;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> >& operator/=(std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   v = v / w;
//   return v;
// };
/////////////////////////////////////////////////
//@ vector - scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(std::vector<T> v, const T din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const T tmp)
				   { return tmp - din; });
	return v;
};
//! matrix - scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> v, const T &din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp)
				   { return tmp - din; });
	return v;
};
//@ scaler - vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(const T din, std::vector<T> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const T tmp)
				   { return din - tmp; });
	return v;
};
//! scaler - matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(const T din, std::vector<std::vector<T>> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp)
				   { return din - tmp; });
	return v;
};
//@ vector - vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator-(std::vector<T> v, const std::vector<T> &w)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const T a, const T b)
				   { return a - b; });
	return v;
};
//! matrix - matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> v, const std::vector<std::vector<T>> &w)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const std::vector<T> &b)
				   { return a - b; });
	return v;
};
//%matrix - vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(std::vector<std::vector<T>> v, const std::vector<T> &w)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b)
				   { return a - b; });
	return v;
};
//%vector - matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator-(const std::vector<T> &w, std::vector<std::vector<T>> v)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b)
				   { return b - a; }); //逆にする
	return v;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator-=(std::vector<T> &v, const T &w) { return (v = v - w); };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator-=(std::vector<T> &v, const std::vector<T> &w) { return (v = v - w); };
//! matrix -= matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator-=(std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w)
{
	return (v = v - w);
};
//%vector -= matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator-=(std::vector<T> &v, const std::vector<std::vector<T>> &w)
{
	return (v = v - w);
};
//%matrix -= vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator-=(std::vector<std::vector<T>> &v, const std::vector<T> &w)
{
	return (v = v - w);
};
// //vector x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator-(const std::vector<T>& w, const std::vector< std::vector<T> >& v){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return b-a; });
//   return ret;
// };
// //matrix x vector
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator-(const std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return a-b; });
//   return ret;
// };
// //matrix x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator-(const std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, std::vector<T> b){ return a-b; });
//   return ret;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> >& operator-=(std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   v = v - w;
//   return v;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> >& operator-=(std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   v = v - w;
//   return v;
// };
/////////////////////////////////////////////////
//@vector + scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator+(std::vector<T> v, const T &din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const T &tmp)
				   { return tmp + din; });
	return v;
};
//! matrix + scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(std::vector<std::vector<T>> v, const T &din)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp)
				   { return tmp + din; });
	return v;
};
//@scaler + vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator+(const T &din, std::vector<T> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const T &tmp)
				   { return tmp + din; });
	return v;
};
//! scaler + matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(const T &din, std::vector<std::vector<T>> v)
{
	std::transform(v.begin(), v.end(), v.begin(), [&din](const std::vector<T> &tmp)
				   { return tmp + din; });
	return v;
};
//@vector + vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> operator+(std::vector<T> v, const std::vector<T> &w)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const T &a, const T &b)
				   { return a + b; });
	return v;
};
//! matrix + matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(std::vector<std::vector<T>> v, const std::vector<std::vector<T>> &w)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const std::vector<T> &b)
				   { return a + b; });
	return v;
};
//%matrix + vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(std::vector<std::vector<T>> v, const std::vector<T> &w)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b)
				   { return a + b; });
	return v;
};
//%vector + matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> operator+(const std::vector<T> &w, std::vector<std::vector<T>> v)
{
	std::transform(v.begin(), v.end(), w.cbegin(), v.begin(), [](const std::vector<T> &a, const T &b)
				   { return a + b; });
	return v;
};
//@vector += scaler
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator+=(std::vector<T> &v, const T &w)
{
	return (v = v + w);
};
//@vector += vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> &operator+=(std::vector<T> &v, const std::vector<T> &w)
{
	return (v = v + w);
};
//! matrix += matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator+=(std::vector<std::vector<T>> &v, const std::vector<std::vector<T>> &w)
{
	return (v = v + w);
};
//%vector += matrix
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator+=(std::vector<T> &v, const std::vector<std::vector<T>> &w)
{
	return (v = v + w);
};
//%matrix += vector
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> &operator+=(std::vector<std::vector<T>> &v, const std::vector<T> &w)
{
	return (v = v + w);
};
#endif