#ifndef fundamental_vectors_H
#define fundamental_vectors_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <type_traits>
#include <vector>

#include "basic_arithmetic_vector_operations.hpp"
#include "fundamental_exception.hpp"
/* ------------------------------------------------------ */
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
/* ------------------------------------------------------ */
using Tdd = std::tuple<double, double>;
using Tddd = std::tuple<double, double, double>;
using T4d = std::tuple<double, double, double, double>;
using T6d = std::tuple<double, double, double, double, double, double>;
using T7d = std::tuple<double, double, double, double, double, double, double>;
using Tiii = std::tuple<int, int, int>;
using T2Tdd = std::tuple<Tdd, Tdd>;
using T3Tdd = std::tuple<Tdd, Tdd, Tdd>;
using T3Tddd = std::tuple<Tddd, Tddd, Tddd>;
using T4T4d = std::tuple<T4d, T4d, T4d, T4d>;
using T6Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T6T6d = std::tuple<T6d, T6d, T6d, T6d, T6d, T6d>;
using T7T7d = std::tuple<T7d, T7d, T7d, T7d, T7d, T7d, T7d>;
/* ------------------------------------------------------ */
// //vector x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>> operator+(const std::vector<T>& w, const std::vector<std::vector<T>>& v){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return b+a; });
//   return ret;
// };
// //matrix x vector
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator+(const std::vector< std::vector<T> >& v, const std::vector<T>& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, T b){ return a+b; });
//   return ret;
// };
// //matrix x matrix
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector< std::vector<T> > operator+(const std::vector< std::vector<T> >& v, const std::vector< std::vector<T> >& w){
//   std::vector< std::vector<T> > ret;
//   std::transform(v.cbegin(),v.cend(),w.cbegin(),std::back_inserter(ret),[](std::vector<T> a, std::vector<T> b){ return a+b; });
//   return ret;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>>& operator+=(std::vector<std::vector<T>>& v, const std::vector<T>& w){
//   v = v + w;
//   return v;
// };
// template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
// std::vector<std::vector<T>>& operator+=(std::vector<std::vector<T>>& v, const std::vector< std::vector<T> >& w){
//   v = v + w;
//   return v;
// };
/////////////////////
std::ostream &operator<<(std::ostream &stream, const std::vector<std::string> &v)
{
	stream << "{";
	if (!v.empty())
	{
		for (size_t i = 0; i < v.size() - 1; i++)
			stream << v[i] << ",";
		stream << *v.rbegin();
	}
	stream << "}";
	return stream;
};
template <typename T>
std::ostream &operator<<(std::ostream &stream, const std::vector<T *> &v)
{
	stream << "{";
	if (!v.empty())
	{
		for (size_t i = 0; i < v.size() - 1; i++)
			stream << v[i] << ",";
		stream << *v.rbegin();
	}
	stream << "}";
	return stream;
};

std::ostream &operator<<(std::ostream &stream, const std::vector<bool> &v)
{
	stream << "{";
	if (!v.empty())
	{
		for (size_t i = 0; i < v.size() - 1; i++)
			stream << v[i] << ",";
		stream << *v.rbegin();
	}
	stream << "}";
	return stream;
};

std::ostream &operator<<(std::ostream &stream, const std::vector<int> &v)
{
	stream << "{";
	if (!v.empty())
	{
		for (size_t i = 0; i < v.size() - 1; i++)
			stream << v[i] << ",";
		stream << *v.rbegin();
	}
	stream << "}";
	return stream;
};

std::ostream &operator<<(std::ostream &stream, const std::vector<double> &v)
{
	stream << "{";
	if (!v.empty())
	{
		for (size_t i = 0; i < v.size() - 1; i++)
			stream << v[i] << ",";
		stream << *v.rbegin();
	}
	stream << "}";
	return stream;
};

std::ostream &operator<<(std::ostream &stream, const std::vector<std::vector<int>> &v)
{
	stream << "{";
	if (!v.empty())
	{
		for (size_t i = 0; i < v.size() - 1; i++)
			stream << v[i] << ",";
		stream << *v.rbegin();
	}
	stream << "}";
	return stream;
};
std::ostream &operator<<(std::ostream &stream, const std::vector<std::vector<double>> &v)
{
	stream << "{";
	if (!v.empty())
	{
		for (size_t i = 0; i < v.size() - 1; i++)
			stream << v[i] << ",";
		stream << *v.rbegin();
	}
	stream << "}";
	return stream;
};
// std::stringstream &operator<<(std::stringstream &stream, const std::tuple<double, double, double> &v)
// {
// 	stream << "{" << std::get<0>(v) << "," << std::get<1>(v) << "," << std::get<2>(v) << "}";
// 	return stream;
// };
// std::stringstream &operator<<(std::stringstream &stream, const std::vector<std::tuple<double, double, double>> &v)
// {
// 	for (const auto &u : v)
// 		stream << "{" << u << "}";
// 	return stream;
// };
std::ostream &operator<<(std::ostream &stream, const Tdd &v)
{
	stream << "{" << std::get<0>(v) << "," << std::get<1>(v) << "}";
	return stream;
};
std::ostream &operator<<(std::ostream &stream, const Tddd &v)
{
	stream << "{" << std::get<0>(v) << "," << std::get<1>(v) << "," << std::get<2>(v) << "}";
	return stream;
};
std::ostream &operator<<(std::ostream &stream, const T4d &v)
{
	stream << "{" << std::get<0>(v) << "," << std::get<1>(v) << "," << std::get<2>(v) << "," << std::get<3>(v) << "}";
	return stream;
};
std::ostream &operator<<(std::ostream &stream, const T6d &v)
{
	stream << "{" << std::get<0>(v) << "," << std::get<1>(v) << "," << std::get<2>(v) << "," << std::get<3>(v) << "," << std::get<4>(v) << "," << std::get<5>(v) << "}";
	return stream;
};
std::ostream &operator<<(std::ostream &stream, const T3Tdd &v)
{
	auto [v0, v1, v2] = v;
	stream << "{{" << std::get<0>(v0) << "," << std::get<1>(v0) << "},";
	stream << "{" << std::get<0>(v1) << "," << std::get<1>(v1) << "},";
	stream << "{" << std::get<0>(v2) << "," << std::get<1>(v2) << "}}";
	return stream;
};
std::ostream &operator<<(std::ostream &stream, const T3Tddd &v)
{
	auto [v0, v1, v2] = v;
	stream << "{{" << std::get<0>(v0) << "," << std::get<1>(v0) << "," << std::get<2>(v0) << "},";
	stream << "{" << std::get<0>(v1) << "," << std::get<1>(v1) << "," << std::get<2>(v1) << "},";
	stream << "{" << std::get<0>(v2) << "," << std::get<1>(v2) << "," << std::get<2>(v2) << "}}";
	return stream;
};
// template <typename T>
// std::ostream &operator<<(std::ostream &stream, const std::vector<T *> &v)
// {
//   stream << "{";
//   if (!v.empty())
//   {
//     for (size_t i = 0; i < v.size() - 1; i++)
//       stream << *v[i] << ",";
//     stream << *(*v.rbegin());
//   }
//   stream << "}";
//   return stream;
// };
///////////////////////////////////////
template <typename T>
class Chain
{
public:
	int isComplete;
	int isReversed;
	int addedBack;
	int addedFront;
	std::vector<T *> chained;
	Chain() : chained({})
	{
		addedBack = false;
		addedFront = false;
		isComplete = false;
		isReversed = false;
	};
	////////////////
	void display()
	{
		std::cout << "--------------------------------------" << std::endl;
		std::cout << "isComplete = " << isComplete << std::endl;
		std::cout << "isReversed = " << isReversed << std::endl;
		std::cout << " addedBack = " << addedBack << std::endl;
		std::cout << " addedFront = " << addedFront << std::endl;
		std::cout << "   chained = " << chained << std::endl;
	};
	///////
	void checkComplete(bool oneside_match)
	{
		if (oneside_match)
		{
			// A--B--Aでは鎖とは言えない
			// A--B--C--AならOK
			if ((*chained.begin()) == (*chained.rbegin()) && chained.size() > 3)
			{
				isComplete = 2;
			}
			else
				isComplete = 1;
		}
		else
		{
			this->chained = {};
			isComplete = 0;
			isReversed = 0;
		}
	};
	/////
	void join_front(const std::vector<T *> &base, const std::vector<T *> &pc)
	{
		this->chained = base;
		bool oneside_match = false;
		if ((*chained.rbegin()) == (*pc.begin()))
		{ // right direction
			chained.insert(chained.end(), /*後尾に挿入*/ pc.begin() + 1, pc.end());
			oneside_match = true;
			isReversed = false;
			addedBack = true;
			addedFront = false;
		}
		else if (*(chained.rbegin()) == *pc.rbegin())
		{ // oppsite direction
			chained.insert(chained.end(), /*後尾に反転して挿入*/ pc.rbegin() + 1, pc.rend());
			oneside_match = true;
			isReversed = true;
			addedBack = true;
			addedFront = false;
		}
		// else if (*(chained.begin()) == *pc.rbegin())
		// { //right direction
		//   chained.insert(chained.begin(), /*前方に挿入*/ pc.begin(), pc.end() - 1);
		//   oneside_match = true;
		//   isReversed = false;
		//   addedBack = false;
		//   addedFront = true;
		// }
		// else if (*(chained.begin()) == *pc.begin())
		// { //oppsite direction
		//   chained.insert(chained.begin(), /*前方に反転して挿入*/ pc.rbegin(), pc.rend() - 1);
		//   oneside_match = true;
		//   isReversed = true;
		//   addedBack = false;
		//   addedFront = true;
		// }
		checkComplete(oneside_match);
	};
	////////////////
	void join_front_fix_order(const std::vector<T *> &base, const std::vector<T *> &pc)
	{
		this->chained = base;
		bool oneside_match = false;
		if ((*chained.rbegin()) == (*pc.begin()))
		{ // right direction
			chained.insert(chained.end(), /*後尾に挿入*/ pc.begin() + 1, pc.end());
			oneside_match = true;
			isReversed = false;
			addedBack = true;
			addedFront = false;
		}
		// else if (*(chained.rbegin()) == *pc.rbegin())
		// { //oppsite direction
		//   chained.insert(chained.end(), /*後尾に反転して挿入*/ pc.rbegin() + 1, pc.rend());
		//   oneside_match = true;
		//   isReversed = true;
		//   addedBack = true;
		//   addedFront = false;
		// }
		// else if (*(chained.begin()) == *pc.rbegin())
		// { //right direction
		//   chained.insert(chained.begin(), /*前方に挿入*/ pc.begin(), pc.end() - 1);
		//   oneside_match = true;
		//   isReversed = false;
		//   addedBack = false;
		//   addedFront = true;
		// }
		// else if (*(chained.begin()) == *pc.begin())
		// { //oppsite direction
		//   chained.insert(chained.begin(), /*前方に反転して挿入*/ pc.rbegin(), pc.rend() - 1);
		//   oneside_match = true;
		//   isReversed = true;
		//   addedBack = false;
		//   addedFront = true;
		// }
		checkComplete(oneside_match);
	};
	//前後ろに適当にくっつけていくのは，うまくいかないことがわかった，
	//なぜなら，cutlineどうしがくっついたりする場合がでるためだ，
	void join(const std::vector<T *> &base, const std::vector<T *> &pc)
	{
		this->chained = base;
		bool oneside_match = false;
		if ((*chained.rbegin()) == (*pc.begin()))
		{ // right direction
			chained.insert(chained.end(), /*後尾に挿入*/ pc.begin() + 1, pc.end());
			oneside_match = true;
			isReversed = false;
			addedBack = true;
			addedFront = false;
		}
		else if (*(chained.rbegin()) == *pc.rbegin())
		{ // oppsite direction
			chained.insert(chained.end(), /*後尾に反転して挿入*/ pc.rbegin() + 1, pc.rend());
			oneside_match = true;
			isReversed = true;
			addedBack = true;
			addedFront = false;
		}
		else if (*(chained.begin()) == *pc.rbegin())
		{ // right direction
			chained.insert(chained.begin(), /*前方に挿入*/ pc.begin(), pc.end() - 1);
			oneside_match = true;
			isReversed = false;
			addedBack = false;
			addedFront = true;
		}
		else if (*(chained.begin()) == *pc.begin())
		{ // oppsite direction
			chained.insert(chained.begin(), /*前方に反転して挿入*/ pc.rbegin(), pc.rend() - 1);
			oneside_match = true;
			isReversed = true;
			addedBack = false;
			addedFront = true;
		}
		checkComplete(oneside_match);
	};
};

///////////////////////////////////////
template <typename T>
std::vector<T *> Flatten(const std::vector<std::vector<T *>> &mat)
{
	std::vector<T *> ret;
	for (const auto &m : mat)
		for (const auto &n : m)
			ret.emplace_back(n);
	return ret;
};
template <typename T>
std::vector<T> Flatten(const std::vector<std::vector<T>> &mat)
{
	std::vector<T> ret;
	for (const auto &m : mat)
		for (const auto &n : m)
			ret.emplace_back(n);
	return ret;
};
template <typename T>
std::vector<T> Flatten(const std::vector<std::unordered_set<T>> &mat)
{
	std::vector<T> ret(0);
	ret.reserve(mat.size() * mat[0].size());
	for (const auto &m : mat)
		for (const auto &n : m)
			ret.emplace_back(n);
	return ret;
};
std::vector<Tddd> Flatten(const std::vector<std::vector<Tddd>> &mat)
{
	std::vector<Tddd> ret(0);
	ret.reserve(mat.size() * mat[0].size());
	for (const auto &m : mat)
		for (const auto &n : m)
			ret.emplace_back(n);
	return ret;
};
// 鎖状につなげていく．
// 1つ目のベクトルの順番を入れ替えることもあり得る．
// 2つ目のベクトルは，同じ方向を向いていなければならない．
template <typename T>
std::vector<T> FlattenAsChain(/*not const*/ std::vector<std::vector<T>> Vps)
{
	int count = 0; //無限ループをさける
	std::vector<T> ret = Vps[0];
	Vps.erase(Vps.begin());
	while (!Vps.empty())
	{
		for (auto it = std::begin(Vps); it < std::end(Vps); it++)
		{
			//前か後につぎつぎにくっつけていく
			if (*ret.rbegin() == *std::begin(*it))
			{
				ret.emplace_back(*std::rbegin(*it) /*next*/);
				Vps.erase(it);
				break;
			}
			else if (*ret.begin() == *std::rbegin(*it))
			{
				ret.insert(ret.begin(), *std::begin(*it) /*next*/);
				Vps.erase(it);
				break;
			}
			if (it == std::end(Vps) - 1)
				return ret;
		}

		if (count++ > 10000)
			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "loop infinity　多分，点の周りの面が重複したりして，数珠繋ぎにできない"));
	};
	return ret;
};
//////////////////////////////////////
template <class T>
std::vector<std::vector<T>> Transpose(const std::vector<std::vector<T>> &mat)
{
	if (mat.empty())
		return mat;
	std::vector<std::vector<T>> ans(mat[0].size(), std::vector<T>(mat.size()));
	for (size_t i = 0; i < mat.size(); i++)
		for (size_t j = 0; j < mat[i].size(); j++)
			ans[j][i] = mat[i][j];
	return ans;
};

VVV_d Transpose(const std::vector<VV_d> &mat)
{
	VVV_d ans(mat[0][0].size(), VV_d(mat[0].size(), V_d(mat.size())));
	for (size_t i = 0; i < mat.size(); i++)
		for (size_t j = 0; j < mat[i].size(); j++)
			for (size_t k = 0; k < mat[i][j].size(); k++)
				ans[k][j][i] = mat[i][j][k];
	return ans;
};

VV_d TensorProduct(const V_d &vec1, const V_d &vec2)
{
	VV_d ret(vec1.size(), V_d(vec2.size()));
	for (auto m = 0; m < vec1.size(); m++)
		for (auto j = 0; j < vec2.size(); j++)
			ret[m][j] = vec1[m] * vec2[j];
	return ret;
};

VVV_d TensorProductSet(const V_d &vec1, const V_d &vec2)
{
	VVV_d ret(vec1.size(), VV_d(vec2.size(), V_d(2, 0)));
	for (size_t m = 0; m < vec1.size(); m++)
		for (size_t j = 0; j < vec2.size(); j++)
			ret[m][j] = {vec1[m], vec2[j]};
	return ret;
};

double Dot3d(const std::vector<double> &vec1, const std::vector<double> &vec2) { return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]; };

template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Dot(const std::vector<T> &vec1, const std::vector<T> &vec2) { return std::inner_product(vec1.cbegin(), vec1.cend(), vec2.cbegin(), 0.); };
//// for pointer
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Dot(const std::vector<T *> &vec1, const std::vector<T *> &vec2)
{
	T ans(0.);
	for (size_t i = 0; i < vec1.size(); i++)
		ans += *vec1[i] * *vec2[i];
	return ans;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> Dot(const std::vector<std::vector<T>> &mat, const std::vector<T> &vec)
{
	std::vector<T> ans(mat.size());
	for (auto i = 0; i < mat.size(); i++)
		ans[i] = Dot(mat[i], vec);
	return ans;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> Dot(const std::vector<T> &vec, const std::vector<std::vector<T>> &mat)
{
	// return Dot(Transpose(mat), vec);
	std::vector<T> ans(mat[0].size(), 0.);
	for (auto j = 0; j < mat.size(); j++)
		for (auto i = 0; i < mat[j].size(); i++)
			ans[i] += mat[j][i] * vec[j];
	return ans;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<std::vector<T>> Dot(const std::vector<std::vector<T>> &mat1, const std::vector<std::vector<T>> &mat2)
{
	if (mat1[0].size() != mat2.size())
	{
		std::stringstream ss;
		ss << "passed variables have dimensions that can not be computed" << std::endl;
		ss << "mat1:" << mat1.size() << "x" << mat1[0].size() << std::endl;
		ss << "mat2:" << mat2.size() << "x" << mat2[0].size() << std::endl;
		error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
	}
	std::vector<std::vector<T>> ans(mat1.size(), std::vector<T>(mat2[0].size(), 0.));
	for (size_t x = 0; x < mat1.size(); x++)
		for (size_t y = 0; y < mat2[0].size(); y++)
			for (size_t j = 0; j < mat2.size(); j++)
				ans[x][y] += mat1[x][j] * mat2[j][y];
	return ans;
};
/* ------------------------------------------------------ */
template <typename T>
std::vector<T *> ToVector(const std::unordered_set<T *> &uo)
{
	std::vector<T *> ret(uo.size());
	int i = 0;
	for (const auto &p : uo)
		ret[i++] = p;
	return ret;
};
/* ------------------------------------------------------ */
template <typename T>
std::unordered_set<T *> ToUnorderedSet(const std::vector<T *> &v)
{
	return std::unordered_set<T *>(v.begin(), v.end());
};
/* ------------------------------------------------------ */
/*                     タプルのベクトル演算                  */
/* ------------------------------------------------------ */
std::vector<double> ToVector(const Tddd &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v)}; };
std::vector<double> ToVector(const T4d &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v)}; };
std::vector<double> ToVector(const T6d &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v),
													 std::get<3>(v), std::get<4>(v), std::get<5>(v)}; };
std::vector<double> ToVector(const T7d &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v),
													 std::get<3>(v), std::get<4>(v), std::get<5>(v),
													 std::get<6>(v)}; };
VV_d ToVector(const std::vector<Tddd> &v)
{
	VV_d ret(v.size(), {0, 0, 0});
	for (auto i = 0; i < v.size(); ++i)
		ret[i] = ToVector(v[i]);
	return ret;
};
VV_d ToVector(const T4T4d &v)
{
	return {ToVector(std::get<0>(v)),
			ToVector(std::get<1>(v)),
			ToVector(std::get<2>(v)),
			ToVector(std::get<3>(v))};
};
VV_d ToVector(const T6T6d &v)
{
	return {ToVector(std::get<0>(v)),
			ToVector(std::get<1>(v)),
			ToVector(std::get<2>(v)),
			ToVector(std::get<3>(v)),
			ToVector(std::get<4>(v)),
			ToVector(std::get<5>(v))};
};
VV_d ToVector(const T7T7d &v)
{
	return {ToVector(std::get<0>(v)),
			ToVector(std::get<1>(v)),
			ToVector(std::get<2>(v)),
			ToVector(std::get<3>(v)),
			ToVector(std::get<4>(v)),
			ToVector(std::get<5>(v)),
			ToVector(std::get<6>(v))};
};
T6d ToT6d(const Tddd tmp)
{
	return {std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), 0., 0., 0.};
};
Tddd ToTddd(const V_d &v) { return {v[0], v[1], v[2]}; };
Tddd ToTddd(const T6d &v) { return {std::get<0>(v), std::get<1>(v), std::get<2>(v)}; };
Tdd ToTdd(const V_d &v) { return {v[0], v[1]}; };
// std::vector<double> ToVector(const Tdd &v) { return {std::get<0>(v), std::get<1>(v)}; };
double Norm(const T4d &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) + std::get<1>(t) * std::get<1>(t) + std::get<2>(t) * std::get<2>(t) + std::get<3>(t) * std::get<3>(t)); };
double Norm(const T6d &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) +
											 std::get<1>(t) * std::get<1>(t) +
											 std::get<2>(t) * std::get<2>(t) +
											 std::get<3>(t) * std::get<3>(t) +
											 std::get<4>(t) * std::get<4>(t) +
											 std::get<5>(t) * std::get<5>(t)); };
double Norm(const T7d &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) +
											 std::get<1>(t) * std::get<1>(t) +
											 std::get<2>(t) * std::get<2>(t) +
											 std::get<3>(t) * std::get<3>(t) +
											 std::get<4>(t) * std::get<4>(t) +
											 std::get<5>(t) * std::get<5>(t) +
											 std::get<6>(t) * std::get<6>(t)); };
double Norm(const Tddd &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) + std::get<1>(t) * std::get<1>(t) + std::get<2>(t) * std::get<2>(t)); };
double Norm(const Tdd &t) { return std::sqrt(std::get<0>(t) * std::get<0>(t) + std::get<1>(t) * std::get<1>(t)); };
T4d Normalize(const T4d &X) { return X / Norm(X); };
T7d Normalize(const T7d &X) { return X / Norm(X); };
Tddd Normalize(const Tddd &X) { return X / Norm(X); };
Tdd Normalize(const Tdd &X) { return X / Norm(X); };
double Dot(const T6d &v, const T6d &u)
{
	return (std::get<0>(v) * std::get<0>(u) +
			std::get<1>(v) * std::get<1>(u) +
			std::get<2>(v) * std::get<2>(u) +
			std::get<3>(v) * std::get<3>(u) +
			std::get<4>(v) * std::get<4>(u) +
			std::get<5>(v) * std::get<5>(u));
};
Tddd Dot(const T6d &v, const T6Tddd &A)
{
	return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<0>(std::get<1>(A)) * std::get<1>(v) + std::get<0>(std::get<2>(A)) * std::get<2>(v) + std::get<0>(std::get<3>(A)) * std::get<3>(v) + std::get<0>(std::get<4>(A)) * std::get<4>(v) + std::get<0>(std::get<5>(A)) * std::get<5>(v),
			std::get<1>(std::get<0>(A)) * std::get<0>(v) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<1>(std::get<2>(A)) * std::get<2>(v) + std::get<1>(std::get<3>(A)) * std::get<3>(v) + std::get<1>(std::get<4>(A)) * std::get<4>(v) + std::get<1>(std::get<5>(A)) * std::get<5>(v),
			std::get<2>(std::get<0>(A)) * std::get<0>(v) + std::get<2>(std::get<1>(A)) * std::get<1>(v) + std::get<2>(std::get<2>(A)) * std::get<2>(v) + std::get<2>(std::get<3>(A)) * std::get<3>(v) + std::get<2>(std::get<4>(A)) * std::get<4>(v) + std::get<2>(std::get<5>(A)) * std::get<5>(v)};
};
double Dot(const T4d &v, const T4d &u)
{
	return (std::get<0>(v) * std::get<0>(u) +
			std::get<1>(v) * std::get<1>(u) +
			std::get<2>(v) * std::get<2>(u) +
			std::get<3>(v) * std::get<3>(u));
};
double Dot(const Tddd &v, const Tddd &u)
{
	return (std::get<0>(v) * std::get<0>(u) +
			std::get<1>(v) * std::get<1>(u) +
			std::get<2>(v) * std::get<2>(u));
};
double Dot(const Tdd &v, const Tdd &u)
{
	return (std::get<0>(v) * std::get<0>(u) + std::get<1>(v) * std::get<1>(u));
};
Tddd Dot(const T3Tddd &A, const Tddd &v)
{
	return {Dot(std::get<0>(A), v),
			Dot(std::get<1>(A), v),
			Dot(std::get<2>(A), v)};
};
Tddd Dot(const Tddd &v, const T3Tddd &A)
{
	return {std::get<0>(std::get<0>(A)) * std::get<0>(v) + std::get<0>(std::get<1>(A)) * std::get<1>(v) + std::get<0>(std::get<2>(A)) * std::get<2>(v),
			std::get<1>(std::get<0>(A)) * std::get<0>(v) + std::get<1>(std::get<1>(A)) * std::get<1>(v) + std::get<1>(std::get<2>(A)) * std::get<2>(v),
			std::get<2>(std::get<0>(A)) * std::get<0>(v) + std::get<2>(std::get<1>(A)) * std::get<1>(v) + std::get<2>(std::get<2>(A)) * std::get<2>(v)};
};

Tddd Cross(const Tddd &A, const Tddd &X)
{
	return {std::get<1>(A) * std::get<2>(X) - std::get<2>(A) * std::get<1>(X),
			std::get<2>(A) * std::get<0>(X) - std::get<0>(A) * std::get<2>(X),
			std::get<0>(A) * std::get<1>(X) - std::get<1>(A) * std::get<0>(X)};
};
T2Tdd Transpose(const T2Tdd &A)
{
	return {{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A))},
			{std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A))}};
};
T3Tdd Transpose(const T2Tddd &A)
{
	return {{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A))},
			{std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A))},
			{std::get<2>(std::get<0>(A)), std::get<2>(std::get<1>(A))}};
};
T3Tddd Transpose(const T3Tddd &A)
{
	return {{std::get<0>(std::get<0>(A)), std::get<0>(std::get<1>(A)), std::get<0>(std::get<2>(A))},
			{std::get<1>(std::get<0>(A)), std::get<1>(std::get<1>(A)), std::get<1>(std::get<2>(A))},
			{std::get<2>(std::get<0>(A)), std::get<2>(std::get<1>(A)), std::get<2>(std::get<2>(A))}};
};
/* ------------------------------------------------------ */

// bool isfinite(const V_d &v_IN)
// {
//   for (const auto &v : v_IN)
//     if (!std::isfinite(v))
//       return false;
//   return true;
// };

bool isfinite(const V_d &v_IN)
{
	for (const auto &v : v_IN)
		if (v < -1E+50 || v > 1E+50)
			return false;
	return true;
};

bool myIsfinite(const double v)
{
	if (v < -1E+50 || v > 1E+50)
		return false;
	else
		return true;
};

bool isFinite(const double v, const double eps = 1E+50)
{
	if (v < -eps || v > eps || v != v)
		return false;
	else
		return true;
};

bool isFinite(const V_d &v_IN)
{
	for (const auto &v : v_IN)
		if (!isFinite(v))
			return false;
	return true;
};

bool isFinite(const VV_d &vv_IN)
{
	for (const auto &v : vv_IN)
		if (!isFinite(v))
			return false;
	return true;
};

bool isFinite(const Tdd &v)
{
	if (isFinite(std::get<0>(v)) && isFinite(std::get<1>(v)))
		return true;
	return false;
};
bool isFinite(const T4d &v)
{
	if (isFinite(std::get<0>(v)) && isFinite(std::get<1>(v)) && isFinite(std::get<2>(v)) && isFinite(std::get<3>(v)))
		return true;
	return false;
};
bool isFinite(const Tddd &v)
{
	if (isFinite(std::get<0>(v)) && isFinite(std::get<1>(v)) && isFinite(std::get<2>(v)))
		return true;
	return false;
};
bool isFinite(const T3Tddd &V)
{
	auto [X, Y, Z] = V;
	if (isFinite(X) && isFinite(Y) && isFinite(Z))
		return true;
	return false;
};
/* ------------------------------------------------------ */
double Mean(const Tdd &X) { return (std::get<0>(X) + std::get<1>(X)) / 2.; };
double Mean(const Tddd &X) { return (std::get<0>(X) + std::get<1>(X) + std::get<2>(X)) / 3.; };
Tddd Mean(const T3Tddd &X)
{
	return {(std::get<0>(std::get<0>(X)) + std::get<0>(std::get<1>(X)) + std::get<0>(std::get<2>(X))) / 3.,
			(std::get<1>(std::get<0>(X)) + std::get<1>(std::get<1>(X)) + std::get<1>(std::get<2>(X))) / 3.,
			(std::get<2>(std::get<0>(X)) + std::get<2>(std::get<1>(X)) + std::get<2>(std::get<2>(X))) / 3.};
};
Tddd Mean(const std::vector<Tddd> &X)
{
	Tddd ret = {0., 0., 0.};
	for (const auto &x : X)
	{
		std::get<0>(ret) += std::get<0>(x);
		std::get<1>(ret) += std::get<1>(x);
		std::get<2>(ret) += std::get<2>(x);
	}
	return ret / ((double)X.size());
};
double Min(const Tdd &A) { return ((std::get<0>(A) < std::get<1>(A)) ? std::get<0>(A) : std::get<1>(A)); };
double Max(const Tdd &A) { return ((std::get<0>(A) > std::get<1>(A)) ? std::get<0>(A) : std::get<1>(A)); };
double Min(const Tddd &A)
{
	//イコールが重要
	auto [X, Y, Z] = A;
	if (X <= Y && X <= Z)
		return X;
	else if (Y <= X && Y <= Z)
		return Y;
	else
		return Z;
};
double FiniteMin(const Tddd &A)
{
	//イコールが重要
	auto [X, Y, Z] = A;
	if (isFinite(X) && !isFinite(Y) && !isFinite(Z))
	{
		Y = X;
		Z = X;
	}
	else if (!isFinite(X) && !isFinite(Y) && isFinite(Z))
	{
		Y = Z;
		X = Z;
	}
	else if (!isFinite(X) && isFinite(Y) && !isFinite(Z))
	{
		Z = Y;
		X = Y;
	}
	else if (isFinite(X) && isFinite(Y) && !isFinite(Z))
		Z = X;
	else if (isFinite(X) && !isFinite(Y) && isFinite(Z))
		Y = Z;
	else if (!isFinite(X) && isFinite(Y) && isFinite(Z))
		X = Z;
	//
	return Min(Tddd{X, Y, Z});
};
double Max(const Tddd &A)
{
	auto [X, Y, Z] = A;
	if (X >= Y && X >= Z)
		return X;
	else if (Y >= X && Y >= Z)
		return Y;
	else
		return Z;
};
double Max(const std::vector<Tddd> &A)
{
	double ret = -1E-10;
	for (const auto &a : A)
		if (ret < Max(a))
			ret = Max(a);
	return ret;
};
// double Max(const std::unordered_set<Tddd> &A)
// {
// 	double ret = -1E-10;
// 	for (const auto &a : A)
// 		if (ret < Max(a))
// 			ret = Max(a);
// 	return ret;
// };
/* ------------------------------------------------------ */
double Rot(const V_d vec1, const V_d vec2)
{
	return vec1[0] * vec2[1] - vec1[1] * vec2[0];
};
std::vector<V_d> Inv(const std::vector<V_d> &mat)
{
	std::vector<V_d> ans(mat.size(), V_d(mat[0].size(), 0.));
	double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
	ans[1][1] = mat[0][0] / det;
	ans[0][1] = -mat[0][1] / det;
	ans[1][0] = -mat[1][0] / det;
	ans[0][0] = mat[1][1] / det;
	return ans;
};
//==========================================================
// vector operators
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> Cross(const std::vector<T> &A)
{
	return {-A[1], A[0]};
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
std::vector<T> Cross(const std::vector<T> &A, const std::vector<T> &X)
{
	if (A.size() == 3)
		return {A[1] * X[2] - A[2] * X[1],
				A[2] * X[0] - A[0] * X[2],
				A[0] * X[1] - A[1] * X[0]};
	else if (A.size() == 2)
		return Cross(std::vector<T>{A[0], A[1], 0.}, std::vector<T>{X[0], X[1], 0.});
	else
	{
		std::stringstream ss;
		ss << "Invalid vectors are passed\n";
		ss << "A " << A << " X " << X;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
	}
};
//==========================================================
V_d log10(const V_d &vec)
{
	V_d ret(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
		ret[i] = std::log10(vec[i]);
	return ret;
};
V_d log(const V_d &vec)
{
	V_d ret(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
		ret[i] = std::log(vec[i]);
	return ret;
};
/* ------------------------------------------------------ */
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Norm(const T x) { return std::abs(x); };
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Norm(const std::vector<T> &vec) { return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.)); };
/* ------------------------------------------------------ */
V_d Normalize(const V_d X) { return X / Norm(X); };
double Norm3d(const V_d &vec)
{
	if (vec.size() != 3)
	{
		std::stringstream ss;
		ss << vec;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
	}
	return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T Abs(const V_d &vec)
{
	T tmp(0);
	for (size_t i = 0; i < vec.size(); i++)
		tmp += vec[i] * vec[i];
	return sqrt(tmp);
};
/* ------------------------------------------------------ */
// 2021/06/14
struct Quaternion
{
	double a, b, c, d;
	Tddd v;
	T4d q;

	// cos(q/2) + (ux*i + uy*j+ uz*k) * sin(q/2)
	Quaternion() : a(1), b(0), c(0), d(0), v({0, 0, 0}), q({1, 0, 0, 0}){};
	Quaternion(const double aIN, const double bIN, const double cIN, const double dIN) : a(aIN),
																						 b(bIN),
																						 c(cIN),
																						 d(dIN),
																						 v({b, c, d}),
																						 q({a, b, c, d}){};
	Quaternion(const T4d &qIN) : a(std::get<0>(qIN)),
								 b(std::get<1>(qIN)),
								 c(std::get<2>(qIN)),
								 d(std::get<3>(qIN)),
								 v({std::get<1>(qIN), std::get<2>(qIN), std::get<3>(qIN)}),
								 q(qIN){};
	Quaternion(const Tddd &axis, const double angle)
	{
		//空間回転　q = cos(theta/2) + n*sin(theta/2)
		this->v = Normalize(axis) * sin(-angle / 2.);
		this->a = cos(-angle / 2.);
		this->b = std::get<0>(v);
		this->c = std::get<1>(v);
		this->d = std::get<2>(v);
		this->q = {this->a, this->b, this->c, this->d};
	};
	/* ------------------------------------------------------ */
	T3Tddd Rv() const
	{
		//%固定した座標系(global座標)における．位置ベクトルの回転をするために使う．
		//ノーマライズされていなくていい
		double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
		return {{a2 + b2 - c2 - d2, 2 * b * c - 2 * a * d, 2 * a * c + 2 * b * d},
				{2 * b * c + 2 * a * d, a2 - b2 + c2 - d2, -2 * a * b + 2 * c * d},
				{-2 * a * c + 2 * b * d, 2 * a * b + 2 * c * d, a2 - b2 - c2 + d2}};
		// OK
	}
	T3Tddd Rs() const
	{
		//%物体座標系を回転することで，global座標が物体座標にとってどのように移動するかを計算するために使う．
		//ノーマライズされていなくていい
		double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
		return {{a2 + b2 - c2 - d2, 2 * b * c + 2 * a * d, -2 * a * c + 2 * b * d},
				{2 * b * c - 2 * a * d, a2 - b2 + c2 - d2, 2 * a * b + 2 * c * d},
				{2 * a * c + 2 * b * d, -2 * a * b + 2 * c * d, a2 - b2 - c2 + d2}};
		// OK
	}
	Tddd Rv(const Tddd &uIN) const { return Dot(this->Rv(), uIN); }
	Tddd Rs(const Tddd &uIN) const { return Dot(this->Rs(), uIN); }

	//ノーマライズされていなくていい
	T3Tddd R() const { return this->Rv(); }
	Tddd R(const Tddd &uIN) const { return Dot(this->Rv(), uIN); }
	/* ------------------------------------------------------ */
	//ドローン工学入門(1.74) or https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
	double yaw() const { return std::atan2(2 * (a * d + b * c), -1. + 2. * (c * c + d * d)); }
	double pitch() const { return std::asin(2 * (a * c - b * d)); }
	double roll() const { return std::atan2(2 * (a * b + c * d), -1. + 2. * (a * a + d * d)); }

	// yaw，pitch，rollは，回転の順序でもある
	Tddd YPR() const { return {this->yaw(), this->pitch(), this->roll()}; }

	T3Tddd Ryaw() const
	{
		double t = this->yaw();
		return T3Tddd{{cos(t), sin(t), 0.}, {-sin(t), cos(t), 0.}, {0., 0., 1.}};
	}

	T3Tddd Rpitch() const
	{
		double t = this->pitch();
		return T3Tddd{{cos(t), 0., -sin(t)}, {0., 1., 0.}, {sin(t), 0., cos(t)}};
	}

	T3Tddd Rroll() const
	{
		double t = this->roll();
		return T3Tddd{{1., 0., 0.}, {0., cos(t), sin(t)}, {0., -sin(t), cos(t)}};
	}
	/* ------------------------------------------------------ */

	const T4d &operator()() const { return q; }

	Quaternion conjugate() const
	{
		return Quaternion(T4d{a, -b, -c, -d});
	};

	Quaternion approxNextQuaternion(const Tddd &w, const double dt) const
	{
		return Quaternion(this->q + this->d_dt(w * dt)());
	};

	Quaternion d_dt(const Tddd &w /*angular velocity*/) const
	{
		//回転変化率w=dtheta/dtは，クォータニオン変化はQ.d_dt(w)と同じである．
		/*
		W = {{0, -w0, -w1, -w2},
			{w0, 0, w2, -w1},
			{w1, -w2, 0, w0},
			{w2, w1, -w0, 0}}

		inverse[W] = {{0, w0, w1, w2},
					{-w0, 0, -w2, w1},
					{-w1, w2, 0, -w0},
					{-w2, -w1, w0,0}}/(w0^2 + w1^2 + w2^2)
					= -W/(w0^2 + w1^2 + w2^2)
		*/

		return Quaternion((-this->b * std::get<0>(w) - this->c * std::get<1>(w) - this->d * std::get<2>(w)) / 2., /*{0,-w0,-w1,-w2}.{a,b,c,d}/2*/
						  (this->a * std::get<0>(w) - this->d * std::get<1>(w) + this->c * std::get<2>(w)) / 2.,  /*{w0,0,w2,-w1}.{a,b,c,d}/2*/
						  (this->d * std::get<0>(w) + this->a * std::get<1>(w) - this->b * std::get<2>(w)) / 2.,  /*{w1,-w2,0,w0}.{a,b,c,d}/2*/
						  (-this->c * std::get<0>(w) + this->b * std::get<1>(w) + this->a * std::get<2>(w)) / 2. /*{w2,w1,-w0,0}.{a,b,c,d}/2*/);
	};

	void set(const T4d &qIN)
	{
		double n = Norm(qIN);
		this->a = std::get<0>(qIN) / n;
		std::get<0>(this->v) = this->b = std::get<1>(qIN) / n;
		std::get<1>(this->v) = this->c = std::get<2>(qIN) / n;
		std::get<2>(this->v) = this->d = std::get<3>(qIN) / n;
		this->q = {this->a, this->b, this->c, this->d};
	};

	void set(const Quaternion &qIN)
	{
		this->set(qIN());
	};
};

Quaternion operator*(const Quaternion &A, const Quaternion &B)
{
	// https://en.wikipedia.org/wiki/Quaternion
	Tddd v = A.a * B.v + B.a * A.v + Cross(A.v, B.v);
	return Quaternion(T4d{(A.a) * (B.a) - Dot(A.v, B.v), std::get<0>(v), std::get<1>(v), std::get<2>(v)}); // ok
};
Quaternion &operator*=(Quaternion &A, const Quaternion &B)
{
	return A = (A * B); // ok
};
Quaternion operator*(Quaternion A, const double dt)
{
	A.set(A() * dt); // ok
	return A;
};
/* ------------------------- 足し算 ------------------------ */
Quaternion operator+(Quaternion B, const double A)
{
	//()はただのT4d
	B.set(A + B());
	return B;
};
Quaternion operator+(const double A, Quaternion B)
{
	//()はただのT4d
	B.set(A + B());
	return B;
};
Quaternion operator+(Quaternion A, const Quaternion &B)
{
	//()はただのT4d
	A.set(A() + B());
	return A;
};
Quaternion operator+(Quaternion A, const T4d &B)
{
	//()はただのT4d
	A.set(A() + B);
	return A;
};
Quaternion operator+(const T4d &B, Quaternion A)
{
	//()はただのT4d
	A.set(A() + B);
	return A;
};
Quaternion &operator+=(Quaternion &A, const Quaternion &B)
{
	A.set(A() + B());
	return A;
};
double Norm(const Quaternion &A)
{
	return Norm(A());
};
/* ------------------------------------------------------ */
double MyVectorAngle(const Tddd &v0, const Tddd &v1, const Tddd &frontdir /*右手系のzとなる*/)
{
	//      /|
	//     / |
	//    /  |y
	//   /   |
	//  /q___|
	//     x
	// q = atan2(double y, double x), q = [-pi, pi]
	//
	//
	// this can distingish ccw(positive) or cw(negative)

	// auto Y = Cross(frontdir, v0); //右手系
	// return atan2(Dot(v1, Y / Norm(Y)), Dot(v1, v0 / Norm(v0)));

	// V_d z = Cross(v0/*x*/, v1);
	// V_d z = frontdir;
	auto z = frontdir / Norm(frontdir);
	auto y = Cross(z, v0 /*x*/);
	y = y / Norm(y);
	auto x = v0 / Norm(v0);
	// {x,y,z}右手系

	//      /|
	//     / |
	//  v0/  |Dot(v1,y)
	//   /   |
	//  /q___|
	//     Dot(v1,v0/Norm(v0)

	return std::atan2(Dot(v1, y), Dot(v1, x));
};
double MyVectorAngle(const V_d &v0, const V_d &v1, const V_d &frontdir /*右手系のzとなる*/)
{
	//      /|
	//     / |
	//    /  |y
	//   /   |
	//  /q___|
	//     x
	// q = atan2(double y, double x), q = [-pi, pi]
	//
	//
	// this can distingish ccw(positive) or cw(negative)

	// auto Y = Cross(frontdir, v0); //右手系
	// return atan2(Dot(v1, Y / Norm(Y)), Dot(v1, v0 / Norm(v0)));

	// V_d z = Cross(v0/*x*/, v1);
	// V_d z = frontdir;
	V_d z = frontdir / Norm(frontdir);
	V_d y = Cross(z, v0 /*x*/);
	y = y / Norm(y);
	V_d x = v0 / Norm(v0);
	// {x,y,z}右手系

	//      /|
	//     / |
	//  v0/  |Dot(v1,y)
	//   /   |
	//  /q___|
	//     Dot(v1,v0/Norm(v0)

	return std::atan2(Dot(v1, y), Dot(v1, x));
};

double MyVectorAngle(const V_d &v0, const V_d &v1)
{
	// cannot distingish ccw(positive) or cw(negative)
	// return MyVectorAngle(v0, v1, Cross(v0, v1));

	return std::acos(Dot(v0, v1) / (Norm(v0) * Norm(v1)));
	//これはMathematicaの定義と同じ:
	// a = {a0, a1, a2}
	// b = {b0, b1, b2}
	// Refine[VectorAngle[a, b],
	//  Assumptions -> {a0 \[Element] Reals && a1 \[Element] Reals &&
	//     a2 \[Element] Reals && b0 \[Element] Reals &&
	//     b1 \[Element] Reals && b2 \[Element] Reals}]
};
double MyVectorAngle(const Tddd &v0, const Tddd &v1)
{
	return std::acos(Dot(v0, v1) / (Norm(v0) * Norm(v1)));
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T MathematicaVectorAngle(const std::vector<T> &V1, const std::vector<T> &V2)
{
	if (V1.size() > 1)
	{
		return std::atan2(Norm(Cross(V1, V2)), Dot(V1, V2));
	}
	else
	{
		throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, " size is invalid"));
		return 0.;
	}
};

// point base
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T VectorAngle(const std::vector<T> &X1, const std::vector<T> &X2, const std::vector<T> &X0)
{
	if (X1.size() > 1)
	{
		return std::atan2(Norm(Cross(X1 - X0, X2 - X0)), Dot(X1 - X0, X2 - X0));
	}
	else
	{
		throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, " size is invalid"));
		return 0.;
	}
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T VectorAngle(const std::vector<T> &X1, const std::vector<T> &X2)
{
	return VectorAngle(X1, X2, {0., 0.});
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T VectorAngleDirected(const std::vector<T> &X1, const std::vector<T> &X2)
{
	T a = VectorAngle(X1, X2, {0., 0.});
	return (a < 0) ? (2 * M_PI - a) : a;
};
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T DirectedArea(const std::vector<T> &vec1,
			   const std::vector<T> &vec2,
			   const std::vector<T> &n /*the direction*/)
{
	//  vec1
	// --->*
	//     | vec2
	//     v
	return Dot(Cross(vec1, vec2), n);
};

template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T DirectedArea(const std::vector<T> &a,
			   const std::vector<T> &b,
			   const std::vector<T> &c,
			   const std::vector<T> &n)
{
	// a --> b
	//       |
	//       v
	//	   c
	return Dot(Cross(b - a, c - b), n);
};

double TriangleArea(const V_d &a, const V_d &b, const V_d &c)
{
	// a --> b
	//       |
	//       v
	//	     c

	// old
	// std::vector<T> n = Cross(b - a, c - b);
	// return 0.5 * Dot(n, n / Norm(n));

	auto A = Norm(a - c);
	auto B = Norm(b - a);
	auto C = Norm(c - b);
	// this->area = std::abs(Norm(c) / Dot(a, b));
	auto s = 0.5 * (A + B + C);
	return std::sqrt(s * (s - A) * (s - B) * (s - C));
};
double TriangleArea(const T3Tddd &abc)
{
	auto [a, b, c] = abc;
	auto A = Norm(a - c);
	auto B = Norm(b - a);
	auto C = Norm(c - b);
	// this->area = std::abs(Norm(c) / Dot(a, b));
	auto s = 0.5 * (A + B + C);
	return std::sqrt(s * (s - A) * (s - B) * (s - C));
};
//
V_d TriangleAngles(const V_d &a, const V_d &b, const V_d &c)
{
	return {MyVectorAngle(b - a, c - a),
			MyVectorAngle(c - b, a - b),
			MyVectorAngle(a - c, b - c)};
};
Tddd TriangleAngles(const T3Tddd &abc)
{
	auto [a, b, c] = abc;
	return {MyVectorAngle(b - a, c - a),
			MyVectorAngle(c - b, a - b),
			MyVectorAngle(a - c, b - c)};
};
//
V_d TriangleNormal(const V_d &a, const V_d &b, const V_d &c)
{
	return Normalize(Cross((b - a), (c - a)));

	// if(isfinite(n))
	//   return n;
	// else
	// {
	//   std::stringstream ss;
	//   ss << "{a,b,c} = " << VV_d{a,b,c} << std::endl;
	//   ss << "n = " << n << std::endl;
	//   ss << "TriangleAngles = " << TriangleAngles(a,b,c) << std::endl;
	//   ss << "TriangleArea = " << TriangleArea(a,b,c) << std::endl;
	//   throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str()));
	// }
};

Tddd TriangleNormal(const Tddd &a, const Tddd &b, const Tddd &c)
{
	return Normalize(Cross((b - a), (c - a)));
};
Tddd TriangleNormal(const T3Tddd &abc)
{
	auto [a, b, c] = abc;
	return Normalize(Cross((b - a), (c - a)));
};

//多角形の頂点における隣り合う辺が作る面積を計算する．法線と逆の場合はマイナス．範囲は[pi,-pi]
V_d DirectedArea(const VV_d &abcd,
				 const V_d &n)
{
	V_d ret(abcd.size());
	ret[0] = DirectedArea(abcd[abcd.size() - 1], abcd[0], abcd[1], n);
	for (size_t i = 0; i < abcd.size() - 2; i++)
		ret[i + 1] = DirectedArea(abcd[i], abcd[i + 1], abcd[i + 2], n);
	ret[abcd.size() - 1] = DirectedArea(abcd[abcd.size() - 2], abcd[abcd.size() - 1], abcd[0], n);
	return ret;
};

double InteriorAngle(const V_d &x, const V_d &b, const V_d &z)
{
	// this can distingish ccw(positive) or cw(negative)
	auto Y = Cross(z, x); //右手系
	return std::atan2(Dot(b, Y / Norm(Y)), Dot(b, x / Norm(x)));
};

V_d RotateVector(const double angle, const V_d &axis, const V_d &v)
{
	V_d k = axis / Norm(axis);
	return v * std::cos(angle) + Cross(k, v) * std::sin(angle) + k * (Dot(k, v)) * (1. - std::cos(angle));
};

V_d RotateVectorPI2(const V_d &axis, const V_d &v)
{
	V_d k = axis / Norm(axis);
	return Cross(k, v) + k * (Dot(k, v));
};

double SphericalInteriorAngle(const V_d &org, const V_d &p0, const V_d &p1, const V_d &p2)
{
	auto v0 = (p0 - org);
	auto v1 = (p1 - org);
	auto v2 = (p2 - org);
	auto s01 = RotateVectorPI2(Cross(v1, v0), v1);
	auto s02 = RotateVectorPI2(Cross(v1, v2), v1);
	return MyVectorAngle(s01, s02);
};

V_d SphericalInteriorAngles(const V_d &org, const VV_d &ps)
{
	int s = ps.size();
	V_d ret(s);
	int i0, i1, i2;
	for (int i = 0; i < s; i++)
	{
		i0 = (s + i - 1) % s;
		i1 = (s + i) % s; // middle
		i2 = (s + i + 1) % s;
		ret[i] = SphericalInteriorAngle(org, ps[i0], ps[i1] /*middle*/, ps[i2]);
	}
	return ret;
};

double SphericalInteriorArea(const V_d &org, const V_d &p0, const V_d &p1, const V_d &p2)
{
	return 0.5 * Norm(p0 - p1) * Norm(p2 - p1) * sin(SphericalInteriorAngle(org, p0, p1, p2));
};

V_d SphericalInteriorArea(const V_d &org, const VV_d &ps)
{
	int s = ps.size();
	V_d ret(s);
	int i0, i1, i2;
	for (auto i = 0; i < s; i++)
	{
		i0 = (s + i - 1) % s;
		i1 = (s + i) % s; // middle
		i2 = (s + i + 1) % s;
		ret[i] = SphericalInteriorArea(org, ps[i0], ps[i1], ps[i2]);
	};
	return ret;
};

/*vector from center of sphere*/
double SphericalVectorAngle(const V_d v0, const V_d v1, const V_d v2)
{
	// ccw
	auto va = v0 / Norm(v0);
	auto s01 = RotateVectorPI2(Cross(va, v1 / Norm(v1)) /*axis*/, va);
	auto s02 = RotateVectorPI2(Cross(va, v2 / Norm(v2)) /*axis*/, va);
	return MyVectorAngle(s01, s02, Cross(s01, s02));
};

double MeanDistance(const VV_d &data)
{
	double ret(0.);
	int N(0);
	for (auto i = 0; i < data.size(); i++)
	{
		for (auto j = i + 1; j < data.size(); j++)
		{
			ret += Norm(data[i] - data[j]);
			N++;
		}
	}
	return ret / N;
};

double MeanDistance(const VV_d &data, const V_d &p)
{
	double ret(0.);
	for (auto i = 0; i < data.size(); i++)
		ret += Norm(data[i] - p);
	return ret / ((int)data.size());
};

#endif
