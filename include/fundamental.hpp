#ifndef fundamental_H
#define fundamental_H
#pragma once

// g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o fundamental_exception.hpp.gch fundamental_exception.hpp
// g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o fundamental_vectors.hpp.gch fundamental_vectors.hpp
// g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o basic_arithmetic_vector_operations.hpp.gch basic_arithmetic_vector_operations.hpp
// g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o fundamental_constants.hpp.gch fundamental_constants.hpp

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <variant>
#include <vector>
#include <functional>
#include <unordered_set>

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

using V_i = std::vector<int>;
using VV_i = std::vector<std::vector<int>>;
using VVV_i = std::vector<std::vector<std::vector<int>>>;

using V_s = std::vector<std::string>;
using VV_s = std::vector<std::vector<std::string>>;

using Tdd = std::tuple<double, double>;
using Tddd = std::tuple<double, double, double>;
using T4d = std::tuple<double, double, double, double>;
using T5d = std::tuple<double, double, double, double, double>;
using T6d = std::tuple<double, double, double, double, double, double>;
using T7d = std::tuple<double, double, double, double, double, double, double>;

double real_sph(const int l, const int m, const double theta, const double psi)
{
	return std::sph_legendre(std::abs(l), std::abs(l), theta) * std::cos(psi);
};
double real_sph_scale_ommited(const int l, const int m, const double theta, const double psi)
{
	return real_sph(l, m, theta, psi) / std::sqrt((2. * std::abs(l) + 1.) / (4. * M_PI));
};
/* ------------------------------------------------------ */
/*              表示方法に関するライブラリ                     */
/* ------------------------------------------------------ */
// 2021/08/02追加
using Var_SID = std::variant<std::string, int, size_t, double, std::vector<double>, std::vector<int>, Tddd, T6d>;
using V_Var_SID = std::vector<Var_SID>;
using VV_Var_SID = std::vector<std::vector<Var_SID>>;

std::string GridVector(const std::vector<int> &V, const int i = 20)
{
	std::stringstream ss;
	for (const auto &v : V)
	{
		ss << std::setw(i) << std::right;
		ss << v;
	}
	return ss.str();
};
std::string GridVector(const std::vector<double> &V, const int i = 20)
{
	std::stringstream ss;
	for (const auto &v : V)
	{
		ss << std::setw(i) << std::right;
		ss << v;
	}
	return ss.str();
};

std::string Grid(const V_Var_SID &V, const int i = 20)
{
	std::stringstream ss;
	for (const auto &var : V)
	{
		ss << std::setw(i) << std::right;
		if (std::holds_alternative<int>(var))
			ss << std::get<int>(var);
		else if (std::holds_alternative<size_t>(var))
			ss << std::get<size_t>(var);
		else if (std::holds_alternative<double>(var))
			ss << std::get<double>(var);
		else if (std::holds_alternative<std::string>(var))
			ss << std::get<std::string>(var);
		else if (std::holds_alternative<Tddd>(var))
		{
			auto v = std::get<Tddd>(var);
			std::stringstream s;
			s << "{"
			  << std::get<0>(v) << ","
			  << std::get<1>(v) << ","
			  << std::get<2>(v) << ","
			  << "}";
			ss << s.str();
		}
		else if (std::holds_alternative<T6d>(var))
		{
			auto v = std::get<T6d>(var);
			std::stringstream s;
			s << "{"
			  << std::get<0>(v) << ","
			  << std::get<1>(v) << ","
			  << std::get<2>(v) << ","
			  << std::get<3>(v) << ","
			  << std::get<4>(v) << ","
			  << std::get<5>(v) << ","
			  << "}";
			ss << s.str();
		}
		else if (std::holds_alternative<std::vector<double>>(var))
		{
			std::vector<double> v = std::get<std::vector<double>>(var);
			std::stringstream s;
			s << "{";
			if (!v.empty())
			{
				for (size_t i = 0; i < v.size() - 1; i++)
					s << v[i] << ",";
				s << *v.rbegin();
			}
			s << "}";
			ss << s.str();
		}
		else if (std::holds_alternative<std::vector<int>>(var))
		{
			std::vector<int> v = std::get<std::vector<int>>(var);
			std::stringstream s;
			s << "{";
			if (!v.empty())
			{
				for (size_t i = 0; i < v.size() - 1; i++)
					s << v[i] << ",";
				s << *v.rbegin();
			}
			s << "}";
			ss << s.str();
		}
		// else if (std::holds_alternative<std::vector<Tddd>>(var))
		// {
		// 	auto v = std::get<std::vector<Tddd>>(var);
		// 	std::stringstream s;
		// 	s << "{";
		// 	if (!v.empty())
		// 	{
		// 		for (size_t i = 0; i < v.size() - 1; i++)
		// 			s << (Tddd)v[i] << ",";
		// 		s << *v.rbegin();
		// 	}
		// 	s << "}";
		// 	ss << s.str();
		// }
	}
	return ss.str();
};

std::string GridGrid(const VV_Var_SID &VV, const int i = 10)
{
	std::stringstream ss;
	for (const auto &V : VV)
	{
		ss << Grid(V, i);
		ss << std::endl;
	}
	return ss.str();
};

/* ------------------------------------------------------ */
template <typename T, typename U>
std::map<T, U> KeyTake(std::map<T, U> &m, const std::vector<T> &V)
{
	std::map<T, U> ret;
	for (const auto &v : V)
		if (m.find(v) != m.end())
			ret[v] = m[v];
	return ret;
};
template <typename T, typename U>
std::vector<T> TakeFirst(const std::map<T, U> &m)
{
	std::vector<T> ret(m.size());
	int i = 0;
	for (auto it = m.begin(); it != m.end(); ++it)
		ret[i++] = it->first;
	return ret;
};
template <typename T, typename U>
std::vector<U> TakeSecond(const std::map<T, U> &m)
{
	std::vector<U> ret(m.size());
	int i = 0;
	for (auto it = m.begin(); it != m.end(); ++it)
		ret[i++] = it->second;
	return ret;
};
template <typename T, typename U>
std::vector<std::vector<U>> TakeSecond(const std::map<T, std::vector<U>> &m)
{
	std::vector<std::vector<U>> ret(m.size());
	int i = 0;
	for (auto it = m.begin(); it != m.end(); ++it)
		ret[i++] = it->second;
	return ret;
};

///////////////////////////////////////////////////////
// 20200526
V_i stoi(const V_s &vec)
{
	V_i ret(vec.size());
	std::transform(vec.begin(), vec.end(), ret.begin(), [](const std::string &str)
				   { return std::stoi(str); });
	return ret;
};
VV_i stoi(const std::vector<V_s> &vec)
{
	VV_i ret(0);
	for (const auto &v : vec)
		ret.emplace_back(stoi(v));
	return ret;
};
VVV_i stoi(const std::vector<std::vector<V_s>> &vec)
{
	VVV_i ret(0);
	for (const auto &v : vec)
		ret.emplace_back(stoi(v));
	return ret;
};
std::vector<bool> stob(const V_s &vec)
{
	std::vector<bool> ret(vec.size());
	std::transform(vec.begin(), vec.end(), ret.begin(), [](const std::string &str)
				   { 
					   bool b;
					   std::istringstream(str) >> std::boolalpha >> b;
					   return b; });
	return ret;
};
V_d stod(const V_s &vec)
{
	V_d ret(vec.size());
	std::transform(vec.begin(), vec.end(), ret.begin(), [](const std::string &str)
				   { return std::stod(str); });
	return ret;
};
VV_d stod(const std::vector<V_s> &vec)
{
	VV_d ret(0);
	for (const auto &v : vec)
		ret.emplace_back(stod(v));
	return ret;
};
VVV_d stod(const std::vector<std::vector<V_s>> &vec)
{
	VVV_d ret(0);
	for (const auto &v : vec)
		ret.emplace_back(stod(v));
	return ret;
};
///////////////////////////////////////////////////////
int InverseQuotientMod(const V_i &row_col, const int b)
{
	return row_col[0] * b + row_col[1];
};
V_i QuotientMod(const int a, const int b)
{
	return {(int)a / b, a % b};
};
//==========================================================
template <class T>
std::vector<T> Take(std::vector<T> &vec, const V_i &beg_end_skip)
{
	std::vector<T> ret;
	switch (beg_end_skip.size())
	{
	case 2:
		if (beg_end_skip[1] < 0)
		{
			for (auto i = beg_end_skip[0]; i < vec.size() + beg_end_skip[1]; i++)
				ret.emplace_back(vec[i]);
			return ret;
		}
		else
		{
			for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
				ret.emplace_back(vec[i]);
			return ret;
		}
	case 3:
		if (beg_end_skip[1] < 0)
		{
			for (auto i = beg_end_skip[0]; i < vec.size() + beg_end_skip[1]; i++)
				if (i % beg_end_skip[2] == 0)
					ret.emplace_back(vec[i]);
			return ret;
		}
		else
		{
			for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
				if (i % beg_end_skip[2] == 0)
					ret.emplace_back(vec[i]);
			return ret;
		}
	default:
		abort();
		return vec;
	};
};
template <class T>
std::vector<std::vector<T>> Take(std::vector<std::vector<T>> &mat, const V_i &beg_end_skip)
{
	std::vector<std::vector<T>> ret;
	switch (beg_end_skip.size())
	{
	case 2:
		if (beg_end_skip[1] < 0)
		{
			for (auto i = beg_end_skip[0]; i < mat.size() + beg_end_skip[1]; i++)
			{
				std::vector<T> vec;
				for (size_t j = 0; j < mat[i].size(); j++)
					vec.emplace_back(mat[i][j]);
				ret.emplace_back(vec);
			};
			return ret;
		}
		else
		{
			for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
			{
				std::vector<T> vec;
				for (size_t j = 0; j < mat[i].size(); j++)
					vec.emplace_back(mat[i][j]);
				ret.emplace_back(vec);
			};
			return ret;
		}
	case 3:
		if (beg_end_skip[1] < 0)
		{
			for (auto i = beg_end_skip[0]; i < mat.size().beg_end_skip[1]; i++)
			{
				if (i % beg_end_skip[2] == 0)
				{
					std::vector<T> vec;
					for (size_t j = 0; j < mat[i].size(); j++)
						vec.emplace_back(mat[i][j]);
					ret.emplace_back(vec);
				}
			};
			return ret;
		}
		else
		{
			for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
			{
				if (i % beg_end_skip[2] == 0)
				{
					std::vector<T> vec;
					for (size_t j = 0; j < mat[i].size(); j++)
						vec.emplace_back(mat[i][j]);
					ret.emplace_back(vec);
				}
			};
			return ret;
		}
	default:
		abort();
		return mat;
	};
};
template <class T>
std::vector<std::vector<T>> Take(std::vector<std::vector<T>> &mat, const V_i &beg_end_skip, const V_i &beg_end_skip2)
{
	std::vector<std::vector<T>> ret;
	switch (beg_end_skip.size())
	{
	case 2:
		if (beg_end_skip[1] < 0)
		{
			for (auto i = beg_end_skip[0]; i < mat.size() + beg_end_skip[1]; i++)
			{
				std::vector<T> vec;
				for (auto j = beg_end_skip[0]; j < beg_end_skip2[1]; j++)
					vec.emplace_back(mat[i][j]);
				ret.emplace_back(vec);
			};
			return ret;
		}
		else
		{
			for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
			{
				std::vector<T> vec;
				for (auto j = beg_end_skip[0]; j < beg_end_skip2[1]; j++)
					vec.emplace_back(mat[i][j]);
				ret.emplace_back(vec);
			};
			return ret;
		}
	case 3:
		if (beg_end_skip[1] < 0)
		{
			for (auto i = beg_end_skip[0]; i < mat.size() + beg_end_skip[1]; i++)
				if (i % beg_end_skip[2] == 0)
					ret.emplace_back(Take(mat[i], beg_end_skip2));
			return ret;
		}
		else
		{
			for (auto i = beg_end_skip[0]; i < beg_end_skip[1]; i++)
				if (i % beg_end_skip[2] == 0)
					ret.emplace_back(Take(mat[i], beg_end_skip2));
			return ret;
		}
	default:
		abort();
		return mat;
	};
};
//==========================================================
VV_d Import(const std::string &fname)
{
	std::ifstream fin(fname);
	VV_d mat;
	std::string str;
	double tmp;
	while (getline(fin, str))
	{
		std::stringstream ss;
		ss << str;
		V_d row;
		while (ss >> tmp)
			row.emplace_back(tmp);
		mat.emplace_back(row);
	}
	fin.close();
	return mat;
};
	// example
	// std::vector<std::vector<double>> mat=Import("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");

#include "fundamental_constants.hpp"
#include "fundamental_exception.hpp"
//==========================================================
// 20200716
template <typename T>
void Print(const T *const v_IN)
{
	std::cout << v_IN << std::endl;
};
template <typename T>
void Print(const T &v_IN)
{
	std::cout << v_IN << std::endl;
};
// void Print(const std::string &v_IN)
// {
// 	std::cout << v_IN << std::endl;
// };
// void Print(const std::vector<double> &v_IN)
// {
// 	std::cout << v_IN << std::endl;
// };
// void Print(const std::vector<int> &v_IN)
// {
// 	std::cout << v_IN << std::endl;
// };
template <typename T>
void Print(const T &v_IN, const std::string &color)
{
	std::cout << color << v_IN << reset << std::endl;
};
template <typename T>
void Print(const std::vector<T> &v_IN, const std::string &color)
{
	std::cout << color << v_IN << reset << std::endl;
};
///////////////////////////////////////////////////////////////
template <class T>
T Min(const std::vector<T> &v)
{
	return *std::min_element(std::begin(v), std::end(v));
};
template <class T>
T Min(const std::vector<std::vector<T>> &v)
{
	T tmp, ret(0.);
	for (size_t i = 0; i < v.size(); i++)
	{
		tmp = Min(v[i]);
		if (ret > tmp)
			ret = tmp;
	}
	return ret;
};

Tdd MinMax(const std::vector<double> &v)
{
	auto it = std::minmax_element(std::begin(v), std::end(v));
	return {*it.first, *it.second};
};

template <class T>
T Max(const std::vector<T> &v) { return *std::max_element(std::begin(v), std::end(v)); };
template <class T>
T Max(const std::vector<std::vector<T>> &v)
{
	T tmp, ret(0.);
	for (size_t i = 0; i < v.size(); i++)
	{
		tmp = Max(v[i]);
		if (ret < tmp)
			ret = tmp;
	}
	return ret;
};

template <class T, class U>
bool MemberQ(const std::map<T, U> &list, const T &form) { return list.find(form) != list.end(); };
template <class T, class U>
bool MemberQ(const std::map<T *, U> &list, const T *form) { return list.find(form) != list.end(); };
template <class T, class U>
bool AllMemberQ(const std::map<T *, U> &list, const std::vector<T *> &form)
{
	if (form.empty())
		throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "empty"));
	for (const auto &f : form)
		if (!MemberQ(list, f))
			return false;
	return true;
};
template <class T, class U>
bool AnyMemberQ(const std::map<T *, U> &list, const std::vector<T *> &form)
{
	if (form.empty())
		throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "empty"));
	for (const auto &f : form)
		if (MemberQ(list, f))
			return true;
	return false;
};
/* ------------------------------------------------------ */
template <class T>
int Count(const std::vector<T *> &list, const T *const v)
{
	return std::count(list.cbegin(), list.cend(), v);
};
/* ------------------------------------------------------ */
template <class T>
bool MemberQ(const std::tuple<T *, T *, T *, T *, T *, T *> &Ps, const T *const p)
{
	return (std::get<0>(Ps) == p ||
			std::get<1>(Ps) == p ||
			std::get<2>(Ps) == p ||
			std::get<3>(Ps) == p ||
			std::get<4>(Ps) == p ||
			std::get<5>(Ps) == p);
};

template <class T>
bool MemberQ(const std::vector<T> &list, const T &form) { return (std::find(list.cbegin(), list.cend(), form) != list.cend()); };
template <class T>
bool MemberQ(const std::vector<std::vector<T>> &list, const std::vector<T> &form) { return (std::find(list.cbegin(), list.cend(), form) != list.cend()); };
template <class T>
bool MemberQ(const std::vector<T *> &list, const T *const form) { return (std::find(list.cbegin(), list.cend(), form) != list.cend()); };
template <class T>
bool MemberQ(const std::vector<std::vector<T *>> &list, const std::vector<T *> &form) { return (std::find(list.cbegin(), list.cend(), form) != list.cend()); };

template <class T>
bool MemberQ(const std::tuple<T, T, T> &list, const T &form)
{
	if (std::get<0>(list) == form || std::get<1>(list) == form || std::get<2>(list) == form)
		return true;
	else
		return false;
};

template <class T>
bool MemberQ(const std::tuple<T *, T *, T *> &list, const T *const form)
{
	if (std::get<0>(list) == form || std::get<1>(list) == form || std::get<2>(list) == form)
		return true;
	else
		return false;
};

template <class T>
bool MemberQ(const std::tuple<T, T> &list, const T &form)
{
	if (std::get<0>(list) == form || std::get<1>(list) == form)
		return true;
	else
		return false;
};

template <class T>
bool MemberQ(const std::unordered_set<T> &list, const T &form)
{
	return list.contains(form);
};

// template <class T>
// bool MemberQ(const std::vector<T> &list, const T &form) { return std::binary_search(list.cbegin(), list.cend(), form); };
// template <class T>
// bool MemberQ(const std::vector<std::vector<T>> &list, const std::vector<T> &form) { return std::binary_search(list.cbegin(), list.cend(), form); };
// template <class T>
// bool MemberQ(const std::vector<T *> &list, const T *form) { return std::binary_search(list.cbegin(), list.cend(), form); };
// template <class T>
// bool MemberQ(const std::vector<std::vector<T *>> &list, const std::vector<T *> &form) { return std::binary_search(list.cbegin(), list.cend(), form); };

////////
template <typename T>
std::vector<T> TakeExcept(const std::vector<T> &list, const T &form)
{
	std::vector<T> ret;
	ret.reserve(list.size());
	for (const auto &l : list)
		if (l != form)
			ret.emplace_back(l);
	return ret;
};
template <typename T>
std::vector<T> TakeExcept(const std::vector<T> &list, const std::vector<T> &form)
{
	std::vector<T> ret;
	ret.reserve(list.size());
	for (const auto &l : list)
		if (!MemberQ(form, l))
			ret.emplace_back(l);
	return ret;
};
//
template <typename T>
std::unordered_set<T> TakeExcept(const std::unordered_set<T> &list, const T &form)
{
	std::unordered_set<T> ret;
	ret.reserve(list.size());
	for (const auto &l : list)
		if (l != form)
			ret.emplace_back(l);
	return ret;
};
template <typename T>
std::unordered_set<T> TakeExcept(const std::unordered_set<T> &list, const std::vector<T> &form)
{
	std::unordered_set<T> ret;
	ret.reserve(list.size());
	for (const auto &l : list)
		if (!MemberQ(form, l))
			ret.emplace(l);
	return ret;
};
template <typename T>
std::vector<T> TakeExceptIndex(const std::vector<T> &list, const int index)
{
	std::vector<T> ret;
	ret.reserve(list.size());
	// for (auto i = 0; (i < list.size() && i < index); i++)
	//   ret.emplace_back(list[i]);
	// for (auto i = index + 1; i < list.size(); i++)
	//   ret.emplace_back(list[i]);
	for (auto i = 0; i < list.size(); i++)
		if (i != index)
			ret.emplace_back(list[i]);
	return ret;
};
template <typename T>
std::vector<T> TakeExceptIndecies(const std::vector<T> &list, const V_i &form)
{
	std::vector<T> ret;
	ret.reserve(list.size());
	for (auto i = 0; i < list.size(); i++)
	{
		if (!MemberQ(form, i))
			ret.emplace_back(list[i]);
	}
	return ret;
};
template <class T>
bool AnyMemberQ(const std::vector<T> &list, const std::vector<T> &form)
{ // Mathematica like
	for (const auto &f : form)
		if (MemberQ(list, f))
			return true;
	return false;
};
template <class T>
bool AllMemberQ(const std::vector<T> &list, const std::vector<T> &form)
{ // Mathematica like
	for (const auto &f : form)
		if (!MemberQ(list, f))
			return false;
	return true;
};
// double Sum(const V_d::iterator first, const V_d::iterator second) { return std::reduce(first, second, 0.); };
// double Sum(const V_d &v) { return std::reduce(v.cbegin(), v.cend(), 0.); };
// double Total(const V_d &v) { return std::reduce(v.cbegin(), v.cend(), 0.); };

double Sum(const V_d::iterator first, const V_d::iterator second) { return std::accumulate(first, second, 0.); };
double Sum(const V_d &v) { return std::accumulate(v.cbegin(), v.cend(), 0.); };
double Total(const V_d &v) { return std::accumulate(v.cbegin(), v.cend(), 0.); };

#include "basic_arithmetic_vector_operations.hpp"

double Mean(const V_d &v)
{
	if (v.empty())
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "An empty vector is passed");
	else if (v.size() == 1)
		return v[0];
	else
		return std::accumulate(v.cbegin(), v.cend(), 0.) / v.size();
};

V_d Sum(const VV_d &v)
{
	if (v.empty())
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "An empty vector is passed");
	else if (v.size() == 1)
		return v[0];

	V_d ret(v[0].size(), 0.);
	for (const auto &u : v)
		ret += u;
	return ret;
};

V_d Total(const VV_d &v)
{
	if (v.empty())
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "An empty vector is passed");
	else if (v.size() == 1)
		return v[0];

	V_d ret(v[0].size(), 0.);
	for (const auto &u : v)
		ret += u;
	return ret;
};

double RMS(const V_d &v)
{
	double ret(0);
	for (const auto &u : v)
		ret += u * u;
	return std::sqrt(ret);
};
template <typename T>
std::vector<T> RandomSample(std::vector<T> ret /*copy may change*/, const int n = 0)
{
	std::shuffle(std::begin(ret), std::end(ret), std::default_random_engine());
	if (!n /*if n==0*/)
	{
		return ret;
	}
	else
	{
		return Take(ret, {0, n});
	}
};
template <typename T>
std::vector<T *> RandomSample(std::vector<T *> ret /*copy may change*/, const int n = 0)
{
	std::shuffle(std::begin(ret), std::end(ret), std::default_random_engine());
	if (!n /*if n==0*/)
	{
		return ret;
	}
	else
	{
		return Take(ret, {0, n});
	}
};

#include "fundamental_vectors.hpp"
//////////////////////////////////////////////////////////
template <typename T>
bool Between(const T x, std::tuple<T, T> minmax)
{
	if (!isFinite(x))
		return false;
	else if (std::get<0>(minmax) <= x && x <= std::get<1>(minmax))
		return true;
	else
		return false;
};
//////////////////////////////////////////////////////////
//==========================================================
struct Timer
{
	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	Timer()
	{
		start = std::chrono::high_resolution_clock::now();
	};
	std::chrono::duration<double> elapsed()
	{
		end = std::chrono::high_resolution_clock::now();
		return end - start;
	};
};
/* ------------------------------------------------------ */
class TimeWatch
{
public:
	// std::chrono::time_point<std::chrono::high_resolution_clock> start;
	std::chrono::time_point<std::chrono::system_clock> start;
	std::chrono::time_point<std::chrono::system_clock> lasttime;
	TimeWatch()
	{
		this->start = std::chrono::system_clock::now();
		this->lasttime = start;
	};
	//!一つ目は，前回からの経過時間，二つ目は作成した時間からの経過時間
	// std::vector<std::chrono::duration<double>> get()
	std::vector<double> get()
	{
		auto tmp = this->lasttime;
		this->lasttime = std::chrono::system_clock::now();
		// return {this->lasttime - tmp, this->lasttime - this->start};
		return {(std::chrono::duration<double>(this->lasttime - tmp)).count(),
				(std::chrono::duration<double>(this->lasttime - this->start)).count()};
	};
	std::vector<double> operator()()
	{
		auto tmp = this->lasttime;
		this->lasttime = std::chrono::system_clock::now();
		// return {this->lasttime - tmp, this->lasttime - this->start};
		return {(std::chrono::duration<double>(this->lasttime - tmp)).count(),
				(std::chrono::duration<double>(this->lasttime - this->start)).count()};
	};
};

//==========================================================
double sgn(const double x) { return (x > 0.) ? 1. : (x < 0. ? -1. : 0.); };
//==========================================================
// template <class T>
// double Mean(const std::vector<T> &v)
// {
//   if (v.size() == 1)
//     return v[0];
//   return std::accumulate(v.cbegin(), v.cend(), 0.) / v.size();
// };
// V_d Mean(const VV_d &v)
// {
//   if (v.size() == 1)
//   {
//     return v[0];
//   }
//   VV_d W = Transpose(v);
//   V_d ret(W.size());
//   for (size_t i = 0; i < ret.size(); i++)
//     ret[i] = Mean(W[i]);
//   ;
//   return ret;
// };
// double Variance(const V_d &v)
// {
//   int N = v.size();
//   double ret(0);
//   for (const auto &w : v)
//     ret = ret + w * w;
//   return ret / N - pow(Mean(v), 2.);
// };
//==========================================================
double SgLog(const double u, const double hsig, const double maxmin)
{
	return hsig + sgn(u) * exp(maxmin * (std::abs(u) - 1.));
};
double InvSgLog(const double h, const double hsig, const double maxmin)
{
	return sgn(h - hsig) * std::log(std::abs(h - hsig)) / maxmin + 1.;
};
V_d SgLog(const V_d &h, const double hsig, const double maxmin)
{
	V_d ret;
	std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [hsig, maxmin](double tmp)
				   { return SgLog(tmp, hsig, maxmin); });
	return ret;
};
double DSgLog(const double u, const double hsig, const double maxmin)
{
	return maxmin * sgn(u) * exp(maxmin * (std::abs(u) - 1.));
}; //単調増加
V_d DSgLog(const V_d &h, const double hsig, const double maxmin)
{
	V_d ret;
	std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [hsig, maxmin](double tmp)
				   { return DSgLog(tmp, hsig, maxmin); });
	return ret;
};
//============================================================
double Sg(const double h, const double h_sig, const double beta)
{
	return sgn(h) * pow(std::abs(h), beta) + h_sig;
};
double InvSg(const double h, const double h_sig, const double beta)
{
	return sgn(h - h_sig) * pow(std::abs(h - h_sig), 1. / beta);
};
V_d Sg(const V_d &h, const double h_sig, const double beta)
{
	V_d ret;
	std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [h_sig, beta](double tmp)
				   { return Sg(tmp, h_sig, beta); });
	return ret;
};
double DSg(const double h, const double h_sig, const double beta)
{
	return beta * pow(std::abs(h), beta - 1.);
}; //単調増加
V_d DSg(const V_d &h, const double h_sig, const double beta)
{
	V_d ret;
	std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [h_sig, beta](double tmp)
				   { return DSg(tmp, h_sig, beta); });
	return ret;
};
//==========================================================
double linspace(const V_d &v, int size, int i)
{
#if defined FULL_DEBUG
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
	return v[0] + (v[1] - v[0]) * i / (double)(size - 1.);
};
/* double init_position(int size, int i, V_d v, double c, double beta){ */
/*   cout << "max = " << v[1] << ", min = " << v[0] << ", d = " << v[1]-v[0] << endl; */
/*   return Sg( v[0] + (v[1] - v[0]) * i/(double)(size - 1.), c, beta);; */
/* }; */
double linspace(const V_d &v, int size, int i, double c, double beta)
{
#if defined FULL_DEBUG
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
	V_d w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
	return Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
};
void linspace(const V_d &v, int size, V_d &O_X)
{
#if defined FULL_DEBUG
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
	O_X.resize(size);
	for (int i = 0; i < size; i++)
		O_X[i] = v[0] + (v[1] - v[0]) * i / ((double)size - 1.);
};
void linspace(const V_d &v, int size, V_d &O_X, double c, double beta)
{
#if defined FULL_DEBUG
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
	V_d w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
	O_X.resize(size);
	for (int i = 0; i < size; i++)
		O_X[i] = Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
};
V_d linspace(const V_d &v, int size)
{
#if defined FULL_DEBUG
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
	V_d O_X(size);
	for (int i = 0; i < size; i++)
		O_X[i] = v[0] + (v[1] - v[0]) * i / (double)(size - 1.);
	return O_X;
};
V_d linspace(const V_d &v, int size, double c, double beta)
{
#if defined FULL_DEBUG
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
	V_d w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
	V_d O_X(size);
	for (int i = 0; i < size; i++)
		O_X[i] = Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
	return O_X;
};
void linspace(const std::vector<V_d> &mat, const V_i size, V_d &O_X)
{
#if defined FULL_DEBUG
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
	V_d v;
	for (size_t n = 0; n < mat.size(); n++)
	{
		linspace(mat[n], size[n], v);
		O_X.insert(std::end(O_X), std::begin(v), std::end(v));
	}
};
////////// Join has been modified on 20200411
template <typename T>
std::unordered_set<T> Join(std::unordered_set<T> ret, const std::unordered_set<T> &b)
{
	ret.insert(b.begin(), b.end());
	return ret;
};
template <typename T>
std::vector<T> Join(std::vector<T> ret, const std::vector<T> &b)
{
	ret.reserve(ret.size() + b.size());
	ret.insert(ret.end(), b.begin(), b.end());
	return ret;
};
template <typename T>
std::vector<std::vector<T>> Join(std::vector<std::vector<T>> ret, const std::vector<std::vector<T>> &b)
{
	ret.reserve(ret.size() + b.size());
	ret.insert(ret.end(), b.begin(), b.end());
	return ret;
};
template <typename T>
std::vector<T> Append(const std::vector<T> &a, const T &b)
{
	std::vector<T> ret = a;
	ret.insert(ret.end(), b);
	return ret;
};
/////////////////////////////////////////////////////
template <typename T>
std::vector<T> Range(const T xmin, const T xmax, const T di)
{
	std::vector<T> ret;
	T i = 0.;
	while (di * i + xmin <= xmax)
	{
		ret.push_back(di * i + xmin);
		i = i + 1.;
	}
	return ret;
};
// template <typename T>
// std::vector<T> Subdivide(const T xmin, const T xmax, const int n)
// {
// 	T dx = (xmax - xmin) / n;
// 	std::vector<T> ret(n + 1);
// 	for (int i = 0; i < n + 1; i++)
// 		ret[i] = i * dx + xmin;
// 	return ret;
// };
std::vector<double> Subdivide(const double xmin, const double xmax, const int n)
{
	if (n <= 0)
	{
		std::stringstream ss;
		ss << "Subdivide(" << xmin << "," << xmax << "," << n << ")の位置3に与えられた分割数は正の整数でなければなりません．";
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
	}
	double dx = (xmax - xmin) / n;
	std::vector<double> ret(n + 1);
	for (int i = 0; i < n + 1; i++)
		ret[i] = i * dx + xmin;
	return ret;
};
template <class T>
std::vector<std::vector<T>> Subdivide(const std::vector<T> &xmin, const std::vector<T> &xmax, const int n)
{
	std::vector<T> dx = (xmax - xmin) / n;
	std::vector<std::vector<T>> ret(xmin.size(), std::vector<T>(n + 1, 0.));
	for (size_t i = 0; i < xmin.size(); i++)
		ret[i] = Subdivide(xmin[i], xmax[i], n);
	return Transpose(ret);
};
template <typename T>
std::vector<T> SubdivideByStep(const T xmin, const T xmax, const T di)
{
	return Subdivide(xmin, xmax, (int)((xmax - xmin) / di));
};
template <typename T>
std::vector<std::vector<T>> SubdivideByStep(const std::vector<T> &xmin, const std::vector<T> &xmax, const T di)
{
	return Subdivide(xmin, xmax, (int)((Norm(xmax - xmin)) / di));
};
template <typename T>
std::vector<std::vector<T>> SubdivideByStep(const std::vector<std::vector<T>> &x, const T di)
{
	std::vector<std::vector<T>> ret;
	for (size_t i = 0; i < x.size() - 1; i++)
	{
		if (!ret.empty())
			ret.pop_back();
		ret = Join(ret, SubdivideByStep(x[i], x[i + 1], di));
	}
	return ret;
};
template <typename T>
std::vector<std::vector<T>> SubdivideByStepExclude(const std::vector<std::vector<T>> &x_in, const T di)
{
	std::vector<std::vector<T>> x = x_in;
	std::vector<std::vector<T>> ret;
	T n;
	for (size_t i = 0; i < x.size() - 1; i++)
	{
		if (x.size() == 1)
		{
			std::cout << "invalid length" << std::endl;
			std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
			return x;
		}
		else if (x.size() == 2)
		{
			return SubdivideByStep(x[i], x[i + 1], di);
		}
		else if (i == 0)
		{
			n = Norm(x[i + 1] - x[i]) / di;
			std::vector<std::vector<T>> tmp = SubdivideByStep(x[i], x[i + 1] - (x[i + 1] - x[i]) / (2. * n), di);
			ret = Join(ret, tmp);
		}
		else if (i == x.size() - 2)
		{
			n = Norm(x[i + 1] - x[i]) / di;
			std::vector<std::vector<T>> tmp = SubdivideByStep(x[i] + (x[i + 1] - x[i]) / (2. * n), x[i + 1], di);
			ret = Join(ret, tmp);
		}
		else if (i > 0 && i < x.size() - 2)
		{
			n = Norm(x[i + 1] - x[i]) / di;
			std::vector<std::vector<T>> tmp = SubdivideByStep(x[i] + (x[i + 1] - x[i]) / (2. * n), x[i + 1] - (x[i + 1] - x[i]) / (2 * n), di);
			ret = Join(ret, tmp);
		}
	}
	return ret;
};
template <typename T>
std::vector<std::vector<T>> Subdivide(const std::vector<std::vector<T>> &x, const int n)
{
	T len(0.);
	for (size_t i = 0; i < x.size() - 1; i++)
		len += Norm(x[i + 1] - x[i]);
	T di = len / n;
	int counter(0);
	while (counter < 100000)
	{
		std::vector<std::vector<T>> tmp = SubdivideByStep(x, di);
		if (tmp.size() == n)
			return SubdivideByStep(x, di);
		else if (tmp.size() > n)
		{
			di = di * 1.01;
		}
		else if (tmp.size() < n)
		{
			di = di * 0.99;
		};
		counter++;
	}
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
	std::cout << Red << "can not subdivide by " << n << " counter = " << counter << reset << std::endl;
	abort();
	return SubdivideByStep(x, di);
};
template <typename T>
std::vector<std::vector<T>> SubdivideExclude(const std::vector<std::vector<T>> &x, const int n)
{
	T len(0.);
	for (size_t i = 0; i < x.size() - 1; i++)
		len += Norm(x[i + 1] - x[i]);
	T di = len / n;
	int counter(0);
	while (counter < 100000)
	{
		std::vector<std::vector<T>> tmp = SubdivideByStepExclude(x, di);
		if (tmp.size() == n)
			return tmp;
		else if (tmp.size() > n)
		{
			di = di * 1.01;
		}
		else if (tmp.size() < n)
		{
			di = di * 0.999;
		};
		counter++;
	}
	std::cout << Red << __func__ << ":" << reset << std::endl;
	std::cout << Red << "can not subdivide by " << n << " counter = " << counter << reset << std::endl;
	abort();
	return SubdivideByStepExclude(x, di);
};
//==========================================================
std::vector<std::tuple<double, double>> GaussianQuadratureWeightsTuple(const int n,
																	   const double x1,
																	   const double x2)
{
	double EPS = 1E-14;
	std::tuple<double, double> zeros = {0., 0.};
	std::vector<std::tuple<double, double>> XW(n, zeros);
	double z1, z, xm, xl, pp, p3, p2, p1;
	int m = (n + 1) / 2;
	xm = 0.5 * (x2 + x1);
	xl = 0.5 * (x2 - x1);
	for (auto i = 0; i < m; i++)
	{
		z = cos(M_PI * (i + 0.75) / (n + 0.5));
		do
		{
			p1 = 1.0;
			p2 = 0.0;
			for (auto j = 0; j < n; j++)
			{
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
			}
			pp = n * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp;
		} while (std::abs(z - z1) > EPS);
		std::get<0>(XW[i]) = xm - xl * z;
		std::get<0>(XW[n - 1 - i]) = xm + xl * z;
		std::get<1>(XW[i]) = 2.0 * xl / ((1.0 - z * z) * pp * pp);
		std::get<1>(XW[n - 1 - i]) = std::get<1>(XW[i]);
	}
	return XW;
};
/* ------------------------------------------------------ */
VV_d GaussianQuadratureWeights(const int n,
							   const double x1,
							   const double x2)
{
	double EPS = 1E-14;
	VV_d XW(n, V_d(2, 0));
	double z1, z, xm, xl, pp, p3, p2, p1;
	int m = (n + 1) / 2;
	xm = 0.5 * (x2 + x1);
	xl = 0.5 * (x2 - x1);
	for (auto i = 0; i < m; i++)
	{
		z = cos(M_PI * (i + 0.75) / (n + 0.5));
		do
		{
			p1 = 1.0;
			p2 = 0.0;
			for (auto j = 0; j < n; j++)
			{
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
			}
			pp = n * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp;
		} while (std::abs(z - z1) > EPS);
		XW[i][0] = xm - xl * z;
		XW[n - 1 - i][0] = xm + xl * z;
		XW[i][1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
		XW[n - 1 - i][1] = XW[i][1];
	}
	return XW;
};

VV_d SingularGaussianQuadratureWeights(const int n,
									   const double x1,
									   const double x2,
									   const double xi_sing,
									   const double b /*magnitude of singularity*/)
{
	VV_d gw_T = Transpose(GaussianQuadratureWeights(n,
													InvSg(x1, xi_sing, b),
													InvSg(x2, xi_sing, b)));

	return Transpose(VV_d{Sg(gw_T[0], xi_sing, b),
						  gw_T[1] * DSg(gw_T[0], xi_sing, b)});
};
//==========================================================
void gauleg(const double x1, const double x2,
			V_d &x,
			V_d &w)
{
#if defined FULL_DEBUG
	std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;
#endif
	const double EPS = 1.0e-14;
	double z1, z, xm, xl, pp, p3, p2, p1;
	int n = x.size();
	int m = (n + 1) / 2;
	xm = 0.5 * (x2 + x1);
	xl = 0.5 * (x2 - x1);
	for (int i = 0; i < m; i++)
	{
		z = cos(M_PI * (i + 0.75) / (n + 0.5));
		do
		{
			p1 = 1.0;
			p2 = 0.0;
			for (int j = 0; j < n; j++)
			{
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
			}
			pp = n * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp;
		} while (std::abs(z - z1) > EPS);
		x[i] = xm - xl * z;
		x[n - 1 - i] = xm + xl * z;
		w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
		w[n - 1 - i] = w[i];
	}
};
//===========================================================
// iはNコのデータがある場合，0<=i<=N-1となる
// i=N-1の場合，q[N-1]=1.なのでh=1.はq[N-1]=1.<= h && h < q[N]

//============= 20191222 Interpolation =============
double Bspline(const double h, const V_d &q, const int i, const int K)
{
	switch (K)
	{
	case 1:
		return (q[i] <= h && h < q[i + 1]) ? 1. : 0.;
	case 2:
		return ((q[i] <= h && h < q[i + 1]) ? (h - q[i]) / (q[i + 1] - q[i]) : 0.) + ((q[i + 1] <= h && h < q[i + 2]) ? (q[i + 2] - h) / (q[i + 2] - q[i + 1]) : 0.);
	case 3:
		return ((q[i] <= h && h < q[i + 1]) ? (h - q[i]) / (q[i + 2] - q[i]) * (h - q[i]) / (q[i + 1] - q[i]) : 0.) + ((q[i + 1] <= h && h < q[i + 2]) ? ((h - q[i]) / (q[i + 2] - q[i]) * (q[i + 2] - h) + (q[i + 3] - h) / (q[i + 3] - q[i + 1]) * (h - q[i + 1])) / (q[i + 2] - q[i + 1]) : 0.) + ((q[i + 2] <= h && h < q[i + 3]) ? (q[i + 3] - h) / (q[i + 3] - q[i + 1]) * (q[i + 3] - h) / (q[i + 3] - q[i + 2]) : 0.);
	default:
		return (h - q[i]) / (q[i + K - 1] - q[i]) * Bspline(h, q, i, K - 1) + (q[i + K] - h) / (q[i + K] - q[i + 1]) * Bspline(h, q, i + 1, K - 1);
	};
};

V_d Bspline(const V_d &h, const V_d &q, const int i, const int K)
{
	// {B_0(x),B_1(x),B_2(x),B_3(x),.......}
	V_d ret(h.size());
	for (auto k = 0; k < h.size(); k++)
		ret[k] = Bspline(h[k], q, i, K);
	return ret;
};

V_d Bspline(const double h, const V_d &q, const int K)
{
	// {B_0(h),B_1(h),B_2(h),B_3(h),.......}
	V_d ret(q.size() - K);
	for (auto i = 0; i < q.size() - K; i++)
		ret[i] = Bspline(h, q, i, K);
	return ret;
};

VV_d Bspline(const V_d &H, const V_d &q, const int K)
{
	// {{B_0(x0),B_1(x0),B_2(x0),B_3(x0),.......},
	//  {B_0(x1),B_1(x1),B_2(x1),B_3(x1),.......},
	//  {B_0(x2),B_1(x2),B_2(x2),B_3(x2),.......}}
	VV_d ret(H.size());
	std::transform(H.begin(), H.end(), ret.begin(),
				   [q, K](const auto h)
				   { return Bspline(h, q, K); });
	return ret;
};

double D_Bspline(const double h, const V_d &q, const int i, const int K)
{
	return (K - 1.) * (1. / (q[i + K - 1] - q[i]) * Bspline(h, q, i, K - 1) - 1. / (q[i + K] - q[i + 1]) * Bspline(h, q, i + 1, K - 1));
};

double D_Bspline(const double h, const V_d &q, const int i, const int K, const int n)
{
	if (K == 1 && n > 0)
	{
		return 0.;
	}
	else
	{
		switch (n)
		{
		case 0:
			return Bspline(h, q, i, K);
		default:
			return (K - 1) * (1. / (q[i + K - 1] - q[i]) * D_Bspline(h, q, i, K - 1, n - 1) - 1. / (q[i + K] - q[i + 1]) * D_Bspline(h, q, i + 1, K - 1, n - 1));
		};
	};
};

V_d D_Bspline(const double h, const V_d &q, const int K)
{
	// {B_0(h),B_1(h),B_2(h),B_3(h),.......}
	V_d ret(q.size() - K);
	for (auto i = 0; i < q.size() - K; i++)
		ret[i] = D_Bspline(h, q, i, K);
	return ret;
};

/* ------------------------------------------------------ */
struct ludcmp_parallel
{
	int n;
	VV_d lu;
	V_i indx;
	double d;
	V_d vv;
	ludcmp_parallel(const VV_d &a) : n(a.size()), lu(a), indx(a.size()), vv(a.size(), 1.)
	{
		const double TINY = 1.0e-40;
		int i, imax, j;
		double big, temp;
		d = 1.0;
		for (auto i = 0; i < n; ++i)
		{
			big = 0.0;
			for (const auto &lu_i_j : lu[i])
				if ((temp = std::abs(lu_i_j)) > big)
					big = temp;
			if (big == 0.0)
				throw("Singular matrix in LUdcmp");
			vv[i] /= big;
		}
		/* ------------------------------------------------------ */
		int k = 0;
		for (auto &lu_k : lu)
		{
			big = 0.0;
			for (i = k; i < n; ++i)
				if ((temp = vv[i] * std::abs(lu[i][k])) > big)
				{
					big = temp;
					imax = i;
				}
			if (k != imax)
			{
				lu_k.swap(lu[imax]);
				d = -d;
				vv[imax] = vv[k];
			}
			indx[k] = imax;
			if (lu_k[k] == 0.0)
				lu_k[k] = TINY;

#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (auto i = k + 1; i < n; ++i)
			{
				lu[i][k] /= lu_k[k];
				for (auto j = k + 1; j < n; ++j)
					lu[i][j] -= lu[i][k] * lu_k[j];
			}

			// #ifdef _OPENMP
			// #pragma omp parallel for
			// #endif
			// 			for (auto i = k + 1; i < n; ++i)
			// 			{
			// 				auto temp = lu[i][k] /= lu_k[k];
			// 				for (auto j = k + 1; j < n; ++j)
			// 					lu[i][j] -= temp * lu_k[j];
			// 			}

			k++;
		}
	};

	void solve(const V_d &b, V_d &x)
	{
		int i, ii = 0, ip, j;
		double sum;
		if (b.size() != n || x.size() != n)
			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));

		// for (i=0;i<n;i++)
		//   x[i] = b[i];

		x = b;

		for (i = 0; i < n; i++)
		{
			ip = indx[i];
			sum = x[ip];
			x[ip] = x[i];
			if (ii != 0)
				for (j = ii - 1; j < i; j++)
					sum -= lu[i][j] * x[j];
			else if (sum != 0.0)
				ii = i + 1;
			x[i] = sum;
		}
		for (i = n - 1; i >= 0; i--)
		{
			sum = x[i];
			for (j = i + 1; j < n; j++)
				sum -= lu[i][j] * x[j];
			x[i] = sum / lu[i][i];
		}
	};

	void solve(const VV_d &b, VV_d &x)
	{
		int i, j, m = b[0].size();
		if (b.size() != n || x.size() != n || b[0].size() != x.size())
			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
		V_d xx(n);
		for (j = 0; j < m; j++)
		{
			for (i = 0; i < n; i++)
				xx[i] = b[i][j];
			solve(xx, xx);
			for (i = 0; i < n; i++)
				x[i][j] = xx[i];
		}
	};

	void inverse(VV_d &ainv)
	{
		int i, j;
		ainv.resize(n, V_d(n, 0));
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
				ainv[i][j] = 0.;
			ainv[i][i] = 1.;
		}
		solve(ainv, ainv);
	};

	VV_d Inverse()
	{
		VV_d ainv(n, V_d(n, 0));
		int i, j;
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
				ainv[i][j] = 0.;
			ainv[i][i] = 1.;
		}
		solve(ainv, ainv);
		return ainv;
	};

	double det()
	{
		double dd = d;
		for (int i = 0; i < n; i++)
			dd *= lu[i][i];
		return dd;
	};

	// void mprove(V_d &b, V_d &x)
	// {
	// 	int i, j;
	// 	V_d r(n);
	// 	this->aref.resize(n);
	// 	for (i = 0; i < n; i++)
	// 	{
	// 		double sdp = -b[i];
	// 		for (j = 0; j < n; j++)
	// 			sdp += (double)aref[i][j] * (double)x[j];
	// 		r[i] = sdp;
	// 	}
	// 	solve(r, r);
	// 	for (i = 0; i < n; i++)
	// 		x[i] -= r[i];
	// };
};

// #ifdef use_lapack

extern "C" void dgetrf_(const int *dim1,
						const int *dim2,
						double *a,
						const int *lda,
						int *ipiv,
						int *info);
extern "C" void dgetrs_(const char *TRANS,
						const int *N,
						const int *NRHS,
						const double *A,
						const int *LDA,
						const int *IPIV,
						double *B,
						const int *LDB,
						int *INFO);

struct lapack_lu
{
	const int dim;
	std::vector<double> a;
	const char trans = 'T';
	const int nrhs = 1, LDA, LDB;
	int info;
	std::vector<int> ipiv;
	~lapack_lu(){};
	lapack_lu(const std::vector<std::vector<double>> &aIN)
		: a(Flatten(aIN)),
		  dim(aIN.size()),
		  LDB(dim),
		  LDA(dim),
		  ipiv(dim)
	{
		dgetrf_(&dim, &dim, &*a.begin(), &LDA, &*ipiv.begin(), &info);
	};

	void solve(const std::vector<double> &rhd, std::vector<double> &b)
	{
		b = rhd;
		dgetrs_(&trans, &dim, &nrhs, &*a.begin(), &LDA, &*ipiv.begin(), &*b.begin(), &LDB, &info);
		if (info)
		{
			// std::cerr << e.what() << reset << std::endl;
			std::stringstream ss;
			ss << "LDB:" << LDB;
			ss << "\nLDA:" << LDA;
			ss << "\nipiv:" << ipiv;
			// ss << "\na:" << a;
			ss << "\nb:" << b;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
		};
	};
};

// #endif

struct ludcmp
{
	int n;
	VV_d lu;
	V_i indx;
	double d;
	V_d vv;
	ludcmp(const VV_d &a) : n(a.size()), lu(a), indx(a.size()), vv(a.size(), 1.)
	{
		const double TINY = 1.0e-40;
		int i, imax, j;
		double big, temp;
		// V_d vv(n);
		d = 1.0;
		// #ifdef _OPENMP
		// #pragma omp parallel for
		// #endif
		for (auto i = 0; i < n; i++)
		{
			big = 0.0;
			for (const auto &lu_i_j : lu[i])
				if ((temp = std::abs(lu_i_j)) > big)
					big = temp;
			if (big == 0.0)
				throw("Singular matrix in LUdcmp");
			vv[i] /= big;
		}
		/* ------------------------------------------------------ */
		// for (auto k = 0; k < n; k++)
		// {
		int k = 0;
		for (auto &lu_k : lu)
		{
			// auto &lu_k = lu[k];
			big = 0.0;
			for (i = k; i < n; ++i)
			{
				temp = vv[i] * std::abs(lu[i][k]);
				if (temp > big)
				{
					big = temp;
					imax = i;
				}
			}
			if (k != imax)
			{
				// for (auto j = 0; j < n; ++j)
				// 	std::swap(lu[imax][j], lu_k[j]);
				lu_k.swap(lu[imax]);
				d = -d;
				vv[imax] = vv[k];
			}
			indx[k] = imax;
			if (lu_k[k] == 0.0)
				lu_k[k] = TINY;

			for (auto i = k + 1; i < n; ++i)
			{
				temp = lu[i][k] /= lu_k[k];
				for (auto j = k + 1; j < n; ++j)
					lu[i][j] -= temp * lu_k[j];
			}
			k++;
		}
	};

	void solve(const V_d &b, V_d &x) const
	{
		int i, ii = 0, ip, j;
		double sum;
		if (b.size() != n || x.size() != n)
			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "サイズが違います"));

		x = b;

		for (i = 0; i < n; ++i)
		{
			ip = indx[i];
			sum = x[ip];
			x[ip] = x[i];
			if (ii != 0)
				for (j = ii - 1; j < i; j++)
					sum -= lu[i][j] * x[j];
			else if (sum != 0.0)
				ii = i + 1;
			x[i] = sum;
		}
		for (i = n - 1; i >= 0; i--)
		{
			sum = x[i];
			for (j = i + 1; j < n; j++)
				sum -= lu[i][j] * x[j];
			x[i] = sum / lu[i][i];
		}
	};

	void solve(const VV_d &b, VV_d &x)
	{
		int i, j, m = b[0].size();
		if (b.size() != n || x.size() != n || b[0].size() != x.size())
			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
		V_d xx(n);
		for (j = 0; j < m; j++)
		{
			for (i = 0; i < n; i++)
				xx[i] = b[i][j];
			solve(xx, xx);
			for (i = 0; i < n; i++)
				x[i][j] = xx[i];
		}
	};

	void inverse(VV_d &ainv)
	{
		int i, j;
		ainv.resize(n, V_d(n, 0));
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
				ainv[i][j] = 0.;
			ainv[i][i] = 1.;
		}
		solve(ainv, ainv);
	};

	VV_d Inverse()
	{
		VV_d ainv(n, V_d(n, 0));
		int i, j;
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
				ainv[i][j] = 0.;
			ainv[i][i] = 1.;
		}
		solve(ainv, ainv);
		return ainv;
	};

	double det()
	{
		double dd = d;
		for (int i = 0; i < n; i++)
			dd *= lu[i][i];
		return dd;
	};

	// void mprove(V_d &b, V_d &x)
	// {
	// 	int i, j;
	// 	V_d r(n);
	// 	this->aref.resize(n);
	// 	for (i = 0; i < n; i++)
	// 	{
	// 		double sdp = -b[i];
	// 		for (j = 0; j < n; j++)
	// 			sdp += (double)aref[i][j] * (double)x[j];
	// 		r[i] = sdp;
	// 	}
	// 	solve(r, r);
	// 	for (i = 0; i < n; i++)
	// 		x[i] -= r[i];
	// };
};

T3Tddd Inverse(const T3Tddd &mat)
{
	auto [x00, x01, x02] = std::get<0>(mat);
	auto [x10, x11, x12] = std::get<1>(mat);
	auto [x20, x21, x22] = std::get<2>(mat);
	double det = -x02 * x11 * x20 + x01 * x12 * x20 + x02 * x10 * x21 - x00 * x12 * x21 - x01 * x10 * x22 + x00 * x11 * x22;
	return {Tddd{-x12 * x21 + x11 * x22, x02 * x21 - x01 * x22, -x02 * x11 + x01 * x12} / det,
			Tddd{x12 * x20 - x10 * x22, -x02 * x20 + x00 * x22, x02 * x10 - x00 * x12} / det,
			Tddd{-x11 * x20 + x10 * x21, x01 * x20 - x00 * x21, -x01 * x10 + x00 * x11} / det};
};

VV_d Inverse(const VV_d &mat)
{
	ludcmp LU(mat);
	VV_d ret;
	LU.inverse(ret);
	return ret;
};

// V_d PeriodicKnot(const V_d& h, const int K){
//   int s(h.size());
//   V_d q(s + K);
//   double EPS = 1.E-15;
//   for(auto i=-K; i <= -1; i++)
//     q[i] = h[i+s-1]-(s-1);
//   for(auto i=s; i <= s+K-2; i++)
//     q[i] = h[i-s+1]+(s-1);
//   for(auto i=0; i < s; i++)
//     q[i + K] = h[i];
//   return q;
// };
// センサー値の微分を計算する方法として使えないか9月25日(土)
V_d UniformKnots(V_d h, const int K)
{
	double scale = *h.rbegin() - *h.begin();
	*h.begin() -= scale * 1E-14;
	*h.rbegin() += scale * 1E-14;
	return Subdivide(*h.begin(), *h.rbegin(), h.size() + K - 1);
};
V_d UniformKnots(const int size, const int K)
{
	return UniformKnots(Subdivide(double(-1), double(1), size - 1), K);
};
V_d OpenUniformKnots(V_d h, const int K)
{
	int s(h.size());
	V_d q(s + K);
	double EPS = 1E-15;
	*h.begin() -= EPS;
	*h.rbegin() += EPS;
	for (auto i = 0; i < K; i++)
	{
		q[i] = h[0] - EPS * (i - K + 1.);
		q[i + s] = h[s - 1] + EPS * i;
	}
	for (auto i = 0; i < s - K; i++)
		q[i + K] = (h[i] + h[(i + K)]) / 2.;

	return q;
};

V_d OpenUniformKnots(const int size, const int K)
{
	return OpenUniformKnots(Subdivide(double(-1), double(1), size - 1), K);
};

V_d PeriodicUniformKnots(const V_d &h, const int K)
{
	// 2021/09/25
	int s(h.size());
	V_d q(s + K);
	for (auto i = 0; i < s; i++)
		q[i] = h[i];
	for (auto i = -2 * K; i < 0; i++)
		q[i] = h[i + s - 1] - (s - 1);
	for (auto i = s; i < s + K; i++)
		q[i] = h[i - s + 1] + (s - 1);
	return q;
};

V_d Bspline_knot(const int size, const int K)
{
	V_d h = Subdivide(double(-1), double(1), size - 1), q(size + K);
	double EPS = 1.E-15;
	for (auto i = 0; i < K; i++)
		q[i] = h[0] - EPS * (i - K + 1.);
	for (auto i = 0; i < size - K; i++)
		q[i + K] = (h[i] + h[(i + K)]) / 2.;
	for (auto i = size; i < size + K; i++)
		q[i] = h[size - 1] + EPS * (i - size);
	// for (auto i = 0; i < K; i++)
	// {
	//   q[i] = h[0] - EPS * (i - K + 1.);
	//   q[i + size] = h[size - 1] + EPS * i;
	// }
	// for (auto i = 0; i < size - K; i++)
	//   q[i + K] = (h[i] + h[i + K]) / 2.;
	return q;
};

V_d Bspline_knot_periodic(const int size, const int K)
{
	V_d h = Subdivide(double(-1), double(1), size - 1), q(size + 4 * K - 2);
	for (auto i = 0; i < size; i++)
		q[i] = h[i];
	for (auto i = -2 * K; i < 0; i++)
		q[i] = h[i + size - 1] - (size - 1);
	for (auto i = size; i < size + 2 * K - 1; i++)
		q[i] = h[i - size + 1] + (size - 1);
	return q;
};

V_d Bspline_knot(const V_d &h, const int K)
{
	int s(h.size());
	V_d q(s + K);
	double EPS = 1E-15;
	for (auto i = 0; i < K; i++)
		q[i] = *h.begin() - EPS * (i - K + 1.);
	for (auto i = 0; i < s - K; i++)
		q[i + K] = (h[i] + h[(i + K)]) / 2.;
	for (auto i = s; i < s + K; i++)
		q[i] = *h.rbegin() + EPS * (i - s);
	return q;
};
VV_d Bspline_knot(const VV_d &h, const int K)
{
	VV_d q(h.size());
	for (size_t i = 0; i < h.size(); i++)
		q[i] = Bspline_knot(h[i], K);
	return q;
};
V_d Bspline_vector(const double h, const int s, const int K)
{
	//  --------------  s --------------->
	// {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)}
	V_d ret(s, 0.), q = OpenUniformKnots(s, K);
	for (auto j = 0; j < s; j++)
		ret[j] = Bspline(h, q, j, K);
	return ret;
}; // parametric -1 to 1
V_d Bspline_vector(const double h, const int s, const int K, const V_d &q)
{
	V_d ret(s, 0.);
	for (auto j = 0; j < s; j++)
		ret[j] = Bspline(h, q, j, K);
	return ret;
};
V_d D_Bspline_vector(const double h, const int s, const int K, int n)
{
	//  --------------  s --------------->
	// {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)}
	V_d ret(s, 0.), q = OpenUniformKnots(s, K);
	for (auto j = 0; j < s; j++)
		ret[j] = D_Bspline(h, q, j, K, n);
	return ret;
};
V_d D_Bspline_vector(const double h, const int s, const int K, int n, const V_d &q)
{
	V_d ret(s, 0.);
	for (auto j = 0; j < s; j++)
		ret[j] = D_Bspline(h, q, j, K, n);
	return ret;
};
//           --------------  s  --------------->
// |        {
// |        {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)},
// |        {B0(h1),B1(h1),B2(h1),B3(h1),B4(h1)},
// |        {B0(h2),B1(h2),B2(h2),B3(h2),B4(h2)},
// h.size() {B0(h3),B1(h3),B2(h3),B3(h3),B4(h3)},
// |        {B0(h4),B1(h4),B2(h4),B3(h4),B4(h4)},
// |        {B0(h5),B1(h5),B2(h5),B3(h5),B4(h5)},
// |        {B0(h6),B1(h6),B2(h6),B3(h6),B4(h6)},
// V        }
VV_d Bspline_matrix(const V_d &h, const int s, const int K)
{
	VV_d ret(h.size(), V_d(s, 0.));
	for (size_t i = 0; i < h.size(); i++)
		ret[i] = Bspline_vector(h[i], s, K);
	return ret;
};
VV_d Bspline_matrix(const V_d &h, const int s, const int K, const V_d &q)
{
	VV_d ret(h.size(), V_d(s, 0.));
	for (size_t i = 0; i < h.size(); i++)
		ret[i] = Bspline_vector(h[i], s, K, q);
	return ret;
};
VV_d Bspline_matrix(const VV_d &h, const int s, const int K, const VV_d &q)
{
	VV_d ret(h.size(), V_d(s, 0.));
	for (size_t i = 0; i < h.size(); i++)
		for (auto j = 0; j < s; j++)
			ret[i][j] = Bspline(h[i][j], q[i], s, K);
	return ret;
};
VV_d Bspline_matrix(const VV_d &h, const int s, const int K, const V_d &q)
{
	VV_d ret(h.size(), V_d(s, 0.));
	for (size_t i = 0; i < h.size(); i++)
		for (auto j = 0; j < s; j++)
			ret[i][j] = Bspline(h[i][j], q, i, K);
	return ret;
};
VV_d D_Bspline_matrix(const V_d &h, const int s, const int K, const int n)
{
	VV_d ret(h.size(), V_d(s, 0.));
	for (size_t i = 0; i < h.size(); i++)
		ret[i] = D_Bspline_vector(h[i], s, K, n);
	return ret;
};
VV_d D_Bspline_matrix(const V_d &h, const int s, const int K, const int n, const V_d &q)
{
	VV_d ret(h.size(), V_d(s, 0.));
	for (size_t i = 0; i < h.size(); i++)
		ret[i] = D_Bspline_vector(h[i], s, K, n, q);
	return ret;
};
//    ---------------  s --------------->
// | {
// | {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)},
// | {B0(h1),B1(h1),B2(h1),B3(h1),B4(h1)},
// s {B0(h2),B1(h2),B2(h2),B3(h2),B4(h2)},
// | {B0(h3),B1(h3),B2(h3),B3(h3),B4(h3)},
// | {B0(h4),B1(h4),B2(h4),B3(h4),B4(h4)}
// V }
VV_d Bspline_matrix(const int s, const int K)
{
	return Bspline_matrix(Subdivide(-1., 1., s - 1) /*h*/, s, K);
};
VV_d D_Bspline_matrix(const int s, const int K, const int n)
{
	return D_Bspline_matrix(Subdivide(-1., 1., s - 1) /*h*/, s, K, n);
};
///////////////////////////////////////////////////////////
struct ParametricInterpolation
{
public:
	VV_d samp3;
	V_d samp2;
	V_d q, xi, eta;
	V_d q1, q2;
	VV_d InvMatB, InvMatBT;
	int s, K;
	int s1, s2;
	///////////// 1 parm  ///////////
	ParametricInterpolation(const V_d &samp2_IN, const int K_IN) : samp2(samp2_IN),
																   s(samp2_IN.size()),
																   xi(Subdivide(double(-1), double(1), samp2_IN.size() - 1)),
																   q(OpenUniformKnots(Subdivide(double(-1), double(1), samp2_IN.size() - 1), K_IN)),
																   K(K_IN),
																   InvMatBT(Inverse(Transpose(Bspline_matrix(samp2_IN.size(), K_IN)))),
																   InvMatB(Inverse(Bspline_matrix(samp2_IN.size(), K_IN))){};
	V_d N(const double a) { return Dot(InvMatBT, Bspline_vector(a, s, K)); };
	double operator()(const double a) { return Dot(N(a), samp2); };
	//===========================
	V_d N(const double a, const int n) { return Dot(InvMatBT, D_Bspline_vector(a, s, K, n)); };
	double operator()(const double a, const int n) { return Dot(Dot(InvMatB, D_Bspline_vector(a, s, K, n)), samp2); };
	//========= specific ========
	double Basis(const double h, const int i) { return Bspline(h, q, i, K); };
	V_d Basis(const V_d &h, const int i) { return Bspline(h, q, i, K); };
	double Shape(const double h, const int i) { return Dot(InvMatBT[i], Bspline_vector(h, s, K)); };
	double DShape(const double h, const int i, const int n) { return Dot(InvMatBT[i], D_Bspline_vector(h, s, K, n)); };
	/////////// 2 parm  ////////////
	/*
	s1 is the size of a column vector : points in y direction
	s2 is the size of a row vector : points in x direction
  */
	ParametricInterpolation(const VV_d &samp3_IN, const int K_IN) : samp3(samp3_IN),
																	s1(samp3_IN.size()) /*y*/,
																	s2(samp3_IN[0].size()) /*x*/,
																	q1(OpenUniformKnots(samp3_IN.size(), K_IN)),
																	q2(OpenUniformKnots(samp3_IN[0].size(), K_IN)),
																	K(K_IN),
																	InvMatB(Inverse(Bspline_matrix(samp3_IN.size() /*y*/, K_IN))),
																	InvMatBT(Inverse(Transpose(Bspline_matrix(samp3_IN[0].size() /*x*/, K_IN)))){};
	ParametricInterpolation(){};
	void reset(const VV_d &samp3_IN, const int K_IN)
	{
		q1.clear();
		q2.clear();
		samp3.clear();
		InvMatB.clear();
		InvMatBT.clear();
		// insert
		samp3 = samp3_IN;
		s1 = samp3_IN.size() /*y*/;
		s2 = samp3_IN[0].size() /*x*/;
		q1 = OpenUniformKnots(samp3_IN.size(), K_IN);
		q2 = OpenUniformKnots(samp3_IN[0].size(), K_IN);
		K = K_IN;
		InvMatB = Inverse(Bspline_matrix(samp3_IN.size(), K_IN));
		InvMatBT = Inverse(Transpose(Bspline_matrix(samp3_IN[0].size(), K_IN)));
	};
	////////////////////////////////////////////////////
	V_d Nx(const double b /*y*/)
	{
		return Dot(InvMatBT, Bspline_vector(b, s2, K, q2));
	};
	V_d Ny(const double a /*x*/)
	{
		return Dot(Bspline_vector(a, s1, K, q1), InvMatB);
	};
	VV_d N(const V_d &ba)
	{
		return TensorProduct(Ny(ba[1]), Nx(ba[0]));
	};
	double operator()(const double b /*x*/, const double a /*y*/)
	{
		return Dot(Ny(a), Dot(samp3, Nx(b)));
	};
	double operator()(const V_d &ba)
	{
		return Dot(Ny(ba[1]), Dot(samp3, Nx(ba[0])));
	};
	////////////////// DERIVATIVES /////////////////
	V_d DNx(const double b /*y*/, const int n)
	{
		return Dot(InvMatBT, D_Bspline_vector(b, s2, K, n, q2));
	};
	V_d DNy(const double a /*x*/, const int n)
	{
		return Dot(D_Bspline_vector(a, s1, K, n, q1), InvMatB);
	};
	VV_d DN(const V_d &ba, const V_i &n)
	{
		return TensorProduct(DNy(ba[1], n[1]), DNx(ba[0], n[0]));
	};
	double operator()(const double b /*x*/, const double a /*y*/, const V_i &n)
	{
		return Dot(DNy(a, n[1]), Dot(samp3, DNx(b, n[0])));
	};
	double operator()(const V_d &ba, const V_i &n)
	{
		return Dot(DNy(ba[1], n[1]), Dot(samp3, DNx(ba[0], n[0])));
	};
	///////////////// specific shape function ////////////
	double Nx(const double b /*y*/, const int j)
	{
		return Dot(InvMatBT[j], Bspline_vector(b, s2, K, q2));
	};
	double Ny(const double a /*x*/, const int m)
	{
		V_d By = Bspline_vector(a, s1, K, q1);
		double ret(0.);
		for (auto l = 0; l < s1; l++)
			ret += InvMatB[l][m] * By[l];
		return ret;
	};
	double N(const V_d &ba, const V_i &jm)
	{ /*2D vector access*/
		return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
	};
	double N(const V_d &ba, const int mj)
	{ /*1D access*/
		return Ny(ba[1], (int)(mj / s2) /*row*/) * Nx(ba[0], mj % s2 /*col*/);
	};
	// double Basis(const V_d& ba, const V_i &jm){
	//   return samp3[jm[1]][jm[0]] * Nx(ba[0], jm[0]) * Ny(ba[1], jm[1]);
	// };
	double Basis(const V_d &ba, const V_i &jm)
	{
		return Bspline(ba[0], q2, jm[0], K) * Bspline(ba[1], q1, jm[1], K);
	};
	double Shape(const V_d &ba, const V_i &jm)
	{
		return Nx(ba[0], jm[0]) * Ny(ba[1], jm[1]);
	};
};
///////////////////////////////////////////////////////////
class node_complex
{
	//    *
	//    |
	// *--*--*
	//    |
	//    *
public:
	VV_d v_complex;
	std::vector<node_complex *> nodes;
	std::vector<bool> hits;
	node_complex(const VV_d &v_complex_IN) : v_complex(v_complex_IN){};

	void connect(node_complex &node_IN)
	{
		node_complex *p = std::addressof(node_IN);
		for (auto &n : nodes)
			if (p == n)
				return;
		nodes.emplace_back(p);
	};
};

///////////////////////////////////////////////////////////
class glLINES
{
public:
	VV_d v_complex;
	std::vector<std::vector<float>> s_complex;
	VV_i f_v_complex; //これのstd::vector<int>が１面 -> 3頂点を指定
	std::vector<bool> hits;

	VV_d samp3X, samp3Y, samp3Z;
	// sample点がない場合
	glLINES(int row = 20, int col = 20, double scale = 5.)
	{
		/* sampleを作成 */
		VV_d samp3X(row, V_d(col, 0));
		VV_d samp3Y(row, V_d(col, 0));
		VV_d samp3Z(row, V_d(col, 0));
		for (auto i = 0; i < row; i++)
		{
			for (auto j = 0; j < col; j++)
			{
				double x = (-1. + 2. / (col - 1) * j);
				double y = (-1. + 2. / (row - 1) * i);
				samp3X[i][j] = x * scale;
				samp3Y[i][j] = y * scale;
				samp3Z[i][j] = sin(2. * M_PI * x) * sin(2. * M_PI * y);
			}
		}
		/* gl.LINESでメッシュがかけるような形式で格納 */
		set(samp3X, samp3Y, samp3Z);
	};
	// sample点がある場合
	glLINES(const VV_d &samp3X_IN,
			const VV_d &samp3Y_IN,
			const VV_d &samp3Z_IN) : samp3X(samp3X_IN), samp3Y(samp3Y_IN), samp3Z(samp3Z_IN)
	{
		/* gl.LINESでメッシュがかけるような形式で格納 */
		set(samp3X_IN, samp3Y_IN, samp3Z_IN);
	};
	// verticesの格納
	void set(const VV_d &samp3X_IN,
			 const VV_d &samp3Y_IN,
			 const VV_d &samp3Z_IN)
	{
		auto col = samp3X_IN[0].size();
		auto row = samp3X_IN.size();
		v_complex.clear();
		s_complex.clear();
		f_v_complex.clear();
		/* gl.LINESでメッシュがかけるような形式で格納 */
		int l = 0;
		//    8 points
		//    0   1
		//   7*---*2
		//    |   |
		//    |   |
		//   6*---*3
		//    5   4
		for (size_t i = 0; i < row - 1; i++)
			for (size_t j = 0; j < col - 1; j++)
			{
				v_complex.push_back({samp3X_IN[i][j], samp3Z_IN[i][j], samp3Y_IN[i][j]});			  // 0
				v_complex.push_back({samp3X_IN[i][j + 1], samp3Z_IN[i][j + 1], samp3Y_IN[i][j + 1]}); // 1
				f_v_complex.push_back({l++, l++});
				v_complex.push_back({samp3X_IN[i][j + 1], samp3Z_IN[i][j + 1], samp3Y_IN[i][j + 1]});			  // 2
				v_complex.push_back({samp3X_IN[i + 1][j + 1], samp3Z_IN[i + 1][j + 1], samp3Y_IN[i + 1][j + 1]}); // 3
				f_v_complex.push_back({l++, l++});
				v_complex.push_back({samp3X_IN[i + 1][j + 1], samp3Z_IN[i + 1][j + 1], samp3Y_IN[i + 1][j + 1]}); // 4
				v_complex.push_back({samp3X_IN[i + 1][j], samp3Z_IN[i + 1][j], samp3Y_IN[i + 1][j]});			  // 5
				f_v_complex.push_back({l++, l++});
				v_complex.push_back({samp3X_IN[i + 1][j], samp3Z_IN[i + 1][j], samp3Y_IN[i + 1][j]}); // 6
				v_complex.push_back({samp3X_IN[i][j], samp3Z_IN[i][j], samp3Y_IN[i][j]});			  // 7
				f_v_complex.push_back({l++, l++});
			}

		s_complex.resize(v_complex.size(), std::vector<float>(4, 0)); // RGBS
		for (auto &tmp : s_complex)
			tmp = {.9, .0, .0, 1.};

		this->hits.resize(f_v_complex.size());
	};
};
///////////////////////////////////////////////////////////
struct ParametricInterpolation3D
{
public:
	VV_d v_complex;
	std::vector<std::vector<float>> s_complex;
	VV_i f_v_complex; //これのstd::vector<int>が１面 -> 3頂点を指定

	VV_d samp3X, samp3Y, samp3Z;
	V_d q1, q2;
	VV_d inv_coefmat, inv_coefmatT;
	int s1, s2, K;
	////////// 2 parm  //////////
	/*
	s1 is the size of a column vector : points in y direction
	s2 is the size of a row vector : points in x direction
  */
	ParametricInterpolation3D(const VV_d &samp3X_IN,
							  const VV_d &samp3Y_IN,
							  const VV_d &samp3Z_IN,
							  const int K_IN) : K(K_IN),
												s1(samp3X_IN.size()) /*y*/,
												s2(samp3X_IN[0].size()) /*x*/,
												q1(OpenUniformKnots(samp3X_IN.size(), K_IN)),
												q2(OpenUniformKnots(samp3X_IN[0].size(), K_IN)),
												inv_coefmat(Inverse(Bspline_matrix(samp3X_IN.size() /*y*/, K_IN))),
												inv_coefmatT(Inverse(Transpose(Bspline_matrix(samp3X_IN[0].size() /*x*/, K_IN)))),
												samp3X(samp3X_IN),
												samp3Y(samp3Y_IN),
												samp3Z(samp3Z_IN)
	{
		if (samp3X_IN.size() != samp3Y_IN.size() || samp3X_IN.size() != samp3Z_IN.size())
		{
			Print("sizes of input samples are differen", Red);
			Print(samp3X_IN.size(), Red);
			Print(samp3Y_IN.size(), Red);
			Print(samp3Z_IN.size(), Red);
			abort();
		}
	};
	ParametricInterpolation3D(){};
	void reset(const VV_d &samp3X_IN,
			   const VV_d &samp3Y_IN,
			   const VV_d &samp3Z_IN,
			   const int K_IN)
	{
		if (samp3X_IN.size() != samp3Y_IN.size() || samp3X_IN.size() != samp3Z_IN.size())
		{
			Print("sizes of input samples are differen", Red);
			Print(samp3X_IN.size(), Red);
			Print(samp3Y_IN.size(), Red);
			Print(samp3Z_IN.size(), Red);
			abort();
		}
		q1.clear();
		q2.clear();
		inv_coefmat.clear();
		inv_coefmatT.clear();
		samp3X.clear();
		samp3Y.clear();
		samp3Z.clear();
		// insert
		K = K_IN;
		s1 = samp3X_IN.size() /*y*/;
		s2 = samp3X_IN[0].size() /*x*/;
		q1 = OpenUniformKnots(samp3X_IN.size(), K_IN);
		q2 = OpenUniformKnots(samp3X_IN[0].size(), K_IN);
		inv_coefmat = Inverse(Bspline_matrix(samp3X_IN.size(), K_IN));
		inv_coefmatT = Inverse(Transpose(Bspline_matrix(samp3X_IN[0].size(), K_IN)));
		samp3X = samp3X_IN;
		samp3Y = samp3Y_IN;
		samp3Z = samp3Z_IN;
	};
	////////////////////////////////////////////////////
	V_d Nx(const double b /*y*/)
	{
		return Dot(inv_coefmatT, Bspline_vector(b, s2, K, q2));
	};
	V_d Ny(const double a /*x*/)
	{
		return Dot(Bspline_vector(a, s1, K, q1), inv_coefmat);
	};
	VV_d N(const V_d &ba)
	{
		return TensorProduct(Ny(ba[1]), Nx(ba[0]));
	};
	V_d operator()(const double b /*x*/, const double a /*y*/)
	{
		V_d shapey = Ny(a), shapex = Nx(b), ret(3, 0.);
		for (auto m = 0; m < s1; m++)
		{
			ret[0] += shapey[m] * std::inner_product(samp3X[m].begin(), samp3X[m].end(), shapex.begin(), 0.);
			ret[1] += shapey[m] * std::inner_product(samp3Y[m].begin(), samp3Y[m].end(), shapex.begin(), 0.);
			ret[2] += shapey[m] * std::inner_product(samp3Z[m].begin(), samp3Z[m].end(), shapex.begin(), 0.);
		}
		return ret;
	};
	V_d operator()(const V_d &ba)
	{
		V_d shapey = Ny(ba[1]), shapex = Nx(ba[0]), ret(3, 0.);
		for (auto m = 0; m < s1; m++)
		{
			ret[0] += shapey[m] * std::inner_product(samp3X[m].begin(), samp3X[m].end(), shapex.begin(), 0.);
			ret[1] += shapey[m] * std::inner_product(samp3Y[m].begin(), samp3Y[m].end(), shapex.begin(), 0.);
			ret[2] += shapey[m] * std::inner_product(samp3Z[m].begin(), samp3Z[m].end(), shapex.begin(), 0.);
		}
		return ret;
	};
	////////////////// DERIVATIVES /////////////////
	V_d DNx(const double b /*y*/, const int n)
	{
		return Dot(inv_coefmatT, D_Bspline_vector(b, s2, K, n, q2));
	};
	V_d DNy(const double a /*x*/, const int n)
	{
		V_d By = D_Bspline_vector(a, s1, K, n, q1);
		V_d ret(s1, 0.);
		for (auto m = 0; m < s1; m++)
			for (auto l = 0; l < s1; l++)
				ret[m] += inv_coefmat[l][m] * By[l];
		return ret;
	};
	VV_d DN(const V_d &ba, const V_i &n)
	{
		V_d shapey = DNy(ba[1], n[1]), shapex = DNx(ba[0], n[0]);
		VV_d ret(s1, V_d(s2, 0.));
		for (auto m = 0; m < s1; m++)
			for (auto j = 0; j < s2; j++)
				ret[m][j] = shapey[m] * shapex[j];
		return ret;
	};
	V_d operator()(const double b /*x*/, const double a /*y*/, const V_i &n)
	{
		V_d shapey = DNy(a, n[1]), shapex = DNx(b, n[0]), ret(3, 0.);
		for (auto m = 0; m < s1; m++)
			for (auto j = 0; j < s2; j++)
			{
				ret[0] += samp3X[m][j] * shapey[m] * shapex[j];
				ret[1] += samp3Y[m][j] * shapey[m] * shapex[j];
				ret[2] += samp3Z[m][j] * shapey[m] * shapex[j];
			}
		return ret;
	};
	V_d operator()(const V_d &ba, const std::vector<int> &n)
	{
		V_d shapey = DNy(ba[1], n[1]), shapex = DNx(ba[0], n[0]), ret(3, 0.);
		for (auto m = 0; m < s1; m++)
			for (auto j = 0; j < s2; j++)
			{
				ret[0] += samp3X[m][j] * shapey[m] * shapex[j];
				ret[1] += samp3Y[m][j] * shapey[m] * shapex[j];
				ret[2] += samp3Z[m][j] * shapey[m] * shapex[j];
			}
		return ret;
	};
	///////////////// specific shape function ////////////
	double Nx(const double b /*y*/, const int j)
	{
		return std::inner_product(inv_coefmatT[j].begin(), inv_coefmatT[j].end(), Bspline_vector(b, s2, K, q2).begin(), 0.);
	};
	double Ny(const double a /*x*/, const int m)
	{
		V_d By = Bspline_vector(a, s1, K, q1);
		double ret(0.);
		for (auto l = 0; l < s1; l++)
			ret += inv_coefmat[l][m] * By[l];
		return ret;
	};
	double N(const V_d &ba, const std::vector<int> &jm)
	{ /*2D vector access*/
		return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
	};
	double N(const V_d &ba, const int mj)
	{ /*1D access*/
		return Ny(ba[1], (int)(mj / s2) /*row*/) * Nx(ba[0], mj % s2 /*col*/);
	};
	V_d basis(const V_d &ba, const std::vector<int> &jm)
	{
		return V_d{samp3X[jm[1]][jm[0]], samp3Y[jm[1]][jm[0]], samp3Z[jm[1]][jm[0]]} * Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
	};
	V_d cross(const V_d &ba)
	{
		V_d Dshapey = DNy(ba[1], 1), Dshapex = DNx(ba[0], 1), shapey = DNy(ba[1], 0), shapex = DNx(ba[0], 0), ret_dxi(3, 0.), ret_deta(3, 0.);
		for (auto m = 0; m < s1; m++)
			for (auto j = 0; j < s2; j++)
			{
				ret_dxi[0] += samp3X[m][j] * shapey[m] * Dshapex[j];
				ret_dxi[1] += samp3Y[m][j] * shapey[m] * Dshapex[j];
				ret_dxi[2] += samp3Z[m][j] * shapey[m] * Dshapex[j];
				ret_deta[0] += samp3X[m][j] * Dshapey[m] * shapex[j];
				ret_deta[1] += samp3Y[m][j] * Dshapey[m] * shapex[j];
				ret_deta[2] += samp3Z[m][j] * Dshapey[m] * shapex[j];
			}
		return Cross(ret_dxi, ret_deta);
	};
};
////////////////////////////////////////////////////////
struct RBF_fn
{
	virtual double rbf(double r) = 0;
};
struct RBF_interp
{
	int dim, n;
	VV_d pts;
	V_d vals;
	V_d w;
	RBF_fn &fn;
	bool norm;

	RBF_interp(const VV_d &ptss /* x for f(x) */,
			   const V_d &valss /*f(x)*/,
			   RBF_fn &func,
			   bool nrbf = false)
		: dim(ptss[0].size()), n(ptss.size()), pts(ptss), vals(valss), w(n), fn(func), norm(nrbf)
	{
		double sum;
		VV_d rbf(n, V_d(n, 0.));
		V_d rhs(n);
		for (auto i = 0; i < n; i++)
		{
			sum = 0.;
			for (auto j = 0; j < n; j++)
				sum += (rbf[i][j] = fn.rbf(rad(&pts[i][0], &pts[j][0])));
			rhs[i] = norm ? sum * vals[i] : vals[i];
		}
		ludcmp lu(rbf);
		lu.solve(rhs, w);
	}
	double operator()(const V_d &pt)
	{
		double fval, sum = 0., sumw = 0.;
		if (pt.size() != dim)
			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
		for (auto i = 0; i < n; i++)
		{
			fval = fn.rbf(rad(&pt[0], &pts[i][0]));
			sumw += w[i] * fval;
			sum += fval;
		}
		return norm ? sumw / sum : sumw;
	}

	double rad(const double *p1, const double *p2)
	{
		double sum = 0.;
		for (int i = 0; i < dim; i++)
			sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
		return sqrt(sum);
	}
};

struct RBF_multiquadric : RBF_fn
{
	double r02;
	RBF_multiquadric(double scale = 1.) : r02(scale * scale) {}
	double rbf(double r) { return sqrt(r * r + r02); }
};
struct RBF_thinplate : RBF_fn
{
	double r0;
	RBF_thinplate(double scale = 1.) : r0(scale) {}
	double rbf(double r) { return r <= 0. ? 0. : r * r * std::log(r / r0); }
};
struct RBF_gauss : RBF_fn
{
	double r0;
	RBF_gauss(double scale = 1.) : r0(scale) {}
	double rbf(double r) { return exp(-0.5 * r * r / (r0 * r0)); }
};
struct RBF_inversemultiquadric : RBF_fn
{
	double r02;
	RBF_inversemultiquadric(double scale = 1.) : r02(scale * scale) {}
	double rbf(double r) { return 1. / sqrt(r * r + r02); }
};
struct Shep_interp
{
	int dim, n;
	const VV_d &pts;
	const V_d &vals;
	double pneg;

	Shep_interp(VV_d &ptss, V_d &valss, double p = 2.)
		: dim(ptss[0].size()), n(ptss.size()), pts(ptss), vals(valss), pneg(-p) {}

	double interp(V_d &pt)
	{
		double r, w, sum = 0., sumw = 0.;
		if (pt.size() != dim)
			throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
		for (auto i = 0; i < n; i++)
		{
			if ((r = rad(&pt[0], &pts[i][0])) == 0.)
				return vals[i];
			sum += (w = pow(r, pneg));
			sumw += w * vals[i];
		}
		return sumw / sum;
	}

	double rad(const double *p1, const double *p2)
	{
		double sum = 0.;
		for (auto i = 0; i < dim; i++)
			sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
		return sqrt(sum);
	}
};

// #ifdef INCL_lapacke
// struct LU_LAPACK{
// public:
//   char TRANS = 'N';
//   int INFO;
//   int LDA;
//   int LDB;
//   int N;
//   int NRHS = 1;
//   std::vector<int> IPIV;
//   V_d Flatten_A;

//   LU_LAPACK(const std::vector<V_d>& A/* must be square matrix */):
//     Flatten_A(Flatten(A)),
//     INFO(A.size()),
//     LDA(A.size()),
//     LDB(A.size()),
//     N(A.size()),
//     IPIV(A.size())
//   {
//     std::cout << "compute the LU factorization..." << std::endl;
//     //void LAPACK_dgetrf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda, lapack_int* ipiv, lapack_int *info );
//     LAPACK_dgetrf(&N, &N, &Flatten_A[0], &LDA, &IPIV[0], &INFO);
//     if(INFO)
//       std::cout << "an error occured : "<< INFO << std::endl;
//   };
//   /////////////////////////////////////////////////////////////
//   V_d solve(const V_d& B){
//     V_d ret(B);
//     // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
//     std::cout << "solving the system..."<< std::endl;
//     // void LAPACK_dgetrs( char* trans, lapack_int* n, lapack_int* nrhs, const double* a, lapack_int* lda, const lapack_int* ipiv,double* b, lapack_int* ldb, lapack_int *info );
//     dgetrs_(&TRANS, &N, &NRHS, &Flatten_A[0], &LDA, &IPIV[0], &ret[0], &LDB,&INFO);
//     if(INFO)
//       {
// 	std::cout << "an error occured : "<< INFO << std::endl;
// 	abort();
//       }
//     return ret;
//   };
//   /////////////////////////////////////////////////////////////
//   void solve(const V_d& B, V_d& ret){
//     // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
//     std::cout << "solving the system..."<< std::endl;
//     // void LAPACK_dgetrs( char* trans, lapack_int* n, lapack_int* nrhs, const double* a, lapack_int* lda, const lapack_int* ipiv,double* b, lapack_int* ldb, lapack_int *info );
//     dgetrs_(&TRANS, &N, &NRHS, &Flatten_A[0], &LDA,&IPIV[0], &ret[0], &LDB,&INFO);
//     if(INFO)
//       {
// 	std::cout << "an error occured : "<< INFO << std::endl;
// 	abort();
//       }
//   };
// };
// #endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//=========================================================
// void Load(const std::string& filename, const V_s& sep, std::vector<std::vector<V_d>>& mat){
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     std::cout << Red << "The file can not be opened" << reset << std::endl;
//   std::string str, S, cell;
//   double v;
//   std::size_t found, s(sep.size()), min_len, left, right;
//   bool loop;
//   int row_counter(0), col_counter(0);
//   while(!strm.eof( ) && mat.size() > row_counter)
//     {
//       std::getline(strm,str);
//       if(!str.empty())
// 	{
// 	  if((left = str.find("{"))!=std::string::npos);
// 	  if((right = str.find("}"))!=std::string::npos);
// 	  cell = str.substr(left,right);
// 	  str = str.substr(right+1);

// 	  V_d row_vector;
// 	  loop = true;
// 	  col_counter = 0;

// 	  while(loop && mat[row_counter].size() > col_counter)
// 	    {
// 	      min_len = 999999;
// 	      for(auto i=0; i<s; i++)
// 		if((found = cell.find(sep[i]))!=std::string::npos)
// 		  if(min_len > found)
// 		    {
// 		      min_len = found;
// 		      S = sep[i];
// 		    }
// 	      col_counter++;
// 	    }
// 	  row[row_counter].emplace_back(row_vector);
// 	  row_counter++;
// 	}
//     }
//   strm.close();
// };
//=========================================================
template <class Ttype>
void load(std::ifstream &in, std::vector<Ttype> &vec, const int row_size)
{
	vec.clear();
	if (!in)
	{
		std::cout << "入力ファイルが開けません" << std::endl;
	}
	std::string str;
	int i;
	Ttype tmp;
	i = 0;
	getline(in, str);
	std::cout << std::string(40, '=') << std::endl;
	std::cout << "この行は，" << str << std::endl;
	std::cout << std::string(40, '-') << std::endl;
	std::cout << "読み込み結果：" << std::endl;
	for (size_t j = 0; j < row_size; j++)
	{
		if (str.find(',') != std::string::npos)
		{
			i = str.find(',');
			//	std::cout << "\',\'は" << i << "番目にあります" << endl;
			//	std::cout << str << "　を" << i << "まで読み込みます -> " << str.substr(0,i) << " -> これを数値に変換すると ->";
			tmp = (Ttype)stod(str.substr(0, i));
			//	std::cout << std::setprecision(10) << tmp << endl;
			str = str.substr(i + 1);
			//	std::cout << i+1 << "番目以降を抜き出し文字列を置き換えます -> " << str << endl;
		}
		else
		{
			//	std::cout << "\',\'は" << "見つかりませんでした -> ";
			tmp = (Ttype)stod(str);
			//	std::cout << str << "　を" << "最後まで読み込みます -> " << std::setprecision(10) << tmp << endl;
		}
		vec.emplace_back(tmp);
		//      std::cout << tmp << "を代入 -> ";
		std::cout << "vec[" << j << "] = " << std::setprecision(5) << vec[j] << std::endl;
	}
	std::cout << std::string(40, '-') << std::endl;
};
//===========================================================
template <class Ttype>
void load(std::ifstream &in, std::vector<std::vector<Ttype>> &vec, const int row_size)
{
	vec.clear();
	vec.resize(row_size, std::vector<Ttype>(0));
	if (!in)
	{
		std::cout << "入力ファイルが開けません" << std::endl;
	}
	std::string str;
	int i;
	Ttype tmp;
	for (int line = 0;; line++)
	{
		i = 0;
		getline(in, str);
		if (in.eof())
			break;
		std::cout << std::string(40, '=') << std::endl;
		std::cout << "この行は，" << str << std::endl;
		std::cout << std::string(40, '-') << std::endl;
		std::cout << "読み込み結果：" << std::endl;
		for (int j = 0; j < row_size; j++)
		{
			if (str.find(',') != std::string::npos)
			{
				i = str.find(',');
				//	std::cout << "\',\'は" << i << "番目にあります" << endl;
				//	std::cout << str << "　を" << i << "まで読み込みます -> " << str.substr(0,i) << " -> これを数値に変換すると ->";
				tmp = (Ttype)stod(str.substr(0, i));
				//	std::cout << std::setprecision(10) << tmp << endl;
				str = str.substr(i + 1);
				//	std::cout << i+1 << "番目以降を抜き出し文字列を置き換えます -> " << str << endl;
			}
			else
			{
				//	std::cout << "\',\'は" << "見つかりませんでした -> ";
				tmp = (Ttype)stod(str);
				//	std::cout << str << "　を" << "最後まで読み込みます -> " << std::setprecision(10) << tmp << endl;
			}
			vec[j].emplace_back(tmp);
			//      std::cout << tmp << "を代入 -> ";
			std::cout << "vec[" << j << "][" << line << "] = " << std::setprecision(5) << vec[j][line] << std::endl;
		}
		std::cout << std::string(40, '-') << std::endl;
	};
};
//=====================================================================
//=================== Mathematica output loader =======================
//=====================================================================
V_s StringSplit(const std::string &strIN, const V_s &SEP)
{
	std::string str(strIN), foundFirstSep; /* expexting "1,2,3,4,5" or "1, 2,3, 4,5"*/
	V_s ret(0), tmp(0);
	while (!str.empty())
	{
		size_t foundFirst(999999);
		bool bfound(false);
		for (const auto &sep : SEP)
		{
			size_t found = str.find(sep);
			if (!sep.empty()									  /* ignore sep, "" */
				&& found != std::string::npos					  /* if sep is found */
				&& (found < foundFirst)) /* choose closest one */ // do not use "<=" instead
			{
				foundFirst = found;
				foundFirstSep = sep;
				bfound = true;
			}
		}
		if (bfound)
		{
			tmp.emplace_back(str.substr(0, foundFirst));
			str = str.substr(foundFirst + foundFirstSep.length());
		}
		else
		{
			tmp.emplace_back(str);

			// remove first and last cells if they are empty strings
			while (tmp.size() > 0)
				if ((*tmp.begin()).empty())
					tmp.erase(tmp.begin());
				else if ((*tmp.rbegin()).empty() ||
						 (*tmp.rbegin()).find_first_not_of(" ") == std::string::npos /* there is no string except spaces */ ||
						 (*tmp.rbegin()).find_first_not_of("\r") == std::string::npos /* there is no string except \r */ ||
						 (*tmp.rbegin()).find_first_not_of("\n") == std::string::npos /* there is no string except \n */ ||
						 (*tmp.rbegin()).find_first_not_of("\t") == std::string::npos /* there is no string except \t */)
					tmp.pop_back();
				else
					break;

			return tmp;
		}
	}

	// remove first and last cells if they are empty strings

	while (tmp.size() > 0)
		if ((*tmp.begin()).empty())
			tmp.erase(tmp.begin());
		else if ((*tmp.rbegin()).empty() ||
				 (*tmp.rbegin()).find_first_not_of(" ") == std::string::npos /* there is no string except spaces */ ||
				 (*tmp.rbegin()).find_first_not_of("\r") == std::string::npos /* there is no string except \r */ ||
				 (*tmp.rbegin()).find_first_not_of("\n") == std::string::npos /* there is no string except \n */ ||
				 (*tmp.rbegin()).find_first_not_of("\t") == std::string::npos /* there is no string except \t */)
			tmp.pop_back();
		else
			break;

	return tmp;
};
std::string StringTrim(const std::string &strIN, const V_s &SEP)
{
	V_s str = StringSplit(strIN, SEP);
	std::string ret;
	for (const auto &s : str)
		ret += s;
	return ret;
};
V_s StringTrim(const V_s &strIN, const V_s &SEP)
{
	V_s strOUT(strIN.size());
	for (size_t i = 0; i < strOUT.size(); i++)
		strOUT[i] = StringTrim(strIN[i], SEP);
	return strOUT;
};
std::string StringJoin(const V_s &strIN, const std::string sep = "")
{
	std::string ret(*strIN.begin());
	for (size_t i = 1; i < strIN.size(); i++)
		ret += sep + strIN[i];
	return ret;
};
std::string Directory(const std::string &f)
{
	auto dir = StringSplit(f, {"/"});
	dir.pop_back();
	return "/" + StringJoin(dir, "/") + "/";
};
////////////////////////////////////////////////////////////////////////
// 2021/03/31導入
using V_S = std::vector<std::string>;

std::map<std::string, V_S> parseJSON(const std::string &str_IN)
{
	std::map<std::string, V_S> map_S_S;
	//
	V_S SEP = {" ", "\t", "\n", "\r"};
	std::string trimed = StringTrim(str_IN, SEP);
	std::string L = "[", R = "]", Col = ":", Cam = ",", Lb = "{", Rb = "}";
	bool is_array = false;
	std::string array = "", value = "";
	V_S vstr(0);
	for (auto it = trimed.begin(); it != trimed.end(); it++)
	{
		if (is_array || L.find(*it) != std::string::npos)
		{
			//[から始まるarrayの場合は，ここでキャッチする
			if ((R.find(*it) != std::string::npos))
			{
				is_array = false;
				vstr.emplace_back(array);
				array.clear();
			}
			else
			{
				if (is_array)
					array.push_back(*it);
				is_array = true;
			}
		}
		else
		{
			// arrayでない場合で，{},:などで始まる場合はここでキャッチする
			if ((Lb.find(*it) != std::string::npos) ||
				(Cam.find(*it) != std::string::npos) ||
				(Col.find(*it) != std::string::npos) ||
				(Rb.find(*it) != std::string::npos))
			{
				if (!value.empty())
				{
					// std::cout << red << value << std::endl;
					vstr.emplace_back(value);
					value.clear();
				}
			}
			else
			{
				value.push_back(*it);
			}
		}
	}
	// 得られたstd::vector<std::string>をmapに保存する．
	for (auto i = 0; i < vstr.size() - 1; i += 2)
		map_S_S[StringTrim(vstr[i], {"\""})] = StringSplit(StringTrim(vstr[i + 1], {"\""}), {","});

	return map_S_S;
};

/* ------------------------------------------------------ */

struct JSON
{
	using V_S = std::vector<std::string>;
	std::map<std::string, V_S> map_S_S;

	// JSON(const std::string &str_IN) : map_S_S() { this->parse(str_IN); };

	JSON(){};
	JSON &operator=(const JSON &other)
	{
		this->map_S_S = other.map_S_S;
		return *this;
	};
	JSON(const std::ifstream &istrm) : map_S_S()
	{
		if (!istrm.is_open())
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not open");
		else
		{
			std::stringstream ss;
			ss << istrm.rdbuf();
			map_S_S = parseJSON(ss.str());
		}
	};

	JSON(const std::string &filename) : map_S_S()
	{
		std::ifstream istrm(filename);
		if (!istrm.is_open())
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not open");
		else
		{
			std::stringstream ss;
			ss << istrm.rdbuf();
			map_S_S = parseJSON(ss.str());
		}
	};

	std::map<std::string, V_S> operator()() const
	{
		return this->map_S_S;
	};

	V_S &operator[](const std::string &key) /*変更を許す*/
	{
		return this->map_S_S[key];
	};

	bool find(const std::string &key) const
	{
		return (this->map_S_S.find(key) != this->map_S_S.end());
	};

	// V_s operator[](const std::string &key) const
	// {
	// 	if (this->map_S_S.find(key) != this->map_S_S.end())
	// 		return this->map_S_S.at(key);
	// 	else
	// 		return {};
	// };
};

/* ------------------------------------------------------ */

std::ofstream &operator<<(std::ofstream &stream, const JSON &json)
{
	stream << "{\n";
	auto count = 0;
	for (const auto &[key, v] : json.map_S_S)
	{
		count++;
		if (v.size() == 1)
		{
			stream << "    "
				   << "\"" << key << "\":"
				   << " "
				   << "\"" << v[0] << "\"";
		}
		else
		{
			stream << "    "
				   << "\"" << key << "\":"
				   << " "
				   << "[";
			for (auto i = 0; i < v.size(); i++)
			{
				stream << "\"" << v[i] << "\"";
				if (i != v.size() - 1)
					stream << ",";
			}
			stream << "]";
		}
		if (json.map_S_S.size() != count)
			stream << ",\n";
		else
			stream << "\n";
	}
	stream << "}";
	return stream;
};
/* ------------------------------------------------------ */

struct JSONoutput
{
	std::map<std::string, T6d> map_S_T6d;
	std::map<std::string, Tddd> map_S_Tddd;
	std::map<std::string, std::vector<double>> map_S_D;
	std::map<std::string, std::vector<int>> map_S_I;
	std::map<std::string, std::vector<std::string>> map_S_S;
	JSONoutput(){};

	void push(std::string s, const T6d D)
	{
		if (this->map_S_T6d.find(s) != this->map_S_T6d.end())
			this->map_S_T6d[s] = D;
		else
			this->map_S_T6d[s] = {D};
	};
	void push(std::string s, const Tddd D)
	{
		if (this->map_S_Tddd.find(s) != this->map_S_Tddd.end())
			this->map_S_Tddd[s] = D;
		else
			this->map_S_Tddd[s] = {D};
	};
	void push(std::string s, double D)
	{
		if (this->map_S_D.find(s) != this->map_S_D.end())
			this->map_S_D[s].emplace_back(D);
		else
			this->map_S_D[s] = {D};
	};
	void push(std::string s, int I)
	{
		if (this->map_S_I.find(s) != this->map_S_I.end())
			this->map_S_I[s].emplace_back(I);
		else
			this->map_S_I[s] = {I};
	};
	void push(std::string s, std::string S)
	{
		if (this->map_S_S.find(s) != this->map_S_S.end())
			this->map_S_S[s].emplace_back(S);
		else
			this->map_S_S[s] = {S};
	};

	void output(std::ofstream &os)
	{
		os << "{" << std::endl;
		{
			int size = this->map_S_T6d.size();
			int i = 0;
			for (const auto &[key, V] : this->map_S_T6d)
			{
				auto [V0, V1, V2, V3, V4, V5] = V;
				os << "\"" << key << "\":[";
				os << std::setprecision(15) << V0 << ", "
				   << std::setprecision(15) << V1 << ", "
				   << std::setprecision(15) << V2 << ", "
				   << std::setprecision(15) << V3 << ", "
				   << std::setprecision(15) << V4 << ", "
				   << std::setprecision(15) << V5 << "]";
				if (i++ != size - 1)
					os << ",\n";
			}
		}
		{
			int size = this->map_S_Tddd.size();
			int i = 0;
			for (const auto &[key, V] : this->map_S_Tddd)
			{
				auto [V0, V1, V2] = V;
				os << "\"" << key << "\":[";
				os << std::setprecision(15) << V0 << ", "
				   << std::setprecision(15) << V1 << ", "
				   << std::setprecision(15) << V2 << "]";
				if (i++ != size - 1)
					os << ",\n";
			}
		}
		{
			int size = this->map_S_D.size();
			int i = 0;
			for (const auto &[key, V] : this->map_S_D)
			{
				os << "\"" << key << "\":[";
				for (auto j = 0; j < V.size(); ++j)
					os << std::setprecision(15) << V[j] << (j != V.size() - 1 ? "," : "]");
				if (i++ != size - 1)
					os << ",\n";
			}
		}
		{
			int size = this->map_S_I.size();
			int i = 0;
			for (const auto &[key, V] : this->map_S_I)
			{
				os << "\"" << key << "\":[";
				for (auto j = 0; j < V.size(); ++j)
					os << V[j] << (j != V.size() - 1 ? "," : "]");
				if (i++ != size - 1)
					os << ",\n";
			}
		}
		os << "\n}" << std::endl;
	};
};

//============================================================
//========================= Load  ============================
//============================================================
// std::vector<V_d> Load(const std::string& filename, const V_s& sep){
//   std::vector<V_d> mat;
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     std::cout << Red << filename << " can not be opened" << reset << std::endl;
//   else
//     std::cout << Blue << filename << " is opened" << reset << std::endl;

//   std::string str, S;
//   double v;
//   std::size_t found, s(sep.size()), min_len;
//   bool loop;
//   while(!strm.eof( ))
//     {
//       std::getline(strm,str);
//       if(!str.empty())
// 	{
// 	  V_d row_vector;
// 	  loop = true;
// 	  while(loop)
// 	    {
// 	      min_len = 999999;
// 	      for(auto i=0; i<s; i++)
// 		if((found = str.find(sep[i]))!=std::string::npos)
// 		  if(min_len > found)
// 		    {
// 		      min_len = found;
// 		      S = sep[i];
// 		    }

// 	      if(min_len!=999999/*if not 0*/)
// 		{
// 		  row_vector.push_back((double)stod(str.substr(0,min_len)));
// 		  str = str.substr(min_len + S.length());
// 		}
// 	      else
// 		{
// 		  row_vector.push_back((double)stod(str));
// 		  loop = false;
// 		}
// 	    }
// 	  mat.emplace_back(row_vector);
// 	}
//     }
//   strm.close();
//   return mat;
// };
//////////////////////////////
void Load(const std::string &filename, std::vector<V_s> &ret_mat, const V_s &SEP)
{
	ret_mat.clear();
	std::ifstream strm(filename, std::ios::in);
	if (!strm)
	{
		std::stringstream ss;
		ss << "The file can not be opened: " << filename;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
	}
	else
		std::cout << Blue << filename << " is opened" << reset << std::endl;

	std::string read;
	int row_counter(0);
	while (!strm.eof()) // loop until the end of file
	{
		read.clear();
		std::getline(strm, read);
		if (!StringTrim(read, {" "}).empty())
		{
			ret_mat.emplace_back(StringTrim(StringSplit(read, SEP), {" ", "\t", "\n", "\r"}));
		}
	}
	strm.close();
	std::cout << Blue << filename << " is closed" << reset << std::endl;
};

std::vector<V_s> Load(const std::string &filename, const V_s &SEP)
{
	std::vector<V_s> ret;
	Load(filename, ret, SEP);
	return ret;
};
////////////////////////////////
V_d string_to_vector_double(const std::string &strIN, const V_s &SEP)
{
	V_s str = StringSplit(strIN, SEP);
	V_d ret;
	for (size_t i = 0; i < str.size(); i++)
		if (!str[i].empty())
			ret.emplace_back(stod(str[i]));

	return ret;
	// std::string str(strIN);/* expexting "1,2,3,4,5" or "1, 2,3, 4,5"*/
	// V_d ret(0);
	// size_t found;
	// int i(0);
	// while(!str.empty())
	//   {
	//     for(auto sep: SEP)
	//       if((found = str.find(sep))!=std::string::npos)
	// 	  {
	// 	    // std::cout << str << std::endl;
	// 	    // std::cout << found << std::endl;
	// 	    // std::cout << sep << std::endl;
	// 	    // std::cout << str << std::endl;
	// 	    ret.emplace_back(stod(
	// 			       str.substr((i==0 ? 0 : found+sep.length()), found)
	// 			       ));
	// 	    str = str.substr(found+sep.length());
	// 	    break;
	// 	  }
	// 	else
	// 	  {
	// 	    ret.emplace_back(stod(
	// 			       str.substr((i==0 ? 0 : found+sep.length()))
	// 			       ));
	// 	    return ret;
	// 	  }
	//     i++;
	//   }
	// return ret;
};
///////////////
V_s to_cell(const std::string &strIN)
{
	std::string str(strIN); /* expexting "{1,2,3,4,5} {1,2,3,4,5} {1,2,3,4,5}"*/
	V_s ret(0);
	size_t found, right, left;
	int i(0);
	while (!str.empty())
	{
		if ((left = str.find("{")) != std::string::npos)
		{
			if ((right = str.find("}")) != std::string::npos)
			{
				std::string cell = str.substr(left, right - left + 1);
				ret.emplace_back(cell);
				str = str.substr(right + 1);
				// std::cout << str << std::endl;
			}
			else
				return ret;
		}
		else
			return ret;
	}
	return ret;
};
///////////////////////////

/////////////////////////////////
void Load(const std::string &filename, VVV_d &vvv)
{
	std::ifstream strm(filename, std::ios::in);
	if (!strm)
	{
		std::stringstream ss;
		ss << "The file can not be opened: " << filename;
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
	}
	std::string str;
	while (!strm.eof())
	{
		std::getline(strm, str);

		// std::cout <<  str << std::endl;
		// str = str.substr(str.find_last_of("{")+1,  str.find_last_of("}")-1);
		// std::cout <<  str << std::endl;

		V_s str_vec = to_cell(str); /* expecting {"{1,2,3,5,5}","{1,2,3,5,5}","{1,2,3,5,5}"} */
		VV_d row_vec(0);
		for (const auto &v : str_vec)
		{
			std::string tmp = v.substr(v.find_first_of("{") + 1, v.find_last_of("}") - 1);
			tmp.find("");
			V_d vec_doub = string_to_vector_double(/* expecting {1,2,3,5,5} */ tmp, {",", ", "});
			row_vec.emplace_back(vec_doub);
		}
		vvv.emplace_back(row_vec);
	}
	strm.close();
};
// /////////////////////////////////
// void Load(const std::string& filename, std::vector<V_d>& vv)
// {
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     {
//       std::cout << Red << "The file can not be opened: " << filename << reset << std::endl;
//       std::cin.ignore();
//       abort();
//     }
//   std::string str;
//   while(!strm.eof( ))
//     {
//       std::getline(strm,str);

//       // std::cout <<  str << std::endl;
//       // str = str.substr(str.find_last_of("{")+1,  str.find_last_of("}")-1);
//       // std::cout <<  str << std::endl;

//       V_d vec_doub = string_to_vector_double(/* expecting {1,2,3,5,5} */str,{",",", "});
//       vv.emplace_back(vec_doub);
//     }
//   strm.close();
// };
//============================================================
//============================================================
template <class T>
std::vector<T> PointsToSurface(T t0, T t1, const std::vector<std::vector<T>> &a)
{
	T EPS = 1E-15;

	if (a.size() == 3)
	{
		T t2 = 1. - t1 - t0;
		T x0 = 1. - t0 - 2. * t1, x1 = 1. - t1 - 2. * t2, x2 = 1. - t2 - 2. * t0;

		if (t0 < EPS)
		{
			return a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
		}
		else if (t1 < EPS)
		{
			return a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
		}
		else if (t2 < EPS)
		{
			return a[2] * (1. - x2) / 2. + a[0] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
		}
		else
		{
			std::vector<T> l0 = a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
			std::vector<T> l1 = a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
			std::vector<T> l2 = a[2] * (1. - x2) / 2. + a[0] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
			std::vector<T> tmp = {1 / t0, 1 / t1, 1 / t2};
			return Dot(tmp / (1 / t0 + 1 / t1 + 1 / t2), {l0, l1, l2});
		}
	}
	else
	{
		T t2 = 1. - t0, t3 = 1. - t1;
		T x0 = 2. * t3 - 1., x1 = 2. * t0 - 1., x2 = 2. * t1 - 1., x3 = 2. * t2 - 1.;

		if (t0 < EPS)
		{
			return a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
		}
		else if (t1 < EPS)
		{
			return a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
		}
		else if (t2 < EPS)
		{
			return a[2] * (1. - x2) / 2. + a[3] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
		}
		else if (t3 < EPS)
		{
			return a[3] * (1. - x3) / 2. + a[0] * (1. + x3) / 2. /*+p2A3*(x3^2-1)*/;
		}
		else
		{
			std::vector<T> l0 = a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
			std::vector<T> l1 = a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
			std::vector<T> l2 = a[2] * (1. - x2) / 2. + a[3] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
			std::vector<T> l3 = a[3] * (1. - x3) / 2. + a[0] * (1. + x3) / 2. /*+p2A3*(x3^2-1)*/;
			std::vector<T> tmp = {1 / t0, 1 / t1, 1 / t2, 1 / t3};
			return Dot(tmp / (1 / t0 + 1 / t1 + 1 / t2 + 1 / t3), {l0, l1, l2, l3});
		}
	}
};
//============================================================
template <class Type>
std::vector<Type> RotateLeft(const std::vector<Type> &vecs, const int n)
{ // 2020/03/22
	std::vector<Type> ret(vecs);
	std::rotate(ret.begin(), ret.begin() + n, ret.end());
	return ret;
}
template <class T>
std::tuple<T, T, T> Reverse(const std::tuple<T, T, T> &vecs)
{ // 2022年3月21日
	return {std::get<2>(vecs), std::get<1>(vecs), std::get<0>(vecs)};
};
Tddd RotateLeft(const Tddd &vecs, int n)
{ // 2021/12/05
	n = n % 3;
	if (n == 0)
		return vecs;
	else if (n == 1)
		return {std::get<1>(vecs), std::get<2>(vecs), std::get<0>(vecs)};
	else
		return {std::get<2>(vecs), std::get<0>(vecs), std::get<1>(vecs)};
}
T3Tddd RotateLeft(const T3Tddd &vecs, int n)
{ // 2021/12/05
	n = n % 3;
	if (n == 0)
		return vecs;
	else if (n == 1)
		return {std::get<1>(vecs), std::get<2>(vecs), std::get<0>(vecs)};
	else
		return {std::get<2>(vecs), std::get<0>(vecs), std::get<1>(vecs)};
}
T3Tddd RotateRight(const T3Tddd &vecs, int n)
{ // 2021/12/05
	n = n % 3;
	if (n == 0)
		return vecs;
	else if (n == 1)
		return {std::get<2>(vecs), std::get<0>(vecs), std::get<1>(vecs)};
	else
		return {std::get<1>(vecs), std::get<2>(vecs), std::get<0>(vecs)};
}
Tddd RotateRight(const Tddd &vecs, int n)
{ // 2021/12/05
	n = n % 3;
	if (n == 0)
		return vecs;
	else if (n == 1)
		return {std::get<2>(vecs), std::get<0>(vecs), std::get<1>(vecs)};
	else
		return {std::get<1>(vecs), std::get<2>(vecs), std::get<0>(vecs)};
}
template <class Type>
std::vector<Type> RotateRight(const std::vector<Type> &vecs, const int n)
{ // 2021/04/07
	std::vector<Type> ret(vecs);
	std::rotate(ret.rbegin(), ret.rbegin() + n, ret.rend());
	return ret;
}

template <typename T>
std::vector<T> Drop(std::vector<T> ret /*copy*/, std::vector<int> vec_n)
{ // 2020/03/22
	std::sort(vec_n.begin(), vec_n.end());
	for (size_t i = 0; i < vec_n.size(); i++)
		ret.erase(ret.begin() + (vec_n[i] - i));
	return ret;
};

int Position(std::vector<int> vecIN, const int n)
{ // 2020/03/22
	std::vector<int>::iterator it = std::find(vecIN.begin(), vecIN.end(), n);
	if (it == vecIN.end())
		std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << reset << std::endl;

	return std::distance(vecIN.begin(), it);
};

template <class T>
std::vector<T> Union(const std::vector<T> &vecs, const std::vector<T> &vecs2)
{ // Union without sorting!
	std::vector<T> ret = vecs;
	for (const auto &v : vecs2)
		if (!MemberQ(ret, v))
			ret.emplace_back(v);
	return ret;
};

template <class T>
std::vector<T> Union(const std::vector<std::vector<T>> &vecs)
{ // Union without sorting!
	std::vector<T> ret(0);
	for (const auto &vec : vecs)
		for (const auto &v : vec)
			if (!MemberQ(ret, v))
				ret.emplace_back(v);

	return ret;
};

template <class Type>
std::vector<std::vector<Type>> SortVectorChain(const std::vector<std::vector<Type>> &vecs)
{ // 2020/03/22
	std::vector<Type> v = *vecs.begin();
	std::vector<std::vector<Type>> rest = Drop(vecs, {0});
	Type knot = *v.rbegin();
	int n = 0;
	std::vector<std::vector<Type>> ret = {v};
	while (rest.size() != 0)
	{
		std::vector<Type> tmp = rest[n];
		// std::cout << "knot= " << knot << std::endl;
		// std::cout << "tmp= " << tmp << std::endl;
		// std::cout << "MemberQ("<< tmp << "," << knot << ") = " << MemberQ(tmp, knot) << std::endl;
		if (MemberQ(tmp, knot))
		{
			// std::cout << tmp << std::endl;
			// std::cout << Green << "RotateLeft(tmp, Position(tmp, knot))= " << RotateLeft(tmp, Position(tmp, knot)) << reset << std::endl;
			//	  v.clear();
			v = RotateLeft(tmp, Position(tmp, knot));
			ret.emplace_back(v);
			knot = *v.rbegin();

			rest = Drop(rest, {n});

			n = 0;
		}
		else
			n++;

		// std::cout << Blue << "ret= " << ret << reset << std::endl;
	}
	return ret;
}

template <class T>
std::vector<T> Reverse(std::vector<T> vec)
{
	std::reverse(vec.begin(), vec.end());
	return vec;
};
//===========================================================
// double normalDirDistanceFromTriangle(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a)
// {
// 	return Dot(TriangleNormal(p0, p1, p2), p0 - a);
// };
// V_d vectorToTriangle(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a)
// {
// 	auto n = TriangleNormal(p0, p1, p2);
// 	return n * Dot(n, p0 - a);
// };
// // line a to b
// double factorOfVectorToReachTriangle(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a, const V_d &b)
// {
// 	//オーダーが匹敵する物を選ぶ
// 	double log_b_a = log10(Norm(b - a));
// 	double diff0 = std::abs(log10(Norm(p0 - a) - log_b_a));
// 	double diff1 = std::abs(log10(Norm(p1 - a) - log_b_a));
// 	double diff2 = std::abs(log10(Norm(p2 - a) - log_b_a));
// 	V_d n = TriangleNormal(p0, p1, p2);
// 	if (diff0 < diff1 && diff0 < diff2)
// 		return Dot(p0 - a, n) / Dot(b - a, n);
// 	else if (diff1 < diff0 && diff1 < diff2)
// 		return Dot(p1 - a, n) / Dot(b - a, n);
// 	else
// 		return Dot(p2 - a, n) / Dot(b - a, n);
// };
// ///////////////////////////////////////////////////////////
// int isPointingTriangle(const V_d &p0, const V_d &p1, const V_d &p2,
// 					   const V_d &a, const V_d &b)
// {
// 	auto d = factorOfVectorToReachTriangle(p0, p1, p2, a, b);

// 	if (!std::isfinite(d))
// 		return false; //nan

// 	V_d ps = a + (b - a) * d;

// 	// double e = 1E-20;
// 	// //ポリゴン頂点の最大最小でチェック
// 	// if (((p0[0] - e > ps[0] && p1[0] - e > ps[0] && p2[0] - e > ps[0]) || (p0[0] + e < ps[0] && p1[0] + e < ps[0] && p2[0] + e < ps[0])) ||
// 	//     ((p0[1] - e > ps[1] && p1[1] - e > ps[1] && p2[1] - e > ps[1]) || (p0[1] + e < ps[1] && p1[1] + e < ps[1] && p2[1] + e < ps[1])) ||
// 	//     ((p0[2] - e > ps[2] && p1[2] - e > ps[2] && p2[2] - e > ps[2]) || (p0[2] + e < ps[2] && p1[2] + e < ps[2] && p2[2] + e < ps[2])))
// 	//   return false;/*面の最大最小範囲にすら入れていない*/

// 	auto ps_p0 = p0 - ps;
// 	auto ps_p1 = p1 - ps;
// 	auto ps_p2 = p2 - ps;

// 	ps_p0 = ps_p0 / Norm(ps_p0);
// 	ps_p1 = ps_p1 / Norm(ps_p1);
// 	ps_p2 = ps_p2 / Norm(ps_p2);

// 	V_d n = TriangleNormal(p0, p1, p2);
// 	if (Dot(Cross(ps_p0, ps_p1), n) >= 0. &&
// 		Dot(Cross(ps_p1, ps_p2), n) >= 0. &&
// 		Dot(Cross(ps_p2, ps_p0), n) >= 0.)
// 		return true;
// 	else
// 		return false;
// };
// ///////////////////////////////////////////////////////////
// int isIntersectingSurface(const V_d &p0, const V_d &p1, const V_d &p2,
// 						  const V_d &a, const V_d &b)
// {
// 	/* 0:頂点の最大最小の範囲の外で，片方にa,bgがある */
// 	/* 1:拡大した面には入れているが，多角形の頂点の最大最小範囲にすら入れていない */
// 	/* 2:a,bは多角形の面と交差していないが，かなり惜しい */
// 	/* 3:a,bは多角形の面と交差 */
// 	// double e = 1E-11; //これを1E-14とすることで，干渉のチェックが行われるようになる場合があった．
// 	// if ((p0[0] - e > a[0] && p1[0] - e > a[0] && p2[0] - e > a[0] && p0[0] - e > b[0] && p1[0] - e > b[0] && p2[0] - e > b[0]) ||
// 	// 	(p0[0] + e < a[0] && p1[0] + e < a[0] && p2[0] + e < a[0] && p0[0] + e < b[0] && p1[0] + e < b[0] && p2[0] + e < b[0]) ||
// 	// 	(p0[1] - e > a[1] && p1[1] - e > a[1] && p2[1] - e > a[1] && p0[1] - e > b[1] && p1[1] - e > b[1] && p2[1] - e > b[1]) ||
// 	// 	(p0[1] + e < a[1] && p1[1] + e < a[1] && p2[1] + e < a[1] && p0[1] + e < b[1] && p1[1] + e < b[1] && p2[1] + e < b[1]) ||
// 	// 	(p0[2] - e > a[2] && p1[2] - e > a[2] && p2[2] - e > a[2] && p0[2] - e > b[2] && p1[2] - e > b[2] && p2[2] - e > b[2]) ||
// 	// 	(p0[2] + e < a[2] && p1[2] + e < a[2] && p2[2] + e < a[2] && p0[2] + e < b[2] && p1[2] + e < b[2] && p2[2] + e < b[2]))
// 	// 	return 0;

// 	if ((p0[0] > a[0] && p1[0] > a[0] && p2[0] > a[0] && p0[0] > b[0] && p1[0] > b[0] && p2[0] > b[0]) ||
// 		(p0[0] < a[0] && p1[0] < a[0] && p2[0] < a[0] && p0[0] < b[0] && p1[0] < b[0] && p2[0] < b[0]) ||
// 		(p0[1] > a[1] && p1[1] > a[1] && p2[1] > a[1] && p0[1] > b[1] && p1[1] > b[1] && p2[1] > b[1]) ||
// 		(p0[1] < a[1] && p1[1] < a[1] && p2[1] < a[1] && p0[1] < b[1] && p1[1] < b[1] && p2[1] < b[1]) ||
// 		(p0[2] > a[2] && p1[2] > a[2] && p2[2] > a[2] && p0[2] > b[2] && p1[2] > b[2] && p2[2] > b[2]) ||
// 		(p0[2] < a[2] && p1[2] < a[2] && p2[2] < a[2] && p0[2] < b[2] && p1[2] < b[2] && p2[2] < b[2]))
// 		return 0;
// 	double e = 1E-11; //これを1E-14とすることで，干渉のチェックが行われるようになる場合があった．

// 	auto d = factorOfVectorToReachTriangle(p0, p1, p2, a, b);
// 	if (d < 0. || d > 1.)
// 		return 0; /*面に到達できていない*/

// 	V_d b_a = b - a;
// 	V_d n = TriangleNormal(p0, p1, p2);
// 	V_d ps = a + (b_a)*d;

// 	//ポリゴン頂点の最大最小でチェック
// 	if (((p0[0] - e > ps[0] && p1[0] - e > ps[0] && p2[0] - e > ps[0]) || (p0[0] + e < ps[0] && p1[0] + e < ps[0] && p2[0] + e < ps[0])) ||
// 		((p0[1] - e > ps[1] && p1[1] - e > ps[1] && p2[1] - e > ps[1]) || (p0[1] + e < ps[1] && p1[1] + e < ps[1] && p2[1] + e < ps[1])) ||
// 		((p0[2] - e > ps[2] && p1[2] - e > ps[2] && p2[2] - e > ps[2]) || (p0[2] + e < ps[2] && p1[2] + e < ps[2] && p2[2] + e < ps[2])))
// 		return 1; /*面の最大最小範囲にすら入れていない*/

// 	auto ps_p0 = p0 - ps;
// 	auto ps_p1 = p1 - ps;
// 	auto ps_p2 = p2 - ps;

// 	ps_p0 = ps_p0 / Norm(ps_p0);
// 	ps_p1 = ps_p1 / Norm(ps_p1);
// 	ps_p2 = ps_p2 / Norm(ps_p2);

// 	if (Dot(Cross(ps_p0, ps_p1), n) >= 0. &&
// 		Dot(Cross(ps_p1, ps_p2), n) >= 0. &&
// 		Dot(Cross(ps_p2, ps_p0), n) >= 0.)
// 		return 3; /*a,bは面と交差*/
// 	else
// 		return 2; /*a,bは面と交差していないが，かなり惜しい*/
// };

// /////////////////////////////////
// // V_d pOnSurface(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a, const V_d &b)
// // {
// //   V_d n = TriangleNormal(p0, p1, p2);
// //   return a + (b - a) * Dot(/*tangential vector*/ Norm(p0 - a) < Norm(p1 - a) ? (p1 - a) : (p0 - a),
// //                            /*normal vector*/ n) /
// //                  Dot((b - a), n);
// //   //return a + (b-a)*Dot(p0-a,n)/Dot(b-a,n);
// // }
// V_d pOnSurface(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a, const V_d &b)
// {
// 	V_d n = TriangleNormal(p0, p1, p2), b_a = b - a;
// 	return a + b_a * Dot(p0 - a, n) / Dot(b_a, n); //分母が0の場合はあり得る
// }

// V_d pOnSurface(const VV_d &p0p1p2, const VV_d &ab)
// {
// 	if (p0p1p2.size() != 3)
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point size = " + std::to_string(p0p1p2.size()));

// 	return pOnSurface(p0p1p2[0],
// 					  p0p1p2[1],
// 					  p0p1p2[2],
// 					  ab[0], ab[1]);
// }

// int isIntersectingSurface(const VV_d &p0p1p2, const VV_d &ab)
// {
// 	if (p0p1p2.size() != 3)
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point size = " + std::to_string(p0p1p2.size()));

// 	return isIntersectingSurface(p0p1p2[0],
// 								 p0p1p2[1],
// 								 p0p1p2[2],
// 								 ab[0], ab[1]);
// };
//===========================================================
#include "fundamental_statistics.hpp"
#include "fundamental_geometry.hpp"
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
using Tdd = std::tuple<double, double>;
using Tddd = std::tuple<double, double, double>;
using T4d = std::tuple<double, double, double, double>;
using T5d = std::tuple<double, double, double, double, double>;
using T6d = std::tuple<double, double, double, double, double, double>;
using T7d = std::tuple<double, double, double, double, double, double, double>;
using Tiii = std::tuple<int, int, int>;
using T2Tdd = std::tuple<Tdd, Tdd>;
using T3Tdd = std::tuple<Tdd, Tdd, Tdd>;
using T3Tddd = std::tuple<Tddd, Tddd, Tddd>;
using T4Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd>;
using T5Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd>;
using T6Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T7Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
using T8Tddd = std::tuple<Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd, Tddd>;
/* ------------------------------------------------------ */
double Total(const Tddd &v) { return std::get<0>(v) + std::get<1>(v) + std::get<2>(v); };
double Total(const Tdd &v) { return std::get<0>(v) + std::get<1>(v); };
Tddd Total(const T3Tddd &v)
{
	return {std::get<0>(std::get<0>(v)) + std::get<0>(std::get<1>(v)) + std::get<0>(std::get<2>(v)),
			std::get<1>(std::get<0>(v)) + std::get<1>(std::get<1>(v)) + std::get<1>(std::get<2>(v)),
			std::get<2>(std::get<0>(v)) + std::get<2>(std::get<1>(v)) + std::get<2>(std::get<2>(v))};
};
Tdd Total(const T3Tdd &v)
{
	return {std::get<0>(std::get<0>(v)) + std::get<0>(std::get<1>(v)) + std::get<0>(std::get<2>(v)),
			std::get<1>(std::get<0>(v)) + std::get<1>(std::get<1>(v)) + std::get<1>(std::get<2>(v))};
};
Tddd Total(const std::vector<Tddd> &V)
{
	Tddd ret = {0, 0, 0};
	for (const auto &v : V)
		ret += v;
	return ret;
};
T4d Total(const std::vector<T4d> &V)
{
	T4d ret = {0, 0, 0, 0};
	for (const auto &v : V)
		ret += v;
	return ret;
};
T6d Total(const std::vector<T6d> &V)
{
	T6d ret = {0, 0, 0, 0, 0, 0};
	for (const auto &v : V)
		ret += v;
	return ret;
};
/* ------------------------------------------------------ */
// b% ------------------- タプルからparticlize ------------------ */
std::vector<std::tuple<Tddd /*実際の座標*/, Tdd /*パラメタt0t1*/>>
triangleIntoPoints(const T3Tddd &X0X1X2, const double dx)
{
	/*
	|-o--o--o--o-| n = 4
	x,yはパラメタ
	*/
	auto intp = [&X0X1X2](double t0, double t1)
	{
		return Dot(Tddd{t0, t1, 1 - t0 - t1}, X0X1X2);
	};
	auto y_list = [&intp, &dx](const double x)
	{
		int n = std::round(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
		// int n = std::ceil(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
		V_d ret(n);
		for (auto i = 0; i < n; ++i)
			ret[i] = ((1. - x) / n) * (i + 0.5);
		return ret;
	};
	/* ------------------------------------------------------ */
	Tddd v = intp(0., 1.) - intp(0., 0.);
	Tddd u = intp(1., 0.) - intp(0., 0.);
	double height = Norm(u - Dot(Normalize(v), u) * Normalize(v));

	int n = std::round(height / dx);
	// int n = std::ceil(height / dx);
	// std::cout << Grid({"n", n, "height", height}) << std::endl;
	double dH = 1. / n;
	V_d x_list(n);
	for (auto i = 0; i < n; i++)
		x_list[i] = dH * (i + 0.5);
	/* ------------------------------------------------------ */
	std::vector<std::tuple<Tddd /*実際の座標*/, Tdd /*パラメタt0t1*/>> ret;
	ret.reserve(x_list.size() * x_list.size());
	for (const auto &x : x_list)
		for (const auto &y : y_list(x))
			ret.push_back({intp(x, y), Tdd{x, y}});
	return ret;
};

std::vector<std::tuple<Tddd, Tdd>> particlize(const T3Tddd &X0X1X2, const double dx)
{
	return triangleIntoPoints(X0X1X2, dx);
};
/* ------------------------------------------------------ */
/*         近傍のオブジェクト探査のための，空間分割バケツ         */
/* ------------------------------------------------------ */
//テンプレートどんなオブジェクトでもできるはず
// #define debug_BaseBuckets
template <typename T>
struct BaseBuckets
{
	//! getX()でxyz座標を取得できるオブジェクトTのための，バケツ
	Tdd xbounds;
	Tdd ybounds;
	Tdd zbounds;
	int xsize, ysize, zsize;
	T3Tdd bounds;
	Tddd center;
	Tiii dn;
	// std::vector</*x*/ std::vector</*y*/ std::vector</*z*/ std::vector<T *>>>> buckets;
	std::vector</*x*/ std::vector</*y*/ std::vector</*z*/ std::unordered_set<T *>>>> buckets;
	std::unordered_set<T *> all_stored_objects;
	/*
	bucketsベクトルは，unordered_setで重複を認めないが，他のベクトルには，同じオブジェクトが入っていてもよい．
	例えば，鏡像関係を作りたい場合，他のバケット位置に，同じオブジェクトを入れておけばいい．
	*/
	double dL;
	double bucketVolume() const
	{
		return std::pow(this->dL, 3.);
	};
	double bucketVolume(const int depth) const
	{
		//    depth  1  2  3  4
		// edge len  1  3  5  7   = 2*d-1
		//  buckets  1  9  25  49 = (2*d-1)^2
		auto vol = std::pow(this->dL, 3.);
		return pow(2. * depth - 1., 2.) * vol;
	};
	BaseBuckets(const geometry::CoordinateBounds &c_bounds, const double dL_IN) : bounds(c_bounds.bounds), dL(dL_IN), dn({0, 0, 0})
	{
		set(this->bounds, dL_IN);
	};
	BaseBuckets(const T3Tdd &boundingboxIN, const double dL_IN) : bounds(boundingboxIN), dL(dL_IN), dn({0, 0, 0})
	{
		set(this->bounds, dL_IN);
	};
	void set(const T3Tdd &boundingboxIN, const double dL_IN)
	{
		try
		{
			this->bounds = boundingboxIN;
			this->dL = dL_IN;
			this->xbounds = std::get<0>(this->bounds);
			this->ybounds = std::get<1>(this->bounds);
			this->zbounds = std::get<2>(this->bounds);
			this->center = {Mean(this->xbounds), Mean(this->ybounds), Mean(this->zbounds)};
			this->xsize = std::ceil((std::get<1>(this->xbounds) - std::get<0>(this->xbounds)) / this->dL);
			this->ysize = std::ceil((std::get<1>(this->ybounds) - std::get<0>(this->ybounds)) / this->dL);
			this->zsize = std::ceil((std::get<1>(this->zbounds) - std::get<0>(this->zbounds)) / this->dL);
			this->dn = {xsize, ysize, zsize};
			this->buckets.clear();
			this->all_stored_objects.clear();
			this->buckets.resize(xsize,
								 std::vector</*y*/ std::vector</*z*/ std::unordered_set<T *>>>(ysize,
																							   std::vector</*z*/ std::unordered_set<T *>>(zsize, std::unordered_set<T *>{})));
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};
	void resize(const T3Tdd &boundingboxIN, const double dL_IN)
	{
		try
		{
			this->bounds = boundingboxIN;
			this->dL = dL_IN;
			this->xbounds = std::get<0>(this->bounds);
			this->ybounds = std::get<1>(this->bounds);
			this->zbounds = std::get<2>(this->bounds);
			this->xsize = std::ceil((std::get<1>(this->xbounds) - std::get<0>(this->xbounds)) / this->dL);
			this->ysize = std::ceil((std::get<1>(this->ybounds) - std::get<0>(this->ybounds)) / this->dL);
			this->zsize = std::ceil((std::get<1>(this->zbounds) - std::get<0>(this->zbounds)) / this->dL);
			this->dn = {xsize, ysize, zsize};

			std::cout << "１バケット幅 = " << this->dL << std::endl;
			std::cout << "バウンド = " << this->bounds << std::endl;
			std::cout << "サイズ = " << this->dn << std::endl;

			this->buckets.clear();
			this->all_stored_objects.clear();
			this->buckets.resize(xsize,
								 std::vector</*y*/ std::vector</*z*/ std::unordered_set<T *>>>(ysize,
																							   std::vector</*z*/ std::unordered_set<T *>>(zsize, std::unordered_set<T *>{})));
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		};
	};

	const std::unordered_set<T *> &getAll() const
	{
		return this->all_stored_objects;
	};

	// x座標を内包するバケツのインデックスを返す
	void indices(const Tddd &x, int &i, int &j, int &k) const
	{
		i = (int)((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL);
		j = (int)((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL);
		k = (int)((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL);
	};
	Tiii indices(const Tddd &x) const
	{
		// intのキャストはゼロ方向へ実数を切り捨てた結果を返す
		/*
		 0.**    1.**  2.**   3.**  <= X-minX というわけで，(int)によって正しくセルのインデックスに変換できる
		<-dL-> <-dL-> <-dL-> <-dL->
		*-----*------*------*------*
		|  0  |   1  |   2  |   3  |
		*-----*------*------*------*
		*/
		return {(int)((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL),
				(int)((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL),
				(int)((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL)};
	};
	Tddd indices2X(const int i, const int j, const int k) const
	{
		return {i * this->dL + std::get<0>(this->xbounds) /*min*/,
				j * this->dL + std::get<0>(this->ybounds) /*min*/,
				k * this->dL + std::get<0>(this->zbounds) /*min*/};
	};
	Tddd indices2X_center(const int i, const int j, const int k) const
	{
		return {this->dL * (i + 0.5) + std::get<0>(this->xbounds) /*min*/,
				this->dL * (j + 0.5) + std::get<0>(this->ybounds) /*min*/,
				this->dL * (k + 0.5) + std::get<0>(this->zbounds) /*min*/};
	};
	//インデックスがboundsに収まっているかどうか
	bool isInside(const int i, const int j, const int k) const
	{
		return (!(i < 0 || j < 0 || k < 0 || i >= this->xsize || j >= this->ysize || k >= this->zsize));
	};
	bool isInside(const Tiii &i) const
	{
		return (!(std::get<0>(i) < 0 || std::get<1>(i) < 0 || std::get<2>(i) < 0 || std::get<0>(i) >= this->xsize || std::get<1>(i) >= this->ysize || std::get<2>(i) >= this->zsize));
		// return isInside(std::get<0>(i), std::get<1>(i), std::get<2>(i));
	};
	void add(const Tddd &x, T *const p)
	{
		auto [i, j, k] = indices(x);
		try
		{
			if (!(i < 0 || j < 0 || k < 0 || i >= this->xsize || j >= this->ysize || k >= this->zsize))
			{
				this->buckets[i][j][k].emplace(p);
				this->all_stored_objects.emplace(p);
			}
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			std::stringstream ss;
			ss << "[i,j,k] = "
			   << "[" << i << "," << j << "," << k << "]" << std::endl;
			ss << "[xsize,ysize,zsize] = "
			   << "[" << xsize << "," << ysize << "," << zsize << "]" << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
		};
	};
	//座標を入力し，バケツを指定する．
	std::vector<T *> getAllBuckets() const
	{
		return this->buckets;
	};
	std::vector<T4Tddd> getT4Tddd() const
	{
		std::vector<T4Tddd> ret;
		for (auto k = 0; k < this->zsize; ++k)
			for (auto j = 0; j < this->ysize; ++j)
				for (auto i = 0; i < this->xsize; ++i)
				{
					Tddd X0 = {(i * this->dL + std::get<0>(this->xbounds)),
							   (j * this->dL + std::get<0>(this->ybounds)),
							   (k * this->dL + std::get<0>(this->zbounds))};
					Tddd X1 = {((i + 1) * this->dL + std::get<0>(this->xbounds)),
							   (j * this->dL + std::get<0>(this->ybounds)),
							   (k * this->dL + std::get<0>(this->zbounds))};
					Tddd X2 = {((i + 1) * this->dL + std::get<0>(this->xbounds)),
							   ((j + 1) * this->dL + std::get<0>(this->ybounds)),
							   (k * this->dL + std::get<0>(this->zbounds))};
					Tddd X3 = {(i * this->dL + std::get<0>(this->xbounds)),
							   ((j + 1) * this->dL + std::get<0>(this->ybounds)),
							   (k * this->dL + std::get<0>(this->zbounds))};
					Tddd X4 = {(i * this->dL + std::get<0>(this->xbounds)),
							   (j * this->dL + std::get<0>(this->ybounds)),
							   ((k + 1) * this->dL + std::get<0>(this->zbounds))};
					Tddd X5 = {((i + 1) * this->dL + std::get<0>(this->xbounds)),
							   (j * this->dL + std::get<0>(this->ybounds)),
							   ((k + 1) * this->dL + std::get<0>(this->zbounds))};
					Tddd X6 = {((i + 1) * this->dL + std::get<0>(this->xbounds)),
							   ((j + 1) * this->dL + std::get<0>(this->ybounds)),
							   ((k + 1) * this->dL + std::get<0>(this->zbounds))};
					Tddd X7 = {(i * this->dL + std::get<0>(this->xbounds)),
							   ((j + 1) * this->dL + std::get<0>(this->ybounds)),
							   ((k + 1) * this->dL + std::get<0>(this->zbounds))};
					ret.push_back({X0, X1, X2, X3});
					ret.push_back({X0, X1, X5, X4});
					ret.push_back({X1, X2, X6, X5});
					ret.push_back({X0, X4, X7, X3});
					ret.push_back({X2, X3, X7, X6});
					ret.push_back({X4, X5, X6, X7});
				}
		return ret;
	};
	std::vector<T *> getBucket(const Tddd &x) const
	{
		Tiii i = this->indices(x);
		return this->buckets[std::get<0>(i)][std::get<1>(i)][std::get<2>(i)];
	};
	/* ------------------------------------------------------ */
	void clear()
	{
		this->buckets.clear();
		this->all_stored_objects.clear();
		// for (auto &vvu_b : this->buckets)
		// 	for (auto &vu_b : vvu_b)
		// 		for (auto &u_b : vu_b)
		// 			u_b.clear();
	};
	//@ 2021年11月9日さらに早くするために，
	//@ std::vector<std::vector<T *>> getObjects(const Tddd &x, const int limit_depth /*limit depth*/, const int limit_number = 100000) const
	//@ を改良
	/* -------------- デフォルトのバケツは，深さ毎に粒子を保存していく -------------- */

	int min_index_x(const int i, const int d) const { return (i - d < 0 ? 0 : i - d); };
	int max_index_x(const int i, const int d) const { return (i + d > this->xsize - 1 ? this->xsize - 1 : i + d); };
	int min_index_y(const int i, const int d) const { return (i - d < 0 ? 0 : i - d); };
	int max_index_y(const int i, const int d) const { return (i + d > this->ysize - 1 ? this->ysize - 1 : i + d); };
	int min_index_z(const int i, const int d) const { return (i - d < 0 ? 0 : i - d); };
	int max_index_z(const int i, const int d) const { return (i + d > this->zsize - 1 ? this->zsize - 1 : i + d); };

	std::unordered_set<T *> getObjects_unorderedset(const Tddd &x,
													const int limit_depth /*limit depth*/,
													const int limit_number = 100000) const
	{
		std::unordered_set<T *> ret;
		// ret.reserve(limit_depth);
		// std::vector<T *> ret_last(0);
		int i0, j0, k0;
		this->indices(x, i0, j0, k0);
		//* ------------------------ depth=0 ------------------------ */
		// auto &r = *ret.rbegin();
		if (isInside(i0, j0, k0))
			ret.insert(this->buckets[i0][j0][k0].begin(), this->buckets[i0][j0][k0].end());
		if (ret.size() >= limit_number)
			return ret;
		/* ------------------------------------------------------ */
		int i, j, k, j_min, j_max, k_min, k_max, i_min_1, i_max_1, j_min_1, j_max_1, i0_m_d, i0_p_d, j0_m_d, j0_p_d, k0_m_d, k0_p_d;
		for (auto d = 1; d <= limit_depth; ++d)
		{
			j_min = min_index_y(j0, d);
			j_max = max_index_y(j0, d);
			k_min = min_index_z(k0, d);
			k_max = max_index_z(k0, d);
			j_min_1 = min_index_y(j0, d - 1);
			j_max_1 = max_index_y(j0, d - 1);
			i0_m_d = i0 - d;
			i0_p_d = i0 + d;
			for (auto j = j_min; j <= j_max /*need equal (=)*/; ++j)
				for (auto k = k_min; k <= k_max; ++k)
				{
					if (i0_m_d >= 0 && i0_m_d < this->xsize)
						ret.insert(std::begin(this->buckets[i0_m_d][j][k]), std::end(this->buckets[i0_m_d][j][k]));
					if (i0_p_d >= 0 && i0_p_d < this->xsize)
						ret.insert(std::begin(this->buckets[i0_p_d][j][k]), std::end(this->buckets[i0_p_d][j][k]));
				}
			i_min_1 = min_index_x(i0, d - 1);
			i_max_1 = max_index_x(i0, d - 1);
			//
			j0_m_d = j0 - d;
			j0_p_d = j0 + d;
			//
			k0_m_d = k0 - d;
			k0_p_d = k0 + d;
			//
			for (auto i = i_min_1; i <= i_max_1; ++i)
			{
				for (auto k = k_min; k <= k_max; ++k)
				{
					if (j0_m_d >= 0 && j0_m_d < this->ysize)
						ret.insert(std::begin(this->buckets[i][j0_m_d][k]), std::end(this->buckets[i][j0_m_d][k]));
					if (j0_p_d >= 0 && j0_p_d < this->ysize)
						ret.insert(std::begin(this->buckets[i][j0_p_d][k]), std::end(this->buckets[i][j0_p_d][k]));
				}
				for (auto j = j_min_1; j <= j_max_1; ++j)
				{
					if (k0_m_d >= 0 && k0_m_d < this->zsize)
						ret.insert(std::begin(this->buckets[i][j][k0_m_d]), std::end(this->buckets[i][j][k0_m_d]));
					if (k0_p_d >= 0 && k0_p_d < this->zsize)
						ret.insert(std::begin(this->buckets[i][j][k0_p_d]), std::end(this->buckets[i][j][k0_p_d]));
				}
			}
			if (ret.size() >= limit_number)
				return ret;
		}
		return ret;
	};
	std::unordered_set<T *> getObjects_unorderedset(const Tddd &x, const double smoothing_length, const int limit_number = 100000) const
	{
		int depth = std::ceil(smoothing_length / this->dL);
		return getObjects_unorderedset(x, depth, limit_number);
	};
	/* -------------- デフォルトのバケツは，深さ毎に粒子を保存していく -------------- */
	std::vector<std::vector<T *>> getObjects(const Tddd &x,
											 const int limit_depth /*limit depth*/,
											 const int limit_number = 100000) const
	{
		std::vector<std::vector<T *>> ret(1);
		// ret.reserve(limit_depth);
		// std::vector<T *> ret_last(0);
		int i0, j0, k0;
		this->indices(x, i0, j0, k0);
		//* ------------------------ depth=0 ------------------------ */
		auto &r = *ret.rbegin();
		if (isInside(i0, j0, k0))
			r.insert(r.end(), this->buckets[i0][j0][k0].begin(), this->buckets[i0][j0][k0].end());
		int tot = r.size();
		if (tot >= limit_number)
			return ret;
		/* ------------------------------------------------------ */
		int i, j, k;
		for (auto d = 1; d < limit_depth; d++)
		{
			ret.push_back({});
			auto &r = *ret.rbegin();
			for (auto j = (j0 - d < 0 ? 0 : j0 - d); j <= j0 + d && j < this->ysize; j++)
				for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; k++)
				{
					i = i0 + d;
					if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
					i = i0 - d;
					if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
				}
			for (auto i = (i0 - d + 1 < 0 ? 0 : i0 - d + 1); i < i0 + d && i < this->xsize; i++)
			{
				for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; k++)
				{
					j = j0 + d;
					if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
					j = j0 - d;
					if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
				}
				for (auto j = (j0 - d + 1 < 0 ? 0 : j0 - d + 1); j < j0 + d && j < this->ysize; j++)
				{
					k = k0 + d;
					if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
					k = k0 - d;
					if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
				}
			}
			tot += r.size();
			if (tot >= limit_number)
				return ret;
		}
		return ret;
	};
	std::vector<std::vector<T *>> getObjects(const Tddd &x, const double smoothing_length, const int limit_number = 100000) const
	{
		int depth = std::ceil(smoothing_length / this->dL);
		return getObjects(x, depth, limit_number);
	};
	/* ------------------------------------------------------ */
	std::vector<T *> getObjectsFlattened(const Tddd &x, const int d /*limit depth*/) const
	{
#if defined(debug_BaseBuckets)
		std::cout << "d = " << d << std::endl;
#endif
		std::vector<T *> ret;
		int i0, j0, k0;
		this->indices(x, i0, j0, k0);
#if defined(debug_BaseBuckets)
		std::cout << "[i0,j0,k0] = [" << i0 << "," << j0 << "," << k0 << "]" << std::endl;
#endif
		int i_beg = ((i0 - d) >= 0 ? (i0 - d) : 0);
		int j_beg = ((j0 - d) >= 0 ? (j0 - d) : 0);
		int k_beg = ((k0 - d) >= 0 ? (k0 - d) : 0);

		int i_end = ((i0 + d) <= this->xsize ? (i0 + d) : this->xsize);
		int j_end = ((j0 + d) <= this->ysize ? (j0 + d) : this->ysize);
		int k_end = ((k0 + d) <= this->zsize ? (k0 + d) : this->zsize);

		i_beg = (i_beg <= this->xsize ? i_beg : this->xsize);
		j_beg = (j_beg <= this->ysize ? j_beg : this->ysize);
		k_beg = (k_beg <= this->zsize ? k_beg : this->zsize);

		i_end = (i_end >= 0 ? i_end : 0);
		j_end = (j_end >= 0 ? j_end : 0);
		k_end = (k_end >= 0 ? k_end : 0);
#if defined(debug_BaseBuckets)
		std::cout << "[i_beg,j_beg,k_beg] = "
				  << "[" << i_beg << "," << j_beg << "," << k_beg << "]" << std::endl;
		std::cout << "[i_end,j_end,k_end] = "
				  << "[" << i_end << "," << j_end << "," << k_end << "]" << std::endl;
		std::cout << "[xsize,ysize,zsize] = "
				  << "[" << xsize << "," << ysize << "," << zsize << "]" << std::endl;
		try
		{
#endif
			for (auto it = this->buckets.begin() + i_beg; it != this->buckets.begin() + i_end; ++it)
			{
				for (auto jt = it->begin() + j_beg; jt != it->begin() + j_end; ++jt)
					for (auto kt = jt->begin() + k_beg; kt != jt->begin() + k_end; ++kt)
						ret.insert(ret.end(), kt->begin(), kt->end());
			}
			return ret;
#if defined(debug_BaseBuckets)
		}
		catch (std::exception &e)
		{
			std::cerr << e.what() << reset << std::endl;
			std::stringstream ss;
			ss << "[i_beg,j_beg,k_beg] = "
			   << "[" << i_beg << "," << j_beg << "," << k_beg << "]" << std::endl;
			ss << "[i_end,j_end,k_end] = "
			   << "[" << i_end << "," << j_end << "," << k_end << "]" << std::endl;
			ss << "[xsize,ysize,zsize] = "
			   << "[" << xsize << "," << ysize << "," << zsize << "]" << std::endl;
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
		};
#endif
	};
	std::vector<T *> getObjectsFlattened(const Tddd &x, const double smoothing_length) const
	{
		int depth = std::round(std::ceil(smoothing_length / this->dL));
		return getObjectsFlattened(x, depth);
	};
	/* ------------------------------------------------------ */
	void apply(const std::function<void(T *)> &fun, const Tddd &x, const int d /*limit depth*/, const int limit_number = 100000)
	{
		int i0, j0, k0;
		this->indices(x, i0, j0, k0);
		int i_beg = ((i0 - d) >= 0 ? (i0 - d) : 0);
		int j_beg = ((j0 - d) >= 0 ? (j0 - d) : 0);
		int k_beg = ((k0 - d) >= 0 ? (k0 - d) : 0);
		int i_end = ((i0 + d) <= this->xsize ? (i0 + d) : this->xsize);
		int j_end = ((j0 + d) <= this->ysize ? (j0 + d) : this->ysize);
		int k_end = ((k0 + d) <= this->zsize ? (k0 + d) : this->zsize);
		for (auto it = this->buckets.begin() + i_beg; it != this->buckets.begin() + i_end; ++it)
			for (auto jt = it->begin() + j_beg; jt != it->begin() + j_end; ++jt)
				for (auto kt = jt->begin() + k_beg; kt != jt->begin() + k_end; ++kt)
					for (auto objT = kt->begin(); objT != kt->end(); ++objT)
						fun(*objT);
	};
};

template <typename T>
struct Buckets : public BaseBuckets<T>
{
	Buckets(const geometry::CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<T>(c_bounds, dL_IN){};
	Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<T>(boundingboxIN, dL_IN){};
};

/* ------------------------------------------------------ */
// template <typename T>
// class Buckets
// {
// public:
// 	//!getX()でxyz座標を取得できるオブジェクトTのための，バケツ
// 	Tdd xbounds;
// 	Tdd ybounds;
// 	Tdd zbounds;
// 	int xsize, ysize, zsize;
// 	T3Tdd bounds;
// 	Tiii dn;
// 	std::vector</*x*/ std::vector</*y*/ std::vector</*z*/ std::vector<T *>>>> buckets;
// 	double dL;

// 	Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : bounds(boundingboxIN), dL(dL_IN), dn({0, 0, 0})
// 	{
// 		set(boundingboxIN, dL_IN);
// 	};
// 	void set(const T3Tdd &boundingboxIN, const double dL_IN)
// 	{
// 		this->bounds = boundingboxIN;
// 		this->dL = dL_IN;
// 		this->xbounds = std::get<0>(this->bounds);
// 		this->ybounds = std::get<1>(this->bounds);
// 		this->zbounds = std::get<2>(this->bounds);
// 		this->xsize = (int)((std::get<1>(this->xbounds) - std::get<0>(this->xbounds)) / this->dL);
// 		this->ysize = (int)((std::get<1>(this->ybounds) - std::get<0>(this->ybounds)) / this->dL);
// 		this->zsize = (int)((std::get<1>(this->zbounds) - std::get<0>(this->zbounds)) / this->dL);
// 		this->dn = {xsize, ysize, zsize};
// 		this->buckets.resize(xsize,
// 							 std::vector</*y*/ std::vector</*z*/ std::vector<T *>>>(ysize,
// 																					std::vector</*z*/ std::vector<T *>>(zsize, std::vector<T *>(0))));
// 	};
// 	//x座標を内包するバケツのインデックスを返す
// 	void indices(const Tddd &x, int &i, int &j, int &k) const
// 	{
// 		i = (int)((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL);
// 		j = (int)((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL);
// 		k = (int)((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL);
// 	};
// 	Tiii indices(const Tddd &x) const
// 	{
// 		return {(int)((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL),
// 				(int)((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL),
// 				(int)((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL)};
// 	};
// 	Tddd indices2X(const int i, const int j, const int k) const
// 	{
// 		return {i * this->dL + std::get<0>(this->xbounds) /*min*/,
// 				j * this->dL + std::get<0>(this->ybounds) /*min*/,
// 				k * this->dL + std::get<0>(this->zbounds) /*min*/};
// 	};
// 	//インデックスがboundsに収まっているかどうか
// 	bool isInside(const int i, const int j, const int k) const
// 	{
// 		return (!(i < 0 || j < 0 || k < 0 || i >= this->xsize || j >= this->ysize || k >= this->zsize));
// 	};
// 	bool isInside(const Tiii &i) const
// 	{
// 		return isInside(std::get<0>(i), std::get<1>(i), std::get<2>(i));
// 	};
// 	//
// 	// void add(const Tddd &x, T *p)
// 	// {
// 	// 	try
// 	// 	{
// 	// 		int i, j, k;
// 	// 		this->indices(x, i, j, k);
// 	// 		if (isInside(i, j, k))
// 	// 			this->buckets[i][j][k].emplace_back(p);
// 	// 	}
// 	// 	catch (error_message &e)
// 	// 	{
// 	// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	// 	}
// 	// };
// 	void add(const Tddd &x, T *p)
// 	{
// 		try
// 		{
// 			auto [i, j, k] = this->indices(x);
// 			if (isInside(i, j, k) && !std::binary_search(this->buckets[i][j][k].begin(), this->buckets[i][j][k].end(), p))
// 				this->buckets[i][j][k].insert(std::lower_bound(this->buckets[i][j][k].begin(), this->buckets[i][j][k].end(), p), p);
// 		}
// 		catch (error_message &e)
// 		{
// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 		}
// 	};
// 	//座標を入力し，バケツを指定する．
// 	std::vector<T *> getAllBuckets() const { return this->buckets; };
// 	std::vector<T *> getBucket(const Tddd &x) const
// 	{
// 		Tiii i = this->indices(x);
// 		return this->buckets[std::get<0>(i)][std::get<1>(i)][std::get<2>(i)];
// 	};
// 	/* -------------- デフォルトのバケツは，深さ毎に粒子を保存していく -------------- */
// 	std::vector<std::vector<T *>> getObjects(const Tddd &x, const int depth /*limit depth*/, const int limit_number = 100000) const
// 	{
// 		std::vector<std::vector<T *>> ret(depth);
// 		// std::vector<T *> ret_last(0);
// 		int i0, j0, k0;
// 		this->indices(x, i0, j0, k0);
// 		//* ------------------------ depth=0 ------------------------ */
// 		auto &r = ret[0];
// 		if (isInside(i0, j0, k0))
// 			r.insert(r.end(), this->buckets[i0][j0][k0].begin(), this->buckets[i0][j0][k0].end());
// 		int tot = r.size();
// 		int i, j, k;
// 		for (auto d = 1; d < depth; d++)
// 		{
// 			r = ret[d];
// 			for (auto j = (j0 - d < 0 ? 0 : j0 - d); j <= j0 + d && j < this->ysize; j++)
// 				for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; k++)
// 				{
// 					i = i0 + d;
// 					if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 					i = i0 - d;
// 					if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 				}
// 			for (auto i = (i0 - d + 1 < 0 ? 0 : i0 - d + 1); i < i0 + d && i < this->xsize; i++)
// 			{
// 				for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; k++)
// 				{
// 					j = j0 + d;
// 					if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 					j = j0 - d;
// 					if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 				}
// 				for (auto j = (j0 - d + 1 < 0 ? 0 : j0 - d + 1); j < j0 + d && j < this->ysize; j++)
// 				{
// 					k = k0 + d;
// 					if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 					k = k0 - d;
// 					if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
// 						r.insert(r.begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
// 				}
// 			}
// 			tot += r.size();
// 			if (tot >= limit_number)
// 				return ret;
// 		}
// 		return ret;
// 	};
// 	std::vector<std::vector<T *>> getObjects(const Tddd &x, const double smoothing_length, const int limit_number = 100000) const
// 	{
// 		int depth = std::ceil(smoothing_length / this->dL);
// 		return getObjects(x, depth, limit_number);
// 	};
// 	/* ------------------------------------------------------ */
// 	std::vector<T *> getObjectsFlattened(const Tddd &x, const int d /*limit depth*/, const int limit_number = 100000) const
// 	{
// 		std::vector<T *> ret;
// 		int i0, j0, k0;
// 		this->indices(x, i0, j0, k0);
// 		int i_beg = ((i0 - d) >= 0 ? (i0 - d) : 0);
// 		int j_beg = ((j0 - d) >= 0 ? (j0 - d) : 0);
// 		int k_beg = ((k0 - d) >= 0 ? (k0 - d) : 0);
// 		int i_end = ((i0 + d) <= this->xsize ? (i0 + d) : this->xsize);
// 		int j_end = ((j0 + d) <= this->ysize ? (j0 + d) : this->ysize);
// 		int k_end = ((k0 + d) <= this->zsize ? (k0 + d) : this->zsize);
// 		for (auto it = this->buckets.begin() + i_beg; it != this->buckets.begin() + i_end; ++it)
// 			for (auto jt = it->begin() + j_beg; jt != it->begin() + j_end; ++jt)
// 				for (auto kt = jt->begin() + k_beg; kt != jt->begin() + k_end; ++kt)
// 					ret.insert(ret.end(), kt->begin(), kt->end());
// 		return ret;
// 	};
// 	std::vector<T *> getObjectsFlattened(const Tddd &x, const double smoothing_length, const int limit_number = 100000) const
// 	{
// 		int depth = std::ceil(smoothing_length / this->dL);
// 		return getObjectsFlattened(x, depth, limit_number);
// 	};
// 	/* ------------------------------------------------------ */
// 	void apply(const std::function<void(T *)> &fun, const Tddd &x, const int d /*limit depth*/, const int limit_number = 100000)
// 	{
// 		int i0, j0, k0;
// 		this->indices(x, i0, j0, k0);
// 		int i_beg = ((i0 - d) >= 0 ? (i0 - d) : 0);
// 		int j_beg = ((j0 - d) >= 0 ? (j0 - d) : 0);
// 		int k_beg = ((k0 - d) >= 0 ? (k0 - d) : 0);
// 		int i_end = ((i0 + d) <= this->xsize ? (i0 + d) : this->xsize);
// 		int j_end = ((j0 + d) <= this->ysize ? (j0 + d) : this->ysize);
// 		int k_end = ((k0 + d) <= this->zsize ? (k0 + d) : this->zsize);
// 		for (auto it = this->buckets.begin() + i_beg; it != this->buckets.begin() + i_end; ++it)
// 			for (auto jt = it->begin() + j_beg; jt != it->begin() + j_end; ++jt)
// 				for (auto kt = jt->begin() + k_beg; kt != jt->begin() + k_end; ++kt)
// 					for (auto objT = kt->begin(); objT != kt->end(); ++objT)
// 						fun(*objT);
// 	};
// };

//===========================================================
struct LoadObj
{
	V_d eachMax;
	V_d eachMin;

	VV_d v_complex;
	std::vector<Tddd> v_complex_tuple;
	std::vector<std::vector<float>> s_complex;
	VV_i f_v_complex; //これのstd::vector<int>が１面 -> 3頂点を指定

	std::vector<std::vector<float>> s;
	VV_d v, vn;

	VV_i f_v; // start from 0
	VV_i f_t;
	VV_i f_vn;

	std::vector<int> checkAllIntersection(const V_d &a, const V_d &b)
	{
		//物体全体の最大最小範囲の外で偏っている場合をチェック
		for (auto k = 0; k < 3; k++)
			if ((eachMax[k] < a[k] && eachMax[k] < b[k]) || (eachMin[k] > a[k] && eachMin[k] > b[k]))
				return {};
		std::vector<int> check0123(f_v_complex.size(), 0);
		int i = 0;
		for (const auto &f : f_v_complex)
		{
			//上以外
			auto mean = v_complex_tuple[f[0]] + v_complex_tuple[f[1]] + v_complex_tuple[f[2]];
			mean /= 3.;
			check0123[i] = isIntersectingSurface(v_complex_tuple[f[0]] + (v_complex_tuple[f[0]] - mean) * 1E-10,
												 v_complex_tuple[f[1]] + (v_complex_tuple[f[1]] - mean) * 1E-10,
												 v_complex_tuple[f[2]] + (v_complex_tuple[f[2]] - mean) * 1E-10,
												 Tddd{a[0], a[1], a[2]},
												 Tddd{b[0], b[1], b[2]});

			if (check0123[i] == 3)
			{
				s_complex[f[0]] = {0., 1., 0., 1.};
				s_complex[f[1]] = {0., 1., 0., 1.};
				s_complex[f[2]] = {0., 1., 0., 1.};
				return check0123;
			}
			// if(check0123[i] == 2){
			// 	s_complex[f[0]] += {0.,1.,0.,1.};
			// 	s_complex[f[1]] += {0.,1.,0.,1.};
			// 	s_complex[f[2]] += {0.,1.,0.,1.};
			// }

			i++;
		}
		return check0123;
	}

	void calculateMaxMin()
	{
		eachMax.resize(3);
		eachMin.resize(3);
		for (const auto &V : v_complex)
			for (auto i = 0; i < 3; i++)
			{
				if (eachMax[i] < V[i])
					eachMax[i] = V[i];
				else if (eachMin[i] > V[i])
					eachMin[i] = V[i];
			}
	};

	std::vector<int> checkAllIntersection(const VV_d &ab) { return checkAllIntersection(ab[0], ab[1]); }

	V_d surface(int i, double t0, double t1)
	{
		std::vector<int> indexes = f_v[i];

		VV_d a(indexes.size());
		for (auto n = 0; n < indexes.size(); n++)
			a[n] = v[indexes[n] - 1];

		return PointsToSurface(t0, t1, a);
	};

	LoadObj() : eachMax(3, 0.), eachMin(3, 0.){};
	LoadObj(const std::string &filename) : eachMax(3, 0.), eachMin(3, 0.)
	{
		std::vector<V_s> read_line;
		Load(filename, read_line, {"    ", "   ", "  ", " "});
		this->load(read_line);
		generateComplex();
		calculateMaxMin();
	};

	void generateComplex()
	{
		int j = 0;
		std::vector<int> tri_int{0, 0, 0};
		if (s.empty())
		{
			v_complex.clear();
			v_complex_tuple.clear();
			f_v_complex.clear();
			for (const auto &Ind : f_v)
			{
				for (const auto &i : Ind)
				{
					v_complex.emplace_back(v[i]);
					v_complex_tuple.emplace_back(Tddd{v[i][0], v[i][0], v[i][0]});
				}
				tri_int = {j, j + 1, j + 2};
				f_v_complex.emplace_back(tri_int);
				j = j + 3;
			}
		}
		else
		{
			v_complex.clear();
			v_complex_tuple.clear();
			s_complex.clear();
			f_v_complex.clear();
			for (const auto &Ind : f_v)
			{
				//	std::vector<float> tmp = {std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),1.};
				std::vector<float> tmp = {.4, .4, .4, 1.};
				for (const auto &i : Ind)
				{
					v_complex.emplace_back(v[i]);
					v_complex_tuple.emplace_back(Tddd{v[i][0], v[i][0], v[i][0]});
					s_complex.emplace_back(tmp);
				}
				tri_int = {j, j + 1, j + 2};
				f_v_complex.emplace_back(tri_int);
				j = j + 3;
			}
		}
		// }

		// s_complex.resize(f_v.size(),std::vector<float>(4,0));
		// int j=0;
		// for(const auto& Ind:f_v)
		//   for(const auto& i:Ind){
		//     s_complex[j] = v[i];
		// 	j++;
		//   }
	};

	// stringを変換しv, f_nに格納
	void load(VV_s &read_line)
	{
		Print(__PRETTY_FUNCTION__, Red);
		for (auto &line : read_line)
			if (line[0] == "v")
			{
				line.erase(line.begin());
				v.emplace_back(stod(line));
			}
			else if (line[0] == "f")
			{
				line.erase(line.begin());
				std::vector<int> f_v_tmp, f_t_tmp, f_vn_tmp;
				for (auto &l : line)
				{
					V_s v_t_vn = StringSplit(l, {"/"});
					if (v_t_vn.size() == 1)
						f_v_tmp.emplace_back(stoi(v_t_vn[0]) - 1);
					else if (v_t_vn.size() == 2)
					{
						f_v_tmp.emplace_back(stoi(v_t_vn[0]) - 1);
						if (!v_t_vn[1].empty())
							f_t_tmp.emplace_back(stoi(v_t_vn[1]) - 1);
					}
					else if (v_t_vn.size() == 3)
					{
						f_v_tmp.emplace_back(stoi(v_t_vn[0]) - 1);
						if (!v_t_vn[1].empty())
							f_t_tmp.emplace_back(stoi((v_t_vn[1])) - 1);
						if (!v_t_vn[2].empty())
							f_vn_tmp.emplace_back(stoi(v_t_vn[2]) - 1);
					}
				}
				if (f_v_tmp.size() > 0)
					f_v.emplace_back(f_v_tmp);
				if (f_t_tmp.size() > 0)
					f_t.emplace_back(f_t_tmp);
				if (f_vn_tmp.size() > 0)
					f_vn.emplace_back(f_vn_tmp);
			}
			else if (line[0] == "vn")
			{
				line.erase(line.begin());
				vn.emplace_back(stod(line));
			}
		s.resize(v.size(), std::vector<float>(4, 0)); // RGBS
		for (auto &tmp : s)
			tmp = {.9, .9, .9, 1.};

		// for(auto& tmp:s)
		//   tmp = {std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),1.};
	};

	void ConstructFromString(const std::string &strIN)
	{
		std::stringstream strm = std::stringstream(strIN);
		std::vector<V_s> read_line;
		std::string read;
		while (std::getline(strm, read)) // loop until the end of file
		{
			if (!StringTrim(read, {" "}).empty())
			{
				read_line.emplace_back(StringTrim(StringSplit(read, {"    ", "   ", "  ", " "}), {"    ", "   ", "  ", " ", "\t", "\n", "\r"}));
			}
		}
		load(read_line);
		generateComplex();
		calculateMaxMin();
		glLINES lines;
		for (auto i = 0; i < lines.f_v_complex.size(); i++)
		{
			checkAllIntersection({lines.v_complex[lines.f_v_complex[i][0]],
								  lines.v_complex[lines.f_v_complex[i][1]]});
		}
	};

	// std::string JSON(){
	//   std::string vertices = "{\n\"vertices\": [\n";
	//   for(auto i = 0; i<v.size(); i++){
	//     for(auto j = 0; j<v[i].size(); j++)
	//   	vertices += std::to_string(float(v[i][j])) + ", ";
	//     if(i != v.size()-1)
	// 	vertices += "\n";
	//     if(i == v.size()-1)
	//       vertices.pop_back();
	//   }
	//   vertices.pop_back();
	//   vertices += "],\n";
	//   std::string indices = "\"indices\": [\n";
	//   for(auto i = 0; i<f_v.size(); i++){
	//     for(auto j = 0; j<f_v[i].size(); j++)
	//       indices += std::to_string(f_v[i][j]-1) + ", ";
	//     if(i != f_v.size()-1)
	//       indices += "\n";
	//     if(i == f_v.size()-1)
	//       indices.pop_back();
	//   }
	//   indices.pop_back();
	//   indices += "]\n}";
	//   return vertices + indices;
	// };

	////////// 自動的に三角形にする
	// objindex をなおそう
	std::string JSON()
	{
		///////////////// vertices
		std::string str = "{\n\"vertices\": [\n";
		for (size_t i = 0; i < v.size(); i++)
		{
			for (size_t j = 0; j < v[i].size(); j++)
				str += std::to_string(float(v[i][j])) + ", ";
			if (i != v.size() - 1)
				str += "\n";
			if (i == v.size() - 1)
				str.pop_back();
		}
		str.pop_back();
		str += "],\n";
		/////////////////// scalars
		str += "\"scalars\": [\n";
		for (size_t i = 0; i < s.size(); i++)
		{
			for (size_t j = 0; j < s[i].size(); j++)
				str += std::to_string(float(s[i][j])) + ", ";
			if (i != s.size() - 1)
				str += "\n";
			if (i == s.size() - 1)
				str.pop_back();
		}
		str.pop_back();
		str += "],\n";
		///////////////// vertices_complex
		str += "\n\"vertices_complex\": [\n";
		for (size_t i = 0; i < v_complex.size(); i++)
		{
			for (size_t j = 0; j < v_complex[i].size(); j++)
				str += std::to_string(float(v_complex[i][j])) + ", ";
			if (i != v_complex.size() - 1)
				str += "\n";
			if (i == v_complex.size() - 1)
				str.pop_back();
		}
		str.pop_back();
		str += "],\n";
		/////////////////// scalars_complex
		str += "\"scalars_complex\": [\n";
		for (size_t i = 0; i < s_complex.size(); i++)
		{
			for (size_t j = 0; j < s_complex[i].size(); j++)
				str += std::to_string(float(s_complex[i][j])) + ", ";
			if (i != s_complex.size() - 1)
				str += "\n";
			if (i == s_complex.size() - 1)
				str.pop_back();
		}
		str.pop_back();
		str += "],\n";
		////////////////// index_complex
		str += "\"indices_complex\": [\n";
		for (size_t i = 0; i < f_v_complex.size(); i++)
		{
			if (f_v_complex[i].size() == 3)
			{
				for (size_t j = 0; j < f_v_complex[i].size(); j++)
					str += std::to_string(f_v_complex[i][j]) + ", ";
				if (i != f_v_complex.size() - 1)
					str += "\n";
				if (i == f_v_complex.size() - 1)
					str.pop_back();
			}
			else if (f_v_complex[i].size() == 4)
			{
				str += std::to_string(f_v_complex[i][0]) + ", " + std::to_string(f_v_complex[i][1]) + ", " + std::to_string(f_v_complex[i][2]) + ", ";
				str += std::to_string(f_v_complex[i][2]) + ", " + std::to_string(f_v_complex[i][3]) + ", " + std::to_string(f_v_complex[i][0]) + ", ";
				;
				if (i != f_v_complex.size() - 1)
					str += "\n";
				if (i == f_v_complex.size() - 1)
					str.pop_back();
			}
		}
		str.pop_back();
		str += "],\n";
		////////////////// indices
		str += "\"indices\": [\n";
		for (size_t i = 0; i < f_v.size(); i++)
		{
			if (f_v[i].size() == 3)
			{
				for (size_t j = 0; j < f_v[i].size(); j++)
					str += std::to_string(f_v[i][j]) + ", ";
				if (i != f_v.size() - 1)
					str += "\n";
				if (i == f_v.size() - 1)
					str.pop_back();
			}
			else if (f_v[i].size() == 4)
			{
				str += std::to_string(f_v[i][0]) + ", " + std::to_string(f_v[i][1]) + ", " + std::to_string(f_v[i][2]) + ", ";
				str += std::to_string(f_v[i][2]) + ", " + std::to_string(f_v[i][3]) + ", " + std::to_string(f_v[i][0]) + ", ";
				;
				if (i != f_v.size() - 1)
					str += "\n";
				if (i == f_v.size() - 1)
					str.pop_back();
			}
		}
		str.pop_back();
		str += "]\n}";
		////////////////////
		return str;
	};
	std::string JSON(const std::string &strIN)
	{
		ConstructFromString(strIN);
		return JSON();
	};
};
std::string obj2json(const std::string &filename)
{
	LoadObj loaded(filename);
	return loaded.JSON();
};
//============================================================
// class MyArray{
// public:
//   V_d vec;

//   MyArray(){
//   }

//   V_d getvec(){
//     return vec;
//   }

//   void setDoubleVec(double& address, int len){
//     vec.resize(len);
//     for(auto i=0; i<len; i++)
//       vec[i] = *(address + i);
//   }

//   V_d show(){
//     return vec;
//   }
// };
//======================================
// VV_d RotationMatrix(const double theta, const V_d &V)
// {
// 	// Euler Rodrigues
// 	double c = cos(theta), s = sin(theta), e = (1. - c);
// 	return {{V[0] * V[0] * e + c, V[0] * V[1] * e - V[2] * s, V[0] * V[2] * e + V[1] * s},
// 			{V[0] * V[1] * e + V[2] * s, V[1] * V[1] * e + c, V[1] * V[2] * e - V[0] * s},
// 			{V[0] * V[2] * e - V[1] * s, V[2] * V[1] * e + V[0] * s, V[2] * V[2] * e + c}};
// };

T3Tddd RotationMatrix(const double theta, const Tddd &V)
{
	// Euler Rodrigues
	double c = cos(theta), s = sin(theta), e = (1. - c);
	return {{std::get<0>(V) * std::get<0>(V) * e + c, std::get<0>(V) * std::get<1>(V) * e - std::get<2>(V) * s, std::get<0>(V) * std::get<2>(V) * e + std::get<1>(V) * s},
			{std::get<0>(V) * std::get<1>(V) * e + std::get<2>(V) * s, std::get<1>(V) * std::get<1>(V) * e + c, std::get<1>(V) * std::get<2>(V) * e - std::get<0>(V) * s},
			{std::get<0>(V) * std::get<2>(V) * e - std::get<1>(V) * s, std::get<2>(V) * std::get<1>(V) * e + std::get<0>(V) * s, std::get<2>(V) * std::get<2>(V) * e + c}};
};

template <class T>
std::vector<T> BarycentriCoordinate(const std::vector<std::vector<T>> &a, const T t0, const T t1)
{
	switch (a.size())
	{
	case 3:
		return Dot({t0, 1. - t0 - t1, t1}, a);
	case 4:
		return Dot({t0, 1. - t0, t1, 1. - t1}, a);
	}
};
//=====================================

// V_d TriangleArea(const VV_d& abcd){
//   V_d ret(abcd.size());
//   ret[0] = TriangleArea(abcd[abcd.size()-1], abcd[0], abcd[1]);
//   for(size_t i=0; i<abcd.size()-2; i++)
//     ret[i+1] = TriangleArea(abcd[i],abcd[i+1],abcd[i+2]);
//   ret[abcd.size()-1] = TriangleArea(abcd[abcd.size()-2],abcd[abcd.size()-1],abcd[0]);
//   return ret;
// };

std::vector<int> removePositiveMinArea(const V_d &area,
									   VV_d &vertices,
									   std::vector<int> &indices)
{
	V_d min_area;
	std::vector<int> min_ind, k_at_min;
	for (size_t k = 0; k < indices.size(); k++)
	{
		if (area[k] > 0)
		{
			if ((min_area.empty() ? 1E+30 : min_area[0]) > area[k])
			{
				min_area.insert(min_area.begin(), area[k]);
				min_ind.insert(min_ind.begin(), indices[k]);
				k_at_min.insert(k_at_min.begin(), k);
			}
			else
			{
				min_area.emplace_back(area[k]);
				min_ind.emplace_back(indices[k]);
				k_at_min.emplace_back(k);
			}
		}
	}

	int select_num = k_at_min[0]; // 20200716のところ最小のものから埋めていく他ない

	for (size_t k = 0; k < k_at_min.size(); k++)
	{
		if (area[k_at_min[k]] > Mean(area) * 0.01 /*too small*/)
		{
			select_num = k_at_min[k];
			break;
		}
	}

	if (area[select_num] < 0.01)
		Print(area[select_num], Red);

	// eraseの前にretを準備
	std::vector<int> ret = {indices[(select_num + 1) % indices.size()],
							indices[select_num - 1 == -1 ? indices.size() - 1 : select_num - 1]};
	// area.erase(area.begin()+select_num);
	vertices.erase(vertices.begin() + select_num);
	indices.erase(indices.begin() + select_num);
	//結んだ点の番号

	return ret;
};

std::vector<int> removeMinArea(const V_d &area,
							   VV_d &vertices,
							   std::vector<int> &indices)
{
	V_d min_area;
	std::vector<int> min_ind, k_at_min;
	for (size_t k = 0; k < indices.size(); k++)
	{
		if (area[k] > 0 && (min_area.empty() ? 1E+30 : min_area[0]) > std::abs(area[k]))
		{
			min_area.insert(min_area.begin(), area[k]);
			min_ind.insert(min_ind.begin(), indices[k]);
			k_at_min.insert(k_at_min.begin(), k);
		}
		else if (std::abs(area[k]) > 0)
		{
			min_area.emplace_back(area[k]);
			min_ind.emplace_back(indices[k]);
			k_at_min.emplace_back(k);
		}
	}

	int select_num = k_at_min[0]; // 20200716のところ最小のものから埋めていく他なない

	// eraseの前にretを準備
	std::vector<int> ret = {indices[select_num + 1 == indices.size() ? 0 : select_num + 1],
							indices[select_num - 1 == -1 ? indices.size() - 1 : select_num - 1]};
	// area.erase(area.begin()+select_num);
	vertices.erase(vertices.begin() + select_num);
	indices.erase(indices.begin() + select_num);
	//結んだ点の番号

	return ret;
};

// VV_i PolygonTriangulation(const VV_d& abcd,
// 						   const V_d& n){
//   int s = abcd.size(), i=0;
//   VV_d vertices = abcd;
//   std::vector<int> indices(s);
//   for(auto& ind:indices)
//     ind = i++;
//   VV_i connectInd;
//   while(true){
//     connectInd.emplace_back(removePositiveMinArea(DirectedArea(vertices,n),
// 					       vertices,
// 					       indices));
//     if(vertices.size()==3)
//       break;
//   }
//   return connectInd;
// };

// VV_i PolygonTriangulation(const VV_d& abcd){
//   int s = abcd.size(), i=0;
//   VV_d vertices = abcd;//順番に並んでいる
//   std::vector<int> indices(s);
//   for(auto& ind:indices)
//     ind = i++;
//   VV_i connectInd;
//   while(true){
//     connectInd.emplace_back(removeMinArea(TriangleArea(vertices),
// 				       vertices,
// 				       indices));
//     if(vertices.size()==3)
//       break;
//   }
//   return connectInd;
// };

// template <class T>
// std::vector<T> DeleteDuplicates(const std::vector<T> &vec)
// {
// 	std::vector<T> ret(0);
// 	ret.reserve(vec.size());
// 	for (const auto &v : vec)
// 	{
// 		if (std::find(ret.begin(), ret.end(), v) == ret.end())
// 			ret.emplace_back(v);
// 	}
// 	return ret;
// };

template <class T>
std::vector<T> DeleteDuplicates(const std::vector<T> &V)
{
	std::vector<T> ret(0);
	for (const auto &v : V)
		if (!std::binary_search(ret.begin(), ret.end(), v))
			ret.insert(std::lower_bound(ret.begin(), ret.end(), v), v);
	return ret;
};
template <class T>
std::vector<T *> DeleteDuplicates(const std::vector<T *> &V)
{
	std::vector<T *> ret(0);
	for (const auto &v : V)
		if (!std::binary_search(ret.begin(), ret.end(), v))
			ret.insert(std::lower_bound(ret.begin(), ret.end(), v), v);
	return ret;
};
template <class T>
bool DuplicateFreeQ(const std::vector<T> &vec)
{
	for (auto i = 0; i < vec.size(); i++)
		for (auto j = i + 1; j < vec.size(); j++)
			if (vec[i] == vec[j])
				return false;
	return true;
};
template <class T>
bool DuplicateFreeQ(const std::vector<T *> &vec)
{
	for (auto i = 0; i < vec.size(); i++)
		for (auto j = i + 1; j < vec.size(); j++)
			if (vec[i] == vec[j])
				return false;
	return true;
};
//===============
template <class T>
T Product(const std::vector<T> &vec)
{
	T ret = 1;
	for (const auto &v : vec)
		ret = ret * v;
	return ret;
};
//===============
VV_d nearestPointsOfLines(const VV_d &a0a1,
						  const VV_d &b0b1)
{
	V_d a1 = a0a1[1], a0 = a0a1[0];
	V_d b1 = b0b1[1], b0 = b0b1[0];
	V_d va = a1 - a0;
	V_d vb = b1 - b0;
	V_d a0Tob0 = b0b1[0] - a0a1[0];
	V_d n = Cross(va, vb), n1 = Cross(va, n), n2 = Cross(vb, n);
	return {a0 + va * Dot(a0Tob0, n2) / Dot(va, n2), b0 + vb * Dot(-(a0Tob0), n1) / Dot(vb, n1)};
};
V_d midPointOfLines(const VV_d &a0a1,
					const VV_d &b0b1)
{
	VV_d points = nearestPointsOfLines(a0a1, b0b1);
	return (points[0] + points[1]) / 2.;
};
//=============
// 20200817
V_d midPointOfLines(const VVV_d &vecs)
{
	// vecs is {{a0,a1},{b0,b1},{c0,c1},...}
	int n = vecs.size() - 1;
	VV_d A(vecs.size(), V_d(3));
	VV_d V(vecs.size(), V_d(3));

	for (size_t i = 0; i < vecs.size(); i++)
	{
		A[i] = vecs[i][0];
		V[i] = vecs[i][1] - vecs[i][0];
	}

	V_d rhs(vecs.size());

	for (size_t j = 0; j < vecs.size(); j++)
		for (size_t i = 0; i < vecs.size(); i++)
			rhs[j] -= (double)((i == j) ? -n : 1.) * Dot(V[j], A[i]);

	VV_d tmp(vecs.size(), V_d(vecs.size()));

	for (size_t j = 0; j < vecs.size(); j++)
		for (size_t i = 0; i < vecs.size(); i++)
			tmp[i][j] = (double)((i == j) ? -n : 1.) * Dot(V[j], V[i]);

	auto param = Dot(rhs, Inverse(tmp));

	V_d ret{0, 0, 0};
	for (size_t i = 0; i < vecs.size(); i++)
		ret += A[i] + V[i] * param[i];

	return ret / ((double)(vecs.size()));
};

// 20200818

double RandomReal(const Tdd &minmax)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(std::get<0>(minmax), std::get<1>(minmax));
	return dis(gen);
};

Tdd RandomRealTdd(double min, double max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(min, max);
	return {dis(gen), dis(gen)};
};
Tddd RandomRealTddd(double min, double max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(min, max);
	return {dis(gen), dis(gen), dis(gen)};
};
T3Tdd RandomRealT3Tdd(double min, double max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(min, max);
	return {{dis(gen), dis(gen)}, {dis(gen), dis(gen)}, {dis(gen), dis(gen)}};
};
//-----------------------
template <class T>
std::vector<T *> Intersection(const std::vector<T *> &A, const std::vector<T *> &B)
{
	std::vector<T *> ret({});
	ret.reserve((A.size() > B.size()) ? B.size() : A.size());
	for (const auto &b : B)
	{
		if (std::find_if(A.cbegin(), A.cend(), [&](const auto &a)
						 { return a == b; }) != A.end())
			ret.emplace_back(b);
	}
	return ret;
};
template <class T>
std::vector<T> Intersection(const std::vector<T> &A, const std::vector<T> &B)
{
	std::vector<T> ret({});
	ret.reserve((A.size() > B.size()) ? B.size() : A.size());
	for (const auto &b : B)
	{
		if (std::find_if(A.cbegin(), A.cend(), [&](const auto &a)
						 { return a == b; }) != A.end())
			ret.emplace_back(b);
	}
	return ret;
};
//------------------------
VV_d CoordinateBounds(const VV_d &v)
{
	// v = {{x0,y0,z0},{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}}
	VV_d xyzT = Transpose(v);
	V_d r(v.size(), 0);
	V_d mean = Mean(v); // mean = {xc,yc,zc}
	std::transform(v.cbegin(), v.cend(), r.begin(), [&mean](const auto &xyz)
				   { return Norm(xyz - mean); });
	VV_d bs = {{Min(xyzT[0]), Max(xyzT[0])},
			   {Min(xyzT[1]), Max(xyzT[1])},
			   {Min(xyzT[2]), Max(xyzT[2])},
			   {Min(r), Max(r)}}; // r:中心から各点までの距離
	return bs;
};

// class BoundingBox
// {
// public:
//   V_d center, xbounds, ybounds, zbounds, rbounds;
//   VV_d bounds;
//   BoundingBox(const VV_d &v) { set(v); };
//   VV_d getBounds() const { return {xbounds, ybounds, zbounds, rbounds}; };
//   double getScale() const
//   {
//     return Norm(V_d{std::abs(this->xbounds[1] - this->xbounds[0]),
//                     std::abs(this->ybounds[1] - this->ybounds[0]),
//                     std::abs(this->zbounds[1] - this->zbounds[0])});
//   };
//   void setBounds(const VV_d &bs)
//   {
//     double rand = (RandomReal({-1, 1.}) > 0 ? 1. : -1) * RandomReal({0.1, 1.}) * 1E-13;
//     this->xbounds = {bs[0][0], bs[0][1]};
//     this->ybounds = {bs[1][0], bs[1][1]};
//     this->zbounds = {bs[2][0], bs[2][1]};
//     this->rbounds = {bs[3][0], bs[3][1]};
//     this->xbounds += rand;
//     this->ybounds += rand;
//     this->zbounds += rand;
//     this->center = V_d{Mean(this->xbounds), Mean(this->ybounds), Mean(this->zbounds)};
//   };
//   void setBoundsFromPoints(const VV_d &v)
//   {
//     VV_d xyzT = Transpose(v);
//     V_d r(v.size(), 0.);
//     V_d center = Mean(v); //mean = {xc,yc,zc}
//     std::transform(v.cbegin(), v.cend(), r.begin(), [&center](const auto &xyz) { return Norm(xyz - center); });
//     this->xbounds = {Min(xyzT[0]), Max(xyzT[0])};
//     this->ybounds = {Min(xyzT[1]), Max(xyzT[1])};
//     this->zbounds = {Min(xyzT[2]), Max(xyzT[2])};
//     this->rbounds = {Min(r), Max(r)}; //r:中心から各点までの距離
//     this->bounds = getBounds();
//   };
// };

// //境界面の重なりを検知するのみ，完全に内部に入った場合はfalse
// bool isBoundariesOverlap(const BoundingBox &boxA, const BoundingBox &boxB)
// {
//   auto boundsA = boxA.bounds;
//   auto boundsB = boxB.bounds;
//   if ((boundsA[0][1] /*max*/ < boundsB[0][0] /*min*/ || boundsB[0][1] /*max*/ < boundsA[0][0] /*min*/) ||
//       (boundsA[1][1] /*max*/ < boundsB[1][0] /*min*/ || boundsB[1][1] /*max*/ < boundsA[1][0] /*min*/) ||
//       (boundsA[2][1] /*max*/ < boundsB[2][0] /*min*/ || boundsB[2][1] /*max*/ < boundsA[2][0] /*min*/))
//     return false;

//   return true;
// };

template <class T>
bool AnyTrue(const std::vector<T> &vec)
{
	for (const auto &o : vec)
		if (o)
			return true;
	return false;
}

template <class T>
bool AllTrue(const std::vector<T> &vec)
{
	for (const auto &o : vec)
		if (!o)
			return false;
	return true;
}

/* ------------------------------------------------------ */
/*                       連立1次方程式の解法                 */
/* ------------------------------------------------------ */

template <typename U, typename T>
class SLE
{
	//マップから行列を生成する
	// mapと違い，行列は，行列と，それにかかるベクトルの並びが正しくなくてはならない．
	//正しい並びで生成する方法を考える．
	//*1)indexをもつ新たなmapを生成し，それを使って，行列に割り振っていく
	//*2)pointerをソートして，そのソートした順番に従う．バイナリーサーチ
public:
	VV_d A;
	V_d b;
	std::vector<T *> keys;
	std::vector<U *> eqs_keys;
	//! ------------------------------------------------------ */
	SLE(const std::map<T *, std::map<T *, double>> &A_IN,
		const std::map<T *, double> &b_IN)
		: keys(A_IN.size(), nullptr),
		  A(A_IN.size(), V_d(A_IN.size(), 0.)),
		  b(b_IN.size(), 0.)
	{
		if (A_IN.size() != b_IN.size())
		{
			std::stringstream ss;
			ss << "A.size() = " << A_IN.size() << ", b.size() = " << b.size();
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
		}
		//* -------------------- 使用するキーを生成，ソート ------------------- */
		//このキーは初期に作成したあとは，固定．
		this->keys = TakeFirst(A_IN);
		std::sort(this->keys.begin(), this->keys.end());
		//この順番がマトリックスの添字となる
		//これ以降index()が使える
		//* -------------------- 行列とベクトルを生成する -------------------- */
		int i = 0;
		for (const auto &[p, m] : A_IN) //必ず全ての点がある
		{
			auto &row = this->A[index(p)];
			for (const auto &[q, v] : m)
				if ((i = index(q)) >= 0)
					row[i] = v;
		}
		for (const auto &[p, v] : b_IN)
			this->b[index(p)] = v;
		Print("eq has made", Red);
		Print("A.size() = " + std::to_string(A.size()), Red);
		Print("b.size() = " + std::to_string(b.size()), Red);
		Print("keys.size() = " + std::to_string(keys.size()), Red);
	};
	//! ------------------------------------------------------ */
	int index(const T *key) const
	{
		if (std::binary_search(this->keys.cbegin(), this->keys.cend(), key))
		{
			return std::distance(this->keys.cbegin(), std::lower_bound(this->keys.cbegin(), this->keys.cend(), key));
		}
		else
			return -1;
	};
	//! ------------------------------------------------------ */
	int eqs_index(const T *eq_key) const
	{
		if (std::binary_search(this->eqs_keys.cbegin(), this->eqs_keys.cend(), eq_key))
		{
			return std::distance(this->eqs_keys.cbegin(), std::lower_bound(this->eqs_keys.cbegin(), this->eqs_keys.cend(), eq_key));
		}
		else
			return -1;
	};
	//! ------------------------------------------------------ */
	SLE(const std::vector<T *> &keys_IN /*variables*/) : keys(keys_IN)
	{
		//* -------------------- 使用するキーを生成，ソート ------------------- */
		std::sort(this->keys.begin(), this->keys.end());
	};
	//
	SLE(const std::vector<U *> &eqs_keys_IN, const std::vector<T *> &keys_IN) : eqs_keys(eqs_keys_IN), keys(keys_IN)
	{
		std::sort(this->keys.begin(), this->keys.end());
		std::sort(this->eqs_keys.begin(), this->eqs_keys.end());
	};
	//! ------------------------------------------------------ */
	//単純に，方程式を追加したい場合，方程式に対してインデックスを割り振る必要はない．
	void add(const std::map<U *, std::map<T *, double>> &A_IN,
			 const std::map<U *, double> &b_IN = {})
	{
		//* ------------------------------------------------------ */
		auto eqn = TakeFirst(A_IN);
		std::sort(eqn.begin(), eqn.end());
		auto index_eq = [&eqn](const T *key) -> int
		{
			if (std::binary_search(eqn.cbegin(), eqn.cend(), key))
			{
				return std::distance(eqn.cbegin(), std::lower_bound(eqn.cbegin(), eqn.cend(), key));
			}
			else
				return -1;
		};
		//* -------------------- 行列とベクトルを生成する -------------------- */
		int i = 0;
		VV_d newA(eqn.size(), V_d(this->keys.size(), 0.));
		V_d newb(eqn.size(), 0.);
		for (const auto &[p, m] : A_IN) //必ず全ての点がある
		{
			auto &row = newA[index_eq(p)];
			for (const auto &[q, v] : m)
				if ((i = index(q)) >= 0)
					row[i] = v;
		}
		for (const auto &[p, v] : b_IN)
			newb[index_eq(p)] = v;
		//
		this->A = Join(this->A, newA);
		this->b = Join(this->b, newb);
	};
	//! ------------------------------------------------------ */
	//方程式を修正する場合は，修正する方程式を識別するインデックスが必要
	//インデックスは初期に作っておく必要がある．
	void modify(const std::map<U *, std::map<T *, double>> &A_IN,
				const std::map<U *, double> &b_IN = {})
	{
		try
		{
			int i = 0, j = 0;
			for (const auto &[p, m] : A_IN) //必ず全ての点がある
			{
				j = this->eqs_index(p);
				if (j >= 0)
				{
					auto &row = this->A[this->eqs_index(p)];
					for (const auto &[q, v] : m)
						if ((i = this->index(q)) >= 0)
							row[i] += v; //! ここでいう修正とは，行列成分への足し算のこと
				}
			}
			for (const auto &[p, v] : b_IN)
			{
				j = this->eqs_index(p);
				if (j >= 0)
				{
					this->b[eqs_index(p)] += v;
				}
			}
		}
		catch (const error_message &e)
		{
			e.print();
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
		}
	};
	//! ------------------------------------------------------ */
};

#endif
