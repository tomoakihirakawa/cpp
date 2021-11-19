#ifndef fundamental_statistics_H
#define fundamental_statistics_H

#include <cmath>
#include <vector>

#include "fundamental_vectors.hpp"
using V_d = std::vector<double>;
using V_i = std::vector<int>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

template <class T>
double Mean(const std::vector<T> &v)
{
	if (v.empty())
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "An empty vector is passed");
	else if (v.size() == 1)
		return v[0];
	return std::accumulate(v.cbegin(), v.cend(), 0.) / v.size();
};
V_d Mean(const VV_d &v)
{
	if (v.empty())
		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "An empty vector is passed");
	else if (v.size() == 1)
		return v[0];
	VV_d W = Transpose(v);
	V_d ret(W.size());
	for (size_t i = 0; i < ret.size(); i++)
		ret[i] = Mean(W[i]);
	return ret;
};
double Variance(const V_d &v)
{
	int N = v.size();
	double ret(0);
	for (const auto &w : v)
		ret = ret + w * w;
	return ret / N - pow(Mean(v), 2.);
};
/* ------------------------------------------------------ */
// 2021/08/06追加
// #define debug_histogram
class Histogram
{
public:
	V_d data;			  //ソートした生データ
	V_i count;			  //角瓶の中にあるデータの数
	V_i cumulative_count; //角瓶の中にあるデータの数を初めから累積してった結果
	int data_size;		  //全dataの長さ
	int sturges;		  //dataの数に応じてbinの幅を決める
	double bin_width;
	VV_d bins;
	V_i diff;
	V_d interval;
	V_d mid_interval;
	Histogram(const V_d &dataIN)
	{
		this->data = dataIN;
		std::sort(this->data.begin(), this->data.end(), [](double a, double b)
				  { return a < b; });
		this->data_size = dataIN.size();

#if defined(debug_histogram)
		std::cout << "data size: " << this->data_size << std::endl;
#endif
		this->sturges = (int)(std::log2(this->data_size) + 3);

#if defined(debug_histogram)
		std::cout << "sturges: " << this->sturges << std::endl;
#endif
		this->bins.resize(this->sturges, {});
		auto min_data = this->data[0];
		this->bin_width = (*this->data.rbegin() - min_data) / this->sturges;
#if defined(debug_histogram)
		std::cout << "bin_width: " << this->bin_width << std::endl;
#endif
		int n;
		for (auto i = 0; i < this->data.size(); i++)
		{
#if defined(debug_histogram)
			std::cout << "bins: " << this->bins << ", data: " << this->data[i] << ", floor((this->data[i] - min_data) / this->bin_width)=" << floor((this->data[i] - min_data) / this->bin_width) << std::endl;
#endif
			n = floor((this->data[i] - min_data) / this->bin_width);
			this->bins[n >= bins.size() ? bins.size() - 1 : n].emplace_back(this->data[i]);
		}
		//
#if defined(debug_histogram)
		std::cout << "bins: " << this->bins << std::endl;
#endif
		//
		this->diff = {};
		for (auto i = 0; i < this->bins.size() - 1; i++)
			this->diff.emplace_back((this->bins[i + 1].size() - this->bins[i].size()));
			//
#if defined(debug_histogram)
		std::cout << "diff: " << this->diff << std::endl;
#endif
		this->count.resize(this->bins.size(), 0);
		for (auto i = 0; i < this->count.size(); i++)
			this->count[i] = this->bins[i].size();
			//
#if defined(debug_histogram)
		std::cout << "count: " << this->count << std::endl;
#endif
		this->interval.resize(this->bins.size() + 1, 0);
		for (auto i = 0; i < this->interval.size(); i++)
			this->interval[i] = i * bin_width + min_data;
			//
#if defined(debug_histogram)
		std::cout << "interval: " << this->interval << std::endl;
#endif
		this->mid_interval.resize(this->interval.size() - 1, 0);
		for (auto i = 0; i < this->mid_interval.size(); i++)
			this->mid_interval[i] = (this->interval[i] + this->interval[i + 1]) / 2.;
		//
		this->cumulative_count.resize(this->bins.size(), 0);
		this->cumulative_count[0] = this->count[0];
		for (auto i = 1; i < this->cumulative_count.size(); i++)
			this->cumulative_count[i] = this->cumulative_count[i - 1] + this->count[i];
	};
};

#endif