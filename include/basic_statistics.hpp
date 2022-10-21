#ifndef basic_statistics_H
#define basic_statistics_H

#include <cmath>
#include <vector>
#include "basic_vectors.hpp"
/* ------------------------------------------------------ */
V_d Mean(const VV_d &v) {
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
double Variance(const V_d &v) {
   double ret(0);
   for (const auto &w : v)
      ret += w * w;
   return ret / (v.size()) - std::pow(Mean(v), 2.);
};
template <typename T>
double RootMeanSquare(const T &TUPLE) {
   double ret = 0, div = 0;
   for_each(TUPLE, [&](const auto &t) { ret += t * t; div += 1.; });
   return std::sqrt(ret / div);
};
/* ------------------------------------------------------ */
// 2021/08/06追加
// #define debug_histogram
class Histogram {
  public:
   V_d data;       //ソートした生データ
   int data_size;  //全dataの長さ
   int sturges;    // dataの数に応じてbinの幅を決める
   VV_d bins;
   double bin_width;
   V_i count;             //角瓶の中にあるデータの数
   V_i cumulative_count;  //角瓶の中にあるデータの数を初めから累積してった結果
   V_i diff;
   V_d interval;
   V_d mid_interval;

   V_d sort(V_d dataIN) const {
      std::sort(dataIN.begin(), dataIN.end(), [](const double a, const double b) { return a < b; });
      return dataIN;
   };

   Histogram(const V_d &dataIN)
       : data(this->sort(dataIN)),
         data_size(dataIN.size()),
         sturges((int)(std::log2(this->data_size) + 3)),
         bins(VV_d(this->sturges, V_d{})),
         bin_width((*this->data.rbegin() - data[0]) / sturges) {
      auto min_data = this->data[0];
#if defined(debug_histogram)
      std::cout << "bin_width: " << this->bin_width << std::endl;
#endif
      int n;
      for (const auto &d : this->data) {
         n = floor((d - min_data) / this->bin_width);
         this->bins[n >= bins.size() ? bins.size() - 1 : n].emplace_back(d);
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