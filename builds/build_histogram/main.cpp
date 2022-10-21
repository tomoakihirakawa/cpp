#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "fundamental.hpp"

using V_d = std::vector<double>;
using V_i = std::vector<int>;
using VV_d = std::vector<std::vector<double>>;

int main()
{
  /* ------------------------------------------------------ */
  class Histogram
  {
  public:
    V_d data;             //ソートした生データ
    V_i count;            //角瓶の中にあるデータの数
    V_i cumulative_count; //角瓶の中にあるデータの数を初めから累積してった結果
    int data_size;        //全dataの長さ
    int sturges;          //dataの数に応じてbinの幅を決める
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
      this->sturges = (int)(std::log2(this->data_size) + 1);
      this->bins.resize(this->sturges, {});
      auto min_data = this->data[0];
      this->bin_width = (*this->data.rbegin() - min_data) / this->sturges;
      for (auto i = 0; i < this->data.size() - 1; i++)
        this->bins[floor((this->data[i] - min_data) / this->bin_width)].emplace_back(this->data[i]);
      this->bins.rbegin()->emplace_back(*this->data.rbegin());
      //
      this->diff = {};
      for (auto i = 0; i < this->bins.size() - 1; i++)
        this->diff.emplace_back((this->bins[i + 1].size() - this->bins[i].size()));
      //
      this->count.resize(this->bins.size(), 0);
      for (auto i = 0; i < this->count.size(); i++)
        this->count[i] = this->bins[i].size();
      //
      this->interval.resize(this->bins.size() + 1, 0);
      for (auto i = 0; i < this->interval.size(); i++)
        this->interval[i] = i * bin_width + min_data;
      //
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

  /* ------------------------------------------------------ */
  // random data
  // h=RandomVariate[NormalDistribution[0, 1], 100]
  // Histogram[h]
  V_d data = {-1.45114, 0.0175009, -0.477811, -0.000444694, 0.34841, 1.25773, -0.0843609, 1.01718, 1.93293, -0.711627, -0.50765, 0.343237, 1.34123, 0.256944, -1.76313, -2.03918, 1.76659, 1.5619, -1.81814, -0.162943, -0.221754, 1.94962, 0.621411, -0.778915, -0.0465956, 0.671615, -0.396592, -0.241664, -1.49107, -0.791977, 0.14373, 1.40084, -0.877371, 0.445416, -0.0222211, -0.969357, -0.777598, 0.307552, 2.10962, 0.968876, 0.615126, 0.515048, 0.228766, -1.72247, 0.190768, -0.616721, 1.07264, 2.13789, -0.802044, 0.749316, 0.681747, -0.141591, -0.079582, -1.37473, -0.821742, 0.870849, -1.32781, -0.640009, 1.50405, -0.208336, 1.3714, 0.332523, 0.123145, 1.05069, -0.750535, -0.181517, 0.867686, -0.0403632, -0.86369, -0.21725, -1.11271, -0.68791, -0.804072, -0.688389, -0.193886, 0.970914, 2.20058, -0.115096, -0.00476043, -1.75345, 2.16769, -0.85975, -0.704658, -1.17081, -0.0716246, 2.169, 0.244943, -1.2936, -2.26912, -1.08769, 0.399829, 0.54558, 2.0042, 0.0316217, 2.25, -1.37579, 1.37137, -2.27198, 0.165188, -0.742774};
  auto h = Histogram(data);
  std::cout << h.bins << std::endl;
  std::cout << h.count << std::endl;
  std::cout << h.cumulative_count << std::endl;
};
