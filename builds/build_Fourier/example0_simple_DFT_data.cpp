#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include "interpolations.hpp"

// 窓関数を生成する関数（例：ハニング窓）
std::vector<double> generateHanningWindow(int N) {
   std::vector<double> window(N);
   for (int i = 0; i < N; ++i) {
      window[i] = 0.5 * (1 - std::cos(2 * M_PI * i / (N - 1)));
   }
   return window;
}

std::complex<double> coeff(const std::vector<double>& sample, const int n) {
   int N = sample.size();
   std::complex<double> sum = 0;
   for (int k = 0; k <= N - 1; ++k) {
      sum += std::polar(sample[k], -n * 2 * M_PI / N * k);
   }
   return sum / static_cast<double>(N);
}

std::vector<std::complex<double>> DFT(const std::vector<double>& sample) {
   int N = sample.size();
   std::vector<std::complex<double>> result(N);
   for (int n = 0; n < N; ++n)
      result[n] = coeff(sample, n);
   return result;
}

std::vector<double> InverseDFT(const std::vector<std::complex<double>>& cn) {
   int N = cn.size();
   double dt;
   std::complex<double> sum = 0;
   std::vector<double> result(N);
   for (int k = 0; k < N; ++k) {
      sum = 0;
      dt = 2. * M_PI / N * k;
      for (int n = 0; n < N; ++n)
         sum += cn[n] * std::polar(1., n * dt);
      result[k] = sum.real();
   }
   return result;
}

int main() {
   double use_from_time = 2;  //! [s]
   std::string file_name = "data4";

   // ファイルの読み込み
   std::ifstream ifs("./input_data/" + file_name + ".csv");
   if (!ifs.is_open()) {
      std::cerr << "Error: file not found" << std::endl;
      return 1;
   }

   std::vector<double> value, time;
   std::string line;
   std::getline(ifs, line);  // 1行目を読むが使わない

   double last_time = 0;
   double first_time = -999.;
   double mean = 0;
   for (; std::getline(ifs, line);) {
      std::array<double, 2> d;
      sscanf(line.c_str(), "%lf,%lf", &d[0], &d[1]);
      if (line[0] == '#' || d[0] <= use_from_time)
         continue;
      time.push_back(d[0]);
      value.push_back(d[1]);
      last_time = d[0];
      if (first_time <= -999.)
         first_time = d[0];
   }

   value -= Mean(value);

   InterpolationBspline intB(4, time, value);

   time.clear();
   value.clear();
   double dt = 0.001;
   for (auto t = first_time; t < last_time; t += dt) {
      time.push_back(t);
      value.push_back(intB(t));
   }

   {
      std::ofstream ofs("./input_data/" + file_name + "_interpolated.csv");
      for (int i = 0; i < time.size(); ++i)
         ofs << time[i] << ", " << value[i] << std::endl;
      ofs.close();
   }

   // 窓関数の生成と適用
   std::vector<double> window = generateHanningWindow(value.size());
   for (size_t i = 0; i < value.size(); ++i) {
      value[i] *= window[i];
   }

   std::ofstream ofs("./input_data/" + file_name + "_DFT.csv");
   auto cn = DFT(value);

   auto n2freq = [&](const int i) {
      return i / (last_time - use_from_time);
   };

   auto n2time = [&](const int i) {
      return i * (last_time - use_from_time) / value.size();
   };

   // 30Hz以上の周波数成分を0にする
   int i = 0;
   for (auto& c : cn) {
      double freq = 0;
      if (i > (value.size() / 2.))
         freq = (value.size() - i) / (last_time - use_from_time);
      else
         freq = i / (last_time - use_from_time);
      if (freq > 30)
         c = 0.;
      i++;
   }

   for (i = 0; auto&& c : cn)
      ofs << n2freq(i++) << ", " << c.real() << "," << c.imag() << std::endl;
   ofs.close();

   std::vector<double> inv = InverseDFT(cn);
   std::ofstream ofs2("./input_data/" + file_name + "_invDFT.csv");
   for (i = 0; auto&& c : inv)
      ofs2 << n2time(i++) << ", " << c << std::endl;
   ofs2.close();
}
