#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

int main() {

   double N = 10.;
   std::vector<std::array<double, 3>> k1;
   for (int j = 1; j < 10; ++j) {
      double r = j * 0.1;
      for (int i = 0; i < N; ++i) {
         double t = i / N * 2 * M_PI;
         k1.push_back({r, std::cos(t), std::sin(t)});
      }
   }
   std::vector<std::array<double, 3>> k2 = k1;
   std::vector<std::array<double, 3>> k3 = k1;
   std::vector<std::array<double, 3>> k4 = k1;
}