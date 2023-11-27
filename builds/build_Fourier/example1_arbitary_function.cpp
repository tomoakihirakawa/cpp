#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include "rootFinding.hpp"

std::complex<double> coeff(const std::vector<double>& sample, const int n) {
   int N = sample.size();
   std::complex<double> sum = 0;
   for (int k = 0; k <= N - 2; ++k) {
      sum += std::polar(sample[k], -n * 2 * M_PI / N * k);
   }
   // sample[k] means f(k * dt) or f(k * T / N). T is the maximum time.
   return sum / static_cast<double>(N);
};

std::vector<std::complex<double>> DFT(const std::vector<double>& sample) {
   int N = sample.size();
   std::vector<std::complex<double>> result(N);
   for (int n = 0; n < N; ++n)
      result[n] = coeff(sample, n);
   return result;
}

// LighthillRobot class

double L = 0.71;
double T = 2.0;
double w = 2. * M_PI / T;
double k = 2. * M_PI / L;
double c1 = 0.1;
double c2 = 0.1;
int nodes = 5;
LighthillRobot lhr(L, w, k, c1, c2, nodes);

// const auto f = [](double t) { return lhr.getAngles(t); };
std::vector<double> f(double t) {
   return {cos(w * t), sin(w * t)};
};

int main() {

   const int vec_size = f(0.).size();
   const int steps = 1000;
   const double fundamental_T = 5 * T;
   double t = 0., dt = fundamental_T / steps;
   //! --------------------------- make dataOrg vector -------------------------- */
   std::vector<std::vector<double>> dataOrg(vec_size);
   std::ofstream dataOrg_file("./example1_dataOrg.dat");
   for (auto i = 0; i < steps; i++) {
      dataOrg_file << t << " ";
      int j = 0;
      for (const auto& q : f(t)) {
         dataOrg_file << q << " ";
         dataOrg[j++].push_back(q);  // head node is fixed
      }
      dataOrg_file << std::endl;
      t += dt;
   }
   dataOrg_file.close();
   //! ----------------------------- make DFT vector ---------------------------- */
   std::vector<std::vector<std::complex<double>>> dataDFT(vec_size, std::vector<std::complex<double>>(0));
   for (auto i = 0; i < vec_size; i++) {
      dataDFT[i] = DFT(dataOrg[i]);
   }

   std::ofstream output_dft_Re("./example1_dataDFT_Re.dat");
   std::ofstream output_dft_Im("./example1_dataDFT_Im.dat");
   // output row: freq, column: nodes

   for (auto j = 0; j < dataDFT[0].size(); j++) {
      output_dft_Re << j / fundamental_T << " ";
      output_dft_Im << j / fundamental_T << " ";
      for (auto i = 0; i < vec_size; i++) {
         output_dft_Re << dataDFT[i][j].real() << " ";
         output_dft_Im << dataDFT[i][j].imag() << " ";
      }
      output_dft_Re << std::endl;
      output_dft_Im << std::endl;
   }

   return 0;
}
