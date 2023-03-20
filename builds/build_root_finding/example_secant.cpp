/*
ニュートン法で使うヤコビアンなどを別のものに置き換えた場合，準ニュートン法と呼ぶ
*/
#include "rootFinding.hpp"

double f(const double x) {
   return x * x + 20 * cos(x) - 10;
}

int main() {
   double x1 = RandomReal({-10, 10}), f_x1 = f(x1);
   double x2 = RandomReal({-10, 10}), f_x2 = f(x2);
   double X, F, X_, F_;
   double dX;
   double xacc = 1E-5;
   if (std::abs(f_x1) < std::abs(f_x2)) {
      X = x1;
      F = f_x1;
      X_ = x2;
      F_ = f_x2;
   } else {
      X = x2;
      F = f_x2;
      X_ = x1;
      F_ = f_x1;
   }
   std::cout << X << ", " << F << std::endl;
   double Xnext;
   for (auto j = 0; j < 20; j++) {
      double k = (F_ - F) / (X_ - X);
      dX = -F / k;
      X += dX;
      Xnext = X;
      //
      F = F_;
      X = X_;
      F_ = f(Xnext);
      X_ = Xnext;
      //
      std::cout << X << ", " << F << std::endl;
      if (std::abs(F) < xacc || F == 0.0)
         break;
   }
}