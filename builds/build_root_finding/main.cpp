#include "rootFinding.hpp"

// #define simple_test
#ifdef simple_test
double f(const double x) {
   return sin(x);
}
double dfdx(const double x) {
   return cos(x);
}

V_d F(const V_d &X) {
   return {sin(X[0]), cos(X[1]), sin(X[2])};
}
VV_d dFdx(const V_d &X) {
   return {{cos(X[0]), 0, 0},
           {0, -sin(X[1]), 0},
           {0, 0, cos(X[2])}};
}

int main() {
   /* ----------------------- 一次元の場合 ----------------------- */
   double xn = 1.;
   double xn1 = 0.;
   for (auto i = 0; i < 10; i++) {
      Print(xn);
      xn1 = xn - f(xn) / dfdx(xn);
      xn = xn1;
   }
   /* ----------------------- 多次元の場合 ----------------------- */
   V_d Xn = {1., 1., 1.};
   V_d Xn1(3), dX(3);
   for (auto i = 0; i < 10; i++) {
      Print(Xn);
      ludcmp lu(dFdx(Xn));
      lu.solve(-F(Xn), dX);
      Xn1 = Xn + dX;
      Xn = Xn1;
   }
};
#else

// class NewtonRaphson
// {
// public:
// 	V_d X;
// 	V_d dX; // tmp
// 	NewtonRaphson(const V_d &Xinit) : X(Xinit), dX(Xinit){};
// 	void update(const V_d &F, const VV_d &dFdx)
// 	{
// 		ludcmp lu(dFdx);
// 		lu.solve(-F, dX);
// 		X += dX;
// 	};
// };

double f(const double x) {
   return sin(x);
}
double dfdx(const double x) {
   return cos(x);
}

Tdd F2(const Tdd &X) {
   return {sin(std::get<0>(X)), cos(std::get<1>(X))};
}
T2Tdd dF2dx(const Tdd &X) {
   return {{cos(std::get<0>(X)), 0}, {0, -sin(std::get<1>(X))}};
}

Tddd F3(const Tddd &X) {
   return {sin(std::get<0>(X)), cos(std::get<1>(X)), sin(std::get<2>(X))};
}
T3Tddd dF3dx(const Tddd &X) {
   return {{cos(std::get<0>(X)), 0, 0},
           {0, -sin(std::get<1>(X)), 0},
           {0, 0, cos(std::get<2>(X))}};
};

T3Tddd X0X1X2 = {{10, 5, 1}, {-4, 5, 5}, {5, 10, 1}};

double G(const double t0, const double t1) {
   return Norm(Dot({t0, t1, 1 - t0 - t1}, X0X1X2));
};

double dGdt0(const double t0, const double t1) {
   auto X = Dot({t0, t1, 1 - t0 - t1}, X0X1X2);
   auto dXdt = Dot({1, 0, -1}, X0X1X2);
   return Dot(X, X) * Dot(X, dXdt);
};

double dGdt1(const double t0, const double t1) {
   auto X = Dot({t0, t1, 1 - t0 - t1}, X0X1X2);
   auto dXdt = Dot({0, 1, -1}, X0X1X2);
   return Dot(X, X) * Dot(X, dXdt);
};
//
double w(const double k, const double h) {
   return Sqrt(_GRAVITY_ * k * Tanh(h * k));
};

double dwdk(const double k, const double h) {
   return (_GRAVITY_ * (h * k * Power(Sech(h * k), 2) + Tanh(h * k))) / (2. * Sqrt(_GRAVITY_ * k * Tanh(h * k)));
};

int main() {
   /* ----------------------- 1次元の場合 ----------------------- */
   std::cout << "1次元の場合" << std::endl;
   {
      NewtonRaphson nr(4.2);
      for (auto i = 0; i < 10; i++) {
         nr.update(f(nr.X), dfdx(nr.X));
         std::cout << "nr.dX = " << nr.dX << ", f(nr.X) = " << f(nr.X) << std::endl;
      }
   }
   /* ----------------------- 1次元の場合 ----------------------- */
   std::cout << "1次元の場合" << std::endl;
   {
      NewtonRaphson nr0(1.);
      NewtonRaphson nr1(1.);
      for (auto i = 0; i < 500; i++) {
         nr0.update(G(nr0.X, nr1.X), dGdt0(nr0.X, nr1.X));
         nr1.update(G(nr0.X, nr1.X), dGdt1(nr0.X, nr1.X));
         std::cout << "{nr0.dX, nr1.dX} = " << Tdd{nr0.dX, nr1.dX} << ", G(nr0.X, nr1.X) = " << G(nr0.X, nr1.X) << std::endl;
      }
   }
   /* ----------------------- 2次元の場合 ----------------------- */
   std::cout << "2次元の場合" << std::endl;
   {
      NewtonRaphson nr(Tdd{1., 1.});
      for (auto i = 0; i < 10; i++) {
         nr.update(F2(nr.X), dF2dx(nr.X));
         std::cout << "nr.dX = " << nr.dX << ", F2(nr.X) = " << F2(nr.X) << std::endl;
      }
   }
   /* ----------------------- 3次元の場合 ----------------------- */
   std::cout << "3次元の場合" << std::endl;
   {
      NewtonRaphson nr(Tddd{1., 1., 1.});
      for (auto i = 0; i < 10; i++) {
         nr.update(F3(nr.X), dF3dx(nr.X));
         std::cout << "nr.dX = " << nr.dX << ", F3(nr.X) = " << F3(nr.X) << std::endl;
      }
   }
   /* -------------------------------------------------------------------------- */
   std::cout << "分散関係の lambda -> omega" << std::endl;
   for (const auto &T : Subdivide({0.1, 20}, 50)) {
      double h = 1000, omega = 2 * M_PI / T;
      NewtonRaphson nr(1.);
      for (auto i = 0; i < 10; i++)
         nr.update(w(nr.X, h) - omega, dwdk(nr.X, h));
      // std::cout << "分散関係を満たす組：{w,k,h} = " << Tddd{omega, nr.X, h} << ", w(nr.X, h) - omega = " << w(nr.X, h) - omega << std::endl;
      auto ds = DispersionRelation(omega, h);
      std::cout << "{水深h, 角周波数w, 波数k, 周期T, 波長L} = " << T5d{ds.h, ds.w, ds.k, ds.T, ds.L} << std::endl;
   }
};
#endif