#ifndef minMaxOfFunctions_H
#define minMaxOfFunctions_H
#pragma once

#include <functional>
#include "basic_linear_systems.hpp"
// #include "svd.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

/* -------------------------------------------------------------------------- */
/*                               Gradient Method                              */
/* -------------------------------------------------------------------------- */

struct GradientMethod {
  public:
   VV_d A;
   VV_d org_A;
   int count;

   bool use_diagonal_scaling;
   V_d diagonal_scaling_vector;
   GradientMethod(const VV_d &A_IN) : org_A(A_IN), A(A_IN), use_diagonal_scaling(false){};

   void diagonal_scaling() {
      this->use_diagonal_scaling = true;
      this->diagonal_scaling_vector = V_d(org_A.size());
      //
      int s = org_A.size();
      for (auto i = 0; i < s; ++i)
         this->diagonal_scaling_vector[i] = 1. / org_A[i][i];
      for (auto i = 0; i < s; ++i)
         this->A[i] = this->org_A[i] * this->diagonal_scaling_vector[i];
   };

   V_d solve(V_d b, const V_d &x_init = {}, double eps = 1E-10) {
      V_d x(b.size());
      if (use_diagonal_scaling)
         b *= this->diagonal_scaling_vector;

      for (auto i = 0; i < b.size(); ++i) {
         if (i >= x_init.size())
            x[i] = 0.;
         else if (!isFinite(x_init[i]))
            x[i] = 0.;
         else
            x[i] = x_init[i];
      }
      count = 0;
      double norm;
      V_d p = b - Dot(A, x);  // 修正ベクトルp
      norm = Norm(p);
      while (!(norm < eps)) {
         //  std::cout << "count = " << count << ", norm = " << norm << std::endl;
         p = b - Dot(A, x += Dot(p, p) / Dot(Dot(A, p), p) * p);
         if (!isFinite(norm = Norm(p)))
            return x;
         else if (count++ > 10000)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
      //   std::cout << "count = " << count << ", norm = " << norm << std::endl;
      return x;
   };
   V_d solve(V_d b, double eps) {
      return solve(b, V_d{}, eps);
   };
   /* ------------------------------------------------------ */
   V_d solveCG(V_d b, const V_d &x_init = {}, double eps = 1E-10) {
      if (use_diagonal_scaling)
         b *= this->diagonal_scaling_vector;
      V_d x(b.size());
      for (auto i = 0; i < b.size(); i++) {
         if (i >= x_init.size())
            x[i] = 0.;
         else if (!isFinite(x_init[i]))
            x[i] = 0.;
         else
            x[i] = x_init[i];
      }
      count = 0;
      V_d r = b - Dot(A, x);  // 残差ベクトル
      double norm = Norm(r);
      //   std::cout << "count = " << count << ", norm = " << norm << std::endl;
      if (norm < eps)
         return x;
      //
      V_d p = r;
      //
      double alpha = Dot(r, p) / Dot(Dot(A, p), p);
      x += alpha * p;  // xを修正
      r = b - Dot(this->A, x /*x=x0+alpha*p0*/);
      // check
      double beta;
      norm = Norm(r);
      while (!(norm < eps)) {
         //  std::cout << "count = " << count << ", norm = " << norm << std::endl;
         beta = -Dot(r, Dot(A, p)) / Dot(Dot(A, p), p);
         p = r + beta * p;  // 修正ベクトルは，最急降下方向を少し変更した物
         alpha = Dot(r, p) / Dot(Dot(A, p), p);
         x += alpha * p;                       // xを修正
         r = b - Dot(A, x /*x=x0+alpha*p0*/);  // 残差ベクトルr（最急降下方向）
         if (!isFinite(norm = Norm(r)))
            return x;
         else if (count++ > 1000)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
      //   std::cout << "count = " << count << ", norm = " << norm << std::endl;
      return x;
   };
   V_d solveCG(V_d b, double eps) {
      return solveCG(b, V_d{}, eps);
   };

   /* ------------------------------------------------------ */
   GradientMethod(){};
   V_d minimize(const std::function<double(V_d)> &F,
                const std::function<V_d(V_d)> &dF,
                V_d Xinit,
                double initial_step = 1E-12,
                int max_loop = 100) {
      V_d X = Xinit, X_last = Xinit, m(X.size(), 0.);
      double s = initial_step;
      for (auto i = 0; i < max_loop; ++i) {
         X_last = X;
         double min = 1E+10;
         int size = 50;
         double d, f, S, extend = 2;
         for (auto k = 0; k < size; ++k) {
            d = 2. * extend * s / (size - 1);
            S = d * k - extend * s;
            f = F(X_last + S * m);
            if (k == 0 || min >= f) {
               min = f;
               s = S;
            }
         }
         X = X_last + s * m;
         /* ------------------------------------------------------ */
         auto dF_X = dF(X);
         auto dF_Xlast = dF(X_last);
         double a = -Dot(dF_X, dF_X - dF_Xlast) / Dot(dF_Xlast, dF_Xlast);
         // double a = -Dot(dF_X, dF_X - dF_Xlast) / Dot(m, dF_X - dF_Xlast);
         /* ------------------------------------------------------ */
         m = dF(X) + a * m;
         // std::cout << "s = " << s << std::endl;
         std::cout << "i = " << i << std::endl;
         std::cout << "X = " << X << std::endl;
         std::cout << "F(X) = " << F(X) << std::endl;
         std::cout << "Norm(m) = " << Norm(m) << std::endl;
         if (Norm(m) < 1E-8)
            break;
      };
      return X;
   };
};

/* -------------------------------------------------------------------------- */
/*                               Broyden Method                               */
/* -------------------------------------------------------------------------- */
template <typename T>
struct BroydenMethod;

template <typename T>
   requires is_std_array<T>::value
struct BroydenMethod<T> {
   T X, dX;
   std::array<T, std::tuple_size<T>::value> J;

   BroydenMethod(const T &Xin, const T &Xin_) : X(Xin), dX(Xin_ - Xin), J(TensorProduct(Xin, Xin)) {
      IdentityMatrix(J);
   }

   void initialize(const T &Xin) { X = Xin; }

   void update(const T &F, const T &F_, const double alpha = 1.) {
      auto dot = Dot(dX, dX);
      if (dot != static_cast<double>(0.))
         J += TensorProduct((F - F_ - Dot(J, dX)), dX) / dot;
      X += (dX = -alpha * Dot(Inverse(J), F));  // inverseが計算できるかは未検証
   }
};

// Deduction guide
template <typename T1, typename T2>
BroydenMethod(T1, T2) -> BroydenMethod<std::decay_t<T1>>;

using Vd = std::vector<double>;
using VVd = std::vector<std::vector<double>>;

template <>
struct BroydenMethod<Vd> {
   Vd X, dX;
   VVd J;

   BroydenMethod(const Vd &Xin, const Vd &Xin_) : X(Xin), dX(Xin_ - Xin), J(VVd(Xin.size(), Vd(Xin.size()))) {
      IdentityMatrix(J);
   }

   void initialize(const Vd &Xin) { X = Xin; }

   void update(const Vd &F, const Vd &F_, const double alpha = 1.) {
      auto dot = Dot(dX, dX);
      if (dot != static_cast<double>(0.))
         J += TensorProduct((F - F_ - Dot(J, dX)), dX) / dot;
      X += (dX = -alpha * Dot(Inverse(J), F));
   }
};

#endif