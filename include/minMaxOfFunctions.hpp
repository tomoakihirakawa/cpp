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
   GradientMethod(const VV_d &A_IN) : org_A(A_IN), A(A_IN), use_diagonal_scaling(false) {};

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
   GradientMethod() {}
   void initialize(const V_d &Xin) {
      X = Xin;
      isFirst = true;
   }
   V_d X, dX, dFdX_last, s, s_last;
   bool isFirst;

   V_d update(const V_d &dFdX /*direction*/, const double a = 1.) {
      return X -= (dX = a * dFdX);
   };

   V_d getFletcherReevesCD(const V_d &dFdX) {
      auto dX = -dFdX;
      auto dX_last = -dFdX_last;
      double beta = Dot(dX, dX) / Dot(dX_last, dX_last);
      this->dFdX_last = dFdX;
      if (isFirst) {
         s_last = V_d(dFdX.size(), 0.);
         isFirst = false;
      }
      return s = s_last = -dFdX + beta * s_last;
   };

   V_d getPolakRibiereCD(const V_d &dFdX) {
      auto dX = -dFdX;
      auto dX_last = -dFdX_last;
      double beta = Dot(dX, dX - dX_last) / Dot(dX_last, dX_last);
      this->dFdX_last = dFdX;
      if (isFirst) {
         s_last = V_d(dFdX.size(), 0.);
         isFirst = false;
      }
      return s = s_last = -dFdX + beta * s_last;
   };

   // V_d update(const V_d &F, const V_d &dF_X) {
   //    V_d X_last = X, m(X.size(), 0.);
   //    double min = 1E+10;
   //    int size = 50;
   //    double d, f, S, extend = 2;
   //    for (auto k = 0; k < size; ++k) {
   //       d = 2. * extend * step / (size - 1);
   //       S = d * k - extend * step;
   //       f = F(X_last + S * m);
   //       if (k == 0 || min >= f) {
   //          min = f;
   //          step = S;
   //       }
   //    }
   //    X = X_last + step * m;
   //    if (this->isFirst) {
   //       dF_X_last = dF_X;
   //       m = dF_X;
   //    } else {
   //       double a = -Dot(dF_X, dF_X - dF_X_last) / Dot(dF_X_last, dF_X_last);
   //       m = dF_X + a * m;
   //    }
   //    dF_X_last = dF_X;
   //    return X;
   // };

   V_d minimize(const std::function<double(V_d)> &F, const std::function<V_d(V_d)> &dF, V_d Xinit, double initial_step = 1E-12, int max_loop = 100) {
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

/*DOC_EXTRACT Broydens_Method

準ニュートン法は，
ヤコビ行列を近似し，
ニュートン法を適用する．
近似ヤコビ行列を計算する指針として，セカント条件を満たすようにする．

```math
{\bf J}_{k} \cdot \Delta {\bf x}_k = \Delta {\bf f}_k
```


```math
{\bf J}_{k} = {\bf J}_{k-1} + \frac{(\Delta {\bf f}_k - {\bf J}_{k-1} \cdot \Delta {\bf x}_k) \otimes \Delta  {\bf x}_k}{\Delta  {\bf x}_k \cdot \Delta {\bf x}_k},\quad \Delta {\bf f}_k = {\bf f}_{k} - {\bf f}_{k-1}
```

*/

template <typename T>
struct BroydenMethod;

template <typename T>
   requires is_std_array<T>::value
struct BroydenMethod<T> {
   T X, dX;
   std::array<T, std::tuple_size<T>::value> J, Inv_J;

   // Default constructor
   BroydenMethod() = default;

   BroydenMethod(const T &Xin, const T &Xin_) : X(Xin), dX(Xin_ - Xin), J(TensorProduct(Xin, Xin)) {
      IdentityMatrix(J);
      Inv_J = J;
   }

   void initialize(const T &Xin, const T &dXin) {
      X = Xin;
      dX = dXin;
      J = TensorProduct(Xin, Xin);
      IdentityMatrix(J);
      Inv_J = J;
   }

   void initialize(const T &Xin) {
      X = Xin;
      J = TensorProduct(Xin, Xin);
      IdentityMatrix(J);
      Inv_J = J;
   }

   void updateGoodBroyden(const T &F, const T &F_, const double alpha = 1.) {
      auto dot = Dot(dX, dX);
      auto dF = F - F_;
      /* -------------------------------------------------------------------------- */
      // if (dot != static_cast<double>(0.))
      //    J += TensorProduct((dF - Dot(J, dX)), dX) / dot;
      // // use lapack_lu
      // lapack_lu(J, dX, -alpha * F);
      // X += dX;
      // X += (dX = -alpha * Dot(Inverse(J), F));  // inverseが計算できるかは未検証
      /* -------------------------------------------------------------------------- */
      // good Broyden's method
      if (Dot(dF, dF) != static_cast<double>(0.))
         Inv_J += TensorProduct((dX - Dot(Inv_J, dF)), Dot(dX, Inv_J)) / Dot(dX, Dot(Inv_J, dF));
      X += (dX = -alpha * Dot(Inv_J, F));
   }

   void updateBFGS(const T &F, const T &F_, const double alpha = 1.) {
      auto dot = Dot(dX, dX);
      auto dF = F - F_;
      auto s = dX;
      auto y = dF;
      double sy = Dot(s, y);
      if (sy != static_cast<double>(0.))
         J += TensorProduct(y, y) / Dot(y, s) - Dot(TensorProduct(Dot(J, s), s), Transpose(J)) / Dot(s, Dot(J, s));
      // X += (dX = -alpha * Dot(Inv_J, F));
      // Solve(J, dX, -alpha * F);
      // lapack_lu(J, dX, -alpha * F);
      lapack_svd(J, dX, -alpha * F);

      X += dX;
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
   VVd J, Inv_J;

   // Default constructor
   BroydenMethod() = default;

   BroydenMethod(const Vd &Xin, const Vd &Xin_) : X(Xin), dX(Xin_ - Xin), J(VVd(Xin.size(), Vd(Xin.size()))) {
      IdentityMatrix(J);
      Inv_J = J;
   }

   void initialize(const Vd &Xin) { X = Xin; }

   // void update(const Vd &F, const Vd &F_, const double alpha = 1.) {
   //    auto dot = Dot(dX, dX);
   //    if (dot != static_cast<double>(0.))
   //       J += TensorProduct((F - F_ - Dot(J, dX)), dX) / dot;

   //    auto dF = F - F_;
   //    // good Broyden's method
   //    if (Dot(dF, dF) != static_cast<double>(0.))
   //       Inv_J += TensorProduct((dX - Dot(Inv_J, dF)), Dot(dX, Inv_J)) / Dot(dX, Dot(Inv_J, dF));

   //    X += (dX = -alpha * Dot(Inv_J, F));
   //    // X += (dX = -alpha * Dot(Inverse(J), F));
   // }

   void update(const Vd &F, const Vd &F_, const double alpha = 1.) {
      // auto dot = Dot(dX, dX);
      // auto dF = F - F_;

      // auto s = dX;
      // auto y = dF;
      // double sy = Dot(s, y);
      // if (sy != static_cast<double>(0.))
      //    J += TensorProduct(y, y) / Dot(y, s) - Dot(TensorProduct(Dot(J, s), s), Transpose(J)) / Dot(s, Dot(J, s));

      // if (sy != static_cast<double>(0.))
      //    Inv_J += Dot(Dot(s, y) + Dot(Inv_J, y), TensorProduct(s, s)) / std::pow(Dot(s, y), 2) - TensorProduct(Dot(Inv_J, y), s) + Dot(TensorProduct(s, y), Inv_J) / Dot(s, y);

      // // X += (dX = -alpha * Dot(Inv_J, F));
      // X += (dX = -alpha * Dot(Inverse(J), F));
      /* -------------------------------------------------------------------------- */
      auto dot = Dot(dX, dX);
      auto dF = F - F_;
      // good Broyden's method
      if (Dot(dF, dF) != static_cast<double>(0.))
         Inv_J += TensorProduct((dX - Dot(Inv_J, dF)), Dot(dX, Inv_J)) / Dot(dX, Dot(Inv_J, dF));
      X += (dX = -alpha * Dot(Inv_J, F));
   }
};

#endif