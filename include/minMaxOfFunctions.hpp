#ifndef minMaxOfFunctions_H
#define minMaxOfFunctions_H
#pragma once

#include <functional>
#include "basic_vectors.hpp"
#include "svd.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

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

/* ------------------------------------------------------ */
/*                          GMRES                         */
/* ------------------------------------------------------ */

V_d forward_substitution(const VV_d &mat, const V_d &b) {
   int s = b.size();
   V_d x(s);
   double tmp = 0;
   // for (auto i = 0; i < s; ++i)
   // {
   // 	tmp = 0;
   // 	for (auto j = 0; j < i; ++j)
   // 		tmp += mat[i][j] * x[j];
   // 	x[i] = (b[i] - tmp) / mat[i][i];
   // }

   int i = 0;
   for (const auto &a : mat) {
      tmp = 0;
      for (auto j = 0; j < i; ++j)
         tmp += a[j] * x[j];
      x[i] = (b[i] - tmp) / a[i];
      i++;
   }
   return x;
};

V_d back_substitution(const VV_d &mat, V_d b, const Tii &mat_size) {
   auto [row, col] = mat_size;
   for (auto i = row - 1; i >= 0; --i) { /*　0~row-1まで，長さは，row　*/
      for (auto j = col - 1; j > i; --j) /*　長さは，col-i-1　*/
         b[i] -= mat[i][j] * b[j];
      b[i] /= mat[i][i];
   }
   b.erase(std::next(b.begin(), row + 1), b.end());
   return b;
};

V_d back_substitution(const VV_d &mat, V_d b, const int &mat_size) {
   return back_substitution(mat, b, Tii{mat_size, mat_size});
};

struct ArnoldiProcess {
   // ヘッセンベルグ行列H[0:k-1]は，Aと相似なベクトルであり，同じ固有値を持つ
   // GMRESで使う場合，V0にはNormalize(b-A.x0)を与える．
   // x0は初期値
   // Therefore V is an orthonormal basis of the Krylov subspace Km(A,r0)
   //
   // https://en.wikipedia.org/wiki/Arnoldi_iteration
   // アーノルディ法は固有値問題の数値解法であり反復解法である．
   // 一般的な行列の固有ベクトルと固有値を
   // クリロフ空間の直行基底によって近似する方法計算する方法である．
   /* ------------------------------------------------------ */
   int n;  // the number of interation
   double beta;
   V_d v0;
   VV_d H;  // ヘッセンベルグ行列
   VV_d V;  // an orthonormal basis of the Krylov subspace like {v0,A.v0,A^2.v0,...}
   V_d w;
   int i, j;
   ArnoldiProcess(const VV_d &A, const V_d &v0IN /*the first direction*/, const int nIN)
       : n(nIN),
         beta(Norm(v0IN)),
         v0(Normalize(v0IN)),
         H(VV_d(nIN + 1, V_d(nIN, 0.))),
         V(VV_d(nIN + 1, v0 /*V[0]=v0であればいい．ここではv0=v1=v2=..としている*/)),
         w(A.size()) {
      for (j = 0; j < n /*展開項数*/; ++j) {
         w = Dot(A, V[j]);  //@ 行列-ベクトル積
         for (i = 0; i < j + 1; ++i)
            w -= (H[i][j] = Dot(V[i], w)) * V[i];  //@ ベクトル内積
         V[j + 1] = w / (H[j + 1][j] = Norm(w));
         // MatrixForm(H);
      }
      // Print("done");
   };
};

struct gmres : public ArnoldiProcess {
   int n;  // th number of interation
   V_d e1, x, y;
   double err;
   /* NOTE:
   r0 = Normalize(b - A.x0)
   to find {r0,A.r0,A^2.r0,...}
   Therefore V is an orthonormal basis of the Krylov subspace Km(A,r0)
   */
   gmres(const VV_d &A, const V_d &b, const V_d &x0, const int nIN)
       : ArnoldiProcess(A, b - Dot(A, x0) /*行列-ベクトル積*/, nIN), n(nIN), e1(V_d(nIN + 1, 0.)) {
      e1[0] = ArnoldiProcess::beta;
      QR qr(ArnoldiProcess::H);
      // Print("Q");
      // MatrixForm(qr.Q);
      // Print("H");
      // MatrixForm(ArnoldiProcess::H);
      int i = 0;
      V_d g(qr.Q.size());
      for (const auto &q : qr.Q)
         g[i++] = q[0] * beta;
      // Print("予想される誤差", err = *g.rbegin());
      err = g[i];
      // std::cout << "予想される誤差" << err << std::endl;
      // std::cout << "n = " << nIN << std::endl;
      // std::cout << "g = " << g << std::endl;
      g.pop_back();
      y = back_substitution(qr.R, g, g.size());
      this->x = x0;
      for (auto i = 0; i < n; i++)
         this->x += y[i] * ArnoldiProcess::V[i];
   };

   ~gmres(){
       // delete ap;
   };
};
#endif