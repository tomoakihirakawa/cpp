#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>
#include "basic_vectors.hpp"  // あなたの lapack_svd / lapack_lu を含む想定

// ↓ U(), S(), Vt() を提供しているなら有効化
// #define LAPACK_SVD_HAS_USV

using Mat = std::vector<std::vector<double>>;
using Vec = std::vector<double>;

// ---------- 小ユーティリティ ----------
static Mat zeros(int m, int n) { return Mat(m, std::vector<double>(n, 0.0)); }
static Mat eye(int n) {
   Mat I = zeros(n, n);
   for (int i = 0; i < n; ++i) I[i][i] = 1.0;
   return I;
}
static Mat transpose(const Mat& A) {
   int m = A.size(), n = A[0].size();
   Mat At = zeros(n, m);
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) At[j][i] = A[i][j];
   return At;
}
static Mat matmul(const Mat& A, const Mat& B) {
   int m = A.size(), k = A[0].size(), n = B[0].size();
   Mat C = zeros(m, n);
   for (int i = 0; i < m; ++i)
      for (int s = 0; s < k; ++s) {
         double a = A[i][s];
         if (a == 0.0) continue;
         for (int j = 0; j < n; ++j) C[i][j] += a * B[s][j];
      }
   return C;
}
static Vec matvec(const Mat& A, const Vec& x) {
   int m = A.size(), n = A[0].size();
   Vec y(m, 0.0);
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) y[i] += A[i][j] * x[j];
   return y;
}
static Mat subtract(const Mat& A, const Mat& B) {
   int m = A.size(), n = A[0].size();
   Mat C = zeros(m, n);
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) C[i][j] = A[i][j] - B[i][j];
   return C;
}
static Vec subtract(const Vec& a, const Vec& b) {
   Vec c(a.size());
   for (size_t i = 0; i < a.size(); ++i) c[i] = a[i] - b[i];
   return c;
}
static double frob_norm(const Mat& A) {
   double s = 0.0;
   for (auto& r : A)
      for (double v : r) s += v * v;
   return std::sqrt(s);
}
static double vec_norm2(const Vec& x) {
   double s = 0.0;
   for (double v : x) s += v * v;
   return std::sqrt(s);
}
static Mat diag_from(const Vec& s) {
   Mat D = zeros((int)s.size(), (int)s.size());
   for (int i = 0; i < (int)s.size(); ++i) D[i][i] = s[i];
   return D;
}
static Mat crop_cols(const Mat& A, int k) {
   int m = A.size(), n = A[0].size();
   k = std::min(k, n);
   Mat B = zeros(m, k);
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < k; ++j) B[i][j] = A[i][j];
   return B;
}
static Mat crop_rows(const Mat& A, int k) {
   int m = A.size(), n = A[0].size();
   k = std::min(k, m);
   Mat B = zeros(k, n);
   for (int i = 0; i < k; ++i)
      for (int j = 0; j < n; ++j) B[i][j] = A[i][j];
   return B;
}
static void print_mat(const std::string& name, const Mat& A, int maxr = 6, int maxc = 6) {
   std::cout << name << " =\n";
   int m = A.size(), n = A[0].size();
   for (int i = 0; i < std::min(m, maxr); ++i) {
      std::cout << "  ";
      for (int j = 0; j < std::min(n, maxc); ++j)
         std::cout << std::setw(10) << std::setprecision(6) << A[i][j] << " ";
      if (n > maxc) std::cout << "...";
      std::cout << "\n";
   }
   if (m > maxr) std::cout << "  ...\n";
}
static void print_vec(const std::string& name, const Vec& v) {
   std::cout << name << " = [ ";
   for (double x : v) std::cout << std::setprecision(6) << x << " ";
   std::cout << "]\n";
}
static double projector_error_idempotent(const Mat& P) {
   // ||P^2 - P||_F
   return frob_norm(subtract(matmul(P, P), P));
}
static double projector_error_symmetric(const Mat& P) {
   // ||P - P^T||_F
   return frob_norm(subtract(P, transpose(P)));
}

int main() {
   std::cout << std::fixed;

   // ======================================================
   // 例1：SVDの本質（A^+ を使った列空間/行空間への射影と安定解）
   // ======================================================
   {
      std::cout << "==== 例1：SVDの本質（射影と安定解） ====\n";
      Mat A = {
          {1.0, 2.0, 3.0},
          {2.0, 3.0, 4.0},
          {7.0, 8.0, 9.0}};
      Vec b = {1.0, 1.0, 1.0};
      print_mat("A", A);

      // あなたのSVDラッパ：solve は A^+ b を返す実装を想定
      lapack_svd svd(A);

      // 最小二乗/擬似逆解
      Vec x_svd((int)A[0].size(), 0.0);
      svd.solve(b, x_svd);
      print_vec("x_svd (A^+ b)", x_svd);
      std::cout << "||Ax - b||_2 = " << vec_norm2(subtract(matvec(A, x_svd), b)) << "\n";

      // 擬似逆 A^+ を取得（あなたの inverse() を使用）
      Mat Aplus = svd.inverse();  // A^+ （rank欠損時は小特異値を抑制したものでもOK）
      print_mat("A^+ (pseudo-inverse)", Aplus);

      // 列空間への直交射影 P_col = A A^+
      Mat Pcol = matmul(A, Aplus);
      // 行空間への直交射影 P_row = A^+ A
      Mat Prow = matmul(Aplus, A);

      std::cout << "列空間射影 P_col の性質: "
                << " idempotent ||P^2-P||_F=" << projector_error_idempotent(Pcol)
                << ", symmetric ||P-P^T||_F=" << projector_error_symmetric(Pcol) << "\n";
      std::cout << "行空間射影 P_row の性質: "
                << " idempotent ||P^2-P||_F=" << projector_error_idempotent(Prow)
                << ", symmetric ||P-P^T||_F=" << projector_error_symmetric(Prow) << "\n";

      // 射影の幾何：任意の y を列空間に直交射影
      Vec y = {2.0, -1.0, 0.5};
      Vec y_proj = matvec(Pcol, y);
      print_vec("y", y);
      print_vec("P_col * y (yのA列空間への直交射影)", y_proj);

      // 再構成：A ≈ (A A^+) A = A（数学的には等号．数値的に誤差のみ）
      Mat Arec = matmul(Pcol, A);
      std::cout << "reconstruction error ||A - (A A^+)A||_F = "
                << frob_norm(subtract(A, Arec)) << "\n";
   }
   return 0;
}