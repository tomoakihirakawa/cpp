#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>
#include "basic_vectors.hpp"  // あなたの lapack_svd / lapack_lu を含む

using Mat = std::vector<std::vector<double>>;
using Vec = std::vector<double>;

// --- 小ユーティリティ ---
static Mat zeros(int m, int n) { return Mat(m, std::vector<double>(n, 0.0)); }
static Mat matmul(const Mat& A, const Mat& B) {
   int m = A.size(), k = A[0].size(), n = B[0].size();
   Mat C = zeros(m, n);
   for (int i = 0; i < m; ++i) {
      for (int s = 0; s < k; ++s) {
         double a = A[i][s];
         if (a == 0.0) continue;
         for (int j = 0; j < n; ++j) C[i][j] += a * B[s][j];
      }
   }
   return C;
}
static Mat transpose(const Mat& A) {
   int m = A.size(), n = A[0].size();
   Mat At = zeros(n, m);
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) At[j][i] = A[i][j];
   return At;
}
static Mat subtract(const Mat& A, const Mat& B) {
   int m = A.size(), n = A[0].size();
   Mat C = zeros(m, n);
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) C[i][j] = A[i][j] - B[i][j];
   return C;
}
static double frob_norm(const Mat& A) {
   double s = 0.0;
   for (auto& r : A)
      for (double v : r) s += v * v;
   return std::sqrt(s);
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
static Mat diag_from_topK(const std::vector<double>& S, int k) {
   Mat D = zeros(k, k);
   for (int i = 0; i < k && i < (int)S.size(); ++i) D[i][i] = S[i];
   return D;
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

// --- 99%到達の打ち切りランクを決める ---
// 既定: エネルギ基準 sum σ^2（Frobenius）。sum σ にしたい場合は USE_ENERGY=false。
static int rank_by_cumulative(const std::vector<double>& S, double target = 0.99, bool USE_ENERGY = true) {
   if (S.empty()) return 0;
   std::vector<double> w(S.size());
   if (USE_ENERGY) {
      for (size_t i = 0; i < S.size(); ++i) w[i] = S[i] * S[i];
   } else {
      for (size_t i = 0; i < S.size(); ++i) w[i] = S[i];
   }
   double total = std::accumulate(w.begin(), w.end(), 0.0);
   if (total <= 0.0) return 0;

   double acc = 0.0;
   for (size_t i = 0; i < S.size(); ++i) {
      acc += w[i];
      if (acc / total >= target) return (int)(i + 1);
   }
   return (int)S.size();
}

// 要素ごとの最大絶対誤差・最大相対誤差（相対はAの要素基準）
static void max_errors(const Mat& A, const Mat& B, double& max_abs, double& max_rel) {
   max_abs = 0.0;
   max_rel = 0.0;
   const double eps = 1e-15;
   int m = A.size(), n = A[0].size();
   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         double diff = std::abs(A[i][j] - B[i][j]);
         max_abs = std::max(max_abs, diff);
         double denom = std::max(eps, std::abs(A[i][j]));
         max_rel = std::max(max_rel, diff / denom);
      }
   }
}

int main() {
   std::cout << std::fixed;

   // =========================================
   // 例2：情報圧縮（99%エネルギ到達で打ち切り ＋ 再構築度の評価）
   // =========================================
   {
      std::cout << "\n==== 例2：情報圧縮（99%基準 + 再構築度） ====\n";
      const int m = 8, n = 8;
      Mat B(m, std::vector<double>(n));
      for (int i = 0; i < m; ++i)
         for (int j = 0; j < n; ++j)
            B[i][j] = 2.0 + 0.6 * i + 0.4 * j + 0.05 * std::sin(0.4 * i * j);

      lapack_svd svdB(B);
      const Mat& U = svdB.U;                  // m×m
      const std::vector<double>& S = svdB.S;  // σ降順
      const Mat& Vt = svdB.VT;                // n×n

      const double percent = 0.999;
      // 99% 到達ランク（エネルギ=Σσ^2 基準）
      int k = rank_by_cumulative(S, percent, /*USE_ENERGY=*/true);
      std::cout << "selected rank k = " << k << " (to reach " << percent * 100 << "% energy)\n";

      // 低ランク近似 B_k = U_k Σ_k V_k^T
      Mat Uk = crop_cols(U, k);              // m×k
      Mat Vtk = crop_rows(Vt, k);            // k×n
      Mat Sk = diag_from_topK(S, k);         // k×k
      Mat Bk = matmul(matmul(Uk, Sk), Vtk);  // m×n

      // --- 再構築の評価 ---
      // 絶対誤差（Frobenius）
      double errF = frob_norm(subtract(B, Bk));
      // 相対誤差（Frobenius）
      double relErrF = errF / std::max(1e-15, frob_norm(B));
      // 再構築率（保持エネルギー率） = Σ_{i<=k} σ_i^2 / Σ σ_i^2
      double sumS2 = 0.0, sumS2k = 0.0;
      for (double v : S) sumS2 += v * v;
      for (int i = 0; i < k && i < (int)S.size(); ++i) sumS2k += S[i] * S[i];
      double retention = sumS2k / std::max(1e-15, sumS2);

      std::cout << "||B - B_k||_F (abs) = " << errF << "\n";
      std::cout << "||B - B_k||_F / ||B||_F (relative) = " << relErrF << "\n";
      std::cout << "reconstruction rate (energy kept) = " << retention << "\n";

      // 参考：σの和で 99% にしたい場合（“強度”基準）
      int k_sum = rank_by_cumulative(S, percent, /*USE_ENERGY=*/false);
      std::cout << "(sum σ 基準) rank = " << k_sum << "\n";

      // --- スペクトルの可視化と理論検証を追加 ---
      auto print_spectrum = [](const std::vector<double>& S) {
         std::cout << "singular values (σ): ";
         for (double v : S) std::cout << std::setprecision(6) << v << " ";
         std::cout << "\n";
      };
      print_spectrum(S);

      // 累積エネルギー表（上位iまでの Σσ^2 / Σσ^2）
      double totalS2 = 0.0;
      for (double v : S) totalS2 += v * v;
      std::cout << "cumulative energy (Σσ^2 ratio):\n  i   ratio\n";
      double accS2 = 0.0;
      for (int i = 0; i < (int)S.size(); ++i) {
         accS2 += S[i] * S[i];
         std::cout << "  " << (i + 1) << "   " << (accS2 / std::max(1e-15, totalS2)) << "\n";
      }

      // Eckart–Young の等式の数値検証： ||B-Bk||_F^2 vs. Σ_{i>k} σ_i^2
      double tailS2 = 0.0;
      for (int i = k; i < (int)S.size(); ++i) tailS2 += S[i] * S[i];
      double errF_sq = errF * errF;
      std::cout << "check: ||B-B_k||_F^2 = " << errF_sq
                << " , sum_{i>k} σ_i^2 = " << tailS2
                << " , diff = " << std::abs(errF_sq - tailS2) << "\n";

      // 行列の中身を表示
      print_mat("Original B", B, /*maxr=*/m, /*maxc=*/n);
      print_mat("Reconstructed B_k", Bk, /*maxr=*/m, /*maxc=*/n);

      // 差分を表示
      Mat Diff = subtract(B, Bk);
      print_mat("Difference (B - B_k)", Diff, /*maxr=*/m, /*maxc=*/n);

      // 要素ごとの最大誤差（目安）
      double max_abs = 0.0, max_rel = 0.0;
      max_errors(B, Bk, max_abs, max_rel);
      std::cout << "max |B - B_k| = " << max_abs << "\n";
      std::cout << "max relative error = " << max_rel << "\n";
   }

   return 0;
}