#ifndef basic_linear_systems_H
#define basic_linear_systems_H

#include "basic_arithmetic_vector_operations.hpp"
#include "basic_exception.hpp"

struct ludcmp_parallel {
   int n;
   VV_d lu;
   V_i indx;
   double d;
   V_d vv;
   ludcmp_parallel(const VV_d &a) : n(a.size()), lu(a), indx(a.size()), vv(a.size(), 1.) {
      const double TINY = 1.0e-40;
      int i, imax = 0;
      double big, temp;
      d = 1.0;
      for (auto i = 0; i < n; ++i) {
         big = 0.0;
         for (const auto &lu_i_j : lu[i])
            if ((temp = std::abs(lu_i_j)) > big)
               big = temp;
         if (big == 0.0)
            throw("Singular matrix in LUdcmp");
         vv[i] /= big;
      }
      /* ------------------------------------------------------ */
      int k = 0;
      for (auto &lu_k : lu) {
         big = 0.0;
         for (i = k; i < n; ++i)
            if ((temp = vv[i] * std::abs(lu[i][k])) > big) {
               big = temp;
               imax = i;
            }
         if (k != imax) {
            lu_k.swap(lu[imax]);
            d = -d;
            vv[imax] = vv[k];
         }
         indx[k] = imax;
         if (lu_k[k] == 0.0)
            lu_k[k] = TINY;

#ifdef _OPENMP
#pragma omp parallel for
#endif
         for (auto i = k + 1; i < n; ++i) {
            lu[i][k] /= lu_k[k];
            for (auto j = k + 1; j < n; ++j)
               lu[i][j] -= lu[i][k] * lu_k[j];
         }

         // #ifdef _OPENMP
         // #pragma omp parallel for
         // #endif
         // 			for (auto i = k + 1; i < n; ++i)
         // 			{
         // 				auto temp = lu[i][k] /= lu_k[k];
         // 				for (auto j = k + 1; j < n; ++j)
         // 					lu[i][j] -= temp * lu_k[j];
         // 			}

         k++;
      }
   };

   void solve(const V_d &b, V_d &x) {
      int i, ii = 0, ip, j;
      double sum;
      if (b.size() != n || x.size() != n)
         throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));

      // for (i=0;i<n;i++)
      //   x[i] = b[i];

      x = b;

      for (i = 0; i < n; i++) {
         ip = indx[i];
         sum = x[ip];
         x[ip] = x[i];
         if (ii != 0)
            for (j = ii - 1; j < i; j++)
               sum -= lu[i][j] * x[j];
         else if (sum != 0.0)
            ii = i + 1;
         x[i] = sum;
      }
      for (i = n - 1; i >= 0; i--) {
         sum = x[i];
         for (j = i + 1; j < n; j++)
            sum -= lu[i][j] * x[j];
         x[i] = sum / lu[i][i];
      }
   };

   void solve(const VV_d &b, VV_d &x) {
      int i, j, m = b[0].size();
      if (b.size() != n || x.size() != n || b[0].size() != x.size())
         throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
      V_d xx(n);
      for (j = 0; j < m; j++) {
         for (i = 0; i < n; i++)
            xx[i] = b[i][j];
         solve(xx, xx);
         for (i = 0; i < n; i++)
            x[i][j] = xx[i];
      }
   };

   void inverse(VV_d &ainv) {
      int i, j;
      ainv.resize(n, V_d(n, 0));
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++)
            ainv[i][j] = 0.;
         ainv[i][i] = 1.;
      }
      solve(ainv, ainv);
   };

   VV_d Inverse() {
      VV_d ainv(n, V_d(n, 0));
      int i, j;
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++)
            ainv[i][j] = 0.;
         ainv[i][i] = 1.;
      }
      solve(ainv, ainv);
      return ainv;
   };

   double det() {
      double dd = d;
      for (int i = 0; i < n; i++)
         dd *= lu[i][i];
      return dd;
   };

   // void mprove(V_d &b, V_d &x)
   // {
   // 	int i, j;
   // 	V_d r(n);
   // 	this->aref.resize(n);
   // 	for (i = 0; i < n; i++)
   // 	{
   // 		double sdp = -b[i];
   // 		for (j = 0; j < n; j++)
   // 			sdp += (double)aref[i][j] * (double)x[j];
   // 		r[i] = sdp;
   // 	}
   // 	solve(r, r);
   // 	for (i = 0; i < n; i++)
   // 		x[i] -= r[i];
   // };
};

// #ifdef use_lapack

extern "C" void dgetrf_(const int *dim1,
                        const int *dim2,
                        double *a,
                        const int *lda,
                        int *ipiv,
                        int *info);
extern "C" void dgetrs_(const char *TRANS,
                        const int *N,
                        const int *NRHS,
                        const double *A,
                        const int *LDA,
                        const int *IPIV,
                        double *B,
                        const int *LDB,
                        int *INFO);

struct lapack_lu {
   std::vector<double> a;
   const int dim;
   const char trans = 'T';
   const int nrhs = 1, LDA, LDB;
   int info;
   std::vector<int> ipiv;
   ~lapack_lu(){};
   lapack_lu(const std::vector<std::vector<double>> &aIN)
       : a(this->flatten(aIN)),
         dim(aIN.size()),
         LDB(dim),
         LDA(dim),
         ipiv(dim) {
      dgetrf_(&dim, &dim, &*a.begin(), &LDA, &*ipiv.begin(), &info);
   };

   lapack_lu(const std::vector<std::vector<double>> &aIN,
             const std::vector<double> &rhd,
             std::vector<double> &b)
       : a(this->flatten(aIN)),
         dim(aIN.size()),
         LDB(dim),
         LDA(dim),
         ipiv(dim) {
      b = rhd;
      dgetrf_(&dim, &dim, &*a.begin(), &LDA, &*ipiv.begin(), &info);
      dgetrs_(&trans, &dim, &nrhs, &*a.begin(), &LDA, &*ipiv.begin(), &*b.begin(), &LDB, &info);
      if (info) {
         // std::cerr << e.what() << colorOff << std::endl;
         std::stringstream ss;
         ss << "LDB:" << LDB;
         ss << "\nLDA:" << LDA;
         ss << "\nipiv:" << ipiv;
         // ss << "\na:" << a;
         ss << "\nb:" << b;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };

   std::vector<double> flatten(const std::vector<std::vector<double>> &mat) const {
      int i = 0;
      std::vector<double> ret(mat.size() * mat[0].size());
      for (const auto &part : mat)
         for (const auto &value : part)
            ret[i++] = value;
      return ret;
   };

   void solve(const std::vector<double> &rhd, std::vector<double> &b) {
      b = rhd;
      dgetrs_(&trans, &dim, &nrhs, &*a.begin(), &LDA, &*ipiv.begin(), &*b.begin(), &LDB, &info);
      if (info) {
         // std::cerr << e.what() << colorOff << std::endl;
         std::stringstream ss;
         ss << "LDB:" << LDB;
         ss << "\nLDA:" << LDA;
         ss << "\nipiv:" << ipiv;
         // ss << "\na:" << a;
         ss << "\nb:" << b;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };
};

// #endif

struct ludcmp {
   int n;
   VV_d lu;
   V_i indx;
   double d;
   V_d vv;
   ludcmp(const VV_d &a) : n(a.size()), lu(a), indx(a.size()), vv(a.size(), 1.) {
      const double TINY = 1.0e-40;
      int i, imax = 0, j;
      double big, temp;
      // V_d vv(n);
      d = 1.0;
      // #ifdef _OPENMP
      // #pragma omp parallel for
      // #endif
      for (auto i = 0; i < n; i++) {
         big = 0.0;
         for (const auto &lu_i_j : lu[i])
            if ((temp = std::abs(lu_i_j)) > big)
               big = temp;
         if (big == 0.0)
            throw("Singular matrix in LUdcmp");
         vv[i] /= big;
      }
      /* ------------------------------------------------------ */
      // for (auto k = 0; k < n; k++)
      // {
      int k = 0;
      for (auto &lu_k : lu) {
         // auto &lu_k = lu[k];
         big = 0.0;
         for (i = k; i < n; ++i) {
            temp = vv[i] * std::abs(lu[i][k]);
            if (temp > big) {
               big = temp;
               imax = i;
            }
         }
         if (k != imax) {
            // for (auto j = 0; j < n; ++j)
            // 	std::swap(lu[imax][j], lu_k[j]);
            lu_k.swap(lu[imax]);
            d = -d;
            vv[imax] = vv[k];
         }
         indx[k] = imax;
         if (lu_k[k] == 0.0)
            lu_k[k] = TINY;

         for (auto i = k + 1; i < n; ++i) {
            temp = lu[i][k] /= lu_k[k];
            for (auto j = k + 1; j < n; ++j)
               lu[i][j] -= temp * lu_k[j];
         }
         k++;
      }
   };

   void solve(const V_d &b, V_d &x) const {
      int i, ii = 0, ip, j;
      double sum;
      if (b.size() != n || x.size() != n)
         throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "サイズが違います"));

      x = b;

      for (i = 0; i < n; ++i) {
         ip = indx[i];
         sum = x[ip];
         x[ip] = x[i];
         if (ii != 0)
            for (j = ii - 1; j < i; j++)
               sum -= lu[i][j] * x[j];
         else if (sum != 0.0)
            ii = i + 1;
         x[i] = sum;
      }
      for (i = n - 1; i >= 0; i--) {
         sum = x[i];
         for (j = i + 1; j < n; j++)
            sum -= lu[i][j] * x[j];
         x[i] = sum / lu[i][i];
      }
   };

   void solve(const VV_d &b, VV_d &x) {
      int i, j, m = b[0].size();
      if (b.size() != n || x.size() != n || b[0].size() != x.size())
         throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
      V_d xx(n);
      for (j = 0; j < m; j++) {
         for (i = 0; i < n; i++)
            xx[i] = b[i][j];
         solve(xx, xx);
         for (i = 0; i < n; i++)
            x[i][j] = xx[i];
      }
   };

   void inverse(VV_d &ainv) {
      int i, j;
      ainv.resize(n, V_d(n, 0));
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++)
            ainv[i][j] = 0.;
         ainv[i][i] = 1.;
      }
      solve(ainv, ainv);
   };

   VV_d Inverse() {
      VV_d ainv(n, V_d(n, 0));
      int i, j;
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++)
            ainv[i][j] = 0.;
         ainv[i][i] = 1.;
      }
      solve(ainv, ainv);
      return ainv;
   };

   double det() {
      double dd = d;
      for (int i = 0; i < n; i++)
         dd *= lu[i][i];
      return dd;
   };

   // void mprove(V_d &b, V_d &x)
   // {
   // 	int i, j;
   // 	V_d r(n);
   // 	this->aref.resize(n);
   // 	for (i = 0; i < n; i++)
   // 	{
   // 		double sdp = -b[i];
   // 		for (j = 0; j < n; j++)
   // 			sdp += (double)aref[i][j] * (double)x[j];
   // 		r[i] = sdp;
   // 	}
   // 	solve(r, r);
   // 	for (i = 0; i < n; i++)
   // 		x[i] -= r[i];
   // };
};

VV_d Inverse(const VV_d &mat) {
   ludcmp LU(mat);
   VV_d ret;
   LU.inverse(ret);
   return ret;
};

// LAPACKにもQRがあるが．
#include <iomanip>
#include "basic_vectors.hpp"
struct QR {
   VV_d Q, R, A;
   QR(const VV_d &AIN)
       : R(AIN),
         A(AIN),
         Q(AIN.size(), V_d(AIN.size(), 0.)) {
      IdentityMatrix(Q);
      // MatrixForm(A, std::setw(10));
      int n = AIN.size();
      int m = AIN[0].size();
      double r, c, s, a, b;
      double Q0, Q1, R0, R1;  // Rのためのtmp
      double eps = 1e-12;
      for (auto j = 0; j < m; ++j) {
         // for (auto i = n - 2; i >= j; --i) {  // givensの位置
         for (auto i = j; i <= n - 2; ++i) {
            if (!Between(R[i + 1][j], {-eps, eps})) {  // givensの位置
               //! {i+1,j}がゼロとしたい成分

               a = R[i][j];
               b = R[i + 1][j];
               // s*a+c*b=0
               // -> s*a=-c*b
               // -> tan(theta)=-b/a
               // ----------------------------
               // r = atan(-b / a);
               // s = sin(r);  //! F_{i+k,i}
               // c = cos(r);  //! F_{i+k,i+k}
               // ----------------------------
               r = std::sqrt(a * a + b * b);
               s = -b / r;
               c = a / r;
               // ----------------------------
               // std::cout << "R{" << i << "," << j << "} = " << a << std::endl;
               // std::cout << "R{" << i + 1 << "," << j << "} = " << b << std::endl;
               if (s != s || c != c) {
                  s = 0.;
                  c = 1.;
               }

               // Dot(F,R)の省略版
               for (auto k = 0; k < m; ++k) {          // # column direction
                  R0 = c * R[i][k] - s * R[i + 1][k];  //$ row i
                  R1 = s * R[i][k] + c * R[i + 1][k];  //$ row i+k
                  R[i][k] = R0;
                  R[i + 1][k] = R1;
               };
               for (auto k = 0; k < n; ++k) {          // # column direction
                  Q0 = c * Q[i][k] - s * Q[i + 1][k];  //$ row i
                  Q1 = s * Q[i][k] + c * Q[i + 1][k];  //$ row i+k
                  Q[i][k] = Q0;
                  Q[i + 1][k] = Q1;
               };
               // Print("Q");
               // MatrixForm(Q, std::setw(10));
               // Print("R");
               // MatrixForm(R, std::setw(15));
            }
         }
      }
   };
   void IdentityMatrix(VV_d &mat) {
      int i = 0, j = 0;
      for (auto &m : mat) {
         j = 0;
         for (auto &n : m) {
            if (i == j)
               n = 1.;
            else
               n = 0.;
            j++;
         }
         i++;
      }
   };
};

#endif