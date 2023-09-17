#ifndef basic_linear_systems_H
#define basic_linear_systems_H

#include <concepts>
#include <execution>
#include <numeric>  // for std::transform_reduce
//
#include "basic_IO.hpp"
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

/* -------------------------------------------------------------------------- */
/*                                   LAPACK                                   */
/* -------------------------------------------------------------------------- */

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
   const int nrhs = 1, LDB, LDA;
   int info;
   std::vector<int> ipiv;
   ~lapack_lu(){};

   lapack_lu(const std::vector<std::vector<double>> &aIN) : a(this->flatten(aIN)), dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim) {
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
   };

   template <size_t N>
   lapack_lu(const std::array<std::array<double, N>, N> &aIN) : a(this->flatten(aIN)), dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim) {
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
   };

   lapack_lu(const std::vector<std::vector<double>> &aIN, const std::vector<double> &rhd, std::vector<double> &b)
       : a(this->flatten(aIN)), dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim) {
      b = rhd;
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
      dgetrs_("T", &dim, &nrhs, a.data(), &LDA, ipiv.data(), b.data(), &LDB, &info);
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

   template <typename Container>
   std::vector<typename Container::value_type::value_type> flatten(const Container &mat) {
      using ValueType = typename Container::value_type::value_type;
      std::vector<ValueType> flattened;
      flattened.reserve(mat.size() * mat.size());
      for (const auto &part : mat)
         flattened.insert(flattened.end(), part.begin(), part.end());
      return flattened;
   }

   void solve(const std::vector<double> &rhd, std::vector<double> &b) {
      if ((dim != rhd.size()) || (dim != b.size()) || (rhd.size() != b.size()))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "dim != rhd.size() || dim != b.size() || rhd.size() != b.size() ");
      b = rhd;
      dgetrs_("T", &dim, &nrhs, a.data(), &LDA, ipiv.data(), b.data(), &LDB, &info);
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

   template <size_t N>
   void solve(const std::array<double, N> &rhd, std::array<double, N> &b) {
      b = rhd;
      dgetrs_("T", &dim, &nrhs, a.data(), &LDA, ipiv.data(), b.data(), &LDB, &info);
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

#include <algorithm>

extern "C" void dgesvd_(const char *jobu, const char *jobvt,
                        const int *m, const int *n,
                        double *a, const int *lda,
                        double *s,
                        double *u, const int *ldu,
                        double *vt, const int *ldvt,
                        double *work, const int *lwork, int *info);

struct lapack_svd {
   std::vector<double> a;
   const int m, n;
   std::vector<double> s, u, vt, work;
   int lwork;  // remove const
   int info;

   lapack_svd(const std::vector<std::vector<double>> &aIN)
       : a(flatten(aIN)), m(aIN.size()), n(aIN[0].size()), s(std::min(m, n)), u(m * m), vt(n * n) {
      char jobu = 'A', jobvt = 'A';
      double work_query;
      int lwork_query = -1;
      dgesvd_(&jobu, &jobvt, &m, &n, a.data(), &m,
              s.data(), u.data(), &m, vt.data(), &n,
              &work_query, &lwork_query, &info);
      lwork = static_cast<int>(work_query);
      work.resize(lwork);

      dgesvd_(&jobu, &jobvt, &m, &n, a.data(), &m,
              s.data(), u.data(), &m, vt.data(), &n,
              work.data(), &lwork, &info);
      if (info != 0) {
         throw std::runtime_error("Error in SVD computation");
      }
   };

   template <typename Container>
   std::vector<typename Container::value_type::value_type> flatten(const Container &mat) {
      using ValueType = typename Container::value_type::value_type;
      std::vector<ValueType> flattened;
      flattened.reserve(mat.size() * mat[0].size());
      for (const auto &part : mat)
         flattened.insert(flattened.end(), part.begin(), part.end());
      return flattened;
   }

   void solve(const std::vector<double> &b, std::vector<double> &x) {
      if (m != b.size())
         throw std::runtime_error("dimension mismatch");

      // Compute S_inv * Ut * b
      std::vector<double> tmp(m);
      for (int i = 0; i < m; ++i) {
         double inv_s = (s[i] > 1e-9) ? (1 / s[i]) : 0.0;
         for (int j = 0; j < m; ++j) {
            tmp[j] += u[j * m + i] * inv_s * b[j];
         }
      }

      // Compute V * (S_inv * Ut * b)
      x.resize(n);
      for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j) {
            x[i] += vt[j * n + i] * tmp[j];
         }
      }
   }
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

// V_d forward_substitution(const VV_d &mat, const V_d &b) {
//    int s = b.size();
//    V_d x(s);
//    double tmp = 0;
//    // for (auto i = 0; i < s; ++i)
//    // {
//    // 	tmp = 0;
//    // 	for (auto j = 0; j < i; ++j)
//    // 		tmp += mat[i][j] * x[j];
//    // 	x[i] = (b[i] - tmp) / mat[i][i];
//    // }

//    int i = 0, j;
//    for (const auto &a : mat) {
//       tmp = 0;
//       for (j = 0; j < i; ++j)
//          tmp += a[j] * x[j];
//       x[i] = (b[i] - tmp) / a[i];
//       i++;
//    }
//    return x;
// };

// // faster version
// V_d forward_substitution(const VV_d &mat, V_d b /*copy*/) {
//    double tmp = 0;
//    int i = 0, j;
//    for (const auto &a : mat) {
//       tmp = 0;
//       for (j = 0; j < i; ++j)
//          tmp += a[j] * b[j];
//       b[i] = (b[i] - tmp) / a[i];
//       i++;
//    }
//    return b;
// };

// V_d back_substitution(const VV_d &mat, V_d b, const Tii &mat_size) {
//    auto [row, col] = mat_size;
//    int i, j;
//    for (i = row - 1; i >= 0; --i) { /*　0~row-1まで，長さは，row　*/
//       for (j = col - 1; j > i; --j) /*　長さは，col-i-1　*/
//          b[i] -= mat[i][j] * b[j];
//       b[i] /= mat[i][i];
//    }
//    b.erase(std::next(b.begin(), row + 1), b.end());
//    return b;
// };

V_d forward_substitution(const VV_d &mat, V_d b /*copy*/) {
   int i = 0;
   for (const auto &a : mat) {
      double tmp = 0;
      for (int j = 0; j < i; ++j) {
         tmp = std::fma(a[j], b[j], tmp);
      }
      b[i] = std::fma(tmp, -1.0 / a[i], b[i]);
      i++;
   }
   return b;
};

V_d back_substitution(const VV_d &mat, V_d b, const Tii &mat_size) {
   auto [row, col] = mat_size;
   for (int i = row - 1; i >= 0; --i) {
      for (int j = col - 1; j > i; --j) {
         b[i] = std::fma(-mat[i][j], b[j], b[i]);
      }
      b[i] /= mat[i][i];
   }
   b.erase(std::next(b.begin(), row + 1), b.end());
   return b;
};

V_d back_substitution(const VV_d &mat, V_d b, const int mat_size) {
   return back_substitution(mat, b, Tii{mat_size, mat_size});
};

/* -------------------------------------------------------------------------- */
/*                              QR decomposition                              */
/* -------------------------------------------------------------------------- */

// struct QR {
//    VV_d Q, R, A;

//    // Copy constructor
//    QR(const QR &other)
//        : Q(other.Q),
//          R(other.R),
//          A(other.A) {
//       // No need to repeat computation; copy the computed Q, R, and A
//    }

//    ~QR() {
//       std::cout << "QR destructor" << std::endl;
//    };

//    QR(const VV_d &AIN) : R(AIN), A(AIN), Q(AIN.size(), V_d(AIN.size(), 0.)) { Initialize(AIN, true); };

//    void Initialize(const VV_d &AIN, const bool constractor = false) {
//       if (!constractor) {
//          A = R = AIN;
//          Q.resize(AIN.size(), V_d(AIN.size(), 0.));
//       }
//       IdentityMatrix(Q);
//       // MatrixForm(A, std::setw(10));
//       int n = AIN.size();
//       int m = AIN[0].size();
//       double r, c, s, a, b;
//       double Q0, Q1, R0, R1;  // Rのためのtmp
//       double eps = 1e-12;
//       for (auto j = 0; j < m; ++j) {
//          // for (auto i = n - 2; i >= j; --i) {  // givensの位置
//          for (auto i = j; i <= n - 2; ++i) {
//             if (!Between(R[i + 1][j], {-eps, eps})) {  // givensの位置
//                //! {i+1,j}がゼロとしたい成分

//                a = R[i][j];
//                b = R[i + 1][j];
//                // s*a+c*b=0
//                // -> s*a=-c*b
//                // -> tan(theta)=-b/a
//                // ----------------------------
//                // r = atan(-b / a);
//                // s = sin(r);  //! F_{i+k,i}
//                // c = cos(r);  //! F_{i+k,i+k}
//                // ----------------------------
//                r = std::sqrt(a * a + b * b);
//                s = -b / r;
//                c = a / r;
//                // ----------------------------
//                // std::cout << "R{" << i << "," << j << "} = " << a << std::endl;
//                // std::cout << "R{" << i + 1 << "," << j << "} = " << b << std::endl;
//                if (s != s || c != c) {
//                   s = 0.;
//                   c = 1.;
//                }

//                // Dot(F,R)の省略版
//                for (auto k = 0; k < m; ++k) {          // # column direction
//                   R0 = c * R[i][k] - s * R[i + 1][k];  //$ row i
//                   R1 = s * R[i][k] + c * R[i + 1][k];  //$ row i+k
//                   R[i][k] = R0;
//                   R[i + 1][k] = R1;
//                };
//                for (auto k = 0; k < n; ++k) {          // # column direction
//                   Q0 = c * Q[i][k] - s * Q[i + 1][k];  //$ row i
//                   Q1 = s * Q[i][k] + c * Q[i + 1][k];  //$ row i+k
//                   Q[i][k] = Q0;
//                   Q[i + 1][k] = Q1;
//                };
//                // Print("Q");
//                // MatrixForm(Q, std::setw(10));
//                // Print("R");
//                // MatrixForm(R, std::setw(15));
//             }
//          }
//       }
//    };

//    void IdentityMatrix(VV_d &mat) {
//       size_t i = 0;
//       for (auto &m : mat) {
//          m.assign(m.size(), 0.0);
//          m[i++] = 1.0;
//       }
//    }
// };

/* -------------------------------------------------------------------------- */

struct QR {
   VV_d Q, R, A, QT;

   // Copy constructor
   // No need to repeat computation; copy the computed Q, R, and A
   QR(const QR &other) : Q(other.Q), QT(other.QT), R(other.R), A(other.A) {}
   //  std::cout << "QR destructor" << std::endl;
   ~QR(){};
   QR(const VV_d &AIN) : R(AIN), A(AIN), Q(AIN.size(), V_d(AIN.size(), 0.)) { Initialize(AIN, true); };

   void Initialize(const VV_d &AIN, const bool constractor = false) {
      int N_ROW = AIN.size();
      int N_COL = AIN[0].size();
      int nR = AIN.size();
      int mR = AIN[0].size();
      if (!constractor) {
         A = R = AIN;
         Q.resize(AIN.size(), V_d(AIN.size(), 0.));
      }
      IdentityMatrix(Q);
      QT = Q;
      double r, c, s, a, b;
      double Q0, Q1, R0, R1;  // Rのためのtmp
      double eps = 1e-15;
      int max = std::max(N_ROW, N_COL);
      for (auto j = 0; j < max; ++j) {
         for (auto i = nR - 2; i >= j; --i) {
            // 下から
            // if (!Between(R[i + 1][j], {-eps, eps}))
            if (R[i + 1][j] != 0.) {  // givensの位置
               a = R[i][j];
               b = R[i + 1][j];  //! {i+1,j}がゼロとしたい成分
               r = std::sqrt(a * a + b * b);
               s = -b / r;
               c = a / r;
               if (s != s || c != c) {
                  s = 0.;
                  c = 1.;
               }

               //
               // A = F_1^T F_2^T ... F_2 F_1 A = QR
               // R = ... F_2 F_1 A
               // Q = F_1^T F_2^T...
               //
               // Rの更新
               for (auto col = 0; col < mR; ++col) {
                  R0 = c * R[i][col] - s * R[i + 1][col];  //$ row i
                  R1 = s * R[i][col] + c * R[i + 1][col];  //$ row i+k
                  R[i][col] = R0;
                  R[i + 1][col] = R1;
               };

               for (auto col = 0; col < QT[i].size(); ++col) {
                  Q0 = c * QT[i][col] - s * QT[i + 1][col];  //$ row i
                  Q1 = s * QT[i][col] + c * QT[i + 1][col];  //$ row i+k
                  QT[i][col] = Q0;
                  QT[i + 1][col] = Q1;
               };

               // Qの更新 * 左から書けるのでRとは逆順
               for (auto row = 0; row < Q.size(); ++row) {  // # row direction
                  Q0 = c * Q[row][i] - s * Q[row][i + 1];   //$ row i
                  Q1 = s * Q[row][i] + c * Q[row][i + 1];   //$ row i+k
                  Q[row][i] = Q0;
                  Q[row][i + 1] = Q1;
               };

               // std::cout << "-----------------" << std::endl;
               // std::cout << "i: " << i << ", j: " << j << std::endl;
               // std::cout << "-----------------" << std::endl;
               // Print("Q");
               // std::cout << MatrixForm(Q, std::setw(10)) << std::endl;
               // std::cout << "-----------------" << std::endl;
               // Print("R");
               // std::cout << "{i+1,j} = " << i + 1 << "," << j << std::endl;
               // auto func = [&](auto i_in, auto j_in) { return (i + 1 == i_in && j_in == j); };
               // std::cout << MatrixForm(R, func, 5, 20) << std::endl;
            }
         }
         // ここで，j列目の下半分はゼロになっているはず
      }
      // ここで，Rは上三角行列になっているはず
   };

   void IdentityMatrix(VV_d &mat) {
      size_t i = 0;
      for (auto &m : mat) {
         m.assign(m.size(), 0.0);
         m[i++] = 1.0;
      }
   }
};

/* -------------------------------------------------------------------------- */
struct CSR {
   std::unordered_map<CSR *, double> column_value;
   void clearColumnValue() {
      this->column_value.clear();
      this->canUseVector = false;
   };
   double value;
   double diagonal_value;
   double tmp_value;
   bool canUseVector;
   std::array<double, 3> value3d;
   size_t __index__;
   void setIndexCSR(size_t i) { this->__index__ = i; };
   size_t getIndexCSR() const { return __index__; };
   CSR() : canUseVector(false){};
   void clear() { this->column_value.clear(); }
   double at(CSR *const p) const { return column_value.at(p); };
   bool contains(CSR *const p) const { return column_value.contains(p); };
   void increment(CSR *const p, const double v) {
      auto [it, inserted] = this->column_value.insert({p, v});
      if (!inserted)
         it->second += v;
      this->canUseVector = false;
   };

   // 高速化のために，vectorに変換する．
   std::vector<std::tuple<int, double>> column_value_vector;
   void setVectorCSR() {
      column_value_vector.clear();
      column_value_vector.reserve(column_value.size());
      for (const auto &[crs, value] : column_value)
         column_value_vector.push_back({crs->__index__, value});
      this->canUseVector = true;
   };
};

template <typename T>
   requires std::derived_from<T, CSR>
double Dot(const std::unordered_map<T *, double> &column_value, const V_d &V) {
   double ret = 0.;
   for (const auto &[crs, value] : column_value)
      ret += value * V[crs->__index__];
   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
double Dot(const std::unordered_map<T *, double> &column_value, const std::unordered_set<T *> &V_crs) {
   double ret = 0.;

   for (const auto &[crs, value] : column_value)
      if (V_crs.contains(crs))
         ret += crs->value * value;

   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
V_d Dot(const std::unordered_set<T *> &V_crs) {
   V_d ret(V_crs.size(), 0.);
#pragma omp parallel
   for (const auto &crs0 : V_crs)
#pragma omp single nowait
   {
      double tmp = 0.;
      for (const auto &[crs1, value] : crs0->column_value)
         tmp += crs1->value * value;
      ret[crs0->__index__] = tmp;
   }
   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
V_d Dot(const std::unordered_set<T *> &A, const V_d &V) {
   V_d ret = V;
   // \label{CSR:parrallel}
#pragma omp parallel
   for (const auto &crs : A)
#pragma omp single nowait
   {
      double tmp = 0.;
      if (crs->canUseVector) {
         for (const auto &[i, value] : crs->column_value_vector)
            tmp += value * V[i];
      } else {
         for (const auto &[crs_local, value] : crs->column_value)
            tmp += value * V[crs_local->__index__];
      }
      ret[crs->__index__] = tmp;
   }
   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
V_d Dot(const std::vector<T *> &A, const V_d &V) {
   V_d ret = V;
   // \label{CSR:parrallel}
#pragma omp parallel
   for (const auto &crs : A)
#pragma omp single nowait
   {
      double tmp = 0.;
      if (crs->canUseVector) {
         for (const auto &[i, value] : crs->column_value_vector)
            tmp += value * V[i];
      } else {
         for (const auto &[crs_local, value] : crs->column_value)
            tmp += value * V[crs_local->__index__];
      }
      ret[crs->__index__] = tmp;
   }
   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
void DotOutput(const std::unordered_set<T *> &A, const V_d &V, V_d &w) {
   // \label{CSR:parrallel}
#pragma omp parallel
   for (const auto &crs : A)
#pragma omp single nowait
   {
      double tmp = 0.;
      if (crs->canUseVector) {
         for (const auto &[i, value] : crs->column_value_vector)
            tmp += value * V[i];
      } else {
         for (const auto &[crs_local, value] : crs->column_value)
            tmp += value * V[crs_local->__index__];
      }
      w[crs->__index__] = tmp;
   }
};

template <typename T>
   requires std::derived_from<T, CSR>
void DotOutput(const std::vector<T *> &A, const V_d &V, V_d &w) {
   // \label{CSR:parrallel}
#pragma omp parallel
   for (const auto &crs : A)
#pragma omp single nowait
   {
      double tmp = 0.;
      if (crs->canUseVector) {
         for (const auto &[i, value] : crs->column_value_vector)
            tmp += value * V[i];
      } else {
         for (const auto &[crs_local, value] : crs->column_value)
            tmp += value * V[crs_local->__index__];
      }
      w[crs->__index__] = tmp;
   }
};

/* -------------------------------------------------------------------------- */
struct ILU {
   std::vector<CSR *> LU;

   ILU(const std::unordered_set<CSR *> &A) {
      LU.reserve(A.size());

      for (const auto &a : A) {
         LU.push_back(new CSR(*a));  // assuming CSR has a copy constructor
      }

      size_t n = LU.size();

      for (size_t k = 0; k < n; ++k) {
         if (!LU[k]->contains(LU[k])) {
            throw std::runtime_error("Zero diagonal element in ILU factorization.");
         }
         double diag = LU[k]->at(LU[k]);

         for (auto &[i, lij] : LU[k]->column_value) {
            if (i->getIndexCSR() > k) {
               lij /= diag;

               for (auto &[j, ljk] : LU[i->getIndexCSR()]->column_value) {
                  if (j->getIndexCSR() > i->getIndexCSR()) {
                     LU[k]->increment(j, -lij * ljk);
                  }
               }
            }
         }
      }
   }

   ~ILU() {
      for (auto lu : LU) {
         delete lu;
      }
   }

   V_d solve(const V_d &b) {
      size_t n = LU.size();
      V_d x(n), y(n);

      // Forward solve Ly = b
      for (size_t i = 0; i < n; ++i) {
         y[i] = b[i];
         for (auto &[j, lij] : LU[i]->column_value) {
            if (j->getIndexCSR() < i) {
               y[i] -= lij * y[j->getIndexCSR()];
            }
         }
      }

      // Backward solve Ux = y
      for (int i = n - 1; i >= 0; --i) {
         x[i] = y[i];
         for (auto &[j, uij] : LU[i]->column_value) {
            if (j->getIndexCSR() > i) {
               x[i] -= uij * x[j->getIndexCSR()];
            }
         }
         x[i] /= LU[i]->at(LU[i]);
      }

      return x;
   }
};

/* ------------------------------------------------------ */
/*                          GMRES                         */
/* ------------------------------------------------------ */
/*DOC_EXTRACT ArnoldiProcess

## Arnoldi過程

1. 正規化した$`{\bf v}_1`$を与えておく．
2. $`{\bf v}_2 = {\rm Normalize}(\,\,\,\quad\quad\quad\quad\quad A{\bf v}_1 - ((A{\bf v}_1) \cdot {\bf v}_1){\bf v}_1\,\,\qquad\qquad\qquad\qquad\qquad\qquad)`$を計算する．
3. $`{\bf v}_3 = {\rm Normalize}(\quad\quad\quad({\bf w}=A{\bf v}_2 - ((A{\bf v}_2) \cdot {\bf v}_1){\bf v}_1)) - ({\bf w} \cdot {\bf v}_2){\bf v}_2\quad\quad\quad\quad\quad\quad)`$を計算する．
4. $`{\bf v}_4 = {\rm Normalize}(({\bf w}=(({\bf w}=A{\bf v}_3 - ((A{\bf v}_3) \cdot {\bf v}_1){\bf v}_1)) - ({\bf w} \cdot {\bf v}_2){\bf v}_2)) - ({\bf w} \cdot {\bf v}_3){\bf v}_3)`$を計算する．

言い換えると，

1. 正規化した$`{\bf v}_1`$を与えておく．
2. $`{\bf w}=A{\bf v}_1, {\bf v}_2 = {\rm Normalize}({\rm Chop}({\bf w},{\bf v}_1))`$を計算する．
3. $`{\bf w}=A{\bf v}_2, {\bf v}_3 = {\rm Normalize}({\rm Chop}({\rm Chop}({\bf w}, {\bf v}_1), {\bf v}_2))`$を計算する．
4. $`{\bf w}=A{\bf v}_3, {\bf v}_4 = {\rm Normalize}({\rm Chop}({\rm Chop}({\rm Chop}({\bf w}, {\bf v}_1), {\bf v}_2), {\bf v}_3))`$を計算する．

NOTE: ここで最も計算コストがかかるのは，$`{\bf w}=A{\bf v}_i`$の行列-ベクトル積である．

$`A{\bf v}_i`$の直交化の際に，
それに含まれる各基底$`{\bf v}_0,{\bf v}_1,...,{\bf v}_i`$の成分を計算している．
この成分からなる行列が，Hessenberg行列$`H`$である．

```math
\begin{align*}
A{\bf v}_1 & = h_{1,1} {\bf v}_1 + h_{2,1} {\bf v}_2\\
A{\bf v}_2 & = h_{1,2} {\bf v}_1 + h_{2,2} {\bf v}_2 + h_{3,2} {\bf v}_3\\
& \dots\\
A{\bf v}_{n} & = h_{1,n} {\bf v}_1 + h_{2,n} {\bf v}_2 + \cdots + h_{n,n+1} {\bf v}_{n+1}
\end{align*}
```

行列を使ってまとめると，

```math
A V_n = V_{n+1} \tilde H_n, \quad V_n = [v_1|v_2|...|v_n],
\quad \tilde H_n = \begin{bmatrix} h_{1,1} & h_{1,2} & \cdots & h_{1,n} & h_{1,n+1} \\ h_{2,1} & h_{2,2} & \cdots & h_{2,n} & h_{2,n+1} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & h_{n,n} & h_{n,n+1} \\ 0 & 0 & \cdots & 0 & h_{n+1,n+1} \end{bmatrix}
```

これをArnoldi分解という．ここで，$`[v_1|v_2|...|v_n]`$の$`|`$は列ベクトルを連結して行列を形成することを示している．

*/

// #define DEBUG_GMRES

template <typename Matrix>
struct ArnoldiProcess {

   int n;  // the number of interation
   double beta;
   V_d v0;
   VV_d H;  // ((n+1) x n) Hessenberg matrix
   VV_d V;  // ((n+1) x n) an orthonormal basis of the Krylov subspace like {v0,A.v0,A^2.v0,...}
   V_d w;
   ~ArnoldiProcess() { std::cout << "destucting ArnoldiProcess" << std::endl; };
   ArnoldiProcess(const Matrix &A, const V_d &v0IN /*the first direction*/, const int nIN)
       : n(nIN),
         beta(Norm(v0IN)),
         v0(v0IN / beta),
         H(VV_d(nIN + 1, V_d(nIN, 0.))),
         V(VV_d(nIN + 1, v0 /*V[0]=v0であればいい．ここではv0=v1=v2=..としている*/)),
         w(A.size()) {
      Initialize(A, v0IN, nIN, false);
   };

   void Initialize(const Matrix &A, const V_d &v0IN /*the first direction*/, const int nIN, const bool do_constract = true) {
      if (do_constract) {
         n = nIN;
         beta = Norm(v0IN);
         v0 = (v0IN / beta);
         H = (VV_d(nIN + 1, V_d(nIN, 0.)));
         V = (VV_d(nIN + 1, v0 /*V[0]=v0であればいい．ここではv0=v1=v2=..としている*/));
         w.resize(A.size());
      }
#if defined(DEBUG_GMRES)
      TimeWatch watch;
      std::cout << "ArnoldiProcess::Initialize" << std::endl;
#endif
#if defined(DEBUG_GMRES)
      TimeWatch watch;
      std::array<double, 2> tmp;
      double mean_elapsed_time_for_AV = 0., mean_elapsed_time_for_orthogonalization = 0., mean_elapsed_time_for_normalizing_w = 0.;
      int mean_elapsed_time_for_AV_count = 0, mean_elapsed_time_for_orthogonalization_count = 0, mean_elapsed_time_for_normalizing_w_count = 0;
      std::cout << "ArnoldiProcess" << std::endl;
#endif
      size_t i, j;
      for (j = 0; j < n /*展開項数*/; ++j) {
         DotOutput(A, V[j], w);  //@ 行列-ベクトル積\label{ArnoldiProcess:matrix-vector}
                                 // w = Dot(A, V[j]);  //@ 行列-ベクトル積\label{ArnoldiProcess:matrix-vector}
#if defined(DEBUG_GMRES)
         // std::cout << "Elapsed time for Dot(A, V[j]) " << tmp = watch() << std::endl;
         tmp = watch();
         mean_elapsed_time_for_AV += tmp[0];
         mean_elapsed_time_for_AV_count++;
#endif
         // orthogonalization
         for (i = 0; i <= j; ++i)
            w -= (H[i][j] = Dot(V[i], w)) * V[i];
#if defined(DEBUG_GMRES)
         // std::cout << "Elapsed time for orthogonalization " << watch() << std::endl;
         tmp = watch();
         mean_elapsed_time_for_orthogonalization += tmp[0];
         mean_elapsed_time_for_orthogonalization_count++;
#endif
         // normalize w
         V[j + 1] = w / (H[j + 1][j] = Norm(w));
#if defined(DEBUG_GMRES)
         // std::cout << "Elapsed time for normalizing w " << watch() << std::endl;
         tmp = watch();
         mean_elapsed_time_for_normalizing_w += tmp[0];
         mean_elapsed_time_for_normalizing_w_count++;
#endif
      }
#if defined(DEBUG_GMRES)
      std::cout << "Elapsed time" << watch() << std::endl;
      std::cout << "mean_elapsed_time_for_AV = " << mean_elapsed_time_for_AV / mean_elapsed_time_for_AV_count << std::endl;
      std::cout << "mean_elapsed_time_for_orthogonalization = " << mean_elapsed_time_for_orthogonalization / mean_elapsed_time_for_orthogonalization_count << std::endl;
      std::cout << "mean_elapsed_time_for_normalizing_w = " << mean_elapsed_time_for_normalizing_w / mean_elapsed_time_for_normalizing_w_count << std::endl;

#endif
#if defined(DEBUG_GMRES)
      std::cout << "Elapsed time" << watch() << std::endl;
#endif
   };

   void AddBasisVector(const Matrix &A) {
      // std::cout << "AddBasisVector" << std::endl;
      size_t i;
      // update n, V, H
      this->n++;
      H.push_back(V_d(n, 0.));
      w = Dot(A, V[n - 1]);  //@ 行列-ベクトル積\label{ArnoldiProcess:matrix-vector}
      // orthogonalization
      // std::cout << "orthogonalization" << std::endl;
      for (i = 0; i < n; ++i) {
         H[i].push_back(0.);
         w -= (H[i][n - 1] = Dot(V[i], w)) * V[i];
      }
      V.push_back(w / (H[n][n - 1] = Norm(w)));
      // std::cout << "done" << std::endl;
   };
};

/*DOC_EXTRACT GMRES

## 一般化最小残差法/GMRES

残差$`\|{\bf b} - A{\bf x}\|`$を最小とするような$`{\bf x}`$を求めたい．
そのような$`{\bf x}`$を，クリロフ部分空間の正規直交基底を用いた，$`{\bf x}_n = V_n {\bf y}_n`$の形で近似し，追い求めていく．
$`n`$はこの表現での展開項数である．$`V_n = \{{\bf v}_1,{\bf v}_2,...,{\bf v}_n\}`$は，アーノルディ過程によって計算する，クリロフ部分空間の正規直交基底である．

1. クリロフ部分空間法の考えから，$`\|{\bf b} - A V_n {\bf y}_n\|`$を最小とするような，$`{\bf y}_n`$を求める問題に書き換える．
2. $`A V_n = V_{n+1} \tilde H_n`$（アーノルディ分解）と書き換える．
3. $`V_{n+1}`$でくくる．
4. QR分解を使って，$`{\bf y}_n`$に関する最小二乗問題を$`{\bf y}_n`$について解く．

```math
\begin{align*}
\|{\bf b} - A{\bf x}_n\| & = \|{\bf b} - A V_n {\bf y}_n\|\\
& = \|{\bf b} - V_{n+1} \tilde H_n {\bf y}_n\|\quad \text{(use Arnoldi decomposition)}\\
& = \|V_{n+1} (\|{\bf b}\| {\bf e}_1 - \tilde H_n {\bf y}_n)\|\\
& = \|\|{\bf b}\| {\bf e}_1 - \tilde H_n {\bf y}_n\|\quad \text{(the dimension has been reduced!)}\\
& = \|\|{\bf b}\| {\bf e}_1 - QR {\bf y}_n\|\quad \text{(use QR decomposition)}\\
\end{align*}
```

<details>
<summary>なぜアーノルディ分解をするのか？</summary>

* $`A`$は$`m \times m`$とすると
* $`{\bf x}`$と$`{\bf b}`$は，$`m \times 1`$ベクトル（列ベクトル）.
* $`V_n`$は，$`m \times n`$行列で，$`A`$のクリロフ部分空間の基底ベクトルを列に持つ行列．
* $`{\bf y}_n`$は$`n \times 1`$ベクトル．
* $`\tilde H_n`$は$`(n+1) \times n`$行列．

従って，$`n`$が$`m`$よりも大幅に小さい場合，
アーノルディ分解によって作られた問題$`\min\|{\bf b} - V_{n+1}{\tilde H}_n {\bf y}_n\|`$は，
元の問題$`\min\|{\bf b}-A{\bf x}\|`$より計算量が少ない問題となる．

$`A{\bf x} = {\bf b}`$の問題を解くよりも，
$`{\tilde H}_n {\bf y}_n = {\bf b}`$という問題を解く方が計算量が少ない．

</details>

NOTE: アーノルディ過程が逐次的に計算できるため，展開項数$`n`$を$`n+1`$へと大きくしようとする際に（精度が$`n`$では十分でない場合），GMRESで近似解$`{\bf x}_{n+1}`$を始めから計算しなおす必要はない．$`V_{n+1}`$と$`{\tilde H}_{n+1}`$は，$`V_n`$と$`{\tilde H}_n`$を再利用するようにして計算でき，従って，比較的安く，得られている$`{\bf x}_n`$から$`{\bf x}_{n+1}`$へと更新できる．

*/

template <typename Matrix>
struct gmres : public ArnoldiProcess<Matrix> {
   /* NOTE:
   r0 = Normalize(b - A.x0)
   to find {r0,A.r0,A^2.r0,...}
   Therefore V is an orthonormal basis of the Krylov subspace Km(A,r0)
   */
   // int n;  // th number of interation
   V_d x, y;
   double err;
   QR qr;
   V_d g;
   ~gmres() { std::cout << "destucting gmres" << std::endl; };
   gmres(const Matrix &A, const V_d &b, const V_d &x0, const int nIN)
       : ArnoldiProcess<Matrix>(A, b - Dot(A, x0) /*行列-ベクトル積*/, nIN),
         // n(nIN),
         x(x0),
         y(b.size()),
         qr(ArnoldiProcess<Matrix>::H),
         g(qr.Q.size()),
         err(0.) {
#if defined(DEBUG_GMRES)
      TimeWatch watch;
      std::cout << "gmres" << std::endl;
#endif
      if (this->beta /*initial error*/ == static_cast<double>(0.))
         return;
      //
      size_t i = 0;
      for (const auto &q : qr.QT)
         g[i++] = q[0] * this->beta;

      err = g.back();  // 予想される誤差
      g.pop_back();
      this->y = back_substitution(qr.R, g, g.size());
      i = 0;
      for (i = 0; i < this->n; ++i)
         this->x += this->y[i] * this->V[i];
#if defined(DEBUG_GMRES)
      std::cout << "Elapsed time" << watch() << std::endl;
#endif
   };

   void Restart(const Matrix &A, const V_d &b, const V_d &x0, const int nIN) {
      // std::cout << "Restart" << std::endl;
      this->Initialize(A, b - Dot(A, x0), nIN);
      this->x = x0;
      // this->y.resize(b.size());
      this->qr.Initialize(this->H);
      this->g.resize(this->qr.Q.size());
      if (this->beta /*initial error*/ == static_cast<double>(0.))
         return;
      for (size_t i = 0; const auto &q : qr.QT)
         g[i++] = q[0] * this->beta;

      err = g.back();  // 予想される誤差
      g.pop_back();
      this->y = back_substitution(qr.R, g, g.size());
      for (size_t i = 0; i < this->n; ++i)
         this->x += this->y[i] * this->V[i];
      std::cout << "done" << std::endl;
   }

   void Iterate(const Matrix &A) {
      // std::cout << "Iterate" << std::endl;
      // do not change v0
      this->AddBasisVector(A);
      this->qr = QR(this->H);
      g.resize(qr.Q.size());
      err = 0.;
      if (this->beta /*initial error*/ == static_cast<double>(0.))
         return;
      // std::cout << "g.size() = " << g.size() << std::endl;
      size_t i = 0;
      for (const auto &q : this->qr.QT)
         g[i++] = q[0] * this->beta;
      err = g.back();  // 予想される誤差
      g.pop_back();
      // std::cout << "back_substitution" << std::endl;
      // this->y = back_substitution(this->qr.R, g, g.size());
      // i = this->n - 1;
      // this->x += this->y[i] * this->V[i];

      std::fill(this->x.begin(), this->x.end(), 0);
      this->y = back_substitution(this->qr.R, g, g.size());
      for (i = 0; i < this->n; ++i)
         this->x += this->y[i] * this->V[i];
      // std::cout << "done" << std::endl;
   };
};

/* -------------------------------------------------------------------------- */
V_d Eigenvalues(const VV_d &A, const double tol = 1e-9, const int maxIter = 1000) {
   VV_d Ak = A, I = A;
   IdentityMatrix(I);
   V_d eigenvalues;
   for (int i = 0; i < maxIter; ++i) {
      QR qr(Ak);
      eigenvalues = Diagonal(Ak = Dot(qr.R, qr.Q));
      if (std::ranges::all_of(eigenvalues, [&](const auto lambda) { return std::abs(Det(A - lambda * I)) < tol; })) {
         return eigenvalues;
      }
   }
   return eigenvalues;
}

std::pair<V_d, VV_d> Eigensystem(const VV_d &A, const double tol = 1e-9, const int maxIter = 1000) {
   VV_d Ak = A, I = A, Qk = A;
   IdentityMatrix(Qk);
   IdentityMatrix(I);
   V_d eigenvalues;
   for (int i = 0; i < maxIter; ++i) {
      QR qr(Ak);
      Ak = Dot(qr.R, qr.Q);        // A(k+1) = R * Q
      eigenvalues = Diagonal(Ak);  // Extract diagonal as eigenvalues
      Qk = Dot(Qk, qr.Q);          // Update the transformation matrix
      if (std::ranges::all_of(eigenvalues, [&](const auto lambda) { return std::abs(Det(A - lambda * I)) < tol; })) {
         return {eigenvalues, Qk};
      }
   }
   return {eigenvalues, Qk};  // Return both eigenvalues and eigenvectors
}

template <std::size_t N>
std::array<double, N> Eigenvalues(const std::array<std::array<double, N>, N> &A, const double tol = 1e-9, const int maxIter = 1000) {
   std::array<double, N> eigenvalues;
   int i = 0;
   for (const auto &a : Eigenvalues(ToVector(A), tol, maxIter))
      eigenvalues[i++] = a;
   return eigenvalues;
}

template <std::size_t N>
std::pair<std::array<double, N>, std::array<std::array<double, N>, N>> Eigensystem(const std::array<std::array<double, N>, N> &A, const double tol = 1e-9, const int maxIter = 1000) {
   auto [E, Qk] = Eigensystem(ToVector(A), tol, maxIter);
   std::array<std::array<double, N>, N> Q;
   std::array<double, N> eigenvalues;
   int i = 0;
   for (const auto &e : E)
      eigenvalues[i++] = e;

   i = 0;
   int j = 0;
   for (const auto &q : Qk) {
      for (const auto &qj : q)
         Q[i][j++] = qj;
      i++;
      j = 0;
   }
   // normalize vectors and modify eigenvalues
   for (auto i = 0; i < N; ++i) {
      double norm = Norm(Q[i]);
      Q[i] /= norm;
      eigenvalues[i] *= norm;
   }

   return {eigenvalues, Q};
}

#endif