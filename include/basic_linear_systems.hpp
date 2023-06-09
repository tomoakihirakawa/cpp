#ifndef basic_linear_systems_H
#define basic_linear_systems_H

#include <concepts>
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
      size_t i = 0;
      for (auto &m : mat) {
         m.assign(m.size(), 0.0);
         m[i++] = 1.0;
      }
   }
};

/* -------------------------------------------------------------------------- */
struct CSR {
   double value;
   double diagonal_value;
   double tmp_value;
   std::array<double, 3> value3d;
   size_t __index__;
   void setIndexCSR(size_t i) {
      this->__index__ = i;
   };
   size_t getIndexCSR() const { return __index__; };
   std::unordered_map<CSR *, double> column_value;
   CSR(){};
   void clear() { this->column_value.clear(); }
   double at(CSR *const p) const { return column_value.at(p); };
   bool contains(CSR *const p) const { return column_value.contains(p); };
   void increment(CSR *const p, const double v) {
      auto [it, inserted] = this->column_value.insert({p, v});
      if (!inserted)
         it->second += v;
   };
};
template <typename T>
   requires std::derived_from<T, CSR>
double Dot(const std::unordered_map<T *, double> &column_value, const V_d &V) {
   double ret = 0.;
   for (const auto &[crs, value] : column_value) {
      ret += value * V[crs->__index__];
   }
   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
double Dot(const std::unordered_map<T *, double> &column_value, const std::unordered_set<T *> &V_crs) {
   double ret = 0.;

   if (V_crs.size() <= column_value.size()) {
      for (const auto &crs : V_crs) {
         auto it = column_value.find(crs);
         if (it != column_value.end()) {
            ret += crs->value * it->second;
         }
      }
   } else {
      for (const auto &[crs, value] : column_value) {
         if (V_crs.find(crs) != V_crs.end()) {
            ret += value * crs->value;
         }
      }
   }

   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
V_d Dot(const std::unordered_set<T *> &V_crs) {
   V_d ret(V_crs.size(), 0.);
   for (const auto &crs0 : V_crs) {
      auto &tmp = ret[crs0->__index__];
      for (const auto &[crs1, value] : crs0->column_value) {
         tmp += value * crs1->value;
      }
   }
   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
V_d Dot(const std::unordered_set<T *> &A, const V_d &V) {
   V_d ret(V.size(), 0.);
   // for (const auto &crs : A) {
   //    ret[crs->__index__] += Dot(crs->column_value, V);
   // }

   // \label{CSR:parrallel}
   for (const auto &crs : A)
      crs->tmp_value = 0;
#pragma omp parallel
   for (const auto &crs : A)
#pragma omp single nowait
      crs->tmp_value += Dot(crs->column_value, V);

#pragma omp parallel
   for (const auto &crs : A)
#pragma omp single nowait
      ret[crs->__index__] = crs->tmp_value;

   return ret;
};

template <typename T>
   requires std::derived_from<T, CSR>
V_d Dot(const std::unordered_set<T *> &A, const std::unordered_set<T *> &V_crs) {
   V_d ret(V_crs.size(), 0.);
   // for (const auto &crs : A) {
   //    ret[crs->__index__] += Dot(crs->column_value, V_crs);
   // }

   // \label{CSR:parrallel}
   for (const auto &crs : A)
      crs->tmp_value = 0;

#pragma omp parallel
   for (const auto &crs : A)
#pragma omp single nowait
      crs->tmp_value += Dot(crs->column_value, V_crs);

#pragma omp parallel
   for (const auto &crs : A)
#pragma omp single nowait
      ret[crs->__index__] = crs->tmp_value;

   return ret;
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

## ArnoldiProcess

ヘッセンベルグ行列$H[0:k-1]$は，Aと相似なベクトルであり，同じ固有値を持つ
GMRESで使う場合，$V0$にはNormalize(b-A.x0)を与える．
x0は初期値

アーノルディ法は固有値問題の数値解法であり反復解法．
一般的な行列の固有ベクトルと固有値をクリロフ空間の直行基底によって近似する方法計算する方法．

   1. 正規化した${\bf v}_0$を与えておく．
   2. $\quad\quad\quad\quad\quad{\bf v}_1 = {\rm Normalize}(A{\bf v}_0 - ((A{\bf v}_0) \cdot {\bf v}_0){\bf v}_0)$を計算する．
   3. $\quad\quad\quad{\bf v}_2 = {\rm Normalize}((w=A{\bf v}_1 - ((A{\bf v}_1) \cdot {\bf v}_0){\bf v}_0)) - (w \cdot {\bf v}_1){\bf v}_1)$を計算する．
   4. ${\bf v}_3 = {\rm Normalize}((w=((w=A{\bf v}_2 - ((A{\bf v}_2) \cdot {\bf v}_0){\bf v}_0)) - (w \cdot {\bf v}_1){\bf v}_1)) - (w \cdot {\bf v}_2){\bf v}_2)$を計算する．

*/
template <typename Matrix>
struct ArnoldiProcess {

   int n;  // the number of interation
   double beta;
   V_d v0;
   VV_d H;  // ((n+1) x n) Hessenberg matrix
   VV_d V;  // ((n+1) x n) an orthonormal basis of the Krylov subspace like {v0,A.v0,A^2.v0,...}
   V_d w;
   ArnoldiProcess(const Matrix &A, const V_d &v0IN /*the first direction*/, const int nIN)
       : n(nIN),
         beta(Norm(v0IN)),
         v0(v0IN / beta),
         H(VV_d(nIN + 1, V_d(nIN, 0.))),
         V(VV_d(nIN + 1, v0 /*V[0]=v0であればいい．ここではv0=v1=v2=..としている*/)),
         w(A.size()) {
      Print("ArnoldiProcess");
      size_t i, j;
      for (j = 0; j < n /*展開項数*/; ++j) {
         w = Dot(A, V[j]);  //@ 行列-ベクトル積\label{ArnoldiProcess:matrix-vector}

         // orthogonalization
         for (i = 0; i <= j; ++i)
            w -= (H[i][j] = Dot(V[i], w)) * V[i];
         V[j + 1] = w / (H[j + 1][j] = Norm(w));
      }
      Print("done");
   };
};

template <typename Matrix>
struct gmres : public ArnoldiProcess<Matrix> {
   /* NOTE:
   r0 = Normalize(b - A.x0)
   to find {r0,A.r0,A^2.r0,...}
   Therefore V is an orthonormal basis of the Krylov subspace Km(A,r0)
   */
   int n;  // th number of interation
   V_d x, y;
   double err;
   const QR qr;
   V_d g;
   ~gmres(){};
   gmres(const Matrix &A, const V_d &b, const V_d &x0, const int nIN)
       : ArnoldiProcess<Matrix>(A, b - Dot(A, x0) /*行列-ベクトル積*/, nIN),
         n(nIN),
         x(x0),
         y(b.size()),
         qr(ArnoldiProcess<Matrix>::H),
         g(qr.Q.size()),
         err(0.) {
      if (ArnoldiProcess<Matrix>::beta /*initial error*/ == static_cast<double>(0.))
         return;
      //
      size_t i = 0;
      for (const auto &q : qr.Q)
         g[i++] = q[0] * ArnoldiProcess<Matrix>::beta;

      err = g.back();  // 予想される誤差
      g.pop_back();
      this->y = back_substitution(qr.R, g, g.size());
      i = 0;
      for (i = 0; i < n; ++i)
         this->x += this->y[i] * ArnoldiProcess<Matrix>::V[i];
   };
   // gmres(const Matrix &A, const V_d &b, const V_d &x0, const int nIN, const auto ILU)
   //     : ArnoldiProcess<Matrix>(A, ILU(A).solve(b - Dot(A, x0)) /*行列-ベクトル積*/, nIN),
   //       n(nIN),
   //       x(x0),
   //       y(b.size()),
   //       qr(ArnoldiProcess<Matrix>::H),
   //       g(qr.Q.size()),
   //       err(0.) {
   //    if (ArnoldiProcess<Matrix>::beta /*initial error*/ == static_cast<double>(0.))
   //       return;
   //    //
   //    size_t i = 0;
   //    for (const auto &q : qr.Q)
   //       g[i++] = q[0] * ArnoldiProcess<Matrix>::beta;

   //    err = g.back();  // 予想される誤差
   //    g.pop_back();
   //    this->y = back_substitution(qr.R, g, g.size());
   //    i = 0;
   //    for (i = 0; i < n; ++i)
   //       this->x += this->y[i] * ArnoldiProcess<Matrix>::V[i];
   // };
};

#endif