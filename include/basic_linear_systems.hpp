#ifndef basic_linear_systems_H
#define basic_linear_systems_H

#include <cmath>
#include "basic_vectors.hpp"
// basic_linear_systems2.hpp

//! A.x = b
void Solve(const std::array<std::array<double, 2>, 2> &A, std::array<double, 2> &x, const std::array<double, 2> &y) {
   const double inv_det = 1. / std::fma(A[0][0], A[1][1], -A[0][1] * A[1][0]);
   x[0] = std::fma(A[1][1], y[0], -A[0][1] * y[1]) * inv_det;
   x[1] = std::fma(-A[1][0], y[0], A[0][0] * y[1]) * inv_det;
};

//! x.A = b
void Solve(std::array<double, 2> &x, const std::array<std::array<double, 2>, 2> &A, const std::array<double, 2> &y) {
   const double inv_det = 1. / std::fma(A[0][0], A[1][1], -A[0][1] * A[1][0]);
   x[0] = std::fma(A[1][1], y[0], -A[1][0] * y[1]) * inv_det;
   x[1] = std::fma(-A[0][1], y[0], A[0][0] * y[1]) * inv_det;
};

// ! x.A = b
void Solve(std::array<double, 3> &x, const std::array<std::array<double, 3>, 3> &A, const std::array<double, 3> &y) {
   const double inv_det = 1. / std::fma(A[0][2], std::fma(A[1][1], A[2][0], -A[1][0] * A[2][1]), std::fma(A[0][1], std::fma(-A[1][2], A[2][0], A[1][0] * A[2][2]), A[0][0] * std::fma(A[1][2], A[2][1], -A[1][1] * A[2][2])));
   x[0] = inv_det * std::fma(y[2], A[1][1] * A[2][0], std::fma(-y[1], A[1][2] * A[2][0], std::fma(-y[2], A[1][0] * A[2][1], std::fma(y[0], A[1][2] * A[2][1], std::fma(y[1], A[1][0], -y[0] * A[1][1]) * A[2][2]))));
   x[1] = -inv_det * std::fma(y[2], A[0][1] * A[2][0], std::fma(-y[1], A[0][2] * A[2][0], std::fma(-y[2], A[0][0] * A[2][1], std::fma(y[0], A[0][2] * A[2][1], std::fma(y[1], A[0][0], -y[0] * A[0][1]) * A[2][2]))));
   x[2] = inv_det * std::fma(y[2], A[0][1] * A[1][0], std::fma(-y[1], A[0][2] * A[1][0], std::fma(-y[2], A[0][0] * A[1][1], std::fma(y[0], A[0][2] * A[1][1], std::fma(y[1], A[0][0], -y[0] * A[0][1]) * A[1][2]))));
};

// ! A.x = b
void Solve(std::array<std::array<double, 3>, 3> A, std::array<double, 3> &x, const std::array<double, 3> &y) {
   std::swap(A[0][1], A[1][0]);
   std::swap(A[0][2], A[2][0]);
   std::swap(A[1][2], A[2][1]);
   Solve(x, A, y);
};

#include <concepts>
#include <execution>
#include <numeric>  // for std::transform_reduce
//
#include "basic_IO.hpp"
#include "basic_arithmetic_vector_operations.hpp"
#include "basic_exception.hpp"
#include "lib_measurement.hpp"

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

extern "C" void dgetrf_(const int *dim1, const int *dim2, double *a, const int *lda, int *ipiv, int *info);
extern "C" void dgetrs_(const char *TRANS, const int *N, const int *NRHS, const double *A, const int *LDA, const int *IPIV, double *B, const int *LDB, int *INFO);
extern "C" void dgetri_(const int *n, double *a, const int *lda, const int *ipiv, double *work, int *lwork, int *info);

struct lapack_lu {
   const int dim, nrhs = 1, LDB, LDA;
   int info;
   std::vector<int> ipiv;
   std::vector<double> a;
   char TRANS = 'T';
   ~lapack_lu() {};

   lapack_lu(const std::vector<std::vector<double>> &aIN) : dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim), a(dim * dim) {
      std::size_t i, j, k = 0;
      for (i = 0; i < dim; ++i)
         for (j = 0; j < dim; ++j)
            a[k++] = aIN[i][j];
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
   };

   template <std::size_t N>
   lapack_lu(const std::array<std::array<double, N>, N> &aIN) : dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim), a(dim * dim) {
      std::size_t i, j, k = 0;
      for (i = 0; i < dim; ++i)
         for (j = 0; j < dim; ++j)
            a[k++] = aIN[i][j];
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
   };

   //!  solve at initialization
   //! solve A.x = b
   lapack_lu(const std::vector<std::vector<double>> &aIN, std::vector<double> &x, const std::vector<double> &rhd) : dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim), a(dim * dim) {
      std::size_t i, j, k = 0;
      for (i = 0; i < dim; ++i)
         for (j = 0; j < dim; ++j)
            a[k++] = aIN[i][j];
      x = rhd;
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
      dgetrs_(&TRANS, &dim, &nrhs, a.data(), &LDA, ipiv.data(), x.data(), &LDB, &info);
      if (info) {
         std::stringstream ss;
         ss << "LDB:" << LDB;
         ss << "\nLDA:" << LDA;
         ss << "\nipiv:" << ipiv;
         // ss << "\na:" << a;
         ss << "\nb:" << x;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };

   //! solve x.A = b
   lapack_lu(std::vector<double> &x, const std::vector<std::vector<double>> &aIN, const std::vector<double> &rhd) : dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim), a(dim * dim), TRANS('N') {
      std::size_t i, j, k = 0;
      for (i = 0; i < dim; ++i)
         for (j = 0; j < dim; ++j)
            a[k++] = aIN[i][j];
      x = rhd;
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
      dgetrs_(&TRANS, &dim, &nrhs, a.data(), &LDA, ipiv.data(), x.data(), &LDB, &info);
      if (info) {
         std::stringstream ss;
         ss << "LDB:" << LDB;
         ss << "\nLDA:" << LDA;
         ss << "\nipiv:" << ipiv;
         // ss << "\na:" << a;
         ss << "\nb:" << x;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };

   //! solve A.x = b
   template <std::size_t N>
   lapack_lu(const std::array<std::array<double, N>, N> &aIN, std::array<double, N> &x, const std::array<double, N> &rhd) : dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim), a(dim * dim) {
      std::size_t i, j, k = 0;
      for (i = 0; i < dim; ++i)
         for (j = 0; j < dim; ++j)
            a[k++] = aIN[i][j];
      x = rhd;
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
      dgetrs_(&TRANS, &dim, &nrhs, a.data(), &LDA, ipiv.data(), x.data(), &LDB, &info);
      if (info) {
         std::stringstream ss;
         ss << "LDB:" << LDB;
         ss << "\nLDA:" << LDA;
         ss << "\nipiv:" << ipiv;
         // ss << "\na:" << a;
         ss << "\nb:" << x;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };

   //! solve x.A = b
   template <std::size_t N>
   lapack_lu(std::array<double, N> &x, const std::array<std::array<double, N>, N> &aIN, const std::array<double, N> &rhd) : dim(aIN.size()), LDB(dim), LDA(dim), ipiv(dim), a(dim * dim), TRANS('N') {
      std::size_t i, j, k = 0;
      for (i = 0; i < dim; ++i)
         for (j = 0; j < dim; ++j)
            a[k++] = aIN[i][j];
      x = rhd;
      dgetrf_(&dim, &dim, a.data(), &LDA, ipiv.data(), &info);
      dgetrs_(&TRANS, &dim, &nrhs, a.data(), &LDA, ipiv.data(), x.data(), &LDB, &info);
      if (info) {
         std::stringstream ss;
         ss << "LDB:" << LDB;
         ss << "\nLDA:" << LDA;
         ss << "\nipiv:" << ipiv;
         // ss << "\na:" << a;
         ss << "\nb:" << x;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };

   //! solve A.x = b or x.A = b depending on the TRANS = 'T' or 'N'
   void solve(const std::vector<double> &rhd, std::vector<double> &x) {
      if ((dim != rhd.size()) || (dim != x.size()) || (rhd.size() != x.size()))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "dim != rhd.size() || dim != x.size() || rhd.size() != x.size() ");
      x = rhd;
      dgetrs_(&TRANS, &dim, &nrhs, a.data(), &LDA, ipiv.data(), x.data(), &LDB, &info);
      if (info) {
         std::stringstream ss;
         ss << "LDB:" << LDB;
         ss << "\nLDA:" << LDA;
         ss << "\nipiv:" << ipiv;
         // ss << "\na:" << a;
         ss << "\nb:" << x;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };

   template <std::size_t N>
   void solve(const std::array<double, N> &rhd, std::array<double, N> &x) {
      x = rhd;
      dgetrs_(&TRANS, &dim, &nrhs, a.data(), &LDA, ipiv.data(), x.data(), &LDB, &info);
      if (info) {
         std::stringstream ss;
         ss << "LDB:" << LDB;
         ss << "\nLDA:" << LDA;
         ss << "\nipiv:" << ipiv;
         // ss << "\na:" << a;
         ss << "\nb:" << x;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };

   /* --------------------------------------------- */

   template <typename ARRAY>
   void inverse_helper(ARRAY &ret) {
      if (dim * dim != a.size())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "dim * dim != a.size()");

      //* Query and allocate the optimal workspace
      int lwork = -1;
      double wkopt = 0.0;
      auto a_copy = a;
      dgetri_(&dim, a_copy.data(), &LDA, ipiv.data(), &wkopt, &lwork, &info);
      lwork = static_cast<int>(wkopt);
      std::vector<double> work(lwork);
      //* Compute the inverse

      dgetri_(&dim, a_copy.data(), &LDA, ipiv.data(), work.data(), &lwork, &info);
      if (info != 0) {
         std::stringstream ss;
         ss << "Error in dgetri_ with info = " << info;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }

      for (std::size_t i = 0; i < dim; ++i)
         for (std::size_t j = 0; j < dim; ++j)
            ret[i][j] = a_copy[i * dim + j];  //! do not Transpose back because LAPACK uses column-major order
                                              //! because this comput A^-1
   }

   std::vector<std::vector<double>> inverse() {
      std::vector<std::vector<double>> ret(dim, std::vector<double>(dim));
      inverse_helper(ret);
      return ret;
   }

   // std::vector<std::vector<double>> inverse() {
   //    //! Check if the matrix is square
   //    if (dim * dim != a.size())
   //       throw std::runtime_error("The matrix must be square to compute its inverse.");
   //    //! Create an identity matrix
   //    std::vector<double> id_mat(dim * dim, 0.);
   //    for (std::size_t i = 0; i < dim; ++i)
   //       id_mat[i * dim + i] = 1;
   //    // Use dgetrs_ to find each column of the inverse
   //    dgetrs_(&TRANS, &dim, &dim, a.data(), &LDA, ipiv.data(), id_mat.data(), &LDB, &info);
   //    // Check for errors
   //    if (info) {
   //       std::stringstream ss;
   //       ss << "Error in dgetrs_ with info = " << info;
   //       throw std::runtime_error(ss.str());
   //    }

   //    // Convert the flat array back to a 2D array
   //    std::vector<std::vector<double>> inv_mat(dim, std::vector<double>(dim));
   //    for (std::size_t i = 0; i < dim; ++i) {
   //       for (std::size_t j = 0; j < dim; ++j) {
   //          inv_mat[i][j] = id_mat[j * dim + i];
   //       }
   //    }
   //    return inv_mat;
   // }
};

template <std::size_t N>
std::array<std::array<double, N>, N> Inverse(const std::array<std::array<double, N>, N> &A) {
   lapack_lu lu(A);
   std::array<std::array<double, N>, N> inv_array;

   lu.inverse_helper(inv_array);
   return inv_array;

   // auto inv = lu.inverse();
   // for (std::size_t i = 0; i < N; ++i)
   //    for (std::size_t j = 0; j < N; ++j)
   //       inv_array[i][j] = inv[i][j];
   // return inv_array;
}

// for three dimensional array
// template <>
// std::array<std::array<double, 3>, 3> Inverse(const std::array<std::array<double, 3>, 3> &a) {
//    const double inv_det = 1. / std::fma(-a[0][2], a[1][1] * a[2][0], std::fma(a[0][1], a[1][2] * a[2][0], std::fma(a[0][2], a[1][0] * a[2][1], std::fma(-a[0][0], a[1][2] * a[2][1], std::fma(-a[0][1], a[1][0] * a[2][2], a[0][0] * a[1][1] * a[2][2])))));
//    return {{{std::fma(-a[1][2], a[2][1], a[1][1] * a[2][2]) * inv_det, std::fma(a[0][2], a[2][1], -a[0][1] * a[2][2]) * inv_det, std::fma(-a[0][2], a[1][1], a[0][1] * a[1][2]) * inv_det},
//             {std::fma(a[1][2], a[2][0], -a[1][0] * a[2][2]) * inv_det, std::fma(-a[0][2], a[2][0], a[0][0] * a[2][2]) * inv_det, std::fma(a[0][2], a[1][0], -a[0][0] * a[1][2]) * inv_det},
//             {std::fma(-a[1][1], a[2][0], a[1][0] * a[2][1]) * inv_det, std::fma(a[0][1], a[2][0], -a[0][0] * a[2][1]) * inv_det, std::fma(-a[0][1], a[1][0], a[0][0] * a[1][1]) * inv_det}}};
// }

template <>
std::array<std::array<double, 3>, 3> Inverse(const std::array<std::array<double, 3>, 3> &a) {
   // lapack_lu lu(a);
   // auto ret = lu.inverse();
   // return {{{ret[0][0], ret[0][1], ret[0][2]}, {ret[1][0], ret[1][1], ret[1][2]}, {ret[2][0], ret[2][1], ret[2][2]}}};

   const double inv_det = 1.0 / std::fma(
                                    -std::get<2>(std::get<0>(a)), std::get<1>(std::get<1>(a)) * std::get<0>(std::get<2>(a)),
                                    std::fma(
                                        std::get<1>(std::get<0>(a)), std::get<2>(std::get<1>(a)) * std::get<0>(std::get<2>(a)),
                                        std::fma(
                                            std::get<2>(std::get<0>(a)), std::get<0>(std::get<1>(a)) * std::get<1>(std::get<2>(a)),
                                            std::fma(
                                                -std::get<0>(std::get<0>(a)), std::get<2>(std::get<1>(a)) * std::get<1>(std::get<2>(a)),
                                                std::fma(
                                                    -std::get<1>(std::get<0>(a)), std::get<0>(std::get<1>(a)) * std::get<2>(std::get<2>(a)),
                                                    std::get<0>(std::get<0>(a)) * std::get<1>(std::get<1>(a)) * std::get<2>(std::get<2>(a)))))));

   return {{{std::fma(-std::get<2>(std::get<1>(a)), std::get<1>(std::get<2>(a)), std::get<1>(std::get<1>(a)) * std::get<2>(std::get<2>(a))) * inv_det,
             std::fma(std::get<2>(std::get<0>(a)), std::get<1>(std::get<2>(a)), -std::get<1>(std::get<0>(a)) * std::get<2>(std::get<2>(a))) * inv_det,
             std::fma(-std::get<2>(std::get<0>(a)), std::get<1>(std::get<1>(a)), std::get<1>(std::get<0>(a)) * std::get<2>(std::get<1>(a))) * inv_det},
            {std::fma(std::get<2>(std::get<1>(a)), std::get<0>(std::get<2>(a)), -std::get<0>(std::get<1>(a)) * std::get<2>(std::get<2>(a))) * inv_det,
             std::fma(-std::get<2>(std::get<0>(a)), std::get<0>(std::get<2>(a)), std::get<0>(std::get<0>(a)) * std::get<2>(std::get<2>(a))) * inv_det,
             std::fma(std::get<2>(std::get<0>(a)), std::get<0>(std::get<1>(a)), -std::get<0>(std::get<0>(a)) * std::get<2>(std::get<1>(a))) * inv_det},
            {std::fma(-std::get<1>(std::get<1>(a)), std::get<0>(std::get<2>(a)), std::get<0>(std::get<1>(a)) * std::get<1>(std::get<2>(a))) * inv_det,
             std::fma(std::get<1>(std::get<0>(a)), std::get<0>(std::get<2>(a)), -std::get<0>(std::get<0>(a)) * std::get<1>(std::get<2>(a))) * inv_det,
             std::fma(-std::get<1>(std::get<0>(a)), std::get<0>(std::get<1>(a)), std::get<0>(std::get<0>(a)) * std::get<1>(std::get<1>(a))) * inv_det}}};
}

VV_d Inverse(const VV_d &mat) {
   lapack_lu lu(mat);
   return lu.inverse();
};

// Solve Ax = b using lapack
// vector case
// void SolveLinearSystem(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x) {
//    lapack_lu(A, b, x);
// };
// void SolveLinearSystem(const std::vector<double> &b, const std::vector<std::vector<double>> &A, std::vector<double> &x) {
//    lapack_lu(A, b, x);
// };

// array case
// template <size_t N>
// void SolveLinearSystem(const std::array<std::array<double, N>, N> &A, std::array<double, N> &x, const std::array<double, N> &b) {
//    lapack_lu(A, b, x);
// };

// // Solve xA = b using lapack
// // vector case
// void SolveLinearSystem(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x) {
//    lapack_lu(b, A, x);
// };

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
   int lwork;
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
   }

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

      std::vector<double> tmp(std::min(m, n), 0.0);
      for (std::size_t i = 0; i < std::min(m, n); ++i) {
         double inv_s = (s[i] > 1e-9) ? (1 / s[i]) : 0.0;
         for (std::size_t j = 0; j < m; ++j) {
            tmp[i] += u[j * m + i] * b[j];
         }
         tmp[i] *= inv_s;
      }

      x.resize(n, 0.0);
      for (std::size_t i = 0; i < n; ++i) {
         for (std::size_t j = 0; j < std::min(m, n); ++j) {
            x[i] += vt[j * n + i] * tmp[j];
         }
      }
   }
};

std::array<double, 3> optimalVectorSVD(std::vector<double> Vsample,
                                       std::vector<Tddd> Directions,
                                       const std::vector<double> &weights) {
   if (Vsample.size() == 1)
      return Vsample[0] * Directions[0];

   double mean = 0.;
   for (const auto &v : Vsample)
      mean += std::abs(v);
   mean /= Vsample.size();
   if (mean == 0)
      return {0., 0., 0.};

   for (auto &d : Directions)
      d = Normalize(d);

   const double tolerance = 1E-12 * mean;
   Tddd Vinit = {0., 0., 0.};
   for (std::size_t i = 0; i < Vsample.size(); ++i)
      Vinit += Vsample[i] * Directions[i];
   Vinit /= Vsample.size();

   auto diff = [&](const Tddd &U, const std::size_t i) -> double { return Dot(U, Directions[i]) - Vsample[i]; };

   auto optimizing_function = [&](const Tddd &U) -> double {
      double S = 0;
      for (std::size_t i = 0; i < Vsample.size(); ++i) {
         // S += weights[i] * std::pow(diff(U, i), 2);
         S += weights[i] * std::pow(Dot(U, Directions[i]) - Vsample[i], 2);
      }
      return 0.5 * S;
   };

   if (optimizing_function(Vinit) < tolerance) {
      return Vinit;
   }

   std::vector<double> X(3);
   std::vector<std::vector<double>> A(Directions.size(), std::vector<double>(3));
   for (std::size_t i = 0; i < Directions.size(); ++i) {
      Vsample[i] *= weights[i];
      for (std::size_t j = 0; j < 3; ++j)
         A[i][j] = Directions[i][j] * weights[i];
   }

   lapack_svd svd(A);
   svd.solve(Vsample, X);

   return {X[0], X[1], X[2]};
}

std::array<double, 3> optimalVectorSVD(std::vector<double> Vsample,
                                       std::vector<Tddd> Directions) {
   return optimalVectorSVD(Vsample, Directions, std::vector<double>(Vsample.size(), 1.));
}

/* -------------------------------------------------------------------------- */
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
   std::size_t i = 0, j;
   double tmp;
   for (const auto &a : mat) {
      tmp = 0;
      for (j = 0; j < i; ++j) {
         tmp = std::fma(a[j], b[j], tmp);
      }
      b[i] = std::fma(tmp, -1.0 / a[i], b[i]);
      i++;
   }
   return b;
};

// V_d back_substitution(const VV_d &mat, V_d b, const Tii &mat_size) {
//    auto [row, col] = mat_size;
//    for (int i = row - 1; i >= 0; --i) {
//       for (int j = col - 1; j > i; --j) {
//          b[i] = std::fma(-mat[i][j], b[j], b[i]);
//       }
//       b[i] /= mat[i][i];
//    }
//    b.erase(std::next(b.begin(), row + 1), b.end());
//    return b;
// };

// V_d back_substitution(const VV_d &mat, V_d &b, const Tii &mat_size) {
//    const auto [row, col] = mat_size;
//    double bi;
//    int i, j;
//    for (i = row - 1; i >= 0; --i) {
//       auto &mat_i = mat[i];
//       bi = b[i];  // Cache b[i] for better cache locality
//       for (j = col - 1; j > i; --j)
//          bi = std::fma(-mat_i[j], b[j], bi);
//       b[i] = bi / mat_i[i];
//    }
//    b.resize(row);  // Resize vector to `row`, this will automatically remove extra elements
//    return b;
// };

V_d back_substitution(const VV_d &mat, V_d b, const std::size_t row, const std::size_t col) {
   double bi;
   int i, j;
   for (i = row - 1; i >= 0; --i) {
      auto &mat_i = mat[i];
      bi = b[i];  // Cache b[i] for better cache locality
      for (j = col - 1; j > i; --j)
         bi = std::fma(-mat_i[j], b[j], bi);
      b[i] = bi / mat_i[i];
   }
   b.resize(row);  // Resize vector to `row`, this will automatically remove extra elements
   return b;
};

// V_d back_substitution(const VV_d &mat, V_d b, const std::size_t row_col) {
//    const auto row = row_col;
//    const auto col = row_col;
//    // const auto [row, col] = mat_size;
//    double bi;
//    int i, j;
//    for (i = row - 1; i >= 0; --i) {
//       auto &mat_i = mat[i];
//       bi = b[i];  // Cache b[i] for better cache locality
//       for (j = col - 1; j > i; --j)
//          bi = std::fma(-mat_i[j], b[j], bi);
//       b[i] = bi / mat_i[i];
//    }
//    b.resize(row);  // Resize vector to `row`, this will automatically remove extra elements
//    return b;
// };

V_d back_substitution(const VV_d &mat, const V_d &b, const std::size_t row_col) {
   const auto row = row_col;
   const auto col = row_col;
   int i, j;
   V_d result(b);
   double bi;
   for (i = row - 1; i >= 0; --i) {
      const auto &mat_i = mat[i];
      bi = result[i];  // Use result instead of modifying b directly
      for (j = col - 1; j > i; --j)
         // bi -= mat_i[j] * result[j];
         bi = std::fma(-mat_i[j], result[j], bi);
      result[i] = bi / mat_i[i];
   }
   return result;
}

/* -------------------------------------------------------------------------- */
/*                              QR decomposition                              */
/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT QR_decomposition

## QR分解

QR分解は，行列を直交行列と上三角行列の積に分解することである．
このQR分解では，ギブンズ回転を用いてQR分解を行っている．

$`I = F_1^{-1} F_1`$となるような$`F_1`$があるとする．これを使えば，

```math
\begin{align*}
A&=IA\\
&=F_1^{-1} F_1 A \quad {\text{where}}\quad I = F_1^{-1} F_1
\end{align*}
```

となる．これを繰り返した結果
$`A = F_1^{-1} F_2^{-1} ...F_n^{-1} F_n ... F_2 F_1 A`$の内，
$`R = F_n ... F_2 F_1 A`$が上三角行列となってくれれば，これはQR分解である．
$`Q`$は$`Q = F_1^{-1} F_2^{-1} ...F_n^{-1}`$．

### ギブンズ回転

そんな都合がいい$`F_1`$はどうやって作るのか？

少なくとも
$`F_1 = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}`$と決めると，$`F_1`$は直交行列となり，

```math
\begin{align*}
F_1^{-1} &= F_1^{T}\\
\rightarrow I&=F_1^{T} F_1=F_1^{-1} F_1
\end{align*}
```

あとは，$`R = F_n ... F_2 F_1 A`$が上三角行列となってくれるような，$`\cos\theta`$と$`\sin\theta`$を決めればよい．

そうするには，$`F_1 A`$を計算する際に，左下の成分をゼロにするように$`\theta`$を決めればよい．
次に，$`F_2`$を決める際に，$`F_2 F_1 A`$の左下の一つ上の成分をゼロにするように$`\theta`$を決めればよい．

というふうに繰り返す．

後で具体的にどのような値をかけているかここに示しておく．

#### ヘッセンベルグ行列に対するQR分解

$`A`$がヘッセンベルグ行列の場合，
ゼロとする必要がある成分は，$`A_{i+1,i}`$のみである．
そのため，普通のQR分解の計算量が$`O(n^2)`$であるのに対し，ヘッセンベルグ行列の場合は$`O(n)`$である．

*/

template <typename T, bool IsHessenberg = false>
struct QR {
   T Q;
   T QT;
   T R;
   T A;

   // Copy constructor
   // No need to repeat computation; copy the computed Q, R, and A
   QR(const QR &other) : Q(other.Q), QT(other.QT), R(other.R), A(other.A) {
      // DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
   }
   //  std::cout << "QR destructor" << std::endl;
   ~QR() {};
   QR(const T &AIN) : Q(AIN.size(), typename T::value_type(AIN.size(), 0.)), A(AIN), R(AIN) {
      // DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
      Initialize(AIN, true);
   };

   QR &operator=(const QR &other) {
      if (this != &other) {
         this->Q = other.Q;
         this->QT = other.QT;
         this->R = other.R;
         this->A = other.A;
      }
      return *this;
   }

   void Initialize(const T &AIN, const bool constractor = false) {
      /*DOC_EXTRACT

      もし，`A`が`N+1\times N`の行列である場合，
      `Q`は`N+1\times N+1`の行列，
      `R`は`N+1\times N`の行列となる．
      */
      const int N_ROW = AIN.size();
      const int N_COL = AIN[0].size();
      const auto N_ROW_Q = AIN.size();
      const auto N_COL_Q = N_ROW_Q;
      const auto N_ROW_R = AIN.size();
      const auto N_COL_R = N_COL;
      this->A = this->R = AIN;
      // this->Q = std::vector(N_ROW_Q, std::vector<double>(N_COL_Q, 0.));

      this->Q = T(N_ROW_Q, typename T::value_type(N_COL_Q, 0.));

      IdentityMatrix(this->Q);
      this->QT = this->Q;
      // double r, c, s, a, b;
      // std::vector<double> V_COL_i(N_COL), V_ROW_i(Q.size()), V_COL_j(N_COL), V_ROW_j(Q.size());

      typename T::value_type::value_type r, c, s, a, b;
      std::vector<typename T::value_type::value_type> V_COL_i(N_COL), V_ROW_i(Q.size()), V_COL_j(N_COL), V_ROW_j(Q.size());

      for (int j = 0; j < N_COL; ++j) {         // Process up to the second-to-last row
         for (int i = j + 1; i < N_ROW; ++i) {  // Start from the row below j
            if (this->R[i][j] != 0.) {          // Check if the current element is non-zero for rotation
               a = R[j][j];                     // Pivot element in row j, column j
               b = R[i][j];                     // Element to be annihilated in row i, column j
               r = std::hypot(a, b);            // Compute hypotenuse
               if (r == 0.0) {
                  c = 1.0;
                  s = 0.0;
               } else {
                  c = a / r;   // Cosine for Givens rotation
                  s = -b / r;  // Sine for Givens rotation
               }

               // Update R matrix
               for (int col = j; col < N_COL_R; ++col) {
                  V_COL_i[col] = s * R[j][col] + c * R[i][col];  // Rotate row i
                  V_COL_j[col] = c * R[j][col] - s * R[i][col];  // Rotate row j
               }
               for (int col = j; col < N_COL_R; ++col) {
                  R[j][col] = V_COL_j[col];
                  R[i][col] = V_COL_i[col];
               }

               // Update Q matrix
               for (int row = 0; row < N_ROW_Q; ++row) {
                  V_ROW_i[row] = c * Q[row][i] + s * Q[row][j];   // Rotate row j in Q
                  V_ROW_j[row] = -s * Q[row][i] + c * Q[row][j];  // Rotate row i in Q
               }

               for (int row = 0; row < N_ROW_Q; ++row) {
                  Q[row][i] = V_ROW_i[row];
                  Q[row][j] = V_ROW_j[row];
               }

               //
               // std::cout << "Q = " << MatrixForm(Q, 5, 10) << std::endl;
               // std::cout << "QT = " << MatrixForm(Transpose(Q), 5, 10) << std::endl;
               // std::cout << "R = " << MatrixForm(R, 5, 10) << std::endl;
               // std::cin.ignore();
            }
            if (IsHessenberg) {
               /*DOC_EXTRACT
               Hessenberg行列の場合，既に`R`の`j+2`行目以降はゼロになっているので，
               `j+1`行目の`j`列目をゼロにすることで，`R`は上三角行列になる．
               */
               break;
            }
         }
      }
      QT = Transpose(Q);
   };

   void IdentityMatrix(T &mat) {
      std::size_t i = 0;
      for (auto &m : mat) {
         m.assign(m.size(), 0.0);
         m[i++] = 1.0;
      }
   }
};

template <std::size_t N, std::size_t M>
struct QR<std::array<std::array<double, N>, M>> {
   std::array<std::array<double, N>, M> R, A;
   std::array<std::array<double, M>, M> Q, QT;

   // Copy constructor
   // No need to repeat computation; copy the computed Q, R, and A
   QR(const QR &other) : Q(other.Q), QT(other.QT), R(other.R), A(other.A) {
      // DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
   }
   //  std::cout << "QR destructor" << std::endl;
   ~QR() {};
   QR(const std::array<std::array<double, N>, M> &AIN) : R(AIN), A(AIN) {
      // DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
      Initialize(AIN, true);
   };

   QR &operator=(const QR &other) {
      if (this != &other) {
         Q = other.Q;
         QT = other.QT;
         R = other.R;
         A = other.A;
      }
      return *this;
   }

   void Initialize(const std::array<std::array<double, N>, M> &AIN, const bool constractor = false) {
      // DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
      const int N_ROW = M;
      const int N_COL = N;
      const auto N_ROW_Q = M;
      const auto N_COL_Q = N_ROW_Q;
      const auto N_ROW_R = M;
      const auto N_COL_R = N_COL;
      this->A = this->R = AIN;
      IdentityMatrix(this->Q);
      this->QT = this->Q;
      double r, c, s, a, b;
      std::vector<double> V_COL_i(N_COL), V_ROW_i(Q.size()), V_COL_j(N_COL), V_ROW_j(Q.size());
      for (int j = 0; j < N_COL; ++j) {         // Process up to the second-to-last row
         for (int i = j + 1; i < N_ROW; ++i) {  // Start from the row below j
            if (this->R[i][j] != 0.) {          // Check if the current element is non-zero for rotation
               a = R[j][j];                     // Pivot element in row j, column j
               b = R[i][j];                     // Element to be annihilated in row i, column j
               r = std::hypot(a, b);            // Compute hypotenuse
               if (r == 0.0) {
                  c = 1.0;
                  s = 0.0;
               } else {
                  c = a / r;   // Cosine for Givens rotation
                  s = -b / r;  // Sine for Givens rotation
               }

               // Update R matrix
               for (int col = j; col < N_COL_R; ++col) {
                  V_COL_i[col] = s * R[j][col] + c * R[i][col];  // Rotate row i
                  V_COL_j[col] = c * R[j][col] - s * R[i][col];  // Rotate row j
               }
               for (int col = j; col < N_COL_R; ++col) {
                  R[j][col] = V_COL_j[col];
                  R[i][col] = V_COL_i[col];
               }

               // Update Q matrix
               for (int row = 0; row < N_ROW_Q; ++row) {
                  V_ROW_i[row] = c * Q[row][i] + s * Q[row][j];   // Rotate row j in Q
                  V_ROW_j[row] = -s * Q[row][i] + c * Q[row][j];  // Rotate row i in Q
               }
               for (int row = 0; row < N_ROW_Q; ++row) {
                  Q[row][i] = V_ROW_i[row];
                  Q[row][j] = V_ROW_j[row];
               }
            }
         }
      }
      QT = Transpose(Q);
   };

   void IdentityMatrix(auto &mat) {
      std::size_t i = 0;
      for (auto &m : mat) {
         m.fill(0.0);
         m[i++] = 1.0;
      }
   }
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT compressed_row_storage

# 圧縮行格納法 (Compressed Row Storage, CRS)

CRSは疎行列を表現する一つの手法．
このクラスでは，行列-ベクトル積の高速化と，その行列とベクトルの管理とそれらへのアクセスの容易さを目的とし，
次のような考えで，CRSを実装している．

## 実装方法

多くの数値計算で最も高速化したいのは，行列-ベクトル積である．
整数のインデックスを使って，行列-ベクトル積を計算すると，インデックスの管理が大変になる．
例えば，インデックスが消えたり増えたりする場合が大変だ（例えば，要素の節点番号や，粒子法の粒子番号）．

そこで，ポインタをキーとした連装配列として，行列-ベクトル積を計算することが一つの解決策である．
行ベクトルの成分がどの節点や粒子に対応しているかを，ポインタで管理するものである．

また，CRSクラスに，できるだけ`Dot(A,V)`で利用する情報を保存しておくことで，管理とアクセスの容易さを向上させたい．

多くの場合，行ベクトルはある節点や粒子に対して成り立つ方程式を表しているので，
行のインデックスは，その節点や粒子のインデックスと一致する．
強く関連するので，この行ベクトル自体を，それが成り立つ節点や粒子クラスに保存しておくのは，自然なことである．

さらに，そのCRSには方程式`A[i]`だけでなく，節点上または粒子上の`V[i]`も保存しておくことも，自然なことである．

つまり，`Dot(A,V)`の計算を考えて，
CRSには次のような情報を保存しておくことにする．

- 方程式：`A[i]`は，`CRS->column_value`に保存されている．`column_value`は，`std::unordered_map<CRS *, double>`
- 値：`V[i]`は，`CRS->value`に保存されている
- 当然CRSとしてのポインタ
- 行番号(`i`)

特に、Arnoldiプロセスにおける行列-ベクトル積の計算は計算コストが高いです（参照: ArnoldiProcessの行列-ベクトル積）。CRSのDot積を並列化することで、大幅な高速化が可能です（参照: CRSのDot積を並列化）。

## CRS構造体の仕様

### 概要

CRS（Compressed Row Storage）構造体は、疎行列の一部を効率的に格納するためのデータ構造であり、高速な線形代数の計算を実現します。

### メンバ変数

| 変数名 | 型 | 説明 |
|:------:|:--:|:----:|
| `column_value` | `std::unordered_map<CRS *, double>` | 行の非ゼロ要素を格納する連想配列 |
| `value` | `double` | 一般的な値を格納 |
| `diagonal_value` | `double` | 対角要素の値 |
| `tmp_value` | `double` | 一時的な値の格納用 |
| `canUseVector` | `bool` | ベクタが使用可能かのフラグ |
| `value3d` | `std::array<double, 3>` | 3次元空間の値 |
| `__index__` | `std::size_t` | インデックス |

### メンバ関数

| 関数名 | 引数 | 戻り値 | 説明 |
|:------:|:----:|:------:|:----:|
| `clearColumnValue` | なし | `void` | `column_value`をクリアし、`canUseVector`を`false`に設定する |
| `setIndexCRS` | `std::size_t i` | `void` | インデックス`__index__`を設定する |
| `getIndexCRS` | なし | `std::size_t` | インデックス`__index__`を取得する |
| `at` | `CRS *const p` | `double` | 指定された`p`に対応する`column_value`の値を取得する |
| `contains` | `CRS *const p` | `bool` | 指定された`p`が`column_value`に含まれているかを確認する |
| `increment` | `CRS *const p, const double v` | `void` | 指定された`p`に対する`column_value`の値に`v`を加算、または新規挿入する |
| `setVectorCRS` | なし | `void` | `column_value`を`std::vector`形式に変換し、`canUseVector`を`true`に設定する |

*/

// #define GMRES_USING_FLOAT

struct CRS {
   std::unordered_map<CRS *, double> column_value;
   void clearColumnValue() {
      this->column_value.clear();
      this->canUseVector = false;
   };
   double value = 0.;
   double diagonal_value = 0.;
#ifdef GMRES_USING_FLOAT
   float tmp_value = 0.;
   float tmp_tmp_value = 0.;
#else
   double tmp_value = 0.;
   double tmp_tmp_value = 0.;
#endif
   double tmp_b = 0.;
   bool canUseVector = false;
   std::array<double, 3> value3d = {0., 0., 0.};
   std::size_t __index__;
   double *value_ptr = nullptr;
   float *value_ptr_float = nullptr;
   void setIndexCRS(const std::size_t i) { this->__index__ = i; };
   std::size_t getIndexCRS() const { return __index__; };
   CRS() {
      column_value.reserve(50000);
      column_value_vector.reserve(50000);
   };

   // double self_dot_tmp() const {
   //    double ret = 0.;
   //    if (this->canUseVector) {
   //       for (const auto &[crs, value, i] : this->column_value_vector)
   //          ret = std::fma(value, crs->tmp_value, ret);
   //    } else {
   //       for (const auto &[crs, value] : this->column_value)
   //          ret = std::fma(value, crs->tmp_value, ret);
   //    }
   //    return ret;
   // }

   // double self_dot_tmp_store() {
   //    this->tmp_b = 0.;
   //    if (this->canUseVector) {
   //       std::for_each(std::execution::unseq, this->column_value_vector.begin(), this->column_value_vector.end(),
   //                     [&](const auto &c_v) { this->tmp_b = std::fma(std::get<1>(c_v), std::get<0>(c_v)->tmp_value, this->tmp_b); });
   //    } else {
   //       std::for_each(std::execution::unseq, this->column_value.begin(), this->column_value.end(),
   //                     [&](const auto &c_v) { this->tmp_b = std::fma(std::get<1>(c_v), std::get<0>(c_v)->tmp_value, this->tmp_b); });
   //    }
   //    return this->tmp_b;
   // }

   void clear() { this->column_value.clear(); }
   double at(CRS *const p) const { return column_value.at(p); };
   bool contains(CRS *const p) const { return column_value.contains(p); };

   void set(CRS *const p, const double v) {
      if (v == 0.) return;
      auto [it, inserted] = this->column_value.try_emplace(p, v);
      if (!inserted)
         it->second = v;
      this->canUseVector = false;
   };

   void increment(CRS *const p, const double v) {
      if (v == 0.) return;
      auto [it, inserted] = this->column_value.try_emplace(p, v);
      if (!inserted)
         it->second += v;
      this->canUseVector = false;
   };

   // 高速化のために，vectorに変換する．

#ifdef GMRES_USING_FLOAT
   std::vector<std::tuple<CRS *, float, std::size_t>> column_value_vector;
#else
   std::vector<std::tuple<CRS *, double, std::size_t>> column_value_vector;
#endif

   void setVectorCRS() {
      column_value_vector.resize(column_value.size());
      std::size_t i = 0;
      for (const auto &[crs, value] : column_value) {
         // column_value_vector.push_back({crs, value, crs->__index__});
         column_value_vector[i++] = {crs, value, crs->__index__};
      }
      this->canUseVector = true;
   };

#ifdef GMRES_USING_FLOAT
   float selfDotTmpValue() const {
      if (this->canUseVector) {
         return std::transform_reduce(std::execution::unseq,
                                      this->column_value_vector.cbegin(),
                                      this->column_value_vector.cend(),
                                      0.0,
                                      std::plus<>(),
                                      [](const auto &c_v) {
                                         return std::get<1>(c_v) * std::get<0>(c_v)->tmp_value;
                                      });
      } else {
         return std::transform_reduce(std::execution::unseq,
                                      this->column_value.cbegin(),
                                      this->column_value.cend(),
                                      0.0,
                                      std::plus<>(),
                                      [](const auto &c_v) {
                                         return std::get<1>(c_v) * std::get<0>(c_v)->tmp_value;
                                      });
      }
   };
#else
   // double selfDotTmpValue() const {
   //    double ret = 0.;
   //    if (this->canUseVector) {
   //       // std::ranges::for_each(this->column_value_vector, [&](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->tmp_value, ret); });
   //       // std::for_each(std::execution::unseq, this->column_value_vector.begin(), this->column_value_vector.end(),
   //       //               [&ret](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->tmp_value, ret); });

   //       std::for_each(std::execution::unseq, this->column_value_vector.cbegin(), this->column_value_vector.cend(),
   //                     [&ret](const auto &c_v) { ret += std::get<1>(c_v) * std::get<0>(c_v)->tmp_value; });

   //    } else {
   //       // std::ranges::for_each(this->column_value, [&](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->tmp_value, ret); });
   //       std::for_each(std::execution::unseq, this->column_value.begin(), this->column_value.end(),
   //                     [&ret](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->tmp_value, ret); });
   //    }
   //    return ret;
   // };

   void selfDotTmpValue_Save() {
      if (this->canUseVector) {
         this->tmp_tmp_value = std::transform_reduce(std::execution::unseq,
                                                     this->column_value_vector.cbegin(),
                                                     this->column_value_vector.cend(),
                                                     0.0,
                                                     std::plus<>(),
                                                     [](const auto &c_v) {
                                                        return std::get<1>(c_v) * std::get<0>(c_v)->tmp_value;
                                                     });
      } else {
         this->tmp_tmp_value = std::transform_reduce(std::execution::unseq,
                                                     this->column_value.cbegin(),
                                                     this->column_value.cend(),
                                                     0.0,
                                                     std::plus<>(),
                                                     [](const auto &c_v) {
                                                        return std::get<1>(c_v) * std::get<0>(c_v)->tmp_value;
                                                     });
      }
   }

   double selfDotTmpValue() const {
      if (this->canUseVector) {
         return std::transform_reduce(std::execution::unseq,
                                      this->column_value_vector.cbegin(),
                                      this->column_value_vector.cend(),
                                      0.0,
                                      std::plus<>(),
                                      [](const auto &c_v) {
                                         return std::get<1>(c_v) * std::get<0>(c_v)->tmp_value;
                                      });
      } else {
         return std::transform_reduce(std::execution::unseq,
                                      this->column_value.cbegin(),
                                      this->column_value.cend(),
                                      0.0,
                                      std::plus<>(),
                                      [](const auto &c_v) {
                                         return std::get<1>(c_v) * std::get<0>(c_v)->tmp_value;
                                      });
      }
   }
#endif
   double selfDot() const {
      double ret = 0.;
      if (this->canUseVector) {
         // std::ranges::for_each(this->column_value_vector, [&](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->value, ret); });
         std::for_each(std::execution::unseq, this->column_value_vector.begin(), this->column_value_vector.end(),
                       [&ret](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->value, ret); });
      } else {
         // std::ranges::for_each(this->column_value, [&](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->value, ret); });
         std::for_each(std::execution::unseq, this->column_value.begin(), this->column_value.end(),
                       [&ret](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->value, ret); });
      }
      return ret;
   };

   double selfDot(const V_d &V) const {
      double ret = 0.;
      if (this->canUseVector) {
         // std::ranges::for_each(this->column_value_vector, [&](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->value, ret); });
         // for (std::size_t i = 0; const auto &__index : column_value_vector)
         //    ret = std::fma(std::get<1>(column_value_vector[i++]), V[std::get<2>(__index)], ret);
         // for (const auto &[crs, v, i] : column_value_vector)
         //    ret = std::fma(v, V[i], ret);
         std::for_each(std::execution::unseq, this->column_value_vector.begin(), this->column_value_vector.end(),
                       [&ret, &V](const auto &c_v) { ret = std::fma(std::get<1>(c_v), V[std::get<2>(c_v)], ret); });
      } else {
         // std::ranges::for_each(this->column_value, [&](const auto &c_v) { ret = std::fma(std::get<1>(c_v), std::get<0>(c_v)->value, ret); });
         // for (const auto &[crs_local, value] : column_value)
         //    ret = std::fma(value, V[crs_local->__index__], ret);
         std::for_each(std::execution::unseq, this->column_value.begin(), this->column_value.end(),
                       [&ret, &V](const auto &c_v) { ret = std::fma(std::get<1>(c_v), V[std::get<0>(c_v)->__index__], ret); });
      }
      return ret;
   };
};

template <template <typename, typename...> class Container, typename T>
   requires std::derived_from<T, CRS>
void DotOutput(const Container<T *> &A, const V_d &V, V_d &w) {
   w.resize(A.size());
   // \label{CRS:parrallel}
   std::for_each(std::execution::par_unseq, A.begin(), A.end(), [&w, &V](T *crs) {
      w[crs->__index__] = crs->selfDot(V);
   });
}

template <template <typename, typename...> class Container, typename T>
   requires std::derived_from<T, CRS>
V_d Dot(const Container<T *> &A, const V_d &V) {
   V_d ret(A.size());
   std::for_each(std::execution::par_unseq, A.begin(), A.end(), [&ret, &V](T *crs) {
      ret[crs->__index__] = crs->selfDot(V);
   });
   return ret;
};

template <typename T>
   requires std::derived_from<T, CRS>
V_d Dot(const std::vector<T *> &A) {
   V_d ret(A.size(), 0.);
   std::for_each(std::execution::par_unseq, A.begin(), A.end(), [&ret](T *crs) {
      ret[crs->__index__] = crs->selfDot();
   });
   return ret;
};

template <template <typename, typename...> class Container, typename T>
   requires std::derived_from<T, CRS>
void PreparedDotOutput(const Container<T *> &A, V_d &V) {

   std::for_each(std::execution::par_unseq, A.cbegin(), A.cend(), [&](T *crs) {
      double &value = V[crs->__index__];  // V[crs->__index__]に一度だけアクセス
      crs->tmp_value = value;             // tmp_valueに保存
      crs->value_ptr = &value;            // V[index]のポインタを設定
   });

   //! ポインタを使って計算し、Vに直接書き戻す
   std::for_each(std::execution::par_unseq, A.cbegin(), A.cend(), [&](T *crs) {
      *crs->value_ptr = crs->selfDotTmpValue();  // 計算結果をポインタ経由でVに書き戻す
   });
}

// template <template <typename, typename...> class Container, typename T>
//    requires std::derived_from<T, CRS>
// V_d PreparedDot(const Container<T *> &A, V_d V) {
//    //! Vに対するランダムアクセスの回数を減らすために，Vの値をCRSに保存しておく
//    std::for_each(std::execution::par, A.begin(), A.end(), [&](T *crs) {
//       crs->tmp_value = V[crs->__index__];
//    });

//    std::for_each(std::execution::par, A.begin(), A.end(), [&](T *crs) {
//       V[crs->__index__] = crs->selfDotTmpValue();
//    });

//    return V;
// }

/*
PreparedDotは，GMRESのArnoldi過程において，行列-ベクトル積を計算するために使われている．
現在の実装では，

double &value = V[crs->__index__];

VがどのAの行に対応しているかがわからないため，遅くなっている．
もしはじめからわかっていれば，

for(auto &[crs, value] : V) {
   crs->tmp_value = value;
}

std::for_each(std::execution::par_unseq, A.cbegin(), A.cend(), [&](T *crs) {
   *crs->value_ptr = crs->selfDotTmpValue();  // 計算結果をポインタ経由でVに書き戻す
});

そう考えると，AにVの機能を組み込むことで，高速化が可能である．

*/

template <template <typename, typename...> class Container, typename T>
   requires std::derived_from<T, CRS>
V_d PreparedDot(const Container<T *> &A, V_d V) {
   //! V_dの値を取得しつつポインタを設定
   std::for_each(std::execution::par_unseq, A.cbegin(), A.cend(), [&](T *crs) {
      double &value = V[crs->__index__];  // V[crs->__index__]に一度だけアクセス
      crs->tmp_value = value;             // tmp_valueに保存
      crs->value_ptr = &value;            // V[index]のポインタを設定
   });
   //! ポインタを使って計算し、Vに直接書き戻す
   std::for_each(std::execution::par_unseq, A.cbegin(), A.cend(), [&](T *crs) {
      *crs->value_ptr = crs->selfDotTmpValue();  // 計算結果をポインタ経由でVに書き戻す
   });
   return V;
}

template <template <typename, typename...> class Container, typename T>
   requires std::derived_from<T, CRS>
V_f PreparedDot(const Container<T *> &A, V_f V) {
   //! V_dの値を取得しつつポインタを設定
   std::for_each(std::execution::par_unseq, A.cbegin(), A.cend(), [&](T *crs) {
      float &value = V[crs->__index__];  // V[crs->__index__]に一度だけアクセス
      crs->tmp_value = value;            // tmp_valueに保存
      crs->value_ptr_float = &value;     // V[index]のポインタを設定
   });
   //! ポインタを使って計算し、Vに直接書き戻す
   std::for_each(std::execution::par_unseq, A.cbegin(), A.cend(), [&](T *crs) {
      *crs->value_ptr_float = crs->selfDotTmpValue();  // 計算結果をポインタ経由でVに書き戻す
   });
   return V;
}

/* -------------------------------------------------------------------------- */
struct ILU {
   std::vector<CRS *> LU;

   ILU(const std::unordered_set<CRS *> &A) {
      LU.reserve(A.size());

      for (const auto &a : A) {
         LU.push_back(new CRS(*a));  // assuming CRS has a copy constructor
      }

      std::size_t n = LU.size();

      for (std::size_t k = 0; k < n; ++k) {
         if (!LU[k]->contains(LU[k])) {
            throw std::runtime_error("Zero diagonal element in ILU factorization.");
         }
         double diag = LU[k]->at(LU[k]);

         for (auto &[i, lij] : LU[k]->column_value) {
            if (i->getIndexCRS() > k) {
               lij /= diag;

               for (auto &[j, ljk] : LU[i->getIndexCRS()]->column_value) {
                  if (j->getIndexCRS() > i->getIndexCRS()) {
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
      std::size_t n = LU.size();
      V_d x(n), y(n);

      // Forward solve Ly = b
      for (std::size_t i = 0; i < n; ++i) {
         y[i] = b[i];
         for (auto &[j, lij] : LU[i]->column_value) {
            if (j->getIndexCRS() < i) {
               y[i] -= lij * y[j->getIndexCRS()];
            }
         }
      }

      // Backward solve Ux = y
      for (int i = n - 1; i >= 0; --i) {
         x[i] = y[i];
         for (auto &[j, uij] : LU[i]->column_value) {
            if (j->getIndexCRS() > i) {
               x[i] -= uij * x[j->getIndexCRS()];
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

アーノルディ分解は，Krylov部分空間の生成のために使われる．
GMRESは，Krylov部分空間と呼ばれる線形空間内で反復解を探すアルゴリズム．
この部分空間は，行列Aと初期ベクトルから生成される．
GMRESにとってアーノルディ分解は，ヘッセンベルグ行列を生成するための手段．
得られたヘッセンベルグ行列はQR分解され，直交行列と上三角行列に分解される．

1. 正規化した$`{\bf v}_1`$を与えておく．
2. $`{\bf v}_2 = {\rm Normalize}(\,\,\,\quad\quad\quad\quad\quad A{\bf v}_1 - ((A{\bf v}_1) \cdot {\bf v}_1){\bf v}_1\,\,\qquad\qquad\qquad\qquad\qquad\qquad)`$を計算する．
3. $`{\bf v}_3 = {\rm Normalize}(\quad\quad\quad({\bf w}=A{\bf v}_2 - ((A{\bf v}_2) \cdot {\bf v}_1){\bf v}_1)) - ({\bf w} \cdot {\bf v}_2){\bf v}_2\quad\quad\quad\quad\quad\quad)`$を計算する．
4. $`{\bf v}_4 = {\rm Normalize}(({\bf w}=(({\bf w}=A{\bf v}_3 - ((A{\bf v}_3) \cdot {\bf v}_1){\bf v}_1)) - ({\bf w} \cdot {\bf v}_2){\bf v}_2)) - ({\bf w} \cdot {\bf v}_3){\bf v}_3)`$を計算する．

言い換えると，

1. 正規化した$`{\bf v}_1`$を与えておく．
2. $`{\bf w}=A{\bf v}_1, {\bf v}_2 = {\rm Normalize}({\rm Chop}({\bf w},{\bf v}_1))`$を計算する．
3. $`{\bf w}=A{\bf v}_2, {\bf v}_3 = {\rm Normalize}({\rm Chop}({\rm Chop}({\bf w}, {\bf v}_1), {\bf v}_2))`$を計算する．
4. $`{\bf w}=A{\bf v}_3, {\bf v}_4 = {\rm Normalize}({\rm Chop}({\rm Chop}({\rm Chop}({\bf w}, {\bf v}_1), {\bf v}_2), {\bf v}_3))`$を計算する．

これは，既存のベクトルの成分を削り落として，新しいベクトルを作っているに過ぎない．

NOTE: ここで最も計算コストがかかるのは，$`{\bf w}=A{\bf v}_i`$の行列-ベクトル積である．

$`A{\bf v}_i`$の直交化の際に，
それに含まれる各基底$`{\bf v}_0,{\bf v}_1,...,{\bf v}_i`$の成分を計算している．
この成分からなる行列が，Hessenberg行列$`H`$である（ほとんど上三角行列のことを指す）．

```math
\begin{align*}
A{\bf v}_1 & = h_{1,1} {\bf v}_1 + h_{2,1} {\bf v}_2\\
A{\bf v}_2 & = h_{1,2} {\bf v}_1 + h_{2,2} {\bf v}_2 + h_{3,2} {\bf v}_3\\
& \dots\\
A{\bf v}_{n} & = h_{1,n} {\bf v}_1 + h_{2,n} {\bf v}_2 + \cdots + h_{n,n+1} {\bf v}_{n+1}
\end{align*}
```

NOTE: ここで，行数よりも項数が1多いことに注目しよう．

行列を使ってまとめると，

```math
A V_n = V_{n+1} \tilde H_n, \quad V_n = [v_1|v_2|...|v_n],
\quad \tilde H_n = \begin{bmatrix} h_{1,1} & h_{1,2} & \cdots & h_{1,n} & h_{1,n+1} \\ h_{2,1} & h_{2,2} & \cdots & h_{2,n} & h_{2,n+1} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & h_{n,n} & h_{n,n+1} \\ 0 & 0 & \cdots & 0 & h_{n+1,n+1} \end{bmatrix}
```

NOTE: $`\tilde H_n`$は，Hessenberg行列が１行長くなった行列になっている．これは，前の式において，行数よりも項数が1多いことによる．

これをArnoldi分解という．ここで，$`[v_1|v_2|...|v_n]`$の$`|`$は列ベクトルを連結して行列を形成することを示している．

### 基底ベクトルの追加

基底ベクトルを追加したい場合にどのような操作が必要となるか整理しておこう．
これは，GMRES法の繰り返し計算の中で必要となる．

*/

// #define DEBUG_GMRES
//
struct ArnoldiProcess {

   std::size_t n;  // the number of interation
   double beta;
   V_d v0;
   VV_d H;  // ((n+1) x n) Hessenberg matrix
   VV_d V;  // ((n+1) x n) an orthonormal basis of the Krylov subspace like {v0,A.v0,A^2.v0,...}
   V_d w;
   //
   float beta_float;
   V_f v0_float;
   VV_f H_float;  // ((n+1) x n) Hessenberg matrix
   VV_f V_float;  // ((n+1) x n) an orthonormal basis of the Krylov subspace like {v0,A.v0,A^2.v0,...}
   V_f w_float;

   virtual ~ArnoldiProcess() = default;

   ArnoldiProcess(const std::function<V_d(const V_d &v)> &return_A_dot_v,
                  const V_d &v0IN,
                  const std::size_t nIN)
       : n(nIN),
         beta(Norm(v0IN)),
         v0(v0IN / beta),
         H(nIN + 1, V_d(nIN, 0.)),
         V(nIN + 1, v0 /*V[0]=v0であればいい．ここではv0=v1=v2=..としている*/),
         w(v0IN.size()) {
      H.reserve(300);
      V.reserve(300);
      v0.reserve(300);
      w.reserve(300);
      DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
      this->Initialize(return_A_dot_v, v0IN, nIN, false);
   };

   /* -------------------------------------------------------------------------- */

   void Initialize_float(const std::function<V_f(const V_f &v)> &return_A_dot_v, const V_f &v0IN, const int nIN, const bool do_constract = true) {
      DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
      if (do_constract) {
         n = nIN;
         beta = beta_float = Norm(v0IN);
         v0_float = (v0IN / beta_float);
         // 既存のHとVのサイズを変更し、メモリを再利用する
         H_float.resize(nIN + 1);  // Hは (nIN + 1) x nIN のサイズになる
         H.resize(nIN + 1);
         V_float.resize(nIN + 1);  // Vは (nIN + 1) x v0_float.size() のサイズになる
         V.resize(nIN + 1);
         // Hの各行をゼロで初期化
         for (std::size_t i = 0; i < nIN + 1; ++i) {
            H_float[i].assign(nIN, 0.0);  // 既存のメモリを再利用しつつゼロで初期化
            H[i].assign(nIN, 0.0);        // 既存のメモリを再利用しつつゼロで初期化
         }
         // Vの各行をv0で初期化（V[0] = v0_float, 他の行は0にする場合は別途初期化）
         V_float[0] = v0_float;  // V_float[0] は v0_float に初期化される
         for (std::size_t i = 1; i <= nIN; ++i) {
            V_float[i].assign(v0_float.size(), 0.0);  // 他の行はゼロで初期化
            V[i].assign(v0_float.size(), 0.0);        // 他の行はゼロで初期化
         }
         w_float.resize(v0IN.size());
         w.resize(v0IN.size());
      }
      /* --------------------------------------- */
      std::size_t i, j;
      for (j = 0; j < n /*展開項数*/; ++j) {
         w_float = return_A_dot_v(V_float[j]);
         for (i = 0; i <= j; ++i)
            FusedMultiplyIncrement(-(H_float[i][j] = Dot(V_float[i], w_float)), V_float[i], w_float);
         V_float[j + 1] = w_float / (H_float[j + 1][j] = Norm(w_float));
      }
      // copy to double H, V, w, v0, beta
      for (std::size_t i = 0; i < nIN + 1; ++i) {
         for (std::size_t j = 0; j < H_float[i].size(); ++j)
            H[i][j] = H_float[i][j];
         for (std::size_t j = 0; j < V_float[i].size(); ++j)
            V[i][j] = V_float[i][j];
      }
      for (std::size_t i = 0; i < w_float.size(); ++i)
         w[i] = w_float[i];
      beta = beta_float;
   };

   /* -------------------------------------------------------------------------- */

   void Initialize(const std::function<V_d(const V_d &v)> &return_A_dot_v, const V_d &v0IN, const int nIN, const bool do_constract = true) {
      DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
      if (do_constract) {
         n = nIN;
         beta = Norm(v0IN);
         v0 = (v0IN / beta);
         // 既存のHとVのサイズを変更し、メモリを再利用する
         H.resize(nIN + 1);  // Hは (nIN + 1) x nIN のサイズになる
         V.resize(nIN + 1);  // Vは (nIN + 1) x v0.size() のサイズになる
         // Hの各行をゼロで初期化
         for (auto &row : H) {
            row.assign(nIN, 0.0);  // 既存のメモリを再利用しつつゼロで初期化
         }
         // Vの各行をv0で初期化（V[0] = v0, 他の行は0にする場合は別途初期化）
         V[0] = v0;  // V[0] は v0 に初期化される
         for (std::size_t i = 1; i <= nIN; ++i)
            V[i].assign(v0.size(), 0.0);  // 他の行はゼロで初期化
         w.resize(v0IN.size());
      }
      std::size_t i, j;
      double tmp = 0;
      for (j = 0; j < n /*展開項数*/; ++j) {
         w = return_A_dot_v(V[j]);
         for (i = 0; i <= j; ++i)
            FusedMultiplyIncrement(-(H[i][j] = Dot(V[i], w)), V[i], w);
         V[j + 1] = w / (H[j + 1][j] = Norm(w));
      }
   };

   /* -------------------------------------------------------------------------- */

   void AddBasisVector(const std::function<V_d(const V_d &v)> &return_A_dot_v) {
      // std::cout << "AddBasisVector" << std::endl;
      std::size_t i;
      // update n, V, H
      this->n++;
      // w = Dot(A, *V.rbegin());  //@ 行列-ベクトル積\label{ArnoldiProcess:matrix-vector}
      w = return_A_dot_v(*V.rbegin());
      for (auto &h : H)
         h.push_back(0.);
      H.push_back(V_d(n, 0.));
      // orthogonalization
      // std::cout << "orthogonalization" << std::endl;
      for (i = 0; i < n; ++i)
         w -= (*H[i].rbegin() = Dot(w, V[i])) * V[i];
      V.push_back(w / (H[n][n - 1] = Norm(w)));
      // std::cout << "done" << std::endl;
   };
};

/*DOC_EXTRACT GMRES

## 一般化最小残差法 (Generalized Minimal Residual Method, GMRES)

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
& = \|Q^\intercal \|{\bf b}\| {\bf e}_1 - R {\bf y}_n\|\quad \text{(use QR decomposition)}\\
\end{align*}
```

ただし，これは理論の話であって，
簡単には，アーノルディ分解で得られたヘッセンベルグ行列をQR分解して，$`{\bf y}_n`$を求める，ということである．

NOTE: アーノルディ過程が逐次的に計算できるため，展開項数$`n`$を$`n+1`$へと大きくしようとする際に（精度が$`n`$では十分でない場合），GMRESで近似解$`{\bf x}_{n+1}`$を始めから計算しなおす必要はない．$`V_{n+1}`$と$`{\tilde H}_{n+1}`$は，$`V_n`$と$`{\tilde H}_n`$を再利用するようにして計算でき，従って，比較的安く，得られている$`{\bf x}_n`$から$`{\bf x}_{n+1}`$へと更新できる．

### GMRESの計算複雑性

単純に考えて，$`A`$が$`m \times m`$行列であるとすると，
GMRESの行列ベクトル積の計算量は，$`O(m^2)`$．
もし，クリロフ部分空間の基底を$`n`$個まで展開するとすると，$`O(m^2 n)`$の計算量が必要となる．
さらに，収束するまでの反復回数を$`k`$とすると，$`O(m^2 n k)`$の計算量が必要となる．

LU分解の場合は，$`O(m^3)`$の計算量が必要となる．
従って，$`m`$が大きい場合は，GMRESの方が計算量が少なくて済む．

GMRESと多重極展開法（もし$`m`$が$`m/d`$になったとすると）を組み合わせれば，GMRESは$`O(knm^2/d^2)`$で計算できる．

*/

/*
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
*/

struct gmres : public ArnoldiProcess {

   /* NOTE:
   r0 = Normalize(b - A.x0)
   to find {r0,A.r0,A^2.r0,...}
   Therefore V is an orthonormal basis of the Krylov subspace Km(A,r0)
   */
   // int n;  // th number of interation

   V_d x;
   V_d y;
   QR<VV_d, true> qr;
   V_d g;
   double err;
   ~gmres() { std::cout << "destructing gmres" << std::endl; };
   gmres(const std::function<V_d(const V_d &v)> &return_A_dot_v, const V_d &b, const V_d x0, const std::size_t nIN)
       : ArnoldiProcess(return_A_dot_v, b - return_A_dot_v(x0), nIN),
         x(x0),
         y(b.size()),
         qr(ArnoldiProcess::H),
         g(qr.Q.size()),
         err(0.) {
      DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
      if (this->beta /*initial error*/ == static_cast<double>(0.))
         return;
      std::size_t i = 0;
      for (const auto &q : qr.QT)
         this->g[i++] = q[0] * this->beta;
      this->err = std::abs(this->g.back());  // 予想される誤差
      this->g.pop_back();
      this->y = back_substitution(this->qr.R, this->g, this->g.size());
      for (auto i = 0; i < this->n; ++i)
         FusedMultiplyIncrement(this->y[i], this->V[i], this->x);
   };

   // std::vector<double>からstd::vector<float>への変換コンストラクタ
   std::vector<float> convertToFloat(const std::vector<double> &vec) {
      return std::vector<float>(vec.begin(), vec.end());
   }

   /* -------------------------------------------------------------------------- */

   void Restart(const std::function<V_d(const V_d &v)> &return_A_dot_v, const V_d &b, const V_d x0, const int nIN) {
      DebugPrint(Yellow, __FILE__, " ", __PRETTY_FUNCTION__, " ", __LINE__);
      std::cout << "Restart:Restarting" << std::endl;
      this->ArnoldiProcess::Initialize(return_A_dot_v, b - return_A_dot_v(x0), nIN);
      std::cout << "Restart:Initialize done" << std::endl;
      this->x = x0;
      std::cout << "QR" << std::endl;
      this->qr.Initialize(this->H);
      std::cout << "QR done" << std::endl;
      if (this->beta /*initial error*/ == static_cast<double>(0.)) return;
      this->g = this->qr.Q[0] * this->beta;
      err = std::abs(this->g.back());  // 予想される誤差
      this->g.pop_back();
      this->y = back_substitution(this->qr.R, this->g, this->g.size());
      for (auto i = 0; i < this->n; ++i)
         FusedMultiplyIncrement(this->y[i], this->V[i], this->x);
      std::cout << "Restart done" << std::endl;
   }

   /* -------------------------------------------------------------------------- */
   //@今の所このIterateは正しく実装されていない．
   void Iterate(const std::function<V_d(const V_d &v)> &return_A_dot_v) {
      this->AddBasisVector(return_A_dot_v);
      this->qr = QR<VV_d, true>(this->H);
      err = 0.;
      if (this->beta /*initial error*/ == static_cast<double>(0.))
         return;
      std::size_t i = 0;
      this->g = this->qr.Q[0] * this->beta;
      err = std::abs(g.back());  // 予想される誤差
      g.pop_back();
      std::fill(this->x.begin(), this->x.end(), 0);
      this->y = back_substitution(this->qr.R, g, g.size());
      for (i = 0; i < this->n; ++i) {
         // this->x += this->y[i] * this->V[i];
         //! use std::fma
         FusedMultiplyIncrement(this->y[i], this->V[i], this->x);
      }
      // std::cout << "done" << std::endl;
   };
};

/* -------------------------------------------------------------------------- */
V_d Eigenvalues(const VV_d &A, const double tol = 1e-9, const std::size_t maxIter = 1000) {
   VV_d Ak = A, I = A;
   IdentityMatrix(I);
   V_d eigenvalues;
   QR qr(Ak);
   for (std::size_t i = 0; i < maxIter; ++i) {
      qr.Initialize(Ak);
      eigenvalues = Diagonal(Ak = Dot(qr.R, qr.Q));
      if (std::ranges::all_of(eigenvalues, [&](const auto lambda) { return std::abs(Det(A - lambda * I)) < tol; })) {
         return eigenvalues;
      }
   }
   return eigenvalues;
}

std::pair<V_d, VV_d> Eigensystem(const VV_d &A, const double tol = 1e-9, const std::size_t maxIter = 1000) {
   VV_d Ak = A, I = A, Qk = A;
   IdentityMatrix(Qk);
   IdentityMatrix(I);
   V_d eigenvalues;
   QR qr(Ak);
   for (std::size_t i = 0; i < maxIter; ++i) {
      qr.Initialize(Ak);
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
std::array<double, N> Eigenvalues(const std::array<std::array<double, N>, N> &A, const double tol = 1e-9, const std::size_t maxIter = 1000) {
   std::array<double, N> eigenvalues;
   std::size_t i = 0;
   for (const auto &a : Eigenvalues(ToVector(A), tol, maxIter))
      eigenvalues[i++] = a;
   return eigenvalues;
}

template <std::size_t N>
std::pair<std::array<double, N>, std::array<std::array<double, N>, N>> Eigensystem(const std::array<std::array<double, N>, N> &A, const double tol = 1e-9, const std::size_t maxIter = 1000) {
   std::array<std::array<double, N>, N> Ak = A, I = A, Qk = A;
   IdentityMatrix(Qk);
   IdentityMatrix(I);
   std::array<double, N> eigenvalues;
   QR qr(Ak);
   for (std::size_t i = 0; i < maxIter; ++i) {
      qr.Initialize(Ak);
      Ak = Dot(qr.R, qr.Q);        // A(k+1) = R * Q
      eigenvalues = Diagonal(Ak);  // Extract diagonal as eigenvalues
      Qk = Dot(Qk, qr.Q);          // Update the transformation matrix
      if (std::ranges::all_of(eigenvalues, [&](const auto lambda) { return std::abs(Det(A - lambda * I)) < tol; })) {
         return {eigenvalues, Qk};
      }
   }
   return {eigenvalues, Qk};  // Return both eigenvalues and eigenvectors
}

// template <std::size_t N>
// std::pair<std::array<double, N>, std::array<std::array<double, N>, N>> Eigensystem(const std::array<std::array<double, N>, N> &A, const double tol = 1e-9, const std::size_t maxIter = 1000) {
//    auto [E, Qk] = Eigensystem(ToVector(A), tol, maxIter);
//    std::array<std::array<double, N>, N> Q;
//    std::array<double, N> eigenvalues;
//    std::size_t i = 0;
//    for (const auto &e : E)
//       eigenvalues[i++] = e;

//    i = 0;
//    std::size_t j = 0;
//    for (const auto &q : Qk) {
//       for (const auto &qj : q)
//          Q[i][j++] = qj;
//       i++;
//       j = 0;
//    }
//    // normalize vectors and modify eigenvalues
//    for (auto i = 0; i < N; ++i) {
//       double norm = Norm(Q[i]);
//       Q[i] /= norm;
//       eigenvalues[i] *= norm;
//    }

//    return {eigenvalues, Q};
// }

#endif