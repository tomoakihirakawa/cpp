#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

extern "C" void dgesvd_(const char *jobu, const char *jobvt,
                        const int *m, const int *n,
                        double *a, const int *lda,
                        double *s,
                        double *u, const int *ldu,
                        double *vt, const int *ldvt,
                        double *work, const int *lwork, int *info);

void computeSVD(const std::vector<std::vector<double>> &A,
                std::vector<std::vector<double>> &U,
                std::vector<double> &S,
                std::vector<std::vector<double>> &VT) {
   int m = A.size();
   int n = A[0].size();
   int lda = m;
   int ldu = m;
   int ldvt = n;

   std::vector<double> a(m * n);
   std::vector<double> u(m * m);
   std::vector<double> vt(n * n);
   S.resize(std::min(m, n));

   // Flatten A into a 1D array
   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         a[i + j * m] = A[i][j];
      }
   }

   // Workspace query
   double work_query;
   int lwork = -1;
   int info;

   dgesvd_("A", "A", &m, &n, a.data(), &lda, S.data(), u.data(), &ldu, vt.data(), &ldvt, &work_query, &lwork, &info);

   if (info < 0) {
      throw std::runtime_error("Illegal value in argument during workspace query.");
   }

   // Allocate workspace
   lwork = static_cast<int>(work_query);
   std::vector<double> work(lwork);

   // Perform SVD
   dgesvd_("A", "A", &m, &n, a.data(), &lda, S.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork, &info);

   if (info < 0) {
      throw std::runtime_error("Illegal value in argument during SVD computation.");
   } else if (info > 0) {
      throw std::runtime_error("SVD did not converge.");
   }

   // Reshape U and VT into 2D arrays
   U.resize(m, std::vector<double>(m));
   VT.resize(n, std::vector<double>(n));

   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
         U[i][j] = u[i + j * m];
      }
   }

   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
         VT[i][j] = vt[i + j * n];
      }
   }
}

std::vector<double> solveWithSVD(const std::vector<std::vector<double>> &A, const std::vector<double> &b) {
   int m = A.size();
   int n = A[0].size();

   if (b.size() != static_cast<size_t>(m)) {
      throw std::runtime_error("Dimension mismatch between A and b.");
   }

   std::vector<std::vector<double>> U, VT;
   std::vector<double> S;

   computeSVD(A, U, S, VT);

   // Compute U^T * b
   std::vector<double> c(m, 0.0);
   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
         c[i] += U[j][i] * b[j];
      }
   }

   // Solve for y = S^(-1) * c
   std::vector<double> y(std::min(m, n), 0.0);
   for (int i = 0; i < std::min(m, n); ++i) {
      if (S[i] > 1e-9) {
         y[i] = c[i] / S[i];
      }
   }

   // Compute x = V * y
   std::vector<double> x(n, 0.0);
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < std::min(m, n); ++j) {
         x[i] += VT[j][i] * y[j];
      }
   }

   return x;
}

std::vector<std::vector<double>> computePseudoInverse(const std::vector<std::vector<double>> &A) {
   int m = A.size();
   int n = A[0].size();

   std::vector<std::vector<double>> U, VT;
   std::vector<double> S;

   computeSVD(A, U, S, VT);

   // Compute pseudoinverse
   std::vector<std::vector<double>> pseudoInverse(n, std::vector<double>(m, 0.0));
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
         for (int k = 0; k < std::min(m, n); ++k) {
            if (S[k] > 1e-9) {
               pseudoInverse[i][j] += VT[k][i] * (1.0 / S[k]) * U[j][k];
            }
         }
      }
   }

   return pseudoInverse;
}

int main() {
   // Define a 3x3 matrix
   std::vector<std::vector<double>> A = {
       {1.0, 2.0, 3.0},
       {0.0, 1.0, 6.0},
       {7.0, 2.0, 9.0}};

   // Containers for SVD results
   std::vector<std::vector<double>> U;
   std::vector<double> S;
   std::vector<std::vector<double>> VT;

   try {
      computeSVD(A, U, S, VT);

      // Print the results
      std::cout << "Matrix U:" << std::endl;
      for (const auto &row : U) {
         for (double val : row) {
            std::cout << std::setw(10) << val << " ";
         }
         std::cout << std::endl;
      }

      std::cout << "\nSingular values (S):" << std::endl;
      for (double val : S) {
         std::cout << std::setw(10) << val << " ";
      }
      std::cout << std::endl;

      std::cout << "\nMatrix VT:" << std::endl;
      for (const auto &row : VT) {
         for (double val : row) {
            std::cout << std::setw(10) << val << " ";
         }
         std::cout << std::endl;
      }

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }

   std::vector<double> b = {1.0, 1.0, 1.0};

   try {
      std::vector<double> x = solveWithSVD(A, b);

      std::cout << "Solution x:" << std::endl;
      for (double val : x) {
         std::cout << val << " ";
      }
      std::cout << std::endl;

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }

   try {
      std::vector<std::vector<double>> A_pseudoInverse = computePseudoInverse(A);

      // Print the pseudoinverse
      std::cout << "Pseudoinverse of A:" << std::endl;
      for (const auto &row : A_pseudoInverse) {
         for (double val : row) {
            std::cout << std::setw(10) << val << " ";
         }
         std::cout << std::endl;
      }

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }

   return 0;
}
