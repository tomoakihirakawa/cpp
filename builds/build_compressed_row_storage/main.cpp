#include <iostream>
#include <map>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include "basic_linear_systems.hpp"
#include "basic_vectors.hpp"
#include "lib_measurement.hpp"

class CRSMatrix {
  private:
   std::vector<double> values;
   std::vector<int> col_indices;
   std::vector<int> row_pointers;
   int numRows, numCols;

  public:
   CRSMatrix(int numRows, int numCols) : numRows(numRows), numCols(numCols) {
      row_pointers.assign(numRows + 1, 0);
   }

   double operator()(int row, int col) const {
      int index = row_pointers[row];
      while (index < row_pointers[row + 1] && col_indices[index] != col) {
         index++;
      }

      if (index < row_pointers[row + 1]) {
         return values[index];
      }

      return 0;
   }

   double& operator()(int row, int col) {
      int index = row_pointers[row];
      while (index < row_pointers[row + 1] && col_indices[index] != col) {
         index++;
      }

      if (index < row_pointers[row + 1]) {
         return values[index];
      }

      col_indices.insert(col_indices.begin() + row_pointers[row + 1], col);
      values.insert(values.begin() + row_pointers[row + 1], 0);
      for (int i = row + 1; i <= numRows; ++i) {
         row_pointers[i]++;
      }

      return values[index];
   }

   int size() const {
      return numRows;
   }
   int numOmittedElements() const {
      int totalElements = numRows * numCols;
      int nonZeroElements = values.size();
      return totalElements - nonZeroElements;
   }
   std::vector<double> CRS_Dot(const std::vector<double>& vec) const {
      std::vector<double> result(numRows, 0);

      for (int i = 0; i < numRows; ++i) {
         for (int j = row_pointers[i]; j < row_pointers[i + 1]; ++j) {
            result[i] += values[j] * vec[col_indices[j]];
         }
      }

      return result;
   }
};

class CRS_QR {
  private:
   CRSMatrix Q, R;
   int n, m;

   void IdentityMatrix(CRSMatrix& mat) {
      for (int i = 0; i < mat.size(); ++i) {
         mat(i, i) = 1.0;
      }
   }

  public:
   CRS_QR(const CRSMatrix& AIN)
       : R(AIN),
         Q(AIN.size(), AIN.size()) {
      IdentityMatrix(Q);
      n = AIN.size();
      m = AIN.size();  // Assuming square matrix
      double r, c, s, a, b;
      double Q0, Q1, R0, R1;
      double eps = 1e-12;

      for (auto j = 0; j < m; ++j) {
         for (auto i = j; i <= n - 2; ++i) {
            if (fabs(R(i + 1, j)) > eps) {
               a = R(i, j);
               b = R(i + 1, j);

               r = std::sqrt(a * a + b * b);
               s = -b / r;
               c = a / r;

               if (s != s || c != c) {
                  s = 0.;
                  c = 1.;
               }

               for (auto k = 0; k < m; ++k) {
                  R0 = c * R(i, k) - s * R(i + 1, k);
                  R1 = s * R(i, k) + c * R(i + 1, k);
                  R(i, k) = R0;
                  R(i + 1, k) = R1;
               }

               for (auto k = 0; k < n; ++k) {
                  Q0 = c * Q(i, k) - s * Q(i + 1, k);
                  Q1 = s * Q(i, k) + c * Q(i + 1, k);
                  Q(i, k) = Q0;
                  Q(i + 1, k) = Q1;
               }
            }
         }
      }
   }

   CRSMatrix getQ() const { return Q; }
   CRSMatrix getR() const { return R; }
};

int main() {
   CRSMatrix matrix(4, 4);

   matrix(0, 1) += 1;
   matrix(1, 2) = 2;
   matrix(2, 0) = 3;
   matrix(2, 0) = 4;
   matrix(3, 3) = 4;
   std::cout << "numOmittedElements: " << matrix.numOmittedElements() << "\n";

   std::vector<double> vec = {1, 2, 3, 4};

   std::vector<double> result = matrix.CRS_Dot(vec);

   for (auto i = 0; i < 4; i++) {
      for (auto j = 0; j < 4; j++) {
         std::cout << matrix(i, j) << " ";
      }
      std::cout << std::endl;
   }

   for (double val : result) {
      std::cout << val << " ";
   }

   CRS_QR qr(matrix);
   CRSMatrix Q = qr.getQ();
   CRSMatrix R = qr.getR();
   std::cout << "Q numOmittedElements: " << Q.numOmittedElements() << "\n";
   std::cout << "R numOmittedElements: " << R.numOmittedElements() << "\n";

   std::cout << "size: " << Q.size() << std::endl;
   std::cout << "size: " << R.size() << std::endl;

   for (auto i = 0; i < 4; i++) {
      for (auto j = 0; j < 4; j++) {
         std::cout << Q(i, j) << " ";
      }
      std::cout << std::endl;
   }

   for (auto i = 0; i < 4; i++) {
      for (auto j = 0; j < 4; j++) {
         std::cout << R(i, j) << " ";
      }
      std::cout << std::endl;
   }

   return 0;
}
