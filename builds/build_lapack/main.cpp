// LAPACK test code
// compile with: g++ main.cpp -llapack -lblas -o testprog

// #include <iostream>
// #include <vector>

// using namespace std;

// extern "C" void dgetrf_(int *dim1, int *dim2, double *a, int *lda, int *ipiv, int *info);
// extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

// int main()
// {
//   char trans = 'N';
//   int dim = 2;
//   int nrhs = 1;
//   int LDA = dim;
//   int LDB = dim;
//   int info;

//   vector<double> a, b;

//   a.push_back(1);
//   a.push_back(0);
//   a.push_back(0);
//   a.push_back(1);

//   b.push_back(2);
//   b.push_back(0);

//   int ipiv[3];

//   dgetrf_(&dim, &dim, &*a.begin(), &LDA, ipiv, &info);
//   dgetrs_(&trans, &dim, &nrhs, &*a.begin(), &LDA, ipiv, &*b.begin(), &LDB, &info);

//   std::cout << "solution is:";
//   std::cout << "[" << b[0] << ", " << b[1] << ", "
//             << "]" << std::endl;
//   std::cout << "Info = " << info << std::endl;

//   return (0);
// }
/* ------------------------------------------------------ */
// LAPACK test code
// compile with: g++ main.cpp -llapack -lblas -o testprog

#include <fundamental.hpp>

int main()
{
  std::vector<std::vector<double>> a({{1.4, 1.2, 3.0, 4.3},
                                      {1.5, 4.5, 2.0, 1.2},
                                      {5.2, 0.5, 7.0, 4.2},
                                      {2.0, 7.0, 1.0, 2.2}});
  std::vector<double> b({2, 1.4, 4, 1}), ans(3);
  /* ------------------------------------------------------ */
  lapack_lu la(a);
  la.solve(b, ans);
  std::cout << "solution is:";
  std::cout << ans << std::endl;
  /* ------------------------------------------------------ */
  ludcmp lu(a);
  lu.solve(b, ans);
  std::cout << "solution is:";
  std::cout << ans << std::endl;
  /* ------------------------------------------------------ */
};