#include <cmath>
#include <iostream>

int main() {

   auto f = [](auto x) { return std::sin(x); };
   auto analitical_integral = [](auto a, auto b) { return -(std::cos(b) - std::cos(a)); };

   double a = 0, b = M_PI;
   for (int N = 1; N < 100; ++N) {
      double integral = (f(a) + f(b)) / 2.;
      double dx = (b - a) / N;
      for (int i = 1; i < N; ++i)
         integral += f(a + i * dx);
      integral *= dx;
      std::cout << "N = " << N << ", numerical = " << integral << ", analitical = " << analitical_integral(a, b) << std::endl;
   }
   return 0;
};

// #include <iostream>
// #include <cmath>
// #include <vector>
// #include <Eigen/Dense>

// using namespace std;
// using namespace Eigen;

// double analyticalIntegral(const Vector3d &p1, const Vector3d &p2, const Vector3d &p3) {
//     const double epsilon = 1e-12;

//     // Compute edge vectors
//     Vector3d e1 = p2 - p1;
//     Vector3d e2 = p3 - p1;
//     Vector3d e3 = p3 - p2;

//     // Compute normal vector
//     Vector3d n = e1.cross(e2).normalized();

//     // Check if the origin is on the triangle plane
//     double d = -n.dot(p1);
//     bool onPlane = abs(d) < epsilon;

//     // Check if the origin is on the triangle
//     Vector3d v1 = -p1;
//     Vector3d v2 = -p2;
//     Vector3d v3 = -p3;
//     double area = e1.cross(e2).norm() / 2;
//     double s = v1.cross(v2).norm() / (2 * area);
//     double t = v2.cross(v3).norm() / (2 * area);
//     bool insideTriangle = s >= 0 && t >= 0 && s + t <= 1;

//     if (insideTriangle) {
//         // If the origin is inside the triangle, return the special case value
//         double L1 = e1.norm();
//         double L2 = e2.norm();
//         double L3 = e3.norm();
//         double beta = acos(-e1.dot(e3) / (L1 * L3));
//         double alpha = asin(sin(beta) * L1 / L2);

//         return -2 * M_PI * (1 - cos(alpha));
//     } else {
//         // If the origin is not inside the triangle, compute the integral using the general formula
//         Vector3d f = (e1 - e3).cross(e2 - e3);
//         Vector3d g = (v1 - v3).cross(v2 - v3);
//         double rho = g.norm() / f.norm();

//         double a = e1.dot(e2) / (e1.norm() * e2.norm());
//         double b = e1.dot(e3) / (e1.norm() * e3.norm());
//         double c = e2.dot(e3) / (e2.norm() * e3.norm());

//         double I = (acos(a) - acos(b) + acos(c)) / (1 - rho * rho);

//         if (onPlane) {
//             return I / (4 * M_PI);
//         } else {
//             return (2 * M_PI - I) / (4 * M_PI);
//         }
//     }
// }

// int main() {
//     Vector3d p1(0, 0, 0);
//     Vector3d p2(1, 0, 0);
//     Vector3d p3(0, 1, 0);

//     double result = analyticalIntegral(p1, p2, p3);
//     cout << "Analytical Integral Result: " << result << endl;

//     return 0;
// }
