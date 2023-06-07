/*DOC_EXTRACT pybind11

# pybind11の例

この例は，c++のNewton法を利用して作ったLight Hill Robotをpythonで使うためのもの.
以下をターミナルで実行して`make`すると，Macだと`LightHillRobot_pybind.cpython-311-darwin.so`が作られる.

```
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./
make
```

`run_robot.py`にあるように`import`して利用できる．

*/

#define NOMINMAX
#define _CRT_SECURE_NO_WARNINGS

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../include/minMaxOfFunctions.hpp"
#include "../../include/rootFinding.hpp"

struct LightHillRobot {
   double L;
   double w;
   double k;
   double c1;
   double c2;
   int n;  // node + 1 (head node is dummy)

   LightHillRobot(double L, double w, double k, double c1, double c2, int n)
       : L(L), w(w), k(k), c1(c1), c2(c2), n(n + 1){};

   auto yLH(const double x, const double t) { return (c1 * x / L + c2 * std::pow(x / L, 2)) * sin(k * (x / L) - w * t); };

   auto X_RB(const std::array<double, 2>& a, const double q) {
      double r = L / n;
      return a + r * std::array<double, 2>{cos(q), sin(q)};
   };

   auto f(const std::array<double, 2>& a, const double q, const double t) {
      auto [x, y] = X_RB(a, q);
      return yLH(x, t) - y;
   };

   auto ddx_yLH(const double x, const double t) {
      return (c1 / L + 2 * c2 * x / std::pow(L, 2)) * sin(k * (x / L) - w * t) +
             (c1 / L * x + c2 * std::pow(x / L, 2)) * cos(k * (x / L) - w * t) * k / L;
   };

   auto ddq_f(const double q, const double t) {
      double r = L / n;
      double x = r * cos(q);
      return -r * sin(q) * ddx_yLH(x, t) - r * cos(q);
   };

   V_d getAngles(const double t) {
      std::vector<double> Q(n, 0.);  // thetas
      std::array<double, 2> a{{0., 0.}};
      double error = 0, F;
      double scale = 0.25, q0;
      for (auto i = 0; i < Q.size(); i++) {
         q0 = atan(ddx_yLH(a[0], t));
         NewtonRaphson nr(q0);
         error = 0;
         for (auto k = 0; k < 100; k++) {
            F = f(a, nr.X, t);
            nr.update(F * F / 2., F * ddq_f(nr.X, t), scale);
            if ((error = std::abs(F)) < 1E-10)
               break;
         }
         Q[i] = nr.X;
         a = X_RB(a, Q[i]);
      }
      return Q;
   };

   std::vector<std::array<double, 2>> anglesToX(const V_d& Q) {
      std::array<double, 2> a = {0., 0.};
      std::vector<std::array<double, 2>> ret;
      ret.reserve(Q.size() + 1);
      ret.push_back(a);
      for (auto i = 0; i < Q.size(); i++)
         ret.push_back(a = X_RB(a, Q[i]));
      return ret;
   };
};

namespace py = pybind11;

PYBIND11_MODULE(LightHillRobot_pybind, m) {
   py::class_<LightHillRobot>(m, "LightHillRobot")
       .def(py::init<double, double, double, double, double, int>())
       .def_readwrite("c1", &LightHillRobot::c1)
       .def_readwrite("c2", &LightHillRobot::c2)
       .def("yLH", &LightHillRobot::yLH)
       .def("X_RB", &LightHillRobot::X_RB)
       .def("f", &LightHillRobot::f)
       .def("ddx_yLH", &LightHillRobot::ddx_yLH)
       .def("ddq_f", &LightHillRobot::ddq_f)
       .def("getAngles", &LightHillRobot::getAngles)
       .def("anglesToX", &LightHillRobot::anglesToX);
}
