#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"

auto w = std::setw(20);

// struct that put nodes of a robot on the LightHill curve

struct LightHillRobot {
   double L;
   double w;
   double k;
   double c1;
   double c2;
   int n;

   LightHillRobot(double L, double w, double k, double c1, double c2, int n)
       : L(L), w(w), k(k), c1(c1), c2(c2), n(n){};

   auto yLH(const double x, const double t) { return (c1 * x / L + c2 * std::pow(x / L, 2)) * sin(k * (x / L) - w * t); };

   auto X_RB(const std::array<double, 2> a, const double q) {
      double r = L / n;
      return a + r * std::array<double, 2>{cos(q), sin(q)};
   };

   auto f(const std::array<double, 2> a, const double q, const double t) {
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
      double error = 0;
      for (auto i = 1; i < Q.size(); i++) {
         NewtonRaphson nr(Q[i - 1]);
         error = 0;
         for (auto k = 0; k < 50; k++) {
            auto F = f(a, nr.X, t);
            nr.update(F * F / 2., F * ddq_f(nr.X, t));
            if ((error = std::abs(F)) < 1E-5)
               break;
         }
         Q[i] = nr.X;
         a = X_RB(a, Q[i]);
      }
      return Q;
   };

   std::vector<std::array<double, 2>> anglesToX(const V_d Q) {
      std::array<double, 2> a = {0., 0.};
      std::vector<std::array<double, 2>> ret(Q.size());
      ret[0] = a;
      for (auto i = 1; i < Q.size(); i++)
         ret[i] = (a = X_RB(a, Q[i]));
      return ret;
   };
};

int main() {

   /*DOC_EXTRACT newton

    ## ロボットの節の位置をLightHillの曲線上に乗せる

    LightHillの式は以下のようになる．

    $$
    {\bf x}^{\rm LH}(x,t) = (x,y^{\rm LH}(x,t)),\quad
    y^{\rm LH}(x,t) = \left( \frac{c_1}{L} x + {c_2} \left(\frac{x}{L}\right)^2 \right) \sin \left( \frac{2 \pi}{L} x - \omega t \right)
    $$

    ここで，$c_1, c_2, L, \omega$は定数である．

   | variable | meaning |
   |:---:|:---:|
   | $L$ | 全長 |
   | $\omega$ | 角周波数 |
   | $k$ | 波数 |
   | $c_1$ | 振幅1 |
   | $c_2$ | 振幅2 |
   | $n$ | number of nodes of the robot |
   | $r$ | length of a node of the robot |
   | $\theta_i$ | angle of the $i$-th node of the robot |

    ロボットの$i$番目の節の位置は，${\bf x}_{i}^{\rm rb} = {\bf x}_{i-1}^{\rm rb} + r \left( \cos \theta_i, \sin \theta_i \right)$である．
    次の関数を使って表すことにする．

    $$
    {\bf x}_{i}^{\rm rb} = X^{\rm rb}({\bf x}_{i-1}^{\rm rb},r,\theta_i),
    \quad X^{\rm rb}({\bf a},r,\theta) = {\bf a} + r \left( \cos \theta, \sin \theta \right)
    $$

    ここで，$r$はロボットの節の長さ，$\theta_i$はロボットの節の角度である．
    頭の位置を$X_{0}^{\rm rb}=(0,0)$とする．
    次の節の位置は，$X_{1}^{\rm rb} = r (\cos \theta_1, \sin \theta_1)$である．
    さらに次の節の位置は，$X_{2}^{\rm rb} = X_{1}^{\rm rb} + r (\cos \theta_2, \sin \theta_2)$となる．

    目的関数$f$の微分は，

    $$
    \frac{df}{d\theta} = -r \sin\theta\frac{d y^{\rm LH} }{dx}-r\cos\theta
    $$

    NOTE: この目的関数$f$には，前の節の位置を与える必要がある．節の位置は，後ろの節の位置によって変わらないので，この目的関数を先頭から順番に最適化することは問題ない．

   ![sample.gif](sample.gif)

    */

   TimeWatch time;

   std::cout << time() << std::endl;

   double L = 0.71;
   double w = 2. * M_PI * 1.0;
   double k = 2. * M_PI * 2.0;
   double c1 = 0.05;
   double c2 = 0.05;
   int nodes = 10;
   int steps = 20;

   LightHillRobot lhr(L, w, k, c1, c2, nodes);

   for (auto i = 0; i < steps; i++) {
      std::ofstream outFile("./output_lighthill/lighthill" + std::to_string(i) + ".txt");
      double t = (double)i / steps;
      auto Q = lhr.getAngles(t);
      auto xy = lhr.anglesToX(Q);
      for (auto j = 0; j < xy.size(); j++) {
         auto [x, y] = xy[j];
         outFile << x << " " << y << " " << Q[j] << std::endl;
         std::cout << x << " " << y << " " << Q[j] << std::endl;
      }
      outFile.close();
   }
}
