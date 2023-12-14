#ifndef rootFinding_H
#define rootFinding_H
#pragma once

#include "basic.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

template <typename T>
struct NewtonRaphson_Common {
   T X, dX;  // tmp
   NewtonRaphson_Common(const T &Xinit) : X(Xinit), dX(Xinit){};
   void initialize(const T &Xin) { X = Xin; };
};

template <typename T>
struct NewtonRaphson : public NewtonRaphson_Common<T> {
   NewtonRaphson(const T &Xinit) : NewtonRaphson_Common<T>(Xinit){};
};

template <>
struct NewtonRaphson<V_d> : public NewtonRaphson_Common<V_d> {
   NewtonRaphson(const V_d &Xinit) : NewtonRaphson_Common<V_d>(Xinit){};
   void update(const V_d &F, const VV_d &dFdx) {
      lapack_lu lu(dFdx);
      lu.solve(-F, dX);
      X += dX;
   };
};
/* ------------------------------------------------------ */
template <>
struct NewtonRaphson<double> : public NewtonRaphson_Common<double> {
   NewtonRaphson(const double Xinit = 0) : NewtonRaphson_Common<double>(Xinit){};
   void update(const double F, const double dFdx) {
      if (std::abs(dFdx) > 1E-20)
         X += (dX = F / (-dFdx));
   };
   void update(const double F, const double dFdx, const double a) {
      if (std::abs(dFdx) > 1E-20)
         X += a * (dX = F / (-dFdx));
   };
};
template <>
struct NewtonRaphson<Tdd> : public NewtonRaphson_Common<Tdd> {
   NewtonRaphson(const Tdd &Xinit) : NewtonRaphson_Common<Tdd>(Xinit){};
   void update(const Tdd &F, const T2Tdd &dFdx) {
      lapack_lu(dX, dFdx, -F);
      X += dX;
   };
   void update(const Tdd &F, const T2Tdd &dFdx, const double a) {
      lapack_lu(dX, dFdx, -F);
      X += a * dX;
   };
};
template <>
struct NewtonRaphson<Tddd> : public NewtonRaphson_Common<Tddd> {
   NewtonRaphson(const Tddd &Xinit) : NewtonRaphson_Common<Tddd>(Xinit){};
   void update(const Tddd &F, const T3Tddd &dFdx) {
      lapack_lu(dX, dFdx, -F);
      X += dX;
   };
};
template <>
struct NewtonRaphson<T4d> : public NewtonRaphson_Common<T4d> {
   NewtonRaphson(const T4d &Xinit) : NewtonRaphson_Common<T4d>(Xinit) {
      this->ans = V_d(4, 0.);
   };
   V_d ans;
   void update(const T4d &F, const T4T4d &dFdx) {
      lapack_lu(dX, dFdx, -F);
      X += dX;
      // ludcmp lu(ToVector(dFdx));
      // lu.solve(ToVector(-F), ans);
      // std::get<0>(X) += (std::get<0>(dX) = ans[0]);
      // std::get<1>(X) += (std::get<1>(dX) = ans[1]);
      // std::get<2>(X) += (std::get<2>(dX) = ans[2]);
      // std::get<3>(X) += (std::get<3>(dX) = ans[3]);
   };
};

/* ------------------------------------------------------ */
/*                           利用例                        */
/* ------------------------------------------------------ */

template <std::size_t N>
std::array<double, N> optimumVector(const std::vector<std::array<double, N>> &sample_vectors,
                                    const std::array<double, N> &init_vector,
                                    const double tolerance = 1E-12) {
   std::array<NewtonRaphson<double>, N> NRs;
   for (std::size_t i = 0; i < N; ++i)
      NRs[i].X = init_vector[i];
   std::array<double, N> Fs, dFs;
   double w, drdx;
   for (auto j = 0; j < 500; ++j) {
      Fs.fill(0);
      dFs.fill(0);
      for (const auto &vec : sample_vectors) {
         w = 1;
         drdx = -w;
         for (std::size_t i = 0; i < N; ++i) {
            // Fs[i] += (w * (vec[i] - NRs[i].X)) * drdx;  //<- d/dx (d*d)
            // dFs[i] += drdx * drdx;
            // use std::fma
            Fs[i] = std::fma(w * (vec[i] - NRs[i].X), drdx, Fs[i]);
            dFs[i] = std::fma(drdx, drdx, dFs[i]);
         }
      }
      bool converged = true;
      for (std::size_t i = 0; i < N; ++i) {
         NRs[i].update(Fs[i], dFs[i]);
         if (std::abs(NRs[i].dX) >= tolerance) converged = false;
      }
      if (converged) break;
   }
   std::array<double, N> result;
   for (std::size_t i = 0; i < N; ++i) result[i] = NRs[i].X;
   return result;
}

template <std::size_t N>
std::array<double, N> optimumVector(const std::vector<std::array<double, N>> &sample_vectors,
                                    const std::array<double, N> &init_vector,
                                    const std::vector<double> &weights,
                                    const double tolerance = 1E-12) {
   if (weights.size() != sample_vectors.size())
      throw std::runtime_error("The size of the weights vector must match the size of the sample_vectors vector.");

   double max_weight = *std::max_element(weights.begin(), weights.end());
   std::vector<double> normalized_weights(weights.size());
   std::transform(weights.begin(), weights.end(), normalized_weights.begin(), [&](double w) { return w / max_weight; });

   std::array<NewtonRaphson<double>, N> NRs;
   for (std::size_t i = 0; i < N; ++i)
      NRs[i].X = init_vector[i];

   std::array<double, N> Fs, dFs;
   double w, drdx;
   for (int iteration = 0; iteration < 500; ++iteration) {
      Fs.fill(0);
      dFs.fill(0);
      for (size_t i = 0; i < sample_vectors.size(); ++i) {
         w = normalized_weights[i];
         drdx = -w;
         for (size_t j = 0; j < N; ++j) {
            // Fs[j] += w * (sample_vectors[i][j] - NRs[j].X) * drdx;
            // dFs[j] += drdx * drdx;
            // use std::fma
            Fs[j] = std::fma(w * (sample_vectors[i][j] - NRs[j].X), drdx, Fs[j]);
            dFs[j] = std::fma(drdx, drdx, dFs[j]);
         }
      }

      bool converged = true;
      for (std::size_t i = 0; i < N; ++i) {
         NRs[i].update(Fs[i], dFs[i]);
         if (std::abs(NRs[i].dX) >= tolerance) converged = false;
      }
      if (converged) break;
   }

   std::array<double, N> result;
   for (std::size_t i = 0; i < N; ++i) result[i] = NRs[i].X;
   return result;
}

double optimumValue(const std::vector<double> &sample_values,
                    const double init_value,
                    std::vector<double> weights,
                    const double tolerance = 1E-12) {
   if (weights.size() != sample_values.size())
      throw std::runtime_error("The size of the weights vector must match the size of the sample_vectors vector.");

   double m = 0;
   for (const auto &max : weights)
      if (m < std::abs(max))
         m = std::abs(max);
   weights /= m;
   NewtonRaphson<double> NRs;
   NRs.X = init_value;

   double Fs, dFs;
   double w, drdx;
   for (int iteration = 0; iteration < 500; ++iteration) {
      Fs = 0;
      dFs = 0;
      for (size_t i = 0; i < sample_values.size(); ++i) {
         w = weights[i];
         drdx = -w;
         Fs += w * (sample_values[i] - NRs.X) * drdx;
         dFs += drdx * drdx;
      }

      bool converged = true;
      NRs.update(Fs, dFs);
      if (std::abs(NRs.dX) >= tolerance)
         converged = false;
      if (converged) break;
   }

   return NRs.X;
}

/* -------------------------------------------------------------------------- */

struct DispersionRelation {
   double w;
   double T;
   double k;
   double L;
   double h;
   DispersionRelation(const double wIN, const double hIN) : w(wIN), h(hIN), T(2 * M_PI / wIN) {
      NewtonRaphson nr(1.);
      for (auto i = 0; i < 10; i++)
         nr.update(omega(nr.X, h) - w, domegadk(nr.X, h));
      this->k = std::abs(nr.X);
      this->L = 2 * M_PI / this->k;
   };
   double omega(const double k, const double h) {
      return Sqrt(_GRAVITY_ * k * Tanh(h * k));
   };

   double domegadk(const double k, const double h) {
      return (_GRAVITY_ * (h * k * Power(Sech(h * k), 2) + Tanh(h * k))) / (2. * Sqrt(_GRAVITY_ * k * Tanh(h * k)));
   };
};

/* -------------------------------------------------------------------------- */

// \label{newton:LighthillRobot}
struct LighthillRobot {
   double L;
   double w;
   double k;
   double c1;
   double c2;
   int n;  // node + 1 (head node is dummy)

   LighthillRobot(double L, double w, double k, double c1, double c2, int n)
       : L(L), w(w), k(k), c1(c1), c2(c2), n(n + 1){};

   auto yLH(const double x, const double t) { return (c1 * x / L + c2 * std::pow(x / L, 2)) * std::sin(k * (x / L) - w * t); };

   auto X_RB(const std::array<double, 2> &a, const double q) {
      double r = L / n;
      return a + r * std::array<double, 2>{std::cos(q), std::sin(q)};
   };

   auto f(const std::array<double, 2> &a, const double q, const double t) {
      auto [x, y] = X_RB(a, q);
      return yLH(x, t) - y;
   };

   auto ddx_yLH(const double x, const double t) {
      const double a = k * (x / L) - w * t;
      return (c1 / L + 2 * c2 * x / std::pow(L, 2)) * std::sin(a) + (c1 / L * x + c2 * std::pow(x / L, 2)) * std::cos(a) * k / L;
   };

   auto ddq_f(const double q, const double t) {
      const double r = L / n;
      const double x = r * std::cos(q);
      return -r * std::sin(q) * ddx_yLH(x, t) - r * std::cos(q);
   };

   V_d getAngles(const double t) {
      std::vector<double> Q(n, 0.);
      std::array<double, 2> a{{0., 0.}};
      double error = 0, F, scale = 0.3;  //\label{LighthillRobot:scale}
      NewtonRaphson nr(0.);
      for (auto i = 0; i < Q.size(); i++) {
         nr.initialize(std::atan(ddx_yLH(std::get<0>(a), t)));
         error = 0;
         for (auto k = 0; k < 100; k++) {
            F = f(a, nr.X, t);
            nr.update(F * F * 0.5, F * ddq_f(nr.X, t), scale);
            if ((error = std::abs(F)) < 1E-10) break;
         }
         Q[i] = nr.X;
         a = X_RB(a, Q[i]);
      }
      return Q;
   };

   std::vector<std::array<double, 2>> anglesToX(const V_d &Q) {
      std::array<double, 2> a = {0., 0.};
      std::vector<std::array<double, 2>> ret;
      ret.reserve(Q.size() + 1);
      ret.push_back(a);
      for (auto i = 0; i < Q.size(); i++)
         ret.push_back(a = X_RB(a, Q[i]));
      return ret;
   };
};

#endif