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
      ludcmp lu(dFdx);
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
   void update(const Tdd &F, const T2Tdd &dFdx) { X += (dX = -Dot(Inverse(dFdx), F)); };
   void update(const Tdd &F, const T2Tdd &dFdx, const double a) { X += a * (dX = -Dot(Inverse(dFdx), F)); };
};
template <>
struct NewtonRaphson<Tddd> : public NewtonRaphson_Common<Tddd> {
   NewtonRaphson(const Tddd &Xinit) : NewtonRaphson_Common<Tddd>(Xinit){};
   void update(const Tddd &F, const T3Tddd &dFdx) { X += (dX = -Dot(Inverse(dFdx), F)); };
};
template <>
struct NewtonRaphson<T4d> : public NewtonRaphson_Common<T4d> {
   NewtonRaphson(const T4d &Xinit) : NewtonRaphson_Common<T4d>(Xinit) {
      this->ans = V_d(4, 0.);
   };
   V_d ans;
   void update(const T4d &F, const T4T4d &dFdx) {
      ludcmp lu(ToVector(dFdx));
      lu.solve(ToVector(-F), ans);
      std::get<0>(X) += (std::get<0>(dX) = ans[0]);
      std::get<1>(X) += (std::get<1>(dX) = ans[1]);
      std::get<2>(X) += (std::get<2>(dX) = ans[2]);
      std::get<3>(X) += (std::get<3>(dX) = ans[3]);
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
            Fs[i] += (w * (vec[i] - NRs[i].X)) * drdx;  //<- d/dx (d*d)
            dFs[i] += drdx * drdx;
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
            Fs[j] += w * (sample_vectors[i][j] - NRs[j].X) * drdx;
            dFs[j] += drdx * drdx;
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

#endif