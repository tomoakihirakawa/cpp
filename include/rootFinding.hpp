#ifndef rootFinding_H
#define rootFinding_H
#pragma once

#include "basic.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

template <typename T>
struct NewtonRaphson_Common {
   T X;
   T dX;  // tmp
   NewtonRaphson_Common(const T &Xinit) : X(Xinit), dX(Xinit){};
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
   NewtonRaphson(const double Xinit) : NewtonRaphson_Common<double>(Xinit){};
   void update(const double F, const double dFdx) { X += (dX = F / (-dFdx)); };
   void update(const double F, const double dFdx, const double a) { X += a * (dX = F / (-dFdx)); };
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
template <>
struct NewtonRaphson<T7d> : public NewtonRaphson_Common<T7d> {
   NewtonRaphson(const T7d &Xinit) : NewtonRaphson_Common<T7d>(Xinit) {
      this->ans = V_d(7, 0.);
   };
   V_d ans;
   void update(const T7d &F, const T7T7d &dFdx) {
      ludcmp lu(ToVector(dFdx));
      lu.solve(ToVector(-F), ans);
      std::get<0>(X) += (std::get<0>(dX) = ans[0]);
      std::get<1>(X) += (std::get<1>(dX) = ans[1]);
      std::get<2>(X) += (std::get<2>(dX) = ans[2]);
      std::get<3>(X) += (std::get<3>(dX) = ans[3]);
      std::get<4>(X) += (std::get<4>(dX) = ans[4]);
      std::get<5>(X) += (std::get<5>(dX) = ans[5]);
      std::get<6>(X) += (std::get<6>(dX) = ans[6]);
   };
};

/* ------------------------------------------------------ */
/*                           利用例                        */
/* ------------------------------------------------------ */
Tddd optimumVector(const std::vector<Tddd> &sample_vectors, const Tddd &init_vector) {
   NewtonRaphson NR0(std::get<0>(init_vector)), NR1(std::get<1>(init_vector)), NR2(std::get<2>(init_vector));
   double r, F0, F1, F2, dF0dx, dF1dx, dF2dx, drdx, w;
   for (auto i = 0; i < 100; ++i) {
      F0 = F1 = F2 = dF0dx = dF1dx = dF2dx = 0;
      for (const auto &vec : sample_vectors) {
         w = 1;  // kernel_Bspline3(Norm(X - p->getXtuple()), p->radius);
         r = w * (std::get<0>(vec) - NR0.X);
         drdx = -w;
         F0 += r * r;  //<- d/dx (d*d)
         dF0dx += 2. * r * drdx;
         r = w * (std::get<1>(vec) - NR1.X);
         drdx = -w;
         F1 += r * r;  //<- d/dx (d*d)
         dF1dx += 2. * r * drdx;
         r = w * (std::get<2>(vec) - NR2.X);
         drdx = -w;
         F2 += r * r;  //<- d/dx (d*d)
         dF2dx += 2. * r * drdx;
      }
      NR0.update(F0, dF0dx);
      NR1.update(F1, dF1dx);
      NR2.update(F2, dF2dx);
      if (std::abs(NR0.dX) < 1E-12 && std::abs(NR1.dX) < 1E-12 && std::abs(NR2.dX) < 1E-12)
         break;
   }
   return {NR0.X, NR1.X, NR2.X};
};

double optimumVector_(const std::vector<double> &sample_vectors, const double init_vector) {
   NewtonRaphson NR0(init_vector);
   double r, F0, F1, F2, dF0dx, dF1dx, dF2dx, drdx, w;
   for (auto i = 0; i < 100; ++i) {
      F0 = F1 = F2 = dF0dx = dF1dx = dF2dx = 0;
      for (const auto &vec : sample_vectors) {
         w = 1;  // kernel_Bspline3(Norm(X - p->getXtuple()), p->radius);
         r = w * (vec - NR0.X);
         drdx = -w;
         F0 += 2. * r * drdx;  //<- d/dx (d*d)
         dF0dx += 2. * drdx * drdx;
      }
      NR0.update(F0, dF0dx);
      if (std::abs(NR0.dX) < 1E-12)
         break;
   }
   return NR0.X;
};

Tddd optimumVector_(const std::vector<Tddd> &sample_vectors, const Tddd &init_vector) {
   NewtonRaphson NR0(std::get<0>(init_vector)), NR1(std::get<1>(init_vector)), NR2(std::get<2>(init_vector));
   double r, F0, F1, F2, dF0dx, dF1dx, dF2dx, drdx, w;
   for (auto i = 0; i < 100; ++i) {
      F0 = F1 = F2 = dF0dx = dF1dx = dF2dx = 0;
      for (const auto &vec : sample_vectors) {
         w = 1;  // kernel_Bspline3(Norm(X - p->getXtuple()), p->radius);
         r = w * (std::get<0>(vec) - NR0.X);
         drdx = -w;
         F0 += 2. * r * drdx;  //<- d/dx (d*d)
         dF0dx += 2. * drdx * drdx;
         r = w * (std::get<1>(vec) - NR1.X);
         drdx = -w;
         F1 += 2. * r * drdx;  //<- d/dx (d*d)
         dF1dx += 2. * drdx * drdx;
         r = w * (std::get<2>(vec) - NR2.X);
         drdx = -w;
         F2 += 2. * r * drdx;  //<- d/dx (d*d)
         dF2dx += 2. * drdx * drdx;
      }
      NR0.update(F0, dF0dx);
      NR1.update(F1, dF1dx);
      NR2.update(F2, dF2dx);
      if (std::abs(NR0.dX) < 1E-12 && std::abs(NR1.dX) < 1E-12 && std::abs(NR2.dX) < 1E-12)
         break;
   }
   return {NR0.X, NR1.X, NR2.X};
};

Tddd optimumVector_(const std::vector<Tddd> &sample_vectors, const Tddd &init_vector, std::vector<double> weight) {
   if (weight.size() != sample_vectors.size())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   weight /= Max(weight);
   NewtonRaphson NR0(std::get<0>(init_vector)), NR1(std::get<1>(init_vector)), NR2(std::get<2>(init_vector));
   double r, F0, F1, F2, dF0dx, dF1dx, dF2dx, drdx, w;
   Tddd vec;
   for (auto i = 0; i < 200; ++i) {
      F0 = F1 = F2 = dF0dx = dF1dx = dF2dx = 0;
      for (auto i = 0; i < sample_vectors.size(); ++i) {
         vec = sample_vectors[i];
         w = weight[i];
         r = w * (std::get<0>(vec) - NR0.X);
         drdx = -w;
         F0 += 2. * r * drdx;  //<- d/dx (d*d)
         dF0dx += 2. * drdx * drdx;
         r = w * (std::get<1>(vec) - NR1.X);
         drdx = -w;
         F1 += 2. * r * drdx;  //<- d/dx (d*d)
         dF1dx += 2. * drdx * drdx;
         r = w * (std::get<2>(vec) - NR2.X);
         drdx = -w;
         F2 += 2. * r * drdx;  //<- d/dx (d*d)
         dF2dx += 2. * drdx * drdx;
      }
      NR0.update(F0, dF0dx);
      NR1.update(F1, dF1dx);
      NR2.update(F2, dF2dx);
      if (std::abs(NR0.dX) < 1E-13 && std::abs(NR1.dX) < 1E-13 && std::abs(NR2.dX) < 1E-13)
         break;
   }
   return {NR0.X, NR1.X, NR2.X};
};

double optimumVector_(const std::vector<double> &sample_vectors, const double init_vector, std::vector<double> weight) {
   if (weight.size() != sample_vectors.size())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   weight /= Max(weight);
   NewtonRaphson NR0(init_vector);
   double r, F0, F1, F2, dF0dx, dF1dx, dF2dx, drdx, w;
   double vec;
   for (auto i = 0; i < 200; ++i) {
      F0 = F1 = F2 = dF0dx = dF1dx = dF2dx = 0;
      for (auto i = 0; i < sample_vectors.size(); ++i) {
         vec = sample_vectors[i];
         w = weight[i];
         r = w * (vec - NR0.X);
         drdx = -w;
         F0 += 2. * r * drdx;  //<- d/dx (d*d)
         dF0dx += 2. * drdx * drdx;
      }
      NR0.update(F0, dF0dx);
      if (std::abs(NR0.dX) < 1E-13)
         break;
   }
   return NR0.X;
};

Tddd optimumVector_(const std::vector<Tddd> &sample_vectors,
                    const Tddd &init_vector,
                    const std::vector<Tddd> &N,
                    const std::vector<double> &weight) {
   NewtonRaphson NR0(std::get<0>(init_vector)), NR1(std::get<1>(init_vector)), NR2(std::get<2>(init_vector));
   double r, F0, F1, F2, dF0dx, dF1dx, dF2dx, drdx, w;
   Tddd vec;
   for (auto i = 0; i < 100; ++i) {
      F0 = F1 = F2 = dF0dx = dF1dx = dF2dx = 0;
      for (auto i = 0; i < sample_vectors.size(); ++i) {
         vec = sample_vectors[i];
         auto n = N[i];
         w = weight[i];
         Tddd X = {NR0.X, NR1.X, NR2.X};
         r = w * Dot(vec - X, n);
         drdx = -w * std::get<0>(n);
         F0 += 2. * r * drdx;  //<- d/dx (d*d)
         dF0dx += 2. * drdx * drdx;
         r = w * Dot(vec - X, n);
         drdx = -w * std::get<1>(n);
         F1 += 2. * r * drdx;  //<- d/dx (d*d)
         dF1dx += 2. * drdx * drdx;
         r = w * Dot(vec - X, n);
         drdx = -w * std::get<2>(n);
         F2 += 2. * r * drdx;  //<- d/dx (d*d)
         dF2dx += 2. * drdx * drdx;
      }
      if (isFinite(dF0dx))
         NR0.update(F0, dF0dx);
      if (isFinite(dF1dx))
         NR1.update(F1, dF1dx);
      if (isFinite(dF2dx))
         NR2.update(F2, dF2dx);
      if (isFinite(NR0.dX) && isFinite(NR1.dX) && isFinite(NR2.dX))
         if (std::abs(NR0.dX) < 1E-12 && std::abs(NR1.dX) < 1E-12 && std::abs(NR2.dX) < 1E-12)
            break;
   }
   Tddd ret = {NR0.X, NR1.X, NR2.X};
   if (isFinite(ret))
      return ret;
   else
      return {0., 0., 0.};
};

Tddd optimumVector(const std::vector<Tddd> &sample_vectors,
                   const Tddd &init_vector,
                   const std::vector<Tddd> &N,
                   std::vector<double> weight) {
   NewtonRaphson NR0(std::get<0>(init_vector)), NR1(std::get<1>(init_vector)), NR2(std::get<2>(init_vector));
   double r, F0, F1, F2, dF0dx, dF1dx, dF2dx, drdx, w;
   Tddd vec;
   for (auto i = 0; i < 100; ++i) {
      F0 = F1 = F2 = dF0dx = dF1dx = dF2dx = 0;
      for (auto i = 0; i < sample_vectors.size(); ++i) {
         vec = sample_vectors[i];
         auto n = N[i];
         w = weight[i];
         Tddd X = {NR0.X, NR1.X, NR2.X};
         r = w * Dot(vec - X, n);
         drdx = -w * std::get<0>(n);
         F0 += r * r;  //<- d/dx (d*d)
         dF0dx += 2. * r * drdx;
         drdx = -w * std::get<1>(n);
         F1 += r * r;  //<- d/dx (d*d)
         dF1dx += 2. * r * drdx;
         drdx = -w * std::get<2>(n);
         F2 += r * r;  //<- d/dx (d*d)
         dF2dx += 2. * r * drdx;
      }
      if (F0 + F1 + F2 < 1E-5)
         break;
      if (isFinite(dF0dx))
         NR0.update(F0, dF0dx);
      if (isFinite(dF1dx))
         NR1.update(F1, dF1dx);
      if (isFinite(dF2dx))
         NR2.update(F2, dF2dx);
      if (isFinite(NR0.dX) && isFinite(NR1.dX) && isFinite(NR2.dX))
         if (std::abs(NR0.dX) < 1E-12 && std::abs(NR1.dX) < 1E-12 && std::abs(NR2.dX) < 1E-12)
            break;
   }
   return {NR0.X, NR1.X, NR2.X};
};
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
      this->k = nr.X;
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