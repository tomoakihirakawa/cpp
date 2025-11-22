#pragma once

// Extracted Moments struct from lib_spatial_partitioning.hpp (nested Buckets<T,N>::Moments)
// This header is preparatory: after integration, remove the nested struct in lib_spatial_partitioning.hpp
// and include this header instead. Until that removal, DO NOT include this file in a TU that also
// includes the old nested definition to avoid ODR / duplicate symbol issues.

#include <array>
#include <complex>
#include <cstddef>
#include <iostream>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Project specific dependencies (provide SphericalCoordinates, Fourier2DConvolution, constants, helpers)
#include "lib_Fourier.hpp"             // Fourier2DConvolution
#include "lib_multipole_expansion.hpp" // make_zero_MM, make_nm_set, computeSizeM2M, AAA_* tables, i_absk_A_FMM

// Forward declaration of Buckets to allow pointer members without full definition
template <typename T, int N> struct Buckets;

// Moments class extracted from Buckets<T,N>
// T: underlying object type used by Buckets (only appears indirectly via Buckets<T,N>* caches)
// N: truncation order of multipole/local expansions

//@ -------------------------------------------------------------------------- */
//@              モーメント:　ある観測点で，ソース点を近似するための係数配列             */
//@ -------------------------------------------------------------------------- */
template <typename T, int N> struct Moments {
  /*
  N = 0 [ 0,  ,  ,  ,  ] n=0   total 1 = (N+1)^2             center : 0, in 1D 0
  N = 1 [-1, 0, 1,  ,  ] n=1   total 1 + 3 = 4 = (N+1)^2     center : 1, in 1D 2
  N = 2 [-2,-1, 0, 1, 2] n=2=N total 4 + 5 = 9 = (N+1)^2     center : 2, in 1D 6
        <-- 2*N+1=5 -->
  */
  std::array<double, 3> X;
  using type_MM = std::array<std::array<cmplx, 2>, (N + 1) * (N + 1)>;
  type_MM MM_ = make_zero_MM<N>();
  //! Mは，sum {q_i*\rho_i^n*Y_n^-m}
  //! Greengard & Rokhlin(1997)の式(3.19)
  constexpr static std::array<std::array<int, 2>, (N + 1) * (N + 1)> nm_set = make_nm_set<N>();

  Moments() : X{0.0, 0.0, 0.0} {}
  Moments(const std::array<double, 3> &XIN) : X(XIN) {}

  const cmplx zero = {0., 0.};

  std::size_t index(int n, int m) const { return n * (n + 1) + m; }

  /* -------------------------------------------------------------------------- */

  void initialize() { this->MM_ = make_zero_MM<N>(); }

  void initialize(const std::array<double, 3> &XIN) {
    this->X = XIN;
    this->MM_ = make_zero_MM<N>();
  }

  /* -------------------------------------------------------------------------- */

  using TupleType = std::tuple<std::array<cmplx, 2>, std::array<cmplx, 2> *>;
  using ArrayType = std::array<TupleType, (N + 1) * (N + 1)>;
  std::vector<std::tuple<Tdd *, ArrayType>> vector_weightedSourceDensities_RRn_MM_;

  //^OK
  template <typename TYPE> void increment_M(const std::vector<TYPE> &sources) {
    std::array<int, 2> n_m;
    this->vector_weightedSourceDensities_RRn_MM_.clear();
    this->vector_weightedSourceDensities_RRn_MM_.reserve(sources.size());
    for (const auto &source : sources) {
      auto [w_phi, w_phin] = source->weighted_source_densities;
      // ! 重みとは，数値積分w0w1と変数変換ヤコビアンの積
      // ! w_phiとw_phinは，q_iに相当する in G&R 1997 (3.19)
      // ! BIEは，Rにはw_phinが，Rnにはw_phiがかかる形になっている
      // ! R = {rho_i^n Y_n^-m}
      SphericalCoordinates kernel(source->X - this->X);
      ArrayType RRn_MM_;
      std::array<cmplx, 2> RRn;
      for (std::size_t ind = nm_set.size(); ind-- > 0;) {
        auto &[n, m] = nm_set[ind];
        RRn = kernel.p2mFunction(n, m, source->normal);
        std::get<0>(MM_[ind]) += w_phin * std::get<0>(RRn) /*３重和の１項 G&R 1997 (3.19)*/;
        std::get<1>(MM_[ind]) += w_phi * std::get<1>(RRn) /*３重和の１項*/;
        RRn_MM_[ind] = {RRn, &MM_[ind]}; //! {R, grad_Rnm1_dot_normal}
      }; //! Wは{R, grad_Rnm1_dot_normal}
      this->vector_weightedSourceDensities_RRn_MM_.emplace_back(&source->weighted_source_densities, RRn_MM_);
    }
  }

  //^OK
  void increment_M_reuse() {
    std::array<double, 2> w_phi_w_phin;
    for (const auto &[weighted_source_densities, RRn_MM_] : this->vector_weightedSourceDensities_RRn_MM_) {
      w_phi_w_phin = *weighted_source_densities;
      for (auto &[RRn, mm_] : RRn_MM_) {
        std::get<0>(*mm_) += std::get<1>(w_phi_w_phin) * std::get<0>(RRn);
        std::get<1>(*mm_) += std::get<0>(w_phi_w_phin) * std::get<1>(RRn);
      }
    }
  }

  //^OK
  std::array<double, 2> L2P(const std::array<double, 3> &a) const {
    SphericalCoordinates P(a - this->X);
    std::array<cmplx, 2> ret = {0, 0};
    cmplx Ynm_rhon1;
    for (std::size_t ind = nm_set.size(); ind-- > 0;) {
      auto &[n, m] = nm_set[ind];
      Ynm_rhon1 = P.l2pFunction(n, m);
      std::get<0>(ret) += Ynm_rhon1 * std::get<0>(this->MM_[ind]);
      std::get<1>(ret) += Ynm_rhon1 * std::get<1>(this->MM_[ind]);
    }
    //   return {std::get<0>(ret).real(), std::get<1>(ret).real()};
    return {std::get<0>(ret).real(), std::get<1>(ret).real()};
  }

  //$ ----------------------------------- M2M ---------------------------------- */

  /*DOC_EXTRACT M

  ### 定数の読み込み

  `AAA_M2M_FMM`，`AAA_M2L_FMM`，`AAA_L2L_FMM`は．`定数.nb`あらかじめ計算しておいたものを読み込む．

  */

  static constexpr std::size_t ARRAY_SIZE = computeSizeM2M(N);

  using Type_IICCIIIC = std::tuple<int, int, std::array<cmplx, 2> *, int, int, int, cmplx>;
  std::vector<Type_IICCIIIC> index_map_for_M2M = index_map_for_M2M_(); // i,j,index, n, m, index
  constexpr std::vector<Type_IICCIIIC> index_map_for_M2M_() {
    std::vector<Type_IICCIIIC> index_map_for_M2M;
    index_map_for_M2M.reserve(ARRAY_SIZE);
    std::size_t index = 0;
    int j, k, n, m;
    //   for (j = 0; j <= N; ++j)
    for (j = N; j >= 0; --j)
      for (k = -j; k <= j; ++k) {
        auto &AAA_M2M_FMM_j_k = AAA_M2M_FMM[j][k + N_AAA_M2M_FMM];
        // for (n = 0; n <= j; ++n)
        for (n = j; n >= 0; --n)
          for (m = -n; m <= n; ++m) {
            auto AAA = AAA_M2M_FMM_j_k[n][m + N_AAA_M2M_FMM];
            if (AAA.real() != 0.0 || AAA.imag() != 0.0)
              index_map_for_M2M.emplace_back(j, k, &this->MM_[this->index(j, k)], n, m, this->index(j - n, k - m), AAA);
          }
      }
    return index_map_for_M2M;
  }

  std::vector<std::tuple<std::unique_ptr<SphericalCoordinates>, Moments *>> m2m_cache;

  //! これは計算毎に変わる
  template <typename TYPE> void set_m2m(const std::vector<TYPE> &children) {
    this->m2m_cache.resize(children.size()); // これでOK
    int i = 0;
    for (const auto &b : children) {
      //  auto sph = std::make_unique<SphericalCoordinates>(b->X - this->X);
      if (std::get<0>(this->m2m_cache[i]) == nullptr) {
        auto sph = std::make_unique<SphericalCoordinates>((b->X - this->X));
        sph->precompute_sph(2 * N);
        this->m2m_cache[i] = {std::move(sph), &b->MomentsMultipoleExpansion};
      } else {
        std::get<0>(this->m2m_cache[i])->initialize((b->X - this->X));
        std::get<0>(this->m2m_cache[i])->precompute_sph(2 * N);
        std::get<1>(this->m2m_cache[i]) = &b->MomentsMultipoleExpansion;
      }
      i++;
    }
  }

  void m2m() {
    cmplx R;
    for (const auto &[sph, moments_class] : this->m2m_cache) {
      for (const auto &[j, k, MMjk, n, m, index, AAA] : this->index_map_for_M2M) {
        // R = AAA * sph->sph_harmonics_rho(n, -m);
        R = sph->m2mFunction(j, k, n, m);
        std::get<0>(*MMjk) += R * std::get<0>(moments_class->MM_[index]);
        std::get<1>(*MMjk) += R * std::get<1>(moments_class->MM_[index]);
      }
    }
  }

  //! 4重和のインデックスとそれがさす自身Moment係数へのポインターを，あらかじめ計算しておく．
  std::vector<Type_IICCIIIC> index_map_for_M2L = index_map_for_M2L_();
  constexpr std::vector<Type_IICCIIIC> index_map_for_M2L_() {
    std::vector<Type_IICCIIIC> ret;
    ret.reserve((N + 1) * (2 * N + 1) * (N + 1) * (2 * N + 1));
    int j, k, n, m;
    //   for (j = 0; j <= N; ++j)
    for (j = N; j >= 0; --j)
      for (k = -j; k <= j; ++k) {
        auto &AAA_M2L_FMM_j_k = AAA_M2L_FMM[j][k + N_AAA_M2L_FMM];
        // for (n = 0; n <= N; ++n)
        for (n = N; n >= 0; --n)
          for (m = -n; m <= n; ++m) {
            auto AAA = AAA_M2L_FMM_j_k[n][m + N_AAA_M2L_FMM];
            if (AAA.real() != 0.0 || AAA.imag() != 0.0)
              ret.emplace_back(j, k, &this->MM_[this->index(j, k)], n, m, this->index(n, m), AAA);
          }
      }
    return ret;
  }

  //% DFT畳み込み用

  // #define SimpleM2L
#define FourierM2L_double
  // #define FourierM2L_DoubleDouble
  // #define FourierM2L_BlockDecomposition

#if defined(SimpleM2L)
  using TREATED = double;
#elif defined(FourierM2L_double)
  using TREATED = double;
#elif defined(FourierM2L_DoubleDouble)
  using TREATED = DoubleDouble;
#elif defined(FourierM2L_BlockDecomposition)
  using TREATED = double;
  //% --------------------------------
#define N_block 4
  int getBlockIndex(const int n) const {
    if (n <= 2)
      return 0;
    else if (n <= 5)
      return 1;
    else if (n <= 8)
      return 2;
    else
      return 3;
  }
#endif

  // using TREATED = std::float128_t;
  using INPUT = std::complex<double>; //! これは固定とする．なぜならstd::float128_tならこれで問題なく収束したから
  using INPUT2D = std::array<std::array<INPUT, 2 * N + 1>, N + 1>;
  using fourier2DConvolutionINPUT = Fourier2DConvolution<INPUT, double, N + 1, 2 * N + 1, N + 1, 2 * N + 1>;
  using fourier2DConvolution = Fourier2DConvolution<INPUT, TREATED, N + 1, 2 * N + 1, N + 1, 2 * N + 1>;

#if defined(FourierM2L_BlockDecomposition)
  std::array<std::vector<std::tuple<std::shared_ptr<fourier2DConvolution>, Buckets<T, N> *>>, N_block> Yqp_Onm0_Onm1_cache_blocked;
#endif
  std::vector<std::tuple<std::shared_ptr<fourier2DConvolutionINPUT>, Buckets<T, N> *>> Yqp_Onm0_Onm1_cache;
  fourier2DConvolutionINPUT Onm0_DFT_M2L, Onm1_DFT_M2L;
  using cd = std::complex<double>;
  std::vector<std::tuple<std::array<cd, 2> *, std::vector<std::tuple<cd, std::array<cd, 2> *>>>> m2l_cache;

  /* --------------------------------------------------------------------------- */
  // \label{Moments::set_m2l}

  template <typename TYPE> void set_m2l(std::vector<TYPE> children) {
    std::unordered_map<std::array<cd, 2> *, std::vector<std::tuple<cd, std::array<cd, 2> *>>> uo_m2l_cache;
    uo_m2l_cache.reserve(m2l_cache.size());

    this->Yqp_Onm0_Onm1_cache.resize(children.size(), {nullptr, nullptr});

    INPUT2D Ypq{};
    // Yqp_Onm0_Onm1_cache_blocked
#if defined(FourierM2L_BlockDecomposition)
    for (auto &v : this->Yqp_Onm0_Onm1_cache_blocked)
      v.clear();
    for (auto &v : this->Yqp_Onm0_Onm1_cache_blocked)
      v.reserve(children.size() / N_block + 1);
#endif
    this->m2l_cache.clear();
    // 削減率を出そう．
    int total_operations = 0;
    int total_reduced_operations = 0;
    const double threshold_AAAY = 0.; // 係数の閾値,スレッショルド
    std::sort(children.begin(), children.end(), [X = this->X](const auto &lhs, const auto &rhs) { return Norm(lhs->X - X) > Norm(rhs->X - X); });
    //
    int i = 0;
    for (auto &b : children) {
      auto sph = SphericalCoordinates(b->MomentsMultipoleExpansion.X - this->X);
      sph.precompute_sph(2 * N);
      //! 4重和に対応
      for (auto &[j, k, MMjk_local, n, m, index, AAA] : this->index_map_for_M2L) {
        auto AAAY = sph.m2lFunction(j, k, n, m);
        if (threshold_AAAY < std::abs(AAAY)) {
          uo_m2l_cache[MMjk_local].emplace_back(AAAY, &b->MomentsMultipoleExpansion.MM_[index]);
          total_reduced_operations++;
        }
        total_operations++;
      }
      //$ ===================================== */
#if defined(FourierM2L_double) || defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition)
      for (int p = 0; p <= N; ++p)
        for (int q = -p; q <= p; ++q)
          Ypq[p][q + N] = sph.Ypq(p, -q);
#endif
#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
      if (std::get<0>(this->Yqp_Onm0_Onm1_cache[i]) == nullptr)
        this->Yqp_Onm0_Onm1_cache[i] = {std::make_shared<fourier2DConvolutionINPUT>(Ypq, true, false), b};
      else {
        std::get<0>(this->Yqp_Onm0_Onm1_cache[i])->initialize(Ypq, true, false);
        std::get<1>(this->Yqp_Onm0_Onm1_cache[i]) = b;
      }
#elif defined(FourierM2L_BlockDecomposition)
      std::array<INPUT2D, N_block> Yqp_blocked{};
      for (int p = 0; p <= N; ++p)
        for (int q = -p; q <= p; ++q)
          Yqp_blocked[getBlockIndex(p)][p][q + N] = Ypq[p][q + N];
      for (int i = 0; i < N_block; ++i)
        this->Yqp_Onm0_Onm1_cache_blocked[i].emplace_back(std::make_shared<fourier2DConvolutionINPUT>(Yqp_blocked[i], true, false), b);
#endif
      //$ ===================================== */
      i++;
    }

    this->m2l_cache.clear();
    this->m2l_cache.reserve(uo_m2l_cache.size());
    for (auto &[MMjk_local, Vec_AAAY_MM] : uo_m2l_cache) {
      std::sort(Vec_AAAY_MM.begin(), Vec_AAAY_MM.end(), [](const auto &lhs, const auto &rhs) { return std::abs(std::get<0>(lhs)) < std::abs(std::get<0>(rhs)); });
      this->m2l_cache.emplace_back(MMjk_local, std::move(Vec_AAAY_MM));
    }
  }

  /* --------------------------------------------------------------------------- */
#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
  Fourier2DConvolution<INPUT, double, N + 1, 2 * N + 1, N + 1, 2 * N + 1> accum_Onm0Yqp, accum_Onm1Yqp;
#elif defined(FourierM2L_BlockDecomposition)
  Fourier2DConvolution<INPUT, double, N + 1, 2 * N + 1, N + 1, 2 * N + 1> accum_Onm0Yqp, accum_Onm1Yqp;
  std::array<fourier2DConvolution, N_block> Onm0_DFT_M2L_Blocked, Onm1_DFT_M2L_Blocked;
#endif

  //% -------------------------------------------------------------------------- */
  //% -------------------------------------------------------------------------- */
  //% -------------------------------------------------------------------------- */

  void update_m2l_cache() {
#if defined(FourierM2L_double) || defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition)
    INPUT2D Onm0, Onm1;
    std::complex<TREATED> A;
    for (int n = 0; n <= N; ++n)
      for (int m = -n; m <= n; ++m) {
        A = static_cast<std::complex<TREATED>>(i_absk_A_FMM[n][m + N_i_absk_A_FMM]) * static_cast<TREATED>(1 - ((n & 1) << 1));
        auto mm = this->MM_[this->index(n, m)];
        Onm0[n][m + N] = static_cast<INPUT>(A * static_cast<std::complex<TREATED>>(std::get<0>(mm)));
        Onm1[n][m + N] = static_cast<INPUT>(A * static_cast<std::complex<TREATED>>(std::get<1>(mm)));
      }
#endif
#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
    this->Onm0_DFT_M2L.clear();
    this->Onm1_DFT_M2L.clear();
    this->Onm0_DFT_M2L.add(Onm0);
    this->Onm1_DFT_M2L.add(Onm1);
#elif defined(FourierM2L_BlockDecomposition)
    /* ------------------------ distribute into 4 blocks ------------------------ */
    std::array<INPUT2D, N_block> Onm0Blocked, Onm1Blocked;
    for (int n = 0; n <= N; ++n) {
      int index = getBlockIndex(n);
      for (int m = -n; m <= n; ++m) {
        Onm0Blocked[index][n][m + N] = Onm0[n][m + N];
        Onm1Blocked[index][n][m + N] = Onm1[n][m + N];
      }
    }
    for (int i = 0; i < N_block; ++i) {
      this->Onm0_DFT_M2L_Blocked[i].clear();
      this->Onm1_DFT_M2L_Blocked[i].clear();
      this->Onm0_DFT_M2L_Blocked[i].add(Onm0Blocked[i]);
      this->Onm1_DFT_M2L_Blocked[i].add(Onm1Blocked[i]);
    }
    /* --------------------------------------------------------------------------- */
#endif
  }

  //% -------------------------------------------------------------------------- */
  //% -------------------------------------------------------------------------- */
  //% -------------------------------------------------------------------------- */

  void m2l() {
#ifdef SimpleM2L
    std::complex<double> sum0 = 0.0, sum1 = 0.0;
    for (const auto &[MMjk_local, Vec_AAAY_MM] : this->m2l_cache) {
      sum1 = sum0 = 0.0;
      for (const auto &[AAAY, MM] : Vec_AAAY_MM) {
        sum0 += AAAY * std::get<0>(*MM);
        sum1 += AAAY * std::get<1>(*MM);
      }
      std::get<0>(*MMjk_local) += sum0;
      std::get<1>(*MMjk_local) += sum1;
    }

#elif defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition) || defined(FourierM2L_double)
    //$ ===================================== */
    //  clear this->MM_
    this->MM_.fill({});

    // \label{DFT2D_Onm0_Onm1_Yqp}
    accum_Onm0Yqp.clear();
    accum_Onm1Yqp.clear();

    /* -------------------------------------------------------------------------- */
#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
    for (const auto &[Ypq_ptr, b] : this->Yqp_Onm0_Onm1_cache) {
      accum_Onm0Yqp.Hadamard_product_add(Ypq_ptr->MatrixDFT, b->MomentsMultipoleExpansion.Onm0_DFT_M2L.MatrixDFT);
      accum_Onm1Yqp.Hadamard_product_add(Ypq_ptr->MatrixDFT, b->MomentsMultipoleExpansion.Onm1_DFT_M2L.MatrixDFT);
    }

    if constexpr (std::same_as<INPUT, std::complex<DoubleDouble>>) {
      accum_Onm0Yqp.copyDoubleDoubleToDouble();
      accum_Onm1Yqp.copyDoubleDoubleToDouble();
    }

    accum_Onm0Yqp.convolve();
    accum_Onm1Yqp.convolve();
    /* ------------------------------- block版あくまでテスト ----------------------- */
#elif defined(FourierM2L_BlockDecomposition)
    for (int k = 0; k < N_block; ++k)
      for (int i = 0; i < N_block; ++i)
        for (const auto &[Ypq_ptr, b] : this->Yqp_Onm0_Onm1_cache_blocked[i])
          for (int j = 0; j < N_block; ++j)
            if (k == i + j) {
              accum_Onm0Yqp.Hadamard_product_add(Ypq_ptr->MatrixDFT, b->MomentsMultipoleExpansion.Onm0_DFT_M2L_Blocked[j].MatrixDFT);
              accum_Onm1Yqp.Hadamard_product_add(Ypq_ptr->MatrixDFT, b->MomentsMultipoleExpansion.Onm1_DFT_M2L_Blocked[j].MatrixDFT);
            }

    accum_Onm0Yqp.convolve();
    accum_Onm1Yqp.convolve();
#endif
    /* -------------------------------------------------------------------------- */
    std::complex<double> A; //! ここはdoubleでも良いようだ．

    for (int j = N; j >= 0; --j) {
      for (int k = -j; k <= j; ++k) {
        A = i_absk_A_FMM[j][k + N_i_absk_A_FMM];
        auto &MM = this->MM_[this->index(j, k)];
        std::get<0>(MM) += static_cast<cd>(A * accum_Onm0Yqp.convolution_T[k + 2 * N][N - j]);
        std::get<1>(MM) += static_cast<cd>(A * accum_Onm1Yqp.convolution_T[k + 2 * N][N - j]);
      }
    }

    //$ ===================================== */
#endif
  }

  //^ ----------------------------------- L2L ---------------------------------- */

  std::vector<Type_IICCIIIC> index_map_for_L2L = index_map_for_L2L_(); // i,j,index, n, m, index

  constexpr std::vector<Type_IICCIIIC> index_map_for_L2L_() {
    std::vector<Type_IICCIIIC> index_map_for_L2L;
    index_map_for_L2L.reserve((N + 1) * (2 * N + 1) * (N + 1) * (2 * N + 1));
    int j, k, n, m;
    //   for (j = 0; j <= N; ++j)
    for (j = N; j >= 0; --j)
      for (k = -j; k <= j; ++k) {
        auto &AAA_L2L_FMM_j_k = AAA_L2L_FMM[j][k + N_AAA_L2L_FMM];
        // for (n = j; n <= N; ++n)
        for (n = N; n >= j; --n)
          for (m = -n; m <= n; ++m) {
            auto AAA = AAA_L2L_FMM_j_k[n][m + N_AAA_L2L_FMM];
            if (AAA.real() != 0.0 || AAA.imag() != 0.0)
              index_map_for_L2L.emplace_back(j, k, &this->MM_[this->index(j, k)], n, m, this->index(n, m), AAA);
          }
      }
    return index_map_for_L2L;
  }

  std::tuple<std::unique_ptr<SphericalCoordinates>, Moments *> l2l_cache;

  //! これは計算毎に変わる
  void set_l2l(Moments *MomentsLocalExpansion) {
    //   auto sph = std::make_unique<SphericalCoordinates>(MomentsLocalExpansion->X - this->X);
    auto sph = std::make_unique<SphericalCoordinates>((MomentsLocalExpansion->X - this->X));
    sph->precompute_sph(2 * N);
    this->l2l_cache = {std::move(sph), MomentsLocalExpansion};
  }

  void l2l() {
    cmplx R;
    for (const auto &[j, k, MMjk, n, m, index, AAA] : this->index_map_for_L2L) {
      // R = AAA * std::get<0>(l2l_cache)->sph_harmonics_rho(n - j, m - k);
      R = std::get<0>(l2l_cache)->l2lFunction(j, k, n, m);
      std::get<0>(*MMjk) += R * std::get<0>(std::get<1>(l2l_cache)->MM_[index]);
      std::get<1>(*MMjk) += R * std::get<1>(std::get<1>(l2l_cache)->MM_[index]);
    }
  }
};

//! ========================================================================== */
//! ========================================================================== */
//! ========================================================================== */

#include "lib_spatial_partitioning.hpp"

//! ========================================================================== */
//! ========================================================================== */
//! ========================================================================== */

//% -------------------------------------------------------------------------- */
//%        Functions Related to Multipole Expansion Utilizing Buckets          */
//% -------------------------------------------------------------------------- */
// #define _DEBUG_FMM_
void setM2M(auto &B_poles) {
  TimeWatch tw, tw_all;
  for (int level = B_poles.max_level - 1; level >= 0; level--) {
    B_poles.forEachAtLevel({level}, [&](auto *B) { B->MomentsMultipoleExpansion.set_m2m(B->getAllChildren()); });
#if defined(_DEBUG_FMM_)
    std::cout << yellow << "setM2M" << ", level=" << level << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
#endif
  }
  std::cout << Yellow << "setM2M" << Green << ", Elapsed time : " << tw_all() << colorReset << std::endl;
}

void M2M(auto &B_poles) {
  for (int level = B_poles.max_level - 1; level >= 0; level--)
    B_poles.forEachAtLevel({level}, [&](auto *B) { B->MomentsMultipoleExpansion.m2m(); });
};

/* -------------------------------------------------------------------------- */

void setL2L(auto &B_poles) {
  for (auto &buckets_from_top_level : B_poles.level_buckets)
    for (auto &parent : buckets_from_top_level)
      parent->traverseChildren([&](auto *child) { child->MomentsLocalExpansion.set_l2l(&parent->MomentsLocalExpansion); });
}

void L2L(auto &B_poles) {
  for (auto &buckets_from_top_level : B_poles.level_buckets)
    for (auto &parent : buckets_from_top_level)
      parent->traverseChildren([&](auto *child) { child->MomentsLocalExpansion.l2l(); });
}

/* -------------------------------------------------------------------------- */

const double scale = 3.;
/*

scaledBoundsが返すのは，
xmin <-- L -->xmax　として，
xmin <-- L/2 -- xcenter -- L/2 -->xmax
xmin <-- scale*L/2 -- xcenter -- scale*L/2 -->xmax
xmin <-- scale*L -->xmax と結果的になる．
*/
bool isFar(auto *A, auto *B) { return !InsideQ(B->X, A->scaledBounds(scale)); };
bool isNear(auto *A, auto *B) { return InsideQ(B->X, A->scaledBounds(scale)); };
bool isFar(const auto *A, const auto *B);

/*
M2Lの対象の条件
1. 同じレベル
2. 十分離れているisFar==true, また直接近傍ではないisNear==false
3. 親同士はisNear==true
*/

template <typename T> void checkAndAddBuckets(T A, T B) {
  // isNear(A,B) && isFar(A->c,B->c)となるA->cとB->cを探し，A->c->buckets_for_M2LにB->cを追加する．
  if (isNear(A, B)) {
    if (A != B)
      A->buckets_near.emplace_back(B);
    A->traverseChildren([&](auto *A_c) {
      B->traverseChildren([&](auto *B_c) {
        if (A_c != B_c && isFar(A_c, B_c))
          A_c->buckets_for_M2L.emplace_back(B_c);
        else
          checkAndAddBuckets(A_c, B_c);
      });
    });
  }
}

// 極が少ないバケツは，M2Lの時に省略する．

void setBucketsForM2L(auto &B_poles) {

  B_poles.traverseTree([&](auto *A) {
    A->buckets_for_M2L.clear();
    A->buckets_near.clear();
    A->buckets_for_L2M.clear();
  });

  B_poles.traverseChildren([&](auto *A) { B_poles.traverseChildren([&](auto *B) { checkAndAddBuckets(A, B); }); });

  for (auto &buckets_from_top_level : B_poles.level_buckets) {
    for (auto &A : buckets_from_top_level)
      for (auto &B : A->buckets_for_M2L)
        B->buckets_for_L2M.emplace_back(A);
  }
}

void setM2L(auto &B_poles) {
  TimeWatch tw, tw_all;
  // std::cout << "各レベルの各セルのM2Lの相手を保存する" << std::endl;
  setBucketsForM2L(B_poles);
  // A -> M2L -> B
  int level = 0;
  for (auto &buckets_at_a_level : B_poles.level_buckets) {
#pragma omp parallel for
    for (auto &A : buckets_at_a_level)
      A->MomentsLocalExpansion.set_m2l(A->buckets_for_L2M);
    std::cout << yellow << "setM2L" << ", level=" << level << ", buckets_at_a_level.size()=" << buckets_at_a_level.size() << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
    level++;
  }
  std::cout << Yellow << "setM2L" << Green << ", Elapsed time : " << tw_all() << colorReset << std::endl;
};

void M2L(auto &B_poles) {
  for (auto &buckets_at_a_level : B_poles.level_buckets) {
#pragma omp parallel for
    for (auto &A : buckets_at_a_level)
      A->MomentsMultipoleExpansion.update_m2l_cache();
  }

  for (auto &buckets_at_a_level : B_poles.level_buckets) {
#pragma omp parallel for
    for (auto &A : buckets_at_a_level)
      A->MomentsLocalExpansion.m2l();
  }
};

/* -------------------------------------------------------------------------- */

/*
source4FMMは，
Tddd X;
Tdd weights;
Tddd normal;
std::function<Tdd()> getValues;
を持つ．これらを使うと，直接数値面積分ができる．
*/

/*

BIE：alpha_phi = w_Gn_phin - w_Gn_phiの，w_Gn_phin，w_Gn_phiをまず計算する．

ここで，w_Gn_phi，w_Gn_phiはそれぞれ，

$$
\begin{align*}
\int\int_{\Gamma} G \nabla \phi \cdot {\mathbfit n} dS\\
\int\int_{\Gamma} \phi \nabla G \cdot {\mathbfit n} dS
\end{align*}
$$

w_Gn_phiのincrementにおいて，符号が負になっている．
これは，勘違いして修正してしまいそうだが，w_Gn_phi自体が負になっているためで，
BIEの符号と混同しないように注意する．

BIEは次の形である．

alpha_phi = w_Gn_phin - w_Gn_phi

多重極展開（see L2P）も，w_Gn_phin，w_Gn_phiの近似を計算しており，直接積分しているものと対応している．

*/

/* -------------------------------------------------------------------------- */

void initializeFMM(auto &B_poles, const auto &targets) {

  std::cout << "initializeFMM" << std::endl;

// M2Lの方法を表示
#if defined(SimpleM2L)
  std::cout << "SimpleM2L" << std::endl;
#elif defined(FourierM2L_double)
  std::cout << "FourierM2L_double" << std::endl;
#elif defined(FourierM2L_DoubleDouble)
  std::cout << "FourierM2L_DoubleDouble" << std::endl;
#elif defined(FourierM2L_BlockDecomposition)
  std::cout << "FourierM2L_BlockDecomposition" << std::endl;
#endif

  std::cout << "[FMM:init] Step0: update poles start, #poles=" << B_poles.data1D_vector.size() << std::endl;

  // ------------------------------------------------------------------
  // Step 0: 事前準備
  //  - data1D_vector に入っている全 source の値を最新化
  //  - ここでの source->update() は，各 source が保持する (phi, phin) 等の動的量を再計算し，多重極係数生成に必要な weighted_source_densities を一貫させる目的。
  // ------------------------------------------------------------------
  for (auto source : B_poles.data1D_vector)
    source->updateDensity();
  std::cout << "[FMM:init] Step0: update poles done" << std::endl;

  // ------------------------------------------------------------------
  // Step 1: 全ノードのモーメントをゼロ初期化
  //  - Multipole(Local) Moments の係数配列 MM_ をクリア
  //  - FMM の各フェーズ (M2M, M2L, L2L) を再構築する前に必須
  // ------------------------------------------------------------------
  size_t total_nodes_before_reset = 0;
  for (auto &lv : B_poles.level_buckets)
    total_nodes_before_reset += lv.size();
  std::cout << "[FMM:init] Step1: reset moments for total nodes=" << total_nodes_before_reset << ", levels=" << B_poles.level_buckets.size() << std::endl;
  B_poles.forEachAll([&](auto *b) {
    b->MomentsMultipoleExpansion.initialize(); //! MとM_の初期化
    b->MomentsLocalExpansion.initialize();     //! MとM_の初期化
  });
  std::cout << "[FMM:init] Step1: reset moments done" << std::endl;

  // ------------------------------------------------------------------
  // Step 2: 葉ノード(最深バケツ)で P2M (increment_M) を実行
  //  - 各葉が保持するソース点集合 data1D_vector から多重極係数を構築
  //  - increment_M は再利用キャッシュ(vector_weightedSourceDensities_RRn_MM_) を再生成
  //  - 並列 (par_unseq/OMP) で高速化
  // ------------------------------------------------------------------
  size_t n_leaves = B_poles.deepest_level_buckets.size();
  size_t sources_in_leaves = 0;
  for (auto *leaf : B_poles.deepest_level_buckets)
    sources_in_leaves += leaf->data1D_vector.size();
  std::cout << "[FMM:init] Step2: P2M begin, #leaves=" << n_leaves << ", total leaf sources=" << sources_in_leaves << std::endl;
  B_poles.forEachAtDeepestParallel([&](auto *b) {
    b->MomentsMultipoleExpansion.increment_M(b->data1D_vector); // !vector_weightedSourceDensities_RRn_MMの初期化
  });
  std::cout << "[FMM:init] Step2: P2M (increment_M) done" << std::endl;

  // ------------------------------------------------------------------
  // Step 3: M2M 準備 (setM2M) → 親方向への集約で必要な幾何キャッシュ生成
  //         M2L 準備 (setM2L) → 相互作用(遠方)リストを構築し M2L キャッシュ生成
  // ------------------------------------------------------------------
  std::cout << "[FMM:init] Step3: setM2M start" << std::endl;
  setM2M(B_poles);
  std::cout << "[FMM:init] Step3: setM2M done, setM2L start" << std::endl;
  setM2L(B_poles);
  std::cout << "[FMM:init] Step3: setM2L done" << std::endl;

  // ------------------------------------------------------------------
  // Debug: 各レベルの M2L 接続数統計
  //  - level ごとの bucket 数とその平均 M2L 相手数を出力し，分割条件や遠方判定 (isFar / isNear) が期待どおりかを確認する。
  // ------------------------------------------------------------------
  //! check number of couping
  int level = 0;
  for (auto &buckets_at_a_level : B_poles.level_buckets) {
    std::cout << "level=" << level++ << std::endl;
    if (buckets_at_a_level.size() == 0)
      continue;
    double size = 0.;
    for (auto &A : buckets_at_a_level)
      size += A->buckets_for_M2L.size();
    std::cout << "mean A->buckets_for_M2L=" << size / buckets_at_a_level.size() << std::endl;
  }
  std::cout << "[FMM:init] Debug: coupling stats printed" << std::endl;

  // ------------------------------------------------------------------
  // Step 4: L2L 準備 (setL2L)
  //  - 親 Local Moment を子へ伝搬するための幾何キャッシュ (spherical harmonics)を設定
  // ------------------------------------------------------------------
  std::cout << "[FMM:init] Step4: setL2L start" << std::endl;
  setL2L(B_poles);
  std::cout << "[FMM:init] Step4: setL2L done" << std::endl;

  // ------------------------------------------------------------------
  // Step 5: L2P 用ターゲット前処理
  //  - 各ターゲット点 t に対し，必要な (Y_nm ρ^{n+1}) * M_local のペアを前計算
  //  - integrateFMM() が O(N^2) 反復で足し合わせるだけになる
  // ------------------------------------------------------------------
  TimeWatch tw;
  std::cout << "各レベルの各セルのL2Pの相手を保存する (#targets=" << targets.size() << ")" << std::endl;
#pragma omp parallel for
  for (auto &t : targets)
    t->setL2P(B_poles); // !L2Pの初期化
  std::cout << Yellow << "setL2P" << Green << ", Elapsed time : " << tw() << colorReset << " [FMM:init Step5 done]" << std::endl;

  // ------------------------------------------------------------------
  // Step 6: 直接積分 (near field) 用前処理
  //  - 近傍バケツ (buckets_near) とそのソースを列挙し，ターゲットごとにソース寄与をキャッシュ (vec_phiphin_WGNWGnN)
  //  - GMRES の反復で再構築不要としコスト削減
  // ------------------------------------------------------------------
  std::cout << "各ターゲットの近傍ソースを保存する (#targets=" << targets.size() << ")" << std::endl;
#pragma omp parallel for
  for (auto &t : targets)
    t->setDirectIntegration(B_poles); // !Direct integrationの初期化
  std::cout << Yellow << "setDirectIntegration" << Green << ", Elapsed time : " << tw() << colorReset << " [FMM:init Step6 done]" << std::endl;

  // ------------------------------------------------------------------
  // Step 7: 近傍ソースキャッシュの平均サイズを出力 (疎密バランスチェック)
  //  - 過剰に大きい場合は分割閾値 (grow_condition) の調整を検討
  // ------------------------------------------------------------------
  double size = 0.;
  for (auto &t : targets) {
    if (t->vec_phiphin_WGNWGnN.size() == 0)
      continue;
    size = size + t->vec_phiphin_WGNWGnN.size();
  }
  std::cout << "mean vec_phiphin_WGNWGnN=" << size / targets.size() << " [FMM:init Step7 done]" << std::endl;
};

/* -------------------------------------------------------------------------- */

void updateFMM(auto &B_poles, std::array<double, 6> &elapsed_time) {
  TimeWatch tw;
  for (auto source : B_poles.data1D_vector)
    source->updateDensity();

  elapsed_time[0] = tw()[0];

  // !この操作は，Momentsのvector_weightedSourceDensities_RRn_MM_を削除しない．
  // !vector_weightedSourceDensities_RRn_MM_の保存するポインタ変数を初期化するだけで，配列内のポインタは生きたまま残る
  // !increment_Mを使うと，vector_weightedSourceDensities_RRn_MM_は初期化され再計算される
  B_poles.forEachAllParallel([&](auto *b) {
    b->MomentsMultipoleExpansion.initialize();
    b->MomentsLocalExpansion.initialize();
  });
  elapsed_time[1] = tw()[0];

  B_poles.forEachAtDeepestParallel([&](auto *b) { b->MomentsMultipoleExpansion.increment_M_reuse(); });

  elapsed_time[2] = tw()[0];
  M2M(B_poles);

  elapsed_time[3] = tw()[0];
  M2L(B_poles);

  elapsed_time[4] = tw()[0];
  L2L(B_poles);

  elapsed_time[5] = tw()[0];
};

void updateFMM(auto &B_poles) {
  TimeWatch tw;
  for (auto source : B_poles.data1D_vector)
    source->updateDensity();

  std::cout << Yellow << "update poles" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

  // !この操作は，Momentsのvector_weightedSourceDensities_RRn_MM_を削除しない．
  // !vector_weightedSourceDensities_RRn_MM_の保存するポインタ変数を初期化するだけで，配列内のポインタは生きたまま残る
  // !increment_Mを使うと，vector_weightedSourceDensities_RRn_MM_は初期化され再計算される
  B_poles.forEachAllParallel([&](auto *b) {
    b->MomentsMultipoleExpansion.initialize();
    b->MomentsLocalExpansion.initialize();
  });

  std::cout << Yellow << "reset Moments" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

  B_poles.forEachAtDeepestParallel([&](auto *b) { b->MomentsMultipoleExpansion.increment_M_reuse(); });
  std::cout << Yellow << "increment_M_reuse" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

  M2M(B_poles);
  std::cout << Magenta << "M2M" << ", Elapsed time : " << tw() << colorReset << std::endl;

  M2L(B_poles);
  std::cout << Magenta << "M2L" << ", Elapsed time : " << tw() << colorReset << std::endl;

  L2L(B_poles);
  std::cout << Magenta << "L2L" << ", Elapsed time : " << tw() << colorReset << std::endl;

  std::cout << Yellow << "updateFMM" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
};

/* -------------------------------------------------------------------------- */

void write2Dcsv(const std::string &filename, const std::vector<std::vector<double>> &data) {
  std::ofstream ofs(filename);
  for (const auto &row : data) {
    for (size_t j = 0; j < row.size(); ++j) {
      ofs << std::setprecision(15) << row[j];
      if (j < row.size() - 1)
        ofs << ",";
    }
    ofs << std::endl;
  }
  std::cout << "write2Dcsv: " << filename << ", size = {" << data.size() << ", " << data[0].size() << "}" << std::endl;
  ofs.close();
};

void write2Dcsv(const std::string &filename, const std::vector<std::vector<std::complex<double>>> &data) {
  std::ofstream ofs(filename);
  for (const auto &row : data) {
    for (size_t j = 0; j < row.size(); ++j) {
      ofs << std::setprecision(15) << row[j].real();
      if (j < row.size() - 1)
        ofs << ",";
    }
    ofs << std::endl;
  }
  std::cout << "write2Dcsv: " << filename << ", size = {" << data.size() << ", " << data[0].size() << "}" << std::endl;
  ofs.close();
};

std::vector<std::vector<double>> Re(const auto &data) {
  std::vector<std::vector<double>> result(data.size(), std::vector<double>(data[0].size()));
  for (size_t i = 0; i < data.size(); ++i)
    for (size_t j = 0; j < data[0].size(); ++j)
      result[i][j] = data[i][j].real();
  return result;
}
std::vector<std::vector<double>> Im(const auto &data) {
  std::vector<std::vector<double>> result(data.size(), std::vector<double>(data[0].size()));
  for (size_t i = 0; i < data.size(); ++i)
    for (size_t j = 0; j < data[0].size(); ++j)
      result[i][j] = data[i][j].imag();
  return result;
}
