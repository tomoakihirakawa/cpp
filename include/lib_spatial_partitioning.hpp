#ifndef lib_spatial_partitioning_H
#define lib_spatial_partitioning_H

template <typename T>
struct BaseBuckets {
   using sizeType = int;
   using ST = sizeType;
   using ST2 = std::array<sizeType, 2>;
   using ST3 = std::array<sizeType, 3>;
   using ST6 = std::array<sizeType, 6>;
   //! getX()でxyz座標を取得できるオブジェクトTのための，バケツ
   Tdd xbounds, ybounds, zbounds;
   ST xsize, ysize, zsize;
   T3Tdd bounds;
   Tddd center;
   ST3 dn;
   std::vector<std::vector<std::vector<std::unordered_set<T>>>> buckets;
   std::unordered_set<T> all_stored_objects;
   std::unordered_map<T, ST3> map_to_ijk;
   double dL;
   double bucketVolume() const { return std::pow(this->dL, 3.); };
   BaseBuckets(const CoordinateBounds &c_bounds, const double dL_IN) : bounds(c_bounds.bounds), dL(dL_IN), dn(ST3{0, 0, 0}) { initialize(this->bounds, dL_IN); };
   BaseBuckets(const T3Tdd &boundingboxIN, const double dL_IN) : bounds(boundingboxIN), dL(dL_IN), dn(ST3{0, 0, 0}) { initialize(this->bounds, dL_IN); };
   void initialize(const T3Tdd &boundingboxIN, const double dL_IN) {
      this->buckets.clear();
      this->all_stored_objects.clear();
      //
      this->bounds = boundingboxIN;
      this->dL = dL_IN;
      this->xbounds = std::get<0>(this->bounds);
      this->ybounds = std::get<1>(this->bounds);
      this->zbounds = std::get<2>(this->bounds);
      this->center = compute_center(this->xbounds, this->ybounds, this->zbounds);
      this->xsize = std::ceil((std::get<1>(this->xbounds) - std::get<0>(this->xbounds)) / this->dL);
      this->ysize = std::ceil((std::get<1>(this->ybounds) - std::get<0>(this->ybounds)) / this->dL);
      this->zsize = std::ceil((std::get<1>(this->zbounds) - std::get<0>(this->zbounds)) / this->dL);
      this->dn = {xsize, ysize, zsize};
      this->buckets.resize(xsize, std::vector</*y*/ std::vector</*z*/ std::unordered_set<T>>>(ysize, std::vector</*z*/ std::unordered_set<T>>(zsize, std::unordered_set<T>{})));

      std::cout << "１バケット幅 = " << this->dL << std::endl;
      std::cout << "バウンド = " << this->bounds << std::endl;
      std::cout << "サイズ = " << this->dn << std::endl;
   }
   const std::unordered_set<T> &getAll() const { return this->all_stored_objects; };
   constexpr Tddd compute_center(const Tdd &xbounds, const Tdd &ybounds, const Tdd &zbounds) {
      return {(std::get<1>(xbounds) + std::get<0>(xbounds)) / 2,
              (std::get<1>(ybounds) + std::get<0>(ybounds)) / 2,
              (std::get<1>(zbounds) + std::get<0>(zbounds)) / 2};
   }
   //@ -------------------------------- インデックス変換 -------------------------------- */
   // x座標を内包するバケツのインデックスを返す
   constexpr Tddd itox(const ST i, const ST j, const ST k) const {
      return {this->dL * 0.5 + this->dL * i + std::get<0>(this->xbounds),
              this->dL * 0.5 + this->dL * j + std::get<0>(this->ybounds),
              this->dL * 0.5 + this->dL * k + std::get<0>(this->zbounds)};
   };
   constexpr Tddd itox(const ST3 &ijk) const { return itox(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk)); };
   constexpr ST3 indices_no_clamp(const Tddd &x) const {
      return {static_cast<ST>(std::floor((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL)),
              static_cast<ST>(std::floor((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL)),
              static_cast<ST>(std::floor((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL))};
   };
   constexpr ST3 indices(const Tddd &x) const {
      /*
      intキャストはゼロ方向へ実数を切り捨てた結果を返すので，static_cast<int>によって正しくセルのインデックスに変換できる．
       0.**    1.**  2.**   3.**
      <-dL-> <-dL-> <-dL-> <-dL->
      *-----*------*------*------*
      |  0  |   1  |   2  |   3  |
      *-----*------*------*------*
      */
      const auto ijk = indices_no_clamp(x);
      return {std::clamp(std::get<0>(ijk), static_cast<ST>(0), this->xsize - 1),
              std::clamp(std::get<1>(ijk), static_cast<ST>(0), this->ysize - 1),
              std::clamp(std::get<2>(ijk), static_cast<ST>(0), this->zsize - 1)};
   };
   constexpr ST6 indices_ranges(const Tddd &x, const double d) const {
      auto [i_min, j_min, k_min] = indices(x - d);
      auto [i_max, j_max, k_max] = indices(x + d);
      return {i_min, i_max, j_min, j_max, k_min, k_max};
   };
   //@ ------------------------ インデックスがboundsに収まっているかどうか ------------------------ */
   bool isInside(const ST i, const ST j, const ST k) const { return (i >= 0 && j >= 0 && k >= 0 && i < this->xsize && j < this->ysize && k < this->zsize); };
   bool isInside(const ST3 &ijk) const { return isInside(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk)); };
   bool isInside(const Tddd &x) const { return isInside(indices_no_clamp(x)); };
   void erase(T const p) {
      auto it = this->map_to_ijk.find(p);
      if (it != this->map_to_ijk.end()) {
         auto [i, j, k] = it->second;
         this->all_stored_objects.erase(p);
         this->buckets[i][j][k].erase(p);
         this->map_to_ijk.erase(it);
      }
   };
   //@ -------------------------------------------------------------------------- */
   //! explicitly specify coordinates
   bool add_force(const ST3 &ijk, const T p) {
      const auto [i, j, k] = ijk;
      this->map_to_ijk[p] = ijk;
      const bool bucket_inserted = this->buckets[i][j][k].emplace(p).second;
      const bool all_objects_inserted = this->all_stored_objects.emplace(p).second;
      return bucket_inserted && all_objects_inserted;
   };
   bool add(const Tddd &x, const T p) { return add_force(indices(x), p); };
   //! automatically specify coordinates
   bool add(const std::unordered_set<T> &P) {
      bool ret = true;
      for (const auto &p : P) {
         auto ijk = indices(ToX(p));
         this->map_to_ijk[p] = ijk;
         ret = ret && add_force(ijk, p);
      }
      return ret;
   };
   /* -------------------------------------------------------------------------- */
   auto getBucket(const Tddd &x) const {
      const auto ijk = this->indices(x);
      return this->buckets[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
   };

   auto getBucket(const Tddd &x, const double &range) const {
      const auto ijk = this->indices(x);
      std::unordered_set<T> ret;
      this->apply(x, range, [&](auto &a) { ret.emplace(a); });
      return ret;
   };

   void clear() {
      this->buckets.clear();
      this->all_stored_objects.clear();
   };
   // b@ -------------------------------------------------------------------------- */
   // b@                             STL like functions                             */
   // b@ -------------------------------------------------------------------------- */

   //! none_of
   bool none_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (this->buckets.empty()) {
         return true;
      }
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);

      return std::none_of(this->buckets.cbegin() + i_min, this->buckets.cbegin() + i_max + 1, [&](const auto &Bi) {
         return std::none_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
            return std::none_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
               return std::none_of(Bijk.cbegin(), Bijk.cend(), func);
            });
         });
      });
   }

   //! all_of
   bool all_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (this->buckets.empty()) {
         return true;
      }
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);

      return std::all_of(this->buckets.cbegin() + i_min, this->buckets.cbegin() + i_max + 1, [&](const auto &Bi) {
         return std::all_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
            return std::all_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
               return std::all_of(Bijk.cbegin(), Bijk.cend(), func);
            });
         });
      });
   }

   //! any_of
   bool any_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (this->buckets.empty()) {
         return false;
      }
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);

      return std::any_of(this->buckets.cbegin() + i_min, this->buckets.cbegin() + i_max + 1, [&](const auto &Bi) {
         return std::any_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
            return std::any_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
               return std::any_of(Bijk.cbegin(), Bijk.cend(), func);
            });
         });
      });
   }

   //! apply
   void apply(const Tddd &x, const double d, const std::function<void(const T &)> &func) const {
      if (this->buckets.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D buckets is empty");
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
      std::for_each(std::execution::unseq, this->buckets.cbegin() + i_min, this->buckets.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
         std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
            std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
               for (const auto &p : Bijk) func(p);
            });
         });
      });
   };
};

template <typename T>
struct Buckets : public BaseBuckets<T> {
   Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<T>(c_bounds, dL_IN){};
   Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<T>(boundingboxIN, dL_IN){};
};

#endif